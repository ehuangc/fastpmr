"""Gaussian process model for spatiotemporal pairwise mismatch rate (PMR).

Fits a variational GP where inputs are sample-pair coordinates
(x_i, y_i, t_i, x_j, y_j, t_j) in an equal-area projection (km) and the output is PMR.
Write s_i = (x_i, y_i, t_i) for the spacetime coordinates of sample i.

The pair kernel is symmetric under (i, j) swap and has two additive components, each
capturing a structurally different aspect of a pair:
  1. A latent per-sample function f(x, y, t) ~ GP(0, k_f) with a separable
     Matern-5/2(x, y) x Matern-5/2(time) kernel, lifted to pairs as
     K_pair((s_i1, s_j1), (s_i2, s_j2)) = sum over {a in pair1, b in pair2} k_f(a, b),
     which by bilinearity of covariance equals Cov(f(s_i1) + f(s_j1), f(s_i2) + f(s_j2)).
     This captures sample-level structure: f's contribution to a pair depends on
     *where* each sample is (in absolute coords and time), not on the within-pair gap.
  2. A within-pair function g(d_space, d_time) ~ GP(0, k_g) with a separable
     Matern-5/2(d_space) x Matern-5/2(d_time) kernel, where
     d_space = ||(x_i, y_i) - (x_j, y_j)|| and d_time = |t_i - t_j|.
     This captures the direct dependence of PMR on *how separated* two samples are in
     space and time.

As such, the model is
    r_ij = c + f(s_i) + f(s_j) + g(d_space_ij, d_time_ij) + epsilon,
where c is a constant mean and epsilon is i.i.d. Gaussian likelihood noise.

Predicting the "latent relatedness" l_ii for coincident samples s_i is then
    l_ii = c + 2 * f(s_i) + g(0, 0).
"""

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import gpytorch
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import torch

# gpytorch.variational.VariationalStrategy._cholesky_factor unconditionally upcasts
# to gpytorch.settings._linalg_dtype_cholesky.value() (default float64). On MPS,
# float64 is unsupported, so we must force float32 via this (private) setting. There
# is no public alternative in this gpytorch version.
from gpytorch.settings import _linalg_dtype_cholesky
from torch.utils.data import DataLoader, TensorDataset

from evaluation_utils import (
    AADR_DIR,
    DATE_MEAN_BP_FIELD,
    INDIVIDUAL_ID_FIELD,
    LAT_FIELD,
    LON_FIELD,
    ensure_aadr_npz_present,
    is_archaic_or_reference_sample,
    load_aadr_metadata,
    load_aadr_npz_arrays,
)

OUTPUT_DIR = AADR_DIR / "results" / "gp"
CHECKPOINT_PATH = OUTPUT_DIR / "gp_model.pt"
MAP_PATH = OUTPUT_DIR / "gp_pmr_map.pdf"

OVERLAP_THRESHOLD = 30_000
MAX_PAIRS = 1_000_000
LOW_PMR_Z_CUT = 5.0
SEED = 42

N_INDUCING = 512
N_EPOCHS = 20
BATCH_SIZE = 4096
LEARNING_RATE = 0.01

# Equal-area projection used for both training inputs and map plotting
PROJECTION = ccrs.LambertAzimuthalEqualArea(central_longitude=15, central_latitude=50)

TIME_SLICES_BP = [9_000, 6_000, 3_000, 1_000]
GRID_RESOLUTION_DEG = 1.0

# Lat/lon bounding box used for sample filtering, as (lon_min, lon_max, lat_min, lat_max).
# This is broader than Europe proper: samples in adjacent regions (e.g., North Africa, Anatolia,
# the Caucasus) that fall inside the box are included to provide boundary support
MAP_EXTENT_DEG = (-15.0, 45.0, 34.0, 72.0)


def project_coords(lats: np.ndarray, lons: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Project geographic lat/lon to equal-area (x, y) in km."""
    pts = PROJECTION.transform_points(ccrs.PlateCarree(), np.asarray(lons), np.asarray(lats))
    return pts[..., 0] / 1000.0, pts[..., 1] / 1000.0


def lat_lon_box_path(extent_deg: tuple[float, float, float, float], n: int = 200) -> mpath.Path:
    """Project the boundary of a lat/lon rectangle into PROJECTION coords as a closed Path.

    The rectangle's edges become curves under LAEA, so we densely sample each edge
    before projecting. The resulting path is used as the axes boundary so the
    displayed region matches the lat/lon prediction region exactly.
    """
    lon_min, lon_max, lat_min, lat_max = extent_deg
    lons = np.concatenate(
        [
            np.linspace(lon_min, lon_max, n),
            np.full(n, lon_max),
            np.linspace(lon_max, lon_min, n),
            np.full(n, lon_min),
        ]
    )
    lats = np.concatenate(
        [
            np.full(n, lat_min),
            np.linspace(lat_min, lat_max, n),
            np.full(n, lat_max),
            np.linspace(lat_max, lat_min, n),
        ]
    )
    pts = PROJECTION.transform_points(ccrs.PlateCarree(), lons, lats)
    return mpath.Path(pts[:, :2])


def get_device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def extract_pair_data(
    overlap_threshold: int = OVERLAP_THRESHOLD,
    max_pairs: int = MAX_PAIRS,
    low_pmr_z_cut: float = LOW_PMR_Z_CUT,
    seed: int = SEED,
) -> tuple[np.ndarray, np.ndarray]:
    """Load AADR data and extract pairwise observations for GP fitting.

    Returns (pair_coords, pmr), where pair_coords has columns
    (x_km_i, y_km_i, t_i, x_km_j, y_km_j, t_j) in the equal-area projection.
    """
    ensure_aadr_npz_present()
    samples, site_overlaps, mismatch_rates, _, _, covered_snps = load_aadr_npz_arrays()
    metadata = load_aadr_metadata()

    # Collect samples with valid coordinates, dates, and Individual IDs
    candidates: list[tuple[int, float, float, float, str, int]] = []
    for sample_idx, sample in enumerate(samples):
        sample_metadata = metadata.get(sample)
        if sample_metadata is None or is_archaic_or_reference_sample(sample, sample_metadata):
            continue
        try:
            lat = float(sample_metadata[LAT_FIELD])
            lon = float(sample_metadata[LON_FIELD])
            date_bp = float(sample_metadata[DATE_MEAN_BP_FIELD])
        except (KeyError, TypeError, ValueError):
            continue
        if not (np.isfinite(lat) and np.isfinite(lon) and np.isfinite(date_bp) and date_bp >= 0):
            continue
        individual_id = sample_metadata.get(INDIVIDUAL_ID_FIELD, "").strip()
        candidates.append((sample_idx, lat, lon, date_bp, individual_id, int(covered_snps[sample_idx])))

    # Deduplicate by Individual ID, keeping the sample with the most covered SNPs;
    # samples with no / ".." ID are kept as-is since they can't be matched
    best_pos_by_individual: dict[str, int] = {}
    unkeyed_positions: list[int] = []
    for pos, (_, _, _, _, individual_id, n_snps) in enumerate(candidates):
        if not individual_id or individual_id == "..":
            unkeyed_positions.append(pos)
            continue
        prev = best_pos_by_individual.get(individual_id)
        if prev is None or n_snps > candidates[prev][5]:
            best_pos_by_individual[individual_id] = pos
    selected_positions = sorted(unkeyed_positions + list(best_pos_by_individual.values()))

    original_indices = np.asarray([candidates[p][0] for p in selected_positions], dtype=int)
    lats = np.asarray([candidates[p][1] for p in selected_positions], dtype=float)
    lons = np.asarray([candidates[p][2] for p in selected_positions], dtype=float)
    dates = np.asarray([candidates[p][3] for p in selected_positions], dtype=float)

    # Restrict to the fixed MAP_EXTENT_DEG bounding box. This includes Europe proper
    # plus adjacent samples (North Africa, Anatolia, Caucasus) that happen to fall
    # inside the box, which provide boundary support for the GP
    lon_min, lon_max, lat_min, lat_max = MAP_EXTENT_DEG
    box_mask = (lats >= lat_min) & (lats <= lat_max) & (lons >= lon_min) & (lons <= lon_max)
    original_indices = original_indices[box_mask]
    lats = lats[box_mask]
    lons = lons[box_mask]
    dates = dates[box_mask]
    n_samples = len(original_indices)
    print(f"Samples inside lat/lon box after filtering: {n_samples}.")

    triu_i, triu_j = np.triu_indices(n_samples, k=1)
    row_idx = original_indices[triu_i]
    col_idx = original_indices[triu_j]
    pair_overlaps = site_overlaps[row_idx, col_idx]
    pair_rates = mismatch_rates[row_idx, col_idx]
    valid_mask = (pair_overlaps >= overlap_threshold) & np.isfinite(pair_rates) & (pair_rates > 0)
    (pair_positions,) = np.nonzero(valid_mask)
    print(f"Valid pairs (overlap >= {overlap_threshold}): {len(pair_positions)}.")

    # Drop relatives / cryptic kinship: pairs more than `low_pmr_z_cut` sigma below the mean
    rates_valid = pair_rates[pair_positions]
    low_cut = rates_valid.mean() - low_pmr_z_cut * rates_valid.std()
    keep_mask = rates_valid >= low_cut
    n_dropped = int((~keep_mask).sum())
    pair_positions = pair_positions[keep_mask]
    print(f"Dropped {n_dropped} pairs with PMR < {low_cut:.6f} (z < -{low_pmr_z_cut}); {len(pair_positions)} remain.")

    rng = np.random.default_rng(seed)
    if len(pair_positions) > max_pairs:
        pair_positions = rng.choice(pair_positions, size=max_pairs, replace=False)
        pair_positions.sort()
        print(f"Subsampled to {max_pairs} pairs.")

    pi, pj = triu_i[pair_positions], triu_j[pair_positions]
    # Project to equal-area km coords so kernels operate on true Euclidean distances
    x_km, y_km = project_coords(lats, lons)
    pair_coords = np.column_stack([x_km[pi], y_km[pi], dates[pi], x_km[pj], y_km[pj], dates[pj]])
    pmr = pair_rates[pair_positions].astype(np.float32)
    return pair_coords, pmr


class SwapSymmetricKernel(gpytorch.kernels.Kernel):
    """Swap-symmetric wrapper for unordered sample pairs.

    K(p1, p2) = sum over (a in {s_i1, s_j1}, b in {s_i2, s_j2}) k(a, b),
    which is the covariance of f(s_i) + f(s_j) when f ~ GP(0, k).
    """

    has_lengthscale = False

    def __init__(self, base_kernel: gpytorch.kernels.Kernel, **kwargs: object) -> None:
        super().__init__(**kwargs)
        self.base_kernel = base_kernel

    def forward(  # type: ignore[override]
        self,
        x1: torch.Tensor,
        x2: torch.Tensor,
        diag: bool = False,
        **params: object,
    ) -> torch.Tensor:
        # Split each 6-d pair encoding into its two 3-d points: (..., n, 3) each
        si1, sj1 = x1[..., :3], x1[..., 3:]
        si2, sj2 = x2[..., :3], x2[..., 3:]
        # Stack along a new batch dim of size 4 so that a[t], b[t] enumerate
        # all 2x2 point combinations: (i1, i2), (i1, j2), (j1, i2), (j1, j2)
        a = torch.stack([si1, si1, sj1, sj1])
        b = torch.stack([si2, sj2, si2, sj2])
        # One batched kernel call: K[t] = k(a[t], b[t]), shape (4, ..., n, m)
        # (or (4, ..., n) if diag=True). Hyperparameters broadcast over the 4.
        K = self.base_kernel(a, b, diag=diag)
        if not isinstance(K, torch.Tensor):
            K = K.to_dense()
        # By bilinearity of covariance,
        # Cov(f(si1) + f(sj1), f(si2) + f(sj2))
        #   = k(si1, si2) + k(si1, sj2) + k(sj1, si2) + k(sj1, sj2),
        # which is exactly K[0] + K[1] + K[2] + K[3]
        return K.sum(0)


class WithinPairDistanceKernel(gpytorch.kernels.Kernel):
    """Kernel acting on (d_space, d_time) gaps between the two samples in a pair.

    Features are d_space = ||(lat_i, lon_i) - (lat_j, lon_j)|| and d_time = |t_i - t_j|,
    both of which are invariant under (i, j) swap.
    """

    has_lengthscale = False

    def __init__(self, base_kernel: gpytorch.kernels.Kernel, **kwargs: object) -> None:
        super().__init__(**kwargs)
        self.base_kernel = base_kernel

    @staticmethod
    def _features(x: torch.Tensor) -> torch.Tensor:
        # eps keeps d/dx sqrt(...) finite when two samples in a pair are co-located
        sq = (x[..., 0] - x[..., 3]) ** 2 + (x[..., 1] - x[..., 4]) ** 2
        d_space = torch.sqrt(sq + 1e-12)
        d_time = torch.abs(x[..., 2] - x[..., 5])
        return torch.stack([d_space, d_time], dim=-1)

    def forward(  # type: ignore[override]
        self,
        x1: torch.Tensor,
        x2: torch.Tensor,
        diag: bool = False,
        **params: object,
    ) -> torch.Tensor:
        f1 = self._features(x1)
        f2 = self._features(x2)
        K = self.base_kernel(f1, f2, diag=diag)
        if not isinstance(K, torch.Tensor):
            K = K.to_dense()
        return K


class PairwisePMRModel(gpytorch.models.ApproximateGP):
    """Variational GP for pairwise mismatch rate prediction."""

    def __init__(self, inducing_points: torch.Tensor) -> None:
        variational_distribution = gpytorch.variational.CholeskyVariationalDistribution(inducing_points.size(0))
        variational_strategy = gpytorch.variational.VariationalStrategy(
            self, inducing_points, variational_distribution, learn_inducing_locations=True
        )
        super().__init__(variational_strategy)
        self.mean_module = gpytorch.means.ConstantMean()

        # Per-sample latent f: separable Matern-5/2(x, y) * Matern-5/2(time). Default
        # lengthscale of 1 is too small for km / years, so initialise to reasonable
        # scales to keep the initial kernel matrix well-conditioned.
        self.space_kernel = gpytorch.kernels.MaternKernel(nu=2.5, active_dims=[0, 1])
        self.time_kernel = gpytorch.kernels.MaternKernel(nu=2.5, active_dims=[2])
        self.space_kernel.lengthscale = 1000.0
        self.time_kernel.lengthscale = 500.0
        self.swap_kernel = SwapSymmetricKernel(gpytorch.kernels.ScaleKernel(self.space_kernel * self.time_kernel))

        # Within-pair g: separable Matern-5/2 on (d_space, d_time) features
        self.within_space_kernel = gpytorch.kernels.MaternKernel(nu=2.5, active_dims=[0])
        self.within_time_kernel = gpytorch.kernels.MaternKernel(nu=2.5, active_dims=[1])
        self.within_space_kernel.lengthscale = 1000.0
        self.within_time_kernel.lengthscale = 500.0
        self.within_kernel = WithinPairDistanceKernel(
            gpytorch.kernels.ScaleKernel(self.within_space_kernel * self.within_time_kernel)
        )

        self.covar_module = self.swap_kernel + self.within_kernel

    def forward(self, x: torch.Tensor) -> gpytorch.distributions.MultivariateNormal:
        return gpytorch.distributions.MultivariateNormal(self.mean_module(x), self.covar_module(x))


def train_model(
    model: PairwisePMRModel,
    likelihood: gpytorch.likelihoods.GaussianLikelihood,
    train_x: torch.Tensor,
    train_y: torch.Tensor,
    n_epochs: int = N_EPOCHS,
    batch_size: int = BATCH_SIZE,
    lr: float = LEARNING_RATE,
) -> None:
    """Train variational GP by maximising the ELBO with Adam."""
    model.train()
    likelihood.train()

    params = list(model.parameters()) + list(likelihood.parameters())
    optimizer = torch.optim.Adam(params, lr=lr)
    mll = gpytorch.mlls.VariationalELBO(likelihood, model, num_data=train_y.size(0))

    loader = DataLoader(
        TensorDataset(train_x, train_y),
        batch_size=batch_size,
        shuffle=True,
        generator=torch.Generator().manual_seed(SEED),
    )

    with _linalg_dtype_cholesky(torch.float32), gpytorch.settings.cholesky_jitter(float_value=1e-4):
        for epoch in range(n_epochs):
            epoch_loss = 0.0
            for xb, yb in loader:
                optimizer.zero_grad()
                loss = -mll(model(xb), yb)
                loss.backward()
                optimizer.step()
                epoch_loss += loss.item()

            avg_loss = epoch_loss / len(loader)

            if (epoch + 1) % 5 == 0 or epoch == 0:
                swap_scale = model.swap_kernel.base_kernel
                within_scale = model.within_kernel.base_kernel
                print(
                    f"  epoch {epoch + 1:4d}/{n_epochs} | loss {avg_loss:.4f}"
                    f" | sigma2_pair {swap_scale.outputscale.item():.4f}"
                    f" | l_space {model.space_kernel.lengthscale.item():.1f} km"
                    f" | l_time {model.time_kernel.lengthscale.item():.0f} yr"
                    f" | sigma2_within {within_scale.outputscale.item():.4f}"
                    f" | l_dspace {model.within_space_kernel.lengthscale.item():.1f} km"
                    f" | l_dtime {model.within_time_kernel.lengthscale.item():.0f} yr"
                )


@torch.no_grad()
def predict_grid(
    model: PairwisePMRModel,
    lat_range: tuple[float, float],
    lon_range: tuple[float, float],
    time_bp: float,
    resolution_deg: float = GRID_RESOLUTION_DEG,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Predict expected PMR on a lat/lon grid at a fixed time."""
    model.eval()

    lat_vec = np.arange(lat_range[0], lat_range[1] + resolution_deg / 2, resolution_deg)
    lon_vec = np.arange(lon_range[0], lon_range[1] + resolution_deg / 2, resolution_deg)
    lon_grid, lat_grid = np.meshgrid(lon_vec, lat_vec)
    # Project to equal-area km coords to match the training input space
    x_km, y_km = project_coords(lat_grid.ravel(), lon_grid.ravel())
    grid_flat = np.column_stack(
        [
            x_km,
            y_km,
            np.full(x_km.size, time_bp),
        ]
    ).astype(np.float32)

    device = next(model.parameters()).device
    # Predict at coincident pairs (s*, s*)
    x_pred = torch.from_numpy(np.column_stack([grid_flat, grid_flat])).to(device)

    means: list[np.ndarray] = []
    chunk = 2048
    with (
        gpytorch.settings.fast_pred_var(),
        _linalg_dtype_cholesky(torch.float32),
        gpytorch.settings.cholesky_jitter(float_value=1e-4),
    ):
        for i in range(0, len(x_pred), chunk):
            pred = model(x_pred[i : i + chunk])
            means.append(pred.mean.cpu().numpy())

    mean_grid = np.concatenate(means).reshape(lat_grid.shape)
    return lat_vec, lon_vec, mean_grid


def plot_predictions(
    lat_vec: np.ndarray,
    lon_vec: np.ndarray,
    mean_grids: dict[float, np.ndarray],
    output_path: Path,
) -> None:
    """Plot predicted PMR maps for several time slices."""
    n = len(mean_grids)
    ncols = min(n, 2)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(10, 10),
        subplot_kw={"projection": PROJECTION},
        constrained_layout=True,
    )
    axes_flat = np.atleast_1d(axes).ravel()
    lon_grid, lat_grid = np.meshgrid(lon_vec, lat_vec)
    boundary = lat_lon_box_path(MAP_EXTENT_DEG)
    # Size the rectangular data window to the boundary's projected bounding box so
    # the curvilinear boundary fits inside without being cropped (cartopy's
    # set_extent in PlateCarree densifies independently and can clip the bulges)
    x_lo, y_lo = boundary.vertices.min(axis=0)
    x_hi, y_hi = boundary.vertices.max(axis=0)

    all_values = np.concatenate([g.ravel() for g in mean_grids.values()])
    vmin = float(np.nanpercentile(all_values, 2))
    vmax = float(np.nanpercentile(all_values, 98))

    mesh = None
    for ax, (t_bp, mean_grid) in zip(axes_flat, mean_grids.items(), strict=True):
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_boundary(boundary)
        ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=0)
        ax.add_feature(cfeature.LAND, facecolor="#f0f0f0", zorder=0)
        mesh = ax.pcolormesh(
            lon_grid,
            lat_grid,
            mean_grid,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            transform=ccrs.PlateCarree(),
            shading="auto",
            zorder=1,
        )
        # Overlay ocean on top of the data layer so color only shows on land
        ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=2)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor="#555555", zorder=3)
        ax.set_title(f"{int(t_bp)} BP", fontsize=14, y=0.97)

    for ax in axes_flat[n:]:
        ax.set_visible(False)

    if mesh is not None:
        cbar = fig.colorbar(
            mesh,
            ax=axes_flat[:n].tolist(),
            orientation="horizontal",
            shrink=0.5,
            aspect=40,
            pad=0.04,
            extend="both",
        )
        cbar.set_label("Predicted PMR at coincident points", fontsize=14)
        cbar.ax.tick_params(labelsize=14)

    fig.suptitle("Gaussian process-predicted PMRs", fontsize=16)
    fig.savefig(output_path, dpi=600)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Seed all torch RNGs so the variational distribution init is reproducible
    torch.manual_seed(SEED)

    print("Loading AADR pairwise data...\n")
    pair_coords, pmr = extract_pair_data()

    train_x = torch.from_numpy(pair_coords.astype(np.float32))
    train_y = torch.from_numpy(pmr)

    print(f"\nTraining pairs : {len(train_x)}")
    print(f"PMR range      : [{pmr.min():.4f}, {pmr.max():.4f}]")
    print(f"PMR mean / std : {pmr.mean():.4f} / {pmr.std():.4f}")

    device = get_device()
    print(f"Device         : {device}")

    rng = np.random.default_rng(SEED)
    inducing_idx = rng.choice(len(train_x), size=N_INDUCING, replace=False)
    inducing_points = train_x[inducing_idx].clone()

    model = PairwisePMRModel(inducing_points)
    # PMR variance is ~0.005, so the default noise lower bound of 1e-4 can clip data-driven
    # initial noise; loosen it so we can initialise from the data
    likelihood = gpytorch.likelihoods.GaussianLikelihood(noise_constraint=gpytorch.constraints.GreaterThan(1e-6))

    # Initialise scales from data so the first ELBO step is not dominated by
    # shrinking GPyTorch's default outputscale=1.0 and noise≈0.693 down to PMR variance
    pmr_var = float(pmr.var())
    with torch.no_grad():
        model.mean_module.constant.copy_(torch.tensor(float(pmr.mean())))
        model.swap_kernel.base_kernel.outputscale = pmr_var / 5.0
        model.within_kernel.base_kernel.outputscale = pmr_var / 5.0
    likelihood.noise = pmr_var * 0.1

    model.to(dtype=torch.float32, device=device)
    likelihood.to(dtype=torch.float32, device=device)
    train_x = train_x.to(device)
    train_y = train_y.to(device)

    print(f"\nTraining ({N_INDUCING} inducing pts, {N_EPOCHS} epochs, batch {BATCH_SIZE})...\n")
    train_model(model, likelihood, train_x, train_y)

    torch.save(
        {"model_state": model.state_dict(), "likelihood_state": likelihood.state_dict()},
        CHECKPOINT_PATH,
    )
    print(f"Wrote model checkpoint to {CHECKPOINT_PATH}.")

    lon_min, lon_max, lat_min, lat_max = MAP_EXTENT_DEG
    lat_range = (lat_min, lat_max)
    lon_range = (lon_min, lon_max)
    print(f"\nPredicting on grid (lat {lat_range}, lon {lon_range})...")

    mean_grids: dict[float, np.ndarray] = {}
    for t_bp in TIME_SLICES_BP:
        lat_vec, lon_vec, mean_grid = predict_grid(model, lat_range, lon_range, float(t_bp))
        mean_grids[float(t_bp)] = mean_grid
        print(f"  {t_bp} BP: mean PMR in [{mean_grid.min():.4f}, {mean_grid.max():.4f}]")

    plot_predictions(lat_vec, lon_vec, mean_grids, MAP_PATH)
    print(f"\nWrote PMR maps to {MAP_PATH}.")


if __name__ == "__main__":
    main()
