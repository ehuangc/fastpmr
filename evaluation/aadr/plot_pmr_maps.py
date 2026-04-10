"""
Plot median within-site pairwise mismatch rate (PMR) over space and time.

For each (locality, time bin) cell with at least MIN_PAIRS_PER_CELL within-site
sample pairs with site overlap > OVERLAP_THRESHOLD, the median PMR is plotted
on a global map stratified by time period.
"""

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from evaluation_utils import (
    AADR_DIR,
    AADR_EXCLUDED_LOCALITY_PREFIXES,
    AADR_METADATA_PATH,
    AADR_NPZ_PATH,
    DATE_MEAN_BP_FIELD,
    LAT_FIELD,
    LOCALITY_FIELD,
    LON_FIELD,
    MASTER_ID_FIELD,
    ensure_aadr_npz_present,
    is_archaic_or_reference_sample,
    load_aadr_metadata,
    load_aadr_npz_arrays,
)

OUTPUT_DIR = AADR_DIR / "results" / "maps"
MAP_PATH = OUTPUT_DIR / "pmr_map_by_time.pdf"
CELLS_CSV_PATH = OUTPUT_DIR / "pmr_map_cells.csv"

OVERLAP_THRESHOLD = 30_000
MIN_PAIRS_PER_CELL = 5

TIME_BINS: list[tuple[float, float, str]] = [
    (10_000.0, np.inf, ">10,000 BP"),
    (5_000.0, 10_000.0, "10,000-5,000 BP"),
    (2_000.0, 5_000.0, "5,000-2,000 BP"),
    (-np.inf, 2_000.0, "<2,000 BP"),
]


def filter_and_extract(
    samples: list[str],
    covered_snps: np.ndarray,
    site_overlaps: np.ndarray,
    mismatch_rates: np.ndarray,
    metadata: dict[str, dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # First pass: collect every sample that passes the per-sample filters,
    # along with its Master ID and covered-SNP count for deduplication.
    candidates: list[tuple[int, float, float, float, str]] = []
    candidate_master_ids: list[str] = []
    candidate_covered_snps: list[int] = []

    def parse_float(value: str) -> float:
        try:
            return float(value)
        except (TypeError, ValueError):
            return float("nan")

    for sample_idx, sample in enumerate(samples):
        sample_metadata = metadata[sample]
        if is_archaic_or_reference_sample(sample, sample_metadata):
            continue
        locality = sample_metadata[LOCALITY_FIELD].strip()
        if not locality or locality == ".." or locality.lower() == "n/a":
            continue
        if locality.startswith(AADR_EXCLUDED_LOCALITY_PREFIXES):
            continue
        lat = parse_float(sample_metadata[LAT_FIELD])
        lon = parse_float(sample_metadata[LON_FIELD])
        date = parse_float(sample_metadata[DATE_MEAN_BP_FIELD])
        if not (np.isfinite(lat) and np.isfinite(lon) and np.isfinite(date)):
            continue

        master_id = sample_metadata[MASTER_ID_FIELD].strip()
        candidates.append((sample_idx, lat, lon, date, locality))
        candidate_master_ids.append(master_id)
        candidate_covered_snps.append(int(covered_snps[sample_idx]))

    # Drop known duplicates by AADR master ID, keeping the sample with the most covered SNPs in each group.
    # Samples without a master ID are kept as-is.
    best_pos_by_master: dict[str, int] = {}
    for pos, master_id in enumerate(candidate_master_ids):
        if not master_id or master_id == "..":
            continue
        prev = best_pos_by_master.get(master_id)
        if prev is None or candidate_covered_snps[pos] > candidate_covered_snps[prev]:
            best_pos_by_master[master_id] = pos
    selected_positions = {pos for pos, mid in enumerate(candidate_master_ids) if not mid or mid == ".."}
    selected_positions.update(best_pos_by_master.values())

    keep: list[int] = []
    lats: list[float] = []
    lons: list[float] = []
    dates: list[float] = []
    localities: list[str] = []
    for pos, (sample_idx, lat, lon, date, locality) in enumerate(candidates):
        if pos not in selected_positions:
            continue
        keep.append(sample_idx)
        lats.append(lat)
        lons.append(lon)
        dates.append(date)
        localities.append(locality)

    keep_arr = np.asarray(keep, dtype=int)
    return (
        site_overlaps[np.ix_(keep_arr, keep_arr)],
        mismatch_rates[np.ix_(keep_arr, keep_arr)],
        np.asarray(lats, dtype=float),
        np.asarray(lons, dtype=float),
        np.asarray(dates, dtype=float),
        np.asarray(localities, dtype=object),
    )


def assign_time_bins(dates: np.ndarray) -> np.ndarray:
    bin_indices = np.full(dates.shape, -1, dtype=int)
    for idx, (low, high, _label) in enumerate(TIME_BINS):
        mask = (dates >= low) & (dates < high)
        bin_indices[mask] = idx
    return bin_indices


def compute_site_cells(
    site_overlaps: np.ndarray,
    mismatch_rates: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
    localities: np.ndarray,
    bin_indices: np.ndarray,
) -> pd.DataFrame:
    """For each (locality, time bin) group with enough valid within-site pairs,
    compute median PMR + mean coordinates."""
    groups: dict[tuple[str, int], list[int]] = {}
    for idx, (locality, bin_idx) in enumerate(zip(localities, bin_indices, strict=True)):
        if bin_idx < 0:
            continue
        groups.setdefault((locality, int(bin_idx)), []).append(idx)

    rows: list[dict[str, object]] = []
    for (locality, bin_idx), sample_indices in groups.items():
        sample_indices_arr = np.asarray(sample_indices, dtype=int)
        sub_rates = mismatch_rates[np.ix_(sample_indices_arr, sample_indices_arr)]
        sub_overlaps = site_overlaps[np.ix_(sample_indices_arr, sample_indices_arr)]

        triu_i, triu_j = np.triu_indices(len(sample_indices_arr), k=1)
        pair_rates = sub_rates[triu_i, triu_j]
        pair_overlaps = sub_overlaps[triu_i, triu_j]
        valid = (pair_overlaps >= OVERLAP_THRESHOLD) & ~np.isnan(pair_rates)
        n_valid = int(valid.sum())
        if n_valid < MIN_PAIRS_PER_CELL:
            continue
        median_pmr = float(np.median(pair_rates[valid]))
        rows.append(
            {
                "locality": locality,
                "bin_idx": bin_idx,
                "bin_label": TIME_BINS[bin_idx][2],
                "n_samples": len(sample_indices),
                "n_pairs": n_valid,
                "median_pmr": median_pmr,
                "lat": float(np.mean(lats[sample_indices_arr])),
                "lon": float(np.mean(lons[sample_indices_arr])),
            }
        )
    columns = ["locality", "bin_idx", "bin_label", "n_samples", "n_pairs", "median_pmr", "lat", "lon"]
    return pd.DataFrame.from_records(rows, columns=columns)


def plot_map(cells: pd.DataFrame, output_path: Path) -> None:
    pmr_values = cells["median_pmr"].to_numpy()
    vmin, vmax = np.percentile(pmr_values, [2, 98])

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(14, 8.75),
        subplot_kw={"projection": ccrs.Robinson(central_longitude=150)},
        constrained_layout=True,
    )
    axes_flat = axes.ravel()

    scatter = None
    for bin_i, (_low, _high, label) in enumerate(TIME_BINS):
        ax = axes_flat[bin_i]
        ax.set_global()
        # Crop Antarctica by adjusting the y-limit in projected coordinates
        _x, y_south = ax.projection.transform_point(0, -60, ccrs.PlateCarree())
        ax.set_ylim(bottom=y_south)
        ax.add_feature(cfeature.LAND, facecolor="#f0f0f0", zorder=0)
        ax.add_feature(cfeature.OCEAN, facecolor="#ffffff", zorder=0)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor="#555555")

        bin_cells = cells[cells["bin_idx"] == bin_i]
        ax.set_title(f"{label}\n(n={len(bin_cells)} sites)", fontsize=12)
        if bin_cells.empty:
            continue

        sizes = 20 + 3 * np.clip(bin_cells["n_pairs"].to_numpy(), 0, 20)
        scatter = ax.scatter(
            bin_cells["lon"].to_numpy(),
            bin_cells["lat"].to_numpy(),
            c=bin_cells["median_pmr"].to_numpy(),
            s=sizes,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            transform=ccrs.PlateCarree(),
            edgecolor="black",
            linewidth=0.3,
            alpha=0.9,
        )

    if scatter is not None:
        cbar = fig.colorbar(
            scatter,
            ax=axes_flat.tolist(),
            orientation="horizontal",
            shrink=0.5,
            aspect=40,
            extend="both",
        )
        cbar.set_label("Median within-site pairwise mismatch rate", fontsize=12)
        cbar.ax.tick_params(labelsize=12)

    fig.suptitle("AADR median pairwise mismatch rates, by site and era", fontsize=14)
    fig.savefig(output_path, dpi=600)
    plt.close(fig)


def main() -> None:
    ensure_aadr_npz_present(AADR_NPZ_PATH, AADR_METADATA_PATH)
    samples, site_overlaps, mismatch_rates, _mismatch_rates_95_ci_lower, _mismatch_rates_95_ci_upper, covered_snps = (
        load_aadr_npz_arrays(AADR_NPZ_PATH)
    )
    metadata = load_aadr_metadata(AADR_METADATA_PATH)
    assert len(samples) == mismatch_rates.shape[0] == site_overlaps.shape[0]

    site_overlaps_filt, mismatch_rates_filt, lats, lons, dates, localities = filter_and_extract(
        samples, covered_snps, site_overlaps, mismatch_rates, metadata
    )
    bin_indices = assign_time_bins(dates)
    cells = compute_site_cells(site_overlaps_filt, mismatch_rates_filt, lats, lons, localities, bin_indices)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cells_sorted = cells.sort_values(["bin_idx", "median_pmr"]).reset_index(drop=True)
    cells_sorted.to_csv(CELLS_CSV_PATH, index=False)

    counts = cells_sorted.groupby("bin_label", sort=False).size().to_dict()
    print(f"Wrote {len(cells_sorted)} (locality, time bin) cells to {CELLS_CSV_PATH}.")
    for _low, _high, label in TIME_BINS:
        print(f"  {label}: {counts.get(label, 0)} cells")

    plot_map(cells_sorted, MAP_PATH)
    print(f"\nWrote PMR map to {MAP_PATH}.")


if __name__ == "__main__":
    main()
