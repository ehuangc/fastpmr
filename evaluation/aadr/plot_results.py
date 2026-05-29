"""
1) Plot median within-locality pairwise mismatch rate (PMR) over space and time.

   For each (locality, time bin) cell with at least MIN_PAIRS_PER_CELL within-locality
   sample pairs with site overlap > OVERLAP_THRESHOLD, the median PMR is plotted
   on a global map stratified by time period.

2) Plot within-locality PMR vs. migratory distance from Addis Ababa, Ethiopia.
"""

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from haversine import haversine
from scipy import stats

from evaluation_utils import (
    AADR_ANOMALOUS_LOCALITY_PREFIXES,
    AADR_DIR,
    AADR_METADATA_PATH,
    AADR_NPZ_PATH,
    DATE_MEAN_BP_FIELD,
    INDIVIDUAL_ID_FIELD,
    LAT_FIELD,
    LOCALITY_FIELD,
    LON_FIELD,
    PUBLICATION_FIELD,
    classify_coords,
    ensure_aadr_npz_present,
    is_archaic_or_reference_sample,
    load_aadr_metadata,
    load_aadr_npz_arrays,
)

OUTPUT_DIR = AADR_DIR / "results" / "plots"
MAP_PATH = OUTPUT_DIR / "pmr_map_by_time.pdf"
REGRESSION_PATH = OUTPUT_DIR / "pmr_vs_migratory_distance.pdf"
BOXPLOT_PATH = OUTPUT_DIR / "pmr_boxplot_by_time.pdf"
CELLS_CSV_PATH = OUTPUT_DIR / "localities.csv"

OVERLAP_THRESHOLD = 30_000
MIN_PAIRS_PER_CELL = 5
TIME_BINS: list[tuple[float, float, str]] = [
    (10_000.0, np.inf, ">10,000 BP"),
    (5_000.0, 10_000.0, "10,000–5,000 BP"),
    (2_000.0, 5_000.0, "5,000–2,000 BP"),
    (500.0, 2_000.0, "2,000–500 BP"),
]
PALETTE = sns.color_palette("colorblind")
BIN_COLORS = {
    ">10,000 BP": PALETTE[0],  # Blue
    "10,000–5,000 BP": PALETTE[1],  # Orange
    "5,000–2,000 BP": PALETTE[4],  # Purple
    "2,000–500 BP": PALETTE[2],  # Teal
}

# Waypoints approximate overland routes around major water bodies, inspired by the
# approach of Ramachandran et al. 2005 (https://doi.org/10.1073/pnas.0507611102)
WAYPOINT_COORDS = {
    "addis_ababa": (9, 38),
    "cairo": (30, 31),
    "istanbul": (41, 28),
    "phnom_penh": (11, 104),
    "anadyr": (64, 177),
    "prince_rupert": (54, -130),
}
REGION_PATHS = {
    "africa": ["addis_ababa"],
    "asia": ["addis_ababa", "cairo"],
    "europe": ["addis_ababa", "cairo", "istanbul"],
    "oceania": ["addis_ababa", "cairo", "phnom_penh"],
    "americas": ["addis_ababa", "cairo", "anadyr", "prince_rupert"],
}


def filter_and_extract(
    samples: list[str],
    covered_snps: np.ndarray,
    site_overlaps: np.ndarray,
    mismatch_rates: np.ndarray,
    metadata: dict[str, dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # First pass: collect every sample that passes the per-sample filters,
    # along with its Individual ID and covered-SNP count for deduplication
    candidates: list[tuple[int, float, float, float, str, str]] = []
    candidate_individual_ids: list[str] = []
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
        if locality.startswith(AADR_ANOMALOUS_LOCALITY_PREFIXES):
            continue
        lat = parse_float(sample_metadata[LAT_FIELD])
        lon = parse_float(sample_metadata[LON_FIELD])
        date = parse_float(sample_metadata[DATE_MEAN_BP_FIELD])
        if not (np.isfinite(lat) and np.isfinite(lon) and np.isfinite(date)):
            continue

        individual_id = sample_metadata[INDIVIDUAL_ID_FIELD].strip()
        publication = sample_metadata[PUBLICATION_FIELD].strip()
        candidates.append((sample_idx, lat, lon, date, locality, publication))
        candidate_individual_ids.append(individual_id)
        candidate_covered_snps.append(int(covered_snps[sample_idx]))

    # Drop known duplicates by AADR individual ID, keeping the sample with the most covered SNPs in each group
    # Samples without a individual ID are kept as-is
    best_pos_by_individual: dict[str, int] = {}
    for pos, individual_id in enumerate(candidate_individual_ids):
        if not individual_id or individual_id == "..":
            continue
        prev = best_pos_by_individual.get(individual_id)
        if prev is None or candidate_covered_snps[pos] > candidate_covered_snps[prev]:
            best_pos_by_individual[individual_id] = pos
    selected_positions = {pos for pos, mid in enumerate(candidate_individual_ids) if not mid or mid == ".."}
    selected_positions.update(best_pos_by_individual.values())

    keep: list[int] = []
    lats: list[float] = []
    lons: list[float] = []
    dates: list[float] = []
    localities: list[str] = []
    publications: list[str] = []
    for pos, (sample_idx, lat, lon, date, locality, publication) in enumerate(candidates):
        if pos not in selected_positions:
            continue
        keep.append(sample_idx)
        lats.append(lat)
        lons.append(lon)
        dates.append(date)
        localities.append(locality)
        publications.append(publication)

    keep_arr = np.asarray(keep, dtype=int)
    return (
        site_overlaps[np.ix_(keep_arr, keep_arr)],
        mismatch_rates[np.ix_(keep_arr, keep_arr)],
        np.asarray(lats, dtype=float),
        np.asarray(lons, dtype=float),
        np.asarray(dates, dtype=float),
        np.asarray(localities, dtype=object),
        np.asarray(publications, dtype=object),
    )


def assign_time_bins(dates: np.ndarray) -> np.ndarray:
    bin_indices = np.full(dates.shape, -1, dtype=int)
    for idx, (low, high, _label) in enumerate(TIME_BINS):
        mask = (dates >= low) & (dates < high)
        bin_indices[mask] = idx
    return bin_indices


def compute_locality_cells(
    bin_indices: np.ndarray,
    site_overlaps: np.ndarray,
    mismatch_rates: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
    localities: np.ndarray,
    publications: np.ndarray,
) -> pd.DataFrame:
    """For each (locality, time bin) group with enough valid within-locality pairs,
    compute median PMR + mean coordinates."""
    groups: dict[tuple[str, int], list[int]] = {}
    for idx, (locality, bin_idx) in enumerate(zip(localities, bin_indices, strict=True)):
        if bin_idx < 0:
            continue
        groups.setdefault((locality, int(bin_idx)), []).append(idx)

    rows: list[dict[str, object]] = []
    for (locality, bin_idx), sample_indices in groups.items():
        sample_indices_arr = np.asarray(sample_indices, dtype=int)
        pubs = sorted({publications[i] for i in sample_indices_arr if publications[i]})
        sub_rates = mismatch_rates[np.ix_(sample_indices_arr, sample_indices_arr)]
        sub_overlaps = site_overlaps[np.ix_(sample_indices_arr, sample_indices_arr)]

        triu_i, triu_j = np.triu_indices(len(sample_indices_arr), k=1)
        pair_rates = sub_rates[triu_i, triu_j]
        pair_overlaps = sub_overlaps[triu_i, triu_j]
        valid = (pair_overlaps >= OVERLAP_THRESHOLD) & ~np.isnan(pair_rates)
        n_valid = int(valid.sum())
        if n_valid < MIN_PAIRS_PER_CELL:
            continue
        median_pmr = round(float(np.median(pair_rates[valid])), 6)
        rows.append(
            {
                "locality": locality,
                "publications": "; ".join(pubs),
                "bin_idx": bin_idx,
                "bin_label": TIME_BINS[bin_idx][2],
                "n_samples": len(sample_indices),
                "n_pairs": n_valid,
                "median_pmr": median_pmr,
                "lat": float(np.mean(lats[sample_indices_arr])),
                "lon": float(np.mean(lons[sample_indices_arr])),
            }
        )
    columns = ["locality", "publications", "bin_idx", "bin_label", "n_samples", "n_pairs", "median_pmr", "lat", "lon"]
    return pd.DataFrame.from_records(rows, columns=columns)


def plot_map(cells: pd.DataFrame, output_path: Path) -> None:
    pmr_values = cells["median_pmr"].to_numpy()
    vmin, vmax = np.percentile(pmr_values, [2, 98])

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(10, 6.75),
        subplot_kw={"projection": ccrs.Robinson(central_longitude=150)},
        constrained_layout=True,
        gridspec_kw={"hspace": 0.05},
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
        ax.set_title(f"{label}\n(n={len(bin_cells)} localities)", fontsize=14, pad=6)
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
            pad=0.04,
        )
        cbar.set_label("Median within-locality PMR", fontsize=14)
        cbar.ax.tick_params(labelsize=14)

    fig.suptitle("AADR median PMRs by locality and time period", fontsize=16)
    fig.savefig(output_path, dpi=600)
    plt.close(fig)


def migratory_distance(lat: float, lon: float, region: str) -> float:
    """Approximate migratory distance (km) from Addis Ababa to (lat, lon).

    Routes through a fixed sequence of waypoints determined by the
    destination's geographic region, then adds the great-circle distance
    from the last waypoint to the destination.
    """
    path = REGION_PATHS[region]
    total = 0.0
    for i in range(len(path) - 1):
        total += haversine(WAYPOINT_COORDS[path[i]], WAYPOINT_COORDS[path[i + 1]])
    total += haversine(WAYPOINT_COORDS[path[-1]], (lat, lon))
    return total


def plot_pmr_vs_migratory_distance(cells: pd.DataFrame, output_path: Path) -> None:
    """Scatter plot + linear regression of within-locality PMR vs. migratory distance."""
    x = cells["migratory_distance_km"].to_numpy() / 1_000  # Thousands of km
    y = cells["median_pmr"].to_numpy()
    slope, intercept, r_value, _p_value, _se = stats.linregress(x, y)

    fig, ax = plt.subplots(figsize=(10, 6.875), constrained_layout=True)
    n_pairs = cells["n_pairs"].to_numpy()
    sizes = 30 + 5 * np.clip(n_pairs, 0, 20)
    for bin_i, (_low, _high, label) in enumerate(TIME_BINS):
        mask = cells["bin_idx"].to_numpy() == bin_i
        ax.scatter(
            x[mask],
            y[mask],
            label=label,
            color=BIN_COLORS[label],
            s=sizes[mask],
            alpha=0.7,
            edgecolor="black",
            linewidth=0.3,
        )
    sns.regplot(
        x=x,
        y=y,
        ci=None,
        scatter=False,
        color="black",
        line_kws={"linewidth": 1.5, "alpha": 0.5},
        label=f"OLS ($R^2 = {r_value**2:.3f}$)",
        ax=ax,
    )

    # Circle outliers
    outliers = [
        (
            "Norwich (England, Norwich, Chapelfield Shopping Center, well shaft)",
            "Medieval Ashkenazi\nJewish individuals",
        ),
        ("(Eastern Settlement)", "Norse Greenland\nindividuals"),
    ]
    localities = cells["locality"].to_numpy()
    for match_str, label in outliers:
        idxs = [i for i, loc in enumerate(localities) if match_str in loc]
        ox = np.mean([x[i] for i in idxs])
        oy = np.mean([y[i] for i in idxs])
        ax.scatter(
            [ox],
            [oy],
            s=400,
            facecolors="none",
            edgecolors="black",
            alpha=0.5,
        )
        ax.annotate(
            label,
            xy=(ox, oy),
            xytext=(18, -18),
            textcoords="offset points",
            fontsize=11,
            bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "black", "alpha": 0.5},
        )

    ax.set_xlabel("Waypoint distance from Addis Ababa (×1000 km)", fontsize=16)
    ax.set_ylabel("Median within-locality PMR", fontsize=16)
    ax.set_title('AADR median within-locality PMR vs. migratory "Out of Africa" distance', fontsize=16)
    ax.tick_params(axis="both", labelsize=14)
    ax.legend(fontsize=11, loc="upper right")

    fig.savefig(output_path, dpi=600)
    plt.close(fig)


def format_p_value(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def draw_significance_brackets(
    ax: plt.Axes, data: pd.DataFrame, bins: list[str], bracket_gap: float, tip_len: float
) -> None:
    """Draw brackets between time bins, using Mann-Whitney U p-values"""
    groups = [data.loc[data["bin_label"] == b, "median_pmr"].to_numpy() for b in bins]
    comparisons = [(0, 1), (1, 2), (0, 2)]

    y_max = max(g.max() for g in groups if len(g) > 0)

    # Sort comparisons by span so narrower brackets are drawn lower
    comparisons.sort(key=lambda pair: pair[1] - pair[0])

    y_cursor = y_max + bracket_gap
    for i, j in comparisons:
        _stat, p = stats.mannwhitneyu(groups[i], groups[j], alternative="two-sided")
        label = format_p_value(p)

        y = y_cursor
        xi = i
        xj = j
        ax.plot([xi, xi, xj, xj], [y - tip_len, y, y, y - tip_len], lw=1.0, color="black")
        ax.text((i + j) / 2, y, label, ha="center", va="bottom", fontsize=12)
        y_cursor = y + bracket_gap + tip_len


def plot_pmr_box_plot(cells: pd.DataFrame, output_path: Path) -> None:
    """Multi-panel box plot of median within-locality PMRs,
    stratified by time period for Europe, Asia, and the Americas."""
    bins = ["10,000–5,000 BP", "5,000–2,000 BP", "2,000–500 BP"]
    panels = [
        ("europe", "Europe"),
        ("asia", "Asia"),
        ("americas", "Americas"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(10, 4.5), constrained_layout=True, sharey=True)

    # Compute shared bracket spacing from the global y-range
    all_pmr = cells[cells["bin_label"].isin(bins)]["median_pmr"]
    y_range = all_pmr.max() - all_pmr.min()
    bracket_gap = 0.06 * y_range
    tip_len = 0.012 * y_range

    for ax, (region, title) in zip(axes, panels, strict=True):
        region_cells = cells[(cells["region"] == region) & (cells["bin_label"].isin(bins))].copy()
        region_cells["bin_label"] = pd.Categorical(region_cells["bin_label"], categories=bins, ordered=True)
        sns.boxplot(
            data=region_cells,
            x="bin_label",
            y="median_pmr",
            hue="bin_label",
            palette=BIN_COLORS,
            dodge=False,
            fliersize=4,
            width=0.6,
            legend=False,
            ax=ax,
        )
        ax.set_xlabel("")
        ax.set_xticks(range(len(bins)))
        ax.set_xticklabels([""] * len(bins))
        ax.tick_params(axis="x", length=0)
        ax.margins(x=0.15)
        ax.set_title(title, fontsize=14)
        ax.tick_params(axis="y", labelsize=14)
        draw_significance_brackets(ax, region_cells, bins, bracket_gap, tip_len)

    # Set shared y-limits: preserve auto lower bound, add modest headroom on top
    y_lo, y_hi = axes[0].get_ylim()
    axes[0].set_ylim(y_lo, y_hi + 0.03 * (y_hi - y_lo))

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=BIN_COLORS[b], edgecolor="black", linewidth=0.8) for b in bins]
    fig.legend(handles, bins, loc="outside lower center", ncol=len(bins), fontsize=14, frameon=False)

    axes[0].set_ylabel("Median within-locality PMR", fontsize=14)
    fig.suptitle("AADR median within-locality PMRs by region and time period", fontsize=16)
    fig.savefig(output_path, dpi=600)
    plt.close(fig)


def main() -> None:
    ensure_aadr_npz_present(AADR_NPZ_PATH, AADR_METADATA_PATH)
    samples, site_overlaps, mismatch_rates, _mismatch_rates_95_ci_lower, _mismatch_rates_95_ci_upper, covered_snps = (
        load_aadr_npz_arrays(AADR_NPZ_PATH)
    )
    metadata = load_aadr_metadata(AADR_METADATA_PATH)
    assert len(samples) == mismatch_rates.shape[0] == site_overlaps.shape[0]

    site_overlaps_filt, mismatch_rates_filt, lats, lons, dates, localities, publications = filter_and_extract(
        samples, covered_snps, site_overlaps, mismatch_rates, metadata
    )
    bin_indices = assign_time_bins(dates)
    cells = compute_locality_cells(
        bin_indices, site_overlaps_filt, mismatch_rates_filt, lats, lons, localities, publications
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cells_sorted = cells.sort_values(["bin_idx", "median_pmr"]).reset_index(drop=True)
    cell_coords = list(zip(cells_sorted["lat"], cells_sorted["lon"], strict=True))
    cell_regions = classify_coords(cell_coords)
    cells_sorted["region"] = cell_regions
    cells_sorted["migratory_distance_km"] = [
        round(migratory_distance(lat, lon, region), 1)
        for lat, lon, region in zip(cells_sorted["lat"], cells_sorted["lon"], cell_regions, strict=True)
    ]
    cells_sorted.to_csv(CELLS_CSV_PATH, index=False)

    counts = cells_sorted.groupby("bin_label", sort=False).size().to_dict()
    print(f"Wrote {len(cells_sorted)} (locality, time bin) cells to {CELLS_CSV_PATH}.")
    for _low, _high, label in TIME_BINS:
        print(f"  {label}: {counts.get(label, 0)} cells")

    plot_map(cells_sorted, MAP_PATH)
    print(f"\nWrote PMR map to {MAP_PATH}.")

    plot_pmr_vs_migratory_distance(cells_sorted, REGRESSION_PATH)
    print(f"\nWrote PMR vs. migratory distance plot to {REGRESSION_PATH}.")

    plot_pmr_box_plot(cells_sorted, BOXPLOT_PATH)
    print(f"\nWrote PMR box plot to {BOXPLOT_PATH}.")


if __name__ == "__main__":
    main()
