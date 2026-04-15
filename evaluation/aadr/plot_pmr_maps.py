"""
1) Plot median within-site pairwise mismatch rate (PMR) over space and time.

   For each (locality, time bin) cell with at least MIN_PAIRS_PER_CELL within-site
   sample pairs with site overlap > OVERLAP_THRESHOLD, the median PMR is plotted
   on a global map stratified by time period.

2) Plot within-site PMR vs. migratory distance from Addis Ababa, Ethiopia.
"""

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import reverse_geocoder as rg
import seaborn as sns
from haversine import haversine
from scipy import stats

from evaluation_utils import (
    AADR_ANOMALOUS_LOCALITY_PREFIXES,
    AADR_DIR,
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

OUTPUT_DIR = AADR_DIR / "results" / "plots"
MAP_PATH = OUTPUT_DIR / "pmr_map_by_time.pdf"
REGRESSION_PATH = OUTPUT_DIR / "pmr_vs_migratory_distance.pdf"
CELLS_CSV_PATH = OUTPUT_DIR / "sites.csv"

OVERLAP_THRESHOLD = 30_000
MIN_PAIRS_PER_CELL = 5
TIME_BINS: list[tuple[float, float, str]] = [
    (10_000.0, np.inf, ">10,000 BP"),
    (5_000.0, 10_000.0, "10,000-5,000 BP"),
    (2_000.0, 5_000.0, "5,000-2,000 BP"),
    (500.0, 2_000.0, "2,000-500 BP"),
]

# Waypoints approximate overland routes around major water bodies, inspired by the
# approach of Ramachandran et al. 2005 (https://doi.org/10.1073/pnas.0507611102).
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
# fmt: off
COUNTRY_TO_REGION: dict[str, str] = {
    # ── Africa ──
    "AO": "africa", "BF": "africa", "BI": "africa", "BJ": "africa", "BW": "africa", "CD": "africa", "CF": "africa",
    "CG": "africa", "CI": "africa", "CM": "africa", "CV": "africa", "DJ": "africa", "DZ": "africa", "EG": "africa",
    "EH": "africa", "ER": "africa", "ET": "africa", "GA": "africa", "GH": "africa", "GM": "africa", "GN": "africa",
    "GQ": "africa", "GW": "africa", "KE": "africa", "KM": "africa", "LR": "africa", "LS": "africa", "LY": "africa",
    "MA": "africa", "MG": "africa", "ML": "africa", "MR": "africa", "MU": "africa", "MW": "africa", "MZ": "africa",
    "NA": "africa", "NE": "africa", "NG": "africa", "RE": "africa", "RW": "africa", "SC": "africa", "SD": "africa",
    "SH": "africa", "SL": "africa", "SN": "africa", "SO": "africa", "SS": "africa", "ST": "africa", "SZ": "africa",
    "TD": "africa", "TG": "africa", "TN": "africa", "TZ": "africa", "UG": "africa", "YT": "africa", "ZA": "africa",
    "ZM": "africa", "ZW": "africa",
    # ── Asia ──
    "AE": "asia", "AF": "asia", "AM": "asia", "AZ": "asia", "BD": "asia", "BH": "asia", "BN": "asia",
    "BT": "asia", "CN": "asia", "CY": "asia", "GE": "asia", "HK": "asia", "ID": "asia", "IL": "asia",
    "IN": "asia", "IQ": "asia", "IR": "asia", "JO": "asia", "JP": "asia", "KG": "asia", "KH": "asia",
    "KP": "asia", "KR": "asia", "KW": "asia", "KZ": "asia", "LA": "asia", "LB": "asia", "LK": "asia",
    "MM": "asia", "MN": "asia", "MO": "asia", "MV": "asia", "MY": "asia", "NP": "asia", "OM": "asia",
    "PH": "asia", "PK": "asia", "PS": "asia", "QA": "asia", "RU": "asia", "SA": "asia", "SG": "asia",
    "SY": "asia", "TH": "asia", "TJ": "asia", "TL": "asia", "TM": "asia", "TR": "asia", "TW": "asia",
    "UZ": "asia", "VN": "asia", "YE": "asia",
    # ── Europe ──
    "AD": "europe", "AL": "europe", "AT": "europe", "AX": "europe", "BA": "europe", "BE": "europe",
    "BG": "europe", "BY": "europe", "CH": "europe", "CZ": "europe", "DE": "europe", "DK": "europe",
    "EE": "europe", "ES": "europe", "FI": "europe", "FO": "europe", "FR": "europe", "GB": "europe",
    "GG": "europe", "GI": "europe", "GR": "europe", "HR": "europe", "HU": "europe", "IE": "europe",
    "IM": "europe", "IS": "europe", "IT": "europe", "JE": "europe", "LI": "europe", "LT": "europe",
    "LU": "europe", "LV": "europe", "MC": "europe", "MD": "europe", "ME": "europe", "MK": "europe",
    "MT": "europe", "NL": "europe", "NO": "europe", "PL": "europe", "PT": "europe", "RO": "europe",
    "RS": "europe", "SE": "europe", "SI": "europe", "SJ": "europe", "SK": "europe", "SM": "europe",
    "UA": "europe", "VA": "europe", "XK": "europe",
    # ── Oceania ──
    "AS": "oceania", "AU": "oceania", "CK": "oceania", "FJ": "oceania", "FM": "oceania", "GU": "oceania",
    "KI": "oceania", "MH": "oceania", "MP": "oceania", "NC": "oceania", "NF": "oceania", "NR": "oceania",
    "NU": "oceania", "NZ": "oceania", "PF": "oceania", "PG": "oceania", "PN": "oceania", "PW": "oceania",
    "SB": "oceania", "TK": "oceania", "TO": "oceania", "TV": "oceania", "VU": "oceania", "WF": "oceania",
    "WS": "oceania",
    # ── Americas ──
    "AG": "americas", "AI": "americas", "AR": "americas", "AW": "americas", "BB": "americas", "BL": "americas",
    "BM": "americas", "BO": "americas", "BQ": "americas", "BR": "americas", "BS": "americas", "BZ": "americas",
    "CA": "americas", "CL": "americas", "CO": "americas", "CR": "americas", "CU": "americas", "CW": "americas",
    "DM": "americas", "DO": "americas", "EC": "americas", "FK": "americas", "GD": "americas", "GF": "americas",
    "GL": "americas", "GP": "americas", "GT": "americas", "GY": "americas", "HN": "americas", "HT": "americas",
    "JM": "americas", "KN": "americas", "KY": "americas", "LC": "americas", "MF": "americas", "MQ": "americas",
    "MS": "americas", "MX": "americas", "NI": "americas", "PA": "americas", "PE": "americas", "PM": "americas",
    "PR": "americas", "PY": "americas", "SR": "americas", "SV": "americas", "SX": "americas", "TC": "americas",
    "TT": "americas", "US": "americas", "UY": "americas", "VC": "americas", "VE": "americas", "VG": "americas",
    "VI": "americas",
    # ── Miscellaneous territories ──
    "CC": "oceania", "CX": "oceania",  # Cocos/Christmas Islands (AU territory)
    "HM": "oceania",  # Heard/McDonald Islands (AU territory)
    "IO": "asia",  # British Indian Ocean Territory (nearest path via Cairo→India)
    "TF": "africa",  # French Southern Territories (least-wrong among five regions)
    "UM": "oceania",  # US Minor Outlying Islands
}
# fmt: on
# US Pacific territories that reverse_geocoder reports as "US" but should classify as Oceania.
US_OCEANIA_ADMIN1 = {"Guam", "Northern Mariana Islands"}


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
        if locality.startswith(AADR_ANOMALOUS_LOCALITY_PREFIXES):
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


def classify_cells(cells: list[tuple[float, float]]) -> list[str]:
    """Classify (lat, lon) cells into regions using reverse_geocoder."""
    results = rg.search(cells)
    regions: list[str] = []
    for (lat, lon), result in zip(cells, results, strict=True):
        cc = result["cc"]
        # US Pacific territories (Guam, CNMI) are returned as "US" by reverse_geocoder
        if cc == "US" and (result["admin1"] in US_OCEANIA_ADMIN1 or (lon > 144 and lat < 25)):
            regions.append("oceania")
            continue
        region = COUNTRY_TO_REGION[cc]
        regions.append(region)
    return regions


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
    """Scatter plot + linear regression of within-site PMR vs. migratory distance."""
    x = cells["migratory_distance_km"].to_numpy() / 1_000  # Thousands of km
    y = cells["median_pmr"].to_numpy()
    slope, intercept, r_value, p_value, _se = stats.linregress(x, y)
    p_exponent = int(np.floor(np.log10(p_value)))
    p_mantissa = p_value / 10**p_exponent

    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    palette = sns.color_palette("colorblind")
    colors = [palette[i] for i in (0, 1, 4, 2)]  # Blue, orange, purple, teal
    n_pairs = cells["n_pairs"].to_numpy()
    sizes = 20 + 3 * np.clip(n_pairs, 0, 20)
    for bin_i, (_low, _high, label) in enumerate(TIME_BINS):
        mask = cells["bin_idx"].to_numpy() == bin_i
        ax.scatter(
            x[mask],
            y[mask],
            label=label,
            color=colors[bin_i],
            s=sizes[mask],
            alpha=0.7,
            edgecolor="black",
            linewidth=0.3,
        )
    sns.regplot(
        x=x,
        y=y,
        ci=95,
        scatter=False,
        color="black",
        line_kws={"linewidth": 1.5, "alpha": 0.7},
        label=f"OLS ($R^2 = {r_value**2:.3f}$, $p = {p_mantissa:.1f} \\times 10^{{{p_exponent}}}$)",
        ax=ax,
    )
    ax.set_xlabel("Waypoint distance from Addis Ababa (×1000 km)", fontsize=12)
    ax.set_ylabel("Median within-site PMR", fontsize=12)
    ax.set_title('Within-site PMR vs. migratory "Out of Africa" distance', fontsize=13)
    ax.legend(fontsize=9, loc="upper right")

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
    cell_coords = list(zip(cells_sorted["lat"], cells_sorted["lon"], strict=True))
    cell_regions = classify_cells(cell_coords)
    cells_sorted["region"] = cell_regions
    cells_sorted["migratory_distance_km"] = [
        migratory_distance(lat, lon, region)
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
    print(f"Wrote PMR vs. migratory distance plot to {REGRESSION_PATH}.")


if __name__ == "__main__":
    main()
