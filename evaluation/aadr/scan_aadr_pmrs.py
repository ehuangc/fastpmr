from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from haversine import haversine

from evaluation_utils import (
    AADR_DIR,
    AADR_EXCLUDED_LOCALITY_PREFIXES,
    AADR_METADATA_PATH,
    AADR_NPZ_PATH,
    ensure_aadr_npz_present,
    is_archaic_or_reference_sample,
    load_aadr_metadata,
    load_aadr_npz_arrays,
)

OUTPUT_DIR = AADR_DIR / "results" / "scans"
SAME_MASTER_OUTPUT_CSV = OUTPUT_DIR / "potential_spurious_duplicates.csv"
DIFF_MASTER_OUTPUT_CSV = OUTPUT_DIR / "potential_missed_duplicates.csv"
DIFF_LOCALITY_OUTPUT_CSV = OUTPUT_DIR / "potential_cross_site_relatives.csv"
SAME_MASTER_HISTOGRAM_PATH = OUTPUT_DIR / "same_master_id_pmrs.pdf"
DIFF_MASTER_HISTOGRAM_PATH = OUTPUT_DIR / "diff_master_id_pmrs.pdf"

IDENTICAL_PMR_THRESHOLD = 0.14
NON_IDENTICAL_PMR_THRESHOLD = 0.17
FIRST_DEGREE_PMR_THRESHOLD = 0.18
OVERLAP_THRESHOLD = 30000

EURASIA_LAT_RANGE = (0.0, 85.0)
EURASIA_LON_RANGE = (-30.0, 180.0)


def filter_samples(
    samples: list[str],
    site_overlaps: np.ndarray,
    mismatch_rates: np.ndarray,
    mismatch_rates_95_ci_lower: np.ndarray,
    mismatch_rates_95_ci_upper: np.ndarray,
    metadata: dict[str, dict[str, str]],
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    keep_indices: list[int] = []
    for idx, sample in enumerate(samples):
        if is_archaic_or_reference_sample(sample, metadata[sample]):
            continue
        try:
            lat = float(metadata[sample]["Lat."])
            lon = float(metadata[sample]["Long."])
        except (TypeError, ValueError):
            continue
        if not (
            EURASIA_LAT_RANGE[0] <= lat <= EURASIA_LAT_RANGE[1] and EURASIA_LON_RANGE[0] <= lon <= EURASIA_LON_RANGE[1]
        ):
            continue
        keep_indices.append(idx)

    keep = np.array(keep_indices, dtype=int)
    return (
        [samples[idx] for idx in keep],
        site_overlaps[np.ix_(keep, keep)],
        mismatch_rates[np.ix_(keep, keep)],
        mismatch_rates_95_ci_lower[np.ix_(keep, keep)],
        mismatch_rates_95_ci_upper[np.ix_(keep, keep)],
    )


def get_pair_metadata(metadata: dict[str, dict[str, str]], sample1: str, sample2: str) -> dict[str, str]:
    date_field = (
        "Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age "
        "(Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). "
        "(Format 2) Archaeological context range, e.g. 2500-1700 BCE"
    )
    return {
        "publication1": metadata[sample1]["Publication abbreviation"],
        "publication2": metadata[sample2]["Publication abbreviation"],
        "skeletal_code1": metadata[sample1]["Skeletal code"],
        "skeletal_code2": metadata[sample2]["Skeletal code"],
        "date1": metadata[sample1][date_field],
        "date2": metadata[sample2][date_field],
        "group_id1": metadata[sample1]["Group ID"],
        "group_id2": metadata[sample2]["Group ID"],
        "locality1": metadata[sample1]["Locality"],
        "locality2": metadata[sample2]["Locality"],
        "political_entity1": metadata[sample1]["Political Entity"],
        "political_entity2": metadata[sample2]["Political Entity"],
        "lat1": metadata[sample1]["Lat."],
        "lon1": metadata[sample1]["Long."],
        "lat2": metadata[sample2]["Lat."],
        "lon2": metadata[sample2]["Long."],
    }


def find_same_master_id_high_pmr_pairs(
    metadata: dict[str, dict[str, str]],
    samples: list[str],
    master_ids: list[str],
    site_overlaps: np.ndarray,
    mismatch_rates: np.ndarray,
    mismatch_rates_95_ci_lower: np.ndarray,
    mismatch_rates_95_ci_upper: np.ndarray,
    threshold: float,
    overlap_threshold: int,
) -> pd.DataFrame:
    master_id_to_indices: dict[str, list[int]] = {}
    for idx, master_id in enumerate(master_ids):
        if not master_id:
            continue
        master_id_to_indices.setdefault(master_id, []).append(idx)

    rows: list[dict[str, str | int | float]] = []
    for master_id, indices in master_id_to_indices.items():
        if len(indices) < 2:
            continue
        for offset, idx_i in enumerate(indices[:-1]):
            for idx_j in indices[offset + 1 :]:
                pmr = mismatch_rates[idx_i, idx_j]
                if np.isnan(pmr):
                    continue
                n_overlap = int(site_overlaps[idx_i, idx_j])
                if n_overlap < overlap_threshold:
                    continue
                if pmr > threshold:
                    row = {
                        "master_id": master_id,
                        "genetic_id1": samples[idx_i],
                        "genetic_id2": samples[idx_j],
                        "site_overlap": n_overlap,
                        "mismatch_rate": pmr,
                        "mismatch_rate_95_ci_lower": mismatch_rates_95_ci_lower[idx_i, idx_j],
                        "mismatch_rate_95_ci_upper": mismatch_rates_95_ci_upper[idx_i, idx_j],
                    }
                    row.update(get_pair_metadata(metadata, samples[idx_i], samples[idx_j]))
                    rows.append(row)
    rows.sort(key=lambda row: row["mismatch_rate"], reverse=True)
    return pd.DataFrame.from_records(rows)


def find_diff_master_id_low_pmr_pairs(
    metadata: dict[str, dict[str, str]],
    samples: list[str],
    master_ids: list[str],
    site_overlaps: np.ndarray,
    mismatch_rates: np.ndarray,
    mismatch_rates_95_ci_lower: np.ndarray,
    mismatch_rates_95_ci_upper: np.ndarray,
    threshold: float,
    overlap_threshold: int,
) -> pd.DataFrame:
    master_ids_array = np.asarray(master_ids, dtype=object)
    rows: list[dict[str, str | int | float]] = []
    for idx_i in range(len(samples) - 1):
        master_id_i = master_ids_array[idx_i]
        if not master_id_i:
            continue
        row_rates = mismatch_rates[idx_i, idx_i + 1 :]
        row_overlaps = site_overlaps[idx_i, idx_i + 1 :]
        row_master_ids = master_ids_array[idx_i + 1 :]
        mask = (
            (row_master_ids != "")
            & (row_master_ids != master_id_i)
            & (row_rates < threshold)
            & (row_overlaps >= overlap_threshold)
        )
        mask &= ~np.isnan(row_rates)
        match_offsets = np.nonzero(mask)[0]
        for offset in match_offsets:
            idx_j = idx_i + 1 + int(offset)
            row = {
                "master_id1": master_id_i,
                "master_id2": row_master_ids[offset],
                "genetic_id1": samples[idx_i],
                "genetic_id2": samples[idx_j],
                "site_overlap": int(row_overlaps[offset]),
                "mismatch_rate": row_rates[offset],
                "mismatch_rate_95_ci_lower": mismatch_rates_95_ci_lower[idx_i, idx_j],
                "mismatch_rate_95_ci_upper": mismatch_rates_95_ci_upper[idx_i, idx_j],
            }
            row.update(get_pair_metadata(metadata, samples[idx_i], samples[idx_j]))
            rows.append(row)
    rows.sort(key=lambda row: row["mismatch_rate"])
    return pd.DataFrame.from_records(rows)


def find_diff_locality_low_pmr_pairs(
    metadata: dict[str, dict[str, str]],
    samples: list[str],
    master_ids: list[str],
    localities: list[str],
    site_overlaps: np.ndarray,
    mismatch_rates: np.ndarray,
    mismatch_rates_95_ci_lower: np.ndarray,
    mismatch_rates_95_ci_upper: np.ndarray,
    lower_threshold: float,
    upper_threshold: float,
    overlap_threshold: int,
) -> pd.DataFrame:
    master_ids_array = np.asarray(master_ids, dtype=object)
    localities_array = np.asarray(localities, dtype=object)
    rows: list[dict[str, str | int | float]] = []
    for idx_i in range(len(samples) - 1):
        locality_i = localities_array[idx_i]
        master_id_i = master_ids_array[idx_i]
        if not locality_i or not master_id_i:
            continue
        row_rates = mismatch_rates[idx_i, idx_i + 1 :]
        row_overlaps = site_overlaps[idx_i, idx_i + 1 :]
        row_localities = localities_array[idx_i + 1 :]
        row_master_ids = master_ids_array[idx_i + 1 :]
        mask = (
            (row_localities != "")
            & (row_localities != locality_i)
            & (row_master_ids != "")
            & (row_master_ids != master_id_i)
            & (row_rates >= lower_threshold)
            & (row_rates < upper_threshold)
            & (row_overlaps >= overlap_threshold)
        )
        mask &= ~np.isnan(row_rates)
        match_offsets = np.nonzero(mask)[0]
        for offset in match_offsets:
            idx_j = idx_i + 1 + int(offset)
            locality_j = row_localities[offset]
            if locality_i.startswith(AADR_EXCLUDED_LOCALITY_PREFIXES) or locality_j.startswith(
                AADR_EXCLUDED_LOCALITY_PREFIXES
            ):
                continue

            row = {
                "master_id1": master_id_i,
                "master_id2": row_master_ids[offset],
                "genetic_id1": samples[idx_i],
                "genetic_id2": samples[idx_j],
                "site_overlap": int(row_overlaps[offset]),
                "mismatch_rate": row_rates[offset],
                "mismatch_rate_95_ci_lower": mismatch_rates_95_ci_lower[idx_i, idx_j],
                "mismatch_rate_95_ci_upper": mismatch_rates_95_ci_upper[idx_i, idx_j],
            }
            row.update(get_pair_metadata(metadata, samples[idx_i], samples[idx_j]))
            row["distance_km"] = round(
                haversine(
                    (float(row["lat1"]), float(row["lon1"])),
                    (float(row["lat2"]), float(row["lon2"])),
                    unit="km",
                )
            )
            rows.append(row)
    rows.sort(key=lambda row: row["distance_km"], reverse=True)
    return pd.DataFrame.from_records(rows)


def collect_pairwise_mismatch_rates(
    master_ids: list[str],
    rates: np.ndarray,
    site_overlaps: np.ndarray,
    overlap_threshold: int,
) -> tuple[np.ndarray, np.ndarray]:
    master_ids_array = np.asarray(master_ids, dtype=object)
    pair_i, pair_j = np.triu_indices(len(master_ids_array), k=1)
    pair_rates = rates[pair_i, pair_j]
    pair_site_overlaps = site_overlaps[pair_i, pair_j]
    master_i = master_ids_array[pair_i]
    master_j = master_ids_array[pair_j]

    valid = (pair_site_overlaps >= overlap_threshold) & ~np.isnan(pair_rates) & (master_i != "") & (master_j != "")
    same_mask = valid & (master_i == master_j)
    diff_mask = valid & (master_i != master_j)
    return pair_rates[same_mask], pair_rates[diff_mask]


def plot_pairwise_mismatch_rate_histograms(
    same_rates: np.ndarray,
    diff_rates: np.ndarray,
    same_output_path: Path,
    diff_output_path: Path,
) -> None:
    Y_AXIS_LIMIT = 100

    def plot_histogram(rates: np.ndarray, title: str, color: str, output_path: Path) -> None:
        bins = np.histogram_bin_edges(rates, bins=50)
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.hist(rates, bins=bins, color=color, alpha=0.85, edgecolor="black", linewidth=0.5)
        ax.set_title(title)
        ax.set_xlabel("Pairwise mismatch rate")
        ax.set_ylabel("Count")
        ax.set_ylim(0, Y_AXIS_LIMIT)
        fig.tight_layout()
        fig.savefig(output_path, dpi=600)
        plt.close(fig)

    plot_histogram(
        same_rates, f"AADR Sample Pairs, Matching Master IDs\n(n={len(same_rates)})", "#4C78A8", same_output_path
    )
    plot_histogram(
        diff_rates,
        f"AADR Sample Pairs, Non-Matching Master IDs\n(n={len(diff_rates)}, y-axis limit={Y_AXIS_LIMIT})",
        "#F58518",
        diff_output_path,
    )


def main() -> None:
    ensure_aadr_npz_present(AADR_NPZ_PATH, AADR_METADATA_PATH)
    samples, site_overlaps, mismatch_rates, mismatch_rates_95_ci_lower, mismatch_rates_95_ci_upper, _covered_snps = (
        load_aadr_npz_arrays(AADR_NPZ_PATH)
    )
    metadata = load_aadr_metadata(AADR_METADATA_PATH)
    matched = sum(sample in metadata for sample in samples)
    assert len(samples) == len(metadata) == matched == mismatch_rates.shape[0] == site_overlaps.shape[0]

    (
        filtered_samples,
        filtered_site_overlaps,
        filtered_mismatch_rates,
        filtered_mismatch_rates_95_ci_lower,
        filtered_mismatch_rates_95_ci_upper,
    ) = filter_samples(
        samples, site_overlaps, mismatch_rates, mismatch_rates_95_ci_lower, mismatch_rates_95_ci_upper, metadata
    )
    filtered_master_ids = [metadata[sample]["Master ID"] for sample in filtered_samples]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    same_rates, diff_rates = collect_pairwise_mismatch_rates(
        filtered_master_ids, filtered_mismatch_rates, filtered_site_overlaps, OVERLAP_THRESHOLD
    )
    plot_pairwise_mismatch_rate_histograms(
        same_rates, diff_rates, SAME_MASTER_HISTOGRAM_PATH, DIFF_MASTER_HISTOGRAM_PATH
    )
    print(
        f"Wrote pairwise mismatch rate histograms to {SAME_MASTER_HISTOGRAM_PATH} and {DIFF_MASTER_HISTOGRAM_PATH}.\n"
    )

    high_pmr_pairs = find_same_master_id_high_pmr_pairs(
        metadata,
        filtered_samples,
        filtered_master_ids,
        filtered_site_overlaps,
        filtered_mismatch_rates,
        filtered_mismatch_rates_95_ci_lower,
        filtered_mismatch_rates_95_ci_upper,
        NON_IDENTICAL_PMR_THRESHOLD,
        OVERLAP_THRESHOLD,
    )
    rate_cols = ["mismatch_rate", "mismatch_rate_95_ci_lower", "mismatch_rate_95_ci_upper"]
    high_pmr_pairs[rate_cols] = high_pmr_pairs[rate_cols].round(6)
    high_pmr_pairs.to_csv(SAME_MASTER_OUTPUT_CSV, index=False)
    print(
        f"Wrote {len(high_pmr_pairs)} same-master-ID pairs with PMR > {NON_IDENTICAL_PMR_THRESHOLD} "
        f"and site overlap >= {OVERLAP_THRESHOLD} to {SAME_MASTER_OUTPUT_CSV}.\n"
    )

    low_pmr_pairs = find_diff_master_id_low_pmr_pairs(
        metadata,
        filtered_samples,
        filtered_master_ids,
        filtered_site_overlaps,
        filtered_mismatch_rates,
        filtered_mismatch_rates_95_ci_lower,
        filtered_mismatch_rates_95_ci_upper,
        IDENTICAL_PMR_THRESHOLD,
        OVERLAP_THRESHOLD,
    )
    low_pmr_pairs[rate_cols] = low_pmr_pairs[rate_cols].round(6)
    low_pmr_pairs.to_csv(DIFF_MASTER_OUTPUT_CSV, index=False)
    print(
        f"Wrote {len(low_pmr_pairs)} different-master-ID pairs with PMR < {IDENTICAL_PMR_THRESHOLD} "
        f"and site overlap >= {OVERLAP_THRESHOLD} to {DIFF_MASTER_OUTPUT_CSV}.\n"
    )

    filtered_localities = [metadata[sample]["Locality"] for sample in filtered_samples]
    diff_locality_pairs = find_diff_locality_low_pmr_pairs(
        metadata,
        filtered_samples,
        filtered_master_ids,
        filtered_localities,
        filtered_site_overlaps,
        filtered_mismatch_rates,
        filtered_mismatch_rates_95_ci_lower,
        filtered_mismatch_rates_95_ci_upper,
        IDENTICAL_PMR_THRESHOLD,
        FIRST_DEGREE_PMR_THRESHOLD,
        OVERLAP_THRESHOLD,
    )
    diff_locality_pairs[rate_cols] = diff_locality_pairs[rate_cols].round(6)
    diff_locality_pairs.to_csv(DIFF_LOCALITY_OUTPUT_CSV, index=False)
    print(
        f"Wrote {len(diff_locality_pairs)} different-locality pairs with "
        f"{IDENTICAL_PMR_THRESHOLD} <= PMR < {FIRST_DEGREE_PMR_THRESHOLD} "
        f"and site overlap >= {OVERLAP_THRESHOLD} to {DIFF_LOCALITY_OUTPUT_CSV}."
    )


if __name__ == "__main__":
    main()
