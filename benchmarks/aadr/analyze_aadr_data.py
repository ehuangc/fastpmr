import csv
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

SCRIPT_DIR = Path(__file__).resolve().parent
COUNTS_PATH = SCRIPT_DIR / "results" / "fastpmr" / "mismatch_counts.npz"
METADATA_PATH = SCRIPT_DIR / "data" / "v62.0_1240k_public.anno"
OUTPUT_DIR = SCRIPT_DIR / "results" / "analysis"
SAME_MASTER_OUTPUT_CSV = OUTPUT_DIR / "same_master_id_high_pmr.csv"
DIFF_MASTER_OUTPUT_CSV = OUTPUT_DIR / "diff_master_id_low_pmr.csv"
IDENTICAL_PMR_THRESHOLD = 0.14
NON_IDENTICAL_PMR_THRESHOLD = 0.17
OVERLAP_THRESHOLD = 30000


def ensure_data_present(counts_path: Path, metadata_path: Path) -> None:
    missing = [path for path in (counts_path, metadata_path) if not path.is_file()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        raise SystemExit(
            f"Missing input files: {missing_str}. Run `pixi run fastpmr-aadr` after `pixi run prepare-aadr-data`."
        )


def load_samples(npz_path: Path) -> list[str]:
    with zipfile.ZipFile(npz_path) as archive:
        with archive.open("samples.json") as handle:
            return json.loads(handle.read().decode("utf-8"))


def load_metadata(metadata_path: Path) -> dict[str, dict[str, str]]:
    with metadata_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        metadata = {}
        for row in reader:
            if not row:
                continue
            metadata[row[0]] = dict(zip(header, row))
    return metadata


def load_mismatch_counts(
    npz_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    with np.load(npz_path, allow_pickle=False) as npz:
        mismatches = npz["mismatches"]
        totals = npz["totals"]
        site_overlaps = npz["site_overlaps"]
    samples = load_samples(npz_path)
    return mismatches, totals, site_overlaps, samples


def compute_mismatch_rates(mismatches: np.ndarray, totals: np.ndarray) -> np.ndarray:
    rates = np.full(mismatches.shape, np.nan, dtype=np.float64)
    np.divide(mismatches, totals, out=rates, where=totals != 0)
    return rates


def compute_mismatch_rate_cis(
    mismatches: np.ndarray, totals: np.ndarray, confidence_level: float = 0.95
) -> tuple[np.ndarray, np.ndarray]:
    # Assume pseudohaploid data, so 1 potential mismatch per site instead of 2
    # Note that this gives us conservative estimates for diploid data
    pseudo_mismatches = mismatches.astype(np.float64) / 2.0
    pseudo_totals = totals.astype(np.float64) / 2.0
    lower = np.full(pseudo_mismatches.shape, np.nan, dtype=np.float64)
    upper = np.full(pseudo_mismatches.shape, np.nan, dtype=np.float64)

    alpha = 1.0 - confidence_level
    valid = pseudo_totals > 0

    zero = valid & (pseudo_mismatches <= 0)
    non_zero = valid & (pseudo_mismatches > 0)
    lower[zero] = 0.0
    lower[non_zero] = beta.ppf(
        alpha / 2, pseudo_mismatches[non_zero], pseudo_totals[non_zero] - pseudo_mismatches[non_zero] + 1
    )

    one = valid & (pseudo_mismatches >= pseudo_totals)
    non_one = valid & (pseudo_mismatches < pseudo_totals)
    upper[one] = 1.0
    upper[non_one] = beta.ppf(
        1 - alpha / 2, pseudo_mismatches[non_one] + 1, pseudo_totals[non_one] - pseudo_mismatches[non_one]
    )
    return lower, upper


def filter_samples(
    mismatches: np.ndarray,
    totals: np.ndarray,
    site_overlaps: np.ndarray,
    samples: list[str],
    metadata: dict[str, dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    keep_indices: list[int] = []
    for idx, sample in enumerate(samples):
        # Exclude archaic humans
        sample_group = metadata[sample]["Group ID"]
        if "neanderthal" in sample_group.lower() or "denisova" in sample_group.lower():
            continue
        # Exclude references
        if ".ref" in sample.lower():
            continue
        keep_indices.append(idx)

    keep = np.array(keep_indices, dtype=int)
    return (
        mismatches[np.ix_(keep, keep)],
        totals[np.ix_(keep, keep)],
        site_overlaps[np.ix_(keep, keep)],
        [samples[idx] for idx in keep],
    )


def get_pair_metadata(metadata: dict[str, dict[str, str]], sample1: str, sample2: str) -> dict[str, str]:
    date_field = (
        "Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age "
        "(Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). "
        "(Format 2) Archaeological context range, e.g. 2500-1700 BCE"
    )
    return {
        "publication_abbreviation1": metadata[sample1]["Publication abbreviation"],
        "publication_abbreviation2": metadata[sample2]["Publication abbreviation"],
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
    }


def find_same_master_id_high_pmr_pairs(
    metadata: dict[str, dict[str, str]],
    samples: list[str],
    master_ids: list[str],
    site_overlaps: np.ndarray,
    rates: np.ndarray,
    rates_ci_lower: np.ndarray,
    rates_ci_upper: np.ndarray,
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
                pmr = rates[idx_i, idx_j]
                lower = rates_ci_lower[idx_i, idx_j]
                upper = rates_ci_upper[idx_i, idx_j]
                if np.isnan(pmr):
                    continue
                n_overlap = int(site_overlaps[idx_i, idx_j])
                if n_overlap <= overlap_threshold:
                    continue
                if pmr > threshold:
                    row = {
                        "master_id": master_id,
                        "genetic_id1": samples[idx_i],
                        "genetic_id2": samples[idx_j],
                        "site_overlap": n_overlap,
                        "mismatch_rate": pmr,
                        "mismatch_rate_95_ci_lower": lower,
                        "mismatch_rate_95_ci_upper": upper,
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
    rates: np.ndarray,
    rates_ci_lower: np.ndarray,
    rates_ci_upper: np.ndarray,
    threshold: float,
    overlap_threshold: int,
) -> pd.DataFrame:
    master_ids_array = np.asarray(master_ids, dtype=object)
    rows: list[dict[str, str | int | float]] = []
    for idx_i in range(len(samples) - 1):
        master_id_i = master_ids_array[idx_i]
        if not master_id_i:
            continue
        row_rates = rates[idx_i, idx_i + 1 :]
        row_overlaps = site_overlaps[idx_i, idx_i + 1 :]
        row_master_ids = master_ids_array[idx_i + 1 :]
        mask = (
            (row_master_ids != "")
            & (row_master_ids != master_id_i)
            & (row_rates < threshold)
            & (row_overlaps > overlap_threshold)
        )
        mask &= ~np.isnan(row_rates)
        match_offsets = np.nonzero(mask)[0]
        for offset in match_offsets:
            idx_j = idx_i + 1 + int(offset)
            lower = rates_ci_lower[idx_i, idx_j]
            upper = rates_ci_upper[idx_i, idx_j]
            row = {
                "master_id1": master_id_i,
                "master_id2": row_master_ids[offset],
                "genetic_id1": samples[idx_i],
                "genetic_id2": samples[idx_j],
                "site_overlap": row_overlaps[offset],
                "mismatch_rate": row_rates[offset],
                "mismatch_rate_95_ci_lower": lower,
                "mismatch_rate_95_ci_upper": upper,
            }
            row.update(get_pair_metadata(metadata, samples[idx_i], samples[idx_j]))
            rows.append(row)
    rows.sort(key=lambda row: row["mismatch_rate"])
    return pd.DataFrame.from_records(rows)


def main() -> None:
    ensure_data_present(COUNTS_PATH, METADATA_PATH)
    mismatches, totals, site_overlaps, samples = load_mismatch_counts(COUNTS_PATH)
    metadata = load_metadata(METADATA_PATH)
    matched = sum(sample in metadata for sample in samples)
    assert len(samples) == len(metadata) == matched == mismatches.shape[0] == totals.shape[0] == site_overlaps.shape[0]

    filtered_mismatches, filtered_totals, filtered_site_overlaps, filtered_samples = filter_samples(
        mismatches, totals, site_overlaps, samples, metadata
    )
    filtered_master_ids = [metadata[sample]["Master ID"] for sample in filtered_samples]
    filtered_rates = compute_mismatch_rates(filtered_mismatches, filtered_totals)
    filtered_rates_ci_lower, filtered_rates_ci_upper = compute_mismatch_rate_cis(filtered_mismatches, filtered_totals)

    high_pmr_pairs = find_same_master_id_high_pmr_pairs(
        metadata,
        filtered_samples,
        filtered_master_ids,
        filtered_site_overlaps,
        filtered_rates,
        filtered_rates_ci_lower,
        filtered_rates_ci_upper,
        NON_IDENTICAL_PMR_THRESHOLD,
        OVERLAP_THRESHOLD,
    )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    high_pmr_pairs.to_csv(SAME_MASTER_OUTPUT_CSV, index=False)
    print(
        f"Wrote {len(high_pmr_pairs)} same-master-ID pairs with PMR > {NON_IDENTICAL_PMR_THRESHOLD} "
        f"and site overlap > {OVERLAP_THRESHOLD} to {SAME_MASTER_OUTPUT_CSV}."
    )

    low_pmr_pairs = find_diff_master_id_low_pmr_pairs(
        metadata,
        filtered_samples,
        filtered_master_ids,
        filtered_site_overlaps,
        filtered_rates,
        filtered_rates_ci_lower,
        filtered_rates_ci_upper,
        IDENTICAL_PMR_THRESHOLD,
        OVERLAP_THRESHOLD,
    )
    low_pmr_pairs.to_csv(DIFF_MASTER_OUTPUT_CSV, index=False)
    print(
        f"Wrote {len(low_pmr_pairs)} different-master-ID pairs with PMR < {IDENTICAL_PMR_THRESHOLD} "
        f"and site overlap > {OVERLAP_THRESHOLD} to {DIFF_MASTER_OUTPUT_CSV}."
    )


if __name__ == "__main__":
    main()
