import csv
import json
import shutil
import subprocess
import zipfile
from pathlib import Path

import numpy as np
from prepare_data import subset_prefix

from evaluation_utils.constants import (
    COMPARISON_RUNS,
    COMPARISON_SAMPLE_SET_SIZES,
    FASTPMR_BIN,
    PLINK_EXTS,
)
from evaluation_utils.core import (
    ensure_data_present,
    quote_path,
    remove_readv2_intermediates,
    run_benchmark,
)

READV2_REPO = "https://github.com/GuntherLab/READv2"
READV2_DIR = Path(__file__).resolve().parent / "READv2"
READV2_SCRIPT = READV2_DIR / "READ2.py"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
OUTPUTS_DIR = RESULTS_DIR / "outputs"
PMR_TOLERANCE = 1e-6
# READv2 relatedness classifications and their fastpmr equivalents
READV2_DEGREES = {
    "IdenticalTwins/SameIndividual": "Identical/Twin",
    "First Degree": "First Degree",
    "Second Degree": "Second Degree",
    "Third Degree": "Third Degree",
    "Unrelated/Consistent with Third Degree": "Unrelated",
    "Unrelated": "Unrelated",
}


def clone_readv2(readv2_dir: Path) -> None:
    if readv2_dir.exists():
        print(f"Removing existing {readv2_dir}...")
        shutil.rmtree(readv2_dir)

    readv2_dir.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["git", "clone", "--depth", "1", READV2_REPO, str(readv2_dir)],
        check=True,
    )


def build_fastpmr_command(prefix: Path, output_dir: Path) -> str:
    parts = [
        FASTPMR_BIN,
        # Ensure we test fastpmr on PLINK dataset, not EIGENSTRAT
        f"--prefix {quote_path(prefix.with_suffix('.bed'))}",
        f"--output-directory {quote_path(output_dir)}",
        "--chromosomes 1-22",
        "--min-covered-snps 0",
        "--ci",
        "--degrees",
        "--npz",
    ]
    return " ".join(parts)


def build_readv2_command(readv2_script: Path, prefix: Path, work_dir: Path) -> str:
    return f"cd {quote_path(work_dir)} && python {quote_path(readv2_script)} -i {quote_path(prefix)}"


def check_results_match(fastpmr_output_dir: Path, readv2_output_dir: Path) -> None:
    """Check that fastpmr and READv2 report the same PMR and degree of relatedness for every pair."""
    npz_path = fastpmr_output_dir / "fastpmr_results.npz"
    with np.load(npz_path, allow_pickle=False) as npz:
        mismatch_rates = npz["mismatch_rates"]
        degrees = npz["degrees"]
    with zipfile.ZipFile(npz_path) as archive:
        samples = json.loads(archive.read("samples.json").decode("utf-8"))
        degree_labels = json.loads(archive.read("degree_labels.json").decode("utf-8"))
    # fastpmr labels samples as "<family ID>:<individual ID>" while READv2 uses the individual ID alone
    sample_indices = {sample.split(":", 1)[-1]: i for i, sample in enumerate(samples)}

    # READv2 drops pairs without overlapping SNPs, so only the pairs it reports are compared
    with (readv2_output_dir / "Read_Results.tsv").open(newline="", encoding="utf-8") as handle:
        readv2_rows = list(csv.DictReader(handle, delimiter="\t"))
    for readv2_row in readv2_rows:
        first, second = readv2_row["PairIndividuals"].split(",")
        pair_indices = sample_indices[first], sample_indices[second]
        readv2_pmr = float(readv2_row["Nonnormalized_P0"])
        fastpmr_pmr = mismatch_rates[pair_indices]
        if abs(fastpmr_pmr - readv2_pmr) > PMR_TOLERANCE:
            raise SystemExit(f"PMR mismatch for {first},{second}: fastpmr={fastpmr_pmr:.10g}, READv2={readv2_pmr:.10g}")
        readv2_degree = READV2_DEGREES[readv2_row["Rel"]]
        fastpmr_degree = degree_labels[degrees[pair_indices]]
        if fastpmr_degree != readv2_degree:
            raise SystemExit(
                f"Degree mismatch for {first},{second}: fastpmr={fastpmr_degree}, READv2={readv2_row['Rel']}"
            )
    print(f"\nAll {len(readv2_rows)} pairwise PMRs match between fastpmr and READv2 to within {PMR_TOLERANCE:g}")
    print(f"All {len(readv2_rows)} pairwise degrees of relatedness match between fastpmr and READv2")


def main() -> None:
    for size in COMPARISON_SAMPLE_SET_SIZES:
        ensure_data_present(
            subset_prefix(size),
            PLINK_EXTS,
            command="pixi run prepare-comparison",
            dataset="Maravall-López et al. dataset",
        )
    clone_readv2(READV2_DIR)

    fastpmr_output_dir = OUTPUTS_DIR / "fastpmr"
    readv2_output_dir = OUTPUTS_DIR / "readv2"
    fastpmr_output_dir.mkdir(parents=True, exist_ok=True)
    readv2_output_dir.mkdir(parents=True, exist_ok=True)

    export_path = RESULTS_DIR / "comparison_benchmark.csv"
    # Both tools write to the same output directories for every sample count (each run overwriting
    # the last), so only the outputs of the final and largest subset, the full dataset, survive for
    # check_results_match
    configs = []
    for size in COMPARISON_SAMPLE_SET_SIZES:
        prefix = subset_prefix(size)
        configs.append((f"tool=fastpmr_samples={size}", build_fastpmr_command(prefix, fastpmr_output_dir)))
        configs.append((f"tool=READv2_samples={size}", build_readv2_command(READV2_SCRIPT, prefix, readv2_output_dir)))
    run_benchmark(configs, export_path, runs=COMPARISON_RUNS)
    for size in COMPARISON_SAMPLE_SET_SIZES:
        remove_readv2_intermediates(subset_prefix(size))

    check_results_match(fastpmr_output_dir, readv2_output_dir)


if __name__ == "__main__":
    main()
