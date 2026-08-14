import csv
import json
import shutil
import subprocess
import zipfile
from pathlib import Path

import numpy as np

from evaluation_utils.constants import (
    COMPARISON_DATA_PREFIX,
    FASTPMR_BIN,
    PERFORMANCE_RUNS,
    PLINK_EXTS,
)
from evaluation_utils.core import (
    ensure_data_present,
    quote_path,
    run_benchmark,
)

READV2_REPO = "https://github.com/GuntherLab/READv2"
READV2_DIR = Path(__file__).resolve().parent / "READv2"
READV2_SCRIPT = READV2_DIR / "READ2.py"
PLINK_PREFIX = COMPARISON_DATA_PREFIX
RESULTS_DIR = Path(__file__).resolve().parent / "results"
OUTPUTS_DIR = RESULTS_DIR / "outputs"
PMR_TOLERANCE = 1e-6


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


def check_pmrs_match(fastpmr_output_dir: Path, readv2_output_dir: Path) -> None:
    """Check that fastpmr and READv2 report the same PMR for every pair."""
    npz_path = fastpmr_output_dir / "fastpmr_results.npz"
    with np.load(npz_path, allow_pickle=False) as npz:
        mismatch_rates = npz["mismatch_rates"]
    with zipfile.ZipFile(npz_path) as archive:
        samples = json.loads(archive.read("samples.json").decode("utf-8"))
    # fastpmr labels samples as "<family ID>:<individual ID>" while READv2 uses the individual ID alone
    sample_indices = {sample.split(":", 1)[-1]: i for i, sample in enumerate(samples)}

    # READv2 drops pairs without overlapping SNPs, so only the pairs it reports are compared
    with (readv2_output_dir / "Read_Results.tsv").open(newline="", encoding="utf-8") as handle:
        readv2_rows = list(csv.DictReader(handle, delimiter="\t"))
    for readv2_row in readv2_rows:
        first, second = readv2_row["PairIndividuals"].split(",")
        readv2_pmr = float(readv2_row["Nonnormalized_P0"])
        fastpmr_pmr = mismatch_rates[sample_indices[first], sample_indices[second]]
        if abs(fastpmr_pmr - readv2_pmr) > PMR_TOLERANCE:
            raise SystemExit(f"PMR mismatch for {first},{second}: fastpmr={fastpmr_pmr:.10g}, READv2={readv2_pmr:.10g}")
    print(f"\nAll {len(readv2_rows)} pairwise PMRs match between fastpmr and READv2 to within {PMR_TOLERANCE:g}")


def main() -> None:
    data_prefix = Path(COMPARISON_DATA_PREFIX)
    ensure_data_present(data_prefix, PLINK_EXTS)
    clone_readv2(READV2_DIR)

    fastpmr_output_dir = OUTPUTS_DIR / "fastpmr"
    readv2_output_dir = OUTPUTS_DIR / "readv2"
    fastpmr_output_dir.mkdir(parents=True, exist_ok=True)
    readv2_output_dir.mkdir(parents=True, exist_ok=True)

    export_path = RESULTS_DIR / "readv2_comparison_benchmark.csv"

    fastpmr_cmd = build_fastpmr_command(PLINK_PREFIX, fastpmr_output_dir)
    readv2_cmd = build_readv2_command(READV2_SCRIPT, PLINK_PREFIX, readv2_output_dir)
    configs = [("fastpmr", fastpmr_cmd), ("READv2", readv2_cmd)]
    run_benchmark(configs, export_path, runs=PERFORMANCE_RUNS)

    check_pmrs_match(fastpmr_output_dir, readv2_output_dir)


if __name__ == "__main__":
    main()
