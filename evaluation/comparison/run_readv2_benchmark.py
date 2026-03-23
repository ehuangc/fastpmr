import shutil
import subprocess
import tempfile
from pathlib import Path

from evaluation_utils import (
    COMPARISON_DATA_PREFIX,
    PERFORMANCE_RUNS,
    PLINK_EXTS,
    ensure_data_present,
    quote_path,
    run_benchmark,
)

READV2_REPO = "https://github.com/GuntherLab/READv2"
READV2_DIR = Path(__file__).resolve().parent / "READv2"
READV2_SCRIPT = READV2_DIR / "READ2.py"
PLINK_PREFIX = COMPARISON_DATA_PREFIX
RESULTS_DIR = Path(__file__).resolve().parent / "results"


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
        "fastpmr",
        # Ensure we test fastpmr on PLINK dataset, not EIGENSTRAT
        f"--prefix {quote_path(prefix.with_suffix('.bed'))}",
        f"--output-directory {quote_path(output_dir)}",
        "--min-covered-snps 0",
        "-n",
    ]
    return " ".join(parts)


def build_readv2_command(readv2_script: Path, prefix: Path, work_dir: Path) -> str:
    return f"cd {quote_path(work_dir)} && python {quote_path(readv2_script)} -i {quote_path(prefix)}"


def main() -> None:
    data_prefix = Path(COMPARISON_DATA_PREFIX)
    ensure_data_present(data_prefix, PLINK_EXTS)
    clone_readv2(READV2_DIR)

    export_path = RESULTS_DIR / "readv2_comparison_benchmark.csv"

    output_dir = tempfile.mkdtemp()
    fastpmr_cmd = build_fastpmr_command(PLINK_PREFIX, Path(output_dir))
    readv2_cmd = build_readv2_command(READV2_SCRIPT, PLINK_PREFIX, Path(output_dir))
    configs = [("fastpmr", fastpmr_cmd), ("READv2", readv2_cmd)]
    run_benchmark(configs, export_path, runs=PERFORMANCE_RUNS)


if __name__ == "__main__":
    main()
