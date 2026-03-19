import subprocess
import tempfile
from pathlib import Path

from benchmark_utils import (
    DATA_PREFIX,
    PLINK_EXTS,
    SCRIPT_DIR,
    ensure_data_present,
    quote_path,
    run_benchmark,
)

READV2_REPO = "https://github.com/GuntherLab/READv2"
READV2_DIR = SCRIPT_DIR / "READv2"
READV2_SCRIPT = READV2_DIR / "READ2.py"
PLINK_PREFIX = DATA_PREFIX


def ensure_readv2(readv2_dir: Path) -> None:
    if (readv2_dir / "READ2.py").exists():
        return
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
    return f"cd {quote_path(work_dir)} && pixi run python {quote_path(readv2_script)} -i {quote_path(prefix)}"


def main() -> None:
    data_prefix = Path(DATA_PREFIX)
    ensure_data_present(data_prefix, PLINK_EXTS)
    ensure_readv2(READV2_DIR)

    results_dir = SCRIPT_DIR / "results"
    export_path = results_dir / "readv2_comparison_benchmark.csv"

    output_dir = tempfile.mkdtemp()
    fastpmr_cmd = build_fastpmr_command(PLINK_PREFIX, output_dir)
    readv2_cmd = build_readv2_command(READV2_SCRIPT, PLINK_PREFIX, output_dir)
    configs = [("fastpmr", fastpmr_cmd), ("READv2", readv2_cmd)]
    run_benchmark(configs, export_path)


if __name__ == "__main__":
    main()
