import subprocess
from pathlib import Path

from benchmark_utils import AADR_DATA_PREFIX, AADR_DIR, AADR_EXTS, AADR_RUNS, ensure_data_present, quote_path


def build_command(prefix: Path, output_dir: Path) -> str:
    parts = [
        "fastpmr",
        f"--prefix {quote_path(prefix)}",
        f"--output-directory {quote_path(output_dir)}",
        "-n",
    ]
    return " ".join(parts)


def main() -> None:
    ensure_data_present(AADR_DATA_PREFIX, AADR_EXTS)

    results_dir = AADR_DIR / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    fastpmr_output_dir = results_dir / "fastpmr"
    fastpmr_output_dir.mkdir(parents=True, exist_ok=True)
    hyperfine_output_dir = results_dir / "runtime"
    hyperfine_output_dir.mkdir(parents=True, exist_ok=True)
    hyperfine_export_path = hyperfine_output_dir / "aadr_benchmark.csv"

    command = build_command(AADR_DATA_PREFIX, fastpmr_output_dir)
    hyperfine_args = [
        "hyperfine",
        "--runs",
        str(AADR_RUNS),
        "--show-output",
        "--export-csv",
        str(hyperfine_export_path),
        "-n",
        "fastpmr-aadr",
        command,
    ]

    subprocess.run(hyperfine_args, check=True)
    print(f"\nHyperfine results written to {hyperfine_export_path}")
    print(f"fastpmr outputs written under {fastpmr_output_dir}")


if __name__ == "__main__":
    main()
