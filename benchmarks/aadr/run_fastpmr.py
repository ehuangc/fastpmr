import subprocess
from pathlib import Path

from aadr_utils import (
    DATA_PREFIX,
    FASTPMR_BIN,
    RUNS,
    SCRIPT_DIR,
    ensure_data_present,
    quote_path,
)


def build_command(fastpmr_bin: Path, prefix: Path, output_dir: Path) -> str:
    parts = [
        quote_path(fastpmr_bin),
        f"--prefix {quote_path(prefix)}",
        f"--output-directory {quote_path(output_dir)}",
        "-n",
    ]
    return " ".join(parts)


def main() -> None:
    data_prefix = Path(DATA_PREFIX)
    fastpmr_bin = Path(FASTPMR_BIN)
    ensure_data_present(data_prefix)

    results_dir = SCRIPT_DIR / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    fastpmr_output_dir = results_dir / "fastpmr"
    fastpmr_output_dir.mkdir(parents=True, exist_ok=True)
    hyperfine_output_dir = results_dir / "runtime"
    hyperfine_output_dir.mkdir(parents=True, exist_ok=True)
    hyperfine_export_path = hyperfine_output_dir / "aadr_benchmark.csv"

    command = build_command(fastpmr_bin, data_prefix, fastpmr_output_dir)
    hyperfine_args = [
        "hyperfine",
        "--runs",
        str(RUNS),
        "--show-output",
        "--export-csv",
        str(hyperfine_export_path),
        "-n",
        "fastpmr_aadr",
        command,
    ]

    subprocess.run(hyperfine_args, check=True)
    print(f"\nHyperfine results written to {hyperfine_export_path}")
    print(f"fastpmr outputs written under {fastpmr_output_dir}")


if __name__ == "__main__":
    main()
