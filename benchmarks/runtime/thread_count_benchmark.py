import subprocess
import tempfile
from pathlib import Path

from benchmark_utils import (
    DATA_PREFIX,
    FASTPMR_BIN,
    RUNS,
    SCRIPT_DIR,
    ensure_data_present,
    quote_path,
)

THREAD_COUNTS = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512)


def build_command(fastpmr_bin: Path, prefix: Path, threads: int, output_dir: Path) -> str:
    parts = [
        quote_path(fastpmr_bin),
        f"--prefix {quote_path(prefix)}",
        f"--output-directory {quote_path(output_dir)}",
        f"--threads {threads}",
        "-n",
    ]
    return " ".join(parts)


def main() -> None:
    data_prefix = Path(DATA_PREFIX)
    fastpmr_bin = Path(FASTPMR_BIN)
    thread_counts = THREAD_COUNTS
    ensure_data_present(data_prefix)

    results_dir = SCRIPT_DIR / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    export_path = results_dir / "thread_count_benchmark.csv"
    hyperfine_args = [
        "hyperfine",
        "--runs",
        str(RUNS),
        "--show-output",
        "--export-csv",
        str(export_path),
    ]

    for count in thread_counts:
        with tempfile.TemporaryDirectory() as output_dir:
            command = build_command(fastpmr_bin, data_prefix, count, output_dir)
            hyperfine_args.extend(["-n", f"threads={count}", command])
    subprocess.run(hyperfine_args, check=True)
    print(f"\nHyperfine results written to {export_path}")


if __name__ == "__main__":
    main()
