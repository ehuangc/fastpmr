import tempfile
from pathlib import Path

from benchmark_utils import (
    RUNTIME_DATA_PREFIX,
    RUNTIME_DIR,
    ensure_data_present,
    quote_path,
    run_benchmark,
)

THREAD_COUNTS = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512)


def build_command(prefix: Path, threads: int, output_dir: Path) -> str:
    parts = [
        "fastpmr",
        f"--prefix {quote_path(prefix)}",
        f"--output-directory {quote_path(output_dir)}",
        f"--threads {threads}",
        "--min-covered-snps 0",
        "-n",
    ]
    return " ".join(parts)


def main() -> None:
    data_prefix = Path(RUNTIME_DATA_PREFIX)
    ensure_data_present(data_prefix)

    results_dir = RUNTIME_DIR / "results"
    export_path = results_dir / "thread_count_benchmark.csv"

    configs = []
    for count in THREAD_COUNTS:
        output_dir = tempfile.mkdtemp()
        command = build_command(data_prefix, count, output_dir)
        configs.append((f"threads={count}", command))
    run_benchmark(configs, export_path)


if __name__ == "__main__":
    main()
