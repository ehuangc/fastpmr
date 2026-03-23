import tempfile
from pathlib import Path

from evaluation_utils import (
    PERFORMANCE_DATA_PREFIX,
    PERFORMANCE_DIR,
    PERFORMANCE_RUNS,
    PERFORMANCE_SAMPLE_SET_DIR,
    PERFORMANCE_SAMPLE_SET_SIZES,
    ensure_data_present,
    quote_path,
    run_benchmark,
)

# Keep threads constant to avoid fastpmr's automatic switch from single- to multi-threaded mode at 500 samples
THREADS = 512


def ensure_sample_set_data_present() -> None:
    missing = [
        PERFORMANCE_SAMPLE_SET_DIR / f"indo_european_samples_{size}.csv"
        for size in PERFORMANCE_SAMPLE_SET_SIZES
        if not (PERFORMANCE_SAMPLE_SET_DIR / f"indo_european_samples_{size}.csv").is_file()
    ]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        raise SystemExit(
            f"Missing sample set files: {missing_str}. Run `pixi run prepare-performance` to generate them."
        )


def build_command(
    prefix: Path,
    csv_path: Path,
    output_dir: Path,
) -> str:
    parts = [
        "fastpmr",
        f"--prefix {quote_path(prefix)}",
        f"--output-directory {quote_path(output_dir)}",
        f"--sample-pairs-csv {quote_path(csv_path)}",
        f"--threads {THREADS}",
        "--min-covered-snps 0",
        "-n",
    ]
    return " ".join(parts)


def main() -> None:
    ensure_data_present(PERFORMANCE_DATA_PREFIX)
    ensure_sample_set_data_present()

    def sample_size(path: Path) -> int:
        return int(path.stem.rsplit("_", 1)[1])

    csv_files = sorted(PERFORMANCE_SAMPLE_SET_DIR.glob("indo_european_samples_*.csv"), key=sample_size)

    results_dir = PERFORMANCE_DIR / "results"
    export_path = results_dir / "pair_count_benchmark.csv"

    def count_lines(path: Path) -> int:
        with path.open("r", encoding="utf-8") as handle:
            return sum(1 for _ in handle)

    configs = []
    for csv_path in csv_files:
        sample_count = count_lines(csv_path)
        pair_count = sample_count * (sample_count - 1) // 2
        output_dir = tempfile.mkdtemp()
        command = build_command(PERFORMANCE_DATA_PREFIX, csv_path, output_dir)
        configs.append((f"samples={sample_count}_pairs={pair_count}", command))
    run_benchmark(configs, export_path, runs=PERFORMANCE_RUNS)


if __name__ == "__main__":
    main()
