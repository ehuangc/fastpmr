from pathlib import Path

from evaluation_utils import (
    AADR_DATA_PREFIX,
    AADR_DIR,
    AADR_EXTS,
    AADR_RUNS,
    FASTPMR_BIN,
    ensure_data_present,
    quote_path,
    run_benchmark,
)


def build_command(prefix: Path, output_dir: Path) -> str:
    parts = [
        FASTPMR_BIN,
        f"--prefix {quote_path(prefix)}",
        f"--output-directory {quote_path(output_dir)}",
        "--chromosomes 1-22",
        "--min-covered-snps 0",
        "--ci",
        "--npz",
    ]
    return " ".join(parts)


def main() -> None:
    ensure_data_present(AADR_DATA_PREFIX, AADR_EXTS)

    results_dir = AADR_DIR / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    fastpmr_output_dir = results_dir / "fastpmr"
    fastpmr_output_dir.mkdir(parents=True, exist_ok=True)
    benchmark_output_path = results_dir / "fastpmr" / "aadr_benchmark.csv"

    command = build_command(AADR_DATA_PREFIX, fastpmr_output_dir)
    run_benchmark([("aadr", command)], benchmark_output_path, AADR_RUNS)


if __name__ == "__main__":
    main()
