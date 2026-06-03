from pathlib import Path

from evaluation_utils import (
    AADR_EXTS,
    AADR_RUNS,
    FASTPMR_BIN,
    quote_path,
    run_benchmark,
)

AADR_V62_DIR = Path(__file__).resolve().parent
AADR_V62_DATA_PREFIX = AADR_V62_DIR / "data" / "v62.0_1240k_public"


def ensure_v62_data_present() -> None:
    data_files = [AADR_V62_DATA_PREFIX.parent / f"{AADR_V62_DATA_PREFIX.name}{ext}" for ext in AADR_EXTS]
    missing = [path for path in data_files if not path.is_file()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        raise SystemExit(
            f"Missing data files: {missing_str}. Run `pixi run prepare-aadr-v62` to download the v62 AADR dataset."
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
    ensure_v62_data_present()

    results_dir = AADR_V62_DIR / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    fastpmr_output_dir = results_dir / "fastpmr"
    fastpmr_output_dir.mkdir(parents=True, exist_ok=True)
    benchmark_output_path = results_dir / "fastpmr" / "aadr_benchmark.csv"

    command = build_command(AADR_V62_DATA_PREFIX, fastpmr_output_dir)
    run_benchmark([("aadr_v62", command)], benchmark_output_path, AADR_RUNS)


if __name__ == "__main__":
    main()
