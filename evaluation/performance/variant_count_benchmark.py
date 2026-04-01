import tempfile
from pathlib import Path

from evaluation_utils import (
    FASTPMR_BIN,
    PERFORMANCE_DATA_PREFIX,
    PERFORMANCE_DIR,
    PERFORMANCE_RUNS,
    ensure_data_present,
    quote_path,
    run_benchmark,
)

VARIANT_SPECS = (
    "1-100000",
    "1-200000",
    "1-300000",
    "1-400000",
    "1-500000",
    "1-600000",
    "1-700000",
    "1-800000",
    "1-900000",
    "1-1000000",
    "1-1100000",
    "1-1150639",  # Run all autosomal variants
)


def build_command(prefix: Path, spec: str, output_dir: Path) -> str:
    parts = [
        FASTPMR_BIN,
        f"--prefix {quote_path(prefix)}",
        f"--output-directory {quote_path(output_dir)}",
        f"--variant-indices {spec}",
        "--chromosomes 1-22",
        "--min-covered-snps 0",
        "--npz",
    ]
    return " ".join(parts)


def main() -> None:
    variant_specs = list(VARIANT_SPECS)
    data_prefix = Path(PERFORMANCE_DATA_PREFIX)
    ensure_data_present(data_prefix)

    results_dir = PERFORMANCE_DIR / "results"
    export_path = results_dir / "variant_count_benchmark.csv"

    configs = []
    for spec in variant_specs:
        output_dir = tempfile.mkdtemp()
        command = build_command(data_prefix, spec, output_dir)
        configs.append((f"variants={spec}", command))
    run_benchmark(configs, export_path, runs=PERFORMANCE_RUNS)


if __name__ == "__main__":
    main()
