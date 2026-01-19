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
    "1-1200000",
)


def build_command(fastpmr_bin: Path, prefix: Path, spec: str, output_dir: Path) -> str:
    parts = [
        quote_path(fastpmr_bin),
        f"--prefix {quote_path(prefix)}",
        f"--output-directory {quote_path(output_dir)}",
        f"--variant-indices {spec}",
        "-n",
    ]
    return " ".join(parts)


def main() -> None:
    variant_specs = list(VARIANT_SPECS)
    data_prefix = Path(DATA_PREFIX)
    fastpmr_bin = Path(FASTPMR_BIN)
    ensure_data_present(data_prefix)

    results_dir = SCRIPT_DIR / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    export_path = results_dir / "variant_count_benchmark.csv"
    hyperfine_args = [
        "hyperfine",
        "--runs",
        str(RUNS),
        "--export-csv",
        str(export_path),
    ]
    for spec in variant_specs:
        with tempfile.TemporaryDirectory() as output_dir:
            command = build_command(fastpmr_bin, data_prefix, spec, output_dir)
            hyperfine_args.extend(["-n", f"variants={spec}", command])

    subprocess.run(hyperfine_args, check=True)
    print(f"\nHyperfine results written to {export_path}")


if __name__ == "__main__":
    main()
