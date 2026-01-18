import subprocess
import tempfile
from pathlib import Path

from benchmark_utils import (
    DATA_PREFIX,
    FASTPMR_BIN,
    RUNS,
    SCRIPT_DIR,
    quote_path,
)

SAMPLE_DIR = SCRIPT_DIR / "data" / "indo_european_sample_sets"


def build_command(
    fastpmr_bin: Path,
    prefix: Path,
    csv_path: Path,
    output_dir: Path,
) -> str:
    parts = [
        quote_path(fastpmr_bin),
        f"--prefix {quote_path(prefix)}",
        f"--output-directory {quote_path(output_dir)}",
        f"--sample-pairs-csv {quote_path(csv_path)}",
        "-n",
    ]
    return " ".join(parts)


def main() -> None:
    data_prefix = Path(DATA_PREFIX)
    fastpmr_bin = Path(FASTPMR_BIN)
    sample_dir = Path(SAMPLE_DIR)

    def sample_size(path: Path) -> int:
        return int(path.stem.rsplit("_", 1)[1])

    csv_files = sorted(sample_dir.glob("indo_european_samples_*.csv"), key=sample_size)

    results_dir = SCRIPT_DIR / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    export_path = results_dir / "pair_count_benchmark.csv"
    hyperfine_args = [
        "hyperfine",
        "--runs",
        str(RUNS),
        "--export-csv",
        str(export_path),
    ]

    def count_lines(path: Path) -> int:
        with path.open("r", encoding="utf-8") as handle:
            return sum(1 for _ in handle)

    for csv_path in csv_files:
        sample_count = count_lines(csv_path)
        pair_count = sample_count * (sample_count - 1) // 2
        with tempfile.TemporaryDirectory() as output_dir:
            command = build_command(fastpmr_bin, data_prefix, csv_path, output_dir)
            hyperfine_args.extend(["-n", f"samples={sample_count}_pairs={pair_count}", command])
    subprocess.run(hyperfine_args, check=True)
    print(f"\nHyperfine results written to {export_path}")


if __name__ == "__main__":
    main()
