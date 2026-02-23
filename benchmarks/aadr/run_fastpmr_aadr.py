import shlex
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_PREFIX = SCRIPT_DIR / "data" / "v62.0_1240k_public"
FASTPMR_BIN = REPO_ROOT / "target" / "release" / "fastpmr"
DATA_EXTS = (".anno", ".ind", ".snp", ".geno")
RUNS = 1


def quote_path(path: Path) -> str:
    return shlex.quote(str(path))


def ensure_data_present(prefix: Path) -> None:
    def data_file(ext: str) -> Path:
        return prefix.parent / f"{prefix.name}{ext}"

    missing = [data_file(ext) for ext in DATA_EXTS if not data_file(ext).is_file()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        raise SystemExit(
            f"Missing data files: {missing_str}. Run `pixi run prepare-aadr-data` to download the AADR dataset."
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
        "fastpmr-aadr",
        command,
    ]

    subprocess.run(hyperfine_args, check=True)
    print(f"\nHyperfine results written to {hyperfine_export_path}")
    print(f"fastpmr outputs written under {fastpmr_output_dir}")


if __name__ == "__main__":
    main()
