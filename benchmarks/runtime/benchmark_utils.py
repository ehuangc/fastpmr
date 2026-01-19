import shlex
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_PREFIX = SCRIPT_DIR / "data" / "IEdata"
FASTPMR_BIN = REPO_ROOT / "target" / "release" / "fastpmr"
DATA_EXTS = (".ind", ".snp", ".geno")
RUNS = 3


def quote_path(path: Path) -> str:
    return shlex.quote(str(path))


def ensure_data_present(prefix: Path) -> None:
    def data_file(ext: str) -> Path:
        return prefix.parent / f"{prefix.name}{ext}"

    missing = [data_file(ext) for ext in DATA_EXTS if not data_file(ext).is_file()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        raise SystemExit(
            f"Missing data files: {missing_str}. "
            "Run `pixi run prepare-runtime-data` to download and unpack the Indo-European dataset."
        )
