import shlex
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
