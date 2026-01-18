from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_PREFIX = SCRIPT_DIR / "data" / "v62.0_1240k_public"
FASTPMR_BIN = REPO_ROOT / "target" / "release" / "fastpmr"
DATA_EXTS = (".anno", ".ind", ".snp", ".geno")


def ensure_data_present(prefix: Path) -> None:
    missing = [prefix.with_suffix(ext) for ext in DATA_EXTS if not prefix.with_suffix(ext).is_file()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        raise SystemExit(
            f"Missing data files: {missing_str}. Run `pixi run prepare-aadr-data` to download the AADR dataset."
        )
