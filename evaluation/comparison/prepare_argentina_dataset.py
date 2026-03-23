import subprocess
import zipfile
from pathlib import Path

from evaluation_utils import (
    COMPARISON_DATA_PREFIX,
    PLINK_EXTS,
)

DATA_DIR = COMPARISON_DATA_PREFIX.parent
ARGENTINA_URL = "https://dataverse.harvard.edu/api/access/datafile/12077020?version=2.0"
ARCHIVE_PATH = DATA_DIR / "argentina_dataset.zip"


def download_archive(url: str, destination: Path) -> None:
    print(f"Downloading {url}...")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["curl", "-L", "-#", "-o", str(destination), url], check=True)
    print(f"Downloaded {url} -> {destination}\n")


def extract_plink_files(archive_path: Path, destination: Path, prefix: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    target_exts = set(PLINK_EXTS)
    extracted = []

    with zipfile.ZipFile(archive_path) as zf:
        for member in zf.namelist():
            path = Path(member)
            if path.suffix in target_exts:
                out_path = prefix.with_suffix(path.suffix)
                out_path.write_bytes(zf.read(member))
                extracted.append(out_path.name)

    print(f"Extracted {', '.join(sorted(extracted))} -> {destination}\n")


def main() -> None:
    download_archive(ARGENTINA_URL, ARCHIVE_PATH)
    extract_plink_files(ARCHIVE_PATH, DATA_DIR, COMPARISON_DATA_PREFIX)


if __name__ == "__main__":
    main()
