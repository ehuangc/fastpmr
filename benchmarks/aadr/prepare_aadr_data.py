import subprocess
from pathlib import Path

from aadr_utils import SCRIPT_DIR

DATA_DIR = SCRIPT_DIR / "data"
REMOTE_FILES = [
    (
        "v62.0_1240k_public.anno",
        "https://dataverse.harvard.edu/api/access/datafile/10537413?version=9.1",
    ),
    (
        "v62.0_1240k_public.ind",
        "https://dataverse.harvard.edu/api/access/datafile/10537414?version=9.1",
    ),
    (
        "v62.0_1240k_public.snp",
        "https://dataverse.harvard.edu/api/access/datafile/10537415?version=9.1",
    ),
    (
        "v62.0_1240k_public.geno",
        "https://dataverse.harvard.edu/api/access/datafile/10537126?version=9.1",
    ),
]


def download_file(url: str, destination: Path) -> None:
    print(f"Downloading {url}...")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["curl", "-L", "-#", "-o", str(destination), url], check=True)


def download_aadr_dataset() -> None:
    for file_name, url in REMOTE_FILES:
        download_file(url, DATA_DIR / file_name)
    print(f"Downloaded AADR data -> {DATA_DIR}")


def main() -> None:
    download_aadr_dataset()


if __name__ == "__main__":
    main()
