from pathlib import Path

from evaluation_utils.core import download_file

DATA_DIR = Path(__file__).resolve().parent / "data"
REMOTE_FILES = [
    (
        "v62.0.p1_1240k_public.anno",
        "https://dataverse.harvard.edu/api/access/datafile/13994492?version=11.0",
    ),
    (
        "v62.0.p1_1240k_public.ind",
        "https://dataverse.harvard.edu/api/access/datafile/13987487?version=11.0",
    ),
    (
        "v62.0.p1_1240k_public.snp",
        "https://dataverse.harvard.edu/api/access/datafile/13987485?version=11.0",
    ),
    (
        "v62.0.p1_1240k_public.geno",
        "https://dataverse.harvard.edu/api/access/datafile/13994086?version=11.0",
    ),
]


def download_aadr_dataset() -> None:
    for file_name, url in REMOTE_FILES:
        download_file(url, DATA_DIR / file_name)


def main() -> None:
    download_aadr_dataset()


if __name__ == "__main__":
    main()
