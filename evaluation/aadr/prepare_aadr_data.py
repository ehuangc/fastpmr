from evaluation_utils import AADR_DIR, download_file

DATA_DIR = AADR_DIR / "data"
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


def download_aadr_dataset() -> None:
    for file_name, url in REMOTE_FILES:
        download_file(url, DATA_DIR / file_name)


def main() -> None:
    download_aadr_dataset()


if __name__ == "__main__":
    main()
