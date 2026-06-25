from evaluation_utils import AADR_DIR, download_file

DATA_DIR = AADR_DIR / "data"
REMOTE_FILES = [
    (
        "v66.1240K.aadr.PUB.anno",
        "https://dataverse.harvard.edu/api/access/datafile/13663706?version=10.0",
    ),
    (
        "v66.1240K.aadr.PUB.ind",
        "https://dataverse.harvard.edu/api/access/datafile/13663698?version=10.0",
    ),
    (
        "v66.1240K.aadr.PUB.snp",
        "https://dataverse.harvard.edu/api/access/datafile/13664260?version=10.0",
    ),
    (
        "v66.1240K.aadr.PUB.geno",
        "https://dataverse.harvard.edu/api/access/datafile/13664080?version=10.0",
    ),
]


def download_aadr_dataset() -> None:
    for file_name, url in REMOTE_FILES:
        download_file(url, DATA_DIR / file_name)


def main() -> None:
    download_aadr_dataset()


if __name__ == "__main__":
    main()
