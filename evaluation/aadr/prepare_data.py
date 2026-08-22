from evaluation_utils.constants import AADR_DIR
from evaluation_utils.core import download_file

DATA_DIR = AADR_DIR / "data"
REMOTE_FILES = [
    (
        "v66.p1_1240K.aadr.PUB.anno",
        "https://dataverse.harvard.edu/api/access/datafile/13994515?version=14.0",
    ),
    (
        "v66.p1_1240K.aadr.PUB.ind",
        "https://dataverse.harvard.edu/api/access/datafile/13994513?version=14.0",
    ),
    (
        "v66.p1_1240K.aadr.PUB.snp",
        "https://dataverse.harvard.edu/api/access/datafile/13994514?version=14.0",
    ),
    (
        "v66.p1_1240K.aadr.PUB.geno",
        "https://dataverse.harvard.edu/api/access/datafile/13994829?version=14.0",
    ),
]


def download_aadr_dataset() -> None:
    for file_name, url in REMOTE_FILES:
        download_file(url, DATA_DIR / file_name)


def main() -> None:
    download_aadr_dataset()


if __name__ == "__main__":
    main()
