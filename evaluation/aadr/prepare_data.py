from evaluation_utils.constants import AADR_DIR
from evaluation_utils.core import download_verified_file

DATA_DIR = AADR_DIR / "data"
# Checksums published by Dataverse
REMOTE_FILES = [
    (
        "v66.p1_1240K.aadr.PUB.anno",
        "https://dataverse.harvard.edu/api/access/datafile/13994515?version=14.0",
        "a2db1ac16f0f3558ed66fb251e1d5c7d",
    ),
    (
        "v66.p1_1240K.aadr.PUB.ind",
        "https://dataverse.harvard.edu/api/access/datafile/13994513?version=14.0",
        "19a434ac954bcd10dbb8dba1d1188a09",
    ),
    (
        "v66.p1_1240K.aadr.PUB.snp",
        "https://dataverse.harvard.edu/api/access/datafile/13994514?version=14.0",
        "50f66178fc81b8aa087cc4b135317e59",
    ),
    (
        "v66.p1_1240K.aadr.PUB.geno",
        "https://dataverse.harvard.edu/api/access/datafile/13994829?version=14.0",
        "5ea1d2675a271c81e55b8f8b08b3ff3b",
    ),
]


def download_aadr_dataset() -> None:
    for file_name, url, md5 in REMOTE_FILES:
        download_verified_file(url, DATA_DIR / file_name, md5)


def main() -> None:
    download_aadr_dataset()


if __name__ == "__main__":
    main()
