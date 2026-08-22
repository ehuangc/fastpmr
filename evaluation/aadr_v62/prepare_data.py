from pathlib import Path

from evaluation_utils.core import download_verified_file

DATA_DIR = Path(__file__).resolve().parent / "data"
# Checksums published by Dataverse
REMOTE_FILES = [
    (
        "v62.0.p1.1240k.aadr.anno",
        "https://dataverse.harvard.edu/api/access/datafile/13994492?version=11.0",
        "6468eb195e6cad5b98dac4a1286a42ae",
    ),
    (
        "v62.0.p1.1240k.aadr.ind",
        "https://dataverse.harvard.edu/api/access/datafile/13987487?version=11.0",
        "3f23dd87521dc1ccf6335c262571ecef",
    ),
    (
        "v62.0.p1.1240k.aadr.snp",
        "https://dataverse.harvard.edu/api/access/datafile/13987485?version=11.0",
        "50f66178fc81b8aa087cc4b135317e59",
    ),
    (
        "v62.0.p1.1240k.aadr.geno",
        "https://dataverse.harvard.edu/api/access/datafile/13994086?version=11.0",
        "24419bbaf673274e21fc9e2dc3dd145f",
    ),
]


def download_aadr_dataset() -> None:
    for file_name, url, md5 in REMOTE_FILES:
        download_verified_file(url, DATA_DIR / file_name, md5)


def main() -> None:
    download_aadr_dataset()


if __name__ == "__main__":
    main()
