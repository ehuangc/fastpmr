from evaluation_utils import (
    COMPARISON_DATA_PREFIX,
    PLINK_EXTS,
    download_file,
    extract_files,
)

DATA_DIR = COMPARISON_DATA_PREFIX.parent
ARGENTINA_URL = "https://dataverse.harvard.edu/api/access/datafile/12077020?version=2.0"
ARCHIVE_PATH = DATA_DIR / "argentina_dataset.zip"


def main() -> None:
    download_file(ARGENTINA_URL, ARCHIVE_PATH)
    extract_files(ARCHIVE_PATH, DATA_DIR, COMPARISON_DATA_PREFIX, PLINK_EXTS)


if __name__ == "__main__":
    main()
