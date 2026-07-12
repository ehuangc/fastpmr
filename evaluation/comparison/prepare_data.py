import shutil
import tempfile
from pathlib import Path

from plinkio import plinkfile

from evaluation_utils.constants import (
    COMPARISON_DATA_PREFIX,
    PLINK_EXTS,
)
from evaluation_utils.core import (
    download_file,
    extract_files,
)

DATA_DIR = COMPARISON_DATA_PREFIX.parent
MARAVALL_LOPEZ_URL = "https://dataverse.harvard.edu/api/access/datafile/12077020?version=2.0"
ARCHIVE_PATH = DATA_DIR / "maravall_lopez_et_al_dataset.zip"
SEX_CHROMOSOMES = {23, 24}


# Filter sex chromosomes at the data level for comparison benchmark because READv2 doesn't
# support filtering by chromosome
def filter_sex_chromosomes(prefix: Path, chrs_to_exclude: set[int]) -> None:
    plink_in = plinkfile.open(str(prefix))
    samples = plink_in.get_samples()
    loci = plink_in.get_loci()

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_prefix = str(Path(tmpdir) / prefix.name)
        plink_out = plinkfile.create(tmp_prefix, samples)
        for locus, row in zip(loci, plink_in, strict=True):
            if locus.chromosome not in chrs_to_exclude:
                plink_out.write_row(locus, row)
        plink_in.close()
        plink_out.close()

        for ext in (".bed", ".bim", ".fam"):
            shutil.move(tmp_prefix + ext, prefix.with_suffix(ext))

    print(f"Filtered chromosomes {sorted(chrs_to_exclude)} from {prefix}")


def main() -> None:
    download_file(MARAVALL_LOPEZ_URL, ARCHIVE_PATH)
    extract_files(ARCHIVE_PATH, DATA_DIR, COMPARISON_DATA_PREFIX, PLINK_EXTS)
    filter_sex_chromosomes(COMPARISON_DATA_PREFIX, SEX_CHROMOSOMES)


if __name__ == "__main__":
    main()
