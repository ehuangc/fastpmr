import random
import shutil
import subprocess
import tempfile
from pathlib import Path

from plinkio import plinkfile

from evaluation_utils.constants import (
    COMPARISON_DATA_PREFIX,
    COMPARISON_SAMPLE_SET_DIR,
    COMPARISON_SAMPLE_SET_SIZES,
    PLINK_EXTS,
)
from evaluation_utils.core import (
    download_verified_file,
    extract_files,
)

DATA_DIR = COMPARISON_DATA_PREFIX.parent
MARAVALL_LOPEZ_URL = "https://dataverse.harvard.edu/api/access/datafile/12077020?version=2.0"
# Checksum published by Dataverse
MARAVALL_LOPEZ_MD5 = "a03bfca18ad51fbffe2c5bac4e27a11f"
ARCHIVE_PATH = DATA_DIR / "maravall_lopez_dataset.zip"
SEX_CHROMOSOMES = {23, 24}
SAMPLE_SHUFFLE_SEED = 42


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


def read_sample_ids(fam_path: Path) -> list[tuple[str, str]]:
    samples: list[tuple[str, str]] = []
    with fam_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if parts:
                samples.append((parts[0], parts[1]))
    return samples


def subset_prefix(size: int) -> Path:
    return COMPARISON_SAMPLE_SET_DIR / f"{COMPARISON_DATA_PREFIX.name}_samples_{size}"


# Subset samples at the data level for the sample count benchmark because READv2, unlike
# fastpmr, doesn't support filtering samples with a command-line flag
def generate_sample_subsets(prefix: Path, sample_set_dir: Path) -> None:
    samples = read_sample_ids(prefix.with_suffix(".fam"))
    rng = random.Random(SAMPLE_SHUFFLE_SEED)
    shuffled = samples.copy()
    rng.shuffle(shuffled)

    sample_set_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        keep_path = Path(tmpdir) / "keep.txt"
        for size in COMPARISON_SAMPLE_SET_SIZES:
            subset = sorted(shuffled[:size])
            keep_path.write_text("".join(f"{family_id}\t{individual_id}\n" for family_id, individual_id in subset))
            subprocess.run(
                [
                    "plink",
                    "--bfile",
                    str(prefix),
                    "--keep",
                    str(keep_path),
                    # Preserve the original A1/A2 assignment so genotype codes match the full dataset
                    "--keep-allele-order",
                    "--allow-no-sex",
                    "--make-bed",
                    "--silent",
                    "--out",
                    str(subset_prefix(size)),
                ],
                check=True,
            )

    print(f"Wrote sample subsets -> {sample_set_dir}")


def main() -> None:
    download_verified_file(MARAVALL_LOPEZ_URL, ARCHIVE_PATH, MARAVALL_LOPEZ_MD5)
    extract_files(ARCHIVE_PATH, DATA_DIR, COMPARISON_DATA_PREFIX, PLINK_EXTS)
    filter_sex_chromosomes(COMPARISON_DATA_PREFIX, SEX_CHROMOSOMES)
    generate_sample_subsets(COMPARISON_DATA_PREFIX, COMPARISON_SAMPLE_SET_DIR)


if __name__ == "__main__":
    main()
