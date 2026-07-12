import subprocess
from pathlib import Path

import pandas as pd

from evaluation_utils.constants import (
    CHICHEN_ITZA_GENO_DIR,
    CHICHEN_ITZA_RESULTS_DIR,
    EIGENSTRAT_EXTS,
    FASTPMR_BIN,
)
from evaluation_utils.core import (
    ensure_data_present,
    quote_path,
)

# eager's pileupcaller genotypes single- and double-stranded libraries into separate
# EIGENSTRAT sets. The two sets hold disjoint individuals (no YCH sample has both library
# types), so mergeit concatenates them into one set with no renaming.
CHICHEN_ITZA_GENO_SINGLE_PREFIX = CHICHEN_ITZA_GENO_DIR / "pileupcaller.single"
CHICHEN_ITZA_GENO_DOUBLE_PREFIX = CHICHEN_ITZA_GENO_DIR / "pileupcaller.double"
CHICHEN_ITZA_GENO_MERGED_PREFIX = CHICHEN_ITZA_GENO_DIR / "pileupcaller.merged"


def merge_eigenstrat(single_prefix: Path, double_prefix: Path, merged_prefix: Path) -> Path:
    par_path = merged_prefix.parent / f"{merged_prefix.name}.par"
    par_path.write_text(
        f"geno1: {single_prefix}.geno\n"
        f"snp1: {single_prefix}.snp\n"
        f"ind1: {single_prefix}.ind\n"
        f"geno2: {double_prefix}.geno\n"
        f"snp2: {double_prefix}.snp\n"
        f"ind2: {double_prefix}.ind\n"
        f"genooutfilename: {merged_prefix}.geno\n"
        f"snpoutfilename: {merged_prefix}.snp\n"
        f"indoutfilename: {merged_prefix}.ind\n"
        "outputformat: EIGENSTRAT\n"
        "docheck: YES\n"
        "hashcheck: NO\n"
        "strandcheck: NO\n"
    )
    subprocess.run(["mergeit", "-p", str(par_path)], check=True)
    return merged_prefix


def build_command(prefix: Path, output_dir: Path) -> str:
    parts = [
        FASTPMR_BIN,
        f"--prefix {quote_path(prefix)}",
        f"--output-directory {quote_path(output_dir)}",
        "--chromosomes 1-22",
        "--min-covered-snps 0",
        "--ci",
    ]
    return " ".join(parts)


def sort_results(results_path: Path) -> Path:
    # Sort fastpmr pair results by increasing PMR
    sorted_path = results_path.with_name(f"{results_path.stem}_sorted{results_path.suffix}")
    df = pd.read_csv(results_path).sort_values("mismatch_rate")
    df.to_csv(sorted_path, index=False)
    return sorted_path


def main() -> None:
    ensure_data_present(CHICHEN_ITZA_GENO_SINGLE_PREFIX, EIGENSTRAT_EXTS)
    ensure_data_present(CHICHEN_ITZA_GENO_DOUBLE_PREFIX, EIGENSTRAT_EXTS)
    merge_eigenstrat(CHICHEN_ITZA_GENO_SINGLE_PREFIX, CHICHEN_ITZA_GENO_DOUBLE_PREFIX, CHICHEN_ITZA_GENO_MERGED_PREFIX)

    fastpmr_output_dir = CHICHEN_ITZA_RESULTS_DIR / "fastpmr"
    fastpmr_output_dir.mkdir(parents=True, exist_ok=True)

    command = build_command(CHICHEN_ITZA_GENO_MERGED_PREFIX, fastpmr_output_dir)
    subprocess.run(command, shell=True, check=True)
    sort_results(fastpmr_output_dir / "fastpmr_pair_results.csv")


if __name__ == "__main__":
    main()
