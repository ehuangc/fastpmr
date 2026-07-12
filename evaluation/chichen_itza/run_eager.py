import csv
import json
import shutil
import subprocess
from pathlib import Path

from evaluation_utils.constants import (
    CHICHEN_ITZA_CUSTOM_CONFIG,
    CHICHEN_ITZA_FILEREPORT,
    CHICHEN_ITZA_RENAMED_BAM_DIR,
    CHICHEN_ITZA_RESULTS_DIR,
    CHICHEN_ITZA_TSV,
    CHICHEN_ITZA_WORK_DIR,
)

# ============================= CONFIG ========================================
# Strandedness is resolved per library from the file report (see get_strandedness)
# since the YCH libraries are a mix of single- and double-stranded
UDG_TREATMENT = "half"  # none | half | full
COLOUR_CHEM = "4"  # 2 = NovaSeq/NextSeq, 4 = HiSeq/MiSeq
SEQTYPE = "SE"  # Merged aDNA BAMs are usually SE

# Reference + cluster paths
# The paths below are specific to the GRACE cluster at MPI-EVA (Leipzig)
PROFILE = "eva_grace"  # Public nf-core institutional config
EAGER_VER = "2.4.0"
REFDIR = Path("/mnt/archgen/Reference_Genomes/Human/hs37d5")
FASTA = REFDIR / "hs37d5.fa"
FASTA_IDX = REFDIR / "hs37d5.fa.fai"
BWA_IDX = REFDIR  # Directory holding the bwa index
SEQ_DICT = REFDIR / "hs37d5.dict"
SNP_BED = REFDIR / "SNPCapBEDs" / "1240K.pos.list_hs37d5.0based.bed"
SNP_FILE = Path("/mnt/archgen/public_data/Datashare_Boston_Jena_June2018.backup/1240K.snp")
# =============================================================================

TSV_HEADER = [
    "Sample_Name",
    "Library_ID",
    "Lane",
    "Colour_Chemistry",
    "SeqType",
    "Organism",
    "Strandedness",
    "UDG_Treatment",
    "R1",
    "R2",
    "BAM",
]


def get_strandedness(record: dict) -> str:
    # The YCH libraries are a mix of single- and double-stranded, and the A/B library code
    # does not encode which. The submitted library_construction_protocol field is the only
    # reliable signal.
    protocol = record.get("library_construction_protocol", "")
    assert "single stranded" in protocol or "double stranded" in protocol
    if "single stranded" in protocol:
        return "single"
    else:
        return "double"


def build_eager_tsv() -> None:
    if not any(CHICHEN_ITZA_RENAMED_BAM_DIR.glob("*.bam")):
        raise SystemExit(f"No reheadered BAMs in {CHICHEN_ITZA_RENAMED_BAM_DIR} — run prepare-bams-chichen-itza first")

    rows = []
    for record in json.loads(CHICHEN_ITZA_FILEREPORT.read_text()):
        bams = [e.strip() for e in record.get("submitted_ftp", "").split(";") if e.strip().endswith(".bam")]
        if not bams:
            continue
        stem = bams[0].rsplit("/", 1)[-1].removesuffix(".bam")  # e.g. YCH003.A0102.1240K
        sample = stem.split(".")[0]  # YCH003 -> merge key
        bam = CHICHEN_ITZA_RENAMED_BAM_DIR / f"{stem}.bam"
        rows.append(
            [
                sample,
                stem,
                "1",
                COLOUR_CHEM,
                SEQTYPE,
                "Homo_sapiens",
                get_strandedness(record),
                UDG_TREATMENT,
                "NA",
                "NA",
                str(bam),
            ]
        )
    rows.sort()  # Group libraries of each sample together (Sample_Name is the first field)

    with CHICHEN_ITZA_TSV.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(TSV_HEADER)
        writer.writerows(rows)
    n_single = sum(1 for row in rows if row[6] == "single")
    print(
        f"Wrote {len(rows)} libraries to {CHICHEN_ITZA_TSV} "
        f"({n_single} single-stranded, {len(rows) - n_single} double-stranded)"
    )


def eager_command() -> list[str]:
    # Resolve nextflow from the pixi `eager` env so a login shell's PATH can't shadow it
    nextflow = shutil.which("nextflow") or "nextflow"
    return [
        nextflow,
        "run",
        "nf-core/eager",
        "-r",
        EAGER_VER,
        "-profile",
        PROFILE,
        "-resume",
        "--input",
        str(CHICHEN_ITZA_TSV),
        "--outdir",
        str(CHICHEN_ITZA_RESULTS_DIR),
        "-c",
        str(CHICHEN_ITZA_CUSTOM_CONFIG),
        "-w",
        str(CHICHEN_ITZA_WORK_DIR),
        "--run_bam_filtering",
        "--bam_mapping_quality_threshold",
        "25",
        "--bam_filter_minreadlength",
        "30",
        "--bam_unmapped_type",
        "discard",
        "--run_trim_bam",
        "--run_sexdeterrmine",
        "--run_nuclear_contamination",
        "--contamination_chrom_name",
        "X",
        "--run_mtnucratio",
        "--mtnucratio_header",
        "MT",
        "--run_genotyping",
        "--genotyping_source",
        "trimmed",
        "--genotyping_tool",
        "pileupcaller",
        "--pileupcaller_bedfile",
        str(SNP_BED),
        "--pileupcaller_snpfile",
        str(SNP_FILE),
        "--fasta",
        str(FASTA),
        "--fasta_index",
        str(FASTA_IDX),
        "--bwa_index",
        str(BWA_IDX),
        "--seq_dict",
        str(SEQ_DICT),
        "--snpcapture_bed",
        str(SNP_BED),
        "--sexdeterrmine_bedfile",
        str(SNP_BED),
        "--anno_file",
        str(SNP_BED),
    ]


def run_eager() -> None:
    CHICHEN_ITZA_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    CHICHEN_ITZA_WORK_DIR.mkdir(parents=True, exist_ok=True)
    # `module` is a shell function from the cluster's env-modules: load apptainer in a
    # login shell, then exec the nextflow argv unchanged so no manual quoting is needed
    subprocess.run(["bash", "-lc", 'module load apptainer && exec "$@"', "eager", *eager_command()], check=True)


def main() -> None:
    build_eager_tsv()
    run_eager()


if __name__ == "__main__":
    main()
