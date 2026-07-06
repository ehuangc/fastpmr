import json
import re
import subprocess

from evaluation_utils import (
    CHICHEN_ITZA_BAM_DIR,
    CHICHEN_ITZA_FILEREPORT,
    CHICHEN_ITZA_RENAMED_BAM_DIR,
    download_file,
)


def collect_urls() -> list[str]:
    records = json.loads(CHICHEN_ITZA_FILEREPORT.read_text())
    return [
        "https://" + entry.strip()
        for record in records
        for entry in record.get("submitted_ftp", "").split(";")
        if entry.strip().endswith((".bam", ".bai"))
    ]


def download_bams() -> None:
    CHICHEN_ITZA_BAM_DIR.mkdir(parents=True, exist_ok=True)
    for url in collect_urls():
        destination = CHICHEN_ITZA_BAM_DIR / url.rsplit("/", 1)[-1]
        download_file(url, destination)


def reheader_line(line: str) -> str:
    # Rename contigs on @SQ lines to match hs37d5: chr<x> -> <x>, chrM/chrMT -> MT
    if not line.startswith("@SQ"):
        return line
    line = re.sub(r"\tSN:chrMT?\t", "\tSN:MT\t", line)
    return re.sub(r"\tSN:chr", "\tSN:", line)


def rename_bams() -> None:
    CHICHEN_ITZA_RENAMED_BAM_DIR.mkdir(parents=True, exist_ok=True)
    for bam in sorted(CHICHEN_ITZA_BAM_DIR.glob("*.bam")):
        out = CHICHEN_ITZA_RENAMED_BAM_DIR / bam.name
        print(f"Reheadering {bam.name}")
        header = subprocess.run(["samtools", "view", "-H", bam], check=True, capture_output=True, text=True).stdout
        new_header = "".join(reheader_line(line) for line in header.splitlines(keepends=True))
        with out.open("wb") as handle:
            subprocess.run(["samtools", "reheader", "-", bam], input=new_header.encode(), stdout=handle, check=True)
        subprocess.run(["samtools", "index", out], check=True)


def main() -> None:
    download_bams()
    rename_bams()


if __name__ == "__main__":
    main()
