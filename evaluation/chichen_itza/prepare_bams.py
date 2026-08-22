import hashlib
import json
import re
import subprocess
from pathlib import Path

from evaluation_utils.constants import (
    CHICHEN_ITZA_BAM_DIR,
    CHICHEN_ITZA_FILEREPORT,
    CHICHEN_ITZA_FILEREPORT_URL,
    CHICHEN_ITZA_RENAMED_BAM_DIR,
)
from evaluation_utils.core import download_file


def collect_downloads() -> list[tuple[str, str]]:
    downloads = []
    for record in json.loads(CHICHEN_ITZA_FILEREPORT.read_text()):
        # submitted_ftp and submitted_md5 are semicolon-separated lists in matching order
        entries = record.get("submitted_ftp", "").split(";")
        checksums = record.get("submitted_md5", "").split(";")
        for entry, md5 in zip(entries, checksums, strict=True):
            if entry.strip().endswith((".bam", ".bai")):
                downloads.append(("https://" + entry.strip(), md5.strip()))
    return downloads


def get_file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def download_bams(attempts: int = 5) -> None:
    """Download the submitted BAMs and BAIs, verifying each against its submitted checksum."""
    CHICHEN_ITZA_BAM_DIR.mkdir(parents=True, exist_ok=True)
    for url, md5 in collect_downloads():
        destination = CHICHEN_ITZA_BAM_DIR / url.rsplit("/", 1)[-1]
        if destination.is_file() and get_file_md5(destination) == md5:
            print(f"Skipping {url} ({destination.name} already matches its checksum)\n")
            continue

        for attempt in range(1, attempts + 1):
            download_file(url, destination)
            downloaded_md5 = get_file_md5(destination)
            if downloaded_md5 == md5:
                break
            print(
                f"Checksum mismatch for {destination.name}: expected {md5}, got {downloaded_md5} ({attempt}/{attempts})"
            )
            destination.unlink()
        else:
            raise SystemExit(f"{url} failed its checksum after {attempts} attempts")


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
    download_file(CHICHEN_ITZA_FILEREPORT_URL, CHICHEN_ITZA_FILEREPORT)
    download_bams()
    rename_bams()


if __name__ == "__main__":
    main()
