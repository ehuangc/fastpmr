import csv
import hashlib
import json
import multiprocessing
import resource
import shlex
import statistics
import subprocess
import sys
import tarfile
import time
import zipfile
from multiprocessing.connection import Connection
from pathlib import Path

import numpy as np
import reverse_geocoder as rg

from evaluation_utils.constants import (
    AADR_METADATA_PATH,
    AADR_NPZ_PATH,
    COUNTRY_TO_REGION,
    CSV_FIELDS,
    EIGENSTRAT_EXTS,
    GROUP_ID_FIELD,
    US_OCEANIA_ADMIN1,
)


def download_file(url: str, destination: Path) -> None:
    print(f"Downloading {url}...")
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Space out requests to lighten server load and avoid the rate-limiting 403s the ENA server
    # returns when many files are fetched back-to-back
    time.sleep(2)
    # Download to a temporary file and rename only on success, so an interrupted download
    # never leaves a partial file at the final path
    partial = destination.with_name(destination.name + ".part")
    # The ENA server drops large transfers mid-download. wget retries and resumes within this single
    # invocation by default, so progress accumulates across drops until the file completes. Omitting
    # -c redownloads from scratch so that every run downloads fresh, and -O writes that fresh download
    # over our .part instead of a new .part.1. The server also returns transient 403/5xx under load,
    # which wget treats as fatal unless the codes are listed in --retry-on-http-error.
    subprocess.run(
        [
            "wget",
            "-nv",
            "--tries=100",
            "--waitretry=10",
            "--retry-connrefused",
            "--retry-on-http-error=403,429,500,502,503,504",
            "--timeout=60",
            "-O",
            str(partial),
            url,
        ],
        check=True,
    )
    partial.replace(destination)
    print(f"Downloaded {url} -> {destination}\n")


def get_file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def download_verified_file(url: str, destination: Path, md5: str, attempts: int = 5) -> None:
    """Download url to destination and verify it against the checksum the source publishes."""
    if destination.is_file() and get_file_md5(destination) == md5:
        print(f"Skipping {url} ({destination.name} already matches its checksum)\n")
        return

    for attempt in range(1, attempts + 1):
        download_file(url, destination)
        downloaded_md5 = get_file_md5(destination)
        if downloaded_md5 == md5:
            return
        print(f"Checksum mismatch for {destination.name}: expected {md5}, got {downloaded_md5} ({attempt}/{attempts})")
        destination.unlink()
    raise SystemExit(f"{url} failed its checksum after {attempts} attempts")


def extract_files(archive_path: Path, destination: Path, prefix: Path, exts: tuple[str, ...]) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    target_exts = set(exts)
    extracted = []

    if archive_path.suffix == ".zip":
        with zipfile.ZipFile(archive_path) as zf:
            for member in zf.namelist():
                path = Path(member)
                if path.suffix in target_exts:
                    out_path = prefix.with_name(prefix.name + path.suffix)
                    out_path.write_bytes(zf.read(member))
                    extracted.append(out_path.name)
    elif "".join(archive_path.suffixes[-2:]) == ".tar.gz":
        with tarfile.open(archive_path, "r:gz") as tar:
            for member in tar.getmembers():
                path = Path(member.name)
                if path.suffix in target_exts:
                    out_path = prefix.with_name(prefix.name + path.suffix)
                    file_obj = tar.extractfile(member)
                    out_path.write_bytes(file_obj.read())
                    extracted.append(out_path.name)
    else:
        raise ValueError(f"Unsupported archive format: {archive_path.suffix}")

    print(f"Extracted {', '.join(sorted(extracted))} -> {destination}\n")


def quote_path(path: Path) -> str:
    return shlex.quote(str(path))


def ensure_aadr_npz_present(npz_path: Path = AADR_NPZ_PATH, metadata_path: Path = AADR_METADATA_PATH) -> None:
    missing = [path for path in (npz_path, metadata_path) if not path.is_file()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        raise SystemExit(
            f"Missing input files: {missing_str}. Run `pixi run fastpmr-aadr` after `pixi run prepare-aadr`."
        )


def load_aadr_samples(npz_path: Path = AADR_NPZ_PATH) -> list[str]:
    with zipfile.ZipFile(npz_path) as archive:
        with archive.open("samples.json") as handle:
            return json.loads(handle.read().decode("utf-8"))


def load_aadr_npz_arrays(
    npz_path: Path = AADR_NPZ_PATH,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with np.load(npz_path, allow_pickle=False) as npz:
        site_overlaps = npz["n_site_overlaps"]
        mismatch_rates = npz["mismatch_rates"]
        mismatch_rates_95_ci_lower = npz["mismatch_rates_95_ci_lower"]
        mismatch_rates_95_ci_upper = npz["mismatch_rates_95_ci_upper"]
        covered_snps = npz["covered_snps"]
    samples = load_aadr_samples(npz_path)
    return samples, site_overlaps, mismatch_rates, mismatch_rates_95_ci_lower, mismatch_rates_95_ci_upper, covered_snps


def load_aadr_metadata(metadata_path: Path = AADR_METADATA_PATH) -> dict[str, dict[str, str]]:
    with metadata_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        metadata = {}
        for row in reader:
            if not row:
                continue
            metadata[row[0]] = dict(zip(header, row, strict=True))
    return metadata


def is_archaic_or_reference_sample(sample: str, sample_metadata: dict[str, str]) -> bool:
    group = sample_metadata[GROUP_ID_FIELD].lower()
    if "neanderthal" in group or "denisova" in group:
        return True
    if ".ref" in sample.lower():
        return True
    return False


def ensure_data_present(prefix: Path, exts: tuple[str, ...] = EIGENSTRAT_EXTS, *, command: str, dataset: str) -> None:
    def data_file(ext: str) -> Path:
        return prefix.parent / f"{prefix.name}{ext}"

    missing = [data_file(ext) for ext in exts if not data_file(ext).is_file()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        raise SystemExit(f"Missing data files: {missing_str}. Run `{command}` to prepare the {dataset}.")


def classify_coords(coords: list[tuple[float, float]]) -> list[str]:
    """Classify (lat, lon) coordinates into regions using reverse_geocoder."""
    results = rg.search(coords)
    regions: list[str] = []
    for (lat, lon), result in zip(coords, results, strict=True):
        cc = result["cc"]
        if cc == "US" and (result["admin1"] in US_OCEANIA_ADMIN1 or (lon > 144 and lat < 25)):
            regions.append("oceania")
        else:
            regions.append(COUNTRY_TO_REGION[cc])
    return regions


def measure_worker(command: str, conn: Connection) -> None:
    start = time.perf_counter()
    subprocess.run(command, shell=True, check=True)
    runtime = time.perf_counter() - start

    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_rss = usage.ru_maxrss
    if sys.platform != "darwin":
        peak_rss *= 1024  # Linux reports KB, macOS reports bytes
    conn.send({"runtime": runtime, "peak_bytes": peak_rss})
    conn.close()


def measure_command(command: str) -> dict[str, float]:
    ctx = multiprocessing.get_context("spawn")
    parent_conn, child_conn = ctx.Pipe(duplex=False)
    proc = ctx.Process(target=measure_worker, args=(command, child_conn))
    proc.start()
    child_conn.close()
    result = parent_conn.recv()
    parent_conn.close()
    proc.join()
    if proc.exitcode != 0:
        raise RuntimeError(f"Measurement process exited with code {proc.exitcode}")
    return result


def run_benchmark(
    configs: list[tuple[str, str]],
    output_path: Path,
    runs: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    total = len(configs)
    for i, (label, command) in enumerate(configs, 1):
        runtimes = []
        peak_bytes_list = []
        for trial in range(1, runs + 1):
            print(f"\n[{i}/{total}] {label} (trial {trial}/{runs})")
            result = measure_command(command)
            runtimes.append(result["runtime"])
            peak_bytes_list.append(result["peak_bytes"])
        rows.append(
            {
                "label": label,
                "mean_s": statistics.mean(runtimes),
                "stddev_s": statistics.stdev(runtimes) if len(runtimes) > 1 else 0.0,
                "min_s": min(runtimes),
                "max_s": max(runtimes),
                "mean_bytes": statistics.mean(peak_bytes_list),
                "stddev_bytes": statistics.stdev(peak_bytes_list) if len(peak_bytes_list) > 1 else 0.0,
                "min_bytes": min(peak_bytes_list),
                "max_bytes": max(peak_bytes_list),
            }
        )
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nResults written to {output_path}")


# READv2 transposes its input into "<prefix>_test.*" beside the input files
# and only cleans those up when given a relative path
def remove_readv2_intermediates(prefix: Path) -> None:
    for path in sorted(prefix.parent.glob(f"{prefix.name}_test.*")):
        path.unlink()
        print(f"Removed READv2 intermediate file {path}")
