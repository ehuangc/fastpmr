import csv
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

EVALUATION_DIR = Path(__file__).resolve().parent.parent
# Local cargo build; switch back to "fastpmr" once the bioconda package is up to date
FASTPMR_BIN = str(EVALUATION_DIR.parent / "target" / "release" / "fastpmr")
AADR_DIR = EVALUATION_DIR / "aadr"
AADR_DATA_PREFIX = AADR_DIR / "data" / "v62.0_1240k_public"
AADR_EXTS = (".anno", ".ind", ".snp", ".geno")
AADR_RUNS = 1

COMPARISON_DIR = EVALUATION_DIR / "comparison"
COMPARISON_DATA_PREFIX = COMPARISON_DIR / "data" / "southern_cone"

PERFORMANCE_DIR = EVALUATION_DIR / "performance"
PERFORMANCE_DATA_PREFIX = PERFORMANCE_DIR / "data" / "IEdata"
PERFORMANCE_SAMPLE_SET_DIR = PERFORMANCE_DIR / "data" / "indo_european_sample_sets"
PERFORMANCE_SAMPLE_SET_SIZES = (
    128,
    256,
    384,
    512,
    640,
    768,
    896,
    1024,
    1152,
    1280,
)
PERFORMANCE_RUNS = 3

EIGENSTRAT_EXTS = (".ind", ".snp", ".geno")
PLINK_EXTS = (".bed", ".bim", ".fam")
CSV_FIELDS = [
    "label",
    "mean_s",
    "stddev_s",
    "min_s",
    "max_s",
    "mean_bytes",
    "stddev_bytes",
    "min_bytes",
    "max_bytes",
]


def download_file(url: str, destination: Path) -> None:
    print(f"Downloading {url}...")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["curl", "-L", "-#", "-o", str(destination), url], check=True)
    print(f"Downloaded {url} -> {destination}\n")


def extract_files(archive_path: Path, destination: Path, prefix: Path, exts: tuple[str, ...]) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    target_exts = set(exts)
    extracted = []

    if archive_path.suffix == ".zip":
        with zipfile.ZipFile(archive_path) as zf:
            for member in zf.namelist():
                path = Path(member)
                if path.suffix in target_exts:
                    out_path = prefix.with_suffix(path.suffix)
                    out_path.write_bytes(zf.read(member))
                    extracted.append(out_path.name)
    elif "".join(archive_path.suffixes[-2:]) == ".tar.gz":
        with tarfile.open(archive_path, "r:gz") as tar:
            for member in tar.getmembers():
                path = Path(member.name)
                if path.suffix in target_exts:
                    out_path = prefix.with_suffix(path.suffix)
                    file_obj = tar.extractfile(member)
                    out_path.write_bytes(file_obj.read())
                    extracted.append(out_path.name)
    else:
        raise ValueError(f"Unsupported archive format: {archive_path.suffix}")

    print(f"Extracted {', '.join(sorted(extracted))} -> {destination}\n")


def quote_path(path: Path) -> str:
    return shlex.quote(str(path))


def ensure_data_present(prefix: Path, exts: tuple[str, ...] = EIGENSTRAT_EXTS) -> None:
    def data_file(ext: str) -> Path:
        return prefix.parent / f"{prefix.name}{ext}"

    missing = [data_file(ext) for ext in exts if not data_file(ext).is_file()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        if prefix == AADR_DATA_PREFIX:
            command = "pixi run prepare-aadr"
            dataset = "AADR dataset"
        elif prefix == COMPARISON_DATA_PREFIX:
            command = "pixi run prepare-comparison"
            dataset = "Southern Cone dataset"
        else:
            command = "pixi run prepare-performance"
            dataset = "Indo-European dataset"
        raise SystemExit(f"Missing data files: {missing_str}. Run `{command}` to download and unpack the {dataset}.")


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
