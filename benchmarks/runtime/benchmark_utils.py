import csv
import multiprocessing
import resource
import shlex
import statistics
import subprocess
import sys
import time
from multiprocessing.connection import Connection
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_PREFIX = SCRIPT_DIR / "data" / "IEdata"
FASTPMR_BIN = REPO_ROOT / "target" / "release" / "fastpmr"
EIGENSTRAT_EXTS = (".ind", ".snp", ".geno")
PLINK_EXTS = (".bed", ".bim", ".fam")
RUNS = 3
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


def quote_path(path: Path) -> str:
    return shlex.quote(str(path))


def ensure_data_present(prefix: Path, exts: tuple[str, ...] = EIGENSTRAT_EXTS) -> None:
    def data_file(ext: str) -> Path:
        return prefix.parent / f"{prefix.name}{ext}"

    missing = [data_file(ext) for ext in exts if not data_file(ext).is_file()]
    if missing:
        missing_str = ", ".join(str(path) for path in missing)
        raise SystemExit(
            f"Missing data files: {missing_str}. "
            "Run `pixi run prepare-runtime-data` to download and unpack the Indo-European dataset."
        )


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
    runs: int = RUNS,
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
