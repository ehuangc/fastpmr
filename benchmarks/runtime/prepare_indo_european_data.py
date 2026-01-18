import random
import subprocess
import tarfile
from pathlib import Path

from benchmark_utils import DATA_PREFIX, DATA_EXTS

DATA_DIR = DATA_PREFIX.parent
SAMPLE_DIR = DATA_DIR / "indo_european_sample_sets"
SAMPLE_SIZES = (
    64,
    128,
    192,
    256,
    320,
    384,
    448,
    512,
    576,
    640,
    704,
    768,
    832,
    896,
    960,
    1024,
    1088,
    1152,
    1216,
    1280,
)
SAMPLE_SHUFFLE_SEED = 42
INDO_EUROPEAN_URL = (
    "https://dataverse.harvard.edu/api/access/datafile/10629469?version=1.2"
)
ARCHIVE_PATH = DATA_DIR / "indo_european_dataset.tar.gz"


def download_archive(url: str, destination: Path) -> None:
    print(f"Downloading {url}...")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["curl", "-L", "-#", "-o", str(destination), url], check=True)
    print(f"Downloaded {url} -> {destination}\n")


def extract_required_files(archive_path: Path, destination: Path, prefix: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    expected = {f"{prefix.name}{ext}": prefix.with_suffix(ext) for ext in DATA_EXTS}

    with tarfile.open(archive_path, "r:gz") as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            if name in expected:
                file_obj = tar.extractfile(member)
                expected[name].write_bytes(file_obj.read())

    print(f"Extracted {', '.join(sorted(expected))} -> {destination}\n")


def read_sample_ids(ind_path: Path) -> list[str]:
    samples: list[str] = []
    with ind_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split()
            if parts:
                samples.append(parts[0])
    return samples


def format_pair_count(size: int) -> str:
    pairs = size * (size - 1) // 2
    return f"{pairs:,}"


def generate_sample_sets(prefix: Path, sample_dir: Path) -> None:
    samples = read_sample_ids(prefix.with_suffix(".ind"))
    rng = random.Random(SAMPLE_SHUFFLE_SEED)
    shuffled = samples.copy()
    rng.shuffle(shuffled)

    sample_dir.mkdir(parents=True, exist_ok=True)
    for size in SAMPLE_SIZES:
        subset = sorted(shuffled[:size])
        path = sample_dir / f"indo_european_samples_{size}.csv"
        path.write_text("\n".join(subset) + "\n")

    print(f"Wrote sample pair lists -> {sample_dir}")


def main() -> None:
    download_archive(INDO_EUROPEAN_URL, ARCHIVE_PATH)
    extract_required_files(ARCHIVE_PATH, DATA_DIR, DATA_PREFIX)
    generate_sample_sets(DATA_PREFIX, SAMPLE_DIR)


if __name__ == "__main__":
    main()
