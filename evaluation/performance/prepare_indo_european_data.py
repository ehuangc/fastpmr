import random
from pathlib import Path

from evaluation_utils import (
    EIGENSTRAT_EXTS,
    PERFORMANCE_DATA_PREFIX,
    PERFORMANCE_SAMPLE_SET_DIR,
    PERFORMANCE_SAMPLE_SET_SIZES,
    download_file,
    extract_files,
)

DATA_DIR = PERFORMANCE_DATA_PREFIX.parent
SAMPLE_SHUFFLE_SEED = 42
INDO_EUROPEAN_URL = "https://dataverse.harvard.edu/api/access/datafile/10629469?version=1.2"
ARCHIVE_PATH = DATA_DIR / "indo_european_dataset.tar.gz"


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
    for size in PERFORMANCE_SAMPLE_SET_SIZES:
        subset = sorted(shuffled[:size])
        path = sample_dir / f"indo_european_samples_{size}.csv"
        path.write_text("\n".join(subset) + "\n")

    print(f"Wrote sample pair lists -> {sample_dir}")


def main() -> None:
    download_file(INDO_EUROPEAN_URL, ARCHIVE_PATH)
    extract_files(ARCHIVE_PATH, DATA_DIR, PERFORMANCE_DATA_PREFIX, EIGENSTRAT_EXTS)
    generate_sample_sets(PERFORMANCE_DATA_PREFIX, PERFORMANCE_SAMPLE_SET_DIR)


if __name__ == "__main__":
    main()
