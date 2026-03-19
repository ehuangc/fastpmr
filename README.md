🚀 <b>`fastpmr`</b> is a command-line tool, written in Rust, that calculates pairwise mismatch rates (PMRs) between ancient DNA sequences.

It supports EIGENSTRAT (unpacked, packed, and transposed packed formats) and PLINK inputs, as well as multithreading and sample/variant filtering.

## Installation
### Pixi (recommended)
```bash
pixi global install -c conda-forge -c bioconda fastpmr
```
We recommend installing `fastpmr` using [`Pixi`](https://pixi.prefix.dev/latest/), a fast and reproducible package manager. To install `Pixi`, see [this page](https://pixi.prefix.dev/latest/installation/).

### Conda
```bash
conda create -n fastpmr -c conda-forge -c bioconda fastpmr
conda activate fastpmr
```

## Usage

```bash
fastpmr -p PREFIX [-o OUTPUT_DIRECTORY] [-n] [-s SAMPLE_PAIRS_CSV] [-v VARIANT_INDICES] [-m MIN_COVERED_SNPS] [-t THREADS]
```

> [!NOTE]
> Minimal command:
> ```bash
> fastpmr -p /path/to/dataset_prefix
> ```

### Inputs
**Prefix** (`-p`, `--prefix`) (*required*): Input file prefix. `fastpmr` auto-detects supported formats from the files on disk. If both EIGENSTRAT and PLINK datasets exist for the same prefix, an explicit extension in `--prefix` can force selection (for example `--prefix data/study.bed`).

**Output Directory** (`-o`, `--output-directory`) (*optional*): Directory for outputs. Defaults to `./fastpmr_output_<timestamp>`.

**NPZ Output** (`-n`, `--npz`) (*flag*): Write compressed count matrices to `mismatch_counts.npz` instead of CSV count outputs. Recommended for large sample sets.

**Sample Pairs CSV** (`-s`, `--sample-pairs-csv`) (*optional*): CSV with no header that controls which sample pairs are computed. This parameter accepts 1-column CSVs and 2-column CSVs. 1-column CSVs specify a sample list, and PMRs for all pairs of these samples are computed. 2-column CSVs specify the exact sample pairs for which PMRs are computed.

**Minimum Covered SNPs** (`--min-covered-snps`) (*optional*): Exclude samples with too-few covered SNPs from PMR calculations. Defaults to `30000`. Set to `0` to disable.

**Variant Indices** (`-v`, `--variant-indices`) (*optional*): 1-based inclusive variant index ranges to keep, e.g., `1,2,3000-5000,10000-20000`.

**Threads** (`-t`, `--threads`) (*optional*): Number of threads to use. Default behavior is:
- `<500` samples: single-threaded
- `>=500` samples: uses logical core count

### Outputs

**Always written**:
- `mismatch_rates.png`: histogram of pairwise mismatch rates (%), including only pairs with `n_site_overlaps >= 30000`.

**When `--npz` is not set**:
- `covered_snps.csv`
- `mismatch_rates.csv`

**When `--npz` is set**:
- `mismatch_counts.npz`
  - Sample list: `samples.json`
  - Arrays: `covered_snps`, `mismatches`, `totals`, `site_overlaps`

## Reproducibility
We recommend using [`Pixi`](https://pixi.sh/) to reproduce benchmarks. `Pixi` benchmark tasks can be found in `benchmarks/pixi.toml`.
