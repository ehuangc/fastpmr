🚀 <b>`fastpmr`</b> is a command-line tool, written in Rust, that computes pairwise mismatch rates (PMRs) and (optionally) degrees of relatedness between ancient DNA sequences.

It supports EIGENSTRAT (unpacked, packed, and transposed packed formats) and PLINK inputs, as well as multithreading and sample/variant/chromosome filtering.

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
fastpmr -p PREFIX [-o OUTPUT_DIRECTORY] [-n] [-d] [-i] [-s SAMPLE_PAIRS_CSV] [-v VARIANT_INDICES] [-c CHROMOSOMES] [-m MIN_COVERED_SNPS] [-t THREADS]
```

> [!NOTE]
> Example command:
> ```bash
> fastpmr -p /path/to/dataset_prefix -c 1-22 
> ```
> This command runs `fastpmr` on the given dataset over all autosomes (chr 1-22).

### Inputs
**Prefix** (`-p`, `--prefix`) (*required*): Input file prefix. `fastpmr` auto-detects supported formats from the files on disk. If both EIGENSTRAT and PLINK datasets exist for the same prefix, an explicit extension in `--prefix` can force selection (for example `--prefix data/study.bed`).

**Output Directory** (`-o`, `--output-directory`) (*optional*): Directory for outputs. Defaults to `./fastpmr_output_<timestamp>`.

**NPZ Output** (`-n`, `--npz`) (*flag*): Write outputs as compressed `numpy` arrays to `fastpmr_results.npz` instead of as CSV files. Recommended for large sample sets.

**Degrees** (`-d`, `--degrees`) (*flag*): Call degrees of relatedness for each sample pair, using the procedure described in [READv2](https://doi.org/10.1186/s13059-024-03350-3) (Kuhn et al. 2018). Each pair is classified as Identical/Twin, First Degree, Second Degree, Third Degree, or Unrelated based on normalized mismatch rates (mismatch rate divided by the median across all pairs). Third degree inference additionally requires an "expected mismatch" count of at least 3000.

**Confidence Intervals** (`-i`, `--ci`) (*flag*): Compute 95% confidence intervals for pairwise mismatch rates using the Wald method.

**Sample Pairs CSV** (`-s`, `--sample-pairs-csv`) (*optional*): CSV with no header that controls which sample pairs are computed. This parameter accepts 1-column CSVs and 2-column CSVs. 1-column CSVs specify a sample list, and PMRs for all pairs of these samples are computed. 2-column CSVs specify the exact sample pairs for which PMRs are computed.

**Minimum Covered SNPs** (`--min-covered-snps`) (*optional*): Exclude samples with too-few covered SNPs from PMR calculations. Defaults to `30000`. Set to `0` to disable.

**Variant Indices** (`-v`, `--variant-indices`) (*optional*): 1-based inclusive variant index ranges to keep, e.g., `1,2,3000-5000,10000-20000`.

**Chromosomes** (`-c`, `--chromosomes`) (*optional*): Chromosomes to include, as comma-separated values or numeric ranges, e.g., `1-22`, `1-22,X,Y`, `X,Y,MT`, `23,24,90`. Chromosome values are matched against the chromosome column in the `.snp` or `.bim` file. When used with `--variant-indices`, only variants matching both filters are kept.

**Threads** (`-t`, `--threads`) (*optional*): Number of threads to use. Defaults to 1 (single-threaded) for datasets with fewer than 500 individuals and the number of logical cores for datasets with 500 or more individuals.

### Outputs

**Always written**:
- `mismatch_rates.png`: histogram of pairwise mismatch rates (%), including only pairs with `n_site_overlaps >= 30000`.

**When `--npz/-n` is not set**:
- `fastpmr_covered_snps.csv`
- `fastpmr_pair_results.csv`
  - When `--ci/-i` is set, includes additional `mismatch_rate_95_ci_lower` and `mismatch_rate_95_ci_upper` columns
  - When `--degrees/-d` is set, includes additional `normalized_mismatch_rate` and `degree` columns

**When `--npz/-n` is set**:
- `fastpmr_results.npz`
  - Sample list: `samples.json`
  - Arrays: `covered_snps`, `mismatches`, `totals`, `site_overlaps`
  - When `--ci/-i` is set, includes additional arrays `mismatch_rate_95_ci_lower` and `mismatch_rate_95_ci_upper`
  - When `--degrees/-d` is set, includes additional arrays `normalized_mismatch_rates`, `degrees`, and a `degree_labels.json` mapping

## Reproducibility
We recommend using [`Pixi`](https://pixi.sh/) to reproduce evaluations. `Pixi` evaluation tasks can be found in `evaluation/pyproject.toml`.
