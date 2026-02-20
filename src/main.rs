mod cli;
mod counts;
mod error;
mod model;
mod output;
mod reader;

use crate::error::Result;
use chrono::Local;
use clap::Parser;
use miette::IntoDiagnostic;

/// Compute pairwise mismatch rates between genetic sequences.
#[derive(Parser, Debug)]
#[command(version, about)]
pub struct Args {
    /// Input file prefix.
    #[arg(short, long)]
    prefix: String,

    /// Output directory. Defaults to "./fastpmr_output_<timestamp>".
    #[arg(
        short,
        long,
        value_hint = clap::ValueHint::DirPath,
        default_value_t = format!("fastpmr_output_{}", Local::now().format("%Y%m%d_%H%M%S"))
    )]
    output_directory: String,

    /// Flag to write outputs in compressed .npz format. If true, count matrices, covered
    /// SNP counts, and a samples.json will be written to mismatch_counts.npz. We recommend
    /// setting this flag to true for large sample sizes. Defaults to false.
    #[arg(short, long, default_value_t = false)]
    npz: bool,

    /// CSV file containing either explicit sample pairs (two columns) or a single
    /// column of samples for which all pairwise mismatch rates will be calculated.
    /// No header row is expected.
    #[arg(short, long, value_hint = clap::ValueHint::FilePath)]
    sample_pairs_csv: Option<String>,

    /// Minimum SNPs covered per sample. Samples with fewer covered SNPs will be excluded
    /// from pairwise mismatch rate calculations. Set to 0 to disable filtering. Defaults
    /// to 30000.
    #[arg(long = "min-covered-snps", default_value_t = 30000)]
    min_covered_snps: u64,

    /// 1-based, inclusive range(s) of variant indices to keep.
    /// Examples: "1-5000,10000-20000", "1,2,3000-4000".
    #[arg(short, long = "variant-indices")]
    variant_indices_spec: Option<String>,

    /// Number of threads to use. When run with fewer than 500 samples, defaults to 1.
    /// When run with 500 or more samples, defaults to the number of logical cores.
    #[arg(short, long)]
    threads: Option<usize>,
}

fn try_main() -> Result<()> {
    let args = Args::parse();
    std::fs::create_dir_all(&args.output_directory)
        .map_err(|e| error::CustomError::OutputDir { source: e })?;

    let input_spec = cli::build_input_spec(&args)?;
    input_spec.print_paths();

    cli::run(&input_spec)?;
    Ok(())
}

fn main() -> miette::Result<()> {
    try_main().into_diagnostic()
}
