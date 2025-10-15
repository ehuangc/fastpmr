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

    /// Output directory.
    #[arg(
        short,
        long,
        value_hint = clap::ValueHint::DirPath,
        default_value_t = format!("fastpmr_output_{}", Local::now().format("%Y%m%d_%H%M%S"))
    )]
    output_directory: String,

    /// 1-based, inclusive range(s) of variant indices to keep.
    /// Examples: "1-5000,10000-20000", "1,2,3000-4000".
    #[arg(short, long = "variant-indices")]
    variant_indices_spec: Option<String>,
}

fn try_main() -> Result<()> {
    let args = Args::parse();
    std::fs::create_dir_all(&args.output_directory)
        .map_err(|e| error::CustomError::OutputDir { source: e })?;

    let input_spec = cli::build_input_spec(&args)?;
    input_spec.print_paths();

    let mut reader = input_spec.open_reader()?;
    cli::run(reader.as_mut(), input_spec.output_dir())?;
    Ok(())
}

fn main() -> miette::Result<()> {
    try_main().into_diagnostic()
}
