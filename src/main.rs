mod cli;
mod counts;
mod error;
mod model;
mod output;
mod reader;

use crate::error::Result;
use clap::Parser;
use miette::IntoDiagnostic;

/// Compute pairwise mismatch rates between genomes
#[derive(Parser, Debug)]
#[command(version, about)]
pub struct Args {
    /// Input file prefix
    #[arg(short, long)]
    prefix: String,
}

fn try_main() -> Result<()> {
    let args = Args::parse();
    let input_spec = cli::build_input_spec(&args)?;
    input_spec.print_paths();

    let mut reader = input_spec.open_reader()?;
    cli::run(reader.as_mut())?;
    Ok(())
}

fn main() -> miette::Result<()> {
    try_main().into_diagnostic()
}
