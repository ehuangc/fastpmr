mod cli;
mod counts;
mod error;
mod model;
mod output;
mod reader;

use crate::error::{CustomError, Result};
use miette::IntoDiagnostic;
use std::env;

fn parse_args() -> Result<Vec<String>> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() {
        return Err(CustomError::Args);
    }
    Ok(args)
}

fn try_main() -> Result<()> {
    let args = parse_args()?;
    let input_spec = cli::build_input_spec(&args)?;
    input_spec.print_paths();

    let mut reader = input_spec.open_reader()?;
    cli::run(&mut reader)?;
    Ok(())
}

fn main() -> miette::Result<()> {
    try_main().into_diagnostic()
}
