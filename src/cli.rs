use crate::Args;
use crate::counts::Counts;
use crate::error::Result;
use crate::output::{plot_mismatch_rates, write_mismatch_rates};
use crate::reader::SiteReader;
use crate::reader::packedancestrymap::PackedAncestryMapReader;
use crate::reader::transposed_packedancestrymap::TransposedPackedAncestryMapReader;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub enum InputSpec {
    PackedAncestryMap {
        ind: PathBuf,
        geno: PathBuf,
        snp: PathBuf,
    },
}

impl InputSpec {
    pub fn from_prefix_packedancestrymap(prefix: &str) -> Self {
        Self::PackedAncestryMap {
            ind: PathBuf::from(prefix.to_string() + ".ind"),
            geno: PathBuf::from(prefix.to_string() + ".geno"),
            snp: PathBuf::from(prefix.to_string() + ".snp"),
        }
    }

    pub fn print_paths(&self) {
        match self {
            InputSpec::PackedAncestryMap { ind, geno, snp } => {
                println!("IND : {}", ind.display());
                println!("GENO: {}", geno.display());
                println!("SNP : {}", snp.display());
                println!();
            }
        }
    }

    // Open the appropriate reader for the given input spec
    pub fn open_reader(&self) -> Result<Box<dyn SiteReader>> {
        match self {
            InputSpec::PackedAncestryMap { ind, geno, snp } => {
                // Check .geno header to determine if it's transposed or not
                let f = File::open(geno).map_err(|e| crate::error::CustomError::ReadWithPath {
                    source: e,
                    path: geno.to_path_buf(),
                })?;
                let mut reader = BufReader::new(f);
                let buffer =
                    reader
                        .fill_buf()
                        .map_err(|e| crate::error::CustomError::ReadWithPath {
                            source: e,
                            path: geno.to_path_buf(),
                        })?;
                let header_prefix = &buffer[..buffer.len().min(5)];

                if header_prefix.starts_with(b"GENO") {
                    let reader = PackedAncestryMapReader::open(ind, geno, snp)?;
                    Ok(Box::new(reader))
                } else if header_prefix.starts_with(b"TGENO") {
                    let reader = TransposedPackedAncestryMapReader::open(ind, geno, snp)?;
                    Ok(Box::new(reader))
                } else {
                    Err(crate::error::CustomError::PackedAncestryMapHeaderPrefix)
                }
            }
        }
    }
}

pub fn build_input_spec(args: &Args) -> Result<InputSpec> {
    Ok(InputSpec::from_prefix_packedancestrymap(&args.prefix))
}

pub fn run(reader: &mut dyn SiteReader) -> Result<()> {
    const PARALLEL_THRESHOLD: usize = 500;
    let samples: Vec<String> = reader.samples().to_vec();

    let mut counts = Counts::new(samples);
    if counts.n_samples() < PARALLEL_THRESHOLD {
        counts = counts.consume_reader(reader)?;
    } else {
        counts = counts.consume_reader_parallel(reader)?;
    }

    let output_path = "mismatch_rates.csv";
    println!("Writing pairwise mismatch rates to {}...", output_path);
    write_mismatch_rates(&counts, output_path)?;

    let plot_path = "mismatch_rates.png";
    println!("Writing pairwise mismatch rate plot to {}...", plot_path);
    plot_mismatch_rates(&counts, plot_path)?;
    Ok(())
}
