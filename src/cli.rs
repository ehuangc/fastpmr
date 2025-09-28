use crate::counts::Counts;
use crate::error::Result;
use crate::output::{plot_mismatch_rates, write_mismatch_rates};
use crate::reader::SiteReader;
use crate::reader::packedancestrymap::PackedAncestryMapReader;
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
                let reader = PackedAncestryMapReader::open(ind, geno, snp)?;
                Ok(Box::new(reader))
            }
        }
    }
}

pub fn build_input_spec(args: &Vec<String>) -> Result<InputSpec> {
    if args.len() != 1 {
        return Err(crate::error::CustomError::Args);
    }
    let prefix = &args[0];
    Ok(InputSpec::from_prefix_packedancestrymap(prefix))
}

pub fn run(reader: &mut dyn SiteReader) -> Result<()> {
    let samples: Vec<String> = reader.samples().to_vec();
    let mut counts = Counts::new(samples);
    counts = counts.consume_reader(reader)?;

    let output_path = "mismatch_rates.csv";
    write_mismatch_rates(&counts, output_path)?;
    println!("Wrote pairwise mismatch rates to {}", output_path);

    let plot_path = "mismatch_rates.png";
    plot_mismatch_rates(&counts, plot_path)?;
    println!("Wrote mismatch rate plot to {}", plot_path);
    Ok(())
}
