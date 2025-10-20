use crate::Args;
use crate::counts::Counts;
use crate::error::{CustomError, Result};
use crate::output::{plot_mismatch_rates, write_mismatch_rates};
use crate::reader::SiteReader;
use crate::reader::packedancestrymap::PackedAncestryMapReader;
use crate::reader::transposed_packedancestrymap::TransposedPackedAncestryMapReader;
use rayon::ThreadPoolBuilder;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub enum InputSpec {
    PackedAncestryMap {
        ind: PathBuf,
        geno: PathBuf,
        snp: PathBuf,
        output_dir: PathBuf,
        // Parsed 0-based indices of variants to keep
        variant_indices: Option<HashSet<usize>>,
        threads: Option<usize>,
    },
}

impl InputSpec {
    pub fn from_prefix_packedancestrymap(
        prefix: &str,
        output_dir: &str,
        variant_indices: Option<HashSet<usize>>,
        threads: Option<usize>,
    ) -> Self {
        Self::PackedAncestryMap {
            ind: PathBuf::from(prefix.to_string() + ".ind"),
            geno: PathBuf::from(prefix.to_string() + ".geno"),
            snp: PathBuf::from(prefix.to_string() + ".snp"),
            output_dir: PathBuf::from(output_dir.to_string()),
            variant_indices,
            threads,
        }
    }

    pub fn print_paths(&self) {
        match self {
            InputSpec::PackedAncestryMap { ind, geno, snp, .. } => {
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
            InputSpec::PackedAncestryMap {
                ind,
                geno,
                snp,
                variant_indices,
                ..
            } => {
                // Check .geno header to determine if it's transposed or not
                let f = File::open(&geno).map_err(|e| crate::error::CustomError::ReadWithPath {
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
                    let reader =
                        PackedAncestryMapReader::open(ind, geno, snp, variant_indices.clone())?;
                    Ok(Box::new(reader))
                } else if header_prefix.starts_with(b"TGENO") {
                    let reader = TransposedPackedAncestryMapReader::open(
                        ind,
                        geno,
                        snp,
                        variant_indices.clone(),
                    )?;
                    Ok(Box::new(reader))
                } else {
                    Err(crate::error::CustomError::PackedAncestryMapHeaderPrefix)
                }
            }
        }
    }

    pub fn output_dir(&self) -> &Path {
        match self {
            InputSpec::PackedAncestryMap { output_dir, .. } => output_dir.as_path(),
        }
    }

    pub fn threads(&self) -> Option<usize> {
        match self {
            InputSpec::PackedAncestryMap { threads, .. } => *threads,
        }
    }
}

pub fn parse_indices(spec: &str) -> Result<HashSet<usize>> {
    let mut indices = HashSet::new();

    for raw in spec.split(',').map(|s| s.trim()).filter(|s| !s.is_empty()) {
        if let Some((a, b)) = raw.split_once('-') {
            let start: usize = a.parse().map_err(|e| CustomError::VariantIndexInt {
                source: e,
                arg: a.to_string(),
            })?;
            let end: usize = b.parse().map_err(|e| CustomError::VariantIndexInt {
                source: e,
                arg: b.to_string(),
            })?;

            if start == 0 || end == 0 {
                return Err(CustomError::VariantIndexLow);
            }

            let (lo, hi) = if start <= end {
                (start, end)
            } else {
                (end, start)
            };
            for i in lo..=hi {
                indices.insert(i - 1); // convert to 0-based index
            }
        } else {
            let idx: usize = raw.parse().map_err(|e| CustomError::VariantIndexInt {
                source: e,
                arg: raw.to_string(),
            })?;
            if idx == 0 {
                return Err(CustomError::VariantIndexLow);
            }
            indices.insert(idx - 1);
        }
    }
    Ok(indices)
}

pub fn build_input_spec(args: &Args) -> Result<InputSpec> {
    let variant_indices = match &args.variant_indices_spec {
        Some(spec) => Some(parse_indices(spec)?),
        None => None,
    };
    Ok(InputSpec::from_prefix_packedancestrymap(
        &args.prefix,
        &args.output_directory,
        variant_indices,
        args.threads,
    ))
}

pub fn run(
    reader: &mut dyn SiteReader,
    output_dir: impl AsRef<Path>,
    threads: Option<usize>,
) -> Result<()> {
    const PARALLEL_THRESHOLD: usize = 500;
    let samples: Vec<String> = reader.samples().to_vec();

    let mut counts = Counts::new(samples);
    if (threads.is_none() && counts.n_samples() < PARALLEL_THRESHOLD) || threads == Some(1) {
        counts = counts.consume_reader(reader)?;
    } else if let Some(n) = threads {
        println!("Using {} threads", n);
        let pool = ThreadPoolBuilder::new().num_threads(n).build()?;
        counts = pool.install(|| counts.consume_reader_parallel(reader))?;
    } else {
        counts = counts.consume_reader_parallel(reader)?;
    }

    let rates_path = output_dir.as_ref().join("mismatch_rates.csv");
    println!(
        "Writing pairwise mismatch rates to {}...",
        rates_path.display()
    );
    write_mismatch_rates(&counts, &rates_path)?;

    let plot_path = output_dir.as_ref().join("mismatch_rates.png");
    println!(
        "Writing pairwise mismatch rate plot to {}...",
        plot_path.display()
    );
    plot_mismatch_rates(&counts, &plot_path)?;
    Ok(())
}
