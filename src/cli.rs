use crate::Args;
use crate::counts::Counts;
use crate::error::{CustomError, Result};
use crate::output::{plot_mismatch_rates, write_counts_npz, write_mismatch_rates};
use crate::reader::SiteReader;
use crate::reader::packedancestrymap::PackedAncestryMapReader;
use crate::reader::transposed_packedancestrymap::TransposedPackedAncestryMapReader;
use rayon::ThreadPoolBuilder;
use std::collections::{HashMap, HashSet};
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
        npz: bool,
        sample_pairs: Option<Vec<(String, String)>>,
        // Parsed 0-based indices of variants to keep
        variant_indices: Option<HashSet<usize>>,
        threads: Option<usize>,
    },
}

impl InputSpec {
    pub fn from_prefix_packedancestrymap(
        prefix: &str,
        output_dir: &str,
        npz: bool,
        sample_pairs: Option<Vec<(String, String)>>,
        variant_indices: Option<HashSet<usize>>,
        threads: Option<usize>,
    ) -> Self {
        Self::PackedAncestryMap {
            ind: PathBuf::from(prefix.to_string() + ".ind"),
            geno: PathBuf::from(prefix.to_string() + ".geno"),
            snp: PathBuf::from(prefix.to_string() + ".snp"),
            output_dir: PathBuf::from(output_dir.to_string()),
            npz,
            sample_pairs,
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

    pub fn npz(&self) -> bool {
        match self {
            InputSpec::PackedAncestryMap { npz, .. } => *npz,
        }
    }

    pub fn sample_pairs(&self) -> Option<&[(String, String)]> {
        match self {
            InputSpec::PackedAncestryMap { sample_pairs, .. } => sample_pairs.as_deref(),
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
    let sample_pairs = match &args.sample_pairs_csv {
        Some(path) => Some(load_sample_pairs_csv(path)?),
        None => None,
    };
    Ok(InputSpec::from_prefix_packedancestrymap(
        &args.prefix,
        &args.output_directory,
        args.npz,
        sample_pairs,
        variant_indices,
        args.threads,
    ))
}

fn load_sample_pairs_csv(path: &str) -> Result<Vec<(String, String)>> {
    let csv_path = PathBuf::from(path);
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(&csv_path)
        .map_err(|source| CustomError::CsvRead {
            source,
            path: csv_path.clone(),
        })?;

    let mut pairs = Vec::new();
    for result in reader.records() {
        let record = result.map_err(|source| CustomError::CsvRead {
            source,
            path: csv_path.clone(),
        })?;

        if record.iter().all(|field| field.trim().is_empty()) {
            continue;
        }

        if record.len() != 2 {
            return Err(CustomError::SamplePairsColumns);
        }

        let id1 = record[0].trim();
        let id2 = record[1].trim();

        if pairs.is_empty() && id1.eq_ignore_ascii_case("id1") && id2.eq_ignore_ascii_case("id2") {
            continue;
        }
        if id1.is_empty() || id2.is_empty() {
            continue;
        }
        pairs.push((id1.to_string(), id2.to_string()));
    }

    if pairs.is_empty() {
        return Err(CustomError::SamplePairsEmpty);
    }
    Ok(pairs)
}

/// Convert sample pairs from a vector of pairs of IDs to a set of pairs of indices
fn resolve_sample_pairs(
    samples: &[String],
    pairs: &[(String, String)],
) -> Result<HashSet<(usize, usize)>> {
    let mut lookup = HashMap::with_capacity(samples.len());
    for (idx, sample) in samples.iter().enumerate() {
        lookup.insert(sample.as_str(), idx);
    }

    let mut to_keep = HashSet::new();
    for (left_raw, right_raw) in pairs {
        let left = left_raw.trim();
        let right = right_raw.trim();
        if left == right {
            return Err(CustomError::SamplePairDuplicate {
                sample: left.to_string(),
            });
        }

        let left_idx =
            lookup
                .get(left)
                .copied()
                .ok_or_else(|| CustomError::SamplePairUnknownSample {
                    sample: left.to_string(),
                })?;
        let right_idx =
            lookup
                .get(right)
                .copied()
                .ok_or_else(|| CustomError::SamplePairUnknownSample {
                    sample: right.to_string(),
                })?;

        let (lo, hi) = if left_idx < right_idx {
            (left_idx, right_idx)
        } else {
            (right_idx, left_idx)
        };
        if lo == hi {
            return Err(CustomError::SamplePairDuplicate {
                sample: samples[lo].clone(),
            });
        }
        to_keep.insert((lo, hi));
    }

    if to_keep.is_empty() {
        return Err(CustomError::SamplePairsEmpty);
    }
    Ok(to_keep)
}

pub fn run(
    reader: &mut dyn SiteReader,
    output_dir: impl AsRef<Path>,
    npz: bool,
    threads: Option<usize>,
    sample_pairs: Option<&[(String, String)]>,
) -> Result<()> {
    const PARALLEL_THRESHOLD: usize = 500;
    let samples: Vec<String> = reader.samples().to_vec();

    let pairs_to_keep = sample_pairs
        .map(|pairs| resolve_sample_pairs(&samples, pairs))
        .transpose()?;
    let mut counts = Counts::new(samples, pairs_to_keep);
    if (threads.is_none() && counts.n_samples() < PARALLEL_THRESHOLD) || threads == Some(1) {
        counts = counts.consume_reader(reader)?;
    } else if let Some(n) = threads {
        let pool = ThreadPoolBuilder::new().num_threads(n).build()?;
        counts = pool.install(|| counts.consume_reader_parallel(reader))?;
    } else {
        counts = counts.consume_reader_parallel(reader)?;
    }

    if npz {
        let npz_path = output_dir.as_ref().join("mismatch_counts.npz");
        println!(
            "Writing pairwise mismatch counts to {}...",
            npz_path.display()
        );
        write_counts_npz(&counts, &npz_path)?;
    } else {
        let rates_path = output_dir.as_ref().join("mismatch_rates.csv");
        println!(
            "Writing pairwise mismatch rates to {}...",
            rates_path.display()
        );
        write_mismatch_rates(&counts, &rates_path)?;
    }

    let plot_path = output_dir.as_ref().join("mismatch_rates.png");
    println!(
        "Writing pairwise mismatch rate plot to {}...",
        plot_path.display()
    );
    plot_mismatch_rates(&counts, &plot_path)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_sample_pairs_succeeds() {
        let samples = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let pairs = vec![
            ("A".to_string(), "C".to_string()),
            ("C".to_string(), "B".to_string()),
        ];
        let keep = resolve_sample_pairs(&samples, &pairs).expect("pairs should resolve");
        assert_eq!(keep.len(), 2);
        assert!(keep.contains(&(0, 2)));
        assert!(keep.contains(&(1, 2)));
    }

    #[test]
    fn resolve_sample_pairs_rejects_unknown_sample() {
        let samples = vec!["A".to_string(), "B".to_string()];
        let pairs = vec![("A".to_string(), "Z".to_string())];
        let err = resolve_sample_pairs(&samples, &pairs).unwrap_err();
        match err {
            CustomError::SamplePairUnknownSample { sample } => assert_eq!(sample, "Z"),
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn resolve_sample_pairs_rejects_duplicate_member() {
        let samples = vec!["A".to_string(), "B".to_string()];
        let pairs = vec![("A".to_string(), "A".to_string())];
        let err = resolve_sample_pairs(&samples, &pairs).unwrap_err();
        match err {
            CustomError::SamplePairDuplicate { sample } => assert_eq!(sample, "A"),
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
