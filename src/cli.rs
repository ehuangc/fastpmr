use crate::Args;
use crate::counts::Counts;
use crate::error::{CustomError, Result};
use crate::output::{plot_mismatch_rates, write_counts_npz, write_mismatch_rates};
use crate::reader::SiteReader;
use crate::reader::eigenstrat::EigenstratReader;
use crate::reader::packedancestrymap::PackedAncestryMapReader;
use crate::reader::plink::PlinkBedReader;
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
        samples_to_keep: Option<HashSet<String>>,
        // Parsed 0-based indices of variants to keep
        variant_indices: Option<HashSet<usize>>,
        threads: Option<usize>,
    },
    Plink {
        bed: PathBuf,
        bim: PathBuf,
        fam: PathBuf,
        output_dir: PathBuf,
        npz: bool,
        sample_pairs: Option<Vec<(String, String)>>,
        samples_to_keep: Option<HashSet<String>>,
        // Parsed 0-based indices of variants to keep
        variant_indices: Option<HashSet<usize>>,
        threads: Option<usize>,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FormatHint {
    PackedAncestryMap,
    Plink,
}

impl InputSpec {
    fn from_packedancestrymap(
        prefix: &Path,
        output_dir: &str,
        npz: bool,
        sample_pairs: Option<Vec<(String, String)>>,
        samples_to_keep: Option<HashSet<String>>,
        variant_indices: Option<HashSet<usize>>,
        threads: Option<usize>,
    ) -> Self {
        Self::PackedAncestryMap {
            ind: path_with_extension(prefix, "ind"),
            geno: path_with_extension(prefix, "geno"),
            snp: path_with_extension(prefix, "snp"),
            output_dir: PathBuf::from(output_dir.to_string()),
            npz,
            sample_pairs,
            samples_to_keep,
            variant_indices,
            threads,
        }
    }

    fn from_plink(
        prefix: &Path,
        output_dir: &str,
        npz: bool,
        sample_pairs: Option<Vec<(String, String)>>,
        samples_to_keep: Option<HashSet<String>>,
        variant_indices: Option<HashSet<usize>>,
        threads: Option<usize>,
    ) -> Self {
        Self::Plink {
            bed: path_with_extension(prefix, "bed"),
            bim: path_with_extension(prefix, "bim"),
            fam: path_with_extension(prefix, "fam"),
            output_dir: PathBuf::from(output_dir.to_string()),
            npz,
            sample_pairs,
            samples_to_keep,
            variant_indices,
            threads,
        }
    }

    pub fn from_prefix(
        prefix: &str,
        output_dir: &str,
        npz: bool,
        sample_pairs: Option<Vec<(String, String)>>,
        samples_to_keep: Option<HashSet<String>>,
        variant_indices: Option<HashSet<usize>>,
        threads: Option<usize>,
    ) -> Result<Self> {
        let (base_prefix, hint) = normalize_prefix(prefix);
        let packed_ind = path_with_extension(&base_prefix, "ind");
        let packed_geno = path_with_extension(&base_prefix, "geno");
        let packed_snp = path_with_extension(&base_prefix, "snp");
        let plink_bed = path_with_extension(&base_prefix, "bed");
        let plink_bim = path_with_extension(&base_prefix, "bim");
        let plink_fam = path_with_extension(&base_prefix, "fam");

        let packed_available =
            packed_ind.is_file() && packed_geno.is_file() && packed_snp.is_file();
        let plink_available = plink_bed.is_file() && plink_bim.is_file() && plink_fam.is_file();

        let selected_format = if packed_available
            && (!plink_available || hint == Some(FormatHint::PackedAncestryMap))
        {
            Some(FormatHint::PackedAncestryMap)
        } else if plink_available && (!packed_available || hint == Some(FormatHint::Plink)) {
            Some(FormatHint::Plink)
        } else if packed_available {
            Some(FormatHint::PackedAncestryMap)
        } else if plink_available {
            Some(FormatHint::Plink)
        } else {
            None
        };

        match selected_format {
            Some(FormatHint::PackedAncestryMap) => Ok(Self::from_packedancestrymap(
                &base_prefix,
                output_dir,
                npz,
                sample_pairs,
                samples_to_keep,
                variant_indices,
                threads,
            )),
            Some(FormatHint::Plink) => Ok(Self::from_plink(
                &base_prefix,
                output_dir,
                npz,
                sample_pairs,
                samples_to_keep,
                variant_indices,
                threads,
            )),
            None => Err(CustomError::InputFilesMissing {
                prefix: base_prefix.display().to_string(),
            }),
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
            InputSpec::Plink { bed, bim, fam, .. } => {
                println!("BED: {}", bed.display());
                println!("BIM: {}", bim.display());
                println!("FAM: {}", fam.display());
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
                samples_to_keep,
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
                    let reader = PackedAncestryMapReader::open(
                        ind,
                        geno,
                        snp,
                        samples_to_keep.clone(),
                        variant_indices.clone(),
                    )?;
                    Ok(Box::new(reader))
                } else if header_prefix.starts_with(b"TGENO") {
                    let reader = TransposedPackedAncestryMapReader::open(
                        ind,
                        geno,
                        snp,
                        samples_to_keep.clone(),
                        variant_indices.clone(),
                    )?;
                    Ok(Box::new(reader))
                } else {
                    let reader = EigenstratReader::open(
                        ind,
                        geno,
                        snp,
                        samples_to_keep.clone(),
                        variant_indices.clone(),
                    )?;
                    Ok(Box::new(reader))
                }
            }
            InputSpec::Plink {
                bed,
                bim,
                fam,
                variant_indices,
                samples_to_keep,
                ..
            } => {
                let reader = PlinkBedReader::open(
                    bed,
                    bim,
                    fam,
                    samples_to_keep.clone(),
                    variant_indices.clone(),
                )?;
                Ok(Box::new(reader))
            }
        }
    }

    pub fn output_dir(&self) -> &Path {
        match self {
            InputSpec::PackedAncestryMap { output_dir, .. }
            | InputSpec::Plink { output_dir, .. } => output_dir.as_path(),
        }
    }

    pub fn npz(&self) -> bool {
        match self {
            InputSpec::PackedAncestryMap { npz, .. } | InputSpec::Plink { npz, .. } => *npz,
        }
    }

    pub fn sample_pairs(&self) -> Option<&[(String, String)]> {
        match self {
            InputSpec::PackedAncestryMap { sample_pairs, .. }
            | InputSpec::Plink { sample_pairs, .. } => sample_pairs.as_deref(),
        }
    }

    pub fn threads(&self) -> Option<usize> {
        match self {
            InputSpec::PackedAncestryMap { threads, .. } | InputSpec::Plink { threads, .. } => {
                *threads
            }
        }
    }
}

fn path_with_extension(prefix: &Path, ext: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(".");
    path.push(ext);
    PathBuf::from(path)
}

fn normalize_prefix(prefix: &str) -> (PathBuf, Option<FormatHint>) {
    let path = PathBuf::from(prefix);
    let hint = match path.extension().and_then(|ext| ext.to_str()) {
        Some("geno") | Some("ind") | Some("snp") => Some(FormatHint::PackedAncestryMap),
        Some("bed") | Some("bim") | Some("fam") => Some(FormatHint::Plink),
        _ => None,
    };
    let base = match hint {
        Some(_) => strip_extension(&path),
        None => path,
    };
    (base, hint)
}

fn strip_extension(path: &Path) -> PathBuf {
    match (path.file_stem(), path.parent()) {
        (Some(stem), Some(parent)) => parent.join(Path::new(stem)),
        (Some(stem), None) => Path::new(stem).to_path_buf(),
        _ => path.to_path_buf(),
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
    let samples_to_keep = sample_pairs
        .as_ref()
        .map(|pairs| samples_to_keep_from_pairs(pairs));
    InputSpec::from_prefix(
        &args.prefix,
        &args.output_directory,
        args.npz,
        sample_pairs,
        samples_to_keep,
        variant_indices,
        args.threads,
    )
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

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum SampleCsvMode {
        Unknown,
        PairList,
        SampleList,
    }

    let mut mode = SampleCsvMode::Unknown;
    let mut pairs = Vec::new();
    let mut sample_ids = Vec::new();
    let mut seen_samples = HashSet::new();

    for result in reader.records() {
        let record = result.map_err(|source| CustomError::CsvRead {
            source,
            path: csv_path.clone(),
        })?;
        if record.iter().all(|field| field.trim().is_empty()) {
            continue;
        }

        match mode {
            SampleCsvMode::Unknown => match record.len() {
                1 => mode = SampleCsvMode::SampleList,
                2 => mode = SampleCsvMode::PairList,
                _ => return Err(CustomError::SamplePairsColumns),
            },
            SampleCsvMode::PairList => {
                if record.len() != 2 {
                    return Err(CustomError::SamplePairsColumns);
                }
            }
            SampleCsvMode::SampleList => {
                if record.len() != 1 {
                    return Err(CustomError::SamplePairsColumns);
                }
            }
        }

        match mode {
            SampleCsvMode::PairList => {
                let id1 = record[0].trim();
                let id2 = record[1].trim();
                if id1.is_empty() || id2.is_empty() {
                    continue;
                }
                pairs.push((id1.to_string(), id2.to_string()));
            }
            SampleCsvMode::SampleList => {
                let id = record[0].trim();
                if id.is_empty() {
                    continue;
                }
                if seen_samples.insert(id.to_string()) {
                    sample_ids.push(id.to_string());
                }
            }
            SampleCsvMode::Unknown => unreachable!(),
        }
    }

    match mode {
        SampleCsvMode::SampleList => {
            if sample_ids.len() < 2 {
                return Err(CustomError::SamplePairsEmpty);
            }
            let mut expanded = Vec::new();
            for i in 0..sample_ids.len() {
                for j in (i + 1)..sample_ids.len() {
                    expanded.push((sample_ids[i].to_string(), sample_ids[j].to_string()));
                }
            }
            Ok(expanded)
        }
        SampleCsvMode::PairList | SampleCsvMode::Unknown => {
            if pairs.is_empty() {
                Err(CustomError::SamplePairsEmpty)
            } else {
                Ok(pairs)
            }
        }
    }
}

fn samples_to_keep_from_pairs(pairs: &[(String, String)]) -> HashSet<String> {
    let mut keep = HashSet::new();
    for (left, right) in pairs {
        keep.insert(left.clone());
        keep.insert(right.clone());
    }
    keep
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
    // This gets the filtered sample set if a sample/sample pair CSV was provided
    let samples: Vec<String> = reader.samples().to_vec();
    // Resolve sample pairs to indices
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
    fn path_with_extension_preserves_multi_dot_prefix() {
        let prefix = std::path::Path::new("/tmp/data/v62.0_1240k_public");
        let expected = std::path::PathBuf::from("/tmp/data/v62.0_1240k_public.geno");
        assert_eq!(path_with_extension(prefix, "geno"), expected);
    }

    #[test]
    fn from_prefix_accepts_multi_dot_paths() {
        let unique_dir = format!(
            "fastpmr-prefix-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        );
        let base_dir = std::env::temp_dir().join(unique_dir);
        std::fs::create_dir_all(&base_dir).unwrap();

        let prefix = base_dir.join("dataset.v62.0");
        for ext in ["ind", "geno", "snp"] {
            let filename = format!("{}.{ext}", prefix.file_name().unwrap().to_string_lossy());
            let path = prefix.with_file_name(filename);
            std::fs::write(path, b"").unwrap();
        }
        let output_dir = base_dir.join("output");

        let spec = InputSpec::from_prefix(
            prefix.to_str().unwrap(),
            output_dir.to_str().unwrap(),
            false,
            None,
            None,
            None,
            None,
        )
        .expect("should detect packed ancestry map inputs");

        match spec {
            InputSpec::PackedAncestryMap { ind, geno, snp, .. } => {
                assert_eq!(
                    ind,
                    prefix.with_file_name(format!(
                        "{}.ind",
                        prefix.file_name().unwrap().to_string_lossy()
                    ))
                );
                assert_eq!(
                    geno,
                    prefix.with_file_name(format!(
                        "{}.geno",
                        prefix.file_name().unwrap().to_string_lossy()
                    ))
                );
                assert_eq!(
                    snp,
                    prefix.with_file_name(format!(
                        "{}.snp",
                        prefix.file_name().unwrap().to_string_lossy()
                    ))
                );
            }
            other => panic!("unexpected input spec: {other:?}"),
        }

        std::fs::remove_dir_all(base_dir).unwrap();
    }

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
