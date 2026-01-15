use itertools::Itertools;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::error::{CustomError, Result};
use crate::model::{Allele, Site};
use crate::reader::SiteReader;
use crate::reader::sample_filter::select_samples;

const BED_MAGIC: [u8; 2] = [0x6c, 0x1b];
const BED_SNP_MAJOR: u8 = 0x01;
const BED_HEADER_LEN: usize = 3;
const FAM_FIELDS: usize = 6;
const BIM_FIELDS: usize = 6;

pub struct PlinkBedReader {
    reader: BufReader<File>,
    n_samples: usize,
    n_variants: usize,
    samples: Vec<String>,
    sample_indices_to_keep: Option<Vec<usize>>,
    variant_indices_to_keep: Option<HashSet<usize>>,
    next_variant_idx: usize,
    block_buf: Vec<u8>,
}

impl PlinkBedReader {
    pub fn open(
        bed_path: &impl AsRef<Path>,
        bim_path: &impl AsRef<Path>,
        fam_path: &impl AsRef<Path>,
        samples_to_keep: Option<HashSet<String>>,
        variant_indices_to_keep: Option<HashSet<usize>>,
    ) -> Result<Self> {
        let samples = read_plink_fam(fam_path)?;
        let n_variants = count_plink_bim(bim_path)?;
        let bytes_per_variant = samples.len().div_ceil(4);

        let f = File::open(bed_path).map_err(|e| CustomError::ReadWithPath {
            source: e,
            path: bed_path.as_ref().to_path_buf(),
        })?;
        let mut reader = BufReader::new(f);

        // Read header block
        let mut header = [0u8; BED_HEADER_LEN];
        reader
            .read_exact(&mut header)
            .map_err(|e| CustomError::ReadWithPath {
                source: e,
                path: bed_path.as_ref().to_path_buf(),
            })?;
        // Validate header bytes (two-byte magic + SNP-major mode flag)
        if header[..2] != BED_MAGIC {
            return Err(CustomError::PlinkBedHeaderMagic);
        }
        if header[2] != BED_SNP_MAJOR {
            return Err(CustomError::PlinkBedMode);
        }

        // Sanity-check samples
        if samples.len() < 2 {
            return Err(CustomError::SampleCount {
                n_samples: samples.len(),
            });
        }

        // Sanity-check variants
        if n_variants < 1 {
            return Err(CustomError::VariantCount { n_variants });
        }
        let expected_size = BED_HEADER_LEN as u64 + (bytes_per_variant as u64 * n_variants as u64);
        let actual_size = reader
            .get_ref()
            .metadata()
            .map_err(|e| CustomError::ReadWithPath {
                source: e,
                path: bed_path.as_ref().to_path_buf(),
            })?
            .len();
        if actual_size != expected_size {
            return Err(CustomError::PlinkBedFileSize {
                expected: expected_size,
                found: actual_size,
            });
        }

        // Sanity-check variant indices to keep
        if let Some(set) = &variant_indices_to_keep
            && let Some(&bad_idx) = set.iter().sorted().find(|&&idx| idx >= n_variants)
        {
            return Err(CustomError::VariantIndexHigh {
                idx: bad_idx + 1,
                n_variants,
            });
        }

        // Overwrite samples
        // Also record indices of the kept samples so we can read only the relevant genotypes later
        let (samples, sample_indices_to_keep) = select_samples(samples, samples_to_keep)?;

        Ok(Self {
            reader,
            n_samples: samples.len(),
            n_variants,
            samples,
            sample_indices_to_keep,
            variant_indices_to_keep,
            next_variant_idx: 0,
            block_buf: vec![0u8; bytes_per_variant],
        })
    }
}

impl SiteReader for PlinkBedReader {
    fn samples(&self) -> &[String] {
        &self.samples
    }

    fn n_sites(&self) -> usize {
        if let Some(set) = &self.variant_indices_to_keep {
            set.len()
        } else {
            self.n_variants
        }
    }
}

impl Iterator for PlinkBedReader {
    type Item = Result<Site>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.next_variant_idx < self.n_variants {
            let keep = match &self.variant_indices_to_keep {
                Some(set) => set.contains(&self.next_variant_idx),
                None => true,
            };

            if let Err(e) = self.reader.read_exact(&mut self.block_buf) {
                // Poison iterator to prevent further reads
                self.next_variant_idx = self.n_variants;
                return Some(Err(CustomError::ReadWithoutPath { source: e }));
            }

            self.next_variant_idx += 1;
            if keep {
                let genotypes = parse_variant_block(
                    &self.block_buf,
                    self.n_samples,
                    self.sample_indices_to_keep.as_deref(),
                );
                return Some(Ok(Site { genotypes }));
            }
        }
        None
    }
}

pub(super) fn read_plink_fam(path: &impl AsRef<Path>) -> Result<Vec<String>> {
    let f = File::open(path).map_err(|e| CustomError::ReadWithPath {
        source: e,
        path: path.as_ref().to_path_buf(),
    })?;
    let mut samples = Vec::new();
    for (line_idx, line) in BufReader::new(f).lines().enumerate() {
        let line = line.map_err(|e| CustomError::ReadWithPath {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() != FAM_FIELDS {
            return Err(CustomError::PlinkFamFields {
                line_num: line_idx + 1,
                n_fields: fields.len(),
                expected: FAM_FIELDS,
            });
        }
        let fid = fields[0];
        let iid = fields[1];
        let sample_id = if fid == "0" {
            iid.to_string()
        } else {
            format!("{fid}:{iid}")
        };
        samples.push(sample_id);
    }
    Ok(samples)
}

pub(super) fn count_plink_bim(path: &impl AsRef<Path>) -> Result<usize> {
    let f = File::open(path).map_err(|e| CustomError::ReadWithPath {
        source: e,
        path: path.as_ref().to_path_buf(),
    })?;
    let mut n_variants = 0usize;
    for (line_idx, line) in BufReader::new(f).lines().enumerate() {
        let line = line.map_err(|e| CustomError::ReadWithPath {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() != BIM_FIELDS {
            return Err(CustomError::PlinkBimFields {
                line_num: line_idx + 1,
                n_fields: fields.len(),
                expected: BIM_FIELDS,
            });
        }
        n_variants += 1;
    }
    Ok(n_variants)
}

fn parse_variant_block(
    block: &[u8],
    n_samples: usize,
    indices_to_keep: Option<&[usize]>,
) -> Vec<Allele> {
    match indices_to_keep {
        Some(indices) => {
            let mut genotypes = Vec::with_capacity(indices.len());
            for &sample_idx in indices {
                genotypes.push(decode_sample(block, sample_idx));
            }
            genotypes
        }
        None => {
            let mut genotypes = Vec::with_capacity(n_samples);
            for sample_idx in 0..n_samples {
                genotypes.push(decode_sample(block, sample_idx));
            }
            genotypes
        }
    }
}

fn decode_sample(bytes: &[u8], sample_idx: usize) -> Allele {
    // PLINK stores genotypes little-endian within the byte; sample 0 uses the lowest two bits
    let byte_idx = sample_idx / 4;
    let shift = (sample_idx % 4) * 2;
    let code = (bytes[byte_idx] >> shift) & 0b11;
    match code {
        0b00 => Allele::Ref,
        0b01 => Allele::Missing,
        0b10 => Allele::Het,
        0b11 => Allele::Alt,
        _ => unreachable!(),
    }
}
