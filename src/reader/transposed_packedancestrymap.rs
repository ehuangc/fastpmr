use itertools::Itertools;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::error::{CustomError, Result};
use crate::model::{Allele, Site};
use crate::reader::SiteReader;

const HEADER_BLOCK_SIZE: usize = 48;
const GENO_HEADER_FIELDS: usize = 5;
const IND_FIELDS: usize = 3;
const SNP_FIELDS: usize = 6;

pub struct TransposedPackedAncestryMapReader {
    header: Header,
    samples: Vec<String>,
    sample_block_size: usize,
    variant_indices_to_keep: Option<HashSet<usize>>,
    next_variant_idx: usize,
    // Entire TGENO matrix (w/o header); length = n_samples * sample_block_size
    genotype_matrix: Vec<u8>,
}

struct Header {
    n_samples: usize,
    n_variants: usize,
}

// See https://www.cog-genomics.org/plink/2.0/formats#geno for format description
impl TransposedPackedAncestryMapReader {
    pub fn open(
        ind_path: &impl AsRef<Path>,
        geno_path: &impl AsRef<Path>,
        snp_path: &impl AsRef<Path>,
        variant_indices_to_keep: Option<HashSet<usize>>,
    ) -> Result<Self> {
        let samples = read_ind(ind_path)?;
        let variants = read_snp(snp_path)?;
        let sample_block_size = HEADER_BLOCK_SIZE.max(variants.len().div_ceil(4));

        let f = File::open(geno_path).map_err(|e| CustomError::ReadWithPath {
            source: e,
            path: geno_path.as_ref().to_path_buf(),
        })?;
        let mut reader = BufReader::new(f);

        // Read header block
        let buffer = reader.fill_buf().map_err(|e| CustomError::ReadWithPath {
            source: e,
            path: geno_path.as_ref().to_path_buf(),
        })?;
        let header_block = &buffer[..HEADER_BLOCK_SIZE];
        if header_block.len() < HEADER_BLOCK_SIZE {
            return Err(CustomError::PackedAncestryMapFileSize);
        }
        let header = parse_header_block(header_block)?;
        // Consume header block
        reader.consume(HEADER_BLOCK_SIZE);

        // Sanity-check samples
        if samples.len() != header.n_samples {
            return Err(CustomError::PackedAncestryMapNAgreement {
                n_header: header.n_samples,
                n_ind: samples.len(),
            });
        }
        if header.n_samples < 2 {
            return Err(CustomError::SampleCount {
                n_samples: header.n_samples,
            });
        }

        // Sanity-check variants
        if variants.len() != header.n_variants {
            return Err(CustomError::PackedAncestryMapVAgreement {
                n_header: header.n_variants,
                n_snp: variants.len(),
            });
        }
        if header.n_variants < 1 {
            return Err(CustomError::VariantCount {
                n_variants: header.n_variants,
            });
        }

        // Sanity-check variant indices to keep
        if let Some(set) = &variant_indices_to_keep {
            if let Some(&bad_idx) = set.iter().sorted().find(|&&idx| idx >= header.n_variants) {
                return Err(CustomError::VariantIndexHigh {
                    idx: bad_idx + 1,
                    n_variants: header.n_variants,
                });
            }
        }

        // Read entire matrix so we can iterate over sites efficiently
        let expected_bytes = header.n_samples * sample_block_size;
        let mut genotype_matrix = vec![0u8; expected_bytes];
        reader
            .read_exact(&mut genotype_matrix)
            .map_err(|e| CustomError::ReadWithPath {
                source: e,
                path: geno_path.as_ref().to_path_buf(),
            })?;

        // Ensure no trailing bytes
        let mut tmp = [0u8; 1];
        match reader.read(&mut tmp) {
            Ok(0) => {}
            Ok(_) => return Err(CustomError::PackedAncestryMapFileSize),
            Err(e) => {
                return Err(CustomError::ReadWithPath {
                    source: e,
                    path: geno_path.as_ref().to_path_buf(),
                });
            }
        }

        Ok(Self {
            header,
            samples,
            sample_block_size,
            variant_indices_to_keep,
            next_variant_idx: 0,
            genotype_matrix,
        })
    }

    fn genotypes_for_variant(&self, variant_idx: usize) -> Vec<Allele> {
        // For a given sample block, byte (variant_idx / 4) holds 4 genotypes (2 bits each)
        // We shift by 6, 4, 2, or 0 to get the relevant 2 right-most bits
        let byte_idx = variant_idx / 4;
        let shift = 6 - 2 * (variant_idx % 4);

        let n_samples = self.header.n_samples;
        let mut genotypes = Vec::with_capacity(n_samples);

        // Decode directly from the in-memory matrix
        for s in 0..n_samples {
            let sample_block_start = s * self.sample_block_size;
            let matrix_idx = sample_block_start + byte_idx;
            let byte = self.genotype_matrix[matrix_idx];
            let code = (byte >> shift) & 0b11;
            genotypes.push(match code {
                0b00 => Allele::Alt,
                0b01 => Allele::Het,
                0b10 => Allele::Ref,
                0b11 => Allele::Missing,
                _ => unreachable!(),
            });
        }
        genotypes
    }
}

impl SiteReader for TransposedPackedAncestryMapReader {
    fn samples(&self) -> &[String] {
        &self.samples
    }
    fn n_sites(&self) -> usize {
        if let Some(set) = &self.variant_indices_to_keep {
            set.len()
        } else {
            self.header.n_variants
        }
    }
}

impl Iterator for TransposedPackedAncestryMapReader {
    type Item = Result<Site>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.next_variant_idx < self.header.n_variants {
            let keep = match &self.variant_indices_to_keep {
                Some(set) => set.contains(&self.next_variant_idx),
                None => true,
            };

            let genotypes = self.genotypes_for_variant(self.next_variant_idx);
            self.next_variant_idx += 1;

            if keep {
                return Some(Ok(Site { genotypes }));
            }
        }
        None
    }
}

fn read_ind(path: &impl AsRef<Path>) -> Result<Vec<String>> {
    let f = File::open(path).map_err(|e| CustomError::ReadWithPath {
        source: e,
        path: path.as_ref().to_path_buf(),
    })?;
    let f = BufReader::new(f);
    let mut sample_ids: Vec<String> = Vec::new();

    for (line_idx, line) in f.lines().enumerate() {
        let line = line.map_err(|e| CustomError::ReadWithPath {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;
        let line = line.trim();
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != IND_FIELDS {
            return Err(CustomError::EigenstratIndFields {
                line_num: line_idx + 1,
                n_fields: fields.len(),
                expected: IND_FIELDS,
            });
        }
        let sample_id = fields[0].to_string();
        sample_ids.push(sample_id);
    }
    Ok(sample_ids)
}

fn read_snp(path: &impl AsRef<Path>) -> Result<Vec<String>> {
    let f = File::open(path).map_err(|e| CustomError::ReadWithPath {
        source: e,
        path: path.as_ref().to_path_buf(),
    })?;
    let f = BufReader::new(f);
    let mut variant_ids: Vec<String> = Vec::new();

    for (line_idx, line) in f.lines().enumerate() {
        let line = line.map_err(|e| CustomError::ReadWithPath {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;
        let line = line.trim();
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != SNP_FIELDS {
            return Err(CustomError::EigenstratSnpFields {
                line_num: line_idx + 1,
                n_fields: fields.len(),
                expected: SNP_FIELDS,
            });
        }
        let snp_id = fields[0].to_string();
        let chr = fields[1].to_string();
        variant_ids.push(format!("{}:{}:", chr, snp_id));
    }
    Ok(variant_ids)
}

fn parse_header_block(block: &[u8]) -> Result<Header> {
    let null_pos = block
        .iter()
        .position(|&b| b == b'\0')
        .ok_or_else(|| CustomError::PackedAncestryMapHeaderNullByte)?;
    let header_str = std::str::from_utf8(&block[..null_pos])
        .map_err(|e| CustomError::PackedAncestryMapHeaderUtf8 { source: e })?;

    let fields: Vec<_> = header_str.split_whitespace().collect();
    if fields.len() != GENO_HEADER_FIELDS {
        return Err(CustomError::PackedAncestryMapHeaderFields {
            n_fields: fields.len(),
            expected: GENO_HEADER_FIELDS,
        });
    }

    if fields[0] != "TGENO" {
        return Err(CustomError::PackedAncestryMapHeaderTgeno);
    }

    let n_samples = fields[1]
        .parse::<usize>()
        .map_err(|e| CustomError::PackedAncestryMapHeaderN { source: e })?;
    let n_variants = fields[2]
        .parse::<usize>()
        .map_err(|e| CustomError::PackedAncestryMapHeaderV { source: e })?;

    // TO-DO: Verify hashes
    // let hash_samples = fields[3].parse::<u32>()
    //     .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid sample hash: {e}")))?;
    // let hash_variants = fields[4].parse::<u32>()
    //     .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid variant hash: {e}")))?;

    Ok(Header {
        n_samples,
        n_variants,
    })
}
