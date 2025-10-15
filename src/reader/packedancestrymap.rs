use itertools::Itertools;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::error::{CustomError, Result};
use crate::model::{Allele, Site};
use crate::reader::SiteReader;
use crate::reader::eigenstrat::{read_eigenstrat_ind, read_eigenstrat_snp};

const MIN_BLOCK_BYTES: usize = 48;
const GENO_HEADER_FIELDS: usize = 5;

pub struct PackedAncestryMapReader {
    reader: BufReader<File>,
    header: Header,
    samples: Vec<String>,
    variant_indices_to_keep: Option<HashSet<usize>>,
    next_variant_idx: usize,
    block_buf: Vec<u8>,
}

struct Header {
    n_samples: usize,
    n_variants: usize,
}

// See https://www.cog-genomics.org/plink/2.0/formats#geno for format description
impl PackedAncestryMapReader {
    pub fn open(
        ind_path: &impl AsRef<Path>,
        geno_path: &impl AsRef<Path>,
        snp_path: &impl AsRef<Path>,
        variant_indices_to_keep: Option<HashSet<usize>>,
    ) -> Result<Self> {
        let samples = read_eigenstrat_ind(ind_path)?;
        let variants = read_eigenstrat_snp(snp_path)?;
        let block_size = MIN_BLOCK_BYTES.max(samples.len().div_ceil(4));

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
        let header_block = &buffer[..block_size];
        if header_block.len() < MIN_BLOCK_BYTES {
            return Err(CustomError::PackedAncestryMapFileSize);
        }
        let header = parse_header_block(header_block)?;
        // Consume header block
        reader.consume(block_size);

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

        Ok(Self {
            reader,
            header,
            samples,
            variant_indices_to_keep,
            next_variant_idx: 0,
            block_buf: vec![0u8; block_size],
        })
    }
}

impl SiteReader for PackedAncestryMapReader {
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

impl Iterator for PackedAncestryMapReader {
    type Item = Result<Site>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.next_variant_idx < self.header.n_variants {
            let keep = match &self.variant_indices_to_keep {
                Some(set) => set.contains(&self.next_variant_idx),
                None => true,
            };

            if let Err(e) = self.reader.read_exact(&mut self.block_buf) {
                // Poison iterator to prevent further reads
                self.next_variant_idx = self.header.n_variants;
                return Some(Err(CustomError::ReadWithoutPath { source: e }));
            }

            self.next_variant_idx += 1;
            if keep {
                let genotypes = parse_variant_block(&self.block_buf, self.header.n_samples);
                return Some(Ok(Site { genotypes }));
            }
        }
        None
    }
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

    if fields[0] != "GENO" {
        return Err(CustomError::PackedAncestryMapHeaderGeno);
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

fn parse_variant_block(block: &[u8], n_samples: usize) -> Vec<Allele> {
    let bytes_needed = n_samples.div_ceil(4);
    let bytes = &block[..bytes_needed];
    let mut genotypes = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let byte = bytes[i / 4];
        // Extract the two bits of the i-th sample
        let shift = 6 - 2 * (i % 4);
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
