use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::error::{CustomError, Result};
use crate::model::{Allele, Site};
use crate::reader::SiteReader;

const MIN_BLOCK_BYTES: usize = 48;
const GENO_HEADER_FIELDS: usize = 5;
const IND_FIELDS: usize = 3;
const SNP_FIELDS: usize = 6;

pub struct PackedAncestryMapReader {
    reader: BufReader<File>,
    header: Header,
    samples: Vec<String>,
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
    ) -> Result<Self> {
        let samples = read_ind(ind_path)?;
        let variants = read_snp(snp_path)?;
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
        let header = parse_header_block(buffer)?;
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

        Ok(Self {
            reader,
            header,
            samples,
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
        self.header.n_variants
    }
}

impl Iterator for PackedAncestryMapReader {
    type Item = Result<Site>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_variant_idx >= self.header.n_variants {
            return None;
        }

        match self.reader.read_exact(&mut self.block_buf) {
            Ok(()) => {
                let genotypes = parse_variant_block(&self.block_buf, self.header.n_samples);
                self.next_variant_idx += 1;
                Some(Ok(Site { genotypes }))
            }
            Err(e) => {
                // Poison iterator to prevent further reads
                self.next_variant_idx = self.header.n_variants;
                Some(Err(CustomError::ReadWithoutPath { source: e }))
            }
        }
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
