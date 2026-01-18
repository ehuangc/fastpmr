use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::{CustomError, Result};
use crate::model::{Allele, Site};
use crate::reader::SiteReader;
use crate::reader::common::{read_eigenstrat_ind, read_eigenstrat_snp, select_samples};

pub struct EigenstratReader {
    reader: BufReader<File>,
    n_samples: usize,
    n_variants: usize,
    samples: Vec<String>,
    sample_indices_to_keep: Option<Vec<usize>>,
    variant_indices_to_keep: Option<HashSet<usize>>,
    next_variant_idx: usize,
}

impl EigenstratReader {
    pub fn open(
        ind_path: &impl AsRef<Path>,
        geno_path: &impl AsRef<Path>,
        snp_path: &impl AsRef<Path>,
        samples_to_keep: Option<HashSet<String>>,
        variant_indices_to_keep: Option<HashSet<usize>>,
    ) -> Result<Self> {
        let samples = read_eigenstrat_ind(ind_path)?;
        let n_samples = samples.len();
        let variants = read_eigenstrat_snp(snp_path)?;
        let n_variants = variants.len();

        // Sanity-check samples
        if n_samples < 2 {
            return Err(CustomError::SampleCount { n_samples });
        }
        // Sanity-check variants
        if n_variants < 1 {
            return Err(CustomError::VariantCount { n_variants });
        }
        // Sanity-check variant indices to keep
        if let Some(set) = &variant_indices_to_keep
            && let Some(&bad_idx) = set.iter().find(|&&idx| idx >= n_variants)
        {
            return Err(CustomError::VariantIndexHigh {
                idx: bad_idx + 1,
                n_variants,
            });
        }

        let f = File::open(geno_path).map_err(|e| CustomError::ReadWithPath {
            source: e,
            path: geno_path.as_ref().to_path_buf(),
        })?;

        // Overwrite samples with filtered set and keep indices so we can slice later
        let (samples, sample_indices_to_keep) = select_samples(samples, samples_to_keep)?;

        Ok(Self {
            reader: BufReader::new(f),
            n_samples,
            n_variants,
            samples,
            sample_indices_to_keep,
            variant_indices_to_keep,
            next_variant_idx: 0,
        })
    }
}

impl SiteReader for EigenstratReader {
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

impl Iterator for EigenstratReader {
    type Item = Result<Site>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.next_variant_idx < self.n_variants {
            let keep = match &self.variant_indices_to_keep {
                Some(set) => set.contains(&self.next_variant_idx),
                None => true,
            };
            let line_num = self.next_variant_idx + 1;

            let mut line = String::new();
            if let Err(e) = self.reader.read_line(&mut line) {
                self.next_variant_idx = self.n_variants;
                return Some(Err(CustomError::ReadWithoutPath { source: e }));
            }
            if line.is_empty() {
                // Poison iterator to prevent further reads
                self.next_variant_idx = self.n_variants;
                return Some(Err(CustomError::EigenstratGenoVariantCount {
                    expected: self.n_variants,
                    found: self.next_variant_idx,
                }));
            }
            self.next_variant_idx += 1;

            if keep {
                let cleaned_line: String = line.chars().filter(|c| !c.is_whitespace()).collect();
                if cleaned_line.len() != self.n_samples {
                    return Some(Err(CustomError::EigenstratGenoFields {
                        line_num,
                        n_fields: cleaned_line.len(),
                        expected: self.n_samples,
                    }));
                }
                let genotypes = parse_variant_row(
                    &cleaned_line,
                    self.n_samples,
                    self.sample_indices_to_keep.as_deref(),
                );
                return Some(Ok(Site { genotypes }));
            };
        }
        None
    }
}

fn parse_variant_row(
    cleaned_line: &str,
    n_samples: usize,
    indices_to_keep: Option<&[usize]>,
) -> Vec<Allele> {
    let mut genotypes: Vec<Allele>;
    match indices_to_keep {
        Some(indices) => {
            genotypes = Vec::with_capacity(indices.len());
        }
        None => {
            genotypes = Vec::with_capacity(n_samples);
        }
    }

    for (idx, c) in cleaned_line.chars().enumerate() {
        let keep = match indices_to_keep {
            Some(indices) => indices.contains(&idx),
            None => true,
        };

        if keep {
            let allele = match c {
                '0' => Allele::Alt,
                '1' => Allele::Het,
                '2' => Allele::Ref,
                '9' => Allele::Missing,
                _ => unreachable!(),
            };
            genotypes.push(allele);
        }
    }
    genotypes
}
