use crate::error::Result;
use crate::model::Allele;
use crate::reader::SiteReader;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use rayon::prelude::*;
use std::collections::HashSet;
use std::sync::atomic::{AtomicU64, Ordering};

pub struct Counts {
    samples: Vec<String>,
    n_samples: usize,
    // Note that there are up to 2 mismatches per site so totals = 2 * n_sites
    mismatches: Vec<u64>, // Flat (n x n) row-major
    totals: Vec<u64>,
    // If Some, only calculate PMRs for pairs where indices_to_count[idx(i, j)] is true
    indices_to_count: Option<Vec<bool>>,
}

impl Counts {
    pub fn new(
        samples: Vec<String>,
        pairs_to_indices_to_count: Option<HashSet<(usize, usize)>>,
    ) -> Self {
        let n_samples = samples.len();
        let indices_to_count = pairs_to_indices_to_count.map(|pairs| {
            let mut mask = vec![false; n_samples * n_samples];
            for &(left, right) in &pairs {
                if left != right {
                    let idx = Self::idx_from_parts(n_samples, left, right);
                    mask[idx] = true;
                }
            }
            mask
        });
        Self {
            samples,
            n_samples,
            mismatches: vec![0; n_samples * n_samples],
            totals: vec![0; n_samples * n_samples],
            indices_to_count,
        }
    }

    pub fn idx(&self, i: usize, j: usize) -> usize {
        Self::idx_from_parts(self.n_samples, i, j)
    }

    fn idx_from_parts(n_samples: usize, i: usize, j: usize) -> usize {
        if i < j {
            n_samples * i + j
        } else {
            n_samples * j + i
        }
    }

    pub fn should_count_pair(&self, i: usize, j: usize) -> bool {
        if i == j {
            return false;
        }
        self.indices_to_count
            .as_ref()
            .map_or(true, |mask| mask[self.idx(i, j)])
    }

    pub fn consume_reader(mut self, reader: &mut dyn SiteReader) -> Result<Self> {
        let pb = ProgressBar::new(reader.n_sites() as u64);
        pb.set_style(
            ProgressStyle::with_template("[{elapsed_precise}] {bar:30} {pos}/{len} sites").unwrap(),
        );

        for site in reader {
            let site = site?;
            let present: Vec<(usize, Allele)> = site
                .genotypes
                .iter()
                .copied()
                .enumerate()
                .filter(|&(_, a)| a != Allele::Missing)
                .collect();

            for (i, &(sample_idx_i, genotype_i)) in present.iter().enumerate() {
                for &(sample_idx_j, genotype_j) in &present[i + 1..] {
                    let counter_idx = self.idx(sample_idx_i, sample_idx_j);
                    if self
                        .indices_to_count
                        .as_ref()
                        .map_or(true, |mask| mask[counter_idx])
                    {
                        self.mismatches[counter_idx] += genotype_i.mismatch(genotype_j) as u64;
                        self.totals[counter_idx] += 2; // Two alleles per site
                    }
                }
            }
            pb.inc(1);
        }
        pb.abandon();
        Ok(self)
    }

    pub fn consume_reader_parallel(mut self, reader: &mut dyn SiteReader) -> Result<Self> {
        let n_sites = reader.n_sites();
        let pb = ProgressBar::new(n_sites as u64);
        pb.set_style(
            ProgressStyle::with_template("[{elapsed_precise}] {bar:30} {pos}/{len} sites").unwrap(),
        );

        let n_samples = self.n_samples;
        let mismatches: Vec<AtomicU64> = (0..n_samples * n_samples)
            .map(|_| AtomicU64::new(0))
            .collect();
        let totals: Vec<AtomicU64> = (0..n_samples * n_samples)
            .map(|_| AtomicU64::new(0))
            .collect();

        reader.par_bridge().try_for_each(|site| -> Result<()> {
            let site = site?;
            let present: Vec<(usize, Allele)> = site
                .genotypes
                .iter()
                .copied()
                .enumerate()
                .filter(|&(_, a)| a != Allele::Missing)
                .collect();

            for (i, &(sample_idx_i, genotype_i)) in present.iter().enumerate() {
                for &(sample_idx_j, genotype_j) in &present[i + 1..] {
                    let counter_idx = self.idx(sample_idx_i, sample_idx_j);
                    if self
                        .indices_to_count
                        .as_ref()
                        .map_or(true, |mask| mask[counter_idx])
                    {
                        mismatches[counter_idx]
                            .fetch_add(genotype_i.mismatch(genotype_j) as u64, Ordering::Relaxed);
                        totals[counter_idx].fetch_add(2, Ordering::Relaxed);
                    }
                }
            }
            pb.inc(1);
            Ok(())
        })?;

        self.mismatches = mismatches
            .into_iter()
            .map(|x| x.load(Ordering::Relaxed))
            .collect();
        self.totals = totals
            .into_iter()
            .map(|x| x.load(Ordering::Relaxed))
            .collect();

        pb.abandon();
        println!();
        Ok(self)
    }

    pub fn samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    pub fn n_samples(&self) -> usize {
        self.n_samples
    }

    pub fn site_overlaps(&self) -> Vec<u64> {
        self.totals.iter().map(|x| x / 2).collect()
    }

    pub fn mismatch_rates(&self) -> (Vec<(String, String)>, Vec<f32>) {
        let mut pairs = vec![(String::new(), String::new()); self.n_samples * self.n_samples];
        let mut rates = vec![0.0; self.n_samples * self.n_samples];
        for i in 0..self.n_samples {
            for j in (i + 1)..self.n_samples {
                let idx = self.idx(i, j);
                pairs[idx] = (self.samples[i].clone(), self.samples[j].clone());
                if self.totals[idx] == 0 {
                    rates[idx] = f32::NAN;
                } else {
                    rates[idx] = self.mismatches[idx] as f32 / self.totals[idx] as f32;
                }
            }
        }
        (pairs, rates)
    }

    pub fn mismatches_2d(&self) -> Array2<u64> {
        let mut matrix = Array2::zeros((self.n_samples, self.n_samples));
        for i in 0..self.n_samples {
            for j in (i + 1)..self.n_samples {
                let idx = self.idx(i, j);
                matrix[(i, j)] = self.mismatches[idx];
                matrix[(j, i)] = self.mismatches[idx];
            }
        }
        matrix
    }

    pub fn totals_2d(&self) -> Array2<u64> {
        let mut matrix = Array2::zeros((self.n_samples, self.n_samples));
        for i in 0..self.n_samples {
            for j in (i + 1)..self.n_samples {
                let idx = self.idx(i, j);
                matrix[(i, j)] = self.totals[idx];
                matrix[(j, i)] = self.totals[idx];
            }
        }
        matrix
    }

    pub fn site_overlaps_2d(&self) -> Array2<u64> {
        let mut matrix = Array2::zeros((self.n_samples, self.n_samples));
        for i in 0..self.n_samples {
            for j in (i + 1)..self.n_samples {
                let idx = self.idx(i, j);
                matrix[(i, j)] = self.totals[idx] / 2;
                matrix[(j, i)] = self.totals[idx] / 2;
            }
        }
        matrix
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn should_count_pair_respects_restrictions() {
        let mut indices_to_count = HashSet::new();
        indices_to_count.insert((0, 2));
        let counts = Counts::new(
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            Some(indices_to_count),
        );
        assert!(counts.should_count_pair(0, 2));
        assert!(counts.should_count_pair(2, 0));
        assert!(!counts.should_count_pair(0, 1));
        assert!(!counts.should_count_pair(1, 2));
    }

    #[test]
    fn should_count_pair_defaults_to_all_pairs() {
        let counts = Counts::new(vec!["A".to_string(), "B".to_string()], None);
        assert!(counts.should_count_pair(0, 1));
        assert!(counts.should_count_pair(1, 0));
        assert!(!counts.should_count_pair(0, 0));
    }
}
