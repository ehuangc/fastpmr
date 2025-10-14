use crate::error::Result;
use crate::model::Allele;
use crate::reader::SiteReader;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};

pub struct Counts {
    samples: Vec<String>,
    n_samples: usize,
    // Note that there are up to 2 mismatches per site so totals = 2 * n_sites
    mismatches: Vec<u64>, // Flat (n x n) row-major
    totals: Vec<u64>,
}

impl Counts {
    pub fn new(samples: Vec<String>) -> Self {
        let n_samples = samples.len();
        Self {
            samples,
            n_samples,
            mismatches: vec![0; n_samples * n_samples],
            totals: vec![0; n_samples * n_samples],
        }
    }

    pub fn idx(&self, i: usize, j: usize) -> usize {
        if i < j {
            return self.n_samples * i + j;
        } else {
            return self.n_samples * j + i;
        }
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
                    self.mismatches[counter_idx] += genotype_i.mismatch(genotype_j) as u64;
                    self.totals[counter_idx] += 2; // Two alleles per site
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

        let mismatches: Vec<AtomicU64> = (0..self.n_samples * self.n_samples)
            .map(|_| AtomicU64::new(0))
            .collect();
        let totals: Vec<AtomicU64> = (0..self.n_samples * self.n_samples)
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
                    let idx = self.idx(sample_idx_i, sample_idx_j);
                    mismatches[idx]
                        .fetch_add(genotype_i.mismatch(genotype_j) as u64, Ordering::Relaxed);
                    totals[idx].fetch_add(2, Ordering::Relaxed);
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

    pub fn n_samples(&self) -> usize {
        self.n_samples
    }

    pub fn overlaps(&self) -> Vec<u64> {
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
}
