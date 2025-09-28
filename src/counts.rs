use crate::error::Result;
use crate::model::{Allele, Site};
use crate::reader::SiteReader;
use indicatif::{ProgressBar, ProgressStyle};

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

    pub fn push_site(&mut self, site: &Site) {
        assert_eq!(
            self.n_samples,
            site.genotypes.len(),
            "site/sample size mismatch"
        );

        // Build a compact list of (index, allele) for non-missing calls.
        // This avoids branching from checking for missingness in the loop
        let present: Vec<(usize, Allele)> = site
            .genotypes
            .iter()
            .copied()
            .enumerate()
            .filter(|&(_, a)| a != Allele::Missing)
            .collect();

        for present_idx in 0..present.len() {
            let (i, genotype_i) = present[present_idx];
            for present_idx in (present_idx + 1)..present.len() {
                let (j, genotype_j) = present[present_idx];
                let counter_idx = self.idx(i, j);
                self.mismatches[counter_idx] += genotype_i.mismatch(genotype_j);
                // 2 alleles per site
                self.totals[counter_idx] += 2;
            }
        }
    }

    pub fn consume_reader<R: SiteReader>(mut self, reader: &mut R) -> Result<Self>
    where
        R: SiteReader + ?Sized, // Allow trait objects
    {
        let n_sites = reader.n_sites();
        let pb = ProgressBar::new(n_sites as u64);
        pb.set_style(
            ProgressStyle::with_template("[{elapsed_precise}] {bar:30} {pos}/{len} sites").unwrap(),
        );

        let mut n_processed_sites: u64 = 0;
        while let Some(site) = reader.next_site()? {
            self.push_site(&site);
            n_processed_sites += 1;
            pb.set_position(n_processed_sites);
        }
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
                    rates[idx] = -1.0;
                } else {
                    rates[idx] = self.mismatches[idx] as f32 / self.totals[idx] as f32;
                }
            }
        }
        (pairs, rates)
    }
}
