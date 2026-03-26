use crate::error::Result;
use crate::model::Allele;
use crate::reader::SiteReader;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use rayon::prelude::*;
use std::collections::HashSet;
use std::sync::atomic::{AtomicU64, Ordering};

pub struct ConfidenceIntervals {
    pub lower: Vec<f32>,
    pub upper: Vec<f32>,
}

pub struct Counts {
    samples: Vec<String>,
    n_samples: usize,
    // Note that there are up to 2 mismatches per site so totals = 2 * n_sites
    mismatches: Vec<u64>,   // Flat (n x n) row-major
    totals: Vec<u64>,       // Flat (n x n) row-major
    covered_snps: Vec<u64>, // Length n
    // If Some, only calculate PMRs for pairs where indices_to_count[idx(i, j)] is true
    indices_to_count: Option<Vec<bool>>,
}

impl Counts {
    pub fn new(
        samples: Vec<String>,
        pair_indices_to_count: Option<HashSet<(usize, usize)>>,
        covered_snps: Vec<u64>,
    ) -> Self {
        let n_samples = samples.len();
        let indices_to_count = pair_indices_to_count.map(|pairs| {
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
            covered_snps,
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
            .is_none_or(|mask| mask[self.idx(i, j)])
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
                        .is_none_or(|mask| mask[counter_idx])
                    {
                        self.mismatches[counter_idx] += genotype_i.mismatch(genotype_j);
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
                        .is_none_or(|mask| mask[counter_idx])
                    {
                        mismatches[counter_idx]
                            .fetch_add(genotype_i.mismatch(genotype_j), Ordering::Relaxed);
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

    pub fn covered_snps(&self) -> &[u64] {
        &self.covered_snps
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

    // Compute 95% confidence intervals for pairwise mismatch rates using the Wald method,
    // assuming pseudohaploid data where each site is one independent observation
    // CI = p̂ ± 1.96 × √(p̂(1 − p̂) / n), where n = number of site overlaps, clamped to [0, 1]
    pub fn confidence_intervals_95(&self) -> ConfidenceIntervals {
        let size = self.n_samples * self.n_samples;
        let site_overlaps = self.site_overlaps();
        let mut lower = vec![0.0f32; size];
        let mut upper = vec![0.0f32; size];

        for i in 0..self.n_samples {
            for j in (i + 1)..self.n_samples {
                let pair_idx = self.idx(i, j);
                if site_overlaps[pair_idx] == 0 {
                    lower[pair_idx] = f32::NAN;
                    upper[pair_idx] = f32::NAN;
                } else {
                    // For n, assume pseudohaploid data where each site is one independent observation
                    let n = site_overlaps[pair_idx] as f64;
                    let p = self.mismatches[pair_idx] as f64 / self.totals[pair_idx] as f64;
                    let se = (p * (1.0 - p) / n).sqrt();
                    lower[pair_idx] = (p - 1.96 * se).max(0.0) as f32;
                    upper[pair_idx] = (p + 1.96 * se).min(1.0) as f32;
                }
            }
        }
        ConfidenceIntervals { lower, upper }
    }
}

impl ConfidenceIntervals {
    pub fn lower_2d(&self, n_samples: usize, counts: &Counts) -> Array2<f32> {
        let mut matrix = Array2::from_elem((n_samples, n_samples), f32::NAN);
        for i in 0..n_samples {
            for j in (i + 1)..n_samples {
                if !counts.should_count_pair(i, j) {
                    continue;
                }
                let pair_idx = counts.idx(i, j);
                matrix[(i, j)] = self.lower[pair_idx];
                matrix[(j, i)] = self.lower[pair_idx];
            }
        }
        matrix
    }

    pub fn upper_2d(&self, n_samples: usize, counts: &Counts) -> Array2<f32> {
        let mut matrix = Array2::from_elem((n_samples, n_samples), f32::NAN);
        for i in 0..n_samples {
            for j in (i + 1)..n_samples {
                if !counts.should_count_pair(i, j) {
                    continue;
                }
                let pair_idx = counts.idx(i, j);
                matrix[(i, j)] = self.upper[pair_idx];
                matrix[(j, i)] = self.upper[pair_idx];
            }
        }
        matrix
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Site;
    use std::collections::HashSet;

    #[test]
    fn should_count_pair_respects_restrictions() {
        let mut indices_to_count = HashSet::new();
        indices_to_count.insert((0, 2));
        let counts = Counts::new(
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            Some(indices_to_count),
            vec![0; 3],
        );
        assert!(counts.should_count_pair(0, 2));
        assert!(counts.should_count_pair(2, 0));
        assert!(!counts.should_count_pair(0, 1));
        assert!(!counts.should_count_pair(1, 2));
    }

    #[test]
    fn should_count_pair_defaults_to_all_pairs() {
        let counts = Counts::new(vec!["A".to_string(), "B".to_string()], None, vec![0; 2]);
        assert!(counts.should_count_pair(0, 1));
        assert!(counts.should_count_pair(1, 0));
        assert!(!counts.should_count_pair(0, 0));
    }

    #[derive(Clone)]
    struct TestReader {
        samples: Vec<String>,
        sites: Vec<Vec<Allele>>,
        index: usize,
    }

    impl TestReader {
        fn new(samples: Vec<String>, sites: Vec<Vec<Allele>>) -> Self {
            Self {
                samples,
                sites,
                index: 0,
            }
        }
    }

    impl Iterator for TestReader {
        type Item = Result<Site>;

        fn next(&mut self) -> Option<Self::Item> {
            if self.index >= self.sites.len() {
                return None;
            }
            let genotypes = self.sites[self.index].clone();
            self.index += 1;
            Some(Ok(Site { genotypes }))
        }
    }

    impl SiteReader for TestReader {
        fn samples(&self) -> &[String] {
            &self.samples
        }

        fn n_sites(&self) -> usize {
            self.sites.len()
        }
    }

    fn build_test_reader() -> TestReader {
        let samples = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let sites = vec![
            vec![Allele::Ref, Allele::Het, Allele::Alt],
            vec![Allele::Missing, Allele::Het, Allele::Alt],
            vec![Allele::Ref, Allele::Ref, Allele::Missing],
        ];
        TestReader::new(samples, sites)
    }

    #[test]
    fn consume_reader_parallel_matches_serial() {
        let samples = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let mut serial_reader = build_test_reader();
        let serial = Counts::new(samples.clone(), None, vec![0; samples.len()])
            .consume_reader(&mut serial_reader)
            .expect("serial counts failed");

        let mut parallel_reader = build_test_reader();
        let parallel = Counts::new(samples.clone(), None, vec![0; samples.len()])
            .consume_reader_parallel(&mut parallel_reader)
            .expect("parallel counts failed");

        assert_eq!(serial.mismatches_2d(), parallel.mismatches_2d());
        assert_eq!(serial.totals_2d(), parallel.totals_2d());

        let mismatches = serial.mismatches_2d();
        let totals = serial.totals_2d();
        assert_eq!(totals[(0, 1)], 4);
        assert_eq!(mismatches[(0, 1)], 1);
        assert_eq!(totals[(0, 2)], 2);
        assert_eq!(mismatches[(0, 2)], 2);
        assert_eq!(totals[(1, 2)], 4);
        assert_eq!(mismatches[(1, 2)], 2);
    }

    #[test]
    fn confidence_intervals_95_basic() {
        let samples = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let mut reader = build_test_reader();
        let counts = Counts::new(samples.clone(), None, vec![0; samples.len()])
            .consume_reader(&mut reader)
            .expect("counts failed");

        let ci = counts.confidence_intervals_95();

        // Pair (A, B): mismatches=1, totals=4, site_overlaps=2, p=0.25
        // se = sqrt(0.25 * 0.75 / 2) = 0.30619
        // lower = max(0, 0.25 - 1.96 * 0.30619) = 0.0 (clamped)
        // upper = 0.25 + 1.96 * 0.30619 = 0.8501
        let pair_idx_ab = counts.idx(0, 1);
        assert!((ci.lower[pair_idx_ab] - 0.0).abs() < 1e-6);
        assert!((ci.upper[pair_idx_ab] - 0.8501).abs() < 1e-3);

        // Pair (A, C): mismatches=2, totals=2, site_overlaps=1, p=1.0
        // se = sqrt(1.0 * 0.0 / 1) = 0.0
        // lower = 1.0, upper = 1.0
        let pair_idx_ac = counts.idx(0, 2);
        assert!((ci.lower[pair_idx_ac] - 1.0).abs() < 1e-6);
        assert!((ci.upper[pair_idx_ac] - 1.0).abs() < 1e-6);

        // Pair (B, C): mismatches=2, totals=4, site_overlaps=2, p=0.5
        // se = sqrt(0.5 * 0.5 / 2) = 0.35355
        // lower = max(0, 0.5 - 1.96 * 0.35355) = 0.0 (clamped)
        // upper = min(1, 0.5 + 1.96 * 0.35355) = 1.0 (clamped)
        let pair_idx_bc = counts.idx(1, 2);
        assert!((ci.lower[pair_idx_bc] - 0.0).abs() < 1e-6);
        assert!((ci.upper[pair_idx_bc] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn confidence_intervals_95_nan_for_zero_totals() {
        let samples = vec!["A".to_string(), "B".to_string()];
        let counts = Counts::new(samples, None, vec![0; 2]);
        let ci = counts.confidence_intervals_95();
        let idx = counts.idx(0, 1);
        assert!(ci.lower[idx].is_nan());
        assert!(ci.upper[idx].is_nan());
    }
}
