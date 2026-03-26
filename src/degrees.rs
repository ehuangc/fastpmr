// Degree-calling procedure replicates that of READv2
// (Kuhn et al. 2018, https://doi.org/10.1186/s13059-024-03350-3)
use crate::counts::Counts;
use ndarray::Array2;

const DEG_THRESHOLDS: [f32; 4] = [0.625, 0.8125, 0.90625, 0.953125];
const THIRD_DEG_INFERENCE_EXPECTED_MISMATCHES_CUTOFF: f32 = 3000.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Degree {
    IdenticalTwin = 0,
    First = 1,
    Second = 2,
    Third = 3,
    Unrelated = 4,
}

impl std::fmt::Display for Degree {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Degree::IdenticalTwin => write!(f, "Identical/Twin"),
            Degree::First => write!(f, "First Degree"),
            Degree::Second => write!(f, "Second Degree"),
            Degree::Third => write!(f, "Third Degree"),
            Degree::Unrelated => write!(f, "Unrelated"),
        }
    }
}

impl Degree {
    pub const LABELS: [&str; 5] = [
        "Identical/Twin",
        "First Degree",
        "Second Degree",
        "Third Degree",
        "Unrelated",
    ];
}

pub struct DegreeResults {
    pub degrees: Vec<Degree>,
    pub normalized_mismatch_rates: Vec<f32>,
}

impl DegreeResults {
    pub fn normalized_mismatch_rates_2d(&self, n_samples: usize, counts: &Counts) -> Array2<f32> {
        let mut matrix = Array2::from_elem((n_samples, n_samples), f32::NAN);
        for i in 0..n_samples {
            for j in (i + 1)..n_samples {
                if !counts.should_count_pair(i, j) {
                    continue;
                }
                let pair_idx = counts.idx(i, j);
                matrix[(i, j)] = self.normalized_mismatch_rates[pair_idx];
                matrix[(j, i)] = self.normalized_mismatch_rates[pair_idx];
            }
        }
        matrix
    }

    pub fn degrees_2d(&self, n_samples: usize, counts: &Counts) -> Array2<u8> {
        let mut matrix = Array2::from_elem((n_samples, n_samples), Degree::Unrelated as u8);
        for i in 0..n_samples {
            for j in (i + 1)..n_samples {
                if !counts.should_count_pair(i, j) {
                    continue;
                }
                let pair_idx = counts.idx(i, j);
                matrix[(i, j)] = self.degrees[pair_idx] as u8;
                matrix[(j, i)] = self.degrees[pair_idx] as u8;
            }
        }
        matrix
    }
}

fn median(sorted: &[f32]) -> Option<f32> {
    if sorted.is_empty() {
        return None;
    }
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        Some((sorted[mid - 1] + sorted[mid]) / 2.0)
    } else {
        Some(sorted[mid])
    }
}

fn classify_pair(normalized_mismatch_rate: f32, expected_mismatches: f32) -> Degree {
    if normalized_mismatch_rate.is_nan() {
        return Degree::Unrelated;
    }
    if normalized_mismatch_rate < DEG_THRESHOLDS[0] {
        Degree::IdenticalTwin
    } else if normalized_mismatch_rate < DEG_THRESHOLDS[1] {
        Degree::First
    } else if normalized_mismatch_rate < DEG_THRESHOLDS[2] {
        Degree::Second
    } else if normalized_mismatch_rate <= DEG_THRESHOLDS[3]
        && expected_mismatches >= THIRD_DEG_INFERENCE_EXPECTED_MISMATCHES_CUTOFF
    {
        Degree::Third
    } else {
        Degree::Unrelated
    }
}

pub fn classify_degrees(counts: &Counts) -> DegreeResults {
    let n_samples = counts.n_samples();
    let size = n_samples * n_samples;
    let rates = counts.mismatch_rates();
    let overlaps = counts.site_overlaps();

    // Collect valid PMR values for median computation
    let mut valid_mismatch_rates: Vec<f32> = Vec::new();
    for i in 0..n_samples {
        for j in (i + 1)..n_samples {
            if !counts.should_count_pair(i, j) {
                continue;
            }
            let pair_idx = counts.idx(i, j);
            let rate = rates[pair_idx];
            if rate.is_finite() {
                valid_mismatch_rates.push(rate);
            }
        }
    }

    valid_mismatch_rates.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_mismatch_rate = median(&valid_mismatch_rates).unwrap_or(0.0);

    let mut degrees = vec![Degree::Unrelated; size];
    let mut normalized_mismatch_rates = vec![f32::NAN; size];

    if median_mismatch_rate == 0.0 {
        return DegreeResults {
            degrees,
            normalized_mismatch_rates,
        };
    }

    for i in 0..n_samples {
        for j in (i + 1)..n_samples {
            if !counts.should_count_pair(i, j) {
                continue;
            }
            let pair_idx = counts.idx(i, j);
            let rate = rates[pair_idx];
            let normalized_mismatch_rate = rate / median_mismatch_rate;
            let expected_mismatches = overlaps[pair_idx] as f32 * median_mismatch_rate;

            normalized_mismatch_rates[pair_idx] = normalized_mismatch_rate;
            degrees[pair_idx] = classify_pair(normalized_mismatch_rate, expected_mismatches);
        }
    }

    DegreeResults {
        degrees,
        normalized_mismatch_rates,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classify_pair_identical_twin() {
        assert_eq!(classify_pair(0.5, 5000.0), Degree::IdenticalTwin);
    }

    #[test]
    fn classify_pair_first_degree_at_lower_boundary() {
        assert_eq!(classify_pair(0.625, 5000.0), Degree::First);
    }

    #[test]
    fn classify_pair_first_degree() {
        assert_eq!(classify_pair(0.7, 5000.0), Degree::First);
    }

    #[test]
    fn classify_pair_second_degree_at_lower_boundary() {
        assert_eq!(classify_pair(0.8125, 5000.0), Degree::Second);
    }

    #[test]
    fn classify_pair_second_degree() {
        assert_eq!(classify_pair(0.85, 5000.0), Degree::Second);
    }

    #[test]
    fn classify_pair_third_degree_at_lower_boundary() {
        assert_eq!(classify_pair(0.90625, 5000.0), Degree::Third);
    }

    #[test]
    fn classify_pair_third_degree() {
        assert_eq!(classify_pair(0.93, 4000.0), Degree::Third);
    }

    #[test]
    fn classify_pair_third_degree_at_boundary() {
        assert_eq!(classify_pair(0.953125, 3000.0), Degree::Third);
    }

    #[test]
    fn classify_pair_third_degree_insufficient_snps() {
        assert_eq!(classify_pair(0.93, 2000.0), Degree::Unrelated);
    }

    #[test]
    fn classify_pair_third_degree_snp_cutoff_boundary() {
        assert_eq!(classify_pair(0.93, 3000.0), Degree::Third);
        assert_eq!(classify_pair(0.93, 2999.9), Degree::Unrelated);
    }

    #[test]
    fn classify_pair_unrelated() {
        assert_eq!(classify_pair(0.97, 5000.0), Degree::Unrelated);
    }

    #[test]
    fn classify_pair_nan() {
        assert_eq!(classify_pair(f32::NAN, 5000.0), Degree::Unrelated);
    }

    #[test]
    fn median_odd() {
        assert_eq!(median(&[1.0, 2.0, 3.0]), Some(2.0));
    }

    #[test]
    fn median_even() {
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), Some(2.5));
    }

    #[test]
    fn median_single() {
        assert_eq!(median(&[42.0]), Some(42.0));
    }

    #[test]
    fn median_empty() {
        assert_eq!(median(&[]), None);
    }

    #[test]
    fn degree_display() {
        assert_eq!(Degree::IdenticalTwin.to_string(), "Identical/Twin");
        assert_eq!(Degree::First.to_string(), "First Degree");
        assert_eq!(Degree::Second.to_string(), "Second Degree");
        assert_eq!(Degree::Third.to_string(), "Third Degree");
        assert_eq!(Degree::Unrelated.to_string(), "Unrelated");
    }

    #[test]
    fn degree_repr() {
        assert_eq!(Degree::IdenticalTwin as u8, 0);
        assert_eq!(Degree::First as u8, 1);
        assert_eq!(Degree::Second as u8, 2);
        assert_eq!(Degree::Third as u8, 3);
        assert_eq!(Degree::Unrelated as u8, 4);
    }
}
