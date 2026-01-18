use std::collections::HashSet;

use crate::error::{CustomError, Result};

/// Returns the subset of samples to keep, along with their original indices.
/// We keep indices so that we can read only the relevant genotypes later.
pub fn select_samples(
    samples: Vec<String>,
    filter: Option<HashSet<String>>,
) -> Result<(Vec<String>, Option<Vec<usize>>)> {
    match filter {
        Some(mut keep) => {
            let mut filtered_samples = Vec::with_capacity(keep.len());
            let mut indices = Vec::with_capacity(keep.len());
            for (idx, sample_id) in samples.into_iter().enumerate() {
                if keep.remove(&sample_id) {
                    indices.push(idx);
                    filtered_samples.push(sample_id);
                }
            }

            if !keep.is_empty() {
                let missing = keep.into_iter().next().unwrap();
                return Err(CustomError::SamplePairUnknownSample { sample: missing });
            }

            Ok((filtered_samples, Some(indices)))
        }
        None => Ok((samples, None)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn filters_requested_samples() {
        let samples = vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ];
        let filter = HashSet::from(["A".to_string(), "C".to_string()]);
        let (kept, indices) =
            select_samples(samples, Some(filter)).expect("filtering should succeed");
        assert_eq!(kept, vec!["A".to_string(), "C".to_string()]);
        assert_eq!(indices.unwrap(), vec![0, 2]);
    }

    #[test]
    fn errors_on_missing_sample() {
        let samples = vec!["A".to_string()];
        let filter = HashSet::from(["Z".to_string()]);
        let err = select_samples(samples, Some(filter)).unwrap_err();
        match err {
            CustomError::SamplePairUnknownSample { sample } => assert_eq!(sample, "Z"),
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
