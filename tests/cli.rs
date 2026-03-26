mod common;

use ndarray::{Array1, Array2};
use ndarray_npy::NpzReader;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Read;
use std::path::Path;
use std::process::Command;
use zip::ZipArchive;

#[test]
fn packedancestrymap_cli_generates_outputs() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "packed").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(&dataset, RunOptions::default());
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records"
    );
}

#[test]
fn transposed_packedancestrymap_cli_generates_outputs() {
    let dataset = common::create_dataset(common::GenoFormat::Transposed, "transposed").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(&dataset, RunOptions::default());
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records"
    );
}

#[test]
fn eigenstrat_cli_generates_outputs() {
    let dataset = common::create_dataset(common::GenoFormat::Eigenstrat, "eigenstrat").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(&dataset, RunOptions::default());
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records"
    );
}

#[test]
fn plink_cli_generates_outputs() {
    let dataset = common::create_dataset(common::GenoFormat::Plink, "plink").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(&dataset, RunOptions::default());
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records"
    );
}

#[test]
fn packedancestrymap_cli_generates_outputs_with_threads() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "packed-threads").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            threads: Some(2),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records"
    );
}

#[test]
fn packedancestrymap_cli_generates_npz_outputs() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "packed-npz").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            npz: true,
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_npz_outputs(
        &dataset.output_dir,
        &expected_pairs,
        &expected_coverage,
        false,
    );
}

#[test]
fn csv_outputs_include_degrees() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "packed-degrees-csv").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            degrees: true,
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records"
    );
    for record in &records {
        assert!(
            record.degree.is_some(),
            "expected degree column when --degrees is set"
        );
        assert!(
            record.normalized_mismatch_rate.is_some(),
            "expected normalized_mismatch_rate column when --degrees is set"
        );
    }
}

#[test]
fn npz_outputs_include_degrees() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "packed-degrees-npz").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            npz: true,
            degrees: true,
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_npz_outputs(
        &dataset.output_dir,
        &expected_pairs,
        &expected_coverage,
        true,
    );
}

fn expected_ci(stats: &common::PairStats) -> (f32, f32) {
    if stats.totals == 0 {
        return (f32::NAN, f32::NAN);
    }
    let n = (stats.totals / 2) as f64;
    let p = stats.mismatches as f64 / stats.totals as f64;
    let se = (p * (1.0 - p) / n).sqrt();
    let lower = (p - 1.96 * se).max(0.0) as f32;
    let upper = (p + 1.96 * se).min(1.0) as f32;
    (lower, upper)
}

#[test]
fn csv_outputs_include_ci() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "packed-ci-csv").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            ci: true,
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    for record in &records {
        let lo = record
            .ci_lower
            .expect("expected ci_lower column when --ci is set");
        let hi = record
            .ci_upper
            .expect("expected ci_upper column when --ci is set");
        let key = (record.id1.clone(), record.id2.clone());
        let stats = &expected_pairs[&key];
        let (expected_lo, expected_hi) = expected_ci(stats);
        assert!(
            (lo - expected_lo).abs() < 1e-6,
            "ci_lower mismatch for {} / {}: got {}, expected {}",
            record.id1,
            record.id2,
            lo,
            expected_lo
        );
        assert!(
            (hi - expected_hi).abs() < 1e-6,
            "ci_upper mismatch for {} / {}: got {}, expected {}",
            record.id1,
            record.id2,
            hi,
            expected_hi
        );
    }
}

#[test]
fn npz_outputs_include_ci() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "packed-ci-npz").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            npz: true,
            ci: true,
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_npz_outputs(
        &dataset.output_dir,
        &expected_pairs,
        &expected_coverage,
        false,
    );

    let npz_path = dataset.output_dir.join("mismatch_counts.npz");
    let mut npz = NpzReader::new(File::open(&npz_path).expect("could not open npz"))
        .expect("invalid npz archive");
    let ci_lower: Array2<f32> = npz
        .by_name("mismatch_rate_95_ci_lower")
        .expect("missing mismatch_rate_95_ci_lower array");
    let ci_upper: Array2<f32> = npz
        .by_name("mismatch_rate_95_ci_upper")
        .expect("missing mismatch_rate_95_ci_upper array");
    let expected_shape = &[common::N_SAMPLES, common::N_SAMPLES];
    assert_eq!(ci_lower.shape(), expected_shape);
    assert_eq!(ci_upper.shape(), expected_shape);

    let samples = common::expected_sample_ids();
    for i in 0..common::N_SAMPLES {
        for j in (i + 1)..common::N_SAMPLES {
            assert!(
                (ci_lower[[i, j]] - ci_lower[[j, i]]).abs() < 1e-10
                    || (ci_lower[[i, j]].is_nan() && ci_lower[[j, i]].is_nan()),
                "ci_lower not symmetric at ({i}, {j})"
            );
            assert!(
                (ci_upper[[i, j]] - ci_upper[[j, i]]).abs() < 1e-10
                    || (ci_upper[[i, j]].is_nan() && ci_upper[[j, i]].is_nan()),
                "ci_upper not symmetric at ({i}, {j})"
            );
            let key = (samples[i].clone(), samples[j].clone());
            let stats = &expected_pairs[&key];
            let (exp_lo, exp_hi) = expected_ci(stats);
            assert!(
                (ci_lower[[i, j]] - exp_lo).abs() < 1e-6,
                "ci_lower mismatch at ({i}, {j}): got {}, expected {}",
                ci_lower[[i, j]],
                exp_lo
            );
            assert!(
                (ci_upper[[i, j]] - exp_hi).abs() < 1e-6,
                "ci_upper mismatch at ({i}, {j}): got {}, expected {}",
                ci_upper[[i, j]],
                exp_hi
            );
        }
    }
}

#[test]
fn variant_indices_limit_sites() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "filtered").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_filtered_variants();
    let expected_coverage = common::expected_covered_snps_filtered_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            variant_spec: Some("1-30000"),
            min_covered_snps: Some(0),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records after filtering"
    );
}

#[test]
fn eigenstrat_variant_indices_limit_sites() {
    let dataset =
        common::create_dataset(common::GenoFormat::Eigenstrat, "eigenstrat-filtered").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_filtered_variants();
    let expected_coverage = common::expected_covered_snps_filtered_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            variant_spec: Some("1-30000"),
            min_covered_snps: Some(0),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records after filtering EIGENSTRAT input"
    );
}

#[test]
fn plink_variant_indices_limit_sites() {
    let dataset = common::create_dataset(common::GenoFormat::Plink, "plink-filtered").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_filtered_variants();
    let expected_coverage = common::expected_covered_snps_filtered_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            variant_spec: Some("1-30000"),
            min_covered_snps: Some(0),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records after filtering PLINK input"
    );
}

#[test]
fn sample_pairs_csv_runs_successfully() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "sample-pairs-ok").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }
    let csv_path = dataset.prefix.with_extension("pairs.csv");
    fs::write(
        &csv_path,
        "Sample1,Sample2\nSample3,Sample1\nSample4,Sample2\n",
    )
    .unwrap();

    let output = run_fastpmr(
        &dataset,
        RunOptions {
            sample_pairs_csv: Some(&csv_path),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let expected_coverage = common::expected_covered_snps_all_variants();
    let expected_subset: BTreeMap<_, _> = common::expected_pair_stats_all_variants()
        .into_iter()
        .filter(|((id1, id2), _)| {
            matches!(
                (id1.as_str(), id2.as_str()),
                ("Sample1", "Sample2") | ("Sample1", "Sample3") | ("Sample2", "Sample4")
            )
        })
        .collect();

    let records = assert_outputs(&dataset.output_dir, &expected_subset, &expected_coverage);
    assert_eq!(records.len(), 3, "expected exactly three requested pairs");
    let pairs: std::collections::HashSet<(String, String)> = records
        .iter()
        .map(|record| (record.id1.clone(), record.id2.clone()))
        .collect();
    let expected_pairs: std::collections::HashSet<(String, String)> = [
        ("Sample1".to_string(), "Sample2".to_string()),
        ("Sample1".to_string(), "Sample3".to_string()),
        ("Sample2".to_string(), "Sample4".to_string()),
    ]
    .into_iter()
    .collect();
    assert_eq!(pairs, expected_pairs);
}

#[test]
fn eigenstrat_sample_pairs_csv_runs_successfully() {
    let dataset =
        common::create_dataset(common::GenoFormat::Eigenstrat, "eigenstrat-sample-pairs-ok")
            .unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }
    let csv_path = dataset.prefix.with_extension("pairs.csv");
    fs::write(
        &csv_path,
        "Sample1,Sample2\nSample3,Sample1\nSample4,Sample2\n",
    )
    .unwrap();

    let output = run_fastpmr(
        &dataset,
        RunOptions {
            sample_pairs_csv: Some(&csv_path),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let expected_coverage = common::expected_covered_snps_all_variants();
    let expected_subset: BTreeMap<_, _> = common::expected_pair_stats_all_variants()
        .into_iter()
        .filter(|((id1, id2), _)| {
            matches!(
                (id1.as_str(), id2.as_str()),
                ("Sample1", "Sample2") | ("Sample1", "Sample3") | ("Sample2", "Sample4")
            )
        })
        .collect();

    let records = assert_outputs(&dataset.output_dir, &expected_subset, &expected_coverage);
    assert_eq!(records.len(), 3, "expected exactly three requested pairs");
    let pairs: std::collections::HashSet<(String, String)> = records
        .iter()
        .map(|record| (record.id1.clone(), record.id2.clone()))
        .collect();
    let expected_pairs: std::collections::HashSet<(String, String)> = [
        ("Sample1".to_string(), "Sample2".to_string()),
        ("Sample1".to_string(), "Sample3".to_string()),
        ("Sample2".to_string(), "Sample4".to_string()),
    ]
    .into_iter()
    .collect();
    assert_eq!(pairs, expected_pairs);
}

#[test]
fn plink_sample_pairs_csv_runs_successfully() {
    let dataset =
        common::create_dataset(common::GenoFormat::Plink, "plink-sample-pairs-ok").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }
    let csv_path = dataset.prefix.with_extension("pairs.csv");
    fs::write(
        &csv_path,
        "Sample1,Sample2\nSample3,Sample1\nSample4,Sample2\n",
    )
    .unwrap();

    let output = run_fastpmr(
        &dataset,
        RunOptions {
            sample_pairs_csv: Some(&csv_path),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let expected_coverage = common::expected_covered_snps_all_variants();
    let expected_subset: BTreeMap<_, _> = common::expected_pair_stats_all_variants()
        .into_iter()
        .filter(|((id1, id2), _)| {
            matches!(
                (id1.as_str(), id2.as_str()),
                ("Sample1", "Sample2") | ("Sample1", "Sample3") | ("Sample2", "Sample4")
            )
        })
        .collect();

    let records = assert_outputs(&dataset.output_dir, &expected_subset, &expected_coverage);
    assert_eq!(records.len(), 3, "expected exactly three requested pairs");
    let pairs: std::collections::HashSet<(String, String)> = records
        .iter()
        .map(|record| (record.id1.clone(), record.id2.clone()))
        .collect();
    let expected_pairs: std::collections::HashSet<(String, String)> = [
        ("Sample1".to_string(), "Sample2".to_string()),
        ("Sample1".to_string(), "Sample3".to_string()),
        ("Sample2".to_string(), "Sample4".to_string()),
    ]
    .into_iter()
    .collect();
    assert_eq!(pairs, expected_pairs);
}

#[test]
fn sample_list_csv_runs_successfully() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "sample-list-ok").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }
    let csv_path = dataset.prefix.with_extension("pairs.csv");
    fs::write(&csv_path, "Sample1\nSample3\nSample4\n").unwrap();

    let output = run_fastpmr(
        &dataset,
        RunOptions {
            sample_pairs_csv: Some(&csv_path),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let mut expected_coverage = common::expected_covered_snps_all_variants();
    expected_coverage.retain(|id, _| matches!(id.as_str(), "Sample1" | "Sample3" | "Sample4"));
    let expected_subset: BTreeMap<_, _> = common::expected_pair_stats_all_variants()
        .into_iter()
        .filter(|((id1, id2), _)| {
            matches!(id1.as_str(), "Sample1" | "Sample3" | "Sample4")
                && matches!(id2.as_str(), "Sample1" | "Sample3" | "Sample4")
        })
        .collect();
    assert_eq!(expected_subset.len(), 3, "expected exactly three pairs");

    let records = assert_outputs(&dataset.output_dir, &expected_subset, &expected_coverage);
    let pairs: std::collections::HashSet<(String, String)> = records
        .iter()
        .map(|record| (record.id1.clone(), record.id2.clone()))
        .collect();
    let expected_pairs: std::collections::HashSet<(String, String)> = [
        ("Sample1".to_string(), "Sample3".to_string()),
        ("Sample1".to_string(), "Sample4".to_string()),
        ("Sample3".to_string(), "Sample4".to_string()),
    ]
    .into_iter()
    .collect();
    assert_eq!(pairs, expected_pairs);
}

#[test]
fn sample_pairs_csv_with_unknown_sample_fails() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "sample-pairs-err").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }
    let csv_path = dataset.prefix.with_extension("pairs.csv");
    fs::write(&csv_path, "Sample1,Unknown\n").unwrap();

    let output = run_fastpmr(
        &dataset,
        RunOptions {
            sample_pairs_csv: Some(&csv_path),
            ..Default::default()
        },
    );
    assert!(
        !output.status.success(),
        "fastpmr unexpectedly succeeded: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Unknown"),
        "stderr did not contain expected sample name: {stderr}"
    );
}

#[test]
fn chromosomes_filter_sites_packed() {
    let dataset =
        common::create_multichrom_dataset(common::GenoFormat::Packed, "packed-chr").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    // Chromosome 1 holds the first CORE_VARIANTS variants; filtering to chr 1 should give the
    // same results as --variant-indices 1-30000.
    let expected_pairs = common::expected_pair_stats_filtered_variants();
    let expected_coverage = common::expected_covered_snps_filtered_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            chromosomes_spec: Some("1"),
            min_covered_snps: Some(0),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records after chromosome filtering"
    );
}

#[test]
fn chromosomes_filter_sites_eigenstrat() {
    let dataset =
        common::create_multichrom_dataset(common::GenoFormat::Eigenstrat, "eigenstrat-chr")
            .unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_filtered_variants();
    let expected_coverage = common::expected_covered_snps_filtered_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            chromosomes_spec: Some("1"),
            min_covered_snps: Some(0),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records after chromosome filtering (EIGENSTRAT)"
    );
}

#[test]
fn chromosomes_filter_sites_plink() {
    let dataset =
        common::create_multichrom_dataset(common::GenoFormat::Plink, "plink-chr").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_filtered_variants();
    let expected_coverage = common::expected_covered_snps_filtered_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            chromosomes_spec: Some("1"),
            min_covered_snps: Some(0),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records after chromosome filtering (PLINK)"
    );
}

#[test]
fn chromosomes_range_includes_all() {
    // --chromosomes 1-2 on a dataset split across chr 1 and chr 2 should include all variants.
    let dataset =
        common::create_multichrom_dataset(common::GenoFormat::Packed, "packed-chr-all").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_all_variants();
    let expected_coverage = common::expected_covered_snps_all_variants();
    let output = run_fastpmr(
        &dataset,
        RunOptions {
            chromosomes_spec: Some("1-2"),
            min_covered_snps: Some(0),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs, &expected_coverage);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records with all chromosomes selected"
    );
}

#[test]
fn chromosomes_intersects_with_variant_indices() {
    // --chromosomes 1 combined with --variant-indices 1-15000 should intersect: only the first
    // 15000 variants on chr 1. We just check the command succeeds and produces fewer overlaps
    // than either filter alone.
    let dataset =
        common::create_multichrom_dataset(common::GenoFormat::Packed, "packed-chr-vi").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let output = run_fastpmr(
        &dataset,
        RunOptions {
            variant_spec: Some("1-15000"),
            chromosomes_spec: Some("1"),
            min_covered_snps: Some(0),
            ..Default::default()
        },
    );
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify that the overlap counts are less than the chr-1-only case (30000 sites).
    let csv_path = dataset.output_dir.join("mismatch_rates.csv");
    let content = fs::read_to_string(&csv_path).expect("could not read mismatch_rates.csv");
    let overlaps: Vec<u64> = content
        .lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let mut fields = l.split(',');
            fields.next(); // id1
            fields.next(); // id2
            fields.next().unwrap().parse::<u64>().unwrap() // n_site_overlap
        })
        .collect();
    assert!(!overlaps.is_empty(), "no pairs in output");
    for overlap in &overlaps {
        assert!(
            *overlap <= 15000,
            "overlap {overlap} exceeds variant-indices bound of 15000"
        );
    }
}

#[derive(Default)]
struct RunOptions<'a> {
    variant_spec: Option<&'a str>,
    chromosomes_spec: Option<&'a str>,
    npz: bool,
    degrees: bool,
    ci: bool,
    sample_pairs_csv: Option<&'a Path>,
    min_covered_snps: Option<u64>,
    threads: Option<usize>,
}

fn run_fastpmr(dataset: &common::Dataset, opts: RunOptions) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_fastpmr"));
    command
        .arg("--prefix")
        .arg(dataset.prefix.as_os_str())
        .arg("--output-directory")
        .arg(dataset.output_dir.as_os_str());
    if let Some(spec) = opts.variant_spec {
        command.arg("--variant-indices").arg(spec);
    }
    if let Some(spec) = opts.chromosomes_spec {
        command.arg("--chromosomes").arg(spec);
    }
    if opts.npz {
        command.arg("--npz");
    }
    if opts.degrees {
        command.arg("--degrees");
    }
    if opts.ci {
        command.arg("--ci");
    }
    if let Some(path) = opts.sample_pairs_csv {
        command.arg("--sample-pairs-csv").arg(path);
    }
    if let Some(min) = opts.min_covered_snps {
        command.arg("--min-covered-snps").arg(min.to_string());
    }
    if let Some(n) = opts.threads {
        command.arg("--threads").arg(n.to_string());
    }
    command.output().expect("failed to run fastpmr")
}

#[derive(Clone, Debug, PartialEq)]
struct OutputRecord {
    id1: String,
    id2: String,
    overlap: u64,
    rate: f32,
    ci_lower: Option<f32>,
    ci_upper: Option<f32>,
    normalized_mismatch_rate: Option<f64>,
    degree: Option<String>,
}

fn assert_outputs(
    output_dir: &Path,
    expected_pairs: &BTreeMap<(String, String), common::PairStats>,
    expected_coverage: &BTreeMap<String, u64>,
) -> Vec<OutputRecord> {
    let csv_path = output_dir.join("mismatch_rates.csv");
    let records = read_records(&csv_path);
    assert!(
        !records.is_empty(),
        "expected at least one mismatch rate record in {}",
        csv_path.display()
    );

    for record in &records {
        let key = (record.id1.clone(), record.id2.clone());
        let expectation = expected_pairs.get(&key).unwrap_or_else(|| {
            panic!(
                "unexpected mismatch rate record for {} / {}",
                record.id1, record.id2
            )
        });
        assert_eq!(
            record.overlap,
            expectation.overlap(),
            "unexpected overlap for {} / {}: got {}, expected {}",
            record.id1,
            record.id2,
            record.overlap,
            expectation.overlap()
        );
        let expected_rate = expectation.mismatch_rate();
        assert!(
            (record.rate - expected_rate).abs() < 1e-6,
            "unexpected mismatch rate for {} / {}: got {}, expected {expected_rate}",
            record.id1,
            record.id2,
            record.rate
        );
    }

    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "mismatch rate file contained {} pairs but expected {}",
        records.len(),
        expected_pairs.len()
    );

    assert_covered_snps(output_dir, expected_coverage);

    let plot_path = output_dir.join("mismatch_rates.png");
    let metadata = fs::metadata(plot_path).expect("missing plot output");
    assert!(metadata.len() > 0, "plot output is empty");
    records
}

fn assert_npz_outputs(
    output_dir: &Path,
    expected_pairs: &BTreeMap<(String, String), common::PairStats>,
    expected_coverage: &BTreeMap<String, u64>,
    degrees: bool,
) {
    let csv_path = output_dir.join("mismatch_rates.csv");
    assert!(
        !csv_path.exists(),
        "unexpected mismatch_rates.csv output alongside npz"
    );

    let coverage_path = output_dir.join("covered_snps.csv");
    assert!(
        !coverage_path.exists(),
        "unexpected covered_snps.csv output alongside npz"
    );

    let npz_path = output_dir.join("mismatch_counts.npz");
    assert!(
        npz_path.exists(),
        "missing mismatch_counts.npz output: looked at {}",
        npz_path.display()
    );

    let mut npz =
        NpzReader::new(File::open(&npz_path).expect("could not open mismatch_counts.npz"))
            .expect("invalid npz archive");
    let overlaps: Array2<u64> = npz
        .by_name("site_overlaps")
        .expect("missing site_overlaps array");
    let mismatch_rates: Array2<f32> = npz
        .by_name("mismatch_rates")
        .expect("missing mismatch_rates array");
    let covered_snps: Array1<u64> = npz
        .by_name("covered_snps")
        .expect("missing covered_snps array");
    let expected_shape = &[common::N_SAMPLES, common::N_SAMPLES];
    assert_eq!(overlaps.shape(), expected_shape);
    assert_eq!(mismatch_rates.shape(), expected_shape);
    assert_eq!(covered_snps.len(), common::N_SAMPLES);

    for idx in 0..common::N_SAMPLES {
        assert_eq!(overlaps[[idx, idx]], 0);
        assert!(mismatch_rates[[idx, idx]].is_nan());
    }

    let samples = common::expected_sample_ids();
    for (idx, sample) in samples.iter().enumerate() {
        let expected = expected_coverage
            .get(sample)
            .unwrap_or_else(|| panic!("missing coverage expectation for {sample}"));
        assert_eq!(
            covered_snps[idx], *expected,
            "unexpected covered SNPs for {sample}"
        );
    }
    let mut validated_pairs = 0usize;
    for i in 0..common::N_SAMPLES {
        for j in (i + 1)..common::N_SAMPLES {
            let key = (samples[i].clone(), samples[j].clone());
            let expectation = expected_pairs.get(&key).unwrap_or_else(|| {
                panic!(
                    "missing expectation for pair {} / {}",
                    samples[i], samples[j]
                )
            });
            assert!(
                (mismatch_rates[[i, j]] - expectation.mismatch_rate()).abs() < 1e-6,
                "unexpected mismatch rate for {} / {}: got {} expected {}",
                samples[i],
                samples[j],
                mismatch_rates[[i, j]],
                expectation.mismatch_rate()
            );
            assert!(
                (mismatch_rates[[j, i]] - expectation.mismatch_rate()).abs() < 1e-6,
                "unexpected symmetric mismatch rate for {} / {}: got {} expected {}",
                samples[j],
                samples[i],
                mismatch_rates[[j, i]],
                expectation.mismatch_rate()
            );
            assert_eq!(
                overlaps[[i, j]],
                expectation.overlap(),
                "unexpected overlaps for {} / {}",
                samples[i],
                samples[j]
            );
            assert_eq!(
                overlaps[[j, i]],
                expectation.overlap(),
                "unexpected symmetric overlaps for {} / {}",
                samples[j],
                samples[i]
            );
            validated_pairs += 1;
        }
    }
    assert_eq!(
        validated_pairs,
        expected_pairs.len(),
        "validated {} pairs but expected {}",
        validated_pairs,
        expected_pairs.len()
    );

    if degrees {
        let normalized: Array2<f32> = npz
            .by_name("normalized_mismatch_rates")
            .expect("missing normalized_mismatch_rates array");
        let degree_codes: Array2<u8> = npz.by_name("degrees").expect("missing degrees array");
        assert_eq!(normalized.shape(), expected_shape);
        assert_eq!(degree_codes.shape(), expected_shape);

        for i in 0..common::N_SAMPLES {
            for j in (i + 1)..common::N_SAMPLES {
                assert_eq!(
                    degree_codes[[i, j]],
                    degree_codes[[j, i]],
                    "degree matrix not symmetric at ({i}, {j})"
                );
                assert!(
                    degree_codes[[i, j]] <= 4,
                    "invalid degree code {} at ({i}, {j})",
                    degree_codes[[i, j]]
                );
                assert!(
                    (normalized[[i, j]] - normalized[[j, i]]).abs() < 1e-10
                        || (normalized[[i, j]].is_nan() && normalized[[j, i]].is_nan()),
                    "normalized_mismatch_rates not symmetric at ({i}, {j})"
                );
            }
        }
    }

    drop(npz);

    let mut archive = ZipArchive::new(File::open(&npz_path).expect("could not reopen npz archive"))
        .expect("failed to read npz as zip");
    {
        let mut samples_file = archive
            .by_name("samples.json")
            .expect("missing samples.json in npz");
        let mut json = String::new();
        samples_file
            .read_to_string(&mut json)
            .expect("failed to read samples.json");
        let samples: Vec<String> =
            serde_json::from_str(&json).expect("invalid JSON in samples.json inside npz");
        assert_eq!(samples, common::expected_sample_ids());
    }

    if degrees {
        let mut labels_file = archive
            .by_name("degree_labels.json")
            .expect("missing degree_labels.json in npz");
        let mut labels_json = String::new();
        labels_file
            .read_to_string(&mut labels_json)
            .expect("failed to read degree_labels.json");
        let labels: Vec<String> =
            serde_json::from_str(&labels_json).expect("invalid JSON in degree_labels.json");
        assert_eq!(
            labels,
            vec![
                "Identical/Twin",
                "First Degree",
                "Second Degree",
                "Third Degree",
                "Unrelated"
            ]
        );
    }

    let plot_path = output_dir.join("mismatch_rates.png");
    let metadata = fs::metadata(plot_path).expect("missing plot output");
    assert!(metadata.len() > 0, "plot output is empty");
}

fn assert_covered_snps(output_dir: &Path, expected: &BTreeMap<String, u64>) {
    let coverage_path = output_dir.join("covered_snps.csv");
    let coverage = read_covered_snps(&coverage_path);
    assert_eq!(
        coverage,
        *expected,
        "unexpected sample covered SNPs in {}",
        coverage_path.display()
    );
}

fn read_records(path: &Path) -> Vec<OutputRecord> {
    let content = fs::read_to_string(path).expect("could not read mismatch rates");
    let mut lines = content.lines();
    let header = lines.next().expect("missing header").trim_end_matches('\r');
    let columns: Vec<&str> = header.split(',').collect();

    let has_ci = columns.contains(&"mismatch_rate_95_ci_lower");
    let has_degrees = columns.contains(&"normalized_mismatch_rate");

    let valid_degrees = [
        "Identical/Twin",
        "First Degree",
        "Second Degree",
        "Third Degree",
        "Unrelated",
    ];

    let expected_cols = 4 + if has_ci { 2 } else { 0 } + if has_degrees { 2 } else { 0 };

    let mut records = Vec::new();
    for line in lines {
        let trimmed = line.trim_end_matches('\r');
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split(',').collect();
        assert_eq!(
            fields.len(),
            expected_cols,
            "unexpected column count in record: {trimmed}"
        );
        let id1 = fields[0].to_string();
        let id2 = fields[1].to_string();
        let overlap: u64 = fields[2].parse().expect("invalid overlap value");
        let rate: f32 = fields[3].parse().expect("invalid mismatch rate value");

        let mut col = 4;

        let (ci_lower, ci_upper) = if has_ci {
            let lo: f32 = fields[col].parse().expect("invalid ci_lower value");
            let hi: f32 = fields[col + 1].parse().expect("invalid ci_upper value");
            col += 2;
            (Some(lo), Some(hi))
        } else {
            (None, None)
        };

        let (normalized_mismatch_rate, degree) = if has_degrees {
            let nmr: f64 = fields[col]
                .parse()
                .expect("invalid normalized_mismatch_rate value");
            let deg = fields[col + 1].to_string();
            assert!(
                valid_degrees.contains(&deg.as_str()),
                "invalid degree: {deg}"
            );
            (Some(nmr), Some(deg))
        } else {
            (None, None)
        };

        records.push(OutputRecord {
            id1,
            id2,
            overlap,
            rate,
            ci_lower,
            ci_upper,
            normalized_mismatch_rate,
            degree,
        });
    }
    records
}

fn read_covered_snps(path: &Path) -> BTreeMap<String, u64> {
    let content = fs::read_to_string(path).expect("could not read sample covered SNPs");
    let mut lines = content.lines();
    let header = lines.next().expect("missing header").trim_end_matches('\r');
    assert_eq!(header, "sample_id,covered_snps");

    let mut coverage = BTreeMap::new();
    for line in lines {
        let trimmed = line.trim_end_matches('\r');
        if trimmed.is_empty() {
            continue;
        }
        let mut fields = trimmed.split(',');
        let id = fields.next().expect("missing sample_id field").to_string();
        let covered: u64 = fields
            .next()
            .expect("missing covered_snps field")
            .parse()
            .expect("invalid covered_snps value");
        assert!(
            fields.next().is_none(),
            "unexpected extra columns in record: {trimmed}"
        );
        coverage.insert(id, covered);
    }
    coverage
}
