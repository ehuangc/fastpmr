mod common;

use ndarray::Array2;
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
    let output = run_fastpmr(&dataset, None, false, None);
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs);
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
    let output = run_fastpmr(&dataset, None, false, None);
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs);
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
    let output = run_fastpmr(&dataset, None, true, None);
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_npz_outputs(&dataset.output_dir, &expected_pairs);
}

#[test]
fn variant_indices_limit_sites() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "filtered").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let expected_pairs = common::expected_pair_stats_filtered_variants();
    let output = run_fastpmr(&dataset, Some("1-30000"), false, None);
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let records = assert_outputs(&dataset.output_dir, &expected_pairs);
    assert_eq!(
        records.len(),
        expected_pairs.len(),
        "unexpected number of pairwise records after filtering"
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

    let output = run_fastpmr(&dataset, None, false, Some(&csv_path));
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let expected_subset: BTreeMap<_, _> = common::expected_pair_stats_all_variants()
        .into_iter()
        .filter(|((id1, id2), _)| {
            matches!(
                (id1.as_str(), id2.as_str()),
                ("Sample1", "Sample2") | ("Sample1", "Sample3") | ("Sample2", "Sample4")
            )
        })
        .collect();

    let records = assert_outputs(&dataset.output_dir, &expected_subset);
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
fn sample_pairs_csv_with_unknown_sample_fails() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "sample-pairs-err").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }
    let csv_path = dataset.prefix.with_extension("pairs.csv");
    fs::write(&csv_path, "Sample1,Unknown\n").unwrap();

    let output = run_fastpmr(&dataset, None, false, Some(&csv_path));
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

fn run_fastpmr(
    dataset: &common::Dataset,
    variant_spec: Option<&str>,
    npz: bool,
    sample_pairs_csv: Option<&Path>,
) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_fastpmr"));
    command
        .arg("--prefix")
        .arg(dataset.prefix.as_os_str())
        .arg("--output-directory")
        .arg(dataset.output_dir.as_os_str());
    if let Some(spec) = variant_spec {
        command.arg("--variant-indices").arg(spec);
    }
    if npz {
        command.arg("--npz");
    }
    if let Some(path) = sample_pairs_csv {
        command.arg("--sample-pairs-csv").arg(path);
    }
    command.output().expect("failed to run fastpmr")
}

#[derive(Clone, Debug, PartialEq)]
struct OutputRecord {
    id1: String,
    id2: String,
    overlap: u64,
    rate: f32,
}

fn assert_outputs(
    output_dir: &Path,
    expected_pairs: &BTreeMap<(String, String), common::PairStats>,
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

    let plot_path = output_dir.join("mismatch_rates.png");
    let metadata = fs::metadata(plot_path).expect("missing plot output");
    assert!(metadata.len() > 0, "plot output is empty");
    records
}

fn assert_npz_outputs(
    output_dir: &Path,
    expected_pairs: &BTreeMap<(String, String), common::PairStats>,
) {
    let csv_path = output_dir.join("mismatch_rates.csv");
    assert!(
        !csv_path.exists(),
        "unexpected mismatch_rates.csv output alongside npz"
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
    let mismatches: Array2<u64> = npz.by_name("mismatches").expect("missing mismatches array");
    let totals: Array2<u64> = npz.by_name("totals").expect("missing totals array");
    let overlaps: Array2<u64> = npz
        .by_name("site_overlaps")
        .expect("missing site_overlaps array");
    let expected_shape = &[common::N_SAMPLES, common::N_SAMPLES];
    assert_eq!(mismatches.shape(), expected_shape);
    assert_eq!(totals.shape(), expected_shape);
    assert_eq!(overlaps.shape(), expected_shape);

    for idx in 0..common::N_SAMPLES {
        assert_eq!(mismatches[[idx, idx]], 0);
        assert_eq!(totals[[idx, idx]], 0);
        assert_eq!(overlaps[[idx, idx]], 0);
    }

    let samples = common::expected_sample_ids();
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
            assert_eq!(
                mismatches[[i, j]],
                expectation.mismatches,
                "unexpected mismatches for {} / {}",
                samples[i],
                samples[j]
            );
            assert_eq!(
                mismatches[[j, i]],
                expectation.mismatches,
                "unexpected symmetric mismatches for {} / {}",
                samples[j],
                samples[i]
            );
            assert_eq!(
                totals[[i, j]],
                expectation.totals,
                "unexpected totals for {} / {}",
                samples[i],
                samples[j]
            );
            assert_eq!(
                totals[[j, i]],
                expectation.totals,
                "unexpected symmetric totals for {} / {}",
                samples[j],
                samples[i]
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

    drop(npz);

    let mut archive = ZipArchive::new(File::open(&npz_path).expect("could not reopen npz archive"))
        .expect("failed to read npz as zip");
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

    let plot_path = output_dir.join("mismatch_rates.png");
    let metadata = fs::metadata(plot_path).expect("missing plot output");
    assert!(metadata.len() > 0, "plot output is empty");
}

fn read_records(path: &Path) -> Vec<OutputRecord> {
    let content = fs::read_to_string(path).expect("could not read mismatch rates");
    let mut lines = content.lines();
    let header = lines.next().expect("missing header").trim_end_matches('\r');
    assert_eq!(header, "id1,id2,n_site_overlaps,mismatch_rate");

    let mut records = Vec::new();
    for line in lines {
        let trimmed = line.trim_end_matches('\r');
        if trimmed.is_empty() {
            continue;
        }
        let mut fields = trimmed.split(',');
        let id1 = fields.next().expect("missing id1 field").to_string();
        let id2 = fields.next().expect("missing id2 field").to_string();
        let overlap: u64 = fields
            .next()
            .expect("missing overlap field")
            .parse()
            .expect("invalid overlap value");
        let rate: f32 = fields
            .next()
            .expect("missing mismatch rate field")
            .parse()
            .expect("invalid mismatch rate value");
        assert!(
            fields.next().is_none(),
            "unexpected extra columns in record: {trimmed}"
        );
        records.push(OutputRecord {
            id1,
            id2,
            overlap,
            rate,
        });
    }
    records
}
