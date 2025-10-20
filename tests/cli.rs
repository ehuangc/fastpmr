mod common;

use std::fs;
use std::path::Path;
use std::process::Command;

#[test]
fn packedancestrymap_cli_generates_outputs() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "packed").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let output = run_fastpmr(&dataset, None);
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_outputs(
        &dataset.output_dir,
        common::expected_overlap_all(),
        common::expected_rate_all(),
    );
}

#[test]
fn transposed_packedancestrymap_cli_generates_outputs() {
    let dataset = common::create_dataset(common::GenoFormat::Transposed, "transposed").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let output = run_fastpmr(&dataset, None);
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_outputs(
        &dataset.output_dir,
        common::expected_overlap_all(),
        common::expected_rate_all(),
    );
}

#[test]
fn variant_indices_limit_sites() {
    let dataset = common::create_dataset(common::GenoFormat::Packed, "filtered").unwrap();
    if dataset.output_dir.exists() {
        fs::remove_dir_all(&dataset.output_dir).unwrap();
    }

    let output = run_fastpmr(&dataset, Some("1-30000"));
    assert!(
        output.status.success(),
        "fastpmr failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_outputs(
        &dataset.output_dir,
        common::expected_overlap_filtered(),
        common::expected_rate_filtered(),
    );
}

fn run_fastpmr(dataset: &common::Dataset, variant_spec: Option<&str>) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_fastpmr"));
    command
        .arg("--prefix")
        .arg(dataset.prefix.as_os_str())
        .arg("--output-directory")
        .arg(dataset.output_dir.as_os_str());
    if let Some(spec) = variant_spec {
        command.arg("--variant-indices").arg(spec);
    }
    command.output().expect("failed to run fastpmr")
}

fn assert_outputs(output_dir: &Path, expected_overlap: u64, expected_rate: f32) {
    let csv_path = output_dir.join("mismatch_rates.txt");
    let record = read_single_record(&csv_path);
    assert_eq!(record[0], "Sample1");
    assert_eq!(record[1], "Sample2");

    let overlap: u64 = record[2].parse().expect("invalid overlap");
    assert_eq!(overlap, expected_overlap);

    let rate: f32 = record[3].parse().expect("invalid mismatch rate");
    assert!(
        (rate - expected_rate).abs() < 1e-6,
        "unexpected mismatch rate: got {rate}, expected {expected_rate}"
    );

    let plot_path = output_dir.join("mismatch_rates.png");
    let metadata = fs::metadata(plot_path).expect("missing plot output");
    assert!(metadata.len() > 0, "plot output is empty");
}

fn read_single_record(path: &Path) -> [String; 4] {
    let content = fs::read_to_string(path).expect("could not read mismatch rates");
    let mut lines = content.lines();

    let header = lines.next().expect("missing header").trim_end_matches('\r');
    assert_eq!(header, "id1,id2,n_site_overlaps,mismatch_rate");

    let record_line = lines
        .next()
        .expect("missing record line")
        .trim_end_matches('\r');
    assert!(lines.next().is_none(), "unexpected extra lines in csv");

    let mut fields = record_line.split(',').map(|s| s.to_string());
    let id1 = fields.next().expect("missing id1");
    let id2 = fields.next().expect("missing id2");
    let overlap = fields.next().expect("missing site overlap");
    let rate = fields.next().expect("missing rate");
    assert!(fields.next().is_none(), "unexpected extra columns");

    [id1, id2, overlap, rate]
}
