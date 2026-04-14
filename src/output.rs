use crate::counts::{ConfidenceIntervals, Counts};
use crate::degrees::{Degree, DegreeResults};
use crate::error::{CustomError, Result};
use ndarray::Array1;
use ndarray_npy::WriteNpyExt;
use plotters::coord::combinators::IntoLogRange;
use plotters::prelude::*;
use plotters::style::{FontStyle, register_font};
use std::io::Write;
use std::path::Path;
use zip::ZipWriter;
use zip::write::SimpleFileOptions;

pub fn write_covered_snps(counts: &Counts, path: &impl AsRef<Path>) -> Result<()> {
    let samples = counts.samples();
    let covered_snps = counts.covered_snps();

    let mut wtr = csv::Writer::from_path(path)?;
    wtr.write_record(["sample_id", "covered_snps"])?;
    for (sample, covered) in samples.iter().zip(covered_snps.iter()) {
        wtr.serialize((sample.as_str(), *covered))?;
    }
    wtr.flush().map_err(|e| CustomError::Write {
        source: e,
        path: path.as_ref().into(),
    })?;
    Ok(())
}

pub fn write_mismatch_rates(
    counts: &Counts,
    degree_results: Option<&DegreeResults>,
    ci_results: Option<&ConfidenceIntervals>,
    path: &impl AsRef<Path>,
) -> Result<()> {
    let n_samples = counts.n_samples();
    let overlaps = counts.site_overlaps();
    let pairs = counts.pairs();
    let rates = counts.mismatch_rates();

    let mut wtr = csv::Writer::from_path(path)?;
    {
        let mut header = vec!["id1", "id2", "n_site_overlap", "mismatch_rate"];
        if ci_results.is_some() {
            header.push("mismatch_rate_95_ci_lower");
            header.push("mismatch_rate_95_ci_upper");
        }
        if degree_results.is_some() {
            header.push("normalized_mismatch_rate");
            header.push("degree");
        }
        wtr.write_record(&header)?;
    }

    for i in 0..n_samples {
        for j in (i + 1)..n_samples {
            let pair_idx = counts.pair_idx(i, j);
            if !counts.should_count_pair(i, j) {
                continue;
            }
            let overlap = overlaps[pair_idx];
            let rate = rates[pair_idx];

            let mut record = vec![
                pairs[pair_idx].0.clone(),
                pairs[pair_idx].1.clone(),
                overlap.to_string(),
                rate.to_string(),
            ];
            if let Some(ci) = ci_results {
                record.push(ci.lower[pair_idx].to_string());
                record.push(ci.upper[pair_idx].to_string());
            }
            if let Some(dr) = degree_results {
                record.push(dr.normalized_mismatch_rates[pair_idx].to_string());
                record.push(dr.degrees[pair_idx].to_string());
            }
            wtr.write_record(&record)?;
        }
    }
    wtr.flush().map_err(|e| CustomError::Write {
        source: e,
        path: path.as_ref().into(),
    })?;
    Ok(())
}

pub fn write_counts_npz(
    counts: &Counts,
    degree_results: Option<&DegreeResults>,
    ci_results: Option<&ConfidenceIntervals>,
    path: &impl AsRef<Path>,
) -> Result<()> {
    let options = SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated)
        .large_file(true);
    let mut zip = ZipWriter::new(std::fs::File::create(path).map_err(|e| CustomError::Write {
        source: e,
        path: path.as_ref().into(),
    })?);
    // Write one array at a time and free it immediately
    {
        let covered_snps = Array1::from_vec(counts.covered_snps().to_vec());
        zip.start_file("covered_snps.npy", options)?;
        covered_snps.write_npy(&mut zip)?;
    }
    {
        let site_overlaps = counts.site_overlaps_2d();
        zip.start_file("n_site_overlaps.npy", options)?;
        site_overlaps.write_npy(&mut zip)?;
    }
    {
        let mismatch_rates = counts.mismatch_rates_2d();
        zip.start_file("mismatch_rates.npy", options)?;
        mismatch_rates.write_npy(&mut zip)?;
    }
    if let Some(dr) = degree_results {
        {
            let normalized = dr.normalized_mismatch_rates_2d(counts.n_samples(), counts);
            zip.start_file("normalized_mismatch_rates.npy", options)?;
            normalized.write_npy(&mut zip)?;
        }
        {
            let degrees = dr.degrees_2d(counts.n_samples(), counts);
            zip.start_file("degrees.npy", options)?;
            degrees.write_npy(&mut zip)?;
        }
    }
    if let Some(ci) = ci_results {
        {
            let lower = counts.ci_95_lower_2d(ci);
            zip.start_file("mismatch_rates_95_ci_lower.npy", options)?;
            lower.write_npy(&mut zip)?;
        }
        {
            let upper = counts.ci_95_upper_2d(ci);
            zip.start_file("mismatch_rates_95_ci_upper.npy", options)?;
            upper.write_npy(&mut zip)?;
        }
    }
    // Write JSON metadata
    let samples = counts.samples();
    zip.start_file("samples.json", options)?;
    zip.write_all(serde_json::to_string(&samples)?.as_bytes())
        .map_err(|e| CustomError::Write {
            source: e,
            path: path.as_ref().into(),
        })?;
    if degree_results.is_some() {
        zip.start_file("degree_labels.json", options)?;
        zip.write_all(serde_json::to_string(&Degree::LABELS)?.as_bytes())
            .map_err(|e| CustomError::Write {
                source: e,
                path: path.as_ref().into(),
            })?;
    }
    zip.finish()?;
    Ok(())
}

pub fn plot_mismatch_rates(counts: &Counts, path: &impl AsRef<Path>) -> Result<()> {
    let rates = counts.mismatch_rates();
    let overlaps = counts.site_overlaps();
    let n_samples = counts.n_samples();
    let mut filtered_percentages = Vec::new();
    for i in 0..n_samples {
        for j in (i + 1)..n_samples {
            if !counts.should_count_pair(i, j) {
                continue;
            }
            let pair_idx = counts.pair_idx(i, j);
            let overlap = overlaps[pair_idx];
            let rate = rates[pair_idx];
            if overlap >= 30000 && rate.is_finite() {
                filtered_percentages.push(rate * 100.0);
            }
        }
    }

    const BIN_SIZE: f32 = 0.5;
    let max_percentage = filtered_percentages.iter().copied().fold(0.0f32, f32::max);
    let n_bins: usize = if max_percentage < 50.0 {
        100
    } else {
        let n_bins_unrounded = (max_percentage / BIN_SIZE).ceil() as usize;
        n_bins_unrounded.div_ceil(10) * 10 // Round up to next multiple of 10
    };

    let mut sorted = filtered_percentages.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = sorted.len() / 2;

    let median: Option<f32> = if sorted.is_empty() {
        None
    } else if sorted.len() % 2 == 0 {
        Some((sorted[mid - 1] + sorted[mid]) / 2.0)
    } else {
        Some(sorted[mid])
    };

    let mut bin_counts = vec![0usize; n_bins];
    for &rate in &filtered_percentages {
        let bin_idx = (rate / BIN_SIZE).floor() as usize;
        if bin_idx < n_bins {
            bin_counts[bin_idx] += 1;
        } else {
            unreachable!("rate {rate} out of histogram x-axis range");
        }
    }

    const IBM_PLEX_MONO: &[u8] =
        include_bytes!("../assets/fonts/ibm-plex-mono/IBMPlexMono-Regular.ttf");
    register_font("ibm-plex-mono", FontStyle::Normal, IBM_PLEX_MONO)
        .map_err(|_| CustomError::Font)?;

    let root_area = BitMapBackend::new(path, (3840, 2160)).into_drawing_area();
    root_area.fill(&WHITE).map_err(|e| CustomError::Plot {
        source: Box::new(e),
    })?;

    let mut chart = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 230)
        .set_label_area_size(LabelAreaPosition::Bottom, 160)
        .margin(20)
        .margin_right(60)
        .caption("Pairwise Mismatch Rate Distribution", ("ibm-plex-mono", 96))
        .build_cartesian_2d(
            0.0..(n_bins as f32) * BIN_SIZE,
            (0usize..filtered_percentages.len()).log_scale(),
        )
        .map_err(|e| CustomError::Plot {
            source: Box::new(e),
        })?;

    fn superscript(n: u32) -> String {
        let mut out = String::new();
        for ch in n.to_string().chars() {
            out.push(match ch {
                '0' => '⁰',
                '1' => '¹',
                '2' => '²',
                '3' => '³',
                '4' => '⁴',
                '5' => '⁵',
                '6' => '⁶',
                '7' => '⁷',
                '8' => '⁸',
                '9' => '⁹',
                _ => unreachable!(),
            });
        }
        out
    }

    fn exponent_label(y: usize) -> String {
        if y == 1 {
            return "1".to_string();
        }
        if y == 10 {
            return "10".to_string();
        }
        let log_y = (y as f64).log10();
        if (log_y - log_y.round()).abs() < 1e-8 {
            return format!("10{}", superscript(log_y.round() as u32));
        }
        String::new()
    }

    chart
        .configure_mesh()
        .label_style(("ibm-plex-mono", 80))
        .y_desc("Sample pairs")
        .y_label_formatter(&|&y| exponent_label(y))
        .x_desc("Pairwise mismatch rate (%)")
        .x_label_formatter(&|x| {
            if *x == 0.0 {
                String::new()
            } else {
                format!("{:.0}", x)
            }
        })
        .draw()
        .map_err(|e| CustomError::Plot {
            source: Box::new(e),
        })?;

    // Draw histogram bars
    chart
        .draw_series((0..n_bins).map(|i| {
            let x0 = i as f32 * BIN_SIZE;
            let x1 = x0 + BIN_SIZE;
            Rectangle::new([(x0, 0usize), (x1, bin_counts[i])], BLUE.mix(0.4).filled())
        }))
        .map_err(|e| CustomError::Plot {
            source: Box::new(e),
        })?;

    // Draw histogram bar outlines: top + sides only (no bottom)
    chart
        .draw_series(
            (0..n_bins)
                // Draw nothing if bar height is zero
                .filter(|&i| bin_counts[i] != 0)
                .flat_map(|i| {
                    let x0 = i as f32 * BIN_SIZE;
                    let x1 = x0 + BIN_SIZE;
                    let y = bin_counts[i];

                    let s = BLACK.stroke_width(1);
                    // Don't draw left side for first bar
                    if i == 0 {
                        return vec![
                            PathElement::new(vec![(x0, y), (x1, y)], s),      // top
                            PathElement::new(vec![(x1, 0usize), (x1, y)], s), // right
                        ]
                        .into_iter();
                    }
                    vec![
                        PathElement::new(vec![(x0, y), (x1, y)], s),      // top
                        PathElement::new(vec![(x0, 0usize), (x0, y)], s), // left
                        PathElement::new(vec![(x1, 0usize), (x1, y)], s), // right
                    ]
                    .into_iter()
                }),
        )
        .map_err(|e| CustomError::Plot {
            source: Box::new(e),
        })?;

    if let Some(median) = median {
        // Draw vertical median line
        chart
            .draw_series(std::iter::once(PathElement::new(
                vec![(median, 0usize), (median, filtered_percentages.len())],
                ShapeStyle {
                    color: RED.mix(0.8).to_rgba(),
                    filled: true,
                    stroke_width: 6,
                },
            )))
            .map_err(|e| CustomError::Plot {
                source: Box::new(e),
            })?;

        // Write median label
        chart
            .draw_series(std::iter::once(Text::new(
                format!("Median: {:.2}%", median),
                (median + 1.0, filtered_percentages.len() / 3),
                ("ibm-plex-mono", 80).into_font().color(&RED.mix(0.8)),
            )))
            .map_err(|e| CustomError::Plot {
                source: Box::new(e),
            })?;
    }

    root_area.present().map_err(|e| CustomError::Plot {
        source: Box::new(e),
    })?;
    Ok(())
}
