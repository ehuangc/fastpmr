use crate::counts::Counts;
use crate::error::{CustomError, Result};
use plotters::coord::combinators::IntoLogRange;
use plotters::prelude::*;
use plotters::style::{FontStyle, register_font};

pub fn write_mismatch_rates(counts: &Counts, path: &str) -> Result<()> {
    let n_samples = counts.n_samples();
    let overlaps = counts.overlaps();
    let (pairs, rates) = counts.mismatch_rates();

    let mut wtr = csv::Writer::from_path(path)?;
    wtr.write_record(&["id1", "id2", "n_overlap", "mismatch_rate"])?;

    for i in 0..n_samples {
        for j in (i + 1)..n_samples {
            let counter_idx = counts.idx(i, j);
            let overlap = overlaps[counter_idx];
            let rate = rates[counter_idx];
            wtr.serialize((
                pairs[counter_idx].0.as_str(),
                pairs[counter_idx].1.as_str(),
                overlap,
                rate,
            ))?;
        }
    }
    wtr.flush().map_err(|e| CustomError::Write {
        source: e,
        path: path.into(),
    })?;
    Ok(())
}

pub fn plot_mismatch_rates(counts: &Counts, path: &str) -> Result<()> {
    let (_pairs, rates) = counts.mismatch_rates();
    let overlaps = counts.overlaps();
    let filtered_percentages: Vec<f32> = rates
        .iter()
        .zip(overlaps.iter())
        .filter_map(|(&rate, &overlap)| {
            if overlap >= 30000 && rate.is_finite() {
                Some(rate * 100.0)
            } else {
                None
            }
        })
        .collect();

    const BIN_SIZE: f32 = 0.5;
    let max_percentage = filtered_percentages.iter().copied().fold(0.0f32, f32::max);
    let n_bins: usize;
    if max_percentage < 50.0 {
        n_bins = 100;
    } else {
        let n_bins_unrounded = (max_percentage / BIN_SIZE).ceil() as usize;
        n_bins = n_bins_unrounded.div_ceil(10) * 10; // Round up to next multiple of 10
    }

    let median: f32;
    let mut sorted = filtered_percentages.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        median = (sorted[mid - 1] + sorted[mid]) / 2.0;
    } else {
        median = sorted[mid];
    }

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
    register_font("ibm-plex-mono", FontStyle::Normal, IBM_PLEX_MONO).map_err(|_| CustomError::Font)?;

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
            BIN_SIZE..(n_bins as f32) * BIN_SIZE,
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
        .x_label_formatter(&|x| format!("{:.0}", x))
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
                    vec![
                        PathElement::new(vec![(x0, y), (x1, y)], s.clone()), // top
                        PathElement::new(vec![(x0, 0usize), (x0, y)], s.clone()), // left
                        PathElement::new(vec![(x1, 0usize), (x1, y)], s),    // right
                    ]
                    .into_iter()
                }),
        )
        .map_err(|e| CustomError::Plot {
            source: Box::new(e),
        })?;

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

    root_area.present().map_err(|e| CustomError::Plot {
        source: Box::new(e),
    })?;
    Ok(())
}
