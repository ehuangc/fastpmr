use crate::counts::Counts;
use crate::error::{CustomError, Result};
use plotters::coord::combinators::IntoLogRange;
use plotters::prelude::*;
use plotters::style::{register_font, FontStyle};

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
            wtr.write_record(&[
                pairs[counter_idx].0.clone(),
                pairs[counter_idx].1.clone(),
                overlap.to_string(),
                rate.to_string(),
            ])?;
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

    const N_BINS: usize = 100;
    const BIN_SIZE: f32 = 0.5;

    let median: f32;
    let mut sorted = filtered_percentages.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        median = (sorted[mid - 1] + sorted[mid]) / 2.0;
    } else {
        median = sorted[mid];
    }

    let mut bin_counts = vec![0usize; N_BINS];
    for &rate in &filtered_percentages {
        let bin_idx = (rate / BIN_SIZE).floor() as usize;
        if bin_idx < N_BINS {
            bin_counts[bin_idx] += 1;
        } else {
            unreachable!("rate {rate} out of histogram x-axis range");
        }
    }

    const ROBOTO_MONO: &[u8] = include_bytes!("../assets/fonts/roboto-mono/RobotoMono-Regular.ttf");
    register_font(
        "roboto-mono",
        FontStyle::Normal,
        ROBOTO_MONO,
    ).map_err(|_| CustomError::Font)?;

    let root_area = BitMapBackend::new(path, (3840, 2160)).into_drawing_area();
    root_area.fill(&WHITE).map_err(|e| CustomError::Plot {
        source: Box::new(e),
    })?;

    let mut chart = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 230)
        .set_label_area_size(LabelAreaPosition::Bottom, 160)
        .margin(20)
        .margin_right(60)
        .caption("Pairwise Mismatch Rate Distribution", ("roboto-mono", 96))
        .build_cartesian_2d(
            BIN_SIZE..(N_BINS as f32) * BIN_SIZE,
            (0usize..filtered_percentages.len()).log_scale(),
        )
        .map_err(|e| CustomError::Plot {
            source: Box::new(e),
        })?;

    chart
        .configure_mesh()
        .label_style(("roboto-mono", 80))
        .y_desc("Sample pairs")
        .x_desc("Pairwise mismatch rate (%)")
        .x_label_formatter(&|x| format!("{:.0}", x))
        .draw()
        .map_err(|e| CustomError::Plot {
            source: Box::new(e),
        })?;

    chart
        .draw_series((0..N_BINS).map(|i| {
            let x0 = i as f32 * BIN_SIZE;
            let x1 = x0 + BIN_SIZE;
            Rectangle::new([(x0, 0usize), (x1, bin_counts[i])], BLUE.mix(0.5).filled())
        }))
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

    // Median label
    chart
        .draw_series(std::iter::once(Text::new(
            format!("Median: {:.2}%", median),
            (median + 1.0, filtered_percentages.len() / 3),
            ("roboto-mono", 80).into_font().color(&RED.mix(0.8)),
        )))
        .map_err(|e| CustomError::Plot {
            source: Box::new(e),
        })?;

    root_area.present().map_err(|e| CustomError::Plot {
        source: Box::new(e),
    })?;
    Ok(())
}
