use std::fs::{self, File};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

const N_SAMPLES: usize = 2;
pub const CORE_VARIANTS: usize = 30_000;
const HIGH_MISMATCH_VARIANTS: usize = 14_999;
const EXTRA_IDENTICAL_VARIANTS: usize = 1;
const EXTRA_MISSING_VARIANTS: usize = 1;
pub const TOTAL_VARIANTS: usize = CORE_VARIANTS + EXTRA_IDENTICAL_VARIANTS + EXTRA_MISSING_VARIANTS;

const ALT: u8 = 0b00;
const HET: u8 = 0b01;
const REF: u8 = 0b10;
const MISSING: u8 = 0b11;

static NEXT_ID: AtomicUsize = AtomicUsize::new(0);

#[derive(Clone, Copy)]
pub enum GenoFormat {
    Packed,
    Transposed,
}

pub struct Dataset {
    pub prefix: PathBuf,
    pub output_dir: PathBuf,
}

pub fn create_dataset(format: GenoFormat, label: &str) -> io::Result<Dataset> {
    let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
    let base_dir = std::env::temp_dir().join("fastpmr-tests").join(format!(
        "{}-{}-{}",
        std::process::id(),
        id,
        label
    ));
    fs::create_dir_all(&base_dir)?;

    let prefix = base_dir.join("dataset");
    let output_dir = base_dir.join("output");

    write_ind(prefix.with_extension("ind"))?;
    write_snp(prefix.with_extension("snp"))?;
    let variants = build_variants();
    match format {
        GenoFormat::Packed => write_geno(prefix.with_extension("geno"), &variants)?,
        GenoFormat::Transposed => write_tgeno(prefix.with_extension("geno"), &variants)?,
    }

    Ok(Dataset { prefix, output_dir })
}

pub fn expected_overlap_all() -> u64 {
    (CORE_VARIANTS + EXTRA_IDENTICAL_VARIANTS) as u64
}

pub fn expected_rate_all() -> f32 {
    let core_mismatches =
        (HIGH_MISMATCH_VARIANTS * 2 + (CORE_VARIANTS - HIGH_MISMATCH_VARIANTS)) as f32;
    let extra_mismatches = 0f32;
    let core_totals = (CORE_VARIANTS * 2) as f32;
    let extra_totals = (EXTRA_IDENTICAL_VARIANTS * 2) as f32;
    (core_mismatches + extra_mismatches) / (core_totals + extra_totals)
}

pub fn expected_overlap_filtered() -> u64 {
    CORE_VARIANTS as u64
}

pub fn expected_rate_filtered() -> f32 {
    let mismatches = (HIGH_MISMATCH_VARIANTS * 2 + (CORE_VARIANTS - HIGH_MISMATCH_VARIANTS)) as f32;
    let totals = (CORE_VARIANTS * 2) as f32;
    mismatches / totals
}

fn build_variants() -> Vec<[u8; N_SAMPLES]> {
    let mut variants = Vec::with_capacity(TOTAL_VARIANTS);
    for _ in 0..HIGH_MISMATCH_VARIANTS {
        variants.push([REF, ALT]);
    }
    for _ in HIGH_MISMATCH_VARIANTS..CORE_VARIANTS {
        variants.push([HET, HET]);
    }
    variants.push([ALT, ALT]);
    variants.push([MISSING, ALT]);
    variants
}

fn write_ind(path: impl AsRef<Path>) -> io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "Sample1 M 0")?;
    writeln!(file, "Sample2 F 0")?;
    Ok(())
}

fn write_snp(path: impl AsRef<Path>) -> io::Result<()> {
    let mut file = File::create(path)?;
    for idx in 0..TOTAL_VARIANTS {
        writeln!(file, "rs{} 1 0.0 {} A G", idx + 1, idx + 1)?;
    }
    Ok(())
}

fn write_geno(path: impl AsRef<Path>, variants: &[[u8; N_SAMPLES]]) -> io::Result<()> {
    let mut file = File::create(path)?;
    let block_size = 48usize.max(N_SAMPLES.div_ceil(4));
    let header_str = format!("GENO {} {} 0 0", N_SAMPLES, variants.len());
    let mut header_block = vec![0u8; block_size];
    header_block[..header_str.len()].copy_from_slice(header_str.as_bytes());
    header_block[header_str.len()] = 0;
    file.write_all(&header_block)?;

    let bytes_needed = N_SAMPLES.div_ceil(4);
    for variant in variants {
        let mut block = vec![0u8; block_size];
        for byte_idx in 0..bytes_needed {
            let mut value = 0u8;
            for within in 0..4 {
                let sample_idx = byte_idx * 4 + within;
                if sample_idx >= N_SAMPLES {
                    continue;
                }
                let code = variant[sample_idx] & 0b11;
                let shift = 6 - 2 * within;
                value |= code << shift;
            }
            block[byte_idx] = value;
        }
        file.write_all(&block)?;
    }
    Ok(())
}

fn write_tgeno(path: impl AsRef<Path>, variants: &[[u8; N_SAMPLES]]) -> io::Result<()> {
    let mut file = File::create(path)?;
    let header_str = format!("TGENO {} {} 0 0", N_SAMPLES, variants.len());
    let mut header_block = vec![0u8; 48];
    header_block[..header_str.len()].copy_from_slice(header_str.as_bytes());
    header_block[header_str.len()] = 0;
    file.write_all(&header_block)?;

    let sample_block_size = 48usize.max(variants.len().div_ceil(4));
    for sample_idx in 0..N_SAMPLES {
        let mut block = vec![0u8; sample_block_size];
        for (variant_idx, variant) in variants.iter().enumerate() {
            let code = variant[sample_idx] & 0b11;
            let byte_idx = variant_idx / 4;
            let offset = variant_idx % 4;
            let shift = 6 - 2 * offset;
            block[byte_idx] |= code << shift;
        }
        file.write_all(&block)?;
    }
    Ok(())
}
