use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

pub const N_SAMPLES: usize = 4;
pub const CORE_VARIANTS: usize = 30_000;
const HIGH_MISMATCH_VARIANTS: usize = 14_999;
const EXTRA_IDENTICAL_VARIANTS: usize = 1;
const EXTRA_MISSING_VARIANTS: usize = 1;
pub const TOTAL_VARIANTS: usize = CORE_VARIANTS + EXTRA_IDENTICAL_VARIANTS + EXTRA_MISSING_VARIANTS;

const ALT: u8 = 0b00;
const HET: u8 = 0b01;
const REF: u8 = 0b10;
const MISSING: u8 = 0b11;
const PLINK_HEADER: [u8; 3] = [0x6c, 0x1b, 0x01];

static NEXT_ID: AtomicUsize = AtomicUsize::new(0);

#[derive(Clone, Copy)]
pub enum GenoFormat {
    Packed,
    Transposed,
    Plink,
    Eigenstrat,
}

pub struct Dataset {
    pub prefix: PathBuf,
    pub output_dir: PathBuf,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PairStats {
    pub mismatches: u64,
    pub totals: u64,
}

impl PairStats {
    pub fn overlap(&self) -> u64 {
        self.totals / 2
    }

    pub fn mismatch_rate(&self) -> f32 {
        if self.totals == 0 {
            f32::NAN
        } else {
            self.mismatches as f32 / self.totals as f32
        }
    }
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

    let variants = build_variants();
    let sample_ids = sample_ids();
    let variant_ids = variant_ids();
    write_ind(prefix.with_extension("ind"), &sample_ids)?;
    write_snp(prefix.with_extension("snp"))?;
    match format {
        GenoFormat::Packed => write_geno(
            prefix.with_extension("geno"),
            &variants,
            &sample_ids,
            &variant_ids,
        )?,
        GenoFormat::Transposed => write_tgeno(
            prefix.with_extension("geno"),
            &variants,
            &sample_ids,
            &variant_ids,
        )?,
        GenoFormat::Plink => {
            write_bed(prefix.with_extension("bed"), &variants)?;
            write_bim(prefix.with_extension("bim"), &variant_ids)?;
            write_fam(prefix.with_extension("fam"), &sample_ids)?;
        }
        GenoFormat::Eigenstrat => write_unpacked_geno(prefix.with_extension("geno"), &variants)?,
    }

    Ok(Dataset { prefix, output_dir })
}

pub fn expected_pair_stats_all_variants() -> BTreeMap<(String, String), PairStats> {
    build_pair_expectations(true)
}

pub fn expected_pair_stats_filtered_variants() -> BTreeMap<(String, String), PairStats> {
    build_pair_expectations(false)
}

pub fn expected_sample_ids() -> Vec<String> {
    sample_ids()
}

fn sample_ids() -> Vec<String> {
    let mut samples = Vec::new();
    for i in 0..N_SAMPLES {
        samples.push(format!("Sample{}", i + 1));
    }
    samples
}

fn variant_ids() -> Vec<String> {
    let mut variants = Vec::new();
    for i in 0..TOTAL_VARIANTS {
        variants.push(format!("rs{}", i + 1));
    }
    variants
}

fn build_variants() -> Vec<[u8; N_SAMPLES]> {
    fn high_mismatch_variant() -> [u8; N_SAMPLES] {
        let mut variant = [ALT; N_SAMPLES];
        variant[0] = REF;
        variant
    }

    fn heterozygous_variant() -> [u8; N_SAMPLES] {
        [HET; N_SAMPLES]
    }

    fn identical_variant() -> [u8; N_SAMPLES] {
        [ALT; N_SAMPLES]
    }

    fn missing_variant() -> [u8; N_SAMPLES] {
        let mut variant = [ALT; N_SAMPLES];
        variant[0] = MISSING;
        variant
    }

    let mut variants = Vec::with_capacity(TOTAL_VARIANTS);
    for _ in 0..HIGH_MISMATCH_VARIANTS {
        variants.push(high_mismatch_variant());
    }
    for _ in HIGH_MISMATCH_VARIANTS..CORE_VARIANTS {
        variants.push(heterozygous_variant());
    }
    variants.push(identical_variant());
    variants.push(missing_variant());
    variants
}

fn write_ind(path: impl AsRef<Path>, samples: &[String]) -> io::Result<()> {
    let mut file = File::create(path)?;
    for sample in samples {
        writeln!(file, "{} M 0", sample)?;
    }
    Ok(())
}

fn write_snp(path: impl AsRef<Path>) -> io::Result<()> {
    let mut file = File::create(path)?;
    for idx in 0..TOTAL_VARIANTS {
        writeln!(file, "rs{} 1 0.0 {} A G", idx + 1, idx + 1)?;
    }
    Ok(())
}

fn write_geno(
    path: impl AsRef<Path>,
    variants: &[[u8; N_SAMPLES]],
    sample_ids: &[String],
    variant_ids: &[String],
) -> io::Result<()> {
    let mut file = File::create(path)?;
    let block_size = 48usize.max(N_SAMPLES.div_ceil(4));

    let (sample_hash, variant_hash) = header_hash(sample_ids, variant_ids);
    let header_str = format!(
        "GENO {} {} {} {}",
        N_SAMPLES,
        variants.len(),
        sample_hash,
        variant_hash
    );
    let mut header_block = vec![0u8; block_size];
    header_block[..header_str.len()].copy_from_slice(header_str.as_bytes());
    header_block[header_str.len()] = 0;
    file.write_all(&header_block)?;

    let bytes_needed = N_SAMPLES.div_ceil(4);
    for variant in variants {
        let mut block = vec![0u8; block_size];
        for (byte_idx, byte) in block.iter_mut().enumerate().take(bytes_needed) {
            let mut value = 0u8;
            for within in 0..4 {
                let sample_idx = byte_idx * 4 + within;
                if sample_idx >= N_SAMPLES {
                    break;
                }
                let code = variant[sample_idx] & 0b11;
                let shift = 6 - 2 * within;
                value |= code << shift;
            }
            *byte = value;
        }
        file.write_all(&block)?;
    }
    Ok(())
}

fn write_unpacked_geno(path: impl AsRef<Path>, variants: &[[u8; N_SAMPLES]]) -> io::Result<()> {
    let mut file = File::create(path)?;
    for variant in variants {
        let mut line = String::with_capacity(N_SAMPLES);
        for &code in variant {
            let allele = match code {
                ALT => '0',
                HET => '1',
                REF => '2',
                MISSING => '9',
                _ => unreachable!(),
            };
            line.push(allele);
        }
        writeln!(file, "{line}")?;
    }
    Ok(())
}

fn write_tgeno(
    path: impl AsRef<Path>,
    variants: &[[u8; N_SAMPLES]],
    sample_ids: &[String],
    variant_ids: &[String],
) -> io::Result<()> {
    let mut file = File::create(path)?;

    let (sample_hash, variant_hash) = header_hash(sample_ids, variant_ids);
    let header_str = format!(
        "TGENO {} {} {} {}",
        N_SAMPLES,
        variants.len(),
        sample_hash,
        variant_hash
    );
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

fn write_bed(path: impl AsRef<Path>, variants: &[[u8; N_SAMPLES]]) -> io::Result<()> {
    let mut file = File::create(path)?;
    file.write_all(&PLINK_HEADER)?;

    let bytes_per_variant = N_SAMPLES.div_ceil(4);
    let mut block = vec![0u8; bytes_per_variant];
    for variant in variants {
        block.fill(0);
        for (sample_idx, &sample) in variant.iter().enumerate() {
            let code = match sample {
                REF => 0b00,
                HET => 0b10,
                ALT => 0b11,
                MISSING => 0b01,
                _ => unreachable!(),
            };
            let byte_idx = sample_idx / 4;
            let shift = (sample_idx % 4) * 2;
            block[byte_idx] |= code << shift;
        }
        file.write_all(&block)?;
    }
    Ok(())
}

fn write_bim(path: impl AsRef<Path>, variant_ids: &[String]) -> io::Result<()> {
    let mut file = File::create(path)?;
    for (idx, variant) in variant_ids.iter().enumerate() {
        writeln!(file, "1 {variant} 0 {} A G", idx + 1)?;
    }
    Ok(())
}

fn write_fam(path: impl AsRef<Path>, samples: &[String]) -> io::Result<()> {
    let mut file = File::create(path)?;
    for sample in samples {
        writeln!(file, "0 {sample} 0 0 0 -9")?;
    }
    Ok(())
}

fn header_hash(sample_ids: &[String], variant_ids: &[String]) -> (String, String) {
    fn hashone(id: &str) -> u32 {
        let mut hash: u32 = 0;
        for &b in id.as_bytes() {
            if b == b'\0' {
                break;
            }
            hash = hash.wrapping_mul(23).wrapping_add(b as u32);
        }
        hash
    }

    fn hasharr(ids: &[String]) -> u32 {
        let mut hash: u32 = 0;
        for id in ids {
            hash = hash.wrapping_mul(17) ^ hashone(id);
        }
        hash
    }

    let sample_hash = hasharr(sample_ids);
    let variant_hash = hasharr(variant_ids);
    (
        format!("{:08x}", sample_hash),
        format!("{:08x}", variant_hash),
    )
}

fn build_pair_expectations(include_extra_variants: bool) -> BTreeMap<(String, String), PairStats> {
    let samples = sample_ids();
    let mut expectations = BTreeMap::new();

    let mismatches_with_sample1 =
        (HIGH_MISMATCH_VARIANTS as u64) * 2 + (CORE_VARIANTS - HIGH_MISMATCH_VARIANTS) as u64;
    let mismatches_without_sample1 = (CORE_VARIANTS - HIGH_MISMATCH_VARIANTS) as u64;

    let totals_with_sample1 = if include_extra_variants {
        ((CORE_VARIANTS + EXTRA_IDENTICAL_VARIANTS) * 2) as u64
    } else {
        (CORE_VARIANTS * 2) as u64
    };

    let totals_without_sample1 = if include_extra_variants {
        ((CORE_VARIANTS + EXTRA_IDENTICAL_VARIANTS + EXTRA_MISSING_VARIANTS) * 2) as u64
    } else {
        (CORE_VARIANTS * 2) as u64
    };

    for i in 0..samples.len() {
        for j in (i + 1)..samples.len() {
            let key = (samples[i].clone(), samples[j].clone());
            let (mismatches, totals) = if i == 0 || j == 0 {
                (mismatches_with_sample1, totals_with_sample1)
            } else {
                (mismatches_without_sample1, totals_without_sample1)
            };
            expectations.insert(key, PairStats { mismatches, totals });
        }
    }

    expectations
}
