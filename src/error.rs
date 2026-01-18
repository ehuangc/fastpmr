use thiserror::Error;

#[derive(Debug, Error)]
pub enum CustomError {
    #[error("could not create output directory")]
    OutputDir {
        #[source]
        source: std::io::Error,
    },

    #[error(
        "could not find supported input files for prefix {prefix} (.geno/.ind/.snp or .bed/.bim/.fam)"
    )]
    InputFilesMissing { prefix: String },

    #[error(
        "sample pairs CSV must contain either one column of sample IDs or two columns of explicit pairs"
    )]
    SamplePairsColumns,

    #[error("sample pairs CSV did not contain any pairs")]
    SamplePairsEmpty,

    #[error("sample pair references unknown sample: {sample}")]
    SamplePairUnknownSample { sample: String },

    #[error("sample pair cannot include the same sample twice: {sample}")]
    SamplePairDuplicate { sample: String },

    #[error("could not parse variant index: {arg}")]
    VariantIndexInt {
        #[source]
        source: std::num::ParseIntError,
        arg: String,
    },

    #[error("variant index out-of-bounds: {idx} > {n_variants}")]
    VariantIndexHigh { idx: usize, n_variants: usize },

    #[error("variant indices are 1-based and must be positive")]
    VariantIndexLow,

    #[error("could not build thread pool")]
    ThreadPoolBuild(#[from] rayon::ThreadPoolBuildError),

    #[error("could not read {path}")]
    ReadWithPath {
        #[source]
        source: std::io::Error,
        path: std::path::PathBuf,
    },

    #[error("could not read input file(s)")]
    ReadWithoutPath {
        #[source]
        source: std::io::Error,
    },

    #[error("could not write to {path}")]
    Write {
        #[source]
        source: std::io::Error,
        path: std::path::PathBuf,
    },

    #[error("could not write to CSV")]
    CsvWrite(#[from] csv::Error),

    #[error("could not read CSV {path}")]
    CsvRead {
        #[source]
        source: csv::Error,
        path: std::path::PathBuf,
    },

    #[error("could not write to NPZ")]
    NpzWrite(#[from] ndarray_npy::WriteNpzError),

    #[error("could not write to ZIP")]
    ZipWrite(#[from] zip::result::ZipError),

    #[error("could not write to sample JSON")]
    JsonWrite(#[from] serde_json::Error),

    #[error("could not plot mismatch distribution")]
    Plot {
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    #[error("could not find font")]
    Font,

    #[error("file too small to contain data")]
    PackedAncestryMapFileSize,

    #[error("no null byte found in header block")]
    PackedAncestryMapHeaderNullByte,

    #[error("header block is not valid UTF-8")]
    PackedAncestryMapHeaderUtf8 {
        #[source]
        source: std::str::Utf8Error,
    },

    #[error("header block does not start with \"GENO\"")]
    PackedAncestryMapHeaderGeno,

    #[error("header block does not start with \"TGENO\"")]
    PackedAncestryMapHeaderTgeno,

    #[error("expected {expected} fields (got {n_fields}) in header")]
    PackedAncestryMapHeaderFields { n_fields: usize, expected: usize },

    #[error("could not parse N in header")]
    PackedAncestryMapHeaderN {
        #[source]
        source: std::num::ParseIntError,
    },

    #[error("could not parse V in header")]
    PackedAncestryMapHeaderV {
        #[source]
        source: std::num::ParseIntError,
    },

    #[error("header and .ind file sample count disagree (header N={n_header}, .ind count={n_ind})")]
    PackedAncestryMapNAgreement { n_header: usize, n_ind: usize },

    #[error(
        "header and .snp file variant count disagree (header V={n_header}, .snp count={n_snp})"
    )]
    PackedAncestryMapVAgreement { n_header: usize, n_snp: usize },

    #[error("header sample hash does not match (expected: {expected}, found: {found})")]
    PackedAncestryMapSampleHash { expected: String, found: String },

    #[error("header variant hash does not match (expected: {expected}, found: {found})")]
    PackedAncestryMapVariantHash { expected: String, found: String },

    #[error("invalid PLINK .bed header magic bytes")]
    PlinkBedHeaderMagic,

    #[error("PLINK .bed file must be SNP-major")]
    PlinkBedMode,

    #[error("PLINK .bed file size mismatch (expected {expected} bytes, found {found} bytes)")]
    PlinkBedFileSize { expected: u64, found: u64 },

    #[error("expected {expected} fields (got {n_fields}) in line {line_num} of .fam file")]
    PlinkFamFields {
        line_num: usize,
        n_fields: usize,
        expected: usize,
    },

    #[error("expected {expected} fields (got {n_fields}) in line {line_num} of .bim file")]
    PlinkBimFields {
        line_num: usize,
        n_fields: usize,
        expected: usize,
    },

    #[error("expected {expected} fields (got {n_fields}) in line {line_num} of .ind file")]
    EigenstratIndFields {
        line_num: usize,
        n_fields: usize,
        expected: usize,
    },

    #[error("expected {expected} fields (got {n_fields}) in line {line_num} of .snp file")]
    EigenstratSnpFields {
        line_num: usize,
        n_fields: usize,
        expected: usize,
    },

    #[error("expected {expected} genotypes (got {n_fields}) in line {line_num} of .geno file")]
    EigenstratGenoFields {
        line_num: usize,
        n_fields: usize,
        expected: usize,
    },

    #[error(".geno file contains {found} variants but .snp lists {expected}")]
    EigenstratGenoVariantCount { expected: usize, found: usize },

    #[error("need at least 2 samples (got {n_samples})")]
    SampleCount { n_samples: usize },

    #[error("need at least 1 variant (got {n_variants})")]
    VariantCount { n_variants: usize },
}

pub type Result<T> = std::result::Result<T, CustomError>;
