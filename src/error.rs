use thiserror::Error;

#[derive(Debug, Error)]
pub enum CustomError {
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

    #[error("header block does not start with \"GENO\" or \"TGENO\"")]
    PackedAncestryMapHeaderPrefix,

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

    #[error("need at least 2 samples (got {n_samples})")]
    SampleCount { n_samples: usize },

    #[error("need at least 1 variant (got {n_variants})")]
    VariantCount { n_variants: usize },
}

pub type Result<T> = std::result::Result<T, CustomError>;
