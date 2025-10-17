use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::{CustomError, Result};

pub(super) const IND_FIELDS: usize = 3;
pub(super) const SNP_FIELDS: usize = 6;

pub(super) fn read_eigenstrat_ind(path: &impl AsRef<Path>) -> Result<Vec<String>> {
    let f = File::open(path).map_err(|e| CustomError::ReadWithPath {
        source: e,
        path: path.as_ref().to_path_buf(),
    })?;
    let f = BufReader::new(f);
    let mut sample_ids: Vec<String> = Vec::new();

    for (line_idx, line) in f.lines().enumerate() {
        let line = line.map_err(|e| CustomError::ReadWithPath {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;
        let line = line.trim();
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != IND_FIELDS {
            return Err(CustomError::EigenstratIndFields {
                line_num: line_idx + 1,
                n_fields: fields.len(),
                expected: IND_FIELDS,
            });
        }
        let sample_id = fields[0].to_string();
        sample_ids.push(sample_id);
    }
    Ok(sample_ids)
}

pub(super) fn read_eigenstrat_snp(path: &impl AsRef<Path>) -> Result<Vec<String>> {
    let f = File::open(path).map_err(|e| CustomError::ReadWithPath {
        source: e,
        path: path.as_ref().to_path_buf(),
    })?;
    let f = BufReader::new(f);
    let mut variant_ids: Vec<String> = Vec::new();

    for (line_idx, line) in f.lines().enumerate() {
        let line = line.map_err(|e| CustomError::ReadWithPath {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;
        let line = line.trim();
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != SNP_FIELDS {
            return Err(CustomError::EigenstratSnpFields {
                line_num: line_idx + 1,
                n_fields: fields.len(),
                expected: SNP_FIELDS,
            });
        }
        let variant_id = fields[0].to_string();
        variant_ids.push(variant_id);
    }
    Ok(variant_ids)
}
