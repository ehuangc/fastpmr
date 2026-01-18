pub mod common;
pub mod packedancestrymap;
pub mod plink;
pub mod transposed_packedancestrymap;
pub mod unpacked_eigenstrat;

use crate::error::Result;
use crate::model::Site;

pub trait SiteReader: Iterator<Item = Result<Site>> + Send {
    fn samples(&self) -> &[String];
    fn n_sites(&self) -> usize;
}
