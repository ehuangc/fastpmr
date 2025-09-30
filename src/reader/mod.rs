pub mod packedancestrymap;

use crate::error::Result;
use crate::model::Site;

pub trait SiteReader: Iterator<Item = Result<Site>> + Send {
    fn samples(&self) -> &[String];
    fn n_sites(&self) -> usize;
}
