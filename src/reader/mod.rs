pub mod packedancestrymap;

use crate::error::Result;
use crate::model::Site;

pub trait SiteReader {
    fn samples(&self) -> &[String];
    fn n_sites(&self) -> usize;
    fn next_site(&mut self) -> Result<Option<Site>>;
}
