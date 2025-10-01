#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Allele {
    Ref = 0,
    Het = 1,
    Alt = 2,
    Missing,
}

impl Allele {
    pub fn mismatch(self, other: Self) -> u64 {
        (self as u8).abs_diff(other as u8) as u64
    }
}

pub struct Site {
    pub genotypes: Vec<Allele>,
}
