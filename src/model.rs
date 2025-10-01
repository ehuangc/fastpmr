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
        const MISMATCH: [[u64; 4]; 4] = [
            [0, 1, 2, 0], // Ref vs. *
            [1, 1, 1, 0], // Het vs. *
            [2, 1, 0, 0], // Alt vs. *
            [0, 0, 0, 0], // Missing vs. *
        ];
        MISMATCH[self as usize][other as usize]
    }
}

pub struct Site {
    pub genotypes: Vec<Allele>,
}
