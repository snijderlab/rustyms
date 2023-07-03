use crate::AminoAcid;

#[derive(Debug, PartialEq, Eq)]
pub enum PlacementRule {
    AminoAcid(AminoAcid, Position),
    Terminal(Position),
}

impl PlacementRule {
    pub fn is_possible(&self, aa: AminoAcid, index: usize, length: usize) -> bool {
        match self {
            Self::AminoAcid(r_aa, r_pos) => *r_aa == aa && r_pos.is_possible(index, length),
            Self::Terminal(r_pos) => {
                r_pos.is_possible(index, length) && (index == length - 1 || index == 0)
            }
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum Position {
    Anywhere,
    AnyNTerm,
    AnyCTerm,
    ProteinNTerm,
    ProteinCTerm,
}

impl Position {
    const fn is_possible(&self, index: usize, length: usize) -> bool {
        match self {
            Self::Anywhere => true,
            Self::AnyNTerm | Self::ProteinNTerm => index == 0,
            Self::AnyCTerm | Self::ProteinCTerm => index == length - 1,
        }
    }
}
