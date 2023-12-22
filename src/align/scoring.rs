use serde::{Deserialize, Serialize};

/// The type of a single match step
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum MatchType {
    /// Aminoacid + Mass identity
    FullIdentity,
    /// Aminoacid + Mass mismatch
    IdentityMassMismatch,
    /// Full mismatch
    #[default]
    Mismatch,
    /// Set of aminoacids + mods with the same mass but different sequence
    Isobaric,
    /// Set of aminoacids + mods in a different order in the two sequences
    Rotation,
    /// A gap
    Gap,
}

pub const MISMATCH: i8 = -1;
pub const MASS_MISMATCH_PENALTY: i8 = -1;
pub const BASE_SPECIAL: i8 = 1;
pub const ROTATED: i8 = 3;
pub const ISOBARIC: i8 = 2;
pub const GAP_START_PENALTY: i8 = -5;
pub const GAP_EXTEND_PENALTY: i8 = -1;

/// Matrices from: <https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/util/tables/> and <https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/>
/// The UO columns are added by me (see top left for the original matrix used by me) (B/J/Z is the rounded down average of the corresponding non ambiguous AAs) (All these are exactly the same for all matrices)
pub mod matrices {
    /// BLOSUM45 matrix
    pub const BLOSUM45: &[[i8; crate::AminoAcid::TOTAL_NUMBER]; crate::AminoAcid::TOTAL_NUMBER] =
        include!("matrices/blosum45.txt");
    /// BLOSUM50 matrix
    pub const BLOSUM50: &[[i8; crate::AminoAcid::TOTAL_NUMBER]; crate::AminoAcid::TOTAL_NUMBER] =
        include!("matrices/blosum50.txt");
    /// BLOSUM62 matrix
    pub const BLOSUM62: &[[i8; crate::AminoAcid::TOTAL_NUMBER]; crate::AminoAcid::TOTAL_NUMBER] =
        include!("matrices/blosum62.txt");
    /// BLOSUM80 matrix
    pub const BLOSUM80: &[[i8; crate::AminoAcid::TOTAL_NUMBER]; crate::AminoAcid::TOTAL_NUMBER] =
        include!("matrices/blosum80.txt");
    /// BLOSUM90 matrix
    pub const BLOSUM90: &[[i8; crate::AminoAcid::TOTAL_NUMBER]; crate::AminoAcid::TOTAL_NUMBER] =
        include!("matrices/blosum90.txt");
    /// Identity matrix (9 for equal, -5 for not equal)
    pub const IDENTITY: &[[i8; crate::AminoAcid::TOTAL_NUMBER]; crate::AminoAcid::TOTAL_NUMBER] =
        include!("matrices/identity.txt");
    /// PAM30 matrix
    pub const PAM30: &[[i8; crate::AminoAcid::TOTAL_NUMBER]; crate::AminoAcid::TOTAL_NUMBER] =
        include!("matrices/pam30.txt");
    /// PAM70 matrix
    pub const PAM70: &[[i8; crate::AminoAcid::TOTAL_NUMBER]; crate::AminoAcid::TOTAL_NUMBER] =
        include!("matrices/pam70.txt");
    /// PAM250 matrix
    pub const PAM250: &[[i8; crate::AminoAcid::TOTAL_NUMBER]; crate::AminoAcid::TOTAL_NUMBER] =
        include!("matrices/pam250.txt");
}
