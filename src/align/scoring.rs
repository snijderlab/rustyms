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

/// The (slightly modified) blosum62 matrix for aminoacid homology scoring
pub const BLOSUM62: &[&[i8]] = include!("blosum62.txt");
