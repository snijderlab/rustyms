//! Defines the different levels of complexity a peptide can be.
//! Used for compile time checking for incorrect use of peptides.
use serde::{Deserialize, Serialize};

/// A [`LinearPeptide`] that (potentially) is linked, either with cross-links or branches
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct Linked;

/// A [`LinearPeptide`] that is not cross-linked or branched, but can use the whole breath of the complexity otherwise
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct Linear;

/// A [`LinearPeptide`] that does not have any of the following:
/// * Labile modifications
/// * Global isotope modifications
/// * Charge carriers, use of charged ions apart from protons
/// * Cyclic structures: inter/intra cross-links or branches
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct SimpleLinear;

/// A [`LinearPeptide`] that does not have any of the following:
/// * Ambiguous modifications
/// * Ambiguous amino acid sequence `(?AA)`
///
/// On top of the outlawed features in [`SimpleLinear`].
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct SemiAmbiguous;

/// A [`LinearPeptide`] that does not have any of the following:
/// * Ambiguous amino acids (B/Z)
///
/// On top of the outlawed features in [`SemiAmbiguous`].
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct UnAmbiguous;

impl From<Linear> for Linked {
    fn from(_val: Linear) -> Self {
        Self
    }
}
impl From<SimpleLinear> for Linked {
    fn from(_val: SimpleLinear) -> Self {
        Self
    }
}
impl From<SemiAmbiguous> for Linked {
    fn from(_val: SemiAmbiguous) -> Self {
        Self
    }
}
impl From<UnAmbiguous> for Linked {
    fn from(_val: UnAmbiguous) -> Self {
        Self
    }
}
impl From<SimpleLinear> for Linear {
    fn from(_val: SimpleLinear) -> Self {
        Self
    }
}
impl From<SemiAmbiguous> for Linear {
    fn from(_val: SemiAmbiguous) -> Self {
        Self
    }
}
impl From<UnAmbiguous> for Linear {
    fn from(_val: UnAmbiguous) -> Self {
        Self
    }
}
impl From<SemiAmbiguous> for SimpleLinear {
    fn from(_val: SemiAmbiguous) -> Self {
        Self
    }
}
impl From<UnAmbiguous> for SimpleLinear {
    fn from(_val: UnAmbiguous) -> Self {
        Self
    }
}
impl From<UnAmbiguous> for SemiAmbiguous {
    fn from(_val: UnAmbiguous) -> Self {
        Self
    }
}

/// Indicate that a peptide has at least this level of complexity, or higher
pub trait AtLeast<T> {
    type HighestLevel;
}
impl<T> AtLeast<T> for T {
    type HighestLevel = T;
}
impl AtLeast<Linear> for Linked {
    type HighestLevel = Linked;
}
impl AtLeast<SimpleLinear> for Linked {
    type HighestLevel = Linked;
}
impl AtLeast<SemiAmbiguous> for Linked {
    type HighestLevel = Linked;
}
impl AtLeast<UnAmbiguous> for Linked {
    type HighestLevel = Linked;
}
impl AtLeast<SimpleLinear> for Linear {
    type HighestLevel = Linear;
}
impl AtLeast<SemiAmbiguous> for Linear {
    type HighestLevel = Linear;
}
impl AtLeast<UnAmbiguous> for Linear {
    type HighestLevel = Linear;
}
impl AtLeast<SemiAmbiguous> for SimpleLinear {
    type HighestLevel = SimpleLinear;
}
impl AtLeast<UnAmbiguous> for SimpleLinear {
    type HighestLevel = SimpleLinear;
}
impl AtLeast<UnAmbiguous> for SemiAmbiguous {
    type HighestLevel = SemiAmbiguous;
}

/// Type level max between two Peptide complexities
pub trait HighestOf<T> {
    type HighestLevel;
}

impl<T> HighestOf<T> for T {
    type HighestLevel = T;
}

impl HighestOf<Linked> for Linear {
    type HighestLevel = Linked;
}
impl HighestOf<Linked> for SimpleLinear {
    type HighestLevel = Linked;
}
impl HighestOf<Linked> for SemiAmbiguous {
    type HighestLevel = Linked;
}
impl HighestOf<Linked> for UnAmbiguous {
    type HighestLevel = Linked;
}
impl HighestOf<Linear> for SimpleLinear {
    type HighestLevel = Linear;
}
impl HighestOf<Linear> for SemiAmbiguous {
    type HighestLevel = Linear;
}
impl HighestOf<Linear> for UnAmbiguous {
    type HighestLevel = Linear;
}
impl HighestOf<SimpleLinear> for SemiAmbiguous {
    type HighestLevel = SimpleLinear;
}
impl HighestOf<SimpleLinear> for UnAmbiguous {
    type HighestLevel = SimpleLinear;
}
impl HighestOf<SemiAmbiguous> for UnAmbiguous {
    type HighestLevel = SemiAmbiguous;
}

impl HighestOf<Linear> for Linked {
    type HighestLevel = Linked;
}
impl HighestOf<SimpleLinear> for Linked {
    type HighestLevel = Linked;
}
impl HighestOf<SemiAmbiguous> for Linked {
    type HighestLevel = Linked;
}
impl HighestOf<UnAmbiguous> for Linked {
    type HighestLevel = Linked;
}
impl HighestOf<SimpleLinear> for Linear {
    type HighestLevel = Linear;
}
impl HighestOf<SemiAmbiguous> for Linear {
    type HighestLevel = Linear;
}
impl HighestOf<UnAmbiguous> for Linear {
    type HighestLevel = Linear;
}
impl HighestOf<SemiAmbiguous> for SimpleLinear {
    type HighestLevel = SimpleLinear;
}
impl HighestOf<UnAmbiguous> for SimpleLinear {
    type HighestLevel = SimpleLinear;
}
impl HighestOf<UnAmbiguous> for SemiAmbiguous {
    type HighestLevel = SemiAmbiguous;
}
