//! Defines the different levels of complexity a peptide can be.
//! Used for compile time checking for incorrect use of peptides.
use serde::{Deserialize, Serialize};

use crate::LinearPeptide;

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
/// * or when the sequence is empty.
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct Simple;

/// A [`LinearPeptide`] that does not have any of the following:
/// * Ambiguous modifications
/// * Ambiguous amino acid sequence `(?AA)`
/// On top of the outlawed features in [`Simple`].
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct VerySimple;

/// A [`LinearPeptide`] that does not have any of the following:
/// * Ambiguous amino acids (B/Z)
/// On top of the outlawed features in [`VerySimple`].
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct ExtremelySimple;

impl Into<Linked> for Linear {
    fn into(self) -> Linked {
        Linked
    }
}
impl Into<Linked> for Simple {
    fn into(self) -> Linked {
        Linked
    }
}
impl Into<Linked> for VerySimple {
    fn into(self) -> Linked {
        Linked
    }
}
impl Into<Linked> for ExtremelySimple {
    fn into(self) -> Linked {
        Linked
    }
}
impl Into<Linear> for Simple {
    fn into(self) -> Linear {
        Linear
    }
}
impl Into<Linear> for VerySimple {
    fn into(self) -> Linear {
        Linear
    }
}
impl Into<Linear> for ExtremelySimple {
    fn into(self) -> Linear {
        Linear
    }
}
impl Into<Simple> for VerySimple {
    fn into(self) -> Simple {
        Simple
    }
}
impl Into<Simple> for ExtremelySimple {
    fn into(self) -> Simple {
        Simple
    }
}
impl Into<VerySimple> for ExtremelySimple {
    fn into(self) -> VerySimple {
        VerySimple
    }
}

impl LinearPeptide<Linked> {
    pub fn linear(self) -> Option<LinearPeptide<Linear>> {
        if self
            .sequence
            .iter()
            .all(|seq| seq.modifications.iter().all(|m| m.simple().is_some()))
        {
            Some(self.mark())
        } else {
            None
        }
    }

    pub fn simple(self) -> Option<LinearPeptide<Simple>> {
        self.linear().and_then(LinearPeptide::<Linear>::simple)
    }

    pub fn very_simple(self) -> Option<LinearPeptide<VerySimple>> {
        self.linear().and_then(LinearPeptide::<Linear>::very_simple)
    }

    pub fn extremely_simple(self) -> Option<LinearPeptide<ExtremelySimple>> {
        self.linear()
            .and_then(LinearPeptide::<Linear>::extremely_simple)
    }
}

impl LinearPeptide<Linear> {
    pub fn simple(self) -> Option<LinearPeptide<Simple>> {
        if self.labile.is_empty()
            && self.get_global().is_empty()
            && self.charge_carriers.is_none()
            && !self.sequence.is_empty()
            && !self
                .sequence
                .iter()
                .any(|s| s.modifications.iter().any(|m| m.simple().is_none()))
        {
            Some(self.mark())
        } else {
            None
        }
    }

    pub fn very_simple(self) -> Option<LinearPeptide<VerySimple>> {
        self.simple().and_then(LinearPeptide::<Simple>::very_simple)
    }

    pub fn extremely_simple(self) -> Option<LinearPeptide<ExtremelySimple>> {
        self.simple()
            .and_then(LinearPeptide::<Simple>::extremely_simple)
    }
}

impl LinearPeptide<Simple> {
    pub fn very_simple(self) -> Option<LinearPeptide<VerySimple>> {
        if self.ambiguous_modifications.is_empty()
            && !self.sequence.iter().any(|seq| seq.ambiguous.is_some())
        {
            Some(self.mark())
        } else {
            None
        }
    }

    pub fn extremely_simple(self) -> Option<LinearPeptide<ExtremelySimple>> {
        self.very_simple()
            .and_then(LinearPeptide::<VerySimple>::extremely_simple)
    }
}

impl LinearPeptide<VerySimple> {
    pub fn extremely_simple(self) -> Option<LinearPeptide<ExtremelySimple>> {
        if !self
            .sequence
            .iter()
            .any(|seq| seq.aminoacid == crate::AminoAcid::B || seq.aminoacid == crate::AminoAcid::Z)
        {
            Some(self.mark())
        } else {
            None
        }
    }
}
