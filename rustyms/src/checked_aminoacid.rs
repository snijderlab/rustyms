use std::marker::PhantomData;

use serde::{Deserialize, Serialize};

use crate::{
    AminoAcid, Chemical, MolecularFormula, Multi, MultiChemical, SemiAmbiguous, UnAmbiguous,
};

#[derive(Ord, PartialOrd, Debug, Serialize, Deserialize)]
pub struct CheckedAminoAcid<T> {
    aminoacid: AminoAcid,
    marker: PhantomData<T>,
}

#[allow(non_upper_case_globals, missing_docs)]
impl CheckedAminoAcid<UnAmbiguous> {
    pub const A: Self = Self::Alanine;
    pub const C: Self = Self::Cysteine;
    pub const D: Self = Self::AsparticAcid;
    pub const E: Self = Self::GlutamicAcid;
    pub const F: Self = Self::Phenylalanine;
    pub const G: Self = Self::Glycine;
    pub const H: Self = Self::Histidine;
    pub const I: Self = Self::Isoleucine;
    pub const J: Self = Self::AmbiguousLeucine;
    pub const K: Self = Self::Lysine;
    pub const L: Self = Self::Leucine;
    pub const M: Self = Self::Methionine;
    pub const N: Self = Self::Asparagine;
    pub const O: Self = Self::Pyrrolysine;
    pub const P: Self = Self::Proline;
    pub const Q: Self = Self::Glutamine;
    pub const R: Self = Self::Arginine;
    pub const S: Self = Self::Serine;
    pub const T: Self = Self::Threonine;
    pub const U: Self = Self::Selenocysteine;
    pub const V: Self = Self::Valine;
    pub const W: Self = Self::Tryptophan;
    pub const X: Self = Self::Unknown;
    pub const Y: Self = Self::Tyrosine;
    pub const Ala: Self = Self::Alanine;
    pub const Cys: Self = Self::Cysteine;
    pub const Asn: Self = Self::Asparagine;
    pub const Asp: Self = Self::AsparticAcid;
    pub const Glu: Self = Self::GlutamicAcid;
    pub const Phe: Self = Self::Phenylalanine;
    pub const Gly: Self = Self::Glycine;
    pub const His: Self = Self::Histidine;
    pub const Ile: Self = Self::Isoleucine;
    pub const Xle: Self = Self::AmbiguousLeucine;
    pub const Lys: Self = Self::Lysine;
    pub const Leu: Self = Self::Leucine;
    pub const Met: Self = Self::Methionine;
    pub const Pyl: Self = Self::Pyrrolysine;
    pub const Pro: Self = Self::Proline;
    pub const Gln: Self = Self::Glutamine;
    pub const Arg: Self = Self::Arginine;
    pub const Ser: Self = Self::Serine;
    pub const Thr: Self = Self::Threonine;
    pub const Sec: Self = Self::Selenocysteine;
    pub const Val: Self = Self::Valine;
    pub const Trp: Self = Self::Tryptophan;
    pub const Tyr: Self = Self::Tyrosine;
    pub const Xaa: Self = Self::Unknown;
    pub const Alanine: Self = Self {
        aminoacid: AminoAcid::Alanine,
        marker: PhantomData,
    };
    pub const Cysteine: Self = Self {
        aminoacid: AminoAcid::Cysteine,
        marker: PhantomData,
    };
    pub const AsparticAcid: Self = Self {
        aminoacid: AminoAcid::AsparticAcid,
        marker: PhantomData,
    };
    pub const GlutamicAcid: Self = Self {
        aminoacid: AminoAcid::GlutamicAcid,
        marker: PhantomData,
    };
    pub const Phenylalanine: Self = Self {
        aminoacid: AminoAcid::Phenylalanine,
        marker: PhantomData,
    };
    pub const Glycine: Self = Self {
        aminoacid: AminoAcid::Glycine,
        marker: PhantomData,
    };
    pub const Histidine: Self = Self {
        aminoacid: AminoAcid::Histidine,
        marker: PhantomData,
    };
    pub const Isoleucine: Self = Self {
        aminoacid: AminoAcid::Isoleucine,
        marker: PhantomData,
    };
    pub const AmbiguousLeucine: Self = Self {
        aminoacid: AminoAcid::AmbiguousLeucine,
        marker: PhantomData,
    };
    pub const Lysine: Self = Self {
        aminoacid: AminoAcid::Lysine,
        marker: PhantomData,
    };
    pub const Leucine: Self = Self {
        aminoacid: AminoAcid::Leucine,
        marker: PhantomData,
    };
    pub const Methionine: Self = Self {
        aminoacid: AminoAcid::Methionine,
        marker: PhantomData,
    };
    pub const Asparagine: Self = Self {
        aminoacid: AminoAcid::Asparagine,
        marker: PhantomData,
    };
    pub const Pyrrolysine: Self = Self {
        aminoacid: AminoAcid::Pyrrolysine,
        marker: PhantomData,
    };
    pub const Proline: Self = Self {
        aminoacid: AminoAcid::Proline,
        marker: PhantomData,
    };
    pub const Glutamine: Self = Self {
        aminoacid: AminoAcid::Glutamine,
        marker: PhantomData,
    };
    pub const Arginine: Self = Self {
        aminoacid: AminoAcid::Arginine,
        marker: PhantomData,
    };
    pub const Serine: Self = Self {
        aminoacid: AminoAcid::Serine,
        marker: PhantomData,
    };
    pub const Threonine: Self = Self {
        aminoacid: AminoAcid::Threonine,
        marker: PhantomData,
    };
    pub const Selenocysteine: Self = Self {
        aminoacid: AminoAcid::Selenocysteine,
        marker: PhantomData,
    };
    pub const Valine: Self = Self {
        aminoacid: AminoAcid::Valine,
        marker: PhantomData,
    };
    pub const Tryptophan: Self = Self {
        aminoacid: AminoAcid::Tryptophan,
        marker: PhantomData,
    };
    pub const Unknown: Self = Self {
        aminoacid: AminoAcid::Unknown,
        marker: PhantomData,
    };
    pub const Tyrosine: Self = Self {
        aminoacid: AminoAcid::Tyrosine,
        marker: PhantomData,
    };

    /// All amino acids with a unique mass (no I/L in favour of J, no B, no Z, and no X)
    pub const UNIQUE_MASS_AMINO_ACIDS: &'static [Self] = &[
        Self::Glycine,
        Self::Alanine,
        Self::Arginine,
        Self::Asparagine,
        Self::AsparticAcid,
        Self::Cysteine,
        Self::Glutamine,
        Self::GlutamicAcid,
        Self::Histidine,
        Self::AmbiguousLeucine,
        Self::Lysine,
        Self::Methionine,
        Self::Phenylalanine,
        Self::Proline,
        Self::Serine,
        Self::Threonine,
        Self::Tryptophan,
        Self::Tyrosine,
        Self::Valine,
        Self::Selenocysteine,
        Self::Pyrrolysine,
    ];

    /// All 20 canonical amino acids
    pub const CANONICAL_AMINO_ACIDS: &'static [Self] = &[
        Self::Glycine,
        Self::Alanine,
        Self::Arginine,
        Self::Asparagine,
        Self::AsparticAcid,
        Self::Cysteine,
        Self::Glutamine,
        Self::GlutamicAcid,
        Self::Histidine,
        Self::Leucine,
        Self::Isoleucine,
        Self::Lysine,
        Self::Methionine,
        Self::Phenylalanine,
        Self::Proline,
        Self::Serine,
        Self::Threonine,
        Self::Tryptophan,
        Self::Tyrosine,
        Self::Valine,
    ];
}

#[allow(non_upper_case_globals, missing_docs)]
impl CheckedAminoAcid<SemiAmbiguous> {
    pub const B: Self = Self::AmbiguousAsparagine;
    pub const Z: Self = Self::AmbiguousGlutamine;
    pub const Asx: Self = Self::AmbiguousAsparagine;
    pub const Glx: Self = Self::AmbiguousGlutamine;
    pub const AmbiguousAsparagine: Self = Self {
        aminoacid: AminoAcid::AmbiguousAsparagine,
        marker: PhantomData,
    };
    pub const AmbiguousGlutamine: Self = Self {
        aminoacid: AminoAcid::AmbiguousGlutamine,
        marker: PhantomData,
    };
}

impl Into<CheckedAminoAcid<SemiAmbiguous>> for AminoAcid {
    fn into(self) -> CheckedAminoAcid<SemiAmbiguous> {
        CheckedAminoAcid::new(self)
    }
}

impl Into<CheckedAminoAcid<SemiAmbiguous>> for &AminoAcid {
    fn into(self) -> CheckedAminoAcid<SemiAmbiguous> {
        CheckedAminoAcid::new(*self)
    }
}

impl CheckedAminoAcid<SemiAmbiguous> {
    pub fn new(aminoacid: AminoAcid) -> CheckedAminoAcid<SemiAmbiguous> {
        CheckedAminoAcid {
            aminoacid,
            marker: PhantomData,
        }
    }
}

impl<T> CheckedAminoAcid<T> {
    pub(super) fn mark<M>(self) -> CheckedAminoAcid<M> {
        CheckedAminoAcid {
            aminoacid: self.aminoacid,
            marker: PhantomData,
        }
    }

    /// Check if this amino acid is an unambiguous amino acid
    pub fn is_unambiguous(self) -> Option<CheckedAminoAcid<UnAmbiguous>> {
        if self.aminoacid != AminoAcid::AmbiguousAsparagine
            && self.aminoacid != AminoAcid::AmbiguousGlutamine
        {
            Some(self.mark())
        } else {
            None
        }
    }

    /// Check if two amino acids are considered identical. X is identical to anything, J to IL, B to ND, Z to EQ.
    pub fn canonical_identical(self, rhs: Self) -> bool {
        self.aminoacid.canonical_identical(rhs.aminoacid)
    }

    pub const fn char(self) -> char {
        self.aminoacid.char()
    }

    pub const fn aminoacid(self) -> AminoAcid {
        self.aminoacid
    }
}

impl Chemical for CheckedAminoAcid<UnAmbiguous> {
    /// Get all possible formula for an unambiguous amino acid (X is defined to be an empty formula)
    fn formula(
        &self,
        _sequence_index: crate::SequencePosition,
        _peptide_index: usize,
    ) -> MolecularFormula {
        match self.aminoacid {
            AminoAcid::Alanine => molecular_formula!(H 5 C 3 O 1 N 1),
            AminoAcid::Arginine => molecular_formula!(H 12 C 6 O 1 N 4),
            AminoAcid::Asparagine => molecular_formula!(H 6 C 4 O 2 N 2),
            AminoAcid::AsparticAcid => molecular_formula!(H 5 C 4 O 3 N 1),
            AminoAcid::Cysteine => molecular_formula!(H 5 C 3 O 1 N 1 S 1),
            AminoAcid::Glutamine => molecular_formula!(H 8 C 5 O 2 N 2),
            AminoAcid::GlutamicAcid => molecular_formula!(H 7 C 5 O 3 N 1),
            AminoAcid::Glycine => molecular_formula!(H 3 C 2 O 1 N 1),
            AminoAcid::Histidine => molecular_formula!(H 7 C 6 O 1 N 3),
            AminoAcid::AmbiguousLeucine | AminoAcid::Isoleucine | AminoAcid::Leucine => {
                molecular_formula!(H 11 C 6 O 1 N 1)
            }
            AminoAcid::Lysine => molecular_formula!(H 12 C 6 O 1 N 2),
            AminoAcid::Methionine => molecular_formula!(H 9 C 5 O 1 N 1 S 1),
            AminoAcid::Phenylalanine => molecular_formula!(H 9 C 9 O 1 N 1),
            AminoAcid::Proline => molecular_formula!(H 7 C 5 O 1 N 1),
            AminoAcid::Pyrrolysine => molecular_formula!(H 19 C 11 O 2 N 3),
            AminoAcid::Selenocysteine => molecular_formula!(H 5 C 3 O 1 N 1 Se 1),
            AminoAcid::Serine => molecular_formula!(H 5 C 3 O 2 N 1),
            AminoAcid::Threonine => molecular_formula!(H 7 C 4 O 2 N 1),
            AminoAcid::Tryptophan => molecular_formula!(H 10 C 11 O 1 N 2),
            AminoAcid::Tyrosine => molecular_formula!(H 9 C 9 O 2 N 1),
            AminoAcid::Valine => molecular_formula!(H 9 C 5 O 1 N 1),
            AminoAcid::Unknown => molecular_formula!(),
            _ => unreachable!(),
        }
    }
}

impl<T> MultiChemical for CheckedAminoAcid<T> {
    /// # Panics
    /// Is the sequence index is a terminal index
    fn formulas(
        &self,
        sequence_index: crate::SequencePosition,
        peptide_index: usize,
    ) -> Multi<MolecularFormula> {
        if let Some(unambiguous) = self.is_unambiguous() {
            unambiguous.formula(sequence_index, peptide_index).into()
        } else {
            let crate::SequencePosition::Index(sequence_index) = sequence_index else {
                panic!("Not allowed to call amino acid formulas with a terminal sequence index")
            };
            match self.aminoacid {
            AminoAcid::AmbiguousAsparagine => vec![
                molecular_formula!(H 6 C 4 O 2 N 2 (crate::AmbiguousLabel::AminoAcid{option: AminoAcid::Asparagine, sequence_index, peptide_index})),
                molecular_formula!(H 5 C 4 O 3 N 1 (crate::AmbiguousLabel::AminoAcid{option: AminoAcid::AsparticAcid, sequence_index, peptide_index})),
            ]
            .into(),
            AminoAcid::AmbiguousGlutamine => vec![
                molecular_formula!(H 8 C 5 O 2 N 2 (crate::AmbiguousLabel::AminoAcid{option: AminoAcid::Glutamine, sequence_index, peptide_index})),
                molecular_formula!(H 7 C 5 O 3 N 1 (crate::AmbiguousLabel::AminoAcid{option: AminoAcid::GlutamicAcid, sequence_index, peptide_index})),
            ]
            .into(),
            _ => unreachable!(),
        }
        }
    }
}

impl<T> Copy for CheckedAminoAcid<T> {}

impl<T> Clone for CheckedAminoAcid<T> {
    fn clone(&self) -> Self {
        Self {
            aminoacid: self.aminoacid,
            marker: PhantomData,
        }
    }
}

impl<T> PartialEq for CheckedAminoAcid<T> {
    fn eq(&self, other: &Self) -> bool {
        self.aminoacid == other.aminoacid
    }
}

impl<T> std::hash::Hash for CheckedAminoAcid<T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.aminoacid.hash(state);
    }
}

impl<T> Eq for CheckedAminoAcid<T> {}

impl<T> Default for CheckedAminoAcid<T> {
    fn default() -> Self {
        Self {
            aminoacid: AminoAcid::default(),
            marker: PhantomData,
        }
    }
}

impl<T> std::fmt::Display for CheckedAminoAcid<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.char())
    }
}
