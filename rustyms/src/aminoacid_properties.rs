//! All amino acid property classes according to IMGT.
//! > IMGT standardized criteria for statistical analysis of immunoglobulin V-REGION amino acid properties
//! >
//! > Christelle Pommié, Séverine Levadoux, Robert Sabatier, Gérard Lefranc, Marie-Paule Lefranc
//! >
//! > <https://doi.org/10.1002/jmr.647>
//!
//! Adapted to include J,B,Z,U,O, and X.
#![allow(missing_docs)]

use serde::{Deserialize, Serialize};

use crate::{
    aminoacids::IsAminoAcid, modification::SimpleModification, AminoAcid, Peptidoform,
    SemiAmbiguous,
};

/// All amino acid property classes according to IMGT.
/// > IMGT standardized criteria for statistical analysis of immunoglobulin V-REGION amino acid properties
/// >
/// > Christelle Pommié, Séverine Levadoux, Robert Sabatier, Gérard Lefranc, Marie-Paule Lefranc
/// >
/// > <https://doi.org/10.1002/jmr.647>
///
/// Adapted to include J,B,Z,U,O, and X.
impl crate::AminoAcid {
    pub const fn physiochemical_class(self) -> PhysiochemicalClass {
        match self {
            Self::Alanine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Valine => PhysiochemicalClass::Aliphatic,
            Self::Proline => PhysiochemicalClass::Proline,
            Self::Tyrosine => PhysiochemicalClass::Tyrosine,
            Self::Tryptophan => PhysiochemicalClass::Tryptophan,
            Self::Glycine => PhysiochemicalClass::Glycine,
            Self::Phenylalanine => PhysiochemicalClass::Phenylalanine,
            Self::Cysteine | Self::Methionine => PhysiochemicalClass::Sulphur,
            Self::Serine | Self::Threonine => PhysiochemicalClass::Hydroxyl,
            Self::AsparticAcid | Self::GlutamicAcid => PhysiochemicalClass::Acidic,
            Self::Asparagine | Self::Glutamine => PhysiochemicalClass::Amide,
            Self::Arginine | Self::Histidine | Self::Lysine => PhysiochemicalClass::Basic,
            Self::AmbiguousAsparagine
            | Self::AmbiguousGlutamine
            | Self::Pyrrolysine
            | Self::Selenocysteine
            | Self::Unknown => PhysiochemicalClass::Unknown,
        }
    }
    pub const fn chemical_class(self) -> ChemicalClass {
        match self {
            Self::Alanine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Proline
            | Self::Glycine
            | Self::Valine => ChemicalClass::Aliphatic,
            Self::Tyrosine | Self::Tryptophan | Self::Phenylalanine => ChemicalClass::Aromatic,
            Self::Cysteine | Self::Methionine => ChemicalClass::Sulphur,
            Self::Serine | Self::Threonine => ChemicalClass::Hydroxyl,
            Self::AsparticAcid | Self::GlutamicAcid => ChemicalClass::Acidic,
            Self::Asparagine | Self::Glutamine => ChemicalClass::Amide,
            Self::Arginine | Self::Histidine | Self::Lysine => ChemicalClass::Basic,
            Self::AmbiguousAsparagine
            | Self::AmbiguousGlutamine
            | Self::Pyrrolysine
            | Self::Selenocysteine
            | Self::Unknown => ChemicalClass::Unknown,
        }
    }
    pub const fn volume_class(self) -> VolumeClass {
        match self {
            Self::Alanine | Self::Glycine | Self::Serine => VolumeClass::VerySmall,
            Self::Proline
            | Self::AsparticAcid
            | Self::Asparagine
            | Self::AmbiguousAsparagine
            | Self::Threonine
            | Self::Cysteine
            | Self::Selenocysteine => VolumeClass::Small,
            Self::GlutamicAcid
            | Self::Glutamine
            | Self::AmbiguousGlutamine
            | Self::Histidine
            | Self::Valine => VolumeClass::Medium,
            Self::Arginine
            | Self::Lysine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Methionine => VolumeClass::Large,
            Self::Tyrosine | Self::Tryptophan | Self::Phenylalanine | Self::Pyrrolysine => {
                VolumeClass::VeryLarge
            }
            Self::Unknown => VolumeClass::Unknown,
        }
    }
    pub const fn hydropathy_class(self) -> HydropathyClass {
        match self {
            Self::Alanine
            | Self::Cysteine
            | Self::Selenocysteine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Methionine
            | Self::Phenylalanine
            | Self::Tryptophan
            | Self::Valine => HydropathyClass::Hydrophobic,
            Self::Glycine
            | Self::Serine
            | Self::Histidine
            | Self::Proline
            | Self::Threonine
            | Self::Tyrosine => HydropathyClass::Neutral,
            Self::Lysine
            | Self::AsparticAcid
            | Self::Asparagine
            | Self::AmbiguousAsparagine
            | Self::GlutamicAcid
            | Self::Glutamine
            | Self::AmbiguousGlutamine
            | Self::Pyrrolysine
            | Self::Arginine => HydropathyClass::Hydrophilic,
            Self::Unknown => HydropathyClass::Unknown,
        }
    }
    pub const fn charge_class(self) -> ChargeClass {
        match self {
            Self::Arginine | Self::Histidine | Self::Lysine => ChargeClass::Positive,
            Self::AsparticAcid | Self::GlutamicAcid => ChargeClass::Negative,
            Self::Alanine
            | Self::Cysteine
            | Self::Selenocysteine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Methionine
            | Self::Phenylalanine
            | Self::Tryptophan
            | Self::Valine
            | Self::Glycine
            | Self::Serine
            | Self::Proline
            | Self::Threonine
            | Self::Tyrosine
            | Self::Asparagine
            | Self::Glutamine
            | Self::Pyrrolysine => ChargeClass::Uncharged,
            Self::Unknown | Self::AmbiguousAsparagine | Self::AmbiguousGlutamine => {
                ChargeClass::Unknown
            }
        }
    }
    pub const fn polarity_class(self) -> PolarityClass {
        match self {
            Self::Arginine
            | Self::Histidine
            | Self::Lysine
            | Self::Pyrrolysine
            | Self::AmbiguousAsparagine
            | Self::AmbiguousGlutamine
            | Self::AsparticAcid
            | Self::Asparagine
            | Self::Glutamine
            | Self::Threonine
            | Self::Tyrosine
            | Self::Serine
            | Self::GlutamicAcid => PolarityClass::Polar,
            Self::Alanine
            | Self::Cysteine
            | Self::Selenocysteine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Methionine
            | Self::Phenylalanine
            | Self::Tryptophan
            | Self::Valine
            | Self::Glycine
            | Self::Proline => PolarityClass::Nonpolar,
            Self::Unknown => PolarityClass::Unknown,
        }
    }
    pub const fn hydrogen_bond_class(self) -> HydrogenBondClass {
        match self {
            Self::Arginine | Self::Lysine | Self::Tryptophan => HydrogenBondClass::Donor,
            Self::AsparticAcid | Self::GlutamicAcid => HydrogenBondClass::Acceptor,
            Self::Histidine
            | Self::Pyrrolysine
            | Self::Asparagine
            | Self::Threonine
            | Self::Tyrosine
            | Self::Serine => HydrogenBondClass::Both,
            Self::Glutamine
            | Self::Alanine
            | Self::Cysteine
            | Self::Selenocysteine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Methionine
            | Self::Phenylalanine
            | Self::Valine
            | Self::Glycine
            | Self::Proline => HydrogenBondClass::None,
            Self::AmbiguousAsparagine | Self::AmbiguousGlutamine | Self::Unknown => {
                HydrogenBondClass::Unknown
            }
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum PhysiochemicalClass {
    Aliphatic,
    Proline,
    Tyrosine,
    Tryptophan,
    Glycine,
    Phenylalanine,
    Sulphur,
    Hydroxyl,
    Acidic,
    Amide,
    Basic,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum ChemicalClass {
    Aliphatic,
    Sulphur,
    Hydroxyl,
    Acidic,
    Amide,
    Basic,
    Aromatic,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum VolumeClass {
    VerySmall,
    Small,
    Medium,
    Large,
    VeryLarge,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum HydropathyClass {
    Hydrophobic,
    Neutral,
    Hydrophilic,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum ChargeClass {
    Uncharged,
    Positive,
    Negative,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum PolarityClass {
    Polar,
    Nonpolar,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum HydrogenBondClass {
    Donor,
    Acceptor,
    Both,
    None,
    Unknown,
}

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
pub struct AminoAcidpKa {
    ammonium: f32,
    sidechain: Option<f32>,
    carboxyl: f32,
}

/// A source for pKa values, which can be used to calculate the pKa for peptidoforms.
pub trait pKaSource<AA: IsAminoAcid> {
    /// Get the pKa values for the given amino acid and modifications.
    fn pKa(amino_acid: AA, modifications: &[SimpleModification]) -> Option<AminoAcidpKa>;

    /// Get the calculated pKa value for the given peptidoform, or None if any of the sequence elements do not have a defined pKa.
    fn peptide_pKa(peptidoform: Peptidoform<SemiAmbiguous>) -> Option<f32> {
        todo!()
    }
}

/// The pKa values for an amino acid
impl AminoAcidpKa {
    const fn new(ammonium: f32, sidechain: Option<f32>, carboxyl: f32) -> Self {
        Self {
            ammonium,
            sidechain,
            carboxyl,
        }
    }
    pub const fn ammonium(self) -> f32 {
        self.ammonium
    }
    pub const fn sidechain(self) -> Option<f32> {
        self.sidechain
    }
    pub const fn carboxyl(self) -> f32 {
        self.carboxyl
    }
}

/// pKa values from Lide, D. R. (1991). Handbook of Chemistry and Physics: A Ready Reference Book of Chemical and Physical Data.
pub struct pKaLide1991;

impl pKaSource<AminoAcid> for pKaLide1991 {
    fn pKa(amino_acid: AminoAcid, modifications: &[SimpleModification]) -> Option<AminoAcidpKa> {
        if !modifications.is_empty() {
            return None;
        }
        match amino_acid {
            AminoAcid::Arginine => Some(AminoAcidpKa::new(9.00, Some(12.10), 2.03)),
            AminoAcid::Histidine => Some(AminoAcidpKa::new(9.09, Some(6.04), 1.70)),
            AminoAcid::Lysine => Some(AminoAcidpKa::new(9.16, Some(10.67), 2.15)),
            AminoAcid::AsparticAcid => Some(AminoAcidpKa::new(9.66, Some(3.71), 1.95)),
            AminoAcid::GlutamicAcid => Some(AminoAcidpKa::new(9.58, Some(4.15), 2.16)),
            AminoAcid::Tyrosine => Some(AminoAcidpKa::new(9.04, Some(10.10), 2.24)),
            AminoAcid::Cysteine => Some(AminoAcidpKa::new(10.28, Some(8.14), 1.91)),
            AminoAcid::Alanine => Some(AminoAcidpKa::new(9.71, None, 2.33)),
            AminoAcid::Glycine => Some(AminoAcidpKa::new(9.58, None, 2.34)),
            AminoAcid::Proline => Some(AminoAcidpKa::new(10.47, None, 1.95)),
            AminoAcid::Serine => Some(AminoAcidpKa::new(9.05, None, 2.13)),
            AminoAcid::Threonine => Some(AminoAcidpKa::new(8.96, None, 2.20)),
            AminoAcid::Methionine => Some(AminoAcidpKa::new(9.08, None, 2.16)),
            AminoAcid::Phenylalanine => Some(AminoAcidpKa::new(9.09, None, 2.18)),
            AminoAcid::Tryptophan => Some(AminoAcidpKa::new(9.34, None, 2.38)),
            AminoAcid::Valine => Some(AminoAcidpKa::new(9.52, None, 2.27)),
            AminoAcid::Isoleucine => Some(AminoAcidpKa::new(9.60, None, 2.26)),
            AminoAcid::Leucine => Some(AminoAcidpKa::new(9.58, None, 2.32)),
            AminoAcid::Glutamine => Some(AminoAcidpKa::new(9.00, None, 2.18)),
            AminoAcid::Asparagine => Some(AminoAcidpKa::new(8.73, None, 2.16)),
            AminoAcid::AmbiguousAsparagine => Some(AminoAcidpKa::new(8.73, None, 2.16)),
            AminoAcid::AmbiguousGlutamine => Some(AminoAcidpKa::new(9.00, None, 2.18)),
            _ => None,
        }
    }
}

/// pKa values from Lehninger, A. L., Nelson, D. L., & Cox, M. M. (2005). Lehninger Principles of Biochemistry. Macmillan.
pub struct pKaLehninger;

impl pKaSource<AminoAcid> for pKaLehninger {
    fn pKa(amino_acid: AminoAcid, modifications: &[SimpleModification]) -> Option<AminoAcidpKa> {
        if !modifications.is_empty() {
            return None;
        }
        match amino_acid {
            AminoAcid::Arginine => Some(AminoAcidpKa::new(9.04, Some(12.48), 2.17)),
            AminoAcid::Histidine => Some(AminoAcidpKa::new(9.17, Some(6.00), 1.82)),
            AminoAcid::Lysine => Some(AminoAcidpKa::new(8.95, Some(10.53), 2.18)),
            AminoAcid::AsparticAcid => Some(AminoAcidpKa::new(9.60, Some(3.65), 1.88)),
            AminoAcid::GlutamicAcid => Some(AminoAcidpKa::new(9.67, Some(4.25), 2.19)),
            AminoAcid::Tyrosine => Some(AminoAcidpKa::new(9.11, Some(10.07), 2.20)),
            AminoAcid::Cysteine => Some(AminoAcidpKa::new(10.28, Some(8.18), 1.96)),
            AminoAcid::Alanine => Some(AminoAcidpKa::new(9.69, None, 2.34)),
            AminoAcid::Glycine => Some(AminoAcidpKa::new(9.60, None, 2.34)),
            AminoAcid::Proline => Some(AminoAcidpKa::new(10.96, None, 1.99)),
            AminoAcid::Serine => Some(AminoAcidpKa::new(9.15, None, 2.21)),
            AminoAcid::Threonine => Some(AminoAcidpKa::new(9.62, None, 2.11)),
            AminoAcid::Methionine => Some(AminoAcidpKa::new(9.21, None, 2.28)),
            AminoAcid::Phenylalanine => Some(AminoAcidpKa::new(9.13, None, 1.83)),
            AminoAcid::Tryptophan => Some(AminoAcidpKa::new(9.39, None, 2.38)),
            AminoAcid::Valine => Some(AminoAcidpKa::new(9.62, None, 2.32)),
            AminoAcid::Isoleucine => Some(AminoAcidpKa::new(9.68, None, 2.36)),
            AminoAcid::Leucine => Some(AminoAcidpKa::new(9.60, None, 2.36)),
            AminoAcid::Glutamine => Some(AminoAcidpKa::new(9.13, None, 2.17)),
            AminoAcid::Asparagine => Some(AminoAcidpKa::new(8.80, None, 2.02)),
            AminoAcid::AmbiguousAsparagine => Some(AminoAcidpKa::new(8.80, None, 2.02)),
            AminoAcid::AmbiguousGlutamine => Some(AminoAcidpKa::new(9.13, None, 2.17)),
            _ => None,
        }
    }
}
