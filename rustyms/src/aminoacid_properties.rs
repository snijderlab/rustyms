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

use crate::AminoAcid;

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

    pub fn pka(&self) -> pKa {
        self.pka_with_source(pksources::Lide1991)
    }

    pub fn pka_with_source(self, source: pksources) -> pKa {
        let table = match source {
            pksources::Lide1991 => PKA_LIDE1991,
            pksources::Lehninger => PKA_LEHNINGER,
        };

        table
            .iter()
            .find(|&&(amino_acid, _)| amino_acid == self)
            .map(|&(_, pka)| pka)
            .expect("Peptide not found in the specified table")
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
pub struct pKa {
    carboxyl: f32,
    ammonium: f32,
    sidechain: Option<f32>,
}

impl pKa {
    const fn from(carboxyl: f32, ammonium: f32, sidechain: Option<f32>) -> Self {
        Self {
            carboxyl,
            ammonium,
            sidechain,
        }
    }
    pub const fn carboxyl(self) -> f32 {
        self.carboxyl
    }
    pub const fn ammonium(self) -> f32 {
        self.ammonium
    }
    pub const fn sidechain(self) -> Option<f32> {
        self.sidechain
    }
}

pub enum pksources {
    Lide1991,
    Lehninger,
}

// pKa values from Lide, D. R. (1991). Handbook of Chemistry and Physics: A Ready Reference Book of Chemical and Physical Data.
const PKA_LIDE1991: &[(AminoAcid, pKa)] = &[
    (AminoAcid::Arginine, pKa::from(2.03, 9.00, Some(12.10))),
    (AminoAcid::Histidine, pKa::from(1.70, 9.09, Some(6.04))),
    (AminoAcid::Lysine, pKa::from(2.15, 9.16, Some(10.67))),
    (AminoAcid::AsparticAcid, pKa::from(1.95, 9.66, Some(3.71))),
    (AminoAcid::GlutamicAcid, pKa::from(2.16, 9.58, Some(4.15))),
    (AminoAcid::Tyrosine, pKa::from(2.24, 9.04, Some(10.10))),
    (AminoAcid::Cysteine, pKa::from(1.91, 10.28, Some(8.14))),
    (AminoAcid::Alanine, pKa::from(2.33, 9.71, None)),
    (AminoAcid::Glycine, pKa::from(2.34, 9.58, None)),
    (AminoAcid::Proline, pKa::from(1.95, 10.47, None)),
    (AminoAcid::Serine, pKa::from(2.13, 9.05, None)),
    (AminoAcid::Threonine, pKa::from(2.20, 8.96, None)),
    (AminoAcid::Methionine, pKa::from(2.16, 9.08, None)),
    (AminoAcid::Phenylalanine, pKa::from(2.18, 9.09, None)),
    (AminoAcid::Tryptophan, pKa::from(2.38, 9.34, None)),
    (AminoAcid::Valine, pKa::from(2.27, 9.52, None)),
    (AminoAcid::Isoleucine, pKa::from(2.26, 9.60, None)),
    (AminoAcid::Leucine, pKa::from(2.32, 9.58, None)),
    (AminoAcid::Glutamine, pKa::from(2.18, 9.00, None)),
    (AminoAcid::Asparagine, pKa::from(2.16, 8.73, None)),
    (AminoAcid::AmbiguousAsparagine, pKa::from(2.16, 8.73, None)),
    (AminoAcid::AmbiguousGlutamine, pKa::from(2.18, 9.00, None)),
    // (AminoAcid::Pyrrolysine, todo!()),
    // (AminoAcid::Unknown, todo!()),
];

// pKa values from Lehninger, A. L., Nelson, D. L., & Cox, M. M. (2005). Lehninger Principles of Biochemistry. Macmillan.
const PKA_LEHNINGER: &[(AminoAcid, pKa)] = &[
    (AminoAcid::Arginine, pKa::from(2.17, 9.04, Some(12.48))),
    (AminoAcid::Histidine, pKa::from(1.82, 9.17, Some(6.00))),
    (AminoAcid::Lysine, pKa::from(2.18, 8.95, Some(10.53))),
    (AminoAcid::AsparticAcid, pKa::from(1.88, 9.60, Some(3.65))),
    (AminoAcid::GlutamicAcid, pKa::from(2.19, 9.67, Some(4.25))),
    (AminoAcid::Tyrosine, pKa::from(2.20, 9.11, Some(10.07))),
    (AminoAcid::Cysteine, pKa::from(1.96, 10.28, Some(8.18))),
    (AminoAcid::Alanine, pKa::from(2.34, 9.69, None)),
    (AminoAcid::Glycine, pKa::from(2.34, 9.60, None)),
    (AminoAcid::Proline, pKa::from(1.99, 10.96, None)),
    (AminoAcid::Serine, pKa::from(2.21, 9.15, None)),
    (AminoAcid::Threonine, pKa::from(2.11, 9.62, None)),
    (AminoAcid::Methionine, pKa::from(2.28, 9.21, None)),
    (AminoAcid::Phenylalanine, pKa::from(1.83, 9.13, None)),
    (AminoAcid::Tryptophan, pKa::from(2.38, 9.39, None)),
    (AminoAcid::Valine, pKa::from(2.32, 9.62, None)),
    (AminoAcid::Isoleucine, pKa::from(2.36, 9.68, None)),
    (AminoAcid::Leucine, pKa::from(2.36, 9.60, None)),
    (AminoAcid::Glutamine, pKa::from(2.17, 9.13, None)),
    (AminoAcid::Asparagine, pKa::from(2.02, 8.80, None)),
    (AminoAcid::AmbiguousAsparagine, pKa::from(2.02, 8.80, None)),
    (AminoAcid::AmbiguousGlutamine, pKa::from(2.17, 9.13, None)),
    // (AminoAcid::Pyrrolysine, todo!()),
    // (AminoAcid::Unknown, todo!()),
];
