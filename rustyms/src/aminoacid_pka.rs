use serde::{Deserialize, Serialize};

use crate::{
    aminoacids::IsAminoAcid, modification::SimpleModification, AminoAcid, Peptidoform,
    SemiAmbiguous,
};

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
pub struct AminoAcidpKa {
    ammonium: f64,
    sidechain: Option<f64>,
    carboxyl: f64,
}

/// A source for pKa values, which can be used to calculate the pKa for peptidoforms.
pub trait pKaSource<AA: IsAminoAcid> {
    /// Get the pKa values for the given amino acid and modifications.
    fn pKa(amino_acid: AA, modifications: &[SimpleModification]) -> Option<AminoAcidpKa>;

    /// Get the calculated pKa value for the given peptidoform, or None if any of the sequence elements do not have a defined pKa.
    fn peptide_pKa(peptidoform: Peptidoform<SemiAmbiguous>) -> Option<f64> {
        todo!()
    }
}

/// The pKa values for an amino acid
impl AminoAcidpKa {
    const fn new(ammonium: f64, sidechain: Option<f64>, carboxyl: f64) -> Self {
        Self {
            ammonium,
            sidechain,
            carboxyl,
        }
    }
    pub const fn ammonium(self) -> f64 {
        self.ammonium
    }
    pub const fn sidechain(self) -> Option<f64> {
        self.sidechain
    }
    pub const fn carboxyl(self) -> f64 {
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
