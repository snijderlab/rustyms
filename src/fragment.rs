use std::fmt::{Debug, Display};

use uom::num_traits::Zero;

use crate::{system::f64::*, Chemical, Element, MolecularFormula, NeutralLoss};

/// A theoretical fragment of a peptide
#[derive(Debug, Clone)]
pub struct Fragment {
    /// The theoretical composition
    pub theoretical_mass: MolecularFormula,
    /// The charge
    pub charge: Charge,
    /// All possible annotations for this fragment saved as a tuple of peptide index and its type
    pub ion: FragmentType,
    pub peptide_index: usize,
    /// Any neutral losses applied
    pub neutral_loss: Option<NeutralLoss>,
    /// Additional description for humans
    pub label: String,
}

impl Fragment {
    /// Get the mz
    pub fn mz(&self) -> Option<MassOverCharge> {
        Some(self.theoretical_mass.monoisotopic_mass()? / self.charge)
    }

    /// Get the ppm difference between two fragments
    pub fn ppm(&self, other: &Self) -> Option<MassOverCharge> {
        Some(MassOverCharge::new::<mz>(self.mz()?.ppm(other.mz()?)))
    }

    /// Create a new fragment
    #[must_use]
    pub fn new(
        theoretical_mass: MolecularFormula,
        charge: Charge,
        peptide_index: usize,
        ion: FragmentType,
        label: String,
    ) -> Self {
        Self {
            theoretical_mass,
            charge,
            ion,
            peptide_index,
            label,
            neutral_loss: None,
        }
    }

    /// Generate a list of possible fragments from the list of possible preceding termini and neutral losses
    #[must_use]
    pub fn generate_all(
        theoretical_mass: &MolecularFormula,
        peptide_index: usize,
        annotation: FragmentType,
        termini: &[(MolecularFormula, String)],
        neutral_losses: &[NeutralLoss],
    ) -> Vec<Self> {
        termini
            .iter()
            .map(|term| {
                Self::new(
                    &term.0 + theoretical_mass,
                    Charge::zero(),
                    peptide_index,
                    annotation,
                    term.1.to_string(),
                )
            })
            .flat_map(|m| m.with_neutral_losses(neutral_losses))
            .collect()
    }

    /// Create a copy of this fragment with the given charge
    #[must_use]
    pub fn with_charge(&self, charge: Charge) -> Self {
        let c = charge.value as i16;
        Self {
            theoretical_mass: &self.theoretical_mass + &molecular_formula!(H c)
                - molecular_formula!(Electron c),
            charge,
            ..self.clone()
        }
    }

    /// Create a copy of this fragment with the given neutral loss
    #[must_use]
    pub fn with_neutral_loss(&self, neutral_loss: &NeutralLoss) -> Self {
        Self {
            theoretical_mass: &self.theoretical_mass - &neutral_loss.formula(),
            neutral_loss: Some(*neutral_loss),
            ..self.clone()
        }
    }

    /// Create copies of this fragment with the given neutral losses (and a copy of this fragment itself)
    #[must_use]
    pub fn with_neutral_losses(&self, neutral_losses: &[NeutralLoss]) -> Vec<Self> {
        let mut output = Vec::with_capacity(neutral_losses.len() + 1);
        output.push(self.clone());
        output.extend(
            neutral_losses
                .iter()
                .map(|loss| self.with_neutral_loss(loss)),
        );
        output
    }
}

impl Display for Fragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:?} {:?} {:+}{} {}",
            self.ion,
            self.mz()
                .map_or("Undefined".to_string(), |m| m.value.to_string()),
            self.charge.value,
            self.neutral_loss
                .map_or(String::new(), |loss| format!(" -{loss}")),
            self.label
        )
    }
}

/// The definition of the position of an ion
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct Position {
    /// The sequence index (0 based into the peptide sequence)
    pub sequence_index: usize,
    /// The series number (1 based from the ion series terminal)
    pub series_number: usize,
}

impl Position {
    /// Generate a position for N terminal ion series
    pub const fn n(sequence_index: usize, _length: usize) -> Self {
        Self {
            sequence_index,
            series_number: sequence_index + 1,
        }
    }
    /// Generate a position for C terminal ion series
    pub const fn c(sequence_index: usize, length: usize) -> Self {
        Self {
            sequence_index,
            series_number: length - sequence_index,
        }
    }
}

/// The possible types of fragments
#[derive(Clone, Copy, Eq, PartialEq, Hash, Debug)]
#[allow(non_camel_case_types)]
pub enum FragmentType {
    /// a
    a(Position),
    /// b
    b(Position),
    /// c
    c(Position),
    /// d
    d(Position),
    /// v
    v(Position),
    /// w
    w(Position),
    /// x
    x(Position),
    /// y
    y(Position),
    /// z
    z(Position),
    /// z·
    z·(Position),
    /// precursor
    precursor,
}

impl FragmentType {
    /// Get the position of this ion (or None if it is a precursor ion)
    pub const fn position(&self) -> Option<&Position> {
        match self {
            Self::a(n)
            | Self::b(n)
            | Self::c(n)
            | Self::d(n)
            | Self::v(n)
            | Self::w(n)
            | Self::x(n)
            | Self::y(n)
            | Self::z(n)
            | Self::z·(n) => Some(n),
            Self::precursor => None,
        }
    }
}

impl Display for FragmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::a(_) => "a",
                Self::b(_) => "b",
                Self::c(_) => "c",
                Self::d(_) => "d",
                Self::v(_) => "v",
                Self::w(_) => "w",
                Self::x(_) => "x",
                Self::y(_) => "y",
                Self::z(_) => "z",
                Self::z·(_) => "z·",
                Self::precursor => "precursor",
            }
        )
    }
}

#[cfg(test)]
mod tests {

    use crate::AminoAcid;

    use super::*;

    #[test]
    fn neutral_loss() {
        let a = Fragment::new(
            AminoAcid::AsparticAcid.formula(),
            Charge::new::<e>(1.0),
            0,
            FragmentType::precursor,
            String::new(),
        );
        let loss = a.with_neutral_losses(&[NeutralLoss::Water]);
        dbg!(&a, &loss);
        assert_eq!(a.theoretical_mass, loss[0].theoretical_mass);
        assert_eq!(
            a.theoretical_mass,
            &loss[1].theoretical_mass + &NeutralLoss::Water.formula()
        );
    }
}
