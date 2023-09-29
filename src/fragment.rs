use std::fmt::{Debug, Display};

use itertools::Itertools;
use uom::num_traits::Zero;

use crate::{
    molecular_charge::MolecularCharge, system::f64::*, AminoAcid, Chemical, Element,
    MolecularFormula, NeutralLoss,
};

/// A theoretical fragment of a peptide
#[derive(Debug, Clone)]
pub struct Fragment {
    /// The theoretical composition
    pub formula: MolecularFormula,
    /// The charge
    pub charge: Charge,
    /// All possible annotations for this fragment saved as a tuple of peptide index and its type
    pub ion: FragmentType,
    /// The peptide this fragment comes from, saved as the index into the list of peptides in the overarching [`ComplexPeptide`] struct
    pub peptide_index: usize,
    /// Any neutral losses applied
    pub neutral_loss: Option<NeutralLoss>,
    /// Additional description for humans
    pub label: String,
}

impl Fragment {
    /// Get the mz
    pub fn mz(&self) -> Option<MassOverCharge> {
        Some(self.formula.monoisotopic_mass()? / self.charge)
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
            formula: theoretical_mass,
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
        annotation: &FragmentType,
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
                    annotation.clone(),
                    term.1.to_string(),
                )
            })
            .flat_map(|m| m.with_neutral_losses(neutral_losses))
            .collect()
    }

    /// Create a copy of this fragment with the given charge
    #[must_use]
    pub fn with_charge(&self, charge: &MolecularCharge) -> Self {
        // TODO: Figure out if labelling these in any way would be nice for later checking when used with adduct ions other than protons
        let formula = charge.formula();
        let c = Charge::new::<e>(f64::from(formula.charge()));
        Self {
            formula: &self.formula + &formula,
            charge: c,
            ..self.clone()
        }
    }

    /// Create a copy of this fragment with the given neutral loss
    #[must_use]
    pub fn with_neutral_loss(&self, neutral_loss: &NeutralLoss) -> Self {
        Self {
            formula: &self.formula + neutral_loss,
            neutral_loss: Some(neutral_loss.clone()),
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
            "{}@{}{:+}{} {}",
            self.ion,
            self.mz()
                .map_or("Undefined".to_string(), |m| m.value.to_string()),
            self.charge.value,
            self.neutral_loss
                .as_ref()
                .map(std::string::ToString::to_string)
                .unwrap_or_default(),
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

/// The definition of the position of an ion inside a glycan
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct GlycanPosition {
    /// The depth starting at the amino acid
    pub inner_depth: usize,
    /// The series number (from the ion series terminal)
    pub series_number: usize,
    /// The branch naming
    pub branch: Vec<usize>,
    /// The aminoacid index where this glycan is attached
    pub attachment: (AminoAcid, usize),
}

impl GlycanPosition {
    /// Get the branch names
    /// # Panics
    /// Panics if the first branch number is outside the range of the greek alphabet (small and caps together).
    pub fn branch_names(&self) -> String {
        self.branch
            .iter()
            .enumerate()
            .map(|(i, b)| {
                if i == 0 {
                    char::from_u32(
                        (0x03B1..=0x03C9)
                            .chain(0x0391..=0x03A9)
                            .nth(*b)
                            .expect("Too many branches in glycan, out of greek letters"),
                    )
                    .unwrap()
                    .to_string()
                } else if i == 1 {
                    "\'".repeat(*b)
                } else {
                    format!(",{b}")
                }
            })
            .collect::<String>()
    }
    /// Generate the label for this glycan position, example: `1α'`
    /// # Panics
    /// Panics if the first branch number is outside the range of the greek alphabet (small and caps together).
    pub fn label(&self) -> String {
        format!("{}{}", self.series_number, self.branch_names())
    }
    /// Generate the label for this glycan attachment eg N1 (1 based numbering)
    pub fn attachment(&self) -> String {
        format!("{}{}", self.attachment.0.char(), self.attachment.1 + 1)
    }
}

/// The possible types of fragments
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
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
    /// glycan A fragment
    A(GlycanPosition),
    /// glycan B fragment
    B(GlycanPosition),
    /// glycan C fragment
    C(GlycanPosition),
    /// glycan X fragment
    X(GlycanPosition),
    /// glycan Y fragment
    Y(GlycanPosition),
    /// glycan Z fragment
    Z(GlycanPosition),
    /// glycan Z fragment
    InternalGlycan(Vec<GlycanBreakPos>),
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
            _ => None,
        }
    }

    /// Get the glycan position of this ion (or None nor applicable)
    pub const fn glycan_position(&self) -> Option<&GlycanPosition> {
        match self {
            Self::A(n) | Self::B(n) | Self::C(n) | Self::X(n) | Self::Y(n) | Self::Z(n) => Some(n),
            _ => None,
        }
    }

    /// Get the position label, unless it is a precursor ion
    pub fn position_label(&self) -> Option<String> {
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
            | Self::z·(n) => Some(n.series_number.to_string()),
            Self::A(n) | Self::B(n) | Self::C(n) | Self::X(n) | Self::Y(n) | Self::Z(n) => {
                Some(n.label())
            }
            Self::InternalGlycan(breakages) => Some(
                breakages
                    .iter()
                    .map(std::string::ToString::to_string)
                    .join(""),
            ),
            Self::precursor => None,
        }
    }

    /// Get the label for this fragment type
    pub const fn label(&self) -> &str {
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
            Self::A(_) => "A",
            Self::B(_) => "B",
            Self::C(_) => "C",
            Self::X(_) => "X",
            Self::Y(_) => "Y",
            Self::Z(_) => "Z",
            Self::InternalGlycan(_) => "internal_glycan",
            Self::precursor => "precursor",
        }
    }
}

impl Display for FragmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::a(pos) => format!("a{}", pos.series_number),
                Self::b(pos) => format!("b{}", pos.series_number),
                Self::c(pos) => format!("c{}", pos.series_number),
                Self::d(pos) => format!("d{}", pos.series_number),
                Self::v(pos) => format!("v{}", pos.series_number),
                Self::w(pos) => format!("w{}", pos.series_number),
                Self::x(pos) => format!("x{}", pos.series_number),
                Self::y(pos) => format!("y{}", pos.series_number),
                Self::z(pos) => format!("z{}", pos.series_number),
                Self::z·(pos) => format!("z·{}", pos.series_number),
                Self::A(pos) => format!("A{}", pos.label()),
                Self::B(pos) => format!("B{}", pos.label()),
                Self::C(pos) => format!("C{}", pos.label()),
                Self::X(pos) => format!("X{}", pos.label()),
                Self::Y(pos) => format!("Y{}", pos.label()),
                Self::Z(pos) => format!("Z{}", pos.label()),
                Self::InternalGlycan(positions) => positions
                    .iter()
                    .map(std::string::ToString::to_string)
                    .collect(),
                Self::precursor => "precursor".to_string(),
            }
        )
    }
}

/// All positions where a glycan can break
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub enum GlycanBreakPos {
    /// No breaks just until the end of a chain
    End(GlycanPosition),
    /// Break at a Y position
    Y(GlycanPosition),
    /// Break at a B position
    B(GlycanPosition),
}

impl GlycanBreakPos {
    /// Get the position of this breaking position
    pub const fn position(&self) -> &GlycanPosition {
        match self {
            Self::B(p) | Self::End(p) | Self::Y(p) => p,
        }
    }

    /// Get the label for this breaking position
    pub const fn label(&self) -> &str {
        match self {
            Self::End(_) => "End",
            Self::Y(_) => "Y",
            Self::B(_) => "B",
        }
    }
}

impl std::fmt::Display for GlycanBreakPos {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.label(), self.position().label())
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
        let loss = a.with_neutral_losses(&[NeutralLoss::Loss(molecular_formula!(H 2 O 1))]);
        dbg!(&a, &loss);
        assert_eq!(a.formula, loss[0].formula);
        assert_eq!(a.formula, &loss[1].formula + &molecular_formula!(H 2 O 1));
    }
}
