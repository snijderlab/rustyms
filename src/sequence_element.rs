#![warn(dead_code)]

use crate::{
    error::{Context, CustomError},
    modification::AmbiguousModification,
    MolecularFormula, MultiChemical, MultiMolecularFormula,
};
use serde::{Deserialize, Serialize};

use crate::{aminoacids::AminoAcid, modification::Modification, Chemical};

/// One block in a sequence meaning an aminoacid and its accompanying modifications
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct SequenceElement {
    /// The aminoacid
    pub aminoacid: AminoAcid,
    /// All present modifications
    pub modifications: Vec<Modification>,
    /// All ambiguous modifications (could be placed here or on another position)
    pub possible_modifications: Vec<AmbiguousModification>,
    /// If this aminoacid is part of an ambiguous sequence group `(QA)?` in pro forma
    pub ambiguous: Option<usize>,
}

impl SequenceElement {
    /// Create a new aminoacid without any modifications
    pub const fn new(aminoacid: AminoAcid, ambiguous: Option<usize>) -> Self {
        Self {
            aminoacid,
            modifications: Vec::new(),
            possible_modifications: Vec::new(),
            ambiguous,
        }
    }

    pub(crate) fn display(
        &self,
        f: &mut std::fmt::Formatter<'_>,
        placed: &[usize],
        last_ambiguous: Option<usize>,
    ) -> Result<Vec<usize>, std::fmt::Error> {
        let mut extra_placed = Vec::new();
        if last_ambiguous.is_some() && last_ambiguous != self.ambiguous {
            write!(f, ")")?;
        }
        if self.ambiguous.is_some() && last_ambiguous != self.ambiguous {
            write!(f, "(?")?;
        }
        write!(f, "{}", self.aminoacid.char())?;
        for m in &self.modifications {
            write!(f, "[{m}]")?;
        }
        for m in &self.possible_modifications {
            write!(
                f,
                "[{}#{}{}]",
                m.group.as_ref().map_or(
                    if placed.contains(&m.id) {
                        String::new()
                    } else {
                        extra_placed.push(m.id);
                        m.modification.to_string()
                    },
                    |group| if group.1 {
                        m.modification.to_string()
                    } else {
                        String::new()
                    }
                ),
                m.group
                    .as_ref()
                    .map_or(m.id.to_string(), |g| g.0.to_string()),
                m.localisation_score
                    .map(|v| format!("({v})"))
                    .unwrap_or_default()
            )?;
        }
        Ok(extra_placed)
    }

    /// Get the molecular formula for this position (unless it is B/Z) with the selected ambiguous modifications, without any global isotype modifications
    pub fn formula(&self, selected_ambiguous: &[usize]) -> Option<MolecularFormula> {
        if self.aminoacid == AminoAcid::B || self.aminoacid == AminoAcid::Z {
            None
        } else {
            Some(
                self.aminoacid.formula()
                    + self
                        .modifications
                        .iter()
                        .map(Chemical::formula)
                        .sum::<MolecularFormula>()
                    + self
                        .possible_modifications
                        .iter()
                        .filter(|&m| selected_ambiguous.contains(&m.id))
                        .map(|m| m.modification.formula())
                        .sum::<MolecularFormula>(),
            )
        }
    }

    /// Get the molecular formula for this position (unless it is B/Z) with the ambiguous modifications placed on the very first placed (and updating this in `placed`), without any global isotype modifications
    pub fn formula_greedy(&self, placed: &mut [bool]) -> Option<MolecularFormula> {
        if self.aminoacid == AminoAcid::B || self.aminoacid == AminoAcid::Z {
            None
        } else {
            #[allow(clippy::filter_map_bool_then)] // otherwise crashes
            Some(
                self.aminoacid.formula()
                    + self
                        .modifications
                        .iter()
                        .map(Chemical::formula)
                        .sum::<MolecularFormula>()
                    + self
                        .possible_modifications
                        .iter()
                        .filter_map(|m| {
                            (!placed[m.id]).then(|| {
                                placed[m.id] = true;
                                m.modification.formula()
                            })
                        })
                        .sum::<MolecularFormula>(),
            )
        }
    }

    /// Get the molecular formula for this position (unless it is B/Z) with all ambiguous modifications, without any global isotype modifications
    pub fn formula_all(&self) -> Option<MolecularFormula> {
        if self.aminoacid == AminoAcid::B || self.aminoacid == AminoAcid::Z {
            None
        } else {
            Some(
                self.aminoacid.formula()
                    + self
                        .modifications
                        .iter()
                        .map(Chemical::formula)
                        .sum::<MolecularFormula>()
                    + self
                        .possible_modifications
                        .iter()
                        .map(|m| m.modification.formula())
                        .sum::<MolecularFormula>(),
            )
        }
    }

    /// Get the molecular formulas for this position with the selected ambiguous modifications, without any global isotype modifications
    pub fn formulas(&self, selected_ambiguous: &[usize]) -> MultiMolecularFormula {
        let modifications = self
            .modifications
            .iter()
            .map(Chemical::formula)
            .sum::<MolecularFormula>()
            + self
                .possible_modifications
                .iter()
                .filter(|&m| selected_ambiguous.contains(&m.id))
                .map(|m| m.modification.formula())
                .sum::<MolecularFormula>();
        self.aminoacid
            .formulas()
            .iter()
            .map(|formula| formula + &modifications)
            .collect()
    }

    /// Get the molecular formulas for this position with the ambiguous modifications placed on the very first placed (and updating this in `placed`), without any global isotype modifications
    pub fn formulas_greedy(&self, placed: &mut [bool]) -> MultiMolecularFormula {
        #[allow(clippy::filter_map_bool_then)] // otherwise crashes
        let modifications = self
            .modifications
            .iter()
            .map(Chemical::formula)
            .sum::<MolecularFormula>()
            + self
                .possible_modifications
                .iter()
                .filter_map(|m| {
                    (!placed[m.id]).then(|| {
                        placed[m.id] = true;
                        m.modification.formula()
                    })
                })
                .sum::<MolecularFormula>();
        self.aminoacid
            .formulas()
            .iter()
            .map(|formula| formula + &modifications)
            .collect()
    }

    /// Get the molecular formulas for this position with all ambiguous modifications, without any global isotype modifications
    pub fn formulas_all(&self) -> MultiMolecularFormula {
        let modifications = self
            .modifications
            .iter()
            .map(Chemical::formula)
            .sum::<MolecularFormula>()
            + self
                .possible_modifications
                .iter()
                .map(|m| m.modification.formula())
                .sum::<MolecularFormula>();
        self.aminoacid
            .formulas()
            .iter()
            .map(|formula| formula + &modifications)
            .collect()
    }

    /// Get the molecular formulas for this position, with all formulas for the amino acids combined with all options for the modifications.
    /// If you have 2 options for amino acid mass (B or Z) and 2 ambiguous modifications that gives you 8 total options for the mass. (2 AA * 2 amb1 * 2 amb2)
    pub fn formulas_all_options(&self) -> MultiMolecularFormula {
        let modifications = self
            .modifications
            .iter()
            .map(Chemical::formula)
            .sum::<MolecularFormula>();
        let mut formulas: MultiMolecularFormula = self
            .aminoacid
            .formulas()
            .iter()
            .map(|f| f + modifications.clone())
            .collect();
        for modification in &self.possible_modifications {
            formulas *= [
                MolecularFormula::default(),
                modification.modification.formula(),
            ]
            .iter()
            .collect::<MultiMolecularFormula>();
        }
        formulas
    }

    /// Enforce the placement rules of predefined modifications.
    pub(crate) fn enforce_modification_rules(
        &self,
        index: usize,
        length: usize,
    ) -> Result<(), CustomError> {
        for modification in &self.modifications {
            if !modification.is_possible(self, index, length) {
                return Err(CustomError::error(
                    "Modification incorrectly placed",
                    format!(
                        "Modification {modification} is not allowed on aminoacid {} index {index}",
                        self.aminoacid.char()
                    ),
                    Context::none(), // TODO: go and give the correct context here
                ));
            }
        }
        Ok(())
    }
}

impl<T> From<T> for SequenceElement
where
    T: Into<AminoAcid>,
{
    fn from(value: T) -> Self {
        Self::new(value.into(), None)
    }
}
