#![warn(dead_code)]

use std::collections::HashSet;

use crate::{
    error::{Context, CustomError},
    fragment::PeptidePosition,
    modification::{
        AmbiguousModification, CrossLikeSide, CrossLinkName, LinkerSpecificity, RulePossible,
        SimpleModification,
    },
    peptide::Linked,
    placement_rule::PlacementRule,
    Chemical, DiagnosticIon, LinearPeptide, MolecularFormula, Multi, MultiChemical,
};
use serde::{Deserialize, Serialize};

use crate::{aminoacids::AminoAcid, modification::Modification};

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

    /// # Errors
    /// If the underlying formatter errors.
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
                if m.preferred && !placed.contains(&m.id) {
                    extra_placed.push(m.id);
                    m.modification.to_string()
                } else {
                    String::new()
                },
                m.group,
                m.localisation_score
                    .map(|v| format!("({v})"))
                    .unwrap_or_default()
            )?;
        }
        Ok(extra_placed)
    }

    /// The base constant molecular formula
    fn base_formula(
        &self,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        let (formula, seen) = self
            .modifications
            .iter()
            .map(|m| m.formula_inner(all_peptides, visited_peptides, applied_cross_links))
            .fold((Multi::default(), HashSet::new()), |(am, av), (m, v)| {
                (am * m, av.union(&v).cloned().collect())
            });
        (self.aminoacid.formulas() * formula, seen)
    }

    /// Get the molecular formulas for this position with the selected ambiguous modifications, without any global isotype modifications
    pub fn formulas(
        &self,
        selected_ambiguous: &[usize],
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        let (formula, seen) =
            self.base_formula(all_peptides, visited_peptides, applied_cross_links);
        (
            formula
                + self
                    .possible_modifications
                    .iter()
                    .filter(|&m| selected_ambiguous.contains(&m.id))
                    .map(Chemical::formula)
                    .sum::<MolecularFormula>(),
            seen,
        )
    }

    /// Get the molecular formulas for this position with the ambiguous modifications placed on the very first placed (and updating this in `placed`), without any global isotype modifications
    #[allow(clippy::filter_map_bool_then)] // has side effects
    pub fn formulas_greedy(
        &self,
        placed: &mut [bool],
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        let (formula, seen) =
            self.base_formula(all_peptides, visited_peptides, applied_cross_links);
        (
            formula
                + self
                    .possible_modifications
                    .iter()
                    .filter_map(|m| {
                        (!placed[m.id]).then(|| {
                            placed[m.id] = true;
                            m.formula()
                        })
                    })
                    .sum::<MolecularFormula>(),
            seen,
        )
    }

    /// Get the molecular formulas for this position with all ambiguous modifications, without any global isotype modifications
    pub fn formulas_all(
        &self,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        let (formula, seen) =
            self.base_formula(all_peptides, visited_peptides, applied_cross_links);
        (
            formula
                + self
                    .possible_modifications
                    .iter()
                    .map(Chemical::formula)
                    .sum::<MolecularFormula>(),
            seen,
        )
    }

    /// Get the molecular formulas for this position, with all formulas for the amino acids combined with all options for the modifications.
    /// If you have 2 options for amino acid mass (B or Z) and 2 ambiguous modifications that gives you 8 total options for the mass. (2 AA * 2 amb1 * 2 amb2)
    pub fn formulas_all_options(
        &self,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        let (mut formulas, seen) =
            self.base_formula(all_peptides, visited_peptides, applied_cross_links);
        for modification in &self.possible_modifications {
            formulas = formulas + modification.formula();
        }
        (formulas, seen)
    }

    /// Enforce the placement rules of predefined modifications.
    /// # Errors
    /// If a rule has been broken.
    pub(crate) fn enforce_modification_rules(
        &self,
        position: &PeptidePosition,
    ) -> Result<(), CustomError> {
        for modification in &self.modifications {
            if modification.is_possible(self, position) == RulePossible::No {
                return Err(CustomError::error(
                    "Modification incorrectly placed",
                    format!(
                        "Modification {modification} is not allowed on aminoacid {} index {}",
                        self.aminoacid.char(),
                        position.sequence_index,
                    ),
                    Context::none(),
                ));
            }
        }
        Ok(())
    }

    /// Get all possible diagnostic ions
    pub(crate) fn diagnostic_ions(&self, position: &PeptidePosition) -> Vec<DiagnosticIon> {
        let mut diagnostic_ions = Vec::new();
        for modification in &self.modifications {
            match modification {
                Modification::CrossLink {
                    linker: SimpleModification::Database { specificities, .. },
                    ..
                }
                | Modification::Simple(SimpleModification::Database { specificities, .. }) => {
                    for (rules, _, ions) in specificities {
                        if PlacementRule::any_possible(rules, self, position) {
                            diagnostic_ions.extend_from_slice(ions);
                        }
                    }
                }
                Modification::CrossLink {
                    linker: SimpleModification::Linker { specificities, .. },
                    side,
                    ..
                } => {
                    for rule in specificities {
                        match rule {
                            LinkerSpecificity::Symmetric(rules, _, ions) => {
                                if PlacementRule::any_possible(rules, self, position) {
                                    diagnostic_ions.extend_from_slice(ions);
                                }
                            }
                            LinkerSpecificity::Asymmetric((rules_left, rules_right), _, ions) => {
                                if *side == CrossLikeSide::Symmetric {
                                    if PlacementRule::any_possible(rules_left, self, position)
                                        || PlacementRule::any_possible(rules_right, self, position)
                                    {
                                        diagnostic_ions.extend_from_slice(ions);
                                    }
                                } else if *side == CrossLikeSide::Left {
                                    if PlacementRule::any_possible(rules_left, self, position) {
                                        diagnostic_ions.extend_from_slice(ions);
                                    }
                                } else if *side == CrossLikeSide::Right {
                                    if PlacementRule::any_possible(rules_right, self, position) {
                                        diagnostic_ions.extend_from_slice(ions);
                                    }
                                }
                            }
                        }
                    }
                }
                Modification::Simple(SimpleModification::Linker { specificities, .. }) => {
                    for rule in specificities {
                        match rule {
                            LinkerSpecificity::Symmetric(rules, _, ions) => {
                                if PlacementRule::any_possible(rules, self, position) {
                                    diagnostic_ions.extend_from_slice(ions);
                                }
                            }
                            LinkerSpecificity::Asymmetric((rules_left, rules_right), _, ions) => {
                                if PlacementRule::any_possible(rules_left, self, position)
                                    || PlacementRule::any_possible(rules_right, self, position)
                                {
                                    diagnostic_ions.extend_from_slice(ions);
                                }
                            }
                        }
                    }
                }
                _ => (),
            }
        }
        diagnostic_ions
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
