#![warn(dead_code)]

use std::{collections::HashSet, fmt::Write, marker::PhantomData, num::NonZeroU32};

use crate::{
    error::{Context, CustomError},
    modification::{
        AmbiguousModification, CrossLinkName, LinkerSpecificity, Modification, RulePossible,
        SimpleModification, SimpleModificationInner,
    },
    peptide::{AtLeast, Linked},
    placement_rule::PlacementRule,
    CheckedAminoAcid, Chemical, DiagnosticIon, LinearPeptide, MolecularFormula, Multi,
    MultiChemical, SequencePosition,
};
use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

/// One block in a sequence meaning an aminoacid and its accompanying modifications
#[derive(Default, PartialOrd, Ord, Debug, Serialize, Deserialize)]
pub struct SequenceElement<T> {
    /// The aminoacid
    pub aminoacid: CheckedAminoAcid<T>,
    /// All present modifications
    pub modifications: ThinVec<Modification>,
    /// All ambiguous modifications (could be placed here or on another position)
    pub possible_modifications: ThinVec<AmbiguousModification>,
    /// If this aminoacid is part of an ambiguous sequence group `(QA)?` in ProForma
    pub ambiguous: Option<NonZeroU32>,
    /// The marker indicating which level of complexity this sequence element uses as higher bound
    marker: PhantomData<T>,
}

impl<T> Clone for SequenceElement<T> {
    fn clone(&self) -> Self {
        Self {
            aminoacid: self.aminoacid,
            modifications: self.modifications.clone(),
            possible_modifications: self.possible_modifications.clone(),
            ambiguous: self.ambiguous,
            marker: PhantomData,
        }
    }
}

impl<A, B> PartialEq<SequenceElement<B>> for SequenceElement<A> {
    fn eq(&self, other: &SequenceElement<B>) -> bool {
        self.aminoacid == other.aminoacid
            && self.modifications == other.modifications
            && self.possible_modifications == other.possible_modifications
            && self.ambiguous == other.ambiguous
    }
}

impl<T> std::hash::Hash for SequenceElement<T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.aminoacid.hash(state);
        self.modifications.hash(state);
        self.possible_modifications.hash(state);
        self.ambiguous.hash(state);
    }
}

impl<T> Eq for SequenceElement<T> {}

impl<T> SequenceElement<T> {
    /// Mark this sequence element as the following complexity level, the level is not validated
    pub(super) fn mark<M>(self) -> SequenceElement<M> {
        SequenceElement {
            aminoacid: self.aminoacid.mark::<M>(),
            modifications: self.modifications,
            possible_modifications: self.possible_modifications,
            ambiguous: self.ambiguous,
            marker: PhantomData,
        }
    }

    /// Cast a sequence element into a more complex sequence element. This does not change the
    /// content of the sequence element. It only allows to pass this as higher complexity if needed.
    pub fn cast<OtherComplexity: AtLeast<T>>(self) -> SequenceElement<OtherComplexity> {
        self.mark()
    }

    /// Create a new aminoacid without any modifications
    pub fn new(aminoacid: CheckedAminoAcid<T>, ambiguous: Option<NonZeroU32>) -> Self {
        Self {
            aminoacid,
            modifications: ThinVec::new(),
            possible_modifications: ThinVec::new(),
            ambiguous,
            marker: PhantomData,
        }
    }

    /// Add a modification to this sequence element
    #[must_use]
    pub fn with_simple_modification(mut self, modification: SimpleModification) -> Self {
        self.modifications.push(Modification::Simple(modification));
        self
    }

    /// Add a modification to this sequence element
    pub fn add_simple_modification(&mut self, modification: SimpleModification) {
        self.modifications.push(Modification::Simple(modification));
    }
}

impl<T> SequenceElement<T> {
    /// # Errors
    /// If the underlying formatter errors.
    pub(crate) fn display(
        &self,
        f: &mut impl Write,
        placed: &[usize],
        last_ambiguous: Option<NonZeroU32>,
        specification_compliant: bool,
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
            write!(f, "[")?;
            m.display(f, specification_compliant)?;
            write!(f, "]")?;
        }
        for m in &self.possible_modifications {
            write!(f, "[",)?;
            if m.preferred && !placed.contains(&m.id) {
                extra_placed.push(m.id);
                m.modification.display(f, specification_compliant)?;
            };
            write!(
                f,
                "\x23{}{}]",
                m.group,
                m.localisation_score
                    .map(|v| format!("({v})"))
                    .unwrap_or_default()
            )?;
        }
        Ok(extra_placed)
    }

    /// The base constant molecular formula
    pub(crate) fn base_formula(
        &self,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
        allow_ms_cleavable: bool,
        sequence_index: SequencePosition,
        peptide_index: usize,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        let (formula, seen) = self
            .modifications
            .iter()
            .map(|m| {
                m.formula_inner(
                    all_peptides,
                    visited_peptides,
                    applied_cross_links,
                    allow_ms_cleavable,
                    sequence_index,
                    peptide_index,
                )
            })
            .fold((Multi::default(), HashSet::new()), |(am, av), (m, v)| {
                (am * m, av.union(&v).cloned().collect())
            });
        (
            self.aminoacid.formulas_inner(sequence_index, peptide_index) * formula,
            seen,
        )
    }

    /// Get the molecular formulas for this position with the ambiguous modifications placed on the very first placed (and updating this in `placed`), without any global isotype modifications
    #[allow(clippy::filter_map_bool_then, clippy::too_many_arguments)] // has side effects
    pub(crate) fn formulas_greedy(
        &self,
        placed: &mut [bool],
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
        allow_ms_cleavable: bool,
        sequence_index: SequencePosition,
        peptide_index: usize,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        let (formula, seen) = self.base_formula(
            all_peptides,
            visited_peptides,
            applied_cross_links,
            allow_ms_cleavable,
            sequence_index,
            peptide_index,
        );
        (
            formula
                + self
                    .possible_modifications
                    .iter()
                    .filter_map(|m| {
                        (!placed[m.id]).then(|| {
                            placed[m.id] = true;
                            m.formula_inner(sequence_index, peptide_index)
                        })
                    })
                    .sum::<MolecularFormula>(),
            seen,
        )
    }

    /// Get the molecular formulas for this position with all ambiguous modifications, without any global isotype modifications
    pub(crate) fn formulas_all(
        &self,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
        allow_ms_cleavable: bool,
        sequence_index: SequencePosition,
        peptide_index: usize,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        let (formula, seen) = self.base_formula(
            all_peptides,
            visited_peptides,
            applied_cross_links,
            allow_ms_cleavable,
            sequence_index,
            peptide_index,
        );
        (
            formula
                + self
                    .possible_modifications
                    .iter()
                    .map(|f| f.formula_inner(sequence_index, peptide_index))
                    .sum::<MolecularFormula>(),
            seen,
        )
    }

    /// Enforce the placement rules of predefined modifications.
    /// # Errors
    /// If a rule has been broken.
    /// # Panics
    /// If any placement rule is placement on a PSI modification that does not exist.
    pub(crate) fn enforce_modification_rules(
        &self,
        position: SequencePosition,
    ) -> Result<(), CustomError> {
        for modification in &self.modifications {
            if modification.is_possible(self, position) == RulePossible::No {
                let rules = modification
                    .simple()
                    .map(|s| s.placement_rules())
                    .unwrap_or_default();
                return Err(CustomError::error(
                    "Modification incorrectly placed",
                    format!(
                        "Modification {modification} is not allowed on {}{}",
                        match position {
                            SequencePosition::NTerm => "the N-terminus".to_string(),
                            SequencePosition::CTerm => "the C-terminus".to_string(),
                            SequencePosition::Index(index) =>
                                format!("the side chain of {} at index {index}", self.aminoacid),
                        },
                        if rules.is_empty() {
                            String::new()
                        } else {
                            format!(", this modification is only allowed at the following locations: {}", rules.join(", "))
                        }
                    ),
                    Context::none(),
                ));
            }
        }
        Ok(())
    }

    /// Get all possible diagnostic ions
    pub(crate) fn diagnostic_ions(&self, position: SequencePosition) -> Vec<DiagnosticIon> {
        let mut diagnostic_ions = Vec::new();
        for modification in &self.modifications {
            match modification {
                Modification::CrossLink { linker, side, .. } => {
                    diagnostic_ions.extend_from_slice(&side.allowed_rules(linker).2);
                }
                Modification::Simple(sim) => match &**sim {
                    SimpleModificationInner::Database { specificities, .. } => {
                        for (rules, _, ions) in specificities {
                            if PlacementRule::any_possible(rules, self, position) {
                                diagnostic_ions.extend_from_slice(ions);
                            }
                        }
                    }
                    SimpleModificationInner::Linker { specificities, .. } => {
                        for rule in specificities {
                            match rule {
                                LinkerSpecificity::Symmetric(rules, _, ions) => {
                                    if PlacementRule::any_possible(rules, self, position) {
                                        diagnostic_ions.extend_from_slice(ions);
                                    }
                                }
                                LinkerSpecificity::Asymmetric(
                                    (rules_left, rules_right),
                                    _,
                                    ions,
                                ) => {
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
                },
            }
        }
        diagnostic_ions
    }
}

impl<C, T> From<T> for SequenceElement<C>
where
    T: Into<CheckedAminoAcid<C>>,
{
    fn from(value: T) -> Self {
        Self::new(value.into(), None)
    }
}
