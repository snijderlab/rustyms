#![warn(dead_code)]

use crate::{
    fragment::{DiagnosticPosition, Fragment, FragmentType, PeptidePosition},
    helper_functions::RangeExtension,
    modification::{
        CrossLinkName, GnoComposition, LinkerSpecificity, Modification, SimpleModification,
    },
    molecular_charge::MolecularCharge,
    peptide::*,
    placement_rule::PlacementRule,
    system::usize::Charge,
    Chemical, DiagnosticIon, Element, Model, MolecularFormula, Multi, MultiChemical, NeutralLoss,
    Protease, SequenceElement,
};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::{
    collections::HashSet,
    fmt::Display,
    marker::PhantomData,
    num::NonZeroU16,
    ops::{Index, RangeBounds},
    slice::SliceIndex,
};
use uom::num_traits::Zero;

/// A peptide with all data as provided by pro forma. Preferably generated by using the [`crate::ComplexPeptide::pro_forma`] function.
#[derive(Default, Clone, PartialOrd, Ord, Debug, Serialize, Deserialize)]
pub struct LinearPeptide<T> {
    /// Global isotope modifications, saved as the element and the species that
    /// all occurrence of that element will consist of. Eg (N, 15) will make
    /// all occurring nitrogens be isotope 15.
    pub(super) global: Vec<(Element, Option<NonZeroU16>)>,
    /// Labile modifications, which will not be found in the actual spectrum.
    pub labile: Vec<SimpleModification>,
    /// N terminal modification
    pub n_term: Option<SimpleModification>,
    /// C terminal modification
    pub c_term: Option<SimpleModification>,
    /// The sequence of this peptide (includes local modifications)
    pub sequence: Vec<SequenceElement>,
    /// For each ambiguous modification list all possible positions it can be placed on.
    /// Indexed by the ambiguous modification id.
    pub ambiguous_modifications: Vec<Vec<usize>>,
    /// The adduct ions, if specified
    pub charge_carriers: Option<MolecularCharge>,
    /// The marker indicating which level of complexity this peptide (potentially) uses
    marker: PhantomData<T>,
}

impl<T> PartialEq for LinearPeptide<T> {
    fn eq(&self, other: &Self) -> bool {
        self.global == other.global
            && self.labile == other.labile
            && self.n_term == other.n_term
            && self.c_term == other.c_term
            && self.sequence == other.sequence
            && self.ambiguous_modifications == other.ambiguous_modifications
            && self.charge_carriers == other.charge_carriers
    }
}

impl<T> std::hash::Hash for LinearPeptide<T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.global.hash(state);
        self.labile.hash(state);
        self.n_term.hash(state);
        self.c_term.hash(state);
        self.sequence.hash(state);
        self.ambiguous_modifications.hash(state);
        self.charge_carriers.hash(state);
    }
}

impl<T> Eq for LinearPeptide<T> {}

/// Builder style methods to create a [`LinearPeptide`]
impl<T> LinearPeptide<T> {
    pub(super) fn mark<M>(self) -> LinearPeptide<M> {
        LinearPeptide {
            global: self.global,
            labile: self.labile,
            n_term: self.n_term,
            c_term: self.c_term,
            sequence: self.sequence,
            ambiguous_modifications: self.ambiguous_modifications,
            charge_carriers: self.charge_carriers,
            marker: PhantomData,
        }
    }

    /// Create a new [`LinearPeptide`], if you want an empty peptide look at [`LinearPeptide::default`].
    /// Potentially the `.collect()` or `.into()` methods can be useful as well.
    #[must_use]
    pub fn new(sequence: impl IntoIterator<Item = SequenceElement>) -> Self {
        sequence.into_iter().collect()
    }

    /// Add global isotope modifications, if any is invalid it returns None
    #[must_use]
    pub fn global(
        mut self,
        global: impl IntoIterator<Item = (Element, Option<NonZeroU16>)>,
    ) -> Option<Self> {
        for modification in global {
            if modification.0.is_valid(modification.1) {
                self.global.push(modification);
            } else {
                return None;
            }
        }
        Some(self)
    }

    /// Add labile modifications
    #[must_use]
    pub fn labile(mut self, labile: impl IntoIterator<Item = SimpleModification>) -> Self {
        self.labile.extend(labile);
        self
    }

    /// Add the N terminal modification
    #[must_use]
    pub fn n_term(mut self, term: Option<SimpleModification>) -> Self {
        self.n_term = term;
        self
    }

    /// Add the C terminal modification
    #[must_use]
    pub fn c_term(mut self, term: Option<SimpleModification>) -> Self {
        self.c_term = term;
        self
    }

    /// Add the charge carriers
    #[must_use]
    pub fn charge_carriers(mut self, charge: Option<MolecularCharge>) -> Self {
        self.charge_carriers = charge;
        self
    }
}

impl<T> LinearPeptide<T> {
    /// Get the number of amino acids making up this peptide
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if there are any amino acids in this peptide
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// The mass of the N terminal modifications. The global isotope modifications are NOT applied.
    pub fn get_n_term(&self) -> MolecularFormula {
        molecular_formula!(H 1)
            + self
                .n_term
                .as_ref()
                .map_or_else(MolecularFormula::default, Chemical::formula)
    }

    /// The mass of the C terminal modifications. The global isotope modifications are NOT applied.
    pub fn get_c_term(&self) -> MolecularFormula {
        molecular_formula!(H 1 O 1)
            + self
                .c_term
                .as_ref()
                .map_or_else(MolecularFormula::default, Chemical::formula)
    }

    /// Get the global isotope modifications
    pub fn get_global(&self) -> &[(Element, Option<NonZeroU16>)] {
        &self.global
    }

    /// Find all neutral losses in the given stretch of peptide
    fn potential_neutral_losses(
        &self,
        range: impl RangeBounds<usize>,
    ) -> Vec<(NeutralLoss, PeptidePosition)> {
        self.iter(range)
            .flat_map(|(pos, aa)| {
                aa.modifications
                    .iter()
                    .filter_map(move |modification| {
                        if let Modification::Simple(SimpleModification::Database {
                            specificities,
                            ..
                        }) = modification
                        {
                            Some(specificities)
                        } else {
                            None
                        }
                    })
                    .flatten()
                    .filter_map(move |(rules, rule_losses, _)| {
                        if PlacementRule::any_possible(rules, aa, &pos) {
                            Some(rule_losses)
                        } else {
                            None
                        }
                    })
                    .flatten()
                    .map(move |loss| (loss.clone(), pos))
            })
            .collect()
    }

    /// Find all diagnostic ions for this full peptide
    fn diagnostic_ions(&self) -> Vec<(DiagnosticIon, DiagnosticPosition)> {
        self.iter(..)
            .flat_map(|(pos, aa)| {
                aa.diagnostic_ions(&pos)
                    .into_iter()
                    .map(move |diag| (diag, DiagnosticPosition::Peptide(pos, aa.aminoacid)))
            })
            .chain(self.labile.iter().flat_map(move |modification| {
                if let SimpleModification::Database { specificities, .. } = modification {
                    specificities
                        .iter()
                        .flat_map(|(_, _, diag)| diag)
                        .map(|diag| {
                            (
                                diag.clone(),
                                DiagnosticPosition::Labile(modification.clone().into()),
                            )
                        })
                        .collect_vec()
                } else if let SimpleModification::Linker { specificities, .. } = modification {
                    specificities
                        .iter()
                        .flat_map(|rule| match rule {
                            LinkerSpecificity::Symmetric(_, _, ions)
                            | LinkerSpecificity::Asymmetric(_, _, ions) => ions,
                        })
                        .map(|diag| {
                            (
                                diag.clone(),
                                DiagnosticPosition::Labile(modification.clone().into()),
                            )
                        })
                        .collect_vec()
                } else {
                    Vec::new()
                }
            }))
            .unique()
            .collect()
    }

    /// Iterate over a range in the peptide and keep track of the position
    pub(super) fn iter(
        &self,
        range: impl RangeBounds<usize>,
    ) -> impl DoubleEndedIterator<Item = (PeptidePosition, &SequenceElement)> + '_ {
        let start = range.start_index();
        self.sequence[(range.start_bound().cloned(), range.end_bound().cloned())]
            .iter()
            .enumerate()
            .map(move |(index, seq)| (PeptidePosition::n(index + start, self.len()), seq))
    }

    /// Iterate over a range in the peptide and keep track of the position
    // pub(super) fn iter_mut(
    //     &mut self,
    //     range: impl RangeBounds<usize>,
    // ) -> impl DoubleEndedIterator<Item = (PeptidePosition, &mut SequenceElement)> + '_ {
    //     let start = match range.start_bound() {
    //         std::ops::Bound::Unbounded => 0,
    //         std::ops::Bound::Included(i) => (*i).max(0),
    //         std::ops::Bound::Excluded(ex) => (ex + 1).max(0),
    //     };
    //     let len = self.len();
    //     self.sequence[(range.start_bound().cloned(), range.end_bound().cloned())]
    //         .iter_mut()
    //         .enumerate()
    //         .map(move |(index, seq)| (PeptidePosition::n(index + start, len), seq))
    // }

    /// Generate all possible patterns for the ambiguous positions (Mass, String:Label).
    /// It always contains at least one pattern (being (base mass, "")).
    /// The global isotope modifications are NOT applied.
    /// Additionally it also returns all peptides present as cross-link.
    fn ambiguous_patterns(
        &self,
        range: impl RangeBounds<usize>,
        aa_range: impl RangeBounds<usize>,
        index: usize,
        base: &MolecularFormula,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        let result = self
            .ambiguous_modifications
            .iter()
            .enumerate()
            .fold(vec![Vec::new()], |acc, (id, possibilities)| {
                acc.into_iter()
                    .flat_map(|path| {
                        let mut path_clone = path.clone();
                        let options = possibilities.iter().filter(|pos| range.contains(pos)).map(
                            move |pos| {
                                let mut new = path.clone();
                                new.push((id, *pos));
                                new
                            },
                        );
                        options.chain(possibilities.iter().find(|pos| !range.contains(pos)).map(
                            move |pos| {
                                path_clone.push((id, *pos));
                                path_clone
                            },
                        ))
                    })
                    .collect()
            })
            .into_iter()
            .fold((Multi::default(), HashSet::new()), |acc, pattern| {
                let ambiguous_local = pattern
                    .iter()
                    .filter_map(|(id, pos)| (*pos == index).then_some(id))
                    .collect::<Vec<_>>();
                let (formulas, seen) = self.sequence[(
                    aa_range.start_bound().cloned(),
                    aa_range.end_bound().cloned(),
                )]
                    .iter()
                    .enumerate()
                    .fold(
                        (Multi::default(), HashSet::new()),
                        |acc, (index, aa)| {
                            let (f, s) = aa.formulas(
                                &pattern
                                    .clone()
                                    .iter()
                                    .copied()
                                    .filter_map(|(id, pos)| (pos == index).then_some(id))
                                    .collect_vec(),
                                all_peptides,
                                visited_peptides,
                                applied_cross_links,
                            );
                            (acc.0 * f, acc.1.union(&s).cloned().collect())
                        },
                    );
                (
                    acc.0
                        * formulas
                            .iter()
                            .map(move |m| {
                                self.sequence[index]
                                    .possible_modifications
                                    .iter()
                                    .filter(|&am| ambiguous_local.contains(&&am.id))
                                    .map(Chemical::formula)
                                    .sum::<MolecularFormula>()
                                    + base
                                    + m
                            })
                            .collect::<Multi<MolecularFormula>>(),
                    acc.1.union(&seen).cloned().collect(),
                )
            });
        if result.0.is_empty() {
            (base.into(), HashSet::new())
        } else {
            result
        }
    }

    /// Generate the theoretical fragments for this peptide, with the given maximal charge of the fragments, and the given model.
    /// With the global isotope modifications applied.
    /// # Panics
    /// Panics if the `max_charge` is bigger than [`isize::MAX`].
    pub(crate) fn generate_theoretical_fragments_inner(
        &self,
        max_charge: Charge,
        model: &Model,
        peptidoform_index: usize,
        peptide_index: usize,
        all_peptides: &[LinearPeptide<Linked>],
    ) -> Vec<Fragment> {
        let default_charge = MolecularCharge::proton(
            isize::try_from(max_charge.value)
                .expect("Charge of the precursor cannot be higher then isize::MAX"),
        );
        let charge_carriers = self.charge_carriers.as_ref().unwrap_or(&default_charge);
        let single_charges = charge_carriers.all_single_charge_options();

        let mut output = Vec::with_capacity(20 * self.sequence.len() + 75); // Empirically derived required size of the buffer (Derived from Hecklib)
        for sequence_index in 0..self.sequence.len() {
            let position = PeptidePosition::n(sequence_index, self.len());
            let mut cross_links = Vec::new();
            let (n_term, n_term_seen) = self.all_masses(
                ..=sequence_index,
                ..sequence_index,
                sequence_index,
                &self.get_n_term(),
                model.modification_specific_neutral_losses,
                all_peptides,
                &[peptide_index],
                &mut cross_links,
            );
            let (c_term, c_term_seen) = self.all_masses(
                sequence_index..,
                sequence_index + 1..,
                sequence_index,
                &self.get_c_term(),
                model.modification_specific_neutral_losses,
                all_peptides,
                &[peptide_index],
                &mut cross_links,
            );
            if !n_term_seen.is_disjoint(&c_term_seen) {
                continue; // There is a link reachable from both sides so there is a loop
            }
            let (modifications_total, modifications_cross_links) = self.sequence[sequence_index]
                .modifications
                .iter()
                .fold((Multi::default(), HashSet::new()), |acc, m| {
                    let (f, s) = m.formula_inner(all_peptides, &[peptide_index], &mut cross_links);
                    (acc.0 * f, acc.1.union(&s).cloned().collect())
                });

            output.append(&mut self.sequence[sequence_index].aminoacid.fragments(
                &n_term,
                &c_term,
                &modifications_total,
                charge_carriers,
                sequence_index,
                self.sequence.len(),
                &model.ions(position),
                peptidoform_index,
                peptide_index,
                (
                    // Allow any N terminal fragment if there is no cross-link to the C terminal side
                    c_term_seen.is_disjoint(&modifications_cross_links),
                    n_term_seen.is_disjoint(&modifications_cross_links),
                ),
            ));

            if model.m {
                // m fragment: precursor amino acid side chain losses
                output.extend(
                    self.formulas_inner(peptide_index, all_peptides, &[], &mut Vec::new())
                        .0
                        .iter()
                        .flat_map(|m| {
                            self.sequence[sequence_index]
                                .aminoacid
                                .formulas()
                                .iter()
                                .flat_map(|aa| {
                                    Fragment::generate_all(
                                        &((-modifications_total.clone()) + m.clone() - aa.clone()
                                            + molecular_formula!(C 2 H 2 N 1 O 1)),
                                        peptidoform_index,
                                        peptide_index,
                                        &FragmentType::m(
                                            position,
                                            self.sequence[sequence_index].aminoacid,
                                        ),
                                        &Multi::default(),
                                        &[],
                                    )
                                })
                                .map(|f| f.with_charge(charge_carriers))
                                .collect_vec()
                        }),
                );
            }
        }
        for fragment in &mut output {
            fragment.formula = fragment
                .formula
                .with_global_isotope_modifications(&self.global)
                .expect("Invalid global isotope modification");
        }

        // Generate precursor peak TODO: generate all masses when cross links break
        let (full_precursor, _all_cross_links) =
            self.formulas_inner(peptide_index, all_peptides, &[], &mut Vec::new());
        output.extend(full_precursor.iter().flat_map(|m| {
            Fragment::new(
                m.clone(),
                Charge::zero(),
                peptidoform_index,
                peptide_index,
                FragmentType::precursor,
                String::new(),
            )
            .with_charge(charge_carriers)
            .with_neutral_losses(&model.precursor)
        }));

        // Generate broken cross-link precursor
        // if model.allow_cross_link_cleavage {
        //     for link in all_cross_links {
        //         let mass = self.formulas_inner(peptide_index, all_peptides, &[], &mut vec![link]).0;

        //     }
        // }

        // Add glycan fragmentation to all peptide fragments
        // Assuming that only one glycan can ever fragment at the same time,
        // and that no peptide fragmentation occurs during glycan fragmentation
        for (sequence_index, position) in self.sequence.iter().enumerate() {
            for modification in &position.modifications {
                if let Modification::Simple(SimpleModification::GlycanStructure(glycan)) =
                    modification
                {
                    output.extend(
                        glycan
                            .clone()
                            .determine_positions()
                            .generate_theoretical_fragments(
                                model,
                                peptidoform_index,
                                peptide_index,
                                charge_carriers,
                                &self
                                    .formulas_inner(
                                        peptide_index,
                                        all_peptides,
                                        &[],
                                        &mut Vec::new(),
                                    )
                                    .0,
                                (position.aminoacid, sequence_index),
                            ),
                    );
                } else if let Modification::Simple(SimpleModification::Gno(
                    GnoComposition::Structure(glycan),
                    _,
                )) = modification
                {
                    output.extend(
                        glycan
                            .clone()
                            .determine_positions()
                            .generate_theoretical_fragments(
                                model,
                                peptidoform_index,
                                peptide_index,
                                charge_carriers,
                                &self
                                    .formulas_inner(
                                        peptide_index,
                                        all_peptides,
                                        &[],
                                        &mut Vec::new(),
                                    )
                                    .0,
                                (position.aminoacid, sequence_index),
                            ),
                    );
                }
            }
        }

        if model.modification_specific_diagnostic_ions {
            // Add all modification diagnostic ions
            output.extend(self.diagnostic_ions().into_iter().flat_map(|(dia, pos)| {
                Fragment {
                    formula: dia.0,
                    charge: Charge::default(),
                    ion: FragmentType::diagnostic(pos),
                    peptidoform_index,
                    peptide_index,
                    neutral_loss: None,
                    label: String::new(),
                    cycles: Vec::new(),
                }
                .with_charges(&single_charges)
            }));
        }

        output
    }

    /// Generate all potential masses for the given stretch of amino acids alongside all peptides seen as part of a cross-link.
    /// Applies ambiguous aminoacids and modifications, and neutral losses (if allowed in the model).
    fn all_masses(
        &self,
        range: impl RangeBounds<usize> + Clone,
        aa_range: impl RangeBounds<usize>,
        index: usize,
        base: &MolecularFormula,
        apply_neutral_losses: bool,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        let (ambiguous_mods_masses, seen) = self.ambiguous_patterns(
            range.clone(),
            aa_range,
            index,
            base,
            all_peptides,
            visited_peptides,
            applied_cross_links,
        );
        if apply_neutral_losses {
            let neutral_losses = self.potential_neutral_losses(range);
            let mut all_masses =
                Vec::with_capacity(ambiguous_mods_masses.len() * (1 + neutral_losses.len()));
            all_masses.extend(ambiguous_mods_masses.iter().cloned());
            for loss in &neutral_losses {
                all_masses.extend((ambiguous_mods_masses.clone() + loss.0.clone()).to_vec());
            }
            (all_masses.into(), seen)
        } else {
            (ambiguous_mods_masses, seen)
        }
    }

    /// Gives all the formulas for the whole peptide with no C and N terminal modifications. With the global isotope modifications applied.
    #[allow(clippy::missing_panics_doc)] // global isotope mods are guaranteed to be correct
    fn bare_formulas_inner(
        &self,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> Multi<MolecularFormula> {
        let mut formulas = Multi::default();
        let mut placed = vec![false; self.ambiguous_modifications.len()];
        for pos in &self.sequence {
            formulas *= pos
                .formulas_greedy(
                    &mut placed,
                    all_peptides,
                    visited_peptides,
                    applied_cross_links,
                )
                .0;
        }

        formulas
            .iter()
            .map(|f| {
                f.with_global_isotope_modifications(&self.global)
                    .expect("Invalid global isotope modification in bare_formulas")
            })
            .collect()
    }

    /// Gives the formulas for the whole peptide. With the global isotope modifications applied. (Any B/Z will result in multiple possible formulas.)
    /// # Panics
    /// When this peptide is already in the set of visited peptides.
    pub(crate) fn formulas_inner(
        &self,
        peptide_index: usize,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        debug_assert!(
            !visited_peptides.contains(&peptide_index),
            "Cannot get the formula for a peptide that is already visited"
        );
        let mut new_visited_peptides = vec![peptide_index];
        new_visited_peptides.extend_from_slice(visited_peptides);
        let mut formulas: Multi<MolecularFormula> =
            vec![self.get_n_term() + self.get_c_term()].into();
        let mut placed = vec![false; self.ambiguous_modifications.len()];
        let mut seen = HashSet::new();
        for pos in &self.sequence {
            let (pos_f, pos_seen) = pos.formulas_greedy(
                &mut placed,
                all_peptides,
                &new_visited_peptides,
                applied_cross_links,
            );
            formulas *= pos_f;
            seen.extend(pos_seen);
        }

        (formulas
            .iter()
            .map(|f| f.with_global_isotope_modifications(&self.global).expect("Global isotope modification invalid in determination of all formulas for a peptide"))
            .collect(), seen)
    }

    /// Display this peptide
    pub(crate) fn display(
        &self,
        f: &mut std::fmt::Formatter<'_>,
        show_global_mods: bool,
    ) -> std::fmt::Result {
        if show_global_mods {
            for (element, isotope) in &self.global {
                write!(
                    f,
                    "<{}{}>",
                    isotope.map(|i| i.to_string()).unwrap_or_default(),
                    element
                )?;
            }
        }
        for labile in &self.labile {
            write!(f, "{{{labile}}}")?;
        }
        // Write any modification of unknown position that has no preferred location at the start of the peptide
        let mut any_ambiguous = false;
        for (id, ambiguous) in self.ambiguous_modifications.iter().enumerate() {
            if ambiguous.iter().all(|i| {
                let m = self.sequence[*i]
                    .possible_modifications
                    .iter()
                    .find(|m| m.id == id)
                    .unwrap();
                !m.preferred
            }) {
                let m = self.sequence[ambiguous[0]]
                    .possible_modifications
                    .iter()
                    .find(|m| m.id == id)
                    .unwrap();
                write!(f, "[{}#{}]", m.modification, m.group)?;
                any_ambiguous = true;
            }
        }
        if any_ambiguous {
            write!(f, "?")?;
        }
        if let Some(m) = &self.n_term {
            write!(f, "[{m}]-")?;
        }
        let mut placed = Vec::new();
        let mut last_ambiguous = None;
        for position in &self.sequence {
            placed.extend(position.display(f, &placed, last_ambiguous)?);
            last_ambiguous = position.ambiguous;
        }
        if last_ambiguous.is_some() {
            // TODO: Does not display ambiguous correctly
            write!(f, ")")?;
        }
        if let Some(m) = &self.c_term {
            write!(f, "-[{m}]")?;
        }
        if let Some(c) = &self.charge_carriers {
            write!(f, "/{c}")?;
        }
        Ok(())
    }
}

impl<T: Clone> LinearPeptide<T> {
    /// Get the reverse of this peptide
    #[must_use]
    pub fn reverse(&self) -> Self {
        Self {
            n_term: self.c_term.clone(),
            c_term: self.n_term.clone(),
            sequence: self.sequence.clone().into_iter().rev().collect(),
            ambiguous_modifications: self
                .ambiguous_modifications
                .clone()
                .into_iter()
                .map(|m| m.into_iter().map(|loc| self.len() - loc).collect())
                .collect(),
            ..self.clone()
        }
    }

    /// Get a region of this peptide as a new peptide (with all terminal/global/ambiguous modifications).
    #[must_use]
    pub fn sub_peptide(&self, index: impl RangeBounds<usize>) -> Self {
        Self {
            n_term: if index.contains(&0) {
                self.n_term.clone()
            } else {
                None
            },
            c_term: if index.contains(&(self.len() - 1)) {
                self.c_term.clone()
            } else {
                None
            },
            sequence: self.sequence[(index.start_bound().cloned(), index.end_bound().cloned())]
                .to_vec(),
            ..self.clone()
        }
    }

    /// Digest this sequence with the given protease and the given maximal number of missed cleavages.
    pub fn digest(&self, protease: &Protease, max_missed_cleavages: usize) -> Vec<Self> {
        let mut sites = vec![0];
        sites.extend_from_slice(&protease.match_locations(&self.sequence));
        sites.push(self.len());

        let mut result = Vec::new();

        for (index, start) in sites.iter().enumerate() {
            for end in sites.iter().skip(index).take(max_missed_cleavages + 1) {
                result.push(self.sub_peptide((*start)..*end));
            }
        }
        result
    }

    /// Concatenate another peptide after this peptide. This will fail if any of these conditions are true:
    /// * The global modifications for the two peptides are not identical
    /// * This peptide has a C terminal modification
    /// * The other peptide has a N terminal modification
    /// * Any of the peptides has ambiguous modifications
    /// * The charge carriers of the two peptides are not identical
    pub fn concatenate(self, other: Self) -> Option<Self> {
        if self.global == other.global
            && self.c_term.is_none()
            && other.n_term.is_none()
            && self.charge_carriers == other.charge_carriers
            && self.ambiguous_modifications.is_empty()
            && other.ambiguous_modifications.is_empty()
        {
            Some(Self {
                global: self.global,
                labile: self.labile.into_iter().chain(other.labile).collect(),
                n_term: self.n_term,
                c_term: other.c_term,
                sequence: self.sequence.into_iter().chain(other.sequence).collect(),
                ambiguous_modifications: Vec::new(),
                charge_carriers: self.charge_carriers,
                marker: self.marker,
            })
        } else {
            None
        }
    }
}

impl LinearPeptide<Linked> {
    /// Gives the formulas for the whole peptide. With the global isotope modifications applied. (Any B/Z will result in multiple possible formulas.)
    pub(crate) fn formulas(
        &self,
        peptide_index: usize,
        all_peptides: &[Self],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> Multi<MolecularFormula> {
        self.formulas_inner(
            peptide_index,
            all_peptides,
            visited_peptides,
            applied_cross_links,
        )
        .0
    }

    /// Gives all the formulas for the whole peptide with no C and N terminal modifications. With the global isotope modifications applied.
    #[allow(clippy::missing_panics_doc, dead_code)] // global isotope mods are guaranteed to be correct
    pub(crate) fn bare_formulas(
        &self,
        all_peptides: &[Self],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> Multi<MolecularFormula> {
        self.bare_formulas_inner(all_peptides, visited_peptides, applied_cross_links)
    }
}

impl<T: Into<Linear>> LinearPeptide<T> {
    /// Gives the formulas for the whole peptide. With the global isotope modifications applied. (Any B/Z will result in multiple possible formulas.)
    #[allow(clippy::missing_panics_doc)] // Can not panic (unless state is already corrupted)
    pub fn formulas(&self) -> Multi<MolecularFormula> {
        let mut formulas: Multi<MolecularFormula> =
            vec![self.get_n_term() + self.get_c_term()].into();
        let mut placed = vec![false; self.ambiguous_modifications.len()];
        for pos in &self.sequence {
            formulas *= pos
                .formulas_greedy(&mut placed, &[], &[], &mut Vec::new())
                .0;
        }

        formulas
            .iter()
            .map(|f| f.with_global_isotope_modifications(&self.global).expect("Global isotope modification invalid in determination of all formulas for a peptide"))
            .collect()
    }

    /// Gives all the formulas for the whole peptide with no C and N terminal modifications. With the global isotope modifications applied.
    pub fn bare_formulas(&self) -> Multi<MolecularFormula> {
        self.bare_formulas_inner(&[], &[], &mut Vec::new())
    }

    /// Generate the theoretical fragments for this peptide, with the given maximal charge of the fragments, and the given model.
    /// With the global isotope modifications applied.
    ///
    /// # Panics
    /// If `max_charge` outside the range `1..=u64::MAX`.
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
    ) -> Vec<Fragment> {
        self.generate_theoretical_fragments_inner(max_charge, model, 0, 0, &[])
    }
}

impl<T> Display for LinearPeptide<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(f, true)
    }
}

impl<Collection, Item, T> From<Collection> for LinearPeptide<T>
where
    Collection: IntoIterator<Item = Item>,
    Item: Into<SequenceElement>,
{
    fn from(value: Collection) -> Self {
        Self {
            global: Vec::new(),
            labile: Vec::new(),
            n_term: None,
            c_term: None,
            sequence: value.into_iter().map(std::convert::Into::into).collect(),
            ambiguous_modifications: Vec::new(),
            charge_carriers: None,
            marker: PhantomData,
        }
    }
}

impl<Item, T> FromIterator<Item> for LinearPeptide<T>
where
    Item: Into<SequenceElement>,
{
    fn from_iter<Iter: IntoIterator<Item = Item>>(iter: Iter) -> Self {
        Self::from(iter)
    }
}

impl<I: SliceIndex<[SequenceElement]>, T> Index<I> for LinearPeptide<T> {
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        &self.sequence[index]
    }
}

/// Make sure that any lower level of peptide can be cast to a higher level
macro_rules! into {
    ($a:tt => $b:ty) => {
        impl From<LinearPeptide<$a>> for LinearPeptide<$b>
        where
            $a: Into<$b>,
        {
            fn from(other: LinearPeptide<$a>) -> Self {
                other.mark()
            }
        }
    };
}

into!(Linear => Linked);
into!(Simple => Linked);
into!(VerySimple => Linked);
into!(ExtremelySimple => Linked);
into!(Simple => Linear);
into!(VerySimple => Linear);
into!(ExtremelySimple => Linear);
into!(VerySimple => Simple);
into!(ExtremelySimple => Simple);
into!(ExtremelySimple => VerySimple);
