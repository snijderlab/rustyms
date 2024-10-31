//! Handle modification related issues, access provided if you want to dive deeply into modifications in your own code.

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use std::{
    cmp::Ordering,
    collections::HashSet,
    fmt::{Display, Write},
};

use crate::{
    glycan::{GlycanStructure, MonoSaccharide},
    molecular_charge::CachedCharge,
    peptide::Linked,
    placement_rule::{PlacementRule, Position},
    system::OrderedMass,
    AmbiguousLabel, AminoAcid, Chemical, DiagnosticIon, Fragment, LinearPeptide, Model,
    MolecularFormula, Multi, NeutralLoss, SequenceElement, SequencePosition,
};

include!("shared/modification.rs");

impl ModificationId {
    /// Get the accession number name for the ontology
    pub fn url(&self) -> Option<String> {
        match self.ontology {
            Ontology::Unimod => Some(format!(
                "https://www.unimod.org/modifications_view.php?editid1={}",
                self.id.unwrap_or_default()
            )),
            Ontology::Psimod => Some(format!(
                "https://ontobee.org/ontology/MOD?iri=http://purl.obolibrary.org/obo/MOD_{:05}",
                self.id.unwrap_or_default()
            )),
            Ontology::Gnome => Some(format!(
                "http://glytoucan.org/Structures/Glycans/{}",
                self.name
            )),
            Ontology::Resid => Some(format!(
                "https://proteininformationresource.org/cgi-bin/resid?id=AA{:04}",
                self.id.unwrap_or_default()
            )),
            Ontology::Xlmod | Ontology::Custom => None,
        }
    }
}

/// The result of checking if a modification can be placed somewhere.
#[derive(Debug, PartialEq, Eq, Serialize, Deserialize, Clone)]
pub enum RulePossible {
    /// This modification cannot be placed
    No,
    /// This modification can be placed and if it is a cross-link it can be placed on both ends
    Symmetric(HashSet<usize>),
    /// This modification can be placed and if it is a cross-link it can only be placed on the 'left' side of the cross-link
    AsymmetricLeft(HashSet<usize>),
    /// This modification can be placed and if it is a cross-link it can only be placed on the 'right' side of the cross-link
    AsymmetricRight(HashSet<usize>),
}

impl RulePossible {
    /// Flatten this into a bool, check if the rule is not [`Self::No`]
    pub fn any_possible(self) -> bool {
        self != Self::No
    }
}

impl std::ops::Add for RulePossible {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Symmetric(a), _) | (_, Self::Symmetric(a)) => Self::Symmetric(a),
            (Self::AsymmetricLeft(l), Self::AsymmetricRight(r))
            | (Self::AsymmetricRight(l), Self::AsymmetricLeft(r)) => {
                let overlap: HashSet<usize> = l.intersection(&r).copied().collect();
                if overlap.is_empty() {
                    Self::No
                } else {
                    Self::Symmetric(overlap)
                }
            }
            (Self::AsymmetricLeft(l), _) | (_, Self::AsymmetricLeft(l)) => Self::AsymmetricLeft(l),
            (Self::AsymmetricRight(r), _) | (_, Self::AsymmetricRight(r)) => {
                Self::AsymmetricRight(r)
            }
            _ => Self::No,
        }
    }
}

impl std::iter::Sum for RulePossible {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::No, |acc, i| acc + i)
    }
}

impl Chemical for SimpleModification {
    /// Get the molecular formula for this modification.
    fn formula_inner(&self, position: SequencePosition, peptide_index: usize) -> MolecularFormula {
        match self {
            Self::Mass(m)
            | Self::Gno {
                composition: GnoComposition::Weight(m),
                ..
            } => MolecularFormula::with_additional_mass(m.value),
            Self::Gno {
                composition: GnoComposition::Composition(monosaccharides),
                ..
            }
            | Self::Glycan(monosaccharides) => monosaccharides
                .iter()
                .fold(MolecularFormula::default(), |acc, i| {
                    acc + i.0.formula_inner(position, peptide_index) * i.1 as i32
                }),
            Self::GlycanStructure(glycan)
            | Self::Gno {
                composition: GnoComposition::Topology(glycan),
                ..
            } => glycan.formula_inner(position, peptide_index),
            Self::Formula(formula)
            | Self::Database { formula, .. }
            | Self::Linker { formula, .. } => formula.clone(),
        }
    }
}

impl SimpleModification {
    /// Get a url for more information on this modification. Only defined for modifications from ontologies.
    #[allow(clippy::missing_panics_doc)]
    pub fn ontology_url(&self) -> Option<String> {
        match self {
            Self::Mass(_) | Self::Formula(_) | Self::Glycan(_) | Self::GlycanStructure(_) => None,
            Self::Database { id, .. } | Self::Linker { id, .. } | Self::Gno { id, .. } => id.url(),
        }
    }

    /// Internal formula code with the logic to make all labels right
    pub(crate) fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptide_index: usize,
    ) -> MolecularFormula {
        match self {
            Self::Mass(m)
            | Self::Gno {
                composition: GnoComposition::Weight(m),
                ..
            } => MolecularFormula::with_additional_mass(m.value),
            Self::Gno {
                composition: GnoComposition::Composition(monosaccharides),
                ..
            }
            | Self::Glycan(monosaccharides) => monosaccharides
                .iter()
                .fold(MolecularFormula::default(), |acc, i| {
                    acc + i.0.formula_inner(sequence_index, peptide_index) * i.1 as i32
                }),
            Self::GlycanStructure(glycan)
            | Self::Gno {
                composition: GnoComposition::Topology(glycan),
                ..
            } => glycan.formula_inner(sequence_index, peptide_index),
            Self::Formula(formula)
            | Self::Database { formula, .. }
            | Self::Linker { formula, .. } => formula.clone(),
        }
    }

    /// Check to see if this modification can be placed on the specified element
    pub fn is_possible<T>(
        &self,
        seq: &SequenceElement<T>,
        position: SequencePosition,
    ) -> RulePossible {
        match self {
            Self::Database { specificities, .. } => {
                // If any of the rules match the current situation then it can be placed
                let matching: HashSet<usize> = specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(index, (rules, _, _))| {
                        PlacementRule::any_possible(rules, seq, position).then_some(index)
                    })
                    .collect();
                if matching.is_empty() {
                    RulePossible::No
                } else {
                    RulePossible::Symmetric(matching)
                }
            }
            Self::Linker { specificities, .. } => specificities
                .iter()
                .enumerate()
                .map(|(index, spec)| match spec {
                    LinkerSpecificity::Symmetric(rules, _, _) => {
                        if PlacementRule::any_possible(rules, seq, position) {
                            RulePossible::Symmetric(HashSet::from([index]))
                        } else {
                            RulePossible::No
                        }
                    }
                    LinkerSpecificity::Asymmetric((rules_left, rules_right), _, _) => {
                        let left = PlacementRule::any_possible(rules_left, seq, position);
                        let right = PlacementRule::any_possible(rules_right, seq, position);
                        if left && right {
                            RulePossible::Symmetric(HashSet::from([index]))
                        } else if left {
                            RulePossible::AsymmetricLeft(HashSet::from([index]))
                        } else if right {
                            RulePossible::AsymmetricRight(HashSet::from([index]))
                        } else {
                            RulePossible::No
                        }
                    }
                })
                .sum::<RulePossible>(),
            _ => RulePossible::Symmetric(HashSet::default()),
        }
    }

    /// Check to see if this modification can be placed on the specified element
    pub fn is_possible_aa(&self, aa: AminoAcid, position: Position) -> RulePossible {
        match self {
            Self::Database { specificities, .. } => {
                // If any of the rules match the current situation then it can be placed
                let matching: HashSet<usize> = specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(index, (rules, _, _))| {
                        PlacementRule::any_possible_aa(rules, aa, position).then_some(index)
                    })
                    .collect();
                if matching.is_empty() {
                    RulePossible::No
                } else {
                    RulePossible::Symmetric(matching)
                }
            }
            Self::Linker { specificities, .. } => specificities
                .iter()
                .enumerate()
                .map(|(index, spec)| match spec {
                    LinkerSpecificity::Symmetric(rules, _, _) => {
                        if PlacementRule::any_possible_aa(rules, aa, position) {
                            RulePossible::Symmetric(HashSet::from([index]))
                        } else {
                            RulePossible::No
                        }
                    }
                    LinkerSpecificity::Asymmetric((rules_left, rules_right), _, _) => {
                        let left = PlacementRule::any_possible_aa(rules_left, aa, position);
                        let right = PlacementRule::any_possible_aa(rules_right, aa, position);
                        if left && right {
                            RulePossible::Symmetric(HashSet::from([index]))
                        } else if left {
                            RulePossible::AsymmetricLeft(HashSet::from([index]))
                        } else if right {
                            RulePossible::AsymmetricRight(HashSet::from([index]))
                        } else {
                            RulePossible::No
                        }
                    }
                })
                .sum::<RulePossible>(),
            _ => RulePossible::Symmetric(HashSet::default()),
        }
    }

    /// Display a modification either normalised to the internal representation or as fully valid ProForma
    /// (no glycan structure or custom modifications).
    /// # Errors
    /// When the given writer errors.
    pub fn display(&self, f: &mut impl Write, specification_compliant: bool) -> std::fmt::Result {
        match self {
            Self::Mass(m) => {
                write!(f, "{:+}", m.value)?;
            }
            Self::Formula(elements) => {
                write!(f, "Formula:{}", elements.hill_notation())?;
            }
            Self::Glycan(monosaccharides) => write!(
                f,
                "Glycan:{}",
                monosaccharides
                    .iter()
                    .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
            )?,
            Self::GlycanStructure(glycan) if specification_compliant => write!(
                f,
                "Glycan:{}|INFO:Structure:{glycan}",
                glycan
                    .composition()
                    .iter()
                    .fold(String::new(), |mut acc, (g, a)| {
                        write!(&mut acc, "{g}{a}").unwrap();
                        acc
                    })
            )?,
            Self::GlycanStructure(glycan) => write!(f, "GlycanStructure:{glycan}")?,
            Self::Database {
                formula,
                id:
                    ModificationId {
                        name,
                        ontology: Ontology::Custom,
                        ..
                    },
                ..
            } if specification_compliant => {
                write!(f, "Formula:{formula}|INFO:Custom:{name}")?;
            }
            Self::Database {
                id:
                    ModificationId {
                        name,
                        ontology: Ontology::Custom,
                        ..
                    },
                ..
            } if specification_compliant => {
                write!(f, "C:{name}")?;
            }
            Self::Database { id, .. } | Self::Gno { id, .. } | Self::Linker { id, .. } => {
                write!(f, "{}:{}", id.ontology.char(), id.name)?;
            }
        }
        Ok(())
    }

    /// Get all placement rules as text
    /// # Panics
    /// When a PSI-MOD modification rule uses an non existing modification
    pub(crate) fn placement_rules(&self) -> Vec<String> {
        match self {
            Self::Database { specificities, .. } => specificities
                .iter()
                .flat_map(|set| &set.0)
                .map(|rule| match rule {
                    PlacementRule::AminoAcid(aa, pos) => {
                        format!("{}@{pos}", aa.iter().join(""))
                    }
                    PlacementRule::Terminal(pos) => pos.to_string(),
                    PlacementRule::Anywhere => "Anywhere".to_string(),
                    PlacementRule::PsiModification(index, pos) => {
                        format!(
                            "{}@{pos}",
                            Ontology::Psimod.find_id(*index, None).unwrap_or_else(|| {
                                panic!(
                    "Invalid PsiMod placement rule, non existing modification {index}"
                )
                            })
                        )
                    }
                })
                .collect_vec(),
            Self::Linker { specificities, .. } => specificities
                .iter()
                .flat_map(|set| match set {
                    LinkerSpecificity::Symmetric(rules, _, _) => rules.clone(),
                    LinkerSpecificity::Asymmetric((rules_a, rules_b), _, _) => rules_a
                        .iter()
                        .cloned()
                        .chain(rules_b.iter().cloned())
                        .collect_vec(),
                })
                .map(|rule| match rule {
                    PlacementRule::AminoAcid(aa, pos) => {
                        format!("{}@{pos}", aa.iter().join(""))
                    }
                    PlacementRule::Terminal(pos) => pos.to_string(),
                    PlacementRule::Anywhere => "Anywhere".to_string(),
                    PlacementRule::PsiModification(index, pos) => {
                        format!(
                            "{}@{pos}",
                            Ontology::Psimod.find_id(index, None).unwrap_or_else(|| {
                                panic!(
                "Invalid PsiMod placement rule, non existing modification {index}"
            )
                            })
                        )
                    }
                })
                .collect_vec(),
            _ => Vec::new(),
        }
    }
}

impl Display for SimpleModification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(f, true)
    }
}

impl From<SimpleModification> for Modification {
    fn from(value: SimpleModification) -> Self {
        Self::Simple(value)
    }
}

impl CrossLinkSide {
    /// Get all allowed placement rules with all applicable neutral losses, stubs, and diagnostic ions.
    pub(crate) fn allowed_rules(
        &self,
        linker: &SimpleModification,
    ) -> (
        Vec<NeutralLoss>,
        Vec<(MolecularFormula, MolecularFormula)>,
        Vec<DiagnosticIon>,
    ) {
        let selected_rules = match self {
            Self::Left(r) | Self::Right(r) | Self::Symmetric(r) => r,
        };
        let mut neutral = Vec::new();
        let mut stubs = Vec::new();
        let mut diagnostic = Vec::new();

        match linker {
            SimpleModification::Linker { specificities, .. } => {
                for rule in specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(i, r)| selected_rules.contains(&i).then_some(r))
                {
                    match rule {
                        LinkerSpecificity::Asymmetric(_, n, d) => {
                            diagnostic.extend_from_slice(d);
                            match self {
                                Self::Left(_) => stubs.extend(n.iter().cloned()),
                                Self::Right(_) => {
                                    stubs.extend(n.iter().map(|(l, r)| (r.clone(), l.clone())));
                                }
                                Self::Symmetric(_) => stubs.extend(n.iter().flat_map(|(l, r)| {
                                    vec![(l.clone(), r.clone()), (r.clone(), l.clone())]
                                })),
                            }
                        }
                        LinkerSpecificity::Symmetric(_, n, d) => {
                            stubs.extend_from_slice(n);
                            diagnostic.extend_from_slice(d);
                        }
                    }
                }
            }
            SimpleModification::Database { specificities, .. } => {
                for rule in specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(i, r)| selected_rules.contains(&i).then_some(r))
                {
                    neutral.extend_from_slice(&rule.1);
                    diagnostic.extend_from_slice(&rule.2);
                }
            }
            _ => (),
        };
        (neutral, stubs, diagnostic)
    }
}

impl Modification {
    /// Get the formula for the whole addition (or subtraction) for this modification
    pub(crate) fn formula_inner(
        &self,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
        allow_ms_cleavable: bool,
        sequence_index: SequencePosition,
        peptide_index: usize,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        match self {
            // A linker that is not cross-linked is hydrolysed
            Self::Simple(SimpleModification::Linker { formula, .. }) => (
                (formula.clone() + molecular_formula!(H 2 O 1)).into(),
                HashSet::new(),
            ),
            Self::Simple(s) => (
                s.formula_inner(sequence_index, peptide_index).into(),
                HashSet::new(),
            ),
            Self::CrossLink {
                peptide: other_peptide,
                linker,
                name,
                side,
                ..
            } => {
                if applied_cross_links.contains(name) {
                    (Multi::default(), HashSet::default())
                } else if visited_peptides.contains(other_peptide) {
                    applied_cross_links.push(name.clone());
                    (
                        linker
                            .formula_inner(sequence_index, peptide_index)
                            .with_label(AmbiguousLabel::CrossLinkBound(name.clone()))
                            .into(),
                        HashSet::from([name.clone()]),
                    )
                } else {
                    applied_cross_links.push(name.clone());
                    let link = linker.formula_inner(sequence_index, peptide_index);
                    let (_, stubs, _) = side.allowed_rules(linker);

                    if allow_ms_cleavable && !stubs.is_empty() {
                        let mut options: Vec<MolecularFormula> = stubs
                            .iter()
                            .map(|s| {
                                s.0.clone().with_label(AmbiguousLabel::CrossLinkBroken(
                                    name.clone(),
                                    s.0.clone(),
                                ))
                            })
                            .unique()
                            .collect();
                        let mut seen_peptides = HashSet::from([name.clone()]);

                        options.extend_from_slice(&{
                            let (f, seen) = all_peptides[*other_peptide].formulas_inner(
                                *other_peptide,
                                all_peptides,
                                visited_peptides,
                                applied_cross_links,
                                false,
                            );
                            seen_peptides.extend(seen);
                            (f + link)
                                .with_label(AmbiguousLabel::CrossLinkBound(name.clone()))
                                .to_vec()
                        });

                        (options.into(), seen_peptides)
                    } else {
                        let (f, mut seen) = all_peptides[*other_peptide].formulas_inner(
                            *other_peptide,
                            all_peptides,
                            visited_peptides,
                            applied_cross_links,
                            false,
                        );
                        seen.insert(name.clone());
                        (
                            (f + link).with_label(AmbiguousLabel::CrossLinkBound(name.clone())),
                            seen,
                        )
                    }
                }
            }
        }
    }

    /// Get the formula for a modification, if it is a cross linked modification only get the cross link
    pub fn formula(&self) -> MolecularFormula {
        match self {
            Self::Simple(s) => s.formula(),
            Self::CrossLink { linker, .. } => linker.formula(),
        }
    }
}

impl Modification {
    /// Check if this is a simple modification
    pub const fn simple(&self) -> Option<&SimpleModification> {
        match self {
            Self::Simple(sim) => Some(sim),
            Self::CrossLink { .. } => None,
        }
    }

    /// Check if this is a simple modification
    pub fn into_simple(self) -> Option<SimpleModification> {
        match self {
            Self::Simple(sim) => Some(sim),
            Self::CrossLink { .. } => None,
        }
    }

    /// Get a url for more information on this modification. Only defined for modifications from ontologies.
    #[allow(clippy::missing_panics_doc)]
    pub fn ontology_url(&self) -> Option<String> {
        match self {
            Self::Simple(modification) => modification.ontology_url(),
            Self::CrossLink { linker, .. } => linker.ontology_url(),
        }
    }

    /// Check to see if this modification can be placed on the specified element
    pub fn is_possible<T>(
        &self,
        seq: &SequenceElement<T>,
        position: SequencePosition,
    ) -> RulePossible {
        self.simple()
            .map_or(RulePossible::Symmetric(HashSet::new()), |s| {
                s.is_possible(seq, position)
            })
    }

    /// Generate theoretical fragments for side chains (glycans)
    pub(crate) fn generate_theoretical_fragments(
        &self,
        model: &Model,
        peptidoform_index: usize,
        peptide_index: usize,
        charge_carriers: &mut CachedCharge,
        full_formula: &Multi<MolecularFormula>,
        attachment: Option<(AminoAcid, usize)>,
    ) -> Vec<Fragment> {
        if let Self::Simple(simple) = self {
            simple.generate_theoretical_fragments(
                model,
                peptidoform_index,
                peptide_index,
                charge_carriers,
                full_formula,
                attachment,
            )
        } else {
            Vec::new()
        }
    }
}

impl SimpleModification {
    /// Generate theoretical fragments for side chains (glycans)
    pub(crate) fn generate_theoretical_fragments(
        &self,
        model: &Model,
        peptidoform_index: usize,
        peptide_index: usize,
        charge_carriers: &mut CachedCharge,
        full_formula: &Multi<MolecularFormula>,
        attachment: Option<(AminoAcid, usize)>,
    ) -> Vec<Fragment> {
        match self {
            Self::GlycanStructure(glycan)
            | Self::Gno {
                composition: GnoComposition::Topology(glycan),
                ..
            } => glycan
                .clone()
                .determine_positions()
                .generate_theoretical_fragments(
                    model,
                    peptidoform_index,
                    peptide_index,
                    charge_carriers,
                    full_formula,
                    attachment,
                ),
            Self::Glycan(composition)
            | Self::Gno {
                composition: GnoComposition::Composition(composition),
                ..
            } => MonoSaccharide::theoretical_fragments(
                composition,
                model,
                peptidoform_index,
                peptide_index,
                charge_carriers,
                full_formula,
                attachment,
            ),
            _ => Vec::new(),
        }
    }
}

/// The structure to lookup ambiguous modifications, with a list of all modifications (the order is fixed) with for each modification their name and the actual modification itself (if already defined)
pub type AmbiguousLookup = Vec<(String, Option<SimpleModification>)>;
/// The structure to lookup cross-links, with a list of all linkers (the order is fixed) with for each linker their name or None if it is a branch and the actual linker itself (if already defined)
pub type CrossLinkLookup = Vec<(CrossLinkName, Option<SimpleModification>)>;

/// An ambiguous modification which could be placed on any of a set of locations
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub struct AmbiguousModification {
    /// The id to compare be able to find the other locations where this modifications can be placed
    pub id: usize,
    /// The modification itself
    pub modification: SimpleModification,
    /// If present the localisation score, meaning the chance/ratio for this modification to show up on this exact spot
    pub localisation_score: Option<OrderedFloat<f64>>,
    /// The name of the group
    pub group: String,
    /// If this is the preferred location or not
    pub preferred: bool,
}

impl Chemical for AmbiguousModification {
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptide_index: usize,
    ) -> MolecularFormula {
        self.modification
            .formula_inner(sequence_index, peptide_index)
    }
}

impl Modification {
    /// Display a modification either normalised to the internal representation or as fully valid ProForma
    /// (no glycan structure or custom modifications).
    /// # Errors
    /// When the given writer errors.
    pub fn display(&self, f: &mut impl Write, specification_compliant: bool) -> std::fmt::Result {
        match self {
            Self::Simple(sim) => sim.display(f, specification_compliant),
            Self::CrossLink { name, linker, .. } => write!(f, "{linker}{name}"),
        }
    }
}

impl Display for Modification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(f, true)
    }
}

impl Display for CrossLinkName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Branch => write!(f, "#BRANCH"),
            Self::Name(n) => write!(f, "#XL{n}"),
        }
    }
}

include!("shared/ontology.rs");

#[test]
#[allow(clippy::missing_panics_doc)]
fn test_reading_custom_modifications_json() {
    use serde_json;
    let data = r#"[ [ 1, "uranium linker", { "Linker": { "specificities": [ { "Asymmetric": [ [ [ { "AminoAcid": [ [ "Selenocysteine" ], "AnyCTerm" ] }, { "AminoAcid": [ [ "GlutamicAcid" ], "Anywhere" ] } ], [ { "AminoAcid": [ [ "Selenocysteine" ], "AnyNTerm" ] } ] ], [ [ { "elements": [ [ "U", null, 1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "U", null, 1 ] ], "additional_mass": 0.0 } ] ], [ { "elements": [ [ "Te", null, 1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "Ne", null, 1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "H", null, 2 ], [ "He", null, 3 ] ], "additional_mass": 0.0 }, { "elements": [ [ "H", null, 1 ], [ "He", null, 2 ] ], "additional_mass": 0.0 }, { "elements": [ [ "I", null, 1 ], [ "Er", null, 1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "H", null, 12 ], [ "C", null, 12 ], [ "O", null, 1 ] ], "additional_mass": 0.0 } ] ] } ], "formula": { "elements": [ [ "U", null, 2 ] ], "additional_mass": 0.0 }, "id": { "ontology": "Custom", "name": "Uranium linker", "id": 1, "description": "Have some uranium, its delicious!", "synonyms": [], "cross_ids": [ [ "Pubmed", "21714143" ] ] }, "length": 23.9 } } ], [ 2, "helium", { "Database": { "specificities": [ [ [ { "AminoAcid": [ [ "Alanine" ], "Anywhere" ] } ], [], [] ], [ [ { "AminoAcid": [ [ "Methionine" ], "Anywhere" ] } ], [ { "Loss": { "elements": [], "additional_mass": 12.0 } } ], [] ] ], "formula": { "elements": [ [ "He", null, 2 ] ], "additional_mass": 0.0 }, "id": { "ontology": "Custom", "name": "Helium", "id": 2, "description": "heeeelium", "synonyms": [ "heeeelium", "funny gas" ], "cross_ids": [ [ "Pubmed", "42" ], [ "Unimod", "12" ], [ "Resid", "A12" ] ] } } } ], [ 3, "db18", { "Database": { "specificities": [ [ [ { "AminoAcid": [ [ "Cysteine" ], "Anywhere" ] } ], [ { "Gain": { "elements": [], "additional_mass": 372.25 } }, { "Gain": { "elements": [], "additional_mass": 373.258 } }, { "Gain": { "elements": [], "additional_mass": 371.242 } }, { "Gain": { "elements": [], "additional_mass": 240.171 } }, { "Gain": { "elements": [], "additional_mass": 239.163 } }, { "Gain": { "elements": [], "additional_mass": 241.179 } }, { "Gain": { "elements": [], "additional_mass": 197.129 } }, { "Gain": { "elements": [], "additional_mass": 198.137 } }, { "Gain": { "elements": [], "additional_mass": 196.121 } }, { "Gain": { "elements": [], "additional_mass": 619.418 } }, { "Gain": { "elements": [], "additional_mass": 621.4343 } }, { "Gain": { "elements": [], "additional_mass": 649.465 } }, { "Gain": { "elements": [], "additional_mass": 677.497 } }, { "Gain": { "elements": [], "additional_mass": 618.41 } }, { "Gain": { "elements": [], "additional_mass": 620.426 } }, { "Gain": { "elements": [], "additional_mass": 648.457 } }, { "Gain": { "elements": [], "additional_mass": 676.489 } }, { "Gain": { "elements": [], "additional_mass": 620.426 } }, { "Gain": { "elements": [], "additional_mass": 622.442 } }, { "Gain": { "elements": [], "additional_mass": 650.473 } }, { "Gain": { "elements": [], "additional_mass": 678.504 } }, { "Gain": { "elements": [], "additional_mass": 28.006 } } ], [ { "elements": [], "additional_mass": 372.25 }, { "elements": [], "additional_mass": 240.171 }, { "elements": [], "additional_mass": 197.129 }, { "elements": [], "additional_mass": 619.418 }, { "elements": [], "additional_mass": 621.434 }, { "elements": [], "additional_mass": 649.465 }, { "elements": [], "additional_mass": 677.497 } ] ] ], "formula": { "elements": [], "additional_mass": 676.489 }, "id": { "ontology": "Custom", "name": "DB18", "id": 3, "description": "", "synonyms": [], "cross_ids": [] } } } ], [ 4, "disulfide", { "Linker": { "specificities": [ { "Symmetric": [ [ { "AminoAcid": [ [ "Cysteine" ], "Anywhere" ] } ], [ [ { "elements": [], "additional_mass": 0.0 }, { "elements": [], "additional_mass": 0.0 } ], [ { "elements": [ [ "H", null, -1 ] ], "additional_mass": 0.0 }, { "elements": [], "additional_mass": 0.0 } ], [ { "elements": [ [ "H", null, -1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "H", null, -1 ] ], "additional_mass": 0.0 } ] ], [] ] } ], "formula": { "elements": [ [ "H", null, -2 ] ], "additional_mass": 0.0 }, "id": { "ontology": "Custom", "name": "Disulfide", "id": 4, "description": "A disulfide with all potential neutral losses", "synonyms": [], "cross_ids": [] }, "length": 0.0 } } ], [ 5, "dsso", { "Linker": { "specificities": [ { "Symmetric": [ [ { "AminoAcid": [ [ "Lysine" ], "Anywhere" ] } ], [ [ { "elements": [ [ "H", null, 1 ], [ "C", null, 3 ], [ "N", null, -1 ], [ "O", null, 3 ], [ "S", null, 1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "H", null, 1 ], [ "C", null, 3 ], [ "N", null, -1 ], [ "O", null, 2 ] ], "additional_mass": 0.0 } ] ], [] ] } ], "formula": { "elements": [ [ "H", null, 2 ], [ "C", null, 6 ], [ "N", null, -2 ], [ "O", null, 5 ], [ "S", null, 1 ] ], "additional_mass": 0.0 }, "id": { "ontology": "Custom", "name": "DSSO", "id": 5, "description": "", "synonyms": [], "cross_ids": [] }, "length": 0.0 } } ]]"#;
    let mods: Vec<(usize, String, SimpleModification)> = serde_json::from_str(data).unwrap();
    assert!(mods.len() > 1);
}
