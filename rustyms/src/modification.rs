//! Handle modification related issues, access provided if you want to dive deeply into modifications in your own code.

use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use std::{
    collections::HashSet,
    fmt::{Display, Write},
};

use crate::{
    fragment::PeptidePosition,
    glycan::{GlycanStructure, MonoSaccharide},
    ontologies::CustomDatabase,
    peptide::Linked,
    placement_rule::PlacementRule,
    system::{Mass, OrderedMass},
    Chemical, DiagnosticIon, LinearPeptide, MolecularFormula, Multi, NeutralLoss, SequenceElement,
    Tolerance, WithinTolerance,
};

include!("shared/modification.rs");

impl ModificationId {
    /// Get the accession number name for the ontology
    pub fn url(&self) -> Option<String> {
        match self.ontology {
            Ontology::Unimod => Some(format!(
                "https://www.unimod.org/modifications_view.php?editid1={}",
                self.id
            )),
            Ontology::Psimod => Some(format!(
                "https://ontobee.org/ontology/MOD?iri=http://purl.obolibrary.org/obo/MOD_{:05}",
                self.id
            )),
            Ontology::Gnome => Some(format!(
                "https://gnome.glyomics.org/StructureBrowser.html?focus={}",
                self.name
            )),
            Ontology::Resid => Some(format!(
                "https://proteininformationresource.org/cgi-bin/resid?id=AA{:04}",
                self.id
            )),
            Ontology::Xlmod | Ontology::Custom => None,
        }
    }
}

impl Chemical for SimpleModification {
    fn formula(&self) -> MolecularFormula {
        match self {
            Self::Mass(m) => MolecularFormula::with_additional_mass(m.value),
            Self::Formula(elements) => elements.clone(),
            Self::Glycan(monosaccharides) => monosaccharides
                .iter()
                .fold(MolecularFormula::default(), |acc, i| {
                    acc + i.0.formula() * i.1 as i32
                }),
            Self::GlycanStructure(glycan) | Self::Gno(GnoComposition::Structure(glycan), _) => {
                glycan.formula()
            }
            Self::Database { formula, .. } | Self::Linker { formula, .. } => formula.clone(),
            Self::Gno(GnoComposition::Mass(m), _) => {
                MolecularFormula::with_additional_mass(m.value)
            }
        }
    }
}

/// The result of checking if a modification can be placed somewhere.
#[derive(Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize, Clone, Copy)]
pub enum RulePossible {
    /// This modification cannot be placed
    No,
    /// This modification can be placed and if it is a cross-link it can be placed on both ends
    Symmetric,
    /// This modification can be placed and if it is a cross-link it can only be placed on the 'left' side of the cross-link
    AsymmetricLeft,
    /// This modification can be placed and if it is a cross-link it can only be placed on the 'right' side of the cross-link
    AsymmetricRight,
}

impl RulePossible {
    /// Flatten this into a bool, check if the rule is not [`Self::No`]
    pub fn possible(self) -> bool {
        self != Self::No
    }
}

impl std::ops::Add for RulePossible {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Symmetric, _)
            | (_, Self::Symmetric)
            | (Self::AsymmetricLeft, Self::AsymmetricRight)
            | (Self::AsymmetricRight, Self::AsymmetricLeft) => Self::Symmetric,
            (Self::AsymmetricLeft, _) | (_, Self::AsymmetricLeft) => Self::AsymmetricLeft,
            (Self::AsymmetricRight, _) | (_, Self::AsymmetricRight) => Self::AsymmetricRight,
            _ => Self::No,
        }
    }
}

impl std::iter::Sum for RulePossible {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::No, |acc, i| acc + i)
    }
}

impl SimpleModification {
    /// Get a url for more information on this modification. Only defined for modifications from ontologies.
    #[allow(clippy::missing_panics_doc)]
    pub fn ontology_url(&self) -> Option<String> {
        match self {
            Self::Mass(_) | Self::Formula(_) | Self::Glycan(_) | Self::GlycanStructure(_) => None,
            Self::Database { id, .. } | Self::Linker { id, .. } => id.url(),
            Self::Gno(_, name) => Some(format!(
                "https://gnome.glyomics.org/StructureBrowser.html?focus={name}",
            )),
        }
    }

    /// Check to see if this modification can be placed on the specified element
    pub fn is_possible(&self, seq: &SequenceElement, position: &PeptidePosition) -> RulePossible {
        match self {
            Self::Database { specificities, .. } => {
                // If any of the rules match the current situation then it can be placed
                if specificities
                    .iter()
                    .any(|(rules, _, _)| PlacementRule::any_possible(rules, seq, position))
                {
                    RulePossible::Symmetric
                } else {
                    RulePossible::No
                }
            }
            Self::Linker { specificities, .. } => specificities
                .iter()
                .map(|spec| match spec {
                    LinkerSpecificity::Symmetric(rules, _, _) => {
                        if PlacementRule::any_possible(rules, seq, position) {
                            RulePossible::Symmetric
                        } else {
                            RulePossible::No
                        }
                    }
                    LinkerSpecificity::Asymmetric((rules_left, rules_right), _, _) => {
                        let left = PlacementRule::any_possible(rules_left, seq, position);
                        let right = PlacementRule::any_possible(rules_right, seq, position);
                        if left && right {
                            RulePossible::Symmetric
                        } else if left {
                            RulePossible::AsymmetricLeft
                        } else if right {
                            RulePossible::AsymmetricRight
                        } else {
                            RulePossible::No
                        }
                    }
                })
                .sum::<RulePossible>(),
            _ => RulePossible::Symmetric,
        }
    }
}

impl Display for SimpleModification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
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
            Self::GlycanStructure(glycan) => write!(
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
            Self::Database {
                formula,
                id:
                    ModificationId {
                        name,
                        ontology: Ontology::Custom,
                        ..
                    },
                ..
            } => {
                write!(f, "Formula:{formula}|INFO:Custom:{name}")?;
            }
            Self::Database { id, .. } => {
                write!(f, "{}:{}", id.ontology.char(), id.name)?;
            }
            Self::Gno(_, name) => write!(f, "{}:{name}", Ontology::Gnome.char())?,
            Self::Linker { id, .. } => write!(f, "{}:{}", id.ontology.char(), id.name)?,
        }
        Ok(())
    }
}

impl From<SimpleModification> for Modification {
    fn from(value: SimpleModification) -> Self {
        Self::Simple(value)
    }
}

impl Modification {
    /// Get the formula for the whole addition (or subtraction) for this modification
    pub(crate) fn formula_inner(
        &self,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
    ) -> (Multi<MolecularFormula>, HashSet<CrossLinkName>) {
        match self {
            // A linker that is not cross-linked is hydrolysed
            Self::Simple(SimpleModification::Linker { formula, .. }) => (
                (formula.clone() + molecular_formula!(H 2 O 1)).into(),
                HashSet::new(),
            ),
            Self::Simple(s) => (s.formula().into(), HashSet::new()),
            Self::CrossLink {
                peptide,
                linker,
                name,
                ..
            } => {
                let link = (!applied_cross_links.contains(name))
                    .then(|| {
                        applied_cross_links.push(name.clone());
                        linker.formula()
                    })
                    .unwrap_or_default();
                if visited_peptides.contains(peptide) {
                    (link.into(), HashSet::from([name.clone()]))
                } else {
                    let (f, mut seen) = all_peptides[*peptide].formulas_inner(
                        *peptide,
                        all_peptides,
                        visited_peptides,
                        applied_cross_links,
                    );
                    seen.insert(name.clone());
                    (f + link, seen)
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
    pub fn is_possible(&self, seq: &SequenceElement, position: &PeptidePosition) -> RulePossible {
        self.simple()
            .map_or(RulePossible::Symmetric, |s| s.is_possible(seq, position))
    }
}

impl SimpleModification {
    /// Search matching modification based on what modification is provided. If a mass modification is provided
    /// it returns all modifications with that mass (within the tolerance). If a formula is provided it returns
    /// all modifications with that formula. If a glycan composition is provided it returns all glycans with
    /// that composition. Otherwise it returns the modification itself.
    pub fn search(
        modification: &Self,
        tolerance: Tolerance<Mass>,
        custom_database: Option<&CustomDatabase>,
    ) -> ModificationSearchResult {
        match modification {
            Self::Mass(mass) => ModificationSearchResult::Mass(
                mass.into_inner(),
                tolerance,
                [
                    Ontology::Unimod,
                    Ontology::Psimod,
                    Ontology::Gnome,
                    Ontology::Xlmod,
                    Ontology::Custom,
                ]
                .iter()
                .flat_map(|o| {
                    o.lookup(custom_database)
                        .iter()
                        .map(|(i, n, m)| (*o, *i, n, m))
                })
                .filter(|(_, _, _, m)| {
                    tolerance.within(&mass.into_inner(), &m.formula().monoisotopic_mass())
                })
                .map(|(o, i, n, m)| (o, i, n.clone(), m.clone()))
                .collect(),
            ),
            Self::Formula(formula) => ModificationSearchResult::Formula(
                formula.clone(),
                [
                    Ontology::Unimod,
                    Ontology::Psimod,
                    Ontology::Gnome,
                    Ontology::Xlmod,
                    Ontology::Custom,
                ]
                .iter()
                .flat_map(|o| {
                    o.lookup(custom_database)
                        .iter()
                        .map(|(i, n, m)| (*o, *i, n, m))
                })
                .filter(|(_, _, _, m)| *formula == m.formula())
                .map(|(o, i, n, m)| (o, i, n.clone(), m.clone()))
                .collect(),
            ),
            Self::Glycan(glycan) => {
                let search = MonoSaccharide::search_composition(glycan.clone());
                ModificationSearchResult::Glycan(
                    glycan.clone(),
                    Ontology::Gnome
                        .lookup(custom_database)
                        .iter()
                        .filter(|(_, _, m)| {
                            if let Self::Gno(GnoComposition::Structure(structure), _) = m {
                                MonoSaccharide::search_composition(structure.composition())
                                    == *search
                            } else {
                                false
                            }
                        })
                        .map(|(i, n, m)| (Ontology::Gnome, *i, n.clone(), m.clone()))
                        .collect(),
                )
            }
            m => ModificationSearchResult::Single(m.clone()),
        }
    }
}

/// The result of a modification search, see [`Modification::search`].
pub enum ModificationSearchResult {
    /// The modification was already defined
    Single(SimpleModification),
    /// All modifications with the same mass, within the tolerance
    Mass(
        Mass,
        Tolerance<Mass>,
        Vec<(Ontology, usize, String, SimpleModification)>,
    ),
    /// All modifications with the same formula
    Formula(
        MolecularFormula,
        Vec<(Ontology, usize, String, SimpleModification)>,
    ),
    /// All modifications with the same glycan composition
    Glycan(
        Vec<(MonoSaccharide, isize)>,
        Vec<(Ontology, usize, String, SimpleModification)>,
    ),
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
    fn formula(&self) -> MolecularFormula {
        self.modification.formula()
    }
}

impl Display for Modification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Simple(sim) => write!(f, "{sim}"),
            Self::CrossLink { name, linker, .. } => write!(f, "{linker}{name}"),
        }
        .unwrap();
        Ok(())
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
