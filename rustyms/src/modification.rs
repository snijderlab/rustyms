//! Handle modification related issues, access provided if you want to dive deeply into modifications in your own code.

use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use std::fmt::{Display, Write};

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
            Self::Predefined(formula, _, _, _, _) | Self::Linker { formula, .. } => formula.clone(),
            Self::Gno(GnoComposition::Mass(m), _) => {
                MolecularFormula::with_additional_mass(m.value)
            }
        }
    }
}

impl SimpleModification {
    /// Get a url for more information on this modification. Only defined for modifications from ontologies.
    #[allow(clippy::missing_panics_doc)]
    pub fn ontology_url(&self) -> Option<String> {
        match self {
            Self::Mass(_) | Self::Formula(_) | Self::Glycan(_) | Self::GlycanStructure(_) => None,
            Self::Predefined(_, _, ontology, name, id)
            | Self::Linker {
                ontology, name, id, ..
            } => ontology.url(*id, name),
            Self::Gno(_, name) => Ontology::Gnome.url(0, name),
        }
    }

    /// Check to see if this modification can be placed on the specified element
    pub fn is_possible(&self, seq: &SequenceElement, position: &PeptidePosition) -> bool {
        if let Self::Predefined(_, positions, _, _, _) = self {
            // If any of the rules match the current situation then it can be placed
            positions
                .iter()
                .any(|(rules, _, _)| PlacementRule::any_possible(rules, seq, position))
        } else {
            true
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
            Self::Predefined(formula, _, Ontology::Custom, name, _) => {
                write!(f, "Formula:{formula}|INFO:Custom:{name}")?;
            }
            Self::Predefined(_, _, ontology, name, _) => {
                write!(f, "{}:{name}", ontology.char())?;
            }
            Self::Gno(_, name) => write!(f, "{}:{name}", Ontology::Gnome.char())?,
            Self::Linker { ontology, name, .. } => write!(f, "{}:{name}", ontology.char())?,
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
        applied_cross_links: &mut Vec<Option<String>>,
    ) -> Multi<MolecularFormula> {
        match self {
            // A linker that is not cross-linked is hydrolysed
            Self::Simple(SimpleModification::Linker { formula, .. }) => {
                (formula.clone() + molecular_formula!(H 2 O 1)).into()
            }
            Self::Simple(s) => s.formula().into(),
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
                    link.into()
                } else {
                    all_peptides[*peptide].formulas_inner(
                        *peptide,
                        all_peptides,
                        visited_peptides,
                        applied_cross_links,
                    ) + link
                }
            } // TODO: impl neutral losses for that other peptide
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
    pub fn is_possible(&self, seq: &SequenceElement, position: &PeptidePosition) -> bool {
        self.simple().map_or(true, |s| s.is_possible(seq, position))
    }

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
            Self::Simple(SimpleModification::Mass(mass)) => ModificationSearchResult::Mass(
                mass.into_inner(),
                tolerance,
                [Ontology::Unimod, Ontology::Psimod, Ontology::Gnome]
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
            Self::Simple(SimpleModification::Formula(formula)) => {
                ModificationSearchResult::Formula(
                    formula.clone(),
                    [Ontology::Unimod, Ontology::Psimod, Ontology::Gnome]
                        .iter()
                        .flat_map(|o| {
                            o.lookup(custom_database)
                                .iter()
                                .map(|(i, n, m)| (*o, *i, n, m))
                        })
                        .filter(|(_, _, _, m)| *formula == m.formula())
                        .map(|(o, i, n, m)| (o, i, n.clone(), m.clone()))
                        .collect(),
                )
            }
            Self::Simple(SimpleModification::Glycan(glycan)) => {
                let search = MonoSaccharide::search_composition(glycan.clone());
                ModificationSearchResult::Glycan(
                    glycan.clone(),
                    Ontology::Gnome
                        .lookup(custom_database)
                        .iter()
                        .filter(|(_, _, m)| {
                            if let SimpleModification::Gno(
                                GnoComposition::Structure(structure),
                                _,
                            ) = m
                            {
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
    Single(Modification),
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
pub type AmbiguousLookup = Vec<(Option<String>, Option<SimpleModification>)>;
/// The structure to lookup cross-links, with a list of all linkers (the order is fixed) with for each linker their name or None if it is a branch and the actual linker itself (if already defined)
pub type CrossLinkLookup = Vec<(Option<String>, Option<SimpleModification>)>;

/// An ambiguous modification which could be placed on any of a set of locations
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub struct AmbiguousModification {
    /// The id to compare be able to find the other locations where this modifications can be placed
    pub id: usize,
    /// The modification itself
    pub modification: SimpleModification,
    /// If present the localisation score, meaning the chance/ratio for this modification to show up on this exact spot
    pub localisation_score: Option<OrderedFloat<f64>>,
    /// If this is a named group contain the name and track if this is the preferred location or not
    pub group: Option<(String, bool)>,
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
            Self::CrossLink { name, linker, .. } => write!(
                f,
                "{linker}#{}",
                name.as_ref()
                    .map_or("BRANCH".to_string(), |n| format!("XL{n}"))
            ),
        }
        .unwrap();
        Ok(())
    }
}

include!("shared/ontology.rs");
