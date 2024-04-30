//! Handle modification related issues, access provided if you want to dive deeply into modifications in your own code.

use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use std::{
    fmt::{Display, Write},
    num::NonZeroU16,
    ops::Range,
};

use regex::Regex;

use crate::{
    error::{Context, CustomError},
    fragment::PeptidePosition,
    glycan::{glycan_parse_list, GlycanStructure, MonoSaccharide},
    helper_functions::*,
    ontologies::CustomDatabase,
    peptide_complexity::Linked,
    placement_rule::{PlacementRule, Position},
    system::{dalton, Mass, OrderedMass},
    AminoAcid, Chemical, DiagnosticIon, Element, LinearPeptide, MolecularFormula, Multi,
    NeutralLoss, SequenceElement, Tolerance, WithinTolerance,
};

include!("shared/modification.rs");

impl Chemical for Linker {
    fn formula(&self) -> MolecularFormula {
        self.formula.clone()
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
            Self::Predefined(formula, _, _, _, _) => formula.clone(),
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
            Self::Mass(_) | Self::Formula(_) | Self::Glycan(_) | Self::GlycanStructure(_)=> None,
            Self::Predefined(_, _, ontology, _, id) => match ontology {
                Ontology::Psimod => Some(format!(
                    "https://ontobee.org/ontology/MOD?iri=http://purl.obolibrary.org/obo/MOD_{id:5}",
                )),
                Ontology::Unimod => Some(format!(
                    "https://www.unimod.org/modifications_view.php?editid1={id}",
                )),
                Ontology::Gnome | Ontology::Xlmod | Ontology::Custom => None,
            },
            Self::Gno(_, name) => Some(format!(
                "https://gnome.glyomics.org/StructureBrowser.html?focus={name}"
            )),
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
            Self::Predefined(_, _, context, name, _) => {
                write!(f, "{}:{name}", context.char())?;
            }
            Self::Gno(_, name) => write!(f, "{}:{name}", Ontology::Gnome.char())?,
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
    pub fn formula(
        &self,
        all_peptides: &[LinearPeptide<Linked>],
        visited_peptides: &[usize],
    ) -> Multi<MolecularFormula> {
        match self {
            Self::Simple(s) => s.formula().into(),
            Self::CrossLink {
                peptide, linker, ..
            } => {
                if visited_peptides.contains(peptide) {
                    linker.formula().into()
                } else {
                    all_peptides[*peptide].formulas_inner(*peptide, all_peptides, visited_peptides)
                        + linker.formula()
                }
            } // TODO: impl neutral losses for that other peptide
            Self::IntraLink { .. } => todo!(), // TODO: impossible, return None?
            Self::Branch { peptide, .. } => {
                all_peptides[*peptide].formulas_inner(*peptide, all_peptides, visited_peptides)
            } // TODO: any linking chemistry? + impl neutral losses for that other peptide
        }
    }

    /// Get the formula for a simple modification
    pub fn simple_formula(&self) -> Option<MolecularFormula> {
        match self {
            Self::Simple(s) => s.formula().into(),
            _ => None,
        }
    }
}

impl Modification {
    /// Check if this is a simple modification
    pub const fn simple(&self) -> Option<&SimpleModification> {
        match self {
            Self::Simple(sim) => Some(sim),
            _ => None,
        }
    }

    /// Check if this is a simple modification
    pub fn into_simple(self) -> Option<SimpleModification> {
        match self {
            Self::Simple(sim) => Some(sim),
            _ => None,
        }
    }

    /// Get a url for more information on this modification. Only defined for modifications from ontologies.
    #[allow(clippy::missing_panics_doc)]
    pub fn ontology_url(&self) -> Option<String> {
        match self {
            Self::Simple(modification) => modification.ontology_url(),
            Self::CrossLink { .. } => todo!(),
            Self::IntraLink { .. } => todo!(),
            Self::Branch { .. } => todo!(),
        }
    }

    /// Try to parse the modification. Any ambiguous modification will be numbered
    /// according to the lookup (which may be added to if necessary). The result
    /// is the modification, with, if applicable, its determined ambiguous group.
    /// # Errors
    /// If it is not a valid modification return a `CustomError` explaining the error.
    pub fn try_from(
        line: &str,
        range: Range<usize>,
        ambiguous_lookup: &mut AmbiguousLookup,
        cross_link_lookup: &mut CrossLinkLookup,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<ReturnModification, CustomError> {
        // Because multiple modifications could be chained with the pipe operator
        // the parsing iterates over all links until it finds one it understands
        // it then returns that one. If no 'understandable' links are found it
        // returns the last link, if this is an info it returns a mass shift of 0,
        // but if any of the links returned an error it returns the last error.
        let mut last_result = Ok(None);
        let mut last_error = None;
        let mut offset = range.start;
        for part in line[range].split('|') {
            last_result = parse_single_modification(
                line,
                part,
                offset,
                ambiguous_lookup,
                cross_link_lookup,
                custom_database,
            );
            if let Ok(Some(m)) = last_result {
                return Ok(m);
            }
            if let Err(er) = &last_result {
                last_error = Some(er.clone());
            }
            offset += part.len() + 1;
        }
        last_error.map_or_else(
            || {
                last_result.map(|m| {
                    m.unwrap_or_else(|| {
                        ReturnModification::Defined(Self::Simple(SimpleModification::Mass(
                            OrderedMass::zero(),
                        )))
                    })
                })
            },
            Err,
        )
    }

    /// Parse a modification defined by sloppy names
    #[allow(clippy::missing_panics_doc)]
    pub fn sloppy_modification(
        line: &str,
        location: std::ops::Range<usize>,
        position: Option<&SequenceElement>,
    ) -> Option<SimpleModification> {
        match line[location.clone()].to_lowercase().as_str() {
            "o" => Ontology::Unimod.find_modification_id(35, None), // oxidation
            "cam" => Ontology::Unimod.find_modification_id(4, None), // carbamidomethyl
            "pyro-glu" => Ontology::Unimod.find_modification_id(
                if position.is_some_and(|p| p.aminoacid == AminoAcid::E) {
                    27
                } else {
                    28
                },
                None,
            ), // pyro Glu with the logic to pick the correct modification based on the amino acid it is placed on
            _ => {
                // Try to detect the Opair format
                Regex::new(r"[^:]+:(.*) on [A-Z]")
                    .unwrap()
                    .captures(&line[location.clone()])
                    .and_then(|capture| {
                        Ontology::Unimod
                            .find_modification_name(&capture[1], None)
                            .ok_or_else(|| {
                                parse_named_counter(
                                    &capture[1].to_ascii_lowercase(),
                                    glycan_parse_list(),
                                    false,
                                )
                                .map(SimpleModification::Glycan)
                            })
                            .flat_err()
                            .or_else(|_| {
                                match &capture[1] {
                                    "Deamidation" => {
                                        Ok(Ontology::Unimod.find_modification_id(7, None).unwrap())
                                    } // deamidated
                                    _ => Err(()),
                                }
                            })
                            .ok()
                    })
                    .or_else(|| {
                        // Common sloppy naming: `modification (AAs)`
                        Regex::new(r"(.*)\s*\([a-zA-Z]+\)")
                            .unwrap()
                            .captures(&line[location])
                            .and_then(|capture| {
                                Ontology::Unimod
                                    .find_modification_name(capture[1].trim(), None)
                                    .or_else(|| {
                                        match capture[1].trim() {
                                            "Deamidation" => Some(
                                                Ontology::Unimod
                                                    .find_modification_id(7, None)
                                                    .unwrap(),
                                            ), // deamidated
                                            _ => None,
                                        }
                                    })
                            })
                    })
            }
        }
    }

    fn sloppy_modification_internal(line: &str) -> Option<SimpleModification> {
        Self::sloppy_modification(line, 0..line.len(), None)
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
                        o.modifications_lookup(custom_database)
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
                            o.modifications_lookup(custom_database)
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
                        .modifications_lookup(custom_database)
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
/// The structure to lookup cross-links, with a list of all linkers (the order is fixed) with for each linker their name and the actual linker itself (if already defined)
pub type CrossLinkLookup = Vec<(String, Option<Linker>)>;

/// # Errors
/// It returns an error when the given line cannot be read as a single modification.
#[allow(clippy::missing_panics_doc)]
fn parse_single_modification(
    line: &str,
    full_modification: &str,
    offset: usize,
    ambiguous_lookup: &mut AmbiguousLookup,
    cross_link_lookup: &mut CrossLinkLookup,
    custom_database: Option<&CustomDatabase>,
) -> Result<Option<ReturnModification>, CustomError> {
    // Parse the whole intricate structure of the single modification (see here in action: https://regex101.com/r/pW5gsj/1)
    let regex =
        Regex::new(r"^(([^:#]*)(?::([^#]+))?)(?:#([0-9A-Za-z]+)(?:\((\d+\.\d+)\))?)?$").unwrap();
    if let Some(groups) = regex.captures(full_modification) {
        // Capture the full mod name (head:tail), head, tail, ambiguous group, and localisation score
        let (full, head, tail, label_group, localisation_score) = (
            groups
                .get(1)
                .map(|m| (m.as_str(), m.start(), m.len()))
                .unwrap_or_default(),
            groups
                .get(2)
                .map(|m| (m.as_str().to_ascii_lowercase(), m.start(), m.len())),
            groups.get(3).map(|m| (m.as_str(), m.start(), m.len())),
            groups.get(4).map(|m| (m.as_str(), m.start(), m.len())),
            groups
                .get(5)
                .map(|m| {
                    m.as_str()
                        .parse::<f64>()
                        .map(OrderedFloat::from)
                        .map_err(|_| {
                            CustomError::error(
                        "Invalid modification localisation score",
                        "The ambiguous modification localisation score needs to be a valid number",
                        Context::line(0, line, offset + m.start(), m.len()),
                    )
                        })
                })
                .invert()?,
        );

        let modification = if let (Some(head), Some(tail)) = (head.as_ref(), tail) {
            let basic_error = CustomError::error(
                "Invalid modification",
                "..",
                Context::line(0, line, offset + tail.1, tail.2),
            );
            match (head.0.as_str(), tail.0) {
                ("unimod", tail) => {
                    let id = tail.parse::<usize>().map_err(|_| {
                        basic_error
                            .with_long_description("Unimod accession number should be a number")
                    })?;
                    Ontology::Unimod.find_modification_id(id, custom_database).map(Some).ok_or_else(|| {
                        basic_error.with_long_description(
                            "The supplied Unimod accession number is not an existing modification",
                        )
                    })
                }
                ("mod", tail) => {
                    let id = tail.parse::<usize>().map_err(|_| {
                        basic_error
                            .with_long_description("PSI-MOD accession number should be a number")
                    })?;
                    Ontology::Psimod.find_modification_id(id, custom_database).map(Some).ok_or_else(|| {
                        basic_error.with_long_description(
                            "The supplied PSI-MOD accession number is not an existing modification",
                        )
                    })
                }
                ("xlmod", tail) => {
                    let id = tail.parse::<usize>().map_err(|_| {
                        basic_error
                            .with_long_description("XLMOD accession number should be a number")
                    })?;
                    Ontology::Xlmod.find_modification_id(id, custom_database).map(Some).ok_or_else(|| {
                        basic_error.with_long_description(
                            "The supplied XLMOD accession number is not an existing modification",
                        )
                    })
                }
                ("u", tail) => Ontology::Unimod
                    .find_modification_name(tail, custom_database)
                    .ok_or_else(|| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| {
                        Ontology::Unimod
                            .find_closest(tail, custom_database)
                            .with_context(basic_error.context().clone())
                    }),
                ("m", tail) => Ontology::Psimod
                    .find_modification_name(tail, custom_database)
                    .ok_or_else(|| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| {
                        Ontology::Psimod
                            .find_closest(tail, custom_database)
                            .with_context(basic_error.context().clone())
                    }),
                ("x", tail) => Ontology::Xlmod
                    .find_modification_name(tail, custom_database)
                    .ok_or_else(|| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| {
                        Ontology::Xlmod
                            .find_closest(tail, custom_database)
                            .with_context(basic_error.context().clone())
                    }),
                ("c", tail) => Ontology::Custom
                    .find_modification_name(tail, custom_database)
                    .map(Some)
                    .ok_or_else(|| {
                        Ontology::Custom
                            .find_closest(tail, custom_database)
                            .with_context(basic_error.context().clone())
                    }),
                ("gno" | "g", tail) => Ontology::Gnome.find_modification_name(tail, custom_database).map(Some).ok_or_else(|| {
                    basic_error
                        .with_long_description("This modification cannot be read as a GNO name")
                }),
                ("formula", tail) => Ok(Some(SimpleModification::Formula(
                    MolecularFormula::from_pro_forma(tail).map_err(|e| {
                        basic_error.with_long_description(format!(
                            "This modification cannot be read as a valid formula: {e}"
                        ))
                    })?,
                ))),
                ("glycan", tail) => Ok(Some(SimpleModification::Glycan(
                    MonoSaccharide::simplify_composition(parse_named_counter(&tail.to_ascii_lowercase(), glycan_parse_list(), false).map_err(|e| {
                        basic_error.with_long_description(format!(
                            "This modification cannot be read as a valid glycan: {e}"
                        ))
                    })?),
                ))),
                ("glycanstructure", _) => {
                    GlycanStructure::parse(&line.to_ascii_lowercase(), offset + tail.1..offset + tail.1 + tail.2)
                        .map(|g| Some(SimpleModification::GlycanStructure(g)))
                }
                ("info", _) => Ok(None),
                ("obs", tail) => numerical_mod(tail).map(Some).map_err(|_| {
                    basic_error.with_long_description(
                        "This modification cannot be read as a numerical modification",
                    )
                }),
                (_, _tail) => Ontology::Unimod
                    .find_modification_name(full.0, custom_database)
                    .or_else(|| Ontology::Psimod.find_modification_name(full.0, custom_database))
                    .map(Some)
                    .ok_or_else(||
                        Ontology::find_closest_many(&[Ontology::Unimod, Ontology::Psimod], full.0, custom_database)
                    .with_long_description("This modification cannot be read as a valid Unimod or PSI-MOD name, or as a numerical modification. Or did you intent to use a ontology that is not yet supported?")
                    .with_context(Context::line(0, line, offset+full.1, full.2)))
            }
        } else if full.0.is_empty() {
            Ok(None)
        } else {
            Ontology::Unimod.find_modification_name(full.0 , custom_database)
                .or_else(|| Ontology::Psimod.find_modification_name(full.0, custom_database))
                .ok_or_else(|| numerical_mod(full.0))
                .flat_err()
                .map(Some)
                .map_err(|_|
                    Ontology::find_closest_many(&[Ontology::Unimod, Ontology::Psimod], full.0, custom_database)
                    .with_long_description("This modification cannot be read as a valid Unimod or PSI-MOD name, or as a numerical modification.")
                    .with_context(Context::line(0, line, offset+full.1, full.2))
                )
        };

        let linker = if modification.as_ref().is_err()
            || modification.as_ref().is_ok_and(Option::is_none)
        {
            if let (Some(head), Some(tail)) = (head.as_ref(), tail) {
                let basic_error = CustomError::error(
                    "Invalid linker",
                    "..",
                    Context::line(0, line, offset + tail.1, tail.2),
                );
                match (head.0.as_str(), tail.0) {
                    ("xlmod", tail) => tail
                        .parse::<usize>()
                        .map_err(|_| {
                            basic_error
                                .with_long_description("XLMOD accession number should be a number")
                        })
                        .map(|id| {
                            Ontology::Xlmod
                                .find_linker_id(id, custom_database)
                                .map(Some)
                                .ok_or_else(|| {
                                    basic_error.with_long_description(
                                    "The supplied XLMOD accession number is not an existing linker",
                                )
                                })
                        })
                        .flat_err(),
                    ("x", tail) => Ontology::Xlmod
                        .find_linker_name(tail, custom_database)
                        .map(Some)
                        .ok_or_else(|| {
                            Ontology::Xlmod
                                .find_closest(tail, custom_database)
                                .with_context(basic_error.context().clone())
                        }),
                    ("c", tail) => Ontology::Custom
                        .find_linker_name(tail, custom_database)
                        .map(Some)
                        .ok_or_else(|| {
                            Ontology::Custom
                                .find_closest(tail, custom_database)
                                .with_context(basic_error.context().clone())
                        }),
                    ("info", _) => Ok(None),
                    (_, _tail) => Err(Ontology::find_closest_many(
                        &[Ontology::Xlmod],
                        full.0,
                        custom_database,
                    )
                    .with_long_description("This linker cannot be read as a valid XLMOD name.")
                    .with_context(Context::line(
                        0,
                        line,
                        offset + full.1,
                        full.2,
                    ))),
                }
            } else {
                Err(CustomError::error(
                    "Invalid linker",
                    "A linker has to be defined with an XLMOD name or accession number.",
                    Context::line(0, line, offset + full.1, full.2),
                ))
            }
        } else {
            Ok(None)
        };

        if let Some(group) = label_group {
            if group.0.to_ascii_lowercase() == "branch" {
                Ok(Some(ReturnModification::Defined(Modification::Branch {
                    peptide: 0,
                    index: 0,
                }))) // TODO: figure out correct location, for now store empty and then fix up in a pass when the whole peptidoform is done?
            } else if let Some(name) = group.0.to_ascii_lowercase().strip_prefix("xl") {
                cross_link_lookup
                    .iter()
                    .position(|c| c.0 == name)
                    .map_or_else(
                        || {
                            let index = cross_link_lookup.len();
                            cross_link_lookup.push((name.to_string(), linker?));
                            Ok(Some(ReturnModification::CrossLinkReferenced(index)))
                        },
                        |index| Ok(Some(ReturnModification::CrossLinkReferenced(index))),
                    )
            } else {
                handle_ambiguous_modification(
                    modification,
                    group,
                    localisation_score,
                    ambiguous_lookup,
                    Context::line(0, line, offset + full.1, full.2),
                )
            }
        } else {
            modification.map(|m| m.map(|m| ReturnModification::Defined(m.into())))
        }
    } else {
        Err(CustomError::error(
            "Invalid modification",
            "It does not match the Pro Forma definition for modifications",
            Context::line(0, line, offset, full_modification.len()),
        ))
    }
}

/// Handle the logic for an ambiguous modification
/// # Errors
/// If the content of the ambiguous modification was already defined
fn handle_ambiguous_modification(
    modification: Result<Option<SimpleModification>, CustomError>,
    group: (&str, usize, usize),
    localisation_score: Option<OrderedFloat<f64>>,
    lookup: &mut AmbiguousLookup,
    context: Context,
) -> Result<Option<ReturnModification>, CustomError> {
    let group_name = group.0.to_ascii_lowercase();
    // Search for a previous definition of this name, store as Some((index, modification_definition_present)) or None if there is no definition in place
    let found_definition = lookup
        .iter()
        .enumerate()
        .find(|(_, (name, _))| name.as_ref().map_or(false, |n| n == &group_name))
        .map(|(index, (_, modification))| (index, modification.is_some()));
    // Handle all possible cases of having a modification found at this position and having a modification defined in the ambiguous lookup
    match (modification, found_definition) {
    // Have a mod defined here and already in the lookup (error)
    (Ok(Some(_)), Some((_, true))) => Err(
        CustomError::error(
            "Invalid ambiguous modification",
            "An ambiguous modification cannot be placed twice (for one of the modifications leave out the modification and only provide the group name)",
            context,
        )),
    // Have a mod defined here, the name present in the lookup but not yet the mod
    (Ok(Some(m)), Some((index, false))) => {
        lookup[index].1 = Some(m);
        Ok(Some(ReturnModification::AmbiguousPreferred(index, localisation_score)))
    },
    // Have a mod defined here which is not present in the lookup
    (Ok(Some(m)), None) => {
        let index = lookup.len();
        lookup.push((Some(group_name), Some(m)));
        Ok(Some(ReturnModification::AmbiguousPreferred(index, localisation_score)))
    },
    // No mod defined, but the name is present in the lookup
    (Ok(None), Some((index, _))) => Ok(Some(ReturnModification::AmbiguousReferenced(index, localisation_score))),
    // No mod defined, and no name present in the lookup
    (Ok(None), None) => {
        let index = lookup.len();
        lookup.push((Some(group_name), None));
        Ok(Some(ReturnModification::AmbiguousReferenced(index, localisation_score)))},
        // Earlier error
        (Err(e), _) => Err(e),
    }
}

include!("shared/ontology.rs");

/// A modification as returned by the parser
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub enum ReturnModification {
    /// A fully self contained modification
    Defined(Modification),
    /// A modification that references an ambiguous modification
    AmbiguousReferenced(usize, Option<OrderedFloat<f64>>),
    /// A modification that references an ambiguous modification and is preferred on this location
    AmbiguousPreferred(usize, Option<OrderedFloat<f64>>),
    /// A modification that references a cross-link
    CrossLinkReferenced(usize),
}

impl ReturnModification {
    /// Force this modification to be defined
    #[must_use]
    pub fn defined(self) -> Option<Modification> {
        match self {
            Self::Defined(modification) => Some(modification),
            _ => None,
        }
    }
}

/// An ambiguous modification which could be placed on any of a set of locations
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub struct AmbiguousModification {
    /// The id to compare be able to find the other locations where this modifications can be placed
    pub id: usize,
    /// The modification itself
    pub modification: Modification,
    /// If present the localisation score, meaning the chance/ratio for this modification to show up on this exact spot
    pub localisation_score: Option<OrderedFloat<f64>>,
    /// If this is a named group contain the name and track if this is the preferred location or not
    pub group: Option<(String, bool)>,
}

/// Intermediate representation of a global modification
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum GlobalModification {
    /// A global isotope modification
    Isotope(Element, Option<NonZeroU16>),
    /// Can be placed on any place it fits, if that is the correct aminoacid and it fits according to the placement rules of the modification itself
    Fixed(Position, Option<AminoAcid>, Modification),
}

/// # Errors
/// It returns an error when the text is not numerical
fn numerical_mod(text: &str) -> Result<SimpleModification, String> {
    text.parse().map_or_else(
        |_| Err("Invalid number".to_string()),
        |n| Ok(SimpleModification::Mass(Mass::new::<dalton>(n).into())),
    )
}

impl Display for Modification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Simple(sim) => write!(f, "{sim}"),
            Self::CrossLink { name, linker, .. } => write!(f, "X:{}#XL{}", linker.name, name),
            Self::IntraLink { name, linker, .. } => write!(f, "X:{}#XL{}", linker.name, name),
            Self::Branch { .. } => write!(f, "#BRANCH"),
        }
        .unwrap();
        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use super::*;

    #[test]
    fn sloppy_names() {
        assert_eq!(
            Modification::sloppy_modification_internal("Deamidation (NQ)"),
            Some(
                Ontology::Unimod
                    .find_modification_name("deamidated", None)
                    .unwrap()
            )
        );
    }
}
