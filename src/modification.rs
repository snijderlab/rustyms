//! Handle modification related issues, access provided if you want to dive deeply into modifications in your own code.

use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use std::{fmt::Display, ops::Range};

use regex::Regex;

use crate::{
    error::Context,
    error::CustomError,
    glycan::{glycan_parse_list, GlycanStructure, MonoSaccharide},
    helper_functions::*,
    ontologies::{gnome_ontology, psimod_ontology, unimod_ontology},
    placement_rule::PlacementRule,
    system::dalton,
    system::Mass,
    system::OrderedMass,
    AminoAcid, Chemical, Element, MolecularFormula, SequenceElement,
};

include!("shared/modification.rs");

impl Chemical for Modification {
    fn formula(&self) -> MolecularFormula {
        match self {
            Self::Mass(m) => MolecularFormula::with_additional_mass(m.value),
            Self::Formula(elements) => elements.clone(),
            Self::Glycan(monosaccharides) => monosaccharides
                .iter()
                .fold(MolecularFormula::default(), |acc, i| {
                    acc + i.0.formula() * i.1 as i16
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

impl Modification {
    /// Try to parse the modification. Any ambiguous modification will be number
    /// according to the lookup (which may be added to if necessary). The result
    /// is the modification, with, if applicable, its determined ambiguous group.
    /// # Errors
    /// If it is not a valid modification return a `CustomError` explaining the error.
    pub fn try_from(
        line: &str,
        range: Range<usize>,
        lookup: &mut Vec<(Option<String>, Option<Self>)>,
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
            last_result = parse_single_modification(line, part, offset, lookup);
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
                        ReturnModification::Defined(Self::Mass(OrderedMass::zero()))
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
    ) -> Option<Self> {
        match line[location.clone()].to_lowercase().as_str() {
            "o" => unimod_ontology().find_id(35),  // oxidation
            "cam" => unimod_ontology().find_id(4), // carbamidomethyl
            "pyro-glu" => unimod_ontology().find_id(
                if position.is_some_and(|p| p.aminoacid == AminoAcid::E) {
                    27
                } else {
                    28
                },
            ), // pyro Glu with the logic to pick the correct modification based on the amino acid it is placed on
            _ => {
                // Try to detect the Opair format
                Regex::new(r"[^:]+:(.*) on [A-Z]")
                    .unwrap()
                    .captures(&line[location.clone()])
                    .and_then(|capture| {
                        unimod_ontology()
                            .find_name(&capture[1])
                            .ok_or_else(|| {
                                parse_named_counter(&capture[1], glycan_parse_list(), false)
                                    .map(Modification::Glycan)
                            })
                            .flat_err()
                            .or_else(|_| {
                                match &capture[1] {
                                    "Deamidation" => Ok(unimod_ontology().find_id(7).unwrap()), // deamidated
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
                                unimod_ontology().find_name(capture[1].trim()).or_else(|| {
                                    match capture[1].trim() {
                                        "Deamidation" => {
                                            Some(unimod_ontology().find_id(7).unwrap())
                                        } // deamidated
                                        _ => None,
                                    }
                                })
                            })
                    })
            }
        }
    }

    fn sloppy_modification_internal(line: &str) -> Option<Self> {
        Self::sloppy_modification(line, 0..line.len(), None)
    }

    /// Check to see if this modification can be placed on the specified element
    pub fn is_possible(&self, seq: &SequenceElement, index: usize, length: usize) -> bool {
        if let Self::Predefined(_, rules, _, _, _) = self {
            // If any of the rules match the current situation then it can be placed
            if !rules
                .iter()
                .any(|rule| rule.is_possible(seq, index, length))
            {
                return false;
            }
        }
        true
    }
}

/// The structure to lookup ambiguous modifications, with a list of all modifications (the order is fixed) with for each modification their name and the actual modification itself (if already defined)
pub type AmbiguousLookup = Vec<(Option<String>, Option<Modification>)>;

fn parse_single_modification(
    line: &str,
    full_modification: &str,
    offset: usize,
    lookup: &mut AmbiguousLookup,
) -> Result<Option<ReturnModification>, CustomError> {
    // Parse the whole intricate structure of the single modification (see here in action: https://regex101.com/r/pW5gsj/1)
    let regex =
        Regex::new(r"^(([^:#]*)(?::([^#]+))?)(?:#([0-9A-Za-z]+)(?:\((\d+\.\d+)\))?)?$").unwrap();
    if let Some(groups) = regex.captures(full_modification) {
        // Capture the full mod name (head:tail), head, tail, ambiguous group, and localisation score
        let (full, head, tail, group, localisation_score) = (
            groups
                .get(1)
                .map(|m| (m.as_str(), m.start(), m.len()))
                .unwrap_or_default(),
            groups
                .get(2)
                .map(|m| (m.as_str().to_ascii_lowercase(), m.start(), m.len())),
            groups.get(3).map(|m| (m.as_str(), m.start(), m.len())),
            groups.get(4).map(|m| (m.as_str(), m.start(), m.len())),
            groups.get(5).map(|m| {
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
            }),
        );
        // Handle localisation score errors (could not do this inside the above closure)
        let localisation_score = if let Some(s) = localisation_score {
            Some(s?)
        } else {
            None
        };

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
                    unimod_ontology().find_id(id)
                        .map(Some)
                        .ok_or_else(||
                            basic_error
                            .with_long_description("The supplied Unimod accession number is not an existing modification"))
                }
                ("mod", tail) => {
                    let id = tail.parse::<usize>().map_err(|_| basic_error
                        .with_long_description("PSI-MOD accession number should be a number"))?;
                    psimod_ontology().find_id(id)
                        .map(Some)
                        .ok_or_else(|| basic_error
                            .with_long_description("The supplied PSI-MOD accession number is not an existing modification"))
                }
                ("u", tail) => unimod_ontology().find_name(tail)
                    .ok_or_else(|| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| basic_error
                        .with_long_description("This modification cannot be read as a Unimod name or numerical modification")),
                ("m", tail) => psimod_ontology().find_name(tail)
                    .ok_or_else(|| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| basic_error
                        .with_long_description("This modification cannot be read as a PSI-MOD name or numerical modification")),
                ("gno" | "g", tail) => gnome_ontology().find_name(tail)
                    .map(Some)
                    .ok_or_else(|| basic_error
                        .with_long_description("This modification cannot be read as a GNO name")),
                ("formula", tail) => Ok(Some(Modification::Formula(
                    parse_molecular_formula_pro_forma(tail).map_err(|e| basic_error
                        .with_long_description(format!("This modification cannot be read as a valid formula: {e}")))?,
                ))),
                ("glycan", tail) => Ok(Some(Modification::Glycan(parse_named_counter(
                    tail,
                    glycan_parse_list(),
                    false,
                ).map_err(|e| basic_error
                    .with_long_description(format!("This modification cannot be read as a valid glycan: {e}")))?))),
                ("glycanstructure", _) => GlycanStructure::parse(line, offset + tail.1..offset+tail.1+tail.2).map(|g| Some(Modification::GlycanStructure(g))),
                ("info", _) => Ok(None),
                ("obs", tail) => numerical_mod(tail).map(Some).map_err(|_| basic_error
                    .with_long_description("This modification cannot be read as a numerical modification")),
                (_, _tail) => unimod_ontology().find_name(full.0)
                    .or_else(|| psimod_ontology().find_name(full.0 ))
                    .map(Some)
                    .ok_or_else(||
                        CustomError::error(
                            "Invalid modification",
                            "Rustyms does not support these types of modification (yet)",
                            Context::line(0, line, offset+head.1, head.2),
                        )),
            }
        } else if full.0.is_empty() {
            Ok(None)
        } else {
            unimod_ontology().find_name(full.0 )
                .or_else(|| psimod_ontology().find_name(full.0))
                .ok_or_else(|| numerical_mod(full.0))
                .flat_err()
                .map(Some)
                .map_err(|_|
                    CustomError::error(
                        "Invalid modification",
                        "This modification cannot be read as a Unimod name, PSI-MOD name, or numerical modification",
                        Context::line(0, line, offset+full.1, full.2),
                    ))
        };
        // Handle ambiguous modifications
        if let Some(group) = group {
            // Search for a previous definition of this name, store as Some((index, modification_definition_present)) or None if there is no definition in place
            let found_definition = lookup
                .iter()
                .enumerate()
                .find(|(_, (name, _))| name.as_ref().map_or(false, |n| n == group.0))
                .map(|(index, (_, modification))| (index, modification.is_some()));
            // Handle all possible cases of having a modification found at this position and having a modification defined in the ambiguous lookup
            match (modification, found_definition) {
                // Have a mod defined here and already in the lookup (error)
                (Ok(Some(_)), Some((_, true))) => Err(
                    CustomError::error(
                        "Invalid ambiguous modification",
                        "An ambiguous modification cannot be placed twice (for one of the modifications leave out the modification and only provide the group name)",
                        Context::line(0, line, offset+full.1, full.2),
                    )),
                // Have a mod defined here, the name present in the lookup but not yet the mod
                (Ok(Some(m)), Some((index, false))) => {
                    lookup[index].1 = Some(m);
                    Ok(Some(ReturnModification::Preferred(index, localisation_score)))
                },
                // Have a mod defined here which is not present in the lookup
                (Ok(Some(m)), None) => {
                    let index = lookup.len();
                    lookup.push((Some(group.0.to_string()), Some(m)));
                    Ok(Some(ReturnModification::Preferred(index, localisation_score)))
                },
                // No mod defined, but the name is present in the lookup
                (Ok(None), Some((index, _))) => Ok(Some(ReturnModification::Referenced(index, localisation_score))),
                // No mod defined, and no name present in the lookup
                (Ok(None), None) => {
                    let index = lookup.len();
                    lookup.push((Some(group.0.to_string()), None));
                    Ok(Some(ReturnModification::Referenced(index, localisation_score)))},
                    // Earlier error
                    (Err(e), _) => Err(e),
            }
        } else {
            modification.map(|m| m.map(ReturnModification::Defined))
        }
    } else {
        Err(CustomError::error(
            "Invalid modification",
            "It does not match the Pro Forma definition for modifications",
            Context::line(0, line, offset, full_modification.len()),
        ))
    }
}

include!("shared/ontology.rs");

/// A modification as returned by the parser
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub enum ReturnModification {
    /// A fully self contained modification
    Defined(Modification),
    /// A modification that references an ambiguous modification
    Referenced(usize, Option<OrderedFloat<f64>>),
    /// A modification that references an ambiguous modification and is preferred on this location
    Preferred(usize, Option<OrderedFloat<f64>>),
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
    Isotope(Element, Option<u16>),
    /// Can be placed on any place it fits, if that is the correct aminoacid and it fits according to the placement rules of the modification itself
    Fixed(AminoAcid, Modification),
    /// Can be placed on any place where it can fit (according to the placement rules of the modification itself)
    Free(Modification),
}

fn numerical_mod(text: &str) -> Result<Modification, String> {
    text.parse().map_or_else(
        |_| Err("Invalid number".to_string()),
        |n| Ok(Modification::Mass(Mass::new::<dalton>(n).into())),
    )
}

impl Display for Modification {
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
            Self::GlycanStructure(glycan) => write!(f, "GlycanStructure:{glycan}",)?,
            Self::Predefined(_, _, context, name, _) => {
                write!(f, "{}:{name}", context.char())?;
            }
            Self::Gno(_, name) => write!(f, "{}:{name}", Ontology::Gnome.char())?,
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sloppy_names() {
        assert_eq!(
            Modification::sloppy_modification_internal("Deamidation (NQ)"),
            Some(unimod_ontology().find_name("deamidated").unwrap())
        );
    }
}
