use std::fmt::Display;

use regex::Regex;
use uom::num_traits::Zero;

use crate::{
    dalton,
    element::{Element, ELEMENT_PARSE_LIST},
    helper_functions::*,
    ontologies::{PSI_MOD_ONTOLOGY, UNIMOD_ONTOLOGY},
    placement_rules::PlacementRule,
    HasMass, Mass, MassSystem, MonoSaccharide,
};

#[derive(Debug, Clone, PartialEq)]
pub enum Modification {
    /// Monoisotopic mass shift
    Mass(Mass),
    #[allow(non_snake_case)]
    Formula(Vec<(Element, isize)>),
    Glycan(Vec<(MonoSaccharide, isize)>),
    Predefined(
        &'static [(Element, isize)],
        &'static [(MonoSaccharide, isize)],
        &'static [PlacementRule],
        &'static str, // Context
        &'static str, // Name
    ),
}

impl HasMass for Modification {
    fn mass<M: MassSystem>(&self) -> Mass {
        match self {
            Self::Mass(m) => *m,
            Self::Formula(elements) => elements.mass::<M>(),
            Self::Glycan(monosaccharides) => monosaccharides.mass::<M>(),
            Self::Predefined(elements, monosaccharides, _, _, _) => {
                elements.mass::<M>() + monosaccharides.mass::<M>()
            }
        }
    }
}

impl Modification {
    /// Try to parse the modification. Any ambiguous modification will be number
    /// according to the lookup (which may be added to if necessary). The result
    /// is the modification, with, if applicable, its determined ambiguous group.
    pub fn try_from(
        value: &str,
        lookup: &mut Vec<(Option<String>, Option<Self>)>,
    ) -> Result<ReturnModification, String> {
        // Because multiple modifications could be chained with the pipe operator
        // the parsing iterates over all links until it finds one it understands
        // it then returns that one. If no 'understandable' links are found it
        // returns the last link, if this is an info it returns a mass shift of 0,
        // but if any of the links returned an error it returns the last error.
        let mut last_result = Ok(None);
        let mut last_error = None;
        for part in value.split('|') {
            last_result = parse_single_modification(part, lookup);
            if let Ok(Some(m)) = last_result {
                return Ok(m);
            }
            if let Err(er) = &last_result {
                last_error = Some(er.clone());
            }
        }
        last_error.map_or_else(
            || {
                last_result.map(|m| {
                    m.unwrap_or_else(|| ReturnModification::Defined(Self::Mass(Mass::zero())))
                })
            },
            Err,
        )
    }
}

fn parse_single_modification(
    full_modification: &str,
    lookup: &mut Vec<(Option<String>, Option<Modification>)>,
) -> Result<Option<ReturnModification>, String> {
    // Parse the whole intricate structure of the single modification (see here in action: https://regex101.com/r/pW5gsj/1)
    let regex =
        Regex::new(r"^(([^:#]*)(?::([^#]+))?)(?:#([0-9A-Za-z]+)(?:\((\d+\.\d+)\))?)?$").unwrap();
    if let Some(groups) = regex.captures(full_modification) {
        // Capture the full mod name (head:tail), head, tail, ambiguous group, and localisation score
        let (full, head, tail, group, localisation_score) = (
            groups.get(1).map(|m| m.as_str()).unwrap_or_default(),
            groups.get(2).map(|m| m.as_str().to_ascii_lowercase()),
            groups.get(3).map(|m| m.as_str()),
            groups.get(4).map(|m| m.as_str()),
            groups.get(5).map(|m| {
                m.as_str().parse::<f64>().map_err(|_| {
                    format!(
                        "Ambiguous modification localisation score is not a valid number: \"{}\"",
                        m.as_str()
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
            match (head.as_str(), tail) {
                ("unimod", tail) => {
                    let id = tail.parse::<usize>().map_err(|e| e.to_string())?;
                    find_id_in_ontology(id, UNIMOD_ONTOLOGY)
                        .map(Some)
                        .map_err(|_| format!("{tail} is not a valid Unimod accession number"))
                }
                ("mod", tail) => {
                    let id = tail.parse::<usize>().map_err(|e| e.to_string())?;
                    find_id_in_ontology(id, PSI_MOD_ONTOLOGY)
                        .map(Some)
                        .map_err(|_| format!("{tail} is not a valid PSI-MOD accession number"))
                }
                ("u", tail) => find_name_in_ontology(tail, UNIMOD_ONTOLOGY)
                    .map_err(|_| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| format!("Not a valid Unimod modification: {tail}")),
                ("m", tail) => find_name_in_ontology(tail, PSI_MOD_ONTOLOGY)
                    .map_err(|_| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| format!("Not a valid PSI-MOD modification: {tail}")),
                ("formula", tail) => Ok(Some(Modification::Formula(parse_named_counter(
                    tail,
                    ELEMENT_PARSE_LIST,
                    true,
                )?))),
                ("glycan", tail) => Ok(Some(Modification::Glycan(parse_named_counter(
                    tail,
                    crate::GLYCAN_PARSE_LIST,
                    false,
                )?))),
                ("info", _) => Ok(None),
                ("obs", tail) => numerical_mod(tail).map(Some),
                (head, _tail) => find_name_in_ontology(full, UNIMOD_ONTOLOGY)
                    .map_err(|_| find_name_in_ontology(full, PSI_MOD_ONTOLOGY))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| format!("Does not support these types yet: {head}")),
            }
        } else if full.is_empty() {
            Ok(None)
        } else {
            find_name_in_ontology(full, UNIMOD_ONTOLOGY)
                .map_err(|_| find_name_in_ontology(full, PSI_MOD_ONTOLOGY))
                .flat_err()
                .map_err(|_| numerical_mod(full))
                .flat_err()
                .map(Some)
                .map_err(|_| {
                    format!("Not a valid delta mass, Unimod, or PSI-MOD modification: {full}")
                })
        };
        // Handle ambiguous modifications
        if let Some(group) = group {
            // Search for a previous definition of this name, store as Some((index, modification_definition_present)) or None if there is no definition in place
            let found_definition = lookup
                .iter()
                .enumerate()
                .find(|(_, (name, _))| name.as_ref().map_or(false, |n| n == group))
                .map(|(index, (_, modification))| (index, modification.is_some()));
            // Handle all possible cases of having a modification found at this position and having a modification defined in the ambiguous lookup
            match (modification, found_definition) {
                // Have a mod defined here and already in the lookup (error)
                (Ok(Some(_)), Some((_, true))) => Err("An ambiguous modification cannot be placed twice (leave out the modification and only provide the group name)".to_string()),
                // Have a mod defined here, the name present in the lookup but not yet the mod
                (Ok(Some(m)), Some((index, false))) => {
                    lookup[index].1 = Some(m);
                    Ok(Some(ReturnModification::Preferred(index, localisation_score)))
                },
                // Have a mod defined here which is not present in the lookup
                (Ok(Some(m)), None) => {
                    let index = lookup.len();
                    lookup.push((Some(group.to_string()), Some(m)));
                    Ok(Some(ReturnModification::Preferred(index, localisation_score)))
                },
                // No mod defined, but the name is present in the lookup
                (Ok(None), Some((index, _))) => Ok(Some(ReturnModification::Referenced(index, localisation_score))),
                // No mod defined, and no name present in the lookup
                (Ok(None), None) => {
                    let index = lookup.len();
                    lookup.push((Some(group.to_string()), None));
                    Ok(Some(ReturnModification::Referenced(index, localisation_score)))},
                    // Earlier error
                    (Err(e), _) => Err(e),
            }
        } else {
            modification.map(|m| m.map(ReturnModification::Defined))
        }
    } else {
        Err(format!(
            "Modification does not match pro forma pattern: \"{full_modification}\""
        ))
    }
}

pub enum ReturnModification {
    Defined(Modification),
    Referenced(usize, Option<f64>),
    Preferred(usize, Option<f64>),
}

fn find_name_in_ontology(
    code: &str,
    ontology: &[(usize, &str, Modification)],
) -> Result<Modification, ()> {
    let code = code.to_ascii_lowercase();
    for option in ontology {
        if option.1 == code {
            return Ok(option.2.clone());
        }
    }
    Err(())
}

fn find_id_in_ontology(
    id: usize,
    ontology: &[(usize, &str, Modification)],
) -> Result<Modification, ()> {
    for option in ontology {
        if option.0 == id {
            return Ok(option.2.clone());
        }
    }
    Err(())
}

fn numerical_mod(text: &str) -> Result<Modification, String> {
    text.parse().map_or_else(
        |_| Err("Invalid number".to_string()),
        |n| Ok(Modification::Mass(Mass::new::<dalton>(n))),
    )
}

impl Display for Modification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Mass(m) => {
                write!(f, "{:+}", m.value).unwrap();
            }
            Self::Formula(elements) => {
                write!(f, "Formula:{}", Element::hill_notation(elements)).unwrap();
            }
            Self::Glycan(monosaccharides) => write!(
                f,
                "Glycan:{}",
                monosaccharides
                    .iter()
                    .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
            )
            .unwrap(),
            Self::Predefined(_, _, _, context, name) => write!(f, "{context}:{name}",).unwrap(),
        }
        Ok(())
    }
}
