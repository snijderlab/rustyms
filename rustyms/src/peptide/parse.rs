use std::{collections::BTreeMap, num::NonZeroU16};

use itertools::Itertools;
use ordered_float::OrderedFloat;

use crate::{
    error::{Context, CustomError},
    helper_functions::*,
    modification::{
        AmbiguousLookup, AmbiguousModification, CrossLinkLookup, Modification, SimpleModification,
    },
    molecular_charge::MolecularCharge,
    ontologies::CustomDatabase,
    peptide::Linked,
    AminoAcid, CheckedAminoAcid, CompoundPeptidoform, Element, LinearPeptide, MolecularFormula,
    Peptidoform, SequenceElement, SequencePosition,
};

use super::{GlobalModification, Linear, ReturnModification, SemiAmbiguous};

#[derive(Debug, PartialEq, Eq)]
enum End {
    Empty,
    CrossLink,
    Chimeric,
}

struct LinearPeptideResult {
    peptide: LinearPeptide<Linear>,
    index: usize,
    ending: End,
    cross_links: Vec<(usize, SequencePosition)>,
}

impl LinearPeptide<Linked> {
    /// Convenience wrapper to parse a linear peptide in ProForma notation, to handle all possible ProForma sequences look at [`CompoundPeptidoform::pro_forma`].
    /// # Errors
    /// It gives an error when the peptide is not correctly formatted. (Also see the `CompoundPeptidoform` main function for this.)
    /// It additionally gives an error if the peptide specified was chimeric (see [`CompoundPeptidoform::singular`] and [`Peptidoform::singular`]).
    pub fn pro_forma(
        value: &str,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, CustomError> {
        CompoundPeptidoform::pro_forma(value, custom_database)?
            .singular()
            .ok_or_else(|| {
                CustomError::error(
                    "Complex peptide found",
                    "A linear peptide was expected but a chimeric peptide was found.",
                    crate::error::Context::Show {
                        line: value.to_string(),
                    },
                )
            })
            .and_then(|p| {
                p.singular().ok_or_else(|| {
                    CustomError::error(
                        "Complex peptide found",
                        "A linear peptide was expected but a cross linked peptidoform was found.",
                        crate::error::Context::Show {
                            line: value.to_string(),
                        },
                    )
                })
            })
    }
}

impl Peptidoform {
    /// Parse a peptidoform in the [ProForma specification](https://github.com/HUPO-PSI/ProForma).
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string. Or when the string has multiple peptidoforms.
    #[allow(clippy::too_many_lines)]
    pub fn pro_forma(
        value: &str,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, CustomError> {
        CompoundPeptidoform::pro_forma(value, custom_database)?
            .singular()
            .ok_or_else(|| {
                CustomError::error(
                    "Complex peptide found",
                    "A linear peptide was expected but a chimeric peptide was found.",
                    crate::error::Context::Show {
                        line: value.to_string(),
                    },
                )
            })
    }
}

impl CompoundPeptidoform {
    /// Parse a compound peptidoform in the [ProForma specification](https://github.com/HUPO-PSI/ProForma).
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string.
    #[allow(clippy::too_many_lines)]
    pub fn pro_forma(
        value: &str,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, CustomError> {
        let mut peptidoforms = Vec::new();
        // Global modification(s)
        let (mut start, global_modifications) = global_modifications(value, 0, custom_database)?;
        let (peptidoform, tail) =
            Self::parse_peptidoform(value, start, &global_modifications, custom_database)?;
        start = tail;
        peptidoforms.push(peptidoform);

        // Parse any following chimeric species
        while start < value.len() {
            let (peptidoform, tail) =
                Self::parse_peptidoform(value, start, &global_modifications, custom_database)?;
            peptidoforms.push(peptidoform);
            start = tail;
        }

        if peptidoforms.is_empty() {
            Err(CustomError::error(
                "No peptide found",
                "The peptide definition is empty",
                Context::full_line(0, value),
            ))
        } else {
            Ok(Self(peptidoforms))
        }
    }

    /// # Errors
    /// It returns an error if the line is not a supported ProForma line.
    fn parse_peptidoform(
        line: &str,
        mut index: usize,
        global_modifications: &[GlobalModification],
        custom_database: Option<&CustomDatabase>,
    ) -> Result<(Peptidoform, usize), CustomError> {
        let mut peptides = Vec::new();
        let mut ending = End::CrossLink;
        let mut cross_link_lookup = Vec::new();
        // Grouped on cross link id and stores peptide id, sequence index
        let mut cross_links_found = BTreeMap::new();

        // Parse any following cross-linked species
        while index < line.len() && ending == End::CrossLink {
            let mut result =
                Self::parse_linear_peptide(line, index, custom_database, &mut cross_link_lookup)?;
            if !result
                .peptide
                .apply_global_modifications(global_modifications)
            {
                return Err(CustomError::error(
                    "Invalid global isotope modification",
                    "There is an invalid global isotope modification",
                    Context::full_line(0, line),
                ));
            } else if result.peptide.is_empty() {
                return Err(CustomError::error(
                    "No peptide found",
                    "The peptide definition is empty",
                    Context::full_line(0, line),
                ));
            }
            peptides.push(result.peptide);
            index = result.index;
            ending = result.ending;
            for cross_link in result.cross_links {
                cross_links_found
                    .entry(cross_link.0)
                    .or_insert(Vec::new())
                    .push((peptides.len() - 1, cross_link.1));
            }
        }

        if peptides.is_empty() {
            Err(CustomError::error(
                "No peptide found",
                "The peptidoform definition is empty",
                Context::full_line(0, line),
            ))
        } else {
            let peptidoform = super::validate::cross_links(
                peptides,
                cross_links_found,
                &cross_link_lookup,
                line,
            )?;
            Ok((peptidoform, index))
        }
    }

    /// # Errors
    /// It returns an error if the line is not a supported ProForma line.
    #[allow(clippy::missing_panics_doc)] // Can not panic
    fn parse_linear_peptide(
        line: &str,
        mut index: usize,
        custom_database: Option<&CustomDatabase>,
        cross_link_lookup: &mut CrossLinkLookup,
    ) -> Result<LinearPeptideResult, CustomError> {
        if line.trim().is_empty() {
            return Err(CustomError::error(
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                Context::line(None, line, index, 1),
            ));
        }
        let mut peptide = LinearPeptide::default();
        let chars: &[u8] = line.as_bytes();
        let mut c_term = false;
        let mut ambiguous_aa_counter = 0;
        let mut ambiguous_aa = None;
        let mut ambiguous_lookup = Vec::new();
        let mut cross_link_found_positions: Vec<(usize, SequencePosition)> = Vec::new();
        let mut ambiguous_found_positions: Vec<(usize, bool, usize, Option<OrderedFloat<f64>>)> =
            Vec::new();
        let mut unknown_position_modifications = Vec::new();
        let mut ranged_unknown_position_modifications = Vec::new();
        let mut ending = End::Empty;

        // Unknown position mods
        if let Some(result) =
            unknown_position_mods(chars, index, line, custom_database, &mut ambiguous_lookup)
        {
            let (buf, mods) = result.map_err(|errors| {
                CustomError::error(
                    "Some unknown position modifications are invalid",
                    "See the underlying errors for more details.",
                    Context::Show {
                        line: line.to_string(),
                    },
                )
                .with_underlying_errors(errors)
            })?;
            index = buf;

            unknown_position_modifications = mods;
        }

        // Labile modification(s)
        let (mut index, labile) = labile_modifications(line, index, custom_database)?;
        peptide = peptide.labile(labile);

        // N term modification
        if chars.get(index) == Some(&b'[') {
            let end_index = end_of_enclosure(line, index+1, b'[', b']').and_then(|i| (chars.get(i+1) == Some(&b'-')).then_some(i+1)).ok_or_else(|| CustomError::error(
                    "Invalid N terminal modification",
                    "No valid closing delimiter, an N terminal modification should be closed by ']-'",
                    Context::line(None, line, index, 1),
                ))?;
            peptide.set_simple_n_term(
                SimpleModification::try_from(
                    line,
                    index + 1..end_index - 1,
                    &mut ambiguous_lookup,
                    cross_link_lookup,
                    custom_database,
                )
                .and_then(|m| match m {
                    ReturnModification::Defined(simple) => Ok(Some(simple)),
                    ReturnModification::CrossLinkReferenced(id) => {
                        cross_link_found_positions.push((id, SequencePosition::NTerm));
                        Ok(None)
                    }
                    ReturnModification::AmbiguousPreferred(_, _)
                    | ReturnModification::AmbiguousReferenced(_, _) => Err(CustomError::error(
                        "Invalid N terminal modification",
                        "An N terminal modification cannot be ambiguous",
                        Context::line(None, line, index + 1, end_index - 2 - index),
                    )),
                })?,
            );
            index = end_index + 1;
        }

        // Rest of the sequence
        let mut braces_start = None; // Sequence index where the last unopened braces started
        while index < chars.len() {
            match (c_term, chars[index]) {
                (false, b'(') if chars.get(index + 1) == Some(&b'?') => {
                    if braces_start.is_some() {
                        return Err(CustomError::error(
                            "Invalid ambiguous amino acid set",
                            "Ambiguous amino acid sets cannot be nested within ranged ambiguous modifications",
                            Context::line(None, line, index, 1),
                        ));
                    }
                    if ambiguous_aa.is_some() {
                        return Err(CustomError::error(
                            "Invalid ambiguous amino acid set",
                            "Ambiguous amino acid sets cannot be nested within ambiguous amino acid sets",
                            Context::line(None, line, index, 1),
                        ));
                    }
                    ambiguous_aa = Some(ambiguous_aa_counter);
                    ambiguous_aa_counter += 1;
                    index += 2;
                }
                (false, b')') if ambiguous_aa.is_some() => {
                    ambiguous_aa = None;
                    index += 1;
                }
                (false, b'(') => {
                    if braces_start.is_some() {
                        return Err(CustomError::error(
                            "Invalid ranged ambiguous modification",
                            "Ranged ambiguous modifications cannot be nested within ranged ambiguous modifications",
                            Context::line(None, line, index, 1),
                        ));
                    }
                    if ambiguous_aa.is_some() {
                        return Err(CustomError::error(
                            "Invalid ranged ambiguous modification",
                            "Ranged ambiguous modifications cannot be nested within ambiguous amino acid sets",
                            Context::line(None, line, index, 1),
                        ));
                    }
                    braces_start = Some(peptide.len());
                    index += 1;
                }
                (false, b')') if braces_start.is_some() => {
                    let start = braces_start.unwrap();
                    braces_start = None;
                    index += 1;
                    while chars.get(index) == Some(&b'[') {
                        let end_index = end_of_enclosure(line, index+1, b'[', b']').ok_or_else(||CustomError::error(
                            "Invalid ranged ambiguous modification",
                            "No valid closing delimiter",
                            Context::line(None, line, index, 1),
                        ))?;
                        let modification = SimpleModification::try_from(
                            line, index + 1..end_index,
                            &mut ambiguous_lookup, cross_link_lookup, custom_database,
                        )?.defined().ok_or_else(|| CustomError::error(
                            "Invalid ranged ambiguous modification",
                            "A ranged ambiguous modification has to be fully defined, so no ambiguous modification is allowed",
                            Context::line(None, line, index, 1),
                        ))?;
                        index = end_index + 1;
                        ranged_unknown_position_modifications.push((
                            start,
                            peptide.len().saturating_sub(1),
                            modification,
                        ));
                    }
                }
                (false, b'/') => {
                    // Chimeric peptide
                    if chars.get(index+1) == Some(&b'/') {
                        index += 2; // Potentially this can be followed by another peptide
                        ending = End::CrossLink;
                    } else {
                        let (buf, charge_carriers) = parse_charge_state(line, index)?;
                        index = buf;
                        peptide = peptide.charge_carriers(Some(charge_carriers));
                        if index < chars.len() && chars[index] == b'+' {
                            index += 1; // Potentially this can be followed by another peptide
                            ending = End::Chimeric;
                        }
                    }
                    break;
                }
                (is_c_term, b'[') => {
                    let end_index = end_of_enclosure(line, index+1, b'[', b']').ok_or_else(||CustomError::error(
                        "Invalid modification",
                        "No valid closing delimiter",
                        Context::line(None, line, index, 1),
                    ))?;
                    let modification = SimpleModification::try_from(
                        line, index + 1..end_index,
                        &mut ambiguous_lookup, cross_link_lookup, custom_database,
                    )?;
                    let start_index = index +1;
                    index = end_index + 1;
                    if is_c_term {
                        peptide = peptide.c_term(
                            match modification {
                                ReturnModification::Defined(simple) => Ok(Some(Modification::Simple(simple))),
                                ReturnModification::CrossLinkReferenced(id) =>
                                    {cross_link_found_positions.push((id, SequencePosition::CTerm)); Ok(None)},
                                ReturnModification::AmbiguousPreferred(_, _) |
                                    ReturnModification::AmbiguousReferenced(_,_) => Err(CustomError::error(
                                        "Invalid C terminal modification",
                                        "A C terminal modification cannot be ambiguous",
                                        Context::line(None, line, index + 1, end_index.saturating_sub(2 + index)),
                                    )),
                            }?);

                        if index + 1 < chars.len() && chars[index] == b'/' && chars[index+1] != b'/' {
                            let (buf, charge_carriers) = parse_charge_state(line, index)?;
                            index = buf;
                            peptide = peptide.charge_carriers(Some(charge_carriers));
                        }
                        if index < chars.len() && chars[index] == b'+' {
                            index += 1; // If a peptide in a chimeric definition contains a C terminal modification
                            ending = End::Chimeric;
                        } else if index + 1 < chars.len() && chars[index..=index+1] == *b"//" {
                            index += 2; // If a peptide in a cross-linked definition contains a C terminal modification
                            ending = End::CrossLink;
                        }
                        c_term = false; // Fix false negative errors on single ending hyphen
                        break;
                    }

                    if let Some((sequence_index, aa)) = peptide.sequence_mut().iter_mut().enumerate().next_back() {
                        match modification {
                            ReturnModification::Defined(m) => aa.modifications.push(Modification::Simple(m)),
                            ReturnModification::AmbiguousPreferred(id, localisation_score) =>
                                ambiguous_found_positions.push((sequence_index, true, id, localisation_score)),
                            ReturnModification::AmbiguousReferenced(id, localisation_score) =>
                                ambiguous_found_positions.push((sequence_index, false, id, localisation_score)),
                            ReturnModification::CrossLinkReferenced(id) =>
                                cross_link_found_positions.push((id, SequencePosition::Index(sequence_index))),
                        }
                    } else {
                        return Err(
                            CustomError::error(
                                "Invalid modification",
                                "A modification cannot be placed before any amino acid, did you want to use an N terminal modification ('[mod]-AA..')? or did you want a modification of unknown position ('[mod]?AA..')?",
                                Context::line(None, line, start_index, index - start_index - 1),
                            )
                        )
                    }
                }
                (false, b'-') => {
                    c_term = true;
                    index += 1;
                }
                (false, b'+') => {
                    // Chimeric spectrum stop for now, remove the plus
                    index += 1;
                    ending = End::Chimeric;
                    break;
                }
                (false, ch) => {
                    peptide.sequence_mut().push(SequenceElement::new(
                        CheckedAminoAcid::<SemiAmbiguous>::try_from(ch).map_err(|()| CustomError::error(
                            "Invalid amino acid",
                            "This character is not a valid amino acid",
                            Context::line(None, line, index, 1),
                        ))?.into(),
                        ambiguous_aa,
                    ));
                    index += 1;
                }
                (true, _) => {
                    return Err(
                        CustomError::error(
                            "Parsing error",
                            "A singular hyphen cannot exist ('-'), if this is part of a c-terminus follow the format 'AA-[modification]'",
                            Context::line(None, line, index, 1),
                        )
                    )
                }
            }
        }
        if c_term {
            return Err(CustomError::error(
                "Invalid peptide",
                "A single hyphen cannot end the definition, if a C terminal modification is intended use 'SEQ-[MOD]'",
                Context::line(None, line, line.len().saturating_sub(2), 1),
            ));
        }
        if let Some(pos) = braces_start {
            return Err(CustomError::error(
                "Invalid peptide",
                format!("Unclosed brace at amino acid position {pos}"),
                Context::full_line(0, line),
            ));
        }
        if ambiguous_aa.is_some() {
            return Err(CustomError::error(
                "Invalid peptide",
                "Unclosed ambiguous amino acid group",
                Context::full_line(0, line),
            ));
        }
        if peptide.is_empty() {
            return Err(CustomError::error(
                "No amino acids found",
                "The peptide definition is empty",
                Context::full_line(0, line),
            ));
        }

        // Fill in ambiguous positions, ambiguous contains (index, preferred, id, localisation_score)
        for (id, ambiguous) in ambiguous_found_positions
            .into_iter()
            .into_group_map_by(|aa| aa.2)
        {
            let modification = AmbiguousModification {
                modification: ambiguous_lookup[id].1.clone().ok_or_else(||
                CustomError::error(
                    "Invalid ambiguous modification",
                    format!("Ambiguous modification {} did not have a definition for the actual modification", ambiguous_lookup[id].0),
                    Context::full_line(0, line),
                )
                )?,
                id,
                group: ambiguous_lookup[id].0.clone(),
                localisation_score: None,
                preferred: false,
            };
            let positions = ambiguous
                .iter()
                .map(|(index, _, _, score)| (*index, *score))
                .collect_vec();
            let preferred = ambiguous.iter().find(|p| p.1).map(|p| p.0);
            peptide.add_ambiguous_modification(&modification, &positions, preferred);
        }

        peptide.apply_unknown_position_modification(&unknown_position_modifications)?;
        peptide.apply_ranged_unknown_position_modification(
            &ranged_unknown_position_modifications,
            ambiguous_lookup.len(),
            unknown_position_modifications.len(),
        )?;
        peptide.enforce_modification_rules()?;

        Ok(LinearPeptideResult {
            peptide,
            index,
            ending,
            cross_links: cross_link_found_positions,
        })
    }
}

/// Parse global modifications
/// # Errors
/// If the global modifications are not defined to the specification
pub(super) fn global_modifications(
    line: &str,
    mut index: usize,
    custom_database: Option<&CustomDatabase>,
) -> Result<(usize, Vec<GlobalModification>), CustomError> {
    let chars = line.as_bytes();
    let mut global_modifications = Vec::new();
    while index < chars.len() && chars[index] == b'<' {
        let end_index =
            end_of_enclosure_with_brackets(line, index + 1, b'<', b'>').ok_or_else(|| {
                CustomError::error(
                    "Global modification not closed",
                    "A global modification should be closed with a closing angle bracket '>'",
                    Context::line(None, line, index, 1),
                )
            })?;
        if let Some(offset) = next_char(chars, index, b'@') {
            let at_index = index + 1 + offset;
            if at_index > end_index {
                return Err(CustomError::error(
                    "Invalid global modification",
                    "A global modification should have an at '@' sign inside the enclosing angle brackets '<>'",
                    Context::line(None, line, index + 1, at_index - index - 2),
                ));
            }
            if chars[index + 1] != b'[' || chars[at_index - 2] != b']' {
                return Err(CustomError::error(
                    "Invalid global modification",
                    "A global modification should always be enclosed in square brackets '[]'",
                    Context::line(None, line, index + 1, at_index - index - 2),
                ));
            }
            let modification = SimpleModification::try_from(
                line,
                index + 2..at_index - 2,
                &mut Vec::new(),
                &mut Vec::new(),
                custom_database,
            )
            .map(|m| {
                m.defined().ok_or_else(|| {
                    CustomError::error(
                        "Invalid global modification",
                        "A global modification cannot be ambiguous or a cross-linker",
                        Context::line(None, line, index + 2, at_index - index - 4),
                    )
                })
            })
            .flat_err()?;
            for aa in line[at_index..end_index].split(',') {
                if aa.to_ascii_lowercase().starts_with("n-term") {
                    if let Some((_, aa)) = aa.split_once(':') {
                        global_modifications.push(GlobalModification::Fixed(
                            crate::placement_rule::Position::AnyNTerm,
                            Some(TryInto::<AminoAcid>::try_into(aa).map_err(|()| {
                                CustomError::error(
                                    "Invalid global modification",
                                    "The location could not be read as an amino acid",
                                    Context::line(
                                        None,
                                        line,
                                        at_index + 7,
                                        end_index - at_index - 7,
                                    ),
                                )
                            })?),
                            modification.clone(),
                        ));
                    } else {
                        global_modifications.push(GlobalModification::Fixed(
                            crate::placement_rule::Position::AnyNTerm,
                            None,
                            modification.clone(),
                        ));
                    }
                } else if aa.to_ascii_lowercase().starts_with("c-term") {
                    if let Some((_, aa)) = aa.split_once(':') {
                        global_modifications.push(GlobalModification::Fixed(
                            crate::placement_rule::Position::AnyCTerm,
                            Some(TryInto::<AminoAcid>::try_into(aa).map_err(|()| {
                                CustomError::error(
                                    "Invalid global modification",
                                    "The location could not be read as an amino acid",
                                    Context::line(
                                        None,
                                        line,
                                        at_index + 7,
                                        end_index - at_index - 7,
                                    ),
                                )
                            })?),
                            modification.clone(),
                        ));
                    } else {
                        global_modifications.push(GlobalModification::Fixed(
                            crate::placement_rule::Position::AnyCTerm,
                            None,
                            modification.clone(),
                        ));
                    }
                } else {
                    global_modifications.push(GlobalModification::Fixed(
                        crate::placement_rule::Position::Anywhere,
                        Some(TryInto::<AminoAcid>::try_into(aa).map_err(|()| {
                            CustomError::error(
                                "Invalid global modification",
                                "The location could not be read as an amino acid",
                                Context::line(None, line, at_index, end_index - at_index),
                            )
                        })?),
                        modification.clone(),
                    ));
                }
            }
        } else if &line[index + 1..end_index].to_ascii_lowercase() == "d" {
            global_modifications.push(GlobalModification::Isotope(Element::H, NonZeroU16::new(2)));
        } else {
            let num = &line[index + 1..end_index]
                .chars()
                .take_while(char::is_ascii_digit)
                .collect::<String>();
            let el = &line[index + 1 + num.len()..end_index];
            let el: Element = el.try_into().map_err(|()| {
                CustomError::error(
                    "Invalid global modification",
                    "Could not determine the element",
                    Context::line(
                        None,
                        line,
                        index + num.len(),
                        end_index - (index + 1 + num.len()),
                    ),
                )
            })?;
            let num = Some(num.parse::<NonZeroU16>().map_err(|err| {
                CustomError::error(
                    "Invalid global modification",
                    format!("The isotope number is {}", explain_number_error(&err)),
                    Context::line(None, line, index + 1, end_index - index),
                )
            })?);
            if !el.is_valid(num) {
                return Err(CustomError::error(
                    "Invalid global modification",
                    format!(
                        "This element {el} does not have a defined weight {}",
                        num.map_or_else(String::new, |num| format!("for isotope {num}"))
                    ),
                    Context::line(None, line, index + 1, end_index - index),
                ));
            }
            global_modifications.push(GlobalModification::Isotope(el, num));
        }

        index = end_index + 1;
    }
    Ok((index, global_modifications))
}

/// A list of found modifications
type UnknownPositionMods = (usize, Vec<SimpleModification>);
/// If the text is recognised as a unknown mods list it is Some(..), if it has errors during parsing Some(Err(..))
/// The returned happy path contains the mods and the index from where to continue parsing.
/// # Errors
/// Give all errors when the text cannot be read as mods of unknown position.
/// # Panics
/// Breaks if the text is not valid UTF-8
pub(super) fn unknown_position_mods(
    chars: &[u8],
    start: usize,
    line: &str,
    custom_database: Option<&CustomDatabase>,
    ambiguous_lookup: &mut AmbiguousLookup,
) -> Option<Result<UnknownPositionMods, Vec<CustomError>>> {
    let mut index = start;
    let mut modifications = Vec::new();
    let mut errs = Vec::new();
    let mut cross_link_lookup = Vec::new();

    // Parse until no new modifications are found
    while chars.get(index) == Some(&b'[') {
        let start_index = index;
        index = next_char(chars, index + 1, b']')? + 1;
        let modification = match SimpleModification::try_from(
            std::str::from_utf8(chars).unwrap(),
            start_index + 1..index - 1,
            ambiguous_lookup,
            &mut cross_link_lookup,
            custom_database,
        ) {
            Ok(ReturnModification::Defined(m)) => m,
            Ok(
                ReturnModification::AmbiguousPreferred(_, _)
                | ReturnModification::AmbiguousReferenced(_, _),
            ) => continue, // Added to list and that is the only thing to do
            Ok(ReturnModification::CrossLinkReferenced(_)) => {
                errs.push(CustomError::error(
                    "Invalid unknown position modification",
                    "A modification of unknown position cannot be a cross-link",
                    Context::line_range(None, line, start_index + 1..index),
                ));
                continue;
            }
            Err(e) => {
                errs.push(e);
                continue;
            }
        };
        let number = if chars.get(index) == Some(&b'^') {
            if let Some((len, num)) = next_num(chars, index + 1, false) {
                index += len + 1;
                if num < 0 {
                    errs.push(
                        CustomError::error("Invalid unknown position modification", "A modification of unknown position with multiple copies cannot have more a negative number of copies", Context::line(None, std::str::from_utf8(chars).unwrap(), index, 1)));
                    0
                } else if num > i16::MAX as isize {
                    errs.push(
                        CustomError::error("Invalid unknown position modification", format!("A modification of unknown position with multiple copies cannot have more then {} copies", i16::MAX), Context::line(None, std::str::from_utf8(chars).unwrap(), index, 1)));
                    0
                } else {
                    num as usize
                }
            } else {
                errs.push(
                    CustomError::error("Invalid unknown position modification", "A modification of unknown position with multiple copies needs the copy number after the caret ('^') symbol", Context::line(None, std::str::from_utf8(chars).unwrap(), index, 1)));
                0
            }
        } else {
            1
        };
        modifications.extend(std::iter::repeat(modification).take(number));
    }
    (chars.get(index) == Some(&b'?')).then_some({
        if errs.is_empty() {
            Ok((index + 1, modifications))
        } else {
            Err(errs)
        }
    })
}

/// Parse labile modifications `{mod}{mod2}`. These are assumed to fall off from the peptide in the MS.
/// # Errors
/// If the mods are not followed by a closing brace. Or if the mods are ambiguous.
fn labile_modifications(
    line: &str,
    mut index: usize,
    custom_database: Option<&CustomDatabase>,
) -> Result<(usize, Vec<SimpleModification>), CustomError> {
    let chars = line.as_bytes();
    let mut labile = Vec::new();
    while chars.get(index) == Some(&b'{') {
        let end_index = end_of_enclosure(line, index + 1, b'{', b'}').ok_or_else(|| {
            CustomError::error(
                "Invalid labile modification",
                "No valid closing delimiter, a labile modification should be closed by '}'",
                Context::line(None, line, index, 1),
            )
        })?;

        labile.push(
            SimpleModification::try_from(
                line,
                index + 1..end_index,
                &mut Vec::new(),
                &mut Vec::new(),
                custom_database,
            )
            .and_then(|m| {
                m.defined().ok_or_else(|| {
                    CustomError::error(
                        "Invalid labile modification",
                        "A labile modification cannot be ambiguous or a cross-linker",
                        Context::line(None, line, index + 1, end_index - 1 - index),
                    )
                })
            })?,
        );
        index = end_index + 1;
    }
    Ok((index, labile))
}

/// Parse a charge state `/2` or more complex ones like `/2[+2Na+]`.
/// Assumes the text starts with `/`.
/// # Errors
/// If the charge state is not following the specification.
/// # Panics
/// Panics if the text is not UTF-8.
pub(super) fn parse_charge_state(
    line: &str,
    index: usize,
) -> Result<(usize, MolecularCharge), CustomError> {
    let chars = line.as_bytes();
    let (charge_len, total_charge) = next_num(chars, index + 1, false).ok_or_else(|| {
        CustomError::error(
            "Invalid peptide charge state",
            "There should be a number dictating the total charge of the peptide",
            Context::line(None, line, index + 1, 1),
        )
    })?;
    if chars.get(index + 1 + charge_len) == Some(&b'[') {
        let end_index =
            end_of_enclosure(line, index + 2 + charge_len, b'[', b']').ok_or_else(|| {
                CustomError::error(
                    "Invalid adduct ion",
                    "No valid closing delimiter",
                    Context::line(None, line, index + 2 + charge_len, 1),
                )
            })?;
        let mut offset = index + 2 + charge_len;
        let mut charge_carriers = Vec::new();
        let mut found_charge: isize = 0;

        for set in chars[index + 2 + charge_len..end_index].split(|c| *c == b',') {
            // num
            let (count_len, count) = next_num(chars, offset, true).ok_or_else(|| {
                CustomError::error(
                    "Invalid adduct ion",
                    "Invalid adduct ion count",
                    Context::line(None, line, offset, 1),
                )
            })?;

            // charge
            let charge_len = set.iter().rev().take_while(|c| c.is_ascii_digit()).count();
            let charge = if charge_len == 0 {
                1
            } else {
                line[offset + set.len() - charge_len..offset + set.len()]
                    .parse::<i32>()
                    .map_err(|err| {
                        CustomError::error(
                            "Invalid adduct ion",
                            format!("The adduct ion number {err}"),
                            Context::line(None, line, offset + set.len() - charge_len, charge_len),
                        )
                    })?
            };
            let (charge_len, charge) = match (set.len() - charge_len)
                .checked_sub(1)
                .and_then(|i| set.get(i))
            {
                Some(b'+') => (charge_len + 1, charge),
                Some(b'-') => (charge_len + 1, -charge),
                _ => {
                    return Err(CustomError::error(
                        "Invalid adduct ion",
                        "The adduct ion number should be preceded by a sign",
                        Context::line(None, line, offset + set.len() - charge_len - 1, 1),
                    ))
                }
            };

            // Check for empty formula
            if count_len + charge_len == set.len() {
                return Err(CustomError::error(
                    "Invalid adduct ion",
                    "The adduct ion should have a formula defined",
                    Context::line(None, line, offset, set.len()),
                ));
            }

            // formula
            let mut formula = MolecularFormula::from_pro_forma(
                line,
                offset + count_len..offset + set.len() - charge_len,
                true,
                false,
            )?;
            let _ = formula.add((
                Element::Electron,
                None,
                formula.charge().value as i32 - charge,
            ));

            // Deduplicate
            if let Some((amount, _)) = charge_carriers.iter_mut().find(|(_, f)| *f == formula) {
                *amount += count;
            } else {
                charge_carriers.push((count, formula));
            }

            offset += set.len() + 1;
            found_charge = found_charge
                .checked_add(count.checked_mul(charge as isize).ok_or_else(|| {
                    CustomError::error(
                        "Invalid peptide charge state",
                        "The peptide charge state is too big to store inside an isize",
                        Context::line(None, line, index, offset),
                    )
                })?)
                .ok_or_else(|| {
                    CustomError::error(
                        "Invalid peptide charge state",
                        "The peptide charge state is too big to store inside an isize",
                        Context::line(None, line, index, offset),
                    )
                })?;
        }
        if total_charge == found_charge {
            Ok((end_index + 1, MolecularCharge::new(&charge_carriers)))
        } else {
            Err(CustomError::error(
                "Invalid peptide charge state",
                "The peptide charge state number has to be equal to the sum of all separate adduct ions",
                Context::line(None, line, index, offset),
            ))
        }
    } else {
        // If no adduct ions are provided assume it is just protons
        Ok((
            index + charge_len + 1,
            MolecularCharge::proton(total_charge),
        ))
    }
}
