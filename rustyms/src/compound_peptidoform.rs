use std::{fmt::Display, num::NonZeroU16};

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    helper_functions::*,
    modification::{
        AmbiguousLookup, AmbiguousModification, GlobalModification, Modification,
        ReturnModification, SimpleModification,
    },
    molecular_charge::MolecularCharge,
    peptide_complexity::Linked,
    system::{usize::Charge, OrderedMass},
    AminoAcid, Element, Fragment, LinearPeptide, Model, MolecularFormula, Multi, MultiChemical,
    Peptidoform, SequenceElement,
};

/// A single pro forma entry, can contain multiple peptidoforms
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct CompoundPeptidoform(Vec<Peptidoform>);

impl MultiChemical for CompoundPeptidoform {
    /// Gives all possible formulas for this compound peptidoform
    fn formulas(&self) -> Multi<MolecularFormula> {
        self.0.iter().flat_map(|p| p.formulas()).collect()
    }
}

impl CompoundPeptidoform {
    /// Get all peptidoforms making up this compound peptidoform
    pub fn peptidoforms(&self) -> &[Peptidoform] {
        &self.0
    }
}

impl Display for CompoundPeptidoform {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // TODO: aggregate global, handle inconsistent global cases?
        if let Some(p) = self.0.first() {
            write!(f, "{p}")?;
        }
        for p in self.peptidoforms().iter().skip(1) {
            write!(f, "+{p}")?;
        }
        Ok(())
    }
}

impl CompoundPeptidoform {
    /// [Pro Forma specification](https://github.com/HUPO-PSI/ProForma)
    /// Supports a subset of the specification (see `proforma_grammar.md` for an overview of what is supported).
    ///
    /// # Errors
    /// It fails when the string is not a valid Pro Forma string or if it uses currently unsupported Pro Forma features.
    #[allow(clippy::too_many_lines)]
    pub fn pro_forma(value: &str) -> Result<Self, CustomError> {
        let mut peptides = Vec::new();
        let mut start = 0;
        let (peptide, tail) = Self::pro_forma_inner(value, start)?;
        start = tail;
        peptides.push(peptide);

        // Parse any following chimeric species
        while start < value.len() {
            let (peptide, tail) = Self::pro_forma_inner(value, start)?;
            peptides.push(peptide);
            start = tail;
        }
        Ok(if peptides.len() > 1 {
            Self::Chimeric(peptides)
        } else {
            Self::Singular(peptides.pop().ok_or_else(|| {
                CustomError::error(
                    "No peptide found",
                    "The peptide definition is empty",
                    Context::full_line(0, value),
                )
            })?)
        })
    }

    /// # Errors
    /// It returns an error if the line is not a supported Pro Forma line.
    #[allow(clippy::missing_panics_doc)] // Can not panic
    fn pro_forma_inner(
        line: &str,
        index: usize,
    ) -> Result<(LinearPeptide<Linked>, usize), CustomError> {
        if line.trim().is_empty() {
            return Err(CustomError::error(
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                Context::line(0, line, index, 1),
            ));
        }
        let mut peptide = LinearPeptide::default();
        let chars: &[u8] = line.as_bytes();
        let mut c_term = false;
        let mut ambiguous_aa_counter = 0;
        let mut ambiguous_aa = None;
        let mut ambiguous_lookup = Vec::new();
        let mut ambiguous_found_positions: Vec<(usize, bool, usize, Option<OrderedFloat<f64>>)> =
            Vec::new();
        let mut unknown_position_modifications = Vec::new();
        let mut ranged_unknown_position_modifications = Vec::new();

        // Global modification(s)
        let (mut index, global_modifications) = global_modifications(chars, index, line)?;

        // Unknown position mods
        if let Some(result) = unknown_position_mods(chars, index, line) {
            let (buf, mods, ambiguous_mods) =
                result.map_err(|mut e| e.pop().unwrap_or_else(|| CustomError::error("Missing error in ambiguous mods", "Ambiguous modifications could not be parsed, but no error was returned, please report this error.", Context::Show { line: line.to_string() })))?; // TODO: at some point be able to return all errors
            index = buf;

            unknown_position_modifications = mods;
            ambiguous_lookup.extend(ambiguous_mods);
        }

        // Labile modification(s)
        let (mut index, labile) = labile_modifications(chars, index, line)?;
        peptide.labile = labile;

        // N term modification
        if chars[index] == b'[' {
            let mut end_index = 0;
            for i in index..line.len() - 1 {
                if chars[i] == b']' && chars[i + 1] == b'-' {
                    end_index = i + 1;
                    break;
                }
            }
            if end_index == 0 {
                return Err(CustomError::error(
                    "Invalid N terminal modification",
                    "No valid closing delimiter, an N terminal modification should be closed by ']-'",
                    Context::line(0, line, index, 1),
                ));
            }
            peptide.n_term = Some(
                Modification::try_from(line, index + 1..end_index - 1, &mut ambiguous_lookup)
                    .map(|m| {
                        m.defined().ok_or_else(|| {
                            CustomError::error(
                                "Invalid N terminal modification",
                                "An N terminal modification cannot be ambiguous",
                                Context::line(0, line, index + 1, end_index - 2 - index),
                            )
                        })
                    })
                    .flat_err()?,
            );
            index = end_index + 1;
        }

        // Rest of the sequence
        let mut braces_start = None; // Sequence index where the last unopened braces started
        while index < chars.len() {
            match (c_term, chars[index]) {
                (false, b'(') if chars[index + 1] == b'?' && ambiguous_aa.is_none() => {
                    ambiguous_aa = Some(ambiguous_aa_counter);
                    ambiguous_aa_counter += 1;
                    index += 2;
                }
                (false, b')') if ambiguous_aa.is_some() => {
                    ambiguous_aa = None;
                    index += 1;
                }
                (false, b'(') => {
                    braces_start = Some(peptide.sequence.len());
                    index += 1;
                }
                (false, b')') if braces_start.is_some()=> {
                    braces_start = Some(peptide.sequence.len());
                    index += 1;
                    while chars[index] == b'[' {
                        let end_index = end_of_enclosure(chars, index+1, b'[', b']').ok_or_else(||CustomError::error(
                            "Invalid ranged ambiguous modification",
                            "No valid closing delimiter",
                            Context::line(0, line, index, 1),
                        ))?;
                        let modification = Modification::try_from(
                            line, index + 1..end_index,
                            &mut ambiguous_lookup,
                        )?;
                        index = end_index + 1;
                        ranged_unknown_position_modifications.push((
                            braces_start.unwrap(),
                            peptide.sequence.len(),
                            modification,
                        ));
                    }
                }
                (false, b'/') => {
                    let (buf, charge_carriers) = parse_charge_state(chars, index, line)?;
                    index = buf;
                    peptide.charge_carriers = Some(charge_carriers);
                    if index < chars.len() && chars[index] == b'+' {
                        index+=1; // Potentially this can be followed by another peptide
                    }
                    break;
                }
                (c_term, b'[') => {
                    let end_index = end_of_enclosure(chars, index+1, b'[', b']').ok_or_else(||CustomError::error(
                        "Invalid modification",
                        "No valid closing delimiter",
                        Context::line(0, line, index, 1),
                    ))?;
                    let modification = Modification::try_from(
                        line, index + 1..end_index,
                        &mut ambiguous_lookup,
                    )?;
                    let start_index = index +1;
                    index = end_index + 1;
                    if c_term {
                        peptide.c_term =
                            Some(modification.defined().ok_or_else(|| CustomError::error(
                                "Invalid C terminal modification",
                                "A C terminal modification cannot be ambiguous",
                                Context::line(0, line, start_index, start_index - index - 1),
                            ))?);
                        dbg!(index,chars.len(),&chars[index..]);
                        if index + 1 < chars.len() && chars[index] == b'/' && chars[index+1] != b'/' {
                            let (buf, charge_carriers) = parse_charge_state(chars, index, line)?;
                            index = buf;
                            peptide.charge_carriers = Some(charge_carriers);
                        }
                        if index < chars.len() && chars[index] == b'+' {
                            index+=1; // If a peptide in a chimeric definition contains a C terminal modification
                        }
                        break;
                    }
                    match peptide.sequence.last_mut() {
                        Some(aa) => match modification {
                            ReturnModification::Defined(m) => aa.modifications.push(m),
                            ReturnModification::Preferred(id, localisation_score) =>
                            ambiguous_found_positions.push(
                                (peptide.sequence.len() -1, true, id, localisation_score)),
                            ReturnModification::Referenced(id, localisation_score) =>
                            ambiguous_found_positions.push(
                                (peptide.sequence.len() -1, false, id, localisation_score)),
                        },
                        None => {
                            return Err(
                                CustomError::error(
                                    "Invalid modification",
                                    "A modification cannot be placed before any amino acid, did you want to use an N terminal modification ('[mod]-AA..')? or did you want a modification of unknown position ('[mod]?AA..')?",
                                    Context::line(0, line, start_index, start_index - index - 1),
                                )
                            )
                        }
                    }
                }
                (false, b'-') => {
                    c_term = true;
                    index += 1;
                }
                (false, b'+') => {
                    // Chimeric spectrum stop for now, remove the plus
                    index += 1;
                    break;
                }
                (false, ch) => {
                    peptide.sequence.push(SequenceElement::new(
                        ch.try_into().map_err(|()| CustomError::error(
                            "Invalid amino acid",
                            "This character is not a valid amino acid",
                            Context::line(0, line, index, 1),
                        ))?,
                        ambiguous_aa,
                    ));
                    index += 1;
                }
                (true, _) => {
                    return Err(
                        CustomError::error(
                            "Parsing error",
                            "A singular hyphen cannot exist ('-'), if this is part of a c-terminus follow the format 'AA-[modification]'",
                            Context::line(0, line, index, 1),
                        )
                    )
                }
            }
        }
        // Fill in ambiguous positions
        for (index, preferred, id, localisation_score) in ambiguous_found_positions.iter().copied()
        {
            peptide.sequence[index].possible_modifications.push(
                AmbiguousModification {
                    id,
                    modification: ambiguous_lookup[id].1.clone().ok_or_else(||
                        CustomError::error(
                            "Invalid ambiguous modification",
                            format!("Ambiguous modification {} did not have a definition for the actual modification", ambiguous_lookup[id].0.as_ref().map_or(id.to_string(), ToString::to_string)),
                            Context::full_line(0, line),
                        )
                        )?,
                    localisation_score,
                    group: ambiguous_lookup[id].0.as_ref().map(|n| (n.to_string(), preferred)) });
        }
        peptide.ambiguous_modifications = ambiguous_found_positions
            .iter()
            .copied()
            .group_by(|p| p.2)
            .into_iter()
            .sorted_by(|(key1, _), (key2, _)| key1.cmp(key2))
            .map(|(_, group)| group.into_iter().map(|p| p.0).collect())
            .collect();

        // Check all placement rules
        if !peptide.apply_global_modifications(&global_modifications) {
            return Err(CustomError::error(
                "Invalid global isotope modification",
                "There is an invalid global isotope modification",
                Context::full_line(0, line),
            ));
        }
        peptide.apply_unknown_position_modification(&unknown_position_modifications);
        peptide.apply_ranged_unknown_position_modification(
            &ranged_unknown_position_modifications,
            &ambiguous_lookup,
        );
        peptide.enforce_modification_rules()?;

        Ok((peptide, index))
    }

    /// Assume there is exactly one peptide in this collection.
    #[doc(alias = "assume_linear")]
    pub fn singular(self) -> Option<LinearPeptide<Linked>> {
        match self {
            Self::Singular(pep) => Some(pep),
            Self::Chimeric(_) => None,
        }
    }

    /// Get all peptides making up this `ComplexPeptide`, regardless of its type
    pub fn peptides(&self) -> &[LinearPeptide<Linked>] {
        match self {
            Self::Singular(pep) => std::slice::from_ref(pep),
            Self::Chimeric(peptides) => peptides,
        }
    }

    /// Generate the theoretical fragments for this peptide collection.
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
    ) -> Vec<Fragment> {
        let mut base = Vec::new();
        for peptidoform in self.peptidoforms() {
            base.extend(peptidoform.generate_theoretical_fragments(max_charge, model));
        }
        base
    }
}

impl<Item> From<Item> for CompoundPeptidoform
where
    Item: Into<LinearPeptide<Linked>>,
{
    fn from(value: Item) -> Self {
        Self::Singular(value.into())
    }
}

impl<Item> FromIterator<Item> for CompoundPeptidoform
where
    Item: Into<SequenceElement>,
{
    fn from_iter<T: IntoIterator<Item = Item>>(iter: T) -> Self {
        Self::Singular(LinearPeptide::from(iter))
    }
}

/// Parse global modifications
/// # Errors
/// If the global modifications are not defined to the specification
pub(crate) fn global_modifications(
    chars: &[u8],
    mut index: usize,
    line: &str,
) -> Result<(usize, Vec<GlobalModification>), CustomError> {
    let mut global_modifications = Vec::new();
    while index < chars.len() && chars[index] == b'<' {
        let end_index = next_char(chars, index, b'>').ok_or_else(|| {
            CustomError::error(
                "Global modification not closed",
                "A global modification should be closed with a closing angle bracket '>'",
                Context::line(0, line, index, 1),
            )
        })?;
        if let Some(offset) = next_char(chars, index, b'@') {
            let at_index = index + 1 + offset;
            dbg!(
                char::from(chars[index + 1]),
                char::from(chars[at_index - 2])
            );
            if chars[index + 1] != b'[' || chars[at_index - 2] != b']' {
                return Err(CustomError::error(
                    "Invalid global modification",
                    "A global modification should always be enclosed in square brackets '[]'",
                    Context::line(0, line, index + 1, at_index - index - 2),
                ));
            }
            let modification =
                Modification::try_from(line, index + 2..at_index - 2, &mut Vec::new())
                    .map(|m| {
                        m.defined().ok_or_else(|| {
                            CustomError::error(
                                "Invalid global modification",
                                "A global modification cannot be ambiguous",
                                Context::line(0, line, index + 2, at_index - index - 4),
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
                                    Context::line(0, line, at_index + 7, end_index - at_index - 7),
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
                                    Context::line(0, line, at_index + 7, end_index - at_index - 7),
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
                                Context::line(0, line, at_index, end_index - at_index),
                            )
                        })?),
                        modification.clone(),
                    ));
                }
            }
        } else if &line[index + 1..end_index] == "D" {
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
                        0,
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
                    Context::line(0, line, index + 1, end_index - index),
                )
            })?);
            if !el.is_valid(num) {
                return Err(CustomError::error(
                    "Invalid global modification",
                    format!(
                        "This element {el} does not have a defined weight {}",
                        num.map_or_else(String::new, |num| format!("for isotope {num}"))
                    ),
                    Context::line(0, line, index + 1, end_index - index),
                ));
            }
            global_modifications.push(GlobalModification::Isotope(el, num));
        }

        index = end_index + 1;
    }
    Ok((index, global_modifications))
}

/// A list of found modifications, with the newly generated ambiguous lookup alongside the index into the chars index from where parsing can be continued
type UnknownPositionMods = (usize, Vec<Modification>, AmbiguousLookup);
/// If the text is recognised as a unknown mods list it is Some(..), if it has errors during parsing Some(Err(..))
/// The returned happy path contains the mods and the index from where to continue parsing.
/// # Errors
/// Give all errors when the text cannot be read as mods of unknown position.
/// # Panics
/// Breaks if the text is not valid UTF-8
fn unknown_position_mods(
    chars: &[u8],
    start: usize,
    line: &str,
) -> Option<Result<UnknownPositionMods, Vec<CustomError>>> {
    let mut index = start;
    let mut modifications = Vec::new();
    let mut errs = Vec::new();
    let mut ambiguous_lookup = Vec::new();

    // Parse until no new modifications are found
    while chars[index] == b'[' {
        let end_index = next_char(chars, index + 1, b']')?;
        #[allow(clippy::map_unwrap_or)]
        // using unwrap_or can not be done because that would have a double mut ref to errs (in the eyes of the compiler)
        let modification = Modification::try_from(
            std::str::from_utf8(chars).unwrap(),
            index + 1..end_index,
            &mut ambiguous_lookup,
        )
        .unwrap_or_else(|e| {
            errs.push(e);
            ReturnModification::Defined(Modification::Simple(SimpleModification::Mass(
                OrderedMass::default(),
            )))
        })
        .defined()
        .map_or_else(
            || {
                errs.push(CustomError::error(
                    "Invalid unknown position modification",
                    "A modification of unknown position cannot be ambiguous",
                    Context::line(0, line, index + 1, end_index - 1 - index),
                ));
                Modification::Simple(SimpleModification::Mass(OrderedMass::default()))
            },
            |m| m,
        );
        index = end_index + 1;
        let number = if chars[index] == b'^' {
            if let Some((len, num)) = next_num(chars, index + 1, false) {
                index += len + 1;
                num as usize
            } else {
                errs.push(
                    CustomError::error("Invalid unknown position modification", "A modification of unknown position with multiple copies needs the copy number after the caret ('^') symbol", Context::line(0, std::str::from_utf8(chars).unwrap(), index, 1)));
                0
            }
        } else {
            1
        };
        modifications.extend(std::iter::repeat(modification).take(number));
    }
    (chars[index] == b'?').then_some({
        if errs.is_empty() {
            Ok((index + 1, modifications, ambiguous_lookup))
        } else {
            Err(errs)
        }
    })
}

/// Parse labile modifications `{mod}{mod2}`. These are assumed to fall off from the peptide in the MS.
/// # Errors
/// If the mods are not followed by a closing brace. Or if the mods are ambiguous.
fn labile_modifications(
    chars: &[u8],
    mut index: usize,
    line: &str,
) -> Result<(usize, Vec<SimpleModification>), CustomError> {
    let mut labile = Vec::new();
    while chars[index] == b'{' {
        let end_index = end_of_enclosure(chars, index + 1, b'{', b'}').ok_or_else(|| {
            CustomError::error(
                "Invalid labile modification",
                "No valid closing delimiter, a labile modification should be closed by '}'",
                Context::line(0, line, index, 1),
            )
        })?;

        labile.push(
            Modification::try_from(line, index + 1..end_index, &mut Vec::new())
                .and_then(|m| {
                    m.defined().ok_or_else(|| {
                        CustomError::error(
                            "Invalid labile modification",
                            "A labile modification cannot be ambiguous",
                            Context::line(0, line, index + 1, end_index - 1 - index),
                        )
                    })
                })
                .and_then(|m| {
                    m.simple().ok_or_else(|| {
                        CustomError::error(
                            "Invalid labile modification",
                            "A labile modification cannot be a cross-link or branch",
                            Context::line(0, line, index + 1, end_index - 1 - index),
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
pub(crate) fn parse_charge_state(
    chars: &[u8],
    index: usize,
    line: &str,
) -> Result<(usize, MolecularCharge), CustomError> {
    let (charge_len, total_charge) = next_num(chars, index + 1, false).ok_or_else(|| {
        CustomError::error(
            "Invalid peptide charge state",
            "There should be a number dictating the total charge of the peptide",
            Context::line(0, line, index + 1, 1),
        )
    })?;
    if index + 1 + charge_len < chars.len() && chars[index + 1 + charge_len] == b'[' {
        let end_index =
            end_of_enclosure(chars, index + 2 + charge_len, b'[', b']').ok_or_else(|| {
                CustomError::error(
                    "Invalid adduct ion",
                    "No valid closing delimiter",
                    Context::line(0, line, index + 2 + charge_len, 1),
                )
            })?;
        let mut offset = index + 2 + charge_len;
        let mut charge_carriers = Vec::new();
        let mut found_charge = 0;

        for set in chars[index + 2 + charge_len..end_index].split(|c| *c == b',') {
            // num
            let (count_len, count) = next_num(chars, offset, true).ok_or_else(|| {
                CustomError::error(
                    "Invalid adduct ion",
                    "Invalid adduct ion count",
                    Context::line(0, line, offset, 1),
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
                            Context::line(0, line, offset + set.len() - charge_len, charge_len),
                        )
                    })?
            };
            let (charge_len, charge) = match set[set.len() - charge_len - 1] {
                b'+' => (charge_len + 1, charge),
                b'-' => (charge_len + 1, -charge),
                _ => {
                    return Err(CustomError::error(
                        "Invalid adduct ion",
                        "The adduct ion number should be preceded by a sign",
                        Context::line(0, line, offset + set.len() - charge_len - 1, 1),
                    ))
                }
            };

            // Check for empty formula
            if count_len + charge_len == set.len() {
                return Err(CustomError::error(
                    "Invalid adduct ion",
                    "The adduct ion should have a formula defined",
                    Context::line(0, line, offset, set.len()),
                ));
            }

            // formula
            let mut formula = MolecularFormula::from_pro_forma_inner(
                line,
                offset + count_len..offset + set.len() - charge_len,
                true,
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
            found_charge += count * charge as isize;
        }
        if total_charge == found_charge {
            Ok((end_index + 1, MolecularCharge::new(&charge_carriers)))
        } else {
            Err(CustomError::error(
                "Invalid peptide charge state",
                "The peptide charge state number has to be equal to the sum of all separate adduct ions",
                Context::line(0, line, index, offset),
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
