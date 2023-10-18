use itertools::Itertools;

use crate::{
    error::Context,
    error::CustomError,
    helper_functions::*,
    modification::{
        AmbiguousLookup, AmbiguousModification, GlobalModification, Modification,
        ReturnModification,
    },
    molecular_charge::MolecularCharge,
    Charge, Element, Fragment, LinearPeptide, Mass, Model, MolecularFormula, SequenceElement,
};

/// A single pro forma entry, can contain multiple peptides
#[derive(Debug, Clone, PartialEq)]
pub enum ComplexPeptide {
    /// A single linear peptide
    Singular(LinearPeptide),
    /// A multimeric spectrum, multiple peptides coexist in a single spectrum indicated with '+' in pro forma
    Multimeric(Vec<LinearPeptide>),
}

impl ComplexPeptide {
    /// [Pro Forma specification](https://github.com/HUPO-PSI/ProForma)
    /// Only supports a subset of the specification (see `proforma_grammar.md` for an overview of what is supported), some functions are not possible to be represented.
    ///
    /// # Errors
    /// It fails when the string is not a valid Pro Forma string, with a minimal error message to help debug the cause.
    #[allow(clippy::too_many_lines)]
    pub fn pro_forma(value: &str) -> Result<Self, CustomError> {
        let mut peptides = Vec::new();
        let mut start = 0;
        let (peptide, tail) = Self::pro_forma_inner(value, start)?;
        start = tail;
        peptides.push(peptide);

        // Parse any following multimeric species
        while start < value.len() {
            let (peptide, tail) = Self::pro_forma_inner(value, start)?;
            peptides.push(peptide);
            start = tail;
        }
        Ok(if peptides.len() > 1 {
            Self::Multimeric(peptides)
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

    /// Read sloppy pro forma like sequences. Defined by the use of square or round braces to indicate
    /// modifications and missing any particular method of defining the N or C terminal modifications.
    /// Additionally any underscores will be ignored both on the ends and inside the sequence.
    ///
    /// All modifications follow the same definitions as the strict pro forma syntax.
    ///
    /// # Errors
    /// If it does not fit the above description.
    pub fn sloppy_pro_forma(
        line: &str,
        location: std::ops::Range<usize>,
    ) -> Result<LinearPeptide, CustomError> {
        if line[location.clone()].trim().is_empty() {
            return Err(CustomError::error(
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                Context::line(0, line, location.start, 1),
            ));
        }
        let mut peptide = LinearPeptide::default();
        let mut ambiguous_lookup = Vec::new();
        let chars: &[u8] = line[location.clone()].as_bytes();
        let mut index = 0;

        while index < chars.len() {
            match chars[index] {
                b'_' => index += 1, //ignore
                b'[' | b'(' => {
                    let (open, close) = if chars[index] == b'[' {
                        (b'[', b']')
                    } else {
                        (b'(', b')')
                    };
                    let end_index =
                        end_of_enclosure(chars, index + 1, open, close).ok_or_else(|| {
                            CustomError::error(
                                "Invalid modification",
                                "No valid closing delimiter",
                                Context::line(0, line, location.start + index, 1),
                            )
                        })?;
                    let modification = Modification::try_from(
                        line,
                        location.start + index + 1..location.start + end_index,
                        &mut ambiguous_lookup,
                    )
                    .map(|m| {
                        if let ReturnModification::Defined(m) = m {
                            Ok(m)
                        } else {
                            Err(CustomError::error(
                                "Invalid modification",
                                "A modification in the sloppy peptide format cannot be ambiguous",
                                Context::line(
                                    0,
                                    line,
                                    location.start + index + 1,
                                    end_index - 1 - index,
                                ),
                            ))
                        }
                    })
                    .flat_err()
                    .map_err(|err| {
                        Modification::sloppy_modification(
                            line,
                            location.start + index + 1..location.start + end_index,
                        )
                        .ok_or(err)
                    })
                    .flat_err()?;
                    index = end_index + 1;

                    match peptide.sequence.last_mut() {
                        Some(aa) => aa.modifications.push(modification),
                        None => peptide.n_term = Some(modification),
                    }
                }
                ch => {
                    peptide.sequence.push(SequenceElement::new(
                        ch.try_into().map_err(|_| {
                            CustomError::error(
                                "Invalid amino acid",
                                "This character is not a valid amino acid",
                                Context::line(0, line, location.start + index, 1),
                            )
                        })?,
                        None,
                    ));
                    index += 1;
                }
            }
        }
        peptide.enforce_modification_rules()?;
        Ok(peptide)
    }

    fn pro_forma_inner(
        line: &str,
        mut index: usize,
    ) -> Result<(LinearPeptide, usize), CustomError> {
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
        let mut ambiguous_found_positions = Vec::new();
        let mut global_modifications = Vec::new();
        let mut unknown_position_modifications = Vec::new();
        let mut ranged_unknown_position_modifications = Vec::new();

        // Global modification(s)
        while chars[index] == b'<' {
            let end_index = next_char(chars, index, b'>').ok_or_else(|| {
                CustomError::error(
                    "Global modification not closed",
                    "A global modification should be closed with a closing angle bracket '>'",
                    Context::line(0, line, index, 1),
                )
            })?;
            if let Some(offset) = next_char(chars, index, b'@') {
                let at_index = index + 1 + offset;
                if !chars[index + 1] == b'[' || !chars[at_index - 1] == b']' {
                    return Err(CustomError::error(
                        "Invalid global modification",
                        "A global modification should always be enclosed in square brackets '[]'",
                        Context::line(0, line, index + 1, at_index - index - 2),
                    ));
                }
                let modification =
                    Modification::try_from(line, index + 2..at_index - 2, &mut ambiguous_lookup)
                        .map(|m| {
                            if let ReturnModification::Defined(m) = m {
                                Ok(m)
                            } else {
                                Err(CustomError::error(
                                    "Invalid global modification",
                                    "A global modification cannot be ambiguous",
                                    Context::line(0, line, index + 2, at_index - index - 4),
                                ))
                            }
                        })
                        .flat_err()?;
                for aa in line[at_index..end_index].split(',') {
                    global_modifications.push(GlobalModification::Fixed(
                        aa.try_into().map_err(|_| {
                            CustomError::error(
                                "Invalid global modification",
                                "The location could not be read as an amino acid",
                                Context::line(0, line, at_index, at_index - end_index),
                            )
                        })?,
                        modification.clone(),
                    ));
                }
            } else if &line[index + 1..end_index] == "D" {
                global_modifications.push(GlobalModification::Isotope(Element::H, 2));
            } else {
                let num = &line[index + 1..end_index]
                    .chars()
                    .take_while(char::is_ascii_digit)
                    .collect::<String>();
                let el = &line[index + 1 + num.len()..end_index];
                global_modifications.push(GlobalModification::Isotope(
                    el.try_into().map_err(|()| {
                        CustomError::error(
                            "Invalid global modification",
                            "Could not determine the element",
                            Context::line(
                                0,
                                line,
                                index + num.len(),
                                index + 1 + num.len() - end_index,
                            ),
                        )
                    })?,
                    num.parse().map_err(|_| {
                        CustomError::error(
                            "Invalid global modification",
                            "Could not read isotope number",
                            Context::line(0, line, index + 1, index - end_index),
                        )
                    })?,
                ));
            }

            index = end_index + 1;
        }

        // Unknown position mods
        if let Some(result) = unknown_position_mods(chars, index) {
            let (buf, mods, ambiguous_mods) = result.map_err(|mut e| e.pop().unwrap())?; // TODO: at some point be able to return all errors
            index = buf;

            unknown_position_modifications = mods
                .into_iter()
                .map(|m| {
                    if let ReturnModification::Defined(def) = m {
                        Ok(def)
                    } else {
                        Err(CustomError::error(
                            "Invalid unknown position modification",
                            "An invalid position modification cannot be ambiguous",
                            Context::full_line(0, line),
                        ))
                    }
                })
                .collect::<Result<_, CustomError>>()?;
            ambiguous_lookup.extend(ambiguous_mods);
        }

        // Labile modification(s)
        while chars[index] == b'{' {
            let end_index = end_of_enclosure(chars, index + 1, b'{', b'}').ok_or_else(|| {
                CustomError::error(
                    "Invalid labile modification",
                    "No valid closing delimiter, a labile modification should be closed by '}'",
                    Context::line(0, line, index, 1),
                )
            })?;

            peptide.labile.push(
                Modification::try_from(line, index + 1..end_index, &mut ambiguous_lookup)
                    .map(|m| {
                        if let ReturnModification::Defined(m) = m {
                            Ok(m)
                        } else {
                            Err(CustomError::error(
                                "Invalid labile modification",
                                "A labile modification cannot be ambiguous",
                                Context::line(0, line, index + 1, end_index - 1 - index),
                            ))
                        }
                    })
                    .flat_err()?,
            );
            index = end_index + 1;
        }
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
                        if let ReturnModification::Defined(m) = m {
                            Ok(m)
                        } else {
                            Err(CustomError::error(
                                "Invalid N terminal modification",
                                "An N terminal modification cannot be ambiguous",
                                Context::line(0, line, index + 1, end_index - 2 - index),
                            ))
                        }
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
                    let (charge_len, charge) = next_num(chars, index+1, false).ok_or_else(||CustomError::error(
                        "Invalid peptide charge state",
                        "There should be a number dictating the total charge of the peptide",
                        Context::line(0, line, index+1, 1),
                    ))?;
                    if index+1+charge_len < chars.len() && chars[index+1+charge_len] == b'[' {
                        let end_index = next_char(chars, index+2+charge_len, b']').ok_or_else(||CustomError::error(
                            "Invalid adduct ion",
                            "No valid closing delimiter",
                            Context::line(0, line, index+2+charge_len, 1),
                        ))?;
                        let mut offset = index+2+charge_len;
                        let mut charge_carriers = Vec::new();
                        for set in chars[index+2+charge_len..end_index].split(|c| *c == b',') {
                            // num
                            let (count_len, count) = next_num(chars, offset, true).ok_or_else(||CustomError::error(
                                "Invalid adduct ion",
                                "Invalid adduct ion count",
                                Context::line(0, line, offset, 1),
                            ))?;
                            // element
                            let element_len =
                                chars[offset+count_len..].iter()
                                .take_while(|c| c.is_ascii_alphabetic())
                                .count();
                            let element: Element = std::str::from_utf8(&chars[offset+count_len..offset+count_len+element_len]).unwrap().try_into().map_err(|_| CustomError::error(
                                "Invalid adduct ion",
                                "Invalid element symbol",
                                Context::line(0, line, offset+count_len, element_len),
                            ))?;
                            // charge
                            let (_, ion_charge) = next_num(chars, offset+count_len+element_len, true).ok_or_else(||CustomError::error(
                                "Invalid adduct ion",
                                "Invalid adduct ion charge",
                                Context::line(0, line, offset+count_len+element_len, 1),
                            ))?;

                            charge_carriers.push((count, MolecularFormula::new(&[(element, 0, 1), (Element::Electron, 0, -ion_charge as i16)])));

                            offset += set.len() + 1;
                        }
                        peptide.charge_carriers = Some(MolecularCharge::new(charge_carriers));
                        index = end_index+1;
                        if index < chars.len() && chars[index] == b'+' {
                            index+=1; // If a peptide in a multimeric definition contains a charge state modification
                        }
                        break; // Nothing else to do, the total charge of the adduct ions should sum up to the charge as provided and this definitively is the last thing in a sequence 
                    }
                    // If no adduct ions are provided assume it is just protons
                    peptide.charge_carriers = Some(MolecularCharge::proton(charge));
                    index += charge_len+1;
                    if index < chars.len() && chars[index] == b'+' {
                        index+=1; // If a peptide in a multimeric definition contains a charge state modification
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
                            Some(if let ReturnModification::Defined(m) = modification {
                                Ok(m)
                            } else {
                                Err(CustomError::error(
                                    "Invalid C terminal modification",
                                    "A C terminal modification cannot be ambiguous",
                                    Context::line(0, line, start_index, start_index - index - 1),
                                ))
                            }?);
                            if index < chars.len() && chars[index] == b'+' {
                                index+=1; // If a peptide in a multimeric definition contains a C terminal modification
                            }
                        break; // TODO: Technically a charge state could follow the C term mods
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
                    // Multimeric spectrum stop for now, remove the plus
                    index += 1;
                    break;
                }
                (false, ch) => {
                    peptide.sequence.push(SequenceElement::new(
                        ch.try_into().map_err(|_| CustomError::error(
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
        peptide.apply_global_modifications(&global_modifications);
        peptide.apply_unknown_position_modification(&unknown_position_modifications);
        peptide.apply_ranged_unknown_position_modification(
            &ranged_unknown_position_modifications,
            &ambiguous_lookup,
        );
        peptide.enforce_modification_rules()?;

        Ok((peptide, index))
    }

    /// Assume there is exactly one peptide in this collection
    /// # Panics
    /// If there are no or multiple peptides.
    pub fn assume_linear(self) -> LinearPeptide {
        match self {
            Self::Singular(pep) => pep,
            Self::Multimeric(_) => panic!("This ComplexPeptide is not a singular linear peptide"),
        }
    }

    /// Get all peptides making up this `ComplexPeptide`, regardless of its type
    pub fn peptides(&self) -> &[LinearPeptide] {
        match self {
            Self::Singular(pep) => std::slice::from_ref(pep),
            Self::Multimeric(peptides) => peptides,
        }
    }

    /// Generate the theoretical fragments for this peptide collection.
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
    ) -> Option<Vec<Fragment>> {
        let mut base = Vec::new();
        for (index, peptide) in self.peptides().iter().enumerate() {
            base.extend(peptide.generate_theoretical_fragments(max_charge, model, index)?);
        }
        Some(base)
    }
}

/// A list of found modifications, with the newly generated ambiguous lookup alongside the index into the chars index from where parsing can be continued
type UnknownPositionMods = (usize, Vec<ReturnModification>, AmbiguousLookup);
/// If the text is recognised as a unknown mods list it is Some(..), if it has errors during parsing Some(Err(..))
/// The returned happy path contains the mods and the index from where to continue parsing.
fn unknown_position_mods(
    chars: &[u8],
    start: usize,
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
            ReturnModification::Defined(Modification::Mass(Mass::default()))
        });
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

#[cfg(test)]
mod tests {
    use crate::ComplexPeptide;
    use crate::Element;
    use crate::MolecularFormula;
    use crate::{e, mz, Location, MassOverCharge};

    use super::*;

    #[test]
    fn parse_glycan() {
        let glycan = ComplexPeptide::pro_forma("A[Glycan:Hex]")
            .unwrap()
            .assume_linear();
        let spaces = ComplexPeptide::pro_forma("A[Glycan:    Hex    ]")
            .unwrap()
            .assume_linear();
        assert_eq!(glycan.sequence.len(), 1);
        assert_eq!(spaces.sequence.len(), 1);
        assert_eq!(glycan, spaces);
        let incorrect = ComplexPeptide::pro_forma("A[Glycan:Hec]");
        assert!(incorrect.is_err());
    }

    #[test]
    fn parse_formula() {
        let peptide = ComplexPeptide::pro_forma("A[Formula:C6H10O5]")
            .unwrap()
            .assume_linear();
        let glycan = ComplexPeptide::pro_forma("A[Glycan:Hex]")
            .unwrap()
            .assume_linear();
        assert_eq!(peptide.sequence.len(), 1);
        assert_eq!(glycan.sequence.len(), 1);
        assert_eq!(glycan.formula(), peptide.formula());
    }

    #[test]
    fn parse_labile() {
        let with = ComplexPeptide::pro_forma("{Formula:C6H10O5}A")
            .unwrap()
            .assume_linear();
        let without = ComplexPeptide::pro_forma("A").unwrap().assume_linear();
        assert_eq!(with.sequence.len(), 1);
        assert_eq!(without.sequence.len(), 1);
        assert_eq!(with.formula(), without.formula());
        assert_eq!(with.labile[0].to_string(), "Formula:C6H10O5".to_string());
    }

    #[test]
    fn parse_ambiguous_modification() {
        let with = ComplexPeptide::pro_forma("A[Phospho#g0]A[#g0]")
            .unwrap()
            .assume_linear();
        let without = ComplexPeptide::pro_forma("AA").unwrap().assume_linear();
        assert_eq!(with.sequence.len(), 2);
        assert_eq!(without.sequence.len(), 2);
        assert_eq!(with.sequence[0].possible_modifications.len(), 1);
        assert_eq!(with.sequence[1].possible_modifications.len(), 1);
        assert!(ComplexPeptide::pro_forma("A[#g0]A[#g0]").is_err());
        assert!(ComplexPeptide::pro_forma("A[Phospho#g0]A[Phospho#g0]").is_err());
        assert!(ComplexPeptide::pro_forma("A[Phospho#g0]A[#g0(0.o1)]").is_err());
        assert_eq!(
            ComplexPeptide::pro_forma("A[+12#g0]A[#g0]")
                .unwrap()
                .assume_linear()
                .to_string(),
            "A[+12#g0]A[#g0]".to_string()
        );
        assert_eq!(
            ComplexPeptide::pro_forma("A[#g0]A[+12#g0]")
                .unwrap()
                .assume_linear()
                .to_string(),
            "A[#g0]A[+12#g0]".to_string()
        );
    }

    #[test]
    fn parse_ambiguous_aminoacid() {
        let with = ComplexPeptide::pro_forma("(?AA)C(?A)(?A)")
            .unwrap()
            .assume_linear();
        let without = ComplexPeptide::pro_forma("AACAA").unwrap().assume_linear();
        assert_eq!(with.sequence.len(), 5);
        assert_eq!(without.sequence.len(), 5);
        assert!(with.sequence[0].ambiguous.is_some());
        assert!(with.sequence[1].ambiguous.is_some());
        assert_eq!(with.formula(), without.formula());
        assert_eq!(with.to_string(), "(?AA)C(?A)(?A)".to_string());
    }

    #[test]
    fn parse_hard_tags() {
        let peptide = ComplexPeptide::pro_forma("A[Formula:C6H10O5|INFO:hello world 🦀]")
            .unwrap()
            .assume_linear();
        let glycan = ComplexPeptide::pro_forma(
            "A[info:you can define a tag multiple times|Glycan:Hex|Formula:C6H10O5]",
        )
        .unwrap()
        .assume_linear();
        assert_eq!(peptide.sequence.len(), 1);
        assert_eq!(glycan.sequence.len(), 1);
        assert_eq!(glycan.formula(), peptide.formula());
    }

    #[test]
    fn parse_global() {
        let deuterium = ComplexPeptide::pro_forma("<D>A").unwrap().assume_linear();
        let nitrogen_15 = ComplexPeptide::pro_forma("<15N>A").unwrap().assume_linear();
        assert_eq!(deuterium.sequence.len(), 1);
        assert_eq!(nitrogen_15.sequence.len(), 1);
        // Formula: A + H2O
        assert_eq!(
            deuterium.formula().unwrap(),
            molecular_formula!((2)H 7 C 3 O 2 N 1)
        );
        assert_eq!(
            nitrogen_15.formula().unwrap(),
            molecular_formula!(H 7 C 3 O 2 (15)N 1)
        );
    }

    #[test]
    fn parse_multimeric() {
        let dimeric = ComplexPeptide::pro_forma("A+AA").unwrap();
        let trimeric = dbg!(ComplexPeptide::pro_forma("A+AA-[+2]+AAA").unwrap());
        assert_eq!(dimeric.peptides().len(), 2);
        assert_eq!(dimeric.peptides()[0].len(), 1);
        assert_eq!(dimeric.peptides()[1].len(), 2);
        assert_eq!(trimeric.peptides().len(), 3);
        assert_eq!(trimeric.peptides()[0].len(), 1);
        assert_eq!(trimeric.peptides()[1].len(), 2);
        assert_eq!(trimeric.peptides()[2].len(), 3);
        assert!(trimeric.peptides()[1].c_term.is_some());
    }

    #[test]
    fn parse_unimod() {
        let peptide = dbg!(ComplexPeptide::pro_forma(
            "Q[U:Gln->pyro-Glu]E[Cation:Na]AA"
        ));
        assert!(peptide.is_ok());
    }

    #[test]
    fn dimeric_peptide() {
        // Only generate a single series, easier to reason about
        let test_model = Model::new(
            (Location::SkipN(1), Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            Vec::new(),
            MassOverCharge::new::<mz>(20.0),
            None,
        );

        // With two different sequences
        let dimeric = ComplexPeptide::pro_forma("AA+CC").unwrap();
        let fragments = dbg!(dimeric
            .generate_theoretical_fragments(Charge::new::<e>(1.0), &test_model)
            .unwrap());
        assert_eq!(fragments.len(), 4); // aA, aC, pAA, pCC

        // With two identical sequences
        let dimeric = ComplexPeptide::pro_forma("AA+AA").unwrap();
        let fragments = dbg!(dimeric
            .generate_theoretical_fragments(Charge::new::<e>(1.0), &test_model)
            .unwrap());
        assert_eq!(fragments.len(), 4); // aA, pAA (both twice once for each peptide)
    }

    #[test]
    fn parse_adduct_ions_01() {
        let peptide = ComplexPeptide::pro_forma("A/2[2Na+]+A").unwrap();
        assert_eq!(peptide.peptides().len(), 2);
        assert_eq!(
            peptide.peptides()[0]
                .charge_carriers
                .clone()
                .unwrap()
                .charge_carriers,
            vec![(2, molecular_formula!(Na 1 Electron -1))]
        );
        assert_eq!(
            peptide.peptides()[0].sequence,
            peptide.peptides()[1].sequence
        );
    }

    #[test]
    fn parse_adduct_ions_02() {
        let peptide = ComplexPeptide::pro_forma("A-[+1]/2[1Na+,+H+]+[+1]-A").unwrap();
        assert_eq!(peptide.peptides().len(), 2);
        assert_eq!(
            peptide.peptides()[0]
                .charge_carriers
                .clone()
                .unwrap()
                .charge_carriers,
            vec![
                (1, molecular_formula!(Na 1 Electron -1)),
                (1, molecular_formula!(H 1 Electron -1))
            ]
        );
        // Check if the C term mod is applied
        assert_eq!(
            peptide.peptides()[0].sequence[0].formula_all(),
            peptide.peptides()[1].sequence[0].formula_all()
        );
    }
}
