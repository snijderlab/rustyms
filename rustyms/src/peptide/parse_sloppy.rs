use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    glycan::glycan_parse_list,
    helper_functions::{end_of_enclosure, parse_named_counter, ResultExtensions},
    modification::{Ontology, SimpleModification},
    ontologies::CustomDatabase,
    peptide::*,
    AminoAcid, SequenceElement,
};

use crate::modification::Modification;

/// Parameters to control the parsing of 'sloppy' ProForma sequences.
#[derive(Debug, Clone, Copy, Default, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct SloppyParsingParameters {
    /// Ignore a prefix lowercase n as in `n[211]GC[779]RQSSEEK` as this indicates an N terminal modification in MSFragger
    pub ignore_prefix_lowercase_n: bool,
}

impl LinearPeptide<VerySimple> {
    /// Read sloppy pro forma like sequences. Defined by the use of square or round braces to indicate
    /// modifications and missing any particular method of defining the N or C terminal modifications.
    /// Additionally any underscores will be ignored both on the ends and inside the sequence.
    ///
    /// All modifications follow the same definitions as the strict pro forma syntax, if it cannot be
    /// parsed as a strict pro forma modification it falls back to [`Modification::sloppy_modification`].
    ///
    /// # Errors
    /// If it does not fit the above description.
    #[allow(clippy::missing_panics_doc)] // Cannot panic
    pub fn sloppy_pro_forma(
        line: &str,
        location: std::ops::Range<usize>,
        custom_database: Option<&CustomDatabase>,
        parameters: SloppyParsingParameters,
    ) -> Result<Self, CustomError> {
        if line[location.clone()].trim().is_empty() {
            return Err(CustomError::error(
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                Context::line(None, line, location.start, 1),
            ));
        }
        let mut peptide = Self::default();
        let mut ambiguous_lookup = Vec::new();
        let mut cross_link_lookup = Vec::new();
        let chars: &[u8] = line[location.clone()].as_bytes();
        let mut index = 0;

        while index < chars.len() {
            match chars[index] {
                b'n' if parameters.ignore_prefix_lowercase_n && index == 0 => index += 1, //ignore
                b'_' => index += 1,                                                       //ignore
                b'[' | b'(' => {
                    let (open, close) = if chars[index] == b'[' {
                        (b'[', b']')
                    } else {
                        (b'(', b')')
                    };
                    let end_index =
                        end_of_enclosure(&line[location.clone()], index + 1, open, close)
                            .ok_or_else(|| {
                                CustomError::error(
                                    "Invalid modification",
                                    "No valid closing delimiter",
                                    Context::line(None, line, location.start + index, 1),
                                )
                            })?;
                    let modification = SimpleModification::try_from(
                        line,
                        location.start + index + 1..location.start + end_index,
                        &mut ambiguous_lookup,
                        &mut cross_link_lookup,
                        custom_database,
                    )
                    .map(|m| {
                        m.defined().ok_or_else(|| {
                            CustomError::error(
                                "Invalid modification",
                                "A modification in the sloppy peptide format cannot be ambiguous",
                                Context::line(
                                    None,
                                    line,
                                    location.start + index + 1,
                                    end_index - 1 - index,
                                ),
                            )
                        })
                    })
                    .flat_err()
                    .map_err(|err| {
                        Modification::sloppy_modification(
                            line,
                            location.start + index + 1..location.start + end_index,
                            peptide.sequence.last(),
                        )
                        .ok_or(err)
                    })
                    .flat_err()
                    .map(Modification::Simple)?;
                    index = end_index + 1;

                    match peptide.sequence.last_mut() {
                        Some(aa) => aa.modifications.push(modification),
                        None => {
                            peptide.n_term = Some(Modification::Simple(
                                modification
                                    .simple()
                                    .expect("Can only put a simple modification on an N terminus.")
                                    .clone(),
                            ));
                        }
                    }
                }
                ch => {
                    peptide.sequence.push(SequenceElement::new(
                        ch.try_into().map_err(|()| {
                            CustomError::error(
                                "Invalid amino acid",
                                "This character is not a valid amino acid",
                                Context::line(None, line, location.start + index, 1),
                            )
                        })?,
                        None,
                    ));
                    index += 1;
                }
            }
        }
        if peptide.is_empty() {
            return Err(CustomError::error(
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                Context::line(None, line, location.start, location.len()),
            ));
        }
        peptide.enforce_modification_rules()?;
        Ok(peptide)
    }
}

impl Modification {
    /// Parse a modification defined by sloppy names
    #[allow(clippy::missing_panics_doc)]
    pub fn sloppy_modification(
        line: &str,
        location: std::ops::Range<usize>,
        position: Option<&SequenceElement>,
    ) -> Option<SimpleModification> {
        match line[location.clone()].to_lowercase().as_str() {
            "o" => Ontology::Unimod.find_id(35, None),  // oxidation
            "cam" => Ontology::Unimod.find_id(4, None), // carbamidomethyl
            "pyro-glu" => Ontology::Unimod.find_id(
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
                            .find_name(&capture[1], None)
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
                                    "Deamidation" => Ok(Ontology::Unimod.find_id(7, None).unwrap()), // deamidated
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
                                    .find_name(capture[1].trim(), None)
                                    .or_else(|| {
                                        match capture[1].trim() {
                                            "Deamidation" => {
                                                Some(Ontology::Unimod.find_id(7, None).unwrap())
                                            } // deamidated
                                            _ => None,
                                        }
                                    })
                            })
                    })
            }
        }
    }

    pub(super) fn sloppy_modification_internal(line: &str) -> Option<SimpleModification> {
        Self::sloppy_modification(line, 0..line.len(), None)
    }
}
