use regex::Regex;

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
    ) -> Result<Self, CustomError> {
        if line[location.clone()].trim().is_empty() {
            return Err(CustomError::error(
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                Context::line(0, line, location.start, 1),
            ));
        }
        let mut peptide = Self::default();
        let mut ambiguous_lookup = Vec::new();
        let mut cross_link_lookup = Vec::new();
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
                        &mut cross_link_lookup,
                        custom_database,
                    )
                    .map(|m| {
                        m.defined().ok_or_else(|| {
                            CustomError::error(
                                "Invalid modification",
                                "A modification in the sloppy peptide format cannot be ambiguous",
                                Context::line(
                                    0,
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
                        .map(Modification::Simple)
                    })
                    .flat_err()?;
                    index = end_index + 1;

                    match peptide.sequence.last_mut() {
                        Some(aa) => aa.modifications.push(modification),
                        None => {
                            peptide.n_term = Some(
                                modification
                                    .simple()
                                    .expect("Can only put a simple modification on an N terminus.")
                                    .clone(),
                            );
                        }
                    }
                }
                ch => {
                    peptide.sequence.push(SequenceElement::new(
                        ch.try_into().map_err(|()| {
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

    fn sloppy_modification_internal(line: &str) -> Option<SimpleModification> {
        Self::sloppy_modification(line, 0..line.len(), None)
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
            Some(Ontology::Unimod.find_name("deamidated", None).unwrap())
        );
    }
}
