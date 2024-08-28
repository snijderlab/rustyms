use std::sync::OnceLock;

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
                    let modification = Modification::sloppy_modification(
                        line,
                        location.start + index + 1..location.start + end_index,
                        peptide.sequence.last(),
                        custom_database,
                    )
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

static SLOPPY_MOD_OPAIR_REGEX: OnceLock<Regex> = OnceLock::new();
static SLOPPY_MOD_ON_REGEX: OnceLock<Regex> = OnceLock::new();
static SLOPPY_MOD_NUMERIC_END_REGEX: OnceLock<Regex> = OnceLock::new();

impl Modification {
    /// Parse a modification defined by sloppy names
    /// # Errors
    /// If the name is not in Unimod, PSI-MOD, the custom database, or the predefined list of common trivial names.
    /// Or if this is the case when the modification follows a known structure (eg `mod (AAs)`).
    #[allow(clippy::missing_panics_doc)]
    pub fn sloppy_modification(
        line: &str,
        location: std::ops::Range<usize>,
        position: Option<&SequenceElement>,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<SimpleModification, CustomError> {
        let full_context = Context::line(None, line, location.start, location.len());
        let name = &line[location];

        Self::find_name(name, position, custom_database)
            .or_else( || {
                match name.trim().to_lowercase().split_once(':') {
                    Some(("u", tail)) => Ontology::Unimod.find_name(tail, None),
                    Some(("m", tail)) => Ontology::Psimod.find_name(tail, None),
                    Some(("c", tail)) => Ontology::Custom.find_name(tail, custom_database),
                    _ => None
                }
            })
            .or_else(|| {SLOPPY_MOD_OPAIR_REGEX.get_or_init(|| {Regex::new(r"(?:[^:]+:)?(.*) (?:(?:on)|(?:from)) ([A-Z])").unwrap()})
                .captures(name)
                .and_then(|capture| {
                    let pos = capture[2].chars().next().and_then(|a| AminoAcid::try_from(a).ok().map(|a| SequenceElement::new(a, None)));
                    Self::find_name(&capture[1], position.or(pos.as_ref()), custom_database)
                        .ok_or_else(|| {
                            parse_named_counter(
                                &capture[1].to_ascii_lowercase(),
                                glycan_parse_list(),
                                false,
                            )
                            .map(SimpleModification::Glycan)
                        })
                        .flat_err()
                        .ok()
                })
                .or_else(|| {
                    // Common sloppy naming: `modification (AAs)` also accepts `modification (Protein N-term)`
                    SLOPPY_MOD_ON_REGEX.get_or_init(|| {Regex::new(r"(.*)\s*\([- @a-zA-Z]+\)").unwrap()})
                        .captures(name)
                        .and_then(|capture| {
                            Self::find_name(&capture[1], position, custom_database)
                        })
                })
                .or_else(|| {
                    // Common sloppy naming: `modification1`
                    SLOPPY_MOD_NUMERIC_END_REGEX.get_or_init(|| {Regex::new(r"(.*)\d+").unwrap()})
                        .captures(name)
                        .and_then(|capture| {
                            Self::find_name(&capture[1], position, custom_database)
                        })
                })
            }).ok_or_else(|| {
                CustomError::error(
                    "Could not interpret modification",
                    "Modifications have to be defined as a number, Unimod, or PSI-MOD name, if this is a custom modification make sure to add it to the database",
                    full_context,
                ).with_suggestions(
                    Ontology::find_closest_many(
                        &[Ontology::Unimod, Ontology::Psimod],
                        &name.trim().to_lowercase(),
                        custom_database).suggestions)
            })
    }

    fn find_name(
        name: &str,
        position: Option<&SequenceElement>,
        custom_database: Option<&CustomDatabase>,
    ) -> Option<SimpleModification> {
        let name = name.trim().to_lowercase();
        match name.as_str() {
            "o" => Ontology::Unimod.find_id(35, None), // oxidation
            "cam" | "carbamidomethylation" => Ontology::Unimod.find_id(4, None), // carbamidomethyl
            "nem" => Ontology::Unimod.find_id(108, None), // Nethylmaleimide
            "deamidation" => Ontology::Unimod.find_id(7, None), // deamidation
            "pyro-glu" => Ontology::Unimod.find_id(
                if position.is_some_and(|p| p.aminoacid == AminoAcid::E) {
                    27
                } else {
                    28
                },
                None,
            ), // pyro Glu with the logic to pick the correct modification based on the amino acid it is placed on
            _ => crate::peptide::parse_modification::numerical_mod(&name)
                .ok()
                .or_else(|| Ontology::Unimod.find_name(&name, custom_database))
                .or_else(|| Ontology::Psimod.find_name(&name, custom_database))
                .or_else(|| Ontology::Custom.find_name(&name, custom_database)),
        }
    }
}
