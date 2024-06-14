//! WIP: mzPAF parser
use std::ops::{Range, RangeBounds};

use crate::{
    error::{Context, CustomError},
    helper_functions::{explain_number_error, next_number},
    modification::{Ontology, SimpleModification},
    system::{mz, MassOverCharge},
    AminoAcid, Fragment, MolecularFormula, NeutralLoss, Tolerance,
};

pub fn parse_mzpaf(line: &str) -> Result<Vec<Fragment>, CustomError> {
    Ok(Vec::new())
}

fn parse_annotation(line: &str, range: Range<usize>) -> Result<Fragment, CustomError> {
    if let Some(index) = line[range.clone()].chars().position(|c| c == '/') {
        // Parse analyte index
        let ion = parse_ion(line, range.start..range.start + index)?;
        // Parse neutral losses
        // Parse isotopes
        // Parse charge
        // Parse adduct type
        let deviation = parse_deviation(line, range.end - index..range.end)?;
        Ok(Fragment::default())
    } else {
        // TODO: not required the deviation
        Err(CustomError::error(
            "Invalid annotation",
            "A deviation, in m/z or ppm, needs to be present",
            Context::line_range(0, line, range),
        ))
    }
}

enum mzPAFIonType {
    Unknown(Option<usize>),
    MainSeries(char, usize),
    Immonium(AminoAcid, Option<SimpleModification>),
    Internal(usize, usize),
    Named(String),
    Precursor,
    Reporter(String), //TODO: store as preloaded name
    Formula(MolecularFormula),
}

/// Parse a mzPAF ion.
/// # Errors
/// When the ion is not formatted correctly.
fn parse_ion(line: &str, range: Range<usize>) -> Result<mzPAFIonType, CustomError> {
    match line[range.clone()].chars().next() {
        Some('?') => {
            if let Some(ordinal) = next_number::<false, usize>(line, range.add_start(1).unwrap()) {
                Ok(mzPAFIonType::Unknown(Some(ordinal.2.map_err(|err| {
                    CustomError::error(
                        "Invalid mzPAF unknown ion ordinal",
                        format!("The ordinal number {}", explain_number_error(&err)),
                        Context::line(0, line, range.start_index() + 1, ordinal.0),
                    )
                })?)))
            } else {
                Ok(mzPAFIonType::Unknown(None))
            }
        }
        Some(c @ ('a' | 'b' | 'c' | 'x' | 'y' | 'z')) => {
            if let Some(ordinal) = next_number::<false, usize>(line, range.add_start(1).unwrap()) {
                Ok(mzPAFIonType::MainSeries(
                    c,
                    ordinal.2.map_err(|err| {
                        CustomError::error(
                            "Invalid mzPAF unknown ion ordinal",
                            format!("The ordinal number {}", explain_number_error(&err)),
                            Context::line(0, line, range.start_index() + 1, ordinal.0),
                        )
                    })?,
                ))
            } else {
                Err(CustomError::error(
                    "Invalid mzPAF main series ion ordinal",
                    "For a main series ion the ordinal should be provided, like 'a12'",
                    Context::line(0, line, range.start_index(), 1),
                ))
            }
        }
        Some('I') => {
            let amino_acid = line[range.clone()].chars().nth(1).ok_or_else(|| {
                CustomError::error(
                    "Invalid mzPAF immonium",
                    "The source amino acid for this immonium ion should be present like 'IA'",
                    Context::line(0, line, range.start_index(), 1),
                )
            })?;
            let modification = if line[range.clone()].chars().nth(2) == Some('[') {
                let first = line[range.clone()].char_indices().nth(3).unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(3)
                    .take_while(|(_, c)| *c != ']')
                    .last()
                    .unwrap();
                Some(
                    Ontology::Unimod
                        .find_name(
                            &line[range.clone()][first..last.0 + last.1.len_utf8()],
                            None,
                        )
                        .ok_or_else(|| {
                            Ontology::Unimod.find_closest(
                                &line[range.clone()][first..last.0 + last.1.len_utf8()],
                                None,
                            )
                        })?,
                )
            } else {
                None
            };
            Ok(mzPAFIonType::Immonium(
                AminoAcid::try_from(amino_acid).map_err(|()| {
                    CustomError::error(
                        "Invalid mzPAF immonium ion",
                        "The provided amino acid is not a known amino acid",
                        Context::line(0, line, range.start_index() + 1, 1),
                    )
                })?,
                modification,
            ))
        }
        Some('m') => {
            let first_ordinal = next_number::<false, usize>(line, range.add_start(1).unwrap())
                .ok_or_else(|| {
                    CustomError::error(
                        "Invalid mzPAF internal ion first ordinal",
                        "The first ordinal for an internal ion should be present",
                        Context::line(0, line, range.start_index(), 1),
                    )
                })?;
            if line[range.clone()].chars().nth(first_ordinal.0) != Some(':') {
                return Err(CustomError::error(
                    "Invalid mzPAF internal ion ordinal separator",
                    "The internal ion ordinal separator should be a colon ':', like 'm4:6'",
                    Context::line(0, line, range.start_index() + 1 + first_ordinal.0, 1),
                ));
            }
            assert!(
                line[range.clone()].chars().nth(first_ordinal.0) == Some(':'),
                "Needs to be separated by colon"
            );
            let second_ordinal = next_number::<false, usize>(
                line,
                range.add_start(2 + first_ordinal.0 as isize).unwrap(),
            )
            .ok_or_else(|| {
                CustomError::error(
                    "Invalid mzPAF internal ion second ordinal",
                    "The second ordinal for an internal ion should be present",
                    Context::line(0, line, range.start_index() + 1 + first_ordinal.0, 1),
                )
            })?;
            let first_location = first_ordinal.2.map_err(|err| {
                CustomError::error(
                    "Invalid mzPAF internal ion first ordinal",
                    format!("The ordinal number {}", explain_number_error(&err)),
                    Context::line(0, line, range.start_index() + 1, first_ordinal.0),
                )
            })?;
            let second_location = second_ordinal.2.map_err(|err| {
                CustomError::error(
                    "Invalid mzPAF internal ion second ordinal",
                    format!("The ordinal number {}", explain_number_error(&err)),
                    Context::line(
                        0,
                        line,
                        range.start_index() + 2 + first_ordinal.0,
                        second_ordinal.0,
                    ),
                )
            })?;
            Ok(mzPAFIonType::Internal(first_location, second_location))
        }
        Some('_') => {
            // Format less strings
            // TODO: Potentially recognise the following as known contaminants:
            // 0@_{y1(R)}
            // 0@_{a2(LP)}
            // 0@_{b2(LP)}

            let name = if line[range.clone()].chars().nth(1) == Some('{') {
                let first = line[range.clone()].char_indices().nth(2).unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(2)
                    .take_while(|(_, c)| *c != '}')
                    .last()
                    .unwrap();
                Ok(&line[range][first..last.0 + last.1.len_utf8()])
            } else {
                Err(CustomError::error(
                    "Invalid mzPAF named compound",
                    "A named compound must be named with curly braces '{}' after the '_'",
                    Context::line(0, line, range.start_index(), 1),
                ))
            }?;
            Ok(mzPAFIonType::Named(name.to_string()))
        }
        Some('p') => Ok(mzPAFIonType::Precursor),
        Some('r') => {
            // Same name as neutral losses
            let name = if line[range.clone()].chars().nth(1) == Some('[') {
                let first = line[range.clone()].char_indices().nth(2).unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(2)
                    .take_while(|(_, c)| *c != ']')
                    .last()
                    .unwrap();
                Ok(&line[range][first..last.0 + last.1.len_utf8()])
            } else {
                Err(CustomError::error(
                    "Invalid mzPAF reporter ion",
                    "A reporter ion must be named with square braces '[]' after the 'r'",
                    Context::line(0, line, range.start_index(), 1),
                ))
            }?;
            Ok(mzPAFIonType::Reporter(name.to_string()))
        }
        Some('f') => {
            // Simple formula
            let formula_range = if line[range.clone()].chars().nth(1) == Some('{') {
                let first = line[range.clone()].char_indices().nth(2).unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(2)
                    .take_while(|(_, c)| *c != '}')
                    .last()
                    .unwrap();
                Ok(range.start_index() + first..range.start_index() + last.0 + last.1.len_utf8())
            } else {
                Err(CustomError::error(
                    "Invalid mzPAF formula",
                    "A formula must have the formula defined with curly braces '{}' after the 'f'",
                    Context::line(0, line, range.start_index(), 1),
                ))
            }?;
            let formula = MolecularFormula::from_mz_paf_inner(line, formula_range)?;

            Ok(mzPAFIonType::Formula(formula))
        }
        Some('s') => todo!(), // TODO: return as Formula
        Some(_) => Err(CustomError::error(
            "Invalid ion",
            "An ion cannot start with this character",
            Context::line(0, line, range.start, 1),
        )),
        None => Err(CustomError::error(
            "Invalid ion",
            "An ion cannot be an empty string",
            Context::line_range(0, line, range),
        )),
    }
}

fn parse_neutral_loss(line: &str, range: Range<usize>) -> Result<Vec<NeutralLoss>, CustomError> {
    let mut neutral_losses = Vec::new();
    while let Some('-' | '+') = line[range.clone()].chars().next() {
        //
    }
    Ok(neutral_losses)
}

/// Parse a mzPAF deviation, either a ppm or mz deviation.
/// # Errors
/// When the deviation is not '<number>' or '<number>ppm'.
fn parse_deviation(
    line: &str,
    range: Range<usize>,
) -> Result<Tolerance<MassOverCharge>, CustomError> {
    if let Some(ppm) = line[range.clone()].strip_suffix("ppm") {
        Ok(Tolerance::new_ppm(ppm.parse::<f64>().map_err(|err| {
            CustomError::error(
                "Invalid deviation",
                format!("ppm deviation should be a valid number {err}",),
                Context::line_range(0, line, range.start..range.end - 3),
            )
        })?))
    } else {
        Ok(Tolerance::new_absolute(MassOverCharge::new::<mz>(
            line[range.clone()].parse::<f64>().map_err(|err| {
                CustomError::error(
                    "Invalid deviation",
                    format!("m/z deviation should be a valid number {err}",),
                    Context::line_range(0, line, range),
                )
            })?,
        )))
    }
}

trait RangeExtension
where
    Self: Sized,
{
    fn add_start(&self, amount: isize) -> Option<Self>;
    fn add_end(&self, amount: isize) -> Option<Self>;
    fn start_index(&self) -> usize;
    fn end_index(&self) -> usize;
}

impl RangeExtension for Range<usize> {
    fn add_start(&self, amount: isize) -> Option<Self> {
        let new_start = usize::try_from(isize::try_from(self.start).ok()? + amount).ok()?;
        (new_start <= self.end).then_some(Self {
            start: new_start,
            end: self.end,
        })
    }
    fn add_end(&self, amount: isize) -> Option<Self> {
        let new_end = usize::try_from(isize::try_from(self.end).ok()? + amount).ok()?;
        (self.start <= new_end).then_some(Self {
            start: self.start,
            end: new_end,
        })
    }
    fn start_index(&self) -> usize {
        match self.start_bound() {
            std::ops::Bound::Unbounded => 0,
            std::ops::Bound::Included(s) => *s,
            std::ops::Bound::Excluded(s) => s + 1,
        }
    }
    fn end_index(&self) -> usize {
        match self.end_bound() {
            std::ops::Bound::Unbounded => 0,
            std::ops::Bound::Included(s) => *s,
            std::ops::Bound::Excluded(s) => s - 1,
        }
    }
}
