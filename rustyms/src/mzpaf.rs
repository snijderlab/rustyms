//! WIP: mzPAF parser
use std::ops::Range;

use crate::{
    error::{Context, CustomError},
    helper_functions::{explain_number_error, next_num, next_number},
    modification::{self, Ontology},
    system::{mz, MassOverCharge},
    Fragment, SequenceElement, Tolerance,
};

pub fn parse_mzpaf(line: &str) -> Result<Vec<Fragment>, CustomError> {
    Ok(Vec::new())
}

fn parse_annotation(line: &str, range: Range<usize>) -> Result<Fragment, CustomError> {
    if let Some(index) = line[range.clone()].chars().position(|c| c == '/') {
        let ion = parse_ion(line, range.start..range.start + index)?;
        let deviation = parse_deviation(line, range.end - index..range.end)?;
        Ok(Fragment::default())
    } else {
        Err(CustomError::error(
            "Invalid annotation",
            "A deviation, in m/z or ppm, needs to be present",
            Context::line_range(0, line, range),
        ))
    }
}

/// Parse a mzPAF ion.
/// # Errors
/// When the ion is not formatted correctly.
fn parse_ion(line: &str, range: Range<usize>) -> Result<(), CustomError> {
    match line[range.clone()].chars().next() {
        Some('?') => {
            let potential_ordinal = next_number::<false, usize>(line, range.add_start(1).unwrap());
            todo!();
        }
        Some('a' | 'b' | 'c' | 'x' | 'y' | 'z') => {
            let ordinal = next_number::<false, usize>(line, range.add_start(1).unwrap());
            todo!();
        }
        Some('I') => {
            let amino_acid = line[range.clone()].chars().skip(1).next();
            let modification = if let Some('[') = line[range.clone()].chars().skip(2).next() {
                let first = line[range.clone()].char_indices().skip(3).next().unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(3)
                    .take_while(|(_, c)| *c != ']')
                    .last()
                    .unwrap();
                Ontology::Unimod.find_name(
                    &line[range.clone()][first..last.0 + last.1.len_utf8()],
                    None,
                ) // TODO: error message if not found
            } else {
                None
            };
            todo!()
        }
        Some('m') => {
            let first_ordinal =
                next_number::<false, usize>(line, range.add_start(1).unwrap()).unwrap();
            assert!(
                line[range.clone()].chars().skip(first_ordinal.0).next() == Some(':'),
                "Needs to be separated by colon"
            );
            let second_ordinal = next_number::<false, usize>(
                line,
                range.add_start(2 + first_ordinal.0 as isize).unwrap(),
            );
            todo!();
        }
        Some('_') => {
            let name = if let Some('{') = line[range.clone()].chars().skip(1).next() {
                let first = line[range.clone()].char_indices().skip(2).next().unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(2)
                    .take_while(|(_, c)| *c != '}')
                    .last()
                    .unwrap();
                Some(&line[range.clone()][first..last.0 + last.1.len_utf8()])
            } else {
                None
            };
            todo!()
        }
        Some('p') => todo!(),
        Some('r') => {
            // Same name as neutral losses
            let name = if let Some('[') = line[range.clone()].chars().skip(1).next() {
                let first = line[range.clone()].char_indices().skip(2).next().unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(2)
                    .take_while(|(_, c)| *c != ']')
                    .last()
                    .unwrap();
                Some(&line[range.clone()][first..last.0 + last.1.len_utf8()])
            } else {
                None
            };
            todo!()
        }
        Some('f') => todo!(),
        Some('s') => todo!(),
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
}
