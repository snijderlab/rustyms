use crate::{error::CustomError, helper_functions::InvertResult};
use std::{ops::Range, str::FromStr};

use super::csv::CsvLine;

/// The base location type to keep track of the location of to be parsed pieces in the monadic parser combinators below
pub struct Location<'a> {
    pub(super) line: &'a CsvLine,
    pub(super) location: Range<usize>,
}

impl<'a> Location<'a> {
    pub fn column(column: usize, source: &'a CsvLine) -> Self {
        Location {
            line: source,
            location: source.fields[column].clone(),
        }
    }

    pub fn optional_column(column: Option<usize>, source: &'a CsvLine) -> Option<Self> {
        column.map(|index| Location {
            line: source,
            location: source.fields[index].clone(),
        })
    }

    #[allow(clippy::needless_pass_by_value)]
    pub fn array(self, sep: char) -> impl Iterator<Item = Location<'a>> {
        let mut offset = 0;
        let mut output = Vec::new();
        for part in self.line.line[self.location.clone()].split(sep) {
            output.push(Location {
                line: self.line,
                location: self.location.start + offset..self.location.start + offset + part.len(),
            });
            offset += part.len() + 1;
        }
        output.into_iter()
    }

    pub fn or_empty(self) -> Option<Self> {
        let text = &self.line.line[self.location.clone()];
        if text.is_empty() || text == "-" {
            None
        } else {
            Some(self)
        }
    }

    pub fn parse<T: FromStr>(self, base_error: &CustomError) -> Result<T, CustomError> {
        self.line.line[self.location.clone()]
            .parse()
            .map_err(|_| base_error.with_context(self.line.range_context(self.location)))
    }

    pub fn parse_with<T>(
        self,
        f: impl Fn(Self) -> Result<T, CustomError>,
    ) -> Result<T, CustomError> {
        f(self)
    }

    pub fn get_id(self, base_error: &CustomError) -> Result<(Option<usize>, usize), CustomError> {
        if let Some((start, end)) = self.line.line[self.location.clone()].split_once(':') {
            Ok((
                Some(
                    Self {
                        line: self.line,
                        location: self.location.start + 1..self.location.start + start.len(),
                    }
                    .parse(base_error)?,
                ),
                Self {
                    line: self.line,
                    location: self.location.start + start.len() + 1
                        ..self.location.start + start.len() + 1 + end.len(),
                }
                .parse(base_error)?,
            ))
        } else {
            Ok((None, self.parse(base_error)?))
        }
    }

    pub fn get_string(self) -> String {
        self.line.line[self.location].to_string()
    }

    pub fn apply(self, f: impl FnOnce(Self) -> Self) -> Self {
        f(self)
    }
}

pub trait OptionalLocation<'a> {
    fn or_empty(self) -> Option<Location<'a>>;
    fn parse<T: FromStr>(self, base_error: &CustomError) -> Result<Option<T>, CustomError>;
    fn parse_with<T>(
        self,
        f: impl Fn(Location<'a>) -> Result<T, CustomError>,
    ) -> Result<Option<T>, CustomError>;
    fn get_id(
        self,
        base_error: &CustomError,
    ) -> Result<Option<(Option<usize>, usize)>, CustomError>;
    fn get_string(self) -> Option<String>;
    fn apply(self, f: impl FnOnce(Location<'a>) -> Location<'a>) -> Option<Location<'a>>;
}

impl<'a> OptionalLocation<'a> for Option<Location<'a>> {
    fn or_empty(self) -> Self {
        self.and_then(Location::or_empty)
    }
    fn parse<T: FromStr>(self, base_error: &CustomError) -> Result<Option<T>, CustomError> {
        self.map(|l| l.parse(base_error)).invert()
    }
    fn parse_with<T>(
        self,
        f: impl Fn(Location<'a>) -> Result<T, CustomError>,
    ) -> Result<Option<T>, CustomError> {
        self.map(f).invert()
    }
    fn get_id(
        self,
        base_error: &CustomError,
    ) -> Result<Option<(Option<usize>, usize)>, CustomError> {
        self.map(|l| l.get_id(base_error)).invert()
    }
    fn get_string(self) -> Option<String> {
        self.map(Location::get_string)
    }
    fn apply(self, f: impl FnOnce(Location<'a>) -> Location<'a>) -> Self {
        self.map(|s| s.apply(f))
    }
}
