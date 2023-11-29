use crate::{error::CustomError, helper_functions::InvertResult};
use std::{ops::Range, str::FromStr};

use super::csv::CsvLine;

// Create:
// * XXFormat
// * XXData
// * XXParser?
macro_rules! format_family {
    ($format:ident, $data:ident, $version:ident, $versions:expr, $separator:expr; required {$($rname:ident: $rtyp:ty, $rf:expr);+; } optional { $($oname:ident: $otyp:ty, $of:expr);+;}) => {
        #[non_exhaustive]
        #[derive(Debug, PartialEq, Eq, Clone)]
        pub struct $format {
            $($rname: &'static str),+
            ,$($oname: Option<&'static str>),+
            ,version: $version,
        }

        #[non_exhaustive]
        #[derive(Debug, PartialEq, Clone)]
        pub struct $data {
            $(pub $rname: $rtyp),+
            ,$(pub $oname: Option<$otyp>),+
            ,version: $version,
        }

        impl IdentifiedPeptideSource for $data {
            type Source = CsvLine;
            type Format = $format;
            fn parse(source: &Self::Source) -> Result<(Self, &'static Self::Format), CustomError> {
                for format in $versions {
                    if let Ok(peptide) = Self::parse_specific(source, format) {
                        return Ok((peptide, format));
                    }
                }
                Err(CustomError::error(
                    format!("Invalid {} line", stringify!($format)),
                    "The correct format could not be determined automatically",
                    source.full_context(),
                ))
            }
            fn parse_file(
                path: impl AsRef<std::path::Path>,
            ) -> Result<BoxedIdentifiedPeptideIter<Self>, String> {
                parse_csv(path, $separator, None).map(|lines| {
                    Self::parse_many::<Box<dyn Iterator<Item = Self::Source>>>(Box::new(
                        lines.skip(1).map(Result::unwrap),
                    ))
                })
            }
            fn parse_specific(source: &Self::Source, format: &$format) -> Result<Self, CustomError> {
                Ok(Self {
                    $($rname: $rf(source.column(format.$rname)?)?),+
                    ,$($oname: format.$oname.and_then(|column| source.column(column).ok().map(|c| $of(c))).invert()?),+
                    ,version: format.version.clone(),
                })
            }
        }
    };
}

macro_rules! file_format {
    ($format:ident, $version:expr; $($name:ident: $column:expr),+,) => {
        $format {
            $($name: $column),+
            ,version: $version,
        }
    };
}

impl CsvLine {
    pub fn column<'a>(&'a self, name: &str) -> Result<Location<'a>, CustomError> {
        self.index_column(name).map(|(_v, c)| Location {
            line: self,
            location: c.clone(),
        })
    }
}

/// The base location type to keep track of the location of to be parsed pieces in the monadic parser combinators below
pub struct Location<'a> {
    pub(super) line: &'a CsvLine,
    pub(super) location: Range<usize>,
}

impl<'a> Location<'a> {
    pub fn column(column: usize, source: &'a CsvLine) -> Self {
        Location {
            line: source,
            location: source.range(column).clone(),
        }
    }

    pub fn optional_column(column: Option<usize>, source: &'a CsvLine) -> Option<Self> {
        column.map(|index| Location {
            line: source,
            location: source.range(index).clone(),
        })
    }

    #[allow(clippy::needless_pass_by_value)]
    pub fn array(self, sep: char) -> std::vec::IntoIter<Location<'a>> {
        let mut offset = 0;
        let mut output = Vec::new();
        for part in self.line.line()[self.location.clone()].split(sep) {
            output.push(Location {
                line: self.line,
                location: self.location.start + offset..self.location.start + offset + part.len(),
            });
            offset += part.len() + 1;
        }
        output.into_iter()
    }

    pub fn or_empty(self) -> Option<Self> {
        let text = &self.line.line()[self.location.clone()];
        if text.is_empty() || text == "-" {
            None
        } else {
            Some(self)
        }
    }

    pub fn ignore(self, pattern: &str) -> Option<Self> {
        let text = &self.line.line()[self.location.clone()];
        if text == pattern {
            None
        } else {
            Some(self)
        }
    }

    pub fn parse<T: FromStr>(self, base_error: (&str, &str)) -> Result<T, CustomError> {
        self.line.line()[self.location.clone()]
            .trim()
            .parse()
            .map_err(|_| {
                CustomError::error(
                    base_error.0,
                    base_error.1,
                    self.line.range_context(self.location),
                )
            })
    }

    pub fn parse_with<T>(
        self,
        f: impl Fn(Self) -> Result<T, CustomError>,
    ) -> Result<T, CustomError> {
        f(self)
    }

    pub fn get_id(self, base_error: (&str, &str)) -> Result<(Option<usize>, usize), CustomError> {
        if let Some((start, end)) = self.line.line()[self.location.clone()].split_once(':') {
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
        self.line.line()[self.location].to_string()
    }

    pub fn as_str(&self) -> &str {
        &self.line.line()[self.location.clone()]
    }

    pub fn full_line(&self) -> &str {
        self.line.line()
    }

    pub fn apply(self, f: impl FnOnce(Self) -> Self) -> Self {
        f(self)
    }
}

pub trait OptionalLocation<'a> {
    fn or_empty(self) -> Option<Location<'a>>;
    fn parse<T: FromStr>(self, base_error: (&str, &str)) -> Result<Option<T>, CustomError>;
    fn parse_with<T>(
        self,
        f: impl Fn(Location<'a>) -> Result<T, CustomError>,
    ) -> Result<Option<T>, CustomError>;
    fn get_id(
        self,
        base_error: (&str, &str),
    ) -> Result<Option<(Option<usize>, usize)>, CustomError>;
    fn get_string(self) -> Option<String>;
    fn apply(self, f: impl FnOnce(Location<'a>) -> Location<'a>) -> Option<Location<'a>>;
    type ArrayIter: Iterator<Item = Location<'a>>;
    fn array(self, sep: char) -> Self::ArrayIter;
}

impl<'a> OptionalLocation<'a> for Option<Location<'a>> {
    fn or_empty(self) -> Self {
        self.and_then(Location::or_empty)
    }
    fn parse<T: FromStr>(self, base_error: (&str, &str)) -> Result<Option<T>, CustomError> {
        self.map(|l| l.parse::<T>(base_error)).invert()
    }
    fn parse_with<T>(
        self,
        f: impl Fn(Location<'a>) -> Result<T, CustomError>,
    ) -> Result<Option<T>, CustomError> {
        self.map(f).invert()
    }
    fn get_id(
        self,
        base_error: (&str, &str),
    ) -> Result<Option<(Option<usize>, usize)>, CustomError> {
        self.map(|l| l.get_id(base_error)).invert()
    }
    fn get_string(self) -> Option<String> {
        self.map(Location::get_string)
    }
    fn apply(self, f: impl FnOnce(Location<'a>) -> Location<'a>) -> Self {
        self.map(|s| s.apply(f))
    }
    type ArrayIter = std::vec::IntoIter<Location<'a>>;
    fn array(self, sep: char) -> Self::ArrayIter {
        self.map(|l| l.array(sep)).unwrap_or_default()
    }
}
