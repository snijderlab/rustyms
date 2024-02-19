use crate::error::CustomError;
use std::{ops::Range, str::FromStr};

use crate::csv::CsvLine;
use crate::helper_functions::InvertResult;

// Create:
// * XXFormat
// * XXData
// * XXParser?
macro_rules! format_family {
    (#[doc = $format_doc:expr]
     $format:ident,
     #[doc = $data_doc:expr]
     $data:ident,
     $version:ident, $versions:expr, $separator:expr;
     required { $($(#[doc = $rdoc:expr])? $rname:ident: $rtyp:ty, $rf:expr;)* }
     optional { $($(#[doc = $odoc:expr])? $oname:ident: $otyp:ty, $of:expr;)*}) => {
        use super::common_parser::HasLocation;
        #[non_exhaustive]
        #[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default, Serialize, Deserialize)]
        #[doc = $format_doc]
        pub struct $format {
            $($rname: &'static str,)*
            $($oname: Option<&'static str>,)*
            version: $version
        }

        #[non_exhaustive]
        #[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
        #[doc = $data_doc]
        #[allow(missing_docs)]
        pub struct $data {
            $($(#[doc = $rdoc])? pub $rname: $rtyp,)*
            $($(#[doc = $odoc])? pub $oname: Option<$otyp>,)*
            /// The version used to read in the data
            pub version: $version
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
                        lines.map(Result::unwrap),
                    ))
                })
            }
            #[allow(clippy::redundant_closure_call)] // Macro magic
            fn parse_specific(source: &Self::Source, format: &$format) -> Result<Self, CustomError> {
                Ok(Self {
                    $($rname: $rf(source.column(format.$rname)?)?,)*
                    $($oname: format.$oname.and_then(|column| source.column(column).ok().map(|c| $of(c))).invert()?,)*
                    version: format.version.clone()
                })
            }
        }
    };
}

pub(crate) trait HasLocation {
    fn column<'a>(&'a self, name: &str) -> Result<Location<'a>, CustomError>;
}

impl HasLocation for CsvLine {
    /// Get the specified column
    /// # Errors
    /// If the given column does not exist
    fn column<'a>(&'a self, name: &str) -> Result<Location<'a>, CustomError> {
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
    #[allow(dead_code)]
    pub fn column(column: usize, source: &'a CsvLine) -> Self {
        Location {
            line: source,
            location: source.range(column).clone(),
        }
    }

    #[allow(dead_code)]
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
        for part in self.as_str().split(sep) {
            output.push(Location {
                line: self.line,
                location: self.location.start + offset..self.location.start + offset + part.len(),
            });
            offset += part.len() + 1;
        }
        output.into_iter()
    }

    pub fn or_empty(self) -> Option<Self> {
        let text = self.as_str();
        if text.is_empty() || text == "-" {
            None
        } else {
            Some(self)
        }
    }

    pub fn ignore(self, pattern: &str) -> Option<Self> {
        let text = self.as_str();
        if text == pattern {
            None
        } else {
            Some(self)
        }
    }

    /// # Errors
    /// If the parse method fails. See [`FromStr::parse`].
    pub fn parse<T: FromStr>(self, base_error: (&str, &str)) -> Result<T, CustomError> {
        self.as_str().trim().parse().map_err(|_| {
            CustomError::error(
                base_error.0,
                base_error.1,
                self.line.range_context(self.location),
            )
        })
    }

    /// # Errors
    /// If the provided parse method fails.
    pub fn parse_with<T>(
        self,
        f: impl Fn(Self) -> Result<T, CustomError>,
    ) -> Result<T, CustomError> {
        f(self)
    }

    /// # Errors
    /// If the text could not be read as a valid id.
    pub fn get_id(self, base_error: (&str, &str)) -> Result<(Option<usize>, usize), CustomError> {
        if let Some((start, end)) = self.as_str().split_once(':') {
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
        self.as_str().to_string()
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
    /// # Errors
    /// If the parse method fails. See [`FromStr::parse`].
    fn parse<T: FromStr>(self, base_error: (&str, &str)) -> Result<Option<T>, CustomError>;
    /// # Errors
    /// If the provided parse method fails.
    fn parse_with<T>(
        self,
        f: impl Fn(Location<'a>) -> Result<T, CustomError>,
    ) -> Result<Option<T>, CustomError>;
    /// # Errors
    /// If the text could not be read as a valid id.
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
