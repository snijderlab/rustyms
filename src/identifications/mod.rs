//! Read in the annotations from peptide identification sources

mod common_parser;
#[path = "../shared/csv.rs"]
mod csv;
mod fasta;
mod identified_peptide;
mod novor;
mod peaks;

use self::csv::CsvLine;
use crate::error::CustomError;
pub use fasta::*;
pub use identified_peptide::*;
pub use novor::*;
pub use peaks::*;

#[cfg(test)]
mod novor_tests;
#[cfg(test)]
mod peaks_tests;

use std::str::FromStr;

impl CsvLine {
    pub fn column_context(&self, column: usize) -> crate::error::Context {
        crate::error::Context::line(
            self.line_index,
            self.line.clone(),
            self.fields[column].start,
            self.fields[column].len(),
        )
    }

    pub fn range_context(&self, range: std::ops::Range<usize>) -> crate::error::Context {
        crate::error::Context::line(self.line_index, self.line.clone(), range.start, range.len())
    }

    pub fn full_context(&self) -> crate::error::Context {
        crate::error::Context::full_line(self.line_index, self.line.clone())
    }
    /// Parse a column into the given format, if erroneous extend the base error with the correct context and return that
    pub fn parse_column<F: FromStr>(
        &self,
        column: usize,
        base_error: &CustomError,
    ) -> Result<F, CustomError> {
        self[column]
            .parse()
            .map_err(|_| base_error.with_context(self.column_context(column)))
    }
    /// Parse a column into the given format, if erroneous extend the base error with the correct context and return that
    pub fn parse_column_or_empty<F: FromStr>(
        &self,
        column: usize,
        base_error: &CustomError,
    ) -> Result<Option<F>, CustomError> {
        let text = &self[column];
        if text.is_empty() || text == "-" {
            Ok(None)
        } else {
            Ok(Some(text.parse().map_err(|_| {
                base_error.with_context(self.column_context(column))
            })?))
        }
    }
}

/// The required methods for any source of identified peptides
pub trait IdentifiedPeptideSource
where
    Self: std::marker::Sized,
{
    /// The source data where the peptides are parsed form eg [`CsvLine`]
    type Source;
    /// The format type eg [`PeaksFormat`]
    type Format;
    /// Parse a single identified peptide from its source and return the detected format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse(source: &Self::Source) -> Result<(Self, &'static Self::Format), CustomError>;
    /// Parse a single identified peptide with the given format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_specific(
        source: &Self::Source,
        format: &'static Self::Format,
    ) -> Result<Self, CustomError>;
}
