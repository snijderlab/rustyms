//! Read in the annotations from peptide identification sources

#[path = "../shared/csv.rs"]
mod csv;
mod fasta;
mod identified_peptide;
mod peaks;
use std::str::FromStr;

pub use fasta::*;
pub use identified_peptide::*;
pub use peaks::*;

use crate::error::{Context, CustomError};

use self::csv::CsvLine;

impl CsvLine {
    pub fn column_context(&self, column: usize) -> crate::error::Context {
        crate::error::Context::line(
            self.line_index,
            self.line.clone(),
            self.fields[column].start,
            self.fields[column].len(),
        )
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
}
