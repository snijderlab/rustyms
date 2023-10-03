use crate::{
    error::{Context, CustomError},
    helper_functions::InvertResult,
    ComplexPeptide, LinearPeptide,
};
use std::{ops::Range, str::FromStr};

use super::csv::CsvLine;

/// The file format for any peaks format, determining the existence and location of all possible columns
pub struct PeaksFormat {
    fraction: Option<usize>,
    source_file: Option<usize>,
    feature: Option<usize>,
    de_novo_score: Option<usize>,
    predicted_rt: Option<usize>,
    accession: Option<usize>,
    scan: usize,
    peptide: usize,
    tag_length: usize,
    alc: usize,
    length: usize,
    mz: usize,
    z: usize,
    rt: usize,
    area: usize,
    mass: usize,
    ppm: usize,
    ptm: usize,
    local_confidence: usize,
    tag: usize,
    mode: usize,
    format: PeaksVersion,
}

impl PeaksFormat {
    const fn number_of_columns(&self) -> usize {
        15 + self.fraction.is_some() as usize
            + self.source_file.is_some() as usize
            + self.feature.is_some() as usize
            + self.de_novo_score.is_some() as usize
            + self.predicted_rt.is_some() as usize
            + self.accession.is_some() as usize
    }
}

/// All possible peaks versions
pub enum PeaksVersion {
    /// A custom version of a PEAKS file forma
    Custom,
    /// An older version of a PEAKS export
    Old,
    /// Version X of PEAKS export (made for build 31 january 2019)
    X,
    /// Version X+ of PEAKS export (made for build 20 november 2019)
    Xplus,
    /// Version Ab of PEAKS export
    Ab,
    /// Version 11
    XI,
}

impl std::fmt::Display for PeaksVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::Custom => "Custom",
                Self::Old => "Old",
                Self::X => "X",
                Self::Xplus => "X+",
                Self::Ab => "Ab",
                Self::XI => "11",
            }
        )
    }
}

/// An older version of a PEAKS export
pub const OLD: PeaksFormat = PeaksFormat {
    scan: 0,
    peptide: 1,
    tag_length: 2,
    alc: 3,
    length: 4,
    mz: 5,
    z: 6,
    rt: 7,
    area: 8,
    mass: 9,
    ppm: 10,
    ptm: 11,
    local_confidence: 12,
    tag: 13,
    mode: 14,
    fraction: None,
    source_file: None,
    feature: None,
    de_novo_score: None,
    predicted_rt: None,
    accession: None,
    format: PeaksVersion::Old,
};

/// Version X of PEAKS export (made for build 31 january 2019)
pub const X: PeaksFormat = PeaksFormat {
    fraction: Some(0),
    source_file: Some(1),
    feature: Some(2),
    peptide: 3,
    scan: 4,
    tag_length: 5,
    alc: 6,
    length: 7,
    mz: 8,
    z: 9,
    rt: 10,
    area: 11,
    mass: 12,
    ppm: 13,
    ptm: 14,
    local_confidence: 15,
    tag: 16,
    mode: 17,
    de_novo_score: None,
    predicted_rt: None,
    accession: None,
    format: PeaksVersion::X,
};

/// Version X+ of PEAKS export (made for build 20 november 2019)
pub const XPLUS: PeaksFormat = PeaksFormat {
    fraction: Some(0),
    source_file: Some(1),
    feature: Some(2),
    peptide: 3,
    scan: 4,
    tag_length: 5,
    de_novo_score: Some(6),
    alc: 7,
    length: 8,
    mz: 9,
    z: 10,
    rt: 11,
    predicted_rt: Some(12),
    area: 13,
    mass: 14,
    ppm: 15,
    ptm: 16,
    local_confidence: 17,
    tag: 18,
    mode: 19,
    accession: None,
    format: PeaksVersion::Xplus,
};

/// Version Ab of PEAKS export
pub const AB: PeaksFormat = PeaksFormat {
    scan: 0,
    peptide: 1,
    tag_length: 2,
    alc: 3,
    length: 4,
    mz: 5,
    z: 6,
    rt: 7,
    area: 8,
    mass: 9,
    ppm: 10,
    accession: Some(11),
    ptm: 12,
    local_confidence: 13,
    tag: 14,
    mode: 15,
    fraction: None,
    predicted_rt: None,
    de_novo_score: None,
    source_file: None,
    feature: None,
    format: PeaksVersion::Ab,
};

/// Version 11 of PEAKS export
pub const XI: PeaksFormat = PeaksFormat {
    source_file: Some(0),
    scan: 1,
    peptide: 2,
    tag_length: 3,
    alc: 4,
    length: 5,
    mz: 6,
    z: 7,
    rt: 8,
    area: 9,
    mass: 10,
    ppm: 11,
    ptm: 12,
    local_confidence: 13,
    mode: 14,
    tag: 15,
    feature: Some(16),
    fraction: None,
    accession: None,
    predicted_rt: None,
    de_novo_score: None,
    format: PeaksVersion::XI,
};

/// A single parsed line of a peaks file
#[allow(missing_docs)]
pub struct PeaksData {
    pub fraction: Option<usize>,
    pub source_file: Option<String>,
    pub feature: Option<(Option<usize>, usize)>,
    pub scan: Vec<(Option<usize>, usize)>,
    pub peptide: LinearPeptide,
    pub tag_length: usize,
    pub de_novo_score: Option<usize>,
    /// Average local confidence [0-1]
    pub alc: f64,
    pub length: usize,
    pub mz: f64,
    pub z: usize,
    pub rt: f64,
    pub predicted_rt: Option<f64>,
    pub area: Option<f64>,
    pub mass: f64,
    pub ppm: f64,
    pub ptm: String,
    /// Local confidence [0-1]
    pub local_confidence: Vec<f64>,
    pub tag: String,
    pub mode: String,
    pub accession: Option<String>,
}

impl PeaksData {
    /// Parse a single line into the given peaks format
    /// # Errors
    /// Returns Err when the parsing could not be performed successfully
    pub fn parse_read(line: &CsvLine, format: &PeaksFormat) -> Result<Self, CustomError> {
        if line.fields.len() != format.number_of_columns() {
            return Err(CustomError::error(
                "Invalid Peaks line", 
                format!("The number of columns ({}) is not equal to the expected number of columns ({})", line.fields.len(), format.number_of_columns()), 
                line.full_context()));
        }
        let number_error = CustomError::error(
            "Invalid Peaks line",
            format!("This column is not a number but it is required to be a number in this peaks format ({})", format.format),
            line.full_context(),
        );
        Ok(Self {
            fraction: Location::optional_column(format.fraction, line).parse(&number_error)?,
            source_file: Location::optional_column(format.source_file, line).get_string(),
            feature: Location::optional_column(format.feature, line)
                .or_empty()
                .get_id(&number_error)?,
            scan: Location::column(format.scan, line)
                .or_empty()
                .map(|l| l.array(';').map(|v| v.get_id(&number_error)).collect())
                .invert()?
                .unwrap_or_default(),
            peptide: ComplexPeptide::sloppy_pro_forma(
                &line.line,
                line.fields[format.peptide].clone(),
            )?,
            tag_length: Location::column(format.tag_length, line).parse(&number_error)?,
            de_novo_score: Location::optional_column(format.de_novo_score, line)
                .parse(&number_error)?,
            alc: Location::column(format.alc, line).parse::<f64>(&number_error)? / 100.0,
            length: Location::column(format.length, line).parse(&number_error)?,
            mz: Location::column(format.mz, line).parse(&number_error)?,
            z: Location::column(format.z, line).parse(&number_error)?,
            rt: Location::column(format.rt, line).parse(&number_error)?,
            predicted_rt: Location::optional_column(format.de_novo_score, line)
                .or_empty()
                .parse(&number_error)?,
            area: Location::column(format.mz, line)
                .or_empty()
                .parse(&number_error)?,
            mass: Location::column(format.mass, line).parse(&number_error)?,
            ppm: Location::column(format.ppm, line).parse(&number_error)?,
            ptm: line[format.ptm].to_string(),
            local_confidence: Location::column(format.local_confidence, line)
                .array(' ')
                .map(|l| l.parse::<f64>(&number_error).map(|v| v / 100.0))
                .collect::<Result<Vec<_>, _>>()?,
            tag: line[format.tag].to_string(),
            mode: line[format.mode].to_string(),
            accession: Location::optional_column(format.accession, line).get_string(),
        })
    }
}

/// The base location type to keep track of the location of to be parsed pieces in the monadic parser combinators below
struct Location<'a> {
    line: &'a CsvLine,
    location: Range<usize>,
}

impl<'a> Location<'a> {
    fn column(column: usize, line: &'a CsvLine) -> Self {
        Location {
            line,
            location: line.fields[column].clone(),
        }
    }

    fn optional_column(column: Option<usize>, line: &'a CsvLine) -> Option<Self> {
        column.map(|index| Location {
            line,
            location: line.fields[index].clone(),
        })
    }

    #[allow(clippy::needless_pass_by_value)]
    fn array(self, sep: char) -> impl Iterator<Item = Location<'a>> {
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

    fn or_empty(self) -> Option<Self> {
        let text = &self.line.line[self.location.clone()];
        if text.is_empty() || text == "-" {
            None
        } else {
            Some(self)
        }
    }

    fn parse<T: FromStr>(self, base_error: &CustomError) -> Result<T, CustomError> {
        self.line.line[self.location.clone()]
            .parse()
            .map_err(|_| base_error.with_context(self.line.range_context(self.location)))
    }

    fn get_id(self, base_error: &CustomError) -> Result<(Option<usize>, usize), CustomError> {
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

    fn get_string(self) -> String {
        self.line.line[self.location].to_string()
    }
}

trait OptionalLocation<'a> {
    fn or_empty(self) -> Option<Location<'a>>;
    fn parse<T: FromStr>(self, base_error: &CustomError) -> Result<Option<T>, CustomError>;
    fn get_id(
        self,
        base_error: &CustomError,
    ) -> Result<Option<(Option<usize>, usize)>, CustomError>;
    fn get_string(self) -> Option<String>;
}

impl<'a> OptionalLocation<'a> for Option<Location<'a>> {
    fn or_empty(self) -> Self {
        self.and_then(Location::or_empty)
    }
    fn parse<T: FromStr>(self, base_error: &CustomError) -> Result<Option<T>, CustomError> {
        self.map(|l| l.parse(base_error)).invert()
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
}
