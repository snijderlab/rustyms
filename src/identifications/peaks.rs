use crate::error::{Context, CustomError};
use std::ops::Range;

use super::csv::CsvLine;

/// The file format for any peaks format, determining the existence and location of all possible columns
pub struct PeaksFormat {
    fraction: Option<usize>,
    source_file: Option<usize>,
    feature: Option<usize>,
    scan: usize,
    peptide: usize,
    tag_length: usize,
    de_novo_score: Option<usize>,
    alc: usize,
    length: usize,
    mz: usize,
    z: usize,
    rt: usize,
    predicted_rt: Option<usize>,
    area: usize,
    mass: usize,
    ppm: usize,
    ptm: usize,
    local_confidence: usize,
    tag: usize,
    mode: usize,
    accession: Option<usize>,
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

/// A single parsed line of a peaks file
#[allow(missing_docs)]
pub struct PeaksData {
    pub fraction: Option<usize>,
    pub source_file: Option<String>,
    pub feature: Option<(usize, usize)>,
    pub scan: (usize, usize),
    pub peptide: String,
    pub tag_length: usize,
    pub de_novo_score: Option<usize>,
    pub alc: usize,
    pub length: usize,
    pub mz: usize,
    pub z: usize,
    pub rt: f64,
    pub predicted_rt: Option<f64>,
    pub area: f64,
    pub mass: f64,
    pub ppm: f64,
    pub ptm: String,
    pub local_confidence: Vec<usize>,
    pub tag: String,
    pub mode: String,
    pub accession: Option<String>,
}

impl PeaksData {
    /// Parse a single line into the given peaks format
    /// # Errors
    /// Returns Err when the parsing could not be performed successfully
    pub fn parse_read(line: &CsvLine, format: &PeaksFormat) -> Result<Self, CustomError> {
        //let (line_number, line, locations) = line;
        if line.fields.len() != format.number_of_columns() {
            return Err(CustomError::error(
                "Invalid Peaks line", 
                format!("The number of columns ({}) is not equal to the expected number of columns ({})", line.fields.len(), format.number_of_columns()), 
                line.full_context()));
        }
        // Unpack an optional column either taking the full cell, parsing the column based on context clues for the type, or applying a function to the cell
        macro_rules! optional_column {
            ($column:expr) => {
                $column.map(|index| line[index].to_string())
            };
            ($column:expr, $error:expr) => {
                match $column {
                    Some(index) => Some(line.parse_column(index, $error)?),
                    None => None,
                }
            };
            ($column:expr => $func:ident) => {
                match $column {
                    Some(index) => Some($func(line.fields[index].clone())?),
                    None => None,
                }
            };
        }
        let get_num_from_slice = |location: Range<usize>| -> Result<usize, CustomError> {
            line.line[location.clone()].parse::<usize>().map_err(|_| CustomError::error(
                "Invalid Peaks line",
                format!("This text is not a number but it is required to be a number in this peaks format ({})", format.format),
                Context::line(line.line_index, line.line.clone(), location.start, location.len()),
            ))
        };
        let get_frac = |location: Range<usize>| -> Result<(usize, usize), CustomError> {
            let (start, end) = line.line[location.clone()].split_once(':').ok_or_else(|| {
                CustomError::error(
                    "Invalid Peaks line",
                    "Invalid scan number, it is missing the required colon",
                    Context::line(
                        line.line_index,
                        line.line.clone(),
                        location.start,
                        location.len(),
                    ),
                )
            })?;
            Ok((
                get_num_from_slice(location.start + 1..start.len() - 1)?,
                get_num_from_slice(location.start + start.len() + 1..end.len())?,
            ))
        };
        let get_array = |location: Range<usize>| -> Result<Vec<usize>, CustomError> {
            let mut offset = 0;
            let mut output = Vec::new();
            for part in line.line[location.clone()].split(' ') {
                output.push(get_num_from_slice(location.start + offset..part.len())?);
                offset += part.len() + 1;
            }
            Ok(output)
        };
        let number_error = CustomError::error(
            "Invalid Peaks line",
            format!("This column is not a number but it is required to be a number in this peaks format ({})", format.format),
            line.full_context(),
        );
        Ok(Self {
            fraction: optional_column!(format.fraction, &number_error),
            source_file: optional_column!(format.source_file),
            feature: optional_column!(format.feature => get_frac),
            scan: get_frac(line.fields[format.scan].clone())?,
            peptide: line[format.peptide].to_string(),
            tag_length: line.parse_column(format.tag_length, &number_error)?,
            de_novo_score: optional_column!(format.de_novo_score, &number_error),
            alc: line.parse_column(format.alc, &number_error)?,
            length: line.parse_column(format.length, &number_error)?,
            mz: line.parse_column(format.mz, &number_error)?,
            z: line.parse_column(format.z, &number_error)?,
            rt: line.parse_column(format.rt, &number_error)?,
            predicted_rt: optional_column!(format.de_novo_score, &number_error),
            area: line.parse_column(format.area, &number_error)?,
            mass: line.parse_column(format.mass, &number_error)?,
            ppm: line.parse_column(format.ppm, &number_error)?,
            ptm: line[format.ptm].to_string(),
            local_confidence: get_array(line.fields[format.local_confidence].clone())?,
            tag: line[format.tag].to_string(),
            mode: line[format.mode].to_string(),
            accession: optional_column!(format.accession),
        })
    }
}
