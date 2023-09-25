use crate::error::{Context, CustomError};
use std::ops::Range;

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

// TODO:
// * hard crashes if you have too little elements per line
// * does not crash if you have too many elements per line
impl PeaksData {
    /// Parse a single line into the given peaks format
    /// # Errors
    /// Returns Err when the parsing could not be performed successfully
    pub fn parse_read(
        line: &(usize, String, [Range<usize>]),
        format: &PeaksFormat,
    ) -> Result<Self, CustomError> {
        let (line_number, line, locations) = line;
        macro_rules! get_column {
            ($column:expr) => {
                match $column {
                    Some(index) => Some(line[locations[index].clone()].to_string()),
                    None => None,
                }
            };
            ($column:expr, $typ:ty) => {
                match $column {
                    Some(index) => Some(line[locations[index].clone()].parse::<$typ>().map_err(|_| CustomError::error(
                        "Invalid Peaks line",
                        format!("Column {} is not of the correct type ({}) for this peaks format ({})", index+1, stringify!($typ), format.format),
                        Context::line(*line_number, &line, locations[index].start, locations[index].len())
                    ))?),
                    None => None,
                }
            };
            ($column:expr => $func:ident) => {
                match $column {
                    Some(index) => Some($func(locations[index].clone())?),
                    None => None,
                }
            };
        }
        // TODO: refactor get_num into macro
        let get_num_from_slice = |location: Range<usize>| -> Result<usize, CustomError> {
            line[location.clone()].parse::<usize>().map_err(|_| CustomError::error(
                "Invalid Peaks line",
                format!("This text is not a number but it is required to be a number in this peaks format ({})", format.format),
                Context::line(*line_number, line, location.start, location.len()),
            ))
        };
        let get_num = |location: Range<usize>| -> Result<usize, CustomError> {
            line[location.clone()].parse::<usize>().map_err(|_| CustomError::error(
                "Invalid Peaks line",
                format!("This text is not a number but it is required to be a number in this peaks format ({})", format.format),
                Context::line(*line_number, line, location.start, location.len()),
            ))
        };
        let get_float = |location: Range<usize>| -> Result<f64, CustomError> {
            line[location.clone()].parse::<f64>().map_err(|_| CustomError::error(
                "Invalid Peaks line",
                format!("This text is not a number but it is required to be a number in this peaks format ({})", format.format),
                Context::line(*line_number, line, location.start, location.len()),
            ))
        };
        let get_frac = |location: Range<usize>| -> Result<(usize, usize), CustomError> {
            let (start, end) = line[location.clone()].split_once(':').ok_or_else(|| {
                CustomError::error(
                    "Invalid Peaks line",
                    "Invalid scan number it is missing the required colon",
                    Context::line(*line_number, line, location.start, location.len()),
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
            for part in line[location.clone()].split(' ') {
                output.push(get_num_from_slice(location.start + offset..part.len())?);
                offset += part.len() + 1;
            }
            Ok(output)
        };
        Ok(Self {
            fraction: get_column!(format.fraction, usize),
            source_file: get_column!(format.source_file),
            feature: get_column!(format.feature => get_frac),
            scan: get_frac(locations[format.scan].clone())?,
            peptide: line[locations[format.peptide].clone()].to_string(),
            tag_length: get_num(locations[format.tag_length].clone())?,
            de_novo_score: get_column!(format.de_novo_score, usize),
            alc: get_num(locations[format.alc].clone())?,
            length: get_num(locations[format.length].clone())?,
            mz: get_num(locations[format.mz].clone())?,
            z: get_num(locations[format.z].clone())?,
            rt: get_float(locations[format.rt].clone())?,
            predicted_rt: get_column!(format.de_novo_score, f64),
            area: get_float(locations[format.area].clone())?,
            mass: get_float(locations[format.mass].clone())?,
            ppm: get_float(locations[format.ppm].clone())?,
            ptm: line[locations[format.ptm].clone()].to_string(),
            local_confidence: get_array(locations[format.local_confidence].clone())?,
            tag: line[locations[format.tag].clone()].to_string(),
            mode: line[locations[format.mode].clone()].to_string(),
            accession: get_column!(format.accession),
        })
    }
}
