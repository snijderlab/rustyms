use std::fmt::Display;

use itertools::Itertools;

use crate::{error::CustomError, helper_functions::InvertResult, ComplexPeptide, LinearPeptide};

use super::{
    common_parser::{Location, OptionalLocation},
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
};

/// The file format for any peaks format, determining the existence and location of all possible columns
#[derive(Debug, Clone)]
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
    version: PeaksVersion,
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
#[derive(Clone, Debug)]
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
    version: PeaksVersion::Old,
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
    version: PeaksVersion::X,
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
    version: PeaksVersion::Xplus,
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
    version: PeaksVersion::Ab,
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
    version: PeaksVersion::XI,
};

/// A single parsed line of a peaks file
#[allow(missing_docs)]
#[derive(Debug)]
pub struct PeaksData {
    pub fraction: Option<usize>,
    pub source_file: Option<String>,
    pub feature: Option<PeaksId>,
    pub scan: Vec<PeaksId>,
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
    pub version: PeaksVersion,
}

/// The scans identifier for a peaks identification
#[derive(Debug)]
pub struct PeaksId {
    pub file: Option<usize>,
    pub scans: Vec<usize>,
}

impl Display for PeaksId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}",
            self.file.map_or(String::new(), |f| format!("F{f}:")),
            self.scans.iter().join(",")
        )
    }
}

impl std::str::FromStr for PeaksId {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((start, end)) = s.split_once(':') {
            Ok(Self {
                file: Some(start[1..].parse().map_err(|_| ())?),
                scans: end
                    .split(' ')
                    .map(str::parse)
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|_| ())?,
            })
        } else {
            Ok(Self {
                file: None,
                scans: vec![s.parse().map_err(|_| ())?],
            })
        }
    }
}

impl IdentifiedPeptideSource for PeaksData {
    type Source = CsvLine;
    type Format = PeaksFormat;
    fn parse(source: &CsvLine) -> Result<(Self, &'static PeaksFormat), CustomError> {
        for format in [&XI, &XPLUS, &X, &AB, &OLD] {
            if let Ok(peptide) = Self::parse_specific(source, format) {
                return Ok((peptide, format));
            }
        }
        Err(CustomError::error(
            "Invalid Peaks line",
            "The correct format could not be determined automatically",
            source.full_context(),
        ))
    }
    fn parse_specific(source: &CsvLine, format: &PeaksFormat) -> Result<Self, CustomError> {
        if source.fields.len() != format.number_of_columns() {
            return Err(CustomError::error(
                "Invalid Peaks line", 
                format!("The number of columns ({}) is not equal to the expected number of columns ({})", source.fields.len(), format.number_of_columns()), 
                source.full_context()));
        }
        let number_error = CustomError::error(
            "Invalid Peaks line",
            format!("This column is not a number but it is required to be a number in this peaks format ({})", format.version),
            source.full_context(),
        );
        let id_error = CustomError::error(
            "Invalid Peaks line",
            format!("This column is not a valid peaks ID but it is required to be in this peaks format ({})\nExamples of valid IDs: '1234', 'F2:1234', 'F2:1234 12345'", format.version),
            source.full_context(),
        );
        Ok(Self {
            fraction: Location::optional_column(format.fraction, source).parse(&number_error)?,
            source_file: Location::optional_column(format.source_file, source).get_string(),
            feature: Location::optional_column(format.feature, source)
                .or_empty()
                .parse(&id_error)?,
            scan: Location::column(format.scan, source)
                .or_empty()
                .map(|l| l.array(';').map(|v| v.parse(&id_error)).collect())
                .invert()?
                .unwrap_or_default(),
            peptide: ComplexPeptide::sloppy_pro_forma(
                &source.line,
                source.fields[format.peptide].clone(),
            )?,
            tag_length: Location::column(format.tag_length, source).parse(&number_error)?,
            de_novo_score: Location::optional_column(format.de_novo_score, source)
                .parse(&number_error)?,
            alc: Location::column(format.alc, source).parse::<f64>(&number_error)? / 100.0,
            length: Location::column(format.length, source).parse(&number_error)?,
            mz: Location::column(format.mz, source).parse(&number_error)?,
            z: Location::column(format.z, source).parse(&number_error)?,
            rt: Location::column(format.rt, source).parse(&number_error)?,
            predicted_rt: Location::optional_column(format.de_novo_score, source)
                .or_empty()
                .parse(&number_error)?,
            area: Location::column(format.mz, source)
                .or_empty()
                .parse(&number_error)?,
            mass: Location::column(format.mass, source).parse(&number_error)?,
            ppm: Location::column(format.ppm, source).parse(&number_error)?,
            ptm: source[format.ptm].to_string(),
            local_confidence: Location::column(format.local_confidence, source)
                .array(' ')
                .map(|l| l.parse::<f64>(&number_error).map(|v| v / 100.0))
                .collect::<Result<Vec<_>, _>>()?,
            tag: source[format.tag].to_string(),
            mode: source[format.mode].to_string(),
            accession: Location::optional_column(format.accession, source).get_string(),
            version: format.version.clone(),
        })
    }
    fn parse_file(
        path: impl AsRef<std::path::Path>,
    ) -> Result<BoxedIdentifiedPeptideIter<Self>, String> {
        parse_csv(path, b',').map(|lines| {
            Self::parse_many::<Box<dyn Iterator<Item = Self::Source>>>(Box::new(
                lines.skip(1).map(Result::unwrap),
            ))
        })
    }
}

impl From<PeaksData> for IdentifiedPeptide {
    fn from(value: PeaksData) -> Self {
        Self {
            peptide: value.peptide.clone(),
            local_confidence: Some(value.local_confidence.clone()),
            score: Some(value.de_novo_score.map_or(value.alc, |s| s as f64 / 100.0)),
            metadata: MetaData::Peaks(value),
        }
    }
}
