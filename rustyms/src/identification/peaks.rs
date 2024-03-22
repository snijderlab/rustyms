use std::fmt::Display;

use crate::{
    error::CustomError,
    helper_functions::InvertResult,
    system::{Charge, Mass, MassOverCharge, Time},
    LinearPeptide,
};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use super::{
    common_parser::{Location, OptionalLocation},
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid Peaks line",
    "This column is not a number but it is required to be a number in this peaks format",
);
static ID_ERROR: (&str, &str) =  ("Invalid Peaks line",
    "This column is not a valid peaks ID but it is required to be in this peaks format\nExamples of valid IDs: '1234', 'F2:1234', 'F2:1234 12345'");

format_family!(
    /// The format for any Peaks file
    PeaksFormat,
    /// The data from any peaks file
    PeaksData,
    PeaksVersion, [&OLD, &X, &XPLUS, &AB, &XI], b',';
    required {
        scan: Vec<PeaksId>, |location: Location| location.or_empty()
                        .map_or(Ok(Vec::new()), |l| l.array(';').map(|v| v.parse(ID_ERROR)).collect::<Result<Vec<_>,_>>());
        peptide: LinearPeptide, |location: Location| LinearPeptide::sloppy_pro_forma(
                            location.full_line(),
                            location.location.clone(),
                        );
        tag_length: usize, |location: Location| location.parse(NUMBER_ERROR);
        /// Range [0-1]
        alc: f64, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(|f| f / 100.0);
        length: usize, |location: Location| location.parse(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        z: Charge, |location: Location| location.parse::<usize>(NUMBER_ERROR).map(|c| Charge::new::<crate::system::e>(c as f64));
        mass: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        rt: Time, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        area: Option<f64>, |location: Location| location.or_empty().parse(NUMBER_ERROR);
        ppm: f64, |location: Location| location.parse(NUMBER_ERROR);
        ptm: String, |location: Location| Ok(location.get_string());
        local_confidence: Vec<f64>, |location: Location| location
                .array(' ')
                .map(|l| l.parse::<f64>(NUMBER_ERROR).map(|v| v / 100.0))
                .collect::<Result<Vec<_>, _>>();
        tag: String, |location: Location| Ok(location.get_string());
        mode: String, |location: Location| Ok(location.get_string());
    }
    optional {
        fraction: usize, |location: Location| location.parse(NUMBER_ERROR).map(Some);
        source_file: String, |location: Location| Ok(Some(location.get_string()));
        feature: PeaksId, |location: Location| location.or_empty().parse(ID_ERROR);
        de_novo_score: f64, |location: Location| location
                .parse::<f64>(NUMBER_ERROR)
                .map(|f| f / 100.0);
        predicted_rt: Time, |location: Location| location.or_empty().parse::<f64>(NUMBER_ERROR).map(|o| o.map(Time::new::<crate::system::time::min>));
        accession: String, |location: Location|  Ok(Some(location.get_string()));
    }
);

impl From<PeaksData> for IdentifiedPeptide {
    fn from(value: PeaksData) -> Self {
        Self {
            peptide: value.peptide.clone(),
            local_confidence: Some(value.local_confidence.clone()),
            score: Some(value.de_novo_score.unwrap_or(value.alc)),
            metadata: MetaData::Peaks(value),
        }
    }
}

/// An older version of a PEAKS export
pub const OLD: PeaksFormat = PeaksFormat {
    version: PeaksVersion::Old,
    scan: "scan",
    peptide: "peptide",
    tag_length: "tag length",
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ppm: "ppm",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag (>=0%)",
    mode: "mode",
    length: "length",
    fraction: None,
    source_file: None,
    feature: None,
    de_novo_score: None,
    predicted_rt: None,
    accession: None,
};
/// Version X of PEAKS export (made for build 31 january 2019)
pub const X: PeaksFormat = PeaksFormat {
    version: PeaksVersion::X,
    scan: "scan",
    peptide: "peptide",
    tag_length: "tag length",
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ppm: "ppm",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag (>=0%)",
    mode: "mode",
    length: "length",
    fraction: Some("fraction"),
    source_file: Some("source file"),
    feature: Some("feature"),
    de_novo_score: None,
    predicted_rt: None,
    accession: None,
};
/// Version X+ of PEAKS export (made for build 20 november 2019)
pub const XPLUS: PeaksFormat = PeaksFormat {
    version: PeaksVersion::Xplus,
    scan: "scan",
    peptide: "peptide",
    tag_length: "tag length",
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ppm: "ppm",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag (>=0%)",
    mode: "mode",
    length: "length",
    fraction: Some("fraction"),
    source_file: Some("source file"),
    feature: Some("feature"),
    de_novo_score: Some("denovo score"),
    predicted_rt: Some("predict rt"),
    accession: None,
};
/// Version 11 of PEAKS export
pub const XI: PeaksFormat = PeaksFormat {
    version: PeaksVersion::XI,
    scan: "scan",
    peptide: "peptide",
    tag_length: "tag length",
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ppm: "ppm",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag(>=0.0%)",
    mode: "mode",
    length: "length",
    fraction: None,
    source_file: None,
    feature: Some("feature"),
    de_novo_score: None,
    predicted_rt: None,
    accession: None,
};
/// Version Ab of PEAKS export
pub const AB: PeaksFormat = PeaksFormat {
    version: PeaksVersion::Ab,
    scan: "scan",
    peptide: "peptide",
    tag_length: "tag length",
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ppm: "ppm",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag (>=0%)",
    mode: "mode",
    length: "length",
    fraction: None,
    source_file: None,
    feature: None,
    de_novo_score: None,
    predicted_rt: None,
    accession: Some("accession"),
};

/// All possible peaks versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
pub enum PeaksVersion {
    /// An older version of a PEAKS export
    Old,
    /// Version X of PEAKS export (made for build 31 january 2019)
    X,
    /// Version X+ of PEAKS export (made for build 20 november 2019)
    Xplus,
    /// Version Ab of PEAKS export
    Ab,
    /// Version 11
    #[default]
    XI,
}

impl std::fmt::Display for PeaksVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::Old => "Old",
                Self::X => "X",
                Self::Xplus => "X+",
                Self::Ab => "Ab",
                Self::XI => "11",
            }
        )
    }
}
/// The scans identifier for a peaks identification
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
pub struct PeaksId {
    /// The file, if defined
    pub file: Option<usize>,
    /// The scan(s)
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
