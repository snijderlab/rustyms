use std::path::{Path, PathBuf};

use crate::{
    error::CustomError,
    helper_functions::InvertResult,
    identification::PeaksRelatedId,
    ontologies::CustomDatabase,
    peptide::{SemiAmbiguous, SloppyParsingParameters},
    system::{usize::Charge, Mass, MassOverCharge, Time},
    LinearPeptide,
};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use super::{
    common_parser::{Location, OptionalLocation},
    csv::{parse_csv, CsvLine},
    modification::SimpleModification,
    peptide::PeptideModificationSearch,
    BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData, Modification,
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
    PeaksVersion, [&V12, &V11, &X, &OLD, &XPLUS, &AB], b',';
    required {
        scan: Vec<PeaksRelatedId>, |location: Location, _| location.or_empty()
                        .map_or(Ok(Vec::new()), |l| l.array(';').map(|v| v.parse(ID_ERROR)).collect::<Result<Vec<_>,_>>());
        peptide: LinearPeptide<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| LinearPeptide::sloppy_pro_forma(
                            location.full_line(),
                            location.location.clone(),
                            custom_database,
                            &SloppyParsingParameters::default()
                        );
        alc: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(|f| f / 100.0);
        length: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        z: Charge, |location: Location, _| location.parse::<usize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        area: Option<f64>, |location: Location, _| location.or_empty().parse(NUMBER_ERROR);
        ptm: Vec<SimpleModification>, |location: Location, custom_database: Option<&CustomDatabase>|
            location.or_empty().array(';').map(|v| {
                let v = v.trim();
                Modification::sloppy_modification(v.full_line(), v.location.clone(), None, custom_database)
            }
            ).unique().collect::<Result<Vec<_>,_>>();
        local_confidence: Vec<f64>, |location: Location, _| location
                .array(' ')
                .map(|l| l.parse::<f64>(NUMBER_ERROR).map(|v| v / 100.0))
                .collect::<Result<Vec<_>, _>>();
        tag: String, |location: Location, _| Ok(location.get_string());
        mode: String, |location: Location, _| Ok(location.get_string());
    }
    optional {
        fraction: usize, |location: Location, _| location.parse(NUMBER_ERROR).map(Some);
        raw_file: PathBuf, |location: Location, _| Ok(Some(Path::new(&location.get_string()).to_owned()));
        feature: PeaksRelatedId, |location: Location, _| location.or_empty().parse(ID_ERROR);
        de_novo_score: f64, |location: Location, _| location
                .parse::<f64>(NUMBER_ERROR)
                .map(|f| f / 100.0);
        predicted_rt: Time, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR).map(|o| o.map(Time::new::<crate::system::time::min>));
        accession: String, |location: Location, _|  Ok(Some(location.get_string()));
    }
);

impl From<PeaksData> for IdentifiedPeptide {
    fn from(mut value: PeaksData) -> Self {
        // Add the meaningful modifications to replace mass modifications
        value.peptide = PeptideModificationSearch::in_modifications(value.ptm.clone())
            .tolerance(super::Tolerance::Absolute(super::system::da(0.05)))
            .search(value.peptide);

        Self {
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
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag (>=0%)",
    mode: "mode",
    length: "length",
    fraction: None,
    raw_file: None,
    feature: None,
    de_novo_score: None,
    predicted_rt: None,
    accession: None,
};
/// Version X of PEAKS export (made for build 31 January 2019)
pub const X: PeaksFormat = PeaksFormat {
    version: PeaksVersion::X,
    scan: "scan",
    peptide: "peptide",
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag (>=0%)",
    mode: "mode",
    length: "length",
    fraction: Some("fraction"),
    raw_file: Some("source file"),
    feature: Some("feature"),
    de_novo_score: None,
    predicted_rt: None,
    accession: None,
};
/// Version X+ of PEAKS export (made for build 20 November 2019)
pub const XPLUS: PeaksFormat = PeaksFormat {
    version: PeaksVersion::Xplus,
    scan: "scan",
    peptide: "peptide",
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag (>=0%)",
    mode: "mode",
    length: "length",
    fraction: Some("fraction"),
    raw_file: Some("source file"),
    feature: Some("feature"),
    de_novo_score: Some("denovo score"),
    predicted_rt: Some("predict rt"),
    accession: None,
};
/// Version 11 of PEAKS export
pub const V11: PeaksFormat = PeaksFormat {
    version: PeaksVersion::V11,
    scan: "scan",
    peptide: "peptide",
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag(>=0.0%)",
    mode: "mode",
    length: "length",
    fraction: None,
    raw_file: Some("source file"),
    feature: Some("feature id"),
    de_novo_score: None,
    predicted_rt: None,
    accession: None,
};
/// Version 12 of PEAKS export
pub const V12: PeaksFormat = PeaksFormat {
    version: PeaksVersion::V12,
    scan: "scan",
    peptide: "peptide",
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag(>=0%)",
    mode: "mode",
    length: "length",
    fraction: None,
    raw_file: Some("source file"),
    feature: Some("feature id"),
    de_novo_score: Some("deep novo score (%)"),
    predicted_rt: None,
    accession: None,
};
/// Version Ab of PEAKS export
pub const AB: PeaksFormat = PeaksFormat {
    version: PeaksVersion::Ab,
    scan: "scan",
    peptide: "peptide",
    alc: "alc (%)",
    mz: "m/z",
    z: "z",
    mass: "mass",
    rt: "rt",
    area: "area",
    ptm: "ptm",
    local_confidence: "local confidence (%)",
    tag: "tag (>=0%)",
    mode: "mode",
    length: "length",
    fraction: None,
    raw_file: None,
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
    /// Version X of PEAKS export (made for build 31 January 2019)
    X,
    /// Version X+ of PEAKS export (made for build 20 November 2019)
    Xplus,
    /// Version Ab of PEAKS export
    Ab,
    /// Version 11
    V11,
    /// Version 12
    #[default]
    V12,
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
                Self::V11 => "11",
                Self::V12 => "12",
            }
        )
    }
}
