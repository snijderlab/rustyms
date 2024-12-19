use crate::{
    error::CustomError,
    identification::{
        common_parser::OptionalColumn, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
    },
    ontologies::CustomDatabase,
    system::{usize::Charge, MassOverCharge, Time},
    LinearPeptide, SemiAmbiguous, SloppyParsingParameters,
};

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use super::{
    common_parser::Location,
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid SpectrumSequenceList line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    /// The format for any SSL file
    SpectrumSequenceListFormat,
    /// The data from any SSL file
    SpectrumSequenceListData,
    SpectrumSequenceListVersion, [&CASCADIA_V0_0_5], b'\t', None;
    required {
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        scan: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        z: Charge, |location: Location, _| location
            .trim_end_matches(".0")
            .parse::<usize>(NUMBER_ERROR)
            .map(Charge::new::<crate::system::e>);
        peptide: LinearPeptide<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| LinearPeptide::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters::default(),
        );
        start_time: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        end_time: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
    }
    optional {
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        score_type: String, |location: Location, _| Ok(location.get_string());
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        adduct: String, |location: Location, _| Ok(location.get_string());
        precursormz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        moleculename: String, |location: Location, _| Ok(location.get_string());
        inchikey: String, |location: Location, _| Ok(location.get_string());
        otherkeys: String, |location: Location, _| Ok(location.get_string());
        ion_mobility: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        ion_mobility_units: String, |location: Location, _| Ok(location.get_string());
        ccs: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
    }
);

impl From<SpectrumSequenceListData> for IdentifiedPeptide {
    fn from(value: SpectrumSequenceListData) -> Self {
        Self {
            score: value.score,
            local_confidence: None,
            metadata: MetaData::SpectrumSequenceList(value),
        }
    }
}

/// The only known version of Cascadia
pub const CASCADIA_V0_0_5: SpectrumSequenceListFormat = SpectrumSequenceListFormat {
    version: SpectrumSequenceListVersion::Cascadia_V0_0_5,
    raw_file: "file",
    scan: "scan",
    z: "charge",
    peptide: "sequence",
    start_time: "start-time",
    end_time: "end-time",
    score: OptionalColumn::Required("score"),
    score_type: OptionalColumn::Required("score-type"),
    rt: OptionalColumn::Required("retention-time"),
    adduct: OptionalColumn::NotAvailable,
    precursormz: OptionalColumn::NotAvailable,
    moleculename: OptionalColumn::NotAvailable,
    inchikey: OptionalColumn::NotAvailable,
    otherkeys: OptionalColumn::NotAvailable,
    ion_mobility: OptionalColumn::NotAvailable,
    ion_mobility_units: OptionalColumn::NotAvailable,
    ccs: OptionalColumn::NotAvailable,
};

/// All possible SpectrumSequenceList versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum SpectrumSequenceListVersion {
    #[default]
    /// Cascadia version 0.0.5
    Cascadia_V0_0_5,
}

impl std::fmt::Display for SpectrumSequenceListVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::Cascadia_V0_0_5 => "Cascadia v0.0.5",
            }
        )
    }
}
