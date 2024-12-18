use crate::{
    error::CustomError,
    identification::{
        common_parser::OptionalColumn, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
    },
    ontologies::CustomDatabase,
    system::{usize::Charge, Time},
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
    "Invalid Cascadia line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    /// The format for any Cascadia file
    CascadiaFormat,
    /// The data from any Cascadia file
    CascadiaData,
    CascadiaVersion, [&CASCADIA_V0_0_5], b'\t', None;
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
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
    }
    optional {
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
    }
);

impl From<CascadiaData> for IdentifiedPeptide {
    fn from(value: CascadiaData) -> Self {
        Self {
            score: Some(value.score),
            local_confidence: None,
            metadata: MetaData::Cascadia(value),
        }
    }
}

/// The only known version of Cascadia
pub const CASCADIA_V0_0_5: CascadiaFormat = CascadiaFormat {
    version: CascadiaVersion::V0_0_5,
    raw_file: "file",
    scan: "scan",
    z: "charge",
    peptide: "sequence",
    score: "score",
    rt: OptionalColumn::Optional("retention-time"),
};

/// All possible Cascadia versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum CascadiaVersion {
    #[default]
    /// Cascadia version 0.0.5
    V0_0_5,
}

impl std::fmt::Display for CascadiaVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::V0_0_5 => "v0.0.5",
            }
        )
    }
}
