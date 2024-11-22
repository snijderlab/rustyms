use crate::{
    error::CustomError,
    identification::{
        common_parser::OptionalColumn, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
        SpectrumId,
    },
    ontologies::CustomDatabase,
    LinearPeptide, SemiAmbiguous, SloppyParsingParameters,
};

use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use super::{
    common_parser::Location,
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid PowerNovo line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    /// The format for any PowerNovo file
    PowerNovoFormat,
    /// The data from any PowerNovo file
    PowerNovoData,
    PowerNovoVersion, [&POWERNOVO_V1_0_1], b',', None;
    required {
        raw_file: PathBuf, |location: Location, _| Ok(location.get_string().into());
        peptide: LinearPeptide<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| LinearPeptide::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters::default(),
        );
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        local_confidence: Vec<f64>, |location: Location, _| location.array(' ')
            .map(|l| l.parse::<f64>(NUMBER_ERROR))
            .collect::<Result<Vec<_>, _>>();
    }
    optional {
        scan: SpectrumId, |location: Location, _| Ok(Some(SpectrumId::Native(location.get_string())));
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _custom_database: Option<&CustomDatabase>) -> Result<Self, CustomError> {
        if let Some(SpectrumId::Native(native)) = parsed.scan.as_ref() {
            if let Some(m) = IDENTIFER_REGEX
                .get_or_init(|| regex::Regex::new(r"([^/]+):index=(\d+)").unwrap())
                .captures(native)
            {
                parsed.raw_file = m.get(1).unwrap().as_str().into();
                parsed.scan = Some(SpectrumId::Index(
                    m.get(2).unwrap().as_str().parse::<usize>().unwrap(),
                ));
            }
        }
        Ok(parsed)
    }
);

/// The Regex to match against PowerNovo scan fields
static IDENTIFER_REGEX: std::sync::OnceLock<regex::Regex> = std::sync::OnceLock::new();

impl From<PowerNovoData> for IdentifiedPeptide {
    fn from(value: PowerNovoData) -> Self {
        Self {
            score: Some(value.score),
            local_confidence: Some(value.local_confidence.iter().map(|v| *v / 100.0).collect()),
            metadata: MetaData::PowerNovo(value),
        }
    }
}

/// The only known version of PowerNovo
pub const POWERNOVO_V1_0_1: PowerNovoFormat = PowerNovoFormat {
    version: PowerNovoVersion::V1_0_1,
    scan: OptionalColumn::NotAvailable,
    raw_file: "spectrum name",
    peptide: "powernovo peptides",
    score: "powernovo score",
    local_confidence: "powernovo aascore",
};

/// All possible PowerNovo versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum PowerNovoVersion {
    #[default]
    /// PowerNovo version 1.0.1
    V1_0_1,
}

impl std::fmt::Display for PowerNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::V1_0_1 => "v1.0.1",
            }
        )
    }
}
