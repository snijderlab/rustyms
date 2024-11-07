use std::sync::OnceLock;

use crate::{
    error::CustomError,
    ontologies::CustomDatabase,
    peptide::{SemiAmbiguous, SloppyParsingParameters},
    LinearPeptide,
};
use serde::{Deserialize, Serialize};

use super::{
    common_parser::Location,
    csv::{parse_csv, CsvLine},
    modification::Ontology,
    AminoAcid, BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid DeepNovo line",
    "This column is not a number but it is required to be a number in this DeepNovo format",
);

static PARAMETERS_LOCK: OnceLock<SloppyParsingParameters> = OnceLock::new();

format_family!(
    /// The format for any DeepNovo file
    DeepNovoFormat,
    /// The data from any DeepNovo file
    DeepNovoData,
    DeepNovoVersion, [&V0_0_1], b'\t';
    required {
        scan: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide: LinearPeptide<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| LinearPeptide::sloppy_pro_forma(
                            location.full_line(),
                            location.location.clone(),
                            custom_database,
                            PARAMETERS_LOCK.get_or_init(|| SloppyParsingParameters{
                                mod_indications: vec![
                                    (AminoAcid::Asparagine, Ontology::Unimod.find_id(7, None).unwrap()),
                                    (AminoAcid::Glutamine, Ontology::Unimod.find_id(7, None).unwrap()),
                                    (AminoAcid::Cysteine, Ontology::Unimod.find_id(6, None).unwrap()),
                                    (AminoAcid::Methionine, Ontology::Unimod.find_id(35, None).unwrap()),
                                ],
                                ..Default::default()
                            })
                        );
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        local_confidence: Vec<f64>, |location: Location, _| location
                .array(',')
                .map(|l| l.parse::<f64>(NUMBER_ERROR).map(|v| 2.0 / (1.0 + (-v).exp())))
                .collect::<Result<Vec<_>, _>>();
    }
    optional {}
);

impl From<DeepNovoData> for IdentifiedPeptide {
    fn from(value: DeepNovoData) -> Self {
        Self {
            score: Some(2.0 / (1.0 + (-value.score).exp())),
            metadata: MetaData::DeepNovo(value),
        }
    }
}

/// The only known version of DeepNovo
pub const V0_0_1: DeepNovoFormat = DeepNovoFormat {
    version: DeepNovoVersion::V0_0_1,
    scan: "scan",
    peptide: "predicted_sequence",
    score: "predicted_score",
    local_confidence: "predicted_position_score",
};

/// All possible deepnovo versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
pub enum DeepNovoVersion {
    /// Version 0.0.1
    #[default]
    V0_0_1,
}

impl std::fmt::Display for DeepNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::V0_0_1 => "v0.0.1",
            }
        )
    }
}
