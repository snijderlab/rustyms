use std::sync::OnceLock;

use crate::{
    error::CustomError,
    helper_functions::*,
    identification::PeaksFamilyId,
    ontologies::CustomDatabase,
    peptide::{SemiAmbiguous, SloppyParsingParameters},
    system::{usize::Charge, MassOverCharge},
    LinearPeptide,
};

use serde::{Deserialize, Serialize};

use super::{
    common_parser::{Location, OptionalLocation},
    csv::{parse_csv, CsvLine},
    modification::Ontology,
    AminoAcid, BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid DeepNovoFamily line",
    "This column is not a number but it is required to be a number in this format",
);
static ID_ERROR: (&str, &str) =  (
    "Invalid DeepNovoFamily line",
    "This column is not a valid ID but it is required to be in this peaks format\nExamples of valid IDs: '1234' & 'F2:1234'");

static PARAMETERS_LOCK: OnceLock<SloppyParsingParameters> = OnceLock::new();

format_family!(
    /// The format for any DeepNovoFamily file
    DeepNovoFamilyFormat,
    /// The data from any DeepNovoFamily file
    DeepNovoFamilyData,
    DeepNovoFamilyVersion, [&DEEPNOVO_V0_0_1, &PGPOINTNOVO_V1_0_6], b'\t';
    required {
        scan: Vec<PeaksFamilyId>, |location: Location, _| location.or_empty()
            .map_or(Ok(Vec::new()), |l| l.array(';').map(|v| v.parse(ID_ERROR)).collect::<Result<Vec<_>,_>>());

        peptide: Option<LinearPeptide<SemiAmbiguous>>, |location: Location, custom_database: Option<&CustomDatabase>|
                location.or_empty().map(|location| LinearPeptide::sloppy_pro_forma(
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
                )).transpose();
        score: Option<f64>, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        local_confidence: Option<Vec<f64>>, |location: Location, _| location.or_empty()
                .optional_array(',').map(|array| array.map(|l| l.parse::<f64>(NUMBER_ERROR).map(|v| 2.0 / (1.0 + (-v).exp()))).collect::<Result<Vec<_>, _>>())
                .transpose();
    }
    optional {
        z: Charge, |location: Location, _| location
            .trim_end_matches(".0")
            .parse::<usize>(NUMBER_ERROR)
            .map(Charge::new::<crate::system::e>);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
    }
);

impl From<DeepNovoFamilyData> for IdentifiedPeptide {
    fn from(mut value: DeepNovoFamilyData) -> Self {
        if value.local_confidence.as_ref().map(Vec::len)
            != value.peptide.as_ref().map(LinearPeptide::len)
        {
            value.local_confidence = value.local_confidence.map(interpolate_lc);
        }

        Self {
            score: value.score.map(|score| (2.0 / (1.0 + (-score).exp()))),
            metadata: MetaData::DeepNovoFamily(value),
        }
    }
}

/// Interpolate the local confidence when the confidence between AAs is used instead of the confidence of a single AA
#[allow(clippy::needless_pass_by_value)] // The return value will replace the given value, so moving is fine
fn interpolate_lc(local_confidence: Vec<f64>) -> Vec<f64> {
    let mut reinterpolated = Vec::with_capacity(local_confidence.len() + 1);

    for i in 0..local_confidence.len() {
        if i == 0 {
            reinterpolated.push(local_confidence[i]);
        } else {
            let average = (local_confidence[i - 1] + local_confidence[i]) / 2.0;
            reinterpolated.push(average);
        }
    }
    reinterpolated.push(local_confidence[local_confidence.len() - 1]);

    reinterpolated
}

/// The only known version of DeepNovo
pub const DEEPNOVO_V0_0_1: DeepNovoFamilyFormat = DeepNovoFamilyFormat {
    version: DeepNovoFamilyVersion::DeepNovoV0_0_1,
    scan: "scan",
    peptide: "predicted_sequence",
    score: "predicted_score",
    local_confidence: "predicted_position_score",
    mz: None,
    z: None,
};

/// The only known version of PGPointNovo
pub const PGPOINTNOVO_V1_0_6: DeepNovoFamilyFormat = DeepNovoFamilyFormat {
    version: DeepNovoFamilyVersion::PGPointNovoV1_0_6,
    scan: "scan_list_original",
    peptide: "predicted_sequence",
    score: "predicted_score",
    local_confidence: "predicted_position_score",
    mz: Some("precursor_mz"),
    z: Some("precursor_charge"),
};

/// All possible DeepNovo versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum DeepNovoFamilyVersion {
    #[default]
    /// DeepNovo version 0.0.1
    DeepNovoV0_0_1,
    /// PGPointNovo version 1.0.6
    PGPointNovoV1_0_6,
}

impl std::fmt::Display for DeepNovoFamilyVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::DeepNovoV0_0_1 => "DeepNovo v0.0.1",
                Self::PGPointNovoV1_0_6 => "PGPointNovo v1.0.6",
            }
        )
    }
}
