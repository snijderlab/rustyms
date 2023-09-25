use super::peaks::PeaksData;
use crate::{error::CustomError, ComplexPeptide};

/// A peptide that is identified by a de novo or database matching program
pub struct IdentifiedPeptide {
    peptide: ComplexPeptide,
    local_confidence: Option<Vec<f64>>,
    score: Option<f64>,
    metadata: MetaData,
}

/// The definition of all special metadata for all types of identified peptides that can be read
pub enum MetaData {
    /// Peaks metadata
    Peaks(PeaksData),
}

impl TryFrom<PeaksData> for IdentifiedPeptide {
    type Error = CustomError;
    fn try_from(value: PeaksData) -> Result<Self, Self::Error> {
        Ok(Self {
            peptide: ComplexPeptide::pro_forma(&value.peptide)?,
            local_confidence: Some(
                value
                    .local_confidence
                    .iter()
                    .map(|v| *v as f64 / 100.0)
                    .collect(),
            ),
            score: Some(value.de_novo_score.unwrap_or(value.alc) as f64 / 100.0),
            metadata: MetaData::Peaks(value),
        })
    }
}
