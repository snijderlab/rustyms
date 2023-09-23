use crate::ComplexPeptide;

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
    Peaks(),
}
