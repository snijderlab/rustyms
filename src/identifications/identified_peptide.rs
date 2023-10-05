use super::novor::NovorData;
use super::peaks::PeaksData;
use crate::LinearPeptide;

/// A peptide that is identified by a de novo or database matching program
pub struct IdentifiedPeptide {
    pub(super) peptide: LinearPeptide,
    pub(super) local_confidence: Option<Vec<f64>>,
    pub(super) score: Option<f64>,
    pub(super) metadata: MetaData,
}

/// The definition of all special metadata for all types of identified peptides that can be read
pub enum MetaData {
    /// Peaks metadata
    Peaks(PeaksData),
    /// Novor metadata
    Novor(NovorData),
}
