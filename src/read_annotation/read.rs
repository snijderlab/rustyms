use crate::ComplexPeptide;

pub struct Read {
    peptide: ComplexPeptide,
    local_confidence: Option<Vec<f64>>,
    score: Option<f64>,
    metadata: MetaData,
}

pub enum MetaData {
    Peaks(),
}
