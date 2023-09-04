use crate::Peptide;

#[derive(Debug, Clone, PartialEq, Default)]
pub struct PeptideCollection {
    pub kind: CollectionKind,
    pub peptides: Vec<Peptide>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub enum CollectionKind {
    #[default]
    Multimeric,
}

impl PeptideCollection {
    /// Assume there is exactly one peptide in this collection
    /// # Panics
    /// If there are no or multiple peptides.
    pub fn assume_singular(mut self) -> Peptide {
        if self.peptides.len() == 1 {
            self.peptides.pop().unwrap()
        } else if self.peptides.len() > 1 {
            panic!("This collection contains multiple spectra, while a single one is assumed")
        } else {
            panic!("This collection does not contain any spectra, while a single one is assumed")
        }
    }
}
