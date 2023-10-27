/// All allowed ontologies for modification names
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum Ontology {
    #[default]
    /// Unimod
    Unimod,
    /// PSI-MOD
    Psimod,
}

impl Ontology {
    /// Get the prefix character for the ontology (TODO: the full name is needed when the name is used right, lets make sure the output is always valid pro forma)
    #[allow(dead_code)]
    pub const fn char(self) -> char {
        match self {
            Self::Unimod => 'U',
            Self::Psimod => 'M',
        }
    }
}
