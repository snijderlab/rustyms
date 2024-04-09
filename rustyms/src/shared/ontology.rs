/// All allowed ontologies for modification names
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum Ontology {
    #[default]
    /// Unimod
    Unimod,
    /// PSI-MOD
    Psimod,
    /// GNOme
    Gnome,
    /// XLMOD
    Xlmod,
}

impl Ontology {
    /// Get the prefix character for the ontology
    #[allow(dead_code)]
    pub const fn char(self) -> char {
        match self {
            Self::Unimod => 'U',
            Self::Psimod => 'M',
            Self::Gnome => 'G',
            Self::Xlmod => 'X',
        }
    }

    /// Get the accession number name for the ontology
    #[allow(dead_code)]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Unimod => "UNIMOD",
            Self::Psimod => "MOD",
            Self::Gnome => "GNO",
            Self::Xlmod => "XLMOD",
        }
    }
}

impl std::fmt::Display for Ontology {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Unimod => "Unimod",
                Self::Psimod => "PSI-MOD",
                Self::Gnome => "GNOme",
                Self::Xlmod => "XLMOD",
            },
        )
    }
}

/// The shared type for contact between the build and compile steps
pub type OntologyList = Vec<(usize, String, Modification)>;
