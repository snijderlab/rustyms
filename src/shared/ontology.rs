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
}

impl Ontology {
    /// Get the prefix character for the ontology
    #[allow(dead_code)]
    pub const fn char(self) -> char {
        match self {
            Self::Unimod => 'U',
            Self::Psimod => 'M',
            Self::Gnome => 'G',
        }
    }

    /// Get the accession number name for the ontology
    #[allow(dead_code)]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Unimod => "UNIMOD",
            Self::Psimod => "MOD",
            Self::Gnome => "GNO",
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
            },
        )
    }
}

/// The shared type for contact between the build and compile steps
pub type OntologyList = Vec<(usize, String, Modification)>;

/// Any ontology you can lookup modification in
pub trait OntologyLookup {
    /// Find the given name in this ontology
    fn find_name(&self, code: &str) -> Option<Modification>;
    /// Find the given id in this ontology
    fn find_id(&self, id: usize) -> Option<Modification>;
}

impl OntologyLookup for OntologyList {
    fn find_name(&self, code: &str) -> Option<Modification> {
        let code = code.to_ascii_lowercase();
        for option in self {
            if option.1 == code {
                return Some(option.2.clone());
            }
        }
        None
    }

    fn find_id(&self, id: usize) -> Option<Modification> {
        for option in self {
            if option.0 == id {
                return Some(option.2.clone());
            }
        }
        None
    }
}
