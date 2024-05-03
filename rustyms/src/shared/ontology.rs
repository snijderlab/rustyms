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
    /// Custom
    Custom,
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
            Self::Custom => 'C',
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
            Self::Custom => "CUSTOM",
        }
    }

    /// Get the accession number name for the ontology
    #[allow(dead_code)]
    pub fn url(self, id: usize, name: &str) -> Option<String> {
        match self {
            Self::Unimod => Some(format!(
                "https://www.unimod.org/modifications_view.php?editid1={id}",
            )),
            Self::Psimod => Some(format!(
                "https://ontobee.org/ontology/MOD?iri=http://purl.obolibrary.org/obo/MOD_{id:5}",
            )),
            Self::Gnome => Some(format!(
                "https://gnome.glyomics.org/StructureBrowser.html?focus={name}"
            )),
            Self::Xlmod | Self::Custom => None,
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
                Self::Custom => "Custom",
            },
        )
    }
}

/// The shared type for contact between the build and compile steps
pub type OntologyModificationList = Vec<(usize, String, SimpleModification)>;
