//! The available ontologies

use std::sync::OnceLock;

use itertools::Itertools;

pub use crate::modification::OntologyModificationList;
use crate::{
    error::{Context, CustomError},
    modification::{Linker, Ontology, OntologyLinkerList, SimpleModification},
};

/// A database of custom modifications
pub type CustomDatabase = (OntologyModificationList, OntologyLinkerList);

/// An empty list of modifications (needed for lifetime reasons)
static EMPTY_LIST: OntologyModificationList = Vec::new();

impl Ontology {
    /// Get the modifications lookup list for this ontology
    pub fn modifications_lookup(
        self,
        custom_database: Option<&CustomDatabase>,
    ) -> &OntologyModificationList {
        match self {
            Self::Gnome => gnome_ontology(),
            Self::Psimod => psimod_ontology(),
            Self::Unimod => unimod_ontology(),
            Self::Xlmod => xlmod_modifications_ontology(),
            Self::Custom => custom_database.map_or(&EMPTY_LIST, |c| &c.0),
        }
    }

    /// Get the linkers lookup list for this ontology
    pub fn linkers_lookup(
        self,
        custom_database: Option<&CustomDatabase>,
    ) -> Option<&OntologyLinkerList> {
        match self {
            Self::Xlmod => Some(xlmod_linker_ontology()),
            Self::Custom => custom_database.map(|c| &c.1),
            _ => None,
        }
    }

    /// Find the closest names in this ontology
    pub fn find_closest(self, code: &str, custom_database: Option<&CustomDatabase>) -> CustomError {
        CustomError::error(
            "Invalid modification",
            format!("The provided name does not exists in {}", self.name()),
            Context::show(code),
        )
        .with_suggestions(Self::similar_names(&[self], code, custom_database))
    }

    /// # Panics
    /// Asserts that the ontology list is not empty.
    pub fn find_closest_many(
        ontologies: &[Self],
        code: &str,
        custom_database: Option<&CustomDatabase>,
    ) -> CustomError {
        assert!(!ontologies.is_empty());
        let names = if ontologies.len() > 1 {
            let mut names = ontologies[..ontologies.len() - 1]
                .iter()
                .map(|o| o.name().to_string())
                .collect_vec();
            let last = names.len() - 1;
            names[last] = format!("{} or {}", names[last], ontologies.last().unwrap().name());
            names.join(", ")
        } else {
            ontologies[0].name().to_string()
        };
        CustomError::error(
            "Invalid modification",
            format!("The provided name does not exists in {names}"),
            Context::show(code),
        )
        .with_suggestions(Self::similar_names(ontologies, code, custom_database))
    }

    /// Get the closest similar names in the given ontologies. Finds both modifications and linkers
    pub fn similar_names(
        ontologies: &[Self],
        code: &str,
        custom_database: Option<&CustomDatabase>,
    ) -> Vec<String> {
        let mut resulting = Vec::new();
        for ontology in ontologies {
            let mut options: Vec<&str> = ontology
                .modifications_lookup(custom_database)
                .iter()
                .map(|option| option.1.as_str())
                .collect();
            options.extend(
                ontology
                    .linkers_lookup(custom_database)
                    .iter()
                    .flat_map(|f| f.iter().map(|option| option.1.as_str())),
            );
            resulting.extend(
                similar::get_close_matches(code, &options, 3, 0.75)
                    .iter()
                    .map(|o| format!("{}:{}", ontology.char(), o)),
            );
        }
        resulting
    }

    /// Find the given name in this ontology.
    pub fn find_modification_name(
        self,
        code: &str,
        custom_database: Option<&CustomDatabase>,
    ) -> Option<SimpleModification> {
        let code = code.to_ascii_lowercase();
        for option in self.modifications_lookup(custom_database) {
            if option.1 == code {
                return Some(option.2.clone());
            }
        }
        None
    }

    /// Find the given name in this ontology.
    pub fn find_linker_name(
        self,
        code: &str,
        custom_database: Option<&CustomDatabase>,
    ) -> Option<Linker> {
        let code = code.to_ascii_lowercase();
        for option in self
            .linkers_lookup(custom_database)
            .iter()
            .flat_map(|f| f.iter())
        {
            if option.1 == code {
                return Some(option.2.clone());
            }
        }
        None
    }

    /// Find the given id in this ontology
    pub fn find_modification_id(
        self,
        id: usize,
        custom_database: Option<&CustomDatabase>,
    ) -> Option<SimpleModification> {
        for option in self.modifications_lookup(custom_database) {
            if option.0 == id {
                return Some(option.2.clone());
            }
        }
        None
    }

    /// Find the given id in this ontology
    pub fn find_linker_id(
        self,
        id: usize,
        custom_database: Option<&CustomDatabase>,
    ) -> Option<Linker> {
        for option in self
            .linkers_lookup(custom_database)
            .iter()
            .flat_map(|f| f.iter())
        {
            if option.0 == id {
                return Some(option.2.clone());
            }
        }
        None
    }
}

/// Get the unimod ontology
/// # Panics
/// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
fn unimod_ontology() -> &'static OntologyModificationList {
    UNIMOD_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/unimod.dat"))).unwrap()
    })
}
/// Get the PSI-MOD ontology
/// # Panics
/// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
fn psimod_ontology() -> &'static OntologyModificationList {
    PSIMOD_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/psimod.dat"))).unwrap()
    })
}
/// Get the Gnome ontology
/// # Panics
/// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
fn gnome_ontology() -> &'static OntologyModificationList {
    GNOME_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/gnome.dat"))).unwrap()
    })
}
/// Get the Xlmod ontology
/// # Panics
/// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
fn xlmod_modifications_ontology() -> &'static OntologyModificationList {
    XLMOD_MODS_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/xlmod_mods.dat"))).unwrap()
    })
}
/// Get the Xlmod ontology
/// # Panics
/// Panics when the linkers are not correctly provided at compile time, always report a panic if it occurs here.
fn xlmod_linker_ontology() -> &'static OntologyLinkerList {
    XLMOD_LINKER_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(
            env!("OUT_DIR"),
            "/xlmod_linkers.dat"
        )))
        .unwrap()
    })
}
static UNIMOD_CELL: OnceLock<OntologyModificationList> = OnceLock::new();
static PSIMOD_CELL: OnceLock<OntologyModificationList> = OnceLock::new();
static GNOME_CELL: OnceLock<OntologyModificationList> = OnceLock::new();
static XLMOD_MODS_CELL: OnceLock<OntologyModificationList> = OnceLock::new();
static XLMOD_LINKER_CELL: OnceLock<OntologyLinkerList> = OnceLock::new();
