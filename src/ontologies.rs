//! The available ontologies

use std::sync::OnceLock;

use itertools::Itertools;

pub use crate::modification::OntologyList;
use crate::{
    error::{Context, CustomError},
    modification::Ontology,
    Modification,
};

impl Ontology {
    /// Get the lookup list for this ontology
    pub fn lookup(self) -> &'static OntologyList {
        match self {
            Self::Gnome => gnome_ontology(),
            Self::Psimod => psimod_ontology(),
            Self::Unimod => unimod_ontology(),
        }
    }

    /// Find the closest names in this ontology
    pub fn find_closest(self, code: &str) -> CustomError {
        CustomError::error(
            "Invalid modification",
            format!("The provided name does not exists in {}", self.name()),
            Context::show(code),
        )
        .with_suggestions(Self::similar_names(&[self], code))
    }

    /// # Panics
    /// Asserts that the ontology list is not empty.
    pub fn find_closest_many(ontologies: &[Self], code: &str) -> CustomError {
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
        .with_suggestions(Self::similar_names(ontologies, code))
    }

    /// Get the closest similar names in the given ontologies.
    pub fn similar_names(ontologies: &[Self], code: &str) -> Vec<String> {
        let mut resulting = Vec::new();
        for ontology in ontologies {
            let options: Vec<&str> = ontology
                .lookup()
                .iter()
                .map(|option| option.1.as_str())
                .collect();
            resulting.extend(
                similar::get_close_matches(code, &options, 3, 0.75)
                    .iter()
                    .map(|o| format!("{}:{}", ontology.char(), o)),
            );
        }
        resulting
    }

    /// Find the given name in this ontology. Return the matching modification, or if no matching modification is found the closest modifications in the used ontology
    pub fn find_name(self, code: &str) -> Option<Modification> {
        let code = code.to_ascii_lowercase();
        for option in self.lookup() {
            if option.1 == code {
                return Some(option.2.clone());
            }
        }
        None
    }

    /// Find the given id in this ontology
    pub fn find_id(self, id: usize) -> Option<Modification> {
        for option in self.lookup() {
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
fn unimod_ontology() -> &'static OntologyList {
    UNIMOD_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/unimod.dat"))).unwrap()
    })
}
/// Get the PSI-MOD ontology
/// # Panics
/// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
fn psimod_ontology() -> &'static OntologyList {
    PSIMOD_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/psimod.dat"))).unwrap()
    })
}
/// Get the Gnome ontology
/// # Panics
/// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
fn gnome_ontology() -> &'static OntologyList {
    GNOME_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/gnome.dat"))).unwrap()
    })
}
static UNIMOD_CELL: OnceLock<OntologyList> = OnceLock::new();
static PSIMOD_CELL: OnceLock<OntologyList> = OnceLock::new();
static GNOME_CELL: OnceLock<OntologyList> = OnceLock::new();
