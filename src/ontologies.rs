//! The available ontologies

use std::sync::OnceLock;

pub use crate::modification::OntologyList;
pub use crate::modification::OntologyLookup;

/// Get the unimod ontology
/// # Panics
/// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
pub fn unimod_ontology() -> &'static OntologyList {
    UNIMOD_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/unimod.dat"))).unwrap()
    })
}
/// Get the PSI-MOD ontology
/// # Panics
/// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
pub fn psimod_ontology() -> &'static OntologyList {
    PSIMOD_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/psimod.dat"))).unwrap()
    })
}
/// Get the Gnome ontology
/// # Panics
/// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
pub fn gnome_ontology() -> &'static OntologyList {
    GNOME_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/gnome.dat"))).unwrap()
    })
}
static UNIMOD_CELL: OnceLock<OntologyList> = OnceLock::new();
static PSIMOD_CELL: OnceLock<OntologyList> = OnceLock::new();
static GNOME_CELL: OnceLock<OntologyList> = OnceLock::new();
