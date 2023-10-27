use std::sync::OnceLock;

use crate::OntologyList;

/// Get the unimod ontology
pub fn unimod_ontology() -> &'static OntologyList {
    UNIMOD_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/unimod.dat"))).unwrap()
    })
}
/// Get the PSI-MOD ontology
pub fn psimod_ontology() -> &'static OntologyList {
    PSIMOD_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/psimod.dat"))).unwrap()
    })
}
/// Get the Gnome ontology
pub fn gnome_ontology() -> &'static OntologyList {
    GNOME_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/gnome.dat"))).unwrap()
    })
}
static UNIMOD_CELL: OnceLock<OntologyList> = OnceLock::new();
static PSIMOD_CELL: OnceLock<OntologyList> = OnceLock::new();
static GNOME_CELL: OnceLock<OntologyList> = OnceLock::new();
