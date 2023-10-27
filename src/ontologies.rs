use std::sync::OnceLock;

use crate::modification::Modification;

/// Get the unimod ontology
pub fn unimod_ontology() -> &'static Vec<(usize, String, Modification)> {
    UNIMOD_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/unimod.dat"))).unwrap()
    })
}
/// Get the PSI-MOD ontology
pub fn psimod_ontology() -> &'static Vec<(usize, String, Modification)> {
    PSIMOD_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/psimod.dat"))).unwrap()
    })
}
/// Get the Gnome ontology
pub fn gnome_ontology() -> &'static Vec<(usize, String, Modification)> {
    GNOME_CELL.get_or_init(|| {
        bincode::deserialize(include_bytes!(concat!(env!("OUT_DIR"), "/gnome.dat"))).unwrap()
    })
}
static UNIMOD_CELL: OnceLock<Vec<(usize, String, Modification)>> = OnceLock::new();
static PSIMOD_CELL: OnceLock<Vec<(usize, String, Modification)>> = OnceLock::new();
static GNOME_CELL: OnceLock<Vec<(usize, String, Modification)>> = OnceLock::new();
