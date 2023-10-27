use std::cell::OnceCell;

use crate::aminoacids::*;
use crate::element::*;
use crate::glycan::*;
use crate::modification::*;
use crate::placement_rules::*;

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
const UNIMOD_CELL: OnceCell<Vec<(usize, String, Modification)>> = OnceCell::new();
const PSIMOD_CELL: OnceCell<Vec<(usize, String, Modification)>> = OnceCell::new();
const GNOME_CELL: OnceCell<Vec<(usize, String, Modification)>> = OnceCell::new();
