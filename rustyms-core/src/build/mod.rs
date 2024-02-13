mod atomic_masses;
#[path = "../shared/csv.rs"]
mod csv;
#[path = "../error/mod.rs"]
pub mod error;
mod gnome;
mod obo;
mod ontology_modification;
mod psi_mod;
mod unimod;

pub use atomic_masses::*;
pub use gnome::*;
use ontology_modification::*;
pub use psi_mod::*;
pub use unimod::*;

use serde::{Deserialize, Serialize};

use crate::{
    formula::MolecularFormula,
    glycan::{GlycanStructure, MonoSaccharide},
    system::OrderedMass,
};

include!("../shared/modification.rs");

impl crate::Element {
    pub fn is_valid(self, _isotope: Option<u16>) -> bool {
        true
    }
}
