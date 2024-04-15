mod atomic_masses;
#[path = "../shared/csv.rs"]
mod csv;
#[path = "../error/mod.rs"]
pub mod error;
mod glycan;
mod gnome;
mod obo;
mod ontology_modification;
mod psi_mod;
mod unimod;
mod xlmod;

pub use atomic_masses::*;
pub use gnome::*;
use ontology_modification::*;
pub use psi_mod::*;
pub use unimod::*;

use serde::{Deserialize, Serialize};

use crate::{
    build::glycan::{GlycanStructure, MonoSaccharide},
    formula::MolecularFormula,
    system::OrderedMass,
};

include!("../shared/neutral_loss.rs");
include!("../shared/modification.rs");

impl crate::Element {
    pub fn is_valid(self, _isotope: Option<std::num::NonZeroU16>) -> bool {
        true
    }
}
