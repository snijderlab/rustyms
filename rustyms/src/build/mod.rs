mod atomic_masses;
#[path = "../shared/csv.rs"]
mod csv;
#[path = "../error/mod.rs"]
pub mod error;
pub mod glycan;
mod gnome;
mod obo;
mod ontology_modification;
mod psi_mod;
mod resid;
mod unimod;
mod xlmod;

pub use atomic_masses::*;
pub use gnome::*;
use ontology_modification::*;
pub use psi_mod::*;
pub use resid::*;
pub use unimod::*;
pub use xlmod::*;

use serde::{Deserialize, Serialize};

use crate::{
    build::glycan::{GlycanStructure, MonoSaccharide},
    formula::MolecularFormula,
    formula::MultiChemical,
    system::OrderedMass,
    Multi,
};

use ordered_float::OrderedFloat;
use std::cmp::Ordering;
use std::collections::HashSet;

include!("../shared/neutral_loss.rs");
include!("../shared/modification.rs");
include!("../shared/aminoacid.rs");

impl crate::Element {
    pub fn is_valid(self, _isotope: Option<std::num::NonZeroU16>) -> bool {
        true
    }
}
