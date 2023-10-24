mod atomic_masses;
#[path = "../shared/csv.rs"]
mod csv;
#[path = "../error/mod.rs"]
mod error;
#[path = "../shared/modification.rs"]
mod modification;
mod gnome;
mod obo;
mod ontology_modification;
mod psi_mod;
mod unimod;

pub use atomic_masses::*;
pub use error::*;
pub use gnome::*;
pub use psi_mod::*;
pub use unimod::*;
use modification::*;

trait ToCode {
    fn to_code(&self) -> String;
}
