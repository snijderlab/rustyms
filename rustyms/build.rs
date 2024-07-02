use std::env;

#[macro_use]
extern crate uom;

#[path = "./src/system.rs"]
mod system;
#[macro_use]
#[path = "./src/helper_functions.rs"]
mod helper_functions;
#[path = "./src/shared/element.rs"]
mod element;
#[macro_use]
#[path = "./src/shared/formula/mod.rs"]
mod formula;
#[path = "./src/build/mod.rs"]
mod build;

#[path = "./src/shared/multi.rs"]
mod multi;

use crate::build::glycan;
use crate::build::*;
pub use crate::element::*;
pub use crate::formula::MolecularFormula;
pub use crate::multi::Multi;

fn main() {
    let debug = env::var("DEBUG_BUILD").map(|v| v == "1").unwrap_or(false);

    let out_dir = env::var_os("OUT_DIR").unwrap();
    build_atomic_masses(&out_dir, debug);
    build_gnome_ontology(&out_dir, debug);
    build_psi_mod_ontology(&out_dir, debug);
    build_resid_ontology(&out_dir, debug);
    build_unimod_ontology(&out_dir, debug);
    build_xlmod_ontology(&out_dir, debug);

    println!("cargo:rerun-if-changed=src/system.rs");
    println!("cargo:rerun-if-changed=src/helper_functions.rs");
    println!("cargo:rerun-if-changed=databases");
    println!("cargo:rerun-if-changed=src/build");
    println!("cargo:rerun-if-changed=build.rs");
    print(out_dir.as_os_str().to_str().unwrap(), debug);
}

fn print(text: impl AsRef<str>, debug: bool) {
    if debug {
        println!("cargo:warning={}", text.as_ref())
    }
}
