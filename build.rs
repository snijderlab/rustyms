use std::env;

#[path = "./src/shared/element.rs"]
mod element;
#[macro_use]
#[path = "./src/shared/formula.rs"]
mod formula;
#[path = "./src/shared/glycan.rs"]
mod glycan;

#[path = "./src/helper_functions.rs"]
mod helper_functions;

#[path = "./src/build/mod.rs"]
mod build;

use crate::build::*;
pub use crate::element::*;
use crate::formula::Chemical;
use crate::formula::*;

fn main() {
    let debug = env::var("DEBUG_BUILD").map(|v| v == "1").unwrap_or(false);

    let out_dir = env::var_os("OUT_DIR").unwrap();
    build_unimod_ontology(&out_dir, debug);
    build_psi_mod_ontology(&out_dir, debug);
    build_atomic_masses(&out_dir, debug).unwrap();

    println!("cargo:rerun-if-changed=databases/unimod.obo");
    println!("cargo:rerun-if-changed=databases/PSI-MOD-newstyle.obo");
    println!("cargo:rerun-if-changed=databases/IUPAC-atomic-masses.csv");
    println!("cargo:rerun-if-changed=databases/CIAAW-atomic-weights.csv");
    println!("cargo:rerun-if-changed=databases/CIAAW-isotopic-abundances.csv");
    println!("cargo:rerun-if-changed=src/build/*");
    println!("cargo:rerun-if-changed=build.rs");
    print(out_dir.as_os_str().to_str().unwrap(), debug);
}

fn print(text: impl AsRef<str>, debug: bool) {
    if debug {
        println!("cargo:warning={}", text.as_ref())
    }
}
