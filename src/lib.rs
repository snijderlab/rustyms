#![doc = include_str!("../README.md")]
#![allow(dead_code)]
#![warn(clippy::all, clippy::pedantic, clippy::nursery, missing_docs)]
#![allow(
    clippy::must_use_candidate,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::wildcard_imports,
    clippy::module_name_repetitions,
    clippy::suboptimal_flops,
    clippy::too_many_lines
)]

#[cfg(test)]
mod fragmentation_tests;
#[cfg(test)]
mod pro_forma_parse_tests;
#[macro_use]
mod helper_functions;
#[macro_use]
mod formula;

pub mod align;
mod aminoacids;
mod complex_peptide;
mod element;
mod error;
mod fragment;
mod glycan;
mod helper_functions;
pub mod identifications;
mod isobaric_sets;
mod isotopes;
mod linear_peptide;
mod model;
mod modification;
mod molecular_charge;
mod neutral_loss;
mod ontologies;
mod placement_rules;
pub mod rawfile;
mod spectrum;
mod system;

pub use crate::complex_peptide::*;
pub use crate::element::*;
pub use crate::error::*;
pub use crate::formula::*;
pub use crate::fragment::*;
pub use crate::glycan::*;
pub use crate::isobaric_sets::*;
pub use crate::linear_peptide::*;
pub use crate::model::*;
pub use crate::modification::*;
pub use crate::neutral_loss::*;
pub use crate::spectrum::*;
pub use crate::system::f64::*;
pub use aminoacids::AminoAcid;
pub use fragment::Fragment;
pub use linear_peptide::LinearPeptide;
pub use model::Model;
pub use uom::num_traits::Zero;

#[macro_use]
extern crate uom;

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn simple_fragments() {
        let peptide = ComplexPeptide::pro_forma("WFWF").unwrap().assume_linear();
        let fragments = peptide
            .generate_theoretical_fragments(Charge::new::<e>(1.0), &Model::all(), 0)
            .unwrap();
        println!("{}", fragments.len());
        println!("{fragments:?}");
    }

    #[test]
    fn simple_matching() {
        let model = Model::all();
        let spectrum = rawfile::mgf::open("data/example.mgf").unwrap();
        let peptide = ComplexPeptide::pro_forma("WFWF").unwrap();
        let fragments = peptide
            .generate_theoretical_fragments(Charge::new::<e>(1.0), &model)
            .unwrap();
        let annotated = spectrum[0].annotate(peptide, &fragments, &model);
        println!("{annotated:?}");
    }
}
