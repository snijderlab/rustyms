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
pub mod aminoacid_properties;
mod aminoacids;
mod complex_peptide;
mod element;
pub mod error;
pub mod fragment;
pub mod glycan;
pub mod identifications;
mod isobaric_sets;
mod isotopes;
mod linear_peptide;
pub mod model;
pub mod modification;
mod molecular_charge;
mod neutral_loss;
pub mod ontologies;
pub mod placement_rule;
pub mod rawfile;
mod sequence_element;
pub mod spectrum;
pub mod system;

pub use crate::complex_peptide::ComplexPeptide;
pub use crate::element::*;
pub use crate::formula::*;
pub use crate::isobaric_sets::{building_blocks, find_isobaric_sets, MassTolerance};
pub use crate::linear_peptide::LinearPeptide;
pub use crate::model::Model;
pub use crate::modification::Modification;
pub use crate::neutral_loss::*;
pub use crate::sequence_element::SequenceElement;
pub use crate::spectrum::{AnnotatedSpectrum, MassMode, RawSpectrum};
pub use aminoacids::AminoAcid;
pub use fragment::Fragment;

#[macro_use]
extern crate uom;

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn simple_fragments() {
        let peptide = ComplexPeptide::pro_forma("WFWF")
            .unwrap()
            .singular()
            .unwrap();
        let fragments = peptide
            .generate_theoretical_fragments(system::Charge::new::<system::e>(1.0), &Model::all(), 0)
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
            .generate_theoretical_fragments(system::Charge::new::<system::e>(1.0), &model)
            .unwrap();
        let annotated = spectrum[0].annotate(peptide, &fragments, &model, MassMode::Monoisotopic);
        println!("{annotated:?}");
    }
}
