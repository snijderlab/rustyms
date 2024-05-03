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

#[cfg(feature = "align")]
/// Only available with feature `align`.
pub mod align;

#[cfg(feature = "identification")]
/// Only available with feature `identification`.
pub mod identification;

#[cfg(feature = "imgt")]
/// Only available with feature `imgt`.
pub mod imgt;

#[cfg(test)]
mod fragmentation_tests;
#[cfg(test)]
mod pro_forma_parse_tests;
#[macro_use]
mod helper_functions;
#[macro_use]
mod formula;

#[doc(hidden)]
#[path = "shared/csv.rs"]
pub mod csv;

pub mod aminoacid_properties;
mod aminoacids;
mod element;
pub mod error;
pub mod fragment;
mod fragmentation;
pub mod glycan;
mod isobaric_sets;
mod isotopes;
mod itertools_extension;
mod mass_mode;
pub mod model;
pub mod modification;
mod molecular_charge;
mod multi;
mod multi_formula;
mod neutral_loss;
pub mod ontologies;
mod peptide;
pub mod placement_rule;
mod protease;
pub mod rawfile;
mod sequence_element;
pub mod spectrum;
pub mod system;
mod tolerance;

pub use crate::element::*;
pub use crate::formula::*;
pub use crate::isobaric_sets::{building_blocks, find_isobaric_sets};
pub use crate::mass_mode::MassMode;
pub use crate::model::Model;
pub use crate::modification::Modification;
pub use crate::molecular_charge::MolecularCharge;
pub use crate::multi::*;
pub use crate::multi_formula::*;
pub use crate::neutral_loss::*;
pub use crate::peptide::{
    CompoundPeptidoform, ExtremelySimple, Linear, LinearPeptide, Linked, Peptidoform, Simple,
    VerySimple,
};
pub use crate::protease::*;
pub use crate::sequence_element::SequenceElement;
pub use crate::spectrum::{AnnotatedSpectrum, RawSpectrum};
pub use crate::tolerance::*;
pub use aminoacids::AminoAcid;
pub use fragment::Fragment;

#[macro_use]
extern crate uom;

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod test {
    use super::*;

    #[test]
    fn simple_fragments() {
        let peptide = LinearPeptide::pro_forma("WFWF", None)
            .unwrap()
            .linear()
            .unwrap();
        let fragments = peptide.generate_theoretical_fragments(
            system::usize::Charge::new::<system::e>(1),
            &Model::all(),
        );
        println!("{}", fragments.len());
        println!("{fragments:?}");
    }

    #[test]
    fn simple_matching() {
        let model = Model::all();
        let spectrum = rawfile::mgf::open("data/example.mgf").unwrap();
        let peptide = CompoundPeptidoform::pro_forma("WFWF", None).unwrap();
        let fragments = peptide
            .generate_theoretical_fragments(system::usize::Charge::new::<system::e>(1), &model);
        let annotated = spectrum[0].annotate(peptide, &fragments, &model, MassMode::Monoisotopic);
        println!("{annotated:?}");
    }
}
