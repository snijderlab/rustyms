//! Handle mass spectrometry data in Rust. This crate is set up to handle very complex peptides with
//! loads of ambiguity and complexity. It pivots around the [`ComplexPeptide`] and [`LinearPeptide`]
//! which encode the [ProForma](https://github.com/HUPO-PSI/ProForma) specification. Additionally
//! this crate enables the reading of [mgf](rawfile::mgf), doing [spectrum annotation](RawSpectrum::annotate)
//! (BU/MD/TD), finding [isobaric sequences](find_isobaric_sets), doing [alignments of peptides](align::align)
//! , accessing the [IMGT germline database](imgt), and [reading identified peptide files](identification).
//!
//! ```
//! # fn main() -> Result<(), rustyms::error::CustomError> {
//! # let raw_file_path = "../rustyms-core/data/annotated_example.mgf";
//! // Open some data and see if the given peptide is a valid match
//! use rustyms::{*, system::{Charge, e}};
//! let peptide = ComplexPeptide::pro_forma("Q[Gln->pyro-Glu]VQEVSERTHGGNFD")?;
//! let spectrum = rawfile::mgf::open(raw_file_path)?;
//! let model = Model::ethcd();
//! let fragments = peptide.generate_theoretical_fragments(Charge::new::<e>(2.0), &model);
//! let annotated = spectrum[0].annotate(peptide, &fragments, &model, MassMode::Monoisotopic);
//! let fdr = annotated.fdr(&fragments, &model);
//! // This is the incorrect sequence for this spectrum so the FDR will indicate this
//! # dbg!(&fdr, fdr.sigma(), fdr.fdr(), fdr.score());
//! assert!(fdr.sigma() < 2.0);
//! # Ok(()) }
//! ```
//!
//! ```
//! # fn main() -> Result<(), rustyms::error::CustomError> {
//! // Check how this peptide compares to a similar peptide (using `align`)
//! // (same sequence, repeated for easy reference)
//! use rustyms::{*, align::*};
//! let first_peptide = LinearPeptide::pro_forma("Q[Gln->pyro-Glu]VQEVS")?;
//! let second_peptide = LinearPeptide::pro_forma("E[Glu->pyro-Glu]VQVES")?;
//! let alignment = align::<4>(&first_peptide, &second_peptide,
//!                  align::BLOSUM62, Tolerance::new_ppm(10.0), AlignType::GLOBAL);
//! # dbg!(&alignment);
//! let stats = alignment.stats();
//! # //assert_eq!(stats.identical, 3); // Only three positions are identical
//! assert_eq!(stats.mass_similar, 6); // All positions are mass similar
//! # Ok(()) }
//! ```
//!
//! Rustyms is the main crate tying together multiple smaller crates into one cohesive structure.
//! It has multiple features which allow you to slim it down if needed (all are enabled by default).
//! * `identification` - gives access to methods reading many different identified peptide formats.
//! * `align` - gives access to mass based alignment of peptides.
//! * `imgt` - enables access to the IMGT database of antibodies germline sequences, with annotations.
//! * `rayon` - enables parallel iterators using rayon, mostly for `imgt` but also in consecutive
//!   align.

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
mod complex_peptide;
mod element;
pub mod error;
pub mod fragment;
mod fragmentation;
pub mod glycan;
mod isobaric_sets;
mod isotopes;
mod itertools_extension;
mod linear_peptide;
pub mod model;
pub mod modification;
mod molecular_charge;
mod multi;
mod multi_formula;
mod neutral_loss;
pub mod ontologies;
pub mod placement_rule;
mod protease;
pub mod rawfile;
mod sequence_element;
pub mod spectrum;
pub mod system;
mod tolerance;

pub use crate::complex_peptide::ComplexPeptide;
pub use crate::element::*;
pub use crate::formula::*;
pub use crate::isobaric_sets::{building_blocks, find_isobaric_sets};
pub use crate::linear_peptide::LinearPeptide;
pub use crate::model::Model;
pub use crate::modification::Modification;
pub use crate::multi::*;
pub use crate::multi_formula::*;
pub use crate::neutral_loss::*;
pub use crate::protease::*;
pub use crate::sequence_element::SequenceElement;
pub use crate::spectrum::{AnnotatedSpectrum, MassMode, RawSpectrum};
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
        let peptide = ComplexPeptide::pro_forma("WFWF")
            .unwrap()
            .singular()
            .unwrap();
        let fragments = peptide
            .generate_theoretical_fragments(system::Charge::new::<system::e>(1.0), &Model::all());
        println!("{}", fragments.len());
        println!("{fragments:?}");
    }

    #[test]
    fn simple_matching() {
        let model = Model::all();
        let spectrum = rawfile::mgf::open("data/example.mgf").unwrap();
        let peptide = ComplexPeptide::pro_forma("WFWF").unwrap();
        let fragments =
            peptide.generate_theoretical_fragments(system::Charge::new::<system::e>(1.0), &model);
        let annotated = spectrum[0].annotate(peptide, &fragments, &model, MassMode::Monoisotopic);
        println!("{annotated:?}");
    }
}
