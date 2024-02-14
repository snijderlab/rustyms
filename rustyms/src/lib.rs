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
//! * `consecutive-align` - enables consecutive align, where the mass alignment is done on multiple
//!   IMGT germlines consecutively (enables both `imgt` and `align`).
//! * `rayon` - enables parallel iterators using rayon, mostly for `imgt` but also in consecutive
//!   align (enables both `imgt` and `align`).

#[cfg(feature = "rustyms-align")]
#[doc(inline)]
/// Only available with feature `align`.
pub use rustyms_align as align;

pub use rustyms_core::*;

#[cfg(feature = "rustyms-identification")]
#[doc(inline)]
/// Only available with feature `identification`.
pub use rustyms_identification as identification;

#[cfg(feature = "rustyms-imgt")]
#[doc(inline)]
/// Only available with feature `imgt`.
pub use rustyms_imgt as imgt;
