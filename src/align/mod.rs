//! Code to make alignments of two peptides based on mass mistakes, and genetic information.
//! A mass based alignment handles the case in which multiple amino acids are wrong, but the total mass
//! of this set of amino acids is equal to the mass of a set of different amino acids on the other peptide.
//! This is quite common in mass spectrometry where mistakes based on mass coincidences are very common.
//! For example `N` has the same mass as `GG`, so if we want to make a mass spectrometry faithful alignment
//! of `ANA` with `AGGA` the result should reflect this fact:
//!
//! ```text
//! Identity: 0.500 (2/4), Similarity: 0.750 (3/4), Gaps: 0.000 (0/4), Score: 0.706 (12/17),
//! Equal mass, Tolerance: 10 ppm, Alignment: global
//! Start: A 0 B 0, Path: 1=1:2i1=
//!
//! AN·A A
//! AGGA B
//!  ╶╴
//! ```
//! _Generated using this algorithm bound to a cli tool: <https://github.com/snijderlab/align-cli>_
//! ```rust
//! use rustyms::*;
//! let a = ComplexPeptide::pro_forma("ANA").unwrap().singular().unwrap();
//! let b = ComplexPeptide::pro_forma("AGGA").unwrap().singular().unwrap();
//! let alignment = align::align::<4>(a, b, &align::BLOSUM62,
//!                    MassTolerance::Ppm(10.0), align::Type::GLOBAL);
//! assert_eq!(alignment.short(), "1=1:2i1=");
//! assert_eq!(alignment.ppm(), 0.0);
//! ```
//!
//! Matrices from: <https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/util/tables/> and <https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/>
//!
//! The UO columns are added, for these the B/J/Z score is the rounded down average of the corresponding non ambiguous AAs. All UO scores are exactly the same for all matrices (except identity).

mod align_type;
mod alignment;
mod diagonal_array;
mod mass_alignment;
mod piece;
mod scoring;

pub use align_type::Type;
pub use alignment::Alignment;
pub use mass_alignment::align;
pub use piece::Piece;
pub use scoring::matrices::*;
pub use scoring::MatchType;
