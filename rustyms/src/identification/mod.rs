//! Read in the annotations from peptide identification sources

#[macro_use]
mod common_parser;

mod fasta;
mod general;
mod identified_peptide;
mod maxquant;
mod msfragger;
mod mztab;
mod novor;
mod opair;
mod peaks;
mod sage;

use crate::*;
pub use fasta::*;
pub use general::*;
pub use identified_peptide::*;
pub use maxquant::*;
pub use msfragger::*;
pub use mztab::*;
pub use novor::*;
pub use opair::*;
pub use peaks::*;
pub use sage::*;

#[cfg(test)]
mod maxquant_tests;
#[cfg(test)]
mod msfragger_tests;
#[cfg(test)]
mod mztab_test;
#[cfg(test)]
mod novor_tests;
#[cfg(test)]
mod opair_tests;
#[cfg(test)]
mod peaks_tests;
#[cfg(test)]
mod sage_tests;
