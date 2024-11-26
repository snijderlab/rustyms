//! Read in the annotations from peptide identification sources

#[macro_use]
mod common_parser;

mod deepnovofamily;
mod fasta;
mod general;
mod identified_peptide;
mod instanovo;
mod maxquant;
mod msfragger;
mod mztab;
mod novob;
mod novor;
mod opair;
mod peaks;
mod pepnet;
mod plink;
mod powernovo;
mod sage;

use crate::*;
pub use deepnovofamily::*;
pub use fasta::*;
pub use general::*;
pub use identified_peptide::*;
pub use instanovo::*;
pub use maxquant::*;
pub use msfragger::*;
pub use mztab::*;
pub use novob::*;
pub use novor::*;
pub use opair::*;
pub use peaks::*;
pub use pepnet::*;
pub use plink::*;
pub use powernovo::*;
pub use sage::*;

#[cfg(test)]
mod deepnovofamily_tests;
#[cfg(test)]
mod instanovo_tests;
#[cfg(test)]
mod maxquant_tests;
#[cfg(test)]
mod msfragger_tests;
#[cfg(test)]
mod mztab_test;
#[cfg(test)]
mod novob_tests;
#[cfg(test)]
mod novor_tests;
#[cfg(test)]
mod opair_tests;
#[cfg(test)]
mod peaks_tests;
#[cfg(test)]
mod pepnet_tests;
#[cfg(test)]
mod plink_tests;
#[cfg(test)]
mod powernovo_tests;
#[cfg(test)]
mod sage_tests;
