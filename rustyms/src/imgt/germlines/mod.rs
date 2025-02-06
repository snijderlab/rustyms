// @generated
#![allow(non_snake_case,non_upper_case_globals)]
use std::sync::OnceLock;
use super::shared::{Germlines, Species};
/// Get the germlines for any of the available species. See the main documentation for which species have which data available.
pub fn germlines(species: Species) -> Option<&'static Germlines> {match species {
_=>None}}
/// Get all germlines in one iterator, see the main documentation for more information about the available germlines
pub fn all_germlines() -> impl std::iter::Iterator<Item = &'static Germlines> {
[
].into_iter()
}
/// Get all germlines in one parallel iterator, see the main documentation for more information about the available germlines
#[cfg(feature = "rayon")]
use rayon::prelude::*;
#[cfg(feature = "rayon")]
pub fn par_germlines() -> impl rayon::prelude::ParallelIterator<Item = &'static Germlines> {
[
].into_par_iter()
}
