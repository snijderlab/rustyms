//! This crate handles parsing the [IMGT LIGM-DB database](https://www.imgt.org/) into structures compatible with rustyms.
//! It additionally stores all regions and annotations. There are two main ways of selecting germline(s), specified by name
//! [`get_germline`] or by building a query over the data [`Selection`].
//!
//! <details><summary>Data present per species</summary>
//!
#![doc = include_str!("germlines/germlines.md")]
//!
//! </details>
//!
//! ```
//! use rustyms_imgt::*;
//! let selection = Selection::default()
//!                           .species([Species::HomoSapiens])
//!                           .chain([ChainType::Heavy])
//!                           .gene([GeneType::V]);
//! let first = selection.germlines().next().unwrap();
//! assert_eq!(first.name(), "IGHV1-2*01");
//! ```

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

mod fancy;
mod germlines;
mod select;
mod shared;

pub use fancy::*;
#[cfg(feature = "rayon")]
use germlines::par_germlines;
use germlines::{all_germlines, germlines};

pub use select::*;
#[allow(unused_imports)]
pub use shared::*;
