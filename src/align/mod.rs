//! A module to contain all alignment related structures and algorithms.

mod align_type;
mod alignment;
mod linear_peptide;
mod mass_alignment;
mod multiple_sequence_alignment;
mod piece;
mod scoring;

pub use align_type::Type;
pub use alignment::Alignment;
pub use linear_peptide::*;
pub use mass_alignment::MassAlignable;
pub use piece::Piece;
pub use scoring::{MatchType, BLOSUM62};
