//! A module to contain all alignment related structures and algorithms.

mod align_type;
mod alignment;
mod mass_alignment;
mod piece;
mod scoring;

pub use align_type::Type;
pub use alignment::Alignment;
pub use mass_alignment::align;
pub use piece::Piece;
pub use scoring::{MatchType, BLOSUM62};
