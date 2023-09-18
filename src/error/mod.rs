//! Contain the definition for errors with all additional data that is needed to generate nice error messages

/// The context of an error
mod context;
/// An error with all its properties
mod error;

pub use context::{Context, FilePosition};
pub use error::CustomError;
