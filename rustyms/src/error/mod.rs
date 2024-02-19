//! Contain the definition for errors with all additional data that is needed to generate nice error messages

/// The context of an error
mod context;
/// An error with all its properties
mod custom_error;

#[allow(unused_imports)]
pub use context::{Context, FilePosition};
pub use custom_error::CustomError;
