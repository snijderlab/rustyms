/// The context of an error
mod context;
/// An error with all its properties
mod error;

pub use context::{Context, FilePosition};
pub use error::CustomError;
