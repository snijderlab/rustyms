use serde::{Deserialize, Serialize};

/// The type of alignment to perform
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum Type {
    /// Global alignment, which tries to find the best alignment to link both sequences fully to each other, like the Needleman Wunsch algorithm
    #[default]
    Global,
    /// Local alignment, which tries to find the best patch of both sequences to align to each other, this could lead to trailing ends on both sides of both sequences, like the Smith Waterman
    Local,
    /// Hybrid alignment, the second sequence will be fully aligned to the first sequence, this could lead to trailing ends on the first sequence but not on the second.
    GlobalForB,
    /// Hybrid alignment, the first sequence will be fully aligned to the second sequence, this could lead to trailing ends on the second sequence but not on the first.
    GlobalForA,
}

impl Type {
    pub(super) const fn global(self) -> bool {
        !matches!(self, Self::Local)
    }
}
