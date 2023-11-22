/// The type of alignment to perform
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum Type {
    /// Global alignment, which tries to find the best alignment to link both sequences fully to each other, like the Needleman Wunsch algorithm
    Global,
    /// Local alignment, which tries to find the best patch of both sequences to align to each other, this could lead to trailing ends on both sides of both sequences, like the Smith Waterman
    Local,
    /// Hybrid alignment, the second sequence will be fully aligned to the first sequence, this could lead to trailing ends on the first sequence but not on the second.
    GlobalForB,
}

impl Type {
    pub(super) const fn global(self) -> bool {
        !matches!(self, Self::Local)
    }
}

#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentType {
    pub bind_start_a: bool,
    pub bind_start_b: bool,
    pub bind_end_a: bool,
    pub bind_end_b: bool,
}

impl AlignmentType {
    pub const fn global() -> Self {
        Self {
            bind_start_a: true,
            bind_start_b: true,
            bind_end_a: true,
            bind_end_b: true,
        }
    }
    pub const fn global_for_b() -> Self {
        Self {
            bind_start_a: false,
            bind_start_b: true,
            bind_end_a: false,
            bind_end_b: true,
        }
    }
    pub const fn local() -> Self {
        Self {
            bind_start_a: false,
            bind_start_b: false,
            bind_end_a: false,
            bind_end_b: false,
        }
    }
    pub const fn is_global(&self) -> bool {
        self.bind_start_a || self.bind_start_b || self.bind_end_a || self.bind_end_b
    }
}
