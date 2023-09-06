use super::scoring::MatchType;

/// A piece in an alignment, determining what step was taken in the alignment and how this impacted the score
#[derive(Clone, Default, Debug)]
pub struct Piece {
    /// The total score of the path up till now
    pub score: isize,
    /// The local contribution to the score of this piece
    pub local_score: i8,
    /// The type of the match
    pub match_type: MatchType,
    /// The number of steps on the first sequence
    pub step_a: u8,
    /// The number of steps on the second sequence
    pub step_b: u8,
}

impl Piece {
    /// Create a new alignment piece
    pub const fn new(
        score: isize,
        local_score: i8,
        match_type: MatchType,
        step_a: u8,
        step_b: u8,
    ) -> Self {
        Self {
            score,
            local_score,
            match_type,
            step_a,
            step_b,
        }
    }
}

impl Piece {
    /// Display this piece very compactly
    pub(super) fn short(&self) -> String {
        match (self.step_a, self.step_b) {
            (0, 1) => "I".to_string(),
            (1, 0) => "D".to_string(),
            (1, 1) => "M".to_string(),
            (a, b) => format!("S[{b},{a}]"),
        }
    }
}
