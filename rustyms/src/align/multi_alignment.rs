use std::borrow::Cow;

use crate::LinearPeptide;

use super::{AlignType, MatchType, Score};

use serde::{Deserialize, Serialize};

type MultiAlignment<'lifetime, Complexity> = Vec<MultiAlignmentLine<'lifetime, Complexity>>;

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
struct MultiAlignmentLine<'lifetime, Complexity> {
    sequence: Cow<'lifetime, LinearPeptide<Complexity>>,
    path: Vec<MultiPiece>,
    score: Score,
    start: usize,
    align_type: AlignType,
    maximal_step: u16,
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
struct MultiPiece {
    score: isize,
    local_score: isize,
    match_type: MatchType,
    step: u16,
}

impl<'lifetime, Complexity> MultiAlignmentLine<'lifetime, Complexity> {
    fn debug_display(&self) {
        for piece in self
            .path
            .iter()
            .zip(self.sequence.sequence().iter().skip(self.start))
        {
            print!(
                "{}{}",
                piece.1.aminoacid.char(),
                "·".repeat(piece.0.step as usize - 1)
            );
        }
        println!();
    }
}