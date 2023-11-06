use std::{borrow::Cow, fmt::Display};

use crate::{
    align::{MatchType, Piece},
    system::Mass,
    LinearPeptide, MassTolerance, SequenceElement,
};

use super::{align_type::AlignmentType, Alignment};

/// An alignment of multiple peptides
#[derive(Debug, Clone, PartialEq)]
pub struct MultipleSequenceAlignment {
    /// The sequences
    pub sequences: Vec<MSAPlacement>,
    /// The alignment type
    pub ty: AlignmentType,
}

/// The placement of a single peptide in a multiple sequence alignment
#[derive(Debug, Clone, PartialEq)]
pub struct MSAPlacement {
    /// The sequence
    pub sequence: LinearPeptide,
    /// The start position in the sequence where the alignment starts
    pub start: usize,
    /// The path the alignment follows
    pub path: Vec<MSAPosition>,
    /// The absolute score
    pub score: usize,
    /// The normalised score
    pub normalised_score: f64,
}

/// A single position in a multiple sequence alignment
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum MSAPosition {
    /// A gap one sequence position spanning one sequence position
    Gap,
    /// A selection of sequence positions with its width, the maximal number of other sequence positions it spans
    Placed(usize, usize),
}

impl MultipleSequenceAlignment {
    /// Normalise the alignment by making sure all steps are expanded (no steps that squish any piece)
    /// and that any location that is expanded for every sequence is squished back.
    fn normalise(&mut self) {
        todo!();
    }

    /// If this multiple sequence alignment consists of only a single pair, create a pairwise mass alignment struct.
    fn assume_single(&self) -> Option<Alignment> {
        todo!();
    }

    /// Determine the scores for all alignments (should be called after the full alignment is done)
    fn determine_scores(&mut self) {
        // For each location find the highest score it has for all other sets
        todo!();
    }
}

impl Display for MultipleSequenceAlignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for sequence in &self.sequences {
            write!(f, "{}", " ".repeat(sequence.start))?;
            let mut seq = sequence.sequence.sequence.iter().skip(sequence.start);
            for step in &sequence.path {
                match step {
                    MSAPosition::Gap => write!(f, "-")?,
                    MSAPosition::Placed(steps, width) => write!(
                        f,
                        "{}{}",
                        (&mut seq)
                            .take(*steps)
                            .map(|s| s.aminoacid.char())
                            .collect::<String>(),
                        "·".repeat(width.saturating_sub(*steps)) // TODO: handle the cases where steps is too big
                    )?,
                }
            }
            writeln!(
                f,
                " [Score: {}, Normalised score: {:.3}]",
                sequence.score, sequence.normalised_score
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::ComplexPeptide;

    use super::*;

    #[test]
    fn print_msma() {
        let msma = MultipleSequenceAlignment {
            sequences: vec![
                MSAPlacement {
                    sequence: ComplexPeptide::pro_forma("WGGD").unwrap().assume_linear(),
                    start: 0,
                    path: vec![
                        MSAPosition::Placed(1, 1),
                        MSAPosition::Placed(1, 1),
                        MSAPosition::Placed(1, 1),
                        MSAPosition::Placed(1, 1),
                    ],
                    score: 0,
                    normalised_score: 0.0,
                },
                MSAPlacement {
                    sequence: ComplexPeptide::pro_forma("WND").unwrap().assume_linear(),
                    start: 0,
                    path: vec![
                        MSAPosition::Placed(1, 1),
                        MSAPosition::Placed(1, 2),
                        MSAPosition::Placed(1, 1),
                    ],
                    score: 0,
                    normalised_score: 0.0,
                },
            ],
            ty: AlignmentType {
                bind_start_a: true,
                bind_start_b: true,
                bind_end_a: true,
                bind_end_b: true,
            },
        };
        println!("{msma}");
        assert_eq!(
            msma.to_string(),
            "WGGD [Score: 0, Normalised score: 0.000]\nWN·D [Score: 0, Normalised score: 0.000]\n"
        );
    }
}
