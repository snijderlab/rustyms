use std::{any::TypeId, fmt::Display};

use super::{align_type::AlignmentType, Alignment, MassAlignable, MatchType};
use crate::LinearPeptide;
use std::fmt::Write;

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

impl MSAPlacement {
    pub fn len(&self) -> usize {
        self.path.iter().map(|p| p.len()).sum()
    }
}

/// A single position in a multiple sequence alignment
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum MSAPosition {
    /// A gap one sequence position spanning one sequence position
    Gap,
    /// A selection of sequence positions with its width, the maximal number of other sequence positions it spans
    Placed(MatchType, usize, usize),
}

impl MultipleSequenceAlignment {
    /// Normalise the alignment by making sure all steps are expanded (no steps that squish any piece)
    /// and that any location that is expanded for every sequence is squished back.
    pub(super) fn normalise(&mut self) {
        fn fix(
            msa: &mut MultipleSequenceAlignment,
            search_index: usize,
            expand: usize,
            except: usize,
        ) {
            for sequence_index in (0..msa.sequences.len()).filter(|i| *i != except) {
                let mut index = msa.sequences[sequence_index].start;
                for path_index in 0..msa.sequences[sequence_index].path.len() {
                    if let MSAPosition::Placed(_, _, b) =
                        &mut msa.sequences[sequence_index].path[path_index]
                    {
                        if index + *b >= search_index {
                            *b += expand;
                            break;
                        }
                        index += *b;
                    } else {
                        if index + 1 >= search_index {
                            for _ in 0..expand {
                                msa.sequences[sequence_index]
                                    .path
                                    .insert(path_index, MSAPosition::Gap);
                            }
                            break;
                        }
                        index += 1;
                    }
                }
            }
        }

        for sequence_index in 0..self.sequences.len() {
            let mut index = self.sequences[sequence_index].start;
            for path_index in 0..self.sequences[sequence_index].path.len() {
                if let MSAPosition::Placed(_, a, b) =
                    self.sequences[sequence_index].path[path_index]
                {
                    if a > b {
                        fix(self, index, a - b, sequence_index);
                    }
                    index += b;
                } else {
                    index += 1;
                }
            }
        }
    }

    /// Generate HTML
    pub fn html(&self) -> String {
        let mut res = String::new();
        write!(
            &mut res,
            "<div class='msa alignment' style='--length:{}'>",
            self.total_length()
        )
        .unwrap();
        for sequence in &self.sequences {
            write!(&mut res, "<div class='peptide' title='path: {} norm: {:.4} abs: {}' data-normalised-score='{1}' style='--start:{};--length:{};--score:{1}'>", 
        sequence
        .path
        .iter()
        .map(ToString::to_string)
        .collect::<String>(),
        sequence.normalised_score,
        sequence.score,
        sequence.start,
        sequence.len(),
    ).unwrap();
            let mut start = sequence.start;
            for piece in &sequence.path {
                match piece {
                    MSAPosition::Gap => write!(&mut res, "<del></del>").unwrap(),
                    MSAPosition::Placed(ty, a, b) => {
                        write!(
                            &mut res,
                            "<span class='{}' style='--wa:{};--wb:{};'><span>{}</span></span>",
                            match ty {
                                MatchType::Switched => "rotation",
                                MatchType::Mismatch => "mismatch",
                                MatchType::Isobaric => "isobaric",
                                MatchType::Gap => "gap",
                                MatchType::FullIdentity => "identity",
                                MatchType::IdentityMassMismatch => "massmismatch",
                            },
                            a,
                            b,
                            render_seq(
                                &sequence.sequence,
                                start..(start + b).min(sequence.sequence.len()),
                                *ty
                            )
                        )
                        .unwrap();
                        start += a;
                    }
                }
            }
            write!(&mut res, "</div>").unwrap();
        }
        write!(&mut res, "</div>").unwrap();
        res
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

fn render_seq(seq: &LinearPeptide, sel: std::ops::Range<usize>, ty: MatchType) -> String {
    use std::fmt::Write;
    let mut res = String::new();
    for s in &seq.sequence[sel] {
        if s.modifications.is_empty() {
            write!(&mut res, "{}", s.aminoacid.char()).unwrap();
        } else {
            write!(
                &mut res,
                "<span class='modified{}'>{}</span>",
                match ty {
                    MatchType::IdentityMassMismatch => " massmismatch",
                    MatchType::Mismatch => " mismatch",
                    _ => "",
                },
                s.aminoacid.char()
            )
            .unwrap();
        }
    }
    res
}

impl Display for MultipleSequenceAlignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for sequence in &self.sequences {
            write!(f, "{}", " ".repeat(sequence.start))?;
            let mut seq = sequence.sequence.sequence.iter().skip(sequence.start);
            for step in &sequence.path {
                match step {
                    MSAPosition::Gap => write!(f, "-")?,
                    MSAPosition::Placed(_, steps, width) => write!(
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
                " [Score: {}, Normalised score: {:.3}, Path: {}]",
                sequence.score,
                sequence.normalised_score,
                sequence
                    .path
                    .iter()
                    .map(ToString::to_string)
                    .collect::<String>()
            )?;
        }
        Ok(())
    }
}

impl MSAPosition {
    pub const fn len(&self) -> usize {
        match self {
            Self::Gap => 1,
            Self::Placed(_, a, _) => *a,
        }
    }
    pub const fn ty(&self) -> MatchType {
        match self {
            Self::Gap => MatchType::Gap,
            Self::Placed(ty, _, _) => *ty,
        }
    }
}

impl Display for MSAPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Gap => "[-]".to_string(),
                Self::Placed(ty, a, b) => format!(
                    "[{}{a},{b}]",
                    match ty {
                        MatchType::FullIdentity => '=',
                        MatchType::IdentityMassMismatch => '≅',
                        MatchType::Isobaric => 'i',
                        MatchType::Switched => 's',
                        MatchType::Mismatch => 'x',
                        MatchType::Gap => 'g',
                    },
                ),
            }
        )
    }
}

impl MassAlignable for MultipleSequenceAlignment {
    fn total_length(&self) -> usize {
        self.sequences
            .iter()
            .map(|s| s.len() + s.start)
            .max()
            .unwrap_or(0)
    }
    fn number_of_sequences(&self) -> usize {
        self.sequences.len()
    }
    fn index(&self, index: usize, sequence_index: usize) -> &crate::SequenceElement {
        &self.sequences[sequence_index].sequence.sequence[index]
    }
    fn index_slice(
        &self,
        index: impl std::ops::RangeBounds<usize>,
        sequence_index: usize,
    ) -> std::borrow::Cow<[crate::SequenceElement]> {
        std::borrow::Cow::Borrowed(
            &self.sequences[sequence_index].sequence.sequence
                [(index.start_bound().cloned(), index.end_bound().cloned())],
        )
    }
    fn sequence_bounds(&self) -> Vec<(usize, usize)> {
        self.sequences
            .iter()
            .map(|s| (s.start, s.len() + s.start))
            .collect()
    }
    fn calculate_masses<const STEPS: usize>(&self) -> Vec<Vec<Vec<Option<crate::system::Mass>>>> {
        todo!();
    }
    fn sequences_with_path(
        &self,
        is_a: bool,
        start: usize,
        path: &[super::Piece],
    ) -> Vec<MSAPlacement> {
        todo!();
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
                        MSAPosition::Placed(MatchType::FullIdentity, 1, 1),
                        MSAPosition::Placed(MatchType::Isobaric, 1, 1),
                        MSAPosition::Placed(MatchType::Isobaric, 1, 1),
                        MSAPosition::Placed(MatchType::FullIdentity, 1, 1),
                    ],
                    score: 0,
                    normalised_score: 0.0,
                },
                MSAPlacement {
                    sequence: ComplexPeptide::pro_forma("WND").unwrap().assume_linear(),
                    start: 0,
                    path: vec![
                        MSAPosition::Placed(MatchType::FullIdentity, 1, 1),
                        MSAPosition::Placed(MatchType::Isobaric, 1, 2),
                        MSAPosition::Placed(MatchType::FullIdentity, 1, 1),
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
