use std::{borrow::Cow, fmt::Display};

use super::align_type::*;
use crate::{
    align::{MatchType, Piece},
    system::Mass,
    LinearPeptide, MassTolerance, SequenceElement,
};

/// An alignment of multiple peptides
#[derive(Debug, Clone, PartialEq)]
pub struct MultipleSequenceAlignment {
    /// The sequences
    pub sequences: Vec<MSAPlacement>,
    /// The alignment type
    pub ty: AlignmentType,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentType {
    pub global: bool,
    pub bind_start_a: bool,
    pub bind_start_b: bool,
    pub bind_end_a: bool,
    pub bind_end_b: bool,
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
                        "·".repeat(width - steps)
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

trait MassAlignable<'a>
where
    Self: std::ops::Index<usize, Output = Cow<'a, [SequenceElement]>>
        + std::ops::Index<std::ops::Range<usize>, Output = Cow<'a, [Cow<'a, [SequenceElement]>]>>,
{
    const GAP_EXTEND: i8 = -5;
    const GAP_START: i8 = -1;
    const MASS_MISMATCH: i8 = -1;
    const ISOMASS: i8 = 2;
    const MISMATCH: i8 = -1;
    const SWITCHED: i8 = 3;

    /// Total length of the structure (number of AA for single, length from start to end for multiple)
    fn len(&self) -> usize;
    /// Give the total number of sequences stored in this structure
    fn number_of_sequences(&self) -> usize;
    /// For all sequences in this structure give their bounds (start index, end index)
    fn sequence_bounds(&self) -> Vec<(usize, usize)>;
    /// Calculate all masses for all steps beforehand.
    /// First dimension is the index in the multiple sequences (or just a single sequence if available).
    /// Second dimension is the index into the sequence (offset by the start and the length is the length of this particular sequence, does not have to span to the very end).
    /// Third dimension is the masses for that number of steps - 1 (a step of 1 is at index 0).
    fn calculate_masses<const STEPS: usize>(&self) -> Vec<Vec<Vec<Mass>>>;

    fn sequences_with_path(&self, path: &[Piece]) -> Vec<MSAPlacement>;

    fn align<
        const STEPS: usize,
        const GLOBAL: bool,
        const BIND_START_A: bool,
        const BIND_START_B: bool,
        const BIND_END_A: bool,
        const BIND_END_B: bool,
        A: MassAlignable<'a>,
        B: MassAlignable<'a>,
    >(
        seq_a: A,
        seq_b: B,
        alphabet: &[&[i8]],
        tolerance: MassTolerance,
    ) -> MultipleSequenceAlignment {
        assert!(isize::try_from(seq_a.len()).is_ok());
        assert!(isize::try_from(seq_b.len()).is_ok());
        let mut matrix = vec![vec![Piece::default(); seq_b.len() + 1]; seq_a.len() + 1];
        let mut high = (0, 0, 0);
        let masses_a = seq_a.calculate_masses::<STEPS>();
        let masses_b = seq_b.calculate_masses::<STEPS>();
        let bounds_a = seq_a.sequence_bounds();
        let bounds_b = seq_b.sequence_bounds();

        if BIND_START_A {
            #[allow(clippy::cast_possible_wrap)]
            // b is always less than seq_b
            for index_b in 0..=seq_b.len() {
                matrix[0][index_b] = Piece::new(
                    (index_b as isize).saturating_sub(1) * Self::GAP_EXTEND as isize
                        + Self::GAP_START as isize,
                    if index_b == 0 {
                        Self::GAP_START
                    } else {
                        Self::GAP_EXTEND
                    },
                    MatchType::Gap,
                    0,
                    u8::from(index_b != 0),
                );
            }
        }
        if BIND_START_B {
            #[allow(clippy::cast_possible_wrap)]
            // a is always less than seq_a
            for (index_a, row) in matrix.iter_mut().enumerate() {
                row[0] = Piece::new(
                    (index_a as isize).saturating_sub(1) * Self::GAP_EXTEND as isize
                        + Self::GAP_START as isize,
                    if index_a == 0 {
                        Self::GAP_START
                    } else {
                        Self::GAP_EXTEND
                    },
                    MatchType::Gap,
                    u8::from(index_a != 0),
                    0,
                );
            }
        }

        // Main loop
        let mut values = Vec::with_capacity(
            STEPS * STEPS * seq_a.number_of_sequences() * seq_b.number_of_sequences() + 2,
        );
        for index_a in 0..seq_a.len() {
            for index_b in 0..seq_b.len() {
                values.clear();
                for sequence_index_a in 0..seq_a.number_of_sequences() {
                    // If this sequence does not have info on this location skip it
                    if index_a < bounds_a[sequence_index_a].0
                        || index_a > bounds_a[sequence_index_a].1
                    {
                        continue;
                    }
                    for sequence_index_b in 0..seq_b.number_of_sequences() {
                        // If this sequence does not have info on this location skip it
                        if index_b < bounds_b[sequence_index_b].0
                            || index_b > bounds_b[sequence_index_b].1
                        {
                            continue;
                        }
                        for len_a in 0..=STEPS {
                            for len_b in 0..=STEPS {
                                if len_a == 0 && len_b != 1
                                    || len_a != 1 && len_b == 0
                                    || len_a > index_a
                                    || len_b > index_b
                                {
                                    continue; // Do not allow double gaps, any double gaps will be counted as two gaps after each other
                                }
                                let base_score =
                                    matrix[index_a - len_a + 1][index_b - len_b + 1].score;
                                let piece = if len_a == 0 || len_b == 0 {
                                    // First check the score to be used for affine gaps
                                    let prev = &matrix[index_a - len_a + 1][index_b - len_b + 1];
                                    let score = if prev.step_a == 0 && len_a == 0
                                        || prev.step_b == 0 && len_b == 0
                                    {
                                        Self::GAP_EXTEND
                                    } else {
                                        Self::GAP_START
                                    };
                                    Some(Piece::new(
                                        base_score + score as isize,
                                        score,
                                        MatchType::Gap,
                                        len_a as u8,
                                        len_b as u8,
                                    ))
                                } else if len_a == 1 && len_b == 1 {
                                    Some(Self::score_pair(
                                        &seq_a[index_a][sequence_index_a],
                                        masses_a[sequence_index_a][index_a][0],
                                        &seq_b[index_b][sequence_index_b],
                                        masses_b[sequence_index_b][index_b][0],
                                        alphabet,
                                        base_score,
                                        tolerance,
                                    ))
                                } else {
                                    Self::score(
                                        &seq_a[index_a - len_a..index_a][sequence_index_a],
                                        masses_a[sequence_index_a][index_a][len_a - 1],
                                        &seq_b[index_b - len_b..index_b][sequence_index_a],
                                        masses_b[sequence_index_b][index_a][len_b - 1],
                                        base_score,
                                        tolerance,
                                    )
                                };
                                if let Some(p) = piece {
                                    values.push(p);
                                }
                            }
                        }
                    }
                }
                // Determine the best step
                let value = values
                    .iter()
                    .max_by(|x, y| x.score.cmp(&y.score))
                    .cloned()
                    .unwrap_or_default();
                // Keep track of the highest scoring cell
                if value.score >= high.0 {
                    high = (value.score, index_a + 1, index_b + 1);
                }
                // If local and the score is too low just do not store the result
                if GLOBAL || value.score > 0 {
                    matrix[index_a + 1][index_b + 1] = value;
                }
            }
        }

        // loop back
        let mut target = high;
        if BIND_END_A && BIND_END_B {
            target = (
                matrix[seq_a.len()][seq_b.len()].score,
                seq_a.len(),
                seq_b.len(),
            );
        } else if BIND_END_A {
            let value = (0..=seq_a.len())
                .map(|v| (v, matrix[v][seq_b.len()].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            target = (value.1, value.0, seq_b.len());
        } else if BIND_END_B {
            let value = (0..=seq_b.len())
                .map(|v| (v, matrix[seq_a.len()][v].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            target = (value.1, seq_a.len(), value.0);
        }
        let mut path = Vec::new();
        while !(target.1 == 0 && target.2 == 0) {
            let value = matrix[target.1][target.2].clone();
            if value.step_a == 0 && value.step_b == 0 {
                break;
            }
            target = (
                0,
                target.1 - value.step_a as usize,
                target.2 - value.step_b as usize,
            );
            path.push(value);
        }

        let path: Vec<Piece> = path.into_iter().rev().collect();
        let mut sequences = seq_a.sequences_with_path(&path);
        sequences.extend(seq_b.sequences_with_path(&path));
        MultipleSequenceAlignment {
            sequences,
            ty: AlignmentType {
                global: GLOBAL,
                bind_start_a: BIND_START_A,
                bind_start_b: BIND_START_B,
                bind_end_a: BIND_END_A,
                bind_end_b: BIND_END_B,
            },
        }
    }

    fn score_pair(
        a: &SequenceElement,
        mass_a: Mass,
        b: &SequenceElement,
        mass_b: Mass,
        alphabet: &[&[i8]],
        score: isize,
        tolerance: MassTolerance,
    ) -> Piece {
        match (a == b, tolerance.within(mass_a, mass_b)) {
            (true, true) => {
                let local = alphabet[a.aminoacid as usize][b.aminoacid as usize];
                Piece::new(score + local as isize, local, MatchType::FullIdentity, 1, 1)
            }
            (true, false) => {
                let local =
                    alphabet[a.aminoacid as usize][b.aminoacid as usize] + Self::MASS_MISMATCH;
                Piece::new(
                    score + local as isize,
                    local,
                    MatchType::IdentityMassMismatch,
                    1,
                    1,
                )
            }
            (false, true) => Piece::new(
                score + Self::ISOMASS as isize,
                Self::ISOMASS,
                MatchType::Isobaric,
                1,
                1,
            ), // TODO: I/L/J is now also scored as isobaric, which is correct but the score is considerably lower then in previous iterations
            (false, false) => Piece::new(
                score + Self::MISMATCH as isize,
                Self::MISMATCH,
                MatchType::Mismatch,
                1,
                1,
            ),
        }
    }

    /// Score two sets of aminoacids (it will only be called when at least one of a and b has len > 1)
    /// Returns none if no sensible explanation can be made
    fn score(
        a: &[SequenceElement],
        mass_a: Mass,
        b: &[SequenceElement],
        mass_b: Mass,
        score: isize,
        tolerance: MassTolerance,
    ) -> Option<Piece> {
        if tolerance.within(mass_a, mass_b) {
            let mut b_copy = b.to_owned();
            let switched = a.len() == b.len()
                && a.iter().all(|el| {
                    b_copy.iter().position(|x| x == el).map_or(false, |pos| {
                        b_copy.remove(pos);
                        true
                    })
                });
            #[allow(clippy::cast_possible_wrap)]
            let local = if switched {
                Self::SWITCHED * a.len() as i8
            } else {
                Self::ISOMASS * (a.len() + b.len()) as i8 / 2
            };
            Some(Piece::new(
                score + local as isize,
                local,
                if switched {
                    MatchType::Switched
                } else {
                    MatchType::Isobaric
                },
                a.len() as u8,
                b.len() as u8,
            ))
        } else {
            None
        }
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
                global: true,
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
