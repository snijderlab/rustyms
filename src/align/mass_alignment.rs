use std::borrow::Cow;

use crate::{
    align::multiple_sequence_alignment::MultipleSequenceAlignment, system::Mass, MassTolerance,
    SequenceElement,
};

use super::{align_type::*, multiple_sequence_alignment::MSAPlacement, piece::*, scoring::*};

pub trait MassAlignable {
    const GAP_EXTEND: i8 = -5;
    const GAP_START: i8 = -1;
    const MASS_MISMATCH: i8 = -1;
    const ISOMASS: i8 = 2;
    const MISMATCH: i8 = -1;
    const SWITCHED: i8 = 3;

    /// Index the underlying structure on one location
    fn index(&self, index: usize) -> Cow<'_, [SequenceElement]>;
    /// Get a slice of the underlying structure
    fn index_slice(&self, index: std::ops::Range<usize>) -> Cow<[Cow<[SequenceElement]>]>;

    /// Total length of the structure (number of AA for single, length from start to end for multiple)
    fn total_length(&self) -> usize;
    /// Give the total number of sequences stored in this structure
    fn number_of_sequences(&self) -> usize;
    /// For all sequences in this structure give their bounds (start index, end index)
    fn sequence_bounds(&self) -> Vec<(usize, usize)>;
    /// Calculate all masses for all steps beforehand.
    /// First dimension is the index in the multiple sequences (or just a single sequence if available).
    /// Second dimension is the index into the sequence (offset by the start and the length is the length of this particular sequence, does not have to span to the very end).
    /// Third dimension is the masses for that number of steps - 1 (a step of 1 is at index 0).
    /// Finally the result is None when that selection given you an invalid step (like getting the mass halfway in an aligned aminoacid)
    fn calculate_masses<const STEPS: usize>(&self) -> Vec<Vec<Vec<Option<Mass>>>>;

    /// Get the [`MSAPlacement`] of the underlying sequences given the path resulting from the alignment
    fn sequences_with_path(&self, is_a: bool, start: usize, path: &[Piece]) -> Vec<MSAPlacement>;

    fn align<const STEPS: usize, A: MassAlignable, B: MassAlignable>(
        seq_a: &A,
        seq_b: &B,
        alphabet: &[&[i8]],
        tolerance: MassTolerance,
        ty: AlignmentType,
    ) -> MultipleSequenceAlignment {
        assert!(isize::try_from(seq_a.total_length()).is_ok());
        assert!(isize::try_from(seq_b.total_length()).is_ok());
        let mut matrix =
            vec![vec![Piece::default(); seq_b.total_length() + 1]; seq_a.total_length() + 1];
        let mut high = (0, 0, 0);
        let masses_a = seq_a.calculate_masses::<STEPS>();
        let masses_b = seq_b.calculate_masses::<STEPS>();
        let bounds_a = seq_a.sequence_bounds();
        let bounds_b = seq_b.sequence_bounds();

        if ty.bind_start_a {
            #[allow(clippy::cast_possible_wrap)]
            // b is always less than seq_b
            for index_b in 0..=seq_b.total_length() {
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
        if ty.bind_start_b {
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
        for index_a in 0..seq_a.total_length() {
            for index_b in 0..seq_b.total_length() {
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
                                    masses_a[sequence_index_a][index_a][0].and_then(|mass_a| {
                                        masses_b[sequence_index_b][index_b][0].map(|mass_b| {
                                            Self::score_pair(
                                                &seq_a.index(index_a)[sequence_index_a],
                                                mass_a,
                                                &seq_b.index(index_b)[sequence_index_b],
                                                mass_b,
                                                alphabet,
                                                base_score,
                                                tolerance,
                                            )
                                        })
                                    })
                                } else {
                                    masses_a[sequence_index_a][index_a][len_a - 1].and_then(
                                        |mass_a| {
                                            masses_b[sequence_index_b][index_a][len_b - 1].and_then(
                                                |mass_b| {
                                                    Self::score(
                                                        &seq_a
                                                            .index_slice(index_a - len_a..index_a)
                                                            [sequence_index_a],
                                                        mass_a,
                                                        &seq_b
                                                            .index_slice(index_b - len_b..index_b)
                                                            [sequence_index_a],
                                                        mass_b,
                                                        base_score,
                                                        tolerance,
                                                    )
                                                },
                                            )
                                        },
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
                if ty.is_global() || value.score > 0 {
                    matrix[index_a + 1][index_b + 1] = value;
                }
            }
        }

        // Find the end cell
        let mut target = high;
        if ty.bind_end_a && ty.bind_end_b {
            target = (
                matrix[seq_a.total_length()][seq_b.total_length()].score,
                seq_a.total_length(),
                seq_b.total_length(),
            );
        } else if ty.bind_end_a {
            let value = (0..=seq_a.total_length())
                .map(|v| (v, matrix[v][seq_b.total_length()].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            target = (value.1, value.0, seq_b.total_length());
        } else if ty.bind_end_b {
            let value = (0..=seq_b.total_length())
                .map(|v| (v, matrix[seq_a.total_length()][v].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            target = (value.1, seq_a.total_length(), value.0);
        }

        // Walk the path back
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
        let mut sequences = seq_a.sequences_with_path(true, target.1, &path);
        sequences.extend(seq_b.sequences_with_path(false, target.2, &path));
        MultipleSequenceAlignment { sequences, ty }
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
    use crate::align::align_type::AlignmentType;
    use crate::align::mass_alignment::MassAlignable;
    use crate::align::BLOSUM62;
    use crate::aminoacids::AminoAcid;
    use crate::{ComplexPeptide, LinearPeptide, MolecularFormula, SequenceElement};

    #[test]
    fn pair() {
        let a = [SequenceElement::new(AminoAcid::N, None)];
        let b = [
            SequenceElement::new(AminoAcid::G, None),
            SequenceElement::new(AminoAcid::G, None),
        ];
        let pair = dbg!(LinearPeptide::score(
            &a,
            a.iter()
                .map(|s| s.formula_all().unwrap())
                .sum::<MolecularFormula>()
                .monoisotopic_mass()
                .unwrap(),
            &b,
            b.iter()
                .map(|s| s.formula_all().unwrap())
                .sum::<MolecularFormula>()
                .monoisotopic_mass()
                .unwrap(),
            0,
            crate::MassTolerance::Ppm(10.0)
        ));
        assert!(pair.is_some());
    }

    #[test]
    fn simple_single_case() {
        let a = ComplexPeptide::pro_forma("EVQLN").unwrap().assume_linear();
        let b = ComplexPeptide::pro_forma("VEQLGG").unwrap().assume_linear();
        let alignment = LinearPeptide::align::<3, LinearPeptide, LinearPeptide>(
            &a,
            &b,
            BLOSUM62,
            crate::MassTolerance::Ppm(10.0),
            AlignmentType::global(),
        );
        println!("{alignment}");
        todo!();
    }
}
