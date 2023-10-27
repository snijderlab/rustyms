use crate::{system::Mass, LinearPeptide, MassTolerance, MolecularFormula, SequenceElement};

use super::{align_type::*, piece::*, scoring::*, Alignment};
use crate::uom::num_traits::Zero;

/// Create an alignment of two peptides based on mass and homology.
/// # Panics
/// It panics when the length of `seq_a` or `seq_b` is bigger then [`isize::MAX`].
#[allow(clippy::too_many_lines)]
pub fn align(
    seq_a: LinearPeptide,
    seq_b: LinearPeptide,
    alphabet: &[&[i8]],
    tolerance: MassTolerance,
    ty: Type,
) -> Alignment {
    const STEPS: usize = 3; // can at max be i8::MAX / 2 => 64
    assert!(isize::try_from(seq_a.len()).is_ok());
    assert!(isize::try_from(seq_b.len()).is_ok());
    let mut matrix = vec![vec![Piece::default(); seq_b.len() + 1]; seq_a.len() + 1];
    let mut high = (0, 0, 0);
    let masses_a = calculate_masses(STEPS, &seq_a);
    let masses_b = calculate_masses(STEPS, &seq_b);

    if ty.global() {
        #[allow(clippy::cast_possible_wrap)]
        // b is always less than seq_b
        for index_b in 0..=seq_b.len() {
            matrix[0][index_b] = Piece::new(
                (index_b as isize) * GAP_EXTEND_PENALTY as isize,
                GAP_EXTEND_PENALTY,
                MatchType::Gap,
                0,
                u8::from(index_b != 0),
            );
        }
    }
    if ty == Type::Global {
        #[allow(clippy::cast_possible_wrap)]
        // a is always less than seq_a
        for (index_a, row) in matrix.iter_mut().enumerate() {
            row[0] = Piece::new(
                (index_a as isize) * GAP_EXTEND_PENALTY as isize,
                GAP_EXTEND_PENALTY,
                MatchType::Gap,
                u8::from(index_a != 0),
                0,
            );
        }
    }

    let mut values = Vec::with_capacity(STEPS * STEPS + 2);
    for index_a in 1..=seq_a.len() {
        for index_b in 1..=seq_b.len() {
            values.clear();
            for len_a in 0..=STEPS {
                for len_b in 0..=STEPS {
                    if len_a == 0 && len_b != 1
                        || len_a != 1 && len_b == 0
                        || len_a > index_a
                        || len_b > index_b
                    {
                        continue; // Do not allow double gaps, any double gaps will be counted as two gaps after each other
                    }
                    let base_score = matrix[index_a - len_a][index_b - len_b].score;
                    // len_a and b are always <= STEPS
                    let piece = if len_a == 0 || len_b == 0 {
                        // First check the score to be used for affine gaps
                        let prev = &matrix[index_a - len_a][index_b - len_b];
                        let score =
                            if prev.step_a == 0 && len_a == 0 || prev.step_b == 0 && len_b == 0 {
                                GAP_EXTEND_PENALTY
                            } else {
                                GAP_START_PENALTY
                            };
                        Some(Piece::new(
                            base_score + score as isize,
                            score,
                            MatchType::Gap,
                            len_a as u8,
                            len_b as u8,
                        ))
                    } else if len_a == 1 && len_b == 1 {
                        Some(score_pair(
                            &seq_a.sequence[index_a - 1],
                            masses_a[0][index_a],
                            &seq_b.sequence[index_b - 1],
                            masses_b[0][index_b],
                            alphabet,
                            base_score,
                            tolerance,
                        ))
                    } else {
                        score(
                            &seq_a.sequence[index_a - len_a..index_a],
                            masses_a[len_a - 1][index_a],
                            &seq_b.sequence[index_b - len_b..index_b],
                            masses_b[len_b - 1][index_b],
                            base_score,
                            tolerance,
                        )
                    };
                    if let Some(p) = piece {
                        values.push(p);
                    }
                }
            }
            let value = values
                .iter()
                .max_by(|x, y| x.score.cmp(&y.score))
                .cloned()
                .unwrap_or_default();
            if value.score >= high.0 {
                high = (value.score, index_a, index_b);
            }
            matrix[index_a][index_b] = value;
        }
    }

    // loop back
    if ty == Type::Global {
        high = (
            matrix[seq_a.len()][seq_b.len()].score,
            seq_a.len(),
            seq_b.len(),
        );
    } else if ty == Type::GlobalForB {
        let value = (0..=seq_a.len())
            .map(|v| (v, matrix[v][seq_b.len()].score))
            .max_by(|a, b| a.1.cmp(&b.1))
            .unwrap_or_default();
        high = (value.1, value.0, seq_b.len());
    }
    let mut path = Vec::new();
    let high_score = high.0;
    while ty == Type::Global || !(high.1 == 0 && high.2 == 0) {
        let value = matrix[high.1][high.2].clone();
        if value.step_a == 0 && value.step_b == 0 {
            break;
        }
        high = (
            0,
            high.1 - value.step_a as usize,
            high.2 - value.step_b as usize,
        );
        path.push(value);
    }

    let path: Vec<Piece> = path.into_iter().rev().collect();
    let max_score = seq_a.sequence
        [high.1..high.1 + path.iter().map(|p| p.step_a as usize).sum::<usize>()]
        .iter()
        .map(|a| alphabet[a.aminoacid as usize][a.aminoacid as usize] as isize)
        .sum::<isize>()
        + seq_b.sequence[high.2..high.2 + path.iter().map(|p| p.step_b as usize).sum::<usize>()]
            .iter()
            .map(|a| alphabet[a.aminoacid as usize][a.aminoacid as usize] as isize)
            .sum::<isize>();
    Alignment {
        absolute_score: high_score,
        normalised_score: high_score as f64 / max_score as f64 * 2.0,
        path,
        start_a: high.1,
        start_b: high.2,
        seq_a,
        seq_b,
        ty,
    }
}

/// Score a pair of sequence elements (AA + mods)
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
                alphabet[a.aminoacid as usize][b.aminoacid as usize] + MASS_MISMATCH_PENALTY;
            Piece::new(
                score + local as isize,
                local,
                MatchType::IdentityMassMismatch,
                1,
                1,
            )
        }
        (false, true) => Piece::new(score + ISOMASS as isize, ISOMASS, MatchType::Isobaric, 1, 1), // TODO: I/L/J is now also scored as isobaric, which is correct but the score is considerably lower then in previous iterations
        (false, false) => Piece::new(
            score + MISMATCH as isize,
            MISMATCH,
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
            SWITCHED * a.len() as i8
        } else {
            ISOMASS * (a.len() + b.len()) as i8 / 2
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

/// Get the masses of all subsets of up to the given number of steps as a lookup table.
/// The result should be is index by [steps-1][index]
fn calculate_masses(steps: usize, sequence: &LinearPeptide) -> Vec<Vec<Mass>> {
    (1..=steps)
        .map(|size| {
            (0..=sequence.len())
                .map(|index| {
                    if index < size {
                        Mass::zero()
                    } else {
                        sequence.sequence[index - size..index]
                            .iter()
                            .map(|s| s.formula_all().unwrap())
                            .sum::<MolecularFormula>()
                            .monoisotopic_mass()
                            .unwrap()
                    }
                })
                .collect::<Vec<Mass>>()
        })
        .collect::<Vec<_>>()
}

#[cfg(test)]
mod tests {
    use super::score;
    use crate::aminoacids::AminoAcid;
    use crate::{MolecularFormula, SequenceElement};

    #[test]
    fn pair() {
        let a = [SequenceElement::new(AminoAcid::N, None)];
        let b = [
            SequenceElement::new(AminoAcid::G, None),
            SequenceElement::new(AminoAcid::G, None),
        ];
        let pair = dbg!(score(
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
}
