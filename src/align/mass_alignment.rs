use std::fmt::Debug;

use crate::{
    system::Mass, AminoAcid, LinearPeptide, MassComparable, MolecularFormula, Multi,
    SequenceElement, Tolerance,
};

use super::{align_type::*, diagonal_array::DiagonalArray, piece::*, scoring::*, Alignment};

/// Create an alignment of two peptides based on mass and homology.
/// The substitution matrix is in the exact same order as the definition of [`AminoAcid`].
/// The [`MassTolerance`] sets the tolerance for two sets of amino acids to be regarded as the same mass.
/// The [`Type`] controls the alignment behaviour, global/local or anything in between.
/// # Panics
/// It panics when the length of `seq_a` or `seq_b` is bigger than [`isize::MAX`].
/// The peptides are assumed to be simple (see [`LinearPeptide::assume_simple`]).
#[allow(clippy::too_many_lines)]
pub fn align<'a, const STEPS: u16>(
    seq_a: &'a LinearPeptide,
    seq_b: &'a LinearPeptide,
    scoring_matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    tolerance: Tolerance,
    ty: AlignType,
) -> Alignment<'a> {
    // Enforce some assumptions
    seq_a.assume_simple();
    seq_b.assume_simple();
    assert!(isize::try_from(seq_a.len()).is_ok());
    assert!(isize::try_from(seq_b.len()).is_ok());

    let mut matrix = Matrix::new(seq_a.len(), seq_b.len());
    let mut global_highest = (0, 0, 0);
    let masses_a: DiagonalArray<Multi<Mass>> = calculate_masses::<STEPS>(seq_a);
    let masses_b: DiagonalArray<Multi<Mass>> = calculate_masses::<STEPS>(seq_b);
    let zero: Multi<Mass> = Multi::default();

    if ty.left.global_a() {
        matrix.global_start(true);
    }
    if ty.left.global_b() {
        matrix.global_start(false);
    }

    for index_a in 1..=seq_a.len() {
        for index_b in 1..=seq_b.len() {
            let mut highest = None;
            for len_a in 0..=index_a.min(STEPS as usize) {
                for len_b in 0..=index_b.min(STEPS as usize) {
                    if len_a == 0 && len_b != 1
                        || len_a != 1 && len_b == 0
                        || len_a == 0 && len_b == 0
                    {
                        continue; // Do not allow double gaps, any double gaps will be counted as two gaps after each other
                    }
                    let prev = unsafe { matrix.get_unchecked([index_a - len_a, index_b - len_b]) };
                    let base_score = prev.score;

                    // len_a and b are always <= STEPS
                    let piece = if len_a == 0 || len_b == 0 {
                        // First check the score to be used for affine gaps
                        let score = GAP_EXTEND_PENALTY
                            + GAP_START_PENALTY
                                * isize::from(
                                    prev.step_a == 0 && len_a == 0
                                        || prev.step_b == 0 && len_b == 0,
                                );
                        Some(Piece::new(
                            base_score + score,
                            score,
                            MatchType::Gap,
                            len_a as u16,
                            len_b as u16,
                        ))
                    } else if len_a == 1 && len_b == 1 {
                        Some(score_pair(
                            unsafe {
                                (
                                    &seq_a.sequence.get_unchecked(index_a - 1),
                                    &masses_a.get_unchecked([index_a - 1, 0]),
                                )
                            },
                            unsafe {
                                (
                                    &seq_b.sequence.get_unchecked(index_b - 1),
                                    &masses_b.get_unchecked([index_b - 1, 0]),
                                )
                            },
                            scoring_matrix,
                            base_score,
                            tolerance,
                        ))
                    } else {
                        score(
                            unsafe {
                                (
                                    &seq_a.sequence.get_unchecked(
                                        (index_a - len_a).saturating_sub(1)..index_a - 1,
                                    ),
                                    if len_a == 0 {
                                        &zero
                                    } else {
                                        &masses_a.get_unchecked([index_a - 1, len_a - 1])
                                    },
                                )
                            },
                            unsafe {
                                (
                                    &seq_b.sequence.get_unchecked(
                                        (index_b - len_b).saturating_sub(1)..index_b - 1,
                                    ),
                                    if len_b == 0 {
                                        &zero
                                    } else {
                                        &masses_b.get_unchecked([index_b - 1, len_b - 1])
                                    },
                                )
                            },
                            base_score,
                            tolerance,
                        )
                    };
                    if let Some(p) = piece {
                        if highest.is_none()
                            || highest.as_ref().is_some_and(|h: &Piece| h.score < p.score)
                        {
                            highest = Some(p);
                        }
                    }
                }
            }
            if let Some(highest) = highest {
                if highest.score >= global_highest.0 {
                    global_highest = (highest.score, index_a, index_b);
                }
                if ty.left.global() || highest.score > 0 {
                    unsafe {
                        *matrix.get_unchecked_mut([index_a, index_b]) = highest;
                    }
                }
            } else if ty.left.global() {
                unsafe {
                    *matrix.get_unchecked_mut([index_a, index_b]) = score_pair(
                        (
                            &seq_a.sequence.get_unchecked(index_a - 1),
                            &masses_a.get_unchecked([index_a - 1, 0]),
                        ),
                        (
                            &seq_b.sequence.get_unchecked(index_b - 1),
                            &masses_b.get_unchecked([index_b - 1, 0]),
                        ),
                        scoring_matrix,
                        matrix.get_unchecked([index_a - 1, index_b - 1]).score,
                        tolerance,
                    );
                }
            }
        }
    }
    let (absolute_score, start_a, start_b, path) = matrix.trace_path(ty, global_highest);

    let maximal_score = (seq_a.sequence
        [start_a..start_a + path.iter().map(|p| p.step_a as usize).sum::<usize>()]
        .iter()
        .map(|a| scoring_matrix[a.aminoacid as usize][a.aminoacid as usize] as isize)
        .sum::<isize>()
        + seq_b.sequence[start_b..start_b + path.iter().map(|p| p.step_b as usize).sum::<usize>()]
            .iter()
            .map(|a| scoring_matrix[a.aminoacid as usize][a.aminoacid as usize] as isize)
            .sum::<isize>())
        / 2;

    Alignment {
        absolute_score,
        normalised_score: absolute_score as f64 / maximal_score as f64,
        maximal_score,
        path,
        start_a,
        start_b,
        seq_a,
        seq_b,
        ty,
        maximal_step: STEPS as usize,
    }
}

/// Score a pair of sequence elements (AA + mods)
fn score_pair(
    a: (&SequenceElement, &Multi<Mass>),
    b: (&SequenceElement, &Multi<Mass>),
    alphabet: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    score: isize,
    tolerance: Tolerance,
) -> Piece {
    match (a.0 == b.0, tolerance.within(a.1, b.1)) {
        (true, true) => {
            let local = alphabet[a.0.aminoacid as usize][b.0.aminoacid as usize] as isize;
            Piece::new(score + local, local, MatchType::FullIdentity, 1, 1)
        }
        (true, false) => {
            let local = alphabet[a.0.aminoacid as usize][b.0.aminoacid as usize] as isize
                + MASS_MISMATCH_PENALTY;
            Piece::new(score + local, local, MatchType::IdentityMassMismatch, 1, 1)
        }
        (false, true) => {
            // println!("isobaric: {:?} vs {:?}", a.1, b.1);
            Piece::new(score + ISOBARIC, ISOBARIC, MatchType::Isobaric, 1, 1)
        }
        (false, false) => Piece::new(score + MISMATCH, MISMATCH, MatchType::Mismatch, 1, 1),
    }
}

/// Score two sets of aminoacids (it will only be called when at least one of a and b has len > 1)
/// Returns none if no sensible explanation can be made
fn score(
    a: (&[SequenceElement], &Multi<Mass>),
    b: (&[SequenceElement], &Multi<Mass>),
    score: isize,
    tolerance: Tolerance,
) -> Option<Piece> {
    if tolerance.within(a.1, b.1) {
        let rotated = {
            a.0.len() == b.0.len() && {
                let mut b_copy = vec![false; b.0.len()];
                a.0.iter().all(|el| {
                    b_copy
                        .iter()
                        .enumerate()
                        .position(|(index, used)| !used && b.0[index] == *el)
                        .map_or(false, |pos| {
                            b_copy[pos] = true;
                            true
                        })
                })
            }
        };
        #[allow(clippy::cast_possible_wrap)]
        let local = if rotated {
            // println!("rotated: {:?} vs {:?}", a.1, b.1);
            BASE_SPECIAL + ROTATED * a.0.len() as isize
        } else {
            // println!("isobaric: {:?} vs {:?}", a.1, b.1);
            BASE_SPECIAL + ISOBARIC * (a.0.len() + b.0.len()) as isize / 2
        };
        Some(Piece::new(
            score + local as isize,
            local,
            if rotated {
                MatchType::Rotation
            } else {
                MatchType::Isobaric
            },
            a.0.len() as u16,
            b.0.len() as u16,
        ))
    } else {
        None
    }
}

/// Get the masses of all sequence elements
fn calculate_masses<const STEPS: u16>(sequence: &LinearPeptide) -> DiagonalArray<Multi<Mass>> {
    let mut array = DiagonalArray::new(sequence.len(), STEPS);
    // dbg!(&array, format!("{sequence}"));
    for i in 0..sequence.len() {
        // dbg!(i, 0..=i.min(max_depth));
        for j in 0..=i.min(STEPS as usize) {
            array[[i, j]] = sequence.sequence[i - j..=i]
                .iter()
                .map(SequenceElement::formulas_all)
                .sum::<Multi<MolecularFormula>>()
                .iter()
                .map(MolecularFormula::monoisotopic_mass)
                .collect();
        }
    }
    array
}

struct Matrix {
    value: Vec<Vec<Piece>>,
    a: usize,
    b: usize,
}

impl Debug for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use std::fmt::Write;
        for column in &self.value {
            let mut line_0 = String::new();
            let mut line_1 = String::new();
            for cell in column {
                let top = format!("{}/{} {:2}", cell.step_a, cell.step_b, cell.local_score);
                let bottom = format!(
                    "{} {:3}",
                    match cell.match_type {
                        MatchType::FullIdentity => "FI",
                        MatchType::Gap => "G ",
                        MatchType::IdentityMassMismatch => "IM",
                        MatchType::Isobaric => "I ",
                        MatchType::Rotation => "R ",
                        MatchType::Mismatch => "M ",
                    },
                    cell.score
                );
                write!(&mut line_0, "⎡{top:0$}⎤", top.len().max(bottom.len()))?;
                write!(&mut line_1, "⎣{bottom:0$}⎦", top.len().max(bottom.len()))?;
            }
            writeln!(f, "{line_0}")?;
            writeln!(f, "{line_1}")?;
        }
        Ok(())
    }
}

impl Matrix {
    pub fn new(a: usize, b: usize) -> Self {
        Self {
            value: vec![vec![Piece::default(); b + 1]; a + 1],
            a,
            b,
        }
    }

    #[allow(clippy::cast_possible_wrap)]
    pub fn global_start(&mut self, is_a: bool) {
        let max = if is_a { self.a } else { self.b };
        for index in 0..=max {
            self.value[if is_a { index } else { 0 }][if is_a { 0 } else { index }] = Piece::new(
                (index as isize) * GAP_EXTEND_PENALTY,
                GAP_EXTEND_PENALTY,
                MatchType::Gap,
                if is_a { u16::from(index != 0) } else { 0 },
                if is_a { 0 } else { u16::from(index != 0) },
            );
        }
    }

    pub fn trace_path(
        &self,
        ty: AlignType,
        high: (isize, usize, usize),
    ) -> (isize, usize, usize, Vec<Piece>) {
        let mut path = Vec::new();
        let mut high = self.find_end(ty, high);
        let high_score = high.0;

        // Loop back to left side
        while ty.left.global() || !(high.1 == 0 && high.2 == 0) {
            let value = self.value[high.1][high.2].clone();
            if value.step_a == 0 && value.step_b == 0 || !ty.left.global() && value.score < 0 {
                break;
            }
            high = (
                0,
                high.1 - value.step_a as usize,
                high.2 - value.step_b as usize,
            );
            path.push(value);
        }
        (high_score, high.1, high.2, path.into_iter().rev().collect())
    }

    fn find_end(&self, ty: AlignType, high: (isize, usize, usize)) -> (isize, usize, usize) {
        if ty.right.global_a() && ty.right.global_a() {
            (self.value[self.a][self.b].score, self.a, self.b)
        } else if ty.right.global_b() {
            let value = (0..=self.a)
                .map(|v| (v, self.value[v][self.b].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            (value.1, value.0, self.b)
        } else if ty.right.global_a() {
            let value = (0..=self.b)
                .map(|v| (v, self.value[self.a][v].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            (value.1, self.a, value.0)
        } else if ty.right.global() {
            let value_a = (0..=self.a)
                .map(|v| (v, self.value[v][self.b].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            let value_b = (0..=self.b)
                .map(|v| (v, self.value[self.a][v].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            if value_a.1 >= value_b.1 {
                (value_a.1, value_a.0, self.b)
            } else {
                (value_b.1, self.a, value_b.0)
            }
        } else {
            high
        }
    }

    pub unsafe fn get_unchecked(&self, index: [usize; 2]) -> &Piece {
        self.value.get_unchecked(index[0]).get_unchecked(index[1])
    }

    pub unsafe fn get_unchecked_mut(&mut self, index: [usize; 2]) -> &mut Piece {
        self.value
            .get_unchecked_mut(index[0])
            .get_unchecked_mut(index[1])
    }
}

impl std::ops::Index<[usize; 2]> for Matrix {
    type Output = Piece;
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        assert!(index[0] <= self.a + 1);
        assert!(index[1] <= self.b + 1);
        &self.value[index[0]][index[1]]
    }
}
impl std::ops::IndexMut<[usize; 2]> for Matrix {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        assert!(index[0] <= self.a + 1);
        assert!(index[1] <= self.b + 1);
        &mut self.value[index[0]][index[1]]
    }
}

#[cfg(test)]
mod tests {
    use super::score;
    use crate::aminoacids::AminoAcid;
    use crate::{MolecularFormula, Multi, SequenceElement};

    #[test]
    fn pair() {
        let a = [SequenceElement::new(AminoAcid::N, None)];
        let b = [
            SequenceElement::new(AminoAcid::G, None),
            SequenceElement::new(AminoAcid::G, None),
        ];
        let pair = dbg!(score(
            (
                &a,
                &a.iter()
                    .map(SequenceElement::formulas_all)
                    .sum::<Multi<MolecularFormula>>()[0]
                    .monoisotopic_mass()
                    .into()
            ),
            (
                &b,
                &b.iter()
                    .map(SequenceElement::formulas_all)
                    .sum::<Multi<MolecularFormula>>()[0]
                    .monoisotopic_mass()
                    .into()
            ),
            0,
            crate::Tolerance::new_ppm(10.0)
        ));
        assert!(pair.is_some());
    }
}
