//! Functions to generate alignments of peptides based on homology, while taking mass spec error into account.

use std::fmt::Write;

use crate::uom::num_traits::Zero;
use crate::{Mass, MolecularFormula, Peptide, SequenceElement};

/// An alignment of two reads.
#[derive(Debug, Clone)]
pub struct Alignment {
    /// The score of this alignment
    pub score: isize,
    /// The path or steps taken for the alignment
    pub path: Vec<Piece>,
    /// The position in the first sequence where the alignment starts
    pub start_a: usize,
    /// The position in the second sequence where the alignment starts
    pub start_b: usize,
    /// The first sequence
    pub seq_a: Peptide,
    /// The second sequence
    pub seq_b: Peptide,
    /// The alignment type
    pub ty: Type,
}

impl Alignment {
    /// Get a short representation of the alignment in CIGAR like format. It has one additional class `s[a,b]` denoting any special step with the given a and b step size.
    pub fn short(&self) -> String {
        #[derive(PartialEq, Eq)]
        enum StepType {
            Insertion,
            Deletion,
            Match,
            Mismatch,
            Special(MatchType, u8, u8),
        }
        impl std::fmt::Display for StepType {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(
                    f,
                    "{}",
                    match self {
                        Self::Insertion => String::from("I"),
                        Self::Deletion => String::from("D"),
                        Self::Match => String::from("="),
                        Self::Mismatch => String::from("X"),
                        Self::Special(MatchType::Switched, a, b) => format!("s[{a}, {b}]"),
                        Self::Special(MatchType::Isobaric, a, b) => format!("i[{a}, {b}]"),
                        Self::Special(..) => panic!("A special match cannot be of this match type"),
                    }
                )
            }
        }
        let (_, _, output, last) = self.path.iter().fold(
            (self.start_a, self.start_b, String::new(), None),
            |(a, b, output, last), step| {
                let current_type = match (step.match_type, step.step_a, step.step_b) {
                    (MatchType::Isobaric, a, b) => StepType::Special(MatchType::Isobaric, a, b), // Catch any 1/1 isobaric sets before they are counted as Match/Mismatch
                    (_, 0, 1) => StepType::Insertion,
                    (_, 1, 0) => StepType::Deletion,
                    (_, 1, 1) if self.seq_a.sequence[a] == self.seq_b.sequence[b] => {
                        StepType::Match
                    }
                    (_, 1, 1) => StepType::Mismatch,
                    (m, a, b) => StepType::Special(m, a, b),
                };
                let (str, last) = match last {
                    Some((t @ StepType::Special(..), _)) => {
                        (format!("{output}{t}"), Some((current_type, 1)))
                    }
                    Some((t, n)) if t == current_type => (output, Some((t, n + 1))),
                    Some((t, n)) => (format!("{output}{n}{t}"), Some((current_type, 1))),
                    None => (output, Some((current_type, 1))),
                };
                (
                    a + step.step_a as usize,
                    b + step.step_b as usize,
                    str,
                    last,
                )
            },
        );
        match last {
            Some((t @ StepType::Special(..), _)) => format!("{output}{t}"),
            Some((t, n)) => format!("{output}{n}{t}"),
            _ => output,
        }
    }

    /// Get the error in ppm for this match, if it is a (partial) local match it will only take the matched amino acids into account.
    pub fn ppm(&self) -> Option<f64> {
        Some(
            self.mass_a()?
                .monoisotopic_mass()?
                .ppm(self.mass_b()?.monoisotopic_mass()?),
        )
    }

    /// Get the mass delta for this match, if it is a (partial) local match it will only take the matched amino acids into account.
    pub fn mass_difference(&self) -> Option<crate::Mass> {
        Some(self.mass_a()?.monoisotopic_mass()? - self.mass_b()?.monoisotopic_mass()?)
    }

    fn mass_a(&self) -> Option<MolecularFormula> {
        if self.ty == Type::Global {
            self.seq_a.formula()
        } else {
            let mut placed_a = vec![false; self.seq_a.ambiguous_modifications.len()];
            self.seq_a.sequence[self.start_a..self.start_a + self.len_a()]
                .iter()
                .fold(Some(MolecularFormula::default()), |acc, s| {
                    s.formula_greedy(&mut placed_a)
                        .and_then(|m| acc.map(|a| a + m))
                })
        }
    }

    fn mass_b(&self) -> Option<MolecularFormula> {
        if self.ty == Type::Global {
            self.seq_b.formula()
        } else {
            let mut placed_b = vec![false; self.seq_b.ambiguous_modifications.len()];
            self.seq_b.sequence[self.start_b..self.start_b + self.len_b()]
                .iter()
                .fold(Some(MolecularFormula::default()), |acc, s| {
                    s.formula_greedy(&mut placed_b)
                        .and_then(|m| acc.map(|a| a + m))
                })
        }
    }

    /// Returns statistics for this match. Returns `(identical, gap, length)`. Retrieve the identity (or gap) as percentage
    /// by calculating `identical as f64 / length as f64`.
    pub fn stats(&self) -> (usize, usize, usize) {
        let (identical, gap, _, _) =
            self.path
                .iter()
                .fold((0, 0, self.start_a, self.start_b), |acc, p| {
                    if p.step_a + p.step_b == 1 {
                        (
                            acc.0,
                            acc.1 + 1,
                            acc.2 + p.step_a as usize,
                            acc.3 + p.step_b as usize,
                        )
                    } else if p.step_a == 1
                        && p.step_b == 1
                        && self.seq_a.sequence[acc.2] == self.seq_b.sequence[acc.3]
                    {
                        (
                            acc.0 + 1,
                            acc.1,
                            acc.2 + p.step_a as usize,
                            acc.3 + p.step_b as usize,
                        )
                    } else {
                        (
                            acc.0,
                            acc.1,
                            acc.2 + p.step_a as usize,
                            acc.3 + p.step_b as usize,
                        )
                    }
                });
        (identical, gap, self.len_a().max(self.len_b()))
    }

    /// Generate a summary of this alignment for printing to the command line
    pub fn summary(&self) -> String {
        format!(
            "score: {}\npath: {}\nstart: ({}, {})\naligned:\n{}",
            self.score,
            self.short(),
            self.start_a,
            self.start_b,
            self.aligned()
        )
    }

    /// The total number of residues matched on the first sequence
    pub fn len_a(&self) -> usize {
        self.path.iter().map(|p| p.step_a as usize).sum()
    }

    /// The total number of residues matched on the second sequence
    pub fn len_b(&self) -> usize {
        self.path.iter().map(|p| p.step_b as usize).sum()
    }

    // TODO: find a more graceful way of handling B/Z amino acids
    fn aligned(&self) -> String {
        let blocks: Vec<char> = " ▁▂▃▄▅▆▇█".chars().collect();
        let blocks_neg: Vec<char> = "▔▔▔▔▀▀▀▀█".chars().collect();
        let mut str_a = String::new();
        let mut str_b = String::new();
        let mut str_blocks = String::new();
        let mut str_blocks_neg = String::new();
        let mut loc_a = self.start_a;
        let mut loc_b = self.start_b;
        let max = self
            .path
            .iter()
            .map(|p| p.local_score)
            .max()
            .unwrap_or(i8::MAX);
        let min = self
            .path
            .iter()
            .map(|p| p.local_score)
            .min()
            .unwrap_or(i8::MIN + 1); // +1 to make it also valid as a positive number
        let factor = blocks.len() as f64 / f64::from(min.abs().max(max));
        let index = |n| ((f64::from(n) * factor).floor() as usize).min(blocks.len() - 1);

        for piece in &self.path {
            let l = std::cmp::max(piece.step_b, piece.step_a);
            if piece.step_a == 0 {
                write!(str_a, "{:-<width$}", "", width = l as usize).unwrap();
            } else {
                write!(
                    str_a,
                    "{:·<width$}",
                    self.seq_a.sequence[loc_a..loc_a + piece.step_a as usize]
                        .iter()
                        .map(|a| a.aminoacid.char())
                        .collect::<String>(),
                    width = l as usize
                )
                .unwrap();
            }
            if piece.step_b == 0 {
                write!(str_b, "{:-<width$}", "", width = l as usize).unwrap();
            } else {
                write!(
                    str_b,
                    "{:·<width$}",
                    self.seq_b.sequence[loc_b..loc_b + piece.step_b as usize]
                        .iter()
                        .map(|a| a.aminoacid.char())
                        .collect::<String>(),
                    width = l as usize
                )
                .unwrap();
            }
            write!(
                str_blocks,
                "{}",
                str::repeat(
                    &if piece.local_score < 0 {
                        " ".to_string()
                    } else {
                        #[allow(clippy::cast_sign_loss)] // Checked above
                        blocks[index(piece.local_score)].to_string()
                    },
                    l as usize
                )
            )
            .unwrap();
            write!(
                str_blocks_neg,
                "{}",
                str::repeat(
                    &if piece.local_score > 0 {
                        " ".to_string()
                    } else {
                        #[allow(clippy::cast_sign_loss)] // Checked above
                        blocks_neg[index(-piece.local_score)].to_string()
                    },
                    l as usize
                )
            )
            .unwrap();

            loc_a += piece.step_a as usize;
            loc_b += piece.step_b as usize;
        }

        format!("{str_a}\n{str_b}\n{str_blocks}\n{str_blocks_neg}")
    }
}

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
    pub fn short(&self) -> String {
        match (self.step_a, self.step_b) {
            (0, 1) => "I".to_string(),
            (1, 0) => "D".to_string(),
            (1, 1) => "M".to_string(),
            (a, b) => format!("S[{b},{a}]"),
        }
    }
}

/// The type of alignment to perform
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Type {
    /// Global alignment, which tries to find the best alignment to link both sequences fully to each other, like the Needleman Wunsch algorithm
    Global,
    /// Local alignment, which tries to find the best patch of both sequences to align to each other, this could lead to trailing ends on both sides of both sequences, like the Smith Waterman
    Local,
    /// Hybrid alignment, the second sequence will be fully aligned to the first sequence, this could lead to trailing ends on the first sequence but not on the second.
    GlobalForB,
}

impl Type {
    const fn global(self) -> bool {
        !matches!(self, Self::Local)
    }
}

#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub enum MatchType {
    /// Aminoacid + Mass identity
    FullIdentity,
    /// Aminoacid + Mass mismatch
    IdentityMassMismatch,
    /// Full mismatch
    #[default]
    Mismatch,
    /// Set of aminoacids + mods with the same mass but different sequence
    Isobaric,
    /// Set of aminoacids + mods in a different order in the two sequences
    Switched,
    /// A gap
    Gap,
}

const MISMATCH: i8 = -1;
const MASS_MISMATCH_PENALTY: i8 = -1;
const SWITCHED: i8 = 3;
const ISOMASS: i8 = 2;
const GAP_START_PENALTY: i8 = -5;
const GAP_EXTEND_PENALTY: i8 = -1;

/// Create an alignment of two peptides based on mass and homology.
/// # Panics
/// It panics when the length of `seq_a` or `seq_b` is bigger then [`isize::MAX`].
#[allow(clippy::too_many_lines)]
pub fn align(seq_a: Peptide, seq_b: Peptide, alphabet: &[&[i8]], ty: Type) -> Alignment {
    const STEPS: usize = 3;
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
                        Some(Piece::new(
                            base_score + GAP_EXTEND_PENALTY as isize, // TODO: Check affine gaps
                            GAP_EXTEND_PENALTY,
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
                        ))
                    } else {
                        score(
                            &seq_a.sequence[index_a - len_a..index_a],
                            masses_a[len_a - 1][index_a],
                            &seq_b.sequence[index_b - len_b..index_b],
                            masses_b[len_b - 1][index_b],
                            base_score,
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
    Alignment {
        score: high_score,
        path: path.into_iter().rev().collect(),
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
) -> Piece {
    match (a == b, mass_similar(mass_a, mass_b)) {
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
        (false, true) => Piece::new(score + ISOMASS as isize, ISOMASS, MatchType::Isobaric, 1, 1),
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
) -> Option<Piece> {
    if mass_similar(mass_a, mass_b) {
        let mut b_copy = b.to_owned();
        let switched = a.len() == b.len()
            && a.iter().all(|el| {
                if let Some(pos) = b_copy.iter().position(|x| x == el) {
                    b_copy.remove(pos);
                    true
                } else {
                    false
                }
            });
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

/// Determine if two masses are close enough to be considered similar.
/// This is the case if the two masses are within 10 ppm or 0.1 Da
fn mass_similar(a: Mass, b: Mass) -> bool {
    a.ppm(b) < 10.0 || (a.value - b.value).abs() < 0.1
}

/// Get the masses of all subsets of up to the given number of steps as a lookup table.
/// The result should be is index by [steps-1][index]
fn calculate_masses(steps: usize, sequence: &Peptide) -> Vec<Vec<Mass>> {
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

pub const BLOSUM62: &[&[i8]] = include!("blosum62.txt");

#[cfg(test)]
mod tests {
    use crate::align::score;
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
            0
        ));
        assert!(pair.is_some());
    }
}

//     #[test]
//     fn equal() {
//         let alphabet = Alphabet::default();
//         let a = vec![A, C, C, G, W];
//         let b = vec![A, C, C, G, W];
//         let result = align(&a, &b, &alphabet, Type::Local);
//         dbg!(&result);
//         assert_eq!(40, result.score);
//         assert_eq!("MMMMM", &result.short());
//     }

//     #[test]
//     fn insertion() {
//         let alphabet = Alphabet::default();
//         let a = vec![A, C, G, W];
//         let b = vec![A, C, F, G, W];
//         let result = align(&a, &b, &alphabet, Type::Local);
//         dbg!(&result);
//         assert_eq!(27, result.score);
//         assert_eq!("MMIMM", &result.short());
//     }

//     #[test]
//     fn deletion() {
//         let alphabet = Alphabet::default();
//         let a = vec![A, C, F, G, W];
//         let b = vec![A, C, G, W];
//         let result = align(&a, &b, &alphabet, Type::Local);
//         dbg!(&result);
//         assert_eq!(27, result.score);
//         assert_eq!("MMDMM", &result.short());
//     }

//     #[test]
//     fn iso_mass() {
//         let alphabet = Alphabet::default();
//         let a = vec![A, F, G, G, W];
//         let b = vec![A, F, N, W];
//         let result = align(&a, &b, &alphabet, Type::Local);
//         dbg!(&result);
//         dbg!(result.short());
//         assert_eq!(29, result.score);
//         assert_eq!("MMS[1,2]M", &result.short());
//     }

//     #[test]
//     fn switched() {
//         let alphabet = Alphabet::default();
//         let a = vec![A, F, G, G, W];
//         let b = vec![A, G, F, G, W];
//         let result = align(&a, &b, &alphabet, Type::Local);
//         dbg!(&result);
//         dbg!(result.short());
//         assert_eq!(28, result.score);
//         assert_eq!("MS[2,2]MM", &result.short());
//     }

//     #[test]
//     fn local() {
//         let alphabet = Alphabet::default();
//         let a = vec![A, F, G, G, E, W];
//         let b = vec![F, G, G, D];
//         let result = align(&a, &b, &alphabet, Type::Local);
//         dbg!(&result);
//         dbg!(result.short());
//         assert_eq!(24, result.score);
//         assert_eq!("MMM", &result.short());
//     }

//     #[test]
//     fn global() {
//         let alphabet = Alphabet::default();
//         let a = vec![A, F, G, G, E, W];
//         let b = vec![F, G, G, D];
//         let result = align(&a, &b, &alphabet, Type::Global);
//         dbg!(&result);
//         println!("{}", result.summary());
//         assert_eq!(13, result.score);
//         assert_eq!("DMMMDM", &result.short());
//         assert_eq!(0, result.start_a, "A global alignment should start at 0");
//     }

//     #[test]
//     fn global_for_b() {
//         let alphabet = Alphabet::default();
//         let a = vec![A, F, G, G, E, W];
//         let b = vec![F, G, G, D];
//         let result = align(&a, &b, &alphabet, Type::GlobalForB);
//         dbg!(&result);
//         dbg!(result.short());
//         assert_eq!(23, result.score);
//         assert_eq!("MMMM", &result.short());
//         assert_eq!(0, result.start_b, "A global alignment should start at 0");
//     }
// }
