//! Functions to generate alignments of peptides based on homology, while taking mass spec error into account.

use std::fmt::Write;

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use super::align_type::*;
use super::piece::*;
use crate::align::scoring::*;
use crate::system::Mass;
use crate::LinearPeptide;
use crate::MolecularFormula;
use crate::Multi;
use crate::MultiChemical;

/// An alignment of two reads.
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct Alignment {
    /// The absolute score of this alignment
    pub absolute_score: isize,
    /// The maximal score of this alignment: the average score of the sequence slices on sequence a and b if they were aligned to themself, rounded down.
    /// Think of it like this: `align(sequence_a.sequence[start_a..len_a], sequence_a.sequence[start_a..len_a])`.
    pub maximal_score: isize,
    /// The normalised score, normalised for the alignment length and for the used alphabet.
    /// The normalisation is as follows `absolute_score / max_score`.
    pub normalised_score: f64,
    /// The path or steps taken for the alignment
    pub path: Vec<Piece>,
    /// The position in the first sequence where the alignment starts
    pub start_a: usize,
    /// The position in the second sequence where the alignment starts
    pub start_b: usize,
    /// The first sequence
    pub seq_a: LinearPeptide,
    /// The second sequence
    pub seq_b: LinearPeptide,
    /// The alignment type
    pub ty: AlignType,
    /// The maximal step size (the const generic STEPS)
    pub maximal_step: usize,
}

impl Alignment {
    /// Get a short representation of the alignment in CIGAR like format. It has one additional class `{a}(:{b})?(r|i)` denoting any special step with the given a and b step size, if b is not given it is the same as a.
    pub fn short(&self) -> String {
        #[derive(PartialEq, Eq)]
        enum StepType {
            Insertion,
            Deletion,
            Match,
            Mismatch,
            Special(MatchType, u16, u16),
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
                        Self::Special(MatchType::Rotation, a, b) if a == b => format!("{a}r"),
                        Self::Special(MatchType::Rotation, a, b) => format!("{a}:{b}r"),
                        Self::Special(MatchType::Isobaric, a, b) if a == b => format!("{a}i"),
                        Self::Special(MatchType::Isobaric, a, b) => format!("{a}:{b}i"),
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
    /// If there are multiple possible masses for any of the stretches it returns the smallest difference.
    #[allow(clippy::missing_panics_doc)]
    pub fn ppm(&self) -> f64 {
        self.mass_a()
            .iter()
            .cartesian_product(self.mass_b().iter())
            .map(|(a, b)| a.monoisotopic_mass().ppm(b.monoisotopic_mass()))
            .min_by(f64::total_cmp)
            .expect("An empty Multi<MolecularFormula>  was detected")
    }

    /// Get the mass delta for this match, if it is a (partial) local match it will only take the matched amino acids into account.
    /// If there are multiple possible masses for any of the stretches it returns the smallest difference.
    #[allow(clippy::missing_panics_doc)]
    pub fn mass_difference(&self) -> Mass {
        self.mass_a()
            .iter()
            .cartesian_product(self.mass_b().iter())
            .map(|(a, b)| a.monoisotopic_mass() - b.monoisotopic_mass())
            .min_by(|a, b| a.abs().value.total_cmp(&b.abs().value))
            .expect("An empty Multi<MolecularFormula>  was detected")
    }

    fn mass_a(&self) -> Multi<MolecularFormula> {
        if self.ty.left.global_a() && self.ty.right.global_a() {
            self.seq_a.formulas()
        } else {
            let mut placed_a = vec![false; self.seq_a.ambiguous_modifications.len()];
            self.seq_a.sequence[self.start_a..self.start_a + self.len_a()]
                .iter()
                .fold(Multi::default(), |acc, s| {
                    acc * s.formulas_greedy(&mut placed_a)
                })
        }
    }

    fn mass_b(&self) -> Multi<MolecularFormula> {
        if self.ty.left.global_b() && self.ty.right.global_b() {
            self.seq_b.formulas()
        } else {
            let mut placed_b = vec![false; self.seq_b.ambiguous_modifications.len()];
            self.seq_b.sequence[self.start_b..self.start_b + self.len_b()]
                .iter()
                .fold(Multi::default(), |acc, s| {
                    acc * s.formulas_greedy(&mut placed_b)
                })
        }
    }

    /// Returns statistics for this match. Returns `(identical, mass similar, similar, gap, length)`. Retrieve any stat as percentage
    /// by calculating `stat as f64 / length as f64`. The length is calculated as the max length of `len_a` and `len_b`.
    /// Identical is the number of positions that have identical amino acids, mass similar is the number of positions that have the same mass
    /// (identical, isobaric, rotated), similar is the number of positions that have a positive score in the matrix, and gap is the number
    /// of positions that are gaps.
    pub fn stats(&self) -> (usize, usize, usize, usize, usize) {
        let (identical, mass_similar, similar, gap) =
            self.path.iter().fold((0, 0, 0, 0), |acc, p| {
                let m = p.match_type;
                (
                    acc.0
                        + usize::from(
                            m == MatchType::IdentityMassMismatch || m == MatchType::FullIdentity,
                        ) * p.step_a.max(p.step_b) as usize,
                    acc.1
                        + usize::from(
                            m == MatchType::FullIdentity
                                || m == MatchType::Isobaric
                                || m == MatchType::Rotation,
                        ) * p.step_a.max(p.step_b) as usize,
                    acc.0
                        + usize::from(
                            (m == MatchType::IdentityMassMismatch
                                || m == MatchType::FullIdentity
                                || m == MatchType::Mismatch)
                                && p.local_score >= 0,
                        ) * p.step_a.max(p.step_b) as usize,
                    acc.2 + usize::from(m == MatchType::Gap),
                )
            });
        (
            identical,
            mass_similar,
            similar,
            gap,
            self.len_a().max(self.len_b()),
        )
    }

    /// Generate a summary of this alignment for printing to the command line
    pub fn summary(&self) -> String {
        format!(
            "score: {}\npath: {}\nstart: ({}, {})\naligned:\n{}",
            self.absolute_score,
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
            .unwrap_or(isize::MAX);
        let min = self
            .path
            .iter()
            .map(|p| p.local_score)
            .min()
            .unwrap_or(isize::MIN + 1); // +1 to make it also valid as a positive number
        let factor = blocks.len() as f64 / min.abs().max(max) as f64;
        let index = |n| ((n as f64 * factor).floor() as usize).min(blocks.len() - 1);

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
