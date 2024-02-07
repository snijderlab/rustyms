//! Functions to generate alignments of peptides based on homology, while taking mass spec errors into account.

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::Deserialize;
use serde::Serialize;

use super::align_type::*;
use super::piece::*;
use crate::align::scoring::*;
use crate::system::Mass;
use crate::LinearPeptide;
use crate::MolecularFormula;
use crate::Multi;
use crate::MultiChemical;

/// An alignment of two reads. Which has a reference to the sequences.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct RefAlignment<'a> {
    /// The first sequence
    pub(super) seq_a: &'a LinearPeptide,
    /// The second sequence
    pub(super) seq_b: &'a LinearPeptide,
    pub(super) inner: AlignmentInner,
}

impl<'a> RefAlignment<'a> {
    /// Clone the referenced sequences to make an alignment that owns the sequences.
    /// This can be necessary in some context where the references cannot be guaranteed to stay as long as you need the alignment.
    #[must_use]
    pub fn to_owned(self) -> OwnedAlignment {
        OwnedAlignment {
            seq_a: self.seq_a.clone(),
            seq_b: self.seq_b.clone(),
            inner: self.inner,
        }
    }
}

impl<'a> PartialOrd for RefAlignment<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> Ord for RefAlignment<'a> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.inner.cmp(&other.inner)
    }
}

impl<'a> PrivateAlignment for RefAlignment<'a> {
    fn inner(&self) -> &AlignmentInner {
        &self.inner
    }
}

impl<'a> Alignment for RefAlignment<'a> {
    fn seq_a(&self) -> &LinearPeptide {
        self.seq_a
    }
    fn seq_b(&self) -> &LinearPeptide {
        self.seq_b
    }
}

/// An alignment of two reads. Which owns the sequences.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct OwnedAlignment {
    /// The first sequence
    seq_a: LinearPeptide,
    /// The second sequence
    seq_b: LinearPeptide,
    inner: AlignmentInner,
}

impl PartialOrd for OwnedAlignment {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OwnedAlignment {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.inner.cmp(&other.inner)
    }
}

impl PrivateAlignment for OwnedAlignment {
    fn inner(&self) -> &AlignmentInner {
        &self.inner
    }
}

impl Alignment for OwnedAlignment {
    fn seq_a(&self) -> &LinearPeptide {
        &self.seq_a
    }
    fn seq_b(&self) -> &LinearPeptide {
        &self.seq_b
    }
}

trait PrivateAlignment {
    /// Get the shared inner stuff
    fn inner(&self) -> &AlignmentInner;
}

pub trait Alignment: PrivateAlignment {
    /// The first sequence
    fn seq_a(&self) -> &LinearPeptide;
    /// The second sequence
    fn seq_b(&self) -> &LinearPeptide;

    /// The normalised score, normalised for the alignment length and for the used alphabet.
    /// The normalisation is calculated as follows `absolute_score / max_score`.
    fn normalised_score(&self) -> f64 {
        self.inner().normalised_score.0
    }

    /// All three scores for this alignment (normalised, absolute, max).
    /// 1. The normalised score is normalised for the alignment length and for the used alphabet.
    ///    The normalisation is calculated as follows `absolute_score / max_score`.
    /// 2. The absolute score of this alignment
    /// 3. The maximal score of this alignment: the average score of the sequence slices on sequence a and b if they were aligned to themself, rounded down.
    ///    Think of it like this: `align(sequence_a.sequence[start_a..len_a], sequence_a.sequence[start_a..len_a])`.
    fn scores(&self) -> (f64, isize, isize) {
        (
            self.inner().normalised_score.0,
            self.inner().absolute_score,
            self.inner().maximal_score,
        )
    }

    /// The path or steps taken for the alignment
    fn path(&self) -> &[Piece] {
        &self.inner().path
    }

    /// The position in the sequences where the alignment starts (a, b)
    fn start(&self) -> (usize, usize) {
        (self.inner().start_a, self.inner().start_b)
    }

    /// The position in the first sequence where the alignment starts
    fn start_a(&self) -> usize {
        self.inner().start_a
    }

    /// The position in the second sequence where the alignment starts
    fn start_b(&self) -> usize {
        self.inner().start_b
    }

    /// The alignment type
    fn align_type(&self) -> AlignType {
        self.inner().align_type
    }

    /// The maximal step size (the const generic STEPS)
    fn max_step(&self) -> u16 {
        self.inner().maximal_step
    }

    /// The total number of residues matched on the first sequence
    fn len_a(&self) -> usize {
        self.inner().len_a()
    }

    /// The total number of residues matched on the second sequence
    fn len_b(&self) -> usize {
        self.inner().len_b()
    }

    /// Returns statistics for this match. Returns `(identical, mass similar, similar, gap, length)`. Retrieve any stat as percentage
    /// by calculating `stat as f64 / length as f64`. The length is calculated as the max length of `len_a` and `len_b`.
    /// Identical is the number of positions that have identical amino acids, mass similar is the number of positions that have the same mass
    /// (identical, isobaric, rotated), similar is the number of positions that have a positive score in the matrix, and gap is the number
    /// of positions that are gaps.
    fn stats(&self) -> (usize, usize, usize, usize, usize) {
        self.inner().stats()
    }

    /// The mass(es) for the matched portion of the first sequence
    fn mass_a(&self) -> Multi<MolecularFormula> {
        self.inner().mass_a(self.seq_a())
    }

    /// The mass(es) for the matched portion of the second sequence
    fn mass_b(&self) -> Multi<MolecularFormula> {
        self.inner().mass_b(self.seq_b())
    }

    /// Get the mass delta for this match, if it is a (partial) local match it will only take the matched amino acids into account.
    /// If there are multiple possible masses for any of the stretches it returns the smallest difference.
    fn mass_difference(&self) -> Mass {
        self.inner().mass_difference(self.seq_a(), self.seq_b())
    }

    /// Get the error in ppm for this match, if it is a (partial) local match it will only take the matched amino acids into account.
    /// If there are multiple possible masses for any of the stretches it returns the smallest difference.
    fn ppm(&self) -> f64 {
        self.inner().ppm(self.seq_a(), self.seq_b())
    }

    /// Get a short representation of the alignment in CIGAR like format.
    /// It has one additional class `{a}(:{b})?(r|i)` denoting any special step with the given a and b step size, if b is not given it is the same as a.
    fn short(&self) -> String {
        self.inner().short(self.seq_a(), self.seq_b())
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub(super) struct AlignmentInner {
    /// The absolute score of this alignment
    pub absolute_score: isize,
    /// The maximal score of this alignment: the average score of the sequence slices on sequence a and b if they were aligned to themself, rounded down.
    /// Think of it like this: `align(sequence_a.sequence[start_a..len_a], sequence_a.sequence[start_a..len_a])`.
    pub maximal_score: isize,
    /// The normalised score, normalised for the alignment length and for the used alphabet.
    /// The normalisation is as follows `absolute_score / max_score`.
    pub normalised_score: OrderedFloat<f64>,
    /// The path or steps taken for the alignment
    pub path: Vec<Piece>,
    /// The position in the first sequence where the alignment starts
    pub start_a: usize,
    /// The position in the second sequence where the alignment starts
    pub start_b: usize,
    /// The alignment type
    pub align_type: AlignType,
    /// The maximal step size (the const generic STEPS)
    pub maximal_step: u16,
}

impl AlignmentInner {
    /// Get a short representation of the alignment in CIGAR like format. It has one additional class `{a}(:{b})?(r|i)` denoting any special step with the given a and b step size, if b is not given it is the same as a.
    fn short(&self, seq_a: &LinearPeptide, seq_b: &LinearPeptide) -> String {
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
                    (_, 1, 1) if seq_a.sequence[a] == seq_b.sequence[b] => StepType::Match,
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
    fn ppm(&self, seq_a: &LinearPeptide, seq_b: &LinearPeptide) -> f64 {
        self.mass_a(seq_a)
            .iter()
            .cartesian_product(self.mass_b(seq_b).iter())
            .map(|(a, b)| a.monoisotopic_mass().ppm(b.monoisotopic_mass()))
            .min_by(f64::total_cmp)
            .expect("An empty Multi<MolecularFormula>  was detected")
    }

    /// Get the mass delta for this match, if it is a (partial) local match it will only take the matched amino acids into account.
    /// If there are multiple possible masses for any of the stretches it returns the smallest difference.
    #[allow(clippy::missing_panics_doc)]
    fn mass_difference(&self, seq_a: &LinearPeptide, seq_b: &LinearPeptide) -> Mass {
        self.mass_a(seq_a)
            .iter()
            .cartesian_product(self.mass_b(seq_b).iter())
            .map(|(a, b)| a.monoisotopic_mass() - b.monoisotopic_mass())
            .min_by(|a, b| a.abs().value.total_cmp(&b.abs().value))
            .expect("An empty Multi<MolecularFormula>  was detected")
    }

    fn mass_a(&self, seq_a: &LinearPeptide) -> Multi<MolecularFormula> {
        if self.align_type.left.global_a() && self.align_type.right.global_a() {
            seq_a.formulas()
        } else {
            let mut placed_a = vec![false; seq_a.ambiguous_modifications.len()];
            seq_a[self.start_a..self.start_a + self.len_a()]
                .iter()
                .fold(Multi::default(), |acc, s| {
                    acc * s.formulas_greedy(&mut placed_a)
                })
        }
    }

    fn mass_b(&self, seq_b: &LinearPeptide) -> Multi<MolecularFormula> {
        if self.align_type.left.global_b() && self.align_type.right.global_b() {
            seq_b.formulas()
        } else {
            let mut placed_b = vec![false; seq_b.ambiguous_modifications.len()];
            seq_b[self.start_b..self.start_b + self.len_b()]
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
    fn stats(&self) -> (usize, usize, usize, usize, usize) {
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

    /// The total number of residues matched on the first sequence
    fn len_a(&self) -> usize {
        self.path.iter().map(|p| p.step_a as usize).sum()
    }

    /// The total number of residues matched on the second sequence
    fn len_b(&self) -> usize {
        self.path.iter().map(|p| p.step_b as usize).sum()
    }
}

impl PartialOrd for AlignmentInner {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AlignmentInner {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.normalised_score.total_cmp(&other.normalised_score)
    }
}
