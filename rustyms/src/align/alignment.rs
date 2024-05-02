//! Functions to generate alignments of peptides based on homology, while taking mass spectrometry errors into account.

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::Deserialize;
use serde::Serialize;

use super::align_type::*;
use super::piece::*;
use super::scoring::*;

use crate::peptide_complexity::Linear;
use crate::peptide_complexity::Simple;
use crate::system::Mass;
use crate::system::Ratio;
use crate::LinearPeptide;
use crate::MolecularFormula;
use crate::Multi;

/// An alignment of two reads. Which has a reference to the sequences.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct RefAlignment<'a, A, B> {
    /// The first sequence
    pub(super) seq_a: &'a LinearPeptide<A>,
    /// The second sequence
    pub(super) seq_b: &'a LinearPeptide<B>,
    pub(super) inner: AlignmentInner,
}

impl<'a, A, B> RefAlignment<'a, A, B>
where
    LinearPeptide<A>: Into<LinearPeptide<Simple>> + Clone,
    LinearPeptide<B>: Into<LinearPeptide<Simple>> + Clone,
{
    /// Clone the referenced sequences to make an alignment that owns the sequences.
    /// This can be necessary in some context where the references cannot be guaranteed to stay as long as you need the alignment.
    #[must_use]
    pub fn to_owned(&self) -> OwnedAlignment {
        OwnedAlignment {
            seq_a: self.seq_a.clone().into(),
            seq_b: self.seq_b.clone().into(),
            inner: self.inner.clone(),
        }
    }
}

impl<'a, A: Ord, B: Ord> PartialOrd for RefAlignment<'a, A, B>
where
    LinearPeptide<A>: PartialOrd<LinearPeptide<B>>,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a, A: Ord, B: Ord> Ord for RefAlignment<'a, A, B>
where
    LinearPeptide<A>: PartialOrd<LinearPeptide<B>>,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.inner.cmp(&other.inner)
    }
}

impl<'a, A, B> PrivateAlignment for RefAlignment<'a, A, B> {
    fn inner(&self) -> &AlignmentInner {
        &self.inner
    }
}

impl<'a, A: Into<Simple> + Into<Linear>, B: Into<Simple> + Into<Linear>> Alignment
    for RefAlignment<'a, A, B>
{
    fn seq_a(&self) -> &LinearPeptide<impl Into<Simple> + Into<Linear>> {
        self.seq_a
    }
    fn seq_b(&self) -> &LinearPeptide<impl Into<Simple> + Into<Linear>> {
        self.seq_b
    }
}

/// An alignment of two reads. Which owns the sequences.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct OwnedAlignment {
    /// The first sequence
    seq_a: LinearPeptide<Simple>,
    /// The second sequence
    seq_b: LinearPeptide<Simple>,
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
    fn seq_a(&self) -> &LinearPeptide<impl Into<Simple> + Into<Linear>> {
        &self.seq_a
    }
    fn seq_b(&self) -> &LinearPeptide<impl Into<Simple> + Into<Linear>> {
        &self.seq_b
    }
}

/// The link to the inner data structure shared between all alignments, but private to prevent inner details from leaking to the public
trait PrivateAlignment {
    /// Get the shared inner stuff
    fn inner(&self) -> &AlignmentInner;
}

/// A generalised alignment with all behaviour.
#[allow(private_bounds)] // Intended behaviour no one should build on the inner structure
pub trait Alignment: PrivateAlignment {
    /// The first sequence
    fn seq_a(&self) -> &LinearPeptide<impl Into<Simple> + Into<Linear>>;
    /// The second sequence
    fn seq_b(&self) -> &LinearPeptide<impl Into<Simple> + Into<Linear>>;

    /// The normalised score, normalised for the alignment length and for the used alphabet.
    /// The normalisation is calculated as follows `absolute_score / max_score`.
    fn normalised_score(&self) -> f64 {
        self.inner().score.normalised.0
    }

    /// All three scores for this alignment.
    fn score(&self) -> Score {
        self.inner().score
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
        self.path().iter().map(|p| p.step_a as usize).sum()
    }

    /// The total number of residues matched on the second sequence
    fn len_b(&self) -> usize {
        self.path().iter().map(|p| p.step_b as usize).sum()
    }

    /// Returns statistics for this match.
    fn stats(&self) -> Stats {
        let (identical, mass_similar, similar, gaps, length) =
            self.path().iter().fold((0, 0, 0, 0, 0), |acc, p| {
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
                    acc.2
                        + usize::from(
                            (m == MatchType::IdentityMassMismatch
                                || m == MatchType::FullIdentity
                                || m == MatchType::Mismatch)
                                && p.local_score >= 0,
                        ) * p.step_a.max(p.step_b) as usize,
                    acc.3 + usize::from(m == MatchType::Gap),
                    acc.4 + p.step_a.max(p.step_b) as usize,
                )
            });
        Stats {
            identical,
            mass_similar,
            similar,
            gaps,
            length,
        }
    }

    /// The mass(es) for the matched portion of the first sequence TODO: this assumes no terminal mods
    fn mass_a(&self) -> Multi<MolecularFormula> {
        if self.align_type().left.global_a() && self.align_type().right.global_a() {
            self.seq_a().bare_formulas()
        } else {
            let mut placed_a = vec![false; self.seq_a().ambiguous_modifications.len()];
            self.seq_a()[self.start_a()..self.start_a() + self.len_a()]
                .iter()
                .fold(Multi::default(), |acc, s| {
                    acc * s.formulas_greedy(&mut placed_a, &[], &[], &mut Vec::new())
                })
        }
    }

    /// The mass(es) for the matched portion of the second sequence
    fn mass_b(&self) -> Multi<MolecularFormula> {
        if self.align_type().left.global_b() && self.align_type().right.global_b() {
            self.seq_b().bare_formulas()
        } else {
            let mut placed_b = vec![false; self.seq_b().ambiguous_modifications.len()];
            self.seq_b()[self.start_b()..self.start_b() + self.len_b()]
                .iter()
                .fold(Multi::default(), |acc, s| {
                    acc * s.formulas_greedy(&mut placed_b, &[], &[], &mut Vec::new())
                })
        }
    }

    /// Get the mass delta for this match, if it is a (partial) local match it will only take the matched amino acids into account.
    /// If there are multiple possible masses for any of the stretches it returns the smallest difference.
    #[allow(clippy::missing_panics_doc)]
    fn mass_difference(&self) -> Mass {
        self.mass_a()
            .iter()
            .cartesian_product(self.mass_b().iter())
            .map(|(a, b)| a.monoisotopic_mass() - b.monoisotopic_mass())
            .min_by(|a, b| a.abs().value.total_cmp(&b.abs().value))
            .expect("An empty Multi<MolecularFormula>  was detected")
    }

    /// Get the error in ppm for this match, if it is a (partial) local match it will only take the matched amino acids into account.
    /// If there are multiple possible masses for any of the stretches it returns the smallest difference.
    fn ppm(&self) -> Ratio {
        self.mass_a()
            .iter()
            .cartesian_product(self.mass_b().iter())
            .map(|(a, b)| a.monoisotopic_mass().ppm(b.monoisotopic_mass()))
            .min_by(|a, b| a.value.total_cmp(&b.value))
            .expect("An empty Multi<MolecularFormula>  was detected")
    }

    /// Get a short representation of the alignment in CIGAR like format.
    /// It has one additional class `{a}(:{b})?(r|i)` denoting any special step with the given a and b step size, if b is not given it is the same as a.
    fn short(&self) -> String {
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
        let (_, _, output, last) = self.path().iter().fold(
            (self.start_a(), self.start_b(), String::new(), None),
            |(a, b, output, last), step| {
                let current_type = match (step.match_type, step.step_a, step.step_b) {
                    (MatchType::Isobaric, a, b) => StepType::Special(MatchType::Isobaric, a, b), // Catch any 1/1 isobaric sets before they are counted as Match/Mismatch
                    (_, 0, 1) => StepType::Insertion,
                    (_, 1, 0) => StepType::Deletion,
                    (_, 1, 1) if self.seq_a().sequence[a] == self.seq_b().sequence[b] => {
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
}

#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub(super) struct AlignmentInner {
    /// The scores of this alignment
    pub score: Score,
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

impl PartialOrd for AlignmentInner {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AlignmentInner {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.score.normalised.cmp(&other.score.normalised)
    }
}

/// Statistics for an alignment with some helper functions to easily retrieve the number of interest.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct Stats {
    /// The total number of identical positions
    pub identical: usize,
    /// The total number of mass similar positions, so including isobaric and rotations
    pub mass_similar: usize,
    /// The total number of similar positions, where the scoring matrix scores above 0, and not isobaric
    pub similar: usize,
    /// The total number of gap positions
    pub gaps: usize,
    /// The length of the alignment, the sum of the max of the step for A and B for each position.
    pub length: usize,
}

impl Stats {
    /// Get the identity as fraction.
    pub fn identity(&self) -> f64 {
        self.identical as f64 / self.length as f64
    }

    /// Get the mass similarity as fraction.
    pub fn mass_similarity(&self) -> f64 {
        self.mass_similar as f64 / self.length as f64
    }

    /// Get the similarity as fraction.
    pub fn similarity(&self) -> f64 {
        self.similar as f64 / self.length as f64
    }

    /// Get the gaps as fraction.
    pub fn gaps_fraction(&self) -> f64 {
        self.gaps as f64 / self.length as f64
    }
}

/// The score of an alignment
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct Score {
    /// The normalised score (absolute / max)
    pub normalised: OrderedFloat<f64>,
    /// The absolute score
    pub absolute: isize,
    /// The maximal possible score, the average score of the sequence slices on sequence a and b if they were aligned to themself, rounded down.
    ///    Think of it like this: `align(sequence_a.sequence[start_a..len_a], sequence_a.sequence[start_a..len_a])`.
    pub max: isize,
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use crate::{
        align::{align, matrix::BLOSUM62, AlignType, Alignment},
        peptide_complexity::Simple,
        system::da,
        AminoAcid, LinearPeptide, MultiChemical,
    };

    #[test]
    fn mass_difference() {
        // Test if the mass difference calculation is correct for some harder alignments.
        // A has an ambiguous AA, B and C have the two options, while D has a sub peptide of A.
        let a = LinearPeptide::pro_forma("AABAA", None)
            .unwrap()
            .simple()
            .unwrap();
        let b = LinearPeptide::pro_forma("AANAA", None)
            .unwrap()
            .simple()
            .unwrap();
        let c = LinearPeptide::pro_forma("AADAA", None)
            .unwrap()
            .simple()
            .unwrap();
        let d = LinearPeptide::pro_forma("ADA", None)
            .unwrap()
            .simple()
            .unwrap();

        assert!(
            align::<1, Simple, Simple>(
                &a,
                &b,
                BLOSUM62,
                crate::Tolerance::new_absolute(da(0.1)),
                AlignType::GLOBAL
            )
            .mass_difference()
            .value
            .abs()
                < f64::EPSILON
        );
        assert!(
            align::<1, Simple, Simple>(
                &a,
                &c,
                BLOSUM62,
                crate::Tolerance::new_absolute(da(0.1)),
                AlignType::GLOBAL
            )
            .mass_difference()
            .value
            .abs()
                < f64::EPSILON
        );
        assert!(
            align::<1, Simple, Simple>(
                &a,
                &d,
                BLOSUM62,
                crate::Tolerance::new_absolute(da(0.1)),
                AlignType::GLOBAL_B
            )
            .mass_difference()
            .value
            .abs()
                < f64::EPSILON
        );
        let mass_diff_nd = (AminoAcid::N.formulas()[0].monoisotopic_mass()
            - AminoAcid::D.formulas()[0].monoisotopic_mass())
        .value
        .abs();
        let mass_diff_bc = align::<1, Simple, Simple>(
            &b,
            &c,
            BLOSUM62,
            crate::Tolerance::new_absolute(da(0.1)),
            AlignType::GLOBAL_B,
        )
        .mass_difference()
        .value
        .abs();
        assert!(
            (mass_diff_bc - mass_diff_nd).abs() < 1E-10,
            "{mass_diff_bc} (peptides) should be equal to {mass_diff_nd} (ND)"
        );
    }
}
