use itertools::Itertools;

use crate::{AminoAcid, SequenceElement};

/// A protease defined by it ability to cut at any site identified by the right amino acids at the n and c terminal.
/// Each position is identified by an option, a none means that there is no specificity at this position. If there is
/// a specificity at a certain position any amino acid that is contained in the set is allowed (see [`AminoAcid::canonical_identical`]).
pub struct Protease {
    /// The amino acids n terminal of the cut site.
    pub n_term: Vec<Option<Vec<AminoAcid>>>,
    /// The amino acids c terminal of the cut site.
    pub c_term: Vec<Option<Vec<AminoAcid>>>,
}

impl Protease {
    /// Define a simple protease that cuts exactly between the specified sequences.
    pub fn new(n_term: &[AminoAcid], c_term: &[AminoAcid]) -> Self {
        Self {
            n_term: n_term.iter().map(|aa| Some(vec![*aa])).collect_vec(),
            c_term: c_term.iter().map(|aa| Some(vec![*aa])).collect_vec(),
        }
    }

    /// Define a protease that cuts on the n terminal side of the provided amino acids.
    pub fn n_terminal_of(residues: &[AminoAcid]) -> Self {
        Self {
            n_term: vec![Some(residues.to_vec())],
            c_term: Vec::new(),
        }
    }

    /// Define a protease that cuts on the c terminal side of the provided amino acids.
    pub fn c_terminal_of(residues: &[AminoAcid]) -> Self {
        Self {
            c_term: vec![Some(residues.to_vec())],
            n_term: Vec::new(),
        }
    }

    /// All locations in the given sequence where this protease could cut
    pub fn match_locations<T>(&self, sequence: &[SequenceElement<T>]) -> Vec<usize> {
        (self.n_term.len()..sequence.len() - self.c_term.len())
            .filter(|i| self.matches_at(&sequence[i - self.n_term.len()..i + self.c_term.len()]))
            .collect_vec()
    }

    fn matches_at<T>(&self, slice: &[SequenceElement<T>]) -> bool {
        debug_assert!(slice.len() == self.n_term.len() + self.c_term.len());
        'positions: for (actual, pattern) in slice
            .iter()
            .zip(self.n_term.iter().chain(self.c_term.iter()))
        {
            if let Some(pattern) = pattern {
                for option in pattern {
                    if option.canonical_identical(actual.aminoacid.aminoacid()) {
                        continue 'positions;
                    }
                }
                return false;
            }
        }
        true
    }
}
