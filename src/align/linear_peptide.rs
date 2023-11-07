use std::{borrow::Cow, ops::Range};

use crate::{system::Mass, LinearPeptide, MolecularFormula, SequenceElement};

use super::{
    multiple_sequence_alignment::{MSAPlacement, MSAPosition},
    MassAlignable,
};

impl MassAlignable for LinearPeptide {
    // TODO: any location can have multiple masses, so ambiguous modifications and ambiguous amino acids can be taken into account right?
    fn index(&self, index: usize, sequence_index: usize) -> &SequenceElement {
        debug_assert!(sequence_index == 0);
        &self.sequence[index]
    }
    fn index_slice(
        &self,
        index: Range<usize>,
        sequence_index: usize,
    ) -> Cow<'_, [SequenceElement]> {
        debug_assert!(sequence_index == 0);
        Cow::from(&self.sequence[index])
    }
    fn total_length(&self) -> usize {
        self.sequence.len()
    }
    fn number_of_sequences(&self) -> usize {
        1
    }
    fn sequence_bounds(&self) -> Vec<(usize, usize)> {
        vec![(0, self.sequence.len())]
    }
    fn calculate_masses<const STEPS: usize>(&self) -> Vec<Vec<Vec<Option<Mass>>>> {
        vec![(0..=self.sequence.len())
            .map(|index| {
                (1..=STEPS)
                    .map(|size| {
                        if index < size {
                            None
                        } else {
                            Some(
                                self.sequence[index - size..index]
                                    .iter()
                                    .map(|s| s.formula_all().unwrap())
                                    .sum::<MolecularFormula>()
                                    .monoisotopic_mass()
                                    .unwrap(),
                            )
                        }
                    })
                    .collect::<Vec<Option<Mass>>>()
            })
            .collect::<Vec<_>>()]
    }
    fn sequences_with_path(
        &self,
        is_a: bool,
        start: usize,
        path: &[super::Piece],
    ) -> Vec<super::multiple_sequence_alignment::MSAPlacement> {
        vec![MSAPlacement {
            start,
            sequence: self.clone(),
            path: path
                .iter()
                .map(|piece| {
                    if is_a {
                        if piece.step_a == 0 {
                            MSAPosition::Gap
                        } else {
                            MSAPosition::Placed(piece.step_a as usize, piece.step_b as usize)
                        }
                    } else if piece.step_b == 0 {
                        MSAPosition::Gap
                    } else {
                        MSAPosition::Placed(piece.step_b as usize, piece.step_a as usize)
                    }
                })
                .collect(),
            score: 0,
            normalised_score: 0.0,
        }]
    }
}
