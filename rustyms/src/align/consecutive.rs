use crate::{
    align::*,
    imgt::*,
    peptide::{SimpleLinear, UnAmbiguous},
    system::Mass,
    *,
};
use std::collections::HashSet;

use itertools::Itertools;

/// Only available with if features `align` and `imgt` are turned on.
/// Align one sequence to multiple consecutive genes. Each gene can be controlled to be global to the left or free to allow unmatched residues between it and the previous gene.
/// If the sequence is too short to cover all genes only the genes that could be matched are returned.
/// # Panics
/// If there are not two or more genes listed. If the return number is 0.
#[cfg(feature = "imgt")]
#[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
pub fn consecutive_align<const STEPS: u16, A: Into<SimpleLinear> + Clone + Eq>(
    sequence: &LinearPeptide<A>,
    genes: &[(GeneType, AlignType)],
    species: Option<HashSet<Species, impl std::hash::BuildHasher + Clone + Send + Sync + Default>>,
    chains: Option<HashSet<ChainType, impl std::hash::BuildHasher + Clone + Send + Sync + Default>>,
    allele: AlleleSelection,
    tolerance: Tolerance<Mass>,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    return_number: usize,
) -> Vec<Vec<(Allele<'static>, Alignment<'static, UnAmbiguous, A>)>> {
    assert!(genes.len() >= 2);
    assert!(return_number != 0);

    let mut output: Vec<Vec<(Allele<'static>, Alignment<'static, UnAmbiguous, A>)>> =
        Vec::with_capacity(genes.len());

    let mut prev = 0;
    for gene in genes {
        let (left_sequence, use_species, use_chains) =
            output.last().and_then(|v| v.first()).map_or_else(
                || (sequence.clone(), species.clone(), chains.clone()),
                |last| {
                    prev += last.1.start_b() + last.1.len_b();
                    (
                        sequence.sub_peptide(prev..),
                        Some(std::iter::once(last.0.species).collect()),
                        Some(std::iter::once(last.0.gene.chain).collect()),
                    )
                },
            );

        if left_sequence.is_empty() {
            break;
        }

        output.push(
            Selection {
                species: use_species,
                chains: use_chains,
                allele,
                genes: Some([gene.0].into()),
            }
            .germlines()
            .map(|seq| {
                let alignment = align::<STEPS, UnAmbiguous, A>(
                    seq.sequence,
                    &left_sequence,
                    matrix,
                    tolerance,
                    gene.1,
                )
                .to_owned();
                (seq, alignment)
            })
            .k_largest_by(return_number, |a, b| a.1.cmp(&b.1))
            .collect_vec(),
        );
    }
    output
}

/// Only available with if features `align`, `rayon`, and `imgt` are turned on.
/// Align one sequence to multiple consecutive genes. Each gene can be controlled to be global to the left or free to allow unmatched residues between it and the previous gene.
/// If the sequence is too short to cover all genes only the genes that could be matched are returned.
/// # Panics
/// If there are not two or more genes listed. If the return number is 0.
#[cfg(all(feature = "rayon", feature = "imgt"))]
#[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
pub fn par_consecutive_align<const STEPS: u16, A: Into<SimpleLinear> + Clone + Eq + Send + Sync>(
    sequence: &LinearPeptide<A>,
    genes: &[(GeneType, AlignType)],
    species: Option<HashSet<Species, impl std::hash::BuildHasher + Clone + Send + Sync + Default>>,
    chains: Option<HashSet<ChainType, impl std::hash::BuildHasher + Clone + Send + Sync + Default>>,
    allele: AlleleSelection,
    tolerance: Tolerance<Mass>,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    return_number: usize,
) -> Vec<Vec<(Allele<'static>, Alignment<'static, UnAmbiguous, A>)>> {
    use rayon::iter::ParallelIterator;

    assert!(genes.len() >= 2);
    assert!(return_number != 0);

    let mut output: Vec<Vec<(Allele<'static>, Alignment<'static, UnAmbiguous, A>)>> =
        Vec::with_capacity(genes.len());

    let mut prev = 0;
    for gene in genes {
        let (left_sequence, use_species, use_chains) =
            output.last().and_then(|v| v.first()).map_or_else(
                || (sequence.clone(), species.clone(), chains.clone()),
                |last| {
                    prev += last.1.start_b() + last.1.len_b();
                    (
                        sequence.sub_peptide(prev..),
                        Some(std::iter::once(last.0.species).collect()),
                        Some(std::iter::once(last.0.gene.chain).collect()),
                    )
                },
            );

        if left_sequence.is_empty() {
            break;
        }

        output.push(
            Selection {
                species: use_species,
                chains: use_chains,
                allele,
                genes: Some([gene.0].into()),
            }
            .par_germlines()
            .map(|seq| {
                let alignment = align::<STEPS, UnAmbiguous, A>(
                    seq.sequence,
                    &left_sequence,
                    matrix,
                    tolerance,
                    gene.1,
                );
                (seq, alignment.to_owned())
            })
            .collect::<Vec<_>>()
            .into_iter()
            .k_largest_by(return_number, |a, b| a.1.cmp(&b.1))
            .collect_vec(),
        );
    }
    output
}
