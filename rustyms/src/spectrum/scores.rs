//! Scoring of annotated spectra

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    fragment::{Fragment, FragmentKind},
    peptide::ExtremelySimple,
    system::{
        f64::{MassOverCharge, Ratio},
        mass_over_charge::mz,
    },
    AnnotatedSpectrum, LinearPeptide, MassMode, Model, WithinTolerance,
};

impl AnnotatedSpectrum {
    /// Get the spectrum scores for this annotated spectrum.
    /// The returned tuple has the scores for all peptides combined as first item
    /// and as second item a vector with for each peptide its individual scores.
    pub fn scores(&self, fragments: &[Fragment]) -> (Scores, Vec<Scores>) {
        let mut individual_peptides = Vec::new();
        let total_intensity: f64 = self.spectrum.iter().map(|p| *p.intensity).sum();
        for (peptidoform_index, peptidoform) in self.peptide.peptidoforms().iter().enumerate() {
            for (peptide_index, peptide) in peptidoform.peptides().iter().enumerate() {
                let (recovered_fragments, peaks, intensity_annotated) =
                    self.base_score(fragments, Some(peptide_index), None);
                let positions = self.score_positions(peptide_index, None);
                individual_peptides.push(Scores {
                    score: Score::Position {
                        fragments: recovered_fragments,
                        peaks,
                        intensity: Recovered::new(intensity_annotated, total_intensity),
                        positions: Recovered::new(positions, peptide.len() as u32),
                    },
                    ions: self.score_individual_ions(
                        fragments,
                        Some((peptide_index, peptide)),
                        total_intensity,
                    ),
                });
            }
        }
        // Get the statistics for the combined peptides
        let (recovered_fragments, peaks, intensity_annotated) =
            self.base_score(fragments, None, None);
        let unique_formulas = self.score_unique_formulas(fragments, None, None);
        (
            Scores {
                score: Score::UniqueFormulas {
                    fragments: recovered_fragments,
                    peaks,
                    intensity: Recovered::new(intensity_annotated, total_intensity),
                    unique_formulas,
                },
                ions: self.score_individual_ions::<ExtremelySimple>(
                    fragments,
                    None,
                    total_intensity,
                ),
            },
            individual_peptides,
        )
    }

    /// Get the base score of this spectrum
    fn base_score(
        &self,
        fragments: &[Fragment],
        peptide_index: Option<usize>,
        ion: Option<FragmentKind>,
    ) -> (Recovered<u32>, Recovered<u32>, f64) {
        let (peaks_annotated, fragments_found, intensity_annotated) = self
            .spectrum
            .iter()
            .filter_map(|p| {
                let number = p
                    .annotation
                    .iter()
                    .filter(|(a, _)| {
                        peptide_index.map_or(true, |i| a.peptide_index == i)
                            && ion.map_or(true, |kind| a.ion.kind() == kind)
                    })
                    .count() as u32;
                if number == 0 {
                    None
                } else {
                    Some((number, *p.intensity))
                }
            })
            .fold((0u32, 0u32, 0.0), |(n, f, intensity), p| {
                (n + 1, f + p.0, intensity + p.1)
            });
        let total_fragments = fragments
            .iter()
            .filter(|f| {
                peptide_index.map_or(true, |i| f.peptide_index == i)
                    && ion.map_or(true, |kind| f.ion.kind() == kind)
            })
            .count() as u32;
        (
            Recovered::new(fragments_found, total_fragments),
            Recovered::new(peaks_annotated, self.spectrum.len() as u32),
            intensity_annotated,
        )
    }

    /// Get the total number of positions covered
    fn score_positions(&self, peptide_index: usize, ion: Option<FragmentKind>) -> u32 {
        self.spectrum
            .iter()
            .flat_map(|p| {
                p.annotation
                    .iter()
                    .filter(|(a, _)| {
                        a.peptide_index == peptide_index
                            && ion.map_or(true, |kind| a.ion.kind() == kind)
                    })
                    .filter_map(|(a, _)| a.ion.position())
            })
            .map(|pos| pos.sequence_index)
            .unique()
            .count() as u32
    }

    /// Get the amount of unique formulas recovered
    fn score_unique_formulas(
        &self,
        fragments: &[Fragment],
        peptide_index: Option<usize>,
        ion: Option<FragmentKind>,
    ) -> Recovered<u32> {
        let num_annotated = self
            .spectrum
            .iter()
            .flat_map(|p| {
                p.annotation.iter().filter(|(a, _)| {
                    peptide_index.map_or(true, |i| a.peptide_index == i)
                        && ion.map_or(true, |kind| a.ion.kind() == kind)
                })
            })
            .map(|(f, _)| f.formula.clone())
            .unique()
            .count() as u32;
        let total_fragments = fragments
            .iter()
            .filter(|f| {
                peptide_index.map_or(true, |i| f.peptide_index == i)
                    && ion.map_or(true, |kind| f.ion.kind() == kind)
            })
            .map(|f| f.formula.clone())
            .unique()
            .count() as u32;
        Recovered::new(num_annotated, total_fragments)
    }

    /// Get the scores for the individual ion series
    fn score_individual_ions<T>(
        &self,
        fragments: &[Fragment],
        peptide: Option<(usize, &LinearPeptide<T>)>,
        total_intensity: f64,
    ) -> Vec<(FragmentKind, Score)> {
        [
            FragmentKind::a,
            FragmentKind::b,
            FragmentKind::c,
            FragmentKind::d,
            FragmentKind::v,
            FragmentKind::w,
            FragmentKind::x,
            FragmentKind::y,
            FragmentKind::z,
        ]
        .iter()
        .copied()
        .filter_map(|ion| {
            let (recovered_fragments, peaks, intensity_annotated) =
                self.base_score(fragments, peptide.as_ref().map(|p| p.0), Some(ion));
            if let Some((index, peptide)) = peptide {
                if recovered_fragments.total > 0 {
                    let positions = self.score_positions(index, Some(ion));
                    Some((
                        ion,
                        Score::Position {
                            fragments: recovered_fragments,
                            peaks,
                            intensity: Recovered::new(intensity_annotated, total_intensity),
                            positions: Recovered::new(positions, peptide.len() as u32),
                        },
                    ))
                } else {
                    None
                }
            } else if recovered_fragments.total > 0 {
                let unique_formulas = self.score_unique_formulas(fragments, None, Some(ion));
                Some((
                    ion,
                    Score::UniqueFormulas {
                        fragments: recovered_fragments,
                        peaks,
                        intensity: Recovered::new(intensity_annotated, total_intensity),
                        unique_formulas,
                    },
                ))
            } else {
                None
            }
        })
        .chain(
            [
                FragmentKind::Y,
                FragmentKind::Oxonium,
                FragmentKind::immonium,
                FragmentKind::m,
                FragmentKind::diagnostic,
                FragmentKind::precursor,
            ]
            .iter()
            .copied()
            .filter_map(|ion| {
                let (recovered_fragments, peaks, intensity_annotated) =
                    self.base_score(fragments, peptide.as_ref().map(|p| p.0), Some(ion));
                if recovered_fragments.total > 0 {
                    let unique_formulas = self.score_unique_formulas(
                        fragments,
                        peptide.as_ref().map(|p| p.0),
                        Some(ion),
                    );
                    Some((
                        ion,
                        Score::UniqueFormulas {
                            fragments: recovered_fragments,
                            peaks,
                            intensity: Recovered::new(intensity_annotated, total_intensity),
                            unique_formulas,
                        },
                    ))
                } else {
                    None
                }
            }),
        )
        .collect()
    }

    /// Get a false discovery rate estimation for this annotation. See the [`Fdr`] struct for all statistics that can be retrieved.
    pub fn fdr(&self, fragments: &[Fragment], model: &Model) -> Fdr {
        let masses = fragments
            .iter()
            .map(|f| f.mz(MassMode::Monoisotopic))
            .collect_vec();
        let mut results = Vec::with_capacity(50);

        for offset in -25..=25 {
            let peaks = self
                .spectrum
                .iter()
                .map(|p| {
                    p.experimental_mz
                        + MassOverCharge::new::<mz>(std::f64::consts::PI + f64::from(offset))
                })
                .collect_vec();
            let mut peak_annotated = vec![false; peaks.len()];
            let mut annotated = 0;
            for mass in &masses {
                // Get the index of the element closest to this value (spectrum is defined to always be sorted)
                let index = peaks
                    .binary_search_by(|p| p.value.total_cmp(&mass.value))
                    .map_or_else(|i| i, |i| i);

                // Check index-1, index and index+1 (if existing) to find the one with the lowest ppm
                let mut closest = (0, Ratio::new::<crate::system::ratio::ppm>(f64::INFINITY));
                #[allow(clippy::needless_range_loop)] // I like this better
                for i in if index == 0 { 0 } else { index - 1 }
                    ..=(index + 1).min(self.spectrum.len() - 1)
                {
                    let ppm = peaks[i].ppm(*mass);
                    if ppm < closest.1 {
                        closest = (i, ppm);
                    }
                }

                if model
                    .tolerance
                    .within(&self.spectrum[closest.0].experimental_mz, mass)
                    && !peak_annotated[closest.0]
                {
                    annotated += 1;
                    peak_annotated[closest.0] = true;
                }
            }
            results.push(f64::from(annotated) / self.spectrum.len() as f64);
        }
        let average = results.iter().sum::<f64>() / results.len() as f64;
        let st_dev = results
            .iter()
            .map(|x| (x - average).powi(2))
            .sum::<f64>()
            .sqrt();
        let actual = self
            .spectrum
            .iter()
            .filter(|p| !p.annotation.is_empty())
            .count() as f64
            / self.spectrum.len() as f64;

        Fdr {
            actual,
            average_false: average,
            standard_deviation_false: st_dev,
        }
    }
}

/// The scores for an annotated spectrum
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
#[non_exhaustive]
pub struct Scores {
    /// The scores, based on unique formulas for all peptides combined or based on positions for single peptides.
    pub score: Score,
    /// The scores per [`FragmentKind`], based on unique formulas for all peptides combined or any fragment kind that is not an ion series, or based on positions in the other case.
    pub ions: Vec<(FragmentKind, Score)>,
}

/// The scores for a single fragment series for a single peptide in an annotated spectrum
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub enum Score {
    /// A score for a something that has peptide position coverage
    Position {
        /// The fraction of the total fragments that could be annotated
        fragments: Recovered<u32>,
        /// The fraction of the total peaks that could be annotated
        peaks: Recovered<u32>,
        /// The fraction of the total intensity that could be annotated
        intensity: Recovered<f64>,
        /// The fraction of the total positions that has at least one fragment found
        positions: Recovered<u32>,
    },
    /// A score for something that does not have position coverage, but instead is scored on the number of unique formulas
    UniqueFormulas {
        /// The fraction of the total fragments that could be annotated
        fragments: Recovered<u32>,
        /// The fraction of the total peaks that could be annotated
        peaks: Recovered<u32>,
        /// The fraction of the total intensity that could be annotated
        intensity: Recovered<f64>,
        /// The fraction of with unique formulas that has been found
        unique_formulas: Recovered<u32>,
    },
}
/// A single statistic that has a total number and a subset of that found
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
#[non_exhaustive]
pub struct Recovered<T> {
    /// The number actually found
    pub found: T,
    /// The total number
    pub total: T,
}

impl<T> Recovered<T> {
    /// Create a new recovered statistic
    fn new(found: impl Into<T>, total: impl Into<T>) -> Self {
        Self {
            found: found.into(),
            total: total.into(),
        }
    }
}

impl<T> Recovered<T>
where
    f64: From<T>,
    T: Copy,
{
    /// Get the recovered amount as fraction
    pub fn fraction(&self) -> f64 {
        f64::from(self.found) / f64::from(self.total)
    }
}

/// A false discovery rate for an annotation to a spectrum
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Fdr {
    /// The fraction of the total (assumed to be true) peaks that could be annotated
    pub actual: f64,
    /// The average fraction of the false peaks that could be annotated
    pub average_false: f64,
    /// The standard deviation of the false peaks that could be annotated
    pub standard_deviation_false: f64,
}

impl Fdr {
    /// Get the false discovery rate (as a fraction).
    /// The average number of false peaks annotated divided by the average number of annotated peaks.
    pub fn fdr(&self) -> f64 {
        self.average_false / self.actual
    }

    /// Get the number of standard deviations the number of annotated peaks is from the average number of false annotations.
    pub fn sigma(&self) -> f64 {
        (self.actual - self.average_false) / self.standard_deviation_false
    }

    /// Get the score of this annotation. Defined as the log2 of the sigma.
    pub fn score(&self) -> f64 {
        self.sigma().log2()
    }
}
