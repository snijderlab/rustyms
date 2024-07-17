use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    system::{MassOverCharge, Ratio},
    AnnotatedSpectrum, Fragment, MassMode, Model, WithinTolerance,
};

impl AnnotatedSpectrum {
    /// Get a false discovery rate estimation for this annotation. See the [`Fdr`] struct for all
    /// statistics that can be retrieved. The returned tuple has the FDR for all peptides combined
    /// as first item and as second item a vector with for each peptide its individual scores.
    ///
    /// # Estimation
    /// It estimated the FDR by permutation. It takes all MZs for all fragments (or only the subset
    /// for each separate peptide) and counts how many match with the spectrum with offset
    /// `-25..=25 + π`. The returned result give insight in the average number of matches and
    /// standard deviation. The `π` offset is needed to guarantee non integer offsets preventing
    /// spurious matches from 1 Da isotopes.
    pub fn fdr(
        &self,
        fragments: &[Fragment],
        model: &Model,
        mass_mode: MassMode,
    ) -> (Fdr, Vec<Vec<Fdr>>) {
        let mzs = fragments
            .iter()
            .map(|f| (f.mz(mass_mode), f.peptidoform_index, f.peptide_index))
            .filter(|(mz, _, _)| model.mz_range.contains(mz))
            .collect_vec();

        let individual_peptides = self
            .peptide
            .peptidoforms()
            .iter()
            .enumerate()
            .map(|(peptidoform_index, peptidoform)| {
                peptidoform
                    .peptides()
                    .iter()
                    .enumerate()
                    .map(|(peptide_index, _)| {
                        self.internal_fdr(
                            mzs.iter()
                                .filter_map(|(mz, pi, ppi)| {
                                    (*pi == peptidoform_index && *ppi == peptide_index)
                                        .then_some(*mz)
                                })
                                .collect_vec()
                                .as_slice(),
                            model,
                        )
                    })
                    .collect()
            })
            .collect();
        (
            self.internal_fdr(
                mzs.iter().map(|(mz, _, _)| *mz).collect_vec().as_slice(),
                model,
            ),
            individual_peptides,
        )
    }

    fn internal_fdr(&self, mzs: &[MassOverCharge], model: &Model) -> Fdr {
        let mut results = Vec::with_capacity(51);

        for offset in -25..=25 {
            let peaks = self
                .spectrum
                .iter()
                .map(|p| {
                    p.experimental_mz
                        + MassOverCharge::new::<crate::system::mz>(
                            std::f64::consts::PI + f64::from(offset),
                        )
                })
                .collect_vec();
            let mut peak_annotated = vec![false; peaks.len()];
            let mut annotated = 0;
            for mass in mzs {
                // Get the index of the element closest to this value (spectrum is defined to always be sorted)
                let index = peaks
                    .binary_search_by(|p| p.value.total_cmp(&mass.value))
                    .unwrap_or_else(|i| i);

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
        let st_dev = (results.iter().map(|x| (x - average).powi(2)).sum::<f64>()
            / results.len() as f64)
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
