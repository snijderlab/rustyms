//! Spectrum related code

use std::{cmp::Ordering, iter::FusedIterator};

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use uom::num_traits::Zero;

use crate::{
    fragment::Fragment,
    system::{f64::*, mass_over_charge::mz},
    ComplexPeptide, Model,
};

/// The mode of mass to use
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug, Serialize, Deserialize,
)]
pub enum MassMode {
    /// Monoisotopic mass, use the base isotope to calculate the mass (eg always 12C)
    #[default]
    Monoisotopic,
    /// The average weight, the average between all occurring isotopes (eg something in between 12C and 13C depending on the number of C)
    Average,
    /// The most abundant mass, the most abundant single isotopic species (eg 12C or 13C depending on the number of C)
    MostAbundant,
}

// TODO: Trace Trait to generate the correct time points
// Add optional traces to raw and annotated, plus display nicely in annotator
// Future: add centroiding to build a raw from a trace

/// A trace, generic over the second dimension (eg time (ms1) or mz (ms2))
pub struct Trace<T> {
    data: Vec<f64>,
    step: T,
}

impl<T> Trace<T> {
    /// Create a new trace
    pub fn new(data: &[f64], step: T) -> Self {
        Self {
            data: data.to_owned(),
            step,
        }
    }
}

impl<T> Trace<T>
where
    T: std::ops::Mul<usize, Output = T> + Copy,
{
    /// Get the data of this trace, alongside the value in the second dimension
    pub fn data(&self) -> impl Iterator<Item = (T, f64)> + '_ {
        self.data
            .iter()
            .enumerate()
            .map(|(i, v)| (self.step * (i + 1), *v)) // TODO: Does it start at 1 or at 0?
    }
}

/// The trait for all spectra that contain peaks.
pub trait PeakSpectrum:
    Extend<Self::PeakType>
    + IntoIterator<Item = Self::PeakType>
    + std::ops::Index<usize, Output = Self::PeakType>
{
    /// The type of peaks this spectrum contains
    type PeakType;
    /// The type of spectrum iterator this spectrum generates
    type Iter<'a>: DoubleEndedIterator + ExactSizeIterator + FusedIterator
    where
        Self: 'a;
    /// Return the slice of peaks that is within the given tolerance bounds.
    fn binary_search(&self, low: MassOverCharge, high: MassOverCharge) -> &[Self::PeakType];
    /// Get the full spectrum
    fn spectrum(&self) -> Self::Iter<'_>;
    /// Add a single peak
    fn add_peak(&mut self, item: Self::PeakType);
}

/// A raw spectrum (meaning not annotated yet)
#[derive(Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
pub struct RawSpectrum {
    /// The title (as used in MGF)
    pub title: String,
    /// The number of scans
    pub num_scans: u64,
    /// The retention time
    pub rt: Time,
    /// The found precursor charge
    pub charge: Charge,
    /// The found precursor mass
    pub mass: Mass,
    /// The found precursor intensity
    pub intensity: Option<f64>,
    /// The peaks of which this spectrum consists
    spectrum: Vec<RawPeak>,
    /// MGF: if present the SEQUENCE line
    pub sequence: Option<String>,
    /// MGF TITLE: if present the raw file where this mgf was made from
    pub raw_file: Option<String>,
    /// MGF TITLE: if present the raw file scan number
    pub raw_scan_number: Option<usize>,
    /// MGF TITLE: index number
    pub raw_index: Option<usize>,
    /// MGF TITLE: sample number
    pub sample: Option<usize>,
    /// MGF TITLE: period number
    pub period: Option<usize>,
    /// MGF TITLE: cycle number
    pub cycle: Option<usize>,
    /// MGF TITLE: experiment number
    pub experiment: Option<usize>,
    /// MGF TITLE: controllerType number
    pub controller_type: Option<usize>,
    /// MGF TITLE: controllerNumber number
    pub controller_number: Option<usize>,
}

impl RawSpectrum {
    /// Filter the spectrum to retain all with an intensity above `filter_threshold` times the maximal intensity.
    ///
    /// # Panics
    /// It panics if any peaks has an intensity that is NaN.
    pub fn noise_filter(&mut self, filter_threshold: f64) {
        let max = self
            .spectrum
            .iter()
            .map(|p| *p.intensity)
            .reduce(f64::max)
            .unwrap_or(f64::INFINITY);
        self.spectrum
            .retain(|p| *p.intensity >= max * filter_threshold);
        self.spectrum.shrink_to_fit();
    }

    /// Annotate this spectrum with the given peptide and given fragments see [`crate::ComplexPeptide::generate_theoretical_fragments`].
    ///
    /// # Panics
    /// If any fragment does not have a defined m/z
    pub fn annotate(
        &self,
        peptide: ComplexPeptide,
        theoretical_fragments: &[Fragment],
        model: &Model,
        mode: MassMode,
    ) -> AnnotatedSpectrum {
        let mut annotated = AnnotatedSpectrum {
            title: self.title.clone(),
            num_scans: self.num_scans,
            rt: self.rt,
            charge: self.charge,
            mass: self.mass,
            peptide,
            spectrum: self
                .spectrum
                .iter()
                .map(AnnotatedPeak::background)
                .collect(),
        };

        for fragment in theoretical_fragments {
            // Get the index of the element closest to this value (spectrum is defined to always be sorted)
            let index = self
                .spectrum
                .binary_search_by(|p| p.mz.value.total_cmp(&fragment.mz(mode).value))
                .map_or_else(|i| i, |i| i);

            // Check index-1, index and index+1 (if existing) to find the one with the lowest ppm
            let mut closest = (0, f64::INFINITY);
            for i in
                if index == 0 { 0 } else { index - 1 }..=(index + 1).min(self.spectrum.len() - 1)
            {
                let ppm = self.spectrum[i].ppm(fragment, mode).value;
                if ppm < closest.1 {
                    closest = (i, ppm);
                }
            }

            if closest.1 < model.ppm.value {
                annotated.spectrum[closest.0]
                    .annotation
                    .push(fragment.clone());
            }
        }

        annotated
    }
}

impl Extend<RawPeak> for RawSpectrum {
    fn extend<T: IntoIterator<Item = RawPeak>>(&mut self, iter: T) {
        self.spectrum.extend(iter);
        self.spectrum.sort_unstable();
    }
}

impl IntoIterator for RawSpectrum {
    type Item = RawPeak;
    type IntoIter = std::vec::IntoIter<RawPeak>;
    fn into_iter(self) -> Self::IntoIter {
        self.spectrum.into_iter()
    }
}

impl std::ops::Index<usize> for RawSpectrum {
    type Output = RawPeak;
    fn index(&self, index: usize) -> &Self::Output {
        &self.spectrum[index]
    }
}

impl PeakSpectrum for RawSpectrum {
    type PeakType = RawPeak;
    type Iter<'a> = std::slice::Iter<'a, Self::PeakType>;

    /// Return the slice of peaks that is within the given tolerance bounds.
    fn binary_search(&self, low: MassOverCharge, high: MassOverCharge) -> &[RawPeak] {
        let left_idx = match self
            .spectrum
            .binary_search_by(|a| a.mz.value.total_cmp(&low.value))
        {
            Result::Ok(idx) | Result::Err(idx) => {
                let mut idx = idx.saturating_sub(1);
                while idx > 0 && self.spectrum[idx].mz.value.total_cmp(&low.value) != Ordering::Less
                {
                    idx -= 1;
                }
                idx
            }
        };

        let right_idx = match self.spectrum[left_idx..]
            .binary_search_by(|a| a.mz.value.total_cmp(&high.value))
        {
            Result::Ok(idx) | Err(idx) => {
                let mut idx = idx + left_idx;
                while idx < self.spectrum.len()
                    && self.spectrum[idx].mz.value.total_cmp(&high.value) != Ordering::Greater
                {
                    idx = idx.saturating_add(1);
                }
                idx.min(self.spectrum.len())
            }
        };
        &self.spectrum[left_idx..right_idx]
    }

    fn spectrum(&self) -> Self::Iter<'_> {
        self.spectrum.iter()
    }

    fn add_peak(&mut self, item: Self::PeakType) {
        let index = self.spectrum.binary_search(&item).map_or_else(|i| i, |i| i);
        self.spectrum.insert(index, item);
    }
}

impl Default for RawSpectrum {
    fn default() -> Self {
        Self {
            title: String::new(),
            num_scans: 0,
            rt: Time::zero(),
            charge: Charge::new::<e>(1.0),
            mass: Mass::zero(),
            spectrum: Vec::new(),
            intensity: None,
            sequence: None,
            raw_file: None,
            raw_scan_number: None,
            raw_index: None,
            sample: None,
            period: None,
            cycle: None,
            experiment: None,
            controller_type: None,
            controller_number: None,
        }
    }
}

/// An annotated spectrum
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct AnnotatedSpectrum {
    /// The title (as used in MGF)
    pub title: String,
    /// The number of scans
    pub num_scans: u64,
    /// The retention time
    pub rt: Time,
    /// The found precursor charge
    pub charge: Charge,
    /// The found precursor mass
    pub mass: Mass,
    /// The peptide with which this spectrum was annotated
    pub peptide: ComplexPeptide,
    /// The spectrum
    spectrum: Vec<AnnotatedPeak>,
}

impl Extend<AnnotatedPeak> for AnnotatedSpectrum {
    fn extend<T: IntoIterator<Item = AnnotatedPeak>>(&mut self, iter: T) {
        self.spectrum.extend(iter);
        self.spectrum.sort_unstable();
    }
}

impl IntoIterator for AnnotatedSpectrum {
    type Item = AnnotatedPeak;
    type IntoIter = std::vec::IntoIter<AnnotatedPeak>;
    fn into_iter(self) -> Self::IntoIter {
        self.spectrum.into_iter()
    }
}

impl std::ops::Index<usize> for AnnotatedSpectrum {
    type Output = AnnotatedPeak;
    fn index(&self, index: usize) -> &Self::Output {
        &self.spectrum[index]
    }
}

impl PeakSpectrum for AnnotatedSpectrum {
    type PeakType = AnnotatedPeak;
    type Iter<'a> = std::slice::Iter<'a, Self::PeakType>;

    /// Return the slice of peaks that have experimental mz values within the given tolerance bounds.
    fn binary_search(&self, low: MassOverCharge, high: MassOverCharge) -> &[AnnotatedPeak] {
        let left_idx = match self
            .spectrum
            .binary_search_by(|a| a.experimental_mz.value.total_cmp(&low.value))
        {
            Result::Ok(idx) | Result::Err(idx) => {
                let mut idx = idx.saturating_sub(1);
                while idx > 0
                    && self.spectrum[idx]
                        .experimental_mz
                        .value
                        .total_cmp(&low.value)
                        != Ordering::Less
                {
                    idx -= 1;
                }
                idx
            }
        };

        let right_idx = match self.spectrum[left_idx..]
            .binary_search_by(|a| a.experimental_mz.value.total_cmp(&high.value))
        {
            Result::Ok(idx) | Err(idx) => {
                let mut idx = idx + left_idx;
                while idx < self.spectrum.len()
                    && self.spectrum[idx]
                        .experimental_mz
                        .value
                        .total_cmp(&high.value)
                        != Ordering::Greater
                {
                    idx = idx.saturating_add(1);
                }
                idx.min(self.spectrum.len())
            }
        };
        &self.spectrum[left_idx..right_idx]
    }

    fn spectrum(&self) -> Self::Iter<'_> {
        self.spectrum.iter()
    }

    fn add_peak(&mut self, item: Self::PeakType) {
        let index = self.spectrum.binary_search(&item).map_or_else(|i| i, |i| i);
        self.spectrum.insert(index, item);
    }
}

impl AnnotatedSpectrum {
    /// Get the spectrum scores for this annotated spectrum
    pub fn scores(&self, fragments: &[Fragment]) -> Vec<Scores> {
        let mut results = Vec::new();
        let total_intensity: f64 = self.spectrum.iter().map(|p| *p.intensity).sum();
        for (peptide_index, peptide) in self.peptide.peptides().iter().enumerate() {
            let (num_annotated, intensity_annotated) = self
                .spectrum
                .iter()
                .filter(|p| {
                    p.annotation
                        .iter()
                        .any(|a| a.peptide_index == peptide_index)
                })
                .fold((0, 0.0), |(n, intensity), p| {
                    (n + 1, intensity + *p.intensity)
                });
            let total_fragments = fragments
                .iter()
                .filter(|f| f.peptide_index == peptide_index)
                .count();
            let fragments_found = f64::from(num_annotated) / total_fragments as f64;
            let peaks_annotated = f64::from(num_annotated) / self.spectrum.len() as f64;
            let intensity_annotated = intensity_annotated / total_intensity;
            let positions_covered = self
                .spectrum
                .iter()
                .flat_map(|p| {
                    p.annotation
                        .iter()
                        .filter(|a| a.peptide_index == peptide_index)
                        .filter_map(|a| a.ion.position())
                })
                .map(|pos| pos.sequence_index)
                .unique()
                .count() as f64
                / peptide.len() as f64;
            results.push(Scores {
                peptide_index,
                fragments_found,
                peaks_annotated,
                intensity_annotated,
                positions_covered,
            });
        }
        results
    }

    /// Get a false discovery rate estimation for this annotation. See the [`FDR`] struct for all statistics that can be retrieved.
    pub fn fdr(&self, fragments: &[Fragment], model: &Model) -> FDR {
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
                let mut closest = (0, f64::INFINITY);
                #[allow(clippy::needless_range_loop)] // I like this better
                for i in if index == 0 { 0 } else { index - 1 }
                    ..=(index + 1).min(self.spectrum.len() - 1)
                {
                    let ppm = peaks[i].ppm(*mass);
                    if ppm < closest.1 {
                        closest = (i, ppm);
                    }
                }

                if closest.1 < model.ppm.value && !peak_annotated[closest.0] {
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

        FDR {
            actual_peaks_annotated: actual,
            average_false_peaks_annotated: average,
            standard_deviation_false_peaks_annotated: st_dev,
        }
    }
}

/// The scores for a single peptide in an annotated spectrum
pub struct Scores {
    /// The peptide index
    pub peptide_index: usize,
    /// The fraction of the total fragments that could be annotated
    pub fragments_found: f64,
    /// The fraction of the total peaks that could be annotated
    pub peaks_annotated: f64,
    /// The fraction of the total intensity that could be annotated
    pub intensity_annotated: f64,
    /// The fraction of the total positions that has at least one fragment found
    pub positions_covered: f64,
}

/// A false discovery rate for an annotation to a spectrum
pub struct FDR {
    /// The fraction of the total peaks that could be annotated
    pub actual_peaks_annotated: f64,
    /// The average fraction of the false peaks that could be annotated
    pub average_false_peaks_annotated: f64,
    /// The standard deviation of the false peaks that could be annotated
    pub standard_deviation_false_peaks_annotated: f64,
}

impl FDR {
    /// Get the false discovery rate (as a fraction).
    /// The average number of false peaks annotated divided by the average number of annotated peaks.
    pub fn fdr(&self) -> f64 {
        self.average_false_peaks_annotated / self.actual_peaks_annotated
    }

    /// Get the number of standard deviations the number of annotated peaks is from the average number of false annotations.
    pub fn sigma(&self) -> f64 {
        (self.actual_peaks_annotated - self.average_false_peaks_annotated)
            / self.standard_deviation_false_peaks_annotated
    }

    /// Get the score of this annotation. Defined as the log2 of the sigma.
    pub fn score(&self) -> f64 {
        self.sigma().log2()
    }
}

/// A raw peak
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RawPeak {
    /// The mz value of this peak
    pub mz: MassOverCharge,
    /// The intensity of this peak
    pub intensity: OrderedFloat<f64>,
    /// The charge of this peak
    pub charge: Charge, // TODO: Is this item needed? (mgf has it, not used in rustyms)
}

impl PartialOrd for RawPeak {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RawPeak {
    /// Use `f64::total_cmp` on `self.mz`
    fn cmp(&self, other: &Self) -> Ordering {
        self.mz.value.total_cmp(&other.mz.value)
    }
}

impl PartialEq for RawPeak {
    /// Use `f64::total_cmp` on all fields to detect total equality
    fn eq(&self, other: &Self) -> bool {
        self.mz.value.total_cmp(&other.mz.value) == Ordering::Equal
            && self.intensity.total_cmp(&other.intensity) == Ordering::Equal
            && self.charge.value.total_cmp(&other.charge.value) == Ordering::Equal
    }
}

impl Eq for RawPeak {}

impl RawPeak {
    /// Determine the ppm error for the given fragment
    pub fn ppm(&self, fragment: &Fragment, mode: MassMode) -> MassOverCharge {
        MassOverCharge::new::<mz>(self.mz.ppm(fragment.mz(mode)))
    }
}

/// An annotated peak
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AnnotatedPeak {
    /// The experimental mz
    pub experimental_mz: MassOverCharge,
    /// The experimental intensity
    pub intensity: OrderedFloat<f64>,
    /// The charge
    pub charge: Charge, // TODO: Is this item needed? (mgf has it, not used in rustyms)
    /// The annotation, if present
    pub annotation: Vec<Fragment>,
}

impl AnnotatedPeak {
    /// Make a new annotated peak with the given annotation
    pub fn new(peak: &RawPeak, annotation: Fragment) -> Self {
        Self {
            experimental_mz: peak.mz,
            intensity: peak.intensity,
            charge: peak.charge,
            annotation: vec![annotation],
        }
    }

    /// Make a new annotated peak if no annotation is possible
    pub fn background(peak: &RawPeak) -> Self {
        Self {
            experimental_mz: peak.mz,
            intensity: peak.intensity,
            charge: peak.charge,
            annotation: Vec::new(),
        }
    }
}

impl PartialOrd for AnnotatedPeak {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AnnotatedPeak {
    /// Use `f64::total_cmp` on `self.mz`
    fn cmp(&self, other: &Self) -> Ordering {
        self.experimental_mz
            .value
            .total_cmp(&other.experimental_mz.value)
    }
}

impl PartialEq for AnnotatedPeak {
    /// Use `f64::total_cmp` on all fields to detect total equality
    fn eq(&self, other: &Self) -> bool {
        self.experimental_mz
            .value
            .total_cmp(&other.experimental_mz.value)
            == Ordering::Equal
            && self.intensity.total_cmp(&other.intensity) == Ordering::Equal
            && self.charge.value.total_cmp(&other.charge.value) == Ordering::Equal
            && self.annotation == other.annotation
    }
}

impl Eq for AnnotatedPeak {}
