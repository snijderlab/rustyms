//! Raw spectra (not annotated)

use std::cmp::Ordering;

use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use uom::num_traits::Zero;

use crate::{
    fragment::Fragment,
    itertools_extension::ItertoolsExt,
    system::{
        f64::{Mass, MassOverCharge, Ratio, Time},
        usize::Charge,
    },
    AnnotatedSpectrum, CompoundPeptidoform, MassMode, Model, WithinTolerance,
};

use super::{AnnotatedPeak, PeakSpectrum};

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
    pub fn relative_noise_filter(&mut self, filter_threshold: f64) {
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

    /// Filter the spectrum to retain all with an intensity above `filter_threshold`.
    pub fn absolute_noise_filter(&mut self, filter_threshold: f64) {
        self.spectrum.retain(|p| *p.intensity >= filter_threshold);
        self.spectrum.shrink_to_fit();
    }

    /// Filter a spectrum by dividing it in windows and within each window only retain the `top` number of peaks.
    #[allow(clippy::missing_panics_doc)] // Cannot panic as it checks with peek first
    pub fn top_x_filter(&mut self, window_size: f64, top: usize) {
        let mut new_spectrum = Vec::with_capacity(
            self.spectrum
                .last()
                .and_then(|l| self.spectrum.first().map(|f| (f, l)))
                .map(|(f, l)| ((l.mz.value - f.mz.value) / window_size).round() as usize * top)
                .unwrap_or_default(),
        );
        let mut spectrum = self.spectrum.iter().cloned().peekable();
        let mut window = 1;
        let mut peaks = Vec::new();

        while spectrum.peek().is_some() {
            while let Some(peek) = spectrum.peek() {
                if peek.mz.value <= f64::from(window) * window_size {
                    peaks.push(spectrum.next().unwrap());
                } else {
                    break;
                }
            }
            new_spectrum.extend(
                peaks
                    .iter()
                    .cloned()
                    .k_largest_by(top, |a, b| a.mz.value.total_cmp(&b.mz.value)),
            );
            peaks.clear();
            window += 1;
        }

        self.spectrum = new_spectrum;
    }

    /// Annotate this spectrum with the given peptide and given fragments see [`crate::ComplexPeptide::generate_theoretical_fragments`].
    ///
    /// # Panics
    /// If any fragment does not have a defined m/z
    pub fn annotate(
        &self,
        peptide: CompoundPeptidoform,
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

            if model
                .tolerance
                .within(&self.spectrum[closest.0].mz, &fragment.mz(mode))
            {
                annotated.spectrum[closest.0]
                    .annotation
                    .push((fragment.clone(), Vec::new()));
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
            charge: Charge::new::<crate::system::e>(1),
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
            && self.charge.value.cmp(&other.charge.value) == Ordering::Equal
    }
}

impl Eq for RawPeak {}

impl RawPeak {
    /// Determine the ppm error for the given fragment
    pub fn ppm(&self, fragment: &Fragment, mode: MassMode) -> Ratio {
        self.mz.ppm(fragment.mz(mode))
    }
}
