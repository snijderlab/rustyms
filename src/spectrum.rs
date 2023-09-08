use uom::num_traits::Zero;

use crate::{
    fragment::Fragment,
    system::{f64::*, mass_over_charge::mz},
    ComplexPeptide,
};

/// A raw spectrum (meaning not annotated yet)
#[derive(Clone, Debug)]
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
    pub spectrum: Vec<RawPeak>,
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
            .map(|p| p.intensity)
            .reduce(f64::max)
            .unwrap();
        self.spectrum
            .retain(|p| p.intensity >= max * filter_threshold);
        self.spectrum.shrink_to_fit();
    }

    /// Annotate this spectrum with the given peptide and given fragments see [`Peptide::generate_theoretical_fragments`].
    ///
    /// # Panics
    /// If any fragment does not have a defined m/z
    pub fn annotate(
        &self,
        peptide: ComplexPeptide,
        theoretical_fragments: &[Fragment],
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
                .binary_search_by(|p| p.mz.partial_cmp(&fragment.mz().unwrap()).unwrap())
                .map_or_else(|i| i, |i| i);

            // Check index-1, index and index+1 (if existing) to find the one with the lowest ppm
            let mut closest = (0, f64::INFINITY);
            for i in (index - 1).max(0)..=(index + 1).min(self.spectrum.len() - 1) {
                let ppm = self.spectrum[i].ppm(fragment).unwrap().value;
                if ppm < closest.1 {
                    closest = (i, ppm);
                }
            }

            annotated.spectrum[closest.0]
                .annotation
                .push(fragment.clone());
        }

        annotated
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
        }
    }
}

/// An annotated spectrum
#[derive(Clone, Debug)]
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
    pub spectrum: Vec<AnnotatedPeak>,
}

/// A raw peak
#[derive(Clone, Debug)]
pub struct RawPeak {
    /// The mz value of this peak
    pub mz: MassOverCharge,
    /// The intensity of this peak
    pub intensity: f64,
    /// The charge of this peak
    pub charge: Charge, // #TODO: Is this item needed? (mgf has it, not used in rustyms)
}

impl RawPeak {
    /// Determine the ppm error for the given fragment, optional because the mz of a [Fragment] is optional
    pub fn ppm(&self, fragment: &Fragment) -> Option<MassOverCharge> {
        Some(MassOverCharge::new::<mz>(self.mz.ppm(fragment.mz()?)))
    }
}

/// An annotated peak
#[derive(Clone, Debug)]
pub struct AnnotatedPeak {
    /// The experimental mz
    pub experimental_mz: MassOverCharge,
    /// The experimental intensity
    pub intensity: f64,
    /// The charge
    pub charge: Charge, // #TODO: Is this item needed? (mgf has it, not used in rustyms)
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
