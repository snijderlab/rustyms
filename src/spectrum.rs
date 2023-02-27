use crate::{
    aminoacids::AminoAcid,
    fragment::Fragment,
    model::Model,
    peptide::Peptide,
    system::{f64::*, mass_over_charge::mz},
};

#[derive(Clone, Debug)]
pub struct RawSpectrum {
    pub title: String,
    pub num_scans: u64,
    pub rt: Time,
    pub charge: Charge,
    pub mass: Mass,
    pub spectrum: Vec<RawPeak>,
}

impl RawSpectrum {
    pub fn annotate(
        &self,
        peptide: Peptide,
        theoretical_fragments: Vec<Fragment>,
        model: &Model,
    ) -> AnnotatedSpectrum {
        let mut annotated = AnnotatedSpectrum {
            title: self.title.clone(),
            num_scans: self.num_scans,
            rt: self.rt,
            charge: self.charge,
            mass: self.mass,
            peptide,
            spectrum: Vec::with_capacity(self.spectrum.len()),
        };

        let mut peaks = self.spectrum.clone();

        for fragment in theoretical_fragments {
            let close = peaks
                .iter()
                .enumerate()
                .map(|(i, p)| (i, p.ppm(&fragment)))
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            if let Some((index, ppm)) = close {
                if ppm <= model.ppm {
                    annotated
                        .spectrum
                        .push(AnnotatedPeak::new(&peaks[index], fragment));
                    peaks.remove(index);
                }
            }
        }
        annotated
            .spectrum
            .extend(peaks.iter().map(AnnotatedPeak::background));

        annotated
    }
}

#[derive(Clone, Debug)]
pub struct AnnotatedSpectrum {
    pub title: String,
    pub num_scans: u64,
    pub rt: Time,
    pub charge: Charge,
    pub mass: Mass,
    pub peptide: Peptide,
    pub spectrum: Vec<AnnotatedPeak>,
}

#[derive(Clone, Debug)]
pub struct RawPeak {
    pub mz: MassOverCharge,
    pub intensity: f64,
    pub charge: Charge,
}

impl RawPeak {
    pub fn ppm(&self, fragment: &Fragment) -> MassOverCharge {
        (self.mz - fragment.mz()).abs() / fragment.mz() * MassOverCharge::new::<mz>(1e6)
    }
}

#[derive(Clone, Debug)]
pub struct AnnotatedPeak {
    pub experimental_mz: MassOverCharge,
    pub intensity: f64,
    pub charge: Charge,
    pub annotation: Option<Fragment>,
}

impl AnnotatedPeak {
    pub fn new(peak: &RawPeak, annotation: Fragment) -> Self {
        AnnotatedPeak {
            experimental_mz: peak.mz,
            intensity: peak.intensity,
            charge: peak.charge,
            annotation: Some(annotation),
        }
    }

    pub fn background(peak: &RawPeak) -> Self {
        AnnotatedPeak {
            experimental_mz: peak.mz,
            intensity: peak.intensity,
            charge: peak.charge,
            annotation: None,
        }
    }
}
