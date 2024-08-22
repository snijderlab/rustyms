use mzdata::{prelude::*, spectrum::RefPeakDataLevel};

use crate::{system::MassOverCharge, CompoundPeptidoform};

use super::{AnnotatedPeak, AnnotatedSpectrum, Fragments};

impl<S: SpectrumLike> Fragments for S {
    type Tolerance = Tolerance;

    fn empty_annotated(&self, peptide: CompoundPeptidoform) -> AnnotatedSpectrum {
        AnnotatedSpectrum {
            title: self.description().id.clone(),
            num_scans: self.description().acquisition.scans.len() as u64,
            rt: None,
            charge: None,
            mass: None,
            peptide,
            spectrum: match self.peaks() {
                RefPeakDataLevel::Missing | RefPeakDataLevel::RawData(_) => Vec::new(),
                RefPeakDataLevel::Centroid(data) => data
                    .iter()
                    .map(|p| {
                        AnnotatedPeak::background(&super::RawPeak {
                            mz: MassOverCharge::new::<crate::system::mz>(p.mz),
                            intensity: ordered_float::OrderedFloat(f64::from(p.intensity)),
                        })
                    })
                    .collect(),
                RefPeakDataLevel::Deconvoluted(data) => data
                    .iter()
                    .map(|p| {
                        AnnotatedPeak::background(&super::RawPeak {
                            mz: MassOverCharge::new::<crate::system::mz>(p.neutral_mass), // TODO: Ask if this is M or MH+
                            intensity: ordered_float::OrderedFloat(f64::from(p.intensity)),
                        })
                    })
                    .collect(),
            },
        }
    }

    fn search(
        &self,
        query: crate::system::MassOverCharge,
        tolerance: Self::Tolerance,
    ) -> Option<usize> {
        self.peaks().search(query.value, tolerance)
    }
}

impl From<crate::Tolerance<MassOverCharge>> for Tolerance {
    fn from(value: crate::Tolerance<MassOverCharge>) -> Self {
        match value {
            crate::Tolerance::Absolute(value) => {
                Self::Da(value.get::<crate::system::mass_over_charge::mz>()) // TODO: Ask is this is actually Th or Da
            }
            crate::Tolerance::Relative(value) => {
                Self::PPM(value.get::<crate::system::ratio::ppm>())
            }
        }
    }
}
