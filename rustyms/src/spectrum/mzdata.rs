use mzdata::{prelude::*, spectrum::RefPeakDataLevel};

use crate::{
    spectrum::{AnnotatableSpectrum, AnnotatedPeak, AnnotatedSpectrum},
    system::MassOverCharge,
    CompoundPeptidoform,
};

impl<S: SpectrumLike> AnnotatableSpectrum for S {
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
                            mz: MassOverCharge::new::<crate::system::mz>(p.neutral_mass), // TODO: This is M (not MH+) which is not very well supported in the current matching
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
                Self::Da(value.get::<crate::system::mz>()) // This is in Thompson (verified with crate author)
            }
            crate::Tolerance::Relative(value) => {
                Self::PPM(value.get::<crate::system::ratio::ppm>())
            }
        }
    }
}
