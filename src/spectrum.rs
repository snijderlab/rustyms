use crate::{aminoacids::AminoAcid, fragment::Fragment, system::f64::*};

pub struct RawSpectrum {
    pub title: String,
    pub num_scans: u64,
    pub rt: Time,
    pub charge: Charge,
    pub mass: Mass,
    pub spectrum: Vec<RawPeak>,
}

pub struct AnnotatedSpectrum {
    pub title: String,
    pub num_scans: u64,
    pub rt: Time,
    pub charge: Charge,
    pub mass: Mass,
    pub peptide: Option<Vec<AminoAcid>>,
    pub spectrum: Vec<AnnotatedPeak>,
}

pub struct RawPeak {
    pub mz: MassOverCharge,
    pub intensity: f64,
    pub charge: Charge,
}
pub struct AnnotatedPeak {
    pub experimental_mz: MassOverCharge,
    pub intensity: f64,
    pub charge: Charge,
    pub annotation: Option<Fragment>,
}
