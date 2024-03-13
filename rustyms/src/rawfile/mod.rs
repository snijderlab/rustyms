//! Handling raw files

use std::path::PathBuf;

use ordered_float::OrderedFloat;

use crate::spectrum::{RawPeak, RawSpectrum};

use crate::system::time::min;
use crate::system::{dalton, e, Charge, Mass, MassOverCharge, Time};
pub mod mgf;

pub fn thermo_open(path: impl Into<PathBuf>) -> Result<Vec<RawSpectrum>, String> {
    let data = thermorawfilereader::open(path).map_err(|err| err.to_string())?;
    Ok(data.iter().take(100).map(|s| s.into()).collect())
}

impl From<thermorawfilereader::RawSpectrum> for RawSpectrum {
    fn from(value: thermorawfilereader::RawSpectrum) -> Self {
        let mut spectrum = RawSpectrum::default();
        if let Some(precursor) = value.precursor() {
            spectrum.charge = Charge::new::<e>(precursor.charge() as f64);
            spectrum.mass = Mass::new::<dalton>(precursor.mz() / precursor.charge() as f64);
        }
        spectrum.raw_index = Some(value.index());
        spectrum.rt = Time::new::<min>(value.time());
        if let Some(data) = value.data() {
            if let (Some(mz), Some(intensity)) = (data.mz(), data.intensity()) {
                spectrum.extend(
                    mz.iter()
                        .zip(intensity.iter())
                        .map(|(mz, intensity)| RawPeak {
                            mz: MassOverCharge::new::<crate::system::mz>(mz),
                            intensity: OrderedFloat::<f64>(intensity as f64),
                            charge: Charge::default(),
                        }),
                );
            }
        }
        spectrum
    }
}

#[test]
fn test_thermo_open() {
    dbg!(thermo_open("data/20191211_F1_Ag5_peng0013_SA_HA_tryp.raw"));
    panic!("WHEEE")
}
