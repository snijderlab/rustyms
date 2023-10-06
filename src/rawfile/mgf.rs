//! Handle MGF reader reading
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use uom::num_traits::Zero;

use crate::{
    helper_functions::check_extension,
    spectrum::{RawPeak, RawSpectrum},
    system::{charge::e, f64::*, mass::dalton, mass_over_charge::mz, time::s},
};
use flate2::read::GzDecoder;

/// Open a MGF file and return the contained spectra.
///
/// # Errors
/// It returns an error when:
/// * The file could not be opened
/// * Any line in the file could not be read
/// * When any expected number in the file is not a number
/// * When there is only one column (separated by space or tab) on a data row
pub fn open(path: impl AsRef<Path>) -> Result<Vec<RawSpectrum>, String> {
    let path = path.as_ref();
    let file =
        BufReader::new(File::open(path).map_err(|err| format!("Could not open file: {err}"))?);
    if check_extension(path, "gz") {
        open_raw(BufReader::new(GzDecoder::new(file)))
    } else {
        open_raw(BufReader::new(file))
    }
}

/// Open a MGF file and return the contained spectra. Open it from a raw buffered reader.
///
/// # Errors
/// It returns an error when:
/// * The file could not be opened
/// * Any line in the file could not be read
/// * When any expected number in the file is not a number
/// * When there is only one column (separated by space or tab) on a data row
#[allow(clippy::missing_panics_doc)]
pub fn open_raw<T: std::io::Read>(reader: BufReader<T>) -> Result<Vec<RawSpectrum>, String> {
    let mut current = RawSpectrum::default();
    let mut output = Vec::new();
    for (linenumber, line) in reader.lines().enumerate() {
        let linenumber = linenumber + 1;
        let line = line.map_err(|_| "Error while reading line")?;
        match line.as_str() {
            "BEGIN IONS" | "" => (),
            "END IONS" => {
                output.push(current);
                current = RawSpectrum::default();
            }
            t if t.contains('=') => {
                // THe previous line made sure it will always contain an equals sign
                let (key, value) = t.split_once('=').unwrap();
                match key {
                    "PEPMASS" => match value.split_once(' ') {
                        None => {
                            current.mass = Mass::new::<dalton>(value.parse().map_err(|_| {
                                format!("Not a number {key} for PEPMASS on {linenumber}")
                            })?);
                        }
                        Some((mass, intensity)) => {
                            current.mass = Mass::new::<dalton>(mass.parse().map_err(|_| {
                                format!("Not a number {key} for PEPMASS mass on {linenumber}")
                            })?);
                            current.intensity = Some(intensity.parse().map_err(|_| {
                                format!("Not a number {key} for PEPMASS intensity on {linenumber}")
                            })?);
                        }
                    },
                    "CHARGE" => {
                        current.charge = parse_charge(value).map_err(|_| {
                            format!("Not a number {key} for CHARGE on {linenumber}")
                        })?;
                    }
                    "RT" => {
                        current.rt =
                            Time::new::<s>(value.parse().map_err(|_| {
                                format!("Not a number {key} for RT on {linenumber}")
                            })?);
                    }
                    "RTINSECONDS" => {
                        current.rt =
                            Time::new::<s>(value.parse().map_err(|_| {
                                format!("Not a number {key} for RT on {linenumber}")
                            })?);
                    }
                    "TITLE" => current.title = value.to_owned(),
                    "SEQUENCE" => current.sequence = Some(value.to_owned()),
                    "NUM_SCANS" => {
                        current.num_scans = value.parse().map_err(|_| {
                            format!("Not a number {key} for NUM_SCANS on {linenumber}")
                        })?;
                    }
                    _ => (),
                }
            }
            t if t.contains(' ') || t.contains('\t') => {
                let split = if t.contains(' ') {
                    t.split(' ').collect::<Vec<_>>()
                } else {
                    t.split('\t').collect::<Vec<_>>()
                };
                let mut peak = RawPeak {
                    mz: MassOverCharge::zero(),
                    intensity: 0.0,
                    charge: Charge::new::<e>(1.0),
                };
                if split.len() < 2 {
                    return Err(format!("Not enough columns on line {linenumber}"));
                }
                peak.mz =
                    MassOverCharge::new::<mz>(split[0].parse().map_err(|_| {
                        format!("Not a number {} for MZ on {linenumber}", split[0])
                    })?);
                peak.intensity = split[1].parse().map_err(|_| {
                    format!("Not a number {} for INTENSITY on {linenumber}", split[1])
                })?;
                if split.len() >= 3 {
                    peak.charge = parse_charge(split[2]).map_err(|_| {
                        format!("Not a number {} for CHARGE on {linenumber}", split[2])
                    })?;
                }
                current.spectrum.push(peak);
            }
            _ => {}
        }
    }
    for current in &mut output {
        current
            .spectrum
            .sort_unstable_by(|a, b| a.mz.value.partial_cmp(&b.mz.value).unwrap());
    }
    Ok(output)
}

fn parse_charge(input: &str) -> Result<Charge, ()> {
    if input.ends_with('+') {
        Ok(Charge::new::<e>(
            input.trim_end_matches('+').parse().map_err(|_| ())?,
        ))
    } else if input.ends_with('-') {
        Ok(Charge::new::<e>(
            -input.trim_end_matches('-').parse().map_err(|_| ())?,
        ))
    } else {
        Ok(Charge::new::<e>(input.parse().map_err(|_| ())?))
    }
}

#[test]
fn test_open() {
    let spectra = open(std::env::var("CARGO_MANIFEST_DIR").unwrap() + "/data/example.mgf").unwrap();
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].spectrum.len(), 5);
    assert!(spectra[0].spectrum[0].mz < spectra[0].spectrum[1].mz);
}
