use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use uom::num_traits::Zero;

use crate::{
    spectrum::{RawPeak, RawSpectrum},
    system::{charge::e, f64::*, mass::dalton, mass_over_charge::mz, time::s},
};

pub fn open(path: impl AsRef<Path>) -> Result<Vec<RawSpectrum>, String> {
    let file =
        BufReader::new(File::open(path).map_err(|err| format!("Could not open file: {err}"))?);
    let mut current = RawSpectrum {
        title: String::new(),
        num_scans: 0,
        rt: Time::zero(),
        charge: Charge::zero(),
        mass: Mass::zero(),
        spectrum: Vec::new(),
        intensity: None,
    };
    let mut output = Vec::new();
    for (linenumber, line) in file.lines().enumerate() {
        let linenumber = linenumber + 1;
        let line = line.map_err(|_| "Error while reading line")?;
        match line.as_str() {
            "BEGIN IONS" | "" => (),
            "END IONS" => {
                output.push(current);
                current = RawSpectrum {
                    title: String::new(),
                    num_scans: 0,
                    rt: Time::zero(),
                    charge: Charge::zero(),
                    mass: Mass::zero(),
                    spectrum: Vec::new(),
                    intensity: None,
                }
            }
            t if t.contains('=') => {
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
    let spectra = open("data/example.mgf").unwrap();
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].spectrum.len(), 5);
}

#[test]
fn test_top_down() {
    let spectra = open("data/20201023_L1_VQ1_Tamar002_SA_CI49_T3_S_IgA_prep1_rep1_ETD16_TCEP_strict_LC_combined_extended.mgf").unwrap();
    assert_eq!(spectra.len(), 1);
    dbg!(spectra);
    todo!();
}
