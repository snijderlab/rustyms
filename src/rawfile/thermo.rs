use std::{fs::File, path::Path};

use entab::parsers::thermo::thermo_raw::{ThermoRawReader, ThermoRawRecord};

/// testing code to see if opening Thermo raw files is feasible.
///
/// # Panics
/// * When the file could not be opened
/// * When the file could not be parsed by [entab]
/// * When a record is incorrect
///
/// # Errors
/// None yet, very rudimentary code
pub fn open(path: impl AsRef<Path>) -> Result<(), String> {
    let file = File::open(path).map_err(|e| format!("Could not open thermo file: {e}"))?;
    let mut reader =
        ThermoRawReader::new(file, None).map_err(|e| format!("Could not open thermo file: {e}"))?;
    while let Some(ThermoRawRecord {
        time,
        mz,
        intensity,
    }) = reader
        .next()
        .map_err(|e| format!("Could not get next record: {e}"))?
    {
        println!("{time},{mz},{intensity}");
    }
    Ok(())
}
