use std::path::Path;

use crate::spectrum::RawSpectrum;

pub fn open(path: &Path) -> Result<RawSpectrum, ()> {
    let file = BufReader::new(File::open(path)?);
    for (linenumber, read_line) in file.lines().enumerate() {
        let linenumber = linenumber + 1;

    }
}