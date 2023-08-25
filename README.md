[![rustyms documentation](https://docs.rs/rustyms/badge.svg)](https://docs.rs/rustyms)
[![Crates.io](https://img.shields.io/crates/v/rustyms.svg)](https://crates.io/crates/rustyms)

# Match those fragments!

A work in progress peptide fragmentation matching library for rust.

## Features
 - Read pro forma sequences (very close to 'level 2-ProForma compliant', with the intention to fully support the whole spec)
 - Generate theoretical fragments with control over the fragmentation model from any supported pro forma peptide
 - Generate fragments from satellite ions (w, d, and v)
 - Read mgf files
 - Match spectra to the generated fragments
 - Extensive use of `uom` for compile time unit checking
 - Align peptides based on mass (algorithm will be tweaked extensively over time) (see `Stitch` for more information, but the algorithm has changed)