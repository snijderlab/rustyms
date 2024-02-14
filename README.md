[![rustyms documentation](https://docs.rs/rustyms/badge.svg)](https://docs.rs/rustyms)
[![Crates.io](https://img.shields.io/crates/v/rustyms.svg)](https://crates.io/crates/rustyms)

# Match those fragments!

A peptide fragmentation matching library for rust. Split into multiple smaller crates to help keep things organised.

## Features
 - Read pro forma sequences ('level 2-ProForma + mass spectrum compliant + glycans compliant', with the intention to fully support the whole spec)
 - Generate theoretical fragments with control over the fragmentation model from any supported pro forma peptide
   - Generate fragments from satellite ions (w, d, and v)
   - Generate glycan fragments
   - Generate theoretical fragments for modifications of unknown position
   - Generate theoretical fragments for chimeric spectra
 - Read mgf files
 - Match spectra to the generated fragments
 - Extensive use of `uom` for compile time unit checking
 - Align peptides based on mass (algorithm will be tweaked extensively over time) (see `Stitch` for more information, but the algorithm has been improved)

# Contributing

Any contribution is welcome (especially adding/fixing documentation as that is very hard to do as main developer). 