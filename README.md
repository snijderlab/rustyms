[![rustyms documentation](https://docs.rs/rustyms/badge.svg)](https://docs.rs/rustyms)
[![Crates.io](https://img.shields.io/crates/v/rustyms.svg)](https://crates.io/crates/rustyms)

# Match those fragments!

A peptide fragmentation matching library for rust.

```rust
use rustyms::*;
use rustyms::system::*;
// Parse a peptide
let peptide = ComplexPeptide::pro_forma("VAEINPSNGGTTFNEKFKGGKATJ").unwrap();
// Get the theoretical fragmentation for EThcD
let fragments = peptide.generate_theoretical_fragments(Charge::new::<e>(4.0), &Model::ethcd());
// Load the raw file
let spectra = rawfile::mgf::open("data/annotated_example.mgf").unwrap();
// Annotate the spectrum with this peptide
let matched = spectra[0].annotate(peptide, &fragments, &Model::ethcd(), MassMode::Monoisotopic);
```

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