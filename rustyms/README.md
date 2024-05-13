# Match those fragments!

Handle mass spectrometry data in Rust. This crate is set up to handle very complex peptides with
loads of ambiguity and complexity. It pivots around the [`CompoundPeptidoform`], [`Peptidoform`] and [`LinearPeptide`]
which encode the [ProForma](https://github.com/HUPO-PSI/ProForma) specification. Additionally
this crate enables the reading of [mgf](rawfile::mgf), doing [spectrum annotation](RawSpectrum::annotate)
(BU/MD/TD), finding [isobaric sequences](find_isobaric_sets), doing [alignments of peptides](align::align)
, accessing the [IMGT germline database](imgt), and [reading identified peptide files](identification).

## Library features
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

## Example usage
```rust
# fn main() -> Result<(), rustyms::error::CustomError> {
# let raw_file_path = "data/annotated_example.mgf";
// Open some data and see if the given peptide is a valid match
use rustyms::{*, system::{usize::Charge, e}};
let peptide = CompoundPeptidoform::pro_forma("Q[Gln->pyro-Glu]VQEVSERTHGGNFD", None)?;
let spectrum = rawfile::mgf::open(raw_file_path)?;
let model = Model::ethcd();
let fragments = peptide.generate_theoretical_fragments(Charge::new::<e>(2), &model);
let annotated = spectrum[0].annotate(peptide, &fragments, &model, MassMode::Monoisotopic);
let fdr = annotated.fdr(&fragments, &model);
// This is the incorrect sequence for this spectrum so the FDR will indicate this
# dbg!(&fdr, fdr.sigma(), fdr.fdr(), fdr.score());
assert!(fdr.sigma() < 2.0);
# Ok(()) }
```
```rust
# fn main() -> Result<(), rustyms::error::CustomError> {
// Check how this peptide compares to a similar peptide (using `align`)
// (same sequence, repeated for easy reference)
use rustyms::{*, align::*};
let first_peptide = LinearPeptide::pro_forma("Q[Gln->pyro-Glu]VQEVS", None)?.simple().unwrap();
let second_peptide = LinearPeptide::pro_forma("E[Glu->pyro-Glu]VQVES", None)?.simple().unwrap();
let alignment = align::<4, Simple, Simple>(&first_peptide, &second_peptide,
                 matrix::BLOSUM62, Tolerance::new_ppm(10.0), AlignType::GLOBAL);
# dbg!(&alignment);
let stats = alignment.stats();
# //assert_eq!(stats.identical, 3); // Only three positions are identical
assert_eq!(stats.mass_similar, 6); // All positions are mass similar
# Ok(()) }
```

## Compilation features
Rustyms ties together multiple smaller modules into one cohesive structure.
It has multiple features which allow you to slim it down if needed (all are enabled by default).
* `identification` - gives access to methods reading many different identified peptide formats.
* `align` - gives access to mass based alignment of peptides.
* `imgt` - enables access to the IMGT database of antibodies germline sequences, with annotations.
* `rayon` - enables parallel iterators using rayon, mostly for `imgt` but also in consecutive
  align.
* `isotopes` - gives access to generation of an averagine model for isotopes, also enables two additional dependencies
* `rand` - allows the generation of random peptides