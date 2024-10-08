# Match those fragments!

Handle mass spectrometry data in Rust. This crate is set up to handle very complex peptides with
loads of ambiguity and complexity. It pivots around the [`CompoundPeptidoform`], [`Peptidoform`] and [`LinearPeptide`]
which encode the [ProForma](https://github.com/HUPO-PSI/ProForma) specification. Additionally
this crate enables the reading of [mgf](rawfile::mgf), doing [spectrum annotation](RawSpectrum::annotate)
(BU/MD/TD), finding [isobaric sequences](find_isobaric_sets), doing [alignments of peptides](align::align)
, accessing the [IMGT germline database](imgt), and [reading identified peptide files](identification).

## Library features

 - Read [ProForma](https://github.com/HUPO-PSI/ProForma) sequences (complete specification supported: 'level 2-ProForma + top-down compliant + cross-linking compliant + glycans compliant + mass spectrum compliant')
 - Generate theoretical fragments with control over the fragmentation model from any ProForma peptidoform/proteoform
   - Generate theoretical fragments for chimeric spectra
   - Generate theoretical fragments for cross-links (also disulfides)
   - Generate theoretical fragments for modifications of unknown position
   - Generate peptide backbone (a, b, c, x, y, and z) and satellite ion fragments (w, d, and v)
   - Generate glycan fragments (B, Y, and internal fragments)
 - Integrated with [mzdata](https://crates.io/crates/mzdata) for reading raw data files
 - Match spectra to the generated fragments
 - [Align peptides based on mass](https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00188)
 - Fast access to the IMGT database of antibody germlines
 - Reading of multiple identified peptide file formats (Fasta, MaxQuant, MSFragger, Novor, OPair, Peaks, and Sage)
 - Exhaustively fuzz tested for reliability (using [cargo-afl](https://crates.io/crates/cargo-afl))
 - Extensive use of [uom](https://docs.rs/uom/latest/uom/) for compile time unit checking

## Example usage

```rust
# fn main() -> Result<(), rustyms::error::CustomError> {
# let raw_file_path = "data/annotated_example.mgf";
use rustyms::{*, system::{usize::Charge, e}};
// Open example raw data (this is the built in mgf reader, look into mzdata for more advanced raw file readers)
let spectrum = rawfile::mgf::open(raw_file_path)?;
// Parse the given ProForma definition
let peptide = CompoundPeptidoform::pro_forma("[Gln->pyro-Glu]-QVQEVSERTHGGNFD", None)?;
// Generate theoretical fragments for this peptide given EThcD fragmentation
let model = Model::ethcd();
let fragments = peptide.generate_theoretical_fragments(Charge::new::<e>(2), &model);
// Annotate the raw data with the theoretical fragments
let annotated = spectrum[0].annotate(peptide, &fragments, &model, MassMode::Monoisotopic);
// Calculate a peak false discovery rate for this annotation 
let (fdr, _) = annotated.fdr(&fragments, &model, MassMode::Monoisotopic);
// This is the incorrect sequence for this spectrum so the peak FDR will indicate this
# dbg!(&fdr, fdr.peaks_sigma(), fdr.peaks_fdr(), fdr.peaks_score());
assert!(fdr.peaks_sigma() > 2.0);
# Ok(()) }
```

```rust
# fn main() -> Result<(), rustyms::error::CustomError> {
use rustyms::{*, align::*};
// Check how this peptide compares to a similar peptide (using the feature `align`)
let first_peptide = LinearPeptide::pro_forma("IVQEVS", None)?.into_simple_linear().unwrap();
let second_peptide = LinearPeptide::pro_forma("LEVQVES", None)?.into_simple_linear().unwrap();
// Align the two peptides using mass based alignment
// I-VQEVS A
// LEVQVES B
// ─+  ╶─ 
let alignment = align::<4, SimpleLinear, SimpleLinear>(&first_peptide, &second_peptide,
                 matrix::BLOSUM62, Tolerance::new_ppm(10.0), AlignType::GLOBAL);
# dbg!(&alignment);
// Calculate some more statistics on this alignment
let stats = alignment.stats();
assert_eq!(stats.mass_similar, 6); // 6 out of the 7 positions are mass similar
assert_eq!(stats.gaps, 1); // 1 position is an insertion
# Ok(()) }
```

## Compilation features

Rustyms ties together multiple smaller modules into one cohesive structure.
It has multiple features which allow you to slim it down if needed (all are enabled by default).
* `align` - gives access to mass based alignment of peptides.
* `identification` - gives access to methods reading many different identified peptide formats.
* `imgt` - enables access to the IMGT database of antibodies germline sequences, with annotations.
* `isotopes` - gives access to generation of an averagine model for isotopes, also enables two additional dependencies.
* `rand` - allows the generation of random peptides.
* `rayon` - enables parallel iterators using rayon, mostly for `imgt` but also in consecutive align.
* `mzdata` - enables integration with [mzdata](https://github.com/mobiusklein/mzdata) which has more advanced raw file support.

## Sources for the downloaded files

- PSI-MOD: https://github.com/HUPO-PSI/psi-mod-CV (2021-06-13 v1.031.6)
- Unimod: http://www.unimod.org/obo/unimod.obo (2024-08-12 11:33)
- RESID: ftp://ftp.proteininformationresource.org/pir_databases/other_databases/resid/ (2018-04-31 RESIDUES.XML)
- XL-MOD: https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/cv/XLMOD.obo (2021-03-23 1.1.12)
- GNO: http://purl.obolibrary.org/obo/gno.obo (2024-05-21) structures: https://glycosmos.org/download/ ('List of all GlyCosmos Glycans data.') (downloaded 2024-07-02)
  - To save space (crates.io has a hard limit on crate size) the unused columns of the structures csv are remove (only 0 and 1 are kept) and the `gno.obo` is trimmed using the following regex: `(property_value: GNO:00000(022|023|041|042|101|102) .*$\n)|(def: .*$\n)|(synonym: .*$\n)|(name: [^ ]*$\n)` (any matching line is removed) and the following replacement regex `(is_a: [^ ]*) ! .*\n` with `$1\n`.
  - The structures csv file has only the first two columns kept for the same reason, also remove the two lines starting with `"`
- Isotopic atomic masses: https://ciaaw.org/data/IUPAC-atomic-masses.csv (2021-03-17)

## Ontologies

| Name    | Modifications | Numbered | Rules | Diagnostic ions / neutral losses | Description / synonyms / cross ids |
| ------- | ------------- | -------- | ----- | -------------------------------- | ---------------------------------- |
| Unimod  | Yes           | Yes      | Yes   | Yes                              | Yes                                |
| PSI-MOD | Yes           | Yes      | Yes   | NA                               | Yes                                |
| RESID   | Yes           | Yes      | Yes   | NA                               | Yes                                |
| XL-MOD  | Yes           | Yes      | Yes   | Yes                              | Yes                                |
| GNO     | Yes           | NA       | NA    | NA (solved for all glycans)      | NA                                 |

Note some modifications that do not fit the assumptions of rustyms might be missing from the ontologies. Examples of these are cross-links with more then 2 positions from XL-MOD and RESID, and modifications with different diff_formulas based on which location they bound from RESID. Additionally only the Glycans of a specific mass or with a structure.