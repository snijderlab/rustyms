# ProForma support

## Levels of support

1. Base Level Support (Technical name: Base-ProForma Compliant)
   Represents the lowest level of compliance, this level involves providing support for:

   - [x] 4.1 Amino acid sequences
   - [x] 4.2.1/4.2.2 Protein modifications using two of the supported CVs/ontologies: Unimod and PSI-MOD.
     - [x] Prepare Unimod and parse modifications
     - [x] Prepare PSI-MOD and parse modifications
   - [x] 4.2.6 Protein modifications using delta masses (without prefixes)
   - [x] 4.3 N-terminal, C-terminal and labile modifications.
   - [x] 4.4 Ambiguity in the modification position, including support for localisation scores.
     - [x] Global: `[Phospho]^2?`
     - [x] Local/named: `[Phospho#g1(0.01)]`
     - [x] Stretch: `(AAA)[Phospho]` (note: also handle multiple modifications on a group + scores)
   - [x] 4.7 Ambiguity in the amino acid sequence. `(?DQ)N`
   - [x] 4.8 INFO tag.

2. Additional Separate Support (Technical name: level 2-ProForma compliant)
   These features are independent from each other:

   - [x] 4.1 Unusual amino acids (O and U).
   - [x] 4.1 Ambiguous amino acids (e.g. X, B, Z). This would include support for sequence tags of known mass (using the character X).
   - [x] 4.2.6 Protein modifications using delta masses (using prefixes for the different CVs/ontologies).
   - [x] 4.2.1 Use of prefixes for Unimod (U:) and PSI-MOD (M:) names.
   - [x] 4.9 Support for the joint representation of experimental data and its interpretation. (see 4.9)

3. Top-Down Extensions (Technical name: level 2-ProForma + top-down compliant)

   - [ ] 4.2.1 Additional CV/ontologies for protein modifications: RESID (the prefix R MUST be used for RESID CV/ontology term names)
   - [x] 4.2.8 Chemical formulas (this feature occurs in two places in this list).

4. Cross-Linking Extensions (Technical name: level 2-ProForma + cross-linking compliant)

   - [x] 4.2.1/4.2.3 Cross-linked peptides (using the XL-MOD CV/ontology, the prefix X MUST be used for XL-MOD CV/ontology term names).

5. Glycan Extensions (Technical name: level 2-ProForma + glycans compliant)

   - [x] 4.2.5 Additional CV/ontologies for protein modifications: GNO (the prefix G MUST be used for GNO CV/ontology term names)
   - [x] 4.2.9 Glycan composition.
   - [x] 4.2.8 Chemical formulas (this feature occurs in two places in this list).

6. Spectral Support (Technical name: level 2-ProForma + mass spectrum compliant)

    - [x] 7 Charge and chimeric spectra are special cases (see Appendix II).
        - [x] Parse chimeric spectra
        - [x] Annotate chimeric spectra
        - [x] Parse charge state/adduct ions
        - [x] Handle/annotate charge state/adduct ions
    - [x] 4.6 Global modifications (e.g., every C is C13).
        - [x] 4.6.1 Isotope modifications `<15N>` (applied as a post filter after generating the MolecularFormula for generated fragments and full peptide, not applied on single aminoacids)
        - [x] 4.6.2 Fixed modifications `<[Formula:H-1]@A>`

## Other stuff to implement still

- [x] 4.5 Handle multiple modifications on a single position (or amino acid group)
- [x] 4.2.1.1 Use the Unimod obo file instead of xml
- [x] Keep track of the original mod definition to show it nicer to the user?
- [x] 4.2.7 Gap of known mass
- [x] Show ambiguous amino acids back to the user
- [x] Handle multiple possible backbone masses (fragmenting modifications eg glycans, poorly localised modifications, ..)
- [x] Enforce ontology modification placement rules
- [x] Handle isotopes of elements, amongst others for the missing mods of PSI-MOD and global modifications
- [x] Better error data, allowing the construction of rust-like error messages
- [x] Modification specific neutral losses, phospho: fully deleted
- [x] Modification diagnostic peaks
- [x] 4.2.4 Branched peptides, similar to cross linking

## Sources for the downloaded files

- PSI-MOD: https://github.com/HUPO-PSI/psi-mod-CV (v1.031.6)
- Unimod: http://www.unimod.org/obo/unimod.obo (29:02:2024 10:49)
- RESID: ftp://ftp.proteininformationresource.org/pir_databases/other_databases/resid/ (RESIDUES.XML)
- XL-MOD: https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/cv/XLMOD.obo
- GNO: http://purl.obolibrary.org/obo/gno.obo (2024-03-18) structures: https://glycosmos.org/download/glycosmos_glycans_list.csv (2024-03-18)
  - To save space (crates.io has a hard limit on crate size) the unused columns of the structures csv are remove (only 0 and 1 are kept) and the `gno.obo` is trimmed using the following regex: `(property_value: GNO:00000(022|023|041|042|101|102) .*$\n)|(def: .*$\n)` (any matching line is removed)
  - The structures csv file has only the first two columns kept for the same reason, also remove the two lines starting with `"`
- Isotopic atomic masses: https://ciaaw.org/data/IUPAC-atomic-masses.csv (2021-03-17)

## Ontologies

| Name    | Modifications | Numbered | Rules | Diagnostic ions / neutral losses | Description / synonyms / cross ids |
| ------- | ------------- | -------- | ----- | -------------------------------- | ---------------------------------- |
| Unimod  | Yes           | Yes      | Yes   | Yes                              | -                                  |
| PSI-MOD | Yes           | Yes      | Yes   | NA                               | -                                  |
| RESID   | -             | -        | -     | -                                | -                                  |
| XL-MOD  | YES           | YES      | YES   | YES                              | Yes                                |
| GNO     | Yes           | NA       | NA    | NA (solved for all glycans)      | NA                                 |
