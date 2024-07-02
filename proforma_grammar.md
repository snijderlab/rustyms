# ProForma support

The ProForma specification is fully supported ('level 2-ProForma + top-down compliant + cross-linking compliant + glycans compliant + mass spectrum compliant'). If any valid ProForma sequence is not allowed that is an error of the library. This is also true for any panic (that is not explicitly mentioned in the function documentation).

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
| Unimod  | Yes           | Yes      | Yes   | Yes                              | Yes                                |
| PSI-MOD | Yes           | Yes      | Yes   | NA                               | Yes                                |
| RESID   | Yes           | Yes      | Yes   | NA                               | Yes                                |
| XL-MOD  | Yes           | Yes      | Yes   | Yes                              | Yes                                |
| GNO     | Yes           | NA       | NA    | NA (solved for all glycans)      | NA                                 |

Note some modifications that do not fit the assumptions of rustyms might be missing from the ontologies. Examples of these are cross-links with more then 2 positions from XL-MOD and RESID, and modifications with different diff_formulas based on which location they bound from RESID.