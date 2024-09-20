# rustyms-generate-databases

Generated all databases for inclusion in rustyms. This has to be a separate project to speed up compilation times and because the GNOme database is too big to be submitted to crates.io.

## Sources for the downloaded files

- PSI-MOD: https://github.com/HUPO-PSI/psi-mod-CV (2021-06-13 v1.031.6)
- Unimod: http://www.unimod.org/obo/unimod.obo (2024-08-12 11:33)
- RESID: ftp://ftp.proteininformationresource.org/pir_databases/other_databases/resid/ (2018-04-31 RESIDUES.XML)
- XL-MOD: https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/cv/XLMOD.obo (2021-03-23 1.1.12)
- GNO: http://purl.obolibrary.org/obo/gno.obo (2024-05-21) structures: https://glycosmos.org/download/ ('List of all GlyCosmos Glycans data.') (downloaded 2024-06-16)
- Isotopic atomic masses: https://ciaaw.org/data/IUPAC-atomic-masses.csv (2021-03-17)

## Ontologies

| Name    | Modifications | Numbered | Rules | Diagnostic ions / neutral losses | Description / synonyms / cross ids |
| ------- | ------------- | -------- | ----- | -------------------------------- | ---------------------------------- |
| Unimod  | Yes           | Yes      | Yes   | Yes                              | Yes                                |
| PSI-MOD | Yes           | Yes      | Yes   | NA                               | Yes                                |
| RESID   | Yes           | Yes      | Yes   | NA                               | Yes                                |
| XL-MOD  | Yes           | Yes      | Yes   | Yes                              | Yes                                |
| GNO     | Yes           | NA       | NA    | NA (solved for all glycans)      | Yes                                |

Note some modifications that do not fit the assumptions of rustyms might be missing from the ontologies. Examples of these are cross-links with more then 2 positions from XL-MOD and RESID, and modifications with different diff_formulas based on which location they bound from RESID. Additionally only the Glycans of a specific mass or with a structure.