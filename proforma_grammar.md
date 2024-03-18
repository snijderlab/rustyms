# ProForma support

## Grammar

Here is a full grammar of the pro forma syntax. `\` is used to escape characters. `<name>` is used to refer to another rule. `()` are used to group rules. `|` is used to give alternate options. `*+?` have their common regex meaning. `()x*` indicates a group that is repeated and separated by `x`, any repeating operator can be used.

```
<full_definition> = (<full_sequence>)++
<full_sequence> = <sequence>(//<sequence>)?(\\\\<sequence>)?
<sequence> = <global>*<labile>*<unknown_position_modification>*<n-term>?<aa_sequence><c-term>?(/<charge>([(<charge><atom>(+|-)),*])?)?
<aa_sequence> = (<aa_sequence_inner>|<range>)*
<range> = \(<aa_sequence_inner>\)<modification>
<aa_sequence_inner> = (<aa_sequence_inner_inner>)+|\(\?<aa_sequence_inner_inner>\)
<aa_sequence_inner_inner> = <aminoacid><modification>*
<global> = \<<global_isotope!!>|<modification>@<aminoacids>\>
<global_isotope> = <num><atom>
<aminoacids> = <aminoacid>(,<aminoacid>)*
<unknown_position_modification> = <modification>(^<num>)?\?
<modification> = [<modification_inner>]
<labile> = {<modification_inner>}(^<num>)?
<n-term> = [<modification_inner>]-
<c-term> = -[<modification_inner>]
<modification_inner> = (<modification_inner_inner>|<cross_link>|<branch>)(\|<modification_inner_inner>)*
<modification_inner_inner> = ((<named_modification>|<numbered_modification>|<delta_mass>|<formula>|<glycan>)(#<name>(\(<decimal>\))?)?)|(\|INFO:<info>)
<named_modification> = ((U:)?|(M:)?)<name>|(R:|X:|G:)<name>
<numbered_modification> = (UNIMOD|MOD|RESID|XLMOD|GNO):<num> ?? should be name?
<cross_link> = <modification_inner_inner>?#XL<name>
<branch> = <modification_inner_inner>?#BRANCH
<delta_mass> = ((U|M|R|X|G|Obs):)?<shift>
<formula> = Formula:<pair>+
<pair> = <atom><cardinality>?|[<num><atom><cardinality>?]
<glycan> = Glycan:(<monosaccharide><num>?)+
<charge> =regex= (+|-)?[0-9]+
<cardinality> =regex= (-)?[0-9]+
<shift> =regex= (+|-)[0-9]+\.[0-9]+
<name> =regex= [A-Za-z0-9]+
<num> =regex= [0-9]+
<decimal> =regex= [0-9]+\.[0-9]+
<atom> =human= any valid atom name
<aminoacid> =human= any IUPAC amino acid character or any UniProt ambiguous/unusual
<info> =human= any piece of text, with the exclusion of unpaired brackets
```

## Open questions

- proper order for global/labile/unknown_position and for additional peptides (xl/cys-xl/branched/chimeric) for MS?
  - A: `<GLOBAL_MOD>[UNKNOWN_POS]?{LABILE_MOD}[N_TERM]-PEPTIDE-[C_TERM]` https://github.com/HUPO-PSI/ProForma/issues/3#issuecomment-906448694
- what is the mass for B/Z?
- what defines valid ionic species for ion charge?
- for ionic charge/adduct ions, is there a way to specify higher charged ionic species? (Ca+2 or Ca2+?)

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

1. Additional Separate Support (Technical name: level 2-ProForma compliant)
   These features are independent from each other:

   - [x] 4.1 Unusual amino acids (O and U).
   - [x] 4.1 Ambiguous amino acids (e.g. X, B, Z). This would include support for sequence tags of known mass (using the character X).
   - [x] 4.2.6 Protein modifications using delta masses (using prefixes for the different CVs/ontologies).
   - [x] 4.2.1 Use of prefixes for Unimod (U:) and PSI-MOD (M:) names.
   - [x] 4.9 Support for the joint representation of experimental data and its interpretation. (see 4.9)

1. Top-Down Extensions (Technical name: level 2-ProForma + top-down compliant)

   - [ ] 4.2.1 Additional CV/ontologies for protein modifications: RESID (the prefix R MUST be used for RESID CV/ontology term names)
   - [x] 4.2.8 Chemical formulas (this feature occurs in two places in this list).

1. Cross-Linking Extensions (Technical name: level 2-ProForma + cross-linking compliant)

   - [ ] 4.2.1/4.2.3 Cross-linked peptides (using the XL-MOD CV/ontology, the prefix X MUST be used for XL-MOD CV/ontology term names).

1. Glycan Extensions (Technical name: level 2-ProForma + glycans compliant)

   - [x] 4.2.5 Additional CV/ontologies for protein modifications: GNO (the prefix G MUST be used for GNO CV/ontology term names)
   - [x] 4.2.9 Glycan composition.
   - [x] 4.2.8 Chemical formulas (this feature occurs in two places in this list).

1. Spectral Support (Technical name: level 2-ProForma + mass spectrum compliant)

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
- [ ] 4.2.4 Branched peptides, similar to cross linking
- [ ] Modification specific neutral losses, phospho: fully deleted
- [ ] Modification diagnostic peaks
- [ ] Match isotope patterns in fragmentation matching
- [ ] Custom modifications defined using a CV grammar but loaded at runtime

## Bugs

- [x] Allow isotope definition in formulas
- [ ] Allow for non-named localised ambiguous modifications, using mass or formulas `A[Formula:H-1#g1]AA[#g1]`

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

| Name    | Modifications | Numbered | Rules |
| ------- | ------------- | -------- | ----- |
| Unimod  | Yes           | Yes      | Yes   |
| PSI-MOD | Yes           | Yes      | Yes   |
| RESID   | -             | -        | -     |
| XL-MOD  | -             | -        | -     |
| GNO     | Yes           | NA       | NA    |
