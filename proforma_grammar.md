# ProForma support

## Grammar

Here is a full grammar of the pro forma syntax. `\` is used to escape characters. `<name>` is used to refer to another rule. `()` are used to group rules. `|` is used to give alternate options. `*+?` have their common regex meaning. `()x*` indicates a group that is repeated and separated by `x`, any repeating operator can be used.

```ebnf
full_definition = (full_sequence)++
full_sequence = sequence ("//" sequence)?
sequence = <global>*<labile>*<unknown_position_modification>*<n-term>?<aa_sequence><c-term>?(/<charge>([(<charge><atom>(+|-)),*])?)?
aa_sequence = (<aa_sequence_inner>|<range>)*
range = \(<aa_sequence_inner>\)<modification>
aa_sequence_inner = (<aa_sequence_inner_inner>)+|\(\?<aa_sequence_inner_inner>\)
aa_sequence_inner_inner = <aminoacid><modification>*
global = \<<global_isotope!!>|<modification>@<aminoacids>\>
global_isotope = <num><atom>
aminoacids = <aminoacid>(,<aminoacid>)*
unknown_position_modification = <modification>(^<num>)?\?
modification = [<modification_inner>]
labile = {<modification_inner>}(^<num>)?
n_term = [<modification_inner>]-
c_term = -[<modification_inner>]
modification_inner = (<modification_inner_inner>|<cross_link>|<branch>)(\|<modification_inner_inner>)*
modification_inner_inner = ((<named_modification>|<numbered_modification>|<delta_mass>|<formula>|<glycan>)(#<name>(\(<decimal>\))?)?)|(\|INFO:<info>)
named_modification = ((U:)?|(M:)?)<name>|(R:|X:|G:)<name>
numbered_modification = (UNIMOD|MOD|RESID|XLMOD|GNO):<num> ?? should be name?
cross_link = <modification_inner_inner>?#XL<name>
branch = <modification_inner_inner>?#BRANCH
delta_mass = ((U|M|R|X|G|Obs):)?<shift>
formula = Formula:<pair>+
pair = <atom><cardinality>?|[<num><atom><cardinality>?]
glycan = Glycan:(<monosaccharide><num>?)+
charge =regex= (+|-)?[0-9]+
cardinality =regex= (-)?[0-9]+
shift =regex= (+|-)[0-9]+\.[0-9]+
name =regex= [A-Za-z0-9]+
num =regex= [0-9]+
decimal =regex= [0-9]+\.[0-9]+
atom =human= any valid atom name
aminoacid =human= any IUPAC amino acid character or any UniProt ambiguous/unusual
info =human= any piece of text, with the exclusion of unpaired brackets
```

Here is the grammar https://github.com/bittremieux/spectrum_utils/blob/main/spectrum_utils/proforma.ebnf rewritten to pass this validator: https://thomasgassmann.com/ebnf.
```ebnf
# %import common.DIGIT
# %import common.INT
# %import common.LETTER
# %import common.NUMBER
# %import common.SIGNED_INT
# %import common.SIGNED_NUMBER
# %import common.WS
# %import .monosaccharide.MONOSACCHARIDE

# ProForma specification: https://github.com/HUPO-PSI/ProForma/
# Version: June 29, 2021

<DIGIT> <= 0|1|2|3|4|5|6|7|8|9
<INT> <= <DIGIT>+
<SIGNED_INT> <= [+|-]?<DIGIT>+
<NUMBER> <= <DIGIT>+("."<DIGIT>+)?
<LETTER> <= A|B|C|D|E|F|G|H|I|J|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z|a|b|c|d|e|f|g|h|i|j|k|l|m|n|o|p|q|r|s|t|u|v|w|x|y|z
<WS> <= (" "|"   ")*

<CROSSLINK> <= "//"
<CHIMERIC> <= "+"

<proteoform> <= <peptide> ["/" <charge>]

<peptide> <= <mod_global>* <mod_unknown_pos>? <mod_labile>* <mod_n_term>? (<aa> | <mod_range>)+ <mod_c_term>?
# TODO: Amino acid sequence ambiguity (section 4.7).

<aa> <= <AA> [<mod>+ | (<_MOD_L> <mod_label> <_MOD_R>)]
<AA> <= <LETTER>

<mod_global> <= <_MOD_GLOBAL_L> (<ISOTOPE> | (<mod> "@" (<AA> ",")* <AA>)) <_MOD_GLOBAL_R>
<ISOTOPE> <= <INT>? <LETTER>+ <SIGNED_INT>?

<mod_unknown_pos> <= (<mod> ["^" <MOD_COUNT>])+ "?"

<mod> <=        <_MOD_L>        ((<mod_name> | <mod_accession> | <mod_mass> | <mod_formula> | <mod_glycan> | <info>) <mod_label>? "|")* (<mod_name> | <mod_accession> | <mod_mass> | <mod_formula> | <mod_glycan> | <info>) <mod_label>? <_MOD_R>
<mod_labile> <= <_MOD_LABILE_L> ((<mod_name> | <mod_accession> | <mod_mass> | <mod_formula> | <mod_glycan> | <info>)              "|")* (<mod_name> | <mod_accession> | <mod_mass> | <mod_formula> | <mod_glycan> | <info>)              <_MOD_LABILE_R>
<MOD_COUNT> <= <INT>

<mod_n_term> <= (<mod> | (<_MOD_L> <mod_label> <_MOD_R>)) "-"
<mod_c_term> <= "-" (<mod> | (<_MOD_L> <mod_label> <_MOD_R>))

<mod_range> <= <MOD_RANGE_L> <mod_range_pos> <_MOD_RANGE_R> <mod>+
<mod_range_pos> <= (<aa> | <mod_range>)+

<mod_name> <= ((<CV_ABBREV> ":") | (<CV_ABBREV_OPT> ":")?) <TEXT>
<i> <= ":"
<CV_ABBREV_OPT> <= "U"<i> | "M"<i>
<CV_ABBREV> <= "R"<i> | "X"<i> | "G"<i>

<mod_accession> <= <CV_NAME> ":" <TEXT>
<CV_NAME> <= "UNIMOD"<i> | "MOD"<i> | "RESID"<i> | "XLMOD"<i> | "GNO"<i>

<mod_mass> <= [(<CV_ABBREV_OPT> | <CV_ABBREV> | <MOD_MASS_OBS>) ":"] <MOD_MASS>
<MOD_MASS_OBS> <= "Obs"<i>
<MOD_MASS> <= ("+" | "-") <NUMBER>

<mod_formula> <= "Formula:"<i> (<_MOD_L> <ISOTOPE> <_MOD_R>)* <FORMULA>
<FORMULA> <= (<LETTER>+ <SIGNED_INT>? <WS>?)+

<mod_glycan> <= "Glycan" (<monosaccharide> <WS>?)+
<monosaccharide> <= <MONOSACCHARIDE> <MONOSACCHARIDE_COUNT>?
<MONOSACCHARIDE_COUNT> <= <INT>

<info> <= "Info"<i> <TEXT>

<mod_label> <= "#" (<MOD_LABEL_XL> | <MOD_LABEL_BRANCH> | <MOD_LABEL>) ["(" <MOD_SCORE> ")"]
<MOD_LABEL_XL> <= "XL" <MOD_LABEL>
<MOD_LABEL_BRANCH> <= "BRANCH"
<MOD_LABEL> <= (<LETTER> | <DIGIT>)+
<MOD_SCORE> <= <SIGNED_NUMBER>

<charge> <= <CHARGE> [<_MOD_L> <ion> <_MOD_R>]
<CHARGE> <= <SIGNED_INT>
<ion> <= [<TEXT> ","] <TEXT>

<TEXT> <= [<LETTER>|" "|<DIGIT>]+

<_MOD_L> <= "["
<_MOD_R> <= "]"
<_MOD_LABILE_L> <= "{"
<_MOD_LABILE_R> <= "}"
<_MOD_GLOBAL_L> <= "<"
<_MOD_GLOBAL_R> <= ">"
<MOD_RANGE_L> <= "("
<_MOD_RANGE_R> <= ")"

<proforma> <= (<proteoform> (<CROSSLINK> | <CHIMERIC>))* <proteoform>
```

## Open questions

- proper order for global/labile/unknown_position and for additional peptides (xl/cys-xl/branched/chimeric) for MS?
  - A: `<GLOBAL_MOD>[UNKNOWN_POS]?{LABILE_MOD}[N_TERM]-PEPTIDE-[C_TERM]` https://github.com/HUPO-PSI/ProForma/issues/3#issuecomment-906448694
- what defines valid ionic species for ion charge?
- for ionic charge/adduct ions, is there a way to specify higher charged ionic species? (Ca+2 or Ca2+?)
- separation of multiple linear peptides happens with `+` (chimeric) `//` (xl) and `\\` (branch)

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

   - [ ] 4.2.1/4.2.3 Cross-linked peptides (using the XL-MOD CV/ontology, the prefix X MUST be used for XL-MOD CV/ontology term names).

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
- [ ] 4.2.4 Branched peptides, similar to cross linking

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

| Name    | Modifications | Numbered | Rules | Diagnostic ions / neutral losses |
| ------- | ------------- | -------- | ----- | -------------------------------- |
| Unimod  | Yes           | Yes      | Yes   | Yes (neutral losses)             |
| PSI-MOD | Yes           | Yes      | Yes   | NA                               |
| RESID   | -             | -        | -     | -                                |
| XL-MOD  | -             | -        | -     | -                                |
| GNO     | Yes           | NA       | NA    | NA (solved for all glycans)      |
