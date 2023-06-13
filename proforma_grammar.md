Here is a full grammar of the pro forma syntax. `\` is used to escape characters. `<name>` is used to refer to another rule. `()` are used to group rules. `|` is used to give alternate options. `*+?` have their common regex meaning. `()x*` indicates a group that is repeated and separated by `x`, any repeating operator can be used.

```
<full_definition> = (<full_sequence>)++
<full_sequence> = <sequence>(//<sequence>)?(\\\\<sequence>)?
<sequence> = <global>*<labile>*<unknown_position_modification>*<n-term>?<aa_sequence><c-term>?(/<charge>[(<charge><atom>(+|-)),*])?
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
missing:
proper order for global/labile/unknown_position?

# Levels of support

1) Base Level Support (Technical name: Base-ProForma Compliant)
Represents the lowest level of compliance, this level involves providing support for:
-[x] Amino acid sequences
-[ ] Protein modifications using two of the supported CVs/ontologies: Unimod and PSI-MOD.
-[x] Protein modifications using delta masses (without prefixes)
-[x] N-terminal, C-terminal and labile modifications.
-[ ] Ambiguity in the modification position, including support for localisation scores.
-[ ] Ambiguity in the amino acid sequence. ??
-[ ] INFO tag.
2) Additional Separate Support (Technical name: level 2-ProForma compliant)
These features are independent from each other:
-[x] Unusual amino acids (O and U).
-[x] Ambiguous amino acids (e.g. X, B, Z). This would include support for sequence tags of known mass (using the character X).
-[ ] Protein modifications using delta masses (using prefixes for the different CVs/ontologies).
-[ ] Use of prefixes for Unimod (U:) and PSI-MOD (M:) names.
-[ ] Support for the joint representation of experimental data and its interpretation.
3) Top-Down Extensions (Technical name: level 2-ProForma + top-down compliant)
-[ ] Additional CV/ontologies for protein modifications: RESID (the prefix R MUST be used for RESID CV/ontology term names)
-[x] Chemical formulas (this feature occurs in two places in this list).
4) Cross-Linking Extensions (Technical name: level 2-ProForma + cross-linking
compliant)
-[ ] Cross-linked peptides (using the XL-MOD CV/ontology, the prefix X MUST be used for XL-MOD CV/ontology term names).
5) Glycan Extensions (Technical name: level 2-ProForma + glycans compliant)
-[ ] Additional CV/ontologies for protein modifications: GNO (the prefix G MUST be used for GNO CV/ontology term names)
-[x] Glycan composition.
-[x] Chemical formulas (this feature occurs in two places in this list).
6) Spectral Support (Technical name: level 2-ProForma + mass spectrum compliant)
-[ ] Charge and chimeric spectra are special cases (see Appendix II).
-[ ] Global modifications (e.g., every C is C13).