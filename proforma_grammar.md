Here is a full grammar of the pro forma syntax. `\` is used to escape characters. `<name>` is used to refer to another rule. `()` are used to group rules. `|` is used to give alternate options. `*+?` have their common regex meaning.

```
<full_sequence> = <sequence>(//<sequence>)?(\\\\<sequence>)?
<sequence> = <global>*<labile>*<unknown_position_modification>*<n-term>?<aa_sequence><c-term>?
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