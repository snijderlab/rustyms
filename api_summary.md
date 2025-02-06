# Data types

## Peptidoforms
3 types of peptides defined in mzcore (in order of increasing complexity):

### Peptidoform (former LinearPeptide)
**Description:** A peptidoform with all data as specified by [ProForma](https://github.com/HUPO-PSI/ProForma) (= peptide sequence + modifications).

#### Functionality
*

### PeptidoformIon
**Description:** Peptidoform + charge.

#### Functionality

### CompoundPeptidoformIon
**Description:** A grouping of multiple peptidoform ions that cannot be distinguished from each other by their mass spectrum.

#### Functionality

# Proposed actions

* Rename `LinearPeptide` to `Peptidoform` in order to comply with the ProForma definition and remove any potential cause of ambiguations.


