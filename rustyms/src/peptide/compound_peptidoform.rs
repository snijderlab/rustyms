use std::fmt::{Display, Write};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    peptide::Linked, system::usize::Charge, Fragment, LinearPeptide, Model, MolecularFormula,
    Multi, MultiChemical, Peptidoform, SequenceElement,
};

/// A single pro forma entry, can contain multiple peptidoforms
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct CompoundPeptidoform(pub(super) Vec<Peptidoform>);

impl MultiChemical for CompoundPeptidoform {
    /// Gives all possible formulas for this compound peptidoform
    fn formulas(&self, sequence_index: usize, peptide_index: usize) -> Multi<MolecularFormula> {
        self.0
            .iter()
            .flat_map(|p| p.formulas(sequence_index, peptide_index).to_vec())
            .collect()
    }
}

impl CompoundPeptidoform {
    /// Assume there is exactly one peptidoform in this collection.
    #[doc(alias = "assume_linear")]
    pub fn singular(mut self) -> Option<Peptidoform> {
        if self.0.len() == 1 {
            self.0.pop()
        } else {
            None
        }
    }

    /// Assume there is exactly one peptide in this collection.
    pub fn singular_peptide(self) -> Option<LinearPeptide<Linked>> {
        self.singular().and_then(Peptidoform::singular)
    }

    /// Get all peptidoforms making up this `CompoundPeptidoform`
    pub fn peptidoforms(&self) -> &[Peptidoform] {
        &self.0
    }

    /// Generate the theoretical fragments for this peptide collection.
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
    ) -> Vec<Fragment> {
        let mut base = Vec::new();
        for (index, peptidoform) in self.peptidoforms().iter().enumerate() {
            base.extend(peptidoform.generate_theoretical_fragments(max_charge, model, index));
        }
        base
    }

    /// Display this compound peptidoform.
    /// `specification_compliant` Displays this compound peptidoform either normalised to the internal representation or as fully spec compliant ProForma
    /// (no glycan structure or custom modifications).
    /// # Panics
    /// When some peptides do not have the same global isotope modifications.
    /// # Errors
    /// If the underlying formatter errors.
    pub fn display(&self, f: &mut impl Write, specification_compliant: bool) -> std::fmt::Result {
        let global_equal = self
            .peptidoforms()
            .iter()
            .flat_map(Peptidoform::peptides)
            .map(|p| &p.global)
            .tuple_windows()
            .all(|(a, b)| a == b);
        assert!(global_equal, "Not all global isotope modifications on all peptides on this compound peptidoform are identical");
        let empty = Vec::new();
        let global = self
            .peptidoforms()
            .iter()
            .flat_map(Peptidoform::peptides)
            .next()
            .map_or(&empty, |p| &p.global);
        for (element, isotope) in global {
            write!(
                f,
                "<{}{}>",
                isotope.map(|i| i.to_string()).unwrap_or_default(),
                element
            )?;
        }

        let mut first = true;
        for p in self.peptidoforms() {
            if !first {
                write!(f, "+")?;
            }
            p.display(f, false, specification_compliant)?;
            first = false;
        }
        Ok(())
    }
}

impl Display for CompoundPeptidoform {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(f, true)
    }
}

impl<Item> From<Item> for CompoundPeptidoform
where
    Item: Into<LinearPeptide<Linked>>,
{
    fn from(value: Item) -> Self {
        Self(vec![Peptidoform(vec![value.into()])])
    }
}

impl<Item> FromIterator<Item> for CompoundPeptidoform
where
    Item: Into<SequenceElement>,
{
    fn from_iter<T: IntoIterator<Item = Item>>(iter: T) -> Self {
        Self(vec![Peptidoform(vec![LinearPeptide::from(iter)])])
    }
}
