use std::fmt::Display;

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
    fn formulas(&self) -> Multi<MolecularFormula> {
        self.0.iter().flat_map(|p| p.formulas().to_vec()).collect()
    }
}

impl Display for CompoundPeptidoform {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // TODO: aggregate global, handle inconsistent global cases?
        if let Some(p) = self.0.first() {
            write!(f, "{p}")?;
        }
        for p in self.peptidoforms().iter().skip(1) {
            write!(f, "+{p}")?;
        }
        Ok(())
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
