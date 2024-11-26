use std::fmt::{Display, Write};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    peptide::Linked, system::usize::Charge, Fragment, LinearPeptide, Model, MolecularFormula,
    Multi, Peptidoform,
};

/// A single full ProForma entry. This entry can contain multiple sets of cross-linked peptides.
/// A single set of cross-linked peptides is a [`Peptidoform`]. A ProForma entry with two chimeric
/// peptides will be saved as one [`CompoundPeptidoform`] with two [`Peptidoform`]s that each
/// contain one of the [`LinearPeptide`]s.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct CompoundPeptidoform(pub(super) Vec<Peptidoform>);

impl CompoundPeptidoform {
    /// Create a new [`CompoundPeptidoform`] from many [`Peptidoform`]s. This returns None if the
    /// global isotope modifications of all peptidoforms are not identical.
    pub fn new(iter: impl IntoIterator<Item = Peptidoform>) -> Option<Self> {
        let result = Self(iter.into_iter().collect());
        let global_equal = result
            .peptidoforms()
            .iter()
            .flat_map(Peptidoform::peptides)
            .tuple_windows()
            .all(|(a, b)| a.get_global() == b.get_global());
        global_equal.then_some(result)
    }

    /// Get all possible formulas for this compound peptidoform
    pub fn formulas(&self) -> Multi<MolecularFormula> {
        self.0.iter().flat_map(|p| p.formulas().to_vec()).collect()
    }

    /// Assume there is exactly one peptidoform in this compound peptidoform.
    #[doc(alias = "assume_linear")]
    pub fn singular(mut self) -> Option<Peptidoform> {
        if self.0.len() == 1 {
            self.0.pop()
        } else {
            None
        }
    }

    /// Assume there is exactly one peptide in this compound peptidoform.
    pub fn singular_peptide(self) -> Option<LinearPeptide<Linked>> {
        self.singular().and_then(Peptidoform::singular)
    }

    /// Get all peptidoforms making up this compound peptidoform.
    pub fn peptidoforms(&self) -> &[Peptidoform] {
        &self.0
    }

    /// Get all peptides making up this compound peptidoform.
    pub fn peptides(&self) -> impl Iterator<Item = &LinearPeptide<Linked>> {
        self.0.iter().flat_map(Peptidoform::peptides)
    }

    /// Generate the theoretical fragments for this compound peptidoform.
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
    ) -> Vec<Fragment> {
        let mut base = Vec::new();
        for (index, peptidoform) in self.peptidoforms().iter().enumerate() {
            base.extend(peptidoform.generate_theoretical_fragments_inner(max_charge, model, index));
        }
        base
    }

    /// Display this compound peptidoform.
    /// `specification_compliant` Displays this compound peptidoform either normalised to the
    /// internal representation (with false) or as fully spec compliant ProForma (no glycan
    /// structure or custom modifications) (with true).
    /// # Errors
    /// Only if the underlying formatter (`f`) errors.
    pub fn display(&self, f: &mut impl Write, specification_compliant: bool) -> std::fmt::Result {
        // The global isotope modifications are guaranteed to be identical, so take the first
        let empty = Vec::new();
        let global = self
            .peptidoforms()
            .iter()
            .flat_map(Peptidoform::peptides)
            .next()
            .map_or(empty.as_slice(), |p| p.get_global());
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

impl<Complexity> From<LinearPeptide<Complexity>> for CompoundPeptidoform {
    fn from(value: LinearPeptide<Complexity>) -> Self {
        Self(vec![Peptidoform(vec![value.mark()])])
    }
}

impl From<Peptidoform> for CompoundPeptidoform {
    fn from(value: Peptidoform) -> Self {
        Self(vec![value])
    }
}

impl<Complexity> From<Vec<LinearPeptide<Complexity>>> for CompoundPeptidoform {
    fn from(value: Vec<LinearPeptide<Complexity>>) -> Self {
        Self(value.into_iter().map(Into::into).collect())
    }
}
