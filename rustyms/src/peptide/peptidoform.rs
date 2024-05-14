use serde::{Deserialize, Serialize};

use crate::{
    peptide::Linked, system::usize::Charge, Fragment, LinearPeptide, Model, MolecularFormula,
    Multi, MultiChemical,
};
/// A single peptidoform, can contain multiple linear peptides
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct Peptidoform(pub(crate) Vec<LinearPeptide<Linked>>);

impl MultiChemical for Peptidoform {
    /// Gives all possible formulas for this peptidoform.
    /// Assumes all peptides in this peptidoform are connected.
    /// If there are no peptides in this peptidoform it returns [`Multi::default`].
    fn formulas(&self) -> Multi<MolecularFormula> {
        self.0
            .first()
            .map(|p| p.formulas(0, &self.0, &[], &mut Vec::new()))
            .unwrap_or_default()
    }
}

impl Peptidoform {
    /// Generate the theoretical fragments for this peptide collection.
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
        peptidoform_index: usize,
    ) -> Vec<Fragment> {
        let mut base = Vec::new();
        for (index, peptide) in self.peptides().iter().enumerate() {
            base.extend(peptide.generate_theoretical_fragments_inner(
                max_charge,
                model,
                peptidoform_index,
                index,
                &self.0,
            ));
        }
        base
    }

    /// Assume there is exactly one peptide in this collection.
    pub fn singular(mut self) -> Option<LinearPeptide<Linked>> {
        if self.0.len() == 1 {
            self.0.pop()
        } else {
            None
        }
    }

    /// Get all peptides making up this `Peptidoform`
    pub fn peptides(&self) -> &[LinearPeptide<Linked>] {
        &self.0
    }
}

impl std::fmt::Display for Peptidoform {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(p) = self.0.first() {
            write!(f, "{p}")?;
        }
        for p in self.peptides().iter().skip(1) {
            write!(f, "//{p}")?;
        }
        Ok(())
    }
}
