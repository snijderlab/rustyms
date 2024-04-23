use serde::{Deserialize, Serialize};

use crate::{
    peptide_complexity::Linked, system::usize::Charge, Fragment, LinearPeptide, Model,
    MolecularFormula, Multi, MultiChemical,
};
/// A single peptidoform, can contain multiple linear peptides
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct Peptidoform(pub(crate) Vec<LinearPeptide<Linked>>);

impl MultiChemical for Peptidoform {
    /// Gives all possible formulas for this peptidoform
    fn formulas(&self) -> Multi<MolecularFormula> {
        self.0
            .iter()
            .enumerate()
            .fold(Multi::default(), |acc, (i, p)| {
                acc * p.formulas(i, &self.0, &[])
            })
    }
}

impl Peptidoform {
    /// Generate the theoretical fragments for this peptide collection.
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
    ) -> Vec<Fragment> {
        let mut base = Vec::new();
        for (index, peptide) in self.peptides().iter().enumerate() {
            base.extend(peptide.generate_theoretical_fragments_inner(
                max_charge,
                model,
                index,
                &self.0,
                &[],
            ));
        }
        base
    }

    /// Assume there is exactly one peptide in this collection.
    pub fn singular(mut self) -> Option<LinearPeptide<Linked>> {
        if self.0.len() == 1 {
            Some(self.0.pop().unwrap())
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
