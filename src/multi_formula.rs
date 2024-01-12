use std::ops::{Add, Deref, Mul, MulAssign, Neg, Sub};

use itertools::{Itertools, MinMaxResult};
use serde::{Deserialize, Serialize};

use crate::{system::OrderedMass, MolecularFormula};

/// Any item that has a number of potential chemical formulas
pub trait MultiChemical {
    /// Get all possible molecular formulas
    fn formulas(&self) -> MultiMolecularFormula;
}
/// A set of different molecular formulas resulting from the same entity, if the entity has multiple possible masses, or if the entity is a collection of multiple things.
/// For convenience [`std::ops::Mul`] has been implemented as the cartesian product between two [`MultiMolecularFormula`]s. The [`std::iter::Sum`] implementation does the same.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Hash)]
pub struct MultiMolecularFormula(Vec<MolecularFormula>);

impl MultiMolecularFormula {
    /// Get all possible molecular formulas filtered to only return unique formulas
    #[must_use]
    pub fn unique_formulas(&self) -> Self {
        self.0.iter().unique().cloned().collect()
    }

    /// Get the bounds for the mass.
    pub fn mass_bounds(&self) -> MinMaxResult<&MolecularFormula> {
        self.0
            .iter()
            .minmax_by_key(|f| OrderedMass::from(f.monoisotopic_mass()))
    }

    /// Get the underlying vector
    pub fn as_vec(self) -> Vec<MolecularFormula> {
        self.0
    }
}

impl Deref for MultiMolecularFormula {
    type Target = Vec<MolecularFormula>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Default for MultiMolecularFormula {
    // Default is one empty formula to make the cartesian product with a default return useful results
    fn default() -> Self {
        Self(vec![MolecularFormula::default()])
    }
}

impl Neg for &MultiMolecularFormula {
    type Output = MultiMolecularFormula;
    fn neg(self) -> Self::Output {
        self.0.iter().map(|f| -f).collect()
    }
}

impl Neg for MultiMolecularFormula {
    type Output = Self;
    fn neg(self) -> Self::Output {
        self.0.into_iter().map(|f| -f).collect()
    }
}

impl Add<&MolecularFormula> for &MultiMolecularFormula {
    type Output = MultiMolecularFormula;
    /// Adds this formula to all formulas in the multi formula
    fn add(self, rhs: &MolecularFormula) -> Self::Output {
        MultiMolecularFormula(self.iter().map(|m| m + rhs).collect())
    }
}

impl Sub<&MolecularFormula> for &MultiMolecularFormula {
    type Output = MultiMolecularFormula;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: &MolecularFormula) -> Self::Output {
        MultiMolecularFormula(self.iter().map(|m| m - rhs).collect())
    }
}

impl Mul<&MultiMolecularFormula> for &MultiMolecularFormula {
    type Output = MultiMolecularFormula;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: &MultiMolecularFormula) -> Self::Output {
        MultiMolecularFormula(
            self.iter()
                .cartesian_product(rhs.iter())
                .map(|(a, b)| a + b)
                .collect(),
        )
    }
}

impl_binop_ref_cases!(impl Add, add for MultiMolecularFormula, MolecularFormula, MultiMolecularFormula);
impl_binop_ref_cases!(impl Sub, sub for MultiMolecularFormula, MolecularFormula, MultiMolecularFormula);
impl_binop_ref_cases!(impl Mul, mul for MultiMolecularFormula, MultiMolecularFormula, MultiMolecularFormula);

impl MulAssign<&Self> for MultiMolecularFormula {
    fn mul_assign(&mut self, rhs: &Self) {
        self.0 = self
            .0
            .iter()
            .cartesian_product(rhs.iter())
            .map(|(a, b)| a + b)
            .collect();
    }
}

impl MulAssign<Self> for MultiMolecularFormula {
    fn mul_assign(&mut self, rhs: Self) {
        *self *= &rhs;
    }
}

impl std::iter::Sum<Self> for MultiMolecularFormula {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut res = Self::default();
        iter.for_each(|v| res *= v);
        res
    }
}

impl From<MolecularFormula> for MultiMolecularFormula {
    fn from(value: MolecularFormula) -> Self {
        Self(vec![value])
    }
}

impl From<&MolecularFormula> for MultiMolecularFormula {
    fn from(value: &MolecularFormula) -> Self {
        Self(vec![value.clone()])
    }
}

impl From<Vec<MolecularFormula>> for MultiMolecularFormula {
    fn from(value: Vec<MolecularFormula>) -> Self {
        Self(value)
    }
}

impl From<&[MolecularFormula]> for MultiMolecularFormula {
    fn from(value: &[MolecularFormula]) -> Self {
        Self(value.to_vec())
    }
}

impl std::iter::FromIterator<MolecularFormula> for MultiMolecularFormula {
    fn from_iter<T: IntoIterator<Item = MolecularFormula>>(iter: T) -> Self {
        Self(iter.into_iter().collect_vec())
    }
}

impl<'a> std::iter::FromIterator<&'a MolecularFormula> for MultiMolecularFormula {
    fn from_iter<T: IntoIterator<Item = &'a MolecularFormula>>(iter: T) -> Self {
        Self(iter.into_iter().cloned().collect_vec())
    }
}
