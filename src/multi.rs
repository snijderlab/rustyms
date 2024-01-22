use std::ops::{Add, Deref, Mul, MulAssign, Neg, Sub};

use itertools::{Itertools, MinMaxResult};
use serde::{Deserialize, Serialize};

use crate::system::OrderedMass;

/// A collection of potentially multiple of the generic type, it is used be able to easily
/// combine multiple of this multi struct into all possible combinations.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Hash)]
pub struct Multi<M>(Vec<M>);

impl<M> Multi<M> {
    /// Get the underlying vector
    pub fn as_vec(self) -> Vec<M> {
        self.0
    }
}

impl<M: Eq + std::hash::Hash + Clone> Multi<M> {
    /// Get all unique values
    #[must_use]
    pub fn unique(&self) -> Self {
        self.0.iter().unique().cloned().collect()
    }
}

impl<'a, M> Multi<M>
where
    &'a M: 'a,
    OrderedMass: From<&'a M>,
{
    /// Get the bounds for the mass
    pub fn mass_bounds(&'a self) -> MinMaxResult<&'a M> {
        self.0.iter().minmax_by_key(|f| OrderedMass::from(f))
    }
}

impl<M> Deref for Multi<M> {
    type Target = Vec<M>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<M: Default> Default for Multi<M> {
    // Default is one empty M to make the cartesian product with a default return useful results
    fn default() -> Self {
        Self(vec![M::default()])
    }
}

impl<'a, M> Neg for &'a Multi<M>
where
    &'a M: Neg<Output = M> + 'a,
{
    type Output = Multi<M>;
    fn neg(self) -> Self::Output {
        self.0.iter().map(|f| -f).collect()
    }
}

impl<M> Neg for Multi<M>
where
    M: Neg<Output = M>,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        self.0.into_iter().map(|f| -f).collect()
    }
}

impl<'a, M> Add<&'a M> for &'a Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Multi<M>;
    /// Adds this M to all Ms in the multi
    fn add(self, rhs: &M) -> Self::Output {
        Multi(self.0.iter().map(|m| rhs.clone() + m).collect())
    }
}

impl<'a, M> Add<&'a M> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
{
    type Output = Self;
    /// Adds this formula to all formulas in the multi formula
    fn add(self, rhs: &M) -> Self::Output {
        self.0.into_iter().map(|m| m + rhs.clone()).collect()
    }
}

impl<'a, M> Add<M> for &'a Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Multi<M>;
    /// Adds this formula to all formulas in the multi formula
    fn add(self, rhs: M) -> Self::Output {
        Multi(self.0.iter().map(|m| rhs.clone() + m).collect())
    }
}

impl<M> Add<M> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
{
    type Output = Self;
    /// Adds this formula to all formulas in the multi formula
    fn add(self, rhs: M) -> Self::Output {
        self.0.into_iter().map(|m| m + rhs.clone()).collect()
    }
}

impl<'a, M> Sub<&'a M> for &'a Multi<M>
where
    &'a M: Sub<M, Output = M> + 'a,
    M: Clone,
{
    type Output = Multi<M>;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: &M) -> Self::Output {
        Multi(self.0.iter().map(|m| m - rhs.clone()).collect())
    }
}

impl<'a, M> Sub<&'a M> for Multi<M>
where
    M: Sub<M, Output = M> + Clone,
{
    type Output = Self;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: &M) -> Self::Output {
        self.0.into_iter().map(|m| m - rhs.clone()).collect()
    }
}

impl<'a, M> Sub<M> for Multi<&'a M>
where
    &'a M: Sub<M, Output = M> + 'a,
    M: Clone,
{
    type Output = Multi<M>;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: M) -> Self::Output {
        self.0.into_iter().map(|m| m - rhs.clone()).collect()
    }
}

impl<M> Sub<M> for Multi<M>
where
    M: Sub<M, Output = M> + Clone,
{
    type Output = Self;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: M) -> Self::Output {
        self.0.into_iter().map(|m| m - rhs.clone()).collect()
    }
}

impl<'a, M> Mul<&'a Multi<M>> for &'a Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Multi<M>;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: &Multi<M>) -> Self::Output {
        Multi(
            self.0
                .iter()
                .cartesian_product(rhs.0.iter())
                .map(|(a, b)| b.clone() + a)
                .collect(),
        )
    }
}

impl<'a, M> Mul<&'a Self> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
{
    type Output = Self;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: &Self) -> Self::Output {
        self.0
            .into_iter()
            .cartesian_product(rhs.0.iter())
            .map(|(a, b)| b.clone() + a)
            .collect()
    }
}

impl<'a, M> Mul<Multi<M>> for &'a Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Multi<M>;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: Multi<M>) -> Self::Output {
        self.0
            .iter()
            .cartesian_product(rhs.0.iter())
            .map(|(a, b)| b.clone() + a)
            .collect()
    }
}

impl<M> Mul<Self> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
{
    type Output = Self;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: Self) -> Self::Output {
        self.0
            .into_iter()
            .cartesian_product(rhs.0.iter())
            .map(|(a, b)| b.clone() + a)
            .collect()
    }
}

impl<'a, M> MulAssign<&'a Self> for Multi<M>
where
    M: Add<M, Output = M> + 'a + Clone,
{
    fn mul_assign(&mut self, rhs: &Self) {
        let new_data = self
            .0
            .iter()
            .cartesian_product(rhs.0.iter())
            .map(|(a, b)| b.clone() + a.clone())
            .collect();
        self.0 = new_data;
    }
}

impl<'a, M> MulAssign<Self> for Multi<M>
where
    M: Add<M, Output = M> + 'a + Clone,
{
    fn mul_assign(&mut self, rhs: Self) {
        let new_data = self
            .0
            .iter()
            .cartesian_product(rhs.0.iter())
            .map(|(a, b)| b.clone() + a.clone())
            .collect();
        self.0 = new_data;
    }
}

impl<M: Default> std::iter::Sum<Self> for Multi<M>
where
    Self: MulAssign<Self>,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut res = Self::default();
        iter.for_each(|v| res *= v);
        res
    }
}

impl<M> From<M> for Multi<M> {
    fn from(value: M) -> Self {
        Self(vec![value])
    }
}

impl<M: Clone> From<&M> for Multi<M> {
    fn from(value: &M) -> Self {
        Self(vec![value.clone()])
    }
}

impl<M> From<Vec<M>> for Multi<M> {
    fn from(value: Vec<M>) -> Self {
        Self(value)
    }
}

impl<M: Clone> From<&[M]> for Multi<M> {
    fn from(value: &[M]) -> Self {
        Self(value.to_vec())
    }
}

impl<M> std::iter::FromIterator<M> for Multi<M> {
    fn from_iter<T: IntoIterator<Item = M>>(iter: T) -> Self {
        Self(iter.into_iter().collect_vec())
    }
}

impl<'a, M: Clone + 'a> std::iter::FromIterator<&'a M> for Multi<M> {
    fn from_iter<T: IntoIterator<Item = &'a M>>(iter: T) -> Self {
        Self(iter.into_iter().cloned().collect_vec())
    }
}
