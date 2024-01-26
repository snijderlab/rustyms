use std::{
    ops::{Add, Deref, Mul, MulAssign, Neg, Sub},
    rc::Rc,
};

use itertools::{Itertools, MinMaxResult};
use serde::{Deserialize, Serialize};

use crate::system::OrderedMass;

/// A collection of potentially multiple of the generic type, it is used be able to easily
/// combine multiple of this multi struct into all possible combinations.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Hash)]
pub enum Multi<M> {
    Many(Rc<[M]>),
    Single(M),
}

impl<M: Eq + std::hash::Hash + Clone> Multi<M> {
    /// Get all unique values
    #[must_use]
    pub fn unique(&self) -> Self {
        match self {
            Self::Many(m) => m.iter().unique().cloned().collect(),
            Self::Single(..) => self.clone(),
        }
    }
}

impl<'a, M> Multi<M>
where
    &'a M: 'a,
    OrderedMass: From<&'a M>,
{
    /// Get the bounds for the mass
    pub fn mass_bounds(&'a self) -> MinMaxResult<&'a M> {
        match self {
            Self::Many(m) => m.iter().minmax_by_key(|f| OrderedMass::from(f)),
            Self::Single(m) => MinMaxResult::OneElement(m),
        }
    }
}

// impl<M> IntoIterator for Multi<M> {
//     type IntoIter = ();
//     type Item = &M;
//     fn into_iter() {}
// }

impl<M> Deref for Multi<M> {
    type Target = [M];
    fn deref(&self) -> &Self::Target {
        match self {
            Self::Many(m) => m,
            Self::Single(m) => std::slice::from_ref(m),
        }
    }
}

impl<M: Default> Default for Multi<M> {
    // Default is one empty M to make the cartesian product with a default return useful results
    fn default() -> Self {
        Self::Single(M::default())
    }
}

impl<'a, M> Neg for &'a Multi<M>
where
    &'a M: Neg<Output = M> + 'a,
{
    type Output = Multi<M>;
    fn neg(self) -> Self::Output {
        match self {
            Multi::Many(m) => m.iter().map(|f| -f).collect(),
            Multi::Single(m) => Multi::Single(-m),
        }
    }
}

impl<M> Neg for Multi<M>
where
    M: Neg<Output = M> + Clone,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        match self {
            Self::Many(m) => m.iter().cloned().map(|f| -f).collect(),
            Self::Single(m) => Self::Single(-m),
        }
    }
}

impl<'a, M> Add<&'a M> for &'a Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
    &'a M: Add<&'a M, Output = M> + 'a,
{
    type Output = Multi<M>;
    /// Adds this M to all Ms in the multi
    fn add(self, rhs: &'a M) -> Self::Output {
        match self {
            Multi::Many(m) => m.iter().map(|m| rhs.clone() + m).collect(),
            Multi::Single(m) => Multi::Single(m + rhs),
        }
    }
}

impl<'a, M> Add<&'a M> for Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Self;
    /// Adds this formula to all formulas in the multi formula
    fn add(self, rhs: &'a M) -> Self::Output {
        match self {
            Self::Many(m) => m.iter().cloned().map(|m| m + rhs).collect(),
            Self::Single(m) => Self::Single(m + rhs),
        }
    }
}

impl<'a, M> Add<M> for &'a Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Multi<M>;
    /// Adds this formula to all formulas in the multi formula
    fn add(self, rhs: M) -> Self::Output {
        match self {
            Multi::Many(m) => m.iter().map(|m| rhs.clone() + m).collect(),
            Multi::Single(m) => Multi::Single(rhs + m),
        }
    }
}

impl<'a, M> Add<M> for Multi<M>
where
    M: Add<M, Output = M> + Clone + Add<&'a M, Output = M> + 'a,
{
    type Output = Self;
    /// Adds this formula to all formulas in the multi formula
    fn add(self, rhs: M) -> Self::Output {
        match self {
            Self::Many(m) => Self::Many(m.iter().cloned().map(|m| rhs.clone() + m).collect()),
            Self::Single(m) => Self::Single(m + rhs),
        }
    }
}

impl<'a, M> Sub<&'a M> for &'a Multi<M>
where
    &'a M: Sub<&'a M, Output = M> + 'a,
    M: Clone,
{
    type Output = Multi<M>;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: &'a M) -> Self::Output {
        match self {
            Multi::Many(m) => m.iter().map(|m| m - rhs).collect(),
            Multi::Single(m) => Multi::Single(m - rhs),
        }
    }
}

impl<'a, M> Sub<&'a M> for Multi<M>
where
    M: Sub<&'a M, Output = M> + Clone + 'a,
{
    type Output = Self;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: &'a M) -> Self::Output {
        match self {
            Self::Many(m) => m.iter().cloned().map(|m| m - rhs).collect(),
            Self::Single(m) => Self::Single(m - rhs),
        }
    }
}

impl<'a, M> Sub<M> for &'a Multi<M>
where
    &'a M: Sub<M, Output = M> + 'a,
    M: Clone,
{
    type Output = Multi<M>;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: M) -> Self::Output {
        match self {
            Multi::Many(m) => m.iter().map(|m| m - rhs.clone()).collect(),
            Multi::Single(m) => Multi::Single(m - rhs),
        }
    }
}

impl<'a, M> Sub<M> for Multi<M>
where
    M: Sub<M, Output = M> + Clone + 'a,
{
    type Output = Self;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: M) -> Self::Output {
        match self {
            Self::Many(m) => m.iter().cloned().map(|m| m - rhs.clone()).collect(),
            Self::Single(m) => Self::Single(m - rhs),
        }
    }
}

impl<'a, M> Mul<&'a Multi<M>> for &'a Multi<M>
where
    &'a M: Add<&'a M, Output = M> + 'a,
    &'a Multi<M>: Add<&'a M, Output = Multi<M>> + 'a,
{
    type Output = Multi<M>;
    /// Cartesian product between the two multi formulas, creating all potential combinations of the formula and giving back the sum formula/mass for each combination
    fn mul(self, rhs: &'a Multi<M>) -> Self::Output {
        match (self, rhs) {
            (Multi::Many(s), Multi::Many(o)) => s
                .iter()
                .cartesian_product(o.iter())
                .map(|(a, b)| b + a)
                .collect(),
            (Multi::Many(_), Multi::Single(o)) => self + o,
            (Multi::Single(s), Multi::Many(_)) => rhs + s,
            (Multi::Single(s), Multi::Single(o)) => Multi::Single(s + o),
        }
    }
}

impl<'a, M> Mul<&'a Self> for Multi<M>
where
    &'a M: Add<&'a M, Output = M> + Add<M, Output = M> + 'a,
    M: Add<&'a M, Output = M> + 'a + Clone,
    &'a Self: Add<&'a M, Output = Self> + Add<M, Output = Self> + 'a,
    Self: Add<&'a M, Output = Self> + 'a,
{
    type Output = Self;
    /// Cartesian product between the two multi formulas, creating all potential combinations of the formula and giving back the sum formula/mass for each combination
    fn mul(self, rhs: &'a Self) -> Self::Output {
        match (self, rhs) {
            (Self::Many(s), Self::Many(o)) => s
                .iter()
                .cloned()
                .cartesian_product(o.iter())
                .map(|(a, b)| b + a)
                .collect(),
            (s @ Self::Many(_), Self::Single(o)) => s + o,
            (Self::Single(s), Self::Many(_)) => rhs + s,
            (Self::Single(s), Self::Single(o)) => Self::Single(s + o),
        }
    }
}

impl<'a, M> Mul<Multi<M>> for &'a Multi<M>
where
    &'a M: Add<&'a M, Output = M> + Add<M, Output = M> + 'a,
    M: Add<&'a M, Output = M> + 'a + Clone,
    &'a Multi<M>: Add<&'a M, Output = Multi<M>> + Add<M, Output = Multi<M>> + 'a,
    Multi<M>: Add<&'a M, Output = Multi<M>> + 'a,
{
    type Output = Multi<M>;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: Multi<M>) -> Self::Output {
        match (self, rhs) {
            (Multi::Many(s), Multi::Many(o)) => s
                .iter()
                .cartesian_product(o.iter().cloned())
                .map(|(a, b)| b + a)
                .collect(),
            (s @ Multi::Many(_), Multi::Single(o)) => s + o,
            (Multi::Single(s), o @ Multi::Many(_)) => o + s,
            (Multi::Single(s), Multi::Single(o)) => Multi::Single(s + o),
        }
    }
}

impl<M> Mul<Self> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
    Self: Add<M, Output = Self>,
{
    type Output = Self;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Many(s), Self::Many(o)) => s
                .iter()
                .cloned()
                .cartesian_product(o.iter().cloned())
                .map(|(a, b)| b + a)
                .collect(),
            (s @ Self::Many(_), Self::Single(o)) => s + o,
            (Self::Single(s), o @ Self::Many(_)) => o + s,
            (Self::Single(s), Self::Single(o)) => Self::Single(s + o),
        }
    }
}

impl<M> MulAssign<&Self> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
    Self: Add<M, Output = Self> + Clone,
{
    fn mul_assign(&mut self, rhs: &Self) {
        *self = self.clone() * rhs.clone();
    }
}

impl<M> MulAssign<Self> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
    Self: Add<M, Output = Self> + Clone,
{
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs;
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
        Self::Single(value)
    }
}

impl<M: Clone> From<&M> for Multi<M> {
    fn from(value: &M) -> Self {
        Self::Single(value.clone())
    }
}

impl<M> From<Vec<M>> for Multi<M> {
    fn from(mut value: Vec<M>) -> Self {
        if value.len() == 1 {
            Self::Single(value.pop().unwrap())
        } else if value.is_empty() {
            panic!("Cannot construct an empty Multi")
        } else {
            Self::Many(value.into())
        }
    }
}

impl<M: Clone> From<&[M]> for Multi<M> {
    fn from(value: &[M]) -> Self {
        if value.len() == 1 {
            Self::Single(value[0].clone())
        } else if value.is_empty() {
            panic!("Cannot construct an empty Multi")
        } else {
            Self::Many(value.into())
        }
    }
}

impl<'a, M> From<&'a [Self]> for Multi<M>
where
    Self: Mul<&'a Self, Output = Self> + Add<&'a M, Output = Self> + 'a,
    M: Default + Clone + Add<M, Output = M>,
    &'a M: Add<&'a M, Output = M> + 'a,
    &'a Self: Add<&'a M, Output = Self> + 'a,
{
    /// Get all potential combination from a series of multi elements. If the series is empty it returns the default element.
    fn from(value: &'a [Self]) -> Self {
        value.iter().fold(Self::default(), |acc, a: &Self| match a {
            Self::Many(_) => acc * a,
            Self::Single(m) => acc + m,
        })
    }
}

impl<M> std::iter::FromIterator<M> for Multi<M> {
    fn from_iter<T: IntoIterator<Item = M>>(iter: T) -> Self {
        Self::from(iter.into_iter().collect_vec())
    }
}

impl<'a, M: Clone + 'a> std::iter::FromIterator<&'a M> for Multi<M> {
    fn from_iter<T: IntoIterator<Item = &'a M>>(iter: T) -> Self {
        Self::from(iter.into_iter().cloned().collect_vec())
    }
}
