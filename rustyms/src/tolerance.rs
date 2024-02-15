use std::{fmt::Display, str::FromStr};

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::{
    system::{da, Mass, OrderedMass},
    Multi,
};

/// A tolerance around a given mass for searching purposes
#[allow(non_camel_case_types)]
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub enum Tolerance {
    /// A relative search tolerance in parts per million
    ppm(OrderedFloat<f64>),
    /// An absolute tolerance defined by a constant offset from the mass (bounds are mass - tolerance, mass + tolerance)
    Abs(OrderedMass),
}

impl Tolerance {
    /// Create a new ppm value
    pub fn new_ppm(value: f64) -> Self {
        Self::ppm(value.into())
    }

    /// Create a new absolute value
    pub fn new_absolute(value: Mass) -> Self {
        Self::Abs(value.into())
    }

    /// Find the bounds around a given mass for this tolerance
    pub fn bounds(&self, mass: Mass) -> (Mass, Mass) {
        match self {
            Self::ppm(ppm) => (
                da(mass.value * (1.0 - ppm.into_inner() / 1e6)),
                da(mass.value * (1.0 + ppm.into_inner() / 1e6)),
            ),
            Self::Abs(tolerance) => (mass - tolerance.into_inner(), mass + tolerance.into_inner()),
        }
    }
}

impl Display for Tolerance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Abs(mass) => format!("{} da", mass.value),
                Self::ppm(ppm) => format!("{ppm} ppm"),
            }
        )
    }
}

impl FromStr for Tolerance {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let num_str = String::from_utf8(
            s.bytes()
                .take_while(|c| {
                    c.is_ascii_digit()
                        || *c == b'.'
                        || *c == b'-'
                        || *c == b'+'
                        || *c == b'e'
                        || *c == b'E'
                })
                .collect::<Vec<_>>(),
        )
        .map_err(|_| ())?;
        let num = num_str.parse::<f64>().map_err(|_| ())?;
        match s[num_str.len()..].trim() {
            "ppm" => Ok(Self::ppm(num.into())),
            "da" => Ok(Self::Abs(da(num).into())),
            _ => Err(()),
        }
    }
}

impl TryFrom<&str> for Tolerance {
    type Error = ();
    fn try_from(value: &str) -> Result<Self, ()> {
        value.parse()
    }
}

/// Check if two values are within the specified tolerance from each other.
pub trait MassComparable<A, B> {
    /// Check if two values are within the specified tolerance from each other.
    fn within(&self, a: &A, b: &B) -> bool;
}

impl MassComparable<Mass, Mass> for Tolerance {
    fn within(&self, a: &Mass, b: &Mass) -> bool {
        match self {
            Self::Abs(tol) => (a.value - b.value).abs() <= tol.value,
            Self::ppm(ppm) => a.ppm(*b) <= ppm.into_inner(),
        }
    }
}

impl MassComparable<Multi<Mass>, Multi<Mass>> for Tolerance {
    fn within(&self, a: &Multi<Mass>, b: &Multi<Mass>) -> bool {
        a.iter()
            .cartesian_product(b.iter())
            .any(|(a, b)| self.within(a, b))
    }
}

impl MassComparable<Multi<Mass>, Mass> for Tolerance {
    fn within(&self, a: &Multi<Mass>, b: &Mass) -> bool {
        a.iter().any(|a| self.within(a, b))
    }
}

impl MassComparable<Mass, Multi<Mass>> for Tolerance {
    fn within(&self, a: &Mass, b: &Multi<Mass>) -> bool {
        b.iter().any(|b| self.within(a, b))
    }
}
