//! The measurement system used in this crate.
//! A redefinition of the important SI units for them to be stored in a more sensible base unit for MS purposes.

#![allow(clippy::non_canonical_clone_impl)]
#![allow(clippy::ignored_unit_patterns)]
use std::ops::{Deref, DerefMut};

use num_traits::Zero;
use uom::*;

use serde::{Deserialize, Serialize};

use crate::helper_functions;

pub use self::f64::*;

/// The mass quantity in dalton
#[macro_use]
pub mod mass {
    use uom::*;

    quantity! {
        /// Mass in dalton
        quantity: Mass; "mass";
        /// Mass
        dimension: Q< P1, Z0, Z0>;
        units {
            @millidalton: 0.001; "mDa", "millidalton", "millidaltons";
            @dalton: 1.0; "Da", "dalton", "daltons";
            @kilodalton: 1_000.0; "kDa", "kilodalton", "kilodaltons";
            @megadalton: 1_000_000.0; "MDa", "megadalton", "megadaltons";
        }
    }
}

/// The charge quantity in atomic units of charge aka electrons
#[macro_use]
pub mod charge {
    use uom::*;

    quantity! {
        /// Charge in electrons
        quantity: Charge; "charge";
        /// Charge
        dimension: Q< Z0, P1, Z0>;
        units {
            @e: 1.0; "e", "atomic_unit_of_charge", "atomic_units_of_charge";
        }
    }
}

/// The time quantity in seconds
#[macro_use]
pub mod time {
    use uom::*;

    quantity! {
        /// Time (s)
        quantity: Time; "time";
        /// Time
        dimension: Q< Z0, Z0, P1>;
        units {
            @ns: 0.000_000_001; "ns", "nanosecond", "nanoseconds";
            @μs: 0.000_001; "μs", "microsecond", "microseconds";
            @ms: 0.001; "ms", "millisecond", "milliseconds";
            @s: 1.0; "s", "second", "seconds";
            @min: 60.0; "min", "minute", "minutes";
            @h: 3600.0; "h", "hour", "hours";
        }
    }
}

/// The mass over charge quantity
#[macro_use]
pub mod mass_over_charge {
    use uom::*;

    quantity! {
        /// Mass over charge (da/e)
        quantity: MassOverCharge; "mass_over_charge";
        /// Mass over charge (da/e)
        dimension: Q< P1, N1, Z0>;
        units {
            @mz: 1.0; "mz", "mass_over_charge", "mass_over_charge";
        }
    }
}

/// A unit less quantity for use in general calculations
#[macro_use]
pub mod ratio {
    use uom::*;

    quantity! {
        /// Unit less quantity for general calculations
        quantity: Ratio; "ratio";
        /// Unit less quantity for general calculations
        dimension: Q< Z0, Z0, Z0>;
        units {
            @fraction: 1.0; "⅟", "fraction", "fraction";
            @percent: 0.01; "%", "percent", "percent";
            @promille: 0.01; "‰", "promille", "promille";
            @ppm: 0.000_001; "ppm", "ppm", "ppm";
            @ppb: 0.000_000_001; "ppb", "ppb", "ppb";
            @ppt: 0.000_000_000_001; "ppt", "ppt", "ppt";
            @ppq: 0.000_000_000_000_001; "ppq", "ppq", "ppq";
        }
    }
}

system! {
    /// Quantities
    #[doc(hidden)]
    quantities: Q {
        mass: dalton, M;
        charge: e, C;
        time: s, T;
    }

    /// Units
    units: U {
        mod mass::Mass,
        mod charge::Charge,
        mod time::Time,
        mod mass_over_charge::MassOverCharge,
        mod ratio::Ratio,
    }
}

/// The whole system with f64 as storage type
#[allow(unused_imports)]
pub mod f64 {
    mod mks {
        pub use super::super::*;
    }

    Q!(self::mks, f64);

    pub use super::charge::e;
    pub use super::mass::dalton;
    pub use super::mass_over_charge::mz;
    pub use super::ratio::fraction;
    pub use super::time::s;

    /// Annotate the given number as being in Da
    #[allow(dead_code)]
    pub fn da(v: f64) -> Mass {
        Mass::new::<super::mass::dalton>(v)
    }
}

/// All quantities with usize as underlying type
#[allow(unused_imports)]
pub mod usize {
    mod mks {
        pub use super::super::*;
    }

    Q!(self::mks, usize);

    pub use super::charge::e;
    pub use super::mass::dalton;
    pub use super::mass_over_charge::mz;
    pub use super::ratio::fraction;
    pub use super::time::s;
}

impl MassOverCharge {
    /// Absolute ppm error between this number and the given other
    pub fn ppm(self, b: Self) -> Ratio {
        Ratio::new::<crate::system::ratio::ppm>(((self - b).abs() / self.abs()).value * 1e6)
    }
}
impl Mass {
    /// Absolute ppm error between this number and the given other
    pub fn ppm(self, b: Self) -> Ratio {
        Ratio::new::<crate::system::ratio::ppm>(((self - b).abs() / self.abs()).value * 1e6)
    }
}

/// A wrapper around [`Ratio`] which implements Eq/Ord/Hash to help in auto deriving these on other structs.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct OrderedRatio(Ratio);

impl OrderedRatio {
    /// Use the zero from [`Ratio`] itself
    pub fn zero() -> Self {
        Self(Ratio::zero())
    }

    /// Get a normal [`Ratio`]
    #[allow(dead_code)]
    pub fn into_inner(self) -> Ratio {
        self.0
    }
}

impl Default for OrderedRatio {
    fn default() -> Self {
        Self::zero()
    }
}

impl From<Ratio> for OrderedRatio {
    fn from(value: Ratio) -> Self {
        Self(value)
    }
}

impl Deref for OrderedRatio {
    type Target = Ratio;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for OrderedRatio {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Ord for OrderedRatio {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.value.total_cmp(&other.0.value)
    }
}

impl PartialOrd for OrderedRatio {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for OrderedRatio {}

impl PartialEq for OrderedRatio {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl std::hash::Hash for OrderedRatio {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        helper_functions::f64_bits(self.0.value).hash(state);
    }
}
/// A wrapper around [`Mass`] which implements Eq/Ord/Hash to help in auto deriving these on other structs.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct OrderedMass(Mass);

impl OrderedMass {
    /// Use the zero from [`Mass`] itself
    pub fn zero() -> Self {
        Self(Mass::zero())
    }

    /// Get a normal [`Mass`]
    #[allow(dead_code)]
    pub fn into_inner(self) -> Mass {
        self.0
    }
}

impl Default for OrderedMass {
    fn default() -> Self {
        Self::zero()
    }
}

impl From<Mass> for OrderedMass {
    fn from(value: Mass) -> Self {
        Self(value)
    }
}

impl Deref for OrderedMass {
    type Target = Mass;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for OrderedMass {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Ord for OrderedMass {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.value.total_cmp(&other.0.value)
    }
}

impl PartialOrd for OrderedMass {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for OrderedMass {}

impl PartialEq for OrderedMass {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl std::hash::Hash for OrderedMass {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        helper_functions::f64_bits(self.0.value).hash(state);
    }
}
