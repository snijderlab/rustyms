//! The measurement system used in this crate.
//! A redefinition of the important SI units for them to be stored in a more sensible base unit for MS purposes.

#![allow(clippy::incorrect_clone_impl_on_copy_type)]
#![allow(clippy::ignored_unit_patterns)]
use uom::*;

pub use f64::*;

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
            @r: 1.0; "r", "ratio", "ratios";
        }
    }
}

system! {
    /// Quantities
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
pub mod f64 {
    mod mks {
        pub use super::super::*;
    }

    Q!(self::mks, f64);

    pub use super::charge::e;
    pub use super::mass::dalton;
    pub use super::mass_over_charge::mz;
    pub use super::ratio::r;
    pub use super::time::s;

    /// Annotate the given number as being in Da
    #[allow(dead_code)]
    pub fn da(v: f64) -> Mass {
        Mass::new::<super::mass::dalton>(v)
    }
}

impl MassOverCharge {
    /// Absolute ppm error between this number and the given other
    pub fn ppm(self, b: Self) -> f64 {
        ((self - b).abs() / self).value * 1e6
    }
}
impl Mass {
    /// Absolute ppm error between this number and the given other
    pub fn ppm(self, b: Self) -> f64 {
        ((self - b).abs() / self).value * 1e6
    }
}
