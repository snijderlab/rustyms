//pub use charge::*;
//pub use mass::*;
//pub use mass_over_charge::*;
use uom::*;

#[macro_use]
pub mod mass {
    use uom::*;

    quantity! {
        quantity: Mass; "mass";
        dimension: Q< P1, Z0, Z0>;
        units {
            @dalton: 1.0; "da", "dalton", "daltons";
        }
    }
}

#[macro_use]
pub mod charge {
    use uom::*;

    quantity! {
        quantity: Charge; "charge";
        dimension: Q< Z0, P1, Z0>;
        units {
            @e: 1.0; "e", "atomic_unit_of_charge", "atomic_units_of_charge";
        }
    }
}

#[macro_use]
pub mod time {
    use uom::*;

    quantity! {
        quantity: Time; "time";
        dimension: Q< Z0, Z0, P1>;
        units {
            @s: 1.0; "s", "second", "seconds";
        }
    }
}

#[macro_use]
pub mod mass_over_charge {
    use uom::*;

    quantity! {
        quantity: MassOverCharge; "mass_over_charge";
        dimension: Q< P1, N1, Z0>;
        units {
            @mz: 1.0; "mz", "mass_over_charge", "mass_over_charge";
        }
    }
}

#[macro_use]
pub mod ratio {
    use uom::*;

    quantity! {
        quantity: Ratio; "ratio";
        dimension: Q< Z0, Z0, Z0>;
        units {
            @r: 1.0; "r", "ratio", "ratios";
        }
    }
}

system! {
    quantities: Q {
        mass: dalton, M;
        charge: e, C;
        time: s, T;
    }

    units: U {
        mod mass::Mass,
        mod charge::Charge,
        mod time::Time,
        mod mass_over_charge::MassOverCharge,
        mod ratio::Ratio,
    }
}

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

    pub fn da(v: f64) -> Mass {
        Mass::new::<super::mass::dalton>(v)
    }
}
