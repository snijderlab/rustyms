//pub use charge::*;
//pub use mass::*;
//pub use mass_over_charge::*;
use uom::*;

#[macro_use]
pub mod mass {
    use uom::*;

    quantity! {
        quantity: Mass; "mass";
        dimension: Q< P1, Z0>;
        units {
            @dalton: 1.0; "dalton", "da", "Da";
        }
    }
}

#[macro_use]
pub mod charge {
    use uom::*;

    quantity! {
        quantity: Charge; "charge";
        dimension: Q< Z0, P1>;
        units {
            @e: 1.0; "e", "auc", "atomic_unit_of_charge";
        }
    }
}

#[macro_use]
pub mod mass_over_charge {
    use uom::*;

    quantity! {
        quantity: MassOverCharge; "mass_over_charge";
        dimension: Q< P1, N1>;
        units {
            @mz: 1.0; "mz", "Th", "mass_over_charge";
        }
    }
}

system! {
    quantities: Q {
        mass: dalton, M;
        charge: e, C;
    }

    units: U {
        mod mass::Mass,
        mod charge::Charge,
        mod mass_over_charge::MassOverCharge,
    }
}

pub mod f64 {
    mod mks {
        pub use super::super::*;
    }

    Q!(self::mks, f64);

    pub fn da(v: f64) -> Mass {
        Mass::new::<super::mass::dalton>(v)
    }
}
