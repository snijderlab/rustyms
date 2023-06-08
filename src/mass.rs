#![allow(dead_code)]
#![allow(non_upper_case_globals)]
#![allow(clippy::unreadable_literal)]

use uom::num_traits::Zero;

use crate::Mass;

pub trait MassSystem {
    const H: f64;
    const C: f64;
    const N: f64;
    const O: f64;
    const F: f64;
    const P: f64;
    const S: f64;
    const Se: f64;
    const e: f64 = 5.48579909065e-4;

    // Common combined pieces (for easy reading of calculations)
    const Proton: f64 = Self::H - Self::e;
    const CO: f64 = Self::C + Self::O;
    const OH: f64 = Self::O + Self::H;
    const NH: f64 = Self::N + Self::H;
    const NH2: f64 = Self::N + Self::H * 2.0;
    const NH3: f64 = Self::N + Self::H * 3.0;
    const CH: f64 = Self::C + Self::H;
    const CH2: f64 = Self::C + Self::H * 2.0;
    const CH3: f64 = Self::C + Self::H * 3.0;
    const CHO: f64 = Self::C + Self::H + Self::O;
    const BACKBONE: f64 = Self::CO + Self::CH + Self::NH;
}

/// Source: [CIAAW](https://www.ciaaw.org/atomic-weights.htm)
/// When a range of weight is given the average value of the top and bottom is used.
/// All values are given in dalton.
pub struct AverageWeight {}

impl MassSystem for AverageWeight {
    const H: f64 = 1.007975;
    const C: f64 = 12.0106;
    const N: f64 = 14.006855;
    const O: f64 = 15.9994;
    const F: f64 = 18.9984031625;
    const P: f64 = 30.9737619985;
    const S: f64 = 32.0675;
    const Se: f64 = 78.9718;
}

pub struct MonoIsotopic {}

impl MassSystem for MonoIsotopic {
    const H: f64 = 1.007825031898;
    const C: f64 = 12.0;
    const N: f64 = 14.00307400425;
    const O: f64 = 15.99491461926;
    const F: f64 = 18.99840316207;
    const P: f64 = 30.97376199768;
    const S: f64 = 31.97207117354;
    const Se: f64 = 79.916521761;
}

pub struct Hecklib {}

impl MassSystem for Hecklib {
    const H: f64 = 1.00782503214;
    const C: f64 = 12.0;
    const N: f64 = 14.00307400524;
    const O: f64 = 15.99491462210;
    const F: f64 = 18.99840320500;
    const P: f64 = 30.97376151200;
    const S: f64 = 31.97207069000;
    const Se: f64 = 79.9165218;
}

pub trait HasMass {
    fn mass<M: MassSystem>(&self) -> Mass;
}

impl<T: HasMass> HasMass for Option<T> {
    fn mass<M: MassSystem>(&self) -> Mass {
        self.as_ref().map_or_else(Mass::zero, HasMass::mass::<M>)
    }
}
impl<T: HasMass, U: HasMass> HasMass for (T, U) {
    fn mass<M: MassSystem>(&self) -> Mass {
        self.0.mass::<M>() + self.1.mass::<M>()
    }
}
impl<T: HasMass> HasMass for [T] {
    fn mass<M: MassSystem>(&self) -> Mass {
        self.iter().fold(Mass::zero(), |acc, i| acc + i.mass::<M>())
    }
}
