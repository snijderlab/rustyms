#![allow(dead_code)]
#![allow(non_upper_case_globals)]

pub trait MassSystem {
    const H: f64;
    const C: f64;
    const N: f64;
    const O: f64;
    const F: f64;
    const P: f64;
    const S: f64;
    const Se: f64;
    
    // Common combined pieces (for easy reading of calculations)
    const CO: f64 = Self::C + Self::O;
    const OH: f64 = Self::O + Self::H;
    const NH: f64 = Self::N + Self::H;
    const NH2: f64 = Self::N + Self::H * 2.0;
    const NH3: f64 = Self::N + Self::H * 3.0;
    const CH: f64 = Self::C + Self::H;
    const CH2: f64 = Self::C + Self::H * 2.0;
    const CH3: f64 = Self::C + Self::H * 3.0;
    const BACKBONE: f64 = Self::CO + Self::CH + Self::NH;
}

/// Source: https://www.ciaaw.org/atomic-weights.htm
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
