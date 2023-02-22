//! Source: https://www.ciaaw.org/atomic-weights.htm
//! When a range of weight is given the average value of the top and bottom is used.
//! All values are given in dalton.
#![allow(dead_code)]
#![allow(non_upper_case_globals)]

pub const H: f64 = 1.007975;
pub const C: f64 = 12.0106;
pub const N: f64 = 14.006855;
pub const O: f64 = 15.9994;
pub const F: f64 = 18.9984031625;
pub const P: f64 = 30.9737619985;
pub const S: f64 = 32.0675;
pub const Se: f64 = 78.9718;

// Common combined pieces (for easy reading of calculations)
pub const CO: f64 = C + O;
pub const OH: f64 = O + H;
pub const NH: f64 = N + H;
pub const NH2: f64 = N + H * 2.0;
pub const NH3: f64 = N + H * 3.0;
pub const CH: f64 = C + H;
pub const CH2: f64 = C + H * 2.0;
pub const CH3: f64 = C + H * 3.0;
pub const BACKBONE: f64 = CO + CH + NH;
