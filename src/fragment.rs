use std::fmt::Debug;

use crate::system::f64::*;

pub struct Fragment {
    pub mass: Mass,
    pub sequence_index: usize,
    pub charge: Charge,
    pub ion: FragmentType,
}

impl Fragment {
    pub fn mz(&self) -> MassOverCharge {
        self.mass / self.charge
    }

    pub fn new(mass: Mass, charge: Charge, idx: usize, ion: FragmentType) -> Self {
        Self {
            mass,
            charge,
            sequence_index: idx,
            ion,
        }
    }
}

impl Debug for Fragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:?} {}+{}",
            self.ion,
            self.mz().value,
            self.charge.value
        )
    }
}

#[derive(Debug)]
#[allow(non_camel_case_types)]
pub enum FragmentType {
    a,
    b,
    c,
    d,
    v,
    w,
    x,
    y,
    z,
}
