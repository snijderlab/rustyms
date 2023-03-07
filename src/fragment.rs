use std::fmt::{Debug, Display};

use crate::{system::f64::*, MassSystem};

#[derive(Clone)]
pub struct Fragment {
    pub theoretical_mass: Mass,
    pub charge: Charge,
    pub ion: FragmentType,
}

impl Fragment {
    pub fn mz(&self) -> MassOverCharge {
        self.theoretical_mass / self.charge
    }

    pub fn new(theoretical_mass: Mass, charge: Charge, ion: FragmentType) -> Self {
        Self {
            theoretical_mass,
            charge,
            ion,
        }
    }

    pub fn with_charge<M: MassSystem>(&self, charge: Charge) -> Self {
        Self {
            theoretical_mass: self.theoretical_mass + da(M::Proton * charge.value),
            charge,
            ion: self.ion,
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

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
#[allow(non_camel_case_types)]
pub enum FragmentType {
    a(usize),
    b(usize),
    c(usize),
    d(usize),
    v(usize),
    w(usize),
    x(usize),
    y(usize),
    z(usize),
    z·(usize),
    precursor,
}

impl FragmentType {
    pub fn sequence_index(&self) -> Option<usize> {
        match self {
            FragmentType::a(n) => Some(*n),
            FragmentType::b(n) => Some(*n),
            FragmentType::c(n) => Some(*n),
            FragmentType::d(n) => Some(*n),
            FragmentType::v(n) => Some(*n),
            FragmentType::w(n) => Some(*n),
            FragmentType::x(n) => Some(*n),
            FragmentType::y(n) => Some(*n),
            FragmentType::z(n) => Some(*n),
            FragmentType::z·(n) => Some(*n),
            FragmentType::precursor => None,
        }
    }
}

impl Display for FragmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                FragmentType::a(_) => "a",
                FragmentType::b(_) => "b",
                FragmentType::c(_) => "c",
                FragmentType::d(_) => "d",
                FragmentType::v(_) => "v",
                FragmentType::w(_) => "w",
                FragmentType::x(_) => "x",
                FragmentType::y(_) => "y",
                FragmentType::z(_) => "z",
                FragmentType::z·(_) => "z·",
                FragmentType::precursor => "precursor",
            }
        )
    }
}
