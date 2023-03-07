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

    #[must_use]
    pub fn new(theoretical_mass: Mass, charge: Charge, ion: FragmentType) -> Self {
        Self {
            theoretical_mass,
            charge,
            ion,
        }
    }

    #[must_use]
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
    pub const fn sequence_index(&self) -> Option<usize> {
        match self {
            Self::a(n)
            | Self::b(n)
            | Self::c(n)
            | Self::d(n)
            | Self::v(n)
            | Self::w(n)
            | Self::x(n)
            | Self::y(n)
            | Self::z(n)
            | Self::z·(n) => Some(*n),
            Self::precursor => None,
        }
    }
}

impl Display for FragmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::a(_) => "a",
                Self::b(_) => "b",
                Self::c(_) => "c",
                Self::d(_) => "d",
                Self::v(_) => "v",
                Self::w(_) => "w",
                Self::x(_) => "x",
                Self::y(_) => "y",
                Self::z(_) => "z",
                Self::z·(_) => "z·",
                Self::precursor => "precursor",
            }
        )
    }
}
