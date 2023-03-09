use std::fmt::{Debug, Display};

use crate::{model::NeutralLoss, system::f64::*, HasMass, MassSystem, MonoIsotopic};

#[derive(Clone)]
pub struct Fragment {
    pub theoretical_mass: Mass,
    pub charge: Charge,
    pub ion: FragmentType,
    pub neutral_loss: Option<NeutralLoss>,
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
            neutral_loss: None,
        }
    }

    #[must_use]
    pub fn with_charge<M: MassSystem>(&self, charge: Charge) -> Self {
        Self {
            theoretical_mass: self.theoretical_mass + da(M::Proton * charge.value),
            charge,
            ion: self.ion,
            neutral_loss: self.neutral_loss,
        }
    }

    #[must_use]
    pub fn with_neutral_loss<M: MassSystem>(&self, neutral_loss: &NeutralLoss) -> Self {
        Self {
            theoretical_mass: self.theoretical_mass - neutral_loss.mass::<M>(),
            charge: self.charge,
            ion: self.ion,
            neutral_loss: Some(*neutral_loss),
        }
    }

    #[must_use]
    pub fn with_neutral_losses<M: MassSystem>(&self, neutral_losses: &[NeutralLoss]) -> Vec<Self> {
        let mut output = Vec::with_capacity(neutral_losses.len() + 1);
        output.push(self.clone());
        output.extend(
            neutral_losses
                .iter()
                .map(|loss| self.with_neutral_loss::<M>(loss)),
        );
        output
    }
}

impl Debug for Fragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:?} {} {:+}{}",
            self.ion,
            self.mz().value,
            self.charge.value,
            self.neutral_loss
                .map_or(String::new(), |loss| format!(" -{loss}"))
        )
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct Position {
    pub sequence_index: usize,
    pub series_number: usize,
}

impl Position {
    pub const fn n(sequence_index: usize, _length: usize) -> Self {
        Self {
            sequence_index,
            series_number: sequence_index + 1,
        }
    }
    pub const fn c(sequence_index: usize, length: usize) -> Self {
        Self {
            sequence_index,
            series_number: length - sequence_index,
        }
    }
}

#[derive(Clone, Copy, Eq, PartialEq, Hash)]
#[allow(non_camel_case_types)]
pub enum FragmentType {
    a(Position),
    b(Position),
    c(Position),
    d(Position),
    v(Position),
    w(Position),
    x(Position),
    y(Position),
    z(Position),
    z·(Position),
    precursor,
}

impl FragmentType {
    pub const fn position(&self) -> Option<&Position> {
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
            | Self::z·(n) => Some(n),
            Self::precursor => None,
        }
    }
}

impl Debug for FragmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}",
            self,
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
                | Self::z·(n) => n.series_number.to_string(),
                Self::precursor => String::new(),
            }
        )
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

#[test]
fn neutral_loss() {
    let a = Fragment::new(
        Mass::new::<dalton>(118.0),
        Charge::new::<e>(1.0),
        FragmentType::precursor,
    );
    let loss = a.with_neutral_losses::<MonoIsotopic>(&[NeutralLoss::Water]);
    dbg!(&a, &loss);
    assert_eq!(a.theoretical_mass, loss[0].theoretical_mass);
    assert_eq!(
        a.theoretical_mass,
        loss[1].theoretical_mass + NeutralLoss::Water.mass::<MonoIsotopic>()
    );
}
