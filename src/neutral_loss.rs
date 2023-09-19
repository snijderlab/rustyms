use std::fmt::Display;

use serde::{Deserialize, Serialize};

use crate::{
    formula::{Chemical, MolecularFormula},
    Element,
};

/// All possible neutral losses
#[derive(Debug, Clone, Copy, Eq, PartialEq, Serialize, Deserialize)]
pub enum NeutralLoss {
    /// Loss of water
    Water,
    /// Loss of ammonia (NH3)
    Ammonia,
    /// Loss of carbon monoxide
    CarbonMonoxide,
    /// Loss of hydrogen
    Hydrogen,
    /// Gain of hydrogen
    HydrogenGain,
    /// Gain of water
    WaterGain,
}

impl Chemical for NeutralLoss {
    fn formula(&self) -> crate::formula::MolecularFormula {
        match self {
            Self::Water => molecular_formula!(O 1 H 2),
            Self::Ammonia => molecular_formula!(N 1 H 3),
            Self::CarbonMonoxide => molecular_formula!(C 1 O 1),
            Self::Hydrogen => molecular_formula!(H 1),
            Self::HydrogenGain => molecular_formula!(H - 1),
            Self::WaterGain => molecular_formula!(O -1 H -2),
        }
    }
}
impl Display for NeutralLoss {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Water => "Water",
                Self::Ammonia => "Ammonia",
                Self::CarbonMonoxide => "CarbonMonoxide",
                Self::Hydrogen => "Hydrogen",
                Self::HydrogenGain => "HydrogenGain",
                Self::WaterGain => "WaterGain",
            }
        )
    }
}
