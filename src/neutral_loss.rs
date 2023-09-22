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
    /// Gain of water
    WaterGain,
    /// Loss of 2 water
    DiWater,
    /// Gain of 2 water
    DiWaterGain,
    /// Loss of 3 water
    TriWater,
    /// Gain of 3 water
    TriWaterGain,
    /// Loss of ammonia (NH3)
    Ammonia,
    /// Loss of carbon monoxide
    CarbonMonoxide,
    /// Loss of hydrogen
    Hydrogen,
    /// Gain of hydrogen
    HydrogenGain,
    /// Loss of 2 hydrogens
    Capital,
    /// Gain of 2 hydrogens
    CapitalGain,
    /// Loss of 3 hydrogens
    TriHydrogen,
    /// Gain of 3 hydrogens
    TriHydrogenGain,
}

impl Chemical for NeutralLoss {
    fn formula(&self) -> crate::formula::MolecularFormula {
        match self {
            Self::Water => molecular_formula!(O 1 H 2),
            Self::DiWater => molecular_formula!(O 2 H 4),
            Self::TriWater => molecular_formula!(O 3 H 6),
            Self::Ammonia => molecular_formula!(N 1 H 3),
            Self::CarbonMonoxide => molecular_formula!(C 1 O 1),
            Self::Hydrogen => molecular_formula!(H 1),
            Self::Capital => molecular_formula!(H 2),
            Self::HydrogenGain => molecular_formula!(H - 1),
            Self::CapitalGain => molecular_formula!(H - 2),
            Self::WaterGain => molecular_formula!(O -1 H -2),
            Self::DiWaterGain => molecular_formula!(O -2 H -4),
            Self::TriWaterGain => molecular_formula!(O -3 H -6),
            Self::TriHydrogen => molecular_formula!(H 3),
            Self::TriHydrogenGain => molecular_formula!(H - 3),
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
                Self::DiWater => "DiWater",
                Self::TriWater => "TriWater",
                Self::Ammonia => "Ammonia",
                Self::CarbonMonoxide => "CarbonMonoxide",
                Self::Hydrogen => "Hydrogen",
                Self::Capital => "Capital",
                Self::HydrogenGain => "HydrogenGain",
                Self::CapitalGain => "CapitalGain",
                Self::WaterGain => "WaterGain",
                Self::DiWaterGain => "DiWaterGain",
                Self::TriWaterGain => "TriWaterGain",
                Self::TriHydrogen => "TriHydrogen",
                Self::TriHydrogenGain => "TriHydrogenGain",
            }
        )
    }
}
