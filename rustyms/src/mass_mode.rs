use serde::{Deserialize, Serialize};

/// The mode of mass to use
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug, Serialize, Deserialize,
)]
pub enum MassMode {
    /// Monoisotopic mass, use the base isotope to calculate the mass (eg always 12C)
    #[default]
    Monoisotopic,
    /// The average weight, the average between all occurring isotopes (eg something in between 12C and 13C depending on the number of C)
    Average,
    /// The most abundant mass, the most abundant single isotopic species (eg 12C or 13C depending on the number of C)
    MostAbundant,
}
