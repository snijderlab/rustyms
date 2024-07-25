use crate::{system::isize::Charge, Chemical, Element, MolecularFormula, SequencePosition};
use serde::{Deserialize, Serialize};

/// A selection of ions that together define the charge of a peptide
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct MolecularCharge {
    /// The ions that together define the charge of the peptide.
    /// The first number is the amount of times this adduct ion occurs, the molecular formula is the full formula for the adduct ion.
    /// The charge for each ion is saved as the number of electrons missing or gained in the molecular formula.
    pub charge_carriers: Vec<(isize, MolecularFormula)>,
}

impl MolecularCharge {
    /// Create a default charge state with only protons
    #[allow(clippy::missing_panics_doc)] // Cannot panic
    pub fn proton(charge: isize) -> Self {
        Self {
            charge_carriers: vec![(
                charge,
                MolecularFormula::new(&[(Element::H, None, 1), (Element::Electron, None, -1)], &[])
                    .unwrap(),
            )],
        }
    }

    /// Create a charge state with the given ions
    pub fn new(charge_carriers: &[(isize, MolecularFormula)]) -> Self {
        Self {
            charge_carriers: charge_carriers.to_vec(),
        }
    }

    /// Generate all possible charge carrier options that have a charge of one,
    /// used for small ions (immonium/diagnostic/glycan diagnostic) as these are
    /// not expected to ever have a higher charge.
    pub fn all_single_charge_options(&self) -> Vec<Self> {
        self.charge_carriers
            .iter()
            .filter_map(|(n, c)| {
                (*n > 0 && c.charge().value == 1).then_some(Self::new(&[(1, c.clone())]))
            })
            .collect()
    }

    /// Generate all possible charge carrier options from the selection of ions for use in fragment calculations
    pub fn all_charge_options(&self) -> Vec<Self> {
        let mut options: Vec<Vec<(isize, MolecularFormula)>> = Vec::new();
        for carrier in &self.charge_carriers {
            let mut new_options = Vec::new();
            if options.is_empty() {
                for n in 0..=carrier.0 {
                    new_options.push(vec![(n, carrier.1.clone())]);
                }
            } else {
                for n in 0..=carrier.0 {
                    for o in &options {
                        let mut new = o.clone();
                        new.push((n, carrier.1.clone()));
                        new_options.push(new);
                    }
                }
            }
            options = new_options;
        }
        // Ignore the case where no charge carriers are present at all and combine the cases into single molecular formulas
        options
            .into_iter()
            .filter(|o| o.iter().map(|o| o.0).sum::<isize>() != 0)
            .map(|charge_carriers| Self { charge_carriers })
            .collect()
    }

    /// Get the total charge of these charge carriers
    pub fn charge(&self) -> Charge {
        self.charge_carriers
            .iter()
            .fold(Charge::default(), |acc, (amount, formula)| {
                acc + *amount * formula.charge()
            })
    }
}

impl Chemical for MolecularCharge {
    fn formula(
        &self,
        _sequence_index: SequencePosition,
        _peptide_index: usize,
    ) -> MolecularFormula {
        self.charge_carriers
            .iter()
            .map(|(n, mol)| mol.clone() * *n as i32)
            .sum::<MolecularFormula>()
    }
}

impl std::fmt::Display for MolecularCharge {
    /// Is not guaranteed to fully conform to the Pro Forma standard. Because the data structure accepts more than the standard.
    /// So adducts with other than +1/-1 charge states, or adducts with complex formula (not a single element) will not adhere to the standard.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.charge_carriers
                .iter()
                .map(|c| c.1.charge().value * c.0)
                .sum::<isize>()
        )?;
        if !self.charge_carriers.iter().all(|c| {
            c.1 == MolecularFormula::new(
                &[(Element::H, None, 1), (Element::Electron, None, -1)],
                &[],
            )
            .unwrap()
        }) {
            write!(f, "[")?;
            let mut first = true;
            for (amount, formula) in &self.charge_carriers {
                if first {
                    first = false;
                } else {
                    write!(f, ",")?;
                }
                let charge = formula.charge().value;
                write!(f, "{amount}{formula}{charge:+}")?;
            }
            write!(f, "]")?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use crate::{Chemical, SequencePosition};

    use super::MolecularCharge;

    #[test]
    fn simple_charge_options() {
        let mc = MolecularCharge::new(&[(1, molecular_formula!(H 1 Electron -1))]);
        let options = mc.all_charge_options();
        assert_eq!(options.len(), 1);
        assert_eq!(
            options[0].formula(SequencePosition::default(), 0),
            molecular_formula!(H 1 Electron -1)
        );
    }
}
