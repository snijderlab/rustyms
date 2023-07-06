use crate::{da, r, Mass, Ratio};
use std::fmt::Write;

include!("shared/formula.rs");

impl MolecularFormula {
    pub fn monoisotopic_mass(&self) -> Option<Mass> {
        let mut mass = da(self.additional_mass);
        for (e, i, n) in &self.elements {
            mass += e.mass(*i)? * Ratio::new::<r>(f64::from(*n));
        }
        Some(mass)
    }
    pub fn average_weight(&self) -> Option<Mass> {
        let mut mass = da(self.additional_mass); // Technically this is wrong, the additional mass is defined to be monoisotopic
        for (e, i, n) in &self.elements {
            mass += e.average_weight(*i)? * Ratio::new::<r>(f64::from(*n));
        }
        Some(mass)
    }

    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from this collections of elements
    pub fn hill_notation(&self) -> String {
        // TODO: Handle isotopes
        let mut output = String::new();
        if let Some(carbon) = self.elements.iter().find(|e| e.0 == Element::C) {
            write!(output, "C{}", carbon.2).unwrap();
            if let Some(hydrogen) = self.elements.iter().find(|e| e.0 == Element::H) {
                write!(output, "H{}", hydrogen.2).unwrap();
            }
            for element in self
                .elements
                .iter()
                .filter(|e| e.0 != Element::H && e.0 != Element::C)
            {
                write!(output, "{}{}", element.0, element.2).unwrap();
            }
        } else {
            for element in &self.elements {
                write!(output, "{}{}", element.0, element.2).unwrap();
            }
        }
        output
    }
}

impl std::fmt::Display for MolecularFormula {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.hill_notation())
    }
}

#[cfg(test)]
mod tests {
    use crate::{Element, MolecularFormula};

    #[test]
    fn sorted() {
        assert_eq!(molecular_formula!(H 2 O 2), molecular_formula!(O 2 H 2));
        assert_eq!(
            molecular_formula!(H 6 C 2 O 1),
            molecular_formula!(O 1 C 2 H 6)
        );
        assert_eq!(
            molecular_formula!(H 6 C 2 O 1),
            molecular_formula!(O 1 H 6 C 2)
        );
    }

    #[test]
    fn add() {
        assert_eq!(
            molecular_formula!(H 2 O 2),
            molecular_formula!(H 1 O 1) + molecular_formula!(H 1 O 1)
        );
        assert_eq!(
            molecular_formula!(H 2 O 2),
            molecular_formula!(H 1 O 3) + molecular_formula!(H 1 O -1)
        );
        assert_eq!(
            molecular_formula!(H 2 O 2),
            molecular_formula!(H 1 O -1) + molecular_formula!(H 1 O 3)
        );
        assert_eq!(
            molecular_formula!(H 2 O 2),
            molecular_formula!(H 1 O -1) + molecular_formula!(O 3 H 1)
        );
        assert_eq!(
            molecular_formula!(H 2 O 2),
            molecular_formula!(H 2 O -1) + molecular_formula!(O 3)
        );
    }
}
