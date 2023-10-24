use crate::{da, r, Mass, Ratio};

include!("shared/element.rs");

impl Element {
    /// Get all available isotopes (N, mass, abundance)
    pub const fn isotopes(self) -> &'static [(u16, f64, f64)] {
        ELEMENTAL_DATA[self as usize].2
    }

    /// The mass of the specified isotope of this element (if that isotope exists)
    pub fn mass(&self, isotope: u16) -> Option<Mass> {
        if *self == Self::Electron {
            return Some(da(5.485_799_090_65e-4));
        }
        Some(da(if isotope == 0 {
            crate::element::ELEMENTAL_DATA[*self as usize - 1].0?
        } else {
            // Specific isotope do not change anything
            crate::element::ELEMENTAL_DATA[*self as usize - 1]
                .2
                .iter()
                .find(|(ii, _, _)| *ii == isotope)
                .map(|(_, m, _)| *m)?
        }))
    }

    /// The average weight of the specified isotope of this element (if that isotope exists)
    pub fn average_weight(&self, isotope: u16) -> Option<Mass> {
        if *self == Self::Electron {
            return Some(da(5.485_799_090_65e-4));
        }
        Some(da(if isotope == 0 {
            crate::element::ELEMENTAL_DATA[*self as usize - 1].1?
        } else {
            // Specific isotope do not change anything
            crate::element::ELEMENTAL_DATA[*self as usize - 1]
                .2
                .iter()
                .find(|(ii, _, _)| *ii == isotope)
                .map(|(_, m, _)| *m)?
        }))
    }

    /// Gives the most abundant mass based on the number of this isotope
    pub fn most_abundant_mass(&self, n: i16, isotope: u16) -> Option<Mass> {
        if *self == Self::Electron {
            return Some(da(5.485_799_090_65e-4) * Ratio::new::<r>(f64::from(n)));
        }
        Some(
            da(if isotope == 0 {
                // (mass, chance)
                let mut max = None;
                for iso in crate::element::ELEMENTAL_DATA[*self as usize - 1].2 {
                    let chance = iso.2 * f64::from(n);
                    if max.map_or(true, |m: (f64, f64)| chance > m.1) {
                        max = Some((iso.1, chance));
                    }
                }
                max?.0
            } else {
                // Specific isotope do not change anything
                crate::element::ELEMENTAL_DATA[*self as usize - 1]
                    .2
                    .iter()
                    .find(|(ii, _, _)| *ii == isotope)
                    .map(|(_, m, _)| *m)?
            }) * Ratio::new::<r>(f64::from(n)),
        )
    }
}

include!(concat!(env!("OUT_DIR"), "/elements.rs"));

#[cfg(test)]
mod test {
    use crate::{Element, MolecularFormula};

    #[test]
    fn hill_notation() {
        assert_eq!(
            molecular_formula!(C 6 O 5 H 10).hill_notation(),
            "C6H10O5".to_string()
        );
    }
}
