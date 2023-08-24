use crate::{da, Mass};

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
            crate::element::ELEMENTAL_DATA[*self as usize - 1]
                .2
                .iter()
                .find(|(ii, _, _)| *ii == isotope)
                .map(|(_, m, _)| *m)?
        }))
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
