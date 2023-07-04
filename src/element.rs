use crate::{da, HasMass, Mass, MassSystem};
use std::fmt::Write;

include!("shared/element.rs");

impl Element {
    /// Get all available isotopes (N, mass, abundance)
    pub const fn isotopes(self) -> &'static [(u16, f64, f64)] {
        ELEMENTAL_DATA[self as usize].2
    }

    pub fn mass(&self, isotope: u16) -> Option<Mass> {
        if *self == Element::Electron {
            return Some(da(5.48579909065e-4));
        }
        Some(da(if isotope == 0 {
            crate::element::ELEMENTAL_DATA[*self as usize - 1].0?
        } else {
            crate::element::ELEMENTAL_DATA[*self as usize - 1]
                .2
                .iter()
                .find(|(ii, m, _)| *ii == isotope)
                .map(|(_, m, _)| *m)?
        }))
    }

    pub fn average_weight(&self, isotope: u16) -> Option<Mass> {
        if self == Element::Electron {
            return Some(da(5.48579909065e-4));
        }
        Some(da(if isotope == 0 {
            crate::element::ELEMENTAL_DATA[*self as usize - 1].1?
        } else {
            crate::element::ELEMENTAL_DATA[*self as usize - 1]
                .2
                .iter()
                .find(|(ii, m, _)| *ii == isotope)
                .map(|(_, m, _)| *m)?
        }))
    }
}

include!(concat!(env!("OUT_DIR"), "/elements.rs"));

#[cfg(test)]
mod test {
    use crate::element::Element;

    #[test]
    fn hill_notation() {
        assert_eq!(
            Element::hill_notation(&[(Element::C, 6), (Element::O, 5), (Element::H, 10)]),
            "C6H10O5".to_string()
        );
    }
}
