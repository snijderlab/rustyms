use crate::{da, HasMass, Mass, MassSystem};
use std::fmt::Write;

include!("shared/element.rs");

impl HasMass for Element {
    fn mass<M: MassSystem>(&self) -> Mass {
        match self {
            Self::H => da(M::H),
            Self::C => da(M::C),
            Self::N => da(M::N),
            Self::O => da(M::O),
            Self::F => da(M::F),
            Self::P => da(M::P),
            Self::S => da(M::S),
            Self::Se => da(M::Se),
            _ => panic!("Mass unknown for this element"),
        }
    }
}

impl Element {
    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from this collections of elements
    pub fn hill_notation(elements: &[(Self, isize)]) -> String {
        let mut sorted = elements.iter().map(|(e, c)| (e, c)).collect::<Vec<_>>();
        sorted.sort_unstable_by(|a, b| a.0.cmp(b.0));
        let mut output = String::new();
        if let Some(carbon) = sorted.iter().find(|e| *e.0 == Self::C) {
            write!(output, "C{}", carbon.1).unwrap();
            if let Some(hydrogen) = sorted.iter().find(|e| *e.0 == Self::H) {
                write!(output, "H{}", hydrogen.1).unwrap();
            }
            for element in sorted.iter().filter(|e| *e.0 != Self::H && *e.0 != Self::C) {
                write!(output, "{}{}", element.0, element.1).unwrap();
            }
        } else {
            for element in &sorted {
                write!(output, "{}{}", element.0, element.1).unwrap();
            }
        }
        output
    }

    pub const fn isotopes(self) -> &'static [(f64, f64)] {
        match self {
            Self::H => ISOTOPES[1],
            Self::C => ISOTOPES[6],
            Self::N => ISOTOPES[7],
            Self::O => ISOTOPES[8],
            _ => ISOTOPES[0],
        }
    }
}

#[allow(clippy::unreadable_literal)]
const ISOTOPES: &[&[(f64, f64)]] = &[
    &[],
    &[(1.007825035, 0.999855), (2.014101779, 0.000145)], //H
    &[],                                                 // He
    &[],                                                 // Li
    &[],                                                 // Be
    &[],                                                 // B
    &[(12.0, 0.9894), (13.003354826, 0.0106)],           // C
    &[(14.003074002, 0.996205), (15.000108965, 0.003795)], // N
    &[
        (15.994914630, 0.99757),
        (16.999131220, 0.0003835),
        (17.999160300, 0.0020465000),
    ], // O
];

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
