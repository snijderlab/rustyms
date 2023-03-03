use std::fmt::Display;

use crate::{da, HasMass, Mass, MassSystem};
use std::fmt::Write;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Element {
    H,
    C,
    N,
    O,
    F,
    P,
    S,
    Se,
}

pub const ELEMENT_PARSE_LIST: &[(&str, Element)] = &[
    ("H", Element::H),
    ("C", Element::C),
    ("N", Element::N),
    ("O", Element::O),
    ("F", Element::F),
    ("P", Element::P),
    ("Se", Element::Se),
    ("S", Element::S),
];

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
        }
    }
}

impl Display for Element {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::H => "H",
                Self::C => "C",
                Self::N => "N",
                Self::O => "O",
                Self::F => "F",
                Self::P => "P",
                Self::S => "S",
                Self::Se => "Se",
            }
        )
    }
}

impl Element {
    /// Create a Hill notation from this collections of elements: https://en.wikipedia.org/wiki/Chemical_formula#Hill_system
    pub fn hill_notation(elements: &[(Element, isize)]) -> String {
        let mut sorted = elements
            .iter()
            .map(|(e, c)| (e.to_string(), c))
            .collect::<Vec<_>>();
        sorted.sort_unstable_by(|a, b| a.0.cmp(&b.0));
        if let Some(carbon) = sorted.iter().find(|e| e.0 == "C") {
            let mut output = String::new();
            write!(output, "C{}", carbon.1).unwrap();
            if let Some(hydrogen) = sorted.iter().find(|e| e.0 == "H") {
                write!(output, "H{}", hydrogen.1).unwrap();
            }
            for element in sorted.iter().filter(|e| e.0 == "C" || e.0 == "H") {
                write!(output, "{}{}", element.0, element.1).unwrap();
            }
            output
        } else {
            let mut output = String::new();
            for element in sorted.iter() {
                write!(output, "{}{}", element.0, element.1).unwrap();
            }
            output
        }
    }
}
