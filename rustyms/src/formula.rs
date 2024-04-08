use crate::{
    system::{da, fraction, Mass, OrderedMass, Ratio},
    MassMode,
};
use std::fmt::Write;

include!("shared/formula.rs");

impl From<&MolecularFormula> for OrderedMass {
    /// Create an ordered mass from the monoisotopic mass (needed for [`Multi<MolecularFormula>`](crate::Multi))
    fn from(value: &MolecularFormula) -> Self {
        value.monoisotopic_mass().into()
    }
}

impl MolecularFormula {
    /// The mass of the molecular formula of this element, if all element species (isotopes) exists
    #[allow(clippy::missing_panics_doc)]
    pub fn monoisotopic_mass(&self) -> Mass {
        let mut mass = da(*self.additional_mass);
        for (e, i, n) in &self.elements {
            mass += e
                .mass(*i)
                .expect("An invalid molecular formula was created, please report this crash")
                * Ratio::new::<fraction>(f64::from(*n));
        }
        mass
    }

    /// The average weight of the molecular formula of this element, if all element species (isotopes) exists
    #[allow(clippy::missing_panics_doc)]
    pub fn average_weight(&self) -> Mass {
        let mut mass = da(*self.additional_mass); // Technically this is wrong, the additional mass is defined to be monoisotopic
        for (e, i, n) in &self.elements {
            mass += e
                .average_weight(*i)
                .expect("An invalid molecular formula was created, please report this crash")
                * Ratio::new::<fraction>(f64::from(*n));
        }
        mass
    }

    /// The most abundant mass, meaning the isotope that will have the highest intensity.
    /// It uses an averagine model for the isotopes so the mass will not reflect any isotopomer exact mass
    /// but will be in the form of monoisotopic exact mass + n, where n is the integer dalton offset for that isomer.
    #[allow(clippy::missing_panics_doc)]
    pub fn most_abundant_mass(&self) -> Mass {
        let isotopes = self.isotopic_distribution(0.01);
        let max = isotopes
            .iter()
            .enumerate()
            .max_by_key(|s| OrderedFloat(*s.1));
        self.monoisotopic_mass() + da(max.map_or(0, |f| f.0) as f64)
    }

    /// Get the mass in the given mode
    pub fn mass(&self, mode: MassMode) -> Mass {
        match mode {
            MassMode::Monoisotopic => self.monoisotopic_mass(),
            MassMode::Average => self.average_weight(),
            MassMode::MostAbundant => self.most_abundant_mass(),
        }
    }

    /// Check if the formula is empty
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from this collections of elements merged with the pro forma notation for specific isotopes
    pub fn hill_notation(&self) -> String {
        let mut output = String::new();
        if let Some(carbon) = self
            .elements
            .iter()
            .find(|e| e.0 == Element::C && e.1.is_none())
        {
            write!(output, "C{}", carbon.2).unwrap();
            if let Some(hydrogen) = self
                .elements
                .iter()
                .find(|e| e.0 == Element::H && e.1.is_none())
            {
                write!(output, "H{}", hydrogen.2).unwrap();
            }
            for element in self
                .elements
                .iter()
                .filter(|e| !((e.0 == Element::H || e.0 == Element::C) && e.1.is_none()))
            {
                if let Some(isotope) = element.1 {
                    write!(output, "[{}{}{}]", isotope, element.0, element.2).unwrap();
                } else {
                    write!(output, "{}{}", element.0, element.2).unwrap();
                }
            }
        } else {
            for element in &self.elements {
                if let Some(isotope) = element.1 {
                    write!(output, "[{}{}{}]", isotope, element.0, element.2).unwrap();
                } else {
                    write!(output, "{}{}", element.0, element.2).unwrap();
                }
            }
        }
        if self.additional_mass != 0.0 {
            write!(output, "{:+}", self.additional_mass).unwrap();
        }
        output
    }

    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from this collections of
    /// elements merged with the pro forma notation for specific isotopes. Using fancy unicode characters for subscript
    /// and superscript numbers.
    pub fn hill_notation_fancy(&self) -> String {
        let mut output = String::new();
        if let Some(carbon) = self
            .elements
            .iter()
            .find(|e| e.0 == Element::C && e.1.is_none())
        {
            write!(output, "C{}", to_subscript_num(carbon.2 as isize)).unwrap();
            if let Some(hydrogen) = self
                .elements
                .iter()
                .find(|e| e.0 == Element::H && e.1.is_none())
            {
                write!(output, "H{}", to_subscript_num(hydrogen.2 as isize)).unwrap();
            }
            for element in self
                .elements
                .iter()
                .filter(|e| !((e.0 == Element::H || e.0 == Element::C) && e.1.is_none()))
            {
                if let Some(isotope) = element.1 {
                    write!(
                        output,
                        "{}{}{}",
                        to_superscript_num(isotope),
                        element.0,
                        to_subscript_num(element.2 as isize)
                    )
                    .unwrap();
                } else {
                    write!(
                        output,
                        "{}{}",
                        element.0,
                        to_subscript_num(element.2 as isize)
                    )
                    .unwrap();
                }
            }
        } else {
            for element in &self.elements {
                if let Some(isotope) = element.1 {
                    write!(
                        output,
                        "{}{}{}",
                        to_superscript_num(isotope),
                        element.0,
                        to_subscript_num(element.2 as isize)
                    )
                    .unwrap();
                } else {
                    write!(
                        output,
                        "{}{}",
                        element.0,
                        to_subscript_num(element.2 as isize)
                    )
                    .unwrap();
                }
            }
        }
        if self.additional_mass != 0.0 {
            write!(output, "{:+}", self.additional_mass).unwrap();
        }
        output
    }

    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from this collections of elements encoded in HTML
    pub fn hill_notation_html(&self) -> String {
        let mut output = String::new();
        if let Some(carbon) = self
            .elements
            .iter()
            .find(|e| e.0 == Element::C && e.1.is_none())
        {
            write!(output, "C<sub>{}</sub>", carbon.2).unwrap();
            if let Some(hydrogen) = self
                .elements
                .iter()
                .find(|e| e.0 == Element::H && e.1.is_none())
            {
                write!(output, "H<sub>{}</sub>", hydrogen.2).unwrap();
            }
            for element in self
                .elements
                .iter()
                .filter(|e| !((e.0 == Element::H || e.0 == Element::C) && e.1.is_none()))
                .filter(|e| e.2 != 0)
            {
                write!(
                    output,
                    "{}{}<sub>{}</sub>",
                    element
                        .1
                        .map_or_else(String::new, |isotope| format!("<sup>{isotope}</sup>")),
                    element.0,
                    element.2,
                )
                .unwrap();
            }
        } else {
            for element in self.elements.iter().filter(|e| e.2 != 0) {
                write!(
                    output,
                    "{}{}<sub>{}</sub>",
                    element
                        .1
                        .map_or_else(String::new, |isotope| format!("<sup>{isotope}</sup>")),
                    element.0,
                    element.2,
                )
                .unwrap();
            }
        }
        if self.additional_mass != 0.0 {
            write!(output, "{:+}", self.additional_mass).unwrap();
        }
        output
    }
}

#[allow(clippy::missing_panics_doc)] // Cannot panic
fn to_subscript_num(input: isize) -> String {
    let text = input.to_string();
    let mut output = String::new();
    for c in text.as_bytes() {
        if *c == b'-' {
            output.push('\u{208B}');
        } else {
            output.push(char::from_u32(u32::from(*c) + 0x2080 - 0x30).unwrap());
        }
    }
    output
}

#[allow(clippy::missing_panics_doc)] // Cannot panic
fn to_superscript_num(input: u16) -> String {
    let text = input.to_string();
    let mut output = String::new();
    for c in text.as_bytes() {
        // b'-' could be '\u{207B}' but that is useless when using u16
        if *c == b'1' {
            output.push('\u{00B9}');
        } else if *c == b'2' {
            output.push('\u{00B2}');
        } else if *c == b'3' {
            output.push('\u{00B3}');
        } else {
            output.push(char::from_u32(u32::from(*c) + 0x2070 - 0x30).unwrap());
        }
    }
    output
}

impl std::fmt::Display for MolecularFormula {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.hill_notation())
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
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
    fn simplified() {
        assert_eq!(molecular_formula!(H 2 O 1 O 1), molecular_formula!(O 2 H 2));
        assert_eq!(
            molecular_formula!(H 2 O 1 O 1 H 1 H -2 H 0 H -1 H 2),
            molecular_formula!(O 2 H 2)
        );
        assert_eq!(
            molecular_formula!(H 2 Sb 0 O 1 O 1 H 1 H -2 H 0 H -1 N 0 P 0 Na 0 H 2),
            molecular_formula!(O 2 H 2)
        );
    }

    #[test]
    fn add() {
        assert_eq!(
            molecular_formula!(H 2 O 2).unwrap(),
            molecular_formula!(H 1 O 1).unwrap() + molecular_formula!(H 1 O 1).unwrap()
        );
        assert_eq!(
            molecular_formula!(H 2 O 2).unwrap(),
            molecular_formula!(H 1 O 3).unwrap() + molecular_formula!(H 1 O -1).unwrap()
        );
        assert_eq!(
            molecular_formula!(H 2 O 2).unwrap(),
            molecular_formula!(H 1 O -1).unwrap() + molecular_formula!(H 1 O 3).unwrap()
        );
        assert_eq!(
            molecular_formula!(H 2 O 2).unwrap(),
            molecular_formula!(H 1 O -1).unwrap() + molecular_formula!(O 3 H 1).unwrap()
        );
        assert_eq!(
            molecular_formula!(H 2 O 2).unwrap(),
            molecular_formula!(H 2 O -1).unwrap() + molecular_formula!(O 3).unwrap()
        );
    }
}
