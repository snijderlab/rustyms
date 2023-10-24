use crate::{da, r, Mass, MassMode, Ratio};
use std::fmt::Write;

include!("shared/formula.rs");

impl MolecularFormula {
    /// The mass of the molecular formula of this element, if all element species (isotopes) exists
    pub fn monoisotopic_mass(&self) -> Option<Mass> {
        let mut mass = da(self.additional_mass);
        for (e, i, n) in &self.elements {
            mass += e.mass(*i)? * Ratio::new::<r>(f64::from(*n));
        }
        Some(mass)
    }
    /// The average weight of the molecular formula of this element, if all element species (isotopes) exists
    pub fn average_weight(&self) -> Option<Mass> {
        let mut mass = da(self.additional_mass); // Technically this is wrong, the additional mass is defined to be monoisotopic
        for (e, i, n) in &self.elements {
            mass += e.average_weight(*i)? * Ratio::new::<r>(f64::from(*n));
        }
        Some(mass)
    }
    /// The most abundant mass, meaning the isotope that will have the highest intensity
    pub fn most_abundant_mass(&self) -> Option<Mass> {
        let mut mass = da(self.additional_mass); // Technically this is wrong, the additional mass is defined to be monoisotopic
        for (e, i, n) in &self.elements {
            mass += e.most_abundant_mass(*n, *i)?;
        }
        Some(mass)
    }
    /// Get the mass in the given mode
    pub fn mass(&self, mode: MassMode) -> Option<Mass> {
        match mode {
            MassMode::Monoisotopic => self.monoisotopic_mass(),
            MassMode::Average => self.average_weight(),
            MassMode::MostAbundant => self.most_abundant_mass(),
        }
    }

    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from this collections of elements merged with the pro forma notation for specific isotopes
    pub fn hill_notation(&self) -> String {
        let mut output = String::new();
        if let Some(carbon) = self.elements.iter().find(|e| e.0 == Element::C && e.1 == 0) {
            write!(output, "C{}", carbon.2).unwrap();
            if let Some(hydrogen) = self.elements.iter().find(|e| e.0 == Element::H && e.1 == 0) {
                write!(output, "H{}", hydrogen.2).unwrap();
            }
            for element in self
                .elements
                .iter()
                .filter(|e| !((e.0 == Element::H || e.0 == Element::C) && e.1 == 0))
            {
                if element.1 == 0 {
                    write!(output, "{}{}", element.0, element.2).unwrap();
                } else {
                    write!(output, "[{}{}{}]", element.1, element.0, element.2).unwrap();
                }
            }
        } else {
            for element in &self.elements {
                if element.1 == 0 {
                    write!(output, "{}{}", element.0, element.2).unwrap();
                } else {
                    write!(output, "[{}{}{}]", element.1, element.0, element.2).unwrap();
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
        if let Some(carbon) = self.elements.iter().find(|e| e.0 == Element::C && e.1 == 0) {
            write!(output, "C{}", to_subscript_num(carbon.2 as isize)).unwrap();
            if let Some(hydrogen) = self.elements.iter().find(|e| e.0 == Element::H && e.1 == 0) {
                write!(output, "H{}", to_subscript_num(hydrogen.2 as isize)).unwrap();
            }
            for element in self
                .elements
                .iter()
                .filter(|e| !((e.0 == Element::H || e.0 == Element::C) && e.1 == 0))
            {
                if element.1 == 0 {
                    write!(
                        output,
                        "{}{}",
                        element.0,
                        to_subscript_num(element.2 as isize)
                    )
                    .unwrap();
                } else {
                    write!(
                        output,
                        "{}{}{}",
                        to_superscript_num(element.1 as isize),
                        element.0,
                        to_subscript_num(element.2 as isize)
                    )
                    .unwrap();
                }
            }
        } else {
            for element in &self.elements {
                if element.1 == 0 {
                    write!(
                        output,
                        "{}{}",
                        element.0,
                        to_subscript_num(element.2 as isize)
                    )
                    .unwrap();
                } else {
                    write!(
                        output,
                        "{}{}{}",
                        to_superscript_num(element.1 as isize),
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
        if let Some(carbon) = self.elements.iter().find(|e| e.0 == Element::C && e.1 == 0) {
            write!(
                output,
                "C{}",
                if carbon.2 == 1 {
                    String::new()
                } else {
                    format!("<sub>{}</sub>", carbon.2)
                }
            )
            .unwrap();
            if let Some(hydrogen) = self.elements.iter().find(|e| e.0 == Element::H && e.1 == 0) {
                write!(
                    output,
                    "H{}",
                    if hydrogen.2 == 1 {
                        String::new()
                    } else {
                        format!("<sub>{}</sub>", hydrogen.2)
                    }
                )
                .unwrap();
            }
            for element in self
                .elements
                .iter()
                .filter(|e| !((e.0 == Element::H || e.0 == Element::C) && e.1 == 0))
                .filter(|e| e.2 != 0)
            {
                write!(
                    output,
                    "{}{}{}",
                    if element.1 == 0 {
                        String::new()
                    } else {
                        format!("<sup>{}</sup>", element.1)
                    },
                    element.0,
                    if element.2 == 1 {
                        String::new()
                    } else {
                        format!("<sub>{}</sub>", element.2)
                    }
                )
                .unwrap();
            }
        } else {
            for element in self.elements.iter().filter(|e| e.2 != 0) {
                write!(
                    output,
                    "{}{}{}",
                    if element.1 == 0 {
                        String::new()
                    } else {
                        format!("<sup>{}</sup>", element.1)
                    },
                    element.0,
                    if element.2 == 1 {
                        String::new()
                    } else {
                        format!("<sub>{}</sub>", element.2)
                    }
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

fn to_subscript_num(input: isize) -> String {
    let text = input.to_string();
    let mut output = String::new();
    for c in text.as_bytes() {
        if *c == b'-' {
            output.push('\u{208B}');
        } else {
            output.push(char::from_u32(*c as u32 + 0x2080 - 0x30).unwrap());
        }
    }
    output
}

fn to_superscript_num(input: isize) -> String {
    let text = input.to_string();
    let mut output = String::new();
    for c in text.as_bytes() {
        if *c == b'-' {
            output.push('\u{207B}');
        } else if *c == b'1' {
            output.push('\u{00B9}');
        } else if *c == b'2' {
            output.push('\u{00B2}');
        } else if *c == b'3' {
            output.push('\u{00B3}');
        } else {
            output.push(char::from_u32(*c as u32 + 0x2070 - 0x30).unwrap());
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
