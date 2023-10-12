use crate::Element;
use std::ops::{Add, AddAssign, Mul, Sub};

/// A molecular formula, a selection of elements of specified isotopes together forming a structure
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MolecularFormula {
    /// Save all constituent parts as the element in question, the isotope (or 0 for natural distribution), and the number of this part
    elements: Vec<(crate::Element, u16, i16)>,
    /// Any addition mass, defined to be monoisotopic
    additional_mass: f64,
}

/// Any item that has a clearly defined molecular formula
pub trait Chemical {
    /// Get the molecular formula
    fn formula(&self) -> MolecularFormula;
}

impl<T: Chemical> Chemical for &[T] {
    fn formula(&self) -> MolecularFormula {
        self.iter().map(|i| i.formula()).sum()
    }
}

impl<T: Chemical> Chemical for &Vec<T> {
    fn formula(&self) -> MolecularFormula {
        self.iter().map(|i| i.formula()).sum()
    }
}

impl MolecularFormula {
    /// Create a new molecular formula, the elements will be sorted on element/isotope and deduplicated
    pub fn new(elements: &[(crate::Element, u16, i16)]) -> Self {
        let result = Self {
            elements: elements.to_vec(),
            additional_mass: 0.0,
        };
        result.simplify()
    }

    // The elements will be sorted on element/isotope and deduplicated
    #[must_use]
    fn simplify(mut self) -> Self {
        self.elements.retain(|el| el.2 != 0);
        self.elements.sort_by(|a, b| {
            if a.0 == b.0 {
                // If the elements are the same sort on the isotope number
                a.1.cmp(&b.1)
            } else {
                a.0.cmp(&b.0)
            }
        });
        // Deduplicate
        let mut max = self.elements.len();
        let mut index = 0;
        while index < max {
            let this = self.elements[index];
            let next = self.elements[index + 1];
            if this.0 == next.0 && this.1 == next.1 {
                self.elements[index].2 += next.2;
                self.elements.remove(index + 1);
                max -= 1;
            } else {
                index += 1;
            }
        }
        self.elements.retain(|el| el.2 != 0);
        self
    }

    /// Get an empty molecular formula with only a mass of unspecified origin
    pub const fn with_additional_mass(additional_mass: f64) -> Self {
        Self {
            elements: Vec::new(),
            additional_mass,
        }
    }

    /// Add the given element to this formula (while keeping it ordered and simplified)
    pub fn add(&mut self, element: (crate::Element, u16, i16)) {
        let mut index = 0;
        let mut done = false;
        let (el, i, n) = element;
        while !done {
            let base = self.elements.get(index).copied();
            if let Some((re, ri, _)) = base {
                if el > re || (el == re && i > ri) {
                    index += 1;
                } else if el == re && i == ri {
                    self.elements[index].2 += n;
                    done = true;
                } else {
                    self.elements.insert(index, (el, i, n));
                    done = true;
                }
            } else {
                self.elements.push((el, i, n));
                done = true;
            }
        }
    }

    /// Get the elements making this formula
    pub fn elements(&self) -> &[(Element, u16, i16)] {
        &self.elements
    }

    /// Create a new molecular formula with the given global isotope modifications
    #[must_use]
    pub fn with_global_isotope_modifications(&self, substitutions: &[(Element, u16)]) -> Self {
        let mut new_elements = self.elements.clone();
        for item in &mut new_elements {
            for (substitute_element, substitute_species) in substitutions {
                if item.0 == *substitute_element {
                    item.1 = *substitute_species;
                }
            }
        }
        let result = Self {
            elements: new_elements,
            additional_mass: self.additional_mass,
        };
        result.simplify()
    }

    /// Get the number of electrons (the only charged species, any ionic species is saved as that element +/- the correct number of electrons).
    /// The inverse of that number is given as the charge.
    pub fn charge(&self) -> i16 {
        -self
            .elements
            .iter()
            .find(|el| el.0 == Element::Electron)
            .map_or(0, |el| el.2)
    }
}

impl Add<&MolecularFormula> for &MolecularFormula {
    type Output = MolecularFormula;
    fn add(self, rhs: &MolecularFormula) -> Self::Output {
        let mut result = (*self).clone();
        let mut index_result = 0;
        let mut index_rhs = 0;
        result.additional_mass += rhs.additional_mass;

        while index_rhs < rhs.elements.len() {
            let (el, i, n) = rhs.elements[index_rhs];
            if index_result < result.elements.len() {
                let (re, ri, _) = result.elements[index_result];
                if el > re || (el == re && i > ri) {
                    index_result += 1;
                } else if el == re && i == ri {
                    result.elements[index_result].2 += n;
                    index_rhs += 1;
                } else {
                    result.elements.insert(index_result, (el, i, n));
                    index_rhs += 1;
                }
            } else {
                result.elements.push((el, i, n));
                index_rhs += 1;
            }
        }
        result.elements.retain(|el| el.2 != 0);
        result
    }
}

impl Sub<&MolecularFormula> for &MolecularFormula {
    type Output = MolecularFormula;
    fn sub(self, rhs: &MolecularFormula) -> Self::Output {
        let mut result = (*self).clone();
        let mut index_result = 0;
        let mut index_rhs = 0;
        result.additional_mass -= rhs.additional_mass;
        while index_rhs < rhs.elements.len() {
            let (el, i, n) = rhs.elements[index_rhs];
            if index_result < result.elements.len() {
                let (re, ri, _) = result.elements[index_result];
                if el > re || (el == re && i > ri) {
                    index_result += 1;
                } else if el == re && i == ri {
                    result.elements[index_result].2 -= n;
                    index_rhs += 1;
                } else {
                    result.elements.insert(index_result, (el, i, -n));
                    index_rhs += 1;
                }
            } else {
                result.elements.push((el, i, -n));
                index_rhs += 1;
            }
        }
        result.elements.retain(|el| el.2 != 0);
        result
    }
}

impl Mul<&i16> for &MolecularFormula {
    type Output = MolecularFormula;
    fn mul(self, rhs: &i16) -> Self::Output {
        MolecularFormula {
            additional_mass: self.additional_mass * f64::from(*rhs),
            elements: self
                .elements
                .iter()
                .copied()
                .map(|part| (part.0, part.1, part.2 * rhs))
                .collect(),
        }
    }
}

impl_binop_ref_cases!(impl Add, add for MolecularFormula, MolecularFormula, MolecularFormula);
impl_binop_ref_cases!(impl Sub, sub for MolecularFormula, MolecularFormula, MolecularFormula);
impl_binop_ref_cases!(impl Mul, mul for MolecularFormula, i16, MolecularFormula);

impl AddAssign<&Self> for MolecularFormula {
    fn add_assign(&mut self, rhs: &Self) {
        let mut index_self = 0;
        let mut index_rhs = 0;
        self.additional_mass += rhs.additional_mass;
        while index_rhs < rhs.elements.len() {
            let (el, i, n) = rhs.elements[index_rhs];
            if index_self < self.elements.len() {
                let (re, ri, _) = self.elements[index_self];
                if el > re || (el == re && i > ri) {
                    index_self += 1;
                } else if el == re && i == ri {
                    self.elements[index_self].2 += n;
                    index_rhs += 1;
                } else {
                    self.elements.insert(index_self, (el, i, n));
                    index_rhs += 1;
                }
            } else {
                self.elements.push((el, i, n));
                index_rhs += 1;
            }
        }
    }
}

impl AddAssign<Self> for MolecularFormula {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl std::iter::Sum<Self> for MolecularFormula {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut res = Self::default();
        iter.for_each(|v| res += v);
        res
    }
}

macro_rules! molecular_formula {
    ($($tail:tt)*) => {
        formula_internal!([$($tail)*] -> [])
    };
}

macro_rules! formula_internal {
    ([$e:ident $n:literal $($tail:tt)*] -> [$($output:tt)*]) => {
        formula_internal!([$($tail)*] -> [$($output)*(Element::$e, 0, $n),])
    };
    ([($i:literal)$e:ident $n:literal $($tail:tt)*] -> [$($output:tt)*]) => {
        formula_internal!([$($tail)*] -> [$($output)*(Element::$e, $i, $n),])
    };
    ([$e:ident $n:expr] -> [$($output:tt)*]) =>{
        formula_internal!([] -> [$($output)*(Element::$e, 0, $n),])
    };
    ([($i:literal)$e:ident $n:expr] -> [$($output:tt)*]) =>{
        formula_internal!([] -> [$($output)*(Element::$e, $i, $n),])
    };
    ([] -> [$($output:tt)*]) =>{
        MolecularFormula::new(&[$($output)*])
    };
}
