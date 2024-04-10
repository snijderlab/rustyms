use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    helper_functions::explain_number_error,
    Element, ELEMENT_PARSE_LIST,
};
use std::{
    hash::Hash,
    num::NonZeroU16,
    ops::{Add, AddAssign, Mul, Neg, Sub},
};

/// A molecular formula, a selection of elements of specified isotopes together forming a structure
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Default, Serialize, Deserialize)]
pub struct MolecularFormula {
    /// Save all constituent parts as the element in question, the isotope (or None for natural distribution), and the number of this part
    /// The elements will be sorted on element/isotope and deduplicated, guaranteed to only contain valid isotopes.
    elements: Vec<(crate::Element, Option<NonZeroU16>, i32)>,
    /// Any addition mass, defined to be monoisotopic
    additional_mass: OrderedFloat<f64>,
}

/// Any item that has a clearly defined single molecular formula
pub trait Chemical {
    /// Get the molecular formula
    fn formula(&self) -> MolecularFormula;
}

impl<T: Chemical> Chemical for &[T] {
    fn formula(&self) -> MolecularFormula {
        self.iter().map(Chemical::formula).sum()
    }
}

impl<T: Chemical> Chemical for &Vec<T> {
    fn formula(&self) -> MolecularFormula {
        self.iter().map(Chemical::formula).sum()
    }
}

impl MolecularFormula {
    /// Create a new molecular formula, if the chosen isotopes are not valid it returns None
    pub fn new(elements: &[(crate::Element, Option<NonZeroU16>, i32)]) -> Option<Self> {
        if elements.iter().any(|e| !e.0.is_valid(e.1)) {
            None
        } else {
            let result = Self {
                elements: elements.to_vec(),
                additional_mass: 0.0.into(),
            };
            Some(result.simplify())
        }
    }

    /// Parse Pro Forma formulas: `[13C2][12C-2]H2N`.
    /// # The specification (copied from Pro Forma v2)
    /// As no widely accepted specification exists for expressing elemental formulas, we have adapted a standard with the following rules (taken from <https://github.com/rfellers/chemForma>):
    /// ## Formula Rule 1
    /// A formula will be composed of pairs of atoms and their corresponding cardinality (two Carbon atoms: C2). Pairs SHOULD be separated by spaces but are not required to be.
    /// Atoms and cardinality SHOULD NOT be. Also, the Hill system for ordering (<https://en.wikipedia.org/wiki/Chemical_formula#Hill_system>) is preferred, but not required.
    /// ```text
    /// Example: C12H20O2 or C12 H20 O2
    /// ```
    /// ## Formula Rule 2
    /// Cardinalities must be positive or negative integer values. Zero is not supported. If a cardinality is not included with an atom, it is assumed to be +1.
    /// ```text
    /// Example: HN-1O2
    /// ```
    /// ## Formula Rule 3
    /// Isotopes will be handled by prefixing the atom with its isotopic number in square brackets. If no isotopes are specified, previous rules apply. If no isotope is specified, then it is
    /// assumed the natural isotopic distribution for a given element applies.
    /// ```text
    /// Example: [13C2][12C-2]H2N
    /// Example: [13C2]C-2H2N
    /// ```
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It can panic if the string contains not UTF8 symbols.
    #[allow(dead_code)]
    pub fn from_pro_forma(value: &str) -> Result<Self, CustomError> {
        let mut index = 0;
        let mut element = None;
        let bytes = value.as_bytes();
        let mut result = Self::default();
        while index < value.len() {
            match bytes[index] {
                b'[' => {
                    index += 1; // Skip the open square bracket
                    let len = bytes
                        .iter()
                        .skip(index)
                        .position(|c| *c == b']')
                        .ok_or_else(|| {
                            CustomError::error(
                                "Invalid Pro Forma molecular formula",
                                "No closing square bracket found",
                                Context::line(0, value, index, 1),
                            )
                        })?;
                    let isotope = bytes
                        .iter()
                        .skip(index)
                        .take_while(|c| c.is_ascii_digit())
                        .count();
                    let ele = bytes
                        .iter()
                        .skip(index + isotope)
                        .take_while(|c| c.is_ascii_alphabetic())
                        .count();

                    for possible in ELEMENT_PARSE_LIST {
                        if value[index + isotope..index + isotope + ele] == *possible.0 {
                            element = Some(possible.1);
                            break;
                        }
                    }
                    if let Some(parsed_element) = element {
                        let num = value[index + isotope + ele..index + len]
                            .parse::<i32>()
                            .map_err(|err| {
                                CustomError::error(
                                    "Invalid Pro Forma molecular formula",
                                    format!("The element number {}", explain_number_error(&err)),
                                    Context::line(
                                        0,
                                        value,
                                        index + isotope + ele,
                                        len - isotope - ele,
                                    ),
                                )
                            })?;
                        let isotope = value[index..index + isotope]
                            .parse::<NonZeroU16>()
                            .map_err(|err| {
                                CustomError::error(
                                    "Invalid Pro Forma molecular formula",
                                    format!("The isotope number {}", explain_number_error(&err)),
                                    Context::line(0, value, index, isotope),
                                )
                            })?;

                        if !Self::add(&mut result, (parsed_element, Some(isotope), num)) {
                            return Err(CustomError::error(
                                "Invalid Pro Forma molecular formula",
                                format!("Invalid isotope ({isotope}) added for element ({parsed_element})"),
                                Context::line(0, value, index, len),
                            ),);
                        }
                        element = None;
                        index += len + 1;
                    } else {
                        return Err(CustomError::error(
                            "Invalid Pro Forma molecular formula",
                            "Invalid element",
                            Context::line(0, value, index + isotope, ele),
                        ));
                    }
                }
                b'-' | b'0'..=b'9' if element.is_some() => {
                    let (num, len) = std::str::from_utf8(
                        &bytes
                            .iter()
                            .skip(index)
                            .take_while(|c| c.is_ascii_digit() || **c == b'-')
                            .copied()
                            .collect::<Vec<_>>(),
                    )
                    .map_or_else(
                        |e| panic!("Non UTF8 in Pro Forma molecular formula, error: {e}"),
                        |v| {
                            (
                                v.parse::<i32>().map_err(|err| {
                                    CustomError::error(
                                        "Invalid Pro Forma molecular formula",
                                        format!(
                                            "The element number {}",
                                            explain_number_error(&err)
                                        ),
                                        Context::line(0, value, index, v.len()),
                                    )
                                }),
                                v.len(),
                            )
                        },
                    );
                    let num = num?;
                    if num != 0 && !Self::add(&mut result, (element.unwrap(), None, num)) {
                        return Err(CustomError::error(
                            "Invalid Pro Forma molecular formula",
                            format!(
                                "An element without a defined mass ({}) was used",
                                element.unwrap()
                            ),
                            Context::line(0, value, index - 1, 1),
                        ));
                    }
                    element = None;
                    index += len;
                }
                b' ' => index += 1,
                _ => {
                    if let Some(element) = element {
                        if !Self::add(&mut result, (element, None, 1)) {
                            return Err(CustomError::error(
                                "Invalid Pro Forma molecular formula",
                                format!("An element without a defined mass ({element}) was used"),
                                Context::line(0, value, index - 1, 1),
                            ));
                        }
                    }
                    let mut found = false;
                    for possible in ELEMENT_PARSE_LIST {
                        if value[index..].starts_with(possible.0) {
                            element = Some(possible.1);
                            index += possible.0.len();
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        return Err(CustomError::error(
                            "Invalid Pro Forma molecular formula",
                            "Not a valid character in formula",
                            Context::line(0, value, index, 1),
                        ));
                    }
                }
            }
        }
        if let Some(element) = element {
            if !Self::add(&mut result, (element, None, 1)) {
                return Err(CustomError::error(
                    "Invalid Pro Forma molecular formula",
                    format!("An element without a defined mass ({element}) was used"),
                    Context::line(0, value, index - 1, 1),
                ));
            }
        }
        Ok(result)
    }

    /// PSI-MOD: `(12)C -5 (13)C 5 H 1 N 3 O -1 S 9`
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It can panic if the string contains not UTF8 symbols.
    pub fn from_psi_mod(value: &str) -> Result<Self, CustomError> {
        let mut index = 0;
        let mut isotope = None;
        let mut element = None;
        let bytes = value.as_bytes();
        let mut result = Self::default();
        while index < value.len() {
            match bytes[index] {
                b'(' if isotope.is_none() => {
                    let len = bytes
                        .iter()
                        .skip(index)
                        .position(|c| *c == b')')
                        .ok_or_else(|| {
                            CustomError::error(
                                "Invalid PSI-MOD molecular formula",
                                "No closing round bracket found",
                                Context::line(0, value, index, 1),
                            )
                        })?;
                    isotope = Some(
                        value[index + 1..index + len]
                            .parse::<NonZeroU16>()
                            .map_err(|err| {
                                CustomError::error(
                                    "Invalid PSI-MOD molecular formula",
                                    format!("The isotope number {}", explain_number_error(&err)),
                                    Context::line(0, value, index + 1, len),
                                )
                            })?,
                    );
                    index += len + 1;
                }
                b'-' | b'0'..=b'9' if element.is_some() => {
                    let (num, len) = std::str::from_utf8(
                        &bytes
                            .iter()
                            .skip(index)
                            .take_while(|c| c.is_ascii_digit() || **c == b'-')
                            .copied()
                            .collect::<Vec<_>>(),
                    )
                    .map_or_else(
                        |e| panic!("Non UTF8 in PSI-MOD molecular formula, error: {e}"),
                        |v| {
                            (
                                v.parse::<i32>().map_err(|err| {
                                    CustomError::error(
                                        "Invalid PSI-MOD molecular formula",
                                        format!(
                                            "The isotope number {}",
                                            explain_number_error(&err)
                                        ),
                                        Context::line(0, value, index, v.len()),
                                    )
                                }),
                                v.len(),
                            )
                        },
                    );
                    let num = num?;
                    if num != 0 && !Self::add(&mut result, (element.unwrap(), isotope, num)) {
                        return Err(CustomError::error(
                            "Invalid PSI-MOD molecular formula",
                            format!(
                                "An element without a defined mass ({}) was used",
                                element.unwrap()
                            ),
                            Context::line(0, value, index - 1, 1),
                        ));
                    }
                    element = None;
                    isotope = None;
                    index += len;
                }
                b' ' => index += 1,
                _ => {
                    if let Some(element) = element {
                        if !Self::add(&mut result, (element, None, 1)) {
                            return Err(CustomError::error(
                                "Invalid PSI-MOD molecular formula",
                                format!("An element without a defined mass ({element}) was used"),
                                Context::line(0, value, index - 1, 1),
                            ));
                        }
                    }
                    let mut found = false;
                    for possible in ELEMENT_PARSE_LIST {
                        if value[index..].starts_with(possible.0) {
                            element = Some(possible.1);
                            index += possible.0.len();
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        return Err(CustomError::error(
                            "Invalid PSI-MOD molecular formula",
                            "Not a valid character in formula",
                            Context::line(0, value, index, 1),
                        ));
                    }
                }
            }
        }
        if isotope.is_some() || element.is_some() {
            Err(CustomError::error(
                "Invalid PSI-MOD molecular formula",
                "Last element missed a count",
                Context::line(0, value, index, 1),
            ))
        } else {
            Ok(result)
        }
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
        let mut max = self.elements.len().saturating_sub(1);
        let mut index = 0;
        while index < max {
            let this = self.elements[index];
            let next = self.elements[index + 1];
            if this.0 == next.0 && this.1 == next.1 {
                self.elements[index].2 += next.2;
                self.elements.remove(index + 1);
                max = max.saturating_sub(1);
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
            additional_mass: OrderedFloat(additional_mass),
        }
    }

    /// Add the given element to this formula (while keeping it ordered and simplified).
    /// If the isotope for the added element is not valid it returns `false`.
    #[must_use]
    pub fn add(&mut self, element: (crate::Element, Option<NonZeroU16>, i32)) -> bool {
        if element.0.is_valid(element.1) {
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
            true
        } else {
            false
        }
    }

    /// Get the elements making this formula
    pub fn elements(&self) -> &[(Element, Option<NonZeroU16>, i32)] {
        &self.elements
    }

    /// Create a new molecular formula with the given global isotope modifications. If the given isotope is not valid for this element it returns `None`.
    #[must_use]
    pub fn with_global_isotope_modifications(
        &self,
        substitutions: &[(Element, Option<NonZeroU16>)],
    ) -> Option<Self> {
        if substitutions.iter().all(|e| e.0.is_valid(e.1)) {
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
            Some(result.simplify())
        } else {
            None
        }
    }

    /// Get the number of electrons (the only charged species, any ionic species is saved as that element +/- the correct number of electrons).
    /// The inverse of that number is given as the charge.
    pub fn charge(&self) -> crate::system::isize::Charge {
        -self
            .elements
            .iter()
            .find(|el| el.0 == Element::Electron)
            .map_or_else(crate::system::isize::Charge::default, |el| {
                crate::system::isize::Charge::new::<crate::system::charge::e>(el.2 as isize)
            })
    }
}

impl Neg for &MolecularFormula {
    type Output = MolecularFormula;
    fn neg(self) -> Self::Output {
        let mut res = self.clone();
        for element in &mut res.elements {
            element.2 = -element.2;
        }
        res
    }
}

impl Neg for MolecularFormula {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        for element in &mut self.elements {
            element.2 = -element.2;
        }
        self
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

impl Mul<&i32> for &MolecularFormula {
    type Output = MolecularFormula;
    fn mul(self, rhs: &i32) -> Self::Output {
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
impl_binop_ref_cases!(impl Mul, mul for MolecularFormula, i32, MolecularFormula);

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

#[macro_export]
/// Easily define molecular formulas using the following syntax: `<element> <num>` or `(<isotope>)<element> <num>`.
/// ```
/// # use rustyms::*;
/// molecular_formula!(C 12 (13)C 1 H 24);
/// ```
/// # Panics
/// It panics if the defined molecular formula is not valid. A formula is not valid if not existing isotopes are used
/// or if an element is used that does not have a defined molecular weight (does not have natural abundance).
macro_rules! molecular_formula {
    ($($tail:tt)*) => {
        formula_internal!([$($tail)*] -> [])
    };
}

#[doc(hidden)]
#[macro_export]
/// Internal code for the [`molecular_formula`] macro.
macro_rules! formula_internal {
    ([$e:ident $n:literal $($tail:tt)*] -> [$($output:tt)*]) => {
        formula_internal!([$($tail)*] -> [$($output)*($crate::Element::$e, None, $n),])
    };
    ([($i:literal)$e:ident $n:literal $($tail:tt)*] -> [$($output:tt)*]) => {
        formula_internal!([$($tail)*] -> [$($output)*($crate::Element::$e, Some(std::num::NonZeroU16::new($i).unwrap()), $n),])
    };
    ([$e:ident $n:expr] -> [$($output:tt)*]) =>{
        formula_internal!([] -> [$($output)*($crate::Element::$e, None, $n),])
    };
    ([($i:literal)$e:ident $n:expr] -> [$($output:tt)*]) =>{
        formula_internal!([] -> [$($output)*($crate::Element::$e, Some(std::num::NonZeroU16::new($i).unwrap()), $n),])
    };
    ([] -> [$($output:tt)*]) =>{
        $crate::MolecularFormula::new(&[$($output)*]).unwrap()
    };
}
