#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct MolecularFormula {
    // Save all constituent parts as the element in question, the isotope (or 0 for natural distribution), and the number of this part
    elements: Vec<(crate::Element, u16, i16)>,
}

pub trait Chemical {
    fn formula(&self) -> MolecularFormula;
}

impl MolecularFormula {
    /// Create a new molecular formula, the elements will be sorted on element/isotope but not deduplicated.
    pub fn new(elements: &[(crate::Element, u16, i16)]) -> Self {
        let mut elements = elements.to_vec();
        elements.sort_by(|a, b| {
            if a.0 == b.0 {
                // If the elements are the same sort on the isotope number
                a.1.cmp(&b.1)
            } else {
                a.0.cmp(&b.0)
            }
        });
        MolecularFormula { elements }
    }

    pub fn add(&mut self, element: (crate::Element, u16, i16)) {
        let mut index = 0;
        let mut done = false;
        let (e, i, n) = element;
        while !done {
            let base = self.elements.get(index).copied();
            if let Some((re, ri, _)) = base {
                if e < re || (e == re && i < ri) {
                    index += 1;
                } else if e == re && i == ri {
                    self.elements[index].2 += n;
                    done = true;
                } else {
                    self.elements.insert(index, (e, i, n));
                    done = true;
                }
            } else {
                self.elements.push((e, i, n));
            }
        }
    }
}

impl std::ops::Add<&MolecularFormula> for &MolecularFormula {
    type Output = MolecularFormula;
    fn add(self, rhs: &MolecularFormula) -> Self::Output {
        let mut result = (*self).clone();
        let mut index_result = 0;
        let mut index_rhs = 0;
        while index_rhs < rhs.elements.len() {
            let (e, i, n) = rhs.elements[index_rhs];
            if index_result < result.elements.len() {
                let (re, ri, _) = result.elements[index_result];
                if e < re || (e == re && i < ri) {
                    index_result += 1;
                } else if e == re && i == ri {
                    result.elements[index_result].2 += n;
                    index_rhs += 1;
                } else {
                    result.elements.insert(index_result, (e, i, n));
                    index_rhs += 1;
                }
            } else {
                result.elements.push((e, i, n));
                index_rhs += 1;
            }
        }
        result
    }
}

impl std::ops::Add<MolecularFormula> for MolecularFormula {
    type Output = MolecularFormula;
    fn add(self, rhs: MolecularFormula) -> Self::Output {
        &self + &rhs
    }
}

impl std::ops::Sub<&MolecularFormula> for &MolecularFormula {
    type Output = MolecularFormula;
    fn sub(self, rhs: &MolecularFormula) -> Self::Output {
        let mut result = (*self).clone();
        let mut index_result = 0;
        let mut index_rhs = 0;
        while index_rhs < rhs.elements.len() {
            let (e, i, n) = rhs.elements[index_rhs];
            if index_result < result.elements.len() {
                let (re, ri, _) = result.elements[index_result];
                if e < re || (e == re && i < ri) {
                    index_result += 1;
                } else if e == re && i == ri {
                    result.elements[index_result].2 -= n;
                    index_rhs += 1;
                } else {
                    result.elements.insert(index_result, (e, i, -n));
                    index_rhs += 1;
                }
            } else {
                result.elements.push((e, i, n));
                index_rhs += 1;
            }
        }
        result
    }
}

impl std::ops::Sub<MolecularFormula> for MolecularFormula {
    type Output = MolecularFormula;
    fn sub(self, rhs: MolecularFormula) -> Self::Output {
        &self - &rhs
    }
}

impl std::ops::Mul<i16> for MolecularFormula {
    type Output = MolecularFormula;
    fn mul(mut self, rhs: i16) -> Self::Output {
        self.elements.iter_mut().for_each(|part| part.2 *= rhs);
        self
    }
}

impl std::ops::AddAssign<&MolecularFormula> for MolecularFormula {
    fn add_assign(&mut self, rhs: &MolecularFormula) {
        let mut index_self = 0;
        let mut index_rhs = 0;
        while index_rhs < rhs.elements.len() {
            let (e, i, n) = rhs.elements[index_rhs];
            if index_self < self.elements.len() {
                let (re, ri, _) = self.elements[index_self];
                if e < re || (e == re && i < ri) {
                    index_self += 1;
                } else if e == re && i == ri {
                    self.elements[index_self].2 += n;
                    index_rhs += 1;
                } else {
                    self.elements.insert(index_self, (e, i, n));
                    index_rhs += 1;
                }
            } else {
                self.elements.push((e, i, n));
                index_rhs += 1;
            }
        }
    }
}

impl std::ops::AddAssign<MolecularFormula> for MolecularFormula {
    fn add_assign(&mut self, rhs: MolecularFormula) {
        *self += &rhs;
    }
}
