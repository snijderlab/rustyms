use std::{cmp::Ordering, fmt::Display, str::FromStr};

use crate::{da, modification::Modification, AminoAcid, LinearPeptide, Mass, SequenceElement};

/// A tolerance around a given mass for searching purposes
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum MassTolerance {
    /// A relative search tolerance in parts per million
    Ppm(f64),
    /// An absolute tolerance defined by a constant offset from the mass (bounds are mass - tolerance, mass + tolerance)
    Absolute(Mass),
}

impl MassTolerance {
    /// Find the bounds around a given mass for this tolerance
    pub fn bounds(&self, mass: Mass) -> (Mass, Mass) {
        match self {
            Self::Ppm(ppm) => (
                da(mass.value * (1.0 - ppm / 1e6)),
                da(mass.value * (1.0 + ppm / 1e6)),
            ),
            Self::Absolute(tolerance) => (mass - *tolerance, mass + *tolerance),
        }
    }

    /// See if these two masses are within this tolerance of each other
    pub fn within(&self, a: Mass, b: Mass) -> bool {
        match self {
            Self::Absolute(tol) => (a.value - b.value).abs() <= tol.value,
            Self::Ppm(ppm) => a.ppm(b) <= *ppm,
        }
    }
}

impl Display for MassTolerance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Absolute(mass) => format!("{} da", mass.value),
                Self::Ppm(ppm) => format!("{ppm} ppm"),
            }
        )
    }
}

impl FromStr for MassTolerance {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let num_str = String::from_utf8(
            s.bytes()
                .take_while(|c| {
                    c.is_ascii_digit()
                        || *c == b'.'
                        || *c == b'-'
                        || *c == b'+'
                        || *c == b'e'
                        || *c == b'E'
                })
                .collect::<Vec<_>>(),
        )
        .map_err(|_| ())?;
        let num = num_str.parse::<f64>().map_err(|_| ())?;
        match s[num_str.len()..].trim() {
            "ppm" => Ok(Self::Ppm(num)),
            "da" => Ok(Self::Absolute(da(num))),
            _ => Err(()),
        }
    }
}

impl TryFrom<&str> for MassTolerance {
    type Error = ();
    fn try_from(value: &str) -> Result<Self, ()> {
        value.parse()
    }
}

/// Find the isobaric sets for the given mass with the given modifications and ppm error.
/// The modifications are placed on any location they are allowed based on the given placement
/// rules, so using any modifications which provide those is advised.
/// # Panics
/// Panics if any of the modifications does not have a defined mass.
pub fn find_isobaric_sets(
    mass: Mass,
    tolerance: MassTolerance,
    modifications: &[Modification],
) -> IsobaricSetIterator {
    let bounds = tolerance.bounds(mass);

    // Create the building blocks
    let mut n_term: Vec<(SequenceElement, f64)> = AA
        .iter()
        .flat_map(|aa| {
            let mut options = vec![SequenceElement::new(*aa, None)];
            options.extend(
                modifications
                    .iter()
                    .filter(|&m| can_be_placed(m, *aa, 0, 1))
                    .map(|m| SequenceElement {
                        aminoacid: *aa,
                        ambiguous: None,
                        modifications: vec![m.clone()],
                        possible_modifications: Vec::new(),
                    }),
            );
            options
        })
        .map(|s| {
            (
                s.clone(),
                s.formula_all().unwrap().monoisotopic_mass().unwrap().value,
            )
        })
        .collect();
    n_term.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let mut center: Vec<(SequenceElement, f64)> = AA
        .iter()
        .flat_map(|aa| {
            let mut options = vec![SequenceElement::new(*aa, None)];
            options.extend(
                modifications
                    .iter()
                    .filter(|&m| can_be_placed(m, *aa, 1, 2))
                    .map(|m| SequenceElement {
                        aminoacid: *aa,
                        ambiguous: None,
                        modifications: vec![m.clone()],
                        possible_modifications: Vec::new(),
                    }),
            );
            options
        })
        .map(|s| {
            (
                s.clone(),
                s.formula_all().unwrap().monoisotopic_mass().unwrap().value,
            )
        })
        .collect();
    center.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let mut c_term: Vec<(SequenceElement, f64)> = AA
        .iter()
        .flat_map(|aa| {
            let mut options = vec![SequenceElement::new(*aa, None)];
            options.extend(
                modifications
                    .iter()
                    .filter(|&m| can_be_placed(m, *aa, 1, 1))
                    .map(|m| SequenceElement {
                        aminoacid: *aa,
                        ambiguous: None,
                        modifications: vec![m.clone()],
                        possible_modifications: Vec::new(),
                    }),
            );
            options
        })
        .map(|s| {
            (
                s.clone(),
                s.formula_all().unwrap().monoisotopic_mass().unwrap().value,
            )
        })
        .collect();
    c_term.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    IsobaricSetIterator::new(n_term, c_term, center, bounds)
}

/// Iteratively generate isobaric sets based on the given settings.
#[derive(Debug)]
pub struct IsobaricSetIterator {
    n_term: Vec<(SequenceElement, f64)>,
    c_term: Vec<(SequenceElement, f64)>,
    center: Vec<(SequenceElement, f64)>,
    sizes: (f64, f64),
    bounds: (Mass, Mass),
    state: (Option<usize>, Option<usize>, Vec<usize>),
}

impl IsobaricSetIterator {
    fn new(
        n_term: Vec<(SequenceElement, f64)>,
        c_term: Vec<(SequenceElement, f64)>,
        center: Vec<(SequenceElement, f64)>,
        bounds: (Mass, Mass),
    ) -> Self {
        let sizes = (center.first().unwrap().1, center.last().unwrap().1);
        let mut iter = Self {
            n_term,
            c_term,
            center,
            sizes,
            bounds,
            state: (None, None, Vec::new()),
        };
        while iter.current_mass() < iter.bounds.0.value - iter.sizes.0 {
            iter.state.2.push(0);
        }
        iter
    }

    fn current_mass(&self) -> f64 {
        let mass = self.state.0.map(|i| self.n_term[i].1).unwrap_or_default()
            + self.state.1.map(|i| self.c_term[i].1).unwrap_or_default()
            + self
                .state
                .2
                .iter()
                .copied()
                .map(|i| self.center[i].1)
                .sum::<f64>();
        //println!("{}\t{}", mass.value, self.peptide());
        mass
    }

    fn mass_fits(&self) -> Ordering {
        let mass = self.current_mass();
        if mass < self.bounds.0.value {
            Ordering::Less
        } else if mass > self.bounds.1.value {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }

    fn peptide(&self) -> LinearPeptide {
        let mut sequence = Vec::with_capacity(
            self.state.2.len()
                + usize::from(self.state.0.is_some())
                + usize::from(self.state.1.is_some()),
        );
        if let Some(n) = self.state.0.map(|i| self.n_term[i].clone()) {
            sequence.push(n.0);
        }
        sequence.extend(
            self.state
                .2
                .iter()
                .copied()
                .map(|i| self.center[i].0.clone()),
        );
        if let Some(c) = self.state.1.map(|i| self.c_term[i].clone()) {
            sequence.push(c.0);
        }
        LinearPeptide {
            global: Vec::new(),
            labile: Vec::new(),
            n_term: None,
            c_term: None,
            sequence,
            ambiguous_modifications: Vec::new(),
            charge_carriers: None,
        }
    }
}

impl Iterator for IsobaricSetIterator {
    type Item = LinearPeptide;
    fn next(&mut self) -> Option<Self::Item> {
        // TODO: no check is done for the N and C terminal options
        while !self.state.2.is_empty() {
            // Do state + 1 at the highest level where this is still possible and check if that one fits the bounds
            // Until every level is full then pop and try with one fewer number of aminoacids
            while !self.state.2.iter().all(|s| *s == self.center.len() - 1) {
                let mut level = self.state.2.len() - 1;
                loop {
                    if self.state.2[level] == self.center.len() - 1 {
                        if level == 0 {
                            break;
                        }
                        level -= 1;
                    } else {
                        // Update this level
                        self.state.2[level] += 1;
                        // Reset the levels above, has to start at minimal at the index of this level to prevent 'rotations' of the set to show up
                        for l in level + 1..self.state.2.len() {
                            self.state.2[l] = self.state.2[level];
                        }
                        match self.mass_fits() {
                            Ordering::Greater => {
                                // If the mass is too great the level below will have the be changed, otherwise it could only be getting heavier with the next iteration(s)
                                if level == 0 {
                                    break;
                                }
                                level -= 1;
                            }
                            Ordering::Equal => {
                                return Some(self.peptide());
                            }
                            Ordering::Less => {
                                // If there a way to reach at least the lower limit by having all the heaviest options selected try and reach them.
                                // Otherwise this will increase this level again next iteration.
                                if self.state.2[0..level]
                                    .iter()
                                    .map(|i| self.center[*i].1)
                                    .sum::<f64>()
                                    + (self.state.2.len() - level) as f64 * self.sizes.1
                                    > self.bounds.0.value
                                {
                                    level = self.state.2.len() - 1;
                                }
                            }
                        }
                    }
                }
            }
            self.state.2.pop();
            // Stop the search when there is no possibility for a fitting answer
            if self.state.2.len() as f64 * self.sizes.1 < self.bounds.0.value {
                return None;
            }
            // Reset the levels to be all 0s again
            for level in 0..self.state.2.len() {
                self.state.2[level] = 0;
            }
        }
        None
    }
}

/// Enforce the placement rules of predefined modifications.
fn can_be_placed(modification: &Modification, aa: AminoAcid, index: usize, length: usize) -> bool {
    if let Modification::Predefined(_, rules, _, _, _) = modification {
        rules.is_empty()
            || rules.iter().any(|rule| {
                rule.is_possible(
                    &SequenceElement {
                        aminoacid: aa,
                        modifications: Vec::new(),
                        possible_modifications: Vec::new(),
                        ambiguous: None,
                    },
                    index,
                    length,
                )
            })
    } else {
        true
    }
}

const AA: &[AminoAcid] = &[
    AminoAcid::Glycine,
    AminoAcid::Alanine,
    AminoAcid::Arginine,
    AminoAcid::Asparagine,
    AminoAcid::AsparticAcid,
    AminoAcid::Cysteine,
    AminoAcid::Glutamine,
    AminoAcid::GlutamicAcid,
    AminoAcid::Histidine,
    AminoAcid::AmbiguousLeucine,
    AminoAcid::Lysine,
    AminoAcid::Methionine,
    AminoAcid::Phenylalanine,
    AminoAcid::Proline,
    AminoAcid::Serine,
    AminoAcid::Threonine,
    AminoAcid::Tryptophan,
    AminoAcid::Tyrosine,
    AminoAcid::Valine,
    AminoAcid::Selenocysteine,
    AminoAcid::Pyrrolysine,
];

#[cfg(test)]
mod tests {
    use crate::ComplexPeptide;

    use super::*;
    #[test]
    fn simple_isobaric_sets() {
        let pep = ComplexPeptide::pro_forma("AG").unwrap().assume_linear();
        let sets: Vec<LinearPeptide> = find_isobaric_sets(
            pep.bare_formula().unwrap().monoisotopic_mass().unwrap(),
            MassTolerance::Ppm(10.0),
            &[],
        )
        .collect();
        assert_eq!(
            &sets,
            &[
                ComplexPeptide::pro_forma("GA").unwrap().assume_linear(),
                ComplexPeptide::pro_forma("Q").unwrap().assume_linear(),
            ]
        );
    }
}
