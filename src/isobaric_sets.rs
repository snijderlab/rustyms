use std::{cmp::Ordering, fmt::Display, str::FromStr};

use itertools::Itertools;

use crate::{
    modification::Modification,
    placement_rule::{PlacementRule, Position},
    system::{da, r, Mass, Ratio},
    AminoAcid, Chemical, LinearPeptide, SequenceElement,
};

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

/// A list of building blocks for a sequence defined by its sequence elements and its mass.
pub type BuildingBlocks = Vec<(SequenceElement, Mass)>;
/// A list of all combinations of terminal modifications and their accompanying amino acid
pub type TerminalBuildingBlocks = Vec<(SequenceElement, Modification, Mass)>;
/// Get the possible building blocks for sequences based on the given modifications.
/// Useful for any automated sequence generation, like isobaric set generation or de novo sequencing.
/// The result is for each location (N term, center, C term) the list of all possible building blocks with its mass, sorted on mass.
/// # Panics
/// Panics if any of the modifications does not have a defined mass.
pub fn building_blocks(
    amino_acids: &[AminoAcid],
    fixed: &[(Modification, Option<PlacementRule>)],
    variable: &[(Modification, Option<PlacementRule>)],
) -> (
    TerminalBuildingBlocks,
    BuildingBlocks,
    TerminalBuildingBlocks,
) {
    /// Enforce the placement rules of predefined modifications.
    fn can_be_placed(
        modification: &Modification,
        seq: &SequenceElement,
        index: usize,
        length: usize,
    ) -> bool {
        if let Modification::Predefined(_, rules, _, _, _) = modification {
            rules.is_empty()
                || rules
                    .iter()
                    .any(|rule| rule.is_possible(seq, index, length))
        } else {
            true
        }
    }
    fn n_term_options(amino_acids: &[AminoAcid], rule: &PlacementRule) -> Vec<AminoAcid> {
        match rule {
            PlacementRule::AminoAcid(aa, Position::AnyNTerm | Position::ProteinNTerm) => {
                amino_acids
                    .iter()
                    .filter(|a| aa.contains(a))
                    .copied()
                    .collect_vec()
            }
            PlacementRule::Terminal(Position::AnyNTerm | Position::ProteinNTerm) => {
                amino_acids.iter().copied().collect_vec()
            }
            _ => Vec::new(),
        }
    }
    fn c_term_options(amino_acids: &[AminoAcid], rule: &PlacementRule) -> Vec<AminoAcid> {
        match rule {
            PlacementRule::AminoAcid(aa, Position::AnyCTerm | Position::ProteinCTerm) => {
                amino_acids
                    .iter()
                    .filter(|a| aa.contains(a))
                    .copied()
                    .collect_vec()
            }
            PlacementRule::Terminal(Position::AnyCTerm | Position::ProteinCTerm) => {
                amino_acids.iter().copied().collect_vec()
            }
            _ => Vec::new(),
        }
    }
    fn generate_terminal(
        position: &impl Fn(&PlacementRule) -> Vec<AminoAcid>,
        fixed: &[(Modification, Option<PlacementRule>)],
        variable: &[(Modification, Option<PlacementRule>)],
    ) -> TerminalBuildingBlocks {
        let mut options = fixed
            .iter()
            .chain(variable.iter())
            .flat_map(|(modification, rule)| {
                rule.as_ref().map_or_else(
                    || {
                        if let Modification::Predefined(_, rules, _, _, _) = modification {
                            rules
                                .iter()
                                .flat_map(position)
                                .unique()
                                .map(|a| (SequenceElement::new(a, None), modification.clone()))
                                .collect_vec()
                        } else {
                            Vec::new()
                        }
                    },
                    |rule| {
                        position(rule)
                            .into_iter()
                            .map(|a| (SequenceElement::new(a, None), modification.clone()))
                            .collect_vec()
                    },
                )
            })
            .map(|(a, m)| {
                let mass = a
                    .formula_all()
                    .unwrap()
                    .monoisotopic_mass()
                    .unwrap_or_default()
                    + m.formula().monoisotopic_mass().unwrap_or_default();
                (a, m, mass)
            })
            .collect_vec();
        options.sort_unstable_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
        options
    }

    let generate = |index| {
        let mut options: Vec<(SequenceElement, Mass)> = amino_acids
            .iter()
            .flat_map(|aa| {
                let mut options = Vec::new();
                options.extend(
                    fixed
                        .iter()
                        .filter(|&m| {
                            m.1.as_ref().map_or_else(
                                || can_be_placed(&m.0, &SequenceElement::new(*aa, None), index, 2),
                                |rule| rule.is_possible(&SequenceElement::new(*aa, None), index, 2),
                            )
                        })
                        .map(|m| SequenceElement {
                            aminoacid: *aa,
                            ambiguous: None,
                            modifications: vec![m.0.clone()],
                            possible_modifications: Vec::new(),
                        }),
                );
                if options.is_empty() {
                    vec![SequenceElement::new(*aa, None)]
                } else {
                    options
                }
            })
            .flat_map(|seq| {
                let mut options = vec![seq.clone()];
                options.extend(
                    variable
                        .iter()
                        .filter(|&m| {
                            m.1.as_ref().map_or_else(
                                || can_be_placed(&m.0, &seq, index, 2),
                                |rule| rule.is_possible(&seq, index, 2),
                            )
                        })
                        .map(|m| {
                            let mut modifications = seq.modifications.clone();
                            modifications.push(m.0.clone());
                            SequenceElement {
                                aminoacid: seq.aminoacid,
                                ambiguous: None,
                                modifications,
                                possible_modifications: Vec::new(),
                            }
                        }),
                );
                options
            })
            .map(|s| {
                (
                    s.clone(),
                    s.formula_all().unwrap().monoisotopic_mass().unwrap(),
                )
            })
            .collect();
        options.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        options
    };

    // Create the building blocks
    (
        generate_terminal(&|rule| n_term_options(amino_acids, rule), fixed, variable),
        generate(1),
        generate_terminal(&|rule| c_term_options(amino_acids, rule), fixed, variable),
    )
}

/// Find the isobaric sets for the given mass with the given modifications and ppm error.
/// The modifications are placed on any location they are allowed based on the given placement
/// rules, so using any modifications which provide those is advised.
/// # Panics
/// Panics if any of the modifications does not have a defined mass. Or if the weight of the
/// base selection is already in the tolerance of the given mass.
pub fn find_isobaric_sets(
    mass: Mass,
    tolerance: MassTolerance,
    amino_acids: &[AminoAcid],
    fixed: &[(Modification, Option<PlacementRule>)],
    variable: &[(Modification, Option<PlacementRule>)],
    base: Option<&LinearPeptide>,
) -> IsobaricSetIterator {
    let bounds = tolerance.bounds(mass);
    let base_mass = base
        .and_then(|b| b.formula().and_then(|b| b.monoisotopic_mass()))
        .unwrap_or_default();
    let bounds = (bounds.0 - base_mass, bounds.1 - base_mass);
    assert!(bounds.0.value > 0.0, "Cannot have a base selection that has a weight within the tolerance of the intended final mass for isobaric search.");
    let (n_term, center, c_term) = building_blocks(amino_acids, fixed, variable);

    IsobaricSetIterator::new(n_term, c_term, center, bounds, base)
}

/// Iteratively generate isobaric sets based on the given settings.
#[derive(Debug)]
pub struct IsobaricSetIterator {
    n_term: Vec<(SequenceElement, Modification, Mass)>,
    c_term: Vec<(SequenceElement, Modification, Mass)>,
    center: Vec<(SequenceElement, Mass)>,
    sizes: (Mass, Mass),
    bounds: (Mass, Mass),
    state: (Option<usize>, Option<usize>, Vec<usize>),
    base: Option<LinearPeptide>,
}

impl IsobaricSetIterator {
    /// `n_term` & `c_term` are the possible combinations of terminal modifications with their valid placements and the full mass of this combo
    /// The base sequence is assumed to be simple, see [`LinearPeptide::assume_simple`].
    fn new(
        n_term: Vec<(SequenceElement, Modification, Mass)>,
        c_term: Vec<(SequenceElement, Modification, Mass)>,
        center: Vec<(SequenceElement, Mass)>,
        bounds: (Mass, Mass),
        base: Option<&LinearPeptide>,
    ) -> Self {
        let sizes = (center.first().unwrap().1, center.last().unwrap().1);
        let mut iter = Self {
            n_term,
            c_term,
            center,
            sizes,
            bounds,
            state: (None, None, Vec::new()),
            base: base.map(|b| b.clone().assume_simple()),
        };
        while iter.current_mass() < iter.bounds.0 - iter.sizes.0 {
            iter.state.2.push(0);
        }
        iter
    }

    fn current_mass(&self) -> Mass {
        self.state.0.map(|i| self.n_term[i].2).unwrap_or_default()
            + self.state.1.map(|i| self.c_term[i].2).unwrap_or_default()
            + self
                .state
                .2
                .iter()
                .copied()
                .map(|i| self.center[i].1)
                .sum::<Mass>()
    }

    fn mass_fits(&self) -> Ordering {
        let mass = self.current_mass();
        if mass < self.bounds.0 {
            Ordering::Less
        } else if mass > self.bounds.1 {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }

    fn peptide(&self) -> LinearPeptide {
        let mut sequence = Vec::with_capacity(
            self.base
                .as_ref()
                .map(LinearPeptide::len)
                .unwrap_or_default()
                + self.state.2.len()
                + usize::from(self.state.0.is_some())
                + usize::from(self.state.1.is_some()),
        );
        if let Some((b, _)) = self
            .base
            .as_ref()
            .and_then(|b| b.n_term.as_ref().map(|n| (b, n)))
        {
            sequence.push(b.sequence[0].clone());
        } else if let Some(n) = self.state.0.map(|i| self.n_term[i].clone()) {
            sequence.push(n.0);
        }
        if let Some(base) = &self.base {
            sequence.extend(base.sequence.iter().cloned());
        }
        sequence.extend(
            self.state
                .2
                .iter()
                .copied()
                .map(|i| self.center[i].0.clone()),
        );
        if let Some((b, _)) = self
            .base
            .as_ref()
            .and_then(|b| b.c_term.as_ref().map(|c| (b, c)))
        {
            sequence.push(b.sequence.last().unwrap().clone());
        } else if let Some(c) = self.state.1.map(|i| self.c_term[i].clone()) {
            sequence.push(c.0);
        }
        LinearPeptide {
            global: Vec::new(),
            labile: Vec::new(),
            n_term: self
                .base
                .as_ref()
                .and_then(|b| b.n_term.clone())
                .or_else(|| self.state.0.map(|i| self.n_term[i].1.clone())),
            c_term: self
                .base
                .as_ref()
                .and_then(|b| b.n_term.clone())
                .or_else(|| self.state.1.map(|i| self.c_term[i].1.clone())),
            sequence,
            ambiguous_modifications: Vec::new(),
            charge_carriers: None,
        }
    }

    /// Reset the state for the center selection
    fn reset_center_state(&mut self) {
        self.state.2.clear();
        while self.current_mass() < self.bounds.0 - self.sizes.0 {
            self.state.2.push(0);
        }
    }
}

impl Iterator for IsobaricSetIterator {
    type Item = LinearPeptide;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // N terminal loop
            loop {
                // C terminal loop

                // Main loop
                while !self.state.2.is_empty() {
                    // Do state + 1 at the highest level where this is still possible and check if that one fits the bounds
                    // Until every level is full then pop and try with one fewer number of amino acids
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
                                            .sum::<Mass>()
                                            + Ratio::new::<r>((self.state.2.len() - level) as f64)
                                                * self.sizes.1
                                            > self.bounds.0
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
                    if self.sizes.1 * Ratio::new::<r>(self.state.2.len() as f64) < self.bounds.0 {
                        break;
                    }
                    // Reset the levels to be all 0s again
                    for level in 0..self.state.2.len() {
                        self.state.2[level] = 0;
                    }
                }

                // Try the next C terminal option
                if let Some(c) = self.state.1 {
                    if c + 1 >= self.c_term.len() {
                        break;
                    }
                    self.state.1 = Some(c + 1);
                } else if self.c_term.is_empty()
                    || self
                        .base
                        .as_ref()
                        .map(|b| b.c_term.is_some())
                        .unwrap_or_default()
                {
                    // If the base piece has a defined C term mod do not try possible C term mods in the isobaric generation
                    break;
                } else {
                    self.state.1 = Some(0);
                }
                self.reset_center_state();
            }
            // Try the next N terminal option
            if let Some(n) = self.state.0 {
                if n + 1 >= self.n_term.len() {
                    break;
                }
                self.state.0 = Some(n + 1);
            } else if self.n_term.is_empty()
                || self
                    .base
                    .as_ref()
                    .map(|b| b.n_term.is_some())
                    .unwrap_or_default()
            {
                // If the base piece has a defined N term mod do not try possible N term mods in the isobaric generation
                break;
            } else {
                self.state.1 = Some(0);
            }
            self.reset_center_state();
        }
        None
    }
}

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
            AminoAcid::UNIQUE_MASS_AMINO_ACIDS,
            &[],
            &[],
            None,
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
