use crate::{da, modification::Modification, AminoAcid, Mass, Peptide, SequenceElement};

pub fn find_isobaric_sets(
    mass: Mass,
    tolerance_ppm: f64,
    modifications: &[Modification],
) -> IsobaricSetIterator {
    let bounds = (
        mass * (1.0 - tolerance_ppm / 1e6),
        mass * (1.0 + tolerance_ppm / 1e6),
    );

    // Create the building blocks
    let n_term = AA
        .iter()
        .flat_map(|aa| {
            let mut options = vec![SequenceElement::new(*aa, None)];
            options.extend(modifications.iter().filter_map(|m| {
                can_be_placed(m, *aa, 0, 1)
                    .ok()
                    .map(|_| SequenceElement::modified(*aa, None, &[m.clone()]))
            }));
            options
        })
        .collect();
    let center: Vec<SequenceElement> = AA
        .iter()
        .flat_map(|aa| {
            let mut options = vec![SequenceElement::new(*aa, None)];
            options.extend(modifications.iter().filter_map(|m| {
                can_be_placed(m, *aa, 1, 2)
                    .ok()
                    .map(|_| SequenceElement::modified(*aa, None, &[m.clone()]))
            }));
            options
        })
        .collect();
    let c_term = AA
        .iter()
        .flat_map(|aa| {
            let mut options = vec![SequenceElement::new(*aa, None)];
            options.extend(modifications.iter().filter_map(|m| {
                can_be_placed(m, *aa, 1, 1)
                    .ok()
                    .map(|_| SequenceElement::modified(*aa, None, &[m.clone()]))
            }));
            options
        })
        .collect();
    let lightest = center.iter().fold(da(f64::INFINITY), |acc, s| {
        let mass = s
            .formula_all()
            .and_then(|f| f.monoisotopic_mass())
            .unwrap_or_else(|| da(f64::INFINITY));
        if mass < acc {
            mass
        } else {
            acc
        }
    });

    dbg!(IsobaricSetIterator {
        n_term,
        c_term,
        center,
        lightest,
        bounds,
        state: (None, None, Vec::new()),
    })
}

#[derive(Debug)]
pub struct IsobaricSetIterator {
    n_term: Vec<SequenceElement>,
    c_term: Vec<SequenceElement>,
    center: Vec<SequenceElement>,
    lightest: Mass,
    bounds: (Mass, Mass),
    state: (Option<usize>, Option<usize>, Vec<usize>),
}

impl IsobaricSetIterator {
    fn current_mass(&self) -> Mass {
        let mass = self
            .state
            .0
            .and_then(|i| self.n_term[i].formula_all())
            .and_then(|f| f.monoisotopic_mass())
            .unwrap_or_default()
            + self
                .state
                .1
                .and_then(|i| self.c_term[i].formula_all())
                .and_then(|f| f.monoisotopic_mass())
                .unwrap_or_default()
            + self
                .state
                .2
                .iter()
                .copied()
                .map(|i| {
                    self.center[i]
                        .formula_all()
                        .and_then(|f| f.monoisotopic_mass())
                        .unwrap_or_default()
                })
                .sum();
        //println!("{}\t{}", mass.value, self.peptide());
        mass
    }

    fn mass_fits(&self) -> bool {
        let mass = self.current_mass();
        mass > self.bounds.0 && mass < self.bounds.1
    }

    fn peptide(&self) -> Peptide {
        let mut sequence = Vec::with_capacity(
            self.state.2.len()
                + usize::from(self.state.0.is_some())
                + usize::from(self.state.1.is_some()),
        );
        if let Some(n) = self.state.0.map(|i| self.n_term[i].clone()) {
            sequence.push(n);
        }
        sequence.extend(self.state.2.iter().copied().map(|i| self.center[i].clone()));
        if let Some(c) = self.state.1.map(|i| self.c_term[i].clone()) {
            sequence.push(c);
        }
        Peptide {
            global: Vec::new(),
            labile: Vec::new(),
            n_term: None,
            c_term: None,
            sequence,
            ambiguous_modifications: Vec::new(),
        }
    }

    fn scan(&mut self) -> Option<Peptide> {
        let last = self.state.2.last().copied().unwrap_or(0); // Be sure to not retry combination that where already tried
        while self.current_mass() < self.bounds.0 - self.lightest {
            self.state.2.push(last);
        }
        //println!("Scan added until {} aas", self.state.2.len());
        //dbg!(&self.state.2);
        // See if the naive addition of the first elements worked
        if self.current_mass() <= self.bounds.1 {
            return Some(self.peptide());
        }

        // Now loop over the last elements to see if any SequenceElements fits the mass
        let last = self.state.2.len() - 1;
        let start = self.state.2[last] + 1;
        for n in start..self.center.len() {
            self.state.2[last] = n;
            let mass = self.current_mass();
            if mass < self.bounds.0 - self.lightest {
                return self.scan(); // Too light try to add more AAs
            }

            if mass > self.bounds.0 && mass < self.bounds.1 {
                return Some(self.peptide());
            }
        }
        //println!("Scanned last level");
        None
    }
}

impl Iterator for IsobaricSetIterator {
    type Item = Peptide;
    fn next(&mut self) -> Option<Self::Item> {
        //println!("Whole new element");
        loop {
            // Check the state (a list of selected pieces)
            if let Some(pep) = self.scan() {
                return Some(pep);
            }
            //dbg!(&self.state.2);
            // No match was found do a prune back as many levels as needed and do the scan again
            self.state.2.pop();
            // If we reach rock bottom give up (not sure this works we might need to have tried all options for this level first)
            if self.state.2.is_empty() {
                return None;
            }
            let last = self.state.2.len() - 1;
            self.state.2[last] += 1;
            //println!("Pop");
            while self.state.2[self.state.2.len() - 1] >= self.center.len() {
                self.state.2.pop();
                //println!("Pop one more");
                let last = self.state.2.len() - 1;
                self.state.2[last] += 1;
            }
            // If we reach rock bottom give up (not sure this works we might need to have tried all options for this level first)
            if self.state.2.is_empty() {
                return None;
            }
        }
    }
}

/// Enforce the placement rules of predefined modifications.
fn can_be_placed(
    modification: &Modification,
    aa: AminoAcid,
    index: usize,
    length: usize,
) -> Result<(), String> {
    if let Modification::Predefined(_, rules, _, _) = modification {
        if !rules.is_empty() && !rules.iter().any(|rule| rule.is_possible(aa, index, length)) {
            return Err(format!(
                "Modification {modification} is not allowed on aminoacid {} index {index}",
                aa.char()
            ));
        }
    }
    Ok(())
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
    use super::*;
    #[test]
    fn simple_isobaric_sets() {
        let pep = Peptide::pro_forma("AG").unwrap();
        let sets: Vec<Peptide> = find_isobaric_sets(
            pep.bare_formula().unwrap().monoisotopic_mass().unwrap(),
            10.0,
            &[],
        )
        .collect();
        assert_eq!(
            &sets,
            &[
                Peptide::pro_forma("AG").unwrap(),
                Peptide::pro_forma("Q").unwrap()
            ]
        );
    }
}
