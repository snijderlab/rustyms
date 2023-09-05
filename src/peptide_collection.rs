use crate::{Charge, Fragment, Model, Peptide};

/// A single pro forma entry, can contain multiple peptides
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PeptideCollection {
    /// The reason for having multiple peptides
    pub kind: CollectionKind,
    /// The peptides
    pub peptides: Vec<Peptide>,
}

/// The reason for having multiple peptides
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum CollectionKind {
    /// A multimeric spectrum, multiple peptides coexist in a single spectrum indicated with '+' in pro forma
    #[default]
    Multimeric,
}

impl PeptideCollection {
    /// Assume there is exactly one peptide in this collection
    /// # Panics
    /// If there are no or multiple peptides.
    pub fn assume_singular(mut self) -> Peptide {
        if self.peptides.len() == 1 {
            self.peptides.pop().unwrap()
        } else if self.peptides.len() > 1 {
            panic!("This collection contains multiple spectra, while a single one is assumed")
        } else {
            panic!("This collection does not contain any spectra, while a single one is assumed")
        }
    }

    /// Generate the theoretical fragments for this peptide collection
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
    ) -> Option<Vec<Fragment>> {
        let mut base = Vec::new();
        for (index, peptide) in self.peptides.iter().enumerate() {
            for fragment in peptide.generate_theoretical_fragments(max_charge, model, index)? {
                let (closest_fragment, ppm) =
                    base.iter_mut().fold((None, f64::INFINITY), |acc, i| {
                        let ppm = fragment.ppm(i).map_or(f64::INFINITY, |p| p.value);
                        if acc.1 > ppm {
                            (Some(i), ppm)
                        } else {
                            acc
                        }
                    });
                if ppm < model.ppm.value {
                    // TODO: is this the best combination limit?
                    closest_fragment.unwrap().add_annotation(fragment.ion[0]);
                } else {
                    base.push(fragment);
                }
            }
        }
        Some(base)
    }
}

#[cfg(test)]
mod test {
    use crate::{e, mz, Location, MassOverCharge};

    use super::*;

    #[test]
    fn dimeric_peptide() {
        // Only generate a single series, easier to reason about
        let test_model = Model::new(
            (Location::SkipN(1), Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            Vec::new(),
            MassOverCharge::new::<mz>(20.0),
        );

        // With two different sequences
        let dimeric = Peptide::pro_forma("AA+CC").unwrap();
        let fragments = dbg!(dimeric
            .generate_theoretical_fragments(Charge::new::<e>(1.0), &test_model)
            .unwrap());
        assert_eq!(fragments.len(), 4); // aA, aC, pAA, pCC

        // With two identical sequences
        let dimeric = Peptide::pro_forma("AA+AA").unwrap();
        let fragments = dbg!(dimeric
            .generate_theoretical_fragments(Charge::new::<e>(1.0), &test_model)
            .unwrap());
        assert_eq!(fragments.len(), 2); // aA, pAA
    }
}
