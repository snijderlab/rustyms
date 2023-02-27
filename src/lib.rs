#![allow(dead_code)]

mod aminoacids;
mod fragment;
mod mass;
mod mgf;
mod model;
mod peptide;
mod spectrum;
mod system;
//mod units;

pub use crate::mass::*;

//use crate::units::*;
use aminoacids::AminoAcid;
use fragment::Fragment;
use model::Model;
use peptide::Peptide;
use uom::num_traits::Zero;
//use uom::{num_traits::Zero, si::f64::*};
use crate::system::charge::e;
use crate::system::f64::*;

#[macro_use]
extern crate uom;

pub fn generate_theoretical_fragments<M: MassSystem>(
    peptide: &Peptide,
    max_charge: Charge,
    model: &Model,
) -> Vec<Fragment> {
    assert!(max_charge.value >= 1.0);
    assert!(max_charge.value <= u64::MAX as f64);
    let mut output = Vec::new();
    for index in 0..peptide.sequence.len() {
        for charge in 1..=(max_charge.value as u64) {
            let n_term = peptide
                .n_term
                .as_ref()
                .map_or(Mass::zero(), |m| m.mass::<M>())
                + peptide.sequence[0..index]
                    .iter()
                    .fold(Mass::zero(), |acc, aa| {
                        acc + aa.0.mass::<M>()
                            + aa.1.as_ref().map_or(Mass::zero(), |m| m.mass::<M>())
                    });
            let c_term = peptide
                .c_term
                .as_ref()
                .map_or(Mass::zero(), |m| m.mass::<M>())
                + peptide.sequence[index + 1..peptide.sequence.len()]
                    .iter()
                    .fold(Mass::zero(), |acc, aa| {
                        acc + aa.0.mass::<M>()
                            + aa.1.as_ref().map_or(Mass::zero(), |m| m.mass::<M>())
                    });
            output.append(&mut peptide.sequence[index].0.fragments::<M>(
                n_term,
                c_term,
                Charge::new::<e>(charge as f64),
                index,
                model.ions(index, peptide.sequence.len()),
            ));
        }
    }
    output
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn simple_fragments() {
        let peptide = peptide::Peptide::pro_forma("WFWF").unwrap();
        let fragments = generate_theoretical_fragments::<AverageWeight>(
            &peptide,
            Charge::new::<e>(1.0),
            &Model::all(),
        );
        println!("{}", fragments.len());
        println!("{fragments:?}");
    }

    #[test]
    fn simple_matching() {
        let model = Model::all();
        let spectrum = mgf::open("data/example.mgf").unwrap();
        let peptide = peptide::Peptide::pro_forma("WFWF").unwrap();
        let fragments = generate_theoretical_fragments::<AverageWeight>(
            &peptide,
            Charge::new::<e>(1.0),
            &model,
        );
        let annotated = spectrum[0].annotate(peptide, fragments, &model);
        println!("{annotated:?}")
    }
}
