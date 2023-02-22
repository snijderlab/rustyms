mod aminoacids;
mod atomic_weights;
mod fragment;
mod model;
mod system;
//mod units;

//use crate::units::*;
use aminoacids::AminoAcid;
use fragment::Fragment;
use model::Model;
use uom::num_traits::Zero;
//use uom::{num_traits::Zero, si::f64::*};
use crate::system::charge::e;
use crate::system::f64::*;

#[macro_use]
extern crate uom;

pub fn generate_theoretical_fragments(
    sequence: &[AminoAcid],
    max_charge: Charge,
    model: &Model,
) -> Vec<Fragment> {
    assert!(max_charge.value >= 1.0);
    assert!(max_charge.value <= u64::MAX as f64);
    let mut output = Vec::new();
    for index in 0..sequence.len() {
        for charge in 1..=(max_charge.value as u64) {
            let n_term = sequence[0..index]
                .iter()
                .fold(Mass::zero(), |acc, aa| acc + aa.avg_mass());
            let c_term = sequence[index + 1..sequence.len()]
                .iter()
                .fold(Mass::zero(), |acc, aa| acc + aa.avg_mass());
            output.append(&mut sequence[index].fragments(
                n_term,
                c_term,
                Charge::new::<e>(charge as f64),
                index,
                model.ions(index, sequence.len()),
            ));
        }
    }
    output
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::aminoacids::AminoAcid;

    #[test]
    fn simple_fragments() {
        let sequence = vec![AminoAcid::W, AminoAcid::F, AminoAcid::W, AminoAcid::F];
        let fragments =
            generate_theoretical_fragments(&sequence, Charge::new::<e>(1.0), &Model::all());
        println!("{}", fragments.len());
        println!("{fragments:?}");
    }
}
