mod aminoacids;
mod fragment;
mod mass;
mod mgf;
mod model;
mod spectrum;
mod system;
//mod units;

pub use crate::mass::*;

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

pub fn generate_theoretical_fragments<M: MassSystem>(
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
                .fold(Mass::zero(), |acc, aa| acc + aa.mass::<M>());
            let c_term = sequence[index + 1..sequence.len()]
                .iter()
                .fold(Mass::zero(), |acc, aa| acc + aa.mass::<M>());
            output.append(&mut sequence[index].fragments::<M>(
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
        let fragments = generate_theoretical_fragments::<AverageWeight>(
            &sequence,
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
        let sequence = vec![AminoAcid::W, AminoAcid::F, AminoAcid::W, AminoAcid::F];
        let fragments = generate_theoretical_fragments::<AverageWeight>(
            &sequence,
            Charge::new::<e>(1.0),
            &model,
        );
        let annotated = spectrum[0].annotate(sequence, fragments, &model);
        println!("{:?}", annotated)
    }
}
