#![allow(dead_code)]
#![warn(clippy::all, clippy::pedantic, clippy::nursery)]
#![allow(
    clippy::must_use_candidate,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]

mod aminoacids;
mod element;
mod fragment;
mod glycan;
mod mass;
pub mod mgf;
mod model;
mod peptide;
mod spectrum;
mod system;

pub use crate::mass::*;

pub use crate::fragment::*;
pub use crate::glycan::*;
pub use crate::peptide::*;
pub use crate::spectrum::*;
pub use crate::system::f64::*;
pub use aminoacids::AminoAcid;
pub use fragment::Fragment;
pub use model::Model;
pub use peptide::Peptide;
pub use uom::num_traits::Zero;

#[macro_use]
extern crate uom;

pub fn generate_theoretical_fragments<M: MassSystem>(
    peptide: &Peptide,
    max_charge: Charge,
    model: &Model,
) -> Vec<Fragment> {
    assert!(max_charge.value >= 1.0);
    assert!(max_charge.value <= u64::MAX as f64);
    let mut output = Vec::with_capacity(20 * peptide.sequence.len() + 75); // Empirically derived required size of the buffer (Derived from Hecklib)
    for index in 0..peptide.sequence.len() {
        let n_term = peptide.n_term.mass::<M>()
            + peptide.sequence[0..index]
                .iter()
                .fold(Mass::zero(), |acc, aa| {
                    acc + aa.0.mass::<M>() + aa.1.mass::<M>()
                });
        let c_term = peptide.c_term.mass::<M>()
            + peptide.sequence[index + 1..peptide.sequence.len()]
                .iter()
                .fold(Mass::zero(), |acc, aa| {
                    acc + aa.0.mass::<M>() + aa.1.mass::<M>()
                });
        output.append(&mut peptide.sequence[index].0.fragments::<M>(
            n_term,
            c_term,
            max_charge,
            index,
            &model.ions(index, peptide.sequence.len()),
        ));
    }
    output.push(Fragment::new(
        peptide.mass::<M>(),
        max_charge,
        fragment::FragmentType::precursor,
    ));
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
        let annotated = spectrum[0].annotate(peptide, &fragments, &model);
        println!("{annotated:?}");
    }
}
