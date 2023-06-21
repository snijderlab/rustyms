#![allow(dead_code)]
#![warn(clippy::all, clippy::pedantic, clippy::nursery)]
#![allow(
    clippy::must_use_candidate,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::wildcard_imports,
    clippy::module_name_repetitions,
    clippy::suboptimal_flops
)]

pub mod align;
mod aminoacids;
mod element;
mod formula;
mod fragment;
mod glycan;
mod helper_functions;
mod mass;
mod model;
mod modification;
mod ontologies;
mod peptide;
pub mod rawfile;
mod spectrum;
mod system;

pub use crate::mass::*;

pub use crate::fragment::*;
pub use crate::glycan::*;
pub use crate::model::*;
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

/// Generate the theoretical fragments for the given peptide, with the given maximal charge of the fragments, and the given model.
///
/// # Panics
/// If `max_charge` outside the range `1..=u64::MAX`.
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
                .fold(Mass::zero(), |acc, aa| acc + aa.mass::<M>());
        let c_term = peptide.c_term.mass::<M>()
            + peptide.sequence[index + 1..peptide.sequence.len()]
                .iter()
                .fold(Mass::zero(), |acc, aa| acc + aa.mass::<M>());
        output.append(&mut peptide.sequence[index].aminoacid.fragments::<M>(
            // TODO: does this take the mods on the current position into account, also take possible mods into account
            n_term,
            c_term,
            max_charge,
            index,
            peptide.sequence.len(),
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
        let spectrum = rawfile::mgf::open("data/example.mgf").unwrap();
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
