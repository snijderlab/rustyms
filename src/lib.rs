#![allow(dead_code)]
#![warn(clippy::all, clippy::pedantic, clippy::nursery)]
#![allow(
    clippy::must_use_candidate,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::wildcard_imports,
    clippy::module_name_repetitions,
    clippy::suboptimal_flops,
    clippy::too_many_lines
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn simple_fragments() {
        let peptide = peptide::Peptide::pro_forma("WFWF").unwrap();
        let fragments = peptide
            .generate_theoretical_fragments::<AverageWeight>(Charge::new::<e>(1.0), &Model::all());
        println!("{}", fragments.len());
        println!("{fragments:?}");
    }

    #[test]
    fn simple_matching() {
        let model = Model::all();
        let spectrum = rawfile::mgf::open("data/example.mgf").unwrap();
        let peptide = peptide::Peptide::pro_forma("WFWF").unwrap();
        let fragments =
            peptide.generate_theoretical_fragments::<AverageWeight>(Charge::new::<e>(1.0), &model);
        let annotated = spectrum[0].annotate(peptide, &fragments, &model);
        println!("{annotated:?}");
    }
}
