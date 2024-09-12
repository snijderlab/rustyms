#![allow(clippy::missing_panics_doc)]

use crate::{
    align::{align, matrix, AlignType},
    LinearPeptide, SimpleLinear, Tolerance,
};

#[test]
fn overextended_rotation() {
    test_alignment("IVQEVS", "LEVQVES", "1i1I2=2r1=");
}

fn test_alignment(peptide_one: &str, peptide_two: &str, path: &str) {
    let first_peptide = LinearPeptide::pro_forma(peptide_one, None)
        .unwrap()
        .into_simple_linear()
        .unwrap();
    let second_peptide = LinearPeptide::pro_forma(peptide_two, None)
        .unwrap()
        .into_simple_linear()
        .unwrap();
    let alignment = align::<4, SimpleLinear, SimpleLinear>(
        &first_peptide,
        &second_peptide,
        matrix::BLOSUM62,
        Tolerance::new_ppm(10.0),
        AlignType::GLOBAL,
    );
    assert_eq!(
        alignment.short(),
        path,
        "Alignment of {peptide_one} vs {peptide_two} did not go to plan"
    );
}