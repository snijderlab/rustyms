use rustyms::{find_isobaric_sets, ComplexPeptide, LinearPeptide};

fn main() {
    let pep = ComplexPeptide::pro_forma("AG").unwrap().assume_linear();
    let sets: Vec<LinearPeptide> = find_isobaric_sets(
        pep.bare_formula().unwrap().monoisotopic_mass().unwrap(),
        10.0,
        &[],
    )
    .collect();
    assert_eq!(
        &sets,
        &[
            ComplexPeptide::pro_forma("AG").unwrap().assume_linear(),
            ComplexPeptide::pro_forma("Q").unwrap().assume_linear(),
        ]
    );
}
