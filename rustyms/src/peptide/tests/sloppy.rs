use crate::{
    modification::Ontology, LinearPeptide, Modification, SloppyParsingParameters, VerySimple,
};

#[test]
fn sloppy_names() {
    assert_eq!(
        Modification::sloppy_modification_internal("Deamidation (NQ)"),
        Some(Ontology::Unimod.find_name("deamidated", None).unwrap())
    );
}

#[test]
fn sloppy_msfragger() {
    assert_eq!(
        LinearPeptide::<VerySimple>::sloppy_pro_forma(
            "n[211]GC[779]RQSSEEK",
            0..20,
            None,
            SloppyParsingParameters {
                ignore_prefix_lowercase_n: true
            }
        )
        .unwrap(),
        LinearPeptide::pro_forma("[211]-GC[779]RQSSEEK", None)
            .unwrap()
            .very_simple()
            .unwrap()
    );
}

#[test]
fn sloppy_hang() {
    let _ = LinearPeptide::<VerySimple>::sloppy_pro_forma(
        "ffffffff[gln->|yro-glu]SC2N:iTRAQ4pleeeeeB]",
        0..20,
        None,
        SloppyParsingParameters::default(),
    );
}
