use std::num::NonZeroU16;

use crate::{
    model::Location,
    modification::{self, ModificationId, SimpleModification},
    peptide::{
        parse::{global_modifications, parse_charge_state},
        GlobalModification,
    },
    placement_rule::{self, PlacementRule},
    system::{da, usize::Charge},
    AminoAcid, CompoundPeptidoform, Element, LinearPeptide, Model, MolecularCharge, MultiChemical,
    SequencePosition,
};

#[test]
fn parse_global_modifications() {
    let parse = |str: &str| global_modifications(str, 0, None);
    assert_eq!(
        parse("<[+5]@D>"),
        Ok((
            8,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::Anywhere,
                Some(AminoAcid::D),
                SimpleModification::Mass(da(5.0).into())
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@d>"),
        Ok((
            8,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::Anywhere,
                Some(AminoAcid::D),
                SimpleModification::Mass(da(5.0).into())
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@N-term:D>"),
        Ok((
            15,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::AnyNTerm,
                Some(AminoAcid::D),
                SimpleModification::Mass(da(5.0).into())
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@n-term:D>"),
        Ok((
            15,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::AnyNTerm,
                Some(AminoAcid::D),
                SimpleModification::Mass(da(5.0).into())
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@C-term:D>"),
        Ok((
            15,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::AnyCTerm,
                Some(AminoAcid::D),
                SimpleModification::Mass(da(5.0).into())
            )]
        ))
    );
    assert_eq!(
        parse("<D>"),
        Ok((
            3,
            vec![GlobalModification::Isotope(Element::H, NonZeroU16::new(2))]
        ))
    );
    assert_eq!(
        parse("<12C>"),
        Ok((
            5,
            vec![GlobalModification::Isotope(Element::C, NonZeroU16::new(12))]
        ))
    );
    assert!(parse("<D").is_err());
    assert!(parse("<[+5]>").is_err());
    assert!(parse("<[+5]@DD>").is_err());
    assert!(parse("<[5+]@D>").is_err());
    assert!(parse("<[+5@D>").is_err());
    assert!(parse("<+5]@D>").is_err());
    assert!(parse("<[+5#g1]@D>").is_err());
    assert!(parse("<[+5#g1>").is_err());
    assert!(parse("<C12>").is_err());
    assert!(parse("<>").is_err());
    assert!(parse("<@>").is_err());
    assert!(parse("<@D,E,R,T>").is_err());
    assert!(parse("<[+5]@D,E,R,Te>").is_err());
    assert!(parse("<[+5]@D,E,R,N-term:OO>").is_err());
}

#[test]
fn charge_state_positive() {
    let parse = |str: &str| {
        parse_charge_state(str, 0).map(|(len, res)| {
            assert_eq!(
                len,
                str.len(),
                "Not full parsed: '{str}', amount parsed: {len} as '{res}'"
            );
            res
        })
    };
    assert_eq!(parse("/1"), Ok(MolecularCharge::proton(1)));
    assert_eq!(parse("/5"), Ok(MolecularCharge::proton(5)));
    assert_eq!(parse("/-5"), Ok(MolecularCharge::proton(-5)));
    assert_eq!(parse("/1[+H+]"), Ok(MolecularCharge::proton(1)));
    assert_eq!(parse("/2[+H+,+H+]"), Ok(MolecularCharge::proton(2)));
    assert_eq!(
        parse("/1[+Na+]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Na 1 Electron -1)
        )]))
    );
    assert_eq!(
        parse("/3[2Na+1,1H1+1]"),
        Ok(MolecularCharge::new(&[
            (2, molecular_formula!(Na 1 Electron -1)),
            (1, molecular_formula!(H 1 Electron -1))
        ]))
    );
    assert_eq!(
        parse("/1[-OH-]"),
        Ok(MolecularCharge::new(&[(
            -1,
            molecular_formula!(O 1 H 1 Electron 1)
        ),]))
    );
    assert_eq!(
        parse("/1[+N1H3+]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(N 1 H 3 Electron -1)
        ),]))
    );
    assert_eq!(
        parse("/1[+[15N1]+]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!([15 N 1] Electron -1)
        ),]))
    );
    assert_eq!(
        parse("/3[+Fe+3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
    assert_eq!(
        parse("/3[+ Fe +3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
    assert_eq!(
        parse("/-1[+e-]"),
        Ok(MolecularCharge::new(
            &[(1, molecular_formula!(Electron 1)),]
        ))
    );
    assert_eq!(parse("/1[+H1e-1+]"), Ok(MolecularCharge::proton(1)));
    assert_eq!(
        parse("/3[+Fe1e0+3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
    assert_eq!(
        parse("/3[+Fe1e-1+3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
    assert_eq!(
        parse("/3[+Fe1e-2+3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
    assert_eq!(
        parse("/3[+Fe1e-3+3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
}

#[test]
fn charge_state_negative() {
    let parse = |str: &str| parse_charge_state(str, 0);
    assert!(parse("/3[+Fe+]").is_err());
    assert!(parse("/3[+Fe]").is_err());
    assert!(parse("/3[+Fe 1]").is_err());
    assert!(parse("/3[+[54Fe1+3]").is_err());
    assert!(parse("/3[+54Fe1]+3]").is_err());
    assert!(parse("/1[1H1-1]").is_err());
    assert!(parse("/1[1H1+1").is_err());
    assert!(parse("/1[1+1]").is_err());
    assert!(parse("/1[H+1]").is_err());
    assert!(parse("/1[1H]").is_err());
    assert!(parse("/1[1H1]").is_err());
    assert!(parse("/ 1 [ 1 H 1]").is_err());
}

fn parse_glycan() {
    let glycan = LinearPeptide::pro_forma("A[Glycan:Hex]", None).unwrap();
    let spaces = LinearPeptide::pro_forma("A[Glycan:    Hex    ]", None).unwrap();
    assert_eq!(glycan.sequence.len(), 1);
    assert_eq!(spaces.sequence.len(), 1);
    assert_eq!(glycan, spaces);
    let incorrect = CompoundPeptidoform::pro_forma("A[Glycan:Hec]", None);
    assert!(incorrect.is_err());
}

#[test]
fn parse_formula() {
    let peptide = LinearPeptide::pro_forma("A[Formula:C6H10O5]", None)
        .unwrap()
        .linear()
        .unwrap();
    let glycan = LinearPeptide::pro_forma("A[Glycan:Hex]", None)
        .unwrap()
        .linear()
        .unwrap();
    assert_eq!(peptide.sequence.len(), 1);
    assert_eq!(glycan.sequence.len(), 1);
    assert_eq!(glycan.formulas(), peptide.formulas());
}

#[test]
fn parse_labile() {
    let with = LinearPeptide::pro_forma("{Formula:C6H10O5}A", None)
        .unwrap()
        .linear()
        .unwrap();
    let without = LinearPeptide::pro_forma("A", None)
        .unwrap()
        .linear()
        .unwrap();
    assert_eq!(with.sequence.len(), 1);
    assert_eq!(without.sequence.len(), 1);
    assert_eq!(with.formulas(), without.formulas());
    assert_eq!(with.labile[0].to_string(), "Formula:C6H10O5".to_string());
}

#[test]
fn parse_ambiguous_modification() {
    let with = LinearPeptide::pro_forma("A[Phospho#g0]A[#g0]", None).unwrap();
    let without = LinearPeptide::pro_forma("AA", None).unwrap();
    assert_eq!(with.sequence.len(), 2);
    assert_eq!(without.sequence.len(), 2);
    assert_eq!(with.sequence[0].possible_modifications.len(), 1);
    assert_eq!(with.sequence[1].possible_modifications.len(), 1);
    assert!(CompoundPeptidoform::pro_forma("A[#g0]A[#g0]", None).is_err());
    assert!(CompoundPeptidoform::pro_forma("A[Phospho#g0]A[Phospho#g0]", None).is_err());
    assert!(CompoundPeptidoform::pro_forma("A[Phospho#g0]A[#g0(0.o1)]", None).is_err());
    assert_eq!(
        LinearPeptide::pro_forma("A[+12#g0]A[#g0]", None)
            .unwrap()
            .to_string(),
        "A[+12#g0]A[#g0]".to_string()
    );
    assert_eq!(
        LinearPeptide::pro_forma("A[#g0]A[+12#g0]", None)
            .unwrap()
            .to_string(),
        "A[#g0]A[+12#g0]".to_string()
    );
}

#[test]
fn parse_ambiguous_aminoacid() {
    let with = LinearPeptide::pro_forma("(?AA)C(?A)(?A)", None)
        .unwrap()
        .linear()
        .unwrap();
    let without = LinearPeptide::pro_forma("AACAA", None)
        .unwrap()
        .linear()
        .unwrap();
    assert_eq!(with.sequence.len(), 5);
    assert_eq!(without.sequence.len(), 5);
    assert!(with.sequence[0].ambiguous.is_some());
    assert!(with.sequence[1].ambiguous.is_some());
    assert_eq!(with.formulas(), without.formulas());
    assert_eq!(with.to_string(), "(?AA)C(?A)(?A)".to_string());
}

#[test]
fn parse_hard_tags() {
    let peptide = LinearPeptide::pro_forma("A[Formula:C6H10O5|INFO:hello world 🦀]", None)
        .unwrap()
        .linear()
        .unwrap();
    let glycan = LinearPeptide::pro_forma(
        "A[info:you can define a tag multiple times|Glycan:Hex|Formula:C6H10O5]",
        None,
    )
    .unwrap()
    .linear()
    .unwrap();
    assert_eq!(peptide.sequence.len(), 1);
    assert_eq!(glycan.sequence.len(), 1);
    assert_eq!(glycan.formulas(), peptide.formulas());
}

#[test]
fn parse_global() {
    let deuterium = LinearPeptide::pro_forma("<D>A", None)
        .unwrap()
        .linear()
        .unwrap();
    let nitrogen_15 = LinearPeptide::pro_forma("<15N>A", None)
        .unwrap()
        .linear()
        .unwrap();
    assert_eq!(deuterium.sequence.len(), 1);
    assert_eq!(nitrogen_15.sequence.len(), 1);
    // Formula: A + H2O
    assert_eq!(
        deuterium.formulas(),
        molecular_formula!([2 H 7] C 3 O 2 N 1).into()
    );
    assert_eq!(
        nitrogen_15.formulas(),
        molecular_formula!(H 7 C 3 O 2 [15 N 1]).into()
    );
}

#[test]
fn parse_chimeric() {
    let dimeric = CompoundPeptidoform::pro_forma("A+AA", None).unwrap();
    let trimeric = dbg!(CompoundPeptidoform::pro_forma("A+AA-[+2]+AAA", None).unwrap());
    assert_eq!(dimeric.peptidoforms().len(), 2);
    assert_eq!(dimeric.peptidoforms()[0].peptides()[0].len(), 1);
    assert_eq!(dimeric.peptidoforms()[1].peptides()[0].len(), 2);
    assert_eq!(trimeric.peptidoforms().len(), 3);
    assert_eq!(trimeric.peptidoforms()[0].peptides()[0].len(), 1);
    assert_eq!(trimeric.peptidoforms()[1].peptides()[0].len(), 2);
    assert_eq!(trimeric.peptidoforms()[2].peptides()[0].len(), 3);
    assert!(trimeric.peptidoforms()[1].peptides()[0].c_term.is_some());
}

#[test]
fn parse_unimod() {
    let peptide = dbg!(CompoundPeptidoform::pro_forma(
        "[U:Gln->pyro-Glu]-QE[Cation:Na]AA",
        None
    ));
    assert!(peptide.is_ok());
}

#[test]
fn parse_custom() {
    let peptide = dbg!(CompoundPeptidoform::pro_forma(
        "A[C:WEEE]",
        Some(&vec![(
            0,
            "weee".to_string(),
            SimpleModification::Database {
                formula: molecular_formula!(U 1),
                specificities: vec![(
                    vec![PlacementRule::AminoAcid(
                        AminoAcid::CANONICAL_AMINO_ACIDS.to_vec(),
                        placement_rule::Position::Anywhere
                    )],
                    Vec::new(),
                    Vec::new()
                )],
                id: ModificationId {
                    ontology: modification::Ontology::Custom,
                    name: "WEEE".to_string(),
                    id: 0,
                    ..ModificationId::default()
                }
            }
        )])
    ));
    assert!(peptide.is_ok());
    assert_eq!(
        peptide.as_ref().unwrap().to_string(),
        "A[Formula:U1|INFO:Custom:WEEE]"
    );
    assert_eq!(
        peptide.unwrap().formulas(SequencePosition::default(), 0),
        molecular_formula!(C 3 H 7 N 1 O 2 U 1).into()
    );
}

#[test]
fn parse_xl_intra() {
    let peptide = CompoundPeptidoform::pro_forma("A[XLMOD:02001#XLTEST]A[#XLTEST]", None).unwrap();
    let singular = peptide
        .singular()
        .expect("Peptide is not a singular peptide");
    //dbg!(&singular.sequence[0].modifications);
    assert_eq!(
        singular.formulas(SequencePosition::default(), 0).to_vec()[0].elements(),
        (AminoAcid::Alanine
            .formulas(SequencePosition::default(), 0)
            .to_vec()
            .pop()
            .unwrap()
            * 2
            + molecular_formula!(C 8 H 10 O 2)
            + molecular_formula!(H 2 O 1))
        .elements()
    );
}

#[test]
fn parse_xl_inter() {
    let peptide =
        CompoundPeptidoform::pro_forma("A[XLMOD:02001#XLTEST]//A[#XLTEST]", None).unwrap();
    let peptidoform = peptide.singular();
    assert!(
        peptidoform.is_some(),
        "Peptide is not a singular peptidoform"
    );
    let peptidoform = peptidoform.unwrap();
    //dbg!(&singular.sequence[0].modifications);
    assert_eq!(
        peptidoform
            .formulas(SequencePosition::default(), 0)
            .to_vec()[0]
            .elements(),
        (AminoAcid::Alanine
            .formulas(SequencePosition::default(), 0)
            .to_vec()
            .pop()
            .unwrap()
            * 2
            + molecular_formula!(C 8 H 10 O 2)
            + molecular_formula!(H 2 O 1) * 2)
            .elements()
    );
}

#[test]
fn dimeric_peptide() {
    // Only generate a single series, easier to reason about
    let test_model = Model {
        a: (Location::SkipN(1), Vec::new()),
        ..Model::none()
    };

    // With two different sequences
    let dimeric = CompoundPeptidoform::pro_forma("AA+CC", None).unwrap();
    let fragments = dbg!(dimeric
        .generate_theoretical_fragments(Charge::new::<crate::system::charge::e>(1), &test_model));
    assert_eq!(fragments.len(), 4); // aA, aC, pAA, pCC

    // With two identical sequences
    let dimeric = CompoundPeptidoform::pro_forma("AA+AA", None).unwrap();
    let fragments = dbg!(dimeric
        .generate_theoretical_fragments(Charge::new::<crate::system::charge::e>(1), &test_model));
    assert_eq!(fragments.len(), 4); // aA, pAA (both twice once for each peptide)
}

#[test]
fn parse_adduct_ions_01() {
    let peptide = CompoundPeptidoform::pro_forma("A/2[2Na+]+A", None).unwrap();
    assert_eq!(peptide.peptidoforms().len(), 2);
    assert_eq!(
        peptide.peptidoforms()[0].peptides()[0]
            .charge_carriers
            .clone()
            .unwrap()
            .charge_carriers,
        vec![(2, molecular_formula!(Na 1 Electron -1))]
    );
    assert_eq!(
        peptide.peptidoforms()[0].peptides()[0].sequence,
        peptide.peptidoforms()[1].peptides()[0].sequence
    );
}
