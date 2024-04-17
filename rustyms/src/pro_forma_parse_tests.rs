//! Tests all examples provided in the Pro Forma specification chapter 8 Appendix III
use std::num::NonZeroU16;

use crate::compound_peptidoform::{global_modifications, parse_charge_state};
use crate::modification::{GlobalModification, SimpleModification};
use crate::system::usize::Charge;
use crate::*;
use crate::{model::Location, system::da, CompoundPeptidoform};

macro_rules! parse_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = CompoundPeptidoform::pro_forma($case);
            println!("{}", $case);
            assert!(res.is_ok());
            let back = res.as_ref().unwrap().to_string();
            let res_back = CompoundPeptidoform::pro_forma(&back);
            assert_eq!(res, res_back, "{} != {back}", $case);
        }
    };
    (single $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = CompoundPeptidoform::pro_forma($case);
            println!("{}\n{:?}", $case, res);
            assert!(res.is_ok());
        }
    };
}

parse_test!("EM[Oxidation]EVEES[Phospho]PEK", summary_1_1_01);
// parse_test!(
//     "EM[R: Methionine sulfone]EVEES[O-phospho-L-serine]PEK",
//     summary_1_1_02
// );  TODO: Not implemented yet: resid
// parse_test!("EMEVTK[X:DSS#XL1]SESPEK", summary_1_1_03); // (see Section 4.2.3) TODO: Not implemented yet: cross linking
parse_test!("EM[U:Oxidation]EVEES[U:Phospho]PEK", summary_1_1_04);
parse_test!("EM[+15.9949]EVEES[+79.9663]PEK", summary_1_2_01);
parse_test!("EM[U:+15.995]EVEES[U:+79.966]PEK", summary_1_2_02);
parse_test!("EM[U:+15.995]EVEES[Obs:+79.978]PEK", summary_1_2_03);
parse_test!("RTAAX[+367.0537]WT", summary_1_3);
parse_test!(
    "{Glycan:Hex}EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]",
    summary_1_4
);
parse_test!(
    "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK",
    summary_1_5_01
);
parse_test!(
    "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]",
    summary_1_5_02
);
parse_test!(
    "<[S-carboxamidomethyl-L-cysteine]@C>ATPEILTCNSIGCLK",
    summary_1_6_01
);
parse_test!("<[MOD:01090]@C>ATPEILTCNSIGCLK", summary_1_6_02);
parse_test!(single "[Phospho]?EM[Oxidation]EVTSESPEK", summary_2_1_01);
parse_test!(single
    "[Phospho][Phospho]?[Acetyl]-EM[Oxidation]EVTSESPEK",
    summary_2_1_02
);
parse_test!(
    "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK",
    summary_2_2_01
);
parse_test!(
    "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK",
    summary_2_2_02
);
// parse_test!(
//     "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.90)]ES[#s1(0.90)]PEK",
//     summary_2_3
// );  TODO: Not implemented yet: global ambiguous with localised positions
parse_test!(single "PROT(EOSFORMS)[+19.0523]ISK", summary_2_4_01);
parse_test!(single
    "PROT(EOC[Carbamidomethyl]FORMS)[+19.0523]ISK",
    summary_2_4_02
);
parse_test!("SEQUEN[Formula:C12H20O2]CE", summary_3_1);
parse_test!("SEQUEN[Formula:HN-1O2]CE", summary_3_2);
parse_test!("SEQUEN[Formula:[13C2][12C-2]H2N]CE", summary_3_3);
parse_test!("SEQUEN[Glycan:HexNAc]CE", summary_4);
// parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]", summary_5_01); TODO: Not implemented yet: cross linking
// parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK", summary_5_02); TODO: Not implemented yet: cross linking
// parse_test!(
//     "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
//     summary_5_03
// ); TODO: Not implemented yet: cross linking
// parse_test!("ETFGD[MOD:00093#BRANCH]//R[#BRANCH]ATER", summary_5_04); TODO: Not implemented yet: cross linking/branching
parse_test!("(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K", summary_6);
parse_test!("ELVIS[Phospho|+79.966331]K", summary_7_01);
parse_test!("ELVIS[Phospho|Obs:+79.978]K", summary_7_02);
parse_test!("ELV[INFO:xxxxx]IS", summary_8_01);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:really awesome]K",
    summary_8_02
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created on 2021-06]K",
    summary_8_03
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created by software Tool1]K",
    summary_8_04
);
parse_test!("<13C>ATPEILTVNSIGQLK", summary_9);
parse_test!("EMEVEESPEK/2", summary_10_1);
parse_test!("EMEVEESPEK+ELVISLIVER", summary_10_2_01);
parse_test!("EMEVEESPEK/2+ELVISLIVER/3", summary_10_2_02);

// Personal tests

#[test]
fn parse_glycan() {
    let glycan = CompoundPeptidoform::pro_forma("A[Glycan:Hex]")
        .unwrap()
        .singular()
        .unwrap();
    let spaces = CompoundPeptidoform::pro_forma("A[Glycan:    Hex    ]")
        .unwrap()
        .singular()
        .unwrap();
    assert_eq!(glycan.sequence.len(), 1);
    assert_eq!(spaces.sequence.len(), 1);
    assert_eq!(glycan, spaces);
    let incorrect = CompoundPeptidoform::pro_forma("A[Glycan:Hec]");
    assert!(incorrect.is_err());
}

#[test]
fn parse_formula() {
    let peptide = LinearPeptide::pro_forma("A[Formula:C6H10O5]")
        .unwrap()
        .linear()
        .unwrap();
    let glycan = LinearPeptide::pro_forma("A[Glycan:Hex]")
        .unwrap()
        .linear()
        .unwrap();
    assert_eq!(peptide.sequence.len(), 1);
    assert_eq!(glycan.sequence.len(), 1);
    assert_eq!(glycan.formulas(), peptide.formulas());
}

#[test]
fn parse_labile() {
    let with = LinearPeptide::pro_forma("{Formula:C6H10O5}A")
        .unwrap()
        .linear()
        .unwrap();
    let without = LinearPeptide::pro_forma("A").unwrap().linear().unwrap();
    assert_eq!(with.sequence.len(), 1);
    assert_eq!(without.sequence.len(), 1);
    assert_eq!(with.formulas(), without.formulas());
    assert_eq!(with.labile[0].to_string(), "Formula:C6H10O5".to_string());
}

#[test]
fn parse_ambiguous_modification() {
    let with = CompoundPeptidoform::pro_forma("A[Phospho#g0]A[#g0]")
        .unwrap()
        .singular()
        .unwrap();
    let without = CompoundPeptidoform::pro_forma("AA")
        .unwrap()
        .singular()
        .unwrap();
    assert_eq!(with.sequence.len(), 2);
    assert_eq!(without.sequence.len(), 2);
    assert_eq!(with.sequence[0].possible_modifications.len(), 1);
    assert_eq!(with.sequence[1].possible_modifications.len(), 1);
    assert!(CompoundPeptidoform::pro_forma("A[#g0]A[#g0]").is_err());
    assert!(CompoundPeptidoform::pro_forma("A[Phospho#g0]A[Phospho#g0]").is_err());
    assert!(CompoundPeptidoform::pro_forma("A[Phospho#g0]A[#g0(0.o1)]").is_err());
    assert_eq!(
        CompoundPeptidoform::pro_forma("A[+12#g0]A[#g0]")
            .unwrap()
            .singular()
            .unwrap()
            .to_string(),
        "A[+12#g0]A[#g0]".to_string()
    );
    assert_eq!(
        CompoundPeptidoform::pro_forma("A[#g0]A[+12#g0]")
            .unwrap()
            .singular()
            .unwrap()
            .to_string(),
        "A[#g0]A[+12#g0]".to_string()
    );
}

#[test]
fn parse_ambiguous_aminoacid() {
    let with = LinearPeptide::pro_forma("(?AA)C(?A)(?A)")
        .unwrap()
        .linear()
        .unwrap();
    let without = LinearPeptide::pro_forma("AACAA").unwrap().linear().unwrap();
    assert_eq!(with.sequence.len(), 5);
    assert_eq!(without.sequence.len(), 5);
    assert!(with.sequence[0].ambiguous.is_some());
    assert!(with.sequence[1].ambiguous.is_some());
    assert_eq!(with.formulas(), without.formulas());
    assert_eq!(with.to_string(), "(?AA)C(?A)(?A)".to_string());
}

#[test]
fn parse_hard_tags() {
    let peptide = LinearPeptide::pro_forma("A[Formula:C6H10O5|INFO:hello world 🦀]")
        .unwrap()
        .linear()
        .unwrap();
    let glycan = LinearPeptide::pro_forma(
        "A[info:you can define a tag multiple times|Glycan:Hex|Formula:C6H10O5]",
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
    let deuterium = LinearPeptide::pro_forma("<D>A").unwrap().linear().unwrap();
    let nitrogen_15 = LinearPeptide::pro_forma("<15N>A")
        .unwrap()
        .linear()
        .unwrap();
    assert_eq!(deuterium.sequence.len(), 1);
    assert_eq!(nitrogen_15.sequence.len(), 1);
    // Formula: A + H2O
    assert_eq!(
        deuterium.formulas(),
        molecular_formula!((2)H 7 C 3 O 2 N 1).into()
    );
    assert_eq!(
        nitrogen_15.formulas(),
        molecular_formula!(H 7 C 3 O 2 (15)N 1).into()
    );
}

#[test]
fn parse_chimeric() {
    let dimeric = CompoundPeptidoform::pro_forma("A+AA").unwrap();
    let trimeric = dbg!(CompoundPeptidoform::pro_forma("A+AA-[+2]+AAA").unwrap());
    assert_eq!(dimeric.peptides().len(), 2);
    assert_eq!(dimeric.peptides()[0].len(), 1);
    assert_eq!(dimeric.peptides()[1].len(), 2);
    assert_eq!(trimeric.peptides().len(), 3);
    assert_eq!(trimeric.peptides()[0].len(), 1);
    assert_eq!(trimeric.peptides()[1].len(), 2);
    assert_eq!(trimeric.peptides()[2].len(), 3);
    assert!(trimeric.peptides()[1].c_term.is_some());
}

#[test]
fn parse_unimod() {
    let peptide = dbg!(CompoundPeptidoform::pro_forma(
        "Q[U:Gln->pyro-Glu]E[Cation:Na]AA"
    ));
    assert!(peptide.is_ok());
}

#[test]
fn dimeric_peptide() {
    // Only generate a single series, easier to reason about
    let test_model = Model {
        a: (Location::SkipN(1), Vec::new()),
        ..Model::none()
    };

    // With two different sequences
    let dimeric = CompoundPeptidoform::pro_forma("AA+CC").unwrap();
    let fragments = dbg!(dimeric
        .generate_theoretical_fragments(Charge::new::<crate::system::charge::e>(1), &test_model));
    assert_eq!(fragments.len(), 4); // aA, aC, pAA, pCC

    // With two identical sequences
    let dimeric = CompoundPeptidoform::pro_forma("AA+AA").unwrap();
    let fragments = dbg!(dimeric
        .generate_theoretical_fragments(Charge::new::<crate::system::charge::e>(1), &test_model));
    assert_eq!(fragments.len(), 4); // aA, pAA (both twice once for each peptide)
}

#[test]
fn parse_adduct_ions_01() {
    let peptide = CompoundPeptidoform::pro_forma("A/2[2Na+]+A").unwrap();
    assert_eq!(peptide.peptides().len(), 2);
    assert_eq!(
        peptide.peptides()[0]
            .charge_carriers
            .clone()
            .unwrap()
            .charge_carriers,
        vec![(2, molecular_formula!(Na 1 Electron -1))]
    );
    assert_eq!(
        peptide.peptides()[0].sequence,
        peptide.peptides()[1].sequence
    );
}

#[test]
fn parse_adduct_ions_02() {
    let peptide = dbg!(CompoundPeptidoform::pro_forma("A-[+1]/2[1Na+,+H+]+[+1]-A").unwrap());
    assert_eq!(peptide.peptides().len(), 2);
    assert_eq!(
        peptide.peptides()[0]
            .charge_carriers
            .clone()
            .unwrap()
            .charge_carriers,
        vec![
            (1, molecular_formula!(Na 1 Electron -1)),
            (1, molecular_formula!(H 1 Electron -1))
        ]
    );
    // Check if the C term mod is applied
    assert_eq!(
        peptide.peptides()[0].sequence[0].formulas_all(&[], &[]),
        peptide.peptides()[1].sequence[0].formulas_all(&[], &[])
    );
    assert_eq!(
        peptide.peptides()[0].get_c_term(),
        peptide.peptides()[1].get_n_term()
    );
    assert!(peptide.peptides()[0].get_n_term() != peptide.peptides()[1].get_c_term());
}

#[test]
fn parse_global_modifications() {
    let parse = |str: &str| global_modifications(str.as_bytes(), 0, str);
    assert_eq!(
        parse("<[+5]@D>"),
        Ok((
            8,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::Anywhere,
                Some(AminoAcid::D),
                Modification::Simple(SimpleModification::Mass(da(5.0).into()))
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
                Modification::Simple(SimpleModification::Mass(da(5.0).into()))
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
                Modification::Simple(SimpleModification::Mass(da(5.0).into()))
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
                Modification::Simple(SimpleModification::Mass(da(5.0).into()))
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
                Modification::Simple(SimpleModification::Mass(da(5.0).into()))
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
        parse_charge_state(str.as_bytes(), 0, str).map(|(len, res)| {
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
            molecular_formula!((15)N 1 Electron -1)
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
    let parse = |str: &str| parse_charge_state(str.as_bytes(), 0, str);
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
