//! Tests all examples provided in the Pro Forma specification chapter 8 Appendix III
use crate::*;

macro_rules! parse_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            assert!(Peptide::pro_forma($case).is_ok());
        }
    };
}

parse_test!("EM[Oxidation]EVEES[Phospho]PEK", summary_1_1_01);
parse_test!(
    "EM[R: Methionine sulfone]EVEES[O-phospho-L-serine]PEK",
    summary_1_1_02
);
parse_test!("EMEVTK[X:DSS#XL1]SESPEK", summary_1_1_03); // (see Section 4.2.3)
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
parse_test!("[Phospho]?EM[Oxidation]EVTSESPEK", summary_2_1_01);
parse_test!(
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
parse_test!(
    "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.90)]ES[#s1(0.90)]PEK",
    summary_2_3
);
parse_test!("PROT(EOSFORMS)[+19.0523]ISK", summary_2_4_01);
parse_test!(
    "PROT(EOC[Carbamidomethyl]FORMS)[+19.0523]ISK",
    summary_2_4_02
);
parse_test!("SEQUEN[Formula:C12H20O2]CE", summary_3_1);
parse_test!("SEQUEN[Formula:HN-1O2]CE", summary_3_2);
parse_test!("SEQUEN[Formula:[13C2][12C-2]H2N]CE", summary_3_3);
parse_test!("SEQUEN[Glycan:HexNAc]CE", summary_4);
parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]", summary_5_01);
parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK", summary_5_02);
parse_test!(
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
    summary_5_03
);
parse_test!("ETFGD[MOD:00093#BRANCH]//R[#BRANCH]ATER", summary_5_04);
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
