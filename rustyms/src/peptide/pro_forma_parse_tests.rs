//! Tests all examples provided in the ProForma specification chapter 8 Appendix III,
//! some manual examples, some (past) fuzzer crashes, and some manual tests.
#![allow(clippy::missing_panics_doc)]

use crate::{
    model::Location,
    modification::{ModificationId, SimpleModification},
    placement_rule::PlacementRule,
    system::usize::Charge,
    CompoundPeptidoform, *,
};

macro_rules! parse_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = CompoundPeptidoform::pro_forma($case, None);
            let res_upper = CompoundPeptidoform::pro_forma(&$case.to_ascii_uppercase(), None);
            let res_lower = CompoundPeptidoform::pro_forma(&$case.to_ascii_lowercase(), None);
            println!("{}", $case);
            dbg!(&res);
            assert!(res.is_ok());
            assert_eq!(res, res_upper);
            assert_eq!(res, res_lower);
            let back = res.as_ref().unwrap().to_string();
            let res_back = CompoundPeptidoform::pro_forma(&back, None);
            assert_eq!(res, res_back, "{} != {back}", $case);
        }
    };
    (single $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = CompoundPeptidoform::pro_forma($case, None);
            println!("{}\n{:?}", $case, res);
            assert!(res.is_ok());
        }
    };
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = CompoundPeptidoform::pro_forma($case, None);
            println!("{}\n{:?}", $case, res);
            assert!(res.is_err());
        }
    };
}

// POSITIVE EXAMPLES

parse_test!("AA", positive_example_1);
parse_test!("A[+1]", positive_example_2);
parse_test!("AA[+1]", positive_example_3);
parse_test!("A(AAAA)[+1][+1]", positive_example_4);
parse_test!("UWAKJDNLASNOIJPojkjjdakjn[U:Oxidation]", positive_example_5);
parse_test!("[+1]-A[+1]-[+1]", positive_example_6);
parse_test!("AA+AA", positive_example_7);
parse_test!(
    "EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]SK[#XL1]PEK[#XL2]AR",
    positive_example_8
);
parse_test!(
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
    positive_example_9
);
parse_test!("EM[Oxidation]EVEES[Phospho]PEK", positive_example_10);
parse_test!(
    "EM[R:L-methionine sulfone]EVEES[O-phospho-L-serine]PEK",
    positive_example_11
);
parse_test!("EMEVTK[X:DSS#XL1]SESPEK", positive_example_12);
parse_test!("EM[U:Oxidation]EVEES[U:Phospho]PEK", positive_example_13);
parse_test!("EM[+15.9949]EVEES[+79.9663]PEK", positive_example_14);
parse_test!("EM[U:+15.995]EVEES[U:+79.966]PEK", positive_example_15);
parse_test!("EM[U:+15.995]EVEES[Obs:+79.978]PEK", positive_example_16);
parse_test!("RTAAX[+367.0537]WT", positive_example_17);
parse_test!(
    "{Glycan:Hex}EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]",
    positive_example_18
);
parse_test!(
    "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK",
    positive_example_19
);
parse_test!(
    "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]",
    positive_example_20
);
parse_test!(
    "<[S-carboxamidomethyl-L-cysteine]@C>ATPEILTCNSIGCLK",
    positive_example_21
);
parse_test!("<[MOD:01090]@C>ATPEILTCNSIGCLK", positive_example_22);
parse_test!("[Phospho]?EM[Oxidation]EVTSESPEK", positive_example_23);
parse_test!(
    "[Phospho][Phospho]?[Acetyl]-EM[Oxidation]EVTSESPEK",
    positive_example_24
);
parse_test!(
    "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK",
    positive_example_25
);
parse_test!(
    "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK",
    positive_example_26
);
parse_test!(
    "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.90)]ES[#s1(0.90)]PEK",
    positive_example_27
);
parse_test!("PROT(EOSFORMS)[+19.0523]ISK", positive_example_28);
parse_test!(
    "PROT(EOC[Carbamidomethyl]FORMS)[+19.0523]ISK",
    positive_example_29
);
parse_test!("SEQUEN[Formula:C12H20O2]CE", positive_example_30);
parse_test!("SEQUEN[Formula:HN-1O2]CE", positive_example_31);
parse_test!("SEQUEN[Formula:[13C2][12C-2]H2N]CE", positive_example_32);
parse_test!("SEQUEN[Glycan:HexNAc]CE", positive_example_33);
parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]", positive_example_34);
parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK", positive_example_35);
parse_test!(
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
    positive_example_36
);
// parse_test!(
//     "ETFGD[MOD:00093#BRANCH]//R[#BRANCH]ATER",
//     positive_example_37
// ); This modification cannot be placed on N termini (in the database), so invalid
parse_test!(
    "(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K",
    positive_example_38
);
parse_test!("ELVIS[Phospho|+79.966331]K", positive_example_39);
parse_test!("ELVIS[Phospho|Obs:+79.978]K", positive_example_40);
parse_test!("ELV[INFO:xxxxx]IS", positive_example_41);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:really awesome]K",
    positive_example_42
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created on 2021-06]K",
    positive_example_43
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created by software Tool1]K",
    positive_example_44
);
parse_test!("<13C>ATPEILTVNSIGQLK", positive_example_45);
parse_test!("EMEVEESPEK/2", positive_example_46);
parse_test!("EMEVEESPEK+ELVISLIVER", positive_example_47);
parse_test!("EMEVEESPEK/2+ELVISLIVER/3", positive_example_48);
parse_test!(
    "A[X:DSS#XL1]//B[#XL1]+C[X:DSS#XL1]//D[#XL1]",
    positive_example_49
);
parse_test!("<[Carbamidomethyl]@C>ATPEILTCNSIGCLK", positive_example_50);
parse_test!("<[Oxidation]@C,M>MTPEILTCNSIGCLK", positive_example_51);
parse_test!("<[TMT6plex]@K,N-term>ATPEILTCNSIGCLK", positive_example_52);
parse_test!(
    "<[TMT6plex]@K,N-term:A>ATPEILTCNSIGCLK",
    positive_example_53
);
parse_test!(
    "<[TMT6plex]@K,N-term:A,N-term:B>ATPEILTCNSIGCLK",
    positive_example_54
);
parse_test!("EM[Oxidation]EVEES[Phospho]PEK", positive_example_55);
parse_test!(
    "EM[L-methionine sulfoxide]EVEES[O-phospho-L-serine]PEK",
    positive_example_56
);
parse_test!(
    "EM[R:L-methionine sulfone]EVEES[O-phospho-L-serine]PEK",
    positive_example_57
);
parse_test!("EMEVTK[X:DSS#XL1]SESPEK", positive_example_58);
parse_test!("NEEYN[GNO:G59626AS]K", positive_example_59);
parse_test!("NEEYN[G:G59626AS]K", positive_example_60);
parse_test!("EM[U:Oxidation]EVEES[U:Phospho]PEK", positive_example_61);
parse_test!(
    "EM[M:L-methionine sulfoxide]EVEES[M:O-phospho-L-serine]PEK",
    positive_example_62
);
parse_test!(
    "EM[U:Oxidation]EVEES[M:O-phospho-L-serine]PEK",
    positive_example_63
);
parse_test!(
    "EM[Oxidation]EVEES[O-phospho-L-serine]PEK",
    positive_example_64
);
parse_test!(
    "EM[Oxidation]EVE[Cation:Mg[II]]ES[Phospho]PEK",
    positive_example_65
);
parse_test!("EM[MOD:00719]EVEES[MOD:00046]PEK", positive_example_66);
parse_test!("EM[UNIMOD:35]EVEES[UNIMOD:56]PEK", positive_example_67);
parse_test!(
    "EM[RESID:AA0581]EVEES[RESID:AA0037]PEK",
    positive_example_68
);
parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]", positive_example_69);
parse_test!(
    "EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]SK[#XL1]PEK[#XL2]AR",
    positive_example_70
);
parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK", positive_example_71);
parse_test!("EMEVTK[XLMOD:02001]SESPEK", positive_example_72); // TODO: make sure the hydrolysed status of this linker actually applied
parse_test!(
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
    positive_example_73
);
parse_test!(
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[#XL1]SESPEK",
    positive_example_74
);
parse_test!("EVTSEKC[MOD:00034#XL1]LEMSC[#XL1]EFD", positive_example_75);
parse_test!(
    "EVTSEKC[L-cystine (cross-link)#XL1]LEMSC[#XL1]EFD",
    positive_example_76
);
parse_test!(
    "FVNQHLC[MOD:00034#XL1]GSHLVEALYLVC[MOD:00034#XL2]GERGFFYTPKA//GIVEQC[MOD:00034#XL3]C[#XL1]TSIC[#XL3]SLYQLENYC[#XL2]N",
    positive_example_77
);
parse_test!(
    "EVTSEKC[XLMOD:02009#XL1]LEMSC[#XL1]EFD",
    positive_example_79
);
parse_test!(
    "EVTSEKC[X:Disulfide#XL1]LEMSC[#XL1]EFD",
    positive_example_80
);
parse_test!(
    "EVTSEKC[half cystine]LEMSC[half cystine]EFD",
    positive_example_81
);
parse_test!(
    "EVTSEKC[MOD:00798]LEMSC[MOD:00798]EFDEVTSEKC[MOD:00798]LEMSC[MOD:00798]EFD",
    positive_example_82
);
parse_test!("EVTSEKC[UNIMOD:374#XL1]LEMSC[#XL1]EFD", positive_example_83);
parse_test!("EVTSEKC[Dehydro#XL1]LEMSC[#XL1]EFD", positive_example_84);
// parse_test!(
//     "ETFGD[MOD:00093#BRANCH]//R[#BRANCH]ATER",
//     positive_example_85
// ); This modification cannot be placed on N termini (in the database), so invalid
parse_test!(
    "AVTKYTSSK[MOD:00134#BRANCH]//AGKQLEDGRTLSDYNIQKESTLHLVLRLRG-[#BRANCH]",
    positive_example_86
);
parse_test!("NEEYN[GNO:G59626AS]K", positive_example_87);
parse_test!(
    "YPVLN[GNO:G62765YT]VTMPN[GNO:G02815KT]NSNGKFDK",
    positive_example_88
);
parse_test!("EM[+15.9949]EVEES[+79.9663]PEK", positive_example_89);
parse_test!("EM[+15.995]EVEES[-18.01]PEK", positive_example_90);
parse_test!("EM[U:+15.9949]EVEES[U:+79.9663]PEK", positive_example_91);
parse_test!("EM[U:+15.995]EVEES[U:+79.966]PEK", positive_example_92);
parse_test!("EM[U:+15.995]EVEES[Obs:+79.978]PEK", positive_example_93);
parse_test!("EM[U:+15.995]EVEES[Obs:+79.978]PEK", positive_example_94);
parse_test!("RTAAX[+367.0537]WT", positive_example_95);
parse_test!("SEQUEN[Formula:C12H20O2]CE", positive_example_96);
parse_test!("SEQUEN[Formula:[13C2]CH6N]CE", positive_example_97);
parse_test!("SEQUEN[Formula:[13C2][12C-2]H2N]CE", positive_example_98);
parse_test!("SEQUEN[Glycan:HexNAc1Hex2]CE", positive_example_99);
parse_test!(
    "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK",
    positive_example_100
);
parse_test!(
    "[iTRAQ4plex]-EM[U:Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]",
    positive_example_101
);
parse_test!(
    "{Glycan:Hex}EM[U:Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]",
    positive_example_102
);
parse_test!(
    "{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]",
    positive_example_103
);
parse_test!(
    "{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]",
    positive_example_104
);
parse_test!("{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK", positive_example_105);
parse_test!("[Phospho]?EM[Oxidation]EVTSESPEK", positive_example_106);
parse_test!(
    "[Phospho][Phospho]?[Acetyl]-EM[Oxidation]EVTSESPEK",
    positive_example_107
);
parse_test!(
    "[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK",
    positive_example_108
);
parse_test!(
    "[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK",
    positive_example_109
);
parse_test!(
    "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK",
    positive_example_110
);
parse_test!("PRT(ESFRMS)[+19.0523]ISK", positive_example_111);
parse_test!(
    "PRT(EC[Carbamidomethyl]FRMS)[+19.0523]ISK",
    positive_example_112
);
parse_test!(
    "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK",
    positive_example_113
);
parse_test!(
    "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.09)]ES[#s1(0.90)]PEK",
    positive_example_114
);
parse_test!("MPGLVDSNPAPPESQEKKPLK(PCCACPETKKARDACIIEKGEEHCGHLIEAHKECMRALGFKI)[Oxidation][Oxidation][half cystine][half cystine]", positive_example_115);
parse_test!("<13C>ATPEILTVNSIGQLK", positive_example_116);
parse_test!("<15N>ATPEILTVNSIGQLK", positive_example_117);
parse_test!("<D>ATPEILTVNSIGQLK", positive_example_118);
parse_test!("<13C><15N>ATPEILTVNSIGQLK", positive_example_119);
parse_test!(
    "<[S-carboxamidomethyl-L-cysteine]@C>ATPEILTCNSIGCLK",
    positive_example_120
);
parse_test!("<[MOD:01090]@C>ATPEILTCNSIGCLK", positive_example_121);
parse_test!("<[Oxidation]@C,M>MTPEILTCNSIGCLK", positive_example_122);
parse_test!(
    "<[MOD:01090]@C>[Phospho]?EM[Oxidation]EVTSECSPEK",
    positive_example_123
);
parse_test!(
    "<[MOD:01090]@C>[Acetyl]-EM[Oxidation]EVTSECSPEK",
    positive_example_124
);
parse_test!(
    "(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K",
    positive_example_125
);
parse_test!(
    "(?N)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K",
    positive_example_126
);
parse_test!("ELV[INFO:AnyString]IS", positive_example_127);
parse_test!("ELV[info:AnyString]IS", positive_example_128);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered]K",
    positive_example_129
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:really awesome]K",
    positive_example_130
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created on 2021-06]K",
    positive_example_131
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created by software Tool1]K",
    positive_example_132
);
parse_test!("ELVIS[U:Phospho|+79.966331]K", positive_example_133);
parse_test!("ELVIS[U:Phospho|Obs:+79.978]K", positive_example_134);
parse_test!("ELVIS[Phospho|O-phospho-L-serine]K", positive_example_135);
parse_test!("ELVIS[UNIMOD:21|MOD:00046]K", positive_example_136);
parse_test!("ELVIS[UNIMOD:21|Phospho]K", positive_example_137);
parse_test!(
    "ELVIS[Phospho|O-phospho-L-serine|Obs:+79.966]K",
    positive_example_138
);
parse_test!("ELVIS[Obs:+79.966|Phospho|Sulfo]K", positive_example_139);
parse_test!("EMEVEESPEK/2", positive_example_140);
parse_test!("EM[U:Oxidation]EVEES[U:Phospho]PEK/3", positive_example_141);
parse_test!(
    "[U:iTRAQ4plex]-EM[U:Oxidation]EVNES[U:Phospho]PEK[U:iTRAQ4plex]-[U:Methyl]/3",
    positive_example_142
);
parse_test!("EMEVEESPEK/3[+2Na+,+H+]", positive_example_143);
parse_test!("EMEVEESPEK/1[+2Na+,-H+]", positive_example_144);
parse_test!("EMEVEESPEK/-2[2I-]", positive_example_145);
parse_test!("EMEVEESPEK/-1[+e-]", positive_example_146);
parse_test!("EMEVEESPEK/2+ELVISLIVER/3", positive_example_147);
parse_test!("AA(?AA)", positive_example_148);
parse_test!("AA(?AA)AA", positive_example_149);
parse_test!("[dehydro]^3?Q[gln->pyro-glu]SC", positive_example_150);

// NEGATIVE EXAMPLES

parse_test!(ne "<D>A[UNIMODIFY:+2]+<D>A", negative_example_1);
parse_test!(ne "A[+1]-", negative_example_2);
parse_test!(ne "[Acetyl]-[Phospho]^2?EM[Oxidation]EVTSESPEK", negative_example_3);
parse_test!(ne "PRT(EC[Carbamidomethyl]FRMS)[+19.0523]^2ISK", negative_example_4);
parse_test!(ne "P(RT(ESFRMS)[+19.0523]IS)[+19.0523]K", negative_example_5);
parse_test!(ne "MPGLVDSNPAPPESQEKKPLK(PCCACPETKKARDACIIEKGEEHCGHLIEAHKECMRALGFKI)[Oxidation]^2[half cystine][half cystine]", negative_example_6);
parse_test!(ne "ELVIS[Phospho|INFO:newly]discovered]K", negative_example_7);
parse_test!(ne "<[TMT6plex]>AA", negative_example_8);
parse_test!(ne "<[TMT6plex#g1]@A>AA", negative_example_9);
parse_test!(ne "<[TMT6plex#XL1]@A>AA", negative_example_10);
parse_test!(ne "<[TMT6plex#BRANCH]@A>AA", negative_example_11);
parse_test!(ne "{TMT6plex#g1}AA", negative_example_12);
parse_test!(ne "{TMT6plex#XL1}AA", negative_example_13);
parse_test!(ne "{TMT6plex#BRANCH}AA", negative_example_14);
parse_test!(ne "AA//AA", negative_example_15);
parse_test!(ne "AA(?A(A)[+1])AA", negative_example_16);
parse_test!(ne "AA(A(?A))[+1]AA", negative_example_17);
parse_test!(ne "[dehydro]^3?Q[gln->pyro-glu]S", negative_example_18);
parse_test!(ne "(Q[gln->pyro-glu]S)[Dehydro]", negative_example_19);

// FUZZED CRASHES

parse_test!(ne "ESNCe/", fuzz_1);
parse_test!(ne "/", fuzz_2);
parse_test!(ne "[XLMOd:02001#XL1]", fuzz_3);
parse_test!(ne "[XLOOD:02001#X1]", fuzz_4);
parse_test!(ne "<[TMT6plex]@K,N-term>AtPEILTCNSIGCL/", fuzz_5);
parse_test!(ne "-", fuzz_6);
parse_test!(ne "Y(", fuzz_7);
parse_test!(ne "((", fuzz_8);
parse_test!(ne "()", fuzz_9);
parse_test!(ne "///", fuzz_10);
parse_test!(ne "S///////", fuzz_11);
parse_test!(ne "S//////////(", fuzz_12);
parse_test!(ne "[XLMO01#XL1]", fuzz_13);
parse_test!(ne "SrS/////////", fuzz_14);
parse_test!(ne "S/////////0?", fuzz_15);
parse_test!(ne "[XLMOD:07001#XL1]", fuzz_16);
parse_test!(ne "S//////////ETKvXLSES/0/", fuzz_17);
parse_test!(ne "S///////////////7777777777777777777777777777////E", fuzz_18);
parse_test!(ne "S//Q+++++++/", fuzz_19);
parse_test!(ne "-[Phospho]ddddd(", fuzz_20);
parse_test!(ne "KSPEK/3[1]SE1#XX+2Na", fuzz_21);
parse_test!(ne
    "KgPEKME/3[+222222222222222222222222222222Na1#XX1]SE",
    fuzz_22
);
parse_test!(ne "S//QES/0/0+K/3EK/3EK/3[1]f5[+", fuzz_23);
parse_test!(ne "S//+++/", fuzz_24);
parse_test!(ne "SEK[XLMOD:001#XQ1]SESPE/", fuzz_25);
parse_test!(ne "[[ZPX[[[[[[1]k]]]]]]]", fuzz_26);
parse_test!(ne "[Formula:[09F6]^", fuzz_27);
parse_test!(ne "[Formulaaaaaaaaaaaaaaaaaaaaaa:[13C2][#XL1]SESho]", fuzz_28);
parse_test!(ne "[Formula:[13C2][/EMEVTK[#XL1]", fuzz_29);
parse_test!(ne "[Formul[:\\13C8]^", fuzz_30);
parse_test!(ne "[Fula:[13C23C2]", fuzz_31);
parse_test!(ne "[FoFormula:[13C[EMEVEESPEK/3[+2Na+,+H+]", fuzz_32);
parse_test!(ne "[Formula0m//13CC2]", fuzz_33);
parse_test!(ne "[F][XLMOD:02001#XL1]", fuzz_34);
parse_test!(ne "K[1#Phospho]ddKphoK/", fuzz_35);
parse_test!(ne "K[1#Phospho]tdKSPEK/", fuzz_36);
parse_test!(ne "SEQUEN[Formula:[13C2][12C-2]H2N]CENCE//EME//E/", fuzz_37);
parse_test!(ne "<[TMT6plex]@K,N-term>ddGddd//E/////", fuzz_38);
parse_test!(ne
    "<[TMT6plex]@K,N-term>dSEK[XLMOD:020#XXXXXXLLLLLLLLLL1]UENCE///////////////////3///////",
    fuzz_39
);
parse_test!(ne
    "[NCE//EMEVTK[XLMO[NCE//EMEVTK[XLMOD:02001UUUUUUJ[phlspho]",
    fuzz_40
);
parse_test!(ne "[+2EMEVEENa1#Xa+, 2SPEK/3[+2N+H+]", fuzz_41);
parse_test!(ne "[PhosK[phlspho]d2NaE#Xa+, +H+]", fuzz_42);
parse_test!(ne
    "UE+[Formula:}13Cd [12C-2UE+[Formula:}13Cd [1}13Cd [12C-2BE2NH2N]",
    fuzz_43
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]@K,N-terMT6plex]@K,N-term]@K,N-termS//m>S/",
    fuzz_44
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]@K,N-terMT6plex]@<,N-term]@K,N-term>S//m>S/",
    fuzz_45
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]@K,N-terMT6@K,N-term]@K,N-terMT6plexS//m>S/",
    fuzz_46
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]  KK[Phos -tWrMT6]lex]@K,N-term]@K,N-termY//m>S/",
    fuzz_47
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]@Kpho]erMT6plex<@K,N-term]@K,N-term>S//m>S/",
    fuzz_48
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]@K,N-terMT6plex]ex]@K,C-term>S//Q+++++++//?",
    fuzz_49
);
parse_test!(ne "EMEVEESPEK/3[+5555555555555555555555555Na1#Xa+, ]", fuzz_50);
parse_test!(ne "EMEVEESPEK/3[+26255555555555555555[-hospho]dAddd]", fuzz_51);
parse_test!(ne "DMEVEESPEK/3[+2,2555555555555555555555dNa1#Xa+, ]", fuzz_52);
parse_test!(ne "EMEV+[X://Q++[XS//Q+++++++////////R22222/5NXa+, ]", fuzz_53);
parse_test!(ne "[Formula:   CE[Formula:   C2]", fuzz_54);
parse_test!(ne "EK/3[+2Na+,+H+]PAK/3K/", fuzz_55);
parse_test!(ne "SEK[XLMOD:02001#XL1]UENCE//////////////MEV//E///", fuzz_56);
parse_test!(ne
    "EMEVESEQU//Q++[X://Q++[XS//Q+++++++////////R+/////0//R:/201#X2]",
    fuzz_57
);
parse_test!(ne
    "EMEVESEQU//Q++[j://^++[XS//Q+++++++/3[+2NSRa//////////////R+/R:0201#X2]",
    fuzz_58
);
parse_test!(ne "[TMT6plex]", fuzz_59);
parse_test!(ne "[pphlpho]", fuzz_60);
parse_test!(ne "SSSHSSSSSSSM+S+NPS+S+NPSS+?", fuzz_61);
parse_test!(ne "SESPEK/3[+2N1111111111N1111111111										11111a+,2N]CCE", fuzz_62);
parse_test!(ne "?", fuzz_63);
parse_test!(ne "[TMT6plex]?", fuzz_64);
parse_test!(ne "[]?", fuzz_65);
parse_test!(ne "{}{}", fuzz_66);
parse_test!(ne "EMEVEESH+?", fuzz_67);
parse_test!(ne "EM/3A+H+?", fuzz_68);
parse_test!(ne "{}{}{6}", fuzz_69);
parse_test!(ne "EM/3[+SPEKEVVEE2Na+,+VW( EE3Na+,+HPFKAUBE2Nb+,+HPFKAUBE2`{b,++]", fuzz_70);
parse_test!(ne "EM/3[+SPEKEVEE2Na+,+VEE2Na+,+HPFKAUBE2Nb+,+HPFKAUBW( E2Nb+,++]", fuzz_71);
parse_test!(ne "EM/3[+SPEKEVEE2Na+,+VEE2Na+,+HPFKAUBW( E2Nb+,+HPFKAMOD:0200L11[X]2|||||UBE2Nb++++]", fuzz_72);
parse_test!(ne "{}", fuzz_73);
parse_test!(ne "SEGUEN[Formula:[13C2][12C-2]H W$Ґq]", fuzz_74);
parse_test!(ne "SEGUEN[Formula:[13C2][12C-2]H 555555555W$Ґ555W]", fuzz_75);
parse_test!(ne "SEGUEN[Formula:[13C2][12C-2]Hǂ> 555555555555q]", fuzz_76);
parse_test!(ne "EMEVEESPEK/3[+2W$à+H+]", fuzz_77);
parse_test!(ne "EEMEVEESPEK/3[+2NW$àa+,H+]", fuzz_78);
parse_test!(ne "EMA+HBSPEK/3[+2N<ǂ>A5HAa+,+H+]", fuzz_79);
parse_test!(ne "E/3[+2NE+,+sssss<ǂ>A5HAsssssssssssssssssssssssssssssssssssssssssssssssssssssH+][U:Phossssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssspho]?", fuzz_80);
parse_test!(ne "EMEVEESPEK/3[+2N  皚+]", fuzz_81);
parse_test!(ne "EMEVEVEESPEKrmEESPEJ/3[+ᚚ2/52Na000N00000000000000Na000000Na+,+H+]", fuzz_82);
parse_test!(ne "Z+[13B2]^12777777777777", fuzz_83);
parse_test!(ne "EMEVEESPEK/3[+2Nᛚ-teoedMEVEESPEK/3[+2Nᛚ-teoedAa+,+H+]Aa+,+H+]", fuzz_84);
parse_test!(ne "FGCICNDBXTDOALBAHJFT/3[+2NaSPEKEVE00W$Ґ00000000000000110Na+,+H+]", fuzz_85);
parse_test!(ne "EMEVEESPEKEVEESPEK/3[+[13C2][12C-2]HW$à2Neeedeeeeeeeeeea+,+H+]/", fuzz_86);
parse_test!(ne "HMEVEESPEKEVEESPEK/3Z+[13C2]^-2E-2]H2Na/,+H+]/3[", fuzz_87);
parse_test!(ne "EVEESEVEESPEK/3[+2NaE33W$àEKKB+,+)+]PEK/3[+2NaE33W$àEKKB+,+)+]", fuzz_88);
parse_test!(ne "FJAVEXBDVGDGEYAH/3[+2NaSPSPEKEVE0EKEVE00Na+,+H+01,010010Na+,+H+11,000000100EKEVE00Na+,+H+00,11000001E00Na+,+H+10,010001N<ǂ>a+,+H+11,100000111EKEVE00Na+,*H+00,11000D10]", fuzz_89);
parse_test!(ne "BYNDUCDBCEQUEEeEVEESPEK/3[+2Na+,+H+][i[Oxidatinn^on]^-2]SK[#XL1]PEKCE//M[XLMOR:02001LT////QBT", fuzz_90);
parse_test!(ne "EMECIDWDDMJLOAIAPCFQRBO/3[+2EEEEE0000000000000000+,+HNaSPEKEVBFEEEEEEEEEEEEW$à000000+,+HNaSPEKEVEEEEEEEEEEEEEEEEE0000000000+,+HNaSPEKEVBFEEEEEE8EEEEEEE000000/00+,+H+]0", fuzz_91);
parse_test!(ne "EMEVESPEK/3[+2eW$Ґa+,+HH+]", fuzz_92);
parse_test!(ne "<d><d>NEKEVEESPEMEVEESPEK/3[8P  皚H+]/333333[+2N", fuzz_93);
parse_test!(ne "EMEgEgEE//V//HMEgMEVRAL[Formula:[13C2][12C-2]N4EEMEVEESPEKEEMEVENSPSEM[:Oyidatiom3<cEVENSPSEM[:Pxidatiom3<cvU:Oxidatiom3<cn]W$à{+2Na+,+B+N+K,1,c-term>A+H+MEVSPnU:Oxidatiom3<cn]{+2Na+,+B+N+g,K,c-ter5]WEN[Formula:W1SSSSSSSSSSSSSSS-2EHEVRAL[Formula:[13C2[12C-2]N45]WEN[Formula:W1SSSSSSSSSSSSSSS", fuzz_94);
parse_test!(ne "jTRA[Formula:[1SB2]ymEMK[X:[07C     ---,  򆼡,-----   ]hZyvLMOD:0203[13AU]EVTKSESPgK/31 2Na+    ---,  򆼡,-----   ]hZyvLMeeeee~eeeeeeeeeeeeMeeee eeILTClSIGDE+K+K[U:iTidatioOxidation]nQEVNhxiL1]EVTKSESPrK/", fuzz_95);
parse_test!(ne "EEiTRA[formula:[13C2]KKKKKKKKKKKKKKKKtion]hZA[Formula:o]iTRA[formula:[13333333333333333333333333333333333333333333333333333333333333iTRAQ4plex]-EM[U:Oatdxi33333333333333EMM[U:QdatioɒLdation]nQEVNhS[xxidati]n]EEEM33333333333333333333333333333C2]KKKn]hZ", fuzz_96);
parse_test!(ne "<d>", fuzz_97);
parse_test!(ne "<[TMT2plex]@K,N-Term<Rf>_ZSU\\`R>", fuzz_98);
parse_test!(ne "<[TMT6plex]@K,N-term>", fuzz_99);
parse_test!(ne "<[TMT6plex]@K,N-term>", fuzz_100);
parse_test!(ne "<[TMT6plex]@K,N-termKATPEILTCNSIGoooooooooooooooooooooooCL>", fuzz_101);
parse_test!(ne "<[TMt6plex]@K,N-termK,N-term)A)AqaqIRORXO6plATPEILqq,N-termK,N-term>", fuzz_102);
parse_test!(ne "<111TC>", fuzz_103);
parse_test!(ne "<d>", fuzz_104);
parse_test!(ne "<[TMT6plex]@K,N>", fuzz_105);
parse_test!(ne "<d>", fuzz_106);
parse_test!(ne "<111TC><100Y>", fuzz_107);
parse_test!(ne "<[TMt6plex]@K,N-termK,c-termK,K,N-termK,c-termKNr,c-term<111TC>TEEESSSSSSSSSSSSSSSSSSSSSSSSSSSSSSatLoln->", fuzz_108);
parse_test!(ne "<111TC><111TE><111TE><111TC>", fuzz_109);
parse_test!(ne "<11c>", fuzz_110);
parse_test!(ne "<111nB>", fuzz_111);
parse_test!(ne "<111SB><111SB><111SB><111SB><111SB>", fuzz_112);
parse_test!(ne "<111TC><111TC><111TC><11c><111TC><111TC><111TC><111TC><11c><111TC><111TC><111TC><111TC><111TC><111TC><11c><111TC><111TC><111TC><111TC><10b><111TC><111TC><111TC>", fuzz_113);
parse_test!(ne "<[TMT6plex]@L,N-term<lI><lOlex>TC>", fuzz_114);
parse_test!(ne "<[TMT6plex]@K,k>", fuzz_115);
parse_test!(ne "<[]@K,k>", fuzz_116);
parse_test!(ne "<[TMt6plex]@m>", fuzz_117);
parse_test!(ne "<[TMt6plex]@K,N-termK,N-term,N-term4F:K,r,N-term>", fuzz_118);
parse_test!(ne "[Glycan:SSSFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDEC5FFFFFFF7777777777777777777SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS767777FFFFFDEC5HexS0FFFFFFFFFFFFFFFFSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDEC5FFFFDEC5FFFFFFF7777777777777777777SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS767777FFFFFDEC5HexS0FFFFFFFFFFFFFFFFSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDEC5FFFFFFF7777777777777777777SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSnSSSSSSSSSSSSSSSSSSSSSSSSSSSS767777FFFFFDECFFFFFFF7777777777777777777SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSnSSSSSSSSSSSSSSSSSSSSSSSSSSSS767777FFFFFDEC5FFFFDEC5FFFFFF0]?EMHHHHHHHHFFFFFFFFFNE", fuzz_119);
parse_test!(ne "()[Glycan:H0]AA/LAA)[Glycan:H0](LAA)AA(LAA)[Glycan:H0]AA(LAA)[Glycan:H0](LAA)-[Me[For:[12C-P]4?Q[[Glycan:HexSHHSHgHHHEPP-[Me[", fuzz_120);
parse_test!(ne "<[TMt6plex]@K,N-termK,c-termK,K,N-termK,c-termKNrmG,c-termK,K,N-termK,c-termoK,K,K,N-termKKeq,N-termK,c-termK,K,N-termK,c,K,N-term-tertK,N-term>", fuzz_121);
parse_test!(ne "<[TMt6plex]@K,N-termerm,N-term4F:K,r,N-termcan:K,m,N-term4pxrm,C-term4@111nB><111nB><111nB><111nB><111nB><111nB><111nB><111nB><111nB><111nB>", fuzz_122);
parse_test!(ne "<[TMT2plex]@K,m>", fuzz_123);
parse_test!(ne "<[TMT6plex]@K,C-termn:H>", fuzz_124);
parse_test!(ne "<d><d>", fuzz_125);
parse_test!(ne "<[TMT6plex]@K,N-term>", fuzz_126);

parse_test!(ne "[iT4plex]-EM[Oxidation]EVNES[ES[Phospho]PEK[iTRAQ4plex]-[Methyl]", hang_01);
parse_test!(ne "[TRAQ4plex]-EM[     iTRAQ4plex]-EM[      ", hang_02);
parse_test!(ne "<< II0>", hang_03);
parse_test!(ne "[iTRAQFplex]-EM[Oxidation]EVNES[Phospho]PEKniTREMK[XLMOD:02000#XL1]]EVTKSET:02000#XL1E[XLMOD[XLMOD:02  򆼡", hang_04);
parse_test!(ne "iTRA[Formula:[  |||||||||||||||||||||]hZsqh|||||||    13C      ]E", hang_05);
parse_test!(ne "uHLC[Formula:[13C2][19C2]n12W-20]EO[avqrx|hNiwxy]E", hang_06);
parse_test!(ne "<5BE>[qTR(Qplex]-EM[Oxidation]EVNES[Phospho]EEEEEEEEEEEEEEEEELEEEEEEEEEEEEEMK[XLMOD:021]E+K", hang_07);
parse_test!(ne "[Formula:K13pUKOLM[0s|ho]DLM[0spho]]-", hang_08);
parse_test!(ne "[MOD:020110#EMEspho|O-phosL-serine|OCE//K[][]UENCE//K[]UENCE", hang_09);
parse_test!(ne "[U=hTRAQ4plex]-EM[U:Oxidation]EVNES[U:Phospho]PEK[U:iTRAQ4plex]-[U:Methyl]/3", hang_10);
parse_test!(ne "[L:iTRAQ4plex]-EM[U:Oxidation]EVNES[U:Phospho]PEK[U:iTRAQ4plex]-[U:Methyl]/3", hang_11);
parse_test!(ne "[U=iTRAQ4plex]-EM[U:Oxidation]EVNES[U:Phospho]PEK[U:iTRAQ4plex]-[U:Methyl]/3", hang_12);
parse_test!(ne "[ :iTROQ4plex]-EM[U::::::::::::::::/::::::EK[U:iTRAQ4plexd-[U::::::EK[U:iTRAQ4Methyl]/3", hang_13);
parse_test!(ne "[U TRAQ4plex]-EM[U:Oxidation]EVKpho]PEK[U:iTRAQ7zlex]-[U:Metx]-[U:Methyl]/3", hang_14);
parse_test!(ne "[:Oxidation]K[U:iplex]6dU:MethEVNES[]-EM[U:Oxidation]K[U:iplex]6dU:Methyl]/3", hang_15);
parse_test!(ne "plS[UQ4plex]-[U:Methyl]ex]-EM[U:Oxidation]EVNES[U:Phospho]PEK[U:iTRAQ4plex]-[U:Methyl]/3", hang_16);
parse_test!(ne "[V:iTRAQ4plex]-EM[UOxidation]EVNES[U>Phospho]PEK[U:iTRAQ4plex]-[U:Methyl]/3", hang_17);
parse_test!(ne "[;iTRAQ4plex]-Etiion]EVNES[U:Phthyl]/3", hang_18);
parse_test!(ne "[Uation]EVNES[U:i ospho]PEK[U@iTRAQ4plex]-U:Methyl]/3", hang_19);
parse_test!(ne "[URi:RAQ4plex]-EMM=MMn]EVNlex]-[oOho]PEK[U:iT2AQ4plex]-0U'Methyl]/3", hang_20);
parse_test!(ne "o[	:iTRAQ4Olex]-EM[U:Oxidatioj`DIWES[ K[U:i2[U:iTRAQU:Phospho]PE", hang_21);
parse_test!(ne "[e:iTRAQ<plex]-E:Ph ", hang_22);
parse_test!(ne "[dehyaio]^3?Q[gln->pyro-glu]SC", hang_23);
parse_test!(ne "[deidationhydro]^3?QRAidationo-glu]SC", hang_24);
parse_test!(ne "[dehyehydro]^3?Q[glC->pyro-", hang_25);
parse_test!(ne "  gldro]^3?Q[glC->pyro-glu]Sn", hang_26);
parse_test!(ne "[d))))))))))))))))))))))))ehydro]^3?Q[ghn-=pyro-gluISC", hang_27);
parse_test!(ne "MPGLVDSNPAPPESQEKKPLK(PCCACPETKKARDACIIEKGEEHCGHLIEAHKECMRALGFKI)[Oxdation]^2[half cystine][half cystine]", hang_28);
parse_test!(ne "MPGLVDSNPAPPESQEKKPLK(AHKECMRALGFKI)[Oxidation]^2[half cystine][half cystine]", hang_29);
parse_test!(ne "[dehydro]^6?Q[gln->pyro-gxSEK[XLMOD:02001#?L1[U:U:iTRAQ4plex]-EM[U:Oxidation]	 NiTRA[U:iTRlnWQ4plexgln->pyro-glu]Q4plex]-EM|UFOxidation]E", hang_30);
parse_test!(ne "xSEK[XLMOD:02011#XL1]UENCE//EMEVTK[XLMO:02009#XL1]SESPE", hang_31);
parse_test!(ne "SEQUEN[Formula:B13C2][12C-2]H2N]CE", hang_32);
parse_test!(ne "SEQUTN[Formula:u13C2][12C-2]H2N]CE", hang_33);
parse_test!(ne "SEQUEN[Fo-glu]SC>pyro-gl->pyro-glu]SC> NE				", hang_34);
parse_test!(ne "SEQUEN[Formula:C13C2][1 C-2]H2N]CE", hang_35);
parse_test!(ne "xSETK[xSEK[XLMOD:02001#?L1[U:U:iTRAQ4plex]-EM[U:Oxi   on]	 NiTRAQ4plex]-EM|U:Oxidation]EVNES[U@Phogpho]PEK[", hang_36);
parse_test!(ne "<[TMT2plex]@K,N-term>AWhDRNHA+NHA+[U*iTRAQ4plex]-EM[U:Oxidation]EVNES[U:Phos", hang_37);
parse_test!(ne "SEQWEN[Formula:[13C22N]CE][12C  2]H2N]CE", hang_38);
parse_test!(ne "[1][deh 1    1#XM1][deh 1#XM^3C- 1#XMXLMOD:D:#XM^3C-Q", hang_39);
parse_test!(ne "SEQUEN[Formula:[13C2C-2]H2NNNNNNNNNNNNNNN2][12C-2]H2NNNNNNNNNNNNNNNNN]CE", hang_40);
parse_test!(ne "SEQUEN[Formula:[13C2][12C-2]H2NNN2][12C-2]H2NNNNCE", hang_41);
parse_test!(ne "SEQUEN[Formula:[13C2][12C-2]H2NNNNNNNNNNNN13C2][12C-2]H2NNNNNNNN]CE", hang_42);
parse_test!(ne "SEQUEN[Formula:[13C2][12C-2]H213C2][12C-2]HNNNN]eE", hang_43);
parse_test!(ne "ex[ P:Met(yl]8S][U:O2C-ͨ", hang_44);
parse_test!(ne "SEQUEN[Formula:[05B2][12C2]H2NN0 |NNNNNNNNNNNNNN]CE", hang_45);
parse_test!(ne "SEQUEN[Formula:[16C2][12C-2]H2222][12C-2]22N:NNNNNhONNNNNN2222212:2222NNNNNNNhNNBNNNNNN]]", hang_46);
parse_test!(ne "NFJGDN[Formula:[13C2][12C-2]H2NNNcdrfffffffffffffffffffffffffffNNNNNYAMEVEAXFPKBYAVZBW  2010[C:[pldt]-[LD:i			 NES[thybfmwpxionngzzng[U:hVpl[U:Oxidatio|n]]]]]_][^", hang_47);
parse_test!(ne "[Glycan:HHHNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNPeNNNNGNPNN][U:|TRAQ4plex]-EM", hang_48);
parse_test!(ne "[Formula:C-2][XLMODmula:C-2][XLMODeh 1#XM^3C- 1#XM1][deh 1#XM^(C- ", hang_49);
parse_test!(ne "[12B-2]H2u¨", hang_50);
parse_test!(ne "NFJNNNNNNNN[Glycan:P0][    h][    h 1#XM^3#- 1# M1", hang_51);
parse_test!(ne "[dehydrotio]^<?Q[[1]]SEro dAA]^	 NE	 ", hang_52);
parse_test!(ne "[ͨH2>Z]uU:򆼡Kngthhhhhhhhhh@K,N-serm>Z", hang_53);
parse_test!(ne "[Glycan:HHHHHHHHHHHH0HHHHH0HHHHHNNNHNNNNNN2][12C-2]H22222NNNNNNNhNNNNNNN222222222222NNNNNNfNNHHHHHHH0HHHHH0HHHHHNNNNNNNNNNNNNfNNHHHHHHH0HHHHH0HHHHHNNNNN0]-T", hang_54);
parse_test!(ne "[Glycan:HHHHHTRIHHHH0HHHHH0HHHHHNNNHNNNNNNNfNNHHHHHHH0][12CL2]H2NNNNNENhNNNNHHFYFIFDPTDK/3NENMCHBJQCN[M:06101[C:[plex]-[UC:i			 NES[thydatiotionngting[U:iTRAQ%pl[U:Oxidation]EV  S[U:PhosphoEVEESPEex]-EM[U:O", hang_55);
parse_test!(ne "{Glycan:Hex}[iTRAQ4plex]-EV[Oovered|IRAQ4plex]-EV[Oovered|INFxidatNFxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methylan:]", hang_56);
parse_test!(ne "[iTRAQ4plVx]-EM[O   tion]DVNESHP{M[Ox		hospho]PPPPEK[iTRAQ4plex]-[Met򆼠", hang_57);
parse_test!(ne "[Glycan:HexSHHHHH1HHTET][Hex][dehydro6plfx][Hex][deydro]^ESA", hang_58);
parse_test!(ne "[Glycan: ][Glyc   GS0][Glycan:{GlyexS0][G{GS0x]@K,N>ATQQQLTCNly", hang_59);
parse_test!(ne "[dehydro]^-5[Glycan:HexSHHHH1HH@E3333333fx][d   [de@exSHH^HH1HH@ET6H13333333fx][HexGG ][dehydro]^-2[Glycan:HexSHHHH1HH@E3333333fx][d   [dehydro]^-______________________________2\\H\\H", hang_60);
parse_test!(ne "[Formula:FMTcNu:PHRXRwW2ec|lmZzeqokjzkrTp]-", hang_61);

// Personal tests

#[test]
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
        "Q[U:Gln->pyro-Glu]E[Cation:Na]AA",
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
        peptide.unwrap().formulas(0, 0),
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
        singular.formulas(0, 0).to_vec()[0].elements(),
        (AminoAcid::Alanine.formulas(0, 0).to_vec().pop().unwrap() * 2
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
        peptidoform.formulas(0, 0).to_vec()[0].elements(),
        (AminoAcid::Alanine.formulas(0, 0).to_vec().pop().unwrap() * 2
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

#[test]
fn parse_adduct_ions_02() {
    let peptide = dbg!(CompoundPeptidoform::pro_forma("A-[+1]/2[1Na+,+H+]+[+1]-A", None).unwrap());
    assert_eq!(peptide.peptidoforms().len(), 2);
    assert_eq!(
        peptide.peptidoforms()[0].peptides()[0]
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
        peptide.peptidoforms()[0].peptides()[0].sequence[0].formulas_all(
            &[],
            &[],
            &mut Vec::new(),
            false,
            0,
            0,
        ),
        peptide.peptidoforms()[1].peptides()[0].sequence[0].formulas_all(
            &[],
            &[],
            &mut Vec::new(),
            false,
            0,
            0,
        )
    );
    assert_eq!(
        peptide.peptidoforms()[0].peptides()[0]
            .get_c_term()
            .additional_mass(),
        peptide.peptidoforms()[1].peptides()[0]
            .get_n_term()
            .additional_mass()
    );
    assert!(
        peptide.peptidoforms()[0].peptides()[0].get_n_term()
            != peptide.peptidoforms()[1].peptides()[0].get_c_term()
    );
}

#[test]
fn should_not_crash() {
    // Semantically not valid, but syntactically it is
    let _ = CompoundPeptidoform::pro_forma(
        r"<[Formula:[318  N  005][  7Rf1][93Sb   ]  |5.3|UNIMOD:>]@N-term,N-term,k>[R:-3.70#BRANCH]^012?[Info:[_]=#XL0|R:+734.7|-4.74]^64[Info:[]|GNO:( |Formula: [2  Ca  ] Nd5 [5Sg ]]?{UNIMOD:;}(?B[#BRANCH][Glycan:RES\n1b:x-lgal-HEX-x:x|6:d\nLIN\n235 ?-?-Tetx6 #BRANCH]k)(i[#XLd6]h[#XLW8][#BRANCH])[XLMOD:[=\]-#XL88|UNIMOD:\|[|][?.].|M:=[\][=]_](O[#XLtD]I[Formula: Es Ce  Bi   |Glycan:HexN1sulfate170 ?-?-Trix |Glycan:Neu5Ac903 ][GNO:/?#XL32s|RESID:/[][)-]{#XL35|[{]#mG(-7.5)])[G:\]-[#M(20.5)]",
        None,
    );
}
