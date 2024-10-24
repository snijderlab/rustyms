#![allow(clippy::missing_panics_doc)]
use std::io::{BufReader, Read};

use crate::Modification;
use crate::{modification::SimpleModification, molecular_formula};

use super::error::CustomError;
use super::{IdentifiedPeptideSource, PeaksFormat};

use super::{csv::parse_csv_raw, peaks, IdentifiedPeptide, PeaksData};

#[test]
fn peaks_x() {
    assert_eq!(
        open_file(BufReader::new(DATA_X.as_bytes()), &peaks::X, None).unwrap(),
        19
    );
}

#[test]
fn peaks_x_plus() {
    assert_eq!(
        open_file(BufReader::new(DATA_XPLUS.as_bytes()), &peaks::XPLUS, None).unwrap(),
        19
    );
}

#[test]
fn peaks_11() {
    assert_eq!(
        open_file(BufReader::new(DATA_11.as_bytes()), &peaks::V11, None).unwrap(),
        19
    );
}

#[test]
fn peaks_11_all_candidates() {
    assert_eq!(
        open_file(
            BufReader::new(DATA_11_ALL_CANDIDATES.as_bytes()),
            &peaks::V11,
            None
        )
        .unwrap(),
        19
    );
}

#[test]
fn peaks_11_custom_modification() {
    assert_eq!(
        open_file(
            BufReader::new(DATA_11_CUSTOM_MODIFICATION.as_bytes()),
            &peaks::V11,
            Some(&vec![(
                Some(0),
                "oxidation".to_string(),
                SimpleModification::Formula(molecular_formula!(O 1)),
            )])
        )
        .unwrap(),
        19
    );
}

#[test]
fn peaks_12() {
    assert_eq!(
        open_file(BufReader::new(DATA_12.as_bytes()), &peaks::V12, None).unwrap(),
        19
    );
}

#[test]
fn peaks_ab() {
    assert_eq!(
        open_file(BufReader::new(DATA_AB.as_bytes()), &peaks::AB, None).unwrap(),
        19
    );
}

#[test]
fn full_peaks_file() {
    for pep in PeaksData::parse_file("data/200305_HER_test_04_DENOVO_excerpt.csv", None).unwrap() {
        if let Err(e) = pep {
            panic!("{}", e);
        }
    }
}

#[test]
fn fuzz_crashes() {
    let mut all_passing = true;
    for entry in std::fs::read_dir("src/identification/peaks_fuzz_tests").unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();
        let metadata = std::fs::metadata(&path).unwrap();
        let name = path.file_name().unwrap().to_str().unwrap();
        if metadata.is_file() && name.starts_with("crash") {
            match std::panic::catch_unwind(|| {
                open_file(
                    BufReader::new(std::fs::File::open(&path).unwrap()),
                    &peaks::V11,
                    None,
                )
            }) {
                Ok(_) => {
                    println!("identification::peaks_tests::fuzz_crashes::{name} passed");
                }
                Err(err) => {
                    all_passing = false;
                    println!("---- identification::peaks_tests::fuzz_crashes::{name} crash report ----\n{err:#?}\n");
                }
            }
        }
    }
    assert!(all_passing, "Some fuzz tests did not pass");
}

#[test]
fn fuzz_hangs() {
    let mut all_passing = true;
    for entry in std::fs::read_dir("src/identification/peaks_fuzz_tests").unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();
        let metadata = std::fs::metadata(&path).unwrap();
        let name = path.file_name().unwrap().to_str().unwrap();
        if metadata.is_file() && name.starts_with("hang") {
            match std::panic::catch_unwind(|| {
                open_file(
                    BufReader::new(std::fs::File::open(&path).unwrap()),
                    &peaks::V11,
                    None,
                )
            }) {
                Ok(_) => {
                    println!("identification::peaks_tests::fuzz_hangs::{name} passed");
                }
                Err(err) => {
                    all_passing = false;
                    println!("---- identification::peaks_tests::fuzz_hangs::{name} crash report ----\n{err:#?}\n");
                }
            }
        }
    }
    assert!(all_passing, "Some fuzz tests did not pass");
}

/// Open a Peaks file from the given reader.
/// # Errors
/// If any part of the process errors.
fn open_file(
    reader: impl Read,
    format: &PeaksFormat,
    custom_database: Option<&super::ontologies::CustomDatabase>,
) -> Result<usize, CustomError> {
    let lines = parse_csv_raw(reader, b',', None)?;
    let mut num_lines = 0;
    for line in lines {
        let line = line?;
        let read: IdentifiedPeptide =
            PeaksData::parse_specific(&line, format, custom_database)?.into();
        num_lines += 1;
        assert!(
            read.peptide().unwrap().sequence().iter().all(|s| s
                .modifications
                .iter()
                .all(|m| !matches!(m, Modification::Simple(SimpleModification::Mass(_))))),
            "Peptide contains mass mods: {}",
            read.peptide().unwrap()
        );
    }
    Ok(num_lines)
}

const DATA_AB: &str = r"Scan,Peptide,Tag Length,ALC (%),length,m/z,z,RT,Area,Mass,ppm,Accession,PTM,local confidence (%),tag (>=0%),mode
F2:10351,MHQN(+.98)WLWL,8,98,8,564.7653,2,23.75,5.73E7,1127.5222,-5.5,,Deamidation (NQ),96 98 98 99 99 100 99 99,MHQN(+.98)WLWL,CID
F3:3063,M(+15.99)PHNHHTE,8,98,8,509.7123,2,10.99,4.42E6,1017.4087,1.4,,Oxidation (M),98 98 98 98 99 100 100 96,M(+15.99)PHNHHTE,CID
F3:3534,M(+15.99)PHNHHTE,8,98,8,509.7128,2,11.79,1.82E6,1017.4087,2.4,,Oxidation (M),98 99 99 98 99 99 100 95,M(+15.99)PHNHHTE,CID
F3:2117,HNHHTE,6,97,6,387.6678,2,8.52,2.03E6,773.3205,0.8,,,97 96 98 99 99 96,HNHHTE,CID
F2:13745,STMHWV,6,97,6,380.6751,2,30.62,2.73E6,759.3374,-2.2,,,98 99 95 95 99 99,STMHWV,CID
F1:20191,LFLN(+.98)ESHLTHAF,12,97,12,715.3615,2,37.09,4.4E6,1428.7036,3.3,,Deamidation (NQ),85 97 99 98 98 94 98 99 99 99 100 99,LFLN(+.98)ESHLTHAF,CID
F3:14603,APNTFTCSVLHE,12,97,12,659.8101,2,38.73,9.59E7,1317.6023,2.5,constructed_protein_Heavy,,96 96 98 99 100 99 96 95 96 99 98 94,APNTFTCSVLHE,CID
F3:15521,FTLNLHPVEEE,11,97,11,664.3321,2,41.05,2.65E6,1326.6455,3.1,,,99 99 99 97 99 99 90 92 99 100 94,FTLNLHPVEEE,CID
F4:8195,SNNYATHYAENK(+72.06),12,97,12,742.3466,2,18.77,6.26E7,1482.6792,-0.3,,Carboxyethyl,83 94 98 99 99 99 98 99 99 100 97 99,SNNYATHYAENK(+72.06),CID
F2:13141,WFVDLEEVHTA,11,97,11,673.3289,2,29.59,1.25E7,1344.6350,6.2,,,92 99 99 98 88 95 99 99 100 99 98,WFVDLEEVHTA,CID
F3:2323,HN(+.98)HHTE,6,97,6,388.1600,2,9.21,4.17E4,774.3045,1.2,,Deamidation (NQ),97 97 99 99 99 91,HN(+.98)HHTE,CID
F4:7286,ANN(+.98)HATYYAENK(+72.06),12,96,12,734.8412,2,17.49,7.73E6,1467.6682,-0.2,,Deamidation (NQ); Carboxyethyl,82 94 99 99 100 99 99 99 99 100 95 96,ANN(+.98)HATYYAENK(+72.06),CID
F2:7215,HYLHEV,6,96,6,399.1997,2,18.41,8.67E6,796.3868,-2.5,,,95 95 97 97 99 99,HYLHEV,CID
F1:20356,YLDQTEQWQLY,11,96,11,743.8483,2,37.34,1.04E7,1485.6775,3.1,,,90 97 99 96 99 99 96 98 94 97 96,YLDQTEQWQLY,CID
F3:14820,APNTFTCSVLHE,12,96,12,659.8107,2,39.24,6.09E6,1317.6023,3.5,constructed_protein_Heavy,,92 89 98 99 100 99 99 98 97 98 93 94,APNTFTCSVLHE,CID
F4:9062,SDNYATHYAENK(+72.06),12,96,12,742.8393,2,20.43,5.92E7,1483.6631,0.6,,Carboxyethyl,84 95 98 99 99 99 99 99 99 100 89 97,SDNYATHYAENK(+72.06),CID
F4:6305,NHATYYAENK(+72.06),10,96,10,641.8087,2,16.22,2.99E6,1281.6040,-1.0,,Carboxyethyl,85 96 99 99 99 98 98 100 94 94,NHATYYAENK(+72.06),CID
F2:8831,VCAAVHGV,8,96,8,378.1942,2,21.04,3.21E6,754.3796,-7.5,,,84 92 98 99 99 100 99 99,VCAAVHGV,CID
F4:3135,TPVSEHQK(+72.06),8,96,8,499.2701,2,10.87,3.28E7,996.5292,-3.5,,Carboxyethyl,99 99 99 97 99 94 88 95,TPVSEHQK(+72.06),CID";

const DATA_X: &str = r"Fraction,Source File,Feature,Peptide,Scan,Tag Length,ALC (%),length,m/z,z,RT,Area,Mass,ppm,PTM,local confidence (%),tag (>=0%),mode
1,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_Tryp01.raw,F1:5056,LLYLVSK,F1:6994,7,99,7,418.2689,2,39.59,1.47E6,834.5215,2.2,,100 100 100 100 100 100 100,LLYLVSK,ETHCD
2,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_Chym01.raw,-,MYYSKL,F2:6253,6,99,6,402.7029,2,35.37,,803.3887,3.2,,100 100 100 100 100 100,MYYSKL,ETHCD
4,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_LysN01.raw,F4:9634,KAAVFNHFLSDGVK,F4:6537,14,99,14,511.6133,3,36.67,2.16E6,1531.8147,2.2,,100 100 100 100 100 100 100 100 100 100 100 100 100 100,KAAVFNHFLSDGVK,ETHCD
6,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_TL01.raw,F6:4316,LLSDHRGATYA,F6:4340,11,99,11,401.8766,3,24.21,7.2E5,1202.6042,3.2,,100 100 100 100 100 100 100 100 100 100 100,LLSDHRGATYA,ETHCD
2,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_Chym01.raw,F2:137,KLAKVF,F2:5266,6,99,6,353.2372,2,29.28,6.29E7,704.4584,2.1,,100 100 100 100 100 100,KLAKVF,ETHCD
6,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_TL01.raw,F6:8687,LLKYASQS,F6:4154,8,99,8,455.2574,2,23.10,5.71E5,908.4967,3.9,,100 100 100 100 100 100 100 100,LLKYASQS,ETHCD
2,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_Chym01.raw,F2:137,KLAKVF,F2:5206,6,99,6,353.2372,2,29.28,6.29E7,704.4584,2.1,,100 100 100 100 100 100,KLAKVF,ETHCD
1,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_Tryp01.raw,F1:412,FTLSRDDSK,F1:4190,9,99,9,356.8489,3,23.11,1.72E6,1067.5247,0.3,,100 100 100 100 100 100 100 100 100,FTLSRDDSK,ETHCD
2,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_Chym01.raw,F2:6787,TLVGLVNY,F2:10516,8,99,8,439.7537,2,60.46,1.26E7,877.4909,2.3,,100 100 100 100 100 100 100 100,TLVGLVNY,ETHCD
6,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_TL01.raw,F6:16600,LTGTSTVGVGRGVLGDQKN,F6:5492,19,99,19,620.3372,3,31.03,2.2E7,1857.9907,-0.5,,100 100 99 99 100 100 100 100 100 100 100 100 100 100 100 100 99 100 100,LTGTSTVGVGRGVLGDQKN,ETHCD
1,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_Tryp01.raw,F1:10850,SSTLTLTKDEYER,F1:5255,13,99,13,514.9270,3,29.45,5.87E5,1541.7573,1.2,,100 100 100 100 100 100 100 100 100 100 99 100 100,SSTLTLTKDEYER,ETHCD
6,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_TL01.raw,F6:12891,LRSSVHYSQGYNNA,F6:3949,14,99,14,532.5908,3,21.95,1.03E7,1594.7488,1.2,,99 99 100 100 100 100 100 100 100 100 100 100 100 100,LRSSVHYSQGYNNA,ETHCD
6,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_TL01.raw,F6:1122,LRVEKKNW(+15.99),F6:3772,8,99,8,363.5458,3,20.31,1.99E6,1087.6138,1.8,Oxidation (HW),100 99 100 100 100 100 100 100,LRVEKKNW(+15.99),ETHCD
5,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_Ela01.raw,F5:11029,PYTFGGGTKLELKRA,F5:6127,15,99,15,546.6385,3,34.73,4.04E6,1636.8936,0.2,,99 99 100 100 100 100 100 100 100 100 100 100 100 100 100,PYTFGGGTKLELKRA,ETHCD
6,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_TL01.raw,F6:12937,VTYDYYKGG,F6:4919,9,99,9,533.2510,2,27.58,1.28E7,1064.4814,5.6,,99 99 100 100 100 100 100 100 100,VTYDYYKGG,ETHCD
4,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_LysN01.raw,F4:12679,KASGFTFTDFSLHMK(+58.01),F4:11257,15,99,15,592.2878,3,62.97,4.75E6,1773.8396,1.2,Carboxymethyl (KW  X@N-term),99 100 100 100 100 100 100 100 100 100 100 100 100 100 100,KASGFTFTDFSLHMK(+58.01),ETHCD
6,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_TL01.raw,F6:15331,VFTEQADLSGLTETKK,F6:6459,16,99,16,589.6451,3,36.75,9.07E5,1765.9097,2.2,,100 100 100 100 99 100 100 100 100 100 100 100 100 100 100 100,VFTEQADLSGLTETKK,ETHCD
5,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_Ela01.raw,F5:6860,TFGAGTKLELKRA,F5:5414,13,99,13,464.6069,3,30.47,2.49E6,1390.7932,4.0,,100 100 100 100 100 100 100 100 100 100 100 100 100,TFGAGTKLELKRA,ETHCD
2,20190517_F1_Ag5_3117030_SA_ETHCD_131-2a_Chym01.raw,F2:7876,TLSRDDSKSSVY,F2:3962,12,99,12,453.2252,3,21.77,3.86E6,1356.6521,1.1,,100 100 100 100 100 100 100 100 100 99 100 100,TLSRDDSKSSVY,ETHCD";

const DATA_XPLUS: &str = r"Fraction,Source File,Feature,Peptide,Scan,Tag Length,Denovo Score,ALC (%),length,m/z,z,RT,Predict RT,Area,Mass,ppm,PTM,local confidence (%),tag (>=0%),mode
10,20191211_F1_Ag5_peng0013_SA_her_Asp_N.raw,F10:3434,DYEKHKVYAC(+58.01),F10:3629,10,99,99,10,438.5332,3,19.91,-,2.3176E6,1312.5757,1.6,Carboxymethyl,100 100 100 100 100 100 100 100 100 100,DYEKHKVYAC(+58.01),ETHCD
4,20191211_F1_Ag5_peng0013_SA_her_Ela.raw,F4:4797,SGFGGLKN(+.98)TYLHW,F4:9505,13,99,99,13,494.2459,3,52.43,-,2.4924E7,1479.7146,0.9,Deamidation (NQ),100 100 100 100 100 100 100 100 100 100 100 100 100,SGFGGLKN(+.98)TYLHW,HCD
3,20191211_F1_Ag5_peng0013_SA_her_thermo.raw,F3:12703,LSC(+58.01)AASGFNLKDTY,F3:7983,14,99,99,14,774.3562,2,43.80,-,7.2888E7,1546.6973,0.4,Carboxymethyl,99 100 100 100 100 100 100 100 100 100 100 100 100 100,LSC(+58.01)AASGFNLKDTY,HCD
11,20191211_F1_Ag5_peng0013_SA_her_CB.raw,F11:14212,DSTYSLSSTLTLSK,F11:8043,14,99,99,14,751.8829,2,44.29,-,3.8265E6,1501.7512,0.0,,99 99 100 100 100 100 100 100 100 100 100 100 100 100,DSTYSLSSTLTLSK,ETHCD
2,20191211_F1_Ag5_peng0013_SA_her_chymo.raw,F2:3799,SGFGGLKDTYLHW,F2:9542,13,99,99,13,494.2458,3,52.79,-,5.6351E5,1479.7146,0.6,,100 100 100 100 100 100 100 100 100 100 100 100 100,SGFGGLKDTYLHW,HCD
11,20191211_F1_Ag5_peng0013_SA_her_CB.raw,F11:7245,SSPVTKSFNRGEC(+58.01),F11:4201,13,99,99,13,490.5618,3,22.70,-,3.2875E6,1468.6616,1.3,Carboxymethyl,100 100 100 100 100 100 100 100 100 100 100 100 100,SSPVTKSFNRGEC(+58.01),ETHCD
3,20191211_F1_Ag5_peng0013_SA_her_thermo.raw,-,VGVGRGVLGDQKN,F3:4406,13,99,99,13,433.5779,3,24.14,-,0,1297.7102,1.2,,100 100 100 100 100 100 100 100 100 100 100 100 100,VGVGRGVLGDQKN,ETHCD
7,20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw,F7:14232,YSSDEKVLGEDFSDTR,F7:6514,16,99,99,16,616.6145,3,35.10,-,7.3986E5,1846.8220,-0.2,,100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 98,YSSDEKVLGEDFSDTR,HCD
2,20191211_F1_Ag5_peng0013_SA_her_chymo.raw,F2:4750,TLSKADYEKHKVY,F2:3789,13,99,99,13,527.9471,3,20.26,-,3.2229E9,1580.8198,-0.3,,99 100 100 100 100 100 100 100 100 100 100 100 100,TLSKADYEKHKVY,ETHCD
5,20191211_F1_Ag5_peng0013_SA_her_aLP.raw,F5:10958,TPPVLDSDGSFFLYSKLT,F5:11087,18,99,99,18,663.0070,3,62.15,-,6.1345E6,1985.9985,0.3,,98 99 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100,TPPVLDSDGSFFLYSKLT,ETHCD
4,20191211_F1_Ag5_peng0013_SA_her_Ela.raw,F4:6762,S(+58.01)GFNLKDTYLHWV,F4:10815,13,99,99,13,546.6037,3,60.06,-,2.739E7,1636.7886,0.4,Carboxymethyl (KW  X@N-term),100 100 100 100 100 100 100 100 100 100 100 99 99,S(+58.01)GFNLKDTYLHWV,HCD
8,20191211_F1_Ag5_peng0013_SA_her_Lys_C.raw,F8:5426,FNWYVDGVEVHDAK,F8:8307,14,99,99,14,560.2675,3,46.34,-,3.56E8,1677.7786,1.3,,99 100 100 100 100 100 100 100 100 100 100 100 100 100,FNWYVDGVEVHDAK,HCD
2,20191211_F1_Ag5_peng0013_SA_her_chymo.raw,F2:4750,TLSKADYEKHKVY,F2:3722,13,99,99,13,527.9471,3,20.26,-,3.2229E9,1580.8198,-0.3,,99 100 100 100 100 100 100 100 100 100 100 100 100,TLSKADYEKHKVY,ETHCD
4,20191211_F1_Ag5_peng0013_SA_her_Ela.raw,F4:7864,STLTLSKADYEKHKV,F4:4433,15,99,99,15,573.9813,3,24.22,-,6.8416E6,1718.9202,1.1,,99 99 100 100 100 100 100 100 100 100 100 100 100 100 100,STLTLSKADYEKHKV,ETHCD
11,20191211_F1_Ag5_peng0013_SA_her_CB.raw,F11:3224,LASYLDK,F11:4979,7,99,99,7,405.2239,2,26.98,-,6.041E4,808.4330,0.3,,100 100 100 100 100 100 99,LASYLDK,ETHCD
11,20191211_F1_Ag5_peng0013_SA_her_CB.raw,F11:9144,FNWYVDGVEVHNAK,F11:7980,14,99,99,14,559.9395,3,44.01,-,9.8729E4,1676.7947,1.2,,100 100 100 100 100 100 100 100 100 100 100 99 100 100,FNWYVDGVEVHNAK,HCD
8,20191211_F1_Ag5_peng0013_SA_her_Lys_C.raw,F8:8931,S(+58.01)TSGGTAALGC(+58.01)LVK,F8:7638,14,99,99,14,690.8375,2,42.67,-,2.888E7,1379.6602,0.2,Carboxymethyl (KW  X@N-term); Carboxymethyl,99 99 100 100 100 100 100 100 100 100 100 100 100 100,S(+58.01)TSGGTAALGC(+58.01)LVK,HCD
3,20191211_F1_Ag5_peng0013_SA_her_thermo.raw,F3:3555,LTLSKADYEKHK,F3:3745,12,99,99,12,478.2650,3,20.28,-,1.1906E8,1431.7722,0.7,,100 100 100 100 100 100 100 99 99 100 100 100,LTLSKADYEKHK,HCD
11,20191211_F1_Ag5_peng0013_SA_her_CB.raw,F11:12673,EVQLVESGGGLVQPGGSLRAK,F11:8595,21,99,99,21,694.3840,3,47.38,-,7.0626E7,2080.1274,1.4,,98 98 99 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100,EVQLVESGGGLVQPGGSLRAK,ETHCD";

const DATA_11: &str = r#""Source File","Scan","Peptide","Tag length","ALC (%)","Length","m/z","z","RT","Area","Mass","ppm","PTM","local confidence (%)","mode",tag(>=0.0%),"Feature Id"
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",6514,YSSDEKVLGEDFSDTR,16,99.9,16,616.61450,3,35.1011,190625.72,1846.8220,-0.2,"",100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 98,HCD,YSSDEKVLGEDFSDTR,14223
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",6743,DLQMTQSPSSLSASVGDR,18,99.2,18,626.96686,3,36.1993,794748.8,1877.8789,-0.1,"",99 99 100 100 100 100 99 99 100 100 100 100 100 100 100 100 100 90,ETHCD,DLQMTQSPSSLSASVGDR,14572
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",5644,DLQM(+15.99)TQSPSSLSASVGDR,18,99.0,18,947.94415,2,30.1547,548587.2,1893.8738,0.0,"Oxidation (M)",97 96 100 100 100 100 100 100 100 100 100 100 100 100 100 95 97 97,HCD,DLQM(+15.99)TQSPSSLSASVGDR,22869
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",6616,EVNSQFFR,8,98.8,8,513.75372,2,35.6492,60834.707,1025.4930,-0.1,"",98 97 99 100 99 100 100 97,HCD,EVNSQFFR,10540
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",9428,TYTFDNGTFLLR,12,98.7,12,724.36444,2,51.1627,83370.51,1446.7144,0.0,"",99 99 100 100 100 99 97 98 99 99 99 96,HCD,TYTFDNGTFLLR,17359
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",7805,EVQLVESGGGLVQPGGSLR,19,98.6,19,628.00598,3,42.2104,298160.3,1880.9955,0.3,"",98 97 99 100 100 100 100 100 100 100 100 100 98 99 99 99 100 100 83,ETHCD,EVQLVESGGGLVQPGGSLR,14598
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",5945,VLGEDFSDTR,10,98.0,10,569.77283,2,31.8568,18040.87,1137.5302,0.8,"",100 100 99 99 99 96 95 98 98 95,HCD,VLGEDFSDTR,12617
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",9483,TYTFDDGTFLLR,12,98.0,12,724.85663,2,51.3811,43949.74,1447.6984,0.3,"",96 97 100 100 100 100 96 96 97 99 99 95,HCD,TYTFDDGTFLLR,17378
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",3737,SRSGGGGNGLGSGGSLR,17,97.9,17,492.58188,3,19.8131,29103.092,1474.7236,0.1,"",100 99 100 100 100 99 98 95 100 100 100 100 98 98 99 99 83,ETHCD,SRSGGGGNGLGSGGSLR,9583
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",4286,TSGVLPR,7,97.8,7,365.21625,2,22.6110,744909.94,728.4181,-0.2,"",100 100 100 100 100 100 85,ETHCD,TSGVLPR,1074
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",6256,QAPGKGLEWVAR,12,97.8,12,437.91058,3,33.6266,,1310.7095,0.3,"",94 98 99 100 100 100 100 100 100 100 100 82,ETHCD,QAPGKGLEWVAR,0
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",4163,HEYTHYLQGR,10,97.6,10,435.21124,3,22.1258,218619.0,1302.6105,1.1,"",100 100 100 100 100 100 100 100 100 78,ETHCD,HEYTHYLQGR,6395
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",3647,TKESLSHFR,9,97.6,9,368.86447,3,19.3037,9257.09,1103.5724,-0.7,"",94 96 100 100 100 100 100 100 89,ETHCD,TKESLSHFR,1361
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",11907,LTWFDEGTAEFFAGSTR,17,97.6,17,967.94916,2,64.9324,66995.53,1933.8846,-0.4,"",87 87 100 100 100 100 100 100 100 100 100 100 99 99 97 96 95,HCD,LTWFDEGTAEFFAGSTR,23416
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",4773,DTLM(+15.99)LSR,7,97.6,7,426.21838,2,25.4288,115568.0,850.4218,0.4,"Oxidation (M)",95 95 100 100 100 100 94,HCD,DTLM(+15.99)LSR,5808
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",5396,GRGDSVVYGLR,11,97.6,11,393.54752,3,28.9093,54173.0,1177.6204,0.3,"",96 98 99 100 100 100 100 100 100 100 80,ETHCD,GRGDSVVYGLR,3284
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",8485,DSTYSLSSTLTSRR,14,97.5,14,787.40173,2,45.8671,242945.9,1572.7744,9.3,"",100 100 100 100 100 100 100 100 99 100 100 100 100 67,ETHCD,DSTYSLSSTLTSRR,19182
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",3880,ADSVFKGR,8,97.5,8,440.23810,2,20.6071,27530.41,878.4610,0.7,"",96 99 100 100 99 98 95 93,HCD,ADSVFKGR,6714
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",4217,HEYTHYLQAR,10,97.5,10,439.88297,3,22.4412,1841098.6,1316.6262,0.7,"",94 97 99 100 100 100 100 100 100 84,ETHCD,HEYTHYLQAR,6678"#;

const DATA_11_ALL_CANDIDATES: &str = r#""Source File","Scan","Peptide","Tag length","ALC (%)","Length","m/z","z","RT","Area","Mass","ppm","PTM","local confidence (%)","mode",tag(>=0.0%),"Feature Id"
"20240903_SA_HTA_KR40_500ng_CID_130min.wiff",16281,Q(-17.03)PEDF,5,99.9,5,618.24072,1,48.7702,59523.16,617.2333,0.3,"Pyro-glu from Q",100 100 100 100 100,CID,Q(-17.03)PEDF,16481
"20240903_SA_HTA_KR40_500ng_CID_130min.wiff",16265,Q(-17.03)PEDF,5,99.9,5,618.24072,1,48.7702,59523.16,617.2333,0.3,"Pyro-glu from Q",100 100 100 100 100,CID,Q(-17.03)PEDF,16481
"20240903_SA_HTA_KR40_500ng_CID_130min.wiff",20270,LLGGPSVF,8,99.7,8,789.44788,1,72.3085,5838.5063,788.4432,-3.3,"",100 100 100 99 99 100 100 100,CID,LLGGPSVF,33201
"20240903_SA_HTA_KR40_500ng_CID_130min.wiff",9109,LVKDY,5,99.7,5,637.35510,1,18.7355,4015.6555,636.3483,-0.7,"",100 100 100 100 100,CID,LVKDY,18803
"20240903_SA_HTA_VE20_500ng_CID_60min.wiff",15195,DFATYYCQQHYTTPPTFGQGTKVE,24,99.7,24,1391.62524,2,33.1102,34988.38,2781.2380,-0.8,"",100 98 100 100 100 100 100 99 99 100 100 100 100 100 100 100 100 100 99 100 100 100 100 100,CID,DFATYYCQQHYTTPPTFGQGTKVE,29750
"20240903_SA_HTA_KR40_500ng_CID_90min.wiff",22518,LLGGPSVF,8,99.7,8,789.44580,1,71.3410,6788.68,788.4432,-5.9,"",100 100 100 99 99 100 100 100,CID,LLGGPSVF,32106
"20240903_SA_HTA_VE40_500ng_CID_60min.wiff",15635,APGKGLEWVARLYPTNGY,18,99.6,18,996.52295,2,34.7298,2675.0488,1991.0265,2.5,"",100 99 100 100 100 100 100 100 100 100 99 100 100 99 98 99 100 99,CID,APGKGLEWVARLYPTNGY,25167
"20240903_SA_HTA_KR40_500ng_CID_60min.wiff",11103,PEVKF,5,99.6,5,619.34204,1,23.4924,14385.308,618.3377,-4.7,"",99 100 100 100 99,CID,PEVKF,9777
"20240903_SA_HTA_KR20_500ng_CID_90min.wiff",15449,LLSDF,5,99.6,5,594.31116,1,50.1211,1865.4401,593.3060,-3.7,"",100 99 100 100 99,CID,LLSDF,8217
"20240903_SA_HTA_KR20_500ng_CID_90min.wiff",10862,VDVSHE,6,99.6,6,685.31506,1,25.3097,16523.75,684.3079,-0.1,"",100 100 100 100 98 100,CID,VDVSHE,16688
"20240903_SA_HTA_KR20_500ng_CID_60min.wiff",11421,DDYAGAMHWVR,11,99.5,11,660.79303,2,27.7775,31207.42,1319.5717,-0.1,"",99 98 99 100 100 100 100 100 100 100 100,CID,DDYAGAMHWVR,10245
"20240903_SA_HTA_VE40_500ng_CID_130min.wiff",16631,LLSDF,5,99.5,5,594.31287,1,51.6116,3949.231,593.3060,-0.8,"",100 99 100 100 99,CID,LLSDF,13293
"20240903_SA_HTA_KR40_500ng_CID_90min.wiff",16720,FLYSKL,6,99.5,6,770.44299,1,46.7960,14885.409,769.4374,-2.2,"",99 100 100 100 99 99,CID,FLYSKL,30214
"20240903_SA_HTA_KR40_500ng_CID_60min.wiff",13366,GYGVHWVR,8,99.5,8,487.25134,2,27.3181,,972.4930,-5.0,"",99 99 100 100 99 100 100 100,CID,GYGVHWVR,0
"20240903_SA_HTA_KR20_500ng_CID_130min.wiff",16083,GYGVHWVR,8,99.4,8,487.25189,2,45.3716,,972.4930,-3.8,"",99 99 100 99 99 100 100 100,CID,GYGVHWVR,0
"20240903_SA_HTA_KR40_500ng_CID_90min.wiff",16846,TFSGF,5,99.4,5,558.25446,1,47.0806,6637.823,557.2485,-2.5,"",98 99 100 100 100,CID,TFSGF,8188
"20240903_SA_HTA_KR40_500ng_CID_130min.wiff",15382,GYGVHWVR,8,99.3,8,487.25345,2,44.8251,,972.4930,-0.6,"",98 99 100 99 99 100 100 100,CID,GYGVHWVR,0
"20240903_SA_HTA_KR40_500ng_CID_60min.wiff",9567,VVDVSHED,8,99.3,8,899.41119,1,20.0530,831.1227,898.4032,0.8,"",100 99 100 100 100 98 99 99,CID,VVDVSHED,26075
"20240903_SA_HTA_VE20_500ng_CID_60min.wiff",15912,TKNQVSLTCLVKGF,14,99.3,14,1537.84521,1,34.7840,1830.2261,1536.8334,3.0,"",99 98 97 98 100 100 100 100 100 100 100 100 100 100,CID,TKNQVSLTCLVKGF,29989"#;

const DATA_11_CUSTOM_MODIFICATION: &str = r#""Source File","Scan","Peptide","Tag length","ALC (%)","Length","m/z","z","RT","Area","Mass","ppm","PTM","local confidence (%)","mode",tag(>=0.0%),"Feature Id"
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",6514,YSSDEKVLGEDFSDTR,16,99.9,16,616.61450,3,35.1011,190625.72,1846.8220,-0.2,"",100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 98,HCD,YSSDEKVLGEDFSDTR,14223
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",6743,DLQMTQSPSSLSASVGDR,18,99.2,18,626.96686,3,36.1993,794748.8,1877.8789,-0.1,"",99 99 100 100 100 100 99 99 100 100 100 100 100 100 100 100 100 90,ETHCD,DLQMTQSPSSLSASVGDR,14572
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",5644,DLQM[C:Oxidation]TQSPSSLSASVGDR,18,99.0,18,947.94415,2,30.1547,548587.2,1893.8738,0.0,"Oxidation (M)",97 96 100 100 100 100 100 100 100 100 100 100 100 100 100 95 97 97,HCD,DLQM(+15.99)TQSPSSLSASVGDR,22869
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",6616,EVNSQFFR,8,98.8,8,513.75372,2,35.6492,60834.707,1025.4930,-0.1,"",98 97 99 100 99 100 100 97,HCD,EVNSQFFR,10540
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",9428,TYTFDNGTFLLR,12,98.7,12,724.36444,2,51.1627,83370.51,1446.7144,0.0,"",99 99 100 100 100 99 97 98 99 99 99 96,HCD,TYTFDNGTFLLR,17359
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",7805,EVQLVESGGGLVQPGGSLR,19,98.6,19,628.00598,3,42.2104,298160.3,1880.9955,0.3,"",98 97 99 100 100 100 100 100 100 100 100 100 98 99 99 99 100 100 83,ETHCD,EVQLVESGGGLVQPGGSLR,14598
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",5945,VLGEDFSDTR,10,98.0,10,569.77283,2,31.8568,18040.87,1137.5302,0.8,"",100 100 99 99 99 96 95 98 98 95,HCD,VLGEDFSDTR,12617
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",9483,TYTFDDGTFLLR,12,98.0,12,724.85663,2,51.3811,43949.74,1447.6984,0.3,"",96 97 100 100 100 100 96 96 97 99 99 95,HCD,TYTFDDGTFLLR,17378
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",3737,SRSGGGGNGLGSGGSLR,17,97.9,17,492.58188,3,19.8131,29103.092,1474.7236,0.1,"",100 99 100 100 100 99 98 95 100 100 100 100 98 98 99 99 83,ETHCD,SRSGGGGNGLGSGGSLR,9583
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",4286,TSGVLPR,7,97.8,7,365.21625,2,22.6110,744909.94,728.4181,-0.2,"",100 100 100 100 100 100 85,ETHCD,TSGVLPR,1074
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",6256,QAPGKGLEWVAR,12,97.8,12,437.91058,3,33.6266,,1310.7095,0.3,"",94 98 99 100 100 100 100 100 100 100 100 82,ETHCD,QAPGKGLEWVAR,0
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",4163,HEYTHYLQGR,10,97.6,10,435.21124,3,22.1258,218619.0,1302.6105,1.1,"",100 100 100 100 100 100 100 100 100 78,ETHCD,HEYTHYLQGR,6395
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",3647,TKESLSHFR,9,97.6,9,368.86447,3,19.3037,9257.09,1103.5724,-0.7,"",94 96 100 100 100 100 100 100 89,ETHCD,TKESLSHFR,1361
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",11907,LTWFDEGTAEFFAGSTR,17,97.6,17,967.94916,2,64.9324,66995.53,1933.8846,-0.4,"",87 87 100 100 100 100 100 100 100 100 100 100 99 99 97 96 95,HCD,LTWFDEGTAEFFAGSTR,23416
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",4773,DTLM[C:Oxidation]LSR,7,97.6,7,426.21838,2,25.4288,115568.0,850.4218,0.4,"Oxidation (M)",95 95 100 100 100 100 94,HCD,DTLM(+15.99)LSR,5808
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",5396,GRGDSVVYGLR,11,97.6,11,393.54752,3,28.9093,54173.0,1177.6204,0.3,"",96 98 99 100 100 100 100 100 100 100 80,ETHCD,GRGDSVVYGLR,3284
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",8485,DSTYSLSSTLTSRR,14,97.5,14,787.40173,2,45.8671,242945.9,1572.7744,9.3,"",100 100 100 100 100 100 100 100 99 100 100 100 100 67,ETHCD,DSTYSLSSTLTSRR,19182
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",3880,ADSVFKGR,8,97.5,8,440.23810,2,20.6071,27530.41,878.4610,0.7,"",96 99 100 100 99 98 95 93,HCD,ADSVFKGR,6714
"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw",4217,HEYTHYLQAR,10,97.5,10,439.88297,3,22.4412,1841098.6,1316.6262,0.7,"",94 97 99 100 100 100 100 100 100 84,ETHCD,HEYTHYLQAR,6678"#;

const DATA_12: &str = r#""Source File","Scan","Peptide","Tag length","Deep Novo Score (%)","ALC (%)","Length","m/z","z","RT","Area","Mass","ppm","PTM","local confidence (%)","mode",tag(>=0%),"Feature Id"
"9661_weiwei_mabs_novo_chymo.raw",4245;4251,TLSKADY,7,100.0,100.0,7,399.20514,2,20.9525,4.8573476E7,796.3967,-1.2,"",100 100 100 100 100 100 100,MIX,TLSKADY,1008
"9661_weiwei_mabs_novo_chymo.raw",4478;4483,SSPVTKSF,8,100.0,100.0,8,426.72626,2,27.9639,2.40281728E8,851.4388,-1.1,"",100 100 100 100 100 100 100 100,MIX,SSPVTKSF,1492
"9661_weiwei_mabs_novo_chymo.raw",3553;3556,TLSKADY,7,100.0,100.0,7,399.20514,2,20.9525,4.8573476E7,796.3967,-1.2,"",100 100 100 100 100 100 100,MIX,TLSKADY,1008
"9661_weiwei_mabs_novo_chymo.raw",5170;5177,TFRNY,5,100.0,100.0,5,350.67392,2,30.5363,8.4915336E7,699.3340,-1.0,"",100 100 100 100 100,MIX,TFRNY,5
"9661_weiwei_mabs_novo_chymo.raw",5046;5051,TFRNY,5,100.0,100.0,5,350.67392,2,30.5363,8.4915336E7,699.3340,-1.0,"",100 100 100 100 100,MIX,TFRNY,5
"9661_weiwei_mabs_novo_chymo.raw",6243;6249,SC(+58.01)SVM(+15.99)HEALH,10,100.0,100.0,10,396.49805,3,33.4147,4701863.5,1186.4747,-2.0,"Carboxymethyl; Oxidation (M)",100 100 100 100 100 100 100 100 100 100,MIX,SC(+58.01)SVM(+15.99)HEALH,931
"9661_weiwei_mabs_novo_chymo.raw",2777;2782,SGALTSGVH,9,99.9,99.9,9,414.71368,2,14.7044,6212751.5,827.4137,-1.1,"",100 100 100 100 100 100 100 100 100,MIX,SGALTSGVH,1268
"9661_weiwei_mabs_novo_chymo.raw",6831;6836,KTTPPVLDSDGSF,13,99.9,99.9,13,682.33929,2,36.6478,2.2263128E7,1362.6666,-1.9,"",99 100 100 100 100 100 100 100 100 100 100 100 100,MIX,KTTPPVLDSDGSF,5953
"9661_weiwei_mabs_novo_chymo.raw",2628;2631,YC(+58.01)QHSY,6,99.9,99.9,6,429.65768,2,12.8731,1.723695E7,857.3014,-0.7,"Carboxymethyl",100 100 100 100 100 100,MIX,YC(+58.01)QHSY,1569
"9661_weiwei_mabs_novo_chymo.raw",7151;7157,VSNKALPAPLEK,12,99.9,99.9,12,422.91812,3,37.0323,3034367.8,1265.7343,-1.4,"",100 100 100 100 100 100 100 100 100 100 100 100,MIX,VSNKALPAPLEK,1414
"9661_weiwei_mabs_novo_chymo.raw",8385;8391,SC(+58.01)SVM(+15.99)HEALHNHY,13,99.9,99.9,13,401.16647,4,44.6015,7.0465096E7,1600.6399,-1.9,"Carboxymethyl; Oxidation (M)",100 100 100 100 100 100 100 100 100 100 100 100 100,MIX,SC(+58.01)SVM(+15.99)HEALHNHY,1043
"9661_weiwei_mabs_novo_chymo.raw",10970,LAHLFSDDEKFY,12,99.9,99.9,12,495.57254,3,56.1511,3.4052484E7,1483.6984,-1.7,"",100 100 100 100 100 100 100 100 100 100 100 100,ETHCD,LAHLFSDDEKFY,2942
"9661_weiwei_mabs_novo_chymo.raw",4367;4371,SSPVTKSF,8,99.9,99.9,8,426.72595,2,22.5410,,851.4388,-1.8,"",100 100 100 100 100 100 100 100,MIX,SSPVTKSF,0
"9661_weiwei_mabs_novo_chymo.raw",10971,LAHLFSDDEKFY,12,99.9,99.9,12,495.57254,3,56.1511,3.4052484E7,1483.6984,-1.7,"",100 100 100 100 100 100 100 100 100 100 100 100,ETHCD,LAHLFSDDEKFY,2942
"9661_weiwei_mabs_novo_chymo.raw",10482;10490,LSAHNGDTNYAEKVQDRLTL,20,99.9,99.9,20,562.03442,4,53.8517,,2244.1135,-2.2,"",100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 99 100 100 100,MIX,LSAHNGDTNYAEKVQDRLTL,0
"9661_weiwei_mabs_novo_chymo.raw",10972,LAHLFSDDEKFY,12,99.9,99.9,12,495.57254,3,56.1511,3.4052484E7,1483.6984,-1.7,"",100 100 100 100 100 100 100 100 100 100 100 100,ETHCD,LAHLFSDDEKFY,2942
"9661_weiwei_mabs_novo_chymo.raw",6579;6584,AC(+58.01)EVTHQGLSSPVTK,15,99.9,99.9,15,538.93024,3,35.3031,2.7658062E7,1613.7719,-1.9,"Carboxymethyl",100 100 100 100 100 100 100 100 100 100 100 99 100 100 100,MIX,AC(+58.01)EVTHQGLSSPVTK,3776
"9661_weiwei_mabs_novo_chymo.raw",3222;3229,SGALTSGVH,9,99.9,99.9,9,414.71368,2,14.7044,6212751.5,827.4137,-1.1,"",100 100 100 100 100 100 100 100 100,MIX,SGALTSGVH,1268
"9661_weiwei_mabs_novo_chymo.raw",4947;4951,TFRNY,5,99.9,99.9,5,350.67383,2,25.7710,,699.3340,-1.3,"",100 100 100 100 100,MIX,TFRNY,0"#;
