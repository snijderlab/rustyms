#![allow(clippy::missing_panics_doc)]
mod fuzz_crash;
mod fuzz_hang;
mod parse;
mod pro_forma_negative;
mod pro_forma_positive;
mod sloppy;

#[macro_export]
macro_rules! parse_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::CompoundPeptidoform::pro_forma($case, None);
            let res_upper =
                $crate::CompoundPeptidoform::pro_forma(&$case.to_ascii_uppercase(), None);
            let res_lower =
                $crate::CompoundPeptidoform::pro_forma(&$case.to_ascii_lowercase(), None);
            println!("{}", $case);
            dbg!(&res);
            assert!(res.is_ok());
            assert_eq!(res, res_upper);
            assert_eq!(res, res_lower);
            let back = res.as_ref().unwrap().to_string();
            let res_back = $crate::CompoundPeptidoform::pro_forma(&back, None);
            assert_eq!(res, res_back, "{} != {back}", $case);
        }
    };
    (single $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::CompoundPeptidoform::pro_forma($case, None);
            println!("{}\n{:?}", $case, res);
            assert!(res.is_ok());
        }
    };
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::CompoundPeptidoform::pro_forma($case, None);
            println!("{}\n{:?}", $case, res);
            assert!(res.is_err());
        }
    };
}
