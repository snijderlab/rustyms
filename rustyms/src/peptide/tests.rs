#![allow(clippy::missing_panics_doc)]
use std::num::NonZeroU16;

use tests::parse::{global_modifications, parse_charge_state};

use crate::{modification::SimpleModification, system::da, AminoAcid, Element, MolecularCharge};

use super::*;

#[test]
fn parse_global_modifications() {
    let parse = |str: &str| global_modifications(str.as_bytes(), 0, str, None);
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
