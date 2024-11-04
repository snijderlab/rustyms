#![allow(clippy::missing_panics_doc)]

use crate::{
    align::{align, scoring::AlignScoring, AlignType, Alignment},
    LinearPeptide, SimpleLinear,
};

#[test]
fn global_b() {
    test_alignment(
        "AEDTAVYYC[Carboxymethyl]SRWGGDGFYA",
        "AEDTAVYYM[oxidation]QSW[oxidation]SR",
        AlignScoring::default(),
        AlignType::GLOBAL_B,
        "8=1X3I2=",
    );
}

#[test]
fn global_a() {
    test_alignment(
        "HJKASFHLKJFHAS",
        "JHSLAFKJHFLKSFJ",
        AlignScoring::default(),
        AlignType::GLOBAL_A,
        "1=2X1=2X2i3i3X",
    );
}

#[test]
fn global_a_2() {
    test_alignment(
        "AAASSS",
        "ASSS",
        AlignScoring::default(),
        AlignType::GLOBAL_A,
        "2D4=",
    );
}

#[test]
fn either_global() {
    test_alignment(
        "HHHHHHAA",
        "AAHHHHHHH",
        AlignScoring::default(),
        AlignType::EITHER_GLOBAL,
        "6=",
    );
}

#[test]
fn normal() {
    test_alignment(
        "ANNA",
        "AGGGGA",
        AlignScoring {
            mismatch: 12,
            ..Default::default()
        },
        AlignType::LOCAL,
        "4X",
    );
}

/// Tests that are known to be failing, but will have to be addressed at some point in the future
#[cfg(not(github_action))]
mod known_failing {
    use super::*;

    #[test]
    fn overextended_rotation() {
        test_alignment(
            "IVQEVS",
            "LEVQVES",
            AlignScoring::default(),
            AlignType::GLOBAL,
            "1i1I2=2r1=",
        );
    }

    #[test]
    fn no_detected_rotation() {
        test_alignment(
            "AGGHT",
            "ANTH",
            AlignScoring::default(),
            AlignType::GLOBAL,
            "1=2:1i2r",
        );
    }

    #[test]
    fn no_detected_rotation_2() {
        test_alignment(
            "AGGHTK",
            "ANTHK",
            AlignScoring::default(),
            AlignType::GLOBAL,
            "1=2:1i2r1=",
        );
    }
}

/// Test if the given alignment is as expected and can be recreated
/// # Errors
/// When the alignment is not identical to path and when the alignment cannot be recreated from the path.
fn test_alignment(
    seq_a: &str,
    seq_b: &str,
    scoring: AlignScoring<'_>,
    align_type: AlignType,
    path: &str,
) {
    const MAXIMAL_STEP: u16 = 4;
    let first_peptide = LinearPeptide::pro_forma(seq_a, None)
        .unwrap()
        .into_simple_linear()
        .unwrap();
    let second_peptide = LinearPeptide::pro_forma(seq_b, None)
        .unwrap()
        .into_simple_linear()
        .unwrap();
    let alignment = align::<MAXIMAL_STEP, SimpleLinear, SimpleLinear>(
        &first_peptide,
        &second_peptide,
        scoring,
        align_type,
    );
    assert_eq!(
        alignment.short(),
        path,
        "Alignment of {seq_a} vs {seq_b} did not produce the expected alignment"
    );

    assert_eq!(
        Alignment::create_from_path(
            &first_peptide,
            &second_peptide,
            alignment.start_a,
            alignment.start_b,
            &alignment.short(),
            scoring,
            align_type,
            MAXIMAL_STEP
        )
        .unwrap(),
        alignment
    );
}
