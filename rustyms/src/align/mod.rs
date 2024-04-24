//! Code to make alignments of two peptides based on mass mistakes, and genetic information.
//!
//! A mass based alignment handles the case in which multiple amino acids are wrong, but the total mass
//! of this set of amino acids is equal to the mass of a set of different amino acids on the other peptide.
//! This is quite common in mass spectrometry where mistakes based on mass coincidences are very common.
//! For example `N` has the same mass as `GG`, so if we want to make a mass spectrometry faithful alignment
//! of `ANA` with `AGGA` the result should reflect this fact:
//!
//! ```text
//! Identity: 0.500 (2/4), Similarity: 0.750 (3/4), Gaps: 0.000 (0/4), Score: 0.706 (12/17),
//! Equal mass, Tolerance: 10 ppm, Alignment: global
//! Start: A 0 B 0, Path: 1=1:2i1=
//!
//! AN·A A
//! AGGA B
//!  ╶╴
//! ```
//! _Generated using this algorithm bound to a cli tool: <https://github.com/snijderlab/align-cli>_
//! ```rust
//! use rustyms::{*, align::*, peptide_complexity::*};
//! let a = LinearPeptide::pro_forma("ANA").unwrap().simple().unwrap();
//! let b = LinearPeptide::pro_forma("AGGA").unwrap().simple().unwrap();
//! let alignment = align::<4, Simple, Simple>(&a, &b, &matrix::BLOSUM62,
//!                    Tolerance::new_ppm(10.0), AlignType::GLOBAL);
//! assert_eq!(alignment.short(), "1=1:2i1=");
//! assert_eq!(alignment.ppm().value, 0.0);
//! ```

mod align_type;
mod alignment;
mod diagonal_array;
mod itertools_extension;
mod mass_alignment;
mod piece;
mod scoring;

#[cfg(feature = "imgt")]
mod consecutive;
#[cfg(feature = "imgt")]
pub use consecutive::*;

pub use align_type::{AlignType, Side};
pub use alignment::{Alignment, OwnedAlignment, RefAlignment, Score, Stats};
pub use mass_alignment::align;
pub use piece::Piece;
pub use scoring::MatchType;

/// Different scoring matrices that can be used.
/// Matrices from: <https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/util/tables/> and <https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/>
///
/// The UO columns are added, for these the B/J/Z score is the rounded down average of the corresponding non ambiguous AAs. All UO scores are exactly the same for all matrices (except identity).
pub mod matrix {
    use super::scoring;
    pub use scoring::matrices::*;
}

use alignment::AlignmentInner;

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use crate::{peptide_complexity::Simple, LinearPeptide};

    use super::{AlignType, Alignment, RefAlignment};

    fn align<'a, const STEPS: u16>(
        a: &'a LinearPeptide<Simple>,
        b: &'a LinearPeptide<Simple>,
    ) -> RefAlignment<'a, Simple, Simple> {
        super::align::<STEPS, Simple, Simple>(
            a,
            b,
            super::matrix::BLOSUM62,
            crate::Tolerance::new_ppm(10.0),
            AlignType::GLOBAL,
        )
    }

    fn linear(aa: &str) -> LinearPeptide<Simple> {
        LinearPeptide::pro_forma(aa).unwrap().simple().unwrap()
    }

    #[test]
    fn simple_1() {
        let a = linear("ANGARS");
        let b = linear("AGGQRS");
        let c = dbg!(align::<1>(&a, &b));
        assert_eq!(c.short(), "1=1X1=1X2=");
    }

    #[test]
    fn simple_4() {
        let a = linear("ANGARS");
        let b = linear("AGGQRS");
        let c = dbg!(align::<4>(&a, &b));
        assert_eq!(c.short(), "1=1:2i2:1i2=");
    }

    #[test]
    fn simple_unbounded() {
        let a = linear("ANGARS");
        let b = linear("AGGQRS");
        let c = dbg!(align::<{ u16::MAX }>(&a, &b));
        assert_eq!(c.short(), "1=1:2i2:1i2=");
    }
}
