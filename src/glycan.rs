use crate::{Chemical, MolecularFormula};

include!("shared/glycan.rs");

impl Chemical for MonoSaccharide {
    fn formula(&self) -> MolecularFormula {
        match self {
            Self::Hep => molecular_formula!(H 12 C 7 O 6),
            Self::phosphate => molecular_formula!(H 1 O 3 P 1),
            Self::a_Hex => molecular_formula!(H 8 C 6 O 6),
            Self::Sug => molecular_formula!(H 2 C 2 O 1),
            Self::HexN => molecular_formula!(H 11 C 6 N 1 O 4),
            Self::Pen => molecular_formula!(H 8 C 5 O 4),
            Self::Tet => molecular_formula!(H 6 C 4 O 3),
            Self::HexS => molecular_formula!(H 10 C 6 O 8 S 1),
            Self::HexP => molecular_formula!(H 11 C 6 O 8 P 1),
            Self::Neu5Ac => molecular_formula!(H 17 C 11 N 1 O 8),
            Self::Non => molecular_formula!(H 16 C 9 O 8),
            Self::HexNAcS => molecular_formula!(H 13 C 8 N 1 O 8 S 1),
            Self::Dec => molecular_formula!(H 18 C 10 O 9),
            Self::en_a_Hex => molecular_formula!(H 6 C 6 O 5),
            Self::Neu5Gc => molecular_formula!(H 17 C 11 N 1 O 9),
            Self::Neu => molecular_formula!(H 15 C 9 N 1 O 7),
            Self::HexNAc => molecular_formula!(H 13 C 8 N 1 O 5),
            Self::Fuc => molecular_formula!(H 10 C 6 O 4),
            Self::HexNS => molecular_formula!(H 11 C 6 N 1 O 7 S 1),
            Self::Tri => molecular_formula!(H 4 C 3 O 2),
            Self::Oct => molecular_formula!(H 14 C 8 O 7),
            Self::sulfate => molecular_formula!(O 3 S 1),
            Self::d_Hex | Self::Hex => molecular_formula!(H 10 C 6 O 5),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    #[allow(clippy::float_cmp)] // Already handled in a way clippy does not recognise
    fn mass_glycan() {
        assert_eq!(
            1445.0,
            (MonoSaccharide::Hex.formula().average_weight().unwrap() * 3.0
                + MonoSaccharide::HexNAc.formula().average_weight().unwrap() * 4.0
                + MonoSaccharide::Fuc.formula().average_weight().unwrap())
            .value
            .round()
        );
    }
}
