use crate::{Chemical, MolecularFormula};

include!("shared/glycan.rs");

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
