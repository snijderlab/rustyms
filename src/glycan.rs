use crate::{da, HasMass, Mass, MassSystem};

include!("shared/glycan.rs");

impl HasMass for MonoSaccharide {
    fn mass<M: MassSystem>(&self) -> Mass {
        match self {
            Self::Hep => da(M::C * 7.0 + M::H * 12.0 + M::O * 6.0),
            Self::phosphate => da(M::H + M::O * 3.0 + M::P),
            Self::a_Hex => da(M::C * 6.0 + M::H * 8.0 + M::O * 6.0),
            Self::Sug => da(M::C * 2.0 + M::H * 2.0 + M::O),
            Self::HexN => da(M::C * 6.0 + M::H * 11.0 + M::N + M::O * 4.0),
            Self::Pen => da(M::C * 5.0 + M::H * 8.0 + M::O * 4.0),
            Self::Tet => da(M::C * 4.0 + M::H * 6.0 + M::O * 3.0),
            Self::HexS => da(M::C * 6.0 + M::H * 10.0 + M::O * 8.0 + M::S),
            Self::HexP => da(M::C * 6.0 + M::H * 11.0 + M::O * 8.0 + M::P),
            Self::Neu5Ac => da(M::C * 11.0 + M::H * 17.0 + M::N + M::O * 8.0),
            Self::Non => da(M::C * 9.0 + M::H * 16.0 + M::O * 8.0),
            Self::HexNAcS => da(M::C * 8.0 + M::H * 13.0 + M::N + M::O * 8.0 + M::S),
            Self::Dec => da(M::C * 10.0 + M::H * 18.0 + M::O * 9.0),
            Self::en_a_Hex => da(M::C * 6.0 + M::H * 6.0 + M::O * 5.0),
            Self::Neu5Gc => da(M::C * 11.0 + M::H * 17.0 + M::N + M::O * 9.0),
            Self::Neu => da(M::C * 9.0 + M::H * 15.0 + M::N + M::O * 7.0),
            Self::HexNAc => da(M::C * 8.0 + M::H * 13.0 + M::N + M::O * 5.0),
            Self::Fuc => da(M::C * 6.0 + M::H * 10.0 + M::O * 4.0),
            Self::HexNS => da(M::C * 6.0 + M::H * 11.0 + M::N + M::O * 7.0 + M::S),
            Self::Tri => da(M::C * 3.0 + M::H * 4.0 + M::O * 2.0),
            Self::Oct => da(M::C * 8.0 + M::H * 14.0 + M::O * 7.0),
            Self::sulfate => da(M::O * 3.0 + M::S),
            Self::d_Hex | Self::Hex => da(M::C * 6.0 + M::H * 10.0 + M::O * 5.0),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::mass::AverageWeight;

    #[test]
    #[allow(clippy::float_cmp)] // Already handled in a way clippy does not recognise
    fn mass_glycan() {
        assert_eq!(
            1445.0,
            (MonoSaccharide::Hex.mass::<AverageWeight>() * 3.0
                + MonoSaccharide::HexNAc.mass::<AverageWeight>() * 4.0
                + MonoSaccharide::Fuc.mass::<AverageWeight>())
            .value
            .round()
        );
    }
}
