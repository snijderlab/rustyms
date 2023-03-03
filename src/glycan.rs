use crate::{da, Mass, MassSystem};

pub enum MonoSaccharide {
    Hep,
    phosphate,
    a_Hex,
    Sug,
    d_Hex,
    HexN,
    Pen,
    Tet,
    HexS,
    HexP,
    Neu5Ac,
    Non,
    HexNAcS,
    Dec,
    en_a_Hex,
    Neu5Gc,
    Neu,
    HexNAc,
    Fuc,
    HexNS,
    Tri,
    Oct,
    sulfate,
    Hex,
}

impl MonoSaccharide {
    pub fn mass<M: MassSystem>(&self) -> Mass {
        match self {
            Self::Hep => da(M::C * 7.0 + M::H * 12.0 + M::O * 6.0),
            Self::phosphate => da(M::H + M::O * 3.0 + M::P),
            Self::a_Hex => da(M::C * 6.0 + M::H * 8.0 + M::O * 6.0),
            Self::Sug => da(M::C * 2.0 + M::H * 2.0 + M::O),
            Self::d_Hex => da(M::C * 6.0 + M::H * 10.0 + M::O * 4.0),
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
            Self::Hex => da(M::C * 6.0 + M::H * 10.0 + M::O * 5.0),
        }
    }
}

impl TryFrom<&str> for MonoSaccharide {
    type Error = ();
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "Hep" => Ok(Self::Hep),
            "phosphate" => Ok(Self::phosphate),
            "a-Hex" => Ok(Self::a_Hex),
            "Sug" => Ok(Self::Sug),
            "d-Hex" => Ok(Self::d_Hex),
            "HexN" => Ok(Self::HexN),
            "Pen" => Ok(Self::Pen),
            "Tet" => Ok(Self::Tet),
            "HexS" => Ok(Self::HexS),
            "HexP" => Ok(Self::HexP),
            "Neu5Ac" => Ok(Self::Neu5Ac),
            "Non" => Ok(Self::Non),
            "HexNAc(S)" => Ok(Self::HexNAcS),
            "Dec" => Ok(Self::Dec),
            "en,a-Hex" => Ok(Self::en_a_Hex),
            "Neu5Gc" => Ok(Self::Neu5Gc),
            "Neu" => Ok(Self::Neu),
            "HexNAc" => Ok(Self::HexNAc),
            "Fuc" => Ok(Self::Fuc),
            "HexNS" => Ok(Self::HexNS),
            "Tri" => Ok(Self::Tri),
            "Oct" => Ok(Self::Oct),
            "sulfate" => Ok(Self::sulfate),
            "Hex" => Ok(Self::Hex),
            _ => Err(()),
        }
    }
}
