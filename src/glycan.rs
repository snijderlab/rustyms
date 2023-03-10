use std::fmt::Display;

use crate::{da, HasMass, Mass, MassSystem};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
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

/// All monosaccharides ordered to be able to parse glycans by matching them from the top
pub const GLYCAN_PARSE_LIST: &[(&str, MonoSaccharide)] = &[
    ("phosphate", MonoSaccharide::phosphate),
    ("sulfate", MonoSaccharide::sulfate),
    ("Sug", MonoSaccharide::Sug),
    ("Tri", MonoSaccharide::Tri),
    ("Tet", MonoSaccharide::Tet),
    ("Pen", MonoSaccharide::Pen),
    ("a-Hex", MonoSaccharide::a_Hex),
    ("en,a-Hex", MonoSaccharide::en_a_Hex),
    ("d-Hex", MonoSaccharide::d_Hex),
    ("HexNAc(S)", MonoSaccharide::HexNAcS),
    ("HexNAc", MonoSaccharide::HexNAc),
    ("HexNS", MonoSaccharide::HexNS),
    ("HexN", MonoSaccharide::HexN),
    ("HexS", MonoSaccharide::HexS),
    ("HexP", MonoSaccharide::HexP),
    ("Hex", MonoSaccharide::Hex),
    ("Hep", MonoSaccharide::Hep),
    ("Oct", MonoSaccharide::Oct),
    ("Non", MonoSaccharide::Non),
    ("Dec", MonoSaccharide::Dec),
    ("Neu5Ac", MonoSaccharide::Neu5Ac),
    ("Neu5Gc", MonoSaccharide::Neu5Gc),
    ("Neu", MonoSaccharide::Neu),
    ("Fuc", MonoSaccharide::Fuc),
];

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

impl Display for MonoSaccharide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Hep => "Hep",
                Self::phosphate => "phosphate",
                Self::a_Hex => "a-Hex",
                Self::Sug => "Sug",
                Self::d_Hex => "d-Hex",
                Self::HexN => "HexN",
                Self::Pen => "Pen",
                Self::Tet => "Tet",
                Self::HexS => "HexS",
                Self::HexP => "HexP",
                Self::Neu5Ac => "Neu5Ac",
                Self::Non => "Non",
                Self::HexNAcS => "HexNAc(S)",
                Self::Dec => "Dec",
                Self::en_a_Hex => "en,a-Hex",
                Self::Neu5Gc => "Neu5Gc",
                Self::Neu => "Neu",
                Self::HexNAc => "HexNAc",
                Self::Fuc => "Fuc",
                Self::HexNS => "HexNS",
                Self::Tri => "Tri",
                Self::Oct => "Oct",
                Self::sulfate => "sulfate",
                Self::Hex => "Hex",
            }
        )
    }
}
