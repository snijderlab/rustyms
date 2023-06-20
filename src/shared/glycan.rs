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
#[allow(dead_code)]
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

impl TryFrom<&str> for MonoSaccharide {
    type Error = ();
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "Hep" => Ok(Self::Hep),
            "phosphate" => Ok(Self::phosphate),
            "aHex" | "a-Hex" | "HexA" => Ok(Self::a_Hex),
            "Sug" => Ok(Self::Sug),
            "dHex" | "d-Hex" => Ok(Self::d_Hex),
            "HexN" => Ok(Self::HexN),
            "Pent" | "Pen" => Ok(Self::Pen),
            "Tet" => Ok(Self::Tet),
            "HexS" => Ok(Self::HexS),
            "HexP" => Ok(Self::HexP),
            "NeuAc" | "Neu5Ac" => Ok(Self::Neu5Ac),
            "Non" => Ok(Self::Non),
            "HexNAc(S)" => Ok(Self::HexNAcS),
            "Dec" => Ok(Self::Dec),
            "en,a-Hex" => Ok(Self::en_a_Hex),
            "NeuGc" | "Neu5Gc" => Ok(Self::Neu5Gc),
            "Neu" => Ok(Self::Neu),
            "HexNAc" => Ok(Self::HexNAc),
            "Fuc" => Ok(Self::Fuc),
            "HexNS" => Ok(Self::HexNS),
            "Tri" => Ok(Self::Tri),
            "Oct" => Ok(Self::Oct),
            "Sulf" | "sulfate" => Ok(Self::sulfate),
            "Hex" => Ok(Self::Hex),
            _ => Err(()),
        }
    }
}

impl std::fmt::Display for MonoSaccharide {
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
