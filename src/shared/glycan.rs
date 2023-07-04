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

include!("formula_macro.rs");

impl crate::Chemical for MonoSaccharide {
    fn formula(&self) -> crate::MolecularFormula {
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
