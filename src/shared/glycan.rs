use std::ops::Range;

use crate::{
    build::{Context, CustomError},
    formula::{Chemical, MolecularFormula},
    helper_functions::*,
    Element,
};

/// All monosaccharides as required by pro forma
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(non_camel_case_types, missing_docs)]
pub enum MonoSaccharide {
    /// Any general heptose (seven carbon sugar)
    Heptose,
    phosphate,
    a_Hex,
    Sug,
    d_Hex,
    HexN,
    /// Any general pentose (five carbon sugar)
    Pentose,
    /// Any general tetrose (four carbon sugar)
    Tetrose,
    HexS,
    /// Hexose with two sulfates
    HexS2,
    HexP,
    Neu5Ac,
    Non,
    HexNAcS,
    Dec,
    en_a_Hex,
    Neu5Gc,
    Neu,
    /// Any amino sugar (hexose with N acetylated amino group)
    HexNAc,
    Deoxyhexose,
    HexNS,
    Tri,
    Oct,
    sulfate,
    /// Any general hexose (six carbon sugar)
    Hexose,
}

// TODO: think about a better way to save?
// * Base sugar (Hex, Hep, Oct, etc)
// * 'Flavour' of the base sugar (Man, Glc, etc for Hexoses)
// * Linked amino (N acetylated or not)
// * Added mods: sulfates, pyruvates, deoxy, phosphates(?), Ac, Gc

/// All monosaccharides ordered to be able to parse glycans by matching them from the top
#[allow(dead_code)]
pub const GLYCAN_PARSE_LIST: &[(&str, MonoSaccharide)] = &[
    ("phosphate", MonoSaccharide::phosphate),
    ("sulfate", MonoSaccharide::sulfate),
    ("Sug", MonoSaccharide::Sug),
    ("Tri", MonoSaccharide::Tri),
    ("HexNS", MonoSaccharide::HexNS),
    // hexose amino sugars N acetylated with sulfate TODO: ask to make sure
    ("HexNAc(S)", MonoSaccharide::HexNAcS),
    ("GlcNAc6S", MonoSaccharide::HexNAcS),
    // hexose amino sugars N acetylated
    ("HexNAc", MonoSaccharide::HexNAc),
    ("GlcNAc", MonoSaccharide::HexNAc),
    ("GalNAc", MonoSaccharide::HexNAc),
    ("ManNAc", MonoSaccharide::HexNAc),
    ("AllNAc", MonoSaccharide::HexNAc),
    ("LAltNAc", MonoSaccharide::HexNAc),
    ("GulNAc", MonoSaccharide::HexNAc),
    ("LIdoNAc", MonoSaccharide::HexNAc),
    ("TalNAc", MonoSaccharide::HexNAc),
    // hexose amino sugars
    ("HexN", MonoSaccharide::HexN),
    ("GlcN", MonoSaccharide::HexN),
    ("GalN", MonoSaccharide::HexN),
    ("ManN", MonoSaccharide::HexN),
    ("AllN", MonoSaccharide::HexN),
    ("LAlN", MonoSaccharide::HexN),
    ("GulN", MonoSaccharide::HexN),
    ("LIdN", MonoSaccharide::HexN),
    ("TalN", MonoSaccharide::HexN),
    // tetroses
    ("Tet", MonoSaccharide::Tetrose),
    ("Ery", MonoSaccharide::Tetrose),
    ("Tho", MonoSaccharide::Tetrose),
    // pentoses
    ("Pen", MonoSaccharide::Pentose),
    ("Ribf", MonoSaccharide::Pentose),
    ("Rib", MonoSaccharide::Pentose),
    ("Araf", MonoSaccharide::Pentose),
    ("Ara", MonoSaccharide::Pentose),
    ("LAraf", MonoSaccharide::Pentose),
    ("LAra", MonoSaccharide::Pentose),
    ("Xyl", MonoSaccharide::Pentose),
    ("Lyx", MonoSaccharide::Pentose),
    // hexoses with sulfate TODO: ask to make sure
    ("HexS", MonoSaccharide::HexS),
    ("Gal6S", MonoSaccharide::HexS),
    // hexoses with dual sulfate TODO: ask to make sure
    ("Gal4,6S2", MonoSaccharide::HexS2),
    // hexoses
    ("a-Hex", MonoSaccharide::a_Hex),
    ("en,a-Hex", MonoSaccharide::en_a_Hex),
    ("d-Hex", MonoSaccharide::d_Hex),
    ("HexP", MonoSaccharide::HexP),
    ("Hex", MonoSaccharide::Hexose),
    ("Glc", MonoSaccharide::Hexose),
    ("Galf", MonoSaccharide::Hexose),
    ("Gal", MonoSaccharide::Hexose),
    ("LGal", MonoSaccharide::Hexose),
    ("Man", MonoSaccharide::Hexose),
    ("All", MonoSaccharide::Hexose),
    ("LAlt", MonoSaccharide::Hexose),
    ("Gul", MonoSaccharide::Hexose),
    ("LIdo", MonoSaccharide::Hexose),
    ("Tal", MonoSaccharide::Hexose),
    // deoxyhexose
    ("Fuc", MonoSaccharide::Deoxyhexose),
    ("LFuc", MonoSaccharide::Deoxyhexose),
    ("Rha", MonoSaccharide::Deoxyhexose),
    ("LRha", MonoSaccharide::Deoxyhexose),
    ("Qui", MonoSaccharide::Deoxyhexose),
    ("2dGlc", MonoSaccharide::Deoxyhexose),
    // heptoses
    ("Hep", MonoSaccharide::Heptose),
    ("Lgro-manHep", MonoSaccharide::Heptose),
    ("gro-manHep", MonoSaccharide::Heptose),
    // other
    ("Oct", MonoSaccharide::Oct),
    ("Non", MonoSaccharide::Non),
    ("Dec", MonoSaccharide::Dec),
    ("Neu5Ac", MonoSaccharide::Neu5Ac),
    ("Neu5Gc", MonoSaccharide::Neu5Gc),
    ("Neu", MonoSaccharide::Neu),
];

impl TryFrom<&str> for MonoSaccharide {
    type Error = ();
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "Hep" => Ok(Self::Heptose),
            "phosphate" => Ok(Self::phosphate),
            "aHex" | "a-Hex" | "HexA" => Ok(Self::a_Hex),
            "Sug" => Ok(Self::Sug),
            "dHex" | "d-Hex" => Ok(Self::d_Hex),
            "HexN" => Ok(Self::HexN),
            "Pent" | "Pen" => Ok(Self::Pentose),
            "Tet" => Ok(Self::Tetrose),
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
            "Fuc" => Ok(Self::Deoxyhexose),
            "HexNS" => Ok(Self::HexNS),
            "Tri" => Ok(Self::Tri),
            "Oct" => Ok(Self::Oct),
            "Sulf" | "sulfate" => Ok(Self::sulfate),
            "Hex" => Ok(Self::Hexose),
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
                Self::Heptose => "Hep",
                Self::phosphate => "phosphate",
                Self::a_Hex => "a-Hex",
                Self::Sug => "Sug",
                Self::d_Hex => "d-Hex",
                Self::HexN => "HexN",
                Self::Pentose => "Pen",
                Self::Tetrose => "Tet",
                Self::HexS => "HexS",
                Self::HexS2 => "HexS2",
                Self::HexP => "HexP",
                Self::Neu5Ac => "Neu5Ac",
                Self::Non => "Non",
                Self::HexNAcS => "HexNAc(S)",
                Self::Dec => "Dec",
                Self::en_a_Hex => "en,a-Hex",
                Self::Neu5Gc => "Neu5Gc",
                Self::Neu => "Neu",
                Self::HexNAc => "HexNAc",
                Self::Deoxyhexose => "Fuc",
                Self::HexNS => "HexNS",
                Self::Tri => "Tri",
                Self::Oct => "Oct",
                Self::sulfate => "sulfate",
                Self::Hexose => "Hex",
            }
        )
    }
}

impl MonoSaccharide {
    fn symbol(&self) -> char {
        // ⬠◇♢▭◮⬘
        // ⬟◆♦▬
        match self {
            Self::Hexose => '○',      //●
            Self::HexNAc => '□',      // ■
            Self::Deoxyhexose => '△', // ▲
            Self::Pentose => '☆',     // ★
            Self::HexN => '⬔',
            _ => '⬡', // ⬢
        }
    }
}

impl Chemical for MonoSaccharide {
    fn formula(&self) -> MolecularFormula {
        match self {
            Self::Heptose => molecular_formula!(H 12 C 7 O 6),
            Self::phosphate => molecular_formula!(H 1 O 3 P 1),
            Self::a_Hex => molecular_formula!(H 8 C 6 O 6),
            Self::Sug => molecular_formula!(H 2 C 2 O 1),
            Self::HexN => molecular_formula!(H 11 C 6 N 1 O 4),
            Self::Pentose => molecular_formula!(H 8 C 5 O 4),
            Self::Tetrose => molecular_formula!(H 6 C 4 O 3),
            Self::HexP => molecular_formula!(H 11 C 6 O 8 P 1),
            Self::Neu5Ac => molecular_formula!(H 17 C 11 N 1 O 8),
            Self::Non => molecular_formula!(H 16 C 9 O 8),
            Self::HexNAcS => molecular_formula!(H 13 C 8 N 1 O 8 S 1),
            Self::Dec => molecular_formula!(H 18 C 10 O 9),
            Self::en_a_Hex => molecular_formula!(H 6 C 6 O 5),
            Self::Neu5Gc => molecular_formula!(H 17 C 11 N 1 O 9),
            Self::Neu => molecular_formula!(H 15 C 9 N 1 O 7),
            Self::HexNAc => molecular_formula!(H 13 C 8 N 1 O 5),
            Self::Deoxyhexose => molecular_formula!(H 10 C 6 O 4),
            Self::HexNS => molecular_formula!(H 11 C 6 N 1 O 7 S 1),
            Self::Tri => molecular_formula!(H 4 C 3 O 2),
            Self::Oct => molecular_formula!(H 14 C 8 O 7),
            Self::sulfate => molecular_formula!(O 3 S 1),
            Self::d_Hex | Self::Hexose => molecular_formula!(H 10 C 6 O 5),
            Self::HexS => molecular_formula!(H 10 C 6 O 8 S 1),
            Self::HexS2 => molecular_formula!(H 9 C 6 O 11 S 2), // TODO: ask around
        }
    }
}

/// Rose tree representation of glycan structure
#[derive(Eq, PartialEq, Clone, Hash)]
pub struct GlycanStructure {
    sugar: MonoSaccharide,
    branches: Vec<GlycanStructure>,
}

impl GlycanStructure {
    /// Parse a short IUPAC glycan structure
    /// # Panics
    /// Panics if there is no single sugar found
    /// # Errors
    /// Errors when the format is not correct, could be unknown monosaccharide, or a open brace
    pub fn from_short_iupac(
        line: &str,
        range: Range<usize>,
        line_number: usize,
    ) -> Result<Self, CustomError> {
        // GlcNAc(?1-?)Man(?1-?)[Man(?1-?)Man(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc-ol
        let mut offset = range.start;
        let mut branch: Self = Self {
            sugar: MonoSaccharide::Dec,
            branches: Vec::new(),
        }; // Starting sugar, will be removed
        let mut last_branch: &mut Self = &mut branch;
        let bytes = line.as_bytes();
        while offset < range.end {
            while bytes[offset] == b'[' {
                let end = end_of_enclosure(bytes, offset + 1, b'[', b']').ok_or_else(|| {
                    CustomError::error(
                        "Invalid IUPC short glycan",
                        "No closing brace for branch",
                        Context::line(line_number, line, offset, range.end - offset),
                    )
                })?;
                let offshoot = Self::from_short_iupac(line, offset + 1..end, line_number)?;
                last_branch.branches.push(offshoot);
                offset = end + 1;
            }
            if let Some((name, sugar)) = GLYCAN_PARSE_LIST
                .iter()
                .find(|(name, _)| line[offset..].starts_with(name))
            {
                last_branch.branches.push(Self {
                    sugar: *sugar,
                    branches: Vec::new(),
                });
                last_branch = last_branch.branches.last_mut().unwrap();
                offset += name.len();
                if bytes[offset] == b'(' {
                    if let Some(end) = next_char(bytes, offset + 1, b')') {
                        offset = end + 1; // just ignore all linking stuff I do not care
                    } else {
                        assert!(range.end - offset < 10); // make sure it is the last part
                        offset = range.end; // assume it is the last not closed brace
                    }
                }
            } else {
                // If no monosaccharide could be identified go and report an error
                return Err(CustomError::error(
                    "Invalid IUPC short glycan",
                    "Could not recognise this monosaccharide",
                    Context::line(line_number, line, offset, 1),
                ));
            }
        }
        Ok(branch.branches.pop().unwrap()) // Remove the outer starting sugar
    }
}
