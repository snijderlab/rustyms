use std::ops::Range;

use crate::{
    formula::{Chemical, MolecularFormula},
    helper_functions::*,
    Context, CustomError, Element,
};

/// All monosaccharides as required by pro forma
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(non_camel_case_types, missing_docs)]
pub enum ProFormaMonoSaccharide {
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

/// A monosaccharide with all its complexity
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MonoSaccharide {
    pub base_sugar: BaseSugar,
    pub substituents: Vec<GlycanSubstituent>,
    pub pro_forma_name: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum BaseSugar {
    /// 2 carbon base sugar
    Sugar,
    /// 3 carbon base sugar
    Triose,
    /// 4 carbon base sugar
    Tetrose(Option<TetroseIsomer>),
    /// 5 carbon base sugar
    Pentose(Option<PentoseIsomer>),
    /// 6 carbon base sugar
    Hexose(Option<HexoseIsomer>),
    /// 7 carbon base sugar
    Heptose(Option<HeptoseIsomer>),
    /// 8 carbon base sugar
    Octose,
    /// 9 carbon base sugar
    Nonose,
    /// 10 carbon base sugar
    Decose,
}

impl Chemical for BaseSugar {
    fn formula(&self) -> MolecularFormula {
        match self {
            Self::Sugar => molecular_formula!(H 2 C 2 O 1),
            Self::Triose => molecular_formula!(H 4 C 3 O 2),
            Self::Tetrose(_) => molecular_formula!(H 6 C 4 O 3),
            Self::Pentose(_) => molecular_formula!(H 8 C 5 O 4),
            Self::Hexose(_) => molecular_formula!(H 10 C 6 O 5),
            Self::Heptose(_) => molecular_formula!(H 12 C 7 O 6),
            Self::Octose => molecular_formula!(H 14 C 8 O 7),
            Self::Nonose => molecular_formula!(H 16 C 9 O 8),
            Self::Decose => molecular_formula!(H 18 C 10 O 9),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum TetroseIsomer {
    /// Ery
    Erythrose,
    /// Tho
    Threose,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum PentoseIsomer {
    /// Ribf
    Ribofuranose,
    /// Rib
    Ribopyranose,
    /// Araf
    Arabinofuranose,
    /// Ara
    Arabinopyranose,
    /// Xyl
    Xylopyranose,
    /// Lyx
    Lyxopyranose,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum HexoseIsomer {
    /// glc
    Glucose,
    /// Galf
    Galactofuranose,
    /// Gal
    Galactose,
    /// Man
    Mannose,
    /// All
    Allose,
    /// Alt
    Altrose,
    /// Gul
    Gulose,
    /// Ido
    Idose,
    /// Tal
    Talose,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum HeptoseIsomer {
    /// gro-manHep
    GlyceroMannoHeptopyranose,
}

/// Any substituent on a monosaccharide.
/// Source: https://www.ncbi.nlm.nih.gov/glycans/snfg.html table 3.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum GlycanSubstituent {
    ///Ac acetyl
    Acetyl,
    ///Ala D-alanyl
    Alanyl,
    ///Ala2Ac N-acetyl-D-alanyl
    AcetylAlanyl,
    ///Am N-acetimidoyl
    Acetimidoyl,
    ///AmMe N-(N-methyl-acetimidoyl)
    MethylAcetimidoyl,
    ///AmMe2 N-(N,N-dimethyl-acetimidoyl)
    DiMethylAcetimidoyl,
    ///Fo formyl
    Formyl,
    ///Gc glycolyl
    Glycolyl,
    ///Gln2Ac N-acetyl-glutaminyl
    AcetylGlutaminyl,
    ///5Glu2Me N-methyl-5-glutamyl
    MethylGlutamyl,
    ///Gly glycyl
    Glycyl,
    ///Gr glyceryl
    Glyceryl,
    ///Gr2,3Me2 2,3-di-O-methyl-glyceryl
    DiMethylGlyceryl,
    ///4Hb 4-hydroxybutyryl, 3RHb (R)-3-hydroxybutyryl, 3SHb (S)-3-hydroxybutyryl
    HydroxyButyryl,
    ///3,4Hb 3,4-dihydroxybutyryl
    DiHydroxyButyryl,
    ///Lt lactyl
    Lactyl,
    ///Me methyl
    Methyl,
    ///N amino
    Amino,
    ///NAc N-acetyl
    NAcetyl,
    ///P phosphate
    Phosphate,
    ///Py pyruvyl
    Pyruvyl,
    ///Pyr 1-carboxyethylidene
    CargoxyEthylidene,
    ///S sulfate
    Sulfate,
    ///Tau tauryl
    Tauryl,
}

impl Chemical for GlycanSubstituent {
    fn formula(&self) -> MolecularFormula {
        let side = match self {
            Self::Acetyl => molecular_formula!(H 3 C 2 O 1),
            Self::Alanyl => molecular_formula!(H 6 C 3 N 1 O 1),
            Self::AcetylAlanyl => molecular_formula!(H 8 C 5 N 1 O 2),
            Self::Acetimidoyl => molecular_formula!(H 5 C 2 N 1),
            Self::MethylAcetimidoyl => molecular_formula!(H 7 C 3 N 1),
            Self::DiMethylAcetimidoyl => molecular_formula!(H 9 C 4 N 1),
            Self::Formyl => molecular_formula!(H 1 C 1 O 1),
            Self::Glycolyl => molecular_formula!(H 3 C 2 O 2),
            Self::AcetylGlutaminyl => molecular_formula!(H 11 C 7 N 2 O 3),
            Self::MethylGlutamyl => molecular_formula!(H 10 C 6 N 1 O 3),
            Self::Glycyl => molecular_formula!(H 4 C 2 N 1 O 1),
            Self::Glyceryl => molecular_formula!(H 5 C 3 O 3),
            Self::DiMethylGlyceryl => molecular_formula!(H 9 C 5 O 3),
            Self::HydroxyButyryl => molecular_formula!(H 7 C 4 O 2),
            Self::DiHydroxyButyryl => molecular_formula!(H 7 C 4 O 3),
            Self::Lactyl => molecular_formula!(H 5 C 3 O 2),
            Self::Methyl => molecular_formula!(H 3 C 1),
            Self::Amino => molecular_formula!(H 2 N 1),
            Self::NAcetyl => molecular_formula!(H 4 C 2 N 1 O 1),
            Self::Phosphate => molecular_formula!(H 2 O 4 P 1),
            Self::Pyruvyl => molecular_formula!(H 3 C 3 O 2),
            Self::CargoxyEthylidene => molecular_formula!(H 3 C 3 O 3), // TODO: double substituent?
            Self::Sulfate => molecular_formula!(H 1 O 4 S 1),
            Self::Tauryl => molecular_formula!(H 6 C 2 N 1 O 3 S 1),
        };
        side - molecular_formula!(O 1 H 1) // substituent so replaces a standard oxygen side chain
    }
}

// TODO: think about a better way to save?
// * Base sugar (Hex, Hep, Oct, etc)
// * 'Flavour' of the base sugar (Man, Glc, etc for Hexoses)
// * Linked amino (N acetylated or not)
// * Added mods: sulfates, pyruvates, deoxy, phosphates(?), Ac, Gc

/// All monosaccharides ordered to be able to parse glycans by matching them from the top
#[allow(dead_code)]
pub const GLYCAN_PARSE_LIST: &[(&str, ProFormaMonoSaccharide)] = &[
    ("phosphate", ProFormaMonoSaccharide::phosphate),
    ("sulfate", ProFormaMonoSaccharide::sulfate),
    ("Sug", ProFormaMonoSaccharide::Sug),
    ("Tri", ProFormaMonoSaccharide::Tri),
    ("HexNS", ProFormaMonoSaccharide::HexNS),
    // hexose amino sugars N acetylated with sulfate TODO: ask to make sure
    ("HexNAc(S)", ProFormaMonoSaccharide::HexNAcS),
    ("GlcNAc6S", ProFormaMonoSaccharide::HexNAcS),
    // hexose amino sugars N acetylated
    ("HexNAc", ProFormaMonoSaccharide::HexNAc),
    ("GlcNAc", ProFormaMonoSaccharide::HexNAc),
    ("GalNAc", ProFormaMonoSaccharide::HexNAc),
    ("ManNAc", ProFormaMonoSaccharide::HexNAc),
    ("AllNAc", ProFormaMonoSaccharide::HexNAc),
    ("LAltNAc", ProFormaMonoSaccharide::HexNAc),
    ("GulNAc", ProFormaMonoSaccharide::HexNAc),
    ("LIdoNAc", ProFormaMonoSaccharide::HexNAc),
    ("TalNAc", ProFormaMonoSaccharide::HexNAc),
    // hexose amino sugars
    ("HexN", ProFormaMonoSaccharide::HexN),
    ("GlcN", ProFormaMonoSaccharide::HexN),
    ("GalN", ProFormaMonoSaccharide::HexN),
    ("ManN", ProFormaMonoSaccharide::HexN),
    ("AllN", ProFormaMonoSaccharide::HexN),
    ("LAlN", ProFormaMonoSaccharide::HexN),
    ("GulN", ProFormaMonoSaccharide::HexN),
    ("LIdN", ProFormaMonoSaccharide::HexN),
    ("TalN", ProFormaMonoSaccharide::HexN),
    // tetroses
    ("Tet", ProFormaMonoSaccharide::Tetrose),
    ("Ery", ProFormaMonoSaccharide::Tetrose),
    ("Tho", ProFormaMonoSaccharide::Tetrose),
    // pentoses
    ("Pen", ProFormaMonoSaccharide::Pentose),
    ("Ribf", ProFormaMonoSaccharide::Pentose),
    ("Rib", ProFormaMonoSaccharide::Pentose),
    ("Araf", ProFormaMonoSaccharide::Pentose),
    ("Ara", ProFormaMonoSaccharide::Pentose),
    ("LAraf", ProFormaMonoSaccharide::Pentose),
    ("LAra", ProFormaMonoSaccharide::Pentose),
    ("Xyl", ProFormaMonoSaccharide::Pentose),
    ("Lyx", ProFormaMonoSaccharide::Pentose),
    // hexoses with sulfate TODO: ask to make sure
    ("HexS", ProFormaMonoSaccharide::HexS),
    ("Gal6S", ProFormaMonoSaccharide::HexS),
    // hexoses with dual sulfate TODO: ask to make sure
    ("Gal4,6S2", ProFormaMonoSaccharide::HexS2),
    // hexoses
    ("a-Hex", ProFormaMonoSaccharide::a_Hex),
    ("en,a-Hex", ProFormaMonoSaccharide::en_a_Hex),
    ("d-Hex", ProFormaMonoSaccharide::d_Hex),
    ("HexP", ProFormaMonoSaccharide::HexP),
    ("Hex", ProFormaMonoSaccharide::Hexose),
    ("Glc", ProFormaMonoSaccharide::Hexose),
    ("Galf", ProFormaMonoSaccharide::Hexose),
    ("Gal", ProFormaMonoSaccharide::Hexose),
    ("LGal", ProFormaMonoSaccharide::Hexose),
    ("Man", ProFormaMonoSaccharide::Hexose),
    ("All", ProFormaMonoSaccharide::Hexose),
    ("LAlt", ProFormaMonoSaccharide::Hexose),
    ("Gul", ProFormaMonoSaccharide::Hexose),
    ("LIdo", ProFormaMonoSaccharide::Hexose),
    ("Tal", ProFormaMonoSaccharide::Hexose),
    // deoxyhexose
    ("Fuc", ProFormaMonoSaccharide::Deoxyhexose),
    ("LFuc", ProFormaMonoSaccharide::Deoxyhexose),
    ("Rha", ProFormaMonoSaccharide::Deoxyhexose),
    ("LRha", ProFormaMonoSaccharide::Deoxyhexose),
    ("Qui", ProFormaMonoSaccharide::Deoxyhexose),
    ("2dGlc", ProFormaMonoSaccharide::Deoxyhexose),
    // heptoses
    ("Hep", ProFormaMonoSaccharide::Heptose),
    ("Lgro-manHep", ProFormaMonoSaccharide::Heptose),
    ("gro-manHep", ProFormaMonoSaccharide::Heptose),
    // other
    ("Oct", ProFormaMonoSaccharide::Oct),
    ("Non", ProFormaMonoSaccharide::Non),
    ("Dec", ProFormaMonoSaccharide::Dec),
    ("Neu5Ac", ProFormaMonoSaccharide::Neu5Ac),
    ("Neu5Gc", ProFormaMonoSaccharide::Neu5Gc),
    ("Neu", ProFormaMonoSaccharide::Neu),
];

impl TryFrom<&str> for ProFormaMonoSaccharide {
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

impl std::fmt::Display for ProFormaMonoSaccharide {
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

impl ProFormaMonoSaccharide {
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

impl Chemical for ProFormaMonoSaccharide {
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
    sugar: ProFormaMonoSaccharide,
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
            sugar: ProFormaMonoSaccharide::Dec,
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
