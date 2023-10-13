use std::{fmt::Display, hash::Hash, ops::Range};

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
#[derive(Debug, Clone, Eq)]
pub enum MonoSaccharide {
    Predefined(PredefinedMonosaccharide),
    Allocated(AllocatedMonosaccharide),
}

impl PartialEq for MonoSaccharide {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Allocated(a), Self::Predefined(p)) => {
                p.base_sugar
                    .as_ref()
                    .map(|s| *s == a.base_sugar)
                    .unwrap_or_default()
                    && p.substituents == a.substituents
            }
            (Self::Predefined(p), Self::Allocated(a)) => {
                p.base_sugar
                    .as_ref()
                    .map(|s| *s == a.base_sugar)
                    .unwrap_or_default()
                    && p.substituents == a.substituents
            }
            (Self::Allocated(a), Self::Allocated(b)) => a == b,
            (Self::Predefined(a), Self::Predefined(b)) => a == b,
        }
    }
}

impl Hash for MonoSaccharide {
    /// Hash implementation to try and give a predefined mono sacc with the same composition as a allocated mono sacc the same value
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            Self::Allocated(a) => {
                a.base_sugar.hash(state);
                a.substituents.hash(state);
            }
            Self::Predefined(p) => {
                if let Some(s) = &p.base_sugar {
                    s.hash(state);
                }
                p.substituents.hash(state);
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PredefinedMonosaccharide {
    base_sugar: Option<BaseSugar>,
    substituents: &'static [GlycanSubstituent],
    pro_forma_name: &'static str,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AllocatedMonosaccharide {
    base_sugar: BaseSugar,
    substituents: Vec<GlycanSubstituent>,
}

impl MonoSaccharide {
    pub fn new(sugar: BaseSugar, substituents: &[GlycanSubstituent]) -> Self {
        Self::Allocated(AllocatedMonosaccharide {
            base_sugar: sugar,
            substituents: substituents.to_owned(),
        })
    }
}

impl Chemical for MonoSaccharide {
    fn formula(&self) -> MolecularFormula {
        match self {
            Self::Allocated(AllocatedMonosaccharide {
                base_sugar,
                substituents,
            }) => base_sugar.formula() + substituents.formula(),
            Self::Predefined(PredefinedMonosaccharide {
                base_sugar,
                substituents,
                ..
            }) => {
                base_sugar.as_ref().map(|s| s.formula()).unwrap_or_default()
                    + substituents.formula()
            }
        }
    }
}

impl Display for MonoSaccharide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Allocated(AllocatedMonosaccharide {
                    base_sugar,
                    substituents,
                }) => format!(
                    "{}{}",
                    base_sugar,
                    substituents
                        .iter()
                        .map(ToString::to_string)
                        .collect::<String>()
                ),
                Self::Predefined(PredefinedMonosaccharide { pro_forma_name, .. }) =>
                    pro_forma_name.to_string(),
            }
        )
    }
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

impl Display for BaseSugar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Sugar => "Sug",
                Self::Triose => "Tri",
                Self::Tetrose(_) => "Tet",
                Self::Pentose(_) => "Pen",
                Self::Hexose(_) => "Hex",
                Self::Heptose(_) => "Hep",
                Self::Octose => "Oct",
                Self::Nonose => "Non",
                Self::Decose => "Dec",
            }
        )
    }
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
    ///A acid
    Acid,
    ///Ala D-alanyl
    Alanyl,
    ///Ala2Ac N-acetyl-D-alanyl
    AcetylAlanyl,
    ///Am N-acetimidoyl
    Acetimidoyl,
    ///AmMe N-(N-methyl-acetimidoyl)
    MethylAcetimidoyl,
    ///d Deoxy
    Deoxy,
    ///en didehydro an addition of a double bond
    Didehydro,
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

impl Display for GlycanSubstituent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Acetyl => "Ac",
                Self::Acid => "A",
                Self::Alanyl => "Ala",
                Self::AcetylAlanyl => "Ala2Ac",
                Self::Acetimidoyl => "Am",
                Self::MethylAcetimidoyl => "AmMe",
                Self::Deoxy => "d",
                Self::Didehydro => "en",
                Self::DiMethylAcetimidoyl => "AmMe2",
                Self::Formyl => "Fo",
                Self::Glycolyl => "Gc",
                Self::AcetylGlutaminyl => "Gln2Ac",
                Self::MethylGlutamyl => "5Glu2Me",
                Self::Glycyl => "Gly",
                Self::Glyceryl => "Gr",
                Self::DiMethylGlyceryl => "Gr2,3Me2",
                Self::HydroxyButyryl => "Hb",
                Self::DiHydroxyButyryl => "3,4Hb",
                Self::Lactyl => "Lt",
                Self::Methyl => "Me",
                Self::Amino => "N",
                Self::NAcetyl => "NAc",
                Self::Phosphate => "P",
                Self::Pyruvyl => "Py",
                Self::CargoxyEthylidene => "Pyr",
                Self::Sulfate => "S",
                Self::Tauryl => "Tau",
            }
        )
    }
}

impl Chemical for GlycanSubstituent {
    fn formula(&self) -> MolecularFormula {
        let side = match self {
            Self::Acetyl => molecular_formula!(H 3 C 2 O 1),
            Self::Acid => molecular_formula!(H -1 O 2), // Together with the replacement below this is H-2 O+1
            Self::Alanyl => molecular_formula!(H 6 C 3 N 1 O 1),
            Self::AcetylAlanyl => molecular_formula!(H 8 C 5 N 1 O 2),
            Self::Acetimidoyl => molecular_formula!(H 5 C 2 N 1),
            Self::MethylAcetimidoyl => molecular_formula!(H 7 C 3 N 1),
            Self::Deoxy => molecular_formula!(H 1), // Together with the replacement below this is O-1
            Self::Didehydro => molecular_formula!(H -1 O 1), // Together with the replacement below this is H-2
            Self::DiMethylAcetimidoyl => molecular_formula!(H 9 C 4 N 1),
            Self::Formyl => molecular_formula!(H 1 C 1 O 1),
            Self::Glycolyl => molecular_formula!(H 3 C 2 O 2),
            Self::AcetylGlutaminyl => molecular_formula!(H 11 C 7 N 2 O 3),
            Self::MethylGlutamyl => molecular_formula!(H 10 C 6 N 1 O 3),
            Self::Glycyl | Self::NAcetyl => molecular_formula!(H 4 C 2 N 1 O 1),
            Self::Glyceryl => molecular_formula!(H 5 C 3 O 3),
            Self::DiMethylGlyceryl => molecular_formula!(H 9 C 5 O 3),
            Self::HydroxyButyryl => molecular_formula!(H 7 C 4 O 2),
            Self::DiHydroxyButyryl => molecular_formula!(H 7 C 4 O 3),
            Self::Lactyl => molecular_formula!(H 5 C 3 O 2),
            Self::Methyl => molecular_formula!(H 3 C 1),
            Self::Amino => molecular_formula!(H 2 N 1),
            Self::Phosphate => molecular_formula!(H 2 O 4 P 1),
            Self::Pyruvyl => molecular_formula!(H 3 C 3 O 2),
            Self::CargoxyEthylidene => molecular_formula!(H 3 C 3 O 3), // TODO: double substituent? Calculated to work with the additional side chain deletion
            Self::Sulfate => molecular_formula!(H 1 O 4 S 1),
            Self::Tauryl => molecular_formula!(H 6 C 2 N 1 O 3 S 1),
        };
        side - molecular_formula!(O 1 H 1) // substituent so replaces a standard oxygen side chain
    }
}

/// All monosaccharides ordered to be able to parse glycans by matching them from the top
#[allow(dead_code)]
pub const GLYCAN_PARSE_LIST: &[(&str, MonoSaccharide)] = &[
    (
        "phosphate",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: None,
            substituents: &[GlycanSubstituent::Phosphate],
            pro_forma_name: "phosphate",
        }),
    ),
    (
        "sulfate",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: None,
            substituents: &[GlycanSubstituent::Sulfate],
            pro_forma_name: "sulfate",
        }),
    ),
    (
        "Sug",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Sugar),
            substituents: &[],
            pro_forma_name: "Sug",
        }),
    ),
    (
        "Tri",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Triose),
            substituents: &[],
            pro_forma_name: "Tri",
        }),
    ),
    (
        "Tet",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Tetrose(None)),
            substituents: &[],
            pro_forma_name: "Tet",
        }),
    ),
    (
        "Pen",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Pentose(None)),
            substituents: &[],
            pro_forma_name: "Pen",
        }),
    ),
    (
        "a-Hex",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[GlycanSubstituent::Acid],
            pro_forma_name: "a-Hex",
        }),
    ),
    (
        "en,a-Hex",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[
                GlycanSubstituent::Acid,
                GlycanSubstituent::Didehydro,
                GlycanSubstituent::Deoxy,
            ],
            pro_forma_name: "en,a-Hex",
        }),
    ),
    (
        "d-Hex",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[],
            pro_forma_name: "d-Hex",
        }),
    ),
    (
        "HexNAc(S)",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[GlycanSubstituent::NAcetyl, GlycanSubstituent::Sulfate],
            pro_forma_name: "HexNAc(S)",
        }),
    ),
    (
        "HexNAc",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[GlycanSubstituent::NAcetyl],
            pro_forma_name: "HexNAc",
        }),
    ),
    (
        "HexNS",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[GlycanSubstituent::Amino, GlycanSubstituent::Sulfate],
            pro_forma_name: "HexNS",
        }),
    ),
    (
        "HexN",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[GlycanSubstituent::Amino],
            pro_forma_name: "HexN",
        }),
    ),
    (
        "HexS",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[GlycanSubstituent::Sulfate],
            pro_forma_name: "HexS",
        }),
    ),
    (
        "HexP",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[GlycanSubstituent::Phosphate],
            pro_forma_name: "HexP",
        }),
    ),
    (
        "Hex",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[],
            pro_forma_name: "Hex",
        }),
    ),
    (
        "Hep",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Heptose(None)),
            substituents: &[],
            pro_forma_name: "Hep",
        }),
    ),
    (
        "Oct",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Octose),
            substituents: &[],
            pro_forma_name: "Oct",
        }),
    ),
    (
        "Non",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Nonose),
            substituents: &[],
            pro_forma_name: "Non",
        }),
    ),
    (
        "Dec",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Decose),
            substituents: &[],
            pro_forma_name: "Dec",
        }),
    ),
    (
        "Neu5Ac",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Nonose),
            substituents: &[
                GlycanSubstituent::Amino,
                GlycanSubstituent::Acetyl,
                GlycanSubstituent::Acid,
            ],
            pro_forma_name: "Neu5Ac",
        }),
    ),
    (
        "Neu5Gc",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Nonose),
            substituents: &[
                GlycanSubstituent::Amino,
                GlycanSubstituent::Glycolyl,
                GlycanSubstituent::Acid,
            ],
            pro_forma_name: "Neu5Gc",
        }),
    ),
    (
        "Neu",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Nonose),
            substituents: &[
                GlycanSubstituent::Amino,
                GlycanSubstituent::Deoxy,
                GlycanSubstituent::Acid,
            ],
            pro_forma_name: "Neu",
        }),
    ),
    (
        "Fuc",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(Some(HexoseIsomer::Galactose))),
            substituents: &[GlycanSubstituent::Deoxy],
            pro_forma_name: "Fuc",
        }),
    ),
    // Single letter codes, by defining them like this they will be read but exported to the standard ProForma codes
    (
        "P",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(Some(HexoseIsomer::Mannose))),
            substituents: &[GlycanSubstituent::Phosphate],
            pro_forma_name: "Hexphosphate", // TODO: technically maybe not working when multiple are in there, think it through
        }),
    ),
    (
        "H",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[],
            pro_forma_name: "Hex",
        }),
    ),
    (
        "N",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(None)),
            substituents: &[GlycanSubstituent::NAcetyl],
            pro_forma_name: "HexNAc",
        }),
    ),
    (
        "F",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Hexose(Some(HexoseIsomer::Galactose))),
            substituents: &[GlycanSubstituent::Deoxy],
            pro_forma_name: "Fuc",
        }),
    ),
    (
        "S",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Nonose),
            substituents: &[
                GlycanSubstituent::Amino,
                GlycanSubstituent::Acetyl,
                GlycanSubstituent::Acid,
            ],
            pro_forma_name: "Neu5Ac",
        }),
    ),
    (
        "A",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Nonose),
            substituents: &[
                GlycanSubstituent::Amino,
                GlycanSubstituent::Acetyl,
                GlycanSubstituent::Acid,
            ],
            pro_forma_name: "Neu5Ac",
        }),
    ),
    (
        "G",
        MonoSaccharide::Predefined(PredefinedMonosaccharide {
            base_sugar: Some(BaseSugar::Nonose),
            substituents: &[
                GlycanSubstituent::Amino,
                GlycanSubstituent::Glycolyl,
                GlycanSubstituent::Acid,
            ],
            pro_forma_name: "Neu5Gc",
        }),
    ),
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
            sugar: MonoSaccharide::Allocated(AllocatedMonosaccharide {
                base_sugar: BaseSugar::Decose,
                substituents: Vec::new(),
            }),
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
                    sugar: sugar.clone(),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pro_forma_compliance() {
        let cases = &[
            ("Hep", molecular_formula!(H 12 C 7 O 6)),
            ("phosphate", molecular_formula!(H 1 O 3 P 1)),
            ("a-Hex", molecular_formula!(H 8 C 6 O 6)),
            ("Sug", molecular_formula!(H 2 C 2 O 1)),
            ("HexN", molecular_formula!(H 11 C 6 N 1 O 4)),
            ("Pen", molecular_formula!(H 8 C 5 O 4)),
            ("Tet", molecular_formula!(H 6 C 4 O 3)),
            ("HexP", molecular_formula!(H 11 C 6 O 8 P 1)),
            ("Neu5Ac", molecular_formula!(H 17 C 11 N 1 O 8)),
            ("Non", molecular_formula!(H 16 C 9 O 8)),
            ("HexNAc(S)", molecular_formula!(H 13 C 8 N 1 O 8 S 1)),
            ("Dec", molecular_formula!(H 18 C 10 O 9)),
            ("en,a-Hex", molecular_formula!(H 6 C 6 O 5)),
            ("Neu5Gc", molecular_formula!(H 17 C 11 N 1 O 9)),
            ("Neu", molecular_formula!(H 15 C 9 N 1 O 7)),
            ("HexNAc", molecular_formula!(H 13 C 8 N 1 O 5)),
            ("Fuc", molecular_formula!(H 10 C 6 O 4)),
            ("HexNS", molecular_formula!(H 11 C 6 N 1 O 7 S 1)),
            ("Tri", molecular_formula!(H 4 C 3 O 2)),
            ("Oct", molecular_formula!(H 14 C 8 O 7)),
            ("sulfate", molecular_formula!(O 3 S 1)),
            ("d-Hex", molecular_formula!(H 10 C 6 O 5)),
            ("Hex", molecular_formula!(H 10 C 6 O 5)),
            ("HexS", molecular_formula!(H 10 C 6 O 8 S 1)),
        ];
        for (name, formula) in cases {
            assert_eq!(
                GLYCAN_PARSE_LIST
                    .iter()
                    .find(|p| p.0 == *name)
                    .unwrap_or_else(|| panic!("Assumed {name} would be defined"))
                    .1
                    .formula(),
                *formula,
                "{name}",
            );
        }
    }
}
