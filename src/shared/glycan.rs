use std::{fmt::Display, hash::Hash, ops::Range};

use crate::{
    formula::{Chemical, MolecularFormula},
    helper_functions::*,
    Context, CustomError, Element, ELEMENT_PARSE_LIST,
};

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
    /// Hash implementation to try and give a predefined mono saccharide with the same composition as a allocated mono saccharide the same value
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

/// A predefined monosaccharide which can be used in const contexts
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PredefinedMonosaccharide {
    base_sugar: Option<BaseSugar>,
    substituents: &'static [GlycanSubstituent],
    pro_forma_name: &'static str,
}

/// An allocated monosaccharide which can be made during runtime
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AllocatedMonosaccharide {
    base_sugar: BaseSugar,
    substituents: Vec<GlycanSubstituent>,
    furanose: bool,
}

impl MonoSaccharide {
    pub fn new(sugar: BaseSugar, substituents: &[GlycanSubstituent]) -> Self {
        Self::Allocated(AllocatedMonosaccharide {
            base_sugar: sugar,
            substituents: substituents.to_owned(),
            furanose: false,
        })
    }

    /// Set this saccharide up as to be a furanose
    #[must_use]
    #[allow(dead_code)]
    pub fn furanose(self) -> Self {
        match self {
            Self::Predefined(pre) => Self::Predefined(pre),
            Self::Allocated(alo) => Self::Allocated(AllocatedMonosaccharide {
                furanose: true,
                ..alo
            }),
        }
    }

    /// Parse a short iupac name from this string starting at `start` and returning,
    /// if successful, a monosaccharide and the offset in the string where parsing ended.
    pub fn from_short_iupac(
        line: &str,
        start: usize,
        line_number: usize,
    ) -> Result<(Self, usize), CustomError> {
        let mut index = start;
        let bytes = line.as_bytes();
        let mut substituents = Vec::new();

        // ignore stuff
        if line[index..].starts_with("keto-") {
            index += 5;
        }
        // Isomeric state
        if line[index..].starts_with("D-")
            || line[index..].starts_with("L-")
            || line[index..].starts_with("?-")
        {
            index += 2; // Ignore the isomeric state otherwise
        }
        // Prefix mods
        let mut amount = 1;
        if bytes[index].is_ascii_digit() {
            match bytes[index + 1] {
                b',' if bytes[index + 3] == b':' => {
                    let start_index = index;
                    index += 7;
                    if bytes[index] == b'-' {
                        index += 1;
                    }
                    if !line[index..].starts_with("Anhydro") {
                        return Err(CustomError::error(
                            "Invalid iupac monosaccharide name",
                            "This internally linked glycan could not be parsed, expected Anhydro as modification",
                            Context::Line {
                                linenumber: line_number,
                                line: line.to_string(),
                                offset: start_index,
                                length: index-start_index+5,
                            },
                        ));
                    }
                    index += 7;
                    substituents.extend_from_slice(&[
                        GlycanSubstituent::Didehydro,
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituent::Deoxy,
                    ]);
                }
                b',' => {
                    let num = bytes[index + 1..]
                        .iter()
                        .take_while(|c| c.is_ascii_digit() || **c == b',' || **c == b'?')
                        .count();
                    index += num + 1;
                    amount = num / 2;
                    // X,X{mod} (or 3/4/5/etc mods)
                }
                _ => index += 1, // X{mod}
            }
            if bytes[index] == b'-' {
                index += 1;
            }
        }
        // Detect & ignore epi state
        if bytes[index] == b'e' {
            index += 1;
        }
        for substituent in PREFIX_SUBSTITUENTS {
            if line[index..].starts_with(substituent.0) {
                index += substituent.0.len();
                for _ in 0..amount {
                    substituents.push(substituent.1.clone());
                }
                if bytes[index] == b'-' {
                    index += 1;
                }
                break;
            }
        }
        // Another optional isomeric state
        if line[index..].starts_with("D-")
            || line[index..].starts_with("L-")
            || line[index..].starts_with("?-")
        {
            index += 2; // Ignore the isomeric state otherwise
        }
        // Base sugar
        let mut sugar = None;
        for sug in BASE_SUGARS {
            if line[index..].starts_with(sug.0) {
                index += sug.0.len();
                sugar = Some((sug.1.clone(), sug.2));
                break;
            }
        }
        let mut sugar = sugar
            .map(|(b, s)| {
                let mut alo = AllocatedMonosaccharide {
                    base_sugar: b,
                    substituents,
                    furanose: false,
                };
                alo.substituents.extend(s.iter().cloned());
                alo
            })
            .ok_or_else(|| {
                CustomError::error(
                    "Invalid iupac monosaccharide name",
                    "This name could not be recognised as a standard iupac glycan name",
                    Context::Line {
                        linenumber: line_number,
                        line: line.to_string(),
                        offset: index,
                        length: 3,
                    },
                )
            })?;
        // Furanose
        if index < bytes.len() && bytes[index] == b'f' {
            index += 1;
            sugar.furanose = true;
        }
        // Postfix mods
        while index < bytes.len() {
            if bytes[index] == b'-' {
                index += 1;
            }
            // Location
            let mut amount = 1;
            if bytes[index].is_ascii_digit() || bytes[index] == b'?' {
                match bytes[index + 1] {
                    b',' => {
                        let num = bytes[index + 1..]
                            .iter()
                            .take_while(|c| c.is_ascii_digit() || **c == b',' || **c == b'?')
                            .count();
                        index += num + 1;
                        amount = num / 2 + 1;
                        // X,X{mod} (or 3/4/5/etc mods)
                    }
                    b'-' => {
                        let mut found = false;
                        index += 7;
                        if line[index..].starts_with("(X)")
                            || line[index..].starts_with("(R)")
                            || line[index..].starts_with("(S)")
                        {
                            index += 3; // Ignore isomeric state
                        }
                        for substituent in DOUBLE_LINKED_POSTFIX_SUBSTITUENTS {
                            if line[index..].starts_with(substituent.0) {
                                index += substituent.0.len();
                                for _ in 0..amount {
                                    sugar.substituents.extend(substituent.1.iter().cloned());
                                }
                                found = true;
                                break;
                            }
                        }
                        if !found {
                            return Err(CustomError::error(
                            "Invalid iupac monosaccharide name",
                            "No detected double linked glycan substituent was found, while the pattern for location is for a double linked substituent",
                            Context::Line {
                                linenumber: line_number,
                                line: line.to_string(),
                                offset: index,
                                length: 2,
                            },
                        ));
                        }

                        continue; // Do not accepts an additional mod immediately
                    } // X-X,X-X (Py)
                    _ => index += 1, // X{mod}
                }
            }
            if line[index..].starts_with("(X)")
                || line[index..].starts_with("(R)")
                || line[index..].starts_with("(S)")
            {
                index += 3; // Ignore isomeric state
            }
            // Mod
            let mut found = false;
            for substituent in POSTFIX_SUBSTITUENTS {
                if line[index..].starts_with(substituent.0) {
                    index += substituent.0.len();
                    for _ in 0..amount {
                        sugar.substituents.push(substituent.1.clone());
                    }
                    found = true;
                    break;
                }
            }
            if !found {
                // Parse any element
                for element in ELEMENT_PARSE_LIST {
                    if line[index..].starts_with(element.0) {
                        index += element.0.len();
                        for _ in 0..amount {
                            sugar
                                .substituents
                                .push(GlycanSubstituent::Element(element.1));
                        }
                        found = true;
                        break;
                    }
                }
                if !found {
                    break;
                }
            }
            // Amount
            if amount != 1 {
                // Ignore the amount number, already determined before
                index += 1;
            }
        }
        Ok((Self::Allocated(sugar), index))
    }

    // fn symbol(&self) -> char {
    //     // ⬠◇♢▭◮⬘
    //     // ⬟◆♦▬
    //     match self {
    //         Self::Hexose => '○',      //●
    //         Self::HexNAc => '□',      // ■
    //         Self::Deoxyhexose => '△', // ▲
    //         Self::Pentose => '☆',     // ★
    //         Self::HexN => '⬔',
    //         _ => '⬡', // ⬢
    //     }
    // }
}

impl Chemical for MonoSaccharide {
    fn formula(&self) -> MolecularFormula {
        match self {
            Self::Allocated(AllocatedMonosaccharide {
                base_sugar,
                substituents,
                ..
            }) => base_sugar.formula() + substituents.formula(),
            Self::Predefined(PredefinedMonosaccharide {
                base_sugar,
                substituents,
                ..
            }) => {
                base_sugar
                    .as_ref()
                    .map(Chemical::formula)
                    .unwrap_or_default()
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
                    furanose,
                }) => format!(
                    "{}{}{}",
                    base_sugar,
                    if *furanose { "f" } else { "" },
                    substituents
                        .iter()
                        .map(ToString::to_string)
                        .collect::<String>()
                ),
                Self::Predefined(PredefinedMonosaccharide { pro_forma_name, .. }) =>
                    (*pro_forma_name).to_string(),
            }
        )
    }
}

/// The base sugar of a monosaccharide, optionally with the isomeric state saved as well.
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
                Self::Sugar => "Sug".to_string(),
                Self::Triose => "Tri".to_string(),
                Self::Tetrose(_) => "Tet".to_string(),
                Self::Pentose(_) => "Pen".to_string(),
                Self::Hexose(_) => "Hex".to_string(),
                Self::Heptose(_) => "Hep".to_string(),
                Self::Octose => "Oct".to_string(),
                Self::Nonose => "Non".to_string(),
                Self::Decose => "Dec".to_string(),
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

/// Any 4 carbon glycan
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum TetroseIsomer {
    /// Ery
    Erythrose,
    /// Tho
    Threose,
}

/// Any 5 carbon glycan
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum PentoseIsomer {
    /// Rib
    Ribose,
    /// Ara
    Arabinose,
    /// Xyl
    Xylose,
    /// Lyx
    Lyxose,
    /// Xul
    Xylulose,
}

/// Any 6 carbon glycan
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum HexoseIsomer {
    /// glc
    Glucose,
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
    /// Psi
    Psicose,
    /// Fru
    Fructose,
    /// Sor
    Sorbose,
    /// Tag
    Tagatose,
}

/// Any 7 carbon glycan
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum HeptoseIsomer {
    /// gro-manHep
    GlyceroMannoHeptopyranose, // TODO: Does this indicate some mods?
    /// Sed
    Sedoheptulose,
}

const BASE_SUGARS: &[(&str, BaseSugar, &[GlycanSubstituent])] = &[
    ("Sug", BaseSugar::Sugar, &[]),
    ("Tri", BaseSugar::Triose, &[]),
    ("Tet", BaseSugar::Tetrose(None), &[]),
    (
        "Ery",
        BaseSugar::Tetrose(Some(TetroseIsomer::Erythrose)),
        &[],
    ),
    ("Tho", BaseSugar::Tetrose(Some(TetroseIsomer::Threose)), &[]),
    ("Pen", BaseSugar::Pentose(None), &[]),
    ("Rib", BaseSugar::Pentose(Some(PentoseIsomer::Ribose)), &[]),
    (
        "Ara",
        BaseSugar::Pentose(Some(PentoseIsomer::Arabinose)),
        &[],
    ),
    ("Xyl", BaseSugar::Pentose(Some(PentoseIsomer::Xylose)), &[]),
    ("Lyx", BaseSugar::Pentose(Some(PentoseIsomer::Lyxose)), &[]),
    ("Hex", BaseSugar::Hexose(None), &[]),
    ("Glc", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[]),
    ("Gal", BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[]),
    ("Man", BaseSugar::Hexose(Some(HexoseIsomer::Mannose)), &[]),
    ("All", BaseSugar::Hexose(Some(HexoseIsomer::Allose)), &[]),
    ("Alt", BaseSugar::Hexose(Some(HexoseIsomer::Altrose)), &[]),
    ("Gul", BaseSugar::Hexose(Some(HexoseIsomer::Gulose)), &[]),
    ("Ido", BaseSugar::Hexose(Some(HexoseIsomer::Idose)), &[]),
    ("Tal", BaseSugar::Hexose(Some(HexoseIsomer::Talose)), &[]),
    ("Hep", BaseSugar::Heptose(None), &[]),
    (
        "gro-manHep",
        BaseSugar::Heptose(Some(HeptoseIsomer::GlyceroMannoHeptopyranose)),
        &[],
    ),
    (
        "Neu",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "Sia",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "Kdn",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "Kdo",
        BaseSugar::Octose,
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Acid],
    ),
    (
        "Fuc",
        BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "Rha",
        BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "Qui",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "Oli",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Tyv",
        BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Asc",
        BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Abe",
        BaseSugar::Hexose(Some(HexoseIsomer::Gulose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Par",
        BaseSugar::Hexose(Some(HexoseIsomer::Altrose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Dig",
        BaseSugar::Hexose(Some(HexoseIsomer::Allose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Col",
        BaseSugar::Hexose(Some(HexoseIsomer::Talose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    ("Psi", BaseSugar::Hexose(Some(HexoseIsomer::Psicose)), &[]),
    ("Fru", BaseSugar::Hexose(Some(HexoseIsomer::Fructose)), &[]),
    ("Sor", BaseSugar::Hexose(Some(HexoseIsomer::Sorbose)), &[]),
    ("Tag", BaseSugar::Hexose(Some(HexoseIsomer::Tagatose)), &[]),
    (
        "Xul",
        BaseSugar::Pentose(Some(PentoseIsomer::Xylulose)),
        &[],
    ),
    (
        "Sed",
        BaseSugar::Heptose(Some(HeptoseIsomer::Sedoheptulose)),
        &[],
    ),
    (
        "MurNAc",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::NAcetyl, GlycanSubstituent::OCarboxyEthyl],
    ),
    (
        "MurNGc",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[
            GlycanSubstituent::NGlycolyl,
            GlycanSubstituent::OCarboxyEthyl,
        ],
    ),
    (
        "Mur",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Amino, GlycanSubstituent::OCarboxyEthyl],
    ),
    (
        "Api",
        BaseSugar::Tetrose(Some(TetroseIsomer::Erythrose)),
        &[GlycanSubstituent::HydroxyMethyl],
    ),
    (
        "Dha",
        BaseSugar::Heptose(None),
        &[
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "Bac",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
        ],
    ),
    (
        "Pse",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "Leg",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Acid,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Deoxy,
        ],
    ),
    (
        "Aci",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Acid,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Deoxy,
        ],
    ),
];

/// Any substituent on a monosaccharide.
/// Source: <https://www.ncbi.nlm.nih.gov/glycans/snfg.html> table 3.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum GlycanSubstituent {
    ///Am N-acetimidoyl
    Acetimidoyl,
    ///Ac acetyl
    Acetyl,
    ///Ala2Ac N-acetyl-D-alanyl
    AcetylAlanyl,
    ///Gln2Ac N-acetyl-glutaminyl
    AcetylGlutaminyl,
    ///A acid
    Acid,
    ///Ala D-alanyl
    Alanyl,
    ///ol alcohol
    Alcohol,
    ///N amino
    Amino,
    ///aric ??
    Aric,
    ///Pyr 1-carboxyethylidene
    CargoxyEthylidene,
    ///d Deoxy
    Deoxy,
    ///3,4Hb 3,4-dihydroxybutyryl
    DiHydroxyButyryl,
    ///AmMe2 N-(N,N-dimethyl-acetimidoyl)
    DiMethylAcetimidoyl,
    ///Gr2,3Me2 2,3-di-O-methyl-glyceryl
    DiMethylGlyceryl,
    ///en didehydro an addition of a double bond
    Didehydro,
    ///An element that replaces a side chain
    Element(Element),
    ///Etn Ethanolamine
    Ethanolamine,
    ///EtOH O linked ethanol
    EtOH,
    ///Fo formyl
    Formyl,
    ///Gr glyceryl
    Glyceryl,
    ///Gc glycolyl
    Glycolyl,
    ///Gly glycyl
    Glycyl,
    ///4Hb 4-hydroxybutyryl, 3RHb (R)-3-hydroxybutyryl, 3SHb (S)-3-hydroxybutyryl
    HydroxyButyryl,
    ///HydroxyMethyl
    HydroxyMethyl,
    ///Lac
    Lac,
    ///Lt lactyl
    Lactyl,
    ///Me methyl
    Methyl,
    ///AmMe N-(N-methyl-acetimidoyl)
    MethylAcetimidoyl,
    ///5Glu2Me N-methyl-5-glutamyl
    MethylGlutamyl,
    ///NAc N-acetyl
    NAcetyl,
    ///N2DiMe N linked double methyl
    NDiMe,
    ///NFo N linked formyl
    NFo,
    ///NGc N linked glycolyl
    NGlycolyl,
    ///carboxyethyl used in Mur
    OCarboxyEthyl,
    ///PCho phosphate linked choline
    PCholine,
    ///P phosphate
    Phosphate,
    ///Py pyruvyl
    Pyruvyl,
    ///Suc ??
    Suc,
    ///S sulfate
    Sulfate,
    ///Tau tauryl
    Tauryl,
    ///ulo ??
    Ulo,
    ///ulof ??
    Ulof,
    ///water H2O
    Water,
}

const POSTFIX_SUBSTITUENTS: &[(&str, GlycanSubstituent)] = &[
    ("Ac", GlycanSubstituent::Acetyl),
    ("Ala2Ac", GlycanSubstituent::AcetylAlanyl),
    ("Ala", GlycanSubstituent::Alanyl),
    ("AmMe2", GlycanSubstituent::DiMethylAcetimidoyl),
    ("AmMe", GlycanSubstituent::MethylAcetimidoyl),
    ("Am", GlycanSubstituent::Acetimidoyl),
    ("en", GlycanSubstituent::Didehydro),
    ("Fo", GlycanSubstituent::Formyl),
    ("Gc", GlycanSubstituent::Glycolyl),
    ("Gln2Ac", GlycanSubstituent::AcetylGlutaminyl),
    ("5Glu2Me", GlycanSubstituent::MethylGlutamyl), // TODO: Does this number has to be there? (Not found in the glycosmos data)
    ("Gly", GlycanSubstituent::Glycyl),
    ("Gr", GlycanSubstituent::Glyceryl),
    ("Gr2,3Me2", GlycanSubstituent::DiMethylGlyceryl),
    ("4Hb", GlycanSubstituent::HydroxyButyryl),
    ("3RHb", GlycanSubstituent::HydroxyButyryl),
    ("3SHb", GlycanSubstituent::HydroxyButyryl),
    ("3,4Hb", GlycanSubstituent::DiHydroxyButyryl),
    ("Lt", GlycanSubstituent::Lactyl),
    ("Lac", GlycanSubstituent::Lac),
    ("Me", GlycanSubstituent::Methyl),
    ("NAc", GlycanSubstituent::NAcetyl),
    ("Pyr", GlycanSubstituent::CargoxyEthylidene),
    ("Tau", GlycanSubstituent::Tauryl),
    ("onic", GlycanSubstituent::Acid),
    ("uronic", GlycanSubstituent::Acid),
    ("aric", GlycanSubstituent::Aric),
    ("ol", GlycanSubstituent::Alcohol),
    ("Etn", GlycanSubstituent::Ethanolamine),
    ("EtOH", GlycanSubstituent::Ethanolamine),
    ("ulof", GlycanSubstituent::Ulof),
    ("ulo", GlycanSubstituent::Ulo),
    ("N2DiMe", GlycanSubstituent::NDiMe),
    ("NDiMe", GlycanSubstituent::NDiMe),
    ("PCho", GlycanSubstituent::PCholine),
    ("CE", GlycanSubstituent::Glycyl), // Same molecular formula
    ("Suc", GlycanSubstituent::Suc),
    ("NFo", GlycanSubstituent::NFo),
    ("A", GlycanSubstituent::Acid),
    ("P", GlycanSubstituent::Phosphate),
    ("p", GlycanSubstituent::Phosphate),
    ("S", GlycanSubstituent::Sulfate),
    ("N", GlycanSubstituent::Amino),
];

const DOUBLE_LINKED_POSTFIX_SUBSTITUENTS: &[(&str, &[GlycanSubstituent])] = &[
    ("Py", &[GlycanSubstituent::Pyruvyl]),
    ("N", &[GlycanSubstituent::Water]),
    (
        "P",
        &[GlycanSubstituent::Water, GlycanSubstituent::Phosphate],
    ),
];

const PREFIX_SUBSTITUENTS: &[(&str, GlycanSubstituent)] = &[
    ("deoxy", GlycanSubstituent::Deoxy),
    ("Anhydro", GlycanSubstituent::Deoxy),
    ("d", GlycanSubstituent::Deoxy),
];

impl Display for GlycanSubstituent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Acetimidoyl => "Am".to_string(),
                Self::Acetyl => "Ac".to_string(),
                Self::AcetylAlanyl => "Ala2Ac".to_string(),
                Self::AcetylGlutaminyl => "Gln2Ac".to_string(),
                Self::Acid => "A".to_string(),
                Self::Alanyl => "Ala".to_string(),
                Self::Alcohol => "ol".to_string(),
                Self::Amino => "N".to_string(),
                Self::Aric => "aric".to_string(),
                Self::CargoxyEthylidene => "Pyr".to_string(),
                Self::Deoxy => "d".to_string(),
                Self::Didehydro => "en".to_string(),
                Self::DiHydroxyButyryl => "3,4Hb".to_string(),
                Self::DiMethylAcetimidoyl => "AmMe2".to_string(),
                Self::DiMethylGlyceryl => "Gr2,3Me2".to_string(),
                Self::Ethanolamine => "Etn".to_string(),
                Self::Element(el) => el.to_string(),
                Self::EtOH => "EtOH".to_string(),
                Self::Formyl => "Fo".to_string(),
                Self::Glyceryl => "Gr".to_string(),
                Self::Glycolyl => "Gc".to_string(),
                Self::Glycyl => "Gly".to_string(),
                Self::HydroxyButyryl => "Hb".to_string(),
                Self::HydroxyMethyl => "HMe".to_string(),
                Self::Lac => "Lac".to_string(),
                Self::Lactyl => "Lt".to_string(),
                Self::Methyl => "Me".to_string(),
                Self::MethylAcetimidoyl => "AmMe".to_string(),
                Self::MethylGlutamyl => "5Glu2Me".to_string(),
                Self::NAcetyl => "NAc".to_string(),
                Self::NDiMe => "NDiMe".to_string(),
                Self::NFo => "NFo".to_string(),
                Self::NGlycolyl => "NGc".to_string(),
                Self::OCarboxyEthyl => "carboxyethyl".to_string(),
                Self::PCholine => "PCho".to_string(),
                Self::Phosphate => "P".to_string(),
                Self::Pyruvyl => "Py".to_string(),
                Self::Suc => "Suc".to_string(),
                Self::Sulfate => "S".to_string(),
                Self::Tauryl => "Tau".to_string(),
                Self::Ulo => "ulo".to_string(),
                Self::Ulof => "ulof".to_string(),
                Self::Water => "waterloss".to_string(),
            }
        )
    }
}

impl Chemical for GlycanSubstituent {
    fn formula(&self) -> MolecularFormula {
        let side = match self {
            Self::Acetimidoyl => molecular_formula!(H 5 C 2 N 1),
            Self::Acetyl => molecular_formula!(H 3 C 2 O 1),
            Self::AcetylAlanyl => molecular_formula!(H 8 C 5 N 1 O 2),
            Self::AcetylGlutaminyl => molecular_formula!(H 11 C 7 N 2 O 3),
            Self::Acid => molecular_formula!(H -1 O 2), // Together with the replacement below this is H-2 O+1
            Self::Alanyl => molecular_formula!(H 6 C 3 N 1 O 1),
            Self::Alcohol => molecular_formula!(H 3 O 1), // Together with the replacement below this is H+2
            Self::Amino => molecular_formula!(H 2 N 1),
            Self::Aric => molecular_formula!(H 3 O 3), // Together with replacement below this is H2O2
            Self::CargoxyEthylidene => molecular_formula!(H 3 C 3 O 3), // double substituent, calculated to work with the additional side chain deletion
            Self::Deoxy => molecular_formula!(H 1), // Together with the replacement below this is O-1
            Self::Didehydro => molecular_formula!(H -1 O 1), // Together with the replacement below this is H-2
            Self::DiHydroxyButyryl => molecular_formula!(H 7 C 4 O 3),
            Self::DiMethylAcetimidoyl => molecular_formula!(H 9 C 4 N 1),
            Self::DiMethylGlyceryl => molecular_formula!(H 9 C 5 O 3),
            Self::Ethanolamine => molecular_formula!(H 6 C 2 N 1 O 1),
            Self::EtOH => molecular_formula!(H 5 C 2 O 2),
            Self::Element(el) => MolecularFormula::new(&[(*el, 0, 1)]),
            Self::Formyl => molecular_formula!(H 1 C 1 O 1),
            Self::Glyceryl => molecular_formula!(H 5 C 3 O 3),
            Self::Glycolyl => molecular_formula!(H 3 C 2 O 2),
            Self::Glycyl | Self::NAcetyl => molecular_formula!(H 4 C 2 N 1 O 1),
            Self::HydroxyButyryl => molecular_formula!(H 7 C 4 O 2),
            Self::HydroxyMethyl => molecular_formula!(H 3 C 1 O 2),
            Self::Lac => molecular_formula!(H 5 C 3 O 3),
            Self::Lactyl => molecular_formula!(H 5 C 3 O 2),
            Self::Methyl => molecular_formula!(H 3 C 1),
            Self::MethylAcetimidoyl => molecular_formula!(H 7 C 3 N 1),
            Self::MethylGlutamyl => molecular_formula!(H 10 C 6 N 1 O 3),
            Self::NDiMe => molecular_formula!(H 6 C 2 N 1),
            Self::NFo => molecular_formula!(H 2 C 1 N 1 O 1),
            Self::NGlycolyl => molecular_formula!(H 4 C 2 N 1 O 2),
            Self::OCarboxyEthyl => molecular_formula!(H 6 C 3 O 3), // Replaces H, together with replacement below this is H5C3O2
            Self::PCholine => molecular_formula!(H 14 C 5 N 1 O 4 P 1),
            Self::Phosphate => molecular_formula!(H 2 O 4 P 1),
            Self::Pyruvyl => molecular_formula!(H 3 C 3 O 2),
            Self::Suc => molecular_formula!(H 6 C 4 N 1 O 3),
            Self::Sulfate => molecular_formula!(H 1 O 4 S 1),
            Self::Tauryl => molecular_formula!(H 6 C 2 N 1 O 3 S 1),
            Self::Ulo => molecular_formula!(H 3 C 1 O 2), // Replaces H, together with replacement below this is H2C1O1
            Self::Ulof => molecular_formula!(H 4 C 1 O 2), // Replaces H, together with replacement below this is H3C1O1
            Self::Water => molecular_formula!(H - 1),
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
            pro_forma_name: "Hexphosphate", // TODO: technically maybe not working when multiple are in there, think it through, should be two different elements, both getting counts after them
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
            sugar: MonoSaccharide::new(BaseSugar::Decose, &[]),
            branches: Vec::new(),
        }; // Starting sugar, will be removed
        let mut last_branch: &mut Self = &mut branch;
        let bytes = line.as_bytes();
        while offset < range.end {
            while bytes[offset] == b'[' {
                let end = end_of_enclosure(bytes, offset + 1, b'[', b']').ok_or_else(|| {
                    CustomError::error(
                        "Invalid iupac short glycan",
                        "No closing brace for branch",
                        Context::line(line_number, line, offset, range.end - offset),
                    )
                })?;
                let offshoot = Self::from_short_iupac(line, offset + 1..end, line_number)?;
                last_branch.branches.push(offshoot);
                offset = end + 1;
            }
            let (sugar, new_offset) = MonoSaccharide::from_short_iupac(line, offset, line_number)?;
            offset = new_offset;

            last_branch.branches.push(Self {
                sugar: sugar.clone(),
                branches: Vec::new(),
            });
            last_branch = last_branch.branches.last_mut().unwrap();
            if bytes[offset] == b'(' {
                if let Some(end) = next_char(bytes, offset + 1, b')') {
                    offset = end + 1; // just ignore all linking stuff I do not care
                } else {
                    // This only happens for incomplete branches where the last parts of the branch are unknown.
                    assert!(range.end - offset < 10); // make sure it is the last part
                    offset = range.end; // assume it is the last not closed brace
                }
            }
        }
        if let Some(glycan) = branch.branches.pop() {
            Ok(glycan) // Remove the outer starting sugar
        } else {
            Err(CustomError::error(
                "Invalid iupac short glycan",
                "No glycan found",
                Context::line(line_number, line.to_string(), range.start, range.len()),
            ))
        }
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

    #[test]
    fn iupac_short_names() {
        assert_eq!(
            MonoSaccharide::from_short_iupac("Gal2,3Ac24-1,6-1Py", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    &[
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Pyruvyl,
                    ]
                ),
                18
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("GlcNAc", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                    &[GlycanSubstituent::NAcetyl]
                ),
                6
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("Gal6S", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    &[GlycanSubstituent::Sulfate]
                ),
                5
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("GlcN2Gc", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                    &[GlycanSubstituent::Amino, GlycanSubstituent::Glycolyl,]
                ),
                7
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("GalNAc3S", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    &[GlycanSubstituent::NAcetyl, GlycanSubstituent::Sulfate]
                ),
                8
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("GlcN2,6S2", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                    &[
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Sulfate,
                        GlycanSubstituent::Sulfate
                    ]
                ),
                9
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("Tagf1,6P2", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Tagatose)),
                    &[GlycanSubstituent::Phosphate, GlycanSubstituent::Phosphate]
                )
                .furanose(),
                9
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("Gal2,3Ac24-1,6-1Py", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    &[
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Pyruvyl,
                    ]
                ),
                18
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("D-Araf", 0, 0),
            Ok((
                MonoSaccharide::new(BaseSugar::Pentose(Some(PentoseIsomer::Arabinose)), &[])
                    .furanose(),
                6
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("Xyl-onic", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Pentose(Some(PentoseIsomer::Xylose)),
                    &[GlycanSubstituent::Acid]
                ),
                8
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("Glc2,3,4,6Ac4", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                    &[
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acetyl
                    ]
                ),
                13
            ))
        );
    }
}
