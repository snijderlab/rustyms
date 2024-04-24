const BASE_SUGARS: &[(&str, BaseSugar, &[GlycanSubstituent])] = &[
    ("sug", BaseSugar::Sugar, &[]),
    ("tri", BaseSugar::Triose, &[]),
    ("tet", BaseSugar::Tetrose(None), &[]),
    (
        "ery",
        BaseSugar::Tetrose(Some(TetroseIsomer::Erythrose)),
        &[],
    ),
    ("tho", BaseSugar::Tetrose(Some(TetroseIsomer::Threose)), &[]),
    ("pen", BaseSugar::Pentose(None), &[]),
    ("rib", BaseSugar::Pentose(Some(PentoseIsomer::Ribose)), &[]),
    (
        "ara",
        BaseSugar::Pentose(Some(PentoseIsomer::Arabinose)),
        &[],
    ),
    ("xyl", BaseSugar::Pentose(Some(PentoseIsomer::Xylose)), &[]),
    ("lyx", BaseSugar::Pentose(Some(PentoseIsomer::Lyxose)), &[]),
    ("hex", BaseSugar::Hexose(None), &[]),
    ("glc", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[]),
    ("gal", BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[]),
    ("man", BaseSugar::Hexose(Some(HexoseIsomer::Mannose)), &[]),
    ("all", BaseSugar::Hexose(Some(HexoseIsomer::Allose)), &[]),
    ("alt", BaseSugar::Hexose(Some(HexoseIsomer::Altrose)), &[]),
    ("gul", BaseSugar::Hexose(Some(HexoseIsomer::Gulose)), &[]),
    ("ido", BaseSugar::Hexose(Some(HexoseIsomer::Idose)), &[]),
    ("tal", BaseSugar::Hexose(Some(HexoseIsomer::Talose)), &[]),
    ("hep", BaseSugar::Heptose(None), &[]),
    (
        "gro-manhep",
        BaseSugar::Heptose(Some(HeptoseIsomer::GlyceroMannoHeptopyranose)),
        &[],
    ),
    (
        "neu",
        BaseSugar::Nonose,
        &[GlycanSubstituent::Amino, GlycanSubstituent::Acid],
    ),
    (
        "sia",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "kdn",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "kdo",
        BaseSugar::Octose,
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Acid],
    ),
    (
        "fuc",
        BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "rha",
        BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "qui",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "oli",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "tyv",
        BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "asc",
        BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "abe",
        BaseSugar::Hexose(Some(HexoseIsomer::Gulose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "par",
        BaseSugar::Hexose(Some(HexoseIsomer::Altrose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "dig",
        BaseSugar::Hexose(Some(HexoseIsomer::Allose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "col",
        BaseSugar::Hexose(Some(HexoseIsomer::Talose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    ("psi", BaseSugar::Hexose(Some(HexoseIsomer::Psicose)), &[]),
    ("fru", BaseSugar::Hexose(Some(HexoseIsomer::Fructose)), &[]),
    ("sor", BaseSugar::Hexose(Some(HexoseIsomer::Sorbose)), &[]),
    ("tag", BaseSugar::Hexose(Some(HexoseIsomer::Tagatose)), &[]),
    (
        "xul",
        BaseSugar::Pentose(Some(PentoseIsomer::Xylulose)),
        &[],
    ),
    (
        "sed",
        BaseSugar::Heptose(Some(HeptoseIsomer::Sedoheptulose)),
        &[],
    ),
    (
        "murnac",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::NAcetyl, GlycanSubstituent::OCarboxyEthyl],
    ),
    (
        "murngc",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[
            GlycanSubstituent::NGlycolyl,
            GlycanSubstituent::OCarboxyEthyl,
        ],
    ),
    (
        "mur",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Amino, GlycanSubstituent::OCarboxyEthyl],
    ),
    (
        "api",
        BaseSugar::Tetrose(Some(TetroseIsomer::Erythrose)),
        &[GlycanSubstituent::HydroxyMethyl],
    ),
    (
        "dha",
        BaseSugar::Heptose(None),
        &[
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "bac",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
        ],
    ),
    (
        "pse",
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
        "leg",
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
        "aci",
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

// TODO: Points from mobiusklein in rusteomics/mzcore/pull/2
// * Remove the numbers from the names where already covered by the parsing code
// * Add an additional level which defines the leaving group, to make the chemical formula difference easier
const POSTFIX_SUBSTITUENTS: &[(&str, GlycanSubstituent)] = &[
    ("ac", GlycanSubstituent::Acetyl),
    ("ala2ac", GlycanSubstituent::AcetylAlanyl),
    ("ala", GlycanSubstituent::Alanyl),
    ("amme2", GlycanSubstituent::DiMethylAcetimidoyl),
    ("amme", GlycanSubstituent::MethylAcetimidoyl),
    ("am", GlycanSubstituent::Acetimidoyl),
    ("en", GlycanSubstituent::Didehydro),
    ("fo", GlycanSubstituent::Formyl),
    ("gc", GlycanSubstituent::Glycolyl),
    ("gln2ac", GlycanSubstituent::AcetylGlutaminyl),
    ("5glu2me", GlycanSubstituent::MethylGlutamyl),
    ("gly", GlycanSubstituent::Glycyl),
    ("gr", GlycanSubstituent::Glyceryl),
    ("gr2,3Me2", GlycanSubstituent::DiMethylGlyceryl),
    ("4hb", GlycanSubstituent::HydroxyButyryl),
    ("3rhb", GlycanSubstituent::HydroxyButyryl),
    ("3shb", GlycanSubstituent::HydroxyButyryl),
    ("3,4Hb", GlycanSubstituent::DiHydroxyButyryl),
    ("lt", GlycanSubstituent::Lactyl),
    ("lac", GlycanSubstituent::Lac),
    ("me", GlycanSubstituent::Methyl),
    ("cme", GlycanSubstituent::Methyl), // unsure about the difference with Me
    ("nac", GlycanSubstituent::NAcetyl),
    ("pyr", GlycanSubstituent::CargoxyEthylidene),
    ("tau", GlycanSubstituent::Tauryl),
    ("onic", GlycanSubstituent::Acid),
    ("uronic", GlycanSubstituent::Acid),
    ("aric", GlycanSubstituent::Aric),
    ("ol", GlycanSubstituent::Alcohol),
    ("etn", GlycanSubstituent::Ethanolamine),
    ("etoh", GlycanSubstituent::Ethanolamine),
    ("ulof", GlycanSubstituent::Ulof),
    ("ulo", GlycanSubstituent::Ulo),
    ("n2dime", GlycanSubstituent::NDiMe),
    ("ndime", GlycanSubstituent::NDiMe),
    ("pcho", GlycanSubstituent::PCholine),
    ("ce", GlycanSubstituent::Glycyl), // Same molecular formula
    ("suc", GlycanSubstituent::Suc),
    ("nfo", GlycanSubstituent::NFo),
    ("dime", GlycanSubstituent::DiMethyl),
    ("a", GlycanSubstituent::Acid),
    ("p", GlycanSubstituent::Phosphate),
    ("s", GlycanSubstituent::Sulfate),
    ("n", GlycanSubstituent::Amino),
];

const DOUBLE_LINKED_POSTFIX_SUBSTITUENTS: &[(&str, &[GlycanSubstituent])] = &[
    ("py", &[GlycanSubstituent::Pyruvyl]),
    ("n", &[GlycanSubstituent::Water]),
    (
        "p",
        &[GlycanSubstituent::Water, GlycanSubstituent::Phosphate],
    ),
];

const PREFIX_SUBSTITUENTS: &[(&str, GlycanSubstituent)] = &[
    ("deoxy", GlycanSubstituent::Deoxy),
    ("anhydro", GlycanSubstituent::Deoxy),
    ("d", GlycanSubstituent::Deoxy),
];

/// All monosaccharides ordered to be able to parse glycans by matching them from the top
#[allow(dead_code)]
pub fn glycan_parse_list() -> &'static Vec<(String, MonoSaccharide)> {
    GLYCAN_PARSE_CELL.get_or_init(|| {
        vec![
            (
                "phosphate".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("phosphate".to_string()),
                    furanose: false,
                },
            ),
            (
                "sulfate".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Sulfate],
                    proforma_name: Some("sulfate".to_string()),
                    furanose: false,
                },
            ),
            (
                "sug".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Sugar,
                    substituents: vec![],
                    proforma_name: Some("Sug".to_string()),
                    furanose: false,
                },
            ),
            (
                "tri".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Triose,
                    substituents: vec![],
                    proforma_name: Some("Tri".to_string()),
                    furanose: false,
                },
            ),
            (
                "tet".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Tetrose(None),
                    substituents: vec![],
                    proforma_name: Some("Tet".to_string()),
                    furanose: false,
                },
            ),
            (
                "pen".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(None),
                    substituents: vec![],
                    proforma_name: Some("Pen".to_string()),
                    furanose: false,
                },
            ),
            (
                "a-hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Acid],
                    proforma_name: Some("a-Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "en,a-hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Didehydro,
                        GlycanSubstituent::Deoxy,
                    ],
                    proforma_name: Some("en,a-Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "d-hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("d-Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "hexnac(s)".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl, GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexNAc(S)".to_string()),
                    furanose: false,
                },
            ),
            (
                "hexnac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                },
            ),
            (
                "hexns".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Amino, GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexNS".to_string()),
                    furanose: false,
                },
            ),
            (
                "hexn".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Amino],
                    proforma_name: Some("HexN".to_string()),
                    furanose: false,
                },
            ),
            (
                "hexs".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexS".to_string()),
                    furanose: false,
                },
            ),
            (
                "hexp".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("HexP".to_string()),
                    furanose: false,
                },
            ),
            (
                "hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "hep".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Heptose(None),
                    substituents: vec![],
                    proforma_name: Some("Hep".to_string()),
                    furanose: false,
                },
            ),
            (
                "oct".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Octose,
                    substituents: vec![],
                    proforma_name: Some("Oct".to_string()),
                    furanose: false,
                },
            ),
            (
                "non".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![],
                    proforma_name: Some("Non".to_string()),
                    furanose: false,
                },
            ),
            (
                "dec".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Decose,
                    substituents: vec![],
                    proforma_name: Some("Dec".to_string()),
                    furanose: false,
                },
            ),
            (
                "neu5ac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu5Ac".to_string()),
                    furanose: false,
                },
            ),
            (
                "neu5gc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu5Gc".to_string()),
                    furanose: false,
                },
            ),
            (
                "neu".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu".to_string()),
                    furanose: false,
                },
            ),
            (
                "fuc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::Deoxy],
                    proforma_name: Some("Fuc".to_string()),
                    furanose: false,
                },
            ),
            // Single letter codes, by defining them like this they will be read but exported to the standard ProForma codes
            (
                "p".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("Hexphosphate".to_string()), // TODO: technically maybe not working when multiple are in there, think it through, should be two different elements,  both getting counts after them
                    furanose: false,
                },
            ),
            (
                "h".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "n".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                },
            ),
            (
                "f".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::Deoxy],
                    proforma_name: Some("Fuc".to_string()),
                    furanose: false,
                },
            ),
            (
                "s".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu5Ac".to_string()),
                    furanose: false,
                },
            ),
            (
                "a".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu5Ac".to_string()),
                    furanose: false,
                },
            ),
            (
                "g".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu5Gc".to_string()),
                    furanose: false,
                },
            ),
        ]
    })
}
#[allow(dead_code)]
static GLYCAN_PARSE_CELL: OnceLock<Vec<(String, MonoSaccharide)>> = OnceLock::new();
