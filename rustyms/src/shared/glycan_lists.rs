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
        &[GlycanSubstituent::Amino, GlycanSubstituent::Acid],
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

// TODO: Points from mobiusklein in rusteomics/mzcore/pull/2
// * Remove the numbers from the names where already covered by the parsing code
// * Add an additional level which defines the leaving group, to make the chemical formula difference easier
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
    ("5Glu2Me", GlycanSubstituent::MethylGlutamyl),
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
    ("CMe", GlycanSubstituent::Methyl), // unsure about the difference with Me
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
    ("DiMe", GlycanSubstituent::DiMethyl),
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
                "Sug".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Sugar,
                    substituents: vec![],
                    proforma_name: Some("Sug".to_string()),
                    furanose: false,
                },
            ),
            (
                "Tri".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Triose,
                    substituents: vec![],
                    proforma_name: Some("Tri".to_string()),
                    furanose: false,
                },
            ),
            (
                "Tet".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Tetrose(None),
                    substituents: vec![],
                    proforma_name: Some("Tet".to_string()),
                    furanose: false,
                },
            ),
            (
                "Pen".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(None),
                    substituents: vec![],
                    proforma_name: Some("Pen".to_string()),
                    furanose: false,
                },
            ),
            (
                "a-Hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Acid],
                    proforma_name: Some("a-Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "en,a-Hex".to_string(),
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
                "d-Hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("d-Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexNAc(S)".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl, GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexNAc(S)".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexNAc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexNS".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Amino, GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexNS".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexN".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Amino],
                    proforma_name: Some("HexN".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexS".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexS".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexP".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("HexP".to_string()),
                    furanose: false,
                },
            ),
            (
                "Hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "Hep".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Heptose(None),
                    substituents: vec![],
                    proforma_name: Some("Hep".to_string()),
                    furanose: false,
                },
            ),
            (
                "Oct".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Octose,
                    substituents: vec![],
                    proforma_name: Some("Oct".to_string()),
                    furanose: false,
                },
            ),
            (
                "Non".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![],
                    proforma_name: Some("Non".to_string()),
                    furanose: false,
                },
            ),
            (
                "Dec".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Decose,
                    substituents: vec![],
                    proforma_name: Some("Dec".to_string()),
                    furanose: false,
                },
            ),
            (
                "Neu5Ac".to_string(),
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
                "Neu5Gc".to_string(),
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
                "Neu".to_string(),
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
                "Fuc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::Deoxy],
                    proforma_name: Some("Fuc".to_string()),
                    furanose: false,
                },
            ),
            // Single letter codes, by defining them like this they will be read but exported to the standard ProForma codes
            (
                "P".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("Hexphosphate".to_string()), // TODO: technically maybe not working when multiple are in there, think it through, should be two different elements,  both getting counts after them
                    furanose: false,
                },
            ),
            (
                "H".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "N".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                },
            ),
            (
                "F".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::Deoxy],
                    proforma_name: Some("Fuc".to_string()),
                    furanose: false,
                },
            ),
            (
                "S".to_string(),
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
                "A".to_string(),
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
                "G".to_string(),
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
