use std::{collections::HashMap, ffi::OsString, io::BufWriter, io::Write, path::Path};

use itertools::Itertools;

use crate::{build::csv::parse_csv, glycan::*};

use super::{obo::OboOntology, GnoComposition, Modification, ToCode};

pub fn build_gnome_ontology(out_dir: &OsString, debug: bool) {
    // Get all the basic info
    let mut mods = parse_gnome(debug);
    let read_mods = mods.clone();
    let structures = parse_gnome_structures(debug);

    // Fill all known info points
    for modification in mods.values_mut() {
        if modification.mass.is_none() {
            modification.mass = find_mass(&read_mods, modification.is_a.clone());
        }
        if let Some(id) = &modification.topology {
            modification.structure = structures.get(id).cloned();
        }
    }

    // Write out the data
    let dest_path = Path::new(&out_dir).join("gnome.dat");
    let mut file = std::fs::File::create(dest_path).unwrap();
    let final_mods = mods
        .into_values()
        .filter(|m| m.mass.is_some())
        .take(10)
        .map(|m| (0_usize, m.code_name.clone(), m.to_mod()))
        .collect::<Vec<_>>();
    file.write_all(&bincode::serialize(&final_mods).unwrap())
        .unwrap();
    // let mut writer = BufWriter::new(file);
    // writeln!(
    //     writer,
    //     "pub const GNOME_ONTOLOGY: &[(usize, &str, Modification)] = &["
    // )
    // .unwrap();
    // for modification in mods.values().filter(|m| m.mass.is_some()) {
    //     writeln!(writer, "{},", modification.to_code()).unwrap();
    // }
    // writeln!(writer, "];").unwrap();
}

fn find_mass(mods: &HashMap<String, GNOmeModification>, mut name: String) -> Option<f64> {
    let mut mass = None;
    while mass.is_none() {
        mass = mods.get(&name)?.mass;
        name = mods[&name].is_a.clone();
    }
    mass
}

fn parse_gnome(_debug: bool) -> HashMap<String, GNOmeModification> {
    let obo = OboOntology::from_file("databases/GNOme.obo.gz").expect("Not a valid obo file");
    let mut mods = HashMap::new();

    for obj in obo.objects {
        if obj.name != "Term" || !obj.lines.contains_key("is_a") {
            continue;
        }
        // name: glycan of molecular weight 40.03 Da
        let name = &obj.lines["name"][0];
        let mut modification = GNOmeModification {
            code_name: obj.lines["id"][0][4..].to_lowercase(),
            is_a: obj.lines["is_a"][0]
                .split_once('!')
                .map(|(a, _)| a.trim()[4..].to_lowercase())
                .unwrap(),
            mass: if name.len() > 30 {
                name[27..name.len() - 3].parse::<f64>().ok()
            } else {
                None
            },
            ..Default::default()
        };

        if let Some(values) = obj.lines.get("property_value") {
            for line in values {
                if line.starts_with("GNO:00000035") {
                    modification.topology = Some(line[17..].to_lowercase());
                }
            }
        }
        mods.insert(modification.code_name.clone(), modification);
    }

    mods
}

fn parse_gnome_structures(_debug: bool) -> HashMap<String, GlycanStructure> {
    let mut glycans = HashMap::new();
    let mut errors = 0;
    for line in parse_csv("databases/glycosmos_glycans_list.csv.gz", b',')
        .unwrap()
        .skip(1)
    {
        let line = line.unwrap();
        if !line[1].is_empty() {
            match GlycanStructure::from_short_iupac(
                &line.line,
                line.fields[1].clone(),
                line.line_index + 1,
            ) {
                Ok(glycan) => {
                    glycans.insert(line[0].to_lowercase(), glycan);
                }
                Err(error) => {
                    if errors < 5 {
                        println!("{error}");
                    }
                    errors += 1;
                }
            }
        }
    }
    if errors > 0 {
        panic!(
            "Total glycan structure reading errors: {errors} total read {}",
            glycans.len()
        );
    }
    glycans
}

#[derive(Default, Clone, Debug)]
struct GNOmeModification {
    code_name: String,                  // id of current item
    is_a: String,                       // id to prev in chain
    topology: Option<String>,           // id of topology
    mass: Option<f64>,                  // mass if defined
    structure: Option<GlycanStructure>, // mass if defined
}

impl GNOmeModification {
    fn to_mod(self) -> Modification {
        if let Some(structure) = self.structure {
            Modification::Gno(GnoComposition::Structure(structure), self.code_name)
        } else if let Some(mass) = self.mass {
            Modification::Gno(GnoComposition::Mass(mass), self.code_name)
        } else {
            panic!("unreachable")
        }
    }
}

impl ToCode for GNOmeModification {
    fn to_code(&self) -> String {
        let comp = if let Some(structure) = &self.structure {
            if let GlycanStructure::Allocated(a) = structure {
                format!(
                    "GnoComposition::Structure(GlycanStructure::Predefined({}))",
                    a.to_code()
                )
            } else {
                panic!("Should not have predefined glycans at build time")
            }
        } else {
            format!("GnoComposition::Mass({:.2})", self.mass.unwrap())
        };
        format!(
            "(0, \"{1}\", Modification::Gno({},\"{}\"))",
            comp, self.code_name
        )
    }
}

impl ToCode for AllocatedGlycanStructure {
    fn to_code(&self) -> String {
        format!(
            "PredefinedGlycanStructure{{sugar: {}, branches: &[{}]}}",
            self.sugar.to_code(),
            self.branches.iter().map(ToCode::to_code).join(",")
        )
    }
}

impl ToCode for MonoSaccharide {
    fn to_code(&self) -> String {
        if let MonoSaccharide::Allocated(a) = self {
            format!(
                "PredefinedMonosaccharide {{
            base_sugar: Some({}),
            substituents: &[{}],
            pro_forma_name: None,
            furanose: Some({}),
        }}",
                a.base_sugar.to_code(),
                a.substituents.iter().map(ToCode::to_code).join(","),
                a.furanose
            )
        } else {
            panic!("Should not have predefined glycans at build time")
        }
    }
}

impl ToCode for BaseSugar {
    fn to_code(&self) -> String {
        match self {
            Self::Sugar => "BaseSugar::Sugar".to_string(),
            Self::Triose => "BaseSugar::Triose".to_string(),
            Self::Tetrose(iso) => format!(
                "BaseSugar::Tetrose({})",
                iso.as_ref()
                    .map_or("None".to_string(), |i| format!("Some({})", i.to_code()))
            ),
            Self::Pentose(iso) => format!(
                "BaseSugar::Pentose({})",
                iso.as_ref()
                    .map_or("None".to_string(), |i| format!("Some({})", i.to_code()))
            ),
            Self::Hexose(iso) => format!(
                "BaseSugar::Hexose({})",
                iso.as_ref()
                    .map_or("None".to_string(), |i| format!("Some({})", i.to_code()))
            ),
            Self::Heptose(iso) => format!(
                "BaseSugar::Heptose({})",
                iso.as_ref()
                    .map_or("None".to_string(), |i| format!("Some({})", i.to_code()))
            ),
            Self::Octose => "BaseSugar::Octose".to_string(),
            Self::Nonose => "BaseSugar::Nonose".to_string(),
            Self::Decose => "BaseSugar::Decose".to_string(),
        }
    }
}

impl ToCode for TetroseIsomer {
    fn to_code(&self) -> String {
        match self {
            Self::Erythrose => "TetroseIsomer::Erythrose",
            Self::Threose => "TetroseIsomer::Threose",
        }
        .to_string()
    }
}
impl ToCode for PentoseIsomer {
    fn to_code(&self) -> String {
        match self {
            Self::Ribose => "PentoseIsomer::Ribose",
            Self::Arabinose => "PentoseIsomer::Arabinose",
            Self::Xylose => "PentoseIsomer::Xylose",
            Self::Lyxose => "PentoseIsomer::Lyxose",
            Self::Xylulose => "PentoseIsomer::Xylulose",
        }
        .to_string()
    }
}
impl ToCode for HexoseIsomer {
    fn to_code(&self) -> String {
        match self {
            Self::Glucose => "HexoseIsomer::Glucose",
            Self::Galactose => "HexoseIsomer::Galactose",
            Self::Mannose => "HexoseIsomer::Mannose",
            Self::Allose => "HexoseIsomer::Allose",
            Self::Altrose => "HexoseIsomer::Altrose",
            Self::Gulose => "HexoseIsomer::Gulose",
            Self::Idose => "HexoseIsomer::Idose",
            Self::Talose => "HexoseIsomer::Talose",
            Self::Psicose => "HexoseIsomer::Psicose",
            Self::Fructose => "HexoseIsomer::Fructose",
            Self::Sorbose => "HexoseIsomer::Sorbose",
            Self::Tagatose => "HexoseIsomer::Tagatose",
        }
        .to_string()
    }
}
impl ToCode for HeptoseIsomer {
    fn to_code(&self) -> String {
        match self {
            Self::GlyceroMannoHeptopyranose => "HeptoseIsomer::GlyceroMannoHeptopyranose",
            Self::Sedoheptulose => "HeptoseIsomer::Sedoheptulose",
        }
        .to_string()
    }
}
impl ToCode for GlycanSubstituent {
    fn to_code(&self) -> String {
        match self {
            Self::Acetimidoyl => "GlycanSubstituent::Acetimidoyl".to_string(),
            Self::Acetyl => "GlycanSubstituent::Acetyl".to_string(),
            Self::AcetylAlanyl => "GlycanSubstituent::AcetylAlanyl".to_string(),
            Self::AcetylGlutaminyl => "GlycanSubstituent::AcetylGlutaminyl".to_string(),
            Self::Acid => "GlycanSubstituent::Acid".to_string(),
            Self::Alanyl => "GlycanSubstituent::Alanyl".to_string(),
            Self::Alcohol => "GlycanSubstituent::Alcohol".to_string(),
            Self::Amino => "GlycanSubstituent::Amino".to_string(),
            Self::Aric => "GlycanSubstituent::Aric".to_string(),
            Self::CargoxyEthylidene => "GlycanSubstituent::CargoxyEthylidene".to_string(),
            Self::Deoxy => "GlycanSubstituent::Deoxy".to_string(),
            Self::Didehydro => "GlycanSubstituent::Didehydro".to_string(),
            Self::DiHydroxyButyryl => "GlycanSubstituent::DiHydroxyButyryl".to_string(),
            Self::DiMethyl => "GlycanSubstituent::DiMethyl".to_string(),
            Self::DiMethylAcetimidoyl => "GlycanSubstituent::DiMethylAcetimidoyl".to_string(),
            Self::DiMethylGlyceryl => "GlycanSubstituent::DiMethylGlyceryl".to_string(),
            Self::Element(element) => format!("GlycanSubstituent::Element(Element::{element})"),
            Self::Ethanolamine => "GlycanSubstituent::Ethanolamine".to_string(),
            Self::EtOH => "GlycanSubstituent::EtOH".to_string(),
            Self::Formyl => "GlycanSubstituent::Formyl".to_string(),
            Self::Glyceryl => "GlycanSubstituent::Glyceryl".to_string(),
            Self::Glycolyl => "GlycanSubstituent::Glycolyl".to_string(),
            Self::Glycyl => "GlycanSubstituent::Glycyl".to_string(),
            Self::HydroxyButyryl => "GlycanSubstituent::HydroxyButyryl".to_string(),
            Self::HydroxyMethyl => "GlycanSubstituent::HydroxyMethyl".to_string(),
            Self::Lac => "GlycanSubstituent::Lac".to_string(),
            Self::Lactyl => "GlycanSubstituent::Lactyl".to_string(),
            Self::Methyl => "GlycanSubstituent::Methyl".to_string(),
            Self::MethylAcetimidoyl => "GlycanSubstituent::MethylAcetimidoyl".to_string(),
            Self::MethylGlutamyl => "GlycanSubstituent::MethylGlutamyl".to_string(),
            Self::NAcetyl => "GlycanSubstituent::NAcetyl".to_string(),
            Self::NDiMe => "GlycanSubstituent::NDiMe".to_string(),
            Self::NFo => "GlycanSubstituent::NFo".to_string(),
            Self::NGlycolyl => "GlycanSubstituent::NGlycolyl".to_string(),
            Self::OCarboxyEthyl => "GlycanSubstituent::OCarboxyEthyl".to_string(),
            Self::PCholine => "GlycanSubstituent::PCholine".to_string(),
            Self::Phosphate => "GlycanSubstituent::Phosphate".to_string(),
            Self::Pyruvyl => "GlycanSubstituent::Pyruvyl".to_string(),
            Self::Suc => "GlycanSubstituent::Suc".to_string(),
            Self::Sulfate => "GlycanSubstituent::Sulfate".to_string(),
            Self::Tauryl => "GlycanSubstituent::Tauryl".to_string(),
            Self::Ulo => "GlycanSubstituent::Ulo".to_string(),
            Self::Ulof => "GlycanSubstituent::Ulof".to_string(),
            Self::Water => "GlycanSubstituent::Water".to_string(),
        }
    }
}
