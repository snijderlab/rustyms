use std::{collections::HashMap, ffi::OsString, io::Write, path::Path};

use crate::{build::csv::parse_csv, glycan::*};

use super::{obo::OboOntology, GnoComposition, Modification};

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
