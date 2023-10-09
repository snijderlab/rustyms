use std::{ffi::OsString, fs, path::Path};

use crate::{formula::MolecularFormula, ELEMENT_PARSE_LIST};

use super::{
    obo::OboOntology,
    ontology_modification::{OntologyModification, PlacementRule, Position},
};

pub fn build_gnome_ontology(out_dir: &OsString, debug: bool) {
    let mods = parse_gnome(debug);
    let dest_path = Path::new(&out_dir).join("psi-mod.rs");
    fs::write(
        dest_path,
        format!(
            "pub const GNOME_ONTOLOGY: &[(usize, &str, Modification)] = &[
                {}
                ];",
            mods.iter()
                .fold(String::new(), |acc, m| format!("{}\n{},", acc, m.to_code()))
        ),
    )
    .unwrap();
}

fn parse_gnome(_debug: bool) -> (Vec<GNOmeModification>, Vec<GlycanStructure>) {
    let obo = OboOntology::from_file("databases/GNOme.obo").expect("Not a valid obo file");
    let mut mods = Vec::new();
    let mut glycans = Vec::new();

    for obj in obo.objects {
        if obj.name != "Term" {
            continue;
        }
        let mut modification = GNOmeModification {
            code_name: obj.lines["name"][0].to_string(),
            ..Default::default()
        };

        let mut origins = Vec::new();
        let mut term = None;
        if let Some(values) = obj.lines.get("property_value") {
            for line in values {
                if line.starts_with("GNO:00000021") {
                    modification.mass = line[17..line.len() - 12].to_string();
                } else if line.starts_with("GNO:00000033") {
                    modification.composition = line[17..line.len() - 12].to_string();
                } else if line.starts_with("GNO:00000035") {
                    modification.topology = line[17..line.len() - 12].to_string();
                }
            }
        }
        mods.push(modification);
    }

    mods
}

struct GNOmeModification {
    code_name: String,
    mass: String,
    composition: String,
    topology: String,
}
