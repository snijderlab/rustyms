use std::{ffi::OsString, fs, path::Path};

use crate::{build::csv::parse_csv, glycan::GlycanStructure};

use super::obo::OboOntology;

pub fn build_gnome_ontology(out_dir: &OsString, debug: bool) {
    let mods = parse_gnome(debug);
    let structures = parse_gnome_structures(debug);
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

fn parse_gnome(_debug: bool) -> Vec<GNOmeModification> {
    let obo = OboOntology::from_file("databases/GNOme.obo").expect("Not a valid obo file");
    let mut mods = Vec::new();

    for obj in obo.objects {
        if obj.name != "Term" {
            continue;
        }
        let mut modification = GNOmeModification {
            code_name: obj.lines["name"][0].to_string(),
            ..Default::default()
        };

        //let mut origins = Vec::new();
        //let mut term = None;
        if let Some(values) = obj.lines.get("property_value") {
            for line in values {
                if line.starts_with("GNO:00000021") {
                    modification.mass = line[17..].to_string();
                } else if line.starts_with("GNO:00000033") {
                    modification.composition = line[17..].to_string();
                } else if line.starts_with("GNO:00000035") {
                    modification.topology = line[17..].to_string();
                }
            }
        }
        mods.push(modification);
    }

    mods
}

fn parse_gnome_structures(_debug: bool) -> Vec<(String, GlycanStructure)> {
    let mut glycans = Vec::new();
    for line in parse_csv("databases/glycosmos_glycans_list.csv")
        .unwrap()
        .skip(1)
    {
        let line = line.unwrap();
        glycans.push((
            line[0].to_string(),
            GlycanStructure::from_short_iupac(
                &line.line,
                line.fields[1].clone(),
                line.line_index + 1,
            )
            .unwrap(),
        ));
    }
    glycans
}

#[derive(Default, Debug)]
struct GNOmeModification {
    code_name: String,
    mass: String,
    composition: String,
    topology: String,
}

impl GNOmeModification {
    fn to_code(&self) -> String {
        // composition, placement rules, ontology, name, index (undefined), no place to store glycan composition for now
        format!(
            "Modification::Predefined({},&[],Ontology::Gnome,\"{}\",{})",
            self.composition, self.code_name, 0
        )
    }
}
