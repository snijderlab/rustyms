use std::{ffi::OsString, fs, path::Path};

use crate::helper_functions::parse_molecular_formula_psi_mod;

use super::{obo::OboOntology, ontology_modification::OntologyModification};

pub fn build_psi_mod_ontology(out_dir: &OsString, debug: bool) {
    let mods = parse_psi_mod(debug);
    let dest_path = Path::new(&out_dir).join("psi-mod.rs");
    fs::write(
        dest_path,
        format!(
            "pub const PSI_MOD_ONTOLOGY: &[(usize, &str, Modification)] = &[
                {}
                ];",
            mods.iter()
                .fold(String::new(), |acc, m| format!("{}\n{},", acc, m.to_code()))
        ),
    )
    .unwrap();
}

fn parse_psi_mod(_debug: bool) -> Vec<OntologyModification> {
    let obo =
        OboOntology::from_file("databases/PSI-MOD-newstyle.obo").expect("Not a valid obo file");
    let mut mods = Vec::new();

    for obj in obo.objects {
        if obj.name != "Term" {
            continue;
        }
        let mut modification = OntologyModification {
            id: obj.lines["id"][0]
                .split_once(':')
                .expect("Incorrect psi mod id, should contain a colon")
                .1
                .parse()
                .expect("Incorrect psi mod id, should be numerical"),
            code_name: obj.lines["name"][0].to_string(),
            context: "PSI-MOD".to_string(),
            ..Default::default()
        };

        if let Some(values) = obj.lines.get("property_value") {
            for line in values {
                if line.starts_with("DiffFormula") {
                    modification.elements =
                        parse_molecular_formula_psi_mod(&line[13..line.len() - 12]).unwrap();
                } else if line.starts_with("Origin") {
                    // TODO: parse the rules
                }
            }
        }
        mods.push(modification);
    }

    mods
}
