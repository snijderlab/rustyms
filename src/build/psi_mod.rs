use std::{ffi::OsString, fs, path::Path};

use crate::helper_functions::parse_molecular_formula_psi_mod;

use super::{
    obo::OboOntology,
    ontology_modification::{OntologyModification, PlacementRule, Position},
};

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

        let mut aa = Vec::new();
        let mut term = None;
        if let Some(values) = obj.lines.get("property_value") {
            for line in values {
                if line.starts_with("DiffFormula") {
                    modification.elements =
                        parse_molecular_formula_psi_mod(&line[13..line.len() - 12]).unwrap();
                } else if line.starts_with("Origin") {
                    aa = line[8..line.len() - 12]
                        .split(',')
                        .map(|s| s.trim().chars().next().unwrap())
                        .collect();
                    // TODO: parse the rules
                    //property_value: Origin "E" xsd:string
                    //property_value: TermSpec "N-term" xsd:string
                    // separate
                    //property_value: Origin "D, G" xsd:string
                    //even only on specific mods
                    //property_value: Origin "MOD:01801" xsd:string
                    // c: property_value: TermSpec "C-term" xsd:string
                } else if line.starts_with("TermSpec") {
                    if line[10..].starts_with("N-term") {
                        term = Some(Position::AnyNTerm);
                    } else if line[10..].starts_with("C-term") {
                        term = Some(Position::AnyCTerm);
                    } else {
                        panic!("Invalid TermSpec: {line}")
                    }
                }
            }
        }
        for aminoacid in &aa {
            modification.rules.push(PlacementRule::AminoAcid(
                *aminoacid,
                term.unwrap_or(Position::Anywhere),
            ));
        }
        if aa.is_empty() {
            if let Some(term) = term {
                modification.rules.push(PlacementRule::Terminal(term))
            }
        }
        mods.push(modification);
    }

    mods
}
