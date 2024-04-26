use std::{ffi::OsString, io::Write, path::Path};

use itertools::Itertools;
use ordered_float::OrderedFloat;

use crate::{formula::MolecularFormula, Linker, LinkerSpecificity};

use super::{
    obo::OboOntology,
    ontology_modification::{OntologyList, OntologyModification, PlacementRule, Position},
};

pub fn build_xlmod_ontology(out_dir: &OsString, debug: bool) {
    let mods = parse_xlmod(debug);

    let dest_path = Path::new(&out_dir).join("xlmod.dat");
    let mut file = std::fs::File::create(dest_path).unwrap();
    let final_mods = mods.0.into_iter().map(|m| m.into_mod()).collect::<Vec<_>>();
    file.write_all(&bincode::serialize::<OntologyList>(&final_mods).unwrap())
        .unwrap();
}

fn parse_xlmod(_debug: bool) -> (Vec<OntologyModification>, Vec<Linker>) {
    let obo = OboOntology::from_file("databases/XLMOD.obo.gz").expect("Not a valid obo file");
    let mut mods = Vec::new();
    let mut linkers = Vec::new();

    for obj in obo.objects {
        if obj.name != "Term" {
            continue;
        }
        let id = obj.lines["id"][0]
            .split_once(':')
            .expect("Incorrect XLMOD id, should contain a colon")
            .1
            .parse()
            .expect("Incorrect XLMOD id, should be numerical");
        let name = obj.lines["name"][0].to_string();

        let mut sites = None;
        let mut length = None;
        let mut mass = None;
        let mut formula = None;
        let mut origins = (Vec::new(), Vec::new());
        if let Some(values) = obj.lines.get("property_value") {
            for line in values {
                if line.starts_with("reactionSites") {
                    // reactionSites: "2" xsd:nonNegativeInteger
                    sites = Some((line[16..line.len() - 24]).parse::<u8>().unwrap());
                } else if line.starts_with("spacerLength") {
                    // spacerLength: "5.0" xsd:float
                    if !line[15..line.len() - 11].contains('-') {
                        length = Some(
                            (line[15..line.len() - 11])
                                .parse::<ordered_float::OrderedFloat<f64>>()
                                .unwrap(),
                        );
                    }
                } else if line.starts_with("monoIsotopicMass") {
                    // monoIsotopicMass: "535.15754" xsd:double
                    mass = Some(
                        (line[19..line.len() - 12])
                            .parse::<ordered_float::OrderedFloat<f64>>()
                            .unwrap(),
                    );
                } else if line.starts_with("deadEndFormula") {
                    // deadEndFormula: "C8 H12 O3" xsd:string
                    sites = Some(1);
                    formula =
                        Some(MolecularFormula::from_xlmod(&line[17..line.len() - 12]).unwrap());
                } else if line.starts_with("bridgeFormula") {
                    // bridgeFormula: "C7 D10 H2 N4" xsd:string
                    sites = Some(2);
                    formula =
                        Some(MolecularFormula::from_xlmod(&line[16..line.len() - 12]).unwrap());
                } else if line.starts_with("specificities") {
                    // specificities: "(C,U)" xsd:string
                    // specificities: "(K,N,Q,R,Protein N-term)&(E,D,Protein C-term)" xsd:string
                    if let Some((l, r)) = line[17..line.len() - 13].split_once('&') {
                        origins = (
                            l.trim_matches(['(', ')'])
                                .split(',')
                                .map(|s| s.trim())
                                .collect(),
                            r.trim_matches(['(', ')'])
                                .split(',')
                                .map(|s| s.trim())
                                .collect(),
                        );
                    } else {
                        origins = (
                            line[17..line.len() - 13]
                                .trim_matches(['(', ')'])
                                .split(',')
                                .map(|s| s.trim())
                                .collect(),
                            Vec::new(),
                        );
                    }
                } else if line.starts_with("secondarySpecificities") {
                    // secondarySpecificities: "(S,T,Y)" xsd:string
                    origins.0.extend(
                        line[26..line.len() - 13]
                            .trim_matches(['(', ')'])
                            .split(',')
                            .map(|s| s.trim()),
                    );
                }
            }
        }
        let origins = (
            read_placement_rules(&origins.0),
            read_placement_rules(&origins.1),
        );
        if let Some(mass) = mass {
            if let Some(formula) = &mut formula {
                formula.add_mass(mass);
            } else {
                formula = Some(MolecularFormula::with_additional_mass(mass.0))
            }
        }
        if sites == Some(2) || !origins.1.is_empty() {
            linkers.push(Linker {
                specificities: if !origins.1.is_empty() {
                    LinkerSpecificity::Asymmetric(origins.0, origins.1)
                } else {
                    LinkerSpecificity::Symmetric(origins.0)
                },
                formula: formula.unwrap_or_default(),
                name,
                id,
                length,
                ontology: super::Ontology::Xlmod,
            });
        } else if sites == Some(3) {
            continue; // Ignore
        } else {
            mods.push(OntologyModification {
                diff_formula: formula.unwrap_or_default(),
                monosaccharides: Vec::new(),
                code_name: name.clone(),
                full_name: name,
                ontology: super::Ontology::Xlmod,
                id,
                rules: vec![(origins.0, Vec::new(), Vec::new())],
            });
        }
    }

    (mods, linkers)
}

fn read_placement_rules(bricks: &[&str]) -> Vec<PlacementRule> {
    bricks
        .iter()
        .filter_map(|brick| {
            if brick.len() == 1 {
                Some(PlacementRule::AminoAcid(
                    vec![(*brick).try_into().unwrap()],
                    Position::Anywhere,
                ))
            } else if *brick == "Protein N-term" {
                Some(PlacementRule::Terminal(Position::ProteinNTerm))
            } else if *brick == "Protein C-term" {
                Some(PlacementRule::Terminal(Position::ProteinCTerm))
            } else if ["Thy"].contains(brick) {
                None
            } else {
                panic!("Invalid placement rule: '{brick}'")
            }
        })
        .collect_vec()
}

#[cfg(test)]
mod tests {
    #[test]
    fn parse_molecular_formula() {
        assert_eq!(
            parse_molecular_formula_psi_mod("(12)C -5 (13)C 5 H 0 N 0 O 0 S 0").unwrap(),
            molecular_formula!((12)C -5 (13)C 5 H 0 N 0 O 0 S 0)
        );
        assert_eq!(
            parse_molecular_formula_psi_mod("(12)C -9 (13)C 9").unwrap(),
            molecular_formula!((12)C -9 (13)C 9)
        );
    }
}
