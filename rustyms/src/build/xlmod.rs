use std::{ffi::OsString, io::Write, path::Path};

use crate::formula::MolecularFormula;

use super::{
    obo::OboOntology,
    ontology_modification::{OntologyList, OntologyModification, PlacementRule, Position},
};

pub fn build_xlmod_ontology(out_dir: &OsString, debug: bool) {
    let mods = parse_xlmod(debug);

    let dest_path = Path::new(&out_dir).join("xlmod.dat");
    let mut file = std::fs::File::create(dest_path).unwrap();
    let final_mods = mods.into_iter().map(|m| m.into_mod()).collect::<Vec<_>>();
    file.write_all(&bincode::serialize::<OntologyList>(&final_mods).unwrap())
        .unwrap();
}

fn parse_xlmod(_debug: bool) -> Vec<OntologyModification> {
    let obo = OboOntology::from_file("databases/XLMOD.obo.gz").expect("Not a valid obo file");
    let mut mods = Vec::new();

    for obj in obo.objects {
        if obj.name != "Term" {
            continue;
        }
        let mut modification = OntologyModification {
            id: obj.lines["id"][0]
                .split_once(':')
                .expect("Incorrect XLMOD id, should contain a colon")
                .1
                .parse()
                .expect("Incorrect XLMOD id, should be numerical"),
            code_name: obj.lines["name"][0].to_string(),
            ontology: super::ontology_modification::Ontology::Xlmod,
            ..Default::default()
        };

        let mut sites = None;
        let mut length = None;
        let mut mass = None;
        let mut formula = None;
        let mut origins = Vec::new();
        if let Some(values) = obj.lines.get("property_value") {
            for line in values {
                if line.starts_with("reactionSites") {
                    // reactionSites: "2" xsd:nonNegativeInteger
                    sites = Some((line[16..line.len() - 24]).parse::<u8>().unwrap());
                } else if line.starts_with("spacerLength") {
                    // spacerLength: "5.0" xsd:float
                    length = Some((line[15..line.len() - 11]).parse::<f64>().unwrap());
                } else if line.starts_with("monoIsotopicMass") {
                    // monoIsotopicMass: "535.15754" xsd:double
                    mass = Some((line[19..line.len() - 12]).parse::<f64>().unwrap());
                } else if line.starts_with("deadEndFormula") {
                    // deadEndFormula: "C8 H12 O3" xsd:string
                    // deadEndFormula: "-C1 -H2 O1" xsd:string
                    // TODO: Handle 'D'
                    formula =
                        Some(MolecularFormula::from_xlmod(&line[17..line.len() - 12]).unwrap());
                } else if line.starts_with("bridgeFormula") {
                    // bridgeFormula: "C7 D10 H2 N4" xsd:string
                    // bridgeFormula: "13C6 H6 O2" xsd:string
                    formula =
                        Some(MolecularFormula::from_xlmod(&line[17..line.len() - 12]).unwrap());
                } else if line.starts_with("specificities") {
                    // specificities: "(C,U)" xsd:string
                    // specificities: "(K,N,Q,R,Protein N-term)&(E,D,Protein C-term)" xsd:string
                    origins.extend(line[17..line.len() - 13].split(',').map(|s| s.trim()));
                } else if line.starts_with("secondarySpecificities") {
                    // secondarySpecificities: "(S,T,Y)" xsd:string
                    origins.extend(line[26..line.len() - 13].split(',').map(|s| s.trim()));
                }
            }
        }
        // If the list of possible origins contains "X" than the mod can be placed on any aminoacid
        // But if there is a TermSpec definition that should still be accounted for
        let all_aminoacids = origins.contains(&"X");
        if !all_aminoacids {
            for origin in &origins {
                if origin.len() == 1 {
                    modification.rules.push((
                        vec![PlacementRule::AminoAcid(
                            vec![(*origin).try_into().unwrap()],
                            Position::Anywhere,
                        )],
                        Vec::new(),
                        Vec::new(),
                    ));
                } else if *origin == "Protein N-term" {
                    modification.rules.push((
                        vec![PlacementRule::Terminal(Position::ProteinNTerm)],
                        Vec::new(),
                        Vec::new(),
                    ));
                } else if *origin == "Protein C-term" {
                    modification.rules.push((
                        vec![PlacementRule::Terminal(Position::ProteinCTerm)],
                        Vec::new(),
                        Vec::new(),
                    ));
                }
            }
        }
        mods.push(modification);
    }

    mods
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
