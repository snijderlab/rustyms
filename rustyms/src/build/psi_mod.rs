use std::{ffi::OsString, io::Write, path::Path};

use crate::{formula::MolecularFormula, ELEMENT_PARSE_LIST};

use super::{
    obo::OboOntology,
    ontology_modification::{OntologyList, OntologyModification, PlacementRule, Position},
};

pub fn build_psi_mod_ontology(out_dir: &OsString, debug: bool) {
    let mods = parse_psi_mod(debug);

    let dest_path = Path::new(&out_dir).join("psimod.dat");
    let mut file = std::fs::File::create(dest_path).unwrap();
    let final_mods = mods.into_iter().map(|m| m.into_mod()).collect::<Vec<_>>();
    file.write_all(&bincode::serialize::<OntologyList>(&final_mods).unwrap())
        .unwrap();
}

fn parse_psi_mod(_debug: bool) -> Vec<OntologyModification> {
    let obo =
        OboOntology::from_file("databases/PSI-MOD-newstyle.obo.gz").expect("Not a valid obo file");
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
            ontology: super::ontology_modification::Ontology::Psimod,
            ..Default::default()
        };

        let mut origins = Vec::new();
        let mut term = None;
        if let Some(values) = obj.lines.get("property_value") {
            for line in values {
                if line.starts_with("DiffFormula") {
                    modification.diff_formula =
                        parse_molecular_formula_psi_mod(&line[13..line.len() - 12]).unwrap();
                } else if line.starts_with("Origin") {
                    origins = line[8..line.len() - 12]
                        .split(',')
                        .map(|s| s.trim())
                        .collect();
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
        // If the list of possible origins contains "X" than the mod can be placed on any aminoacid
        // But if there is a TermSpec definition that should still be accounted for
        let all_aminoacids = origins.contains(&"X");
        if !all_aminoacids {
            for origin in &origins {
                if origin.len() == 1 {
                    modification.rules.push(PlacementRule::AminoAcid(
                        vec![(*origin).try_into().unwrap()],
                        term.unwrap_or(Position::Anywhere),
                    ));
                } else {
                    modification.rules.push(PlacementRule::PsiModification(
                        origin
                            .split_once(':')
                            .expect("Incorrect psi mod id, should contain a colon")
                            .1
                            .parse()
                            .expect("Incorrect psi mod id, should be numerical"),
                        term.unwrap_or(Position::Anywhere),
                    ));
                }
            }
        }
        if origins.is_empty() || all_aminoacids {
            if let Some(term) = term {
                modification.rules.push(PlacementRule::Terminal(term))
            }
        }
        mods.push(modification);
    }

    mods
}

// PSI-MOD: (12)C -5 (13)C 5 H 0 N 0 O 0 S 0
pub fn parse_molecular_formula_psi_mod(value: &str) -> Result<MolecularFormula, String> {
    let mut index = 0;
    let mut isotope = None;
    let mut element = None;
    let bytes = value.as_bytes();
    let mut result = MolecularFormula::default();
    while index < value.len() {
        match bytes[index] {
            b'(' if isotope.is_none() => {
                let len = bytes
                    .iter()
                    .skip(index)
                    .position(|c| *c == b')')
                    .ok_or(format!(
                        "No closing round bracket for round bracket at index: {index}"
                    ))?;
                isotope = Some(
                    value[index + 1..index + len]
                        .parse::<u16>()
                        .map_err(|e| e.to_string())?,
                );
                index += len + 1;
            }
            b'-' | b'0'..=b'9' if element.is_some() => {
                let (num, len) = std::str::from_utf8(
                    &bytes
                        .iter()
                        .skip(index)
                        .take_while(|c| c.is_ascii_digit() || **c == b'-')
                        .copied()
                        .collect::<Vec<_>>(),
                )
                .map(|v| (v.parse::<i16>().map_err(|e| e.to_string()), v.len()))
                .map_err(|e| e.to_string())?;
                let num = num?;
                if num != 0 {
                    assert!(result.add((element.unwrap(), isotope, num)));
                }
                element = None;
                isotope = None;
                index += len;
            }
            b' ' => index += 1,
            _ => {
                let mut found = false;
                for possible in ELEMENT_PARSE_LIST {
                    if value[index..].starts_with(possible.0) {
                        element = Some(possible.1);
                        index += possible.0.len();
                        found = true;
                        break;
                    }
                }
                if !found {
                    return Err(format!(
                        "Could not parse PSI-MOD elemental formula, broke down at index: {index}"
                    ));
                }
            }
        }
    }
    if isotope.is_some() || element.is_some() {
        Err("Last element missed a count".to_string())
    } else {
        Ok(result)
    }
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
