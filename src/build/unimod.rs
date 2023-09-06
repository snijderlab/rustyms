use std::{ffi::OsString, fs, iter, path::Path};

use regex::Regex;

use crate::{formula::MolecularFormula, glycan::MonoSaccharide, Element};

use super::{
    obo::OboOntology,
    ontology_modification::{OntologyModification, PlacementRule},
};

pub fn build_unimod_ontology(out_dir: &OsString, debug: bool) {
    let mods = parse_unimod(debug);
    let dest_path = Path::new(&out_dir).join("unimod.rs");
    fs::write(
        dest_path,
        format!(
            "pub const UNIMOD_ONTOLOGY: &[(usize, &str, Modification)] = &[
                {}
                ];",
            mods.iter()
                .fold(String::new(), |acc, m| format!("{}\n{},", acc, m.to_code()))
        ),
    )
    .unwrap();
}

fn parse_unimod(_debug: bool) -> Vec<OntologyModification> {
    let obo = OboOntology::from_file("databases/unimod.obo").expect("Not a valid obo file");
    let mut mods = Vec::new();

    for obj in obo.objects {
        if obj.name != "Term" {
            continue;
        }
        let mut take = false;
        let mut modification = OntologyModification {
            id: obj.lines["id"][0]
                .split_once(':')
                .expect("Incorrect psi mod id, should contain a colon")
                .1
                .parse()
                .expect("Incorrect psi mod id, should be numerical"),
            code_name: obj.lines["name"][0].to_string(),
            context: "Unimod".to_string(),
            ..Default::default()
        };
        if let Some(xref) = obj.lines.get("xref") {
            let re_position = Regex::new("spec_(\\d+)_position \"(.+)\"").unwrap();
            let re_site = Regex::new("spec_(\\d+)_site \"(.+)\"").unwrap();
            let mut rules = Vec::new();
            for line in xref {
                if line.starts_with("delta_composition") {
                    modification
                        .with_unimod_composition(&line[19..line.len() - 1])
                        .expect("Invalid Unimod composition");
                    take = true;
                } else if let Some(groups) = re_position.captures(line) {
                    let index = groups.get(1).unwrap().as_str().parse::<usize>().unwrap() - 1;
                    let position = groups.get(2).unwrap().as_str().to_string();
                    if rules.len() <= index {
                        rules.extend(
                            iter::repeat((String::new(), String::new()))
                                .take(index + 1 - rules.len()),
                        );
                    }
                    rules[index].1 = position;
                } else if let Some(groups) = re_site.captures(line) {
                    let index = groups.get(1).unwrap().as_str().parse::<usize>().unwrap() - 1;
                    let site = groups.get(2).unwrap().as_str().to_string();
                    if rules.len() <= index {
                        rules.extend(
                            iter::repeat((String::new(), String::new()))
                                .take(index + 1 - rules.len()),
                        );
                    }
                    rules[index].0 = site;
                } else {
                    continue;
                }
            }
            modification.rules = rules
                .into_iter()
                .filter_map(|rule| match (rule.0.as_str(), rule.1.as_str()) {
                    ("C-term", pos) => Some(PlacementRule::Terminal(pos.try_into().unwrap())),
                    ("N-term", pos) => Some(PlacementRule::Terminal(pos.try_into().unwrap())),
                    (aa, pos) if aa.len() == 1 => Some(PlacementRule::AminoAcid(
                        aa.chars().next().unwrap(),
                        pos.try_into().unwrap(),
                    )),
                    ("", "") => None,
                    (_, _) => panic!("Invalid rule definition: {rule:?}"),
                })
                .collect();
        }
        if take {
            mods.push(modification);
        }
    }

    mods
}

enum Brick {
    Element(Element),
    Formula(MolecularFormula),
    MonoSaccharide(MonoSaccharide),
}

fn parse_unimod_composition_brick(name: &str) -> Result<Brick, ()> {
    match name.to_lowercase().as_str() {
        "ac" => Ok(Brick::Formula(MolecularFormula::new(&[
            (Element::O, 0, 1),
            (Element::C, 0, 2),
            (Element::H, 0, 2),
        ]))), // TODO: check formula
        "me" => Ok(Brick::Formula(MolecularFormula::new(&[
            (Element::C, 0, 1),
            (Element::H, 0, 2),
        ]))),
        "kdn" => Ok(Brick::Formula(MolecularFormula::new(&[
            (Element::C, 0, 9),
            (Element::H, 0, 14),
            (Element::O, 0, 8),
        ]))),
        "kdo" => Ok(Brick::Formula(MolecularFormula::new(&[
            (Element::C, 0, 8),
            (Element::H, 0, 12),
            (Element::O, 0, 7),
        ]))),
        _ => {
            if let Ok(el) = Element::try_from(name) {
                Ok(Brick::Element(el))
            } else if let Ok(ms) = MonoSaccharide::try_from(name) {
                Ok(Brick::MonoSaccharide(ms))
            } else {
                Err(())
            }
        }
    }
}

impl OntologyModification {
    #[deny(clippy::unwrap_used)]
    fn with_unimod_composition(&mut self, composition: &str) -> Result<(), ()> {
        let mut last_name = String::new();
        let mut last_number = String::new();
        let mut minus = 1;
        for c in composition.bytes() {
            match c {
                b'-' => minus = -1,
                b'(' => (),
                b')' => {
                    let num = last_number.parse::<i16>().map_err(|_| ())? * minus;
                    match parse_unimod_composition_brick(last_name.as_str()) {
                        Ok(Brick::Formula(f)) => self.elements += &(f * num),
                        Ok(Brick::Element(e)) => self.elements.add((e, 0, num)),
                        Ok(Brick::MonoSaccharide(m)) => self.monosaccharides.push((m, num)),
                        Err(()) => return Err(()),
                    }
                    last_name.clear();
                    last_number.clear();
                    minus = 1;
                }
                b' ' => {
                    if !last_name.is_empty() {
                        match parse_unimod_composition_brick(last_name.as_str()) {
                            Ok(Brick::Formula(f)) => self.elements += &f,
                            Ok(Brick::Element(e)) => self.elements.add((e, 0, 1)),
                            Ok(Brick::MonoSaccharide(m)) => self.monosaccharides.push((m, 1)),
                            Err(()) => return Err(()),
                        }
                        last_name.clear();
                    }
                }
                n if n.is_ascii_digit() => last_number.push(n as char),
                n if n.is_ascii_alphabetic() => last_name.push(n as char),
                _ => panic!("Weird formula composition: {composition}"),
            }
        }
        if !last_name.is_empty() {
            match parse_unimod_composition_brick(last_name.as_str()) {
                Ok(Brick::Formula(f)) => self.elements += &f,
                Ok(Brick::Element(e)) => self.elements.add((e, 0, 1)),
                Ok(Brick::MonoSaccharide(m)) => self.monosaccharides.push((m, 1)),
                Err(()) => return Err(()),
            }
        }
        Ok(())
    }
}
