use std::{ffi::OsString, io::Write, iter, path::Path};

use regex::Regex;

use crate::{
    build::glycan::MonoSaccharide,
    formula::{Chemical, MolecularFormula},
    print, Element, NeutralLoss,
};

use super::{
    obo::OboOntology,
    ontology_modification::{OntologyModification, OntologyModificationList, PlacementRule},
    ModData,
};

pub fn build_unimod_ontology(out_dir: &OsString, debug: bool) {
    let mods = parse_unimod(debug);

    let dest_path = Path::new(&out_dir).join("unimod.dat");
    let mut file = std::fs::File::create(dest_path).unwrap();
    let final_mods = mods.into_iter().map(|m| m.into_mod()).collect::<Vec<_>>();
    file.write_all(&bincode::serialize::<OntologyModificationList>(&final_mods).unwrap())
        .unwrap();
}

fn parse_unimod(_debug: bool) -> Vec<OntologyModification> {
    let obo = OboOntology::from_file("databases/unimod.obo.gz").expect("Not a valid obo file");
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
            name: obj.lines["name"][0].to_string(),
            ontology: super::ontology_modification::Ontology::Unimod,
            ..Default::default()
        };
        if let Some(xref) = obj.lines.get("xref") {
            let re_position = Regex::new("spec_(\\d+)_position \"(.+)\"").unwrap();
            let re_site = Regex::new("spec_(\\d+)_site \"(.+)\"").unwrap();
            let re_neutral_loss =
                Regex::new("spec_(\\d+)_neutral_loss_\\d+_composition \"(.+)\"").unwrap();
            let mut mod_rules = Vec::new();
            for line in xref {
                if line.starts_with("delta_composition") {
                    modification
                        .with_unimod_composition(&line[19..line.len() - 1])
                        .expect("Invalid Unimod composition");
                    take = true;
                } else if let Some(groups) = re_position.captures(line) {
                    let index = groups.get(1).unwrap().as_str().parse::<usize>().unwrap() - 1;
                    let position = groups.get(2).unwrap().as_str().to_string();
                    if mod_rules.len() <= index {
                        mod_rules.extend(
                            iter::repeat((String::new(), String::new(), Vec::new()))
                                .take(index + 1 - mod_rules.len()),
                        );
                    }
                    mod_rules[index].1 = position;
                } else if let Some(groups) = re_site.captures(line) {
                    let index = groups.get(1).unwrap().as_str().parse::<usize>().unwrap() - 1;
                    let site = groups.get(2).unwrap().as_str().to_string();
                    if mod_rules.len() <= index {
                        mod_rules.extend(
                            iter::repeat((String::new(), String::new(), Vec::new()))
                                .take(index + 1 - mod_rules.len()),
                        );
                    }
                    mod_rules[index].0.push_str(&site);
                } else if let Some(groups) = re_neutral_loss.captures(line) {
                    let index = groups.get(1).unwrap().as_str().parse::<usize>().unwrap() - 1;
                    if !groups
                        .get(2)
                        .is_some_and(|g| g.is_empty() || g.as_str() == "0")
                    {
                        let formula =
                            parse_unimod_composition(groups.get(2).unwrap().as_str()).unwrap();
                        let loss = NeutralLoss::Loss(
                            formula.0
                                + formula
                                    .1
                                    .iter()
                                    .map(|(m, n)| m.formula() * *n)
                                    .sum::<MolecularFormula>(),
                        );
                        if mod_rules.len() <= index {
                            mod_rules.extend(
                                iter::repeat((String::new(), String::new(), Vec::new()))
                                    .take(index + 1 - mod_rules.len()),
                            );
                        }
                        mod_rules[index].2.push(loss);
                    }
                } else {
                    continue;
                }
            }
            if let ModData::Mod { specificities: rules, .. } = &mut modification.data {
                rules.extend(mod_rules.into_iter().filter_map(|rule| {
                    match (rule.0.as_str(), rule.1.as_str()) {
                        ("C-term", pos) => Some((
                            vec![PlacementRule::Terminal(pos.try_into().unwrap())],
                            rule.2,
                            Vec::new(),
                        )),
                        ("N-term", pos) => Some((
                            vec![PlacementRule::Terminal(pos.try_into().unwrap())],
                            rule.2,
                            Vec::new(),
                        )),
                        ("", "") => None,
                        (aa, pos) => Some((
                            vec![PlacementRule::AminoAcid(
                                aa.chars()
                                    .map(|c| {
                                        c.try_into().unwrap_or_else(|_| panic!("Not an AA: {c}"))
                                    })
                                    .collect(),
                                pos.try_into().unwrap(),
                            )],
                            rule.2,
                            Vec::new(),
                        )),
                    }
                }));
            }
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
        "ac" => Ok(Brick::Formula(molecular_formula!(C 2 H 2 O 1))),
        "me" => Ok(Brick::Formula(molecular_formula!(C 1 H 2))),
        "kdn" => Ok(Brick::Formula(molecular_formula!(C 9 H 14 O 8))),
        "kdo" => Ok(Brick::Formula(molecular_formula!(C 8 H 12 O 7))),
        "sulf" => Ok(Brick::Formula(molecular_formula!(S 1))),
        "water" => Ok(Brick::Formula(molecular_formula!(H 2 O 1))),
        _ => {
            if let Ok(el) = Element::try_from(name) {
                Ok(Brick::Element(el))
            } else if let Ok((ms, _)) = MonoSaccharide::from_short_iupac(name, 0, 0) {
                Ok(Brick::MonoSaccharide(ms))
            } else {
                print(format!("Could not parse unimod brick: `{name}`"), true);
                Err(())
            }
        }
    }
}

fn parse_unimod_composition(
    composition: &str,
) -> Result<(MolecularFormula, Vec<(MonoSaccharide, i32)>), ()> {
    assert!(composition.is_ascii());

    let mut formula = MolecularFormula::default();
    let mut monosaccharides = Vec::new();

    let mut last_name = String::new();
    let mut last_number = String::new();
    let mut sign = 1;
    for (index, c) in composition.bytes().enumerate() {
        match c {
            b'-' => sign = -1,
            b'(' => (),
            b')' => {
                let num = last_number.parse::<i32>().map_err(|_| ())? * sign;
                match parse_unimod_composition_brick(last_name.as_str()) {
                    Ok(Brick::Formula(f)) => formula += &(f * num),
                    Ok(Brick::Element(e)) => assert!(formula.add((e, None, num))),
                    Ok(Brick::MonoSaccharide(m)) => monosaccharides.push((m, num)),
                    Err(()) => return Err(()),
                }
                last_name.clear();
                last_number.clear();
                sign = 1;
            }
            b' ' => {
                if !last_name.is_empty() {
                    match parse_unimod_composition_brick(last_name.as_str()) {
                        Ok(Brick::Formula(f)) => formula += &f,
                        Ok(Brick::Element(e)) => assert!(formula.add((e, None, 1))),
                        Ok(Brick::MonoSaccharide(m)) => monosaccharides.push((m, 1)),
                        Err(()) => return Err(()),
                    }
                    last_name.clear();
                }
            }
            n if n.is_ascii_digit() => last_number.push(n as char),
            n if n.is_ascii_alphabetic() => last_name.push(n as char),
            _ => panic!("Unimod composition parsing broke at: '{composition}' (specifically at '{c}' index {index})"),
        }
    }
    if !last_name.is_empty() {
        match parse_unimod_composition_brick(last_name.as_str()) {
            Ok(Brick::Formula(f)) => formula += &f,
            Ok(Brick::Element(e)) => assert!(formula.add((e, None, 1))),
            Ok(Brick::MonoSaccharide(m)) => monosaccharides.push((m, 1)),
            Err(()) => return Err(()),
        }
    }
    Ok((formula, monosaccharides))
}

impl OntologyModification {
    #[deny(clippy::unwrap_used)]
    fn with_unimod_composition(&mut self, composition: &str) -> Result<(), ()> {
        let (diff_formula, sugars) = parse_unimod_composition(composition)?;
        self.formula = diff_formula;
        if let ModData::Mod {
            monosaccharides, ..
        } = &mut self.data
        {
            monosaccharides.extend(sugars);
        }
        Ok(())
    }
}
