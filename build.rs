use std::collections::HashMap;
use std::env;
use std::ffi::OsString;
use std::fmt::Display;
use std::fs;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::iter;
use std::path::Path;

#[path = "./src/shared/element.rs"]
mod element;
#[macro_use]
#[path = "./src/shared/formula.rs"]
mod formula;
#[path = "./src/shared/glycan.rs"]
mod glycan;

#[path = "./src/helper_functions.rs"]
mod helper_functions;

pub use crate::element::*;
use crate::formula::Chemical;
use crate::formula::*;
use crate::glycan::*;
use crate::helper_functions::*;

use regex::Regex;

fn main() {
    let debug = env::var("DEBUG_BUILD").map(|v| v == "1").unwrap_or(false);

    let out_dir = env::var_os("OUT_DIR").unwrap();
    build_unimod_ontology(&out_dir, debug);
    build_psi_mod_ontology(&out_dir, debug);
    build_atomic_masses(&out_dir, debug).unwrap();

    println!("cargo:rerun-if-changed=databases/unimod.obo");
    println!("cargo:rerun-if-changed=databases/PSI-MOD-newstyle.obo");
    println!("cargo:rerun-if-changed=databases/IUPAC-atomic-masses.csv");
    println!("cargo:rerun-if-changed=databases/CIAAW-atomic-weights.csv");
    println!("cargo:rerun-if-changed=databases/CIAAW-isotopic-abundances.csv");
    println!("cargo:rerun-if-changed=build.rs");
    print(out_dir.as_os_str().to_str().unwrap(), debug);
}

fn build_unimod_ontology(out_dir: &OsString, debug: bool) {
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

fn build_psi_mod_ontology(out_dir: &OsString, debug: bool) {
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

fn parse_psi_mod(debug: bool) -> Vec<OntologyModification> {
    let obo =
        OboOntology::from_file("databases/PSI-MOD-newstyle.obo").expect("Not a valid obo file");
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
            context: "PSI-MOD".to_string(),
            ..Default::default()
        };

        if let Some(values) = obj.lines.get("property_value") {
            for line in values {
                if line.starts_with("DiffFormula") {
                    match parse_molecular_formula_psi_mod(&line[13..line.len() - 12]) {
                        Ok(o) => {
                            modification.elements = o;
                            take = true;
                        }
                        Err(_e) => {
                            print(
                                format!(
                                    "PSI-MOD: could not match formula: \"{}\"",
                                    &line[13..line.len() - 12]
                                ),
                                debug,
                            );
                        }
                    }
                } else if line.starts_with("Origin") {
                    // TODO: parse the rules
                }
            }
        }
        if take {
            mods.push(modification);
        }
    }

    mods
}

fn print(text: impl AsRef<str>, debug: bool) {
    if debug {
        println!("cargo:warning={}", text.as_ref())
    }
}

#[derive(Debug, Default)]
struct OntologyModification {
    elements: MolecularFormula,
    monosaccharides: Vec<(MonoSaccharide, i16)>,
    code_name: String,
    full_name: String,
    context: String,
    id: usize,
    rules: Vec<PlacementRule>,
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

    fn to_code(&self) -> String {
        format!(
            "// {} [code name: {}] rules: {}\n({}, \"{}\", Modification::Predefined(&MolecularFormula::new(&[{}]), &[{}], \"{}\", \"{}\"))",
            self.full_name,
            self.code_name,
            self.rules
                .iter()
                .fold(String::new(), |acc, r| format!("{acc}{r},")),
            self.id,
            self.code_name.to_ascii_lowercase(),
            self.monosaccharides.iter().fold(self.elements.clone(), |acc, m| acc + m.0.formula() * m.1).elements().iter()
            .fold(String::new(), |acc, (e, i, n)| format!(
                "{acc}(Element::{e},{i},{n}),",
            )),
            self.rules
                .iter()
                .fold(String::new(), |acc, e| format!(
                    "{}PlacementRule::{},",
                    acc,
                    e
                )),
                self.context,
                self.code_name,
        )
    }
}

#[derive(Debug)]
enum PlacementRule {
    AminoAcid(char, Position),
    Terminal(Position),
}

impl Display for PlacementRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AminoAcid(aminoacid, position) => {
                write!(
                    f,
                    "AminoAcid(AminoAcid::{}, Position::{:?})",
                    aminoacid, position
                )
            }
            Self::Terminal(t) => write!(f, "Terminal(Position::{:?})", t),
        }
    }
}

#[derive(Debug, Default, PartialEq, Eq)]
enum Position {
    #[default]
    Undefined = 1,
    Anywhere,
    AnyNTerm,
    AnyCTerm,
    ProteinNTerm,
    ProteinCTerm,
}

impl TryInto<Position> for u8 {
    type Error = String;
    fn try_into(self) -> Result<Position, Self::Error> {
        match self {
            b'1' => Ok(Position::Undefined),
            b'2' => Ok(Position::Anywhere),
            b'3' => Ok(Position::AnyNTerm),
            b'4' => Ok(Position::AnyCTerm),
            b'5' => Ok(Position::ProteinNTerm),
            b'6' => Ok(Position::ProteinCTerm),
            n => Err(format!("Outside range: {n}")),
        }
    }
}

impl TryInto<Position> for &str {
    type Error = String;
    fn try_into(self) -> Result<Position, Self::Error> {
        match self {
            "" => Ok(Position::Undefined),
            "Anywhere" => Ok(Position::Anywhere),
            "Any N-term" => Ok(Position::AnyNTerm),
            "Any C-term" => Ok(Position::AnyCTerm),
            "Protein N-term" => Ok(Position::ProteinNTerm),
            "Protein C-term" => Ok(Position::ProteinCTerm),
            n => Err(format!("Not valid position: {n}")),
        }
    }
}

#[derive(Debug, Default)]
struct OboOntology {
    headers: Vec<(String, String)>,
    objects: Vec<OboObject>,
}

#[derive(Debug, Default)]
struct OboObject {
    name: String,
    lines: HashMap<String, Vec<String>>,
}

impl OboOntology {
    fn from_file(path: &str) -> Result<OboOntology, String> {
        let reader = BufReader::new(File::open(path).map_err(|e| e.to_string())?);
        let mut obo = OboOntology::default();
        let mut recent_obj = None;

        for line in reader.lines() {
            let line = line.map_err(|e| e.to_string())?.trim_end().to_string();
            if line.is_empty() {
                continue;
            } else if line.starts_with('[') && line.ends_with(']') {
                if let Some(obj) = recent_obj {
                    obo.objects.push(obj);
                }
                recent_obj = Some(OboObject::new(&line[1..=line.len() - 2]));
            } else if let Some((name, value)) = line.split_once(':') {
                if let Some(obj) = &mut recent_obj {
                    obj.lines
                        .entry(name.trim().to_string())
                        .or_insert(Vec::new())
                        .push(value.trim().to_string());
                } else {
                    obo.headers.push((name.to_string(), value.to_string()));
                }
            } else {
                return Err(format!("Invalid line in obo file: {line}"));
            }
        }
        if let Some(obj) = recent_obj {
            obo.objects.push(obj);
        }
        Ok(obo)
    }
}

impl OboObject {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            ..Self::default()
        }
    }
}

fn parse_csv(path: &str) -> Result<Vec<Vec<String>>, String> {
    let reader = BufReader::new(File::open(path).map_err(|e| e.to_string())?);
    let mut table = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?.trim_end().to_string();
        if line.is_empty() {
            continue;
        }
        let mut enclosed = None;
        let mut row = Vec::new();
        let mut cell = String::new();
        for ch in line.chars() {
            match (ch, enclosed) {
                ('\"', None) | ('\'', None) => enclosed = Some(ch),
                (c, Some(e)) if c == e => enclosed = None,
                (',', None) => {
                    row.push(cell.trim().to_string());
                    cell.clear();
                }
                _ => cell.push(ch),
            }
        }
        row.push(cell.trim().to_string());
        table.push(row);
    }
    Ok(table)
}

fn build_atomic_masses(out_dir: &OsString, _debug: bool) -> Result<(), String> {
    let mut atomic_weights = vec![None; 118];
    let mut isotopic_abundances = vec![Vec::new(); 118];
    let mut atomic_masses = vec![Vec::new(); 118];
    // TODO: write output, read isotopic abundances, average weights, and determine normal monoisotopic mass

    let table = parse_csv("databases/IUPAC-atomic-masses.csv")?;
    for line in table {
        let (nuclide, mass, _uncertainty, year) = (&line[0], &line[1], &line[2], &line[3]);
        if nuclide.starts_with("AME")
            || nuclide.is_empty()
            || nuclide == "nuclide"
            || !year.ends_with("2020</a>")
        {
            continue;
        }
        let isotope = nuclide
            .trim_end_matches(|c: char| c.is_alphabetic())
            .parse::<usize>()
            .map_err(|e| e.to_string())?;
        let element = Element::try_from(nuclide.trim_start_matches(|c: char| c.is_ascii_digit()))
            .map_err(|_| {
            format!("Not a valid isotope+element, could not recognise element: {nuclide}")
        })?;
        let mass = mass.parse::<f64>().map_err(|e| e.to_string())?;
        atomic_masses[element as usize - 1].push((isotope, mass))
    }

    let mut last_element = 0;
    let table = parse_csv("databases/CIAAW-isotopic-abundances.csv")?;
    for line in table {
        let (element, _element, _name, isotope, abundance, _note) =
            (&line[0], &line[1], &line[2], &line[3], &line[4], &line[5]);
        let mut abundance = abundance.to_owned();
        abundance.retain(|c| !c.is_whitespace()); // Remove any whitespace, including any sneaking non breaking spaces
        if element == "Z" || isotope.starts_with("[") || abundance == "-" {
            continue;
        }

        let element = if element.is_empty() {
            last_element
        } else {
            last_element = element
                .parse::<usize>()
                .map_err(|_| format!("Not a valid number for element Z: {element}"))?;
            last_element
        };

        let isotope = isotope.parse::<usize>().unwrap();

        isotopic_abundances[element - 1].push((isotope, get_ciaaw_number(&abundance)?))
    }

    let table = parse_csv("databases/CIAAW-atomic-weights.csv")?;
    for line in table {
        let (element, weight) = (&line[0], &line[3]);
        let mut weight = weight.to_owned();
        weight.retain(|c| !c.is_whitespace()); // Remove any whitespace, including any sneaking non breaking spaces
        if element == "Z" || weight == "—" {
            continue;
        }
        let element = element
            .parse::<usize>()
            .map_err(|_| format!("Not valid element number (not a number): {element}"))?;
        atomic_weights[element - 1] = Some(get_ciaaw_number(&weight)?);
    }

    // Monoisotopic, average weight, isotopes: (num, mass, abundance)
    let combined_data = isotopic_abundances
        .into_iter()
        .zip(atomic_masses)
        .zip(atomic_weights)
        .map(|((isotopic_abundance, isotopic_mass), atomic_weight)| {
            let isotopes = isotopic_mass
                .iter()
                .map(|(i, m)| {
                    (
                        *i,
                        *m,
                        isotopic_abundance
                            .iter()
                            .find(|(ai, _)| ai == i)
                            .map(|(_, a)| *a)
                            .unwrap_or(0.0),
                    )
                })
                .collect::<Vec<_>>();
            (
                isotopes
                    .iter()
                    .fold(
                        (None, 0.0),
                        |acc, (_, m, a)| if *a > acc.1 { (Some(*m), *a) } else { acc },
                    )
                    .0,
                atomic_weight,
                isotopes,
            )
        });

    let dest_path = Path::new(&out_dir).join("elements.rs");
    fs::write(
        dest_path,
        format!(
            "pub const ELEMENTAL_DATA: &[(Option<f64>, Option<f64>, &[(u16, f64, f64)])] = &[
                {}
                ];",
            combined_data.fold(String::new(), |acc, element| format!(
                "{}\n({:?}, {:?}, &[{}]),",
                acc,
                element.0,
                element.1,
                element.2.iter().fold(String::new(), |acc, i| format!(
                    "{}({},{:.1},{:.1}),",
                    acc, i.0, i.1, i.2
                ))
            ))
        ),
    )
    .unwrap();

    Ok(())
}

fn get_ciaaw_number(text: &str) -> Result<f64, String> {
    let parse = |t: &str| {
        t.parse::<f64>()
            .map_err(|_| format!("Not a valid number: {t}"))
    };
    if text.starts_with('[') {
        let (low, high) = text[1..text.len() - 1]
            .split_once(',')
            .ok_or(format!("Not a valid range: {text}"))?;
        Ok((parse(low)? + parse(high)?) / 2.0)
    } else if text.ends_with(')') {
        Ok(parse(
            text.split_once('(')
                .ok_or(format!("Not valid error indication: {text}"))?
                .0,
        )?)
    } else {
        Ok(parse(text)?)
    }
}
