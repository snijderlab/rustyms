use std::collections::HashMap;
use std::env;
use std::ffi::OsString;
use std::fmt::Display;
use std::fs;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

#[path = "./src/shared/element.rs"]
mod element;
#[path = "./src/shared/glycan.rs"]
mod glycan;

#[path = "./src/helper_functions.rs"]
mod helper_functions;

use crate::element::*;
use crate::glycan::*;
use crate::helper_functions::*;

use quick_xml::events::Event;
use quick_xml::reader::Reader;

fn main() {
    let debug = env::var("DEBUG_BUILD").map(|v| v == "1").unwrap_or(false);

    let out_dir = env::var_os("OUT_DIR").unwrap();
    print(out_dir.as_os_str().to_str().unwrap(), debug);
    build_unimod_ontology(&out_dir, debug);
    build_psi_mod_ontology(&out_dir, debug);

    println!("cargo:rerun-if-changed=databases/unimod.xml");
    println!("cargo:rerun-if-changed=databases/PSI-MOD-newstyle.obo");
    println!("cargo:rerun-if-changed=build.rs");
}

fn build_unimod_ontology(out_dir: &OsString, debug: bool) {
    let mods = parse_unimod(debug);
    let dest_path = Path::new(&out_dir).join("unimod.rs");
    fs::write(
        dest_path,
        format!(
            "pub const UNIMOD_ONTOLOGY: &[(usize, &str, &str, Modification)] = &[
                {}
                ];",
            mods.iter()
                .fold(String::new(), |acc, m| format!("{}\n{},", acc, m.to_code()))
        ),
    )
    .unwrap();
}

fn parse_unimod(debug: bool) -> Vec<OntologyModification> {
    let mut reader =
        Reader::from_file("databases/unimod.xml").expect("Unimod definition file not present");
    let mut buf = Vec::new();
    let mut mods = Vec::new();
    reader.trim_text(true);

    // The `Reader` does not implement `Iterator` because it outputs borrowed data (`Cow`s)
    loop {
        // empty
        // alt_names_row -> empty
        // amino_acids_row -> ignore contains all amino acids
        // brick2element -> ignore contains all elements
        // brick contains simple chemicals
        // classifications
        // elements_row
        // fragment_comp_row
        // fragments_row
        // mod2brick_row

        // start+stop
        // modifications_row
        match reader.read_event_into(&mut buf) {
            Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
            // exits the loop when reaching end of file
            Ok(Event::Eof) => break,

            Ok(Event::Start(e)) => match e.name().as_ref() {
                b"modifications_row" => {
                    let mut modification = OntologyModification::default();
                    let mut skip = false;
                    for attr in e.attributes() {
                        match attr {
                            Err(e) => {
                                panic!("Error at position {}: {:?}", reader.buffer_position(), e)
                            }
                            Ok(attribute) => match attribute.key.into_inner() {
                                b"full_name" => {
                                    modification.with_full_name(attribute.value.as_ref())
                                }
                                b"code_name" => {
                                    modification.with_code_name(attribute.value.as_ref())
                                }
                                b"ex_code_name" => {
                                    modification.with_ex_code_name(attribute.value.as_ref())
                                }
                                b"record_id" => {
                                    modification.id =
                                        String::from_utf8(attribute.value.as_ref().to_vec())
                                            .expect("Record_id not valid utf8")
                                            .parse()
                                            .expect("Record_id not valid number");
                                }
                                b"composition" => {
                                    if modification.with_formula(attribute.value.as_ref()).is_err()
                                    {
                                        skip = true;
                                        print(
                                            format!(
                                                "Could not parse: {}",
                                                String::from_utf8(
                                                    attribute.value.as_ref().to_vec()
                                                )
                                                .unwrap(),
                                            ),
                                            debug,
                                        );
                                    }
                                }
                                _ => (),
                            },
                        }
                    }
                    if !skip {
                        mods.push(modification);
                    }
                }
                b"specificity_row" => {
                    let mut place = String::new();
                    let mut position = Position::Undefined;
                    let mut mod_key = None;
                    for attr in e.attributes() {
                        match attr {
                            Err(e) => {
                                panic!("Error at position {}: {:?}", reader.buffer_position(), e)
                            }
                            Ok(attribute) => match attribute.key.into_inner() {
                                b"one_letter" => {
                                    place = String::from_utf8(attribute.value.as_ref().to_vec())
                                        .expect("Invalid utf8");
                                }
                                b"position_key" => {
                                    assert!(attribute.value.as_ref().len() == 1);
                                    position = attribute.value.as_ref()[0]
                                        .try_into()
                                        .expect("Invalid position for placement rule");
                                }
                                b"mod_key" => {
                                    mod_key = Some(
                                        String::from_utf8(attribute.value.as_ref().to_vec())
                                            .expect("Record_id not valid utf8")
                                            .parse()
                                            .expect("Record_id not valid number"),
                                    );
                                }
                                _ => (),
                            },
                        }
                    }
                    #[allow(clippy::option_map_unit_fn)]
                    if let Some(mod_key) = mod_key {
                        let rule = match (place.as_str(), position) {
                            (aa, position) if aa.len() == 1 => {
                                PlacementRule::AminoAcid(aa.chars().next().unwrap(), position)
                            }
                            ("C-term", position)
                                if position == Position::AnyCTerm
                                    || position == Position::ProteinCTerm =>
                            {
                                PlacementRule::Terminal(position)
                            }
                            ("C-term", Position::Anywhere) => {
                                PlacementRule::Terminal(Position::AnyCTerm)
                            }
                            ("N-term", position)
                                if position == Position::AnyNTerm
                                    || position == Position::ProteinNTerm =>
                            {
                                PlacementRule::Terminal(position)
                            }
                            ("N-term", Position::Anywhere) => {
                                PlacementRule::Terminal(Position::AnyNTerm)
                            }
                            (place, position) => {
                                panic!("Invalid placement rule: {place} {position:?} mod_key {mod_key}")
                            }
                        };
                        mods.iter_mut()
                            .find(|m| m.id == mod_key)
                            .map(|m| m.rules.push(rule));
                    }
                }
                _ => (),
            },
            Ok(Event::Text(_e)) => (), //print(format!("{:?}", e), debug),

            // There are several other `Event`s we do not consider here
            _ => (),
        }
        // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
    }

    mods
}

enum Brick {
    Element(Element),
    Formula(Vec<(Element, isize)>),
    MonoSaccharide(MonoSaccharide),
}

fn parse_unimod_composition_brick(name: &str) -> Result<Brick, ()> {
    match name.to_lowercase().as_str() {
        "ac" => Ok(Brick::Formula(vec![
            (Element::O, 1),
            (Element::C, 2),
            (Element::H, 2),
        ])), // TODO: check formula
        "me" => Ok(Brick::Formula(vec![(Element::C, 1), (Element::H, 2)])),
        "kdn" => Ok(Brick::Formula(vec![
            (Element::C, 9),
            (Element::H, 14),
            (Element::O, 8),
        ])),
        "kdo" => Ok(Brick::Formula(vec![
            (Element::C, 8),
            (Element::H, 12),
            (Element::O, 7),
        ])),
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
            "pub const PSI_MOD_ONTOLOGY: &[(usize, &str, &str, Modification)] = &[
                {}
                ];",
            mods.iter()
                .fold(String::new(), |acc, m| format!("{}\n{},", acc, m.to_code()))
        ),
    )
    .unwrap();
}

// TODO: get the name/code out instead of the number, parse the number, provide a way for the user to select modifications based on number (also for unimod) (described in 4.2.2)
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
            ..Default::default()
        };

        if let Some(values) = obj.lines.get("property_value") {
            for line in values {
                if line.starts_with("DiffFormula") {
                    match parse_named_counter(&line[13..line.len() - 12], &ELEMENT_PARSE_LIST, true)
                    {
                        Ok(o) => {
                            modification.elements = o;
                            take = true;
                        }
                        Err(_e) => {
                            print(
                                format!(
                                    "Could not match {obj:?} Formula: \"{}\"",
                                    &line[13..line.len() - 12]
                                ),
                                debug,
                            );
                        }
                    }
                } else if line.starts_with("Origin") {
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
    elements: Vec<(Element, isize)>,
    monosaccharides: Vec<(MonoSaccharide, isize)>,
    code_name: String,
    ex_code_name: String,
    full_name: String,
    id: usize,
    rules: Vec<PlacementRule>,
}

impl OntologyModification {
    #[deny(clippy::unwrap_used)]
    fn with_formula(&mut self, composition: &[u8]) -> Result<(), ()> {
        let mut last_name = String::new();
        let mut last_number = String::new();
        let mut minus = 1;
        for c in composition.iter() {
            match *c {
                b'-' => minus = -1,
                b'(' => (),
                b')' => {
                    let num = last_number.parse::<isize>().map_err(|_| ())? * minus;
                    match parse_unimod_composition_brick(last_name.as_str()) {
                        Ok(Brick::Formula(f)) => self
                            .elements
                            .extend(f.into_iter().map(|pair| (pair.0, pair.1 * num))),
                        Ok(Brick::Element(e)) => self.elements.push((e, num)),
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
                            Ok(Brick::Formula(f)) => self.elements.extend(f),
                            Ok(Brick::Element(e)) => self.elements.push((e, 1)),
                            Ok(Brick::MonoSaccharide(m)) => self.monosaccharides.push((m, 1)),
                            Err(()) => return Err(()),
                        }
                        last_name.clear();
                    }
                }
                n if n.is_ascii_digit() => last_number.push(n as char),
                n if n.is_ascii_alphabetic() => last_name.push(n as char),
                _ => panic!("Weird formula composition"),
            }
        }
        if !last_name.is_empty() {
            match parse_unimod_composition_brick(last_name.as_str()) {
                Ok(Brick::Formula(f)) => self.elements.extend(f),
                Ok(Brick::Element(e)) => self.elements.push((e, 1)),
                Ok(Brick::MonoSaccharide(m)) => self.monosaccharides.push((m, 1)),
                Err(()) => return Err(()),
            }
        }
        Ok(())
    }

    fn with_code_name(&mut self, name: &[u8]) {
        self.code_name = String::from_utf8(name.to_vec())
            .unwrap()
            .replace("&gt;", ">")
    }

    fn with_ex_code_name(&mut self, name: &[u8]) {
        self.ex_code_name = String::from_utf8(name.to_vec())
            .unwrap()
            .replace("&gt;", ">")
    }

    fn with_full_name(&mut self, name: &[u8]) {
        self.full_name = String::from_utf8(name.to_vec())
            .unwrap()
            .replace("&gt;", ">")
    }

    fn to_code(&self) -> String {
        format!(
            "// {} [code name: {}] [ex_code_name: {}] rules: {}\n({}, \"{}\", \"{}\", Modification::Predefined(&[{}], &[{}]))",
            self.full_name,
            self.code_name,
            self.ex_code_name,
            self.rules.iter().fold(String::new(), |acc, r| format!("{acc}{r},")),
            self.id,
            self.code_name.to_ascii_lowercase(),
            self.ex_code_name.to_ascii_lowercase(),
            self.elements.iter().fold(String::new(), |acc, e| format!(
                "{}(Element::{}, {}),",
                acc, e.0, e.1
            )),
            self.monosaccharides.iter().fold(String::new(), |acc, e| format!(
                "{}(MonoSaccharide::{}, {}),",
                acc, e.0.to_string().replace('-', "_"), e.1
            ))
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
                write!(f, "({}, {:?})", aminoacid, position)
            }
            Self::Terminal(t) => write!(f, "{:?}", t),
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
        let mut reader = BufReader::new(File::open(path).map_err(|e| e.to_string())?);
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
