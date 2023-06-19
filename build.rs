use std::fmt::Display;
use std::fs;
use std::path::Path;
use std::{default, env};

use quick_xml::events::Event;
use quick_xml::reader::Reader;

// TODO: This now gives the modifications from unimod in a rustyms accessible rust
// file (unimod.rs, see build output for location). These modifications could be
// read in immediately by rustyms but there are some problems still to be solved:
// * There are unsupported elements in there (eg Na)
// * Some modifications in the unimod.xml uses glycans as there building blocks in the compose
//   (even worse some use both glycans and elements)
// * How to actually recognise a name when parsing a mod
//    - The modifications definitions list should contain the code (and ex code)
//    - Some modification code could contain a colon ':' making life harder

fn main() {
    let debug = env::var("DEBUG_BUILD").map(|v| v == "1").unwrap_or(false);
    let mods = parse_unimod(debug);

    let out_dir = env::var_os("OUT_DIR").unwrap();
    print(out_dir.as_os_str().to_str().unwrap(), debug);
    let dest_path = Path::new(&out_dir).join("unimod.rs");
    fs::write(
        dest_path,
        format!(
            "pub const UNIMOD_ONTOLOGY: &[(&str, &str, Modification)] = &[
                {}
                ];",
            mods.iter()
                .fold(String::new(), |acc, m| format!("{}\n{},", acc, m.to_code()))
        ),
    )
    .unwrap();
    println!("cargo:rerun-if-changed=databases/unimod.xml");
    println!("cargo:rerun-if-changed=build.rs");
}

fn parse_unimod(debug: bool) -> Vec<Modification> {
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
                    let mut modification = Modification::default();
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
                                                .unwrap()
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

fn print(text: impl AsRef<str>, debug: bool) {
    if debug {
        println!("cargo:warning={}", text.as_ref())
    }
}

#[derive(Debug, Default)]
struct Modification {
    formula: Vec<(Element, isize)>,
    code_name: String,
    ex_code_name: String,
    full_name: String,
    id: usize,
    rules: Vec<PlacementRule>,
}

impl Modification {
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
                    self.formula.push((
                        Element::try_from(last_name.as_str()).map_err(|_| ())?,
                        last_number.parse::<isize>().map_err(|_| ())? * minus,
                    ));
                    last_name.clear();
                    last_number.clear();
                    minus = 1;
                }
                b' ' => {
                    if !last_name.is_empty() {
                        self.formula
                            .push((Element::try_from(last_name.as_str()).map_err(|_| ())?, 1));
                        last_name.clear();
                    }
                }
                n if n.is_ascii_digit() => last_number.push(n as char),
                n if n.is_ascii_alphabetic() => last_name.push(n as char),
                _ => panic!("Weird formula composition"),
            }
        }
        if !last_name.is_empty() {
            self.formula
                .push((Element::try_from(last_name.as_str()).map_err(|_| ())?, 1));
            last_name.clear();
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
            "// {} [code name: {}] [ex_code_name: {}] rules: {}\n(\"{}\", \"{}\", Modification::Predefined(&[{}], &[]))",
            self.full_name,
            self.code_name,
            self.ex_code_name,
            self.rules.iter().fold(String::new(), |acc, r| format!("{acc}{r},")),
            self.code_name.to_ascii_lowercase(),
            self.ex_code_name.to_ascii_lowercase(),
            self.formula.iter().fold(String::new(), |acc, e| format!(
                "{}({}, {}),",
                acc, e.0, e.1
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

#[derive(Debug)]
pub enum Element {
    /// Element Hydrogen (H) atomic number: 1
    H = 1,
    /// Element Helium (He) atomic number: 2
    He,
    /// Element Lithium (Li) atomic number: 3
    Li,
    /// Element Beryllium (Be) atomic number: 4
    Be,
    /// Element Boron (B) atomic number: 5
    B,
    /// Element Carbon (C) atomic number: 6
    C,
    /// Element Nitrogen (N) atomic number: 7
    N,
    /// Element Oxygen (O) atomic number: 8
    O,
    /// Element Fluorine (F) atomic number: 9
    F,
    /// Element Neon (Ne) atomic number: 10
    Ne,
    /// Element Sodium (Na) atomic number: 11
    Na,
    /// Element Magnesium (Mg) atomic number: 12
    Mg,
    /// Element Aluminium (Al) atomic number: 13
    Al,
    /// Element Silicon (Si) atomic number: 14
    Si,
    /// Element Phosphorus (P) atomic number: 15
    P,
    /// Element Sulfur (S) atomic number: 16
    S,
    /// Element Chlorine (Cl) atomic number: 17
    Cl,
    /// Element Argon (Ar) atomic number: 18
    Ar,
    /// Element Potassium (K) atomic number: 19
    K,
    /// Element Calcium (Ca) atomic number: 20
    Ca,
    /// Element Scandium (Sc) atomic number: 21
    Sc,
    /// Element Titanium (Ti) atomic number: 22
    Ti,
    /// Element Vanadium (V) atomic number: 23
    V,
    /// Element Chromium (Cr) atomic number: 24
    Cr,
    /// Element Manganese (Mn) atomic number: 25
    Mn,
    /// Element Iron (Fe) atomic number: 26
    Fe,
    /// Element Cobalt (Co) atomic number: 27
    Co,
    /// Element Nickel (Ni) atomic number: 28
    Ni,
    /// Element Copper (Cu) atomic number: 29
    Cu,
    /// Element Zinc (Zn) atomic number: 30
    Zn,
    /// Element Gallium (Ga) atomic number: 31
    Ga,
    /// Element Germanium (Ge) atomic number: 32
    Ge,
    /// Element Arsenic (As) atomic number: 33
    As,
    /// Element Selenium (Se) atomic number: 34
    Se,
    /// Element Bromine (Br) atomic number: 35
    Br,
    /// Element Krypton (Kr) atomic number: 36
    Kr,
    /// Element Rubidium (Rb) atomic number: 37
    Rb,
    /// Element Strontium (Sr) atomic number: 38
    Sr,
    /// Element Yttrium (Y) atomic number: 39
    Y,
    /// Element Zirconium (Zr) atomic number: 40
    Zr,
    /// Element Niobium (Nb) atomic number: 41
    Nb,
    /// Element Molybdenum (Mo) atomic number: 42
    Mo,
    /// Element Technetium (Tc) atomic number: 43
    Tc,
    /// Element Ruthenium (Ru) atomic number: 44
    Ru,
    /// Element Rhodium (Rh) atomic number: 45
    Rh,
    /// Element Palladium (Pd) atomic number: 46
    Pd,
    /// Element Silver (Ag) atomic number: 47
    Ag,
    /// Element Cadmium (Cd) atomic number: 48
    Cd,
    /// Element Indium (In) atomic number: 49
    In,
    /// Element Tin (Sn) atomic number: 50
    Sn,
    /// Element Antimony (Sb) atomic number: 51
    Sb,
    /// Element Tellurium (Te) atomic number: 52
    Te,
    /// Element Iodine (I) atomic number: 53
    I,
    /// Element Xenon (Xe) atomic number: 54
    Xe,
    /// Element Caesium (Cs) atomic number: 55
    Cs,
    /// Element Barium (Ba) atomic number: 56
    Ba,
    /// Element Lanthanum (La) atomic number: 57
    La,
    /// Element Cerium (Ce) atomic number: 58
    Ce,
    /// Element Praseodymium (Pr) atomic number: 59
    Pr,
    /// Element Neodymium (Nd) atomic number: 60
    Nd,
    /// Element Promethium (Pm) atomic number: 61
    Pm,
    /// Element Samarium (Sm) atomic number: 62
    Sm,
    /// Element Europium (Eu) atomic number: 63
    Eu,
    /// Element Gadolinium (Gd) atomic number: 64
    Gd,
    /// Element Terbium (Tb) atomic number: 65
    Tb,
    /// Element Dysprosium (Dy) atomic number: 66
    Dy,
    /// Element Holmium (Ho) atomic number: 67
    Ho,
    /// Element Erbium (Er) atomic number: 68
    Er,
    /// Element Thulium (Tm) atomic number: 69
    Tm,
    /// Element Ytterbium (Yb) atomic number: 70
    Yb,
    /// Element Lutetium (Lu) atomic number: 71
    Lu,
    /// Element Hafnium (Hf) atomic number: 72
    Hf,
    /// Element Tantalum (Ta) atomic number: 73
    Ta,
    /// Element Tungsten (W) atomic number: 74
    W,
    /// Element Rhenium (Re) atomic number: 75
    Re,
    /// Element Osmium (Os) atomic number: 76
    Os,
    /// Element Iridium (Ir) atomic number: 77
    Ir,
    /// Element Platinum (Pt) atomic number: 78
    Pt,
    /// Element Gold (Au) atomic number: 79
    Au,
    /// Element Mercury (Hg) atomic number: 80
    Hg,
    /// Element Thallium (Tl) atomic number: 81
    Tl,
    /// Element Lead (Pb) atomic number: 82
    Pb,
    /// Element Bismuth (Bi) atomic number: 83
    Bi,
    /// Element Polonium (Po) atomic number: 84
    Po,
    /// Element Astatine (At) atomic number: 85
    At,
    /// Element Radon (Rn) atomic number: 86
    Rn,
    /// Element Francium (Fr) atomic number: 87
    Fr,
    /// Element Radium (Ra) atomic number: 88
    Ra,
    /// Element Actinium (Ac) atomic number: 89
    Ac,
    /// Element Thorium (Th) atomic number: 90
    Th,
    /// Element Protactinium (Pa) atomic number: 91
    Pa,
    /// Element Uranium (U) atomic number: 92
    U,
    /// Element Neptunium (Np) atomic number: 93
    Np,
    /// Element Plutonium (Pu) atomic number: 94
    Pu,
    /// Element Americium (Am) atomic number: 95
    Am,
    /// Element Curium (Cm) atomic number: 96
    Cm,
    /// Element Berkelium (Bk) atomic number: 97
    Bk,
    /// Element Californium (Cf) atomic number: 98
    Cf,
    /// Element Einsteinium (Es) atomic number: 99
    Es,
    /// Element Fermium (Fm) atomic number: 100
    Fm,
    /// Element Mendelevium (Md) atomic number: 101
    Md,
    /// Element Nobelium (No) atomic number: 102
    No,
    /// Element Lawrencium (Lr) atomic number: 103
    Lr,
    /// Element Rutherfordium (Rf) atomic number: 104
    Rf,
    /// Element Dubnium (Db) atomic number: 105
    Db,
    /// Element Seaborgium (Sg) atomic number: 106
    Sg,
    /// Element Bohrium (Bh) atomic number: 107
    Bh,
    /// Element Hassium (Hs) atomic number: 108
    Hs,
    /// Element Meitnerium (Mt) atomic number: 109
    Mt,
    /// Element Darmstadtium (Ds) atomic number: 110
    Ds,
    /// Element Roentgenium (Rg) atomic number: 111
    Rg,
    /// Element Copernicium (Cn) atomic number: 112
    Cn,
    /// Element Nihonium (Nh) atomic number: 113
    Nh,
    /// Element Flerovium (Fl) atomic number: 114
    Fl,
    /// Element Moscovium (Mc) atomic number: 115
    Mc,
    /// Element Livermorium (Lv) atomic number: 116
    Lv,
    /// Element Tennessine (Ts) atomic number: 117
    Ts,
    /// Element Oganesson (Og) atomic number: 118
    Og,
}

impl Display for Element {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::H => "Element::H",
                Self::He => "Element::He",
                Self::Li => "Element::Li",
                Self::Be => "Element::Be",
                Self::B => "Element::B",
                Self::C => "Element::C",
                Self::N => "Element::N",
                Self::O => "Element::O",
                Self::F => "Element::F",
                Self::Ne => "Element::Ne",
                Self::Na => "Element::Na",
                Self::Mg => "Element::Mg",
                Self::Al => "Element::Al",
                Self::Si => "Element::Si",
                Self::P => "Element::P",
                Self::S => "Element::S",
                Self::Cl => "Element::Cl",
                Self::Ar => "Element::Ar",
                Self::K => "Element::K",
                Self::Ca => "Element::Ca",
                Self::Sc => "Element::Sc",
                Self::Ti => "Element::Ti",
                Self::V => "Element::V",
                Self::Cr => "Element::Cr",
                Self::Mn => "Element::Mn",
                Self::Fe => "Element::Fe",
                Self::Co => "Element::Co",
                Self::Ni => "Element::Ni",
                Self::Cu => "Element::Cu",
                Self::Zn => "Element::Zn",
                Self::Ga => "Element::Ga",
                Self::Ge => "Element::Ge",
                Self::As => "Element::As",
                Self::Se => "Element::Se",
                Self::Br => "Element::Br",
                Self::Kr => "Element::Kr",
                Self::Rb => "Element::Rb",
                Self::Sr => "Element::Sr",
                Self::Y => "Element::Y",
                Self::Zr => "Element::Zr",
                Self::Nb => "Element::Nb",
                Self::Mo => "Element::Mo",
                Self::Tc => "Element::Tc",
                Self::Ru => "Element::Ru",
                Self::Rh => "Element::Rh",
                Self::Pd => "Element::Pd",
                Self::Ag => "Element::Ag",
                Self::Cd => "Element::Cd",
                Self::In => "Element::In",
                Self::Sn => "Element::Sn",
                Self::Sb => "Element::Sb",
                Self::Te => "Element::Te",
                Self::I => "Element::I",
                Self::Xe => "Element::Xe",
                Self::Cs => "Element::Cs",
                Self::Ba => "Element::Ba",
                Self::La => "Element::La",
                Self::Ce => "Element::Ce",
                Self::Pr => "Element::Pr",
                Self::Nd => "Element::Nd",
                Self::Pm => "Element::Pm",
                Self::Sm => "Element::Sm",
                Self::Eu => "Element::Eu",
                Self::Gd => "Element::Gd",
                Self::Tb => "Element::Tb",
                Self::Dy => "Element::Dy",
                Self::Ho => "Element::Ho",
                Self::Er => "Element::Er",
                Self::Tm => "Element::Tm",
                Self::Yb => "Element::Yb",
                Self::Lu => "Element::Lu",
                Self::Hf => "Element::Hf",
                Self::Ta => "Element::Ta",
                Self::W => "Element::W",
                Self::Re => "Element::Re",
                Self::Os => "Element::Os",
                Self::Ir => "Element::Ir",
                Self::Pt => "Element::Pt",
                Self::Au => "Element::Au",
                Self::Hg => "Element::Hg",
                Self::Tl => "Element::Tl",
                Self::Pb => "Element::Pb",
                Self::Bi => "Element::Bi",
                Self::Po => "Element::Po",
                Self::At => "Element::At",
                Self::Rn => "Element::Rn",
                Self::Fr => "Element::Fr",
                Self::Ra => "Element::Ra",
                Self::Ac => "Element::Ac",
                Self::Th => "Element::Th",
                Self::Pa => "Element::Pa",
                Self::U => "Element::U",
                Self::Np => "Element::Np",
                Self::Pu => "Element::Pu",
                Self::Am => "Element::Am",
                Self::Cm => "Element::Cm",
                Self::Bk => "Element::Bk",
                Self::Cf => "Element::Cf",
                Self::Es => "Element::Es",
                Self::Fm => "Element::Fm",
                Self::Md => "Element::Md",
                Self::No => "Element::No",
                Self::Lr => "Element::Lr",
                Self::Rf => "Element::Rf",
                Self::Db => "Element::Db",
                Self::Sg => "Element::Sg",
                Self::Bh => "Element::Bh",
                Self::Hs => "Element::Hs",
                Self::Mt => "Element::Mt",
                Self::Ds => "Element::Ds",
                Self::Rg => "Element::Rg",
                Self::Cn => "Element::Cn",
                Self::Nh => "Element::Nh",
                Self::Fl => "Element::Fl",
                Self::Mc => "Element::Mc",
                Self::Lv => "Element::Lv",
                Self::Ts => "Element::Ts",
                Self::Og => "Element::Og",
            }
        )
    }
}

impl TryFrom<&str> for Element {
    type Error = ();
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "H" => Ok(Self::H),
            "He" => Ok(Self::He),
            "Li" => Ok(Self::Li),
            "Be" => Ok(Self::Be),
            "B" => Ok(Self::B),
            "C" => Ok(Self::C),
            "N" => Ok(Self::N),
            "O" => Ok(Self::O),
            "F" => Ok(Self::F),
            "Ne" => Ok(Self::Ne),
            "Na" => Ok(Self::Na),
            "Mg" => Ok(Self::Mg),
            "Al" => Ok(Self::Al),
            "Si" => Ok(Self::Si),
            "P" => Ok(Self::P),
            "S" => Ok(Self::S),
            "Cl" => Ok(Self::Cl),
            "Ar" => Ok(Self::Ar),
            "K" => Ok(Self::K),
            "Ca" => Ok(Self::Ca),
            "Sc" => Ok(Self::Sc),
            "Ti" => Ok(Self::Ti),
            "V" => Ok(Self::V),
            "Cr" => Ok(Self::Cr),
            "Mn" => Ok(Self::Mn),
            "Fe" => Ok(Self::Fe),
            "Co" => Ok(Self::Co),
            "Ni" => Ok(Self::Ni),
            "Cu" => Ok(Self::Cu),
            "Zn" => Ok(Self::Zn),
            "Ga" => Ok(Self::Ga),
            "Ge" => Ok(Self::Ge),
            "As" => Ok(Self::As),
            "Se" => Ok(Self::Se),
            "Br" => Ok(Self::Br),
            "Kr" => Ok(Self::Kr),
            "Rb" => Ok(Self::Rb),
            "Sr" => Ok(Self::Sr),
            "Y" => Ok(Self::Y),
            "Zr" => Ok(Self::Zr),
            "Nb" => Ok(Self::Nb),
            "Mo" => Ok(Self::Mo),
            "Tc" => Ok(Self::Tc),
            "Ru" => Ok(Self::Ru),
            "Rh" => Ok(Self::Rh),
            "Pd" => Ok(Self::Pd),
            "Ag" => Ok(Self::Ag),
            "Cd" => Ok(Self::Cd),
            "In" => Ok(Self::In),
            "Sn" => Ok(Self::Sn),
            "Sb" => Ok(Self::Sb),
            "Te" => Ok(Self::Te),
            "I" => Ok(Self::I),
            "Xe" => Ok(Self::Xe),
            "Cs" => Ok(Self::Cs),
            "Ba" => Ok(Self::Ba),
            "La" => Ok(Self::La),
            "Ce" => Ok(Self::Ce),
            "Pr" => Ok(Self::Pr),
            "Nd" => Ok(Self::Nd),
            "Pm" => Ok(Self::Pm),
            "Sm" => Ok(Self::Sm),
            "Eu" => Ok(Self::Eu),
            "Gd" => Ok(Self::Gd),
            "Tb" => Ok(Self::Tb),
            "Dy" => Ok(Self::Dy),
            "Ho" => Ok(Self::Ho),
            "Er" => Ok(Self::Er),
            "Tm" => Ok(Self::Tm),
            "Yb" => Ok(Self::Yb),
            "Lu" => Ok(Self::Lu),
            "Hf" => Ok(Self::Hf),
            "Ta" => Ok(Self::Ta),
            "W" => Ok(Self::W),
            "Re" => Ok(Self::Re),
            "Os" => Ok(Self::Os),
            "Ir" => Ok(Self::Ir),
            "Pt" => Ok(Self::Pt),
            "Au" => Ok(Self::Au),
            "Hg" => Ok(Self::Hg),
            "Tl" => Ok(Self::Tl),
            "Pb" => Ok(Self::Pb),
            "Bi" => Ok(Self::Bi),
            "Po" => Ok(Self::Po),
            "At" => Ok(Self::At),
            "Rn" => Ok(Self::Rn),
            "Fr" => Ok(Self::Fr),
            "Ra" => Ok(Self::Ra),
            "Ac" => Ok(Self::Ac),
            "Th" => Ok(Self::Th),
            "Pa" => Ok(Self::Pa),
            "U" => Ok(Self::U),
            "Np" => Ok(Self::Np),
            "Pu" => Ok(Self::Pu),
            "Am" => Ok(Self::Am),
            "Cm" => Ok(Self::Cm),
            "Bk" => Ok(Self::Bk),
            "Cf" => Ok(Self::Cf),
            "Es" => Ok(Self::Es),
            "Fm" => Ok(Self::Fm),
            "Md" => Ok(Self::Md),
            "No" => Ok(Self::No),
            "Lr" => Ok(Self::Lr),
            "Rf" => Ok(Self::Rf),
            "Db" => Ok(Self::Db),
            "Sg" => Ok(Self::Sg),
            "Bh" => Ok(Self::Bh),
            "Hs" => Ok(Self::Hs),
            "Mt" => Ok(Self::Mt),
            "Ds" => Ok(Self::Ds),
            "Rg" => Ok(Self::Rg),
            "Cn" => Ok(Self::Cn),
            "Nh" => Ok(Self::Nh),
            "Fl" => Ok(Self::Fl),
            "Mc" => Ok(Self::Mc),
            "Lv" => Ok(Self::Lv),
            "Ts" => Ok(Self::Ts),
            "Og" => Ok(Self::Og),
            _ => Err(()),
        }
    }
}
