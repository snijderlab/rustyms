use std::fmt::Display;

use crate::{
    formula::{Chemical, MolecularFormula},
    glycan::MonoSaccharide,
};

#[derive(Debug, Default)]
pub struct OntologyModification {
    pub elements: MolecularFormula,
    pub monosaccharides: Vec<(MonoSaccharide, i16)>,
    pub code_name: String,
    pub full_name: String,
    pub context: String,
    pub id: usize,
    pub rules: Vec<PlacementRule>,
}

impl OntologyModification {
    pub fn to_code(&self) -> String {
        format!(
            "// {} [code name: {}] rules: {}\n({}, \"{}\", Modification::Predefined(&[{}], &[{}], \"{}\", \"{}\"))",
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
pub enum PlacementRule {
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
pub enum Position {
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
