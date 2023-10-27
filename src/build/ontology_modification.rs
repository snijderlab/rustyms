use serde::{Deserialize, Serialize};

use crate::{formula::MolecularFormula, glycan::MonoSaccharide};

use super::Modification;

#[derive(Debug, Default)]
pub struct OntologyModification {
    pub diff_formula: MolecularFormula,
    pub monosaccharides: Vec<(MonoSaccharide, i16)>,
    pub code_name: String,
    pub full_name: String,
    pub ontology: Ontology,
    pub id: usize,
    pub rules: Vec<PlacementRule>,
}

impl OntologyModification {
    pub fn into_mod(self) -> (usize, String, Modification) {
        (
            self.id,
            self.code_name.to_ascii_lowercase(),
            Modification::Predefined(
                self.diff_formula,
                self.rules,
                self.ontology,
                self.code_name,
                self.id,
            ),
        )
    }
}

include!("../shared/placement_rule.rs");
include!("../shared/aminoacid.rs");
include!("../shared/ontology.rs");

impl TryInto<Position> for u8 {
    type Error = String;
    fn try_into(self) -> Result<Position, Self::Error> {
        match self {
            b'1' => Ok(Position::Anywhere),
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
            "" => Ok(Position::Anywhere),
            "Anywhere" => Ok(Position::Anywhere),
            "Any N-term" => Ok(Position::AnyNTerm),
            "Any C-term" => Ok(Position::AnyCTerm),
            "Protein N-term" => Ok(Position::ProteinNTerm),
            "Protein C-term" => Ok(Position::ProteinCTerm),
            n => Err(format!("Not valid position: {n}")),
        }
    }
}
