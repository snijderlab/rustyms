use serde::{Deserialize, Serialize};

use crate::{formula::MolecularFormula, glycan::MonoSaccharide, NeutralLoss};

use super::Modification;

#[derive(Debug, Default)]
pub struct OntologyModification {
    pub diff_formula: MolecularFormula,
    pub monosaccharides: Vec<(MonoSaccharide, i16)>,
    pub code_name: String,
    pub full_name: String,
    pub ontology: Ontology,
    pub id: usize,
    pub rules: Vec<(PlacementRule, Vec<NeutralLoss>)>,
}

impl OntologyModification {
    /// Simplify the placement rules
    pub fn simplify_rules(&mut self) {
        let mut new = Vec::new();
        for rule in &self.rules {
            if new.is_empty() {
                new.push(rule.clone());
            } else {
                let mut found = false;
                for new_rule in &mut new {
                    if let (
                        (PlacementRule::AminoAcid(aa, pos), losses),
                        (PlacementRule::AminoAcid(new_aa, new_pos), new_losses),
                    ) = (rule, new_rule)
                    {
                        if *pos == *new_pos && *losses == *new_losses {
                            new_aa.extend(aa);
                            new_aa.sort_unstable();
                            found = true;
                            break;
                        }
                    }
                }
                if !found {
                    new.push(rule.clone());
                }
            }
        }
        self.rules = new;
    }

    pub fn into_mod(mut self) -> (usize, String, Modification) {
        self.simplify_rules();
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
