use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::{
    build::glycan::MonoSaccharide,
    formula::{Chemical, MolecularFormula},
    DiagnosticIon, LinkerSpecificity, NeutralLoss, SimpleModification,
};

#[derive(Debug, Default)]
pub struct OntologyModification {
    pub formula: MolecularFormula,
    pub code_name: String,
    pub full_name: String,
    pub ontology: Ontology,
    pub id: usize,
    pub data: ModData,
}

#[derive(Debug)]
pub enum ModData {
    Mod {
        rules: Vec<(Vec<PlacementRule>, Vec<NeutralLoss>, Vec<DiagnosticIon>)>,
        monosaccharides: Vec<(MonoSaccharide, i32)>,
    },
    Linker {
        diagnostic_ions: Vec<DiagnosticIon>,
        length: Option<OrderedFloat<f64>>,
        specificities: LinkerSpecificity,
    },
}

impl Default for ModData {
    fn default() -> Self {
        Self::Mod {
            rules: Vec::new(),
            monosaccharides: Vec::new(),
        }
    }
}

impl OntologyModification {
    /// Simplify the placement rules
    pub fn simplify_rules(&mut self) {
        if let ModData::Mod { ref mut rules, .. } = self.data {
            let mut new = Vec::new();
            for rule in rules.iter() {
                if new.is_empty() {
                    new.push(rule.clone());
                } else {
                    let mut found = false;
                    for new_rule in &mut new {
                        // Check if there is a rule with the same neutral loss and diagnostic ions (these can be location specific)
                        if new_rule.1 == rule.1 && new_rule.2 == rule.2 {
                            found = true;
                            // Check if there are other rules in this set of neutral&diagnostic that also use AA placements
                            // If there are, and they are on the same position, merge the AA set
                            for position in &rule.0 {
                                let mut pos_found = false;
                                for new_position in &mut new_rule.0 {
                                    if let (
                                        PlacementRule::AminoAcid(new_aa, new_pos),
                                        PlacementRule::AminoAcid(aa, pos),
                                    ) = (new_position, position)
                                    {
                                        if *new_pos == *pos {
                                            new_aa.extend(aa);
                                            new_aa.sort_unstable();
                                            pos_found = true;
                                            break;
                                        }
                                    }
                                }
                                if !pos_found {
                                    new_rule.0.push(position.clone());
                                }
                            }
                        }
                    }
                    if !found {
                        new.push(rule.clone());
                    }
                }
            }
            rules.clear();
            rules.extend(new);
        }
    }

    pub fn into_mod(mut self) -> (usize, String, SimpleModification) {
        self.simplify_rules();
        match self.data {
            ModData::Mod {
                monosaccharides,
                rules,
            } => (
                self.id,
                self.code_name.to_ascii_lowercase(),
                SimpleModification::Predefined(
                    self.formula
                        + monosaccharides
                            .iter()
                            .map(|(s, n)| s.formula() * n)
                            .sum::<MolecularFormula>(),
                    rules,
                    self.ontology,
                    self.code_name,
                    self.id,
                ),
            ),
            ModData::Linker {
                specificities,
                length,
                diagnostic_ions,
            } => (
                self.id,
                self.code_name.to_ascii_lowercase(),
                SimpleModification::Linker {
                    specificities,
                    formula: self.formula,
                    name: self.code_name,
                    id: self.id,
                    length,
                    ontology: self.ontology,
                    diagnostic_ions,
                },
            ),
        }
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
