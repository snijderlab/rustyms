//! Rules regarding the placement of modifications

use std::str::FromStr;

use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    fragment::PeptidePosition,
    modification::{Modification, ModificationId, Ontology, SimpleModification},
    AminoAcid, SequenceElement,
};

include!("shared/placement_rule.rs");

impl PlacementRule {
    /// Check if this rule fits with the given location
    pub fn is_possible(&self, seq: &SequenceElement, position: &PeptidePosition) -> bool {
        match self {
            Self::AminoAcid(aa, r_pos) => {
                aa.iter().any(|a| *a == seq.aminoacid) && r_pos.is_possible(position)
            }
            Self::PsiModification(mod_index, r_pos) => {
                seq.modifications.iter().any(|m| {
                    if let Modification::Simple(SimpleModification::Database {
                        id:
                            ModificationId {
                                ontology: Ontology::Psimod,
                                id,
                                ..
                            },
                        ..
                    }) = m
                    {
                        id == mod_index
                    } else {
                        false
                    }
                }) && r_pos.is_possible(position)
            }
            Self::Terminal(r_pos) => {
                r_pos.is_possible(position)
                    && (position.is_n_terminal() || position.is_c_terminal())
            }
            Self::Anywhere => true,
        }
    }

    /// Check if any of the given rules are possible
    pub fn any_possible(rules: &[Self], seq: &SequenceElement, position: &PeptidePosition) -> bool {
        rules.iter().any(|r| r.is_possible(seq, position))
    }
}

impl FromStr for PlacementRule {
    type Err = CustomError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((head, tail)) = s.split_once('@') {
            let aa: Vec<AminoAcid> = head
                .chars()
                .enumerate()
                .map(|(i, c)| {
                    AminoAcid::try_from(c).map_err(|()| {
                        CustomError::error(
                            "Invalid amino acid",
                            "Invalid amino acid in specified amino acids in placement rule",
                            Context::line(None, s, i, 1),
                        )
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            tail.parse().map_or_else(
                |()| {
                    Err(CustomError::error(
                        "Invalid position",
                        "Use any of the following for the position: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm",
                        Context::line(None, s, head.len() + 1, tail.len()),
                    ))
                },
                |position| Ok(Self::AminoAcid(aa, position)),
            )
        } else if let Ok(position) = s.parse() {
            Ok(match position {
                Position::Anywhere => Self::Anywhere,
                pos => Self::Terminal(pos),
            })
        } else {
            Err(CustomError::error(
                "Invalid position",
                "Use any of the following for the position: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm",
                Context::full_line(0, s),
            ))
        }
    }
}

impl Position {
    /// See if the given peptide position is a valid position given this [`Position`] as placement rule.
    pub const fn is_possible(self, position: &PeptidePosition) -> bool {
        match self {
            Self::Anywhere => true,
            Self::AnyNTerm | Self::ProteinNTerm => position.is_n_terminal(),
            Self::AnyCTerm | Self::ProteinCTerm => position.is_c_terminal(),
        }
    }
}

impl FromStr for Position {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "anywhere" => Ok(Self::Anywhere),
            "anynterm" => Ok(Self::AnyNTerm),
            "proteinnterm" => Ok(Self::ProteinNTerm),
            "anycterm" => Ok(Self::AnyCTerm),
            "proteincterm" => Ok(Self::ProteinCTerm),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {

    use crate::modification::RulePossible;

    use super::*;
    #[test]
    fn multi_level_rule() {
        assert!(
            !PlacementRule::PsiModification(30, Position::Anywhere).is_possible(
                &SequenceElement {
                    aminoacid: AminoAcid::Alanine,
                    modifications: Vec::new(),
                    possible_modifications: Vec::new(),
                    ambiguous: None
                },
                &PeptidePosition::n(0, 1)
            ),
            "Multi level mod cannot be placed if the dependent mod is not present"
        );
        assert!(
            PlacementRule::PsiModification(30, Position::Anywhere).is_possible(
                &SequenceElement {
                    aminoacid: AminoAcid::Alanine,
                    modifications: vec![Ontology::Psimod.find_id(30, None).unwrap().into()],
                    possible_modifications: Vec::new(),
                    ambiguous: None
                },
                &PeptidePosition::n(0, 1)
            ),
            "Multi level mod can be placed if the dependent mod is present"
        );
    }

    #[test]
    fn place_anywhere() {
        assert!(
            PlacementRule::AminoAcid(vec![AminoAcid::Q], Position::Anywhere).is_possible(
                &SequenceElement {
                    aminoacid: AminoAcid::Q,
                    modifications: Vec::new(),
                    possible_modifications: Vec::new(),
                    ambiguous: None
                },
                &PeptidePosition::n(0, 5)
            ),
            "start"
        );
        assert!(
            PlacementRule::AminoAcid(vec![AminoAcid::Q], Position::Anywhere).is_possible(
                &SequenceElement {
                    aminoacid: AminoAcid::Q,
                    modifications: Vec::new(),
                    possible_modifications: Vec::new(),
                    ambiguous: None
                },
                &PeptidePosition::n(2, 5)
            ),
            "middle"
        );
        assert!(
            PlacementRule::AminoAcid(vec![AminoAcid::Q], Position::Anywhere).is_possible(
                &SequenceElement {
                    aminoacid: AminoAcid::Q,
                    modifications: Vec::new(),
                    possible_modifications: Vec::new(),
                    ambiguous: None
                },
                &PeptidePosition::n(4, 5)
            ),
            "end"
        );
        assert_eq!(
            dbg!(Ontology::Unimod.find_id(7, None).unwrap()).is_possible(
                &SequenceElement {
                    aminoacid: AminoAcid::Q,
                    modifications: Vec::new(),
                    possible_modifications: Vec::new(),
                    ambiguous: None
                },
                &PeptidePosition::n(4, 5)
            ),
            RulePossible::Symmetric,
            "unimod deamidated at end"
        );
    }
}
