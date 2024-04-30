//! Rules regarding the placement of modifications

use serde::{Deserialize, Serialize};

use crate::{
    fragment::PeptidePosition,
    modification::{Modification, Ontology, SimpleModification},
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
                    if let Modification::Simple(SimpleModification::Predefined(
                        _,
                        _,
                        Ontology::Psimod,
                        _,
                        i,
                    )) = m
                    {
                        i == mod_index
                    } else {
                        false
                    }
                }) && r_pos.is_possible(position)
            }
            Self::Terminal(r_pos) => {
                r_pos.is_possible(position)
                    && (position.is_n_terminal() || position.is_c_terminal())
            }
        }
    }

    /// Check if any of the given rules are possible
    pub fn any_possible(rules: &[Self], seq: &SequenceElement, position: &PeptidePosition) -> bool {
        rules.iter().any(|r| r.is_possible(seq, position))
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

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {

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
                    modifications: vec![Ontology::Psimod
                        .find_modification_id(30, None)
                        .unwrap()
                        .into()],
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
        assert!(
            dbg!(Ontology::Unimod.find_modification_id(7, None).unwrap()).is_possible(
                &SequenceElement {
                    aminoacid: AminoAcid::Q,
                    modifications: Vec::new(),
                    possible_modifications: Vec::new(),
                    ambiguous: None
                },
                &PeptidePosition::n(4, 5)
            ),
            "unimod deamidated at end"
        );
    }
}
