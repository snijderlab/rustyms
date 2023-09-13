use crate::{
    modification::{Modification, Ontology},
    AminoAcid, SequenceElement,
};

#[derive(Debug, PartialEq, Eq)]
pub enum PlacementRule {
    AminoAcid(AminoAcid, Position),
    PsiModification(usize, Position),
    Terminal(Position),
}

impl PlacementRule {
    pub fn is_possible(&self, seq: &SequenceElement, index: usize, length: usize) -> bool {
        match self {
            Self::AminoAcid(aa, r_pos) => *aa == seq.aminoacid && r_pos.is_possible(index, length),
            Self::PsiModification(mod_index, r_pos) => {
                seq.modifications.iter().any(|m| {
                    if let Modification::Predefined(_, _, Ontology::PSIMOD, _, i) = m {
                        i == mod_index
                    } else {
                        false
                    }
                }) && r_pos.is_possible(index, length)
            }
            Self::Terminal(r_pos) => {
                r_pos.is_possible(index, length) && (index == length - 1 || index == 0)
            }
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum Position {
    Anywhere,
    AnyNTerm,
    AnyCTerm,
    ProteinNTerm,
    ProteinCTerm,
}

impl Position {
    const fn is_possible(&self, index: usize, length: usize) -> bool {
        match self {
            Self::Anywhere => true,
            Self::AnyNTerm | Self::ProteinNTerm => index == 0,
            Self::AnyCTerm | Self::ProteinCTerm => index == length - 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{modification::find_id_in_ontology, ontologies::PSI_MOD_ONTOLOGY};

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
                0,
                1
            ),
            "Multi level mod cannot be placed if the dependent mod is not present"
        );
        assert!(
            PlacementRule::PsiModification(30, Position::Anywhere).is_possible(
                &SequenceElement {
                    aminoacid: AminoAcid::Alanine,
                    modifications: vec![find_id_in_ontology(30, PSI_MOD_ONTOLOGY).unwrap()],
                    possible_modifications: Vec::new(),
                    ambiguous: None
                },
                0,
                1
            ),
            "Multi level mod can be placed if the dependent mod is present"
        );
    }
}
