/// A rule determining the placement of a modification
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum PlacementRule {
    /// Placed on an aminoacid on the given position
    AminoAcid(Vec<AminoAcid>, Position),
    /// Placed on an another modification on the given position
    PsiModification(usize, Position),
    /// Placed on a terminal position
    Terminal(Position),
    /// Just anywhere
    Anywhere,
}

/// A position where a modification can be placed
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum Position {
    /// At any location
    #[default]
    Anywhere,
    /// At the N term of a peptide or protein
    AnyNTerm,
    /// At the C term of a peptide or protein
    AnyCTerm,
    /// At the N term of a protein
    ProteinNTerm,
    /// At the C term of a protein
    ProteinCTerm,
}

impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Anywhere => "Anywhere",
                Self::AnyNTerm => "AnyNTerm",
                Self::AnyCTerm => "AnyCTerm",
                Self::ProteinNTerm => "ProteinNTerm",
                Self::ProteinCTerm => "ProteinCTerm",
            },
        )
    }
}
