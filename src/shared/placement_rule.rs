#[derive(Debug, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub enum PlacementRule {
    AminoAcid(Vec<AminoAcid>, Position),
    PsiModification(usize, Position),
    Terminal(Position),
}

#[derive(Debug, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub enum Position {
    Anywhere,
    AnyNTerm,
    AnyCTerm,
    ProteinNTerm,
    ProteinCTerm,
}
