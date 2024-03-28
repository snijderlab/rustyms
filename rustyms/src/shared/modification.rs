/// A modification on an amino acid
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum Modification {
    /// A modification defined with a monoisotopic mass shift
    Mass(OrderedMass),
    /// A modification defined with a molecular formula
    #[allow(non_snake_case)]
    Formula(MolecularFormula),
    /// A glycan without a defined structure
    Glycan(Vec<(MonoSaccharide, isize)>),
    /// A glycan with a defined structure
    GlycanStructure(GlycanStructure),
    /// A modification from one of the modification ontologies
    Predefined(
        MolecularFormula,
        Vec<(PlacementRule, Vec<NeutralLoss>)>,
        Ontology, // Context
        String,   // Name
        usize,    // Index
    ),
    /// A modification from one of the modification ontologies
    Gno(
        GnoComposition,
        String, // Name
    ),
}

/// All possible compositions in the GNO ontology
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum GnoComposition {
    /// Only the mass is known
    Mass(OrderedMass),
    /// The (full) structure is known
    Structure(GlycanStructure),
}
