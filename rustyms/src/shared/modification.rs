/// A modification on an amino acid
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum Modification {
    /// Any of the simple modifications
    Simple(SimpleModification),
    /// A cross link to another (or the same) peptide, a branch is also seen as a cross-link but then the name is None.
    CrossLink {
        /// The index of the peptide this cross-link is bound to (can be the index for this peptide if it is an intra link)
        peptide: usize,
        /// The sequence index where this cross-link is bound to
        sequence_index: usize,
        /// The linker that defines the chemical structure that is the actual linker
        linker: SimpleModification,
        /// The name of the cross-linker, if [`CrossLinkName::Branch`] it is a branch instead of cross-link
        name: CrossLinkName,
    },
}

/// A modification on an amino acid
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum SimpleModification {
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
        Vec<(Vec<PlacementRule>, Vec<NeutralLoss>, Vec<DiagnosticIon>)>,
        Ontology, // Context
        String,   // Name
        usize,    // Index
    ),
    /// A modification from the GNOme ontology
    Gno(
        GnoComposition,
        String, // Name
    ),
    /// A cross linker
    Linker {
        specificities: Vec<LinkerSpecificity>,
        formula: MolecularFormula,
        name: String,
        id: usize,
        length: Option<OrderedFloat<f64>>,
        ontology: Ontology,
    },
}

/// The name of a cross-link
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum CrossLinkName {
    /// A branch
    Branch,
    /// A cross-link
    Name(String),
}

/// The linker position specificities for a linker
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum LinkerSpecificity {
    /// A symmetric specificity where both ends have the same specificity.
    /// The first list is all possible positions. The second list is all
    /// stubs that can be left after cleaving or breaking of the cross-link.
    Symmetric(
        Vec<PlacementRule>,
        Vec<(MolecularFormula, MolecularFormula)>,
        Vec<DiagnosticIon>,
    ),
    /// An asymmetric specificity where both ends have a different specificity.
    /// The first list is all possible positions. The second list is all
    /// stubs that can be left after cleaving or breaking of the cross-link.
    Asymmetric(
        (Vec<PlacementRule>, Vec<PlacementRule>),
        Vec<(MolecularFormula, MolecularFormula)>,
        Vec<DiagnosticIon>,
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
