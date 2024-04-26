/// A modification on an amino acid
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum Modification {
    /// Any of the simple modifications
    Simple(SimpleModification),
    /// An intra link, cross linking this position to another position on the same peptide. Has to be present at both locations.
    IntraLink {
        index: usize,
        linker: Linker,
        name: String,
    },
    /// A cross link to another peptide
    CrossLink {
        peptide: usize,
        index: usize,
        linker: Linker,
        name: String,
    },
    /// A branch to another peptide, needed? Or could this be expressed as a special case of a cross link
    Branch {
        peptide: usize,
        index: usize, // Needed?
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
}

/// A database entry for a cross-linker.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct Linker {
    specificities: LinkerSpecificity,
    formula: MolecularFormula, // TODO: what if not defined?, also make sure anything that only has a mass uses that instead
    name: String,
    /// XL:000XXX
    id: usize,
    length: Option<OrderedFloat<f64>>,
    ontology: Ontology,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum LinkerSpecificity {
    Symmetric(Vec<PlacementRule>),
    Asymmetric(Vec<PlacementRule>, Vec<PlacementRule>),
}

/// All possible compositions in the GNO ontology
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum GnoComposition {
    /// Only the mass is known
    Mass(OrderedMass),
    /// The (full) structure is known
    Structure(GlycanStructure),
}
