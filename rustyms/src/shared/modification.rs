/// A modification on an amino acid
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub enum Modification {
    /// Any of the simple modifications
    Simple(SimpleModification),
    /// A cross link to another (or the same) peptide, a branch is also seen as a cross-link but then the name is None.
    CrossLink {
        /// The index of the peptide this cross-link is bound to (can be the index for this peptide if it is an intra link)
        peptide: usize,
        /// The sequence index where this cross-link is bound to
        sequence_index: crate::SequencePosition,
        /// The linker that defines the chemical structure that is the actual linker
        linker: SimpleModification,
        /// The name of the cross-linker, if [`CrossLinkName::Branch`] it is a branch instead of cross-link
        name: CrossLinkName,
        /// To determine if the cross-link is placed symmetrically or if asymmetrically if this is the left or right side
        side: CrossLinkSide,
    },
    /// An ambiguous modification, that can be placed at multiple locations
    Ambiguous {
        /// The name of the group
        group: String,
        /// The id to compare be able to find the other locations where this modifications can be placed
        id: usize,
        /// The modification itself
        modification: SimpleModification,
        /// If present the localisation score, meaning the chance/ratio for this modification to show up on this exact spot
        localisation_score: Option<OrderedFloat<f64>>,
        /// If this is the preferred location or not
        preferred: bool,
    },
}

/// Indicate the cross-link side
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum CrossLinkSide {
    /// The cross-link is symmetric, or if asymmetric it can be placed in both orientations
    Symmetric(std::collections::BTreeSet<usize>),
    /// The cross-link is asymmetric and this is the 'left' side
    Left(std::collections::BTreeSet<usize>),
    /// The cross-link is asymmetric and this is the 'right' side
    Right(std::collections::BTreeSet<usize>),
}

impl PartialOrd for CrossLinkSide {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CrossLinkSide {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Symmetric(_), Self::Symmetric(_)) | (Self::Left(_), Self::Left(_)) => {
                Ordering::Equal
            }
            (Self::Symmetric(_), _) => Ordering::Greater,
            (_, Self::Symmetric(_)) => Ordering::Less,
            (Self::Left(_), _) => Ordering::Greater,
            (_, Self::Left(_)) => Ordering::Less,
            (Self::Right(_), Self::Right(_)) => Ordering::Equal,
        }
    }
}

impl std::hash::Hash for CrossLinkSide {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        use itertools::Itertools;
        let (i, r) = match self {
            Self::Symmetric(r) => (0, r),
            Self::Left(r) => (1, r),
            Self::Right(r) => (2, r),
        };
        state.write_u8(i);
        state.write(
            &r.iter()
                .sorted()
                .flat_map(|r| r.to_ne_bytes())
                .collect_vec(),
        );
    }
}

/// A modification on an amino acid, wrapped in an [`std::sync::Arc`] to not have to clone modifications from databases.
pub type SimpleModification = std::sync::Arc<SimpleModificationInner>;

/// A modification on an amino acid
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum SimpleModificationInner {
    /// A modification defined with a monoisotopic mass shift
    Mass(OrderedMass),
    /// A modification defined with a molecular formula
    #[allow(non_snake_case)]
    Formula(MolecularFormula),
    /// A glycan without a defined structure
    Glycan(Vec<(MonoSaccharide, isize)>),
    /// A glycan with a defined structure
    GlycanStructure(GlycanStructure),
    /// A modification from the GNOme ontology
    Gno {
        /// The composition, weight/composition/topology
        composition: GnoComposition,
        /// The id/name
        id: ModificationId,
        /// The structure score
        structure_score: Option<usize>,
        /// The subsumption level
        subsumption_level: GnoSubsumption,
        /// The underlying glycan motif, first is the human description, the second id the GNOme ID of the motif
        motif: Option<(String, String)>,
        /// Taxonomy of the animals in which this glycan is found, defined as a list of species name with taxonomy ID
        taxonomy: thin_vec::ThinVec<(String, usize)>,
        /// Locations of where the glycan exists
        glycomeatlas: thin_vec::ThinVec<(String, Vec<(String, String)>)>,
    },
    /// A modification from one of the modification ontologies
    Database {
        /// The placement rules, neutral losses, and diagnostic ions
        specificities: Vec<(Vec<PlacementRule>, Vec<NeutralLoss>, Vec<DiagnosticIon>)>,
        /// The chemical formula for this modification (diff formula)
        formula: MolecularFormula,
        /// The id/name
        id: ModificationId,
    },
    /// A cross-linker
    Linker {
        /// All possible specificities for this linker
        specificities: Vec<LinkerSpecificity>,
        /// The chemical formula for this linker (diff formula)
        formula: MolecularFormula,
        /// The id/name
        id: ModificationId,
        /// The length, if known
        length: Option<OrderedFloat<f64>>,
    },
}

/// A modification id/name
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash, Default)]
pub struct ModificationId {
    /// The ontology where this linker is defined
    pub ontology: Ontology,
    /// The name
    pub name: String,
    /// The id
    pub id: Option<usize>,
    /// The description, mostly for search results
    pub description: String,
    /// Any synonyms
    pub synonyms: thin_vec::ThinVec<String>,
    /// Cross reference IDs
    pub cross_ids: thin_vec::ThinVec<(String, String)>,
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
    Weight(OrderedMass),
    /// The composition,
    Composition(Vec<(MonoSaccharide, isize)>),
    /// The (full) structure is known
    Topology(GlycanStructure),
}

/// All possible subsumption levels in the GNOme database indicating different levels of description for a glycan species
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash, Serialize, Deserialize,
)]
pub enum GnoSubsumption {
    /// Indicates only the average weight is defined
    #[default]
    AverageWeight,
    /// Indicates the basic composition, without isomeric information
    BaseComposition,
    /// Indicates the composition, with isomeric information
    Composition,
    /// Indicates the topology, without linkage and anomeric information
    Topology,
    /// Indicates the topology, without reducing end ring and anomeric information
    Saccharide,
}

impl std::fmt::Display for GnoSubsumption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AverageWeight => write!(f, "Average weight"),
            Self::BaseComposition => write!(f, "Base composition (no isomeric information)"),
            Self::Composition => write!(f, "Composition"),
            Self::Topology => write!(f, "Topology (no linkage)"),
            Self::Saccharide => write!(f, "Saccharide"),
        }
    }
}
