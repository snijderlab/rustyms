use serde::{Deserialize, Serialize};

use crate::{glycan::GlycanStructure, *};

/// A modification on an amino acid
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Modification {
    /// A modification defined with a monoisotopic mass shift
    Mass(Mass),
    /// A modification defined with a molecular formula
    #[allow(non_snake_case)]
    Formula(MolecularFormula),
    /// A glycan without a defined structure
    Glycan(Vec<(MonoSaccharide, isize)>),
    /// A glycan with a defined structure
    GlycanStructure(GlycanStructure),
    /// A modification from one of the modification ontologies
    Predefined(
        &'static [(Element, u16, i16)],
        &'static [PlacementRule],
        Ontology,     // Context
        &'static str, // Name
        usize,        // Index
    ),
    /// A modification from one of the modification ontologies
    Gno(
        GnoComposition,
        &'static str, // Name
    ),
}

/// All possible compositions in the GNO ontology
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum GnoComposition {
    /// Only the mass is known
    Mass(f64),
    /// The (full) structure is known
    Structure(GlycanStructure),
}
