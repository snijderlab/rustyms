use serde::{Deserialize, Serialize};

use crate::{fragment::PeptidePosition, AminoAcid, MolecularFormula};

/// Keep track of what ambiguous option is used
/// TODO: Maybe also use the labelling system to keep track of which modification neutral loss is applied?
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
enum AmbiguousLocation {
    /// A ambiguous amino acid, with the actual amino acid used tracked
    AminoAcid {
        /// Which amino acid is used
        option: AminoAcid,
        /// What location in the sequence are we talking about
        location: PeptidePosition,
    },
    /// A ambiguous modification, with the actual position
    Modification {
        /// Which ambiguous modification
        id: usize,
        /// Which location
        location: PeptidePosition,
    },
    /// The actual charge used, when there are multiple charge carriers
    ChargeCarrier(MolecularFormula),
}
