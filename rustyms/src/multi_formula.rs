use crate::{MolecularFormula, Multi};

/// Any item that has a number of potential chemical formulas
pub trait MultiChemical {
    /// Get all possible molecular formulas
    fn formulas(&self) -> Multi<MolecularFormula>;

    /// Get the charge of this chemical, it returns None if no charge is defined.
    fn charge(&self) -> Option<i16> {
        self.formulas()
            .first()
            .map(MolecularFormula::charge)
            .filter(|c| *c != 0)
    }
}

// #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Hash)]
struct AmbiguousLabels<Position, Option>(Vec<AmbiguousLabel<Position, Option>>);

struct AmbiguousLabel<Position, Option> {
    position: Position,
    option: Option,
}

impl<Position, Option> AmbiguousLabel<Position, Option> {
    pub const fn position(&self) -> &Position {
        &self.position
    }
    pub const fn option(&self) -> &Option {
        &self.option
    }
}

impl<Position: Clone, Option: Clone> Clone for AmbiguousLabel<Position, Option> {
    fn clone(&self) -> Self {
        Self {
            position: self.position.clone(),
            option: self.option.clone(),
        }
    }
}
