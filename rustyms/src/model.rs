//! Handle model instantiation.

use std::ops::RangeInclusive;

use serde::{Deserialize, Serialize};

use crate::{
    fragment::PeptidePosition,
    helper_functions::RangeExtension,
    system::{f64::MassOverCharge, mz},
    NeutralLoss, Tolerance,
};

/// A model for the fragmentation, allowing control over what theoretical fragments to generate.
#[non_exhaustive]
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
#[allow(clippy::struct_excessive_bools)]
pub struct Model {
    /// a series ions
    pub a: (Location, Vec<NeutralLoss>),
    /// b series ions
    pub b: (Location, Vec<NeutralLoss>),
    /// c series ions
    pub c: (Location, Vec<NeutralLoss>),
    /// d series ions (side chain fragmentation from a)
    pub d: (Location, Vec<NeutralLoss>),
    /// v series ions (full side chain broken off)
    pub v: (Location, Vec<NeutralLoss>),
    /// w series ions (side chain fragmentation from z)
    pub w: (Location, Vec<NeutralLoss>),
    /// x series ions
    pub x: (Location, Vec<NeutralLoss>),
    /// y series ions
    pub y: (Location, Vec<NeutralLoss>),
    /// z series ions
    pub z: (Location, Vec<NeutralLoss>),
    /// precursor ions
    pub precursor: Vec<NeutralLoss>,
    /// immonium ions
    pub immonium: bool,
    /// m ions, loss of the amino acid side chain from the precursor
    pub m: bool,
    /// If the neutral losses specific for modifications should be generated
    pub modification_specific_neutral_losses: bool,
    /// If the diagnostic ions specific for modifications should be generated
    pub modification_specific_diagnostic_ions: bool,
    /// (allow structural fragments, allow compositional fragments (of this number of mono saccharides (inclusive range)), the allowed neutral losses)
    pub glycan: (bool, (usize, usize), Vec<NeutralLoss>),
    /// Allow any MS cleavable cross-link to be cleaved
    pub allow_cross_link_cleavage: bool,
    /// The matching tolerance
    pub tolerance: Tolerance<MassOverCharge>,
    /// The range in which fragments fall, can be used to limit the theoretical fragments to a known window
    pub mz_range: RangeInclusive<MassOverCharge>,
}

/// A struct to handle all possible fragments that could be generated on a single location
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
#[non_exhaustive]
pub struct PossibleIons<'a> {
    /// a series ions
    pub a: (bool, &'a [NeutralLoss]),
    /// b series ions
    pub b: (bool, &'a [NeutralLoss]),
    /// c series ions
    pub c: (bool, &'a [NeutralLoss]),
    /// d series ions (side chain fragmentation from a)
    pub d: (bool, &'a [NeutralLoss]),
    /// v series ions (full side chain broken off)
    pub v: (bool, &'a [NeutralLoss]),
    /// w series ions (side chain fragmentation from z)
    pub w: (bool, &'a [NeutralLoss]),
    /// x series ions
    pub x: (bool, &'a [NeutralLoss]),
    /// y series ions
    pub y: (bool, &'a [NeutralLoss]),
    /// z series ions
    pub z: (bool, &'a [NeutralLoss]),
    /// precursor ions
    pub precursor: &'a [NeutralLoss],
    /// immonium
    pub immonium: bool,
}

impl<'a> PossibleIons<'a> {
    /// Give an upper bound for the number of theoretical fragment for these possible ions
    pub fn size_upper_bound(&self) -> usize {
        usize::from(self.a.0) * (self.a.1.len() + 1)
            + usize::from(self.b.0) * (self.b.1.len() + 1)
            + usize::from(self.c.0) * (self.c.1.len() + 1)
            + usize::from(self.d.0) * 2 * (self.d.1.len() + 1)
            + usize::from(self.v.0) * (self.v.1.len() + 1)
            + usize::from(self.w.0) * 2 * (self.w.1.len() + 1)
            + usize::from(self.x.0) * (self.x.1.len() + 1)
            + usize::from(self.y.0) * (self.y.1.len() + 1)
            + usize::from(self.z.0) * 2 * (self.z.1.len() + 1)
            + self.precursor.len()
            + 1
    }
}

/// Builder style methods
impl Model {
    /// Set a to the given location and overwrite the neutral losses
    #[must_use]
    pub fn a(self, location: Location, neutral_loss: Vec<NeutralLoss>) -> Self {
        Self {
            a: (location, neutral_loss),
            ..self
        }
    }
    /// Set b to the given location and overwrite the neutral losses
    #[must_use]
    pub fn b(self, location: Location, neutral_loss: Vec<NeutralLoss>) -> Self {
        Self {
            b: (location, neutral_loss),
            ..self
        }
    }
    /// Set c to the given location and overwrite the neutral losses
    #[must_use]
    pub fn c(self, location: Location, neutral_loss: Vec<NeutralLoss>) -> Self {
        Self {
            c: (location, neutral_loss),
            ..self
        }
    }
    /// Set d to the given location and overwrite the neutral losses
    #[must_use]
    pub fn d(self, location: Location, neutral_loss: Vec<NeutralLoss>) -> Self {
        Self {
            d: (location, neutral_loss),
            ..self
        }
    }
    /// Set v to the given location and overwrite the neutral losses
    #[must_use]
    pub fn v(self, location: Location, neutral_loss: Vec<NeutralLoss>) -> Self {
        Self {
            v: (location, neutral_loss),
            ..self
        }
    }
    /// Set w to the given location and overwrite the neutral losses
    #[must_use]
    pub fn w(self, location: Location, neutral_loss: Vec<NeutralLoss>) -> Self {
        Self {
            w: (location, neutral_loss),
            ..self
        }
    }
    /// Set x to the given location and overwrite the neutral losses
    #[must_use]
    pub fn x(self, location: Location, neutral_loss: Vec<NeutralLoss>) -> Self {
        Self {
            x: (location, neutral_loss),
            ..self
        }
    }
    /// Set y to the given location and overwrite the neutral losses
    #[must_use]
    pub fn y(self, location: Location, neutral_loss: Vec<NeutralLoss>) -> Self {
        Self {
            y: (location, neutral_loss),
            ..self
        }
    }
    /// Set z to the given location and overwrite the neutral losses
    #[must_use]
    pub fn z(self, location: Location, neutral_loss: Vec<NeutralLoss>) -> Self {
        Self {
            z: (location, neutral_loss),
            ..self
        }
    }
    /// Overwrite the precursor neutral losses
    #[must_use]
    pub fn precursor(self, neutral_loss: Vec<NeutralLoss>) -> Self {
        Self {
            precursor: neutral_loss,
            ..self
        }
    }
    /// Set immonium
    #[must_use]
    pub fn immonium(self, state: bool) -> Self {
        Self {
            immonium: state,
            ..self
        }
    }
    /// Set m
    #[must_use]
    pub fn m(self, state: bool) -> Self {
        Self { m: state, ..self }
    }
    /// Set modification specific neutral losses
    #[must_use]
    pub fn modification_specific_neutral_losses(self, state: bool) -> Self {
        Self {
            modification_specific_neutral_losses: state,
            ..self
        }
    }
    /// Set modification specific diagnostic ions
    #[must_use]
    pub fn modification_specific_diagnostic_ions(self, state: bool) -> Self {
        Self {
            modification_specific_diagnostic_ions: state,
            ..self
        }
    }
    /// Set glycans, `None` makes no fragments, `Some(_)` makes fragments, with the given neutral losses
    #[must_use]
    pub fn glycan(
        self,
        allow_structural: bool,
        allow_compositional: std::ops::RangeInclusive<usize>,
        neutral_losses: Vec<NeutralLoss>,
    ) -> Self {
        Self {
            glycan: (
                allow_structural,
                (
                    allow_compositional.start_index(),
                    allow_compositional.end_index(usize::MAX),
                ),
                neutral_losses,
            ),
            ..self
        }
    }
    /// Set the tolerance
    #[must_use]
    pub fn allow_cross_link_cleavage(self, state: bool) -> Self {
        Self {
            allow_cross_link_cleavage: state,
            ..self
        }
    }
    /// Set the tolerance
    #[must_use]
    pub fn tolerance(self, tolerance: impl Into<Tolerance<MassOverCharge>>) -> Self {
        Self {
            tolerance: tolerance.into(),
            ..self
        }
    }
    /// Set the mz range
    #[must_use]
    pub fn mz_range(self, mz_range: RangeInclusive<MassOverCharge>) -> Self {
        Self { mz_range, ..self }
    }
}

impl Model {
    /// Give all possible ions for the given N position
    pub fn ions(&self, position: PeptidePosition) -> PossibleIons {
        let c_position = position.flip_terminal();
        PossibleIons {
            a: (self.a.0.possible(position), self.a.1.as_slice()),
            b: (self.b.0.possible(position), self.b.1.as_slice()),
            c: (self.c.0.possible(position), self.c.1.as_slice()),
            d: (self.d.0.possible(position), self.d.1.as_slice()),
            v: (self.v.0.possible(c_position), self.v.1.as_slice()),
            w: (self.w.0.possible(c_position), self.w.1.as_slice()),
            x: (self.x.0.possible(c_position), self.x.1.as_slice()),
            y: (self.y.0.possible(c_position), self.y.1.as_slice()),
            z: (self.z.0.possible(c_position), self.z.1.as_slice()),
            precursor: self.precursor.as_slice(),
            immonium: self.immonium,
        }
    }

    /// Generate all possible fragments
    pub fn all() -> Self {
        Self {
            a: (
                Location::All,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            b: (
                Location::All,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            c: (
                Location::All,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            d: (
                Location::All,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            v: (
                Location::All,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            w: (
                Location::All,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            x: (
                Location::All,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            y: (
                Location::All,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            z: (
                Location::All,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            precursor: vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            immonium: true,
            m: true,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: true,
            glycan: (
                true,
                (1, 3),
                vec![
                    NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 4 O 2)),
                ],
            ),
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// Generate no fragments (except for precursor)
    pub fn none() -> Self {
        Self {
            a: (Location::None, vec![]),
            b: (Location::None, vec![]),
            c: (Location::None, vec![]),
            d: (Location::None, vec![]),
            v: (Location::None, vec![]),
            w: (Location::None, vec![]),
            x: (Location::None, vec![]),
            y: (Location::None, vec![]),
            z: (Location::None, vec![]),
            precursor: vec![],
            immonium: false,
            m: false,
            modification_specific_neutral_losses: false,
            modification_specific_diagnostic_ions: false,
            glycan: (false, (0, 0), Vec::new()),
            allow_cross_link_cleavage: false,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// electron-transfer/higher-energy collisional dissociation
    pub fn ethcd() -> Self {
        Self {
            a: (Location::None, Vec::new()),
            b: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            c: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            d: (Location::None, Vec::new()),
            v: (Location::None, Vec::new()),
            w: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            x: (Location::None, Vec::new()),
            y: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            z: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            precursor: vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            immonium: false,
            m: false,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: true,
            glycan: (
                true,
                (1, 3),
                vec![
                    NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 4 O 2)),
                ],
            ),
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// CID Hcd
    pub fn cid_hcd() -> Self {
        Self {
            a: (
                Location::TakeN { skip: 1, take: 1 },
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            b: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            c: (Location::None, Vec::new()),
            d: (
                Location::TakeN { skip: 1, take: 1 },
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            v: (Location::None, Vec::new()),
            w: (Location::None, Vec::new()),
            x: (Location::None, Vec::new()),
            y: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            z: (Location::None, Vec::new()),
            precursor: vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            immonium: false,
            m: false,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: true,
            glycan: (false, (0, 0), Vec::new()),
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// ETD
    pub fn etd() -> Self {
        Self {
            a: (Location::None, Vec::new()),
            b: (Location::None, Vec::new()),
            c: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            d: (Location::None, Vec::new()),
            v: (Location::None, Vec::new()),
            w: (Location::None, Vec::new()), // TODO: Are w ions also formed here?
            x: (Location::None, Vec::new()),
            y: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            z: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            precursor: vec![
                NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
                NeutralLoss::Loss(molecular_formula!(H 1 O 1)),
                NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 1 O 2)),
                NeutralLoss::Loss(molecular_formula!(C 2 H 3 O 2)),
            ],
            immonium: false,
            m: false,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: true,
            glycan: (false, (0, 0), Vec::new()),
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }
}

/// A location, or range of locations where an ion can be generated
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug, Serialize, Deserialize)]
pub enum Location {
    /// Skip the given number from the N terminal side
    SkipN(usize),
    /// Skip the given number of aminoacids from the N terminal and C terminal side respectively, only using the positions between these two
    SkipNC(usize, usize),
    /// Skip a certain number and then take a certain number of aminoacids
    TakeN {
        /// Skip this number of aminoacids
        skip: usize,
        /// Take this number of aminoacids
        take: usize,
    },
    /// Skip a given number from the C terminal side
    SkipC(usize),
    /// Take a given number of aminoacids from the C terminal side
    TakeC(usize),
    /// All positions (including 0 and len-1)
    All,
    /// Do not allow it anywhere
    #[default]
    None,
}

impl Location {
    /// Determine if an ion is possible on this location
    pub const fn possible(&self, position: PeptidePosition) -> bool {
        match self {
            Self::SkipN(n) => position.sequence_index >= *n,
            Self::SkipNC(n, c) => {
                position.sequence_index >= *n
                    && position.sequence_length - position.sequence_index > *c
            }
            Self::TakeN { skip, take } => {
                position.sequence_index >= *skip && position.sequence_index < *skip + *take
            }
            Self::SkipC(n) => position.sequence_length - position.sequence_index > *n,
            Self::TakeC(n) => position.sequence_length - position.sequence_index <= *n,
            Self::All => position.series_number != position.sequence_length,
            Self::None => false,
        }
    }
}

#[test]
fn location_all() {
    let all = Model::all();
    let ions_n0 = all.ions(PeptidePosition::n(0, 2));
    let ions_c0 = all.ions(PeptidePosition::c(0, 2));
    assert!(ions_n0.a.0);
    assert!(!ions_n0.x.0);
    assert!(!ions_c0.a.0);
    assert!(ions_c0.x.0);
}
