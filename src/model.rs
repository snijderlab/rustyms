use std::fmt::Display;

use crate::{
    da,
    system::{f64::MassOverCharge, mass_over_charge::mz},
    HasMass,
};

pub struct Model {
    pub a: (Location, Vec<NeutralLoss>),
    pub b: (Location, Vec<NeutralLoss>),
    pub c: (Location, Vec<NeutralLoss>),
    pub d: (Location, Vec<NeutralLoss>),
    pub v: (Location, Vec<NeutralLoss>),
    pub w: (Location, Vec<NeutralLoss>),
    pub x: (Location, Vec<NeutralLoss>),
    pub y: (Location, Vec<NeutralLoss>),
    pub z: (Location, Vec<NeutralLoss>),
    pub precursor: Vec<NeutralLoss>,
    pub ppm: MassOverCharge,
}
#[allow(clippy::struct_excessive_bools)]
pub struct PossibleIons<'a> {
    pub a: (bool, &'a [NeutralLoss]),
    pub b: (bool, &'a [NeutralLoss]),
    pub c: (bool, &'a [NeutralLoss]),
    pub d: (bool, &'a [NeutralLoss]),
    pub v: (bool, &'a [NeutralLoss]),
    pub w: (bool, &'a [NeutralLoss]),
    pub x: (bool, &'a [NeutralLoss]),
    pub y: (bool, &'a [NeutralLoss]),
    pub z: (bool, &'a [NeutralLoss]),
    pub precursor: &'a [NeutralLoss],
}

impl<'a> PossibleIons<'a> {
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
    }
}

impl Model {
    pub fn ions<'a>(&'a self, index: usize, length: usize) -> PossibleIons<'a> {
        PossibleIons {
            a: (self.a.0.possible(index, length), self.a.1.as_slice()),
            b: (self.b.0.possible(index, length), self.b.1.as_slice()),
            c: (self.c.0.possible(index, length), self.c.1.as_slice()),
            d: (self.d.0.possible(index, length), self.d.1.as_slice()),
            v: (self.v.0.possible(index, length), self.v.1.as_slice()),
            w: (self.w.0.possible(index, length), self.w.1.as_slice()),
            x: (self.x.0.possible(index, length), self.x.1.as_slice()),
            y: (self.y.0.possible(index, length), self.y.1.as_slice()),
            z: (self.z.0.possible(index, length), self.z.1.as_slice()),
            precursor: self.precursor.as_slice(),
        }
    }

    #[allow(clippy::too_many_arguments, clippy::many_single_char_names)]
    pub fn new(
        a: (Location, Vec<NeutralLoss>),
        b: (Location, Vec<NeutralLoss>),
        c: (Location, Vec<NeutralLoss>),
        d: (Location, Vec<NeutralLoss>),
        v: (Location, Vec<NeutralLoss>),
        w: (Location, Vec<NeutralLoss>),
        x: (Location, Vec<NeutralLoss>),
        y: (Location, Vec<NeutralLoss>),
        z: (Location, Vec<NeutralLoss>),
        precursor: Vec<NeutralLoss>,
        ppm: MassOverCharge,
    ) -> Self {
        Self {
            a,
            b,
            c,
            d,
            v,
            w,
            x,
            y,
            z,
            precursor,
            ppm,
        }
    }

    pub fn all() -> Self {
        Self {
            a: (Location::SkipN(1), vec![NeutralLoss::Water]),
            b: (Location::SkipN(1), vec![NeutralLoss::Water]),
            c: (Location::SkipN(1), vec![NeutralLoss::Water]),
            d: (Location::SkipN(1), vec![NeutralLoss::Water]),
            v: (Location::SkipN(1), vec![NeutralLoss::Water]),
            w: (Location::All, vec![NeutralLoss::Water]),
            x: (Location::All, vec![NeutralLoss::Water]),
            y: (Location::All, vec![NeutralLoss::Water]),
            z: (Location::All, vec![NeutralLoss::Water]),
            precursor: vec![NeutralLoss::Water],
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }

    pub fn ethcd() -> Self {
        Self {
            a: (Location::None, Vec::new()),
            b: (Location::SkipNC(1, 1), vec![NeutralLoss::Water]),
            c: (Location::SkipNC(1, 1), vec![NeutralLoss::Water]),
            d: (Location::None, Vec::new()),
            v: (Location::None, Vec::new()),
            w: (Location::SkipN(1), vec![NeutralLoss::Water]),
            x: (Location::None, Vec::new()),
            y: (Location::SkipN(1), vec![NeutralLoss::Water]),
            z: (Location::SkipN(1), vec![NeutralLoss::Water]),
            precursor: Vec::new(),
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }

    pub fn cid_hcd() -> Self {
        Self {
            a: (
                Location::TakeN { skip: 1, take: 1 },
                vec![NeutralLoss::Water],
            ),
            b: (Location::SkipNC(1, 1), vec![NeutralLoss::Water]),
            c: (Location::None, Vec::new()),
            d: (Location::TakeN { skip: 1, take: 1 }, Vec::new()),
            v: (Location::None, Vec::new()),
            w: (Location::None, Vec::new()),
            x: (Location::None, Vec::new()),
            y: (Location::SkipC(1), vec![NeutralLoss::Water]),
            z: (Location::None, Vec::new()),
            precursor: vec![NeutralLoss::Water],
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }

    pub fn etcid() -> Self {
        Self {
            a: (Location::None, Vec::new()),
            b: (Location::SkipNC(1, 1), vec![NeutralLoss::Water]),
            c: (Location::SkipNC(1, 1), vec![NeutralLoss::Water]),
            d: (Location::None, Vec::new()),
            v: (Location::None, Vec::new()),
            w: (Location::SkipN(1), vec![NeutralLoss::Water]),
            x: (Location::None, Vec::new()),
            y: (Location::SkipN(1), vec![NeutralLoss::Water]),
            z: (Location::SkipN(1), vec![NeutralLoss::Water]),
            precursor: Vec::new(),
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }
}

pub enum Location {
    SkipN(usize),
    SkipNC(usize, usize),
    TakeN { skip: usize, take: usize },
    SkipC(usize),
    TakeC(usize),
    All,
    None,
}

impl Location {
    pub const fn possible(&self, index: usize, length: usize) -> bool {
        match self {
            Self::SkipN(n) => index >= *n,
            Self::SkipNC(n, c) => index >= *n && length - index > *c,
            Self::TakeN { skip, take } => index >= *skip && index < *skip + *take,
            Self::SkipC(n) => length - index > *n,
            Self::TakeC(n) => length - index <= *n,
            Self::All => true,
            Self::None => false,
        }
    }
}
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum NeutralLoss {
    Water,
    Ammonia,
    CarbonMonoxide,
}

impl HasMass for NeutralLoss {
    fn mass<M: crate::MassSystem>(&self) -> crate::Mass {
        match self {
            Self::Water => da(M::OH + M::H),
            Self::Ammonia => da(M::NH3),
            Self::CarbonMonoxide => da(M::CO),
        }
    }
}
impl Display for NeutralLoss {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Water => "Water",
                Self::Ammonia => "Ammonia",
                Self::CarbonMonoxide => "CarbonMonoxide",
            }
        )
    }
}
