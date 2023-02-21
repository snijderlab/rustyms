pub struct Model {
    pub a: Location,
    pub b: Location,
    pub c: Location,
    pub d: Location,
    pub v: Location,
    pub w: Location,
    pub x: Location,
    pub y: Location,
    pub z: Location,
}

pub struct PossibleIons {
    pub a: bool,
    pub b: bool,
    pub c: bool,
    pub d: bool,
    pub v: bool,
    pub w: bool,
    pub x: bool,
    pub y: bool,
    pub z: bool,
}

impl Model {
    pub fn ions(&self, index: usize, length: usize) -> PossibleIons {
        PossibleIons {
            a: self.a.possible(index, length),
            b: self.b.possible(index, length),
            c: self.c.possible(index, length),
            d: self.d.possible(index, length),
            v: self.v.possible(index, length),
            w: self.w.possible(index, length),
            x: self.x.possible(index, length),
            y: self.y.possible(index, length),
            z: self.z.possible(index, length),
        }
    }

    pub fn new(
        a: Location,
        b: Location,
        c: Location,
        d: Location,
        v: Location,
        w: Location,
        x: Location,
        y: Location,
        z: Location,
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
        }
    }

    pub const fn all() -> Self {
        Self {
            a: Location::SkipC(1),
            b: Location::SkipC(1),
            c: Location::SkipC(1),
            d: Location::None,
            v: Location::None,
            w: Location::None,
            x: Location::SkipN(1),
            y: Location::SkipN(1),
            z: Location::SkipN(1),
        }
    }
}

pub enum Location {
    SkipN(usize),
    TakeN(usize),
    SkipC(usize),
    TakeC(usize),
    All,
    None,
}

impl Location {
    pub fn possible(&self, index: usize, length: usize) -> bool {
        match self {
            Location::SkipN(n) => index + 1 > *n,
            Location::TakeN(n) => index < *n,
            Location::SkipC(n) => length - index > *n,
            Location::TakeC(n) => length - index <= *n,
            Location::All => true,
            Location::None => false,
        }
    }
}
