#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Element {
    /// Element Hydrogen (H) atomic number: 1
    H = 1,
    /// Element Helium (He) atomic number: 2
    He,
    /// Element Lithium (Li) atomic number: 3
    Li,
    /// Element Beryllium (Be) atomic number: 4
    Be,
    /// Element Boron (B) atomic number: 5
    B,
    /// Element Carbon (C) atomic number: 6
    C,
    /// Element Nitrogen (N) atomic number: 7
    N,
    /// Element Oxygen (O) atomic number: 8
    O,
    /// Element Fluorine (F) atomic number: 9
    F,
    /// Element Neon (Ne) atomic number: 10
    Ne,
    /// Element Sodium (Na) atomic number: 11
    Na,
    /// Element Magnesium (Mg) atomic number: 12
    Mg,
    /// Element Aluminium (Al) atomic number: 13
    Al,
    /// Element Silicon (Si) atomic number: 14
    Si,
    /// Element Phosphorus (P) atomic number: 15
    P,
    /// Element Sulfur (S) atomic number: 16
    S,
    /// Element Chlorine (Cl) atomic number: 17
    Cl,
    /// Element Argon (Ar) atomic number: 18
    Ar,
    /// Element Potassium (K) atomic number: 19
    K,
    /// Element Calcium (Ca) atomic number: 20
    Ca,
    /// Element Scandium (Sc) atomic number: 21
    Sc,
    /// Element Titanium (Ti) atomic number: 22
    Ti,
    /// Element Vanadium (V) atomic number: 23
    V,
    /// Element Chromium (Cr) atomic number: 24
    Cr,
    /// Element Manganese (Mn) atomic number: 25
    Mn,
    /// Element Iron (Fe) atomic number: 26
    Fe,
    /// Element Cobalt (Co) atomic number: 27
    Co,
    /// Element Nickel (Ni) atomic number: 28
    Ni,
    /// Element Copper (Cu) atomic number: 29
    Cu,
    /// Element Zinc (Zn) atomic number: 30
    Zn,
    /// Element Gallium (Ga) atomic number: 31
    Ga,
    /// Element Germanium (Ge) atomic number: 32
    Ge,
    /// Element Arsenic (As) atomic number: 33
    As,
    /// Element Selenium (Se) atomic number: 34
    Se,
    /// Element Bromine (Br) atomic number: 35
    Br,
    /// Element Krypton (Kr) atomic number: 36
    Kr,
    /// Element Rubidium (Rb) atomic number: 37
    Rb,
    /// Element Strontium (Sr) atomic number: 38
    Sr,
    /// Element Yttrium (Y) atomic number: 39
    Y,
    /// Element Zirconium (Zr) atomic number: 40
    Zr,
    /// Element Niobium (Nb) atomic number: 41
    Nb,
    /// Element Molybdenum (Mo) atomic number: 42
    Mo,
    /// Element Technetium (Tc) atomic number: 43
    Tc,
    /// Element Ruthenium (Ru) atomic number: 44
    Ru,
    /// Element Rhodium (Rh) atomic number: 45
    Rh,
    /// Element Palladium (Pd) atomic number: 46
    Pd,
    /// Element Silver (Ag) atomic number: 47
    Ag,
    /// Element Cadmium (Cd) atomic number: 48
    Cd,
    /// Element Indium (In) atomic number: 49
    In,
    /// Element Tin (Sn) atomic number: 50
    Sn,
    /// Element Antimony (Sb) atomic number: 51
    Sb,
    /// Element Tellurium (Te) atomic number: 52
    Te,
    /// Element Iodine (I) atomic number: 53
    I,
    /// Element Xenon (Xe) atomic number: 54
    Xe,
    /// Element Caesium (Cs) atomic number: 55
    Cs,
    /// Element Barium (Ba) atomic number: 56
    Ba,
    /// Element Lanthanum (La) atomic number: 57
    La,
    /// Element Cerium (Ce) atomic number: 58
    Ce,
    /// Element Praseodymium (Pr) atomic number: 59
    Pr,
    /// Element Neodymium (Nd) atomic number: 60
    Nd,
    /// Element Promethium (Pm) atomic number: 61
    Pm,
    /// Element Samarium (Sm) atomic number: 62
    Sm,
    /// Element Europium (Eu) atomic number: 63
    Eu,
    /// Element Gadolinium (Gd) atomic number: 64
    Gd,
    /// Element Terbium (Tb) atomic number: 65
    Tb,
    /// Element Dysprosium (Dy) atomic number: 66
    Dy,
    /// Element Holmium (Ho) atomic number: 67
    Ho,
    /// Element Erbium (Er) atomic number: 68
    Er,
    /// Element Thulium (Tm) atomic number: 69
    Tm,
    /// Element Ytterbium (Yb) atomic number: 70
    Yb,
    /// Element Lutetium (Lu) atomic number: 71
    Lu,
    /// Element Hafnium (Hf) atomic number: 72
    Hf,
    /// Element Tantalum (Ta) atomic number: 73
    Ta,
    /// Element Tungsten (W) atomic number: 74
    W,
    /// Element Rhenium (Re) atomic number: 75
    Re,
    /// Element Osmium (Os) atomic number: 76
    Os,
    /// Element Iridium (Ir) atomic number: 77
    Ir,
    /// Element Platinum (Pt) atomic number: 78
    Pt,
    /// Element Gold (Au) atomic number: 79
    Au,
    /// Element Mercury (Hg) atomic number: 80
    Hg,
    /// Element Thallium (Tl) atomic number: 81
    Tl,
    /// Element Lead (Pb) atomic number: 82
    Pb,
    /// Element Bismuth (Bi) atomic number: 83
    Bi,
    /// Element Polonium (Po) atomic number: 84
    Po,
    /// Element Astatine (At) atomic number: 85
    At,
    /// Element Radon (Rn) atomic number: 86
    Rn,
    /// Element Francium (Fr) atomic number: 87
    Fr,
    /// Element Radium (Ra) atomic number: 88
    Ra,
    /// Element Actinium (Ac) atomic number: 89
    Ac,
    /// Element Thorium (Th) atomic number: 90
    Th,
    /// Element Protactinium (Pa) atomic number: 91
    Pa,
    /// Element Uranium (U) atomic number: 92
    U,
    /// Element Neptunium (Np) atomic number: 93
    Np,
    /// Element Plutonium (Pu) atomic number: 94
    Pu,
    /// Element Americium (Am) atomic number: 95
    Am,
    /// Element Curium (Cm) atomic number: 96
    Cm,
    /// Element Berkelium (Bk) atomic number: 97
    Bk,
    /// Element Californium (Cf) atomic number: 98
    Cf,
    /// Element Einsteinium (Es) atomic number: 99
    Es,
    /// Element Fermium (Fm) atomic number: 100
    Fm,
    /// Element Mendelevium (Md) atomic number: 101
    Md,
    /// Element Nobelium (No) atomic number: 102
    No,
    /// Element Lawrencium (Lr) atomic number: 103
    Lr,
    /// Element Rutherfordium (Rf) atomic number: 104
    Rf,
    /// Element Dubnium (Db) atomic number: 105
    Db,
    /// Element Seaborgium (Sg) atomic number: 106
    Sg,
    /// Element Bohrium (Bh) atomic number: 107
    Bh,
    /// Element Hassium (Hs) atomic number: 108
    Hs,
    /// Element Meitnerium (Mt) atomic number: 109
    Mt,
    /// Element Darmstadtium (Ds) atomic number: 110
    Ds,
    /// Element Roentgenium (Rg) atomic number: 111
    Rg,
    /// Element Copernicium (Cn) atomic number: 112
    Cn,
    /// Element Nihonium (Nh) atomic number: 113
    Nh,
    /// Element Flerovium (Fl) atomic number: 114
    Fl,
    /// Element Moscovium (Mc) atomic number: 115
    Mc,
    /// Element Livermorium (Lv) atomic number: 116
    Lv,
    /// Element Tennessine (Ts) atomic number: 117
    Ts,
    /// Element Oganesson (Og) atomic number: 118
    Og,
}

impl TryFrom<&str> for Element {
    type Error = ();
    #[allow(clippy::too_many_lines)]
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "H" => Ok(Self::H),
            "He" => Ok(Self::He),
            "Li" => Ok(Self::Li),
            "Be" => Ok(Self::Be),
            "B" => Ok(Self::B),
            "C" => Ok(Self::C),
            "N" => Ok(Self::N),
            "O" => Ok(Self::O),
            "F" => Ok(Self::F),
            "Ne" => Ok(Self::Ne),
            "Na" => Ok(Self::Na),
            "Mg" => Ok(Self::Mg),
            "Al" => Ok(Self::Al),
            "Si" => Ok(Self::Si),
            "P" => Ok(Self::P),
            "S" => Ok(Self::S),
            "Cl" => Ok(Self::Cl),
            "Ar" => Ok(Self::Ar),
            "K" => Ok(Self::K),
            "Ca" => Ok(Self::Ca),
            "Sc" => Ok(Self::Sc),
            "Ti" => Ok(Self::Ti),
            "V" => Ok(Self::V),
            "Cr" => Ok(Self::Cr),
            "Mn" => Ok(Self::Mn),
            "Fe" => Ok(Self::Fe),
            "Co" => Ok(Self::Co),
            "Ni" => Ok(Self::Ni),
            "Cu" => Ok(Self::Cu),
            "Zn" => Ok(Self::Zn),
            "Ga" => Ok(Self::Ga),
            "Ge" => Ok(Self::Ge),
            "As" => Ok(Self::As),
            "Se" => Ok(Self::Se),
            "Br" => Ok(Self::Br),
            "Kr" => Ok(Self::Kr),
            "Rb" => Ok(Self::Rb),
            "Sr" => Ok(Self::Sr),
            "Y" => Ok(Self::Y),
            "Zr" => Ok(Self::Zr),
            "Nb" => Ok(Self::Nb),
            "Mo" => Ok(Self::Mo),
            "Tc" => Ok(Self::Tc),
            "Ru" => Ok(Self::Ru),
            "Rh" => Ok(Self::Rh),
            "Pd" => Ok(Self::Pd),
            "Ag" => Ok(Self::Ag),
            "Cd" => Ok(Self::Cd),
            "In" => Ok(Self::In),
            "Sn" => Ok(Self::Sn),
            "Sb" => Ok(Self::Sb),
            "Te" => Ok(Self::Te),
            "I" => Ok(Self::I),
            "Xe" => Ok(Self::Xe),
            "Cs" => Ok(Self::Cs),
            "Ba" => Ok(Self::Ba),
            "La" => Ok(Self::La),
            "Ce" => Ok(Self::Ce),
            "Pr" => Ok(Self::Pr),
            "Nd" => Ok(Self::Nd),
            "Pm" => Ok(Self::Pm),
            "Sm" => Ok(Self::Sm),
            "Eu" => Ok(Self::Eu),
            "Gd" => Ok(Self::Gd),
            "Tb" => Ok(Self::Tb),
            "Dy" => Ok(Self::Dy),
            "Ho" => Ok(Self::Ho),
            "Er" => Ok(Self::Er),
            "Tm" => Ok(Self::Tm),
            "Yb" => Ok(Self::Yb),
            "Lu" => Ok(Self::Lu),
            "Hf" => Ok(Self::Hf),
            "Ta" => Ok(Self::Ta),
            "W" => Ok(Self::W),
            "Re" => Ok(Self::Re),
            "Os" => Ok(Self::Os),
            "Ir" => Ok(Self::Ir),
            "Pt" => Ok(Self::Pt),
            "Au" => Ok(Self::Au),
            "Hg" => Ok(Self::Hg),
            "Tl" => Ok(Self::Tl),
            "Pb" => Ok(Self::Pb),
            "Bi" => Ok(Self::Bi),
            "Po" => Ok(Self::Po),
            "At" => Ok(Self::At),
            "Rn" => Ok(Self::Rn),
            "Fr" => Ok(Self::Fr),
            "Ra" => Ok(Self::Ra),
            "Ac" => Ok(Self::Ac),
            "Th" => Ok(Self::Th),
            "Pa" => Ok(Self::Pa),
            "U" => Ok(Self::U),
            "Np" => Ok(Self::Np),
            "Pu" => Ok(Self::Pu),
            "Am" => Ok(Self::Am),
            "Cm" => Ok(Self::Cm),
            "Bk" => Ok(Self::Bk),
            "Cf" => Ok(Self::Cf),
            "Es" => Ok(Self::Es),
            "Fm" => Ok(Self::Fm),
            "Md" => Ok(Self::Md),
            "No" => Ok(Self::No),
            "Lr" => Ok(Self::Lr),
            "Rf" => Ok(Self::Rf),
            "Db" => Ok(Self::Db),
            "Sg" => Ok(Self::Sg),
            "Bh" => Ok(Self::Bh),
            "Hs" => Ok(Self::Hs),
            "Mt" => Ok(Self::Mt),
            "Ds" => Ok(Self::Ds),
            "Rg" => Ok(Self::Rg),
            "Cn" => Ok(Self::Cn),
            "Nh" => Ok(Self::Nh),
            "Fl" => Ok(Self::Fl),
            "Mc" => Ok(Self::Mc),
            "Lv" => Ok(Self::Lv),
            "Ts" => Ok(Self::Ts),
            "Og" => Ok(Self::Og),
            _ => Err(()),
        }
    }
}

impl std::fmt::Display for Element {
    #[allow(clippy::too_many_lines)]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::H => "H",
                Self::He => "He",
                Self::Li => "Li",
                Self::Be => "Be",
                Self::B => "B",
                Self::C => "C",
                Self::N => "N",
                Self::O => "O",
                Self::F => "F",
                Self::Ne => "Ne",
                Self::Na => "Na",
                Self::Mg => "Mg",
                Self::Al => "Al",
                Self::Si => "Si",
                Self::P => "P",
                Self::S => "S",
                Self::Cl => "Cl",
                Self::Ar => "Ar",
                Self::K => "K",
                Self::Ca => "Ca",
                Self::Sc => "Sc",
                Self::Ti => "Ti",
                Self::V => "V",
                Self::Cr => "Cr",
                Self::Mn => "Mn",
                Self::Fe => "Fe",
                Self::Co => "Co",
                Self::Ni => "Ni",
                Self::Cu => "Cu",
                Self::Zn => "Zn",
                Self::Ga => "Ga",
                Self::Ge => "Ge",
                Self::As => "As",
                Self::Se => "Se",
                Self::Br => "Br",
                Self::Kr => "Kr",
                Self::Rb => "Rb",
                Self::Sr => "Sr",
                Self::Y => "Y",
                Self::Zr => "Zr",
                Self::Nb => "Nb",
                Self::Mo => "Mo",
                Self::Tc => "Tc",
                Self::Ru => "Ru",
                Self::Rh => "Rh",
                Self::Pd => "Pd",
                Self::Ag => "Ag",
                Self::Cd => "Cd",
                Self::In => "In",
                Self::Sn => "Sn",
                Self::Sb => "Sb",
                Self::Te => "Te",
                Self::I => "I",
                Self::Xe => "Xe",
                Self::Cs => "Cs",
                Self::Ba => "Ba",
                Self::La => "La",
                Self::Ce => "Ce",
                Self::Pr => "Pr",
                Self::Nd => "Nd",
                Self::Pm => "Pm",
                Self::Sm => "Sm",
                Self::Eu => "Eu",
                Self::Gd => "Gd",
                Self::Tb => "Tb",
                Self::Dy => "Dy",
                Self::Ho => "Ho",
                Self::Er => "Er",
                Self::Tm => "Tm",
                Self::Yb => "Yb",
                Self::Lu => "Lu",
                Self::Hf => "Hf",
                Self::Ta => "Ta",
                Self::W => "W",
                Self::Re => "Re",
                Self::Os => "Os",
                Self::Ir => "Ir",
                Self::Pt => "Pt",
                Self::Au => "Au",
                Self::Hg => "Hg",
                Self::Tl => "Tl",
                Self::Pb => "Pb",
                Self::Bi => "Bi",
                Self::Po => "Po",
                Self::At => "At",
                Self::Rn => "Rn",
                Self::Fr => "Fr",
                Self::Ra => "Ra",
                Self::Ac => "Ac",
                Self::Th => "Th",
                Self::Pa => "Pa",
                Self::U => "U",
                Self::Np => "Np",
                Self::Pu => "Pu",
                Self::Am => "Am",
                Self::Cm => "Cm",
                Self::Bk => "Bk",
                Self::Cf => "Cf",
                Self::Es => "Es",
                Self::Fm => "Fm",
                Self::Md => "Md",
                Self::No => "No",
                Self::Lr => "Lr",
                Self::Rf => "Rf",
                Self::Db => "Db",
                Self::Sg => "Sg",
                Self::Bh => "Bh",
                Self::Hs => "Hs",
                Self::Mt => "Mt",
                Self::Ds => "Ds",
                Self::Rg => "Rg",
                Self::Cn => "Cn",
                Self::Nh => "Nh",
                Self::Fl => "Fl",
                Self::Mc => "Mc",
                Self::Lv => "Lv",
                Self::Ts => "Ts",
                Self::Og => "Og",
            }
        )
    }
}

pub const ELEMENT_PARSE_LIST: &[(&str, Element)] = &[
    ("He", Element::He),
    ("Li", Element::Li),
    ("Be", Element::Be),
    ("Ne", Element::Ne),
    ("Na", Element::Na),
    ("Mg", Element::Mg),
    ("Al", Element::Al),
    ("Si", Element::Si),
    ("Cl", Element::Cl),
    ("Ar", Element::Ar),
    ("Ca", Element::Ca),
    ("Sc", Element::Sc),
    ("Ti", Element::Ti),
    ("Cr", Element::Cr),
    ("Mn", Element::Mn),
    ("Fe", Element::Fe),
    ("Co", Element::Co),
    ("Ni", Element::Ni),
    ("Cu", Element::Cu),
    ("Zn", Element::Zn),
    ("Ga", Element::Ga),
    ("Ge", Element::Ge),
    ("As", Element::As),
    ("Se", Element::Se),
    ("Br", Element::Br),
    ("Kr", Element::Kr),
    ("Rb", Element::Rb),
    ("Sr", Element::Sr),
    ("Zr", Element::Zr),
    ("Nb", Element::Nb),
    ("Mo", Element::Mo),
    ("Tc", Element::Tc),
    ("Ru", Element::Ru),
    ("Rh", Element::Rh),
    ("Pd", Element::Pd),
    ("Ag", Element::Ag),
    ("Cd", Element::Cd),
    ("In", Element::In),
    ("Sn", Element::Sn),
    ("Sb", Element::Sb),
    ("Te", Element::Te),
    ("Xe", Element::Xe),
    ("Cs", Element::Cs),
    ("Ba", Element::Ba),
    ("La", Element::La),
    ("Ce", Element::Ce),
    ("Pr", Element::Pr),
    ("Nd", Element::Nd),
    ("Pm", Element::Pm),
    ("Sm", Element::Sm),
    ("Eu", Element::Eu),
    ("Gd", Element::Gd),
    ("Tb", Element::Tb),
    ("Dy", Element::Dy),
    ("Ho", Element::Ho),
    ("Er", Element::Er),
    ("Tm", Element::Tm),
    ("Yb", Element::Yb),
    ("Lu", Element::Lu),
    ("Hf", Element::Hf),
    ("Ta", Element::Ta),
    ("Re", Element::Re),
    ("Os", Element::Os),
    ("Ir", Element::Ir),
    ("Pt", Element::Pt),
    ("Au", Element::Au),
    ("Hg", Element::Hg),
    ("Tl", Element::Tl),
    ("Pb", Element::Pb),
    ("Bi", Element::Bi),
    ("Po", Element::Po),
    ("At", Element::At),
    ("Rn", Element::Rn),
    ("Fr", Element::Fr),
    ("Ra", Element::Ra),
    ("Ac", Element::Ac),
    ("Th", Element::Th),
    ("Pa", Element::Pa),
    ("Np", Element::Np),
    ("Pu", Element::Pu),
    ("Am", Element::Am),
    ("Cm", Element::Cm),
    ("Bk", Element::Bk),
    ("Cf", Element::Cf),
    ("Es", Element::Es),
    ("Fm", Element::Fm),
    ("Md", Element::Md),
    ("No", Element::No),
    ("Lr", Element::Lr),
    ("Rf", Element::Rf),
    ("Db", Element::Db),
    ("Sg", Element::Sg),
    ("Bh", Element::Bh),
    ("Hs", Element::Hs),
    ("Mt", Element::Mt),
    ("Ds", Element::Ds),
    ("Rg", Element::Rg),
    ("Cn", Element::Cn),
    ("Nh", Element::Nh),
    ("Fl", Element::Fl),
    ("Mc", Element::Mc),
    ("Lv", Element::Lv),
    ("Ts", Element::Ts),
    ("Og", Element::Og),
    ("U", Element::U),
    ("W", Element::W),
    ("I", Element::I),
    ("Y", Element::Y),
    ("V", Element::V),
    ("K", Element::K),
    ("S", Element::S),
    ("B", Element::B),
    ("C", Element::C),
    ("N", Element::N),
    ("O", Element::O),
    ("F", Element::F),
    ("H", Element::H),
    ("P", Element::P),
];
