use crate::LinearPeptide;
use serde::{Deserialize, Serialize};

/// An annotated peptide
pub trait AnnotatedPeptide {
    /// The complexity of the peptide
    type Complexity;
    /// Get the peptide
    fn peptide(&self) -> &LinearPeptide<Self::Complexity>;
    /// Get the regions, as a list of the regions in order with th length of each region, these are
    /// required to be as long as the full peptide.
    fn regions(&self) -> &[(Region, usize)];
    /// Get the annotations, specified as the annotation and the into into the peptide where it is
    /// located.
    fn annotations(&self) -> &[(Annotation, usize)];
}

/// A region on an antibody
#[allow(missing_docs)]
#[derive(Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Debug)]
pub enum Region {
    Framework(usize),
    ComplementarityDeterminingRegion(usize),
    Hinge(Option<usize>),
    ConstantHeavy(usize),
    ConstantLight,
    SecratoryTail,
    MembraneTail(Option<usize>),
    Other(String),
    None,
}

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Framework(n) => write!(f, "FR{n}"),
            Self::ComplementarityDeterminingRegion(n) => write!(f, "CDR{n}"),
            Self::Hinge(n) => write!(f, "H{}", n.map_or(String::new(), |n| n.to_string())),
            Self::ConstantHeavy(n) => write!(f, "CH{n}"),
            Self::ConstantLight => write!(f, "CL"),
            Self::SecratoryTail => write!(f, "CHS"),
            Self::MembraneTail(n) => write!(f, "M{}", n.map_or(String::new(), |n| n.to_string())),
            Self::Other(o) => write!(f, "{o}"),
            Self::None => Ok(()),
        }
    }
}

impl std::str::FromStr for Region {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "" => Self::None,
            "CL" => Self::ConstantLight,
            "CHS" => Self::SecratoryTail,
            "H" => Self::Hinge(None),
            "M" => Self::MembraneTail(None),
            cdr if cdr.starts_with("CDR") => cdr[3..]
                .parse::<usize>()
                .map_or(Self::Other(cdr.to_string()), |c| {
                    Self::ComplementarityDeterminingRegion(c)
                }),
            fr if fr.starts_with("FR") => fr[2..]
                .parse::<usize>()
                .map_or(Self::Other(fr.to_string()), Self::Framework),
            ch if ch.starts_with("CH") => ch[2..]
                .parse::<usize>()
                .map_or(Self::Other(ch.to_string()), Self::ConstantHeavy),
            h if h.starts_with('H') => h[1..]
                .parse::<usize>()
                .map_or(Self::Other(h.to_string()), |c| Self::Hinge(Some(c))),
            m if m.starts_with('M') => m[1..]
                .parse::<usize>()
                .map_or(Self::Other(m.to_string()), |c| Self::MembraneTail(Some(c))),
            o => Self::Other(o.to_string()),
        })
    }
}

/// A sequence annotation
#[derive(Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Debug)]
pub enum Annotation {
    /// A conserved residue
    Conserved,
    /// A potential N linked glycan position
    NGlycan,
    /// Any other annotation
    Other(String),
}

impl std::str::FromStr for Annotation {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "C" | "Conserved" => Self::Conserved,
            "N" | "NGlycan" => Self::NGlycan,
            o => Self::Other(o.to_string()),
        })
    }
}
