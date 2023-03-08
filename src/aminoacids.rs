use uom::num_traits::Zero;

use crate::fragment::{Fragment, FragmentType};
use crate::system::f64::*;
use crate::{model::*, HasMass};
use crate::{MassSystem, Position};

#[derive(Clone, Debug, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AminoAcid {
    Alanine,
    Arginine,
    Asparagine,
    AsparticAcid,
    Cysteine,
    Glutamine,
    GlutamicAcid,
    Glycine,
    Histidine,
    Isoleucine,
    Leucine,
    AmbiguousLeucine,
    Lysine,
    Methionine,
    Phenylalanine,
    Proline,
    Serine,
    Threonine,
    Tryptophan,
    Tyrosine,
    Valine,
    Selenocysteine,
    Pyrrolysine,
}

impl TryFrom<char> for AminoAcid {
    type Error = ();
    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'A' => Ok(Self::Alanine),
            'C' => Ok(Self::Cysteine),
            'D' => Ok(Self::AsparticAcid),
            'E' => Ok(Self::GlutamicAcid),
            'F' => Ok(Self::Phenylalanine),
            'G' => Ok(Self::Glycine),
            'H' => Ok(Self::Histidine),
            'I' => Ok(Self::Isoleucine),
            'J' => Ok(Self::AmbiguousLeucine),
            'K' => Ok(Self::Lysine),
            'L' => Ok(Self::Leucine),
            'M' => Ok(Self::Methionine),
            'N' => Ok(Self::Asparagine),
            'O' => Ok(Self::Pyrrolysine),
            'P' => Ok(Self::Proline),
            'Q' => Ok(Self::Glutamine),
            'R' => Ok(Self::Arginine),
            'S' => Ok(Self::Serine),
            'T' => Ok(Self::Threonine),
            'U' => Ok(Self::Selenocysteine),
            'V' => Ok(Self::Valine),
            'W' => Ok(Self::Tryptophan),
            'Y' => Ok(Self::Tyrosine),
            _ => Err(()),
        }
    }
}

impl TryFrom<&u8> for AminoAcid {
    type Error = ();
    fn try_from(value: &u8) -> Result<Self, Self::Error> {
        match value {
            b'A' => Ok(Self::Alanine),
            b'C' => Ok(Self::Cysteine),
            b'D' => Ok(Self::AsparticAcid),
            b'E' => Ok(Self::GlutamicAcid),
            b'F' => Ok(Self::Phenylalanine),
            b'G' => Ok(Self::Glycine),
            b'H' => Ok(Self::Histidine),
            b'I' => Ok(Self::Isoleucine),
            b'J' => Ok(Self::AmbiguousLeucine),
            b'K' => Ok(Self::Lysine),
            b'L' => Ok(Self::Leucine),
            b'M' => Ok(Self::Methionine),
            b'N' => Ok(Self::Asparagine),
            b'O' => Ok(Self::Pyrrolysine),
            b'P' => Ok(Self::Proline),
            b'Q' => Ok(Self::Glutamine),
            b'R' => Ok(Self::Arginine),
            b'S' => Ok(Self::Serine),
            b'T' => Ok(Self::Threonine),
            b'U' => Ok(Self::Selenocysteine),
            b'V' => Ok(Self::Valine),
            b'W' => Ok(Self::Tryptophan),
            b'Y' => Ok(Self::Tyrosine),
            _ => Err(()),
        }
    }
}

impl HasMass for AminoAcid {
    fn mass<M: MassSystem>(&self) -> Mass {
        match self {
            Self::Alanine => da(M::BACKBONE + M::CH3),
            Self::AmbiguousLeucine => da(M::BACKBONE + M::C * 4.0 + M::H * 9.0),
            Self::Arginine => da(M::BACKBONE + M::CH2 * 3.0 + M::NH + M::NH2 * 2.0),
            Self::Asparagine => da(M::BACKBONE + M::CH2 + M::C + M::O + M::NH2),
            Self::AsparticAcid => da(M::BACKBONE + M::CH2 + M::C + M::OH + M::O),
            Self::Cysteine => da(M::BACKBONE + M::CH2 + M::S + M::H),
            Self::GlutamicAcid => da(M::BACKBONE + M::CH2 * 2.0 + M::C + M::OH + M::O),
            Self::Glutamine => da(M::BACKBONE + M::CH2 * 2.0 + M::C + M::O + M::NH2),
            Self::Glycine => da(M::BACKBONE + M::H),
            Self::Histidine => da(M::BACKBONE + M::CH2 + M::C + M::N + M::CH + M::NH + M::CH),
            Self::Isoleucine => da(M::BACKBONE + M::CH + M::CH3 + M::CH2 + M::CH3),
            Self::Leucine => da(M::BACKBONE + M::CH2 + M::CH + M::CH3 + M::CH3),
            Self::Lysine => da(M::BACKBONE + M::CH2 * 4.0 + M::NH2),
            Self::Methionine => da(M::BACKBONE + M::CH2 * 2.0 + M::S + M::CH3),
            Self::Phenylalanine => da(M::BACKBONE + M::CH2 + M::C + M::CH * 5.0),
            Self::Proline => da(M::BACKBONE + M::CH2 * 3.0 - M::H),
            Self::Pyrrolysine => da(M::BACKBONE
                + M::CH2 * 4.0
                + M::NH
                + M::C
                + M::O
                + M::CH
                + M::N
                + M::CH
                + M::CH2
                + M::CH
                + M::CH3),
            Self::Selenocysteine => da(M::BACKBONE + M::CH2 + M::Se),
            Self::Serine => da(M::BACKBONE + M::CH2 + M::OH),
            Self::Threonine => da(M::BACKBONE + M::CH + M::OH + M::CH3),
            Self::Tryptophan => {
                da(M::BACKBONE + M::CH2 + M::C * 2.0 + M::CH * 4.0 + M::C + M::NH + M::CH)
            }
            Self::Tyrosine => da(M::BACKBONE + M::CH2 + M::C * 2.0 + M::CH * 4.0 + M::OH),
            Self::Valine => da(M::BACKBONE + M::CH + M::CH3 * 2.0),
        }
    }
}

#[allow(non_upper_case_globals)]
impl AminoAcid {
    pub const A: Self = Self::Alanine;
    pub const C: Self = Self::Cysteine;
    pub const D: Self = Self::AsparticAcid;
    pub const E: Self = Self::GlutamicAcid;
    pub const F: Self = Self::Phenylalanine;
    pub const G: Self = Self::Glycine;
    pub const H: Self = Self::Histidine;
    pub const I: Self = Self::Isoleucine;
    pub const J: Self = Self::AmbiguousLeucine;
    pub const K: Self = Self::Lysine;
    pub const L: Self = Self::Leucine;
    pub const M: Self = Self::Methionine;
    pub const N: Self = Self::Asparagine;
    pub const O: Self = Self::Pyrrolysine;
    pub const P: Self = Self::Proline;
    pub const Q: Self = Self::Glutamine;
    pub const R: Self = Self::Arginine;
    pub const S: Self = Self::Serine;
    pub const T: Self = Self::Threonine;
    pub const U: Self = Self::Selenocysteine;
    pub const V: Self = Self::Valine;
    pub const W: Self = Self::Tryptophan;
    pub const Ala: Self = Self::Alanine;
    pub const Cys: Self = Self::Cysteine;
    pub const Asp: Self = Self::AsparticAcid;
    pub const Glu: Self = Self::GlutamicAcid;
    pub const Phe: Self = Self::Phenylalanine;
    pub const Gly: Self = Self::Glycine;
    pub const His: Self = Self::Histidine;
    pub const Ile: Self = Self::Isoleucine;
    pub const Xle: Self = Self::AmbiguousLeucine;
    pub const Lys: Self = Self::Lysine;
    pub const Leu: Self = Self::Leucine;
    pub const Met: Self = Self::Methionine;
    pub const Asn: Self = Self::Asparagine;
    pub const Pyl: Self = Self::Pyrrolysine;
    pub const Pro: Self = Self::Proline;
    pub const Gln: Self = Self::Glutamine;
    pub const Arg: Self = Self::Arginine;
    pub const Ser: Self = Self::Serine;
    pub const Thr: Self = Self::Threonine;
    pub const Sec: Self = Self::Selenocysteine;
    pub const Val: Self = Self::Valine;
    pub const Trp: Self = Self::Tryptophan;

    // TODO: Take side chain mutations into account (maybe define pyrrolysine as a mutation)
    pub fn satellite_ion_masses<M: MassSystem>(&self) -> Vec<Mass> {
        match self {
            Self::Alanine => vec![],
            Self::AmbiguousLeucine => vec![],
            Self::Arginine => vec![da(M::CH2 * 2.0 + M::NH + M::NH2 * 2.0)],
            Self::Asparagine => vec![da(M::C + M::O + M::NH2)],
            Self::AsparticAcid => vec![da(M::C + M::O * 2.0)],
            Self::Cysteine => vec![da(M::S + M::H)],
            Self::GlutamicAcid => vec![da(M::CH2 + M::C + M::O * 2.0)],
            Self::Glutamine => vec![da(M::CH2 + M::C + M::O + M::NH2)],
            Self::Glycine => vec![],
            Self::Histidine => vec![], // Aromatic
            Self::Isoleucine => vec![da(M::CH3), da(M::CH2 + M::CH3)],
            Self::Leucine => vec![da(M::CH + M::CH3 * 2.0)],
            Self::Lysine => vec![da(M::CH2 * 3.0 + M::NH2)],
            Self::Methionine => vec![da(M::CH2 + M::S + M::CH3)],
            Self::Phenylalanine => vec![], // Aromatic
            Self::Proline => vec![], // Interesting, TODO: see what other software packages think about this one
            Self::Pyrrolysine => {
                vec![da(M::CH2 * 3.0
                    + M::NH
                    + M::C
                    + M::O
                    + M::CH
                    + M::N
                    + M::CH
                    + M::CH2
                    + M::CH
                    + M::CH3)]
            } // Weird, TODO: figure out what to make of this
            Self::Selenocysteine => vec![da(M::Se)],
            Self::Serine => vec![da(M::OH)],
            Self::Threonine => vec![da(M::OH), da(M::CH3)],
            Self::Tryptophan => vec![],       // Aromatic
            Self::Tyrosine => vec![],         // Aromatic
            Self::Valine => vec![da(M::CH3)], // Technically two options, but both have the same mass TODO: check if the loss of both is an option
        }
    }

    pub fn fragments<M: MassSystem>(
        &self,
        n_term: Mass,
        c_term: Mass,
        max_charge: Charge,
        sequence_index: usize,
        sequence_length: usize,
        ions: &PossibleIons,
    ) -> Vec<Fragment> {
        let mut base_fragments = Vec::with_capacity(ions.size_upper_bound());
        if ions.a {
            base_fragments.push(Fragment::new(
                n_term + self.mass::<M>() + da(-M::CO),
                Charge::zero(),
                FragmentType::a(Position::n(sequence_index, sequence_length)),
            ));
        }
        if ions.b {
            base_fragments.push(Fragment::new(
                n_term + self.mass::<M>(),
                Charge::zero(),
                FragmentType::b(Position::n(sequence_index, sequence_length)),
            ));
        }
        if ions.c {
            base_fragments.push(Fragment::new(
                n_term + self.mass::<M>() + da(M::NH3),
                Charge::zero(),
                FragmentType::c(Position::n(sequence_index, sequence_length)),
            ));
        }
        if ions.d {
            for satellite in self.satellite_ion_masses::<M>() {
                base_fragments.push(Fragment::new(
                    n_term + self.mass::<M>() - satellite + da(-M::CO),
                    Charge::zero(),
                    FragmentType::d(Position::n(sequence_index, sequence_length)),
                ));
            }
        }
        if ions.v {
            base_fragments.push(Fragment::new(
                c_term + da(M::BACKBONE) + da(M::H + M::OH),
                Charge::zero(),
                FragmentType::v(Position::n(sequence_index, sequence_length)),
            ));
        }
        if ions.w {
            for satellite in self.satellite_ion_masses::<M>() {
                base_fragments.push(Fragment::new(
                    c_term + self.mass::<M>() - satellite + da(-M::NH + M::O + M::H),
                    Charge::zero(),
                    FragmentType::w(Position::c(sequence_index, sequence_length)),
                ));
            }
        }
        if ions.x {
            base_fragments.push(Fragment::new(
                c_term + self.mass::<M>() + da(M::CO + M::O),
                Charge::zero(),
                FragmentType::x(Position::c(sequence_index, sequence_length)),
            ));
        }
        if ions.y {
            base_fragments.push(Fragment::new(
                c_term + self.mass::<M>() + da(M::H + M::OH),
                Charge::zero(),
                FragmentType::y(Position::c(sequence_index, sequence_length)),
            ));
        }
        if ions.z {
            base_fragments.push(Fragment::new(
                c_term + self.mass::<M>() + da(-M::NH + M::O),
                Charge::zero(),
                FragmentType::z(Position::c(sequence_index, sequence_length)),
            ));
            base_fragments.push(Fragment::new(
                c_term + self.mass::<M>() + da(-M::NH + M::O + M::H),
                Charge::zero(),
                FragmentType::z·(Position::c(sequence_index, sequence_length)),
            ));
        }
        let mut charged = Vec::with_capacity(base_fragments.len() * max_charge.value as usize);
        for base in base_fragments {
            for charge in 1..=(max_charge.value as u64) {
                charged.push(base.with_charge::<M>(Charge::new::<e>(charge as f64)));
            }
        }
        charged
    }

    pub const fn char(&self) -> char {
        match self {
            Self::Alanine => 'A',
            Self::Cysteine => 'C',
            Self::AsparticAcid => 'D',
            Self::GlutamicAcid => 'E',
            Self::Phenylalanine => 'F',
            Self::Glycine => 'G',
            Self::Histidine => 'H',
            Self::Isoleucine => 'I',
            Self::AmbiguousLeucine => 'J',
            Self::Lysine => 'K',
            Self::Leucine => 'L',
            Self::Methionine => 'M',
            Self::Asparagine => 'N',
            Self::Pyrrolysine => 'O',
            Self::Proline => 'P',
            Self::Glutamine => 'Q',
            Self::Arginine => 'R',
            Self::Serine => 'S',
            Self::Threonine => 'T',
            Self::Selenocysteine => 'U',
            Self::Valine => 'V',
            Self::Tryptophan => 'W',
            Self::Tyrosine => 'Y',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass() {
        let weight_ala = AminoAcid::A.mass::<crate::AverageWeight>();
        let mass_ala = AminoAcid::Ala.mass::<crate::MonoIsotopic>();
        assert_ne!(weight_ala, mass_ala);
        assert_eq!(weight_ala.value, 71.07793);
        assert_eq!(mass_ala.value, 71.037113783);
    }

    #[test]
    fn mass_lysine() {
        let weight_lys = AminoAcid::K.mass::<crate::AverageWeight>();
        let mass_lys = AminoAcid::Lys.mass::<crate::MonoIsotopic>();
        assert_ne!(weight_lys, mass_lys);
        assert_eq!(weight_lys.value, 128.17240999999999);
        assert_eq!(mass_lys.value, 128.094963010536);
    }
}
