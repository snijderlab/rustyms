use crate::fragment::{Fragment, FragmentType};
use crate::model::*;
use crate::system::f64::*;
use crate::MassSystem;

#[derive(Clone, Debug)]
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
            'A' => Ok(AminoAcid::Alanine),
            'C' => Ok(AminoAcid::Cysteine),
            'D' => Ok(AminoAcid::AsparticAcid),
            'E' => Ok(AminoAcid::GlutamicAcid),
            'F' => Ok(AminoAcid::Phenylalanine),
            'G' => Ok(AminoAcid::Glycine),
            'H' => Ok(AminoAcid::Histidine),
            'I' => Ok(AminoAcid::Isoleucine),
            'J' => Ok(AminoAcid::AmbiguousLeucine),
            'K' => Ok(AminoAcid::Lysine),
            'L' => Ok(AminoAcid::Leucine),
            'M' => Ok(AminoAcid::Methionine),
            'N' => Ok(AminoAcid::Asparagine),
            'O' => Ok(AminoAcid::Pyrrolysine),
            'P' => Ok(AminoAcid::Proline),
            'Q' => Ok(AminoAcid::Glutamine),
            'R' => Ok(AminoAcid::Arginine),
            'S' => Ok(AminoAcid::Serine),
            'T' => Ok(AminoAcid::Threonine),
            'U' => Ok(AminoAcid::Selenocysteine),
            'V' => Ok(AminoAcid::Valine),
            'W' => Ok(AminoAcid::Tryptophan),
            'Y' => Ok(AminoAcid::Tyrosine),
            _ => Err(()),
        }
    }
}

impl TryFrom<&u8> for AminoAcid {
    type Error = ();
    fn try_from(value: &u8) -> Result<Self, Self::Error> {
        match value {
            b'A' => Ok(AminoAcid::Alanine),
            b'C' => Ok(AminoAcid::Cysteine),
            b'D' => Ok(AminoAcid::AsparticAcid),
            b'E' => Ok(AminoAcid::GlutamicAcid),
            b'F' => Ok(AminoAcid::Phenylalanine),
            b'G' => Ok(AminoAcid::Glycine),
            b'H' => Ok(AminoAcid::Histidine),
            b'I' => Ok(AminoAcid::Isoleucine),
            b'J' => Ok(AminoAcid::AmbiguousLeucine),
            b'K' => Ok(AminoAcid::Lysine),
            b'L' => Ok(AminoAcid::Leucine),
            b'M' => Ok(AminoAcid::Methionine),
            b'N' => Ok(AminoAcid::Asparagine),
            b'O' => Ok(AminoAcid::Pyrrolysine),
            b'P' => Ok(AminoAcid::Proline),
            b'Q' => Ok(AminoAcid::Glutamine),
            b'R' => Ok(AminoAcid::Arginine),
            b'S' => Ok(AminoAcid::Serine),
            b'T' => Ok(AminoAcid::Threonine),
            b'U' => Ok(AminoAcid::Selenocysteine),
            b'V' => Ok(AminoAcid::Valine),
            b'W' => Ok(AminoAcid::Tryptophan),
            b'Y' => Ok(AminoAcid::Tyrosine),
            _ => Err(()),
        }
    }
}

#[allow(non_upper_case_globals)]
impl AminoAcid {
    pub const A: AminoAcid = AminoAcid::Alanine;
    pub const C: AminoAcid = AminoAcid::Cysteine;
    pub const D: AminoAcid = AminoAcid::AsparticAcid;
    pub const E: AminoAcid = AminoAcid::GlutamicAcid;
    pub const F: AminoAcid = AminoAcid::Phenylalanine;
    pub const G: AminoAcid = AminoAcid::Glycine;
    pub const H: AminoAcid = AminoAcid::Histidine;
    pub const I: AminoAcid = AminoAcid::Isoleucine;
    pub const J: AminoAcid = AminoAcid::AmbiguousLeucine;
    pub const K: AminoAcid = AminoAcid::Lysine;
    pub const L: AminoAcid = AminoAcid::Leucine;
    pub const M: AminoAcid = AminoAcid::Methionine;
    pub const N: AminoAcid = AminoAcid::Asparagine;
    pub const O: AminoAcid = AminoAcid::Pyrrolysine;
    pub const P: AminoAcid = AminoAcid::Proline;
    pub const Q: AminoAcid = AminoAcid::Glutamine;
    pub const R: AminoAcid = AminoAcid::Arginine;
    pub const S: AminoAcid = AminoAcid::Serine;
    pub const T: AminoAcid = AminoAcid::Threonine;
    pub const U: AminoAcid = AminoAcid::Selenocysteine;
    pub const V: AminoAcid = AminoAcid::Valine;
    pub const W: AminoAcid = AminoAcid::Tryptophan;
    pub const Ala: AminoAcid = AminoAcid::Alanine;
    pub const Cys: AminoAcid = AminoAcid::Cysteine;
    pub const Asp: AminoAcid = AminoAcid::AsparticAcid;
    pub const Glu: AminoAcid = AminoAcid::GlutamicAcid;
    pub const Phe: AminoAcid = AminoAcid::Phenylalanine;
    pub const Gly: AminoAcid = AminoAcid::Glycine;
    pub const His: AminoAcid = AminoAcid::Histidine;
    pub const Ile: AminoAcid = AminoAcid::Isoleucine;
    pub const Xle: AminoAcid = AminoAcid::AmbiguousLeucine;
    pub const Lys: AminoAcid = AminoAcid::Lysine;
    pub const Leu: AminoAcid = AminoAcid::Leucine;
    pub const Met: AminoAcid = AminoAcid::Methionine;
    pub const Asn: AminoAcid = AminoAcid::Asparagine;
    pub const Pyl: AminoAcid = AminoAcid::Pyrrolysine;
    pub const Pro: AminoAcid = AminoAcid::Proline;
    pub const Gln: AminoAcid = AminoAcid::Glutamine;
    pub const Arg: AminoAcid = AminoAcid::Arginine;
    pub const Ser: AminoAcid = AminoAcid::Serine;
    pub const Thr: AminoAcid = AminoAcid::Threonine;
    pub const Sec: AminoAcid = AminoAcid::Selenocysteine;
    pub const Val: AminoAcid = AminoAcid::Valine;
    pub const Trp: AminoAcid = AminoAcid::Tryptophan;

    pub fn mass<M: MassSystem>(&self) -> Mass {
        match self {
            AminoAcid::Alanine => da(M::BACKBONE + M::CH3),
            AminoAcid::AmbiguousLeucine => da(M::BACKBONE + M::C * 4.0 + M::H * 9.0),
            AminoAcid::Arginine => da(M::BACKBONE + M::CH2 * 3.0 + M::NH + M::NH2 * 2.0),
            AminoAcid::Asparagine => da(M::BACKBONE + M::CH2 + M::C + M::O + M::NH2),
            AminoAcid::AsparticAcid => da(M::BACKBONE + M::CH2 + M::C + M::O * 2.0),
            AminoAcid::Cysteine => da(M::BACKBONE + M::CH2 + M::S + M::H),
            AminoAcid::GlutamicAcid => da(M::BACKBONE + M::CH2 * 2.0 + M::C + M::O * 2.0),
            AminoAcid::Glutamine => da(M::BACKBONE + M::CH2 * 2.0 + M::C + M::O + M::NH2),
            AminoAcid::Glycine => da(M::BACKBONE + M::H),
            AminoAcid::Histidine => da(M::BACKBONE + M::CH2 + M::C + M::N + M::CH + M::NH + M::CH),
            AminoAcid::Isoleucine => da(M::BACKBONE + M::CH + M::CH3 + M::CH2 + M::CH3),
            AminoAcid::Leucine => da(M::BACKBONE + M::CH2 + M::CH + M::CH3 + M::CH3),
            AminoAcid::Lysine => da(M::BACKBONE + M::CH2 * 4.0 + M::NH2),
            AminoAcid::Methionine => da(M::BACKBONE + M::CH2 * 2.0 + M::S + M::CH3),
            AminoAcid::Phenylalanine => da(M::BACKBONE + M::CH2 + M::C + M::CH * 5.0),
            AminoAcid::Proline => da(M::BACKBONE + M::CH2 * 3.0 - M::H),
            AminoAcid::Pyrrolysine => da(M::BACKBONE
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
            AminoAcid::Selenocysteine => da(M::BACKBONE + M::CH2 + M::Se),
            AminoAcid::Serine => da(M::BACKBONE + M::CH2 + M::OH),
            AminoAcid::Threonine => da(M::BACKBONE + M::CH + M::OH + M::CH3),
            AminoAcid::Tryptophan => {
                da(M::BACKBONE + M::CH2 + M::C * 2.0 + M::CH * 4.0 + M::C + M::NH + M::CH)
            }
            AminoAcid::Tyrosine => da(M::BACKBONE + M::CH2 + M::C * 2.0 + M::CH * 4.0 + M::OH),
            AminoAcid::Valine => da(M::BACKBONE + M::CH + M::CH3 * 2.0),
        }
    }

    // TODO: Take side chain mutations into account (maybe define pyrrolysine as a mutation)
    pub fn satellite_ion_masses<M: MassSystem>(&self) -> Vec<Mass> {
        match self {
            AminoAcid::Alanine => vec![],
            AminoAcid::AmbiguousLeucine => vec![],
            AminoAcid::Arginine => vec![da(M::CH2 * 2.0 + M::NH + M::NH2 * 2.0)],
            AminoAcid::Asparagine => vec![da(M::C + M::O + M::NH2)],
            AminoAcid::AsparticAcid => vec![da(M::C + M::O * 2.0)],
            AminoAcid::Cysteine => vec![da(M::S + M::H)],
            AminoAcid::GlutamicAcid => vec![da(M::CH2 + M::C + M::O * 2.0)],
            AminoAcid::Glutamine => vec![da(M::CH2 + M::C + M::O + M::NH2)],
            AminoAcid::Glycine => vec![],
            AminoAcid::Histidine => vec![], // Aromatic
            AminoAcid::Isoleucine => vec![da(M::CH3), da(M::CH2 + M::CH3)],
            AminoAcid::Leucine => vec![da(M::CH + M::CH3 * 2.0)],
            AminoAcid::Lysine => vec![da(M::CH2 * 3.0 + M::NH2)],
            AminoAcid::Methionine => vec![da(M::CH2 + M::S + M::CH3)],
            AminoAcid::Phenylalanine => vec![], // Aromatic
            AminoAcid::Proline => vec![], // Interesting, TODO: see what other software packages think about this one
            AminoAcid::Pyrrolysine => {
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
            AminoAcid::Selenocysteine => vec![da(M::Se)],
            AminoAcid::Serine => vec![da(M::OH)],
            AminoAcid::Threonine => vec![da(M::OH), da(M::CH3)],
            AminoAcid::Tryptophan => vec![],       // Aromatic
            AminoAcid::Tyrosine => vec![],         // Aromatic
            AminoAcid::Valine => vec![da(M::CH3)], // Technically two options, but both have the same mass TODO: check if the loss of both is an option
        }
    }

    pub fn fragments<M: MassSystem>(
        &self,
        n_term: Mass,
        c_term: Mass,
        charge: Charge,
        idx: usize,
        ions: PossibleIons,
    ) -> Vec<Fragment> {
        let mut output = Vec::new();
        if ions.a {
            output.push(Fragment::new(
                n_term + self.mass::<M>() - da(M::CO + M::H * charge.value),
                charge,
                idx,
                FragmentType::a,
            ))
        }
        if ions.b {
            output.push(Fragment::new(
                n_term + self.mass::<M>() + da(M::H * charge.value),
                charge,
                idx,
                FragmentType::b,
            ))
        }
        if ions.c {
            output.push(Fragment::new(
                n_term + self.mass::<M>() + da(M::NH3 + M::H * charge.value),
                charge,
                idx,
                FragmentType::c,
            ))
        }
        if ions.d {
            for satellite in self.satellite_ion_masses::<M>() {
                output.push(Fragment::new(
                    n_term + self.mass::<M>() - satellite + da(-M::CO + M::H * charge.value),
                    charge,
                    idx,
                    FragmentType::d,
                ))
            }
        }
        if ions.v {
            output.push(Fragment::new(
                c_term + da(M::BACKBONE) + da(M::H * (charge.value + 1.0) + M::OH),
                charge,
                idx,
                FragmentType::v,
            ))
        }
        if ions.w {
            for satellite in self.satellite_ion_masses::<M>() {
                output.push(Fragment::new(
                    c_term + self.mass::<M>() - satellite + da(-M::NH + M::O + M::H * charge.value),
                    charge,
                    idx,
                    FragmentType::w,
                ))
            }
        }
        if ions.x {
            output.push(Fragment::new(
                c_term + self.mass::<M>() + da(M::CO + M::O + M::H * charge.value),
                charge,
                idx,
                FragmentType::x,
            ))
        }
        if ions.y {
            output.push(Fragment::new(
                c_term + self.mass::<M>() + da(M::H * (charge.value + 1.0) + M::OH),
                charge,
                idx,
                FragmentType::y,
            ))
        }
        if ions.z {
            output.push(Fragment::new(
                c_term + self.mass::<M>() + da(-M::NH + M::O + M::H * charge.value),
                charge,
                idx,
                FragmentType::z,
            ))
        }
        output
    }

    pub fn char(&self) -> char {
        match self {
            AminoAcid::Alanine => 'A',
            AminoAcid::Cysteine => 'C',
            AminoAcid::AsparticAcid => 'D',
            AminoAcid::GlutamicAcid => 'E',
            AminoAcid::Phenylalanine => 'F',
            AminoAcid::Glycine => 'G',
            AminoAcid::Histidine => 'H',
            AminoAcid::Isoleucine => 'I',
            AminoAcid::AmbiguousLeucine => 'J',
            AminoAcid::Lysine => 'K',
            AminoAcid::Leucine => 'L',
            AminoAcid::Methionine => 'M',
            AminoAcid::Asparagine => 'N',
            AminoAcid::Pyrrolysine => 'O',
            AminoAcid::Proline => 'P',
            AminoAcid::Glutamine => 'Q',
            AminoAcid::Arginine => 'R',
            AminoAcid::Serine => 'S',
            AminoAcid::Threonine => 'T',
            AminoAcid::Selenocysteine => 'U',
            AminoAcid::Valine => 'V',
            AminoAcid::Tryptophan => 'W',
            AminoAcid::Tyrosine => 'Y',
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
}
