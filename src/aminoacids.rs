use crate::atomic_weights::*;
use crate::fragment::{Fragment, FragmentType};
use crate::model::*;
use crate::system::f64::*;

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

impl TryFrom<u8> for AminoAcid {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
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

impl AminoAcid {
    ///Data source: https://proteomicsresource.washington.edu/protocols06/masses.php
    pub fn avg_mass(&self) -> Mass {
        match self {
            AminoAcid::Alanine => da(BACKBONE + CH3),
            AminoAcid::AmbiguousLeucine => da(BACKBONE + C * 4.0 + H * 9.0),
            AminoAcid::Arginine => da(BACKBONE + CH2 * 3.0 + NH + NH2 * 2.0),
            AminoAcid::Asparagine => da(BACKBONE + CH2 + C + O + NH2),
            AminoAcid::AsparticAcid => da(BACKBONE + CH2 + C + O * 2.0),
            AminoAcid::Cysteine => da(BACKBONE + CH2 + S + H),
            AminoAcid::GlutamicAcid => da(BACKBONE + CH2 * 2.0 + C + O * 2.0),
            AminoAcid::Glutamine => da(BACKBONE + CH2 * 2.0 + C + O + NH2),
            AminoAcid::Glycine => da(BACKBONE + H),
            AminoAcid::Histidine => da(BACKBONE + CH2 + C + N + CH + NH + CH),
            AminoAcid::Isoleucine => da(BACKBONE + CH + CH3 + CH2 + CH3),
            AminoAcid::Leucine => da(BACKBONE + CH2 + CH + CH3 + CH3),
            AminoAcid::Lysine => da(BACKBONE + CH2 * 4.0 + NH2),
            AminoAcid::Methionine => da(BACKBONE + CH2 * 2.0 + S + CH3),
            AminoAcid::Phenylalanine => da(BACKBONE + CH2 + C + CH * 5.0),
            AminoAcid::Proline => da(BACKBONE + CH2 * 3.0 - H),
            AminoAcid::Pyrrolysine => {
                da(BACKBONE + CH2 * 4.0 + NH + C + O + CH + N + CH + CH2 + CH + CH3)
            }
            AminoAcid::Selenocysteine => da(BACKBONE + CH2 + Se),
            AminoAcid::Serine => da(BACKBONE + CH2 + OH),
            AminoAcid::Threonine => da(BACKBONE + CH + OH + CH3),
            AminoAcid::Tryptophan => da(BACKBONE + CH2 + C * 2.0 + CH * 4.0 + C + NH + CH),
            AminoAcid::Tyrosine => da(BACKBONE + CH2 + C * 2.0 + CH * 4.0 + OH),
            AminoAcid::Valine => da(BACKBONE + CH + CH3 * 2.0),
        }
    }

    // TODO: Take side chain mutations into account (maybe define pyrrolysine as a mutation)
    pub fn satellite_ion_masses(&self) -> Vec<Mass> {
        match self {
            AminoAcid::Alanine => vec![],
            AminoAcid::AmbiguousLeucine => vec![],
            AminoAcid::Arginine => vec![da(CH2 * 2.0 + NH + NH2 * 2.0)],
            AminoAcid::Asparagine => vec![da(C + O + NH2)],
            AminoAcid::AsparticAcid => vec![da(C + O * 2.0)],
            AminoAcid::Cysteine => vec![da(S + H)],
            AminoAcid::GlutamicAcid => vec![da(CH2 + C + O * 2.0)],
            AminoAcid::Glutamine => vec![da(CH2 + C + O + NH2)],
            AminoAcid::Glycine => vec![],
            AminoAcid::Histidine => vec![], // Aromatic
            AminoAcid::Isoleucine => vec![da(CH3), da(CH2 + CH3)],
            AminoAcid::Leucine => vec![da(CH + CH3 * 2.0)],
            AminoAcid::Lysine => vec![da(CH2 * 3.0 + NH2)],
            AminoAcid::Methionine => vec![da(CH2 + S + CH3)],
            AminoAcid::Phenylalanine => vec![], // Aromatic
            AminoAcid::Proline => vec![], // Interesting, TODO: see what other software packages think about this one
            AminoAcid::Pyrrolysine => {
                vec![da(CH2 * 3.0 + NH + C + O + CH + N + CH + CH2 + CH + CH3)]
            } // Weird, TODO: figure out what to make of this
            AminoAcid::Selenocysteine => vec![da(Se)],
            AminoAcid::Serine => vec![da(OH)],
            AminoAcid::Threonine => vec![da(OH), da(CH3)],
            AminoAcid::Tryptophan => vec![],    // Aromatic
            AminoAcid::Tyrosine => vec![],      // Aromatic
            AminoAcid::Valine => vec![da(CH3)], // Technically two options, but both have the same mass TODO: check if the loss of both is an option
        }
    }

    pub fn fragments(
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
                n_term + self.avg_mass() - da(CO) + da(H * charge.value),
                charge,
                idx,
                FragmentType::a,
            ))
        }
        if ions.b {
            output.push(Fragment::new(
                n_term + self.avg_mass() + da(H * charge.value),
                charge,
                idx,
                FragmentType::b,
            ))
        }
        if ions.c {
            output.push(Fragment::new(
                n_term + self.avg_mass() + da(NH3) + da(H * charge.value),
                charge,
                idx,
                FragmentType::c,
            ))
        }
        if ions.d {
            for satellite in self.satellite_ion_masses() {
                output.push(Fragment::new(
                    n_term + self.avg_mass() - satellite - da(CO) + da(H * charge.value),
                    charge,
                    idx,
                    FragmentType::d,
                ))
            }
        }
        if ions.v {
            output.push(Fragment::new(
                c_term + da(BACKBONE) + da(H * (charge.value + 1.0)) + da(OH),
                charge,
                idx,
                FragmentType::v,
            ))
        }
        if ions.w {
            for satellite in self.satellite_ion_masses() {
                output.push(Fragment::new(
                    c_term + self.avg_mass() - satellite - da(NH) + da(O) + da(H * charge.value),
                    charge,
                    idx,
                    FragmentType::w,
                ))
            }
        }
        if ions.x {
            output.push(Fragment::new(
                c_term + self.avg_mass() + da(CO) + da(O) + da(H * charge.value),
                charge,
                idx,
                FragmentType::x,
            ))
        }
        if ions.y {
            output.push(Fragment::new(
                c_term + self.avg_mass() + da(H * (charge.value + 1.0)) + da(OH),
                charge,
                idx,
                FragmentType::y,
            ))
        }
        if ions.z {
            output.push(Fragment::new(
                c_term + self.avg_mass() - da(NH) + da(O) + da(H * charge.value),
                charge,
                idx,
                FragmentType::z,
            ))
        }
        output
    }
}
