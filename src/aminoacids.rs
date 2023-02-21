use crate::fragment::{Fragment, FragmentType};
use crate::model::*;
use crate::std_masses;
use crate::system::f64::*;
use crate::system::mass::dalton;

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
            AminoAcid::Alanine => Mass::new::<dalton>(71.0779),
            AminoAcid::AmbiguousLeucine => Mass::new::<dalton>(113.15764),
            AminoAcid::Arginine => Mass::new::<dalton>(156.18568),
            AminoAcid::Asparagine => Mass::new::<dalton>(114.10264),
            AminoAcid::AsparticAcid => Mass::new::<dalton>(115.0874),
            AminoAcid::Cysteine => Mass::new::<dalton>(103.1429),
            AminoAcid::GlutamicAcid => Mass::new::<dalton>(129.11398),
            AminoAcid::Glutamine => Mass::new::<dalton>(128.12922),
            AminoAcid::Glycine => Mass::new::<dalton>(57.05132),
            AminoAcid::Histidine => Mass::new::<dalton>(137.13928),
            AminoAcid::Isoleucine => Mass::new::<dalton>(113.15764),
            AminoAcid::Leucine => Mass::new::<dalton>(113.15764),
            AminoAcid::Lysine => Mass::new::<dalton>(128.17228),
            AminoAcid::Methionine => Mass::new::<dalton>(131.19606),
            AminoAcid::Phenylalanine => Mass::new::<dalton>(147.17386),
            AminoAcid::Proline => Mass::new::<dalton>(97.11518),
            AminoAcid::Pyrrolysine => Mass::new::<dalton>(237.29816),
            AminoAcid::Selenocysteine => Mass::new::<dalton>(150.3079),
            AminoAcid::Serine => Mass::new::<dalton>(87.0773),
            AminoAcid::Threonine => Mass::new::<dalton>(101.10388),
            AminoAcid::Tryptophan => Mass::new::<dalton>(186.2099),
            AminoAcid::Tyrosine => Mass::new::<dalton>(163.17326),
            AminoAcid::Valine => Mass::new::<dalton>(99.13106),
        }
    }

    ///Data source: https://proteomicsresource.washington.edu/protocols06/masses.php
    pub fn monoisotopic_mass(&self) -> Mass {
        match self {
            AminoAcid::Alanine => Mass::new::<dalton>(71.03711381),
            AminoAcid::AmbiguousLeucine => Mass::new::<dalton>(113.084064),
            AminoAcid::Arginine => Mass::new::<dalton>(156.1011111),
            AminoAcid::Asparagine => Mass::new::<dalton>(114.0429275),
            AminoAcid::AsparticAcid => Mass::new::<dalton>(115.0269431),
            AminoAcid::Cysteine => Mass::new::<dalton>(103.0091845),
            AminoAcid::GlutamicAcid => Mass::new::<dalton>(129.0425931),
            AminoAcid::Glutamine => Mass::new::<dalton>(128.0585775),
            AminoAcid::Glycine => Mass::new::<dalton>(57.02146374),
            AminoAcid::Histidine => Mass::new::<dalton>(137.0589119),
            AminoAcid::Isoleucine => Mass::new::<dalton>(113.084064),
            AminoAcid::Leucine => Mass::new::<dalton>(113.084064),
            AminoAcid::Lysine => Mass::new::<dalton>(128.0949631),
            AminoAcid::Methionine => Mass::new::<dalton>(131.0404846),
            AminoAcid::Phenylalanine => Mass::new::<dalton>(147.0684139),
            AminoAcid::Proline => Mass::new::<dalton>(97.05276388),
            AminoAcid::Pyrrolysine => Mass::new::<dalton>(237.1477269),
            AminoAcid::Selenocysteine => Mass::new::<dalton>(150.9536334),
            AminoAcid::Serine => Mass::new::<dalton>(87.03202844),
            AminoAcid::Threonine => Mass::new::<dalton>(101.0476785),
            AminoAcid::Tryptophan => Mass::new::<dalton>(186.079313),
            AminoAcid::Tyrosine => Mass::new::<dalton>(163.0633286),
            AminoAcid::Valine => Mass::new::<dalton>(99.06841395),
        }
    }

    pub fn satellite_ion_mass(&self) -> Option<Mass> {
        unimplemented!();
        //match self {
        //    AminoAcid::Alanine => Mass::new::<dalton>(71.03711381),
        //    AminoAcid::AmbiguousLeucine => Mass::new::<dalton>(113.084064),
        //    AminoAcid::Arginine => Mass::new::<dalton>(156.1011111),
        //    AminoAcid::Asparagine => Mass::new::<dalton>(114.0429275),
        //    AminoAcid::AsparticAcid => Mass::new::<dalton>(115.0269431),
        //    AminoAcid::Cysteine => Mass::new::<dalton>(103.0091845),
        //    AminoAcid::GlutamicAcid => Mass::new::<dalton>(129.0425931),
        //    AminoAcid::Glutamine => Mass::new::<dalton>(128.0585775),
        //    AminoAcid::Glycine => Mass::new::<dalton>(57.02146374),
        //    AminoAcid::Histidine => Mass::new::<dalton>(137.0589119),
        //    AminoAcid::Isoleucine => Mass::new::<dalton>(113.084064),
        //    AminoAcid::Leucine => Mass::new::<dalton>(113.084064),
        //    AminoAcid::Lysine => Mass::new::<dalton>(128.0949631),
        //    AminoAcid::Methionine => Mass::new::<dalton>(131.0404846),
        //    AminoAcid::Phenylalanine => Mass::new::<dalton>(147.0684139),
        //    AminoAcid::Proline => Mass::new::<dalton>(97.05276388),
        //    AminoAcid::Pyrrolysine => Mass::new::<dalton>(237.1477269),
        //    AminoAcid::Selenocysteine => Mass::new::<dalton>(150.9536334),
        //    AminoAcid::Serine => Mass::new::<dalton>(87.03202844),
        //    AminoAcid::Threonine => Mass::new::<dalton>(101.0476785),
        //    AminoAcid::Tryptophan => Mass::new::<dalton>(186.079313),
        //    AminoAcid::Tyrosine => Mass::new::<dalton>(163.0633286),
        //    AminoAcid::Valine => Mass::new::<dalton>(99.06841395),
        //}
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
                n_term + self.avg_mass() - da(std_masses::CO) + da(std_masses::H * charge.value),
                charge,
                idx,
                FragmentType::a,
            ))
        }
        if ions.b {
            output.push(Fragment::new(
                n_term + self.avg_mass() + da(std_masses::H * charge.value),
                charge,
                idx,
                FragmentType::b,
            ))
        }
        if ions.c {
            output.push(Fragment::new(
                n_term + self.avg_mass() + da(std_masses::NH3) + da(std_masses::H * charge.value),
                charge,
                idx,
                FragmentType::c,
            ))
        }
        if ions.x {
            output.push(Fragment::new(
                c_term
                    + self.avg_mass()
                    + da(std_masses::CO)
                    + da(std_masses::O)
                    + da(std_masses::H * charge.value),
                charge,
                idx,
                FragmentType::x,
            ))
        }
        if ions.y {
            output.push(Fragment::new(
                c_term
                    + self.avg_mass()
                    + da(std_masses::H * (charge.value + 1.0))
                    + da(std_masses::OH),
                charge,
                idx,
                FragmentType::y,
            ))
        }
        if ions.z {
            output.push(Fragment::new(
                c_term + self.avg_mass() - da(std_masses::NH)
                    + da(std_masses::O)
                    + da(std_masses::H * charge.value),
                charge,
                idx,
                FragmentType::z,
            ))
        }
        output
    }
}
