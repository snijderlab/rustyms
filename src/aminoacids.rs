use crate::formula::{Chemical, MolecularFormula};
use crate::fragment::{Fragment, FragmentType};
use crate::model::*;
use crate::molecular_charge::MolecularCharge;
use crate::Element;
use crate::Position;

/// An amino acid, alongside the standard ones some ambiguous (J/X) and non-standard (U/O) are included.
/// <https://www.insdc.org/submitting-standards/feature-table/#7.4.3>
#[derive(Clone, Debug, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[allow(missing_docs)]
pub enum AminoAcid {
    Alanine = 0,
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
    Unknown,
    AmbiguousAsparagine,
    AmbiguousGlutamine,
}

impl TryFrom<&str> for AminoAcid {
    type Error = ();
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        if value.is_ascii() && value.len() == 1 {
            let ch = value.chars().next().unwrap();
            ch.try_into()
        } else {
            Err(())
        }
    }
}

impl TryFrom<char> for AminoAcid {
    type Error = ();
    fn try_from(value: char) -> Result<Self, Self::Error> {
        if value.is_ascii() {
            let num = value as u8;
            num.try_into()
        } else {
            Err(())
        }
    }
}

impl TryFrom<u8> for AminoAcid {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'A' => Ok(Self::Alanine),
            b'B' => Ok(Self::AmbiguousAsparagine),
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
            b'X' => Ok(Self::Unknown),
            b'Y' => Ok(Self::Tyrosine),
            b'Z' => Ok(Self::AmbiguousGlutamine),
            _ => Err(()),
        }
    }
}

impl Chemical for AminoAcid {
    fn formula(&self) -> MolecularFormula {
        match self {
            Self::Alanine => molecular_formula!(H 5 C 3 O 1 N 1),
            Self::Arginine => molecular_formula!(H 12 C 6 O 1 N 4), // One of the H's counts as the charge carrier and is added later
            Self::AmbiguousAsparagine => {
                panic!("Molecular formula of Asx/B ambiguous asparagine is not defined.")
            }
            Self::Asparagine => molecular_formula!(H 6 C 4 O 2 N 2),
            Self::AsparticAcid => molecular_formula!(H 5 C 4 O 3 N 1),
            Self::Cysteine => molecular_formula!(H 5 C 3 O 1 N 1 S 1),
            Self::GlutamicAcid => molecular_formula!(H 7 C 5 O 3 N 1),
            Self::Glutamine => molecular_formula!(H 8 C 5 O 2 N 2),
            Self::AmbiguousGlutamine => {
                panic!("Molecular formula of Glx/Z ambiguous glutamine is not defined.")
            }
            Self::Glycine => molecular_formula!(H 3 C 2 O 1 N 1),
            Self::Histidine => molecular_formula!(H 7 C 6 O 1 N 3),
            Self::AmbiguousLeucine | Self::Isoleucine | Self::Leucine => {
                molecular_formula!(H 11 C 6 O 1 N 1)
            }
            Self::Lysine => molecular_formula!(H 12 C 6 O 1 N 2),
            Self::Methionine => molecular_formula!(H 9 C 5 O 1 N 1 S 1),
            Self::Phenylalanine => molecular_formula!(H 9 C 9 O 1 N 1),
            Self::Proline => molecular_formula!(H 7 C 5 O 1 N 1),
            Self::Pyrrolysine => molecular_formula!(H 19 C 11 O 2 N 3),
            Self::Selenocysteine => molecular_formula!(H 4 C 3 O 1 N 1 Se 1),
            Self::Serine => molecular_formula!(H 5 C 3 O 2 N 1),
            Self::Threonine => molecular_formula!(H 7 C 4 O 2 N 1),
            Self::Tryptophan => molecular_formula!(H 10 C 11 O 1 N 2),
            Self::Tyrosine => molecular_formula!(H 9 C 9 O 2 N 1),
            Self::Valine => molecular_formula!(H 9 C 5 O 1 N 1),
            Self::Unknown => molecular_formula!(),
        }
    }
}

#[allow(non_upper_case_globals, missing_docs)]
impl AminoAcid {
    pub const A: Self = Self::Alanine;
    pub const B: Self = Self::AmbiguousAsparagine;
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
    pub const X: Self = Self::Unknown;
    pub const Y: Self = Self::Tyrosine;
    pub const Z: Self = Self::AmbiguousGlutamine;
    pub const Ala: Self = Self::Alanine;
    pub const Cys: Self = Self::Cysteine;
    pub const Asn: Self = Self::Asparagine;
    pub const Asp: Self = Self::AsparticAcid;
    pub const Asx: Self = Self::AmbiguousAsparagine;
    pub const Glu: Self = Self::GlutamicAcid;
    pub const Phe: Self = Self::Phenylalanine;
    pub const Gly: Self = Self::Glycine;
    pub const His: Self = Self::Histidine;
    pub const Ile: Self = Self::Isoleucine;
    pub const Xle: Self = Self::AmbiguousLeucine;
    pub const Lys: Self = Self::Lysine;
    pub const Leu: Self = Self::Leucine;
    pub const Met: Self = Self::Methionine;
    pub const Pyl: Self = Self::Pyrrolysine;
    pub const Pro: Self = Self::Proline;
    pub const Gln: Self = Self::Glutamine;
    pub const Glx: Self = Self::AmbiguousGlutamine;
    pub const Arg: Self = Self::Arginine;
    pub const Ser: Self = Self::Serine;
    pub const Thr: Self = Self::Threonine;
    pub const Sec: Self = Self::Selenocysteine;
    pub const Val: Self = Self::Valine;
    pub const Trp: Self = Self::Tryptophan;
    pub const Tyr: Self = Self::Tyrosine;
    pub const Xaa: Self = Self::Unknown;

    // TODO: Take side chain mutations into account (maybe define pyrrolysine as a mutation)
    pub fn satellite_ion_fragments(&self) -> Vec<MolecularFormula> {
        match self {
            Self::Alanine
            | Self::AmbiguousLeucine
            | Self::Glycine
            | Self::Histidine
            | Self::Phenylalanine
            | Self::Proline
            | Self::Tryptophan
            | Self::Tyrosine
            | Self::AmbiguousAsparagine
            | Self::AmbiguousGlutamine
            | Self::Unknown => vec![],
            Self::Arginine => vec![molecular_formula!(H 9 C 2 N 2)],
            Self::Asparagine => vec![molecular_formula!(H 2 C 1 N 1 O 1)],
            Self::AsparticAcid => vec![molecular_formula!(H 1 C 1 O 2)],
            Self::Cysteine => vec![molecular_formula!(H 1 S 1)],
            Self::GlutamicAcid => vec![molecular_formula!(H 3 C 2 O 2)],
            Self::Glutamine => vec![molecular_formula!(H 4 C 2 N 1 O 1)],
            Self::Isoleucine => vec![molecular_formula!(H 3 C 1), molecular_formula!(H 5 C 2)],
            Self::Leucine => vec![molecular_formula!(H 7 C 3)],
            Self::Lysine => vec![molecular_formula!(H 8 C 3 N 1)],
            Self::Methionine => vec![molecular_formula!(H 5 C 2 S 1)],
            Self::Pyrrolysine => {
                vec![molecular_formula!(H 15 C 9 N 2 O 1)]
            } // Weird, TODO: figure out what to make of this
            Self::Selenocysteine => vec![molecular_formula!(Se 1)],
            Self::Serine => vec![molecular_formula!(H 1 O 1)],
            Self::Threonine => vec![molecular_formula!(H 1 O 1), molecular_formula!(H 3 C 1)],
            Self::Valine => vec![molecular_formula!(H 3 C 1)], // Technically two options, but both have the same mass
        }
    }

    #[allow(clippy::too_many_lines, clippy::too_many_arguments)]
    pub fn fragments(
        &self,
        n_term: &[(MolecularFormula, String)],
        c_term: &[(MolecularFormula, String)],
        modifications: &MolecularFormula,
        charge_carriers: &MolecularCharge,
        sequence_index: usize,
        sequence_length: usize,
        ions: &PossibleIons,
        peptide_index: usize,
    ) -> Vec<Fragment> {
        let mut base_fragments = Vec::with_capacity(ions.size_upper_bound());
        if ions.a.0 {
            base_fragments.extend(Fragment::generate_all(
                &(modifications + &self.formula() - molecular_formula!(H 1 C 1 O 1)),
                peptide_index,
                FragmentType::a(Position::n(sequence_index, sequence_length)),
                n_term,
                ions.a.1,
            ));
        }
        if ions.b.0 {
            base_fragments.extend(Fragment::generate_all(
                &(modifications + &self.formula() - molecular_formula!(H 1)),
                peptide_index,
                FragmentType::b(Position::n(sequence_index, sequence_length)),
                n_term,
                ions.b.1,
            ));
        }
        if ions.c.0 {
            base_fragments.extend(Fragment::generate_all(
                &(modifications + &self.formula() + molecular_formula!(H 2 N 1)),
                peptide_index,
                FragmentType::c(Position::n(sequence_index, sequence_length)),
                n_term,
                ions.c.1,
            ));
        }
        if ions.d.0 {
            for satellite in self.satellite_ion_fragments() {
                base_fragments.extend(Fragment::generate_all(
                    &(modifications + &self.formula()
                        - satellite
                        - molecular_formula!(H 1 C 1 O 1)),
                    peptide_index,
                    FragmentType::d(Position::n(sequence_index, sequence_length)),
                    n_term,
                    ions.d.1,
                ));
            }
        }
        if ions.v.0 {
            base_fragments.extend(Fragment::generate_all(
                &(modifications + &molecular_formula!(H 3 C 2 N 1 O 1)), // TODO: are the modifications needed here?
                peptide_index,
                FragmentType::v(Position::n(sequence_index, sequence_length)),
                c_term,
                ions.v.1,
            ));
        }
        if ions.w.0 {
            for satellite in self.satellite_ion_fragments() {
                base_fragments.extend(Fragment::generate_all(
                    &(modifications + &self.formula() - satellite - molecular_formula!(H 2 N 1)),
                    peptide_index,
                    FragmentType::w(Position::c(sequence_index, sequence_length)),
                    c_term,
                    ions.w.1,
                ));
            }
        }
        if ions.x.0 {
            base_fragments.extend(Fragment::generate_all(
                &(modifications + &self.formula() + molecular_formula!(C 1 O 1)
                    - molecular_formula!(H 1)),
                peptide_index,
                FragmentType::x(Position::c(sequence_index, sequence_length)),
                c_term,
                ions.x.1,
            ));
        }
        if ions.y.0 {
            base_fragments.extend(Fragment::generate_all(
                &(modifications + &self.formula() + molecular_formula!(H 1)),
                peptide_index,
                FragmentType::y(Position::c(sequence_index, sequence_length)),
                c_term,
                ions.y.1,
            ));
        }
        if ions.z.0 {
            base_fragments.extend(Fragment::generate_all(
                &(modifications + &self.formula() - molecular_formula!(H 2 N 1)),
                peptide_index,
                FragmentType::z(Position::c(sequence_index, sequence_length)),
                c_term,
                ions.z.1,
            ));
            base_fragments.extend(Fragment::generate_all(
                &(modifications + &self.formula() - molecular_formula!(H 1 N 1)),
                peptide_index,
                FragmentType::z·(Position::c(sequence_index, sequence_length)),
                c_term,
                ions.z.1,
            ));
        }
        let charge_options = charge_carriers.all_charge_options();
        let mut charged = Vec::with_capacity(base_fragments.len() * charge_options.len());
        for base in base_fragments {
            for charge in &charge_options {
                charged.push(base.with_charge(charge));
            }
        }
        charged
    }

    pub const fn char(&self) -> char {
        match self {
            Self::Alanine => 'A',
            Self::AmbiguousAsparagine => 'B',
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
            Self::Unknown => 'X',
            Self::Tyrosine => 'Y',
            Self::AmbiguousGlutamine => 'Z',
        }
    }
}

#[cfg(test)]
#[allow(clippy::unreadable_literal, clippy::float_cmp)]
mod tests {
    use super::*;

    #[test]
    fn mass() {
        let weight_ala = AminoAcid::A.formula().average_weight().unwrap();
        let mass_ala = AminoAcid::Ala.formula().monoisotopic_mass().unwrap();
        assert_ne!(weight_ala, mass_ala);
        assert!((weight_ala.value - 71.07793).abs() < 1e-5);
        assert!((mass_ala.value - 71.037113783).abs() < 1e-5);
    }

    #[test]
    fn mass_lysine() {
        let weight_lys = AminoAcid::K.formula().average_weight().unwrap();
        let mass_lys = AminoAcid::Lys.formula().monoisotopic_mass().unwrap();
        assert_ne!(weight_lys, mass_lys);
        assert!((weight_lys.value - 128.17240999999999).abs() < 1e-5);
        assert!((mass_lys.value - 128.094963010536).abs() < 1e-5);
    }

    #[test]
    fn masses() {
        let known = &[
            ('A', 71.03711, 71.08),
            ('R', 156.10111, 156.2),
            ('N', 114.04293, 114.1),
            ('D', 115.02694, 115.1),
            ('C', 103.00919, 103.1),
            ('E', 129.04259, 129.1),
            ('Q', 128.05858, 128.1),
            ('G', 57.02146, 57.05),
            ('H', 137.05891, 137.1),
            ('I', 113.08406, 113.2),
            ('L', 113.08406, 113.2),
            ('K', 128.09496, 128.2),
            ('M', 131.04049, 131.2),
            ('F', 147.06841, 147.2),
            ('P', 97.05276, 97.12),
            ('S', 87.03203, 87.08),
            ('T', 101.04768, 101.1),
            ('W', 186.07931, 186.2),
            ('Y', 163.06333, 163.2),
            ('V', 99.06841, 99.13),
        ];

        for (aa, mono_mass, average_weight) in known {
            let aa = AminoAcid::try_from(*aa).unwrap();
            let (mono, weight) = (
                aa.formula().monoisotopic_mass().unwrap().value,
                aa.formula().average_weight().unwrap().value,
            );
            println!(
                "{}: {} {} {} {}",
                aa.char(),
                mono,
                mono_mass,
                weight,
                average_weight
            );
            assert!((mono - *mono_mass).abs() < 1e-5);
            assert!((weight - *average_weight).abs() < 1e-1);
        }
    }

    #[test]
    fn read_aa() {
        assert_eq!(
            AminoAcid::try_from('B').unwrap(),
            AminoAcid::AmbiguousAsparagine
        );
        assert_eq!(
            AminoAcid::try_from(b'B').unwrap(),
            AminoAcid::AmbiguousAsparagine
        );
        assert_eq!(AminoAcid::try_from('c'), Err(()));
        assert_eq!(AminoAcid::try_from('🦀'), Err(()));
    }
}
