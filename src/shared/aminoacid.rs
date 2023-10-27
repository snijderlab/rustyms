/// An amino acid, alongside the standard ones some ambiguous (J/X) and non-standard (U/O) are included.
/// <https://www.insdc.org/submitting-standards/feature-table/#7.4.3>
#[derive(Clone, Debug, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
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
            b'A' | b'a' => Ok(Self::Alanine),
            b'B' | b'b' => Ok(Self::AmbiguousAsparagine),
            b'C' | b'c' => Ok(Self::Cysteine),
            b'D' | b'd' => Ok(Self::AsparticAcid),
            b'E' | b'e' => Ok(Self::GlutamicAcid),
            b'F' | b'f' => Ok(Self::Phenylalanine),
            b'G' | b'g' => Ok(Self::Glycine),
            b'H' | b'h' => Ok(Self::Histidine),
            b'I' | b'i' => Ok(Self::Isoleucine),
            b'J' | b'j' => Ok(Self::AmbiguousLeucine),
            b'K' | b'k' => Ok(Self::Lysine),
            b'L' | b'l' => Ok(Self::Leucine),
            b'M' | b'm' => Ok(Self::Methionine),
            b'N' | b'n' => Ok(Self::Asparagine),
            b'O' | b'o' => Ok(Self::Pyrrolysine),
            b'P' | b'p' => Ok(Self::Proline),
            b'Q' | b'q' => Ok(Self::Glutamine),
            b'R' | b'r' => Ok(Self::Arginine),
            b'S' | b's' => Ok(Self::Serine),
            b'T' | b't' => Ok(Self::Threonine),
            b'U' | b'u' => Ok(Self::Selenocysteine),
            b'V' | b'v' => Ok(Self::Valine),
            b'W' | b'w' => Ok(Self::Tryptophan),
            b'X' | b'x' => Ok(Self::Unknown),
            b'Y' | b'y' => Ok(Self::Tyrosine),
            b'Z' | b'z' => Ok(Self::AmbiguousGlutamine),
            _ => Err(()),
        }
    }
}
