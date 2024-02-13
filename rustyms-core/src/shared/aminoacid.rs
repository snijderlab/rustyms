/// An amino acid, alongside the standard ones some ambiguous (B/J/Z/X) and non-standard (U/O) are included.
/// <https://www.insdc.org/submitting-standards/feature-table/#7.4.3>
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
#[allow(missing_docs)]
pub enum AminoAcid {
    #[default]
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
    Lysine,
    Methionine,
    Phenylalanine,
    Proline,
    Serine,
    Threonine,
    Tryptophan,
    Tyrosine,
    Valine,
    AmbiguousAsparagine,
    AmbiguousLeucine,
    AmbiguousGlutamine,
    Selenocysteine,
    Pyrrolysine,
    Unknown,
}
//ARNDCQEGHILKMFPSTWYVBJZUOX

#[derive(Debug)]
pub struct NotACodon;

impl std::fmt::Display for NotACodon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Not a valid codon")
    }
}

impl std::error::Error for NotACodon {}

#[allow(dead_code)]
impl AminoAcid {
    /// The total number of amino acids
    pub const TOTAL_NUMBER: usize = Self::Unknown as usize + 1;
    /// Translate the dna codon into the corresponding amino acid according to the standard DNA codon table.
    /// It returns None for a stop codon.
    /// <https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables>
    /// # Errors
    /// It returns `Err(NotACodon)` when the given codon is not a valid dna codon.
    #[allow(dead_code)]
    pub fn from_dna(dna: &str) -> Result<Option<Self>, NotACodon> {
        match dna.to_lowercase().as_str() {
            "ttt" | "ttc" => Ok(Some(Self::Phenylalanine)),
            "tta" | "ttg" | "ctt" | "ctc" | "cta" | "ctg" => Ok(Some(Self::Leucine)),
            "att" | "atc" | "ata" => Ok(Some(Self::Isoleucine)),
            "atg" => Ok(Some(Self::Methionine)),
            "gtt" | "gtc" | "gta" | "gtg" => Ok(Some(Self::Valine)),
            "tct" | "tcc" | "tca" | "tcg" | "agt" | "agc" => Ok(Some(Self::Serine)),
            "cct" | "ccc" | "cca" | "ccg" => Ok(Some(Self::Proline)),
            "act" | "acc" | "aca" | "acg" => Ok(Some(Self::Threonine)),
            "gct" | "gcc" | "gca" | "gcg" => Ok(Some(Self::Alanine)),
            "tat" | "tac" => Ok(Some(Self::Tyrosine)),
            "taa" | "tag" | "tga" => Ok(None),
            "cat" | "cac" => Ok(Some(Self::Histidine)),
            "caa" | "cag" => Ok(Some(Self::Glutamine)),
            "aat" | "aac" => Ok(Some(Self::Asparagine)),
            "cgt" | "cgc" | "cga" | "cgg" | "aga" | "agg" => Ok(Some(Self::Arginine)),
            "aaa" | "aag" => Ok(Some(Self::Lysine)),
            "gat" | "gac" => Ok(Some(Self::AsparticAcid)),
            "gaa" | "gag" => Ok(Some(Self::GlutamicAcid)),
            "tgt" | "tgc" => Ok(Some(Self::Cysteine)),
            "tgg" => Ok(Some(Self::Tryptophan)),
            "ggt" | "ggc" | "gga" | "ggg" => Ok(Some(Self::Glycine)),
            _ => Err(NotACodon),
        }
    }
}

impl std::str::FromStr for AminoAcid {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::try_from(s)
    }
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

impl TryFrom<&u8> for AminoAcid {
    type Error = ();
    fn try_from(value: &u8) -> Result<Self, Self::Error> {
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

impl TryFrom<u8> for AminoAcid {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        Self::try_from(&value)
    }
}
