#![warn(dead_code)]

use std::fmt::Display;

use crate::{
    aminoacids::AminoAcid, modification::Modification, system::f64::*, HasMass, MassSystem,
};

#[derive(Debug, Clone, PartialEq)]
pub struct Peptide {
    pub labile: Vec<Modification>,
    pub n_term: Option<Modification>,
    pub c_term: Option<Modification>,
    pub sequence: Vec<(AminoAcid, Vec<Modification>, bool)>,
}

impl Peptide {
    /// Get the number of amino acids making up this peptide
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
    /// Check if there are any amino acids in this peptide
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// [Pro Forma specification](https://github.com/HUPO-PSI/ProForma)
    /// Only supports a subset of the specification, some functions are not possible to be represented.
    ///
    /// # Errors
    /// It fails when the string is not a valid Pro Forma string, with a minimal error message to help debug the cause.
    pub fn pro_forma(value: &str) -> Result<Self, String> {
        let mut peptide = Self {
            labile: Vec::new(),
            n_term: None,
            c_term: None,
            sequence: Vec::new(),
        };
        let chars: &[u8] = value.as_bytes();
        let mut index = 0;
        let mut c_term = false;
        let mut ambiguous_aa = false;

        // N term modification
        while chars[index] == b'{' {
            let mut end_index = 0;
            for (i, ch) in chars[index..].iter().enumerate() {
                if *ch == b'}' {
                    end_index = index + i + 1;
                    break;
                }
            }
            if end_index == 0 {
                return Err("No valid closing delimiter for labile modification".to_string());
            }
            peptide
                .labile
                .push((&value[index + 1..end_index - 1]).try_into()?);
            index = end_index + 1;
        }
        if chars[index] == b'[' {
            let mut end_index = 0;
            for i in index..value.len() - 1 {
                if chars[i] == b']' && chars[i + 1] == b'-' {
                    end_index = i + 1;
                    break;
                }
            }
            if end_index == 0 {
                return Err("No valid closing delimiter for N term modification".to_string());
            }
            peptide.n_term = Some((&value[index + 1..end_index - 1]).try_into()?);
            index = end_index + 1;
        }

        while index < chars.len() {
            match chars[index] {
                b'(' if chars[index + 1] == b'?' => {
                    ambiguous_aa = true;
                    index += 2;
                }
                b')' if ambiguous_aa => {
                    ambiguous_aa = false;
                    index += 1;
                }
                b'[' => {
                    let mut end_index = 0;
                    for (i, ch) in chars[index..].iter().enumerate() {
                        if *ch == b']' {
                            end_index = index + i;
                            break;
                        }
                    }
                    if end_index == 0 {
                        return Err("No valid closing delimiter aminoacid modification".to_string());
                    }
                    let modification = Modification::try_from(&value[index + 1..end_index])?;
                    if c_term {
                        peptide.c_term = Some(modification);
                        if end_index != value.len() - 1 {
                            return Err(
                                "There cannot be any characters after the C terminal modification"
                                    .to_string(),
                            );
                        }
                        break;
                    }
                    match peptide.sequence.last_mut() {
                        Some(aa) => aa.1.push(modification),
                        None => {
                            return Err(
                                "A modification cannot be placed before any amino acid".to_string()
                            )
                        }
                    }
                    index = end_index + 1;
                }
                b'-' => {
                    c_term = true;
                    index += 1;
                }
                ch => {
                    peptide.sequence.push((
                        ch.try_into().map_err(|_| "Invalid Amino Acid code")?,
                        Vec::new(),
                        ambiguous_aa,
                    ));
                    index += 1;
                }
            }
        }
        Ok(peptide)
    }
}

impl HasMass for Peptide {
    fn mass<M: MassSystem>(&self) -> Mass {
        self.n_term
            .as_ref()
            .map_or_else(|| da(M::H), HasMass::mass::<M>)
            + self
                .c_term
                .as_ref()
                .map_or_else(|| da(M::OH), HasMass::mass::<M>)
            + self.sequence.mass::<M>()
    }
}

impl Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(m) = &self.n_term {
            write!(f, "[{m}]-")?;
        }
        for position in &self.sequence {
            write!(f, "{}", position.0.char())?;
            for m in &position.1 {
                write!(f, "[{m}]")?;
            }
        }
        if let Some(m) = &self.c_term {
            write!(f, "-[{m}]")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::{aminoacids::AminoAcid, HasMass, MonoIsotopic};

    use super::Peptide;

    #[test]
    fn test_many_pro_forma() {
        for v in ["AAA", "[+12.2]-AAA-[-12.2]", "[+12.2]-AAA[-10.1]-[-2.1]"] {
            assume_mass_3ala(v);
        }
    }

    fn assume_mass_3ala(value: &str) {
        let peptide = Peptide::pro_forma(value).unwrap();
        assert_eq!(
            peptide.mass::<MonoIsotopic>(),
            AminoAcid::Alanine.mass::<MonoIsotopic>() * 3.0
        );
        assert_eq!(peptide.sequence.len(), 3);
        assert_eq!(peptide.to_string(), value.to_string());
    }

    #[test]
    fn parse_glycan() {
        let glycan = Peptide::pro_forma("A[Glycan:Hex]").unwrap();
        let spaces = Peptide::pro_forma("A[Glycan:    Hex    ]").unwrap();
        assert_eq!(glycan.sequence.len(), 1);
        assert_eq!(spaces.sequence.len(), 1);
        assert_eq!(glycan, spaces);
        let incorrect = Peptide::pro_forma("A[Glycan:Hec]");
        assert!(incorrect.is_err());
        let incorrect = Peptide::pro_forma("A[glycan:Hex]");
        assert!(incorrect.is_err());
    }

    #[test]
    fn parse_formula() {
        let peptide = Peptide::pro_forma("A[Formula:C6H10O5]").unwrap();
        let glycan = Peptide::pro_forma("A[Glycan:Hex]").unwrap();
        assert_eq!(peptide.sequence.len(), 1);
        assert_eq!(glycan.sequence.len(), 1);
        assert_eq!(
            glycan.mass::<MonoIsotopic>(),
            peptide.mass::<MonoIsotopic>()
        );
    }

    #[test]
    fn parse_hard_tags() {
        let peptide = Peptide::pro_forma("A[Formula:C6H10O5|INFO:hello world 🦀]").unwrap();
        let glycan = Peptide::pro_forma(
            "A[info:you can define a tag multiple times|Glycan:Hex|Formula:C6H10O5]",
        )
        .unwrap();
        assert_eq!(peptide.sequence.len(), 1);
        assert_eq!(glycan.sequence.len(), 1);
        assert_eq!(
            glycan.mass::<MonoIsotopic>(),
            peptide.mass::<MonoIsotopic>()
        );
    }

    #[test]
    fn parse_unimod() {
        let peptide = dbg!(Peptide::pro_forma(
            "A[Cation:Na]A[U:Gln->pyro-Glu]A[pyro_Glu]"
        ));
        assert!(peptide.is_ok());
    }
}
