use std::fmt::Display;

use uom::num_traits::Zero;

use crate::{
    aminoacids::AminoAcid, element::*, system::f64::*, HasMass, MassSystem, MonoSaccharide,
};

#[derive(Debug, Clone, PartialEq)]
pub struct Peptide {
    pub labile: Vec<Modification>,
    pub n_term: Option<Modification>,
    pub c_term: Option<Modification>,
    pub sequence: Vec<(AminoAcid, Option<Modification>)>,
}

impl Peptide {
    /// [Pro Forma specification](https://github.com/HUPO-PSI/ProForma)
    /// Only supports a subset of the specification, some functions are not possible to be represented.
    pub fn pro_forma(value: &str) -> Result<Self, String> {
        assert!(value.is_ascii());
        let mut peptide = Self {
            labile: Vec::new(),
            n_term: None,
            c_term: None,
            sequence: Vec::new(),
        };
        let mut last_aa = None;
        let chars: Vec<char> = value.chars().collect();
        let mut index = 0;
        let mut c_term = false;

        // N term modification
        while chars[index] == '{' {
            let mut end_index = 0;
            for i in index..value.len() - 1 {
                if chars[i] == '}' {
                    end_index = i + 1;
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
        if chars[index] == '[' {
            let mut end_index = 0;
            for i in index..value.len() - 1 {
                if chars[i] == ']' && chars[i + 1] == '-' {
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

        while index < value.len() {
            match chars[index] {
                '[' => {
                    let mut end_index = 0;
                    for i in index..value.len() {
                        if chars[i] == ']' {
                            end_index = i;
                            break;
                        }
                    }
                    if end_index == 0 {
                        return Err("No valid closing delimiter aminoacid modification".to_string());
                    }
                    let modification = Some((&value[index + 1..end_index]).try_into()?);
                    if c_term {
                        peptide.c_term = modification;
                        if end_index != value.len() - 1 {
                            return Err(
                                "There cannot be any characters after the C terminal modification"
                                    .to_string(),
                            );
                        }
                        break; // Allow the last aa to be placed
                    }
                    peptide.sequence.push((last_aa.unwrap(), modification));
                    last_aa = None;
                    index = end_index + 1;
                }
                '-' => {
                    c_term = true;
                    index += 1;
                }
                ch => {
                    if let Some(aa) = last_aa {
                        peptide.sequence.push((aa, None));
                    }
                    //dbg!(&ch);
                    last_aa = Some(ch.try_into().map_err(|_| "Invalid Amino Acid code")?);
                    index += 1;
                }
            }
        }
        //dbg!(&value, &peptide, &last_aa);
        if let Some(aa) = last_aa {
            peptide.sequence.push((aa, None));
        }
        //dbg!(&peptide);
        Ok(peptide)
    }
}

impl HasMass for Peptide {
    fn mass<M: MassSystem>(&self) -> Mass {
        let mut mass = self
            .n_term
            .as_ref()
            .map_or_else(Mass::zero, HasMass::mass::<M>)
            + self
                .c_term
                .as_ref()
                .map_or_else(Mass::zero, HasMass::mass::<M>);
        for position in &self.sequence {
            mass += position.0.mass::<M>()
                + position
                    .1
                    .as_ref()
                    .map_or_else(Mass::zero, HasMass::mass::<M>);
        }
        mass
    }
}

impl Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(m) = &self.n_term {
            write!(f, "[{m}]-")?;
        }
        for position in &self.sequence {
            write!(f, "{}", position.0.char())?;
            if let Some(m) = &position.1 {
                write!(f, "[{m}]")?;
            }
        }
        if let Some(m) = &self.c_term {
            write!(f, "-[{m}]")?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Modification {
    /// Monoisotopic mass shift
    Mass(Mass),
    #[allow(non_snake_case)]
    Formula(Vec<(Element, isize)>),
    Glycan(Vec<(MonoSaccharide, isize)>),
}

impl HasMass for Modification {
    fn mass<M: MassSystem>(&self) -> Mass {
        match self {
            Self::Mass(m) => *m,
            Self::Formula(elements) => elements.iter().map(|m| m.1 as f64 * m.0.mass::<M>()).sum(),
            Self::Glycan(monosaccharides) => monosaccharides
                .iter()
                .map(|m| m.1 as f64 * m.0.mass::<M>())
                .sum(),
        }
    }
}

impl TryFrom<&str> for Modification {
    type Error = String;
    /// TODO: support more parts of the Pro Forma spec:
    ///     * Think about a way to enforce knowledge about the monoisotopic nature of the modifications in normal mass shifts
    ///     * Allow zero mass gap (X[+365])
    ///     * Do not crash on other input, provide shift as string with 0 mass?
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        dbg!(&value);
        match value.split_once(':') {
            Some(("Formula", tail)) => Ok(Self::Formula(parse_named_counter(
                tail,
                ELEMENT_PARSE_LIST,
                true,
            )?)),
            Some(("Glycan", tail)) => Ok(Self::Glycan(parse_named_counter(
                tail,
                crate::GLYCAN_PARSE_LIST,
                false,
            )?)),
            Some((head, _tail)) => Err(format!("Does not support these types yet: {head}")),
            None => Ok(Self::Mass(Mass::new::<dalton>(
                value.parse::<f64>().map_err(|_| "Invalid mass shift")?,
            ))),
        }
    }
}

fn parse_named_counter<T: Copy>(
    value: &str,
    names: &[(&str, T)],
    allow_negative: bool,
) -> Result<Vec<(T, isize)>, String> {
    let mut index = 0;
    let mut output = Vec::new();
    while index < value.len() {
        if value[index..].starts_with(' ') {
            index += 1;
        } else {
            let mut found = false;
            for name in names {
                if value[index..].starts_with(name.0) {
                    index += name.0.len();
                    let num = &value[index..]
                        .chars()
                        .take_while(|c| c.is_ascii_digit() || (allow_negative && *c == '-'))
                        .collect::<String>();
                    if num.is_empty() {
                        output.push((name.1, 1));
                    } else {
                        output.push((name.1, dbg!(num).parse().unwrap()));
                        index += num.len();
                    }
                    found = true;
                    break; // Names loop
                }
            }
            if !found {
                return Err(format!("Name not recognised {}", &value[index..]));
            }
        }
    }
    Ok(output)
}

impl Display for Modification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Mass(m) => {
                write!(f, "{:+}", m.value).unwrap();
            }
            Self::Formula(elements) => {
                write!(f, "Formula:{}", Element::hill_notation(elements)).unwrap();
            }
            Self::Glycan(monosaccharides) => write!(
                f,
                "Glycan:{}",
                monosaccharides
                    .iter()
                    .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
            )
            .unwrap(),
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
}
