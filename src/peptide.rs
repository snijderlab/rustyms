use std::fmt::Display;

use uom::num_traits::Zero;

use crate::{aminoacids::AminoAcid, system::f64::*, MassSystem};

#[derive(Debug, Clone)]
pub struct Peptide {
    pub n_term: Option<Modification>,
    pub c_term: Option<Modification>,
    pub sequence: Vec<(AminoAcid, Option<Modification>)>,
}

impl Peptide {
    /// https://github.com/HUPO-PSI/ProForma
    /// Only supports a subset of the specification, some functions are not possible to be represented.
    pub fn pro_forma(value: &str) -> Result<Self, &'static str> {
        assert!(value.is_ascii());
        let mut peptide = Peptide {
            n_term: None,
            c_term: None,
            sequence: Vec::new(),
        };
        let mut last_aa = None;
        let chars: Vec<char> = value.chars().collect();
        let mut index = 0;

        // N term modification
        if chars[index] == '{' {
            let mut end_index = 0;
            for i in index..value.len() {
                if chars[i] == '}' {
                    end_index = i;
                    break;
                }
            }
            if end_index == 0 {
                return Err("No valid closing delimiter for N term modification");
            }
            peptide.n_term = Some((&value[index + 1..end_index]).try_into()?);
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
                        return Err("No valid closing delimiter aminoacid modification");
                    }
                    peptide.sequence.push((
                        last_aa.unwrap(),
                        Some((&value[index + 1..end_index]).try_into()?),
                    ));
                    last_aa = None;
                    index = end_index + 1;
                }
                '{' => {
                    let mut end_index = 0;
                    for i in index..value.len() {
                        if chars[i] == '}' {
                            end_index = i;
                            break;
                        }
                    }
                    if end_index == 0 {
                        return Err("No valid closing delimiter for C term modification");
                    }
                    if end_index != value.len() - 1 {
                        return Err("Characters present after C terminus");
                    }
                    peptide.c_term = Some((&value[index + 1..end_index]).try_into()?);
                    break; // SHOULD be last characters of the string
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
            peptide.sequence.push((aa, None))
        }
        //dbg!(&peptide);
        Ok(peptide)
    }

    pub fn mass<M: MassSystem>(&self) -> Mass {
        let mut mass = self.n_term.as_ref().map_or(Mass::zero(), |m| m.mass::<M>())
            + self.c_term.as_ref().map_or(Mass::zero(), |m| m.mass::<M>());
        for position in &self.sequence {
            mass += position.0.mass::<M>()
                + position.1.as_ref().map_or(Mass::zero(), |m| m.mass::<M>());
        }
        mass
    }
}

impl Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(m) = &self.n_term {
            write!(f, "{{{m}}}")?;
        }
        for position in &self.sequence {
            write!(f, "{}", position.0.char())?;
            if let Some(m) = &position.1 {
                write!(f, "[{m}]")?;
            }
        }
        if let Some(m) = &self.c_term {
            write!(f, "{{{m}}}")?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub enum Modification {
    /// Monoisotopic mass shift
    Mass(Mass),
    Formula {
        H: isize,
        C: isize,
        N: isize,
        O: isize,
        F: isize,
        P: isize,
        S: isize,
        Se: isize,
    },
}

impl Modification {
    pub fn mass<M: MassSystem>(&self) -> Mass {
        match self {
            Self::Mass(m) => *m,
            Self::Formula {
                H,
                C,
                N,
                O,
                F,
                P,
                S,
                Se,
            } => da(M::H * (*H as f64)
                + M::C * (*C as f64)
                + M::N * (*N as f64)
                + M::O * (*O as f64)
                + M::F * (*F as f64)
                + M::P * (*P as f64)
                + M::S * (*S as f64)
                + M::Se * (*Se as f64)),
        }
    }
}

impl TryFrom<&str> for Modification {
    type Error = &'static str;
    /// TODO: support more parts of the spec:
    ///     * Glycans (at least 'formula' based: 'HexNac')
    ///     * Formulas (generic over mass system)
    ///     * Think about a way to enforce knowledge about the monoisotopic nature of the modifications in normal mass shifts
    ///     * Allow zero mass gap (X[+365])
    ///     * Do not crash on other input, provide shift as string with 0 mass
    fn try_from(value: &str) -> Result<Modification, Self::Error> {
        dbg!(&value);
        let mass = value.parse::<f64>().map_err(|_| "Invalid mass shift")?;
        Ok(Modification::Mass(Mass::new::<dalton>(mass)))
    }
}

impl Display for Modification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Modification::Mass(m) => {
                write!(f, "{:+}", m.value).unwrap();
            }
            Self::Formula {
                H,
                C,
                N,
                O,
                F,
                P,
                S,
                Se,
            } => {
                let mut pr = |n: &isize, c: &str| {
                    if *n != 0 {
                        write!(f, "{c}{n}").unwrap();
                    }
                };
                // Display in Hill order (if C: C first followed by H rest alphabetical, otherwise all alphabetical)
                if *C != 0 {
                    pr(C, "C");
                    pr(H, "H");
                    pr(F, "F");
                    pr(N, "N");
                    pr(O, "O");
                    pr(P, "P");
                    pr(S, "S");
                    pr(Se, "Se");
                } else {
                    pr(F, "F");
                    pr(H, "H");
                    pr(N, "N");
                    pr(O, "O");
                    pr(P, "P");
                    pr(S, "S");
                    pr(Se, "Se");
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::{aminoacids::AminoAcid, MonoIsotopic};

    use super::Peptide;

    #[test]
    fn parse_pro_forma_simple() {
        let inp = "AAA";
        let peptide = Peptide::pro_forma(inp).unwrap();
        assert_eq!(
            peptide.mass::<MonoIsotopic>(),
            AminoAcid::Alanine.mass::<MonoIsotopic>() * 3.0
        );
        assert_eq!(peptide.sequence.len(), 3);
        assert_eq!(peptide.to_string(), "AAA".to_string());
    }

    #[test]
    fn parse_pro_forma_termini() {
        let inp = "{+12.2}AAA{-12.2}";
        let peptide = Peptide::pro_forma(inp).unwrap();
        assert_eq!(
            peptide.mass::<MonoIsotopic>(),
            AminoAcid::Alanine.mass::<MonoIsotopic>() * 3.0
        );
        assert_eq!(peptide.sequence.len(), 3);
        assert_eq!(peptide.to_string(), "{+12.2}AAA{-12.2}".to_string());
    }
}
