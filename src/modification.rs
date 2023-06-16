use std::fmt::Display;

use uom::num_traits::Zero;

use crate::{
    dalton,
    element::{Element, ELEMENT_PARSE_LIST},
    HasMass, Mass, MassSystem, MonoSaccharide,
};

#[derive(Debug, Clone, PartialEq)]
pub enum Modification {
    /// Monoisotopic mass shift
    Mass(Mass),
    #[allow(non_snake_case)]
    Formula(Vec<(Element, isize)>),
    Glycan(Vec<(MonoSaccharide, isize)>),
    Predefined(
        &'static [(Element, isize)],
        &'static [(MonoSaccharide, isize)],
    ),
}

impl HasMass for Modification {
    fn mass<M: MassSystem>(&self) -> Mass {
        match self {
            Self::Mass(m) => *m,
            Self::Formula(elements) => elements.mass::<M>(),
            Self::Glycan(monosaccharides) => monosaccharides.mass::<M>(),
            Self::Predefined(elements, monosaccharides) => {
                elements.mass::<M>() + monosaccharides.mass::<M>()
            }
        }
    }
}

fn parse_single_modification(part: &str) -> Result<Option<Modification>, String> {
    let parsed = part
        .split_once(':')
        .map(|v| (v.0.to_ascii_lowercase(), v.1));
    if let Some((head, tail)) = parsed {
        match (head.as_str(), tail) {
            ("formula", tail) => Ok(Some(Modification::Formula(parse_named_counter(
                tail,
                ELEMENT_PARSE_LIST,
                true,
            )?))),
            ("glycan", tail) => Ok(Some(Modification::Glycan(parse_named_counter(
                tail,
                crate::GLYCAN_PARSE_LIST,
                false,
            )?))),
            ("info", _) => Ok(None),
            ("obs", tail) => Ok(Some(Modification::Mass(Mass::new::<dalton>(
                tail.parse::<f64>().map_err(|_| "Invalid mass shift")?,
            )))),
            (head, _tail) => Err(format!("Does not support these types yet: {head}")),
        }
    } else {
        Ok(Some(Modification::Mass(Mass::new::<dalton>(
            part.parse::<f64>().map_err(|_| "Invalid mass shift")?,
        ))))
    }
}

impl TryFrom<&str> for Modification {
    type Error = String;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        // Because multiple modifications could be chained with the pipe operator
        // the parsing iterates over all links until it finds one it understands
        // it then returns that one. If no 'understandable' links are found it
        // returns the last link, if this is an info it returns a mass shift of 0,
        // but if any of the links returned an error it returns the last error.
        let mut last_result = Ok(None);
        let mut last_error = None;
        for part in value.split('|') {
            last_result = parse_single_modification(part);
            if let Ok(Some(m)) = last_result {
                return Ok(m);
            }
            if let Err(er) = &last_result {
                last_error = Some(er.clone());
            }
        }
        last_error.map_or_else(
            || last_result.map(|m| m.unwrap_or_else(|| Self::Mass(Mass::zero()))),
            Err,
        )
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
                        output.push((name.1, num.parse().unwrap()));
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
            Self::Predefined(elements, monosaccharides) => write!(
                f,
                "Predefined:{};{}",
                Element::hill_notation(elements),
                monosaccharides
                    .iter()
                    .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
            )
            .unwrap(),
        }
        Ok(())
    }
}
