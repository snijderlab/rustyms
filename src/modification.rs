use std::fmt::Display;

use uom::num_traits::Zero;

use crate::{
    dalton,
    element::{Element, ELEMENT_PARSE_LIST},
    helper_functions::*,
    ontologies::{PSI_MOD_ONTOLOGY, UNIMOD_ONTOLOGY},
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
            ("unimod", tail) => {
                let id = tail.parse::<usize>().map_err(|e| e.to_string())?;
                find_id_in_ontology(id, UNIMOD_ONTOLOGY)
                    .map(Some)
                    .map_err(|_| format!("{tail} is not a valid Unimod accession number"))
            }
            ("mod", tail) => {
                let id = tail.parse::<usize>().map_err(|e| e.to_string())?;
                find_id_in_ontology(id, PSI_MOD_ONTOLOGY)
                    .map(Some)
                    .map_err(|_| format!("{tail} is not a valid PSI-MOD accession number"))
            }
            ("u", tail) => find_in_ontology(tail, UNIMOD_ONTOLOGY)
                .map_err(|_| numerical_mod(tail))
                .flat_err()
                .map(Some)
                .map_err(|_| format!("Not a valid Unimod modification: {tail}")),
            ("m", tail) => find_in_ontology(tail, PSI_MOD_ONTOLOGY)
                .map_err(|_| numerical_mod(tail))
                .flat_err()
                .map(Some)
                .map_err(|_| format!("Not a valid PSI-MOD modification: {tail}")),
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
            ("obs", tail) => numerical_mod(tail).map(Some),
            (head, _tail) => find_in_ontology(part, UNIMOD_ONTOLOGY)
                .map(Some)
                .map_err(|_| format!("Does not support these types yet: {head}")),
        }
    } else {
        find_in_ontology(part, UNIMOD_ONTOLOGY)
            .map_err(|_| find_in_ontology(part, PSI_MOD_ONTOLOGY))
            .flat_err()
            .map_err(|_| numerical_mod(part))
            .flat_err()
            .map(Some)
            .map_err(|_| format!("Not a valid delta mass, Unimod, or PSI-MOD modification: {part}"))
    }
}

fn find_in_ontology(
    code: &str,
    ontology: &[(usize, &str, &str, Modification)],
) -> Result<Modification, ()> {
    let code = code.to_ascii_lowercase();
    for option in ontology {
        if option.1 == code || option.2 == code {
            return Ok(option.3.clone());
        }
    }
    Err(())
}

fn find_id_in_ontology(
    id: usize,
    ontology: &[(usize, &str, &str, Modification)],
) -> Result<Modification, ()> {
    for option in ontology {
        if option.0 == id {
            return Ok(option.3.clone());
        }
    }
    Err(())
}

fn numerical_mod(text: &str) -> Result<Modification, String> {
    text.parse().map_or_else(
        |_| Err("Invalid number".to_string()),
        |n| Ok(Modification::Mass(Mass::new::<dalton>(n))),
    )
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
