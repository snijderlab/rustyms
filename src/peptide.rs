#![warn(dead_code)]

use std::{
    fmt::Display,
    ops::{Range, RangeBounds},
};

use itertools::Itertools;

use crate::{
    aminoacids::AminoAcid,
    helper_functions::ResultExtensions,
    modification::{Modification, ReturnModification},
    system::f64::*,
    Fragment, FragmentType, HasMass, MassSystem, Model,
};

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Peptide {
    pub labile: Vec<Modification>,
    pub n_term: Option<Modification>,
    pub c_term: Option<Modification>,
    pub sequence: Vec<SequenceElement>,
    total_ambiguous_modifications: usize,
    /// For each ambiguous modification list all possible positions it can be placed on.
    /// Index by the ambiguous modification id.
    ambiguous_modifications: Vec<Vec<usize>>,
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
    #[allow(clippy::too_many_lines)]
    pub fn pro_forma(value: &str) -> Result<Self, String> {
        let mut peptide = Self::default();
        let chars: &[u8] = value.as_bytes();
        let mut index = 0;
        let mut c_term = false;
        let mut ambiguous_aa_counter = 0;
        let mut ambiguous_aa = None;
        let mut ambiguous_lookup = Vec::new();
        let mut ambiguous_found_positions = Vec::new();

        // Labile modification(s)
        while chars[index] == b'{' {
            // TODO: Should I allow for the used of paired curly brackets inside as well?
            let end_index = index
                + 1
                + chars[index..]
                    .iter()
                    .position(|c| *c == b'}')
                    .ok_or(format!(
                        "No valid closing delimiter for labile modification [index: {index}]"
                    ))?;
            peptide.labile.push(
                Modification::try_from(&value[index + 1..end_index - 1], &mut ambiguous_lookup)
                    .map(|m| {
                        if let ReturnModification::Defined(m) = m {
                            Ok(m)
                        } else {
                            Err("A labile modification cannot be ambiguous".to_string())
                        }
                    })
                    .flat_err()?,
            );
            index = end_index;
        }
        // N term modification
        if chars[index] == b'[' {
            let mut end_index = 0;
            for i in index..value.len() - 1 {
                if chars[i] == b']' && chars[i + 1] == b'-' {
                    end_index = i + 1;
                    break;
                }
            }
            if end_index == 0 {
                return Err(format!(
                    "No valid closing delimiter for N term modification [index: {index}]"
                ));
            }
            peptide.n_term = Some(
                Modification::try_from(&value[index + 1..end_index - 1], &mut ambiguous_lookup)
                    .map(|m| {
                        if let ReturnModification::Defined(m) = m {
                            Ok(m)
                        } else {
                            Err("A labile modification cannot be ambiguous".to_string())
                        }
                    })
                    .flat_err()?,
            );
            index = end_index + 1;
        }

        // Rest of the sequence
        while index < chars.len() {
            match chars[index] {
                b'(' if chars[index + 1] == b'?' && ambiguous_aa.is_none() => {
                    ambiguous_aa = Some(ambiguous_aa_counter);
                    ambiguous_aa_counter += 1;
                    index += 2;
                }
                b')' if ambiguous_aa.is_some() => {
                    ambiguous_aa = None;
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
                        return Err(format!(
                            "No valid closing delimiter aminoacid modification [index: {index}]"
                        ));
                    }
                    let modification = Modification::try_from(
                        &value[index + 1..end_index],
                        &mut ambiguous_lookup,
                    )?;
                    if c_term {
                        peptide.c_term =
                            Some(if let ReturnModification::Defined(m) = modification {
                                Ok(m)
                            } else {
                                Err("A labile modification cannot be ambiguous".to_string())
                            }?);
                        if end_index != value.len() - 1 {
                            return Err(
                                format!("There cannot be any characters after the C terminal modification [index: {index}]"),
                            );
                        }
                        break;
                    }
                    match peptide.sequence.last_mut() {
                        Some(aa) => match modification {
                            ReturnModification::Defined(m) => aa.modifications.push(m),
                            ReturnModification::Preferred(id, localisation_score) =>
                            ambiguous_found_positions.push(
                                (peptide.sequence.len() -1, true, id, localisation_score)),
                            ReturnModification::Referenced(id, localisation_score) =>
                            ambiguous_found_positions.push(
                                (peptide.sequence.len() -1, false, id, localisation_score)),
                        },
                        None => {
                            return Err(
                                format!("A modification cannot be placed before any amino acid [index: {index}]")
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
                    peptide.sequence.push(SequenceElement::new(
                        ch.try_into().map_err(|_| "Invalid Amino Acid code")?,
                        ambiguous_aa,
                    ));
                    index += 1;
                }
            }
        }
        // Fill in ambiguous positions
        for (index, preferred, id, localisation_score) in ambiguous_found_positions {
            peptide.sequence[index].possible_modifications.push(
                AmbiguousModification {
                    id,
                    modification: ambiguous_lookup[id].1.as_ref().cloned().ok_or(format!("Ambiguous modification {} did not have a definition for the actual modification", ambiguous_lookup[id].0.as_ref().map_or(id.to_string(), ToString::to_string)))?,
                    localisation_score,
                    group: ambiguous_lookup[id].0.as_ref().map(|n| (n.to_string(), preferred)) });
        }
        peptide.ambiguous_modifications = ambiguous_found_positions
            .iter()
            .group_by(|p| p.2)
            .into_iter()
            .sorted_by(|(key1, _), (key2, _)| key1.cmp(&key2))
            .map(|(key, group)| group.into_iter().map(|p| p.0).collect())
            .collect();

        // Check all placement rules

        Ok(peptide)
    }

    /// Generate all possible patterns for the ambiguous positions (ambiguous id, aa position)
    fn ambiguous_patterns(&self, aa_range: impl RangeBounds<usize>) -> Vec<Vec<(usize, usize)>> {
        self.ambiguous_modifications.iter().enumerate().fold(
            vec![Vec::new()],
            |acc, (id, possibilities)| {
                acc.into_iter()
                    .flat_map(|path| {
                        possibilities
                            .iter()
                            .filter(|pos| aa_range.contains(pos))
                            .map(|pos| {
                                let mut new = path.clone();
                                new.push((id, *pos));
                                new
                            })
                    })
                    .collect()
            },
        )
    }

    pub fn mass<M: MassSystem>(&self) -> Option<Vec<Mass>> {
        let base_mass = self
            .n_term
            .as_ref()
            .map_or_else(|| da(M::H), HasMass::mass::<M>)
            + self
                .c_term
                .as_ref()
                .map_or_else(|| da(M::OH), HasMass::mass::<M>);
        let patterns = self.ambiguous_patterns(..);
        let mut masses = patterns
            .iter()
            .map(|_| base_mass.clone())
            .collect::<Vec<_>>();

        for (mut mass, pattern) in masses.iter_mut().zip(patterns) {
            for (index, pos) in self.sequence.iter().enumerate() {
                *mass += pos.mass::<M>(
                    &pattern
                        .iter()
                        .filter_map(|(id, pos)| (index == *pos).then(|| *id))
                        .collect::<Vec<_>>(),
                )?;
            }
        }
        Some(masses)
    }

    /// Generate the theoretical fragments for this peptide, with the given maximal charge of the fragments, and the given model.
    ///
    /// # Panics
    /// If `max_charge` outside the range `1..=u64::MAX`.
    pub fn generate_theoretical_fragments<M: MassSystem>(
        &self,
        max_charge: Charge,
        model: &Model,
    ) -> Vec<Fragment> {
        assert!(max_charge.value >= 1.0);
        assert!(max_charge.value <= u64::MAX as f64);
        let mut output = Vec::with_capacity(20 * self.sequence.len() + 75); // Empirically derived required size of the buffer (Derived from Hecklib)
        for index in 0..self.sequence.len() {
            let n_term = self.n_term.mass::<M>()
                + self.sequence[0..index]
                    .iter()
                    .fold(Mass::zero(), |acc, aa| acc + aa.mass::<M>());
            let c_term = self.c_term.mass::<M>()
                + self.sequence[index + 1..self.sequence.len()]
                    .iter()
                    .fold(Mass::zero(), |acc, aa| acc + aa.mass::<M>());
            output.append(&mut self.sequence[index].aminoacid.fragments::<M>(
                // TODO: does this take the mods on the current position into account, also take possible mods into account
                n_term,
                c_term,
                max_charge,
                index,
                self.sequence.len(),
                &model.ions(index, self.sequence.len()),
            ));
        }
        output.push(Fragment::new(
            self.mass::<M>(),
            max_charge,
            FragmentType::precursor,
        ));
        output
    }
}

impl Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut placed = Vec::new();
        if let Some(m) = &self.n_term {
            write!(f, "[{m}]-")?;
        }
        let mut last_ambiguous = None;
        for position in &self.sequence {
            placed.extend(position.display(f, &placed, last_ambiguous)?);
            last_ambiguous = position.ambiguous;
        }
        if last_ambiguous.is_some() {
            write!(f, ")")?;
        }
        if let Some(m) = &self.c_term {
            write!(f, "-[{m}]")?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SequenceElement {
    pub aminoacid: AminoAcid,
    pub modifications: Vec<Modification>,
    pub possible_modifications: Vec<AmbiguousModification>,
    pub ambiguous: Option<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AmbiguousModification {
    pub id: usize,
    pub modification: Modification,
    pub localisation_score: Option<f64>,
    /// If this is a named group contain the name and track if this is the preferred location or not
    pub group: Option<(String, bool)>,
}

impl SequenceElement {
    pub const fn new(aminoacid: AminoAcid, ambiguous: Option<usize>) -> Self {
        Self {
            aminoacid,
            modifications: Vec::new(),
            possible_modifications: Vec::new(),
            ambiguous,
        }
    }

    fn display(
        &self,
        f: &mut std::fmt::Formatter<'_>,
        placed: &[usize],
        last_ambiguous: Option<usize>,
    ) -> Result<Vec<usize>, std::fmt::Error> {
        let mut extra_placed = Vec::new();
        if last_ambiguous.is_some() && last_ambiguous != self.ambiguous {
            write!(f, ")")?;
        }
        if self.ambiguous.is_some() && last_ambiguous != self.ambiguous {
            write!(f, "(?")?;
        }
        write!(f, "{}", self.aminoacid.char())?;
        for m in &self.modifications {
            write!(f, "[{m}]")?;
        }
        for m in &self.possible_modifications {
            write!(
                f,
                "[{}#{}{}]",
                m.group.as_ref().map_or(
                    if placed.contains(&m.id) {
                        String::new()
                    } else {
                        extra_placed.push(m.id);
                        m.modification.to_string()
                    },
                    |group| if group.1 {
                        m.modification.to_string()
                    } else {
                        String::new()
                    }
                ),
                m.group
                    .as_ref()
                    .map_or(m.id.to_string(), |g| g.0.to_string()),
                m.localisation_score
                    .map(|v| format!("({v})"))
                    .unwrap_or_default()
            )?;
        }
        Ok(extra_placed)
    }

    fn mass<M: MassSystem>(&self, selected_ambiguous: &[usize]) -> Option<Mass> {
        if self.aminoacid == AminoAcid::B || self.aminoacid == AminoAcid::Z {
            None
        } else {
            Some(
                self.aminoacid.mass::<M>()
                    + self.modifications.mass::<M>()
                    + self
                        .possible_modifications
                        .iter()
                        .filter_map(|m| {
                            selected_ambiguous
                                .contains(&m.id)
                                .then(|| m.modification.mass::<M>())
                        })
                        .sum(),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::MonoIsotopic;

    use super::Peptide;

    #[test]
    fn parse_glycan() {
        let glycan = Peptide::pro_forma("A[Glycan:Hex]").unwrap();
        let spaces = Peptide::pro_forma("A[Glycan:    Hex    ]").unwrap();
        assert_eq!(glycan.sequence.len(), 1);
        assert_eq!(spaces.sequence.len(), 1);
        assert_eq!(glycan, spaces);
        let incorrect = Peptide::pro_forma("A[Glycan:Hec]");
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
    fn parse_labile() {
        let with = Peptide::pro_forma("{Formula:C6H10O5}A").unwrap();
        let without = Peptide::pro_forma("A").unwrap();
        assert_eq!(with.sequence.len(), 1);
        assert_eq!(without.sequence.len(), 1);
        assert_eq!(with.mass::<MonoIsotopic>(), without.mass::<MonoIsotopic>());
        assert_eq!(with.labile[0].to_string(), "Formula:C6H10O5".to_string());
    }

    #[test]
    fn parse_ambiguous_modification() {
        let with = Peptide::pro_forma("A[Phospho#g0]A[#g0]").unwrap();
        let without = Peptide::pro_forma("AA").unwrap();
        assert_eq!(with.sequence.len(), 2);
        assert_eq!(without.sequence.len(), 2);
        assert_eq!(with.sequence[0].possible_modifications.len(), 1);
        assert_eq!(with.sequence[1].possible_modifications.len(), 1);
        assert!(Peptide::pro_forma("A[#g0]A[#g0]").is_err());
        assert!(Peptide::pro_forma("A[Phospho#g0]A[Phospho#g0]").is_err());
        assert!(Peptide::pro_forma("A[Phospho#g0]A[#g0(0.o1)]").is_err());
        assert_eq!(
            Peptide::pro_forma("A[+12#g0]A[#g0]").unwrap().to_string(),
            "A[+12#g0]A[#g0]".to_string()
        );
        assert_eq!(
            Peptide::pro_forma("A[#g0]A[+12#g0]").unwrap().to_string(),
            "A[#g0]A[+12#g0]".to_string()
        );
    }

    #[test]
    fn parse_ambiguous_aminoacid() {
        let with = Peptide::pro_forma("(?AA)C(?A)(?A)").unwrap();
        let without = Peptide::pro_forma("AACAA").unwrap();
        assert_eq!(with.sequence.len(), 5);
        assert_eq!(without.sequence.len(), 5);
        assert!(with.sequence[0].ambiguous.is_some());
        assert!(with.sequence[1].ambiguous.is_some());
        assert_eq!(with.mass::<MonoIsotopic>(), without.mass::<MonoIsotopic>());
        assert_eq!(with.to_string(), "(?AA)C(?A)(?A)".to_string());
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
        let peptide = dbg!(Peptide::pro_forma("A[Cation:Na]A[U:Gln->pyro-Glu]A"));
        assert!(peptide.is_ok());
    }
}
