#![warn(dead_code)]

use std::{fmt::Display, ops::RangeBounds};

use crate::{Element, MolecularFormula};
use itertools::Itertools;
use uom::num_traits::Zero;

use crate::{
    aminoacids::AminoAcid,
    helper_functions::ResultExtensions,
    modification::{Modification, ReturnModification},
    system::f64::*,
    Chemical, Fragment, FragmentType, Model,
};

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Peptide {
    pub labile: Vec<Modification>,
    pub n_term: Option<Modification>,
    pub c_term: Option<Modification>,
    pub sequence: Vec<SequenceElement>,
    /// For each ambiguous modification list all possible positions it can be placed on.
    /// Index by the ambiguous modification id.
    pub ambiguous_modifications: Vec<Vec<usize>>,
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

    pub fn n_term(&self) -> MolecularFormula {
        self.n_term.as_ref().map_or_else(
            || molecular_formula!(H 1),
            |m| molecular_formula!(H 1) + m.formula(),
        )
    }

    pub fn c_term(&self) -> MolecularFormula {
        self.c_term.as_ref().map_or_else(
            || molecular_formula!(H 1 O 1),
            |m| molecular_formula!(H 1 O 1) + m.formula(),
        )
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
        let mut global_modifications = Vec::new();

        // Global modification(s)
        while chars[index] == b'<' {
            let end_index = index
                + 1
                + chars[index..]
                    .iter()
                    .position(|c| *c == b'>')
                    .ok_or(format!(
                        "No valid closing delimiter for global modification [index: {index}]"
                    ))?;
            if let Some(offset) = chars[index..].iter().position(|c| *c == b'@') {
                let at_index = index + 1 + offset;
                if !chars[index + 1] == b'[' || !chars[at_index - 1] == b']' {
                    return Err("A global fixed modification should always be enclosed in square brackets '[]'.".to_string());
                }
                let modification =
                    Modification::try_from(&value[index + 2..at_index - 2], &mut ambiguous_lookup)
                        .map(|m| {
                            if let ReturnModification::Defined(m) = m {
                                Ok(m)
                            } else {
                                Err("A global modification cannot be ambiguous".to_string())
                            }
                        })
                        .flat_err()?;
                for aa in value[at_index..end_index - 1].split(',') {
                    global_modifications.push(GlobalModification::Fixed(
                        aa.try_into().map_err(|_| {
                            format!("Could not read as aminoacid in global modification: {aa}")
                        })?,
                        modification.clone(),
                    ));
                }
            } else if &value[index + 1..end_index - 1] == "D" {
                global_modifications.push(GlobalModification::Isotope(Element::H, 2));
            } else {
                let num = &value[index + 1..end_index - 1]
                    .chars()
                    .take_while(char::is_ascii_digit)
                    .collect::<String>();
                let el = &value[index + 1 + num.len()..end_index - 1];
                global_modifications.push(GlobalModification::Isotope(
                    el.try_into().map_err(|_| {
                        format!("Could not read as element in global modification: {el}")
                    })?,
                    num.parse().map_err(|_| {
                        format!("Could not read as isotope number in global modification: {num}")
                    })?,
                ));
            }

            index = end_index;
        }

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
        for (index, preferred, id, localisation_score) in ambiguous_found_positions.iter().copied()
        {
            peptide.sequence[index].possible_modifications.push(
                AmbiguousModification {
                    id,
                    modification: ambiguous_lookup[id].1.as_ref().cloned().ok_or(format!("Ambiguous modification {} did not have a definition for the actual modification", ambiguous_lookup[id].0.as_ref().map_or(id.to_string(), ToString::to_string)))?,
                    localisation_score,
                    group: ambiguous_lookup[id].0.as_ref().map(|n| (n.to_string(), preferred)) });
        }
        peptide.ambiguous_modifications = ambiguous_found_positions
            .iter()
            .copied()
            .group_by(|p| p.2)
            .into_iter()
            .sorted_by(|(key1, _), (key2, _)| key1.cmp(key2))
            .map(|(_, group)| group.into_iter().map(|p| p.0).collect())
            .collect();

        // Check all placement rules
        peptide.apply_global_modifications(&global_modifications);
        peptide.enforce_modification_rules()?;

        Ok(peptide)
    }

    fn enforce_modification_rules(&self) -> Result<(), String> {
        for (index, element) in self.sequence.iter().enumerate() {
            element.enforce_modification_rules(index, self.sequence.len())?;
        }
        Ok(())
    }

    /// Generate all possible patterns for the ambiguous positions (Mass, String:Label)
    /// It always contains at least one pattern (being (base mass, ""))
    fn ambiguous_patterns(
        &self,
        aa_range: impl RangeBounds<usize>,
        aa: &[SequenceElement],
        index: usize,
        base: MolecularFormula,
    ) -> Option<Vec<(MolecularFormula, String)>> {
        let result = self
            .ambiguous_modifications
            .iter()
            .enumerate()
            .fold(vec![Vec::new()], |acc, (id, possibilities)| {
                acc.into_iter()
                    .flat_map(|path| {
                        let mut path_clone = path.clone();
                        let options = possibilities
                            .iter()
                            .filter(|pos| aa_range.contains(pos))
                            .map(move |pos| {
                                let mut new = path.clone();
                                new.push((id, *pos));
                                new
                            });
                        options.chain(
                            possibilities
                                .iter()
                                .find(|pos| !aa_range.contains(pos))
                                .map(move |pos| {
                                    path_clone.push((id, *pos));
                                    path_clone
                                }),
                        )
                    })
                    .collect()
            })
            .into_iter()
            .map(|pattern| {
                let ambiguous_local = pattern
                    .iter()
                    .filter_map(|(id, pos)| (*pos == index).then_some(id))
                    .collect::<Vec<_>>();
                aa.iter()
                    .enumerate()
                    .fold(Some(MolecularFormula::default()), |acc, (index, aa)| {
                        aa.formula(
                            &pattern
                                .iter()
                                .copied()
                                .filter_map(|(id, pos)| (pos == index).then_some(id))
                                .collect_vec(),
                        )
                        .and_then(|m| acc.map(|a| a + m))
                    })
                    .map(|m| {
                        &base
                            + &m
                            + self.sequence[index]
                                .possible_modifications
                                .iter()
                                .filter_map(|am| {
                                    ambiguous_local
                                        .contains(&&am.id)
                                        .then(|| am.modification.formula())
                                })
                                .sum()
                    })
                    .map(|m| {
                        (
                            m,
                            pattern.iter().fold(String::new(), |acc, (id, pos)| {
                                format!(
                                    "{acc}{}{}@{}",
                                    if acc.is_empty() { "" } else { "," },
                                    &self.sequence[index]
                                        .possible_modifications
                                        .iter()
                                        .find(|am| am.id == *id)
                                        .map_or(String::new(), |v| v
                                            .group
                                            .as_ref()
                                            .map_or(id.to_string(), |g| g.0.clone())),
                                    pos + 1
                                )
                            }),
                        )
                    })
            })
            .collect::<Option<Vec<(MolecularFormula, String)>>>()?;
        if result.is_empty() {
            Some(vec![(base, String::new())])
        } else {
            Some(result)
        }
    }

    pub fn formula(&self) -> Option<MolecularFormula> {
        let mut formula = self.n_term() + self.c_term();
        let mut placed = vec![false; self.ambiguous_modifications.len()];
        for (_, pos) in self.sequence.iter().enumerate() {
            formula += pos.formula_greedy(&mut placed)?;
        }

        Some(formula)
    }

    /// Generate the theoretical fragments for this peptide, with the given maximal charge of the fragments, and the given model.
    ///
    /// # Panics
    /// If `max_charge` outside the range `1..=u64::MAX`.
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
    ) -> Option<Vec<Fragment>> {
        assert!(max_charge.value >= 1.0);
        assert!(max_charge.value <= u64::MAX as f64);

        let mut output = Vec::with_capacity(20 * self.sequence.len() + 75); // Empirically derived required size of the buffer (Derived from Hecklib)
        for index in 0..self.sequence.len() {
            let n_term =
                self.ambiguous_patterns(0..=index, &self.sequence[0..index], index, self.n_term())?;

            let c_term = self.ambiguous_patterns(
                index..self.sequence.len(),
                &self.sequence[index + 1..self.sequence.len()],
                index,
                self.c_term(),
            )?;

            output.append(
                &mut self.sequence[index].aminoacid.fragments(
                    &n_term,
                    &c_term,
                    self.sequence[index]
                        .modifications
                        .iter()
                        .map(Chemical::formula)
                        .sum(),
                    max_charge,
                    index,
                    self.sequence.len(),
                    &model.ions(index, self.sequence.len()),
                ),
            );
        }
        // Generate precursor peak
        output.push(
            Fragment::new(
                self.formula()?,
                Charge::zero(),
                FragmentType::precursor,
                String::new(),
            )
            .with_charge(max_charge),
        );
        Some(output)
    }

    fn apply_global_modifications(&mut self, global_modifications: &[GlobalModification]) {
        for modification in global_modifications {
            match modification {
                GlobalModification::Fixed(aa, modification) => {
                    for seq in self.sequence.iter_mut().filter(|seq| seq.aminoacid == *aa) {
                        seq.modifications.push(modification.clone());
                    }
                }
                GlobalModification::Isotope(_el, _isotope) => (), // TODO: implement
            }
        }
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

    pub fn formula(&self, selected_ambiguous: &[usize]) -> Option<MolecularFormula> {
        if self.aminoacid == AminoAcid::B || self.aminoacid == AminoAcid::Z {
            None
        } else {
            Some(
                self.aminoacid.formula()
                    + self.modifications.iter().map(Chemical::formula).sum()
                    + self
                        .possible_modifications
                        .iter()
                        .filter_map(|m| {
                            selected_ambiguous
                                .contains(&m.id)
                                .then(|| m.modification.formula())
                        })
                        .sum(),
            )
        }
    }

    /// Get the mass of this sequence position while placing each possible modification on the very first place (and updating that fact in `placed`)
    pub fn formula_greedy(&self, placed: &mut [bool]) -> Option<MolecularFormula> {
        if self.aminoacid == AminoAcid::B || self.aminoacid == AminoAcid::Z {
            None
        } else {
            Some(
                self.aminoacid.formula()
                    + self.modifications.iter().map(Chemical::formula).sum()
                    + self
                        .possible_modifications
                        .iter()
                        .filter_map(|m| {
                            (!placed[m.id]).then(|| {
                                placed[m.id] = true;
                                m.modification.formula()
                            })
                        })
                        .sum(),
            )
        }
    }

    /// Get the mass of this sequence position with all possible modifications
    pub fn formula_all(&self) -> Option<MolecularFormula> {
        if self.aminoacid == AminoAcid::B || self.aminoacid == AminoAcid::Z {
            None
        } else {
            Some(
                self.aminoacid.formula()
                    + self.modifications.iter().map(Chemical::formula).sum()
                    + self
                        .possible_modifications
                        .iter()
                        .map(|m| m.modification.formula())
                        .sum(),
            )
        }
    }

    /// Enforce the placement rules of predefined modifications.
    fn enforce_modification_rules(&self, index: usize, length: usize) -> Result<(), String> {
        for modification in &self.modifications {
            if let Modification::Predefined(_, rules, _, _) = modification {
                if !rules.is_empty()
                    && !rules
                        .iter()
                        .any(|rule| rule.is_possible(self.aminoacid, index, length))
                {
                    return Err(format!(
                        "Modification {modification} is not allowed on aminoacid {} index {index}",
                        self.aminoacid.char()
                    ));
                }
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AmbiguousModification {
    pub id: usize,
    pub modification: Modification,
    pub localisation_score: Option<f64>,
    /// If this is a named group contain the name and track if this is the preferred location or not
    pub group: Option<(String, bool)>,
}

pub enum GlobalModification {
    Isotope(Element, isize),
    Fixed(AminoAcid, Modification),
}

#[cfg(test)]
mod tests {
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
        assert_eq!(glycan.formula(), peptide.formula());
    }

    #[test]
    fn parse_labile() {
        let with = Peptide::pro_forma("{Formula:C6H10O5}A").unwrap();
        let without = Peptide::pro_forma("A").unwrap();
        assert_eq!(with.sequence.len(), 1);
        assert_eq!(without.sequence.len(), 1);
        assert_eq!(with.formula(), without.formula());
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
        assert_eq!(with.formula(), without.formula());
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
        assert_eq!(glycan.formula(), peptide.formula());
    }

    #[test]
    fn parse_unimod() {
        let peptide = dbg!(Peptide::pro_forma("Q[U:Gln->pyro-Glu]E[Cation:Na]AA"));
        assert!(peptide.is_ok());
    }
}
