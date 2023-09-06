use itertools::Itertools;

use crate::{
    helper_functions::ResultExtensions,
    modification::{AmbiguousModification, GlobalModification, Modification, ReturnModification},
    Charge, Element, Fragment, LinearPeptide, Model, SequenceElement, MolecularFormula,
};

/// A single pro forma entry, can contain multiple peptides
#[derive(Debug, Clone, PartialEq)]
pub enum ComplexPeptide {
    /// A single linear peptide
    Singular(LinearPeptide),
    /// A multimeric spectrum, multiple peptides coexist in a single spectrum indicated with '+' in pro forma
    Multimeric(Vec<LinearPeptide>),
}

impl ComplexPeptide {
    /// [Pro Forma specification](https://github.com/HUPO-PSI/ProForma)
    /// Only supports a subset of the specification (see `proforma_grammar.md` for an overview of what is supported), some functions are not possible to be represented.
    ///
    /// # Errors
    /// It fails when the string is not a valid Pro Forma string, with a minimal error message to help debug the cause.
    #[allow(clippy::too_many_lines)]
    pub fn pro_forma(value: &str) -> Result<Self, String> {
        let mut peptides = Vec::new();
        let mut start = 0;
        let (peptide, tail) = Self::pro_forma_inner(value, start)?;
        start = tail;
        peptides.push(peptide);

        // Parse any following multimeric species
        while start < value.len() {
            let (peptide, tail) = Self::pro_forma_inner(value, start)?;
            peptides.push(peptide);
            start = tail;
        }
        Ok(if peptides.len() > 1 {
            Self::Multimeric(peptides)
        } else {
            Self::Singular(peptides.pop().ok_or("No peptide found")?)
        })
    }

    fn pro_forma_inner(value: &str, mut index: usize) -> Result<(LinearPeptide, usize), String> {
        if value.trim().is_empty() {
            return Err("Peptide sequence is empty".to_string());
        }
        let mut peptide = LinearPeptide::default();
        let chars: &[u8] = value.as_bytes();
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
            match (c_term, chars[index]) {
                (false, b'(') if chars[index + 1] == b'?' && ambiguous_aa.is_none() => {
                    ambiguous_aa = Some(ambiguous_aa_counter);
                    ambiguous_aa_counter += 1;
                    index += 2;
                }
                (false, b')') if ambiguous_aa.is_some() => {
                    ambiguous_aa = None;
                    index += 1;
                }
                (false, b'/') => {
                    let (charge_len, charge) = next_num(chars, index+1, false).ok_or(format!("Invalid charge state number for peptide [index: {}]", index+1))?;
                    if index+1+charge_len < chars.len() && chars[index+1+charge_len] == b'[' {
                        let end_index = next_char(chars, index+2+charge_len, b']').ok_or(format!(
                            "No valid closing delimiter adduct ion [index: {}]", index+2+charge_len
                        ))?;
                        let mut offset = index+2+charge_len;
                        for set in chars[index+2+charge_len..end_index].split(|c| *c == b',') {
                            // num
                            let (count_len, count) = next_num(chars, offset, true).ok_or(format!("Invalid adduct ion count for peptide [index: {offset}]"))?;
                            // element
                            let element_len = 
                                chars[offset+count_len..].iter()
                                .take_while(|c| c.is_ascii_alphabetic())
                                .count();
                            let element: Element = std::str::from_utf8(&chars[offset+count_len..offset+count_len+element_len]).unwrap().try_into().map_err(|_| format!("Invalid element in adduct ion definition [index: {}]", offset+count_len))?;
                            // charge
                            let (_, ion_charge) = next_num(chars, offset+count_len+element_len, true).ok_or(format!("Invalid adduct ion charge for peptide [index: {offset}]"))?;

                            peptide.adduct_ions.push((count, MolecularFormula::new(&[(element, 0, 1), (Element::Electron, 0, -ion_charge as i16)])));

                            offset += set.len() + 1;
                        }
                        index = end_index+1;
                        if index < chars.len() && chars[index] == b'+' {
                            index+=1; // If a peptide in a multimeric definition contains a charge state modification
                        }
                        break; // Nothing else to do, the total charge of the adduct ions should sum up to the charge as provided and this definitively is the last thing in a sequence 
                    } 
                    // If no adduct ions are provided assume it is just protons
                    peptide.adduct_ions.push((charge, MolecularFormula::new(&[(Element::H, 0, 1), (Element::Electron, 0, -1)])));
                    index += charge_len+1;
                    if index < chars.len() && chars[index] == b'+' {
                        index+=1; // If a peptide in a multimeric definition contains a charge state modification
                    }
                    break;
                }
                (c_term, b'[') => {
                    let end_index = next_char(chars, index, b']').ok_or(format!(
                        "No valid closing delimiter aminoacid modification [index: {index}]"
                    ))?;
                    let modification = Modification::try_from(
                        &value[index + 1..end_index],
                        &mut ambiguous_lookup,
                    )?;
                    index = end_index + 1;
                    if c_term {
                        peptide.c_term =
                            Some(if let ReturnModification::Defined(m) = modification {
                                Ok(m)
                            } else {
                                Err("A labile modification cannot be ambiguous".to_string())
                            }?);
                            if index < chars.len() && chars[index] == b'+' {
                                index+=1; // If a peptide in a multimeric definition contains a C terminal modification
                            }
                        break; // TODO: Technically a charge state could follow the C term mods
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
                }
                (false, b'-') => {
                    c_term = true;
                    index += 1;
                }
                (false, b'+') => {
                    // Multimeric spectrum stop for now, remove the plus
                    index += 1;
                    break;
                }
                (false, ch) => {
                    peptide.sequence.push(SequenceElement::new(
                        ch.try_into().map_err(|_| "Invalid Amino Acid code")?,
                        ambiguous_aa,
                    ));
                    index += 1;
                }
                (true, _) => {
                    return Err(
                        format!("A singular hyphen cannot exist ('-'), if this is part of a c-terminus follow the format 'AA-[modification]' [index: {index}]")
                    )
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

        Ok((peptide, index))
    }

    /// Assume there is exactly one peptide in this collection
    /// # Panics
    /// If there are no or multiple peptides.
    pub fn assume_linear(self) -> LinearPeptide {
        match self {
            Self::Singular(pep) => pep,
            _ => panic!("This ComplexPeptide is not a singular linear peptide"),
        }
    }

    /// Get all peptides making up this `ComplexPeptide`, regardless of its type
    fn peptides(&self) -> &[LinearPeptide] {
        match self {
            Self::Singular(pep) => std::slice::from_ref(pep),
            Self::Multimeric(peptides) => peptides,
        }
    }

    /// Generate the theoretical fragments for this peptide collection.
    /// By iteratively adding every fragment to the set and combining ones that are within the model ppm.
    #[allow(clippy::missing_panics_doc)] // Unwrap is guaranteed to never panic
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
    ) -> Option<Vec<Fragment>> {
        let mut base = Vec::new();
        for (index, peptide) in self.peptides().iter().enumerate() {
            for fragment in peptide.generate_theoretical_fragments(max_charge, model, index)? {
                let (closest_fragment, ppm) =
                    base.iter_mut().fold((None, f64::INFINITY), |acc, i| {
                        let ppm = fragment.ppm(i).map_or(f64::INFINITY, |p| p.value);
                        if acc.1 > ppm {
                            (Some(i), ppm)
                        } else {
                            acc
                        }
                    });
                if ppm < model.ppm.value {
                    // TODO: is this the best combination limit?
                    closest_fragment
                        .unwrap()
                        .add_annotation(fragment.annotations[0]);
                } else {
                    base.push(fragment);
                }
            }
        }
        Some(base)
    }
}

/// Get the index of the next copy of the given char
fn next_char(chars: &[u8], start: usize, char: u8) -> Option<usize> {
    for (i, ch) in chars[start..].iter().enumerate() {
        if *ch == char {
            return Some(start + i);
        }
    }
    None
}

/// Get the next number, returns length in bytes and the number
fn next_num(chars: &[u8], mut start: usize, allow_only_sign: bool) -> Option<(usize, isize)> {
    let mut sign = 1;
    let mut sign_set = false;
    if chars[start] == b'-' {
        sign = -1;
        start += 1;
        sign_set = true;
    } else if chars[start] == b'+' {
        start +=1;
        sign_set = true;
    }
    let len = chars[start..].into_iter().take_while(|c| c.is_ascii_digit()).count();
    if len == 0 {
        if allow_only_sign && sign_set {
            Some((1, sign))
        } else {
            None
        }
    } else {
        let num: isize = std::str::from_utf8(&chars[start..start+len]).unwrap().parse().unwrap();
        Some((sign_set as usize + len, sign * num))
    }
}

#[cfg(test)]
mod tests {
    use crate::ComplexPeptide;
    use crate::Element;
    use crate::MolecularFormula;
    use crate::{e, mz, Location, MassOverCharge};

    use super::*;

    #[test]
    fn parse_glycan() {
        let glycan = ComplexPeptide::pro_forma("A[Glycan:Hex]")
            .unwrap()
            .assume_linear();
        let spaces = ComplexPeptide::pro_forma("A[Glycan:    Hex    ]")
            .unwrap()
            .assume_linear();
        assert_eq!(glycan.sequence.len(), 1);
        assert_eq!(spaces.sequence.len(), 1);
        assert_eq!(glycan, spaces);
        let incorrect = ComplexPeptide::pro_forma("A[Glycan:Hec]");
        assert!(incorrect.is_err());
    }

    #[test]
    fn parse_formula() {
        let peptide = ComplexPeptide::pro_forma("A[Formula:C6H10O5]")
            .unwrap()
            .assume_linear();
        let glycan = ComplexPeptide::pro_forma("A[Glycan:Hex]")
            .unwrap()
            .assume_linear();
        assert_eq!(peptide.sequence.len(), 1);
        assert_eq!(glycan.sequence.len(), 1);
        assert_eq!(glycan.formula(), peptide.formula());
    }

    #[test]
    fn parse_labile() {
        let with = ComplexPeptide::pro_forma("{Formula:C6H10O5}A")
            .unwrap()
            .assume_linear();
        let without = ComplexPeptide::pro_forma("A").unwrap().assume_linear();
        assert_eq!(with.sequence.len(), 1);
        assert_eq!(without.sequence.len(), 1);
        assert_eq!(with.formula(), without.formula());
        assert_eq!(with.labile[0].to_string(), "Formula:C6H10O5".to_string());
    }

    #[test]
    fn parse_ambiguous_modification() {
        let with = ComplexPeptide::pro_forma("A[Phospho#g0]A[#g0]")
            .unwrap()
            .assume_linear();
        let without = ComplexPeptide::pro_forma("AA").unwrap().assume_linear();
        assert_eq!(with.sequence.len(), 2);
        assert_eq!(without.sequence.len(), 2);
        assert_eq!(with.sequence[0].possible_modifications.len(), 1);
        assert_eq!(with.sequence[1].possible_modifications.len(), 1);
        assert!(ComplexPeptide::pro_forma("A[#g0]A[#g0]").is_err());
        assert!(ComplexPeptide::pro_forma("A[Phospho#g0]A[Phospho#g0]").is_err());
        assert!(ComplexPeptide::pro_forma("A[Phospho#g0]A[#g0(0.o1)]").is_err());
        assert_eq!(
            ComplexPeptide::pro_forma("A[+12#g0]A[#g0]")
                .unwrap()
                .assume_linear()
                .to_string(),
            "A[+12#g0]A[#g0]".to_string()
        );
        assert_eq!(
            ComplexPeptide::pro_forma("A[#g0]A[+12#g0]")
                .unwrap()
                .assume_linear()
                .to_string(),
            "A[#g0]A[+12#g0]".to_string()
        );
    }

    #[test]
    fn parse_ambiguous_aminoacid() {
        let with = ComplexPeptide::pro_forma("(?AA)C(?A)(?A)")
            .unwrap()
            .assume_linear();
        let without = ComplexPeptide::pro_forma("AACAA").unwrap().assume_linear();
        assert_eq!(with.sequence.len(), 5);
        assert_eq!(without.sequence.len(), 5);
        assert!(with.sequence[0].ambiguous.is_some());
        assert!(with.sequence[1].ambiguous.is_some());
        assert_eq!(with.formula(), without.formula());
        assert_eq!(with.to_string(), "(?AA)C(?A)(?A)".to_string());
    }

    #[test]
    fn parse_hard_tags() {
        let peptide = ComplexPeptide::pro_forma("A[Formula:C6H10O5|INFO:hello world 🦀]")
            .unwrap()
            .assume_linear();
        let glycan = ComplexPeptide::pro_forma(
            "A[info:you can define a tag multiple times|Glycan:Hex|Formula:C6H10O5]",
        )
        .unwrap()
        .assume_linear();
        assert_eq!(peptide.sequence.len(), 1);
        assert_eq!(glycan.sequence.len(), 1);
        assert_eq!(glycan.formula(), peptide.formula());
    }

    #[test]
    fn parse_global() {
        let deuterium = ComplexPeptide::pro_forma("<D>A").unwrap().assume_linear();
        let nitrogen_15 = ComplexPeptide::pro_forma("<15N>A").unwrap().assume_linear();
        assert_eq!(deuterium.sequence.len(), 1);
        assert_eq!(nitrogen_15.sequence.len(), 1);
        // Formula: A + H2O
        assert_eq!(
            deuterium.formula().unwrap(),
            molecular_formula!((2)H 7 C 3 O 2 N 1)
        );
        assert_eq!(
            nitrogen_15.formula().unwrap(),
            molecular_formula!(H 7 C 3 O 2 (15)N 1)
        );
    }

    #[test]
    fn parse_multimeric() {
        let dimeric = ComplexPeptide::pro_forma("A+AA").unwrap();
        let trimeric = dbg!(ComplexPeptide::pro_forma("A+AA-[+2]+AAA").unwrap());
        assert_eq!(dimeric.peptides().len(), 2);
        assert_eq!(dimeric.peptides()[0].len(), 1);
        assert_eq!(dimeric.peptides()[1].len(), 2);
        assert_eq!(trimeric.peptides().len(), 3);
        assert_eq!(trimeric.peptides()[0].len(), 1);
        assert_eq!(trimeric.peptides()[1].len(), 2);
        assert_eq!(trimeric.peptides()[2].len(), 3);
        assert!(trimeric.peptides()[1].c_term.is_some());
    }

    #[test]
    fn parse_unimod() {
        let peptide = dbg!(ComplexPeptide::pro_forma(
            "Q[U:Gln->pyro-Glu]E[Cation:Na]AA"
        ));
        assert!(peptide.is_ok());
    }

    #[test]
    fn dimeric_peptide() {
        // Only generate a single series, easier to reason about
        let test_model = Model::new(
            (Location::SkipN(1), Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            (Location::None, Vec::new()),
            Vec::new(),
            MassOverCharge::new::<mz>(20.0),
        );

        // With two different sequences
        let dimeric = ComplexPeptide::pro_forma("AA+CC").unwrap();
        let fragments = dbg!(dimeric
            .generate_theoretical_fragments(Charge::new::<e>(1.0), &test_model)
            .unwrap());
        assert_eq!(fragments.len(), 4); // aA, aC, pAA, pCC

        // With two identical sequences
        let dimeric = ComplexPeptide::pro_forma("AA+AA").unwrap();
        let fragments = dbg!(dimeric
            .generate_theoretical_fragments(Charge::new::<e>(1.0), &test_model)
            .unwrap());
        assert_eq!(fragments.len(), 2); // aA, pAA
    }

    #[test]
    fn parse_adduct_ions_01() {
        let peptide = ComplexPeptide::pro_forma("A/2[2Na+]+A").unwrap();
        assert_eq!(peptide.peptides().len(), 2);
        assert_eq!(peptide.peptides()[0].adduct_ions, vec![(2, molecular_formula!(Na 1 Electron -1))]);
        assert_eq!(peptide.peptides()[0].sequence, peptide.peptides()[1].sequence);
    }

    #[test]
    fn parse_adduct_ions_02() {
        let peptide = ComplexPeptide::pro_forma("A-[+1]/2[1Na+,+H+]+[+1]-A").unwrap();
        assert_eq!(peptide.peptides().len(), 2);
        assert_eq!(peptide.peptides()[0].adduct_ions, vec![(1, molecular_formula!(Na 1 Electron -1)), (1, molecular_formula!(H 1 Electron -1))]);
        // Check if the C term mod is applied
        assert_eq!(peptide.peptides()[0].sequence[0].formula_all(), peptide.peptides()[1].sequence[0].formula_all());
    }
}
