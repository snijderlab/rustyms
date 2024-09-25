use std::{cell::RefCell, collections::HashMap};

use crate::{
    glycan::MonoSaccharide,
    modification::{GnoComposition, Ontology, SimpleModification},
    ontologies::CustomDatabase,
    placement_rule::Position,
    system::{ratio::ppm, Mass, Ratio},
    AminoAcid, Chemical, MassMode, Modification, MolecularFormula, Tolerance, WithinTolerance,
};

use super::LinearPeptide;

impl SimpleModification {
    /// Search matching modification based on what modification is provided. If a mass modification is provided
    /// it returns all modifications with that mass (within the tolerance). If a formula is provided it returns
    /// all modifications with that formula. If a glycan composition is provided it returns all glycans with
    /// that composition. Otherwise it returns the modification itself.
    pub fn search(
        modification: &Self,
        tolerance: Tolerance<Mass>,
        mass_mode: MassMode,
        custom_database: Option<&CustomDatabase>,
    ) -> ModificationSearchResult {
        match modification {
            Self::Mass(mass)
            | Self::Gno {
                composition: GnoComposition::Weight(mass),
                ..
            } => ModificationSearchResult::Mass(
                mass.into_inner(),
                tolerance,
                mass_mode,
                [
                    Ontology::Unimod,
                    Ontology::Psimod,
                    Ontology::Gnome,
                    Ontology::Xlmod,
                    Ontology::Custom,
                ]
                .iter()
                .flat_map(|o| {
                    o.lookup(custom_database)
                        .iter()
                        .map(|(i, n, m)| (*o, *i, n, m))
                })
                .filter(|(_, _, _, m)| {
                    tolerance.within(&mass.into_inner(), &m.formula().mass(mass_mode))
                })
                .map(|(o, i, n, m)| (o, i, n.clone(), m.clone()))
                .collect(),
            ),
            Self::Formula(formula) => ModificationSearchResult::Formula(
                formula.clone(),
                [
                    Ontology::Unimod,
                    Ontology::Psimod,
                    Ontology::Gnome,
                    Ontology::Xlmod,
                    Ontology::Custom,
                ]
                .iter()
                .flat_map(|o| {
                    o.lookup(custom_database)
                        .iter()
                        .map(|(i, n, m)| (*o, *i, n, m))
                })
                .filter(|(_, _, _, m)| *formula == m.formula())
                .map(|(o, i, n, m)| (o, i, n.clone(), m.clone()))
                .collect(),
            ),
            Self::Glycan(glycan)
            | Self::Gno {
                composition: GnoComposition::Composition(glycan),
                ..
            } => {
                let search = MonoSaccharide::search_composition(glycan);
                ModificationSearchResult::Glycan(
                    glycan.clone(),
                    Ontology::Gnome
                        .lookup(custom_database)
                        .iter()
                        .filter(|(_, _, m)| {
                            if let Self::Gno {
                                composition: GnoComposition::Topology(structure),
                                ..
                            } = m
                            {
                                MonoSaccharide::search_composition(&structure.composition())
                                    == *search
                            } else if let Self::Gno {
                                composition: GnoComposition::Composition(composition),
                                ..
                            } = m
                            {
                                MonoSaccharide::search_composition(composition) == *search
                            } else {
                                false
                            }
                        })
                        .map(|(i, n, m)| (Ontology::Gnome, *i, n.clone(), m.clone()))
                        .collect(),
                )
            }
            m => ModificationSearchResult::Single(m.clone()),
        }
    }
}

/// The result of a modification search, see [`SimpleModification::search`].
pub enum ModificationSearchResult {
    /// The modification was already defined
    Single(SimpleModification),
    /// All modifications with the same mass, within the tolerance
    Mass(
        Mass,
        Tolerance<Mass>,
        MassMode,
        Vec<(Ontology, Option<usize>, String, SimpleModification)>,
    ),
    /// All modifications with the same formula
    Formula(
        MolecularFormula,
        Vec<(Ontology, Option<usize>, String, SimpleModification)>,
    ),
    /// All modifications with the same glycan composition
    Glycan(
        Vec<(MonoSaccharide, isize)>,
        Vec<(Ontology, Option<usize>, String, SimpleModification)>,
    ),
}

/// Search for named modifications based on mass and/or chemical formula modifications in a peptide.
/// The struct is intended to be reused if multiple peptides need the same replacement strategy.
#[derive(Debug, Clone, PartialEq)]
pub struct PeptideModificationSearch {
    /// If true forces the closest if there are multiple modifications within tolerance, if there are two as close it will still not provide any name
    force_closest: bool,
    /// If true also searches for named modifications for formula modifications
    replace_formulas: bool,
    /// Allow modifications specified as side chain on terminal amino acids to be redefined as terminal modifications
    allow_terminal_redefinition: bool,
    /// The mass mode of the mass modifications (ProForma defines mono isotopic)
    mass_mode: MassMode,
    /// The tolerance for mass based matching
    tolerance: Tolerance<Mass>,
    /// The list of manually given modifications
    modifications: Vec<SimpleModification>,
    /// The list of ontologies to search from
    ontologies: Vec<Ontology>,
    /// The custom modifications, if defined
    custom_database: Option<CustomDatabase>,
    /// The cache to speed up processing from mod + AA to the replacement mod
    cache: RefCell<
        HashMap<(Position, Option<AminoAcid>, SimpleModification), Option<SimpleModification>>,
    >,
}

impl Default for PeptideModificationSearch {
    fn default() -> Self {
        Self {
            force_closest: false,
            replace_formulas: false,
            allow_terminal_redefinition: true,
            mass_mode: MassMode::Monoisotopic,
            tolerance: Tolerance::Relative(Ratio::new::<ppm>(10.0).into()),
            modifications: Vec::new(),
            ontologies: Vec::new(),
            custom_database: None,
            cache: RefCell::new(HashMap::new()),
        }
    }
}

impl PeptideModificationSearch {
    /// Search in the specified list of modifications
    pub fn in_modifications(modifications: Vec<SimpleModification>) -> Self {
        Self {
            modifications,
            ..Self::default()
        }
    }

    /// Search in the given ontologies. Do not forget to add [`Ontology::Custom`] if you want to
    /// allow finding modification in the custom database.
    pub fn in_ontologies(
        ontologies: Vec<Ontology>,
        custom_database: Option<CustomDatabase>,
    ) -> Self {
        Self {
            ontologies,
            custom_database,
            ..Self::default()
        }
    }

    /// Set the tolerance of matches, default is 10 ppm
    pub fn tolerance(self, tolerance: Tolerance<Mass>) -> Self {
        Self {
            tolerance,
            cache: RefCell::new(HashMap::new()),
            ..self
        }
    }

    /// Set the mass mode, all mass modifications will be interpreted as this mode, the default is [`MassMode::MonoIsotopic`]
    pub fn mass_mode(self, mass_mode: MassMode) -> Self {
        Self {
            mass_mode,
            cache: RefCell::new(HashMap::new()),
            ..self
        }
    }

    /// Set the strategy for handling one modification with multiple matches, on default (false) it
    /// will not provide any named modification if multiple are within the tolerance, on true the
    /// closest modification will be picked. If multiple modifications are just as close no
    /// modification will be picked.
    pub fn force_closest(self, force_closest: bool) -> Self {
        Self {
            force_closest,
            cache: RefCell::new(HashMap::new()),
            ..self
        }
    }

    /// Also allow formula modifications to be replaced, defaults to false
    pub fn replace_formulas(self, replace_formulas: bool) -> Self {
        Self {
            replace_formulas,
            cache: RefCell::new(HashMap::new()),
            ..self
        }
    }

    /// Allow modifications on the side chains of terminal amino acids to be redefined as terminal
    /// modifications, default true. It is very common to see such definitions as `Q[-17.02655]AA`
    /// which is a pyroglutamic acid on the Q at the N terminus, these are supposed to be defined
    /// as `[Gln->pyro-glu]-QAA` in ProForma.
    pub fn allow_terminal_redefinition(self, allow_terminal_redefinition: bool) -> Self {
        Self {
            allow_terminal_redefinition,
            cache: RefCell::new(HashMap::new()),
            ..self
        }
    }

    /// Search for modifications that can be replaced by named modifications in this peptide.
    #[allow(
        clippy::similar_names,
        clippy::missing_panics_doc,
        clippy::needless_pass_by_ref_mut
    )] // unwrap is controlled, and the mut is needed for the cache
    pub fn search<Complexity>(
        &mut self,
        mut peptide: LinearPeptide<Complexity>,
    ) -> LinearPeptide<Complexity> {
        let check_matches =
            |in_place: &SimpleModification, provided: &SimpleModification| match in_place {
                SimpleModification::Mass(mass) => self
                    .tolerance
                    .within(&mass.into_inner(), &provided.formula().mass(self.mass_mode)),
                SimpleModification::Formula(formula) if self.replace_formulas => {
                    *formula == provided.formula()
                }
                _ => false,
            };

        let find_replacement_uncached =
            |position: Position, aminoacid: Option<AminoAcid>, in_place: &SimpleModification| {
                let options: Vec<_> = if self.modifications.is_empty() {
                    self.ontologies
                        .iter()
                        .flat_map(|o| o.lookup(self.custom_database.as_ref()))
                        .filter(|modification| {
                            aminoacid.map_or(true, |aa| {
                                modification.2.is_possible_aa(aa, position).any_possible()
                            })
                        })
                        .filter(|modification| check_matches(in_place, &modification.2))
                        .map(|(_, _, modification)| modification)
                        .collect()
                } else {
                    self.modifications
                        .iter()
                        .filter(|modification| {
                            aminoacid.map_or(true, |aa| {
                                modification.is_possible_aa(aa, position).any_possible()
                            })
                        })
                        .filter(|modification| check_matches(in_place, modification))
                        .collect()
                };
                match options.len() {
                    0 => None,
                    1 => Some(options[0].clone()),
                    _ if !self.force_closest => None,
                    _ => {
                        let distances: Vec<_> = options
                            .iter()
                            .map(|m| {
                                in_place
                                    .formula()
                                    .mass(self.mass_mode)
                                    .ppm(m.formula().mass(self.mass_mode))
                                    .value
                            })
                            .collect();
                        let max: f64 = distances
                            .iter()
                            .copied()
                            .max_by(|a, b| (*a).total_cmp(b))
                            .unwrap(); // Guaranteed to always have a value
                        let filtered: Vec<_> = distances
                            .iter()
                            .copied()
                            .enumerate()
                            .filter(|(_, d)| *d == max)
                            .collect();
                        if filtered.len() == 1 {
                            Some(options[filtered[0].0].clone())
                        } else {
                            None
                        }
                    }
                }
            };

        let find_replacement =
            |position: Position, aminoacid: Option<AminoAcid>, in_place: &SimpleModification| {
                if matches!(in_place, SimpleModification::Mass(_))
                    || self.replace_formulas && matches!(in_place, SimpleModification::Formula(_))
                {
                    self.cache
                        .borrow_mut()
                        .entry((position, aminoacid, in_place.clone()))
                        .or_insert_with(|| find_replacement_uncached(position, aminoacid, in_place))
                        .clone()
                } else {
                    None
                }
            };

        let find_replacement_all_positions =
            |n_term: bool,
             c_term: bool,
             aminoacid: Option<AminoAcid>,
             in_place: &SimpleModification| {
                if !self.allow_terminal_redefinition || (!n_term && !c_term) {
                    find_replacement(Position::Anywhere, aminoacid, in_place)
                        .map(|r| (r, Position::Anywhere))
                } else if n_term && c_term {
                    find_replacement(Position::Anywhere, aminoacid, in_place)
                        .map(|r| (r, Position::Anywhere))
                        .or_else(|| {
                            find_replacement(Position::AnyNTerm, aminoacid, in_place)
                                .map(|r| (r, Position::AnyNTerm))
                        })
                        .or_else(|| {
                            find_replacement(Position::AnyCTerm, aminoacid, in_place)
                                .map(|r| (r, Position::AnyCTerm))
                        })
                } else if n_term {
                    find_replacement(Position::Anywhere, aminoacid, in_place)
                        .map(|r| (r, Position::Anywhere))
                        .or_else(|| {
                            find_replacement(Position::AnyNTerm, aminoacid, in_place)
                                .map(|r| (r, Position::AnyNTerm))
                        })
                } else {
                    // The case when c_term
                    find_replacement(Position::Anywhere, aminoacid, in_place)
                        .map(|r| (r, Position::Anywhere))
                        .or_else(|| {
                            find_replacement(Position::AnyCTerm, aminoacid, in_place)
                                .map(|r| (r, Position::AnyCTerm))
                        })
                }
            };

        let find_replacement_modification =
            |position: Position, aminoacid: Option<AminoAcid>, in_place: &Modification| {
                match in_place {
                    Modification::Simple(simple) => {
                        find_replacement(position, aminoacid, simple).map(Modification::Simple)
                    }
                    Modification::CrossLink { .. } => None, // TODO: potentially the cross-linker could be replaced?
                }
            };

        // Start with N and C terminal mods
        let mut n_term = peptide.get_n_term().cloned().and_then(|m| {
            find_replacement_modification(
                Position::AnyNTerm,
                peptide.sequence().first().map(|p| p.aminoacid.aminoacid()),
                &m,
            )
        });
        let mut c_term = peptide.get_c_term().cloned().and_then(|m| {
            find_replacement_modification(
                Position::AnyCTerm,
                peptide.sequence().last().map(|p| p.aminoacid.aminoacid()),
                &m,
            )
        });
        let len = peptide.len();

        // Go over all main stretch mods
        for (index, position) in peptide.sequence_mut().iter_mut().enumerate() {
            let is_n_term = index == 0;
            let is_c_term = index == len;
            let mut remove = None;
            for (i, m) in position.modifications.iter_mut().enumerate() {
                if let Some(simple) = m.simple() {
                    if let Some((replace, location)) = find_replacement_all_positions(
                        is_n_term,
                        is_c_term,
                        Some(position.aminoacid.aminoacid()),
                        simple,
                    ) {
                        if location == Position::AnyNTerm && n_term.is_none() {
                            n_term = Some(Modification::Simple(replace));
                            remove = Some(i);
                        } else if location == Position::AnyCTerm && c_term.is_none() {
                            c_term = Some(Modification::Simple(replace));
                            remove = Some(i);
                        } else if location == Position::Anywhere {
                            *m = Modification::Simple(replace);
                        }
                        // if it can only be a terminal mod but there already is a terminal mod keep it in the original state
                    }
                }
            }
            if let Some(remove) = remove.take() {
                position.modifications.remove(remove);
            }
            for (i, m) in position.possible_modifications.iter_mut().enumerate() {
                if let Some((replace, location)) = find_replacement_all_positions(
                    is_n_term,
                    is_c_term,
                    Some(position.aminoacid.aminoacid()),
                    &m.modification,
                ) {
                    if location == Position::AnyNTerm && n_term.is_none() {
                        n_term = Some(Modification::Simple(replace));
                        remove = Some(i);
                    } else if location == Position::AnyCTerm && c_term.is_none() {
                        c_term = Some(Modification::Simple(replace));
                        remove = Some(i);
                    } else if location == Position::Anywhere {
                        m.modification = replace;
                    }
                    // if it can only be a terminal mod but there already is a terminal mod keep it in the original state
                }
            }
            if let Some(remove) = remove.take() {
                position.possible_modifications.remove(remove);
            }
        }
        for m in peptide.get_labile_mut_inner() {
            if let Some(replace) = find_replacement(Position::Anywhere, None, m) {
                *m = replace;
            }
        }
        peptide.n_term(n_term).c_term(c_term)
    }
}

#[test]
#[allow(clippy::missing_panics_doc)]
fn test_replacement() {
    let mut search = PeptideModificationSearch::in_ontologies(vec![Ontology::Unimod], None)
        .replace_formulas(true);
    let peptide = LinearPeptide::pro_forma("MSFNELT[79.9663]ESNKKSLM[+15.9949]E", None).unwrap();
    let expected = LinearPeptide::pro_forma("MSFNELT[Phospho]ESNKKSLM[Oxidation]E", None).unwrap();
    assert_eq!(search.search(peptide), expected);
    let peptide = LinearPeptide::pro_forma("Q[-17.02655]NKKSLM[+15.9949]E", None).unwrap();
    let expected = LinearPeptide::pro_forma("[Gln->pyro-glu]-QNKKSLM[Oxidation]E", None).unwrap();
    assert_eq!(search.search(peptide), expected);
    let peptide = LinearPeptide::pro_forma("M[Formula:O1]KSLM[+15.9949]E", None).unwrap();
    let expected = LinearPeptide::pro_forma("M[Oxidation]KSLM[Oxidation]E", None).unwrap();
    assert_eq!(search.search(peptide), expected);
}
