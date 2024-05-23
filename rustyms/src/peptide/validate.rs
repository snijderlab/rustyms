use itertools::Itertools;
use std::collections::HashMap;

use crate::{
    error::{Context, CustomError},
    fragment::PeptidePosition,
    modification::{AmbiguousModification, CrossLinkName, SimpleModification},
    placement_rule::Position,
    LinearPeptide, Linked, Modification, Peptidoform,
};

use super::GlobalModification;

/// Validate all cross links
/// # Errors
/// If there is a cross link with more then 2 locations. Or if there never is a definition for this cross link.
pub fn cross_links(
    peptides: Vec<LinearPeptide<Linked>>,
    cross_links_found: HashMap<usize, Vec<(usize, usize)>>,
    cross_link_lookup: &[(CrossLinkName, Option<SimpleModification>)],
    line: &str,
) -> Result<Peptidoform, CustomError> {
    let mut peptidoform = Peptidoform(peptides);
    for (id, locations) in cross_links_found {
        let definition = &cross_link_lookup[id];
        if let Some(linker) = &definition.1 {
            match locations.len() {
                0 => {return Err(CustomError::error(
                    "Invalid cross-link",
                    format!("The cross-link named '{}' has no listed locations, this is an internal error please report this", definition.0),
                    Context::full_line(0, line),
                ))},
                1 => (), // TODO: assumed that the modification is already placed so this works out fine (it is not)
                2 => {
                    if !peptidoform.add_cross_link(locations[0], locations[1], linker.clone(), definition.0.clone()) {
                        return Err(CustomError::error(
                            "Invalid cross-link",
                            format!("The cross-link named '{}' cannot be placed according to its location specificities", definition.0),
                            Context::full_line(0, line),
                        ))
                    }
                },
                _ => {return Err(CustomError::error(
                    "Invalid cross-link",
                    format!("The cross-link named '{}' has more than 2 attachment locations, only cross-links spanning two locations are allowed", definition.0),
                    Context::full_line(0, line),
                ))}
            }
        } else {
            let (c, name, description) = if let CrossLinkName::Branch = definition.0 {
                ("MOD", "00134", " N6-glycyl-L-lysine")
            } else {
                ("X", "DSS", "")
            };
            return Err(CustomError::error(
                "Invalid cross-link",
                format!("The cross-link named '{0}' is never defined, for example for {name}{description} define it like: '[{c}:{name}{0}]'", definition.0),
                Context::full_line(0, line),
            ));
        }
    }

    // TODO: create a function that list all peptides you can reach from a given peptide and check if all are connected
    // let mut connected = peptidoform.formulas_inner().1;
    // if peptidoform.peptides().len() > 1 && connected.len() != peptidoform.peptides().len() {
    //     return Err(CustomError::error(
    //         "Unconnected peptidoform",
    //         "Not all peptides in this peptidoform are connected with cross-links or branches, if separate peptides were intended use the chimeric notation `+` instead of the peptidoform notation `//`.",
    //          Context::full_line(0, line),
    //     ));
    // }

    Ok(peptidoform)
}

impl LinearPeptide<Linked> {
    /// Apply a global modification if this is a global isotope modification with invalid isotopes it returns false
    #[must_use]
    pub(super) fn apply_global_modifications(
        &mut self,
        global_modifications: &[GlobalModification],
    ) -> bool {
        let length = self.len();
        for modification in global_modifications {
            match modification {
                GlobalModification::Fixed(pos, aa, modification) => {
                    for (_, seq) in self.sequence.iter_mut().enumerate().filter(|(index, seq)| {
                        pos.is_possible(&PeptidePosition::n(*index, length))
                            && aa.map_or(true, |aa| aa == seq.aminoacid)
                            && modification
                                .is_possible(seq, &PeptidePosition::n(*index, length))
                                .possible()
                    }) {
                        match pos {
                            Position::Anywhere => seq.modifications.push(modification.clone()),
                            Position::AnyNTerm | Position::ProteinNTerm => {
                                self.n_term = Some(
                                    modification
                                        .simple()
                                        .expect(
                                            "Can only put a simple modification on an N terminus",
                                        )
                                        .clone(),
                                );
                            }
                            Position::AnyCTerm | Position::ProteinCTerm => {
                                self.c_term = Some(
                                    modification
                                        .simple()
                                        .expect(
                                            "Can only put a simple modification on a C terminus",
                                        )
                                        .clone(),
                                );
                            }
                        }
                    }
                }
                GlobalModification::Isotope(el, isotope) if el.is_valid(*isotope) => {
                    self.global.push((*el, *isotope));
                }
                GlobalModification::Isotope(..) => return false,
            }
        }
        true
    }

    /// Place all global unknown positions at all possible locations as ambiguous modifications
    pub(super) fn apply_unknown_position_modification(
        &mut self,
        unknown_position_modifications: &[SimpleModification],
    ) {
        for modification in unknown_position_modifications {
            let id = self.ambiguous_modifications.len();
            let length = self.len();
            #[allow(clippy::unnecessary_filter_map)]
            // Side effects so the lint does not apply here
            self.ambiguous_modifications.push(
                (0..length)
                    .filter_map(|i| {
                        if modification
                            .is_possible(&self.sequence[i], &PeptidePosition::n(i, length))
                            .possible()
                        {
                            self.sequence[i]
                                .possible_modifications
                                .push(AmbiguousModification {
                                    id,
                                    modification: modification.clone(),
                                    localisation_score: None,
                                    group: None,
                                });
                            Some(i)
                        } else {
                            None
                        }
                    })
                    .collect(),
            );
        }
    }

    /// Place all ranged unknown positions at all possible locations as ambiguous modifications
    /// # Panics
    /// It panics when information for an ambiguous modification is missing (name/mod).
    pub(super) fn apply_ranged_unknown_position_modification(
        &mut self,
        ranged_unknown_position_modifications: &[(usize, usize, SimpleModification)],
        ambiguous_lookup: &[(Option<String>, Option<SimpleModification>)],
    ) {
        let mut id = ambiguous_lookup.len();
        for (start, end, modification) in ranged_unknown_position_modifications {
            let length = self.len();
            #[allow(clippy::unnecessary_filter_map)]
            // Side effects so the lint does not apply here
            let positions = (*start..=*end)
                .filter_map(|i| {
                    if modification
                        .is_possible(&self.sequence[i], &PeptidePosition::n(i, length))
                        .possible()
                    {
                        self.sequence[i]
                            .possible_modifications
                            .push(AmbiguousModification {
                                id,
                                modification: modification.clone(),
                                localisation_score: None,
                                group: None,
                            });
                        Some(i)
                    } else {
                        None
                    }
                })
                .collect_vec();
            self.ambiguous_modifications[id].extend(positions);
            id += 1;
        }
    }
}

impl<T> LinearPeptide<T> {
    /// # Errors
    /// If a modification rule is broken it returns an error.
    pub(crate) fn enforce_modification_rules(&self) -> Result<(), CustomError> {
        for (position, seq) in self.iter(..) {
            seq.enforce_modification_rules(&position)?;
        }
        Ok(())
    }
}
