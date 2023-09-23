use std::{fmt::Display, ops::Range, str::FromStr};

use itertools::Itertools;

use crate::{
    error::{Context, CustomError},
    molecular_charge::MolecularCharge,
    Chemical, MolecularFormula,
};

include!("shared/glycan.rs");

/// Rose tree representation of glycan structure
#[derive(Eq, PartialEq, Clone, Hash)]
pub struct GlycanStructure {
    sugar: MonoSaccharide,
    branches: Vec<GlycanStructure>,
}

impl Chemical for GlycanStructure {
    fn formula(&self) -> MolecularFormula {
        self.sugar.formula() + self.branches.iter().map(formula::Chemical::formula).sum()
    }
}

impl FromStr for GlycanStructure {
    type Err = CustomError;
    /// Parse a textual structure representation of a glycan (outside Pro Forma format)
    /// Example: Hex(Hex(HexNAc)) => Hex-Hex-HexNAc (linear)
    /// Example: Hex(Fuc,Hex(HexNAc,Hex(HexNAc)))
    ///          =>  Hex-Hex-HexNAc
    ///              └Fuc  └Hex-HexNAc
    /// # Errors
    /// Return an Err if the format is not correct
    fn from_str(line: &str) -> Result<Self, CustomError> {
        Self::parse(line, 0..line.len())
    }
}

impl GlycanStructure {
    /// Parse a textual structure representation of a glycan (outside Pro Forma format)
    /// Example: Hex(Hex(HexNAc)) => Hex-Hex-HexNAc (linear)
    /// Example: Hex(Fuc,Hex(HexNAc,Hex(HexNAc)))
    ///          =>  Hex-Hex-HexNAc
    ///              └Fuc  └Hex-HexNAc
    /// # Errors
    /// Return an Err if the format is not correct
    pub fn parse(line: &str, range: Range<usize>) -> Result<Self, CustomError> {
        Self::parse_internal(line, range).map(|(g, _)| g)
    }

    fn parse_internal(line: &str, range: Range<usize>) -> Result<(Self, usize), CustomError> {
        // Parse at the start the first recognised glycan name
        for name in GLYCAN_PARSE_LIST {
            if line[range.clone()].starts_with(name.0) {
                // If the name is followed by a bracket parse a list of branches
                return if line.as_bytes()[range.start + name.0.len()] == b'(' {
                    // Find the end of this list
                    let end = end_of_enclosure(
                        line.as_bytes(),
                        range.start + name.0.len() + 1,
                        b'(',
                        b')',
                    )
                    .ok_or_else(|| {
                        CustomError::error(
                            "Invalid glycan branch",
                            "No valid closing delimiter",
                            Context::line(0, line, range.start + name.0.len(), 1),
                        )
                    })?;
                    // Parse the first branch
                    let mut index = range.start + name.0.len() + 1;
                    let mut branches = Vec::new();
                    let (glycan, pos) = Self::parse_internal(line, index..end)?;
                    index = pos;
                    branches.push(glycan);
                    // Keep parsing until the end of this branch level (until ')' is reached)
                    while index < end {
                        if line.as_bytes()[index] != b',' {
                            return Err(CustomError::error(
                                "Invalid glycan structure",
                                "Branches should be separated by commas ','",
                                Context::line(0, line, index, 1),
                            ));
                        }
                        index += 1;
                        let (glycan, pos) = Self::parse_internal(line, index..end)?;
                        branches.push(glycan);
                        index = pos;
                    }
                    Ok((
                        Self {
                            sugar: name.1,
                            branches,
                        },
                        end + 1,
                    ))
                } else {
                    Ok((
                        Self {
                            sugar: name.1,
                            branches: Vec::new(),
                        },
                        range.start + name.0.len(),
                    ))
                };
            }
        }
        Err(CustomError::error(
            "Could not parse glycan structure",
            "Could not parse the following part",
            Context::line(0, line, range.start, range.len()),
        ))
    }

    /// Annotate all positions in this tree with all positions
    pub fn determine_positions(self) -> PositionedGlycanStructure {
        self.internal_pos(0, &[]).0
    }

    /// Given the inner depth determine the correct positions and branch ordering
    /// Return the positioned tree and the outer depth
    fn internal_pos(
        mut self,
        inner_depth: usize,
        branch: &[usize],
    ) -> (PositionedGlycanStructure, usize) {
        // Sort the branches on decreasing molecular weight
        self.branches.sort_unstable_by(|a, b| {
            b.formula()
                .monoisotopic_mass()
                .unwrap()
                .partial_cmp(&a.formula().monoisotopic_mass().unwrap())
                .unwrap()
        });

        // Get the correct branch indices adding a new layer of indices when needed
        let branches: Vec<(PositionedGlycanStructure, usize)> = if self.branches.len() == 1 {
            self.branches
                .into_iter()
                .map(|b| b.internal_pos(inner_depth + 1, branch))
                .collect()
        } else {
            self.branches
                .into_iter()
                .enumerate()
                .map(|(i, b)| {
                    let mut new_branch = branch.to_vec();
                    new_branch.push(i);
                    b.internal_pos(inner_depth + 1, &new_branch)
                })
                .collect()
        };

        let outer_depth = branches.iter().map(|b| b.1).max().unwrap_or(0);
        (
            PositionedGlycanStructure {
                sugar: self.sugar,
                branches: branches.into_iter().map(|b| b.0).collect(),
                branch: branch.to_vec(),
                inner_depth,
                outer_depth,
            },
            outer_depth + 1,
        )
    }

    /// Get the maximal outer depth of all branches for this location
    fn outer_depth(&self) -> usize {
        self.branches
            .iter()
            .map(Self::outer_depth)
            .max()
            .unwrap_or(1)
    }

    /// Recursively show the structure of this glycan
    fn display_tree(&self) -> String {
        if self.branches.is_empty() {
            self.sugar.to_string()
        } else {
            format!(
                "{}({})",
                self.sugar,
                self.branches.iter().map(Self::display_tree).join(",")
            )
        }
    }
}

impl Display for GlycanStructure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.display_tree())
    }
}
impl std::fmt::Debug for GlycanStructure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.display_tree())
    }
}
/// Rose tree representation of glycan structure
#[derive(Debug, Eq, PartialEq, Clone, Hash)]
pub struct PositionedGlycanStructure {
    sugar: MonoSaccharide,
    branches: Vec<PositionedGlycanStructure>,
    inner_depth: usize,
    outer_depth: usize,
    branch: Vec<usize>,
}

impl Chemical for PositionedGlycanStructure {
    fn formula(&self) -> MolecularFormula {
        self.sugar.formula() + self.branches.iter().map(formula::Chemical::formula).sum()
    }
}

impl PositionedGlycanStructure {
    /// Generate all theoretical fragments for this glycan
    /// * `full_formula` the total formula of the whole peptide + glycan
    pub fn generate_theoretical_fragments(
        &self,
        model: &Model,
        peptide_index: usize,
        charge_carriers: &MolecularCharge,
        full_formula: &MolecularFormula,
    ) -> Vec<Fragment> {
        model
            .glycan_fragmentation
            .as_ref()
            .map_or(vec![], |neutral_losses| {
                // Get all base fragments from this node and all its children
                let base_fragments = self.base_theoretical_fragments(peptide_index, full_formula);
                // Apply all neutral losses and all charge options
                let charge_options = charge_carriers.all_charge_options();
                base_fragments
                    .into_iter()
                    .flat_map(|f| f.with_neutral_losses(neutral_losses))
                    .flat_map(|f| charge_options.iter().map(move |c| f.with_charge(c)))
                    .collect()
            })
    }

    /// Generate all fragments without charge and neutral loss options
    fn base_theoretical_fragments(
        &self,
        peptide_index: usize,
        full_formula: &MolecularFormula,
    ) -> Vec<Fragment> {
        // Generate the basic single breakage fragments
        let mut base_fragments = vec![
            Fragment::new(
                self.formula(),
                Charge::zero(),
                peptide_index,
                FragmentType::B(GlycanPosition {
                    inner_depth: self.inner_depth,
                    series_number: self.outer_depth + 1,
                    branch: self.branch.clone(),
                }),
                String::new(),
            ),
            // Fragment::new(
            //     self.formula() + molecular_formula!(O 1 H 2),
            //     Charge::zero(),
            //     peptide_index,
            //     FragmentType::C(GlycanPosition {
            //         inner_depth: self.inner_depth,
            //         series_number: self.outer_depth,
            //         branch: self.branch.clone(),
            //     }),
            //     String::new(),
            // ),
            Fragment::new(
                full_formula - &self.formula(),
                Charge::zero(),
                peptide_index,
                FragmentType::Y(GlycanPosition {
                    inner_depth: self.inner_depth,
                    series_number: self.inner_depth,
                    branch: self.branch.clone(),
                }),
                String::new(),
            ),
            // Fragment::new(
            //     total_glycan - &self.formula() - molecular_formula!(O 1 H 2),
            //     Charge::zero(),
            //     peptide_index,
            //     FragmentType::Z(GlycanPosition {
            //         inner_depth: self.inner_depth,
            //         series_number: self.inner_depth,
            //         branch: self.branch.clone(),
            //     }),
            //     String::new(),
            // ),
            // TODO: for A, X, C, Z just ignore for now?
        ];
        // Extend with all internal fragments, meaning multiple breaking bonds
        base_fragments.extend(
            self.internal_break_points()
                .into_iter()
                .filter(|(_, breakages)| {
                    !breakages
                        .iter()
                        .all(|b| matches!(b, GlycanBreakPos::End(_)))
                })
                .filter(|(m, _)| *m != MolecularFormula::default())
                .map(|(m, b)| {
                    (
                        m,
                        [
                            b,
                            vec![GlycanBreakPos::B(GlycanPosition {
                                inner_depth: self.inner_depth,
                                series_number: self.outer_depth + 1,
                                branch: self.branch.clone(),
                            })],
                        ]
                        .concat(),
                    )
                })
                .map(|(formula, breakages)| {
                    Fragment::new(
                        formula,
                        Charge::zero(),
                        peptide_index,
                        FragmentType::InternalGlycan(breakages),
                        String::new(),
                    )
                }),
        );
        // Extend with the theoretical fragments for all branches of this position
        base_fragments.extend(
            self.branches
                .iter()
                .flat_map(|b| b.base_theoretical_fragments(peptide_index, full_formula)),
        );
        base_fragments
    }

    fn internal_break_points(&self) -> Vec<(MolecularFormula, Vec<GlycanBreakPos>)> {
        // Find every internal fragment ending at this bond (in a B breakage) (all bonds found are Y breakages and endings)
        // Walk through all branches and determine all possible breakages
        if self.branches.is_empty() {
            vec![
                (
                    self.formula(),
                    vec![GlycanBreakPos::End(GlycanPosition {
                        inner_depth: self.inner_depth,
                        series_number: self.inner_depth,
                        branch: self.branch.clone(),
                    })],
                ),
                (
                    MolecularFormula::default(),
                    vec![GlycanBreakPos::Y(GlycanPosition {
                        inner_depth: self.inner_depth,
                        series_number: self.inner_depth,
                        branch: self.branch.clone(),
                    })],
                ),
            ]
        } else {
            self.branches
                .iter()
                .map(Self::internal_break_points) // get all previous options
                .fold(Vec::new(), |accumulator, branch_options| {
                    if accumulator.is_empty() {
                        branch_options
                    } else {
                        let mut new_accumulator = Vec::new();
                        for base in &accumulator {
                            for option in &branch_options {
                                new_accumulator.push((
                                    &option.0 + &base.0,
                                    [option.1.clone(), base.1.clone()].concat(),
                                ));
                            }
                        }
                        new_accumulator
                    }
                })
                .into_iter()
                .map(|(m, b)| (m + self.sugar.formula(), b))
                .chain(std::iter::once((
                    // add the option of it breaking here
                    MolecularFormula::default(),
                    vec![GlycanBreakPos::Y(GlycanPosition {
                        inner_depth: self.inner_depth,
                        series_number: self.inner_depth,
                        branch: self.branch.clone(),
                    })],
                )))
                .collect()
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    #[allow(clippy::float_cmp)] // Already handled in a way clippy does not recognise
    fn mass_glycan() {
        assert_eq!(
            1445.0,
            (MonoSaccharide::Hex.formula().average_weight().unwrap() * 3.0
                + MonoSaccharide::HexNAc.formula().average_weight().unwrap() * 4.0
                + MonoSaccharide::Fuc.formula().average_weight().unwrap())
            .value
            .round()
        );
    }

    #[test]
    fn parse_glycan_structure() {
        assert_eq!(
            GlycanStructure::from_str("HexNAc(Hex)").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::HexNAc,
                branches: vec![GlycanStructure {
                    sugar: MonoSaccharide::Hex,
                    branches: Vec::new()
                }],
            }
        );
        assert_eq!(
            GlycanStructure::from_str("Hex(Hex,HexNAc)").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::Hex,
                branches: vec![
                    GlycanStructure {
                        sugar: MonoSaccharide::Hex,
                        branches: Vec::new()
                    },
                    GlycanStructure {
                        sugar: MonoSaccharide::HexNAc,
                        branches: Vec::new()
                    }
                ],
            }
        );
        assert_eq!(
            GlycanStructure::from_str("Hex(Hex(Hex),HexNAc)").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::Hex,
                branches: vec![
                    GlycanStructure {
                        sugar: MonoSaccharide::Hex,
                        branches: vec![GlycanStructure {
                            sugar: MonoSaccharide::Hex,
                            branches: Vec::new()
                        }]
                    },
                    GlycanStructure {
                        sugar: MonoSaccharide::HexNAc,
                        branches: Vec::new()
                    }
                ],
            }
        );
        assert_eq!(
            GlycanStructure::from_str("HexNAc(Hex(Hex(Hex(HexNAc),Hex)))").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::HexNAc,
                branches: vec![GlycanStructure {
                    sugar: MonoSaccharide::Hex,
                    branches: vec![GlycanStructure {
                        sugar: MonoSaccharide::Hex,
                        branches: vec![
                            GlycanStructure {
                                sugar: MonoSaccharide::Hex,
                                branches: vec![GlycanStructure {
                                    sugar: MonoSaccharide::HexNAc,
                                    branches: Vec::new(),
                                }],
                            },
                            GlycanStructure {
                                sugar: MonoSaccharide::Hex,
                                branches: Vec::new(),
                            },
                        ],
                    }],
                }],
            }
        );
    }

    #[test]
    fn internal_breakages() {
        let glycan = GlycanStructure::from_str("HexNAc(Hex(Hex(Hex(HexNAc),Hex)))")
            .unwrap()
            .determine_positions();
        let fragments = glycan.generate_theoretical_fragments(
            &Model {
                glycan_fragmentation: Some(Vec::new()),
                ..Model::none()
            },
            0,
            &MolecularCharge::proton(1),
            &glycan.formula(),
        );
        for fragment in &fragments {
            println!("{fragment}");
        }
        assert_eq!(fragments.len(), 31);
    }
}
