//! Handle glycan related issues, access provided if you want to work with glycans on your own.

use std::str::FromStr;

use itertools::Itertools;

use crate::{
    fragment::{Fragment, FragmentType, GlycanBreakPos, GlycanPosition},
    molecular_charge::MolecularCharge,
    system::Charge,
    AminoAcid, Model, Multi,
};

use crate::uom::num_traits::Zero;

include!("shared/glycan.rs");

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

    /// # Errors
    /// Return an Err if the format is not correct
    fn parse_internal(line: &str, range: Range<usize>) -> Result<(Self, usize), CustomError> {
        // Parse at the start the first recognised glycan name
        for name in glycan_parse_list() {
            if line[range.clone()].starts_with(&name.0) {
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
                            sugar: name.1.clone(),
                            branches,
                        },
                        end + 1,
                    ))
                } else {
                    Ok((
                        Self {
                            sugar: name.1.clone(),
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
    /// Return the positioned tree and the outer depth.
    /// # Panics
    /// When any of the masses in this glycan cannot be compared see [`f64::partial_cmp`].
    fn internal_pos(
        self,
        inner_depth: usize,
        branch: &[usize],
    ) -> (PositionedGlycanStructure, usize) {
        // Sort the branches on decreasing molecular weight
        let mut branches = self.branches;
        branches.sort_unstable_by(|a, b| {
            b.formula()
                .monoisotopic_mass()
                .partial_cmp(&a.formula().monoisotopic_mass())
                .unwrap()
        });

        // Get the correct branch indices adding a new layer of indices when needed
        let branches: Vec<(PositionedGlycanStructure, usize)> = if branches.len() == 1 {
            branches
                .into_iter()
                .map(|b| b.internal_pos(inner_depth + 1, branch))
                .collect()
        } else {
            branches
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

    // Recursively show the structure of this glycan
    // fn debug_tree(&self) -> String {
    //     if self.branches.is_empty() {
    //         format!("{:?}", self.sugar)
    //     } else {
    //         format!(
    //             "{:?}({})",
    //             self.sugar,
    //             self.branches.iter().map(Self::debug_tree).join(",")
    //         )
    //     }
    // }
}

impl Display for GlycanStructure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.display_tree())
    }
}
// impl std::fmt::Debug for GlycanStructure {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         write!(f, "{}", self.debug_tree())
//     }
// }
/// Rose tree representation of glycan structure
#[derive(Debug, Eq, PartialEq, Clone, Hash, Serialize, Deserialize)]
pub struct PositionedGlycanStructure {
    sugar: MonoSaccharide,
    branches: Vec<PositionedGlycanStructure>,
    inner_depth: usize,
    outer_depth: usize,
    branch: Vec<usize>,
}

impl Chemical for PositionedGlycanStructure {
    fn formula(&self) -> MolecularFormula {
        self.sugar.formula()
            + self
                .branches
                .iter()
                .map(Chemical::formula)
                .sum::<MolecularFormula>()
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
        full_formula: &Multi<MolecularFormula>,
        attachment: (AminoAcid, usize),
    ) -> Vec<Fragment> {
        model
            .glycan_fragmentation
            .as_ref()
            .map_or(vec![], |neutral_losses| {
                // Get all base fragments from this node and all its children
                let base_fragments =
                    self.base_theoretical_fragments(peptide_index, full_formula, attachment);
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
        full_formula: &Multi<MolecularFormula>,
        attachment: (AminoAcid, usize),
    ) -> Vec<Fragment> {
        // Generate the basic single breakage fragments
        let mut base_fragments = vec![Fragment::new(
            self.formula(),
            Charge::zero(),
            peptide_index,
            FragmentType::B(GlycanPosition {
                inner_depth: self.inner_depth,
                series_number: self.outer_depth + 1,
                branch: self.branch.clone(),
                attachment,
            }),
            String::new(),
        )];
        base_fragments.extend(full_formula.iter().map(|f| {
            Fragment::new(
                f - &self.formula(),
                Charge::zero(),
                peptide_index,
                FragmentType::Y(GlycanPosition {
                    inner_depth: self.inner_depth,
                    series_number: self.inner_depth,
                    branch: self.branch.clone(),
                    attachment,
                }),
                String::new(),
            )
        }));
        // Extend with all internal fragments, meaning multiple breaking bonds
        base_fragments.extend(
            self.internal_break_points(attachment)
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
                                attachment,
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
            self.branches.iter().flat_map(|b| {
                b.base_theoretical_fragments(peptide_index, full_formula, attachment)
            }),
        );
        base_fragments
    }

    fn internal_break_points(
        &self,
        attachment: (AminoAcid, usize),
    ) -> Vec<(MolecularFormula, Vec<GlycanBreakPos>)> {
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
                        attachment,
                    })],
                ),
                (
                    MolecularFormula::default(),
                    vec![GlycanBreakPos::Y(GlycanPosition {
                        inner_depth: self.inner_depth,
                        series_number: self.inner_depth,
                        branch: self.branch.clone(),
                        attachment,
                    })],
                ),
            ]
        } else {
            self.branches
                .iter()
                .map(|b| b.internal_break_points(attachment)) // get all previous options
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
                        attachment,
                    })],
                )))
                .collect()
        }
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod test {
    use super::*;

    #[test]
    fn parse_glycan_structure() {
        assert_eq!(
            GlycanStructure::from_str("Hep(Hex)").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::new(BaseSugar::Heptose(None), &[]).with_name("Hep"),
                branches: vec![GlycanStructure {
                    sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                    branches: Vec::new()
                }],
            }
        );
        assert_eq!(
            GlycanStructure::from_str("Hex(Hex,Hep)").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                branches: vec![
                    GlycanStructure {
                        sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                        branches: Vec::new()
                    },
                    GlycanStructure {
                        sugar: MonoSaccharide::new(BaseSugar::Heptose(None), &[]).with_name("Hep"),
                        branches: Vec::new()
                    }
                ],
            }
        );
        assert_eq!(
            GlycanStructure::from_str("Hex(Hex(Hex),Hep)").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                branches: vec![
                    GlycanStructure {
                        sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                        branches: vec![GlycanStructure {
                            sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[])
                                .with_name("Hex"),
                            branches: Vec::new()
                        }]
                    },
                    GlycanStructure {
                        sugar: MonoSaccharide::new(BaseSugar::Heptose(None), &[]).with_name("Hep"),
                        branches: Vec::new()
                    }
                ],
            }
        );
        assert_eq!(
            GlycanStructure::from_str("Hep(Hex(Hex(Hex(Hep),Hex)))").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::new(BaseSugar::Heptose(None), &[]).with_name("Hep"),
                branches: vec![GlycanStructure {
                    sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                    branches: vec![GlycanStructure {
                        sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                        branches: vec![
                            GlycanStructure {
                                sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[])
                                    .with_name("Hex"),
                                branches: vec![GlycanStructure {
                                    sugar: MonoSaccharide::new(BaseSugar::Heptose(None), &[])
                                        .with_name("Hep"),
                                    branches: Vec::new(),
                                }],
                            },
                            GlycanStructure {
                                sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[])
                                    .with_name("Hex"),
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
            &glycan.formula().into(),
            (AminoAcid::N, 0),
        );
        for fragment in &fragments {
            println!("{fragment}");
        }
        assert_eq!(fragments.len(), 31);
    }
}
