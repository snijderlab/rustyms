//! Handle glycan related issues, access provided if you want to work with glycans on your own.

use std::str::FromStr;

use crate::{
    fragment::{Fragment, FragmentType, GlycanBreakPos, GlycanPosition},
    molecular_charge::MolecularCharge,
    system::usize::Charge,
    AminoAcid, Model, Multi, NeutralLoss,
};

use crate::uom::num_traits::Zero;

include!("shared/glycan.rs");

impl MonoSaccharide {
    /// Generate the composition used for searching on glycans
    pub(crate) fn simplify_composition(mut composition: Vec<(Self, isize)>) -> Vec<(Self, isize)> {
        // Sort on monosaccharide
        composition.retain(|el| el.1 != 0);
        composition.sort_unstable_by(|a, b| a.0.cmp(&b.0));

        // Deduplicate
        let mut max = composition.len().saturating_sub(1);
        let mut index = 0;
        while index < max {
            let this = &composition[index];
            let next = &composition[index + 1];
            if this.0 == next.0 {
                composition[index].1 += next.1;
                composition.remove(index + 1);
                max = max.saturating_sub(1);
            } else {
                index += 1;
            }
        }
        composition.retain(|el| el.1 != 0);
        composition
    }

    /// Generate the composition used for searching on glycans
    pub(crate) fn search_composition(
        mut composition: Vec<(Self, isize)>,
    ) -> Vec<(MolecularFormula, isize)> {
        // Sort on monosaccharide
        composition.retain(|el| el.1 != 0);
        let mut composition = composition
            .into_iter()
            .map(|(m, n)| (m.formula(), n))
            .collect_vec();
        composition.sort_unstable_by(|a, b| a.0.cmp(&b.0));

        // Deduplicate
        let mut max = composition.len().saturating_sub(1);
        let mut index = 0;
        while index < max {
            let this = &composition[index];
            let next = &composition[index + 1];
            if this.0 == next.0 {
                composition[index].1 += next.1;
                composition.remove(index + 1);
                max = max.saturating_sub(1);
            } else {
                index += 1;
            }
        }
        composition.retain(|el| el.1 != 0);
        composition
    }

    /// Generate all uncharged diagnostic ions for this monosaccharide.
    /// According to: <https://doi.org/10.1016/j.trac.2018.09.007>.
    fn diagnostic_ions(&self, peptide_index: usize, position: GlycanPosition) -> Vec<Fragment> {
        let base = Fragment::new(
            self.formula(),
            Charge::default(),
            peptide_index,
            FragmentType::diagnostic(crate::fragment::DiagnosticPosition::Glycan(
                position,
                self.clone(),
            )),
            String::new(),
        );
        if matches!(self.base_sugar, BaseSugar::Hexose(_)) && self.substituents.is_empty() {
            vec![
                base.clone(),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(H 2 O 1))),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(H 4 O 2))),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(C 1 H 6 O 3))),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(C 2 H 6 O 3))),
            ]
        } else if matches!(self.base_sugar, BaseSugar::Hexose(_))
            && self.substituents == [GlycanSubstituent::NAcetyl]
        {
            vec![
                base.clone(),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(H 2 O 1))),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(H 4 O 2))),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(C 2 H 4 O 2))),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(C 1 H 6 O 3))),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(C 2 H 6 O 3))),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(C 4 H 8 O 4))),
            ]
        } else if matches!(self.base_sugar, BaseSugar::Nonose)
            && (self.substituents
                == [
                    GlycanSubstituent::Amino,
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acid,
                ]
                || self.substituents
                    == [
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                        GlycanSubstituent::Acid,
                    ])
        {
            // Neu5Ac and Neu5Gc
            vec![
                base.clone(),
                base.with_neutral_loss(&NeutralLoss::Loss(molecular_formula!(H 2 O 1))),
            ]
        } else {
            Vec::new()
        }
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

    /// Get the composition of a `GlycanStructure`. The result is normalised (sorted and deduplicated).
    pub fn composition(&self) -> Vec<(MonoSaccharide, isize)> {
        let composition = self.composition_inner();
        MonoSaccharide::simplify_composition(composition)
    }

    /// Get the composition in monosaccharides of this glycan
    fn composition_inner(&self) -> Vec<(MonoSaccharide, isize)> {
        let mut output = vec![(self.sugar.clone(), 1)];
        output.extend(self.branches.iter().flat_map(Self::composition_inner));
        output
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
        let single_charges = charge_carriers.all_single_charge_options();
        let all_charges = charge_carriers.all_charge_options();
        model.glycan.as_ref().map_or(vec![], |neutral_losses| {
            // Get all base fragments from this node and all its children
            let mut base_fragments = self
                .oxonium_fragments(peptide_index, attachment)
                .into_iter()
                .flat_map(|f| f.with_charges(&single_charges))
                .flat_map(|f| f.with_neutral_losses(neutral_losses))
                .collect_vec();
            // Generate all Y fragments
            base_fragments.extend(
                self.internal_break_points(attachment)
                    .iter()
                    .filter(|(_, bonds)| {
                        bonds.iter().all(|b| !matches!(b, GlycanBreakPos::B(_)))
                            && !bonds.iter().all(|b| matches!(b, GlycanBreakPos::End(_)))
                    })
                    .flat_map(move |(f, bonds)| {
                        full_formula.iter().map(move |full| {
                            Fragment::new(
                                full - self.formula() + f,
                                Charge::zero(),
                                peptide_index,
                                FragmentType::Y(
                                    bonds
                                        .iter()
                                        .filter(|b| !matches!(b, GlycanBreakPos::End(_)))
                                        .map(GlycanBreakPos::position)
                                        .cloned()
                                        .collect(),
                                ),
                                String::new(),
                            )
                        })
                    })
                    .flat_map(|f| f.with_charges(&all_charges))
                    .flat_map(|f| f.with_neutral_losses(neutral_losses)),
            );
            // Generate all diagnostic ions
            base_fragments.extend(
                self.diagnostic_ions(peptide_index, attachment)
                    .into_iter()
                    .flat_map(|f| f.with_charges(&single_charges)),
            );
            base_fragments
        })
    }

    /// Get uncharged diagnostic ions from all positions
    fn diagnostic_ions(
        &self,
        peptide_index: usize,
        attachment: (AminoAcid, usize),
    ) -> Vec<Fragment> {
        let mut output = self
            .sugar
            .diagnostic_ions(peptide_index, self.position(attachment));
        output.extend(
            self.branches
                .iter()
                .flat_map(|b| b.diagnostic_ions(peptide_index, attachment)),
        );

        output
    }

    /// Generate all fragments without charge and neutral loss options
    fn oxonium_fragments(
        &self,
        peptide_index: usize,
        attachment: (AminoAcid, usize),
    ) -> Vec<Fragment> {
        // Generate the basic single breakage B fragments
        let mut base_fragments = vec![Fragment::new(
            self.formula(),
            Charge::zero(),
            peptide_index,
            FragmentType::B(self.position(attachment)),
            String::new(),
        )];
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
                        [b, vec![GlycanBreakPos::B(self.position(attachment))]].concat(),
                    )
                })
                .map(|(formula, breakages)| {
                    Fragment::new(
                        formula,
                        Charge::zero(),
                        peptide_index,
                        FragmentType::Oxonium(breakages),
                        String::new(),
                    )
                }),
        );
        // Extend with the theoretical fragments for all branches of this position
        base_fragments.extend(
            self.branches
                .iter()
                .flat_map(|b| b.oxonium_fragments(peptide_index, attachment)),
        );
        base_fragments
    }

    /// All possible bonds that can be broken and the molecular formula that would be held over if these bonds all broke and the broken off parts are lost.
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
                    vec![GlycanBreakPos::End(self.position(attachment))],
                ),
                (
                    MolecularFormula::default(),
                    vec![GlycanBreakPos::Y(self.position(attachment))],
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
                    vec![GlycanBreakPos::Y(self.position(attachment))],
                )))
                .collect()
        }
    }

    fn position(&self, attachment: (AminoAcid, usize)) -> GlycanPosition {
        GlycanPosition {
            inner_depth: self.inner_depth,
            series_number: self.outer_depth + 1,
            branch: self.branch.clone(),
            attachment,
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
    fn correct_masses() {
        let (sugar, _) = MonoSaccharide::from_short_iupac("Neu5Ac", 0, 0).unwrap();
        dbg!(&sugar);

        assert_eq!(sugar.formula(), molecular_formula!(C 11 H 17 N 1 O 8));
    }

    #[test]
    fn correct_structure_g43728nl() {
        // Furanoses added for error detection
        let structure = GlycanStructure::from_short_iupac(
            "Neu5Ac(?2-?)Galf(?1-?)GlcNAc(?1-?)Man(?1-?)[Galf(?1-?)GlcNAc(?1-?)Man(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc", 
            0..101,
            0
        )
        .unwrap();

        assert_eq!(
            structure.to_string(),
            "HexNAc(HexNAc(Hex(Hex(HexNAc(Hexf(NonNAAc))),Hex(HexNAc(Hexf)))))"
        );
    }

    #[test]
    fn correct_structure_g36564am() {
        let structure = GlycanStructure::from_short_iupac(
            "Gal(?1-?)GlcNAc(?1-?)Man(?1-?)[GlcNAc(?1-?)Man(?1-?)][GlcNAc(?1-?)]Man(?1-?)GlcNAc",
            0..82,
            0,
        )
        .unwrap();

        assert_eq!(
            structure.to_string(),
            "HexNAc(Hex(Hex(HexNAc(Hex)),Hex(HexNAc),HexNAc))"
        );
    }

    #[test]
    fn correct_structure_g67881ee() {
        // Fully specified version of g36564am
        // Furanoses added for error detection
        let structure = GlycanStructure::from_short_iupac(
            "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][Galf(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-",
            0..87,
            0,
        )
        .unwrap();

        assert_eq!(
            structure.to_string(),
            "HexNAc(Hex(Hex(HexNAc),HexNAc,Hex(HexNAc(Hexf))))"
        );
    }

    #[test]
    fn correct_structure_g11771hd() {
        let structure = GlycanStructure::from_short_iupac(
            "GlcNAc(?1-?)[GlcNAc(?1-?)]Man(?1-?)[Man(?1-?)Man(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc(?1-",
            0..86,
            0,
        )
        .unwrap();

        assert_eq!(
            structure.to_string(),
            "HexNAc(HexNAc(Hex(Hex(HexNAc,HexNAc),Hex(Hex))))"
        );
    }
}
