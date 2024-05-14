//! Handle monosaccharides

use crate::{
    fragment::{Fragment, FragmentType, GlycanPosition},
    system::usize::Charge,
    NeutralLoss,
};

include!("../shared/glycan.rs");
include!("../shared/glycan_lists.rs");

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
    pub(crate) fn diagnostic_ions(
        &self,
        peptidoform_index: usize,
        peptide_index: usize,
        position: GlycanPosition,
    ) -> Vec<Fragment> {
        let base = Fragment::new(
            self.formula(),
            Charge::default(),
            peptidoform_index,
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

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use super::*;

    #[test]
    fn pro_forma_compliance() {
        let cases = &[
            ("hep", molecular_formula!(H 12 C 7 O 6)),
            ("phosphate", molecular_formula!(H 1 O 3 P 1)),
            ("a-hex", molecular_formula!(H 8 C 6 O 6)),
            ("sug", molecular_formula!(H 2 C 2 O 1)),
            ("hexn", molecular_formula!(H 11 C 6 N 1 O 4)),
            ("pen", molecular_formula!(H 8 C 5 O 4)),
            ("tet", molecular_formula!(H 6 C 4 O 3)),
            ("hexp", molecular_formula!(H 11 C 6 O 8 P 1)),
            ("neu5ac", molecular_formula!(H 17 C 11 N 1 O 8)),
            ("non", molecular_formula!(H 16 C 9 O 8)),
            ("hexnac(s)", molecular_formula!(H 13 C 8 N 1 O 8 S 1)),
            ("dec", molecular_formula!(H 18 C 10 O 9)),
            ("en,a-hex", molecular_formula!(H 6 C 6 O 5)),
            ("neu5gc", molecular_formula!(H 17 C 11 N 1 O 9)),
            ("neu", molecular_formula!(H 15 C 9 N 1 O 7)),
            ("hexnac", molecular_formula!(H 13 C 8 N 1 O 5)),
            ("fuc", molecular_formula!(H 10 C 6 O 4)),
            ("hexns", molecular_formula!(H 11 C 6 N 1 O 7 S 1)),
            ("tri", molecular_formula!(H 4 C 3 O 2)),
            ("oct", molecular_formula!(H 14 C 8 O 7)),
            ("sulfate", molecular_formula!(O 3 S 1)),
            ("d-hex", molecular_formula!(H 10 C 6 O 5)),
            ("hex", molecular_formula!(H 10 C 6 O 5)),
            ("hexs", molecular_formula!(H 10 C 6 O 8 S 1)),
        ];
        for (name, formula) in cases {
            assert_eq!(
                glycan_parse_list()
                    .iter()
                    .find(|p| p.0 == *name)
                    .unwrap_or_else(|| panic!("Assumed {name} would be defined"))
                    .1
                    .formula(),
                *formula,
                "{name}",
            );
        }
    }

    #[test]
    fn iupac_short_names() {
        let parse = |str: &str| {
            MonoSaccharide::from_short_iupac(str, 0, 0)
                .map(|(res, len)| {
                    assert_eq!(str.len(), len);
                    res
                })
                .unwrap()
        };
        assert_eq!(
            parse("Gal2,3Ac24-1,6-1Py"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                &[
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Pyruvyl,
                ]
            )
        );
        assert_eq!(
            parse("GlcNAc"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                &[GlycanSubstituent::NAcetyl]
            )
        );
        assert_eq!(
            parse("Gal6S"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                &[GlycanSubstituent::Sulfate]
            )
        );
        assert_eq!(
            parse("GlcN2Gc"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                &[GlycanSubstituent::Amino, GlycanSubstituent::Glycolyl,]
            )
        );
        assert_eq!(
            parse("GalNAc3S"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                &[GlycanSubstituent::NAcetyl, GlycanSubstituent::Sulfate]
            )
        );
        assert_eq!(
            parse("GlcN2,6S2"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                &[
                    GlycanSubstituent::Amino,
                    GlycanSubstituent::Sulfate,
                    GlycanSubstituent::Sulfate
                ]
            )
        );
        assert_eq!(
            parse("Tagf1,6P2"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Tagatose)),
                &[GlycanSubstituent::Phosphate, GlycanSubstituent::Phosphate]
            )
            .furanose()
        );
        assert_eq!(
            parse("Gal2,3Ac24-1,6-1Py"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                &[
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Pyruvyl,
                ]
            )
        );
        assert_eq!(
            parse("D-Araf"),
            MonoSaccharide::new(BaseSugar::Pentose(Some(PentoseIsomer::Arabinose)), &[]).furanose()
        );
        assert_eq!(
            parse("Xyl-onic"),
            MonoSaccharide::new(
                BaseSugar::Pentose(Some(PentoseIsomer::Xylose)),
                &[GlycanSubstituent::Acid]
            )
        );
        assert_eq!(
            parse("Glc2,3,4,6Ac4"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                &[
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acetyl
                ]
            )
        );
    }

    #[allow(clippy::float_cmp)] // Handled in a different way
    fn iupac_masses() {
        assert_eq!(
            MonoSaccharide::from_short_iupac("Gal3DiMe(b1-4)GlcNAc(b1-", 0, 0)
                .unwrap()
                .0
                .formula()
                .monoisotopic_mass()
                .value
                .round(),
            411.0
        );
    }
}
