use std::num::NonZeroU16;

use crate::{system::Mass, Element, MolecularFormula, Multi, Tolerance};

/// Find the isobaric sets for the given mass with the given modifications and ppm error.
/// The modifications are placed on any location they are allowed based on the given placement
/// rules, so using any modifications which provide those is advised. If the provided [`LinearPeptide`]
/// has multiple formulas, it uses the formula with the lowest monoisotopic mass.
/// # Panics
/// Panics if any of the modifications does not have a defined mass. Or if the weight of the
/// base selection is already in the tolerance of the given mass.
pub fn find_formulas(
    mass: Mass,
    tolerance: Tolerance<Mass>,
    elements: &[(Element, Option<NonZeroU16>)],
) -> Multi<MolecularFormula> {
    let bounds = tolerance.bounds(mass);
    let mut in_bounds = Vec::new();
    let mut options: Vec<(Mass, MolecularFormula)> = Vec::new();

    for (element, isotope) in elements.iter().copied() {
        if !element.is_valid(isotope) {
            continue;
        }
        let mass = element.mass(isotope).unwrap();
        let mut new_options = Vec::with_capacity(options.len());
        if mass <= bounds.1 {
            new_options.extend((1..=(bounds.1 / mass).value.floor() as i32).map(|n| {
                (
                    mass * f64::from(n),
                    MolecularFormula::new(&[(element, isotope, n)], &[]).unwrap(),
                )
            }));
        }
        for option in &options {
            let rest = bounds.1 - option.0;
            if mass <= rest {
                new_options.extend((1..=(rest / mass).value.floor() as i32).map(|n| {
                    let mut new_formula = option.1.clone();
                    let _ = new_formula.add((element, isotope, n));
                    (option.0 + mass * f64::from(n), new_formula)
                }));
            }
        }
        options.extend_from_slice(&new_options);
    }
    for (mass, option) in options {
        if mass >= bounds.0 && mass <= bounds.1 {
            in_bounds.push(option);
        }
    }

    in_bounds.into()
}
