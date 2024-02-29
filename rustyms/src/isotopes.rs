use crate::{system::da, system::Mass, MolecularFormula};
use itertools::Itertools;
use ndarray::{arr1, concatenate, s, Array1, Axis};
use probability::distribution::{Binomial, Discrete};
use std::{cmp::Ordering, collections::HashMap};

impl MolecularFormula {
    /// Get the isotopic distribution, using the natural distribution as defined by CIAAW.
    /// All elements are considered. The return is an array with the probability per offset.
    /// The first element of the array is the base peak, every consecutive peak is 1 dalton heavier.
    /// The probability os normalized to (approximately) 1 total area.
    fn isotopic_distribution(&self, threshold: f64) -> Array1<f64> {
        let mut result = arr1(&[1.0]);
        for (element, isotope, amount) in self.elements() {
            if isotope.is_some() || *amount <= 0 {
                // TODO: think about negative numbers?
                continue;
            }
            let isotopes = element
                .isotopes()
                .iter()
                .filter(|i| i.2 != 0.0)
                .collect_vec();
            if isotopes.len() < 2 {
                // Only a single species, so no distribution is needed
                continue;
            }
            // Get the probability and base offset (weight) for all non base isotopes
            let base = isotopes[0];
            let isotopes = isotopes
                .into_iter()
                .skip(1)
                .map(|i| (i.0 - base.0, i.2))
                .collect_vec();

            for isotope in isotopes {
                // Generate distribution (take already chosen into account?)
                let binomial = Binomial::new(usize::try_from(*amount).unwrap(), isotope.1);
                let mut last: Option<f64> = None;
                let mut distribution: Array1<f64> = (0..=usize::try_from(*amount).unwrap())
                    .map(|t| binomial.mass(t))
                    .take_while(|a| {
                        // Take all numbers until the threshold is crossed at the tail end of the distribution
                        if let Some(last) = &last {
                            *last < threshold || *a > threshold
                        } else {
                            last = Some(*a);
                            true
                        }
                    })
                    .flat_map(|a| {
                        // Interweave the probability of this isotope with the mass difference to generate the correct distribution
                        std::iter::once(a)
                            .chain(std::iter::repeat(0.0))
                            .take(isotope.0 as usize)
                    })
                    .collect();

                // Make the lengths equal
                match result.len().cmp(&distribution.len()) {
                    Ordering::Less => {
                        let mut zeros = Array1::zeros(distribution.len());
                        for (z, r) in zeros.iter_mut().zip(result) {
                            *z = r;
                        }
                        result = zeros;
                    }
                    Ordering::Greater => {
                        let mut zeros = Array1::zeros(result.len());
                        for (z, d) in zeros.iter_mut().zip(distribution) {
                            *z = d;
                        }
                        distribution = zeros;
                    }
                    Ordering::Equal => (),
                }

                // Combine distribution with previous distribution
                let mut new = Array1::zeros(result.len());
                let shift = 0; // (isotope.0 - 1) as usize;
                for (i, a) in distribution.into_iter().enumerate() {
                    new += &(concatenate(
                        Axis(0),
                        &[
                            Array1::zeros((i + shift).min(result.len())).view(),
                            result.slice(s![0..result.len().saturating_sub(i + shift)]),
                        ],
                    )
                    .unwrap()
                        * a);
                }

                result = new;
            }
        }
        result
    }
}

fn combined_pattern(
    isotopes: &[(MolecularFormula, f64, String)],
) -> Vec<(Mass, f64, usize, String)> {
    let mut combined: Vec<(Mass, f64, usize, String)> = Vec::new();

    for isotope in isotopes {
        let isotope_mass = isotope.0.monoisotopic_mass().value.round();
        if let Some(entry) = combined
            .iter_mut()
            .find(|i| (i.0.value - isotope_mass).abs() < f64::EPSILON)
        {
            entry.1 += isotope.1;
            entry.2 += 1;
            entry.3 += &format!(",{}", isotope.2);
        } else {
            combined.push((da(isotope_mass), isotope.1, 1, isotope.2.clone()));
        }
    }
    combined.sort_unstable_by(|a, b| b.1.total_cmp(&a.1));
    combined
}

fn binom(tries: i16, total: i16, p: f64) -> f64 {
    Binomial::new(total as usize, p).mass(tries as usize)
}

fn memoized_f64_factorial(num: u16, cache: &mut HashMap<u16, f64>) -> f64 {
    *cache
        .entry(num)
        .or_insert_with(|| stupid_f64_factorial(num))
}

fn stupid_f64_factorial(num: u16) -> f64 {
    (2..=num).fold(1.0, |acc, i| acc * f64::from(i))
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use super::*;

    #[test]
    fn distribution() {
        let formula = molecular_formula!(H 10 C 10 O 10).unwrap();
        dbg!(
            formula.monoisotopic_mass(),
            formula.isotopic_distribution(0.0001)
        );
        panic!();
    }

    //#[test]
    // fn simple_isotope_pattern() {
    //     // EVQLVESGGGLVQPGG start 16 AA of herceptin
    //     let peptide = MolecularFormula::new(&[
    //         (Element::H, None, 108),
    //         (Element::C, None, 65),
    //         (Element::N, None, 18),
    //         (Element::O, None, 24),
    //     ])
    //     .unwrap();
    //     let isotopes = peptide.isotopic_distribution(1e-6, true, true, true, true);
    //     save_combinations(
    //         &combined_pattern(&isotopes),
    //         "target/peptide_combined_all.tsv",
    //     );
    //     save_combinations(
    //         &combined_pattern(&peptide.isotopic_distribution(1e-5, true, false, false, false)),
    //         "target/peptide_combined_only_H.tsv",
    //     );
    //     save_combinations(
    //         &combined_pattern(&peptide.isotopic_distribution(1e-5, false, true, false, false)),
    //         "target/peptide_combined_only_C.tsv",
    //     );
    //     save_combinations(
    //         &combined_pattern(&peptide.isotopic_distribution(1e-5, false, false, true, false)),
    //         "target/peptide_combined_only_N.tsv",
    //     );
    //     save_combinations(
    //         &combined_pattern(&peptide.isotopic_distribution(1e-5, false, false, false, true)),
    //         "target/peptide_combined_only_O.tsv",
    //     );
    //     println!(
    //         "MonoIsotopic mass: {} average mass: {}",
    //         peptide.monoisotopic_mass().value,
    //         peptide.average_weight().value
    //     );
    // }

    // //#[test]
    // fn sars_cov_2_spike_isotope_pattern() {
    //     // MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
    //     let spike = MolecularFormula::new(&[
    //         (Element::H, None, 9770),
    //         (Element::C, None, 6336),
    //         (Element::N, None, 1656),
    //         (Element::O, None, 1894),
    //         (Element::S, None, 54),
    //     ])
    //     .unwrap();
    //     let isotopes = spike.isotopic_distribution(1e-5, true, true, true, true);
    //     let file_handler = File::create("target/spike_isotopes.tsv").unwrap();
    //     let mut writer = BufWriter::new(file_handler);
    //     writeln!(writer, "Name\tMass\tProbability").unwrap();
    //     for isotope in &isotopes {
    //         writeln!(writer, "{}\t{}\t{}", isotope.2, isotope.0, isotope.1).unwrap();
    //     }
    //     let combined = combined_pattern(&isotopes);
    //     save_combinations(&combined, "target/spike_combined_all.tsv");
    //     save_combinations(
    //         &combined_pattern(&spike.isotopic_distribution(1e-5, true, false, false, false)),
    //         "target/spike_combined_only_H.tsv",
    //     );
    //     save_combinations(
    //         &combined_pattern(&spike.isotopic_distribution(1e-5, false, true, false, false)),
    //         "target/spike_combined_only_C.tsv",
    //     );
    //     save_combinations(
    //         &combined_pattern(&spike.isotopic_distribution(1e-5, false, false, true, false)),
    //         "target/spike_combined_only_N.tsv",
    //     );
    //     save_combinations(
    //         &combined_pattern(&spike.isotopic_distribution(1e-5, false, false, false, true)),
    //         "target/spike_combined_only_O.tsv",
    //     );
    //     println!(
    //         "MonoIsotopic mass: {} average mass: {}",
    //         spike.monoisotopic_mass().value,
    //         spike.average_weight().value
    //     );
    // }

    // fn save_combinations(data: &[(Mass, f64, usize, String)], name: &str) {
    //     let file_handler = File::create(name).unwrap();
    //     let mut writer = BufWriter::new(file_handler);
    //     writeln!(writer, "Name\tMass\tProbability\tTotalIsotopes").unwrap();
    //     for isotope in data {
    //         writeln!(
    //             writer,
    //             "{}\t{}\t{}\t{}",
    //             isotope.3, isotope.0.value, isotope.1, isotope.2
    //         )
    //         .unwrap();
    //     }
    // }

    // //#[test]
    // fn herceptin_v_heavy_isotope_pattern() {
    //     // EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTK
    //     let herceptin = MolecularFormula::new(&[
    //         (Element::H, None, 915),
    //         (Element::C, None, 602),
    //         (Element::N, None, 165),
    //         (Element::O, None, 185),
    //         (Element::S, None, 4),
    //     ])
    //     .unwrap();
    //     let isotopes = herceptin.isotopic_distribution(1e-5, true, true, true, true);
    //     let file_handler = File::create("target/herceptin_isotopes.tsv").unwrap();
    //     let mut writer = BufWriter::new(file_handler);
    //     writeln!(writer, "Name\tMass\tProbability").unwrap();
    //     for isotope in &isotopes {
    //         writeln!(writer, "{}\t{}\t{}", isotope.2, isotope.0, isotope.1).unwrap();
    //     }
    //     let combined = combined_pattern(&isotopes);
    //     save_combinations(&combined, "target/herceptin_combined_all.tsv");
    //     save_combinations(
    //         &combined_pattern(&herceptin.isotopic_distribution(1e-5, true, false, false, false)),
    //         "target/herceptin_combined_only_H.tsv",
    //     );
    //     save_combinations(
    //         &combined_pattern(&herceptin.isotopic_distribution(1e-5, false, true, false, false)),
    //         "target/herceptin_combined_only_C.tsv",
    //     );
    //     save_combinations(
    //         &combined_pattern(&herceptin.isotopic_distribution(1e-5, false, false, true, false)),
    //         "target/herceptin_combined_only_N.tsv",
    //     );
    //     save_combinations(
    //         &combined_pattern(&herceptin.isotopic_distribution(1e-5, false, false, false, true)),
    //         "target/herceptin_combined_only_O.tsv",
    //     );
    //     println!(
    //         "MonoIsotopic mass: {} average mass: {}",
    //         herceptin.monoisotopic_mass().value,
    //         herceptin.average_weight().value
    //     );
    // }

    // #[test]
    // fn isotope_masses() {
    //     use itertools::Itertools;
    //     fn show(element: Element) {
    //         let isotopes = element
    //             .isotopes()
    //             .iter()
    //             .filter(|i| i.2 != 0.0)
    //             .collect_vec();
    //         if isotopes.len() < 2 {
    //             return;
    //         }
    //         print!("{element}\t");
    //         let first = isotopes[0].1.value;
    //         let base = isotopes[0].0;
    //         for (isotope, mass, fraction) in isotopes {
    //             print!(
    //                 "{}\t{}\t{}\t{}\t",
    //                 isotope,
    //                 mass.value,
    //                 (1.0 - (mass.value - first) / (isotope - base) as f64).abs(),
    //                 fraction,
    //             );
    //         }
    //         println!();
    //     }
    //     for element in (Element::H as usize)..=(Element::Og as usize) {
    //         show(Element::try_from(element).unwrap());
    //     }
    //     panic!()
    // }
}
