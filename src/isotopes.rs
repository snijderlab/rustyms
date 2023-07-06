use crate::{dalton, element::Element, r, Mass, MolecularFormula, Ratio};
use statrs::distribution::{Binomial, Discrete};
use std::collections::HashMap;

impl MolecularFormula {
    #[allow(clippy::fn_params_excessive_bools)]
    fn isotopic_distribution(
        &self,
        threshold: f64,
        use_h: bool,
        use_c: bool,
        use_n: bool,
        use_o: bool,
    ) -> Option<Vec<(MolecularFormula, f64, String)>> {
        let h = Element::H.isotopes();
        let c = Element::C.isotopes();
        let n = Element::N.isotopes();
        let o = Element::O.isotopes();
        let additional_mass: Mass = self
            .elements()
            .iter()
            .filter(|i| {
                i.1 == 0
                    && (i.0 != Element::H
                        && i.0 != Element::C
                        && i.0 != Element::N
                        && i.0 != Element::O)
            })
            .fold(Some(MolecularFormula::default()), |acc, s| {
                s.0.mass(s.1)
                    .zip(acc)
                    .map(|(m, a)| a + m * Ratio::new::<r>(f64::from(s.2)))
            })?;
        let mut isotopes = Vec::new();
        let present_h = self
            .elements()
            .iter()
            .find(|i| i.0 == Element::H && i.1 == 0)
            .map_or(0, |i| i.2);
        let present_c = self
            .elements()
            .iter()
            .find(|i| i.0 == Element::C && i.1 == 0)
            .map_or(0, |i| i.2);
        let present_n = self
            .elements()
            .iter()
            .find(|i| i.0 == Element::N && i.1 == 0)
            .map_or(0, |i| i.2);
        let present_o = self
            .elements()
            .iter()
            .find(|i| i.0 == Element::O && i.1 == 0)
            .map_or(0, |i| i.2);
        let max_h = if use_h { present_h } else { 0 };
        let max_c = if use_c { present_c } else { 0 };
        let max_n = if use_n { present_n } else { 0 };
        let max_o = if use_o { present_o } else { 0 };

        for num_h in 0..=max_h {
            let h_chance = binom(num_h, max_h, h[1].2);
            if h_chance < threshold {
                continue;
            }
            for num_c in 0..=max_c {
                let c_chance = binom(num_c, max_c, c[1].2);
                if h_chance * c_chance < threshold {
                    continue;
                }
                for num_n in 0..=max_n {
                    let n_chance = binom(num_n, max_n, n[1].2);
                    if h_chance * c_chance * n_chance < threshold {
                        continue;
                    }
                    for num_o_1 in 0..=max_o {
                        let o_1_chance = binom(num_o_1, max_o, o[1].2);
                        if h_chance * c_chance * n_chance * o_1_chance < threshold {
                            continue;
                        }
                        for num_o_2 in 0..=max_o - num_o_1 {
                            let o_2_chance = binom(num_o_2, max_o, o[2].2);
                            let chance = h_chance * c_chance * n_chance * o_1_chance * o_2_chance;
                            if chance < threshold {
                                continue;
                            }
                            // Do the calc
                            let mass = Mass::new::<dalton>(
                                f64::from(present_h - num_h) * h[0].1
                                    + f64::from(num_h) * h[1].1
                                    + f64::from(present_c - num_c) * c[0].1
                                    + f64::from(num_c) * c[1].1
                                    + f64::from(present_n - num_n) * n[0].1
                                    + f64::from(num_n) * n[1].1
                                    + f64::from(present_o - num_o_1 - num_o_2) * o[0].1
                                    + f64::from(num_o_1) * o[1].1
                                    + f64::from(num_o_2) * o[2].1,
                            );
                            isotopes.push((
                                mass + additional_mass,
                                chance,
                                format!("[2H{num_h}][13C{num_c}][15N{num_n}][17O{num_o_1}][18O{num_o_2}]"),
                            ));
                        }
                    }
                }
            }
        }
        isotopes.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        Some(isotopes)
    }
}

fn combined_pattern(
    isotopes: &[(MolecularFormula, f64, String)],
) -> Vec<(MolecularFormula, f64, usize, String)> {
    let mut combined: Vec<(MolecularFormula, f64, usize, String)> = Vec::new();

    for isotope in isotopes {
        if let Some(entry) = combined
            .iter_mut()
            .find(|i| (i.0.value - isotope.0.value.round()).abs() < f64::EPSILON)
        {
            entry.1 += isotope.1;
            entry.2 += 1;
            entry.3 += &format!(",{}", isotope.2);
        } else {
            combined.push((
                Mass::new::<dalton>(isotope.0.value.round()),
                isotope.1,
                1,
                isotope.2.clone(),
            ));
        }
    }
    combined.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    combined
}

fn binom(tries: u16, total: u16, p: f64) -> f64 {
    Binomial::new(p, u64::from(total))
        .unwrap()
        .pmf(u64::from(tries))
    //memoized_f64_factorial(total, cache)
    //    / (memoized_f64_factorial(tries, cache) * memoized_f64_factorial(total - tries, cache))
    //    * p.powi(tries as i32)
    //    * (1.0 - p).powi((total - tries) as i32)
    //poisson f.powi(i32::from(n)) * f64::exp(-f) / stupid_f64_factorial(n)
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
mod tests {
    use super::*;
    use std::{
        fs::File,
        io::{BufWriter, Write},
    };

    //#[test]
    fn simple_isotope_pattern() -> Option<()> {
        // EVQLVESGGGLVQPGG start 16 AA of herceptin
        let peptide = MolecularFormula::new(&[
            (Element::H, 0, 108),
            (Element::C, 0, 65),
            (Element::N, 0, 18),
            (Element::O, 0, 24),
        ]);
        let isotopes = peptide.isotopic_distribution(1e-6, true, true, true, true);
        save_combinations(
            &combined_pattern(&isotopes?),
            "target/peptide_combined_all.tsv",
        );
        save_combinations(
            &combined_pattern(&peptide.isotopic_distribution(1e-5, true, false, false, false)?),
            "target/peptide_combined_only_H.tsv",
        );
        save_combinations(
            &combined_pattern(&peptide.isotopic_distribution(1e-5, false, true, false, false)?),
            "target/peptide_combined_only_C.tsv",
        );
        save_combinations(
            &combined_pattern(&peptide.isotopic_distribution(1e-5, false, false, true, false)?),
            "target/peptide_combined_only_N.tsv",
        );
        save_combinations(
            &combined_pattern(&peptide.isotopic_distribution(1e-5, false, false, false, true)?),
            "target/peptide_combined_only_O.tsv",
        );
        println!(
            "MonoIsotopic mass: {} average mass: {}",
            peptide.monoisotopic_mass().unwrap().value,
            peptide.average_weight().unwrap().value
        );
        Some(())
    }

    //#[test]
    fn sars_cov_2_spike_isotope_pattern() -> Option<()> {
        // MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
        let spike = MolecularFormula::new(&[
            (Element::H, 0, 9770),
            (Element::C, 0, 6336),
            (Element::N, 0, 1656),
            (Element::O, 0, 1894),
            (Element::S, 0, 54),
        ]);
        let isotopes = spike.isotopic_distribution(1e-5, true, true, true, true)?;
        let file_handler = File::create("target/spike_isotopes.tsv").unwrap();
        let mut writer = BufWriter::new(file_handler);
        writeln!(writer, "Name\tMass\tProbability").unwrap();
        for isotope in &isotopes {
            writeln!(writer, "{}\t{}\t{}", isotope.2, isotope.0.value, isotope.1).unwrap();
        }
        let combined = combined_pattern(&isotopes);
        save_combinations(&combined, "target/spike_combined_all.tsv");
        save_combinations(
            &combined_pattern(&spike.isotopic_distribution(1e-5, true, false, false, false)?),
            "target/spike_combined_only_H.tsv",
        );
        save_combinations(
            &combined_pattern(&spike.isotopic_distribution(1e-5, false, true, false, false)?),
            "target/spike_combined_only_C.tsv",
        );
        save_combinations(
            &combined_pattern(&spike.isotopic_distribution(1e-5, false, false, true, false)?),
            "target/spike_combined_only_N.tsv",
        );
        save_combinations(
            &combined_pattern(&spike.isotopic_distribution(1e-5, false, false, false, true)?),
            "target/spike_combined_only_O.tsv",
        );
        println!(
            "MonoIsotopic mass: {} average mass: {}",
            spike.monoisotopic_mass().unwrap().value,
            spike.average_weight().unwrap().value
        );
        Some(())
    }

    fn save_combinations(data: &[(Mass, f64, usize, String)], name: &str) {
        let file_handler = File::create(name).unwrap();
        let mut writer = BufWriter::new(file_handler);
        writeln!(writer, "Name\tMass\tProbability\tTotalIsotopes").unwrap();
        for isotope in data {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                isotope.3, isotope.0.value, isotope.1, isotope.2
            )
            .unwrap();
        }
    }

    //#[test]
    fn herceptin_v_heavy_isotope_pattern() -> Option<()> {
        // EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTK
        let herceptin = MolecularFormula::new(&[
            (Element::H, 0, 915),
            (Element::C, 0, 602),
            (Element::N, 0, 165),
            (Element::O, 0, 185),
            (Element::S, 0, 4),
        ]);
        let isotopes = herceptin.isotopic_distribution(1e-5, true, true, true, true)?;
        let file_handler = File::create("target/herceptin_isotopes.tsv").unwrap();
        let mut writer = BufWriter::new(file_handler);
        writeln!(writer, "Name\tMass\tProbability").unwrap();
        for isotope in &isotopes {
            writeln!(writer, "{}\t{}\t{}", isotope.2, isotope.0.value, isotope.1).unwrap();
        }
        let combined = combined_pattern(&isotopes);
        save_combinations(&combined, "target/herceptin_combined_all.tsv");
        save_combinations(
            &combined_pattern(&herceptin.isotopic_distribution(1e-5, true, false, false, false)?),
            "target/herceptin_combined_only_H.tsv",
        );
        save_combinations(
            &combined_pattern(&herceptin.isotopic_distribution(1e-5, false, true, false, false)?),
            "target/herceptin_combined_only_C.tsv",
        );
        save_combinations(
            &combined_pattern(&herceptin.isotopic_distribution(1e-5, false, false, true, false)?),
            "target/herceptin_combined_only_N.tsv",
        );
        save_combinations(
            &combined_pattern(&herceptin.isotopic_distribution(1e-5, false, false, false, true)?),
            "target/herceptin_combined_only_O.tsv",
        );
        println!(
            "MonoIsotopic mass: {} average mass: {}",
            herceptin.monoisotopic_mass().unwrap().value,
            herceptin.average_weight().unwrap().value
        );
        Some(())
    }
}
