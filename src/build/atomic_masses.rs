use std::{ffi::OsString, fs, path::Path};

use crate::Element;

use super::csv::parse_csv;

pub fn build_atomic_masses(out_dir: &OsString, _debug: bool) -> Result<(), String> {
    let mut atomic_weights = vec![None; 118];
    let mut isotopic_abundances = vec![Vec::new(); 118];
    let mut atomic_masses = vec![Vec::new(); 118];
    // TODO: write output, read isotopic abundances, average weights, and determine normal monoisotopic mass

    let table = parse_csv("databases/IUPAC-atomic-masses.csv")?;
    for line in table {
        let line = line?;
        let (nuclide, mass, _uncertainty, year) = (&line[0], &line[1], &line[2], &line[3]);
        if nuclide.starts_with("AME")
            || nuclide.is_empty()
            || nuclide == "nuclide"
            || !year.ends_with("2020</a>")
        {
            continue;
        }
        let isotope = nuclide
            .trim_end_matches(|c: char| c.is_alphabetic())
            .parse::<usize>()
            .map_err(|e| e.to_string())?;
        let element = Element::try_from(nuclide.trim_start_matches(|c: char| c.is_ascii_digit()))
            .map_err(|_| {
            format!("Not a valid isotope+element, could not recognise element: {nuclide}")
        })?;
        let mass = mass.parse::<f64>().map_err(|e| format!("{}@{}", e, mass))?;
        atomic_masses[element as usize - 1].push((isotope, mass))
    }

    let mut last_element = 0;
    let table = parse_csv("databases/CIAAW-isotopic-abundances.csv")?;
    for line in table {
        let line = line?;
        let (element, _element, _name, isotope, abundance, _note) =
            (&line[0], &line[1], &line[2], &line[3], &line[4], &line[5]);
        let mut abundance = abundance.to_owned();
        abundance.retain(|c| !c.is_whitespace()); // Remove any whitespace, including any sneaking non breaking spaces
        if element == "Z" || isotope.starts_with('[') || abundance == "-" {
            continue;
        }

        let element = if element.is_empty() {
            last_element
        } else {
            last_element = element
                .parse::<usize>()
                .map_err(|_| format!("Not a valid number for element Z: {element}"))?;
            last_element
        };

        let isotope = isotope.parse::<usize>().unwrap();

        isotopic_abundances[element - 1].push((isotope, get_ciaaw_number(&abundance)?))
    }

    let table = parse_csv("databases/CIAAW-atomic-weights.csv")?;
    for line in table {
        let line = line?;
        let (element, weight) = (&line[0], &line[3]);
        let mut weight = weight.to_owned();
        weight.retain(|c| !c.is_whitespace()); // Remove any whitespace, including any sneaking non breaking spaces
        if element == "Z" || weight == "—" {
            continue;
        }
        let element = element
            .parse::<usize>()
            .map_err(|_| format!("Not valid element number (not a number): {element}"))?;
        atomic_weights[element - 1] = Some(get_ciaaw_number(&weight)?);
    }

    // Monoisotopic, average weight, isotopes: (num, mass, abundance)
    let combined_data = isotopic_abundances
        .into_iter()
        .zip(atomic_masses)
        .zip(atomic_weights)
        .map(|((isotopic_abundance, isotopic_mass), atomic_weight)| {
            let isotopes = isotopic_mass
                .iter()
                .map(|(i, m)| {
                    (
                        *i,
                        *m,
                        isotopic_abundance
                            .iter()
                            .find(|(ai, _)| ai == i)
                            .map(|(_, a)| *a)
                            .unwrap_or(0.0),
                    )
                })
                .collect::<Vec<_>>();
            (
                isotopes
                    .iter()
                    .fold(
                        (None, 0.0),
                        |acc, (_, m, a)| if *a > acc.1 { (Some(*m), *a) } else { acc },
                    )
                    .0,
                atomic_weight,
                isotopes,
            )
        });

    let dest_path = Path::new(&out_dir).join("elements.rs");
    fs::write(
        dest_path,
        format!(
            "#[allow(clippy::all, clippy::pedantic, clippy::nursery)]\n/// The data for all the elements (atomic_mass/monoisotopic mass, atomic_weight, isotopes: (num, mass, abundance))\npub const ELEMENTAL_DATA: &[(Option<f64>, Option<f64>, &[(u16, f64, f64)])] = &[
                {}
                ];",
            combined_data.fold(String::new(), |acc, element| format!(
                "{}\n({:?}, {:?}, &[{}]),",
                acc,
                element.0,
                element.1,
                element.2.iter().fold(String::new(), |acc, i| format!(
                    "{}({},{:.1},{:.1}),",
                    acc, i.0, i.1, i.2
                ))
            ))
        ),
    )
    .unwrap();

    Ok(())
}

fn get_ciaaw_number(text: &str) -> Result<f64, String> {
    let parse = |t: &str| {
        t.parse::<f64>()
            .map_err(|_| format!("Not a valid number: {t}"))
    };
    if text.starts_with('[') {
        let (low, high) = text[1..text.len() - 1]
            .split_once(',')
            .ok_or(format!("Not a valid range: {text}"))?;
        Ok((parse(low)? + parse(high)?) / 2.0)
    } else if text.ends_with(')') {
        Ok(parse(
            text.split_once('(')
                .ok_or(format!("Not valid error indication: {text}"))?
                .0,
        )?)
    } else {
        Ok(parse(text)?)
    }
}
