use std::collections::HashMap;

use crate::{
    fragment::Fragment,
    model::Model,
    peptide::Peptide,
    system::{f64::*, mass_over_charge::mz},
};

#[derive(Clone, Debug)]
pub struct RawSpectrum {
    pub title: String,
    pub num_scans: u64,
    pub rt: Time,
    pub charge: Charge,
    pub mass: Mass,
    pub intensity: Option<f64>,
    pub spectrum: Vec<RawPeak>,
}

impl RawSpectrum {
    pub fn annotate(
        &self,
        peptide: Peptide,
        theoretical_fragments: &[Fragment],
        model: &Model,
    ) -> AnnotatedSpectrum {
        let mut annotated = AnnotatedSpectrum {
            title: self.title.clone(),
            num_scans: self.num_scans,
            rt: self.rt,
            charge: self.charge,
            mass: self.mass,
            peptide,
            spectrum: Vec::with_capacity(self.spectrum.len()),
        };

        let mut connections = Vec::with_capacity(self.spectrum.len());

        for (fragment_index, fragment) in theoretical_fragments.iter().enumerate() {
            connections.extend(self.spectrum.iter().enumerate().filter_map(|(i, p)| {
                let ppm = p.ppm(fragment);
                if ppm < model.ppm {
                    Some((
                        i,
                        fragment_index,
                        AnnotatedPeak::new(p, fragment.clone()),
                        ppm,
                    ))
                } else {
                    None
                }
            }));
        }
        annotated.spectrum.extend(cluster_matches(
            connections,
            &self.spectrum,
            annotated.peptide.sequence.len(),
            model,
        ));

        annotated
    }
}

type Connection = (usize, usize, AnnotatedPeak, MassOverCharge);

fn cluster_matches(
    matches: Vec<Connection>,
    spectrum: &[RawPeak],
    peptide_length: usize,
    model: &Model,
) -> Vec<AnnotatedPeak> {
    let mut found_peak_indices = HashMap::new();
    let mut found_fragment_indices = HashMap::new();
    for pair in &matches {
        *found_peak_indices.entry(pair.0).or_insert(0) += 1;
        *found_fragment_indices.entry(pair.1).or_insert(0) += 1;
    }
    let mut output = Vec::with_capacity(20 * peptide_length + 75); // Empirically derived required size of the buffer (Derived from Hecklib)
    let mut selected_peaks = Vec::new();
    let mut ambiguous = Vec::new();
    // First get all peaks that are unambiguously matched out of the selection to prevent a lot of computation
    for pair in matches {
        if found_peak_indices.get(&pair.0).map_or(false, |v| *v == 1)
            && found_fragment_indices
                .get(&pair.1)
                .map_or(false, |v| *v == 1)
        {
            output.push(pair.2);
            selected_peaks.push(pair.0);
        } else {
            ambiguous.push(pair);
        }
    }
    dbg!(&output, &ambiguous, output.len(), ambiguous.len());

    // Now find all possible combinations of the ambiguous matches and get the non expandable set with the lowest total ppm error
    let sets = combinations(&ambiguous, &[], 0);
    let max_number_connections =
        (found_peak_indices.len() - output.len()).min(found_fragment_indices.len() - output.len());
    let selected_set = sets.into_iter().fold(
        (
            MassOverCharge::new::<mz>(f64::INFINITY),
            Vec::new(),
            Vec::new(),
        ),
        |acc, item| {
            if acc.0 > item.0 + (max_number_connections - acc.1.len()) as f64 * model.ppm {
                item
            } else {
                acc
            }
        },
    );
    dbg!(&selected_set);
    output.extend(selected_set.1);
    selected_peaks.extend(selected_set.2);
    selected_peaks.sort_unstable();
    output.extend(spectrum.iter().enumerate().filter_map(|(i, p)| {
        if selected_peaks.binary_search(&i).is_err() {
            Some(AnnotatedPeak::background(p))
        } else {
            None
        }
    }));
    output
}

/// Recursively get all possible sets for the connection of a single extra time point
fn combinations(
    connections: &[Connection],
    selected: &[usize],
    skip: usize,
) -> Vec<(MassOverCharge, Vec<AnnotatedPeak>, Vec<usize>)> {
    let mut output = Vec::new();
    let mut found = false;

    let selected_indices0: Vec<usize> = selected.iter().map(|c| connections[*c].0).collect();
    let selected_indices1: Vec<usize> = selected.iter().map(|c| connections[*c].1).collect();
    for (extra_index, connection) in connections.iter().skip(skip).enumerate() {
        if !selected_indices0.contains(&connection.0) && !selected_indices1.contains(&connection.1)
        {
            found = true;
            let mut sel: Vec<_> = selected.into();
            sel.push(skip + extra_index);
            output.append(&mut combinations(connections, &sel, skip + extra_index + 1));
        }
    }
    if !found {
        output.push((
            selected
                .iter()
                .map(|c| connections[*c].3)
                .sum::<MassOverCharge>(),
            selected.iter().map(|c| connections[*c].2.clone()).collect(),
            selected.iter().map(|c| connections[*c].0).collect(),
        ));
    }
    output
}

#[derive(Clone, Debug)]
pub struct AnnotatedSpectrum {
    pub title: String,
    pub num_scans: u64,
    pub rt: Time,
    pub charge: Charge,
    pub mass: Mass,
    pub peptide: Peptide,
    pub spectrum: Vec<AnnotatedPeak>,
}

#[derive(Clone, Debug)]
pub struct RawPeak {
    pub mz: MassOverCharge,
    pub intensity: f64,
    pub charge: Charge,
}

impl RawPeak {
    pub fn ppm(&self, fragment: &Fragment) -> MassOverCharge {
        (self.mz - fragment.mz()).abs() / fragment.mz() * MassOverCharge::new::<mz>(1e6)
    }
}

#[derive(Clone, Debug)]
pub struct AnnotatedPeak {
    pub experimental_mz: MassOverCharge,
    pub intensity: f64,
    pub charge: Charge,
    pub annotation: Option<Fragment>,
}

impl AnnotatedPeak {
    pub fn new(peak: &RawPeak, annotation: Fragment) -> Self {
        Self {
            experimental_mz: peak.mz,
            intensity: peak.intensity,
            charge: peak.charge,
            annotation: Some(annotation),
        }
    }

    pub fn background(peak: &RawPeak) -> Self {
        Self {
            experimental_mz: peak.mz,
            intensity: peak.intensity,
            charge: peak.charge,
            annotation: None,
        }
    }
}
