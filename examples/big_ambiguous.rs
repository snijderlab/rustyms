use std::collections::HashMap;

//use rustyms::model::NeutralLoss;
use rustyms::FragmentType;
use rustyms::Position;
use rustyms::*;

fn main() {
    let mut ambiguous: Vec<_> = include!["../data/data.txt"].into_iter().take(77).collect();
    let mut found_peak_indices = HashMap::new();
    let mut found_fragment_indices = HashMap::new();
    for pair in &ambiguous {
        *found_peak_indices.entry(pair.0).or_insert(0) += 1;
        *found_fragment_indices.entry(pair.1).or_insert(0) += 1;
    }
    let max_number_connections = found_fragment_indices.len().min(found_peak_indices.len());
    ambiguous.sort_unstable_by(|a, b| a.3.partial_cmp(&b.3).unwrap());
    let mut sets = non_recursive_combinations(&ambiguous, MassOverCharge::new::<mz>(20.0));
    dbg!(ambiguous.len(), sets.len(), max_number_connections);
    sets.iter_mut().for_each(|c| {
        c.0 += (max_number_connections - c.1.len()) as f64 * MassOverCharge::new::<mz>(20.0)
    });
    sets.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    dbg!(sets.iter().map(|set| set.0).collect::<Vec<_>>());
    let selected_set = sets.into_iter().fold(
        (MassOverCharge::new::<mz>(f64::INFINITY), Vec::new()),
        |acc, item| {
            if acc.0 > item.0 {
                item
            } else {
                acc
            }
        },
    );
    dbg!(&selected_set);
}
