#![allow(non_snake_case)] // charge_independent_Y needs the capital as it means the glycan fragmentation
use std::{
    fs::File,
    io::{BufReader, BufWriter},
};

use clap::Parser;
use directories::ProjectDirs;
use fragment::FragmentType;
use itertools::Itertools;
use rayon::prelude::*;
use rustyms::{
    spectrum::{Score, Scores},
    system::{e, usize::Charge},
    *,
};
use spectrum::PeakSpectrum;
use std::collections::HashMap;

#[derive(Parser)]
struct Cli {
    /// The input csv file, should have the following columns: 'mgf_path', 'scan_number', 'z', 'sequence', and can have 'fragmentation' (etd/ethcd/hcd/cid/all, default is all)
    #[arg(short, long)]
    in_path: String,
    #[arg(short, long)]
    out_path: String,
    #[arg(long)]
    charge_independent_Y: bool,
    /// To turn off loading the custom modifications database from the Annotator (if installed)
    #[arg(long)]
    no_custom_mods: bool,
}

fn main() {
    let model = Model::all();
    let args = Cli::parse();
    let path = ProjectDirs::from("com", "com.snijderlab.annotator", "")
        .unwrap()
        .config_dir()
        .join("../custom_modifications.json");
    let custom_database = if args.no_custom_mods || !path.exists() {
        None
    } else {
        Some(serde_json::from_reader(BufReader::new(File::open(path).unwrap())).unwrap())
    };
    let files = rustyms::csv::parse_csv(args.in_path, b',', None)
        .unwrap()
        .filter_map(|a| a.ok())
        .into_group_map_by(|l| l.index_column("mgf_path").unwrap().0.to_string());
    let out_file = BufWriter::new(File::create(args.out_path).unwrap());
    let mut out_data = Vec::new();

    for (file_name, lines) in files {
        let file = rustyms::rawfile::mgf::open(&file_name).unwrap();

        let rows = lines
            .par_iter()
            .filter_map(|line| {
                let scan_number = line
                    .index_column("scan_number")
                    .unwrap()
                    .0
                    .parse::<usize>()
                    .unwrap();
                let z = line.index_column("z").unwrap().0.parse::<usize>().unwrap();
                let peptide = CompoundPeptidoform::pro_forma(
                    line.index_column("sequence").unwrap().0,
                    custom_database.as_ref(),
                )
                .unwrap();
                let selected_model = match line
                    .index_column("fragmentation")
                    .map(|v| v.0.to_ascii_lowercase())
                    .as_deref()
                {
                    Ok("etd") => Model::etd(),
                    Ok("ethcd") => Model::ethcd(),
                    Ok("hcd") | Ok("cid") => Model::cid_hcd(),
                    Ok("all") => Model::all(),
                    _ => model.clone(),
                };
                if let Some(spectrum) = file.iter().find(|s| s.raw_scan_number == Some(scan_number))
                {
                    let fragments =
                        peptide.generate_theoretical_fragments(Charge::new::<e>(z), &model);
                    let annotated = spectrum.annotate(
                        peptide,
                        &fragments,
                        &selected_model,
                        MassMode::Monoisotopic,
                    );
                    let scores: &Scores = &annotated
                        .scores(&fragments, &selected_model, MassMode::Monoisotopic)
                        .1[0][0];

                    let mut row: HashMap<_, _> = line.into();

                    if args.charge_independent_Y {
                        let unique_Y = fragments
                            .iter()
                            .filter_map(|fragment| {
                                if let FragmentType::Y(pos) = &fragment.ion {
                                    Some(pos)
                                } else {
                                    None
                                }
                            })
                            .unique()
                            .count();
                        let unique_Y_found = annotated
                            .spectrum()
                            .flat_map(|peak| &peak.annotation)
                            .filter_map(|fragment| {
                                if let FragmentType::Y(pos) = &fragment.ion {
                                    Some(pos)
                                } else {
                                    None
                                }
                            })
                            .unique()
                            .count();
                        row.insert(
                            "ion_Y_charge_independent".to_string(),
                            format!(
                                "{}({}/{})",
                                (unique_Y_found as f64 / unique_Y as f64),
                                unique_Y_found,
                                unique_Y
                            ),
                        );
                    }
                    row.insert(
                        "ion_combined".to_string(),
                        match scores.score {
                            Score::Position {
                                theoretical_positions,
                                ..
                            } => format!(
                                "{}({}/{})",
                                theoretical_positions.fraction(),
                                theoretical_positions.found,
                                theoretical_positions.total
                            ),
                            Score::UniqueFormulas {
                                unique_formulas, ..
                            } => format!(
                                "{}({}/{})",
                                unique_formulas.fraction(),
                                unique_formulas.found,
                                unique_formulas.total
                            ),
                        },
                    );

                    for (ion, score) in &scores.ions {
                        let recovered = match score {
                            Score::Position {
                                theoretical_positions,
                                ..
                            } => theoretical_positions,
                            Score::UniqueFormulas {
                                unique_formulas, ..
                            } => unique_formulas,
                        };
                        row.insert(
                            format!("ion_{ion}"),
                            format!(
                                "{}({}/{})",
                                recovered.fraction(),
                                recovered.found,
                                recovered.total
                            ),
                        );
                    }
                    Some(row)
                } else {
                    eprintln!("Could not find scan number {scan_number} for file {file_name}");
                    None
                }
            })
            .collect::<Vec<_>>();
        out_data.extend_from_slice(&rows);
    }

    rustyms::csv::write_csv(out_file, out_data).unwrap();
}
