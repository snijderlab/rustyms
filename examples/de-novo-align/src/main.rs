use std::{collections::HashMap, fs::File, io::BufWriter};

use clap::Parser;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use rustyms::{
    align::{align, matrix, AlignType},
    csv::write_csv,
    identification::{open_identified_peptides_file, FastaData},
    system::da,
    *,
};

#[derive(Parser)]
struct Cli {
    /// The input identified peptides file
    #[arg(short, long)]
    peptides: String,
    /// The fasta database of known proteins
    #[arg(short, long)]
    database: String,
    /// Where to store the results
    #[arg(long)]
    out_path: String,
}

fn main() {
    let args = Cli::parse();
    let out_file = BufWriter::new(File::create(args.out_path).unwrap());
    let peptides = open_identified_peptides_file(args.peptides, None)
        .unwrap()
        .filter_map(|p| p.ok())
        .collect_vec();
    let database = FastaData::parse_file(args.database).unwrap();

    let alignments: Vec<_> = peptides
        .par_iter()
        .map(|peptide| {
            database
                .iter()
                .map(|db| {
                    (
                        db,
                        peptide,
                        align::<4, SemiAmbiguous, SemiAmbiguous>(
                            &db.peptide,
                            peptide.metadata.peptide().unwrap(),
                            matrix::BLOSUM62,
                            Tolerance::Absolute(da(0.1)),
                            AlignType::EITHER_GLOBAL,
                        ),
                    )
                })
                .max_by_key(|a| OrderedFloat(a.2.normalised_score()))
                .unwrap()
        })
        .map(|(db, peptide, alignment)| {
            HashMap::from([
                ("Peptide".to_string(), alignment.seq_b().to_string()),
                (
                    "Rawfile".to_string(),
                    peptide
                        .metadata
                        .raw_file()
                        .map_or(String::new(), |p| p.to_string_lossy().to_string()),
                ),
                (
                    "De novo score".to_string(),
                    peptide.score.map_or(String::new(), |s| s.to_string()),
                ),
                ("Protein".to_string(), db.id.clone()),
                (
                    "Score".to_string(),
                    alignment.normalised_score().to_string(),
                ),
                ("Start".to_string(), alignment.start_a().to_string()),
                (
                    "End".to_string(),
                    (alignment.start_a() + alignment.len_a()).to_string(),
                ),
                ("Path".to_string(), alignment.short()),
            ])
        })
        .collect();

    write_csv(out_file, alignments).unwrap();
}
