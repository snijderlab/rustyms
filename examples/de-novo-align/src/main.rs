use std::{
    fs::File,
    io::{BufReader, BufWriter},
};

use clap::Parser;
use fragment::FragmentType;
use identification::{open_identified_peptides_file, FastaData};
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
    let peptides = open_identified_peptides_file(args.peptides, None).unwrap();
    let database = FastaData::parse_file(args.database).unwrap();
}
