#![allow(clippy::redundant_pub_crate)]

use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, BufWriter, Write},
};

#[path = "../../rustyms/src/imgt/shared/mod.rs"]
mod shared;

mod combine;
mod imgt_gene;
mod parse;
mod structs;

use crate::shared::*;

use itertools::Itertools;
use rustyms::*;
use structs::{Location, SequenceRegion};

fn main() {
    let file = File::open("data/imgt.dat")
        .expect("Please provide the 'imgt.dat' file in the 'data' directory in the root.");
    let mut output = BufWriter::new(File::create("rustyms/src/imgt/germlines/mod.rs").unwrap());
    let mut docs = BufWriter::new(File::create("rustyms/src/imgt/germlines/germlines.md").unwrap());
    let mut error = BufWriter::new(File::create("errors.dat").unwrap());
    let data = parse::parse_dat(BufReader::new(file));

    let (grouped, errors) = combine::combine(data);

    // Keep track of all errors
    for (species, errors) in errors
        .iter()
        .map(|(species, gene, err)| (species, (gene, err)))
        .into_group_map()
        .into_iter()
        .map(|(species, genes)| {
            (
                species,
                genes
                    .into_iter()
                    .map(|(gene, err)| (gene.key.clone(), (gene, err)))
                    .into_group_map(),
            )
        })
    {
        writeln!(error, "SPECIES: {species}").unwrap();
        for (gene, errors) in errors {
            writeln!(error, "GENE: {gene}").unwrap();
            for (gene, err) in errors {
                writeln!(error, "ERROR FOR GENE:\n{species}\t{gene}\t{err}\n").unwrap();
            }
        }
    }

    writeln!(
        output,
        "// @generated\n#![allow(non_snake_case,non_upper_case_globals)]\nuse std::sync::OnceLock;\nuse super::shared::{{Germlines, Species}};"
    )
    .unwrap();
    writeln!(output, "/// Get the germlines for any of the available species. See the main documentation for which species have which data available.").unwrap();
    let mut found_species = Vec::new();
    let mut found_germlines: Vec<(Species, Germlines)> = grouped.into_iter().collect();
    found_germlines.sort_unstable_by_key(|g| g.0);
    for (species, germlines) in found_germlines {
        writeln!(
            docs,
            "## {} / {}

| Kind | V | J | C |
|------|---|---|---|
|IGHV{}
|IGKV{}
|IGLV{}
|IGIV{}

_Number of genes / number of alleles_
",
            species.scientific_name(),
            species.common_name(),
            germlines.h.doc_row(),
            germlines.k.doc_row(),
            germlines.l.doc_row(),
            germlines.i.doc_row(),
        )
        .unwrap();
        found_species.push(species);

        let mut file =
            std::fs::File::create(format!("rustyms/src/imgt/germlines/{species}.bin")).unwrap();
        file.write_all(&bincode::serialize::<Germlines>(&germlines).unwrap())
            .unwrap();
    }
    // germlines
    writeln!(
        output,
        "pub fn germlines(species: Species) -> Option<&'static Germlines> {{match species {{"
    )
    .unwrap();

    for species in &found_species {
        writeln!(output, "Species::{0} => Some(lock_{0}()),", species.ident()).unwrap();
    }
    writeln!(output, "_=>None}}}}").unwrap();
    // all_germlines
    writeln!(
        output,
"/// Get all germlines in one iterator, see the main documentation for more information about the available germlines
pub fn all_germlines() -> impl std::iter::Iterator<Item = &'static Germlines> {{"
    )
    .unwrap();
    let mut first = true;
    for species in &found_species {
        if first {
            first = false;
            writeln!(output, "std::iter::once(lock_{}())", species.ident()).unwrap();
        } else {
            writeln!(
                output,
                ".chain(std::iter::once(lock_{}()))",
                species.ident()
            )
            .unwrap();
        }
    }
    writeln!(output, "}}").unwrap();
    // par_germlines
    writeln!(
        output,
"/// Get all germlines in one parallel iterator, see the main documentation for more information about the available germlines
#[cfg(feature = \"rayon\")]
use rayon::prelude::*;
#[cfg(feature = \"rayon\")]
pub fn par_germlines() -> impl rayon::prelude::ParallelIterator<Item = &'static Germlines> {{"
    )
    .unwrap();
    let mut first = true;
    for species in &found_species {
        if first {
            first = false;
            writeln!(output, "rayon::iter::once(lock_{}())", species.ident()).unwrap();
        } else {
            writeln!(
                output,
                ".chain(rayon::iter::once(lock_{}()))",
                species.ident()
            )
            .unwrap();
        }
    }
    writeln!(output, "}}").unwrap();

    for species in &found_species {
        writeln!(
            output,
"static LOCK_{0}: OnceLock<Germlines> = OnceLock::new();
fn lock_{0}()->&'static Germlines{{LOCK_{0}.get_or_init(|| {{bincode::deserialize(include_bytes!(\"{species}.bin\")).unwrap()}})}}",
            species.ident(),
        )
        .unwrap();
    }
}

fn complement(s: String) -> String {
    let map = HashMap::from([
        (b'a', b't'),
        (b't', b'a'),
        (b'c', b'g'),
        (b'g', b'c'),
        (b'n', b'n'),
    ]);
    String::from_utf8(
        s.as_bytes()
            .iter()
            .map(|c| {
                *map.get(c)
                    .unwrap_or_else(|| panic!("Invalid sequence: {} in `{s}`", char::from(*c)))
            })
            .rev()
            .collect(),
    )
    .unwrap()
}

fn translate(s: &str) -> Result<(&str, Vec<AminoAcid>), String> {
    if s.len() < 3 {
        Ok((s, Vec::new()))
    } else {
        Ok((
            s,
            (0..=s.len() - 3)
                .step_by(3)
                .filter_map(|chunk| {
                    invert(
                        AminoAcid::from_dna(&s[chunk..chunk + 3])
                            .map_err(|_| format!("Not a codon {}", &s[chunk..chunk + 3])),
                    )
                })
                .collect::<Result<Vec<AminoAcid>, String>>()?,
        ))
    }
}

fn invert<T, E>(x: Result<Option<T>, E>) -> Option<Result<T, E>> {
    match x {
        Ok(None) => None,
        Ok(Some(a)) => Some(Ok(a)),
        Err(e) => Some(Err(e)),
    }
}

fn find_possible_n_glycan_locations(sequence: &[AminoAcid]) -> Vec<usize> {
    let mut result = Vec::new();
    for (index, aa) in sequence.windows(3).enumerate() {
        if let (AminoAcid::N, AminoAcid::S | AminoAcid::T) = (aa[0], aa[2]) {
            if aa[1] != AminoAcid::P {
                result.push(index);
            }
        }
    }
    result
}

fn fix_j(
    j: (Vec<AminoAcid>, Location, String),
    cdr3_length: usize,
) -> (Vec<SequenceRegion>, Vec<(Annotation, usize)>) {
    let (cdr3_loc, fr4_loc) =
        j.1.splice(cdr3_length)
            .expect("CDR3 should fit in full FR4 of J gene");
    let cdr3 = (
        j.0[..cdr3_length].to_vec(),
        cdr3_loc,
        j.2[..cdr3_length].to_owned(),
    );
    let fr4 = (
        j.0[cdr3_length..].to_vec(),
        fr4_loc,
        j.2[cdr3_length..].to_owned(),
    );

    let mut annotations = Vec::new();
    if fr4.0[0] == AminoAcid::W {
        annotations.push((Annotation::Tryptophan, cdr3_length));
    } else if fr4.0[0] == AminoAcid::F {
        annotations.push((Annotation::Phenylalanine, cdr3_length));
    }
    if fr4.0[1] == AminoAcid::G {
        annotations.push((Annotation::Glycine, cdr3_length + 1));
    }
    if fr4.0[3] == AminoAcid::G {
        annotations.push((Annotation::Glycine, cdr3_length + 3));
    }

    (
        vec![(shared::Region::CDR3, cdr3), (shared::Region::FR4, fr4)],
        annotations,
    )
}
