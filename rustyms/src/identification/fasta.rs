use crate::{
    error::{Context, CustomError},
    peptide::SemiAmbiguous,
    CompoundPeptidoform, LinearPeptide, SequenceElement,
};
use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader};

use super::{IdentifiedPeptide, MetaData};

/// A single parsed line of a fasta file
#[allow(missing_docs)]
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct FastaData {
    pub id: String,
    pub full_header: String,
    pub peptide: LinearPeptide<SemiAmbiguous>,
}

impl FastaData {
    /// Parse a single fasta file
    /// # Errors
    /// A custom error when it is not a valid fasta file
    pub fn parse_file(path: &str) -> Result<Vec<Self>, CustomError> {
        let file = std::fs::File::open(path).map_err(|_| {
            CustomError::error(
                "Failed reading fasta file",
                "Error occurred while opening the file",
                Context::show(path),
            )
        })?;
        let reader = BufReader::new(file);
        let mut sequences = Vec::new();
        let mut last_header = None;
        let mut last_sequence: Vec<SequenceElement<SemiAmbiguous>> = Vec::new();

        for (line_index, line) in reader.lines().enumerate() {
            let line = line.map_err(|_| {
                CustomError::error(
                    "Failed reading fasta file",
                    format!("Error occurred while reading line {}", line_index + 1),
                    Context::show(path),
                )
            })?;
            #[allow(clippy::manual_strip)]
            if line.starts_with('>') {
                if let Some((id, full_header)) = last_header {
                    sequences.push(Self {
                        id,
                        full_header,
                        peptide: last_sequence.into(),
                    });
                }
                last_header = Some((
                    line[1..line.find(' ').unwrap_or(line.len())].to_string(),
                    line[1..].to_string(),
                ));
                last_sequence = Vec::new();
            } else {
                let parsed = CompoundPeptidoform::pro_forma(&line, None)
                    .map_err(|e| e.overwrite_line_number(line_index))?
                    .singular()
                    .ok_or_else(||
                        CustomError::error(
                            "Invalid fasta sequence", 
                            "A sequence in a Fasta file is assumed to be a single peptide and not a chimeric compound peptidoform",
                             Context::full_line(line_index, line.clone())))
                    .and_then(|p|
                        p.singular().ok_or_else(||
                            CustomError::error(
                                "Invalid fasta sequence", 
                                "A sequence in a Fasta file is assumed to be a single peptide and not a cross linked peptidoform", 
                                Context::full_line(line_index, line.clone())))
                    )
                    .and_then(|p| p.into_semi_ambiguous().ok_or_else(||
                        CustomError::error(
                            "Invalid fasta sequence", 
                            "A sequence in a Fasta file is assumed to be a simple sequence only consisting of amino acids although this implementation allows simple modifications as well", 
                            Context::full_line(line_index, line.clone())))
                    )?;
                last_sequence.extend(parsed.sequence().iter().cloned());
            }
        }
        if let Some((id, full_header)) = last_header {
            sequences.push(Self {
                id,
                full_header,
                peptide: last_sequence.into(),
            });
        }

        Ok(sequences)
    }
}

impl From<FastaData> for IdentifiedPeptide {
    fn from(value: FastaData) -> Self {
        Self {
            score: None,
            metadata: MetaData::Fasta(value),
        }
    }
}
