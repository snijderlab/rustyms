use std::io::{BufRead, BufReader};

use crate::{
    error::{Context, CustomError},
    ComplexPeptide, LinearPeptide, SequenceElement,
};

use super::{IdentifiedPeptide, MetaData};

/// A single parsed line of a fasta file
#[allow(missing_docs)]
#[derive(Debug)]
pub struct FastaData {
    pub id: String,
    pub full_header: String,
    pub peptide: LinearPeptide,
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
        let mut last_sequence: Vec<SequenceElement> = Vec::new();

        for (line_number, line) in reader.lines().enumerate() {
            let line = line.map_err(|_| {
                CustomError::error(
                    "Failed reading fasta file",
                    format!("Error occurred while reading line {}", line_number + 1),
                    Context::show(path),
                )
            })?;
            #[allow(clippy::manual_strip)]
            if line.starts_with('>') {
                if let Some((id, full_header)) = last_header {
                    sequences.push(Self {
                        id,
                        full_header,
                        peptide: LinearPeptide {
                            global: Vec::new(),
                            labile: Vec::new(),
                            n_term: None,
                            c_term: None,
                            sequence: last_sequence,
                            ambiguous_modifications: Vec::new(),
                            charge_carriers: None,
                        },
                    });
                }
                last_header = Some((
                    line[1..line.find(' ').unwrap_or(line.len())].to_string(),
                    line[1..].to_string(),
                ));
                last_sequence = Vec::new();
            } else {
                let parsed = ComplexPeptide::pro_forma(&line)
                    .map_err(|e| e.overwrite_line_number(line_number))?
                    .assume_linear()
                    .sequence;
                last_sequence.extend(parsed);
            }
        }
        if let Some((id, full_header)) = last_header {
            sequences.push(Self {
                id,
                full_header,
                peptide: LinearPeptide {
                    global: Vec::new(),
                    labile: Vec::new(),
                    n_term: None,
                    c_term: None,
                    sequence: last_sequence,
                    ambiguous_modifications: Vec::new(),
                    charge_carriers: None,
                },
            });
        }

        Ok(sequences)
    }
}

impl From<FastaData> for IdentifiedPeptide {
    fn from(value: FastaData) -> Self {
        Self {
            peptide: value.peptide.clone(),
            local_confidence: None,
            score: None,
            metadata: MetaData::Fasta(value),
        }
    }
}
