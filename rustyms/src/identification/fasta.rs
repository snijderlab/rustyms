use crate::{
    error::{Context, CustomError},
    peptide::SemiAmbiguous,
    LinearPeptide, SequenceElement,
};
use serde::{Deserialize, Serialize};
use std::{
    io::{BufRead, BufReader},
    path::Path,
};

use super::{AminoAcid, IdentifiedPeptide, MetaData};

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
    pub fn parse_file(path: impl AsRef<Path>) -> Result<Vec<Self>, CustomError> {
        let path = path.as_ref();
        let file = std::fs::File::open(path).map_err(|_| {
            CustomError::error(
                "Failed reading fasta file",
                "Error occurred while opening the file",
                Context::show(path.to_string_lossy()),
            )
        })?;
        let reader = BufReader::new(file);
        Self::parse_reader(reader, Some(path))
    }
    /// Parse a single fasta file from a reader
    /// # Errors
    /// A custom error when it is not a valid fasta file
    pub fn parse_reader(
        reader: impl BufRead,
        path: Option<&Path>,
    ) -> Result<Vec<Self>, CustomError> {
        let mut sequences = Vec::new();
        let mut last_header = None;
        let mut last_sequence: Vec<SequenceElement<SemiAmbiguous>> = Vec::new();

        for (line_index, line) in reader.lines().enumerate() {
            let line = line.map_err(|_| {
                CustomError::error(
                    "Failed reading fasta file",
                    format!("Error occurred while reading line {}", line_index + 1),
                    path.map_or(Context::None, |p| Context::show(p.to_string_lossy())),
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
                last_sequence.extend(
                    line.char_indices()
                        .map(|(i, c)| {
                            c.try_into()
                                .map(|aa: AminoAcid| SequenceElement::new(aa.into(), None))
                                .map_err(|()| {
                                    CustomError::error(
                                        "Failed reading fasta file",
                                        "Character is not an amino acid",
                                        Context::line(Some(line_index), &line, i, 1),
                                    )
                                })
                        })
                        .collect::<Result<Vec<SequenceElement<_>>, _>>()?,
                );
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

#[test]
#[allow(clippy::missing_panics_doc)]
fn empty_lines() {
    let file = ">A\naaa\n\naaa";
    let fasta = FastaData::parse_reader(BufReader::new(file.as_bytes()), None).unwrap();
    assert_eq!(fasta.len(), 1);
    assert_eq!(
        fasta[0].peptide,
        LinearPeptide::pro_forma("AAAAAA", None).unwrap()
    );
}
