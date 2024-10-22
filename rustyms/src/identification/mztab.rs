use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
};

use flate2::bufread::GzDecoder;
use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    helper_functions::check_extension,
    identification::{IdentifiedPeptide, MetaData},
    ontologies::CustomDatabase,
    system::{usize::Charge, MassOverCharge, Time},
    AminoAcid, LinearPeptide, SemiAmbiguous,
};

use super::SloppyParsingParameters;

/// Peptide data from a MZTab file
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct MZTabData {
    /// The peptide's sequence corresponding to the PSM
    pub peptide: LinearPeptide<SemiAmbiguous>,
    /// A unique identifier for a PSM within the file. If a PSM can be matched to
    /// multiple proteins, the same PSM should be represented on multiple rows with
    /// different accessions and the same PSM_ID.
    pub psm_id: usize,
    /// The protein's accession the corresponding peptide sequence (coming from the
    /// PSM) is associated with.
    pub accession: Option<String>,
    /// Indicates whether the peptide sequence (coming from the PSM) is unique for
    /// this protein in respect to the searched database.
    pub unique: Option<bool>,
    /// The protein database used for the search (could theoretically come from a
    /// different species) and the peptide sequence comes from.
    pub database: Option<String>,
    /// The protein database's version
    pub database_version: Option<String>,
    /// The search engines that identified this peptide, alongside their identified score
    pub search_engine: Vec<(String, f64)>,
    /// If available the estaimated reliability of the PSM.
    pub reliability: Option<PSMReliability>,
    // TODO: parse modifications columns and inject in peptide
    /// The retention time for this peptide.
    pub rt: Option<Time>,
    /// The charge for this peptide.
    pub z: Charge,
    /// The experimental mz
    pub exp_mz: Option<MassOverCharge>,
    /// A URI pointing to the PSM's entry in the experiment it was identified in (e.g. the peptide’s PRIDE entry).
    pub uri: Option<String>,
    // TODO: fix type
    pub spectra_ref: Vec<String>,
    /// The amino acide before this peptide
    pub preceding_aa: FlankingResidue,
    /// The amino acide after this peptide
    pub following_aa: FlankingResidue,
    /// The start of this peptide in the containing protein (0-based)
    pub start: Option<usize>,
    /// The end of this peptide in the containing protein (0-based)
    pub end: Option<usize>,
    /// Casanovo specific additional metadata with the amino acid confidence
    pub local_confidence: Option<Vec<f64>>,
    /// Any additional metadata
    pub additional: HashMap<String, String>,
}

impl MZTabData {
    /// Parse a MZTab file.
    /// # Errors
    /// If the file is not in the correct format
    pub fn parse_mztab(
        path: impl AsRef<std::path::Path>,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<impl Iterator<Item = Result<MZTabData, CustomError>> + '_, CustomError> {
        let time_unit = 0.0;
        Ok(parse_mztab(path)?.filter_map(move |item| {
            item.transpose().and_then(|item| match item {
                Ok(MZTabLine::MTD(key, value)) => None,
                Ok(MZTabLine::PSM(fields)) => Some(Self::from_line(&fields, custom_database)),
                Err(e) => Some(Err(e)),
            })
        }))
    }

    fn from_line(
        fields: &HashMap<String, String>,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, CustomError> {
        // Ok(Self {
        //     peptide: LinearPeptide::sloppy_pro_forma(
        //         &fields["peptide"],
        //         0..fields["peptide"].len(),
        //         custom_database,
        //         SloppyParsingParameters::default(),
        //     )?,
        // })
        todo!()
    }
}

impl From<MZTabData> for IdentifiedPeptide {
    fn from(value: MZTabData) -> Self {
        Self {
            score: if value.search_engine.is_empty() {
                None
            } else {
                Some(
                    value.search_engine.iter().map(|(_, s)| *s).sum::<f64>()
                        / value.search_engine.len() as f64,
                )
            },
            metadata: MetaData::MZTab(value),
        }
    }
}

/// A flanking residue for a sequence, N or C terminal agnostic
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub enum FlankingResidue {
    /// The flanking residue is unknown (for example in de novo data)
    #[default]
    Unknown,
    /// The residue is terminal
    Terminal,
    /// The flanking residue
    AminoAcid(AminoAcid),
}

/// The reliability of a PSM
#[allow(missing_docs)]
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub enum PSMReliability {
    High,
    Medium,
    #[default]
    Poor,
}

/// A basic structure for an MZTab file line
enum MZTabLine {
    /// Metadata line
    MTD(String, String),
    /// Peptide line, stored as hashmap with the columns names from PSH
    PSM(HashMap<String, String>),
}

/// Parse a MZTab file
/// # Errors
/// If the file is not a valid MZTab file
fn parse_mztab(
    path: impl AsRef<std::path::Path>,
) -> Result<Box<dyn Iterator<Item = Result<Option<MZTabLine>, CustomError>>>, CustomError> {
    let file = File::open(path.as_ref()).map_err(|e| {
        CustomError::error(
            "Could not open file",
            e,
            crate::error::Context::Show {
                line: path.as_ref().to_string_lossy().to_string(),
            },
        )
    })?;
    if check_extension(path, "gz") {
        Ok(Box::new(parse_mztab_raw(BufReader::new(GzDecoder::new(
            BufReader::new(file),
        )))))
    } else {
        Ok(Box::new(parse_mztab_raw(BufReader::new(file))))
    }
}

/// Parse a MZTab file
/// # Errors
/// If the file is not a valid MZTab file
fn parse_mztab_raw<T: BufRead>(
    reader: T,
) -> impl Iterator<Item = Result<Option<MZTabLine>, CustomError>> {
    let mut peptide_header: Option<Vec<String>> = None;
    reader
        .lines()
        .enumerate()
        .map(move |(line_index, line)| {
            line.map_err(|err| {
                CustomError::error(
                    "Could not read line",
                    err,
                    Context::full_line(line_index, "(failed)"),
                )
            })
            .and_then(|line| {
                crate::csv::csv_separate(&line, b'\t').and_then(|fields| {
                    Ok(match &line[fields[0].clone()] {
                        "MTD" => Some(MZTabLine::MTD(
                            line[fields[1].clone()].to_string(),
                            line[fields[2].clone()].to_string(),
                        )),
                        "PSH" => {
                            peptide_header = Some(
                                fields
                                    .iter()
                                    .skip(1)
                                    .map(|r| line[r.clone()].to_string())
                                    .collect(),
                            );
                            None
                        }
                        "PSM" => Some(MZTabLine::PSM(
                            peptide_header
                                .as_ref()
                                .ok_or_else(|| CustomError::error("Missing PSH line", "The PSH line has to be located before the PSM lines in a MZTab file", Context::full_line(line_index, line.to_string())))?
                                .iter()
                                .zip(fields.iter().skip(1))
                                .map(|(header, field)| {
                                    (header.clone(), line[field.clone()].to_string())
                                })
                                .collect(),
                        )),
                        _ => None,
                    })
                })
            })
        })
}
