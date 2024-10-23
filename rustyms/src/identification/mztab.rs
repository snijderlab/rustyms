use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    ops::Range,
    str::FromStr,
};

use flate2::bufread::GzDecoder;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    helper_functions::{check_extension, explain_number_error},
    identification::{IdentifiedPeptide, MetaData},
    modification::SimpleModification,
    ontologies::CustomDatabase,
    system::{usize::Charge, MassOverCharge, Time},
    AminoAcid, LinearPeptide, SemiAmbiguous,
};

use super::{Modification, SequencePosition, SloppyParsingParameters};

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
    /// The search engines that identified this peptide, alongside their identified score and the CV term describing the score
    pub search_engine: Vec<(CVTerm, Option<f64>, String)>,
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
    ) -> Result<Box<dyn Iterator<Item = Result<Self, CustomError>> + '_>, CustomError> {
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
            Ok(Box::new(Self::parse_mztab_raw(
                BufReader::new(GzDecoder::new(BufReader::new(file))),
                custom_database,
            )))
        } else {
            Ok(Box::new(Self::parse_mztab_raw(
                BufReader::new(file),
                custom_database,
            )))
        }
    }

    /// Parse a MZTab file directly from a buf reader
    pub fn parse_mztab_raw<'a, T: std::io::BufRead + 'a>(
        reader: T,
        custom_database: Option<&'a CustomDatabase>,
    ) -> impl Iterator<Item = Result<Self, CustomError>> + 'a {
        let time_unit = 0.0;
        let search_engine_score_type: Vec<CVTerm> = Vec::new();
        let modifications: Vec<SimpleModification> = Vec::new();
        let mut peptide_header: Option<Vec<String>> = None;

        parse_mztab_raw(reader).filter_map(move |item| {
            item.transpose().and_then(|item| match item {
                Ok(MZTabLine::MTD(line_index, line, fields)) => {
                    if fields.len() == 3 {
                        // TODO: finish handling metadata
                        match line[fields[1].clone()].to_ascii_lowercase().as_str() {
                            "columnit_psm" => (),
                            m if m.starts_with("fixed_mod[") && m.ends_with(']') => (),
                            m if m.starts_with("variable_mod[") && m.ends_with(']') => (),
                            m if m.starts_with("psm_search_engine_score[") && m.ends_with(']') => (),
                            m if m.starts_with("ms_run[") && m.ends_with("]-location") => (),
                            _ => ()
                        }
                        None
                    } else {
                        Some(Err(CustomError::error("Invalid MTD line", "MTD lines should contain three columns (the tag, key, and value)", Context::full_line(line_index, line))))
                    }
                },
                Ok(MZTabLine::PSH(line_index, line, fields)) => {
                    let header = fields
                        .into_iter()
                        .map(|field| line[field].to_ascii_lowercase())
                        .collect_vec();
                    // optional: opt_*, reliability, uri,
                    // not checked: search_engine_score[n]
                    for required in [
                        "sequence",
                        "psm_id",
                        "accession",
                        "unique",
                        "database",
                        "database_version",
                        "search_engine",
                        "modifications",
                        "retention_time",
                        "charge",
                        "exp_mass_to_charge",
                        "spectra_ref",
                        "pre",
                        "post",
                        "start",
                        "end",
                    ] {
                        if !header.contains(&required.to_string()) {
                            return Some(Err(CustomError::error(
                                "Invalid peptide table",
                                format!("The required column '{required}' is not present"),
                                Context::full_line(line_index, line),
                            )));
                        }
                    }
                    peptide_header = Some(header);
                    None
                }
                Ok(MZTabLine::PSM(line_index, line, fields)) => Some(
                    PSMLine::new(line_index, peptide_header.as_deref(), &line, &fields)
                        .and_then(|line| Self::from_line(line, custom_database)),
                ),
                Err(e) => Some(Err(e)),
            })
        })
    }

    fn from_line<'a>(
        line: PSMLine<'a>,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, CustomError> {
        let (mod_column, mod_range) =
            line.required_column("modifications", line.line, line.line_index)?;
        let mut mod_index = mod_range.start;
        let modifications: Vec<(usize, SimpleModification)> = mod_column
            .split(',')
            .flat_map(|definition| {
                let pair = definition
                    .split_once('-')
                    .map(|(pos, _)| Ok((
                        pos.parse::<usize>().map_err(|err| CustomError::error(
                        "Invalid modification position",
                        format!("The position {}", explain_number_error(&err)),
                        Context::line_range(Some(line.line_index), line.line, mod_range.clone()),
                    ))?, 
                    SimpleModification::try_from(
                        line.line, 
                        mod_index+1+pos.len()..mod_index+definition.len(), 
                        &mut Vec::new(), 
                        &mut Vec::new(), 
                        custom_database)?.defined()
                        .ok_or_else(
                            || CustomError::error(
                                "Invalid modification", 
                                "A modification should be a fully defined modification, no cross-link or ambiguous modification", 
                                Context::line_range(Some(line.line_index), line.line, mod_range.clone())))?)))
                    .ok_or_else(
                        || CustomError::error(
                            "Invalid modification", 
                            "A modification should be the position followed by a hyphen ('-') followed by the modification", 
                            Context::line_range(Some(line.line_index), line.line, mod_range.clone())));
                            mod_index += definition.len() + 1;
                            pair
            })
            .collect::<Result<Vec<_>, CustomError>>()?;

        Ok(Self {
            peptide: {
                let mut peptide = LinearPeptide::sloppy_pro_forma(
                    line.line,
                    line.required_column("sequence", line.line, line.line_index)?
                        .1,
                    custom_database,
                    SloppyParsingParameters {
                        allow_unwrapped_modifications: true,
                        ..Default::default()
                    },
                )?;
                for (location, modification) in modifications {
                    match location {
                        0 => peptide.set_simple_n_term(Some(modification)),
                        c if c == peptide.len() + 1 => peptide.set_simple_c_term(Some(modification)),
                        i => peptide.sequence_mut()[i - 1].add_simple_modification(modification),
                    }
                }
                peptide
            },
            psm_id: line
                .required_column("psm_id", line.line, line.line_index)?
                .0
                .parse()
                .map_err(|err| {
                    CustomError::error(
                        "Invalid MZTab PSM_ID",
                        format!("The PSM_ID {}", explain_number_error(&err)),
                        Context::line_range(
                            Some(line.line_index),
                            line.line,
                            line.optional_column("psm_id").unwrap().1,
                        ),
                    )
                })?,
            accession: line
                .optional_column("accession")
                .and_then(|(v, _)| (v.to_ascii_lowercase() != "null").then(|| v.to_string())),
            unique: line
                .optional_column("unique")
                .and_then(|(v, _)| (v.to_ascii_lowercase() != "null").then(|| v == "1")),
            database: line
                .optional_column("database")
                .and_then(|(v, _)| (v.to_ascii_lowercase() != "null").then(|| v.to_string())),
            database_version: line
                .optional_column("database_version")
                .and_then(|(v, _)| (v.to_ascii_lowercase() != "null").then(|| v.to_string())),
            search_engine: {
                let (value, range) =
                    line.required_column("search_engine", line.line, line.line_index)?;

                if value.trim().to_ascii_lowercase() == "null" {
                    Vec::new()
                } else {
                    value
                        .split('|')
                        .enumerate()
                        .map(|(i, s)| {
                            line.optional_column(&format!("search_engine_score[{}]", i + 1))
                                .and_then(|(v, _)| {
                                    (v.to_ascii_lowercase() != "null").then(|| {
                                        v.parse::<f64>().map_err(|err| {
                                            CustomError::error(
                                                "Invalid MZTab search engine score",
                                                format!(
                                        "The search engine score can not be parsed as f64: {err}"
                                    ),
                                                Context::line_range(
                                                    Some(line.line_index),
                                                    line.line,
                                                    line.optional_column(&format!(
                                                        "search_engine_score[{}]",
                                                        i + 1
                                                    ))
                                                    .unwrap()
                                                    .1,
                                                ),
                                            )
                                        })
                                    })
                                })
                                .transpose()
                                .and_then(|score| {
                                    CVTerm::from_str(s)
                                        .map_err(|e| {
                                            e.with_context(Context::line_range(
                                                Some(line.line_index),
                                                line.line,
                                                range.clone(),
                                            ))
                                        })
                                        .map(|engine| (engine, score, String::new()))
                                })
                            // TODO: get proper CV term
                        })
                        .collect::<Result<Vec<_>, CustomError>>()?
                }
            },
            reliability: line
                .optional_column("reliability")
                .map(|(v, range)| match v {
                    "1" => Ok(PSMReliability::High),
                    "2" => Ok(PSMReliability::Medium),
                    "3" => Ok(PSMReliability::Poor),
                    _ => Err(CustomError::error(
                        "Invalid PSM reliability",
                        format!("A reliability should be 1, 2, or 3, '{v}' is invalid"),
                        Context::line_range(Some(line.line_index), line.line, range),
                    )),
                })
                .transpose()?,
            rt: line
                .optional_column("retention_time")
                .and_then(|(v, r)| {
                    (v.to_ascii_lowercase() != "null").then(|| {
                        v.parse::<f64>()
                            .map_err(|err| {
                                CustomError::error(
                                    "Invalid MZTab retention time",
                                    format!("The retention time can not be parsed as f64: {err}"),
                                    Context::line_range(Some(line.line_index), line.line, r),
                                )
                            })
                            .map(|v| Time::new::<crate::system::s>(v)) // TODO: use correct unit
                    })
                })
                .transpose()?,
            z: {
                let (value, range) = line.required_column("charge", line.line, line.line_index)?;

                if value.trim().to_ascii_lowercase() == "null" {
                    Charge::new::<crate::system::e>(1)
                } else {
                    value
                        .trim_end_matches(".0")
                        .parse::<usize>()
                        .map_err(|err| {
                            CustomError::error(
                                "Invalid MZTab charge",
                                format!("The charge {}", explain_number_error(&err)),
                                Context::line_range(Some(line.line_index), line.line, range),
                            )
                        })
                        .map(|v| Charge::new::<crate::system::e>(v))?
                }
            },
            exp_mz: line
                .optional_column("exp_mass_to_charge")
                .and_then(|(v, r)| {
                    (v.to_ascii_lowercase() != "null").then(|| {
                        v.parse::<f64>()
                            .map_err(|err| {
                                CustomError::error(
                                    "Invalid MZTab experimental mz",
                                    format!("The experimental mz can not be parsed as f64: {err}"),
                                    Context::line_range(Some(line.line_index), line.line, r),
                                )
                            })
                            .map(|v| MassOverCharge::new::<crate::system::mz>(v))
                    })
                })
                .transpose()?,
            uri: line.optional_column("uri").map(|(v, _)| v.to_string()),
            spectra_ref: line
                .optional_column("spectra_ref")
                .map(|(v, _)| v.to_string())
                .into_iter()
                .collect(),
            preceding_aa: line
                .required_column("pre", line.line, line.line_index)?
                .0
                .parse()
                .map_err(|()| {
                    CustomError::error(
                        "Invalid preceding amino acid",
                        "The pre column should contain null, -, or an aminoacid",
                        Context::line_range(
                            Some(line.line_index),
                            line.line,
                            line.optional_column("pre").unwrap().1,
                        ),
                    )
                })?,
            following_aa: line
                .required_column("post", line.line, line.line_index)?
                .0
                .parse()
                .map_err(|()| {
                    CustomError::error(
                        "Invalid following amino acid",
                        "The post column should contain null, -, or an aminoacid",
                        Context::line_range(
                            Some(line.line_index),
                            line.line,
                            line.optional_column("post").unwrap().1,
                        ),
                    )
                })?,
            start: line
                .optional_column("start")
                .and_then(|(v, r)| {
                    (v.to_ascii_lowercase() != "null").then(|| {
                        v.parse::<usize>().map_err(|err| {
                            CustomError::error(
                                "Invalid MZTab start",
                                format!("The start {}", explain_number_error(&err)),
                                Context::line_range(Some(line.line_index), line.line, r),
                            )
                        })
                    })
                })
                .transpose()?,
            end: line
                .optional_column("end")
                .and_then(|(v, r)| {
                    (v.to_ascii_lowercase() != "null").then(|| {
                        v.parse::<usize>().map_err(|err| {
                            CustomError::error(
                                "Invalid MZTab end",
                                format!("The end {}", explain_number_error(&err)),
                                Context::line_range(Some(line.line_index), line.line, r),
                            )
                        })
                    })
                })
                .transpose()?,
            local_confidence: line
                .optional_column("op_ms_run[1]_aa_scores")
                .map(|(v, r)| {
                    v.split(',')
                        .map(|score| {
                            score.parse::<f64>().map_err(|err| {
                                CustomError::error(
                                    "Invalid MZTab local confidence",
                                    format!("The local confidence can not be parsed: {err}"),
                                    Context::line_range(
                                        Some(line.line_index),
                                        line.line,
                                        r.clone(),
                                    ),
                                )
                            })
                        })
                        .collect()
                })
                .transpose()?,
            additional: HashMap::new(), // TODO: fill
        })
    }
}

struct PSMLine<'a> {
    line_index: usize,
    header: &'a [String],
    pub line: &'a str,
    fields: &'a [Range<usize>],
}

impl<'a> PSMLine<'a> {
    fn new(
        line_index: usize,
        header: Option<&'a [String]>,
        line: &'a str,
        fields: &'a [Range<usize>],
    ) -> Result<Self, CustomError> {
        let header = header.ok_or_else(|| {
            CustomError::error(
                "Missing PSH line",
                "The PSH peptide header line should precede any PSM line",
                Context::full_line(line_index, line),
            )
        })?;
        if header.len() == fields.len() {
            Ok(Self {
                line_index,
                header,
                line,
                fields,
            })
        } else {
            Err(CustomError::error(
                "Invalid PSM line",
                "This PSM line does not have the same number of columns as the PSH line",
                Context::full_line(line_index, line),
            ))
        }
    }

    fn optional_column(&self, column: &str) -> Option<(&str, Range<usize>)> {
        self.header
            .iter()
            .position(|h| h == column)
            .map(|i| (&self.line[self.fields[i].clone()], self.fields[i].clone()))
    }

    fn required_column(
        &self,
        column: &str,
        line: &str,
        line_index: usize,
    ) -> Result<(&str, Range<usize>), CustomError> {
        self.optional_column(column).ok_or(CustomError::error(
            "Missing column",
            format!("The column '{column}' is required but not present"),
            Context::full_line(line_index, line),
        ))
    }
}

impl From<MZTabData> for IdentifiedPeptide {
    fn from(value: MZTabData) -> Self {
        Self {
            score: if value.search_engine.is_empty() {
                None
            } else {
                Some(
                    value
                        .search_engine
                        .iter()
                        .filter_map(|(_, s, _)| *s)
                        .sum::<f64>()
                        / value.search_engine.len() as f64,
                )
            },
            metadata: MetaData::MZTab(value),
        }
    }
}

/// A flanking residue for a sequence, N or C terminal agnostic
#[derive(Clone, PartialEq, Eq, Debug, Default, Serialize, Deserialize)]
pub enum FlankingResidue {
    /// The flanking residue is unknown (for example in de novo data)
    #[default]
    Unknown,
    /// The residue is terminal
    Terminal,
    /// The flanking residue
    AminoAcid(AminoAcid),
}

impl std::str::FromStr for FlankingResidue {
    type Err = <AminoAcid as std::str::FromStr>::Err;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim() {
            "null" => Ok(Self::Unknown),
            "-" => Ok(Self::Terminal),
            _ => AminoAcid::from_str(value).map(Self::AminoAcid),
        }
    }
}

/// A CV term
#[derive(Clone, PartialEq, Eq, Debug, Default, Serialize, Deserialize)]
pub struct CVTerm {
    /// The ontology
    pub ontology: String,
    /// The id within the ontology
    pub id: String,
    /// The human name for the term
    pub term: String,
    /// Additional comments on the term, eg additional specification
    pub comment: String,
}

impl std::str::FromStr for CVTerm {
    type Err = CustomError;
    fn from_str(value: &str) -> Result<Self, CustomError> {
        let value = value.trim();
        if value.starts_with('[') && value.ends_with(']') {
            let value = &value[1..value.len() - 1];
            let mut split = value.splitn(4, ',');
            Ok(Self {
                ontology: split.next().unwrap_or_default().to_string(),
                id: split.next().unwrap_or_default().to_string(),
                term: split.next().unwrap_or_default().to_string(),
                comment: split.next().unwrap_or_default().to_string(),
            })
        } else {
            Err(CustomError::error(
                "Invalid CV term",
                "A CV term should be encolsed by '[]'",
                Context::None,
            ))
        }
    }
}

/// The reliability of a PSM
#[allow(missing_docs)]
#[derive(Clone, PartialEq, Eq, Debug, Default, Serialize, Deserialize)]
pub enum PSMReliability {
    High,
    Medium,
    #[default]
    Poor,
}

/// A basic structure for an MZTab file line
enum MZTabLine {
    /// Metadata line
    MTD(usize, String, Vec<Range<usize>>),
    /// Peptide header line
    PSH(usize, String, Vec<Range<usize>>),
    /// Peptide line, stored as hashmap with the columns names from PSH
    PSM(usize, String, Vec<Range<usize>>),
}

/// Parse a MZTab file
/// # Errors
/// If the file is not a valid MZTab file
fn parse_mztab_raw<T: BufRead>(
    reader: T,
) -> impl Iterator<Item = Result<Option<MZTabLine>, CustomError>> {
    reader.lines().enumerate().map(move |(line_index, line)| {
        line.map_err(|err| {
            CustomError::error(
                "Could not read line",
                err,
                Context::full_line(line_index, "(failed)"),
            )
        })
        .and_then(|line| {
            if line.trim().is_empty() {
                Ok(None)
            } else {
                crate::csv::csv_separate(&line, b'\t').map(|fields| {
                    match &line[fields[0].clone()] {
                        "MTD" => Some(MZTabLine::MTD(line_index, line, fields)),
                        "PSH" => Some(MZTabLine::PSH(line_index, line, fields)),
                        "PSM" => Some(MZTabLine::PSM(line_index, line, fields)),
                        _ => None,
                    }
                })
            }
        })
    })
}
