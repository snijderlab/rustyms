use crate::{error::Context, error::CustomError, AminoAcid, ComplexPeptide, LinearPeptide};

use super::{
    common_parser::Location,
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
};

/// The file format for any opair format, determining the existence and location of all possible columns
#[derive(Debug, Clone)]
pub struct OpairFormat {
    file_name: usize,
    scan_number: usize,
    rt: usize,
    precursor_scan_number: usize,
    mz: usize,
    z: usize,
    mass: usize,
    accession: usize,
    organism: usize,
    protein_name: usize,
    protein_location: usize,
    base_sequence: usize,
    flanking_residues: usize,
    peptide: usize,
    mod_number: usize,
    theoretical_mass: usize,
    score: usize,
    rank: usize,
    matched_ion_series: usize,
    matched_ion_mz_ratios: usize,
    matched_ion_mass_error: usize,
    matched_ion_ppm: usize,
    matched_ion_intensities: usize,
    matched_ion_counts: usize,
    kind: usize,
    q_value: usize,
    pep: usize,
    pep_q_value: usize,
    localisation_score: usize,
    yion_score: usize,
    diagnostic_ion_score: usize,
    plausible_glycan_number: usize,
    total_glycosylation_sites: usize,
    glycan_mass: usize,
    plausible_glycan_composition: usize,
    n_glycan_motif: usize,
    r138_144: usize,
    plausible_glycan_structure: usize,
    glycan_localisation_level: usize,
    glycan_peptide_site_specificity: usize,
    glycan_protein_site_specificity: usize,
    all_potential_glycan_localisations: usize,
    all_site_specific_localisation_probabilities: usize,
    version: OpairVersion,
}

impl OpairFormat {
    const NUMBER_OF_COLUMNS: usize = 43;
}

/// All possible peaks versions
#[derive(Clone, Debug)]
pub enum OpairVersion {
    /// The single known version
    Opair,
}

impl std::fmt::Display for OpairVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::Opair => "Opair",
            }
        )
    }
}

/// The only supported format for Opair data
pub const O_PAIR: OpairFormat = OpairFormat {
    file_name: 0,
    scan_number: 1,
    rt: 2,
    precursor_scan_number: 3,
    mz: 4,
    z: 5,
    mass: 6,
    accession: 7,
    organism: 8,
    protein_name: 9,
    protein_location: 10,
    base_sequence: 11,
    flanking_residues: 12,
    peptide: 13,
    mod_number: 14,
    theoretical_mass: 15,
    score: 16,
    rank: 17,
    matched_ion_series: 18,
    matched_ion_mz_ratios: 19,
    matched_ion_mass_error: 20,
    matched_ion_ppm: 21,
    matched_ion_intensities: 22,
    matched_ion_counts: 23,
    kind: 24,
    q_value: 25,
    pep: 26,
    pep_q_value: 27,
    localisation_score: 28,
    yion_score: 29,
    diagnostic_ion_score: 30,
    plausible_glycan_number: 31,
    total_glycosylation_sites: 32,
    glycan_mass: 33,
    plausible_glycan_composition: 34,
    n_glycan_motif: 35,
    r138_144: 36,
    plausible_glycan_structure: 37,
    glycan_localisation_level: 38,
    glycan_peptide_site_specificity: 39,
    glycan_protein_site_specificity: 40,
    all_potential_glycan_localisations: 41,
    all_site_specific_localisation_probabilities: 42,
    version: OpairVersion::Opair,
};

/// A single parsed line of an Opair file
#[allow(missing_docs)]
#[derive(Debug)]
pub struct OpairData {
    pub file_name: String,
    pub scan: usize,
    pub rt: f64,
    pub precursor_scan_number: usize,
    pub mz: f64,
    pub z: usize,
    pub mass: f64,
    pub accession: String,
    pub organism: String,
    pub protein_name: String,
    pub protein_location: (usize, usize),
    pub base_sequence: String,
    pub flanking_residues: (AminoAcid, AminoAcid),
    pub peptide: LinearPeptide,
    pub mod_number: usize,
    pub theoretical_mass: f64,
    /// [0-1]
    pub score: f64,
    pub rank: usize,
    pub matched_ion_series: String,
    pub matched_ion_mz_ratios: String,
    pub matched_ion_mass_error: String,
    pub matched_ion_ppm: String,
    pub matched_ion_intensities: String,
    pub matched_ion_counts: String,
    pub kind: OpairMatchKind,
    pub q_value: f64,
    pub pep: f64,
    pub pep_q_value: f64,
    pub localisation_score: f64,
    pub yion_score: f64,
    pub diagnostic_ion_score: f64,
    pub plausible_glycan_number: usize,
    pub total_glycosylation_sites: usize,
    pub glycan_mass: f64,
    pub plausible_glycan_composition: String,
    pub n_glycan_motif: bool,
    pub r138_144: f64,
    pub plausible_glycan_structure: String,
    pub glycan_localisation_level: String, // Without the `Level` in from
    pub glycan_peptide_site_specificity: String,
    pub glycan_protein_site_specificity: String,
    pub all_potential_glycan_localisations: String,
    pub all_site_specific_localisation_probabilities: String,
    pub version: OpairVersion,
}

#[derive(Debug)]
#[allow(missing_docs)]
pub enum OpairMatchKind {
    Decoy,
    Contamination,
    Target,
}

impl std::fmt::Display for OpairMatchKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::Decoy => "Decoy",
                Self::Contamination => "Contamination",
                Self::Target => "Target",
            }
        )
    }
}

impl IdentifiedPeptideSource for OpairData {
    type Source = CsvLine;
    type Format = OpairFormat;
    fn parse(source: &CsvLine) -> Result<(Self, &'static OpairFormat), CustomError> {
        for format in [&O_PAIR] {
            if let Ok(peptide) = Self::parse_specific(source, format) {
                return Ok((peptide, format));
            }
        }
        Err(CustomError::error(
            "Invalid Opair line",
            "The correct format could not be determined automatically",
            source.full_context(),
        ))
    }
    fn parse_specific(source: &CsvLine, format: &OpairFormat) -> Result<Self, CustomError> {
        if source.fields.len() != OpairFormat::NUMBER_OF_COLUMNS {
            return Err(CustomError::error(
                "Invalid Opair line", 
                format!("The number of columns ({}) is not equal to the expected number of columns ({})", source.fields.len(), OpairFormat::NUMBER_OF_COLUMNS), 
                source.full_context()));
        }
        let number_error = CustomError::error(
            "Invalid Opair line",
            format!("This column is not a number but it is required to be a number in this Opair format ({})", format.version),
            source.full_context(),
        );
        Ok(Self {
            file_name: Location::column(format.file_name, source).parse(&number_error)?,
            scan: Location::column(format.scan_number, source).parse(&number_error)?,
            rt: Location::column(format.rt, source).parse(&number_error)?,
            precursor_scan_number: Location::column(format.precursor_scan_number, source)
                .parse(&number_error)?,
            mz: Location::column(format.mz, source).parse(&number_error)?,
            z: Location::column(format.z, source).parse(&number_error)?,
            mass: Location::column(format.mass, source).parse(&number_error)?,
            accession: Location::column(format.accession, source).get_string(),
            organism: Location::column(format.organism, source).get_string(),
            protein_name: Location::column(format.protein_name, source).get_string(),
            protein_location: Location::column(format.protein_location, source).parse_with(
                |loc| {
                    if loc.location.len() < 3 {
                        return Err(CustomError::error(
                            "Invalid Opair line",
                            "The location is not defined, it should be defined like this [<start> to <end>]",
                            Context::line(
                                source.line_index+1,
                                source.line.clone(),
                                loc.location.start,
                                loc.location.len(),
                            ),
                        ))
                    }
                    let bytes =
                        loc.line.line[loc.location.start + 1..loc.location.end-1].as_bytes();
                    let start = bytes.iter().take_while(|c| c.is_ascii_digit()).count();
                    let end = bytes
                        .iter()
                        .rev()
                        .take_while(|c| c.is_ascii_digit())
                        .count();
                    Ok((
                        loc.line.line[loc.location.start + 1..loc.location.start + 1 + start]
                            .parse()
                            .map_err(|_| {
                                number_error.with_context(Context::line(
                                    source.line_index+1,
                                    source.line.clone(),
                                    loc.location.start + 1,
                                    start,
                                ))
                            })?,
                        loc.line.line[loc.location.end - 1 - end..loc.location.end - 1]
                            .parse()
                            .map_err(|_| {
                                number_error.with_context(Context::line(
                                    source.line_index+1,
                                    source.line.clone(),
                                    loc.location.end - 1 - end,
                                    end,
                                ))
                            })?,
                    ))
                },
            )?,
            base_sequence: Location::column(format.base_sequence, source).get_string(),
            flanking_residues: Location::column(format.flanking_residues, source).parse_with(
                |loc| {
                    Ok((
                        AminoAcid::try_from(loc.line.line.as_bytes()[loc.location.start]).map_err(
                            |()| {
                                CustomError::error(
                                    "Invalid Opair line",
                                    "The flanking residues could not be parsed as amino acids",
                                    Context::line(
                                        source.line_index+1,
                                        source.line.clone(),
                                        loc.location.start,
                                        1,
                                    ),
                                )
                            },
                        )?,
                        AminoAcid::try_from(loc.line.line.as_bytes()[loc.location.end - 1])
                            .map_err(|()| {
                                CustomError::error(
                                    "Invalid Opair line",
                                    "The flanking residues could not be parsed as amino acids",
                                    Context::line(
                                        source.line_index+1,
                                        source.line.clone(),
                                        loc.location.end - 1,
                                        1,
                                    ),
                                )
                            })?,
                    ))
                },
            )?,
            peptide: ComplexPeptide::sloppy_pro_forma(
                &source.line,
                source.fields[format.peptide].clone(),
            )?,
            mod_number: Location::column(format.mod_number, source).parse(&number_error)?,
            theoretical_mass: Location::column(format.theoretical_mass, source)
                .parse(&number_error)?,
            score: Location::column(format.score, source).parse::<f64>(&number_error)? / 100.0,
            rank: Location::column(format.rank, source).parse(&number_error)?,
            matched_ion_series: Location::column(format.matched_ion_series, source).get_string(),
            matched_ion_mz_ratios: Location::column(format.matched_ion_mz_ratios, source)
                .get_string(),
            matched_ion_mass_error: Location::column(format.matched_ion_mass_error, source)
                .get_string(),
            matched_ion_ppm: Location::column(format.matched_ion_ppm, source).get_string(),
            matched_ion_intensities: Location::column(format.matched_ion_intensities, source)
                .get_string(),
            matched_ion_counts: Location::column(format.matched_ion_counts, source).get_string(),
            kind: Location::column(format.kind, source).parse_with(|loc| {
                match &loc.line.line[loc.location.clone()] {
                    "T" => Ok(OpairMatchKind::Target),
                    "C" => Ok(OpairMatchKind::Contamination),
                    "D" => Ok(OpairMatchKind::Decoy),
                    _ => Err(CustomError::error(
                        "Invalid Opair line",
                        "The kind column does not contain a valid value (T/C/D)",
                        Context::line(
                            source.line_index+1,
                            source.line.clone(),
                            loc.location.start,
                            loc.location.len(),
                        ),
                    )),
                }
            })?,
            q_value: Location::column(format.q_value, source).parse(&number_error)?,
            pep: Location::column(format.pep, source).parse(&number_error)?,
            pep_q_value: Location::column(format.pep_q_value, source).parse(&number_error)?,
            localisation_score: Location::column(format.localisation_score, source)
                .parse(&number_error)?,
            yion_score: Location::column(format.yion_score, source).parse(&number_error)?,
            diagnostic_ion_score: Location::column(format.diagnostic_ion_score, source)
                .parse(&number_error)?,
            plausible_glycan_number: Location::column(format.plausible_glycan_number, source)
                .parse(&number_error)?,
            total_glycosylation_sites: Location::column(format.total_glycosylation_sites, source)
                .parse(&number_error)?,
            glycan_mass: Location::column(format.glycan_mass, source).parse(&number_error)?,
            plausible_glycan_composition: Location::column(
                format.plausible_glycan_composition,
                source,
            )
            .get_string(),
            n_glycan_motif: Location::column(format.n_glycan_motif, source).parse_with(|loc| {
                match &loc.line.line[loc.location.clone()] {
                    "TRUE" => Ok(true),
                    "FALSE" => Ok(false),
                    _ => Err(CustomError::error(
                        "Invalid Opair line",
                        "The N glycan motif check column does not contain a valid value (TRUE/FALSE)",
                        Context::line(
                            source.line_index+1,
                            source.line.clone(),
                            loc.location.start,
                            loc.location.len(),
                        ),
                    )),
                }
            })?,
            r138_144: Location::column(format.r138_144, source).parse(&number_error)?,
            plausible_glycan_structure: Location::column(format.plausible_glycan_structure, source)
                .get_string(),
            glycan_localisation_level: Location::column(format.glycan_localisation_level, source)
                .get_string()
                .trim_start_matches("Level")
                .to_string(),
            glycan_peptide_site_specificity: Location::column(
                format.glycan_peptide_site_specificity,
                source,
            )
            .get_string(),
            glycan_protein_site_specificity: Location::column(
                format.glycan_protein_site_specificity,
                source,
            )
            .get_string(),
            all_potential_glycan_localisations: Location::column(
                format.all_potential_glycan_localisations,
                source,
            )
            .get_string(),
            all_site_specific_localisation_probabilities: Location::column(
                format.all_site_specific_localisation_probabilities,
                source,
            )
            .get_string(),
            version: format.version.clone(),
        })
    }
    fn parse_file(
        path: impl AsRef<std::path::Path>,
    ) -> Result<BoxedIdentifiedPeptideIter<Self>, String> {
        parse_csv(path, b'\t').map(|lines| {
            Self::parse_many::<Box<dyn Iterator<Item = Self::Source>>>(Box::new(
                lines.skip(1).map(Result::unwrap),
            ))
        })
    }
}

impl From<OpairData> for IdentifiedPeptide {
    fn from(value: OpairData) -> Self {
        Self {
            peptide: value.peptide.clone(),
            local_confidence: None,
            score: Some(value.score),
            metadata: MetaData::Opair(value),
        }
    }
}
