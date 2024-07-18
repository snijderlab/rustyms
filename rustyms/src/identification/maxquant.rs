use crate::{
    error::CustomError,
    helper_functions::InvertResult,
    ontologies::CustomDatabase,
    peptide::{SloppyParsingParameters, VerySimple},
    system::{usize::Charge, Mass, MassOverCharge, Ratio, Time},
    LinearPeptide,
};
use serde::{Deserialize, Serialize};

use super::{
    common_parser::{Location, OptionalLocation},
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid MaxQuant line",
    "This column is not a number but it is required to be a number in this MaxQuant format",
);

format_family!(
    /// The format for any MaxQuant file
    MaxQuantFormat,
    /// The data from any MaxQuant file
    MaxQuantData,
    MaxQuantVersion, [&MSMS, &MSMS_SCANS, &NOVO_MSMS_SCANS], b'\t';
    required {
        raw_file: String, |location: Location, _| Ok(location.get_string());
        scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        scan_index: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        modifications: String, |location: Location, _| Ok(location.get_string());
        proteins: String, |location: Location, _| Ok(location.get_string());
        peptide: Option<LinearPeptide<VerySimple>>, |location: Location, custom_database: Option<&CustomDatabase>| location.or_empty().parse_with(|location| LinearPeptide::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            SloppyParsingParameters::default()
        ));
        z: Charge, |location: Location, _| location.parse::<usize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        fragmentation: String, |location: Location, _| Ok(location.get_string());
        mass_analyser: String, |location: Location, _| Ok(location.get_string());
        ty: String, |location: Location, _| Ok(location.get_string());
        scan_event_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        pep: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor: Option<usize>, |location: Location, _| location.ignore("-1").parse::<usize>(NUMBER_ERROR);
        precursor_intensity: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_apex_function: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_apex_offset: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_apex_offset_time: f64, |location: Location, _| location.parse(NUMBER_ERROR);
    }
    optional {
        missed_cleavages: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        isotope_index: isize, |location: Location, _| location.or_empty().parse::<isize>(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        mass_error_da: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        mass_error_ppm: Ratio, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::ppm>);
        simple_mass_error_ppm: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        retention_time:Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        number_of_matches: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        intensity_coverage: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        peak_coverage: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        delta_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        score_diff: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        localisation_probability: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        all_modified_sequences: Vec<LinearPeptide<VerySimple>>,|location: Location, custom_database: Option<&CustomDatabase>| location.array(';')
                .map(|s| LinearPeptide::sloppy_pro_forma(s.line.line(), s.location, custom_database, SloppyParsingParameters::default()))
                .collect::<Result<Vec<LinearPeptide<VerySimple>>, CustomError>>();
        id: usize,|location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        peptide_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        protein_group_ids: Vec<usize>, |location: Location, _| location.array(';').map(|p| p.parse::<usize>(NUMBER_ERROR)).collect::<Result<Vec<_>,_>>();
        modified_peptide_id:usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        evidence_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        base_peak_intensity: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        total_ion_current: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        collision_energy: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        dn_sequence: String, |location: Location, _| Ok(location.get_string());
        dn_combined_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        dn_n_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        dn_c_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        dn_missing_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
    }
);

impl From<MaxQuantData> for IdentifiedPeptide {
    fn from(value: MaxQuantData) -> Self {
        Self {
            peptide: value.peptide.clone().unwrap_or_default(), // TODO: what to do with empty sequences
            local_confidence: None,
            score: Some(value.score),
            metadata: MetaData::MaxQuant(value),
        }
    }
}

/// All possible MaxQuant versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
pub enum MaxQuantVersion {
    /// msms.txt
    #[default]
    #[allow(clippy::upper_case_acronyms)]
    MSMS,
    /// msmsScans.txt
    MSMSScans,
    /// MaxNovo msmsScans.txt
    NovoMSMSScans,
}

impl std::fmt::Display for MaxQuantVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::MSMS => "msms",
                Self::MSMSScans => "msmsScans",
                Self::NovoMSMSScans => "de novo msmsScans",
            }
        )
    }
}

/// msms.txt
pub const MSMS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::MSMS,
    raw_file: "raw file",
    scan_number: "scan number",
    scan_index: "scan index",
    modifications: "modifications",
    proteins: "proteins",
    peptide: "modified sequence",
    z: "charge",
    fragmentation: "fragmentation",
    mass_analyser: "mass analyzer",
    ty: "type",
    scan_event_number: "scan event number",
    pep: "pep",
    score: "score",
    precursor: "precursor full scan number",
    precursor_intensity: "precursor intensity",
    precursor_apex_function: "precursor apex fraction",
    precursor_apex_offset: "precursor apex offset",
    precursor_apex_offset_time: "precursor apex offset time",
    missed_cleavages: Some("missed cleavages"),
    isotope_index: Some("isotope index"),
    mz: Some("m/z"),
    mass: Some("mass"),
    mass_error_da: Some("mass error [da]"),
    mass_error_ppm: Some("mass error [ppm]"),
    simple_mass_error_ppm: Some("simple mass error [ppm]"),
    retention_time: Some("retention time"),
    number_of_matches: Some("number of matches"),
    intensity_coverage: Some("intensity coverage"),
    peak_coverage: Some("peak coverage"),
    delta_score: Some("delta score"),
    score_diff: Some("score diff"),
    localisation_probability: Some("localization prob"),
    all_modified_sequences: Some("all modified sequences"),
    id: Some("id"),
    protein_group_ids: Some("protein group ids"),
    peptide_id: Some("peptide id"),
    modified_peptide_id: Some("mod. peptide id"),
    evidence_id: Some("evidence id"),
    base_peak_intensity: None,
    total_ion_current: None,
    collision_energy: None,
    dn_sequence: None,
    dn_combined_score: None,
    dn_n_mass: None,
    dn_c_mass: None,
    dn_missing_mass: None,
};

/// msmsScans.txt
pub const MSMS_SCANS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::MSMSScans,
    raw_file: "raw file",
    scan_number: "scan number",
    scan_index: "scan index",
    modifications: "modifications",
    proteins: "proteins",
    peptide: "modified sequence",
    z: "charge",
    fragmentation: "fragmentation",
    mass_analyser: "mass analyzer",
    ty: "type",
    scan_event_number: "scan event number",
    pep: "pep",
    score: "score",
    precursor: "precursor full scan number",
    precursor_intensity: "precursor intensity",
    precursor_apex_function: "precursor apex fraction",
    precursor_apex_offset: "precursor apex offset",
    precursor_apex_offset_time: "precursor apex offset time",
    missed_cleavages: None,
    isotope_index: None,
    mz: Some("m/z"),
    mass: Some("mass"),
    mass_error_da: None,
    mass_error_ppm: None,
    simple_mass_error_ppm: None,
    retention_time: Some("retention time"),
    number_of_matches: None,
    intensity_coverage: None,
    peak_coverage: None,
    delta_score: None,
    score_diff: None,
    localisation_probability: None,
    all_modified_sequences: None,
    id: None,
    protein_group_ids: None,
    peptide_id: None,
    modified_peptide_id: None,
    evidence_id: None,
    base_peak_intensity: Some("base peak intensity"),
    total_ion_current: Some("total ion current"),
    collision_energy: Some("collision energy"),
    dn_sequence: None,
    dn_combined_score: None,
    dn_n_mass: None,
    dn_c_mass: None,
    dn_missing_mass: None,
};

/// MaxNovo msmsScans.txt
pub const NOVO_MSMS_SCANS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::NovoMSMSScans,
    raw_file: "raw file",
    scan_number: "scan number",
    scan_index: "scan index",
    modifications: "modifications",
    proteins: "proteins",
    peptide: "modified sequence",
    z: "charge",
    fragmentation: "fragmentation",
    mass_analyser: "mass analyzer",
    ty: "type",
    scan_event_number: "scan event number",
    pep: "pep",
    score: "score",
    precursor: "precursor full scan number",
    precursor_intensity: "precursor intensity",
    precursor_apex_function: "precursor apex fraction",
    precursor_apex_offset: "precursor apex offset",
    precursor_apex_offset_time: "precursor apex offset time",
    missed_cleavages: None,
    isotope_index: None,
    mz: Some("m/z"),
    mass: Some("mass"),
    mass_error_da: None,
    mass_error_ppm: None,
    simple_mass_error_ppm: None,
    retention_time: Some("retention time"),
    number_of_matches: None,
    intensity_coverage: None,
    peak_coverage: None,
    delta_score: None,
    score_diff: None,
    localisation_probability: None,
    all_modified_sequences: None,
    id: None,
    protein_group_ids: None,
    peptide_id: None,
    modified_peptide_id: None,
    evidence_id: None,
    base_peak_intensity: Some("base peak intensity"),
    total_ion_current: Some("total ion current"),
    collision_energy: Some("collision energy"),
    dn_sequence: Some("dn sequence"),
    dn_combined_score: Some("dn combined score"),
    dn_n_mass: Some("dn nterm mass"),
    dn_c_mass: Some("dn cterm mass"),
    dn_missing_mass: Some("dn missing mass"),
};
