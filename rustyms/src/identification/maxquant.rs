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
static BOOL_ERROR: (&str, &str) = (
    "Invalid MaxQuant line",
    "This column is not a boolean ('0' or '1') but it is required to be a boolean in this MaxQuant format",
);

format_family!(
    /// The format for any MaxQuant file
    MaxQuantFormat,
    /// The data from any MaxQuant file
    MaxQuantData,
    MaxQuantVersion, [&MSMS, &MSMS_SCANS, &NOVO_MSMS_SCANS, &SILAC], b'\t';
    required {
        raw_file: String, |location: Location, _| Ok(location.get_string());
        scan_number: Vec<usize>, |location: Location, _| location.or_empty().array(';').map(|s| s.parse(NUMBER_ERROR)).collect::<Result<Vec<usize>, CustomError>>();
        modifications: String, |location: Location, _| Ok(location.get_string());
        proteins: String, |location: Location, _| Ok(location.get_string());
        peptide: Option<LinearPeptide<VerySimple>>, |location: Location, custom_database: Option<&CustomDatabase>| location.or_empty().parse_with(|location| LinearPeptide::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            SloppyParsingParameters::default()
        ));
        z: Charge, |location: Location, _| location.parse::<usize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        ty: String, |location: Location, _| Ok(location.get_string());
        pep: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
    }
    optional {
        all_modified_sequences: Vec<LinearPeptide<VerySimple>>, |location: Location, custom_database: Option<&CustomDatabase>| location.array(';')
                .map(|s| LinearPeptide::sloppy_pro_forma(s.line.line(), s.location, custom_database, SloppyParsingParameters::default()))
                .collect::<Result<Vec<LinearPeptide<VerySimple>>, CustomError>>();
        base_peak_intensity: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        carbamidomethyl_c_probabilities: String, |location: Location, _| Ok(location.get_string());
        carbamidomethyl_c_score_differences: String, |location: Location, _| Ok(location.get_string());
        collision_energy: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        delta_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        dn_c_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        dn_combined_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        dn_missing_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        dn_n_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        dn_sequence: String, |location: Location, _| Ok(location.get_string());
        evidence_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        experiment: String, |location: Location, _| Ok(location.get_string());
        fragmentation: String, |location: Location, _| Ok(location.get_string());
        genes: String, |location: Location, _| Ok(location.get_string());
        id: usize,|location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        intensity_coverage: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        intensity_h: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        intensity_l: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        intensity: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        isotope_index: isize, |location: Location, _| location.or_empty().parse::<isize>(NUMBER_ERROR);
        labeling_state: bool, |location: Location, _| location.or_empty().ignore("-1").parse::<u8>(BOOL_ERROR).map(|n| n.map(|n| n != 0));
        localisation_probability: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        mass_analyser: String, |location: Location, _| Ok(location.get_string());
        mass_error_da: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        mass_error_ppm: Ratio, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::ppm>);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        missed_cleavages: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        modified_peptide_id:usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        nem_probabilities: String, |location: Location, _| Ok(location.get_string());
        nem_score_differences: String, |location: Location, _| Ok(location.get_string());
        number_of_matches: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        oxidation_m_probabilities: String, |location: Location, _| Ok(location.get_string());
        oxidation_m_score_differences: String, |location: Location, _| Ok(location.get_string());
        peak_coverage: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        peptide_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        precursor: usize, |location: Location, _| location.ignore("-1").parse::<usize>(NUMBER_ERROR);
        precursor_intensity: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        precursor_apex_function: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        precursor_apex_offset: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        precursor_apex_offset_time: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        protein_group_ids: Vec<usize>, |location: Location, _| location.array(';').map(|p| p.parse::<usize>(NUMBER_ERROR)).collect::<Result<Vec<_>,_>>();
        ration_h_l_normalised: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        ration_h_l: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        retention_time: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        scan_event_number: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        scan_index: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        score_diff: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        simple_mass_error_ppm: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        total_ion_current: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
    }
);

impl From<MaxQuantData> for IdentifiedPeptide {
    fn from(value: MaxQuantData) -> Self {
        Self {
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
    /// MaxNovo SILAC evidence.txt
    Silac,
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
                Self::Silac => "SILAC evidence",
            }
        )
    }
}

/// msms.txt
pub const MSMS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::MSMS,
    all_modified_sequences: Some("all modified sequences"),
    base_peak_intensity: None,
    carbamidomethyl_c_probabilities: None,
    carbamidomethyl_c_score_differences: None,
    collision_energy: None,
    delta_score: Some("delta score"),
    dn_c_mass: None,
    dn_combined_score: None,
    dn_missing_mass: None,
    dn_n_mass: None,
    dn_sequence: None,
    evidence_id: Some("evidence id"),
    experiment: None,
    fragmentation: Some("fragmentation"),
    genes: None,
    id: Some("id"),
    intensity_coverage: Some("intensity coverage"),
    intensity_h: None,
    intensity_l: None,
    intensity: None,
    isotope_index: Some("isotope index"),
    labeling_state: None,
    localisation_probability: Some("localization prob"),
    mass_analyser: Some("mass analyzer"),
    mass_error_da: Some("mass error [da]"),
    mass_error_ppm: Some("mass error [ppm]"),
    mass: Some("mass"),
    missed_cleavages: Some("missed cleavages"),
    modifications: "modifications",
    modified_peptide_id: Some("mod. peptide id"),
    mz: Some("m/z"),
    nem_probabilities: None,
    nem_score_differences: None,
    number_of_matches: Some("number of matches"),
    oxidation_m_probabilities: None,
    oxidation_m_score_differences: None,
    peak_coverage: Some("peak coverage"),
    pep: "pep",
    peptide_id: Some("peptide id"),
    peptide: "modified sequence",
    precursor_apex_function: Some("precursor apex fraction"),
    precursor_apex_offset_time: Some("precursor apex offset time"),
    precursor_apex_offset: Some("precursor apex offset"),
    precursor_intensity: Some("precursor intensity"),
    precursor: Some("precursor full scan number"),
    protein_group_ids: Some("protein group ids"),
    proteins: "proteins",
    ration_h_l_normalised: None,
    ration_h_l: None,
    raw_file: "raw file",
    retention_time: Some("retention time"),
    scan_event_number: Some("scan event number"),
    scan_index: Some("scan index"),
    scan_number: "scan number",
    score_diff: Some("score diff"),
    score: "score",
    simple_mass_error_ppm: Some("simple mass error [ppm]"),
    total_ion_current: None,
    ty: "type",
    z: "charge",
};

/// msmsScans.txt
pub const MSMS_SCANS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::MSMSScans,
    all_modified_sequences: None,
    base_peak_intensity: Some("base peak intensity"),
    carbamidomethyl_c_probabilities: None,
    carbamidomethyl_c_score_differences: None,
    collision_energy: Some("collision energy"),
    delta_score: None,
    dn_c_mass: None,
    dn_combined_score: None,
    dn_missing_mass: None,
    dn_n_mass: None,
    dn_sequence: None,
    evidence_id: None,
    experiment: None,
    fragmentation: Some("fragmentation"),
    genes: None,
    id: None,
    intensity_coverage: None,
    intensity_h: None,
    intensity_l: None,
    intensity: None,
    isotope_index: None,
    labeling_state: None,
    localisation_probability: None,
    mass_analyser: Some("mass analyzer"),
    mass_error_da: None,
    mass_error_ppm: None,
    mass: Some("mass"),
    missed_cleavages: None,
    modifications: "modifications",
    modified_peptide_id: None,
    mz: Some("m/z"),
    nem_probabilities: None,
    nem_score_differences: None,
    number_of_matches: None,
    oxidation_m_probabilities: None,
    oxidation_m_score_differences: None,
    peak_coverage: None,
    pep: "pep",
    peptide_id: None,
    peptide: "modified sequence",
    precursor_apex_function: Some("precursor apex fraction"),
    precursor_apex_offset_time: Some("precursor apex offset time"),
    precursor_apex_offset: Some("precursor apex offset"),
    precursor_intensity: Some("precursor intensity"),
    precursor: Some("precursor full scan number"),
    protein_group_ids: None,
    proteins: "proteins",
    ration_h_l_normalised: None,
    ration_h_l: None,
    raw_file: "raw file",
    retention_time: Some("retention time"),
    scan_event_number: Some("scan event number"),
    scan_index: Some("scan index"),
    scan_number: "scan number",
    score_diff: None,
    score: "score",
    simple_mass_error_ppm: None,
    total_ion_current: Some("total ion current"),
    ty: "type",
    z: "charge",
};

/// MaxNovo msmsScans.txt
pub const NOVO_MSMS_SCANS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::NovoMSMSScans,
    all_modified_sequences: None,
    base_peak_intensity: Some("base peak intensity"),
    carbamidomethyl_c_probabilities: None,
    carbamidomethyl_c_score_differences: None,
    collision_energy: Some("collision energy"),
    delta_score: None,
    dn_c_mass: Some("dn cterm mass"),
    dn_combined_score: Some("dn combined score"),
    dn_missing_mass: Some("dn missing mass"),
    dn_n_mass: Some("dn nterm mass"),
    dn_sequence: Some("dn sequence"),
    evidence_id: None,
    experiment: None,
    fragmentation: Some("fragmentation"),
    genes: None,
    id: None,
    intensity_coverage: None,
    intensity_h: None,
    intensity_l: None,
    intensity: None,
    isotope_index: None,
    labeling_state: None,
    localisation_probability: None,
    mass_analyser: Some("mass analyzer"),
    mass_error_da: None,
    mass_error_ppm: None,
    mass: Some("mass"),
    missed_cleavages: None,
    modifications: "modifications",
    modified_peptide_id: None,
    mz: Some("m/z"),
    nem_probabilities: None,
    nem_score_differences: None,
    number_of_matches: None,
    oxidation_m_probabilities: None,
    oxidation_m_score_differences: None,
    peak_coverage: None,
    pep: "pep",
    peptide_id: None,
    peptide: "modified sequence",
    precursor_apex_function: Some("precursor apex fraction"),
    precursor_apex_offset_time: Some("precursor apex offset time"),
    precursor_apex_offset: Some("precursor apex offset"),
    precursor_intensity: Some("precursor intensity"),
    precursor: Some("precursor full scan number"),
    protein_group_ids: None,
    proteins: "proteins",
    ration_h_l_normalised: None,
    ration_h_l: None,
    raw_file: "raw file",
    retention_time: Some("retention time"),
    scan_event_number: Some("scan event number"),
    scan_index: Some("scan index"),
    scan_number: "scan number",
    score_diff: None,
    score: "score",
    simple_mass_error_ppm: None,
    total_ion_current: Some("total ion current"),
    ty: "type",
    z: "charge",
};

/// MaxQuant v2.4.14.0 SILAC evidence.txt
pub const SILAC: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::Silac,
    all_modified_sequences: None,
    base_peak_intensity: None,
    carbamidomethyl_c_probabilities: Some("carbamidomethyl (c) probabilities"),
    carbamidomethyl_c_score_differences: Some("carbamidomethyl (c) score diffs"),
    collision_energy: None,
    delta_score: Some("delta score"),
    dn_c_mass: None,
    dn_combined_score: None,
    dn_missing_mass: None,
    dn_n_mass: None,
    dn_sequence: None,
    evidence_id: None,
    experiment: Some("experiment"),
    fragmentation: None,
    genes: Some("gene names"),
    id: Some("id"),
    intensity: Some("intensity"),
    intensity_coverage: None,
    intensity_h: Some("intensity h"),
    intensity_l: Some("intensity l"),
    isotope_index: None,
    labeling_state: Some("labeling state"),
    localisation_probability: None,
    mass_analyser: None,
    mass_error_da: Some("mass error [da]"),
    mass_error_ppm: Some("mass error [ppm]"),
    mass: Some("mass"),
    missed_cleavages: None,
    modifications: "modifications",
    modified_peptide_id: Some("mod. peptide id"),
    mz: Some("m/z"),
    nem_probabilities: Some("nem probabilities"),
    nem_score_differences: Some("nem score diffs"),
    number_of_matches: None,
    oxidation_m_probabilities: Some("oxidation (m) probabilities"),
    oxidation_m_score_differences: Some("oxidation (m) score diffs"),
    peak_coverage: None,
    pep: "pep",
    peptide_id: Some("peptide id"),
    peptide: "modified sequence",
    precursor_apex_function: None,
    precursor_apex_offset_time: None,
    precursor_apex_offset: None,
    precursor_intensity: None,
    precursor: None,
    protein_group_ids: Some("protein group ids"),
    proteins: "proteins",
    ration_h_l: Some("ratio h/l"),
    ration_h_l_normalised: Some("ratio h/l normalized"),
    raw_file: "raw file",
    retention_time: Some("retention time"),
    scan_event_number: None,
    scan_index: None,
    scan_number: "ms/ms scan numbers",
    score_diff: None,
    score: "score",
    simple_mass_error_ppm: None,
    total_ion_current: None,
    ty: "type",
    z: "charge",
};
