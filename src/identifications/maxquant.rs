use crate::{
    error::{Context, CustomError},
    helper_functions::InvertResult,
    system::{Charge, Mass, MassOverCharge, Time},
    ComplexPeptide, LinearPeptide,
};

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
        raw_file: String, |location: Location| Ok(location.get_string());
        scan_number: usize, |location: Location| location.parse(NUMBER_ERROR);
        scan_index: usize, |location: Location| location.parse(NUMBER_ERROR);
        modifications: String, |location: Location| Ok(location.get_string());
        proteins: String, |location: Location| Ok(location.get_string());
        sequence: LinearPeptide, |location: Location| ComplexPeptide::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
        );
        z: Charge, |location: Location| location.parse::<usize>(NUMBER_ERROR).map(|c| Charge::new::<crate::system::e>(c as f64));
        fragmentation: String, |location: Location| Ok(location.get_string());
        mass_analyzer: String, |location: Location| Ok(location.get_string());
        ty: MaxQuantType, |location: Location| location.parse_with(|l| match l.as_str() {
            "PEAK" => Ok(MaxQuantType::Peak),
            "MULTI" => Ok(MaxQuantType::Multi),
            _ => Err(CustomError::error(
                "Invalid MaxQuant line",
                "A MaxQuant type has to be PEAK or MULTI.",
                Context::line(
                    l.line.line_index(),
                    l.line.line(),
                    l.location.start,
                    l.location.len(),
                ),
            )),
        });
        scan_event_number: usize, |location: Location| location.parse(NUMBER_ERROR);
        pep: f64, |location: Location| location.parse(NUMBER_ERROR);
        score: f64, |location: Location| location.parse(NUMBER_ERROR);
        precursor: Option<usize>, |location: Location| location.ignore("-1").parse::<usize>(NUMBER_ERROR);
        precursor_intensity: f64, |location: Location| location.parse(NUMBER_ERROR);
        precursor_apex_function: f64, |location: Location| location.parse(NUMBER_ERROR);
        precursor_apex_offset: f64, |location: Location| location.parse(NUMBER_ERROR);
        precursor_apex_offset_time: f64, |location: Location| location.parse(NUMBER_ERROR);
    }
    optional {
        missed_cleavages: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        isotope_index: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        mass: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        mass_error_da: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        mass_error_ppm: f64, |location: Location| location.parse::<f64>(NUMBER_ERROR);
        simple_mass_error_ppm: f64, |location: Location| location.parse::<f64>(NUMBER_ERROR);
        retention_time:Time, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        number_of_matches: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        intensity_coverage: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        peak_coverage:usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        delta_score: f64, |location: Location| location.parse::<f64>(NUMBER_ERROR);
        score_diff: f64, |location: Location| location.parse::<f64>(NUMBER_ERROR);
        localisation_probability: f64, |location: Location| location.parse::<f64>(NUMBER_ERROR);
        all_modified_sequences: Vec<LinearPeptide>,|location: Location| location.array(';')
                .map(|s| ComplexPeptide::sloppy_pro_forma(&s.line.line(), s.location))
                .collect::<Result<Vec<LinearPeptide>, CustomError>>();
        id: usize,|location: Location| location.parse::<usize>(NUMBER_ERROR);
        protein_group_ids: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        peptide_id: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        modified_peptide_id:usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        evidence_id: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        base_peak_intensity: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        total_ion_current: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        collision_energy: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        dn_sequence: String, |location: Location| Ok(location.get_string());
        dn_combined_score: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        dn_n_mass: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        dn_c_mass: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        dn_missing_mass: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
    }
);

impl From<MaxQuantData> for IdentifiedPeptide {
    fn from(value: MaxQuantData) -> Self {
        Self {
            peptide: value.sequence.clone(),
            local_confidence: None,
            score: Some(value.score),
            metadata: MetaData::MaxQuant(value),
        }
    }
}

/// All possible MaxQuant versions
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MaxQuantVersion {
    /// msms.txt
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
    sequence: "modified sequence",
    z: "charge",
    fragmentation: "fragmentation",
    mass_analyzer: "mass analyzer",
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
    sequence: "modified sequence",
    z: "charge",
    fragmentation: "fragmentation",
    mass_analyzer: "mass analyzer",
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
    sequence: "modified sequence",
    z: "charge",
    fragmentation: "fragmentation",
    mass_analyzer: "mass analyzer",
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
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MaxQuantType {
    Peak,
    Multi,
}
