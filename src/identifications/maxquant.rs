use crate::{
    error::{Context, CustomError},
    helper_functions::InvertResult,
    system::{Mass, MassOverCharge, Time},
    ComplexPeptide, LinearPeptide,
};

use super::{
    common_parser::{Location, OptionalLocation},
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
};

// DB:   Raw file	Scan number	Scan index	Sequence	Length	Missed cleavages	Modifications	Modified sequence	Oxidation (M) Probabilities	oxidation(w) Probabilities	Oxidation (M) Score diffs	oxidation(w) Score diffs	Gln->pyro-Glu	Glu->pyro-Glu	Oxidation (M)	oxidation(w)	Proteins	Charge	Fragmentation	Mass analyzer	Type	Scan event number	Isotope index	m/z	Mass	Mass error [ppm]	Mass error [Da]	Simple mass error [ppm]	Retention time	PEP	Score	Delta score	Score diff	Localization prob	Combinatorics	PIF	Fraction of total spectrum	Base peak fraction	Precursor full scan number	Precursor Intensity	Precursor apex fraction	Precursor apex offset	Precursor apex offset time	Matches	Intensities	Mass deviations [Da]	Mass deviations [ppm]	Masses	Number of matches	Intensity coverage	Peak coverage	Unfragmented precursor intensity	Unfragmented precursor fraction	Neutral loss level	ETD identification type	Reverse	All scores	All sequences	All modified sequences	MS3 scan numbers	Reporter PIF	Reporter fraction	id	Protein group IDs	Peptide ID	Mod. peptide ID	Evidence ID	Oxidation (M) site IDs	oxidation(w) site IDs	Mass deficit
// NOVO: Raw file	Scan number	Scan index	Sequence	Length	Missed cleavages	Modifications	Modified sequence	Oxidation (M) Probabilities	oxidation(w) Probabilities	Oxidation (M) Score diffs	oxidation(w) Score diffs	Gln->pyro-Glu	Glu->pyro-Glu	Oxidation (M)	oxidation(w)	Proteins	Charge	Fragmentation	Mass analyzer	Type	Scan event number	Isotope index	m/z	Mass	Mass error [ppm]	Mass error [Da]	Simple mass error [ppm]	Retention time	PEP	Score	Delta score	Score diff	Localization prob	Combinatorics	PIF	Fraction of total spectrum	Base peak fraction	Precursor full scan number	Precursor Intensity	Precursor apex fraction	Precursor apex offset	Precursor apex offset time	Matches	Intensities	Mass deviations [Da]	Mass deviations [ppm]	Masses	Number of matches	Intensity coverage	Peak coverage	Unfragmented precursor intensity	Unfragmented precursor fraction	Neutral loss level	ETD identification type	Reverse	All scores	All sequences	All modified sequences	MS3 scan numbers	Reporter PIF	Reporter fraction	id	Protein group IDs	Peptide ID	Mod. peptide ID	Evidence ID	Oxidation (M) site IDs	oxidation(w) site IDs	Mass deficit
//

/// The file format for any max quant format, determining the existence and location of all possible columns
#[derive(Debug, Clone)]
pub struct MaxQuantFormat {
    raw_file: usize,
    scan_number: usize,
    scan_index: usize,
    missed_cleavages: Option<usize>,
    modifications: usize,
    modified_sequence: usize,
    proteins: usize,
    charge: usize,
    fragmentation: usize,
    mass_analyzer: usize,
    ty: usize,
    scan_event_number: usize,
    isotope_index: Option<usize>,
    mz: Option<usize>,
    mass: Option<usize>,
    mass_error_da: Option<usize>,
    mass_error_ppm: Option<usize>,
    simple_mass_error_ppm: Option<usize>,
    retention_time: Option<usize>,
    pep: usize,
    score: usize,
    delta_score: Option<usize>,
    score_diff: Option<usize>,
    localisation_probability: Option<usize>,
    precursor_full_scan_number: usize,
    precursor_intensity: usize,
    precursor_apex_function: usize,
    precursor_apex_offset: usize,
    precursor_apex_offset_time: usize,
    number_of_matches: Option<usize>,
    intensity_coverage: Option<usize>,
    peak_coverage: Option<usize>,
    all_modified_sequences: Option<usize>,
    id: Option<usize>,
    protein_group_ids: Option<usize>,
    peptide_id: Option<usize>,
    modified_peptide_id: Option<usize>,
    evidence_id: Option<usize>,
    base_peak_intensity: Option<usize>,
    total_ion_current: Option<usize>,
    collision_energy: Option<usize>,
    dn_sequence: Option<usize>,
    dn_combined_score: Option<usize>,
    dn_n_mass: Option<usize>,
    dn_c_mass: Option<usize>,
    dn_missing_mass: Option<usize>,
    version: MaxQuantVersion,
}

impl MaxQuantFormat {
    const fn number_of_columns(&self) -> usize {
        match self.version {
            MaxQuantVersion::MSMS => 70,
            MaxQuantVersion::MSMSScans => 48,
            MaxQuantVersion::NovoMSMSScans => 80,
        }
    }
}

/// All possible MaxQuant versions
#[derive(Clone, Debug)]
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
    raw_file: 0,
    scan_number: 1,
    scan_index: 2,
    missed_cleavages: Some(5),
    modifications: 6,
    modified_sequence: 7,
    proteins: 16,
    charge: 17,
    fragmentation: 18,
    mass_analyzer: 19,
    ty: 20,
    scan_event_number: 21,
    isotope_index: Some(22),
    mz: Some(23),
    mass: Some(24),
    mass_error_da: Some(26),
    mass_error_ppm: Some(25),
    simple_mass_error_ppm: Some(27),
    retention_time: Some(28),
    pep: 29,
    score: 30,
    delta_score: Some(31),
    score_diff: Some(32),
    localisation_probability: Some(33),
    precursor_full_scan_number: 38,
    precursor_intensity: 39,
    precursor_apex_function: 40,
    precursor_apex_offset: 41,
    precursor_apex_offset_time: 42,
    number_of_matches: Some(48),
    intensity_coverage: Some(49),
    peak_coverage: Some(50),
    all_modified_sequences: Some(58),
    id: Some(62),
    protein_group_ids: Some(63),
    peptide_id: Some(64),
    modified_peptide_id: Some(65),
    evidence_id: Some(66),
    base_peak_intensity: None,
    total_ion_current: None,
    collision_energy: None,
    dn_sequence: None,
    dn_combined_score: None,
    dn_n_mass: None,
    dn_c_mass: None,
    dn_missing_mass: None,
    version: MaxQuantVersion::MSMS,
};

/// msmsScans.txt
pub const MSMS_SCANS: MaxQuantFormat = MaxQuantFormat {
    raw_file: 0,
    scan_number: 1,
    scan_index: 35,
    missed_cleavages: None,
    modifications: 32,
    modified_sequence: 33,
    proteins: 34,
    charge: 18,
    fragmentation: 21,
    mass_analyzer: 22,
    ty: 20,
    scan_event_number: 31,
    isotope_index: None,
    mz: None,
    mass: None,
    mass_error_da: None,
    mass_error_ppm: None,
    simple_mass_error_ppm: None,
    retention_time: None,
    pep: 36,
    score: 35,
    delta_score: None,
    score_diff: None,
    localisation_probability: None,
    precursor_full_scan_number: 26,
    precursor_intensity: 27,
    precursor_apex_function: 28,
    precursor_apex_offset: 29,
    precursor_apex_offset_time: 30,
    number_of_matches: None,
    intensity_coverage: None,
    peak_coverage: None,
    all_modified_sequences: None,
    id: None,
    protein_group_ids: None,
    peptide_id: None,
    modified_peptide_id: None,
    evidence_id: None,
    base_peak_intensity: Some(7),
    total_ion_current: Some(4),
    collision_energy: Some(5),
    dn_sequence: None,
    dn_combined_score: None,
    dn_n_mass: None,
    dn_c_mass: None,
    dn_missing_mass: None,
    version: MaxQuantVersion::MSMSScans,
};

/// MaxNovo msmsScans.txt
pub const NOVO_MSMS_SCANS: MaxQuantFormat = MaxQuantFormat {
    raw_file: 0,
    scan_number: 1,
    scan_index: 35,
    missed_cleavages: None,
    modifications: 32,
    modified_sequence: 33,
    proteins: 34,
    charge: 18,
    fragmentation: 21,
    mass_analyzer: 22,
    ty: 20,
    scan_event_number: 31,
    isotope_index: None,
    mz: None,
    mass: None,
    mass_error_da: None,
    mass_error_ppm: None,
    simple_mass_error_ppm: None,
    retention_time: None,
    pep: 36,
    score: 35,
    delta_score: None,
    score_diff: None,
    localisation_probability: None,
    precursor_full_scan_number: 26,
    precursor_intensity: 27,
    precursor_apex_function: 28,
    precursor_apex_offset: 29,
    precursor_apex_offset_time: 30,
    number_of_matches: None,
    intensity_coverage: None,
    peak_coverage: None,
    all_modified_sequences: None,
    id: None,
    protein_group_ids: None,
    peptide_id: None,
    modified_peptide_id: None,
    evidence_id: None,
    base_peak_intensity: Some(7),
    total_ion_current: Some(4),
    collision_energy: Some(5),
    dn_sequence: Some(41),
    dn_combined_score: Some(51),
    dn_n_mass: Some(53),
    dn_c_mass: Some(54),
    dn_missing_mass: Some(55),
    version: MaxQuantVersion::NovoMSMSScans,
};
#[derive(Debug, Clone)]
pub enum MaxQuantType {
    Peak,
    Multi,
}
#[derive(Debug, Clone)]
pub struct MaxQuantPrecursor {
    pub index: usize,
    pub intensity: f64,
    pub apex_function: f64,
    pub apex_offset: f64,
    pub apex_offset_time: f64,
}

/// A single parsed line of a peaks file
#[allow(missing_docs)]
#[derive(Debug, Clone)]
pub struct MaxQuantData {
    pub raw_file: String,
    pub scan_number: usize,
    pub scan_index: usize,
    pub missed_cleavages: Option<usize>,
    pub modifications: String,
    pub proteins: String,
    pub sequence: LinearPeptide,
    pub charge: usize,
    pub fragmentation: String,
    pub mass_analyzer: String,
    pub ty: MaxQuantType,
    pub scan_event_number: usize,
    pub isotope_index: Option<usize>,
    pub mz: Option<MassOverCharge>,
    pub mass: Option<Mass>,
    pub mass_error_da: Option<Mass>,
    pub mass_error_ppm: Option<f64>,
    pub simple_mass_error_ppm: Option<f64>,
    pub retention_time: Option<Time>,
    pub pep: f64,
    pub score: f64,
    pub delta_score: Option<f64>,
    pub score_diff: Option<f64>,
    pub localisation_probability: Option<f64>,
    pub precursor: Option<MaxQuantPrecursor>,
    pub number_of_matches: Option<usize>,
    pub intensity_coverage: Option<f64>,
    pub peak_coverage: Option<f64>,
    pub all_modified_sequences: Option<Vec<LinearPeptide>>,
    pub id: Option<usize>,
    pub protein_group_ids: Option<usize>,
    pub peptide_id: Option<usize>,
    pub modified_peptide_id: Option<usize>,
    pub evidence_id: Option<usize>,
    pub base_peak_intensity: Option<usize>,
    pub total_ion_current: Option<usize>,
    pub collision_energy: Option<f64>,
    pub dn_sequence: Option<String>,
    pub dn_combined_score: Option<f64>,
    pub dn_n_mass: Option<f64>,
    pub dn_c_mass: Option<f64>,
    pub dn_missing_mass: Option<f64>,
    pub version: MaxQuantVersion,
}

impl IdentifiedPeptideSource for MaxQuantData {
    type Source = CsvLine;
    type Format = MaxQuantFormat;
    fn parse(source: &CsvLine) -> Result<(Self, &'static MaxQuantFormat), CustomError> {
        for format in [&MSMS, &MSMS_SCANS, &NOVO_MSMS_SCANS] {
            if let Ok(peptide) = Self::parse_specific(source, format) {
                return Ok((peptide, format));
            }
        }
        Err(CustomError::error(
            "Invalid MaxQuant line",
            "The correct format could not be determined automatically",
            source.full_context(),
        ))
    }
    fn parse_specific(source: &CsvLine, format: &MaxQuantFormat) -> Result<Self, CustomError> {
        if source.fields.len() != format.number_of_columns() {
            return Err(CustomError::error(
                "Invalid MaxQuant line", 
                format!("The number of columns ({}) is not equal to the expected number of columns ({})", source.fields.len(), format.number_of_columns()), 
                source.full_context()));
        }
        let number_error = CustomError::error(
            "Invalid MaxQuant line",
            format!("This column is not a number but it is required to be a number in this MaxQuant format ({})", format.version),
            source.full_context(),
        );
        Ok(Self {
            raw_file: Location::column(format.raw_file, source).get_string(),
            scan_number: Location::column(format.scan_number, source).parse(&number_error)?,
            scan_index: Location::column(format.scan_index, source).parse(&number_error)?,
            missed_cleavages: Location::optional_column(format.missed_cleavages, source)
                .parse(&number_error)?,
            modifications: Location::column(format.modifications, source).get_string(),
            proteins: Location::column(format.proteins, source).get_string(),
            sequence: ComplexPeptide::sloppy_pro_forma(
                &source.line,
                source.fields[format.modified_sequence].clone(),
            )?,
            charge: Location::column(format.charge, source).parse(&number_error)?,
            fragmentation: Location::column(format.fragmentation, source).get_string(),
            mass_analyzer: Location::column(format.mass_analyzer, source).get_string(),
            ty: Location::column(format.ty, source).parse_with(|l| match l.as_str() {
                "PEAK" => Ok(MaxQuantType::Peak),
                "MULTI" => Ok(MaxQuantType::Multi),
                _ => Err(CustomError::error(
                    "Invalid MaxQuant line",
                    "A MaxQuant type has to be PEAK or MULTI.",
                    Context::line(
                        source.line_index,
                        source.line.clone(),
                        l.location.start,
                        l.location.len(),
                    ),
                )),
            })?,
            scan_event_number: Location::column(format.scan_event_number, source)
                .parse(&number_error)?,
            isotope_index: Location::optional_column(format.isotope_index, source)
                .parse(&number_error)?,
            mz: Location::optional_column(format.mz, source).parse(&number_error)?,
            mass: Location::optional_column(format.mass, source).parse(&number_error)?,
            mass_error_da: Location::optional_column(format.mass_error_da, source)
                .parse(&number_error)?,
            mass_error_ppm: Location::optional_column(format.mass_error_ppm, source)
                .parse(&number_error)?,
            simple_mass_error_ppm: Location::optional_column(format.simple_mass_error_ppm, source)
                .parse(&number_error)?,
            retention_time: Location::optional_column(format.retention_time, source)
                .parse(&number_error)?,
            pep: Location::column(format.pep, source).parse(&number_error)?,
            score: Location::column(format.score, source).parse(&number_error)?,
            delta_score: Location::optional_column(format.delta_score, source)
                .parse(&number_error)?,
            score_diff: Location::optional_column(format.score_diff, source)
                .parse(&number_error)?,
            localisation_probability: Location::optional_column(
                format.localisation_probability,
                source,
            )
            .parse(&number_error)?,
            precursor: Location::column(format.precursor_full_scan_number, source)
                .ignore("-1")
                .parse::<usize>(&number_error)?
                .map(|index| {
                    Ok(MaxQuantPrecursor {
                        index,
                        intensity: Location::column(format.precursor_intensity, source)
                            .parse(&number_error)?,
                        apex_function: Location::column(format.precursor_apex_function, source)
                            .parse(&number_error)?,
                        apex_offset: Location::column(format.precursor_apex_offset, source)
                            .parse(&number_error)?,
                        apex_offset_time: Location::column(
                            format.precursor_apex_offset_time,
                            source,
                        )
                        .parse(&number_error)?,
                    })
                })
                .invert()?,
            number_of_matches: Location::optional_column(format.number_of_matches, source)
                .parse(&number_error)?,
            intensity_coverage: Location::optional_column(format.intensity_coverage, source)
                .parse(&number_error)?,
            peak_coverage: Location::optional_column(format.peak_coverage, source)
                .parse(&number_error)?,
            all_modified_sequences: Location::optional_column(
                format.all_modified_sequences,
                source,
            )
            .map(|l| {
                l.array(';')
                    .map(|s| ComplexPeptide::sloppy_pro_forma(&s.line.line, s.location))
                    .collect::<Result<Vec<LinearPeptide>, CustomError>>()
            })
            .invert()?,
            id: Location::optional_column(format.id, source).parse(&number_error)?,
            protein_group_ids: Location::optional_column(format.protein_group_ids, source)
                .parse(&number_error)?,
            peptide_id: Location::optional_column(format.peptide_id, source)
                .parse(&number_error)?,
            modified_peptide_id: Location::optional_column(format.modified_peptide_id, source)
                .parse(&number_error)?,
            evidence_id: Location::optional_column(format.evidence_id, source)
                .parse(&number_error)?,
            base_peak_intensity: Location::optional_column(format.base_peak_intensity, source)
                .parse(&number_error)?,
            total_ion_current: Location::optional_column(format.total_ion_current, source)
                .parse(&number_error)?,
            collision_energy: Location::optional_column(format.collision_energy, source)
                .parse(&number_error)?,
            dn_sequence: Location::optional_column(format.dn_sequence, source).get_string(),
            dn_combined_score: Location::optional_column(format.dn_combined_score, source)
                .parse(&number_error)?,
            dn_n_mass: Location::optional_column(format.dn_n_mass, source).parse(&number_error)?,
            dn_c_mass: Location::optional_column(format.dn_c_mass, source).parse(&number_error)?,
            dn_missing_mass: Location::optional_column(format.dn_missing_mass, source)
                .parse(&number_error)?,
            version: format.version.clone(),
        })
    }
    fn parse_file(
        path: impl AsRef<std::path::Path>,
    ) -> Result<BoxedIdentifiedPeptideIter<Self>, String> {
        parse_csv(path, b',').map(|lines| {
            Self::parse_many::<Box<dyn Iterator<Item = Self::Source>>>(Box::new(
                lines.skip(1).map(Result::unwrap),
            ))
        })
    }
}

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
