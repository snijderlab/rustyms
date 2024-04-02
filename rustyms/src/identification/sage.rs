use crate::{
    error::CustomError,
    system::{Charge, Mass, Time},
    LinearPeptide,
};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use super::{
    common_parser::Location,
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid Sage line",
    "This column is not a number but it is required to be a number in this Sage format",
);

// TODO: Handle scan number parsing, think about score
format_family!(
    /// The format for any Sage file
    SageFormat,
    /// The data from any Sage file
    SageData,
    SageVersion, [&VERSION_0_14], b'\t';
    required {
        psm_id: usize, |location: Location| location.parse(NUMBER_ERROR);
        peptide: LinearPeptide, |location: Location| LinearPeptide::pro_forma(location.as_str());
        proteins: Vec<String>, |location: Location| Ok(location.get_string().split(';').map(ToString::to_string).collect_vec());
        num_proteins: usize, |location: Location| location.parse(NUMBER_ERROR);
        filename: String, |location: Location| Ok(location.get_string());
        scan_nr: String, |location: Location| Ok(location.get_string());
        rank: usize, |location: Location| location.parse(NUMBER_ERROR);
        label: i8, |location: Location| location.parse(NUMBER_ERROR);
        experimental_mass: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        theoretical_mass: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        z: Charge, |location: Location| location.parse::<usize>(NUMBER_ERROR).map(|c| Charge::new::<crate::system::e>(c as f64));
        missed_cleavages: usize, |location: Location| location.parse(NUMBER_ERROR);
        semi_enzymatic: bool, |location: Location| location.parse::<u8>(NUMBER_ERROR).map(|n| n != 0);
        isotope_error: f64, |location: Location| location.parse(NUMBER_ERROR);
        precursor_ppm: f64, |location: Location| location.parse(NUMBER_ERROR);
        fragment_ppm: f64, |location: Location| location.parse(NUMBER_ERROR);
        hyperscore: f64, |location: Location| location.parse(NUMBER_ERROR);
        delta_next: f64, |location: Location| location.parse(NUMBER_ERROR);
        delta_best: f64, |location: Location| location.parse(NUMBER_ERROR);
        rt: Time, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        aligned_rt: f64, |location: Location| location.parse(NUMBER_ERROR);
        predicted_rt: f64, |location: Location| location.parse(NUMBER_ERROR);
        delta_rt_model: f64, |location: Location| location.parse(NUMBER_ERROR);
        ion_mobility: f64, |location: Location| location.parse(NUMBER_ERROR);
        predicted_mobility: f64, |location: Location| location.parse(NUMBER_ERROR);
        delta_mobility: f64, |location: Location| location.parse(NUMBER_ERROR);
        matched_peaks: usize, |location: Location| location.parse(NUMBER_ERROR);
        longest_b: usize, |location: Location| location.parse(NUMBER_ERROR);
        longest_y: usize, |location: Location| location.parse(NUMBER_ERROR);
        longest_y_pct: f64, |location: Location| location.parse(NUMBER_ERROR);
        matched_intensity_pct: f64, |location: Location| location.parse(NUMBER_ERROR);
        scored_candidates: usize, |location: Location| location.parse(NUMBER_ERROR);
        poisson: f64, |location: Location| location.parse(NUMBER_ERROR);
        sage_discriminant_score: f64, |location: Location| location.parse(NUMBER_ERROR);
        posterior_error: f64, |location: Location| location.parse(NUMBER_ERROR);
        spectrum_q: f64, |location: Location| location.parse(NUMBER_ERROR);
        peptide_q: f64, |location: Location| location.parse(NUMBER_ERROR);
        protein_q: f64, |location: Location| location.parse(NUMBER_ERROR);
        ms2_intensity: f64, |location: Location| location.parse(NUMBER_ERROR);
    }
    optional { }
);

impl From<SageData> for IdentifiedPeptide {
    fn from(value: SageData) -> Self {
        Self {
            peptide: value.peptide.clone(),
            local_confidence: None,
            score: None,
            metadata: MetaData::Sage(value),
        }
    }
}

/// An older version of a Sage export
pub const VERSION_0_14: SageFormat = SageFormat {
    version: SageVersion::Version_0_14,
    psm_id: "psm_id",
    peptide: "peptide",
    proteins: "proteins",
    num_proteins: "num_proteins",
    filename: "filename",
    scan_nr: "scannr",
    rank: "rank",
    label: "label",
    experimental_mass: "expmass",
    theoretical_mass: "calcmass",
    z: "charge",
    missed_cleavages: "missed_cleavages",
    semi_enzymatic: "semi_enzymatic",
    isotope_error: "isotope_error",
    precursor_ppm: "precursor_ppm",
    fragment_ppm: "fragment_ppm",
    hyperscore: "hyperscore",
    delta_next: "delta_next",
    delta_best: "delta_best",
    rt: "rt",
    aligned_rt: "aligned_rt",
    predicted_rt: "predicted_rt",
    delta_rt_model: "delta_rt_model",
    ion_mobility: "ion_mobility",
    predicted_mobility: "predicted_mobility",
    delta_mobility: "delta_mobility",
    matched_peaks: "matched_peaks",
    longest_b: "longest_b",
    longest_y: "longest_y",
    longest_y_pct: "longest_y_pct",
    matched_intensity_pct: "matched_intensity_pct",
    scored_candidates: "scored_candidates",
    poisson: "poisson",
    sage_discriminant_score: "sage_discriminant_score",
    posterior_error: "posterior_error",
    spectrum_q: "spectrum_q",
    peptide_q: "peptide_q",
    protein_q: "protein_q",
    ms2_intensity: "ms2_intensity",
};

/// All possible Sage versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
pub enum SageVersion {
    /// Current sage version
    #[default]
    Version_0_14,
}

impl std::fmt::Display for SageVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::Version_0_14 => "v0.14",
            }
        )
    }
}
