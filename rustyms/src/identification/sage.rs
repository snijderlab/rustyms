use crate::{
    error::CustomError,
    peptide::VerySimple,
    system::{usize::Charge, Mass, Ratio, Time},
    LinearPeptide,
};
use itertools::Itertools;
use regex::Regex;
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

format_family!(
    /// The format for any Sage file
    SageFormat,
    /// The data from any Sage file
    SageData,
    SageVersion, [&VERSION_0_14], b'\t';
    required {
        psm_id: usize, |location: Location| location.parse(NUMBER_ERROR);
        peptide: LinearPeptide<VerySimple>, |location: Location| LinearPeptide::pro_forma(location.as_str(), None).map(|p|p.very_simple().unwrap());
        proteins: Vec<String>, |location: Location| Ok(location.get_string().split(';').map(ToString::to_string).collect_vec());
        filename: String, |location: Location| Ok(location.get_string());
        scan_nr: (usize,usize,usize), |location: Location|
        location.parse_regex(
            &Regex::new("controllerType=(\\d+) controllerNumber=(\\d+) scan=(\\d+)").unwrap(),
             ("Invalid Sage line", "Scan number does not correspond to the expected pattern 'controllerType=0 controllerNumber=1 scan=17666'."))
             .map(|groups| (groups[1].parse().unwrap(), groups[2].parse().unwrap(), groups[3].parse().unwrap()));
        rank: usize, |location: Location| location.parse(NUMBER_ERROR);
        decoy: bool, |location: Location| location.parse::<i8>(NUMBER_ERROR).map(|v| v == -1);
        experimental_mass: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        theoretical_mass: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        z: Charge, |location: Location| location.parse::<usize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        missed_cleavages: usize, |location: Location| location.parse(NUMBER_ERROR);
        semi_enzymatic: bool, |location: Location| location.parse::<u8>(NUMBER_ERROR).map(|n| n != 0);
        isotope_error: f64, |location: Location| location.parse(NUMBER_ERROR);
        precursor_ppm: Ratio, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::ppm>);
        fragment_ppm: Ratio, |location: Location| location.parse(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::ppm>);
        hyperscore: f64, |location: Location| location.parse(NUMBER_ERROR);
        delta_next: f64, |location: Location| location.parse(NUMBER_ERROR);
        delta_best: f64, |location: Location| location.parse(NUMBER_ERROR);
        rt: Time, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        aligned_rt: Ratio, |location: Location| location.parse(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::fraction>);
        predicted_rt: Ratio, |location: Location| location.parse(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::fraction>);
        delta_rt_model: f64, |location: Location| location.parse(NUMBER_ERROR);
        ion_mobility: f64, |location: Location| location.parse(NUMBER_ERROR);
        predicted_mobility: f64, |location: Location| location.parse(NUMBER_ERROR);
        delta_mobility: f64, |location: Location| location.parse(NUMBER_ERROR);
        matched_peaks: usize, |location: Location| location.parse(NUMBER_ERROR);
        longest_b: usize, |location: Location| location.parse(NUMBER_ERROR);
        longest_y: usize, |location: Location| location.parse(NUMBER_ERROR);
        matched_intensity_pct: Ratio, |location: Location| location.parse(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::percent>);
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
            score: Some(value.sage_discriminant_score.clamp(-1.0, 1.0)),
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
    filename: "filename",
    scan_nr: "scannr",
    rank: "rank",
    decoy: "label",
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
#[allow(non_camel_case_types)]
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
