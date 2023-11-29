use crate::{
    error::CustomError,
    helper_functions::InvertResult,
    system::{Charge, Mass, MassOverCharge, Time},
    ComplexPeptide, LinearPeptide,
};

use super::{
    common_parser::Location,
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid Novor line",
    "This column is not a number but it is required to be a number in this Novor format",
);

format_family!(
    NovorFormat, NovorData, NovorVersion, [&OLD_DENOVO, &OLD_PSM, &NEW_DENOVO, &NEW_PSM], b',';
    required {
        scan: usize, |location: Location| location.parse(NUMBER_ERROR);
        mz:MassOverCharge, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        z: Charge, |location: Location| location.parse::<usize>(NUMBER_ERROR).map(|c| Charge::new::<crate::system::e>(c as f64));
        mass:  Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        ppm: f64, |location: Location| location.parse(NUMBER_ERROR);
        score: f64, |location: Location| location        .parse::<f64>(NUMBER_ERROR)        .map(|f| f / 100.0);
        peptide:  LinearPeptide, |location: Location| ComplexPeptide::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
        );
    }
    optional {
        id: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        spectra_id: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        fraction: usize, |location: Location| location
            .apply(|l| Location {
                line: l.line,
                location: l.location.start + 1..l.location.end,
            }) // Skip the F of the F{num} definition
            .parse::<usize>(NUMBER_ERROR);
        rt: Time, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        mass_err: Mass, |location: Location| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        length: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        peptide_no_ptm: String, |location: Location| Ok(Some(location.get_string()));
        protein: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        protein_start: usize, |location: Location| location.parse::<usize>(NUMBER_ERROR);
        protein_origin:String, |location: Location| Ok(Some(location.get_string()));
        protein_all: String, |location: Location| Ok(Some(location.get_string()));
        database_sequence: String, |location: Location| Ok(Some(location.get_string()));
        local_confidence: Vec<f64>, |location: Location| location.array('-')
                    .map(|l| l.parse::<f64>(NUMBER_ERROR).map(|v| v / 100.0))
                    .collect::<Result<Vec<_>, _>>();
    }
);

impl From<NovorData> for IdentifiedPeptide {
    fn from(value: NovorData) -> Self {
        Self {
            peptide: value.peptide.clone(),
            local_confidence: value.local_confidence.clone(),
            score: Some(value.score),
            metadata: MetaData::Novor(value),
        }
    }
}

/// All available Novor versions
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NovorVersion {
    /// An older version for the denovo file
    OldDenovo,
    /// An older version for the psms file
    OldPSM,
    /// Seen since v3.36.893 (not necessarily the time it was rolled out)
    NewDenovo,
    /// Seen since v3.36.893 (not necessarily the time it was rolled out)
    NewPSM,
}
impl std::fmt::Display for NovorVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::OldDenovo => "Older Denovo",
                Self::OldPSM => "Older PSM",
                Self::NewDenovo => "New Denovo",
                Self::NewPSM => "New PSM",
            }
        )
    }
}

/// The older supported format for denovo.csv files from Novor
///# -
/// 1       Fraction
/// 2       Scan #
/// 3       m/z
/// 4       z
/// 5       Score
/// 6       Peptide Mass
/// 7       Error (ppm)
/// 8       Length
/// 9       De Novo Peptide
/// 10      DB Sequence
/// <https://github.com/snijderlab/stitch/issues/156#issuecomment-1097862072>
pub const OLD_DENOVO: NovorFormat = file_format!(NovorFormat, NovorVersion::OldDenovo;
    scan: "scan #",
    mz: "m/z",
    z: "z",
    mass: "peptide mass",
    ppm: "error (ppm)",
    score: "score",
    peptide: "de novo peptide",
    id: None,
    spectra_id: None,
    fraction: Some("fraction"),
    rt: None,
    mass_err: None,
    length: Some("length"),
    peptide_no_ptm: None,
    protein: None,
    protein_start: None,
    protein_origin: None,
    protein_all: None,
    database_sequence: Some("db sequence"),
    local_confidence: None,
);

/// The older supported format for psms.csv files from Novor
///# -
/// 1       ID
/// 2       Fraction
/// 3       Scan
/// 4       m/z
/// 5       z
/// 6       Score
/// 7       Mass
/// 8       Error (ppm)
/// 9       # Proteins
/// 10      Sequence
/// <https://github.com/snijderlab/stitch/issues/156#issuecomment-1097862072>
pub const OLD_PSM: NovorFormat = file_format!(NovorFormat, NovorVersion::OldPSM;
    scan: "scan",
    mz: "m/z",
    z: "z",
    mass: "mass",
    ppm: "error (ppm)",
    score: "score",
    peptide: "sequence",
    id: Some("id"),
    spectra_id: None,
    fraction: Some("fraction"),
    rt: None,
    mass_err: None,
    length: None,
    peptide_no_ptm: None,
    protein: Some("# proteins"),
    protein_start: None,
    protein_origin: None,
    protein_all: None,
    database_sequence: None,
    local_confidence: None,
);

/// denovo: `# id, scanNum, RT, mz(data), z, pepMass(denovo), err(data-denovo), ppm(1e6*err/(mz*z)), score, peptide, aaScore,`
pub const NEW_DENOVO: NovorFormat = file_format!(NovorFormat, NovorVersion::NewDenovo;
    scan: "scannum",
    mz: "mz(data)",
    z: "z",
    mass: "pepmass(denovo)",
    ppm: "ppm(1e6*err/(mz*z))",
    score: "score",
    peptide: "peptide",
    id: Some("# id"),
    spectra_id: None,
    fraction: Some("fraction"),
    rt: None,
    mass_err: Some("err(data-denovo)"),
    length: None,
    peptide_no_ptm: None,
    protein: Some("# proteins"),
    protein_start: None,
    protein_origin: None,
    protein_all: None,
    database_sequence: None,
    local_confidence: Some("aascore"),
);

/// PSM: `#id, spectraId, scanNum, RT, mz, z, pepMass, err, ppm, score, protein, start, length, origin, peptide, noPTMPeptide, aac, allProteins`
pub const NEW_PSM: NovorFormat = file_format!(NovorFormat, NovorVersion::NewPSM;
    scan: "scannum",
    mz: "mz",
    z: "z",
    mass: "pepmass",
    ppm: "ppm",
    score: "score",
    peptide: "peptide",
    id: Some("#id"),
    spectra_id: None,
    fraction: Some("fraction"),
    rt: None,
    mass_err: Some("err(data-denovo)"),
    length: Some("length"),
    peptide_no_ptm: Some("noptmpeptide"),
    protein: Some("protein"),
    protein_start: Some("start"),
    protein_origin: Some("origin"),
    protein_all: Some("allproteins"),
    database_sequence: None,
    local_confidence: Some("aac"),
);
