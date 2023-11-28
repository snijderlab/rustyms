use crate::{error::CustomError, helper_functions::InvertResult, ComplexPeptide, LinearPeptide};

use super::{
    common_parser::{Location, OptionalLocation},
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
};

/// The file format for any peaks format, determining the existence and location of all possible columns
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct NovorFormat {
    scan: usize,
    mz: usize,
    z: usize,
    mass: usize,
    ppm: usize,
    score: usize,
    peptide: usize,
    id: Option<usize>,
    spectra_id: Option<usize>,
    fraction: Option<usize>,
    rt: Option<usize>,
    mass_err: Option<usize>,
    length: Option<usize>,
    peptide_no_ptm: Option<usize>,
    protein: Option<usize>,
    protein_start: Option<usize>,
    protein_origin: Option<usize>,
    protein_all: Option<usize>,
    database_sequence: Option<usize>,
    local_confidence: Option<usize>,
    version: NovorVersion,
}

impl NovorFormat {
    const fn number_of_columns(&self) -> usize {
        7 + self.id.is_some() as usize
            + self.spectra_id.is_some() as usize
            + self.fraction.is_some() as usize
            + self.rt.is_some() as usize
            + self.mass_err.is_some() as usize
            + self.length.is_some() as usize
            + self.peptide_no_ptm.is_some() as usize
            + self.protein.is_some() as usize
            + self.protein_start.is_some() as usize
            + self.protein_origin.is_some() as usize
            + self.protein_all.is_some() as usize
            + self.database_sequence.is_some() as usize
            + self.local_confidence.is_some() as usize
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
pub const OLD_DENOVO: NovorFormat = NovorFormat {
    fraction: Some(0),
    scan: 1,
    mz: 2,
    z: 3,
    score: 4,
    mass: 5,
    ppm: 6,
    length: Some(7),
    peptide: 8,
    database_sequence: Some(9),
    id: None,
    spectra_id: None,
    rt: None,
    mass_err: None,
    peptide_no_ptm: None,
    protein: None,
    protein_all: None,
    protein_start: None,
    protein_origin: None,
    local_confidence: None,
    version: NovorVersion::OldDenovo,
};

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
pub const OLD_PSM: NovorFormat = NovorFormat {
    id: Some(0),
    fraction: Some(1),
    scan: 2,
    mz: 3,
    z: 4,
    score: 5,
    mass: 6,
    ppm: 7,
    protein: Some(8),
    peptide: 9,
    spectra_id: None,
    length: None,
    database_sequence: None,
    rt: None,
    mass_err: None,
    peptide_no_ptm: None,
    protein_all: None,
    protein_start: None,
    protein_origin: None,
    local_confidence: None,
    version: NovorVersion::OldPSM,
};

/// denovo: `# id, scanNum, RT, mz(data), z, pepMass(denovo), err(data-denovo), ppm(1e6*err/(mz*z)), score, peptide, aaScore,`
pub const NEW_DENOVO: NovorFormat = NovorFormat {
    id: Some(0),
    scan: 1,
    rt: Some(2),
    mz: 3,
    z: 4,
    mass: 5,
    mass_err: Some(6),
    ppm: 7,
    score: 8,
    peptide: 9,
    local_confidence: Some(10),
    protein: None,
    spectra_id: None,
    length: None,
    fraction: None,
    database_sequence: None,
    peptide_no_ptm: None,
    protein_all: None,
    protein_start: None,
    protein_origin: None,
    version: NovorVersion::NewDenovo,
};

/// PSM: `#id, spectraId, scanNum, RT, mz, z, pepMass, err, ppm, score, protein, start, length, origin, peptide, noPTMPeptide, aac, allProteins`
pub const NEW_PSM: NovorFormat = NovorFormat {
    id: Some(0),
    spectra_id: Some(1),
    scan: 2,
    rt: Some(3),
    mz: 4,
    z: 5,
    mass: 6,
    mass_err: Some(7),
    ppm: 8,
    score: 9,
    protein: Some(10),
    protein_start: Some(11),
    length: Some(12),
    protein_origin: Some(13),
    peptide: 14,
    peptide_no_ptm: Some(15),
    local_confidence: Some(16),
    protein_all: Some(17),
    fraction: None,
    database_sequence: None,
    version: NovorVersion::NewPSM,
};

/// A single parsed Novor line
#[allow(missing_docs)]
#[derive(Debug)]
pub struct NovorData {
    pub scan: usize,
    pub mz: f64,
    pub z: usize,
    pub mass: f64,
    pub ppm: f64,
    // [0-1]
    pub score: f64,
    pub peptide: LinearPeptide,
    pub id: Option<usize>,
    pub spectra_id: Option<usize>,
    // Indicated as F{num} only the num is saved here
    pub fraction: Option<usize>,
    pub rt: Option<f64>,
    pub mass_err: Option<f64>,
    pub length: Option<usize>,
    pub peptide_no_ptm: Option<String>,
    pub protein: Option<usize>,
    pub protein_start: Option<usize>,
    pub protein_origin: Option<String>,
    pub protein_all: Option<String>,
    pub database_sequence: Option<String>,
    // [0-1]
    pub local_confidence: Option<Vec<f64>>,
    pub version: NovorVersion,
}

impl IdentifiedPeptideSource for NovorData {
    type Source = CsvLine;
    type Format = NovorFormat;
    fn parse(source: &CsvLine) -> Result<(Self, &'static NovorFormat), CustomError> {
        for format in [&NEW_DENOVO, &NEW_PSM, &OLD_DENOVO, &OLD_PSM] {
            if let Ok(peptide) = Self::parse_specific(source, format) {
                return Ok((peptide, format));
            }
        }
        Err(CustomError::error(
            "Invalid Novor line",
            "The correct format could not be determined automatically",
            source.full_context(),
        ))
    }
    fn parse_specific(source: &CsvLine, format: &NovorFormat) -> Result<Self, CustomError> {
        if source.number_of_columns() != format.number_of_columns() {
            return Err(CustomError::error(
                "Invalid Novor line", 
                format!("The number of columns ({}) is not equal to the expected number of columns ({})", source.fields.len(), format.number_of_columns()), 
                source.full_context()));
        }
        let number_error = CustomError::error(
            "Invalid Novor line",
            format!("This column is not a number but it is required to be a number in this novor format ({})", format.version),
            source.full_context(),
        );
        Ok(Self {
            scan: Location::column(format.scan, source).parse(&number_error)?,
            mz: Location::column(format.mz, source).parse(&number_error)?,
            z: Location::column(format.z, source).parse(&number_error)?,
            mass: Location::column(format.mass, source).parse(&number_error)?,
            ppm: Location::column(format.ppm, source).parse(&number_error)?,
            score: Location::column(format.score, source).parse::<f64>(&number_error)? / 100.0,
            peptide: ComplexPeptide::sloppy_pro_forma(
                &source.line,
                source.fields[format.peptide].clone(),
            )?,
            id: Location::optional_column(format.id, source).parse(&number_error)?,
            spectra_id: Location::optional_column(format.spectra_id, source)
                .parse(&number_error)?,
            fraction: Location::optional_column(format.fraction, source)
                .apply(|l| Location {
                    line: l.line,
                    location: l.location.start + 1..l.location.end,
                }) // Skip the F of the F{num} definition
                .parse(&number_error)?,
            rt: Location::optional_column(format.rt, source).parse(&number_error)?,
            mass_err: Location::optional_column(format.mass_err, source).parse(&number_error)?,
            length: Location::optional_column(format.length, source).parse(&number_error)?,
            peptide_no_ptm: Location::optional_column(format.peptide_no_ptm, source).get_string(),
            protein: Location::optional_column(format.protein, source).parse(&number_error)?,
            protein_start: Location::optional_column(format.protein_start, source)
                .parse(&number_error)?,
            protein_origin: Location::optional_column(format.protein_origin, source).get_string(),
            protein_all: Location::optional_column(format.protein_all, source).get_string(),
            database_sequence: Location::optional_column(format.database_sequence, source)
                .get_string(),
            local_confidence: Location::optional_column(format.local_confidence, source)
                .map(|l| {
                    l.array('-')
                        .map(|l| l.parse::<f64>(&number_error).map(|v| v / 100.0))
                        .collect::<Result<Vec<_>, _>>()
                })
                .invert()?,
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
