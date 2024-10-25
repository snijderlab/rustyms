use std::path::PathBuf;

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use super::{
    fasta::FastaData, novor::NovorData, opair::OpairData, peaks::PeaksData, system::MassOverCharge,
    MSFraggerData, MZTabData, MaxQuantData, SageData,
};
use crate::{
    error::CustomError, ontologies::CustomDatabase, peptide::SemiAmbiguous, system::usize::Charge,
    system::Time, LinearPeptide,
};

/// A peptide that is identified by a de novo or database matching program
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct IdentifiedPeptide {
    /// The score -1.0..=1.0 if available in the original format
    pub score: Option<f64>,
    /// The full metadata of this peptide
    pub metadata: MetaData,
}

/// The definition of all special metadata for all types of identified peptides that can be read
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub enum MetaData {
    /// Fasta metadata
    Fasta(FastaData),
    /// MaxQuant metadata
    MaxQuant(MaxQuantData),
    /// MSFragger metadata
    MSFragger(MSFraggerData),
    /// mzTab metadata
    MZTab(MZTabData),
    /// Novor metadata
    Novor(NovorData),
    /// OPair metadata
    Opair(OpairData),
    /// Peaks metadata
    Peaks(PeaksData),
    /// Sage metadata
    Sage(SageData),
}

impl IdentifiedPeptide {
    /// Get the peptide
    pub const fn peptide(&self) -> Option<&LinearPeptide<SemiAmbiguous>> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { peptide, .. })
            | MetaData::Novor(NovorData { peptide, .. })
            | MetaData::Opair(OpairData { peptide, .. })
            | MetaData::Sage(SageData { peptide, .. })
            | MetaData::Fasta(FastaData { peptide, .. })
            | MetaData::MZTab(MZTabData { peptide, .. }) => Some(peptide),
            MetaData::MSFragger(MSFraggerData { peptide, .. })
            | MetaData::MaxQuant(MaxQuantData { peptide, .. }) => peptide.as_ref(),
        }
    }

    /// Get the name of the format
    pub const fn format_name(&self) -> &'static str {
        match &self.metadata {
            MetaData::Fasta(_) => "Fasta",
            MetaData::MaxQuant(_) => "MaxQuant",
            MetaData::MSFragger(_) => "MSFragger",
            MetaData::MZTab(_) => "mzTab",
            MetaData::Novor(_) => "Novor",
            MetaData::Opair(_) => "OPair",
            MetaData::Peaks(_) => "PEAKS",
            MetaData::Sage(_) => "Sage",
        }
    }

    /// Get the identifier
    pub fn id(&self) -> String {
        match &self.metadata {
            MetaData::Peaks(PeaksData { scan, .. }) => scan.iter().join(";"),
            MetaData::Novor(NovorData { id, scan, .. }) => id.unwrap_or(*scan).to_string(),
            MetaData::Opair(OpairData { scan, .. }) => scan.to_string(),
            MetaData::Sage(SageData { id, .. }) | MetaData::MZTab(MZTabData { id, .. }) => {
                id.to_string()
            }
            MetaData::Fasta(FastaData { id, .. }) => id.clone(),
            MetaData::MSFragger(MSFraggerData { scan, .. }) => scan.to_string(),
            MetaData::MaxQuant(MaxQuantData { id, scan, .. }) => {
                id.map_or_else(|| scan.iter().join(";"), |id| id.to_string())
            }
        }
    }

    /// Get the local confidence, it is the same length as the peptide with a local score in 0..=1
    pub fn local_confidence(&self) -> Option<&[f64]> {
        match &self.metadata {
            MetaData::Peaks(PeaksData {
                local_confidence, ..
            }) => Some(local_confidence),
            MetaData::Novor(NovorData {
                local_confidence, ..
            })
            | MetaData::MZTab(MZTabData {
                local_confidence, ..
            }) => local_confidence.as_deref(),
            _ => None,
        }
    }

    /// The charge of the precursor, if known
    pub const fn charge(&self) -> Option<Charge> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { z, .. })
            | MetaData::Novor(NovorData { z, .. })
            | MetaData::Opair(OpairData { z, .. })
            | MetaData::Sage(SageData { z, .. })
            | MetaData::MSFragger(MSFraggerData { z, .. })
            | MetaData::MaxQuant(MaxQuantData { z, .. })
            | MetaData::MZTab(MZTabData { z, .. }) => Some(*z),
            MetaData::Fasta(_) => None,
        }
    }

    /// Which fragmentation mode was used, if known
    pub fn mode(&self) -> Option<&str> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { mode, .. }) => Some(mode),
            MetaData::MaxQuant(MaxQuantData { fragmentation, .. }) => fragmentation.as_deref(),
            _ => None,
        }
    }

    /// The retention time, if known
    pub fn retention_time(&self) -> Option<Time> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { rt, .. })
            | MetaData::Opair(OpairData { rt, .. })
            | MetaData::Sage(SageData { rt, .. })
            | MetaData::MSFragger(MSFraggerData { rt, .. }) => Some(*rt),
            MetaData::MaxQuant(MaxQuantData { rt, .. })
            | MetaData::Novor(NovorData { rt, .. })
            | MetaData::MZTab(MZTabData { rt, .. }) => *rt,
            MetaData::Fasta(_) => None,
        }
    }

    /// The scans per rawfile that are at the basis for this identified peptide, if the rawfile is unknown there will be one
    pub fn scans(&self) -> SpectrumIds {
        match &self.metadata {
            MetaData::Peaks(PeaksData { raw_file, scan, .. }) => raw_file.clone().map_or_else(
                || {
                    SpectrumIds::FileNotKnown(
                        scan.iter()
                            .flat_map(|s| s.scans.clone())
                            .map(SpectrumId::Index)
                            .collect(),
                    )
                },
                |raw_file| {
                    SpectrumIds::FileKnown(vec![(
                        raw_file,
                        scan.iter()
                            .flat_map(|s| s.scans.clone())
                            .map(SpectrumId::Index)
                            .collect(),
                    )])
                },
            ),
            MetaData::Novor(NovorData { scan, .. }) => {
                SpectrumIds::FileNotKnown(vec![SpectrumId::Index(*scan)])
            }
            MetaData::Opair(OpairData { raw_file, scan, .. }) => {
                SpectrumIds::FileKnown(vec![(raw_file.clone(), vec![SpectrumId::Index(*scan)])])
            }
            MetaData::MaxQuant(MaxQuantData { raw_file, scan, .. }) => {
                SpectrumIds::FileKnown(vec![(
                    raw_file.clone(),
                    scan.iter().copied().map(SpectrumId::Index).collect(),
                )])
            }
            MetaData::MZTab(MZTabData { spectra_ref, .. }) => spectra_ref.clone(),
            MetaData::MSFragger(MSFraggerData { raw_file, scan, .. })
            | MetaData::Sage(SageData { raw_file, scan, .. }) => {
                SpectrumIds::FileKnown(vec![(raw_file.clone(), vec![scan.clone()])])
            }
            MetaData::Fasta(_) => SpectrumIds::None,
        }
    }

    /// Get the mz as experimentally determined
    pub fn experimental_mz(&self) -> Option<MassOverCharge> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { mz, .. })
            | MetaData::Novor(NovorData { mz, .. })
            | MetaData::Opair(OpairData { mz, .. })
            | MetaData::MSFragger(MSFraggerData { mz, .. }) => Some(*mz),
            MetaData::MZTab(MZTabData { mz, .. }) | MetaData::MaxQuant(MaxQuantData { mz, .. }) => {
                *mz
            }
            MetaData::Sage(SageData {
                mass: experimental_mass,
                z,
                ..
            }) => Some(MassOverCharge::new::<crate::system::mz>(
                experimental_mass.value / (z.value as f64),
            )),
            MetaData::Fasta(_) => None,
        }
    }

    /// Get the mass as experimentally determined
    pub fn experimental_mass(&self) -> Option<crate::system::Mass> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { mass, .. })
            | MetaData::Novor(NovorData { mass, .. })
            | MetaData::Opair(OpairData { mass, .. })
            | MetaData::MSFragger(MSFraggerData { mass, .. })
            | MetaData::Sage(SageData { mass, .. }) => Some(*mass),
            MetaData::MaxQuant(MaxQuantData { mass, .. }) => *mass,
            MetaData::MZTab(MZTabData { mz, z, .. }) => mz.map(|mz| mz * z.to_float()),
            MetaData::Fasta(_) => None,
        }
    }

    /// Get the absolute ppm error between the experimental and theoretical precursor mass
    pub fn ppm_error(&self) -> Option<crate::system::Ratio> {
        let exp_mass = self.experimental_mass()?;
        let theo_mass = self
            .peptide()
            .and_then(|p| p.formulas().to_vec().pop())
            .map(|f| f.monoisotopic_mass())?;

        Some(theo_mass.ppm(exp_mass))
    }

    /// Get the absolute mass error between the experimental and theoretical precursor mass
    pub fn mass_error(&self) -> Option<crate::system::Mass> {
        let exp_mass = self.experimental_mass()?;
        let theo_mass = self
            .peptide()
            .and_then(|p| p.formulas().to_vec().pop())
            .map(|f| f.monoisotopic_mass())?;

        Some((exp_mass - theo_mass).abs())
    }
}

/// Multiple spectrum identifiers
#[derive(Clone, Default, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum SpectrumIds {
    /// When no spectra references are knwon at all
    #[default]
    None,
    /// When the source file is now known
    FileNotKnown(Vec<SpectrumId>),
    /// When the source file is known, grouped per file
    FileKnown(Vec<(PathBuf, Vec<SpectrumId>)>),
}

/// A spectrum identifier
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum SpectrumId {
    /// A native id, the format differs between vendors
    Native(String),
    /// A spectrum index
    Index(usize),
}

impl Default for SpectrumId {
    fn default() -> Self {
        Self::Index(0)
    }
}

impl std::fmt::Display for SpectrumId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Index(i) => write!(f, "{i}"),
            Self::Native(n) => write!(f, "{n}"),
        }
    }
}

impl SpectrumId {
    /// Get the index if this is an index
    pub const fn index(&self) -> Option<usize> {
        match self {
            Self::Index(i) => Some(*i),
            Self::Native(_) => None,
        }
    }

    /// Get the native ID if this is a native ID
    pub fn native(&self) -> Option<&str> {
        match self {
            Self::Index(_) => None,
            Self::Native(n) => Some(n),
        }
    }
}

/// The required methods for any source of identified peptides
pub trait IdentifiedPeptideSource
where
    Self: std::marker::Sized,
{
    /// The source data where the peptides are parsed form
    type Source;
    /// The format type
    type Format: Clone;
    /// Parse a single identified peptide from its source and return the detected format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse(
        source: &Self::Source,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<(Self, &'static Self::Format), CustomError>;
    /// Parse a single identified peptide with the given format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_specific(
        source: &Self::Source,
        format: &Self::Format,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, CustomError>;
    /// Parse a source of multiple peptides automatically determining the format to use by the first item
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_many<I: Iterator<Item = Result<Self::Source, CustomError>>>(
        iter: I,
        custom_database: Option<&CustomDatabase>,
    ) -> IdentifiedPeptideIter<Self, I> {
        IdentifiedPeptideIter {
            iter: Box::new(iter),
            format: None,
            custom_database,
        }
    }
    /// Parse a file with identified peptides.
    /// # Errors
    /// Returns Err when the file could not be opened
    fn parse_file(
        path: impl AsRef<std::path::Path>,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<BoxedIdentifiedPeptideIter<Self>, CustomError>;
}

/// Convenience type to not have to type out long iterator types
pub type BoxedIdentifiedPeptideIter<'lifetime, T> = IdentifiedPeptideIter<
    'lifetime,
    T,
    Box<dyn Iterator<Item = Result<<T as IdentifiedPeptideSource>::Source, CustomError>>>,
>;

/// An iterator returning parsed identified peptides
pub struct IdentifiedPeptideIter<
    'lifetime,
    R: IdentifiedPeptideSource,
    I: Iterator<Item = Result<R::Source, CustomError>>,
> {
    iter: Box<I>,
    format: Option<R::Format>,
    custom_database: Option<&'lifetime CustomDatabase>,
}

impl<'lifetime, R: IdentifiedPeptideSource, I: Iterator<Item = Result<R::Source, CustomError>>>
    Iterator for IdentifiedPeptideIter<'lifetime, R, I>
where
    R::Format: 'static,
{
    type Item = Result<R, CustomError>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(format) = &self.format {
            self.iter
                .next()
                .map(|source| R::parse_specific(&source?, format, self.custom_database))
        } else {
            match self
                .iter
                .next()
                .map(|source| R::parse(&source?, self.custom_database))
            {
                None => None,
                Some(Ok((pep, format))) => {
                    self.format = Some(format.clone());
                    Some(Ok(pep))
                }
                Some(Err(e)) => Some(Err(e)),
            }
        }
    }
}

impl<'lifetime, R, I> IdentifiedPeptideIter<'lifetime, R, I>
where
    R: IdentifiedPeptideSource + Into<IdentifiedPeptide> + 'lifetime,
    I: Iterator<Item = Result<R::Source, CustomError>> + 'lifetime,
    R::Format: 'static,
{
    pub(super) fn into_box(
        self,
    ) -> Box<dyn Iterator<Item = Result<IdentifiedPeptide, CustomError>> + 'lifetime> {
        Box::new(self.map(|p: Result<R, CustomError>| match p {
            Ok(p) => Ok(p.into()),
            Err(e) => Err(e),
        }))
    }
}
