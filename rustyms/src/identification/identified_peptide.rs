use std::path::Path;

use serde::{Deserialize, Serialize};

use super::{
    fasta::FastaData, novor::NovorData, opair::OpairData, peaks::PeaksData, MSFraggerData,
    MaxQuantData, SageData,
};
use crate::{
    error::CustomError, ontologies::CustomDatabase, peptide::VerySimple, system::usize::Charge,
    system::Time, LinearPeptide,
};

/// A peptide that is identified by a de novo or database matching program
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct IdentifiedPeptide {
    /// The score -1.0..=1.0 if available in the original format
    pub score: Option<f64>,
    /// The full metadata of this peptide
    pub metadata: MetaData,
}

/// The definition of all special metadata for all types of identified peptides that can be read
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub enum MetaData {
    /// No metadata
    #[default]
    None,
    /// Peaks metadata
    Peaks(PeaksData),
    /// Novor metadata
    Novor(NovorData),
    /// OPair metadata
    Opair(OpairData),
    /// Fasta metadata
    Fasta(FastaData),
    /// MaxQuant metadata
    MaxQuant(MaxQuantData),
    /// Sage metadata
    Sage(SageData),
    /// MSFragger metadata
    MSFragger(MSFraggerData),
}

impl MetaData {
    /// Get the peptide
    pub const fn peptide(&self) -> Option<&LinearPeptide<VerySimple>> {
        match self {
            Self::Peaks(PeaksData { peptide, .. })
            | Self::Novor(NovorData { peptide, .. })
            | Self::Opair(OpairData { peptide, .. })
            | Self::Sage(SageData { peptide, .. })
            | Self::Fasta(FastaData { peptide, .. }) => Some(peptide),
            Self::MSFragger(MSFraggerData { peptide, .. })
            | Self::MaxQuant(MaxQuantData { peptide, .. }) => peptide.as_ref(),
            Self::None => None,
        }
    }

    /// Get the local confidence, it is the same lengths as the peptide with a local score in 0..=1
    pub fn local_confidence(&self) -> Option<&[f64]> {
        match self {
            Self::Peaks(PeaksData {
                local_confidence, ..
            }) => Some(local_confidence),
            Self::Novor(NovorData {
                local_confidence, ..
            }) => local_confidence.as_deref(),
            _ => None,
        }
    }

    /// The charge of the precursor, if known
    pub const fn charge(&self) -> Option<Charge> {
        match self {
            Self::Peaks(PeaksData { z, .. })
            | Self::Novor(NovorData { z, .. })
            | Self::Opair(OpairData { z, .. })
            | Self::Sage(SageData { z, .. })
            | Self::MSFragger(MSFraggerData { z, .. })
            | Self::MaxQuant(MaxQuantData { z, .. }) => Some(*z),
            Self::Fasta(_) | Self::None => None,
        }
    }

    /// Which fragmentation mode was used, if known
    pub fn mode(&self) -> Option<&str> {
        match self {
            Self::Peaks(PeaksData { mode, .. }) => Some(mode),
            Self::MaxQuant(MaxQuantData { fragmentation, .. }) => fragmentation.as_deref(),
            _ => None,
        }
    }

    /// The retention time, if known
    pub fn retention_time(&self) -> Option<Time> {
        match self {
            Self::Peaks(PeaksData { rt, .. })
            | Self::Opair(OpairData { rt, .. })
            | Self::Sage(SageData { rt, .. })
            | Self::MSFragger(MSFraggerData { rt, .. }) => Some(*rt),
            Self::MaxQuant(MaxQuantData { rt, .. }) | Self::Novor(NovorData { rt, .. }) => *rt,
            Self::Fasta(_) | Self::None => None,
        }
    }

    /// The scan indices of the spectrum for this identified peptide, if known.
    pub fn scan_indices(&self) -> Option<Vec<usize>> {
        match self {
            Self::Peaks(PeaksData { scan, .. }) => {
                Some(scan.iter().flat_map(|i| &i.scans).copied().collect())
            }
            Self::Novor(NovorData { scan, .. }) | Self::Opair(OpairData { scan, .. }) => {
                Some(vec![*scan])
            }
            Self::MaxQuant(MaxQuantData { scan_number, .. }) => Some(scan_number.clone()),
            Self::MSFragger(MSFraggerData { spectrum, .. }) => Some(vec![spectrum.scan.0]),
            Self::Sage(_) | Self::Fasta(_) | Self::None => None,
        }
    }

    /// The native ids of the spectrum for this identified peptide, if known.
    pub fn spectrum_native_ids(&self) -> Option<Vec<String>> {
        match self {
            Self::Sage(SageData { native_id, .. }) => Some(vec![native_id.clone()]),
            Self::MaxQuant(_)
            | Self::Opair(_)
            | Self::Novor(_)
            | Self::Peaks(_)
            | Self::Fasta(_)
            | Self::MSFragger(_)
            | Self::None => None,
        }
    }

    /// Get the file name for the raw file
    pub fn raw_file(&self) -> Option<&Path> {
        match self {
            Self::Peaks(PeaksData { raw_file, .. }) => raw_file.as_deref(),
            Self::Opair(OpairData { raw_file, .. })
            | Self::MaxQuant(MaxQuantData { raw_file, .. })
            | Self::Sage(SageData { raw_file, .. }) => Some(raw_file),
            Self::MSFragger(MSFraggerData { spectrum, .. }) => Some(&spectrum.file),
            Self::Novor(_) | Self::Fasta(_) | Self::None => None,
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
