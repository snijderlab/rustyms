use super::fasta::FastaData;
use super::novor::NovorData;
use super::opair::OpairData;
use super::peaks::PeaksData;
use crate::error::CustomError;
use crate::LinearPeptide;

/// A peptide that is identified by a de novo or database matching program
#[derive(Debug)]
pub struct IdentifiedPeptide {
    /// The peptide sequence
    pub peptide: LinearPeptide,
    /// The local confidence of this peptide (same length as the peptide)
    pub local_confidence: Option<Vec<f64>>,
    /// The score [0-1] if available in the original format
    pub score: Option<f64>,
    /// The full metadata of this peptide
    pub metadata: MetaData,
}

// impl IdentifiedPeptide {
//     pub fn parse_file(
//         path: impl AsRef<std::path::Path>,
//     ) -> Result<Box<impl Iterator<Item = Result<IdentifiedPeptide, CustomError>>>, String> {
//         let path = path.as_ref();
//         let actual_extension = path
//             .extension()
//             .map(|ex: &std::ffi::OsStr| {
//                 (ex == "gz")
//                     .then(|| path.to_str().and_then(|p| p.rsplitn(3, '.').nth(1)))
//                     .unwrap_or(ex.to_str())
//             })
//             .flatten()
//             .map(|ex| ex.to_lowercase());
//         match actual_extension.as_deref() {
//             Some("csv") => PeaksData::parse_file(path).map(|p| Box::new(p.generalise())),
//             Some("psmtsv") => OpairData::parse_file(path).map(|p| Box::new(p.generalise())),
//             Some("fasta") => {
//                 FastaData::parse_reads(path).map(|p| Box::new(p.into_iter().map(|f| Ok(f.into()))))
//             }
//             _ => Err(
//                 "Could not detect file format for the given identified peptides file".to_string(),
//             ),
//         }
//     }
// }

/// The definition of all special metadata for all types of identified peptides that can be read
#[derive(Debug)]
pub enum MetaData {
    /// Peaks metadata
    Peaks(PeaksData),
    /// Novor metadata
    Novor(NovorData),
    /// Novor metadata
    Opair(OpairData),
    /// Novor metadata
    Fasta(FastaData),
}

impl MetaData {
    /// The charge of the precursor, if known
    pub const fn charge(&self) -> Option<usize> {
        match self {
            Self::Peaks(PeaksData { z, .. })
            | Self::Novor(NovorData { z, .. })
            | Self::Opair(OpairData { z, .. }) => Some(*z),
            Self::Fasta(_) => None,
        }
    }
    /// Which fragmentation mode was used, if known
    pub fn mode(&self) -> Option<&str> {
        match self {
            Self::Peaks(PeaksData { mode, .. }) => Some(mode),
            _ => None,
        }
    }
    /// Which fragmentation mode was used, if known
    pub fn scan_number(&self) -> Option<usize> {
        match self {
            Self::Peaks(PeaksData { scan, .. }) => {
                scan.first().and_then(|i| i.scans.first().copied())
            }
            Self::Novor(NovorData { scan, .. }) | Self::Opair(OpairData { scan, .. }) => {
                Some(*scan)
            }
            Self::Fasta(_) => None,
        }
    }
}

/// The required methods for any source of identified peptides
pub trait IdentifiedPeptideSource
where
    Self: std::marker::Sized,
{
    /// The source data where the peptides are parsed form eg [`CsvLine`]
    type Source;
    /// The format type eg [`PeaksFormat`]
    type Format: Clone;
    /// Parse a single identified peptide from its source and return the detected format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse(source: &Self::Source) -> Result<(Self, &'static Self::Format), CustomError>;
    /// Parse a single identified peptide with the given format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_specific(source: &Self::Source, format: &Self::Format) -> Result<Self, CustomError>;
    /// Parse a source of multiple peptides automatically determining the format to use by the first item
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_many<I: Iterator<Item = Self::Source>>(iter: I) -> IdentifiedPeptideIter<Self, I> {
        IdentifiedPeptideIter {
            iter: Box::new(iter),
            format: None,
        }
    }
    /// Parse a file with identified peptides.
    /// # Errors
    /// Returns Err when the file could not be opened
    fn parse_file(
        path: impl AsRef<std::path::Path>,
    ) -> Result<BoxedIdentifiedPeptideIter<Self>, String>;
}

/// Convenience type to not have to type out long iterator types
pub type BoxedIdentifiedPeptideIter<T> =
    IdentifiedPeptideIter<T, Box<dyn Iterator<Item = <T as IdentifiedPeptideSource>::Source>>>;

// impl<T: Into<IdentifiedPeptide> + IdentifiedPeptideSource> BoxedIdentifiedPeptideIter<T> {
//     fn generalise(self) -> impl Iterator<Item = Result<IdentifiedPeptide, CustomError>> {
//         self.into_iter().map(|element| element.map(|p| p.into()))
//     }
// }

/// An iterator returning parsed identified peptides
pub struct IdentifiedPeptideIter<R: IdentifiedPeptideSource, I: Iterator<Item = R::Source>> {
    iter: Box<I>,
    format: Option<R::Format>,
}

impl<R: IdentifiedPeptideSource, I: Iterator<Item = R::Source>> Iterator
    for IdentifiedPeptideIter<R, I>
where
    R::Format: 'static,
{
    type Item = Result<R, CustomError>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(format) = &self.format {
            self.iter
                .next()
                .map(|source| R::parse_specific(&source, format))
        } else {
            match self.iter.next().map(|source| R::parse(&source)) {
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
