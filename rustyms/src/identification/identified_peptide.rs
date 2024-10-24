use std::path::Path;

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
    /// MSFragger metadata
    MZTab(MZTabData),
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
            MetaData::None => None,
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
            MetaData::Fasta(_) | MetaData::None => None,
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
            MetaData::Fasta(_) | MetaData::None => None,
        }
    }

    /// The scan indices of the spectrum for this identified peptide, if known.
    pub fn scan_indices(&self) -> Option<Vec<usize>> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { scan, .. }) => {
                Some(scan.iter().flat_map(|i| &i.scans).copied().collect())
            }
            MetaData::Novor(NovorData { scan, .. }) | MetaData::Opair(OpairData { scan, .. }) => {
                Some(vec![*scan])
            }
            MetaData::MaxQuant(MaxQuantData { scan_number, .. }) => Some(scan_number.clone()),
            MetaData::MSFragger(MSFraggerData { spectrum, .. }) => Some(vec![spectrum.scan.0]),
            MetaData::MZTab(MZTabData { spectra_ref, .. }) => Some(
                spectra_ref
                    .iter()
                    .filter_map(|(_, _, id, _)| id.index())
                    .collect(),
            ),
            MetaData::Sage(_) | MetaData::Fasta(_) | MetaData::None => None,
        }
    }

    /// The native ids of the spectrum for this identified peptide, if known.
    pub fn spectrum_native_ids(&self) -> Option<Vec<String>> {
        match &self.metadata {
            MetaData::Sage(SageData { native_id, .. }) => Some(vec![native_id.clone()]),
            MetaData::MZTab(MZTabData { spectra_ref, .. }) => Some(
                spectra_ref
                    .iter()
                    .filter_map(|(_, _, id, _)| id.native().map(ToString::to_string))
                    .collect(),
            ),
            MetaData::MaxQuant(_)
            | MetaData::Opair(_)
            | MetaData::Novor(_)
            | MetaData::Peaks(_)
            | MetaData::Fasta(_)
            | MetaData::MSFragger(_)
            | MetaData::None => None,
        }
    }

    /// Get the file name for the raw file
    pub fn raw_file(&self) -> Option<&Path> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { raw_file, .. }) => raw_file.as_deref(),
            MetaData::Opair(OpairData { raw_file, .. })
            | MetaData::MaxQuant(MaxQuantData { raw_file, .. })
            | MetaData::Sage(SageData { raw_file, .. }) => Some(raw_file),
            MetaData::MSFragger(MSFraggerData { spectrum, .. }) => Some(&spectrum.file),
            MetaData::MZTab(MZTabData { spectra_ref, .. }) => {
                spectra_ref.first().map(|r| r.0.as_path()) // TODO: Could contain multiple files
            }
            MetaData::Novor(_) | MetaData::Fasta(_) | MetaData::None => None,
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
            MetaData::Fasta(_) | MetaData::None => None,
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
            MetaData::Fasta(_) | MetaData::None => None,
        }
    }

    /// Get the absolute ppm error between the experimental and theoretical precursor mz
    pub fn ppm_error(&self) -> Option<crate::system::Ratio> {
        let exp_mz = self.experimental_mz()?;
        let z = self.charge()?.to_float();
        let mass = self
            .peptide()
            .and_then(|p| p.formulas().to_vec().pop())
            .map(|f| f.monoisotopic_mass())?;
        let theo_mz = mass / z;

        Some(theo_mz.ppm(exp_mz))
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
