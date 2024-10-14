use std::path::Path;

use super::{
    error::{Context, CustomError},
    ontologies::CustomDatabase,
    FastaData, IdentifiedPeptide, IdentifiedPeptideIter, IdentifiedPeptideSource, MSFraggerData,
    MaxQuantData, NovorData, OpairData, PeaksData, SageData,
};

// TODO:
// * Merge multiple annotations for the same spectrum (e.g. all candidates peaks export, take care not to lose info on chimeric spectra)
// * Merge identical (or similar?) peptide sequences (for faster processing)

/// Open the selected path and automatically determine the file type.
/// # Errors
/// It errors if the file type could not be determined or if opening the file errors.
pub fn open_identified_peptides_file<'a>(
    path: impl AsRef<Path>,
    custom_database: Option<&'a CustomDatabase>,
) -> Result<Box<dyn Iterator<Item = Result<IdentifiedPeptide, CustomError>> + 'a>, CustomError> {
    let path = path.as_ref();
    let actual_extension = path
        .extension()
        .map(|ex| {
            (ex == "gz")
                .then_some(path)
                .and_then(|p| p.file_stem())
                .and_then(|p| Path::new(p).extension())
                .unwrap_or(ex)
        })
        .map(|ex| ex.to_string_lossy().to_lowercase());
    match actual_extension.as_deref() {
        Some("csv") => PeaksData::parse_file(path, custom_database)
            .map(IdentifiedPeptideIter::into_box)
            .or_else(|_| {
                NovorData::parse_file(path, custom_database).map(IdentifiedPeptideIter::into_box)
            })
            .map_err(|_| {
                CustomError::error(
                    "Unknown file",
                    "Could not be recognised as either a Peaks or Novor file",
                    Context::show(path.to_string_lossy()),
                )
            }),
        Some("tsv") => MSFraggerData::parse_file(path, custom_database)
            .map(IdentifiedPeptideIter::into_box)
            .or_else(|_| {
                SageData::parse_file(path, custom_database).map(IdentifiedPeptideIter::into_box)
            })
            .map_err(|_| {
                CustomError::error(
                    "Unknown file",
                    "Could not be recognised as either a MSFragger or Sage file",
                    Context::show(path.to_string_lossy()),
                )
            }),
        Some("psmtsv") => {
            OpairData::parse_file(path, custom_database).map(IdentifiedPeptideIter::into_box)
        }
        Some("fasta") => FastaData::parse_file(path).map(|peptides| {
            Box::new(peptides.into_iter().map(|p| Ok(p.into())))
                as Box<dyn Iterator<Item = Result<IdentifiedPeptide, CustomError>> + 'a>
        }),
        Some("txt") => {
            MaxQuantData::parse_file(path, custom_database).map(IdentifiedPeptideIter::into_box)
        }
        _ => Err(CustomError::error(
            "Unknown extension",
            "Use CSV, TSV, TXT, PSMTSV, or Fasta, or any of these as a gzipped file (eg csv.gz).",
            Context::show(path.to_string_lossy()),
        )),
    }
}
