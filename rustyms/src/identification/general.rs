use std::path::Path;

use super::{
    error::{Context, CustomError},
    ontologies::CustomDatabase,
    DeepNovoData, FastaData, IdentifiedPeptide, IdentifiedPeptideIter, IdentifiedPeptideSource,
    MSFraggerData, MZTabData, MaxQuantData, NovorData, OpairData, PeaksData, SageData,
};

// TODO:
// * Merge multiple annotations for the same spectrum (e.g. all candidates peaks export, take care not to lose info on chimeric spectra)
// * Merge identical (or similar?) peptide sequences (for faster processing)

/// Open the selected path and automatically determine the file type. It will uncompress gzipped
/// files automatically.
///
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
            .or_else(|pe| {
                NovorData::parse_file(path, custom_database)
                    .map(IdentifiedPeptideIter::into_box)
                    .map_err(|ne| (pe, ne))
            })
            .map_err(|(pe, ne)| {
                CustomError::error(
                    "Unknown file",
                    "Could not be recognised as either a Peaks or Novor file",
                    Context::show(path.to_string_lossy()),
                )
                .with_underlying_errors(vec![pe, ne])
            }),
        Some("tsv") => MSFraggerData::parse_file(path, custom_database)
            .map(IdentifiedPeptideIter::into_box)
            .or_else(|me| {
                SageData::parse_file(path, custom_database)
                    .map(IdentifiedPeptideIter::into_box)
                    .map_err(|se| (me, se))
            })
            .map_err(|(me, se)| {
                CustomError::error(
                    "Unknown file",
                    "Could not be recognised as either a MSFragger or Sage file",
                    Context::show(path.to_string_lossy()),
                )
                .with_underlying_errors(vec![me, se])
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
        Some("mztab") => MZTabData::parse_file(path, custom_database).map(|peptides| {
            Box::new(peptides.into_iter().map(|p| p.map(Into::into)))
                as Box<dyn Iterator<Item = Result<IdentifiedPeptide, CustomError>> + 'a>
        }),
        Some("deepnovo_denovo") => {
            DeepNovoData::parse_file(path, custom_database).map(IdentifiedPeptideIter::into_box)
        }
        _ => Err(CustomError::error(
            "Unknown extension",
            "Use CSV, TSV, TXT, PSMTSV, deepnovo_denovo, or Fasta, or any of these as a gzipped file (eg csv.gz).",
            Context::show(path.to_string_lossy()),
        )),
    }
}

#[allow(clippy::missing_panics_doc)]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn open_sage() {
        let peptides =
            open_identified_peptides_file("src/identification/test_files/sage_v0_14.tsv", None)
                .unwrap();
        let mut num_peptides = 0;
        for peptide in peptides {
            let read: IdentifiedPeptide = peptide.unwrap();
            num_peptides += 1;
            let _ = read.peptide().unwrap();
        }
        assert_eq!(num_peptides, 19);
    }

    #[test]
    fn open_msfragger() {
        let peptides =
            open_identified_peptides_file("src/identification/test_files/msfragger_v21.tsv", None)
                .unwrap();
        let mut num_peptides = 0;
        for peptide in peptides {
            let read: IdentifiedPeptide = peptide.unwrap();
            num_peptides += 1;
            let _ = read.peptide().unwrap();
        }
        assert_eq!(num_peptides, 19);
    }
}
