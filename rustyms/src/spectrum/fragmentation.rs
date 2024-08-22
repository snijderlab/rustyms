use crate::{system::MassOverCharge, CompoundPeptidoform, Fragment, MassMode, Model};

use super::AnnotatedSpectrum;

pub trait Fragments {
    type Tolerance: From<crate::Tolerance<MassOverCharge>> + Copy;
    fn empty_annotated(&self, peptide: CompoundPeptidoform) -> AnnotatedSpectrum;
    fn search(&self, query: MassOverCharge, tolerance: Self::Tolerance) -> Option<usize>;

    /// Annotate this spectrum with the given peptidoform and given fragments see [`crate::CompoundPeptidoform::generate_theoretical_fragments`].
    fn annotate(
        &self,
        peptide: CompoundPeptidoform,
        theoretical_fragments: &[Fragment],
        model: &Model,
        mode: MassMode,
    ) -> AnnotatedSpectrum {
        let tolerance = model.tolerance.into();
        let mut annotated = Self::empty_annotated(&self, peptide.clone());

        for fragment in theoretical_fragments {
            // Determine fragment mz and see if it is within the model range.
            let mz = fragment.mz(mode);
            if !model.mz_range.contains(&mz) {
                continue;
            }

            // Get the index of the element closest to this value (spectrum is defined to always be sorted)
            if let Some(index) = Self::search(&self, mz, tolerance) {
                annotated.spectrum[index]
                    .annotation
                    .push((fragment.clone(), Vec::new()));
            }
        }

        annotated
    }
}
