pub(crate) trait Fragments {
    type Model;
    type Fragment;
    type AdditionalInfo;
    fn generate_theoretical_fragments(
        &self,
        model: &Self::Model,
        peptide_index: usize,
        additional_indo: &Self::AdditionalInfo,
    ) -> Vec<Self::Fragment>;

    // One function to calculate all breakages
    // One function to go from breakages to fragments
    // One function to go from fragments to ions
    // One function to go from ions to charged ions
}
