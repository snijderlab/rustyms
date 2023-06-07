# Match those fragments!

A work in progress peptide fragmentation matching library for rust. For now only supports mgf raw files and pro forma sequences.

## Future ideas
 - [ ] Add peaks merging, merge low intensity peaks into high intensity peaks within a small ppm tolerance (look at hecklib for inspiration)
 - [ ] Median correction, determine the median ppm error and shift the whole spectrum with this amount, followed by a new round of matching