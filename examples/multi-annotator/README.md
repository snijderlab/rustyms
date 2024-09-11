# Multi annotator

Usage:
```
cargo run --release --bin multi-annotator -- --in-path .\CIDcurves_file_to_match.csv --out-path out.csv
```
Note: the examples files are not present.

This takes a CSV file as input that contains a peptide and the rawfile it originated from, it then annotates the spectrum with the theoretical fragmentation from rustyms and delivers some statistics on the annotation in a resulting CSV file. This can be used to get a global impression over a whole dataset, so for example see if a certain fragmentation energy increases or decreases the coverage of a particular ion series (peptide or glycan).