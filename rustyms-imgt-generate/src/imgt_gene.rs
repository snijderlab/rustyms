use std::collections::HashMap;
use std::fmt::Display;

use itertools::Itertools;
use rustyms::AminoAcid;

use crate::shared::{AnnotatedSequence, Annotation, Gene, Region};
use crate::structs::{Location, SingleSeq};
use crate::{find_possible_n_glycan_locations, fix_j};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IMGTGene {
    pub acc: String,
    pub key: String,
    pub location: Location,
    pub allele: String,
    pub regions: HashMap<String, crate::structs::Region>,
}

impl IMGTGene {
    pub fn finish(self) -> Result<SingleSeq, String> {
        let (regions, additional_annotations) = self.get_regions()?;

        let sequence: Vec<AminoAcid> = regions.iter().flat_map(|reg| reg.1 .0.clone()).collect();
        let dna: String = regions.iter().map(|reg| reg.1 .2.clone()).collect();
        let region_lengths = regions.iter().map(|reg| (reg.0, reg.1 .0.len())).collect();
        let conserved_map = HashMap::from([
            ("1st-CYS", Annotation::Cysteine1),
            ("2nd-CYS", Annotation::Cysteine2),
            ("CONSERVED-TRP", Annotation::Tryptophan),
            ("J-PHE", Annotation::Phenylalanine),
            ("J-TRP", Annotation::Tryptophan),
        ]);
        let mut conserved = self
            .regions
            .iter()
            .filter(|(key, _)| {
                ["1st-CYS", "2nd-CYS", "CONSERVED-TRP", "J-PHE", "J-TRP"].contains(&key.as_str())
            })
            .map(|(key, region)| {
                region
                    .location
                    .find_aa_location(&regions)
                    .map(|index| (conserved_map[key.as_str()], index))
                    .ok_or(format!("Cannot find location of '{key}' '{region}'"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        conserved.extend(
            find_possible_n_glycan_locations(&sequence)
                .iter()
                .map(|i| (Annotation::NGlycan, *i)),
        );
        conserved.extend(additional_annotations);
        let (name, allele) = Gene::from_imgt_name_with_allele(self.allele.as_str())?;
        Ok(SingleSeq {
            name,
            allele,
            acc: self.acc.clone(),
            sequence: AnnotatedSequence::new(sequence.into(), region_lengths, conserved),
            dna,
        })
    }

    fn get_regions(
        &self,
    ) -> Result<
        (
            Vec<(Region, (Vec<AminoAcid>, Location, String))>,
            Vec<(Annotation, usize)>,
        ),
        String,
    > {
        let mut additional_annotations = Vec::new();

        let regions = match self.key.as_str() {
            "V-GENE" => {
                vec![
                    self.get_region(Region::FR1, "FR1-IMGT")?,
                    self.get_region(Region::CDR1, "CDR1-IMGT")?,
                    self.get_region(Region::FR2, "FR2-IMGT")?,
                    self.get_region(Region::CDR2, "CDR2-IMGT")?,
                    self.get_region(Region::FR3, "FR3-IMGT")?,
                    self.get_region(Region::CDR3, "CDR3-IMGT")?,
                ]
            }
            "J-GENE" => {
                let j = self.get_region(Region::FR4, "J-REGION")?;
                let motif = j.1 .0.iter().tuple_windows().position(|(a, b, _, d)| {
                    (*a == AminoAcid::W || *a == AminoAcid::F)
                        && *b == AminoAcid::G
                        && *d == AminoAcid::G
                });
                if let Some(motif_start) = motif {
                    let j = fix_j(j.1, motif_start);
                    additional_annotations.extend(j.1);
                    j.0
                } else {
                    vec![j] // TODO: not fully correct right, has some CDR3 as well, and has quite some conserved residues
                }
            }
            "D-GENE" => {
                vec![self.get_region(Region::CDR3, "D-REGION")?]
            }
            "C-GENE" => {
                let mut seq = Vec::new();

                // Heavy chain
                seq.extend(self.get_region(Region::CH1, "CH1").ok());
                // Try to detect the best H/CH2
                if self.regions.contains_key("H") && self.regions.contains_key("CH2") {
                    seq.extend(self.get_region(Region::H, "H").ok());
                    seq.extend(self.get_region(Region::CH2, "CH2").ok());
                } else if self.regions.contains_key("H-CH2") {
                    seq.extend(self.get_region(Region::H_CH2, "H-CH2").ok());
                } else {
                    seq.extend(self.get_region(Region::H1, "H1").ok());
                    seq.extend(self.get_region(Region::H2, "H2").ok());
                    seq.extend(self.get_region(Region::H3, "H3").ok());
                    seq.extend(self.get_region(Region::H4, "H4").ok());
                    seq.extend(self.get_region(Region::CH2, "CH2").ok());
                }
                let mut secretory_added = false;

                for (region, with_chs) in [
                    (Region::CH3, Region::CH3_CHS),
                    (Region::CH4, Region::CH4_CHS),
                    (Region::CH5, Region::CH5_CHS),
                    (Region::CH6, Region::CH6_CHS),
                    (Region::CH7, Region::CH7_CHS),
                    (Region::CH8, Region::CH8_CHS),
                    (Region::CH9, Region::CH9_CHS),
                ] {
                    if self.regions.contains_key(&region.to_string())
                        && self.regions.contains_key("CHS")
                    {
                        // If this region is stored separately and the CHS is available
                        seq.extend(self.get_region(region, &region.to_string()).ok());
                    } else if self.regions.contains_key(&with_chs.to_string()) {
                        // If this region is stored together with the CHS
                        seq.extend(self.get_region(with_chs, &with_chs.to_string()).ok());
                        secretory_added = true;
                    } else {
                        // If this region is stored separately but some later region is combined with CHS
                        seq.extend(self.get_region(region, &region.to_string()).ok());
                    }
                }
                if !secretory_added {
                    seq.extend(self.get_region(Region::CHS, "CHS").ok());
                }
                // possibly_add(Region::M, "M")?; // TODO: Potentially find a way to also allow the membrane bound version to be stored
                // possibly_add(Region::M1, "M1")?;
                // possibly_add(Region::M2, "M2")?;

                // Otherwise assume light chain
                if seq.is_empty() {
                    seq.extend(self.get_region(Region::CL, "CL").ok());
                }
                if seq.is_empty() {
                    seq.extend(self.get_region(Region::CL, "C-REGION").ok());
                }

                if seq.is_empty() {
                    return Err("Empty C sequence".to_string());
                }
                seq
            }
            _ => Vec::new(),
        };

        Ok((regions, additional_annotations))
    }

    fn get_region(
        &self,
        region: Region,
        key: &str,
    ) -> Result<(Region, (Vec<AminoAcid>, Location, String)), String> {
        self.regions
            .get(key)
            .ok_or(format!("Could not find {key}"))
            .and_then(|region| {
                region
                    .found_seq
                    .as_ref()
                    .map(|seq| {
                        let mut final_seq = region
                            .splice_aa
                            .map(|aa| vec![aa])
                            .filter(|_| region.shift != 2)
                            .unwrap_or_default();
                        final_seq.extend(seq.1 .0.clone());
                        (final_seq, region.location.clone(), seq.0.clone())
                    })
                    .map_err(|e| e.to_owned())
            })
            .map(|res| (region, res))
    }
}

impl Display for IMGTGene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{}\t{}\t{}", self.key, self.location, self.allele)?;
        for region in self.regions.values().sorted_by_key(|reg| &reg.key) {
            writeln!(f, "  R {region}")?;
        }
        Ok(())
    }
}
