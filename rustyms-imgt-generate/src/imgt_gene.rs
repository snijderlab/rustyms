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
        let get = |key| -> Result<(Vec<AminoAcid>, Location, String), String> {
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
        };
        let mut additional_annotations = Vec::new();
        let regions = if self.key == "V-GENE" {
            vec![
                (Region::FR1, get("FR1-IMGT")?),
                (Region::CDR1, get("CDR1-IMGT")?),
                (Region::FR2, get("FR2-IMGT")?),
                (Region::CDR2, get("CDR2-IMGT")?),
                (Region::FR3, get("FR3-IMGT")?),
                (Region::CDR3, get("CDR3-IMGT")?),
            ]
        } else if self.key == "C-GENE" {
            // if self.allele == "IGHA1*01" {
            //     dbg!(&self);
            // }
            let mut seq = Vec::new();
            let mut possibly_add = |region, key: &str, only_if_empty: bool| -> Result<(), String> {
                if self.regions.contains_key(key)
                    && ((only_if_empty && seq.is_empty()) || !only_if_empty)
                {
                    seq.push((
                        region,
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
                            })?,
                    ))
                }
                Ok(())
            };

            // Heavy chain
            possibly_add(Region::CH1, "CH1", false)?;
            // Try to detect the best H/CH2
            if self.regions.contains_key("H") && self.regions.contains_key("CH2") {
                possibly_add(Region::H, "H", false)?;
                possibly_add(Region::CH2, "CH2", false)?;
            } else if self.regions.contains_key("H-CH2") {
                possibly_add(Region::H_CH2, "H-CH2", false)?;
            } else {
                possibly_add(Region::H1, "H1", false)?;
                possibly_add(Region::H2, "H2", false)?;
                possibly_add(Region::H3, "H3", false)?;
                possibly_add(Region::H4, "H4", false)?;
                possibly_add(Region::CH2, "CH2", false)?;
            }
            let mut secretory = false;
            if self.regions.contains_key("CH3") && self.regions.contains_key("CHS") {
                possibly_add(Region::CH3, "CH3", false)?;
            } else if self.regions.contains_key("CH3-CHS") {
                possibly_add(Region::CH3_CHS, "CH3-CHS", false)?;
                secretory = true;
            }
            if self.regions.contains_key("CH4") && self.regions.contains_key("CHS") {
                possibly_add(Region::CH4, "CH4", false)?;
            } else if self.regions.contains_key("CH4-CHS") {
                possibly_add(Region::CH4_CHS, "CH4-CHS", false)?;
                secretory = true;
            }
            if self.regions.contains_key("CH5") && self.regions.contains_key("CHS") {
                possibly_add(Region::CH5, "CH5", false)?;
            } else if self.regions.contains_key("CH5-CHS") {
                possibly_add(Region::CH5_CHS, "CH5-CHS", false)?;
                secretory = true;
            }
            if self.regions.contains_key("CH6") && self.regions.contains_key("CHS") {
                possibly_add(Region::CH6, "CH6", false)?;
            } else if self.regions.contains_key("CH6-CHS") {
                possibly_add(Region::CH6_CHS, "CH6-CHS", false)?;
                secretory = true;
            }
            if self.regions.contains_key("CH7") && self.regions.contains_key("CHS") {
                possibly_add(Region::CH7, "CH7", false)?;
            } else if self.regions.contains_key("CH7-CHS") {
                possibly_add(Region::CH7_CHS, "CH7-CHS", false)?;
                secretory = true;
            }
            if self.regions.contains_key("CH8") && self.regions.contains_key("CHS") {
                possibly_add(Region::CH8, "CH8", false)?;
            } else if self.regions.contains_key("CH8-CHS") {
                possibly_add(Region::CH8_CHS, "CH8-CHS", false)?;
                secretory = true;
            }
            if self.regions.contains_key("CH9") && self.regions.contains_key("CHS") {
                possibly_add(Region::CH9, "CH9", false)?;
            } else if self.regions.contains_key("CH9-CHS") {
                possibly_add(Region::CH9_CHS, "CH9-CHS", false)?;
                secretory = true;
            }
            if !secretory {
                possibly_add(Region::CHS, "CHS", false)?;
            }
            // possibly_add(Region::M, "M")?; // TODO: Figure out if support for membrane bound is needed, if so provide a way to switch between the two versions
            // possibly_add(Region::M1, "M1")?;
            // possibly_add(Region::M2, "M2")?;

            // Otherwise assume light chain
            possibly_add(Region::CL, "CL", true)?;
            possibly_add(Region::CL, "C-REGION", true)?;

            if seq.is_empty() {
                return Err("Empty C sequence".to_string());
            }
            seq
        } else if self.key == "J-GENE" {
            // if self.regions.contains_key("FR4-IMGT") {
            //     let fr4 = get("FR4-IMGT")?;
            //     let j = get("J-REGION")?;
            //     let cdr3_len = j.0.len() - fr4.0.len();
            //     let j = fix_j(j, cdr3_len);
            //     additional_annotations.extend(j.1);
            //     j.0
            // } else if self.regions.contains_key("J-MOTIF") {
            //     let motif = get("J-MOTIF")?;
            //     let j = get("J-REGION")?;
            //     let loc =
            //         j.1.get_aa_loc(&motif.1)
            //             .ok_or("J-MOTIF does not fall into J-REGION")?;
            //     let j = fix_j(j, *loc.start());
            //     additional_annotations.extend(j.1);
            //     j.0
            // } else {
            // if self.allele == "IGKJ1*01" {
            //     dbg!(&self);
            // }
            let j = get("J-REGION")?;
            let motif = j.0.iter().tuple_windows().position(|(a, b, _, d)| {
                (*a == AminoAcid::W || *a == AminoAcid::F)
                    && *b == AminoAcid::G
                    && *d == AminoAcid::G
            });
            if let Some(motif_start) = motif {
                let j = fix_j(j, motif_start);
                additional_annotations.extend(j.1);
                j.0
            } else {
                vec![(Region::FR4, j)] // TODO: not fully correct right, has some CDR3 as well, and has quite some conserved residues
            }
            // }
        } else if self.key == "D-GENE" {
            vec![(Region::CDR3, get("D-REGION")?)]
        } else {
            Vec::new()
        };
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
