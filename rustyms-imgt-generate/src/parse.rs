use crate::{complement, imgt_gene::IMGTGene, translate};
use crate::{
    shared::Species,
    structs::{AASequence, DataItem, Location, Region},
};
use itertools::Itertools;
use rustyms::AminoAcid;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

pub fn parse_dat<T: std::io::Read>(
    reader: BufReader<T>,
) -> impl Iterator<Item = Result<DataItem, String>> {
    reader
        .lines()
        .batching(|f| {
            let mut data = PreDataItem::default();
            for line in f.filter_map(|i| i.ok()) {
                match &line[..2] {
                    "//" => return Some(data),
                    "ID" => data.id = line,
                    "KW" => data.kw.extend(
                        line[5..]
                            .split(';')
                            .map(|s| s.trim().to_string())
                            .filter(|s| !s.is_empty()),
                    ),
                    "FH" if line.starts_with("FH   Key") => {
                        data.ft_key_width = line.find("Location").expect("Incorrect FH line") - 5
                    }
                    "FT" => data.ft.push(line),
                    "OS" if data.os.is_none() => {
                        data.os = Species::from_imgt(line[5..].trim()).unwrap_or_else(|()| {
                            println!("Not a species name: `{line}`");
                            None
                        })
                    }
                    "  " => {
                        data.sq.extend(line.chars().filter(|c| {
                            *c == 'c' || *c == 'a' || *c == 't' || *c == 'g' || *c == 'n'
                        }))
                    }
                    _ => (),
                }
            }
            None
        })
        .filter(|pre| {
            pre.kw.contains(&"immunoglobulin (IG)".to_string())
                && (pre.kw.contains(&"functional".to_string())
                    || pre.kw.contains(&"germline".to_string()))
                && pre.os.is_some()
        })
        .map(DataItem::new)
}

impl DataItem {
    fn new(data: PreDataItem) -> Result<Self, String> {
        let mut result = Self {
            id: data.id[5..].split(';').next().unwrap().to_string(),
            species: data.os.ok_or("No species found")?,
            sequence: data.sq,
            genes: Vec::new(),
            regions: Vec::new(),
        };
        let mut current: Option<Region> = None;
        let mut sequence = false;
        for line in data.ft {
            let line = &line[5..];
            if !line.starts_with(' ') || current.is_none() {
                if let Some(region) = current {
                    result.add_region(region);
                }
                let (key, location) = (&line[..data.ft_key_width], &line[data.ft_key_width..]);
                if location.contains("join") {
                    return Err("Location is a joined region".to_string());
                }
                if location.contains('^') {
                    return Err("Location is a ^ region".to_string());
                }
                let location = location
                    .trim()
                    .parse()
                    .unwrap_or_else(|_| panic!("`{}` not a valid location", location));
                current = Some(Region {
                    acc: result.id.clone(),
                    key: key.trim().to_string(),
                    location,
                    reported_seq: String::new(),
                    found_seq: Err("Not loaded".to_string()),
                    allele: String::new(),
                    functional: false,
                    partial: false,
                    shift: 0,
                    splice_aa: None,
                });
                continue;
            }
            if let Some(current) = &mut current {
                let trimmed = line.trim();
                let lowercase = trimmed.to_lowercase();
                if sequence {
                    current.reported_seq = trimmed.trim_end_matches('\"').to_string();
                    if trimmed.ends_with('\"') {
                        sequence = false;
                    }
                } else if let Some(tail) = trimmed.strip_prefix("/translation=\"") {
                    current.reported_seq = tail.trim_end_matches('\"').to_string();
                    if !trimmed.ends_with('\"') {
                        sequence = true;
                    }
                } else if let Some(tail) = trimmed.strip_prefix("/IMGT_allele=\"") {
                    current.allele = tail.trim_end_matches('\"').to_string();
                } else if let Some(tail) = trimmed.strip_prefix("/codon_start=") {
                    current.shift = tail
                        .parse::<usize>()
                        .map_err(|_| format!("Not a valid codon_start: '{tail}'"))?
                        - 1;
                } else if let Some(tail) = trimmed.strip_prefix("/splice-expectedcodon=") {
                    if let Some(i) = tail.find(']') {
                        current.splice_aa = AminoAcid::try_from(tail.as_bytes()[i - 1]).ok();
                    }
                } else if lowercase.starts_with("/functional")
                    || lowercase.starts_with("/note=\"functional\"")
                    || lowercase.starts_with("/imgt_note=\"functional\"")
                {
                    current.functional = true;
                } else if trimmed.starts_with("/partial") {
                    current.partial = true;
                }
            }
        }
        if let Some(region) = current {
            result.add_region(region);
        }

        Ok(result)
    }

    fn add_region(&mut self, mut region: Region) {
        // Get the actual sequence
        region.found_seq = self.get_sequence(&region.location, region.shift);

        // Determine if what this region is and if is warrants keeping
        if ["V-GENE", "C-GENE", "J-GENE"].contains(&region.key.as_str()) // , "D-GENE"
            && region.functional
            && !region.partial
            && region.allele.starts_with("IG")
        {
            self.genes.push(IMGTGene {
                acc: region.acc,
                key: region.key,
                location: region.location,
                allele: region.allele,
                regions: HashMap::new(),
            });
        } else if [
            "FR1-IMGT",
            "FR2-IMGT",
            "FR3-IMGT",
            "FR4-IMGT",
            "CDR1-IMGT",
            "CDR2-IMGT",
            "CDR3-IMGT",
            "1st-CYS",
            "2nd-CYS",
            "CONSERVED-TRP",
            "J-REGION",
            // "J-TRP",
            // "J-PHE",
            //"J-MOTIF",
            "CH1",
            "CH2",
            "H-CH2",
            "CH3",
            "CH3-CHS",
            "CH4",
            "CH4-CHS",
            "CH5",
            "CH5-CHS",
            "CH6",
            "CH6-CHS",
            "CH7",
            "CH7-CHS",
            "CH8",
            "CH8-CHS",
            "CH9",
            "CH9-CHS",
            "CHS",
            "CL",
            "C-REGION",
            "H", //"D-REGION",
            "H1",
            "H2",
            "H3",
            "H4",
            "M",
            "M1",
            "M2",
        ]
        .contains(&region.key.as_str())
        {
            if let Some(gene) = self
                .genes
                .iter_mut()
                .find(|g| g.location.contains(&region.location))
            {
                gene.regions.insert(region.key.clone(), region);
            } else {
                self.regions.push(region)
            }
        }
    }

    fn get_sequence(&self, slice: &Location, shift: usize) -> Result<(String, AASequence), String> {
        let (inner_shift, shift) = if shift == 2 { (1, 0) } else { (0, shift) };

        translate(
            &match slice {
                Location::Normal(range) => {
                    if *range.start() < inner_shift {
                        return Err("Shift outside of range".to_string());
                    }
                    self.sequence
                        .get(range.start() - inner_shift..=*range.end())
                        .ok_or("Normal outside of range")?
                        .to_string()
                }
                Location::SingleNormal(index) => char::from(
                    *self
                        .sequence
                        .as_bytes()
                        .get(*index)
                        .ok_or("Single normal outside of range")?,
                )
                .to_string(),
                Location::Complement(range) => complement(
                    self.sequence
                        .get(*range.start()..=*range.end() + inner_shift)
                        .ok_or("Complement outside of range")?
                        .to_string(),
                ),
                Location::SingleComplement(index) => complement(
                    char::from(
                        *self
                            .sequence
                            .as_bytes()
                            .get(*index)
                            .ok_or("Single complement outside of range")?,
                    )
                    .to_string(),
                ),
            }[shift..],
        )
        .map(|(s, v)| (s.to_owned(), AASequence(v)))
    }
}

#[derive(Default, Debug)]
struct PreDataItem {
    id: String,
    kw: Vec<String>,
    ft_key_width: usize,
    ft: Vec<String>,
    os: Option<Species>,
    sq: String,
}
