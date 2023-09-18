/// The file format for any peaks format, determining the existence and location of all possible columns
pub struct PeaksFormat {
    fraction: Option<usize>,
    source_file: Option<usize>,
    feature: Option<usize>,
    scan: usize,
    peptide: usize,
    tag_length: usize,
    de_novo_score: Option<usize>,
    alc: usize,
    length: usize,
    mz: usize,
    z: usize,
    rt: usize,
    predicted_rt: Option<usize>,
    area: usize,
    mass: usize,
    ppm: usize,
    ptm: usize,
    local_confidence: usize,
    tag: usize,
    mode: usize,
    accession: Option<usize>,
    format: PeaksVersion,
}

/// All possible peaks versions
pub enum PeaksVersion {
    /// A custom version of a PEAKS file forma
    Custom,
    /// An older version of a PEAKS export
    Old,
    /// Version X of PEAKS export (made for build 31 january 2019)
    X,
    /// Version X+ of PEAKS export (made for build 20 november 2019)
    Xplus,
    /// Version Ab of PEAKS export
    Ab,
}

/// An older version of a PEAKS export
pub const OLD: PeaksFormat = PeaksFormat {
    scan: 0,
    peptide: 1,
    tag_length: 2,
    alc: 3,
    length: 4,
    mz: 5,
    z: 6,
    rt: 7,
    area: 8,
    mass: 9,
    ppm: 10,
    ptm: 11,
    local_confidence: 12,
    tag: 13,
    mode: 14,
    fraction: None,
    source_file: None,
    feature: None,
    de_novo_score: None,
    predicted_rt: None,
    accession: None,
    format: PeaksVersion::Old,
};

/// Version X of PEAKS export (made for build 31 january 2019)
pub const X: PeaksFormat = PeaksFormat {
    fraction: Some(0),
    source_file: Some(1),
    feature: Some(2),
    peptide: 3,
    scan: 4,
    tag_length: 5,
    alc: 6,
    length: 7,
    mz: 8,
    z: 9,
    rt: 10,
    area: 11,
    mass: 12,
    ppm: 13,
    ptm: 14,
    local_confidence: 15,
    tag: 16,
    mode: 17,
    de_novo_score: None,
    predicted_rt: None,
    accession: None,
    format: PeaksVersion::X,
};

/// Version X+ of PEAKS export (made for build 20 november 2019)
pub const XPLUS: PeaksFormat = PeaksFormat {
    fraction: Some(0),
    source_file: Some(1),
    feature: Some(2),
    peptide: 3,
    scan: 4,
    tag_length: 5,
    de_novo_score: Some(6),
    alc: 7,
    length: 8,
    mz: 9,
    z: 10,
    rt: 11,
    predicted_rt: Some(12),
    area: 13,
    mass: 14,
    ppm: 15,
    ptm: 16,
    local_confidence: 17,
    tag: 18,
    mode: 19,
    accession: None,
    format: PeaksVersion::Xplus,
};

/// Version Ab of PEAKS export
pub const AB: PeaksFormat = PeaksFormat {
    scan: 0,
    peptide: 1,
    tag_length: 2,
    alc: 3,
    length: 4,
    mz: 5,
    z: 6,
    rt: 7,
    area: 8,
    mass: 9,
    ppm: 10,
    accession: Some(11),
    ptm: 12,
    local_confidence: 13,
    tag: 14,
    mode: 15,
    fraction: None,
    predicted_rt: None,
    de_novo_score: None,
    source_file: None,
    feature: None,
    format: PeaksVersion::Ab,
};

/// A single parsed line of a peaks file
pub struct PeaksData {
    fraction: Option<usize>,
    source_file: Option<String>,
    feature: Option<(usize, usize)>,
    scan: (usize, usize),
    peptide: String,
    tag_length: usize,
    de_novo_score: Option<usize>,
    alc: usize,
    length: usize,
    mz: usize,
    z: usize,
    rt: f64,
    predicted_rt: Option<f64>,
    area: f64,
    mass: f64,
    ppm: f64,
    ptm: String,
    local_confidence: Vec<usize>,
    tag: String,
    mode: String,
    accession: Option<String>,
}

// TODO:
// * hard crashes if you have too little elements per line
// * does not crash if you have too many elements per line
// * no nice errors
// * no conversion to common read type yet
// * build peaks metadata struct
impl PeaksData {
    pub fn parse_read(line: &[String], format: PeaksFormat) -> Result<Self, String> {
        macro_rules! get_column {
            ($line:ident, $column:expr) => {
                match $column {
                    Some(index) => Some($line[index].clone()),
                    None => None,
                }
            };
            ($line:ident, $column:expr, $typ:ty) => {
                match $column {
                    Some(index) => Some($line[index].parse::<$typ>().map_err(|e| e.to_string())?),
                    None => None,
                }
            };
            ($line:ident, $column:expr => $func:ident) => {
                match $column {
                    Some(index) => Some($func(&$line[index])?),
                    None => None,
                }
            };
        }
        let get_frac = |text: &str| -> Result<(usize, usize), String> {
            let (start, end) = text
                .split_once(':')
                .ok_or_else(|| "Invalid scan number".to_string())?;
            Ok((
                start[1..].parse::<usize>().map_err(|e| e.to_string())?,
                end.parse::<usize>().map_err(|e| e.to_string())?,
            ))
        };
        Ok(Self {
            fraction: get_column!(line, format.fraction, usize),
            source_file: get_column!(line, format.source_file),
            feature: get_column!(line, format.feature => get_frac),
            scan: get_frac(&line[format.scan])?,
            peptide: line[format.peptide].clone(),
            tag_length: line[format.tag_length]
                .parse::<usize>()
                .map_err(|e| e.to_string())?,
            de_novo_score: get_column!(line, format.de_novo_score, usize),
            alc: line[format.alc]
                .parse::<usize>()
                .map_err(|e| e.to_string())?,
            length: line[format.length]
                .parse::<usize>()
                .map_err(|e| e.to_string())?,
            mz: line[format.mz]
                .parse::<usize>()
                .map_err(|e| e.to_string())?,
            z: line[format.z].parse::<usize>().map_err(|e| e.to_string())?,
            rt: line[format.rt].parse::<f64>().map_err(|e| e.to_string())?,
            predicted_rt: get_column!(line, format.de_novo_score, f64),
            area: line[format.area]
                .parse::<f64>()
                .map_err(|e| e.to_string())?,
            mass: line[format.mass]
                .parse::<f64>()
                .map_err(|e| e.to_string())?,
            ppm: line[format.ppm].parse::<f64>().map_err(|e| e.to_string())?,
            ptm: line[format.ptm].clone(),
            local_confidence: line[format.local_confidence]
                .split(' ')
                .flat_map(|t| t.parse::<usize>().map_err(|e| e.to_string()))
                .collect(),
            tag: line[format.tag].clone(),
            mode: line[format.mode].clone(),
            accession: get_column!(line, format.accession),
        })
    }
}
