use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
};

#[derive(Debug, Default)]
pub struct OboOntology {
    pub headers: Vec<(String, String)>,
    pub objects: Vec<OboObject>,
}

#[derive(Debug, Default)]
pub struct OboObject {
    pub name: String,
    pub lines: HashMap<String, Vec<String>>,
}

impl OboOntology {
    pub fn from_file(path: &str) -> Result<OboOntology, String> {
        let reader = BufReader::new(File::open(path).map_err(|e| e.to_string())?);
        let mut obo = OboOntology::default();
        let mut recent_obj = None;

        for line in reader.lines() {
            let line = line.map_err(|e| e.to_string())?.trim_end().to_string();
            if line.is_empty() {
                continue;
            } else if line.starts_with('[') && line.ends_with(']') {
                if let Some(obj) = recent_obj {
                    obo.objects.push(obj);
                }
                recent_obj = Some(OboObject::new(&line[1..=line.len() - 2]));
            } else if let Some((name, value)) = line.split_once(':') {
                if let Some(obj) = &mut recent_obj {
                    obj.lines
                        .entry(name.trim().to_string())
                        .or_insert(Vec::new())
                        .push(value.trim().to_string());
                } else {
                    obo.headers.push((name.to_string(), value.to_string()));
                }
            } else {
                return Err(format!("Invalid line in obo file: {line}"));
            }
        }
        if let Some(obj) = recent_obj {
            obo.objects.push(obj);
        }
        Ok(obo)
    }
}

impl OboObject {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            ..Self::default()
        }
    }
}
