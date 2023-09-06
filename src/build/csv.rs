use std::{
    fs::File,
    io::{BufRead, BufReader},
};

pub fn parse_csv(path: &str) -> Result<Vec<Vec<String>>, String> {
    let reader = BufReader::new(File::open(path).map_err(|e| e.to_string())?);
    let mut table = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?.trim_end().to_string();
        if line.is_empty() {
            continue;
        }
        let mut enclosed = None;
        let mut row = Vec::new();
        let mut cell = String::new();
        for ch in line.chars() {
            match (ch, enclosed) {
                ('\"', None) | ('\'', None) => enclosed = Some(ch),
                (c, Some(e)) if c == e => enclosed = None,
                (',', None) => {
                    row.push(cell.trim().to_string());
                    cell.clear();
                }
                _ => cell.push(ch),
            }
        }
        row.push(cell.trim().to_string());
        table.push(row);
    }
    Ok(table)
}
