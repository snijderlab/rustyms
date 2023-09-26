use std::{
    fs::File,
    io::{BufRead, BufReader},
};

pub struct CsvLine {
    pub line_index: usize,
    pub line: String,
    pub fields: Vec<std::ops::Range<usize>>,
}

impl std::ops::Index<usize> for CsvLine {
    type Output = str;
    fn index(&self, index: usize) -> &str {
        &self.line[self.fields[index].clone()]
    }
}

pub fn parse_csv(path: &str) -> Result<Vec<CsvLine>, String> {
    let reader = BufReader::new(File::open(path).map_err(|e| e.to_string())?);
    let mut table = Vec::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?.trim_end().to_string();
        if line.is_empty() {
            continue;
        }
        let mut enclosed = None;
        let mut was_enclosed = false;
        let mut row = Vec::new();
        let mut start = None;
        for (index, ch) in line.bytes().enumerate() {
            match (ch, enclosed, start) {
                (b'\"' | b'\'', None, None) => {
                    enclosed = Some(ch);
                    start = Some(index + 1);
                }
                (c, Some(e), Some(s)) if c == e => {
                    enclosed = None;
                    row.push(s..index);
                    start = None;
                    was_enclosed = true;
                }
                (b',', None, Some(s)) => {
                    row.push(s..index);
                    start = None;
                    was_enclosed = false;
                }
                (b',', None, None) => {
                    row.push(index..index);
                    start = None;
                    was_enclosed = false;
                }
                (c, _, _) if c.is_ascii_whitespace() => (), // ignore
                (_, _, None) => start = Some(index),
                _ => (),
            }
        }
        if let Some(s) = start {
            row.push(s..line.len());
        } else if !was_enclosed {
            row.push(line.len()..line.len());
        }
        table.push(CsvLine {
            line_index,
            line,
            fields: row,
        });
    }
    Ok(table)
}
