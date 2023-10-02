use std::{
    fmt::Debug,
    fs::File,
    io::{BufRead, BufReader},
};

#[derive(Debug, PartialEq, Eq, Clone)]
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
    parse_csv_raw(reader)
}

pub fn parse_csv_raw<T: std::io::Read>(reader: BufReader<T>) -> Result<Vec<CsvLine>, String> {
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
                    if !was_enclosed {
                        // ignore any comma directly after an enclosed field
                        row.push(index..index);
                        start = None;
                    }
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

impl std::fmt::Display for CsvLine {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut line = self.line.clone();
        let mut offset = 0;
        let mut display_field = |range: &std::ops::Range<usize>| match range.len() {
            0 => {
                line.insert(range.start + offset, '⌧'); // Insert an empty icon to be able to indicate a field of length 0
                offset += '⌧'.len_utf8();
                "ô".to_string()
            }
            1 => "^".to_string(),
            n => {
                let full = format!("{}-{}", range.start, range.end);
                if full.len() <= n - 2 {
                    format!("└{:·^w$}┘", full, w = n - 2)
                } else {
                    format!("└{}┘", "·".repeat(n - 2))
                }
            }
        };
        let mut fields = String::new();
        let mut index = 0;
        for field in &self.fields {
            if field.start >= index {
                fields = format!(
                    "{fields}{}{}",
                    " ".repeat(field.start - index),
                    display_field(field)
                );
                index = field.end;
            } else {
                fields = format!(
                    "{fields}\n{}{}",
                    " ".repeat(field.start),
                    display_field(field)
                );
            }
        }
        write!(
            f,
            "{}: {}\n{}{}",
            self.line_index + 1,
            line,
            " ".repeat((self.line_index + 1).to_string().len() + 2),
            fields
        )
    }
}
