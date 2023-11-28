use std::{
    fmt::Debug,
    fs::File,
    io::{BufRead, BufReader},
    ops::Range,
};

use flate2::read::GzDecoder;

use crate::helper_functions::check_extension;

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CsvLine {
    line_index: usize,
    line: String,
    fields: Vec<(String, Range<usize>)>,
}

impl CsvLine {
    pub fn line_index(&self) -> usize {
        self.line_index
    }
    pub fn line(&self) -> &str {
        &self.line
    }
    pub fn number_of_columns(&self) -> usize {
        self.fields.len()
    }
    pub fn column_context(&self, column: usize) -> crate::error::Context {
        crate::error::Context::line(
            self.line_index,
            self.line.clone(),
            self.fields[column].1.start,
            self.fields[column].1.len(),
        )
    }
    pub fn range_context(&self, range: std::ops::Range<usize>) -> crate::error::Context {
        crate::error::Context::line(
            self.line_index(),
            self.line.clone(),
            range.start,
            range.len(),
        )
    }
    pub fn full_context(&self) -> crate::error::Context {
        crate::error::Context::full_line(self.line_index, self.line.clone())
    }
    pub fn range(&self, index: usize) -> &Range<usize> {
        &self.fields[index].1
    }
    pub fn column(&self, name: &str) -> Option<&str> {
        self.fields
            .iter()
            .find(|f| f.0 == name)
            .map(|f| &self.line[f.1.clone()])
    }
    pub fn optional_column(&self, name: Option<&str>) -> Option<Option<&str>> {
        name.map(|name| self.column(name))
    }
}

impl std::ops::Index<usize> for CsvLine {
    type Output = str;
    fn index(&self, index: usize) -> &str {
        &self.line[self.fields[index].1.clone()]
    }
}

/// Parse a CSV file into an iterator with the parsed lines.
/// # Errors
/// If the file cannot be opened it returns `Err` with the error.
/// If any single line cannot be read it returns an error for that line.
pub fn parse_csv(
    path: impl AsRef<std::path::Path>,
    separator: u8,
    provided_header: Option<Vec<String>>,
) -> Result<Box<dyn Iterator<Item = Result<CsvLine, String>>>, String> {
    let file = File::open(path.as_ref()).map_err(|e| e.to_string())?;
    if check_extension(path, "gz") {
        Ok(Box::new(parse_csv_raw(
            BufReader::new(GzDecoder::new(file)),
            separator,
            provided_header,
        )?))
    } else {
        Ok(Box::new(parse_csv_raw(
            BufReader::new(file),
            separator,
            provided_header,
        )?))
    }
}

/// Parse a CSV file from a raw `BufReader`
pub fn parse_csv_raw<T: std::io::Read>(
    reader: BufReader<T>,
    separator: u8,
    provided_header: Option<Vec<String>>,
) -> Result<CsvLineIter<T>, String> {
    let mut lines = reader.lines().enumerate();
    let header = if let Some(header) = provided_header {
        header
    } else {
        let (_, header) = lines
            .next()
            .ok_or("Empty CSV file, but it should contain a header")?;
        let header_line = header.map_err(|err| format!("Could not read header line: {err}"))?;
        csv_separate(&header_line, 0, separator)?
            .into_iter()
            .map(|r| header_line[r].to_lowercase())
            .collect()
    };

    Ok(CsvLineIter {
        lines,
        header,
        separator,
    })
}

/// An iterator returning CSV lines
pub struct CsvLineIter<T: std::io::Read> {
    lines: std::iter::Enumerate<std::io::Lines<BufReader<T>>>,
    header: Vec<String>,
    separator: u8,
}

impl<T: std::io::Read> Iterator for CsvLineIter<T> {
    type Item = Result<CsvLine, String>;
    fn next(&mut self) -> Option<Self::Item> {
        self.lines.next().map(|(line_index, line)| {
            let line = if let Ok(line) = line {
                line.trim_end().to_string()
            } else {
                return Err(format!("Could no read line {line_index}"));
            };
            csv_separate(&line, line_index, self.separator).map(|row| CsvLine {
                line_index,
                line,
                fields: self.header.iter().cloned().zip(row).collect(),
            })
        })
    }
}

fn csv_separate(line: &str, line_index: usize, separator: u8) -> Result<Vec<Range<usize>>, String> {
    if line.is_empty() {
        return Err("Empty line".to_string());
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
            (sep, None, Some(s)) if sep == separator => {
                row.push(s..index);
                start = None;
                was_enclosed = false;
            }
            (sep, None, None) if sep == separator => {
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
    Ok(row)
}

impl std::fmt::Display for CsvLine {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut line = self.line.clone();
        let mut offset = 0;
        let mut display_field = |range: &Range<usize>| match range.len() {
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
        for (name, field) in &self.fields {
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
