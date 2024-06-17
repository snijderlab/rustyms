use std::{
    num::{IntErrorKind, ParseFloatError, ParseIntError},
    ops::{Bound, Range, RangeBounds},
    path::Path,
    str::FromStr,
};

use itertools::Itertools;

pub trait ResultExtensions<T, E> {
    /// # Errors
    /// If any of the errors contained within has an error.
    #[allow(dead_code)]
    fn flat_err(self) -> Result<T, E>;
}

impl<T, E> ResultExtensions<T, E> for Result<T, Result<T, E>> {
    fn flat_err(self) -> Result<T, E> {
        match self {
            Ok(o) => Ok(o),
            Err(r) => r,
        }
    }
}

impl<T, E> ResultExtensions<T, E> for Result<Result<T, E>, E> {
    fn flat_err(self) -> Result<T, E> {
        match self {
            Ok(o) => o,
            Err(r) => Err(r),
        }
    }
}

pub trait InvertResult<T, E> {
    /// # Errors
    /// If any of the errors contained within has an error.
    #[allow(dead_code)]
    fn invert(self) -> Result<Option<T>, E>;
}

impl<T, E> InvertResult<T, E> for Option<Result<T, E>> {
    fn invert(self) -> Result<Option<T>, E> {
        self.map_or_else(|| Ok(None), |o| o.map(|v| Some(v)))
    }
}
impl<T, E> InvertResult<T, E> for Option<Result<Option<T>, E>> {
    fn invert(self) -> Result<Option<T>, E> {
        self.map_or_else(|| Ok(None), |o| o)
    }
}

pub trait RangeExtension
where
    Self: Sized,
{
    fn add_start(&self, amount: isize) -> Option<Self>;
    fn add_end(&self, amount: isize) -> Option<Self>;
    fn start_index(&self) -> usize;
    fn end_index(&self) -> usize;
}

impl RangeExtension for Range<usize> {
    fn add_start(&self, amount: isize) -> Option<Self> {
        let new_start = usize::try_from(isize::try_from(self.start).ok()? + amount).ok()?;
        (new_start <= self.end).then_some(Self {
            start: new_start,
            end: self.end,
        })
    }
    fn add_end(&self, amount: isize) -> Option<Self> {
        let new_end = usize::try_from(isize::try_from(self.end).ok()? + amount).ok()?;
        (self.start <= new_end).then_some(Self {
            start: self.start,
            end: new_end,
        })
    }
    fn start_index(&self) -> usize {
        match self.start_bound() {
            std::ops::Bound::Unbounded => 0,
            std::ops::Bound::Included(s) => *s,
            std::ops::Bound::Excluded(s) => s + 1,
        }
    }
    fn end_index(&self) -> usize {
        match self.end_bound() {
            std::ops::Bound::Unbounded => 0,
            std::ops::Bound::Included(s) => *s,
            std::ops::Bound::Excluded(s) => s - 1,
        }
    }
}

/// # Errors
/// If the name cannot be recognised or a number is not valid.
#[allow(dead_code)]
pub fn parse_named_counter<T: Clone>(
    value: &str,
    names: &[(String, T)],
    allow_negative: bool,
) -> Result<Vec<(T, isize)>, String> {
    let mut index = 0;
    let mut output = Vec::new();
    while index < value.len() {
        if value[index..].starts_with(' ') {
            index += 1;
        } else {
            let mut found = false;
            for name in names {
                if value[index..].starts_with(&name.0) {
                    index += name.0.len();
                    let num = &value[index..]
                        .chars()
                        .skip_while(char::is_ascii_whitespace)
                        .take_while(|c| c.is_ascii_digit() || (allow_negative && *c == '-'))
                        .collect::<String>()
                        .trim()
                        .to_string();
                    if num.is_empty() {
                        output.push((name.1.clone(), 1));
                    } else {
                        output.push((
                            name.1.clone(),
                            num.parse()
                                .map_err(|_| format!("Not a valid number '{num}'"))?,
                        ));
                        index += num.len()
                            + value[index..]
                                .chars()
                                .take_while(char::is_ascii_whitespace)
                                .count();
                    }
                    found = true;
                    break; // Names loop
                }
            }
            if !found {
                return Err(format!("Name not recognised {}", &value[index..]));
            }
        }
    }
    Ok(output)
}

/// Split a string into chunks of text separated by whitespace with the offset before each chunk returned for nice error generation.
pub fn split_ascii_whitespace(input: &str) -> Vec<(usize, &str)> {
    let mut index = input.chars().take_while(char::is_ascii_whitespace).count();
    let mut chunks = Vec::new();
    while index < input.len() {
        let chunk_len = input[index..]
            .chars()
            .take_while(|c| !c.is_ascii_whitespace())
            .count();
        chunks.push((index, &input[index..index + chunk_len]));
        index += chunk_len;
        index += input[index..]
            .chars()
            .take_while(char::is_ascii_whitespace)
            .count();
    }
    chunks
}

/// Helper function to check extensions in filenames
pub fn check_extension(filename: impl AsRef<Path>, extension: impl AsRef<Path>) -> bool {
    filename
        .as_ref()
        .extension()
        .map_or(false, |ext| ext.eq_ignore_ascii_case(extension.as_ref()))
}

#[allow(dead_code)]
/// Get the index of the next copy of the given char
pub fn next_char(chars: &[u8], start: usize, char: u8) -> Option<usize> {
    for (i, ch) in chars[start..].iter().enumerate() {
        if *ch == char {
            return Some(start + i);
        }
    }
    None
}

#[allow(dead_code)]
/// Find the enclosed text by the given symbols, assumes a single open is already read just before the start
pub fn end_of_enclosure(text: &str, start: usize, open: u8, close: u8) -> Option<usize> {
    let mut state = 1;
    for (i, ch) in text[start..].as_bytes().iter().enumerate() {
        if !text.is_char_boundary(i) {
            continue;
        }
        if i + 1 < text.len() && !text.is_char_boundary(i + 1) {
            continue;
        }
        if *ch == open {
            state += 1;
        } else if *ch == close {
            state -= 1;
            if state == 0 {
                return Some(start + i);
            }
        }
    }
    None
}

#[allow(dead_code)]
/// Find the enclosed text by the given symbols, assumes a single open is already read just before the start.
/// This also takes brackets '[]' into account and these take precedence over the enclosure searched for.
pub fn end_of_enclosure_with_brackets(
    text: &str,
    start: usize,
    open: u8,
    close: u8,
) -> Option<usize> {
    let mut state = 1;
    let mut index = start;
    while index < text.len() {
        if !text.is_char_boundary(index) {
            continue;
        }
        if index + 1 < text.len() && !text.is_char_boundary(index + 1) {
            continue;
        }
        let ch = text.as_bytes()[index];
        if ch == b'[' {
            index = end_of_enclosure(text, index + 1, b'[', b']')?;
        }
        if ch == open {
            state += 1;
        } else if ch == close {
            state -= 1;
            if state == 0 {
                return Some(index);
            }
        }
        index += 1;
    }
    None
}

#[allow(dead_code)]
/// Get the next number, returns length in bytes and the number.
/// # Panics
/// If the text is not valid UTF-8.
/// # Errors
/// Returns none if the number is too big to fit in a `isize`.
pub fn next_num(chars: &[u8], mut start: usize, allow_only_sign: bool) -> Option<(usize, isize)> {
    let mut sign = 1;
    let mut sign_set = false;
    if chars.get(start) == Some(&b'-') {
        sign = -1;
        start += 1;
        sign_set = true;
    } else if chars.get(start) == Some(&b'+') {
        start += 1;
        sign_set = true;
    }
    let len = chars[start..]
        .iter()
        .take_while(|c| c.is_ascii_digit())
        .count();
    if len == 0 {
        if allow_only_sign && sign_set {
            Some((1, sign))
        } else {
            None
        }
    } else {
        let num: isize = std::str::from_utf8(&chars[start..start + len])
            .unwrap()
            .parse()
            .ok()?;
        Some((usize::from(sign_set) + len, sign * num))
    }
}

/// A number of characters, used as length or index
pub type Characters = usize;

#[allow(dead_code)]
/// Get the next number starting at the character range given, returns length in characters and the number.
/// # Errors
/// Returns none if the number is too big to fit in a `isize`.
pub fn next_number<const ALLOW_SIGN: bool, Number: FromStr>(
    line: &str,
    range: impl RangeBounds<Characters>,
) -> Option<(Characters, bool, Result<Number, Number::Err>)> {
    let start = match range.start_bound() {
        Bound::Unbounded => 0,
        Bound::Excluded(n) => n + 1,
        Bound::Included(n) => *n,
    };
    let end = match range.end_bound() {
        Bound::Unbounded => line.chars().count(),
        Bound::Excluded(n) => n - 1,
        Bound::Included(n) => *n,
    };
    let mut positive = true;
    let mut sign_set = false;
    let mut start_index = 0;
    let mut chars = line.char_indices().skip(start).peekable();
    if ALLOW_SIGN {
        match chars.peek() {
            Some((_, '-')) => {
                positive = false;
                sign_set = true;
            }
            Some((_, '+')) => {
                sign_set = true;
            }
            _ => (),
        }
        if sign_set {
            let _ = chars.next();
        }
    }
    if let Some(index) = chars.peek().map(|(index, _)| index) {
        start_index = *index;
    }

    let mut consumed = usize::from(sign_set);
    chars
        .take_while(|(_, c)| {
            if c.is_ascii_digit() {
                consumed += 1;
                consumed < end - start
            } else {
                false
            }
        })
        .last()
        .map(|(end_index, c)| {
            (
                consumed,
                positive,
                line[start_index..end_index + c.len_utf8()].parse::<Number>(),
            )
        })
}

/// Get a canonicalised u64 for f64 to be able to hash f64, based on the `ordered_float` crate (MIT license)
pub fn f64_bits(value: f64) -> u64 {
    if value.is_nan() {
        0x7ff8_0000_0000_0000_u64 // CANONICAL_NAN_BITS
    } else {
        (value + 0.0).to_bits() // The +0.0 is to guarantee even handling of negative and positive zero
    }
}

/// Implement a binary operator for all ref cases after the implementation for the ref-ref case (assumes deref operator works)
macro_rules! impl_binop_ref_cases {
    (impl $imp:ident, $method:ident for $t:ty, $u:ty, $o:ty) => {
        impl<'a> $imp<$u> for &'a $t {
            type Output = $o;

            #[inline]
            fn $method(self, other: $u) -> $o {
                $imp::$method(self, &other)
            }
        }

        impl<'a> $imp<&'a $u> for $t {
            type Output = $o;

            #[inline]
            fn $method(self, other: &'a $u) -> $o {
                $imp::$method(&self, other)
            }
        }

        impl $imp<$u> for $t {
            type Output = $o;

            #[inline]
            fn $method(self, other: $u) -> $o {
                $imp::$method(&self, &other)
            }
        }
    };
}

/// To be used as `The xx number ` + the explanation from here (does not have a dot).
pub fn explain_number_error(error: &ParseIntError) -> &'static str {
    match error.kind() {
        IntErrorKind::Empty => "is empty",
        IntErrorKind::InvalidDigit => "contains an invalid character",
        IntErrorKind::NegOverflow => "is too small to fit in the internal representation",
        IntErrorKind::PosOverflow => "is too big to fit in the internal representation",
        IntErrorKind::Zero => "is zero, which is not allowed here",
        _ => "is not a valid number",
    }
}
