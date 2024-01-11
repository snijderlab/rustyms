use std::path::Path;

use crate::{formula::MolecularFormula, ELEMENT_PARSE_LIST};
pub trait ResultExtensions<T, E> {
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
                        output.push((name.1.clone(), num.parse().unwrap()));
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

// ProForma: [13C2][12C-2]H2N
#[allow(dead_code)]
pub fn parse_molecular_formula_pro_forma(value: &str) -> Result<MolecularFormula, String> {
    let mut index = 0;
    let mut element = None;
    let bytes = value.as_bytes();
    let mut result = MolecularFormula::default();
    while index < value.len() {
        match bytes[index] {
            b'[' => {
                index += 1; // Skip the open square bracket
                let len = bytes
                    .iter()
                    .skip(index)
                    .position(|c| *c == b']')
                    .ok_or(format!(
                        "No closing square bracket for square bracket at index: {index}"
                    ))?;
                let isotope = bytes
                    .iter()
                    .skip(index)
                    .take_while(|c| c.is_ascii_digit())
                    .count();
                let ele = bytes
                    .iter()
                    .skip(index + isotope)
                    .take_while(|c| c.is_ascii_alphabetic())
                    .count();

                for possible in ELEMENT_PARSE_LIST {
                    if value[index + isotope..index + isotope + ele] == *possible.0 {
                        element = Some(possible.1);
                        break;
                    }
                }
                let num = value[index + isotope + ele..index + len]
                    .parse::<i16>()
                    .map_err(|e| e.to_string())?;
                let isotope = value[index..index + isotope]
                    .parse::<u16>()
                    .map_err(|e| e.to_string())?;

                if !result.add((element.unwrap(), Some(isotope), num)) {
                    return Err(format!("Invalid isotope ({}) added for element ({}) in pro forma molecular formula ({})", isotope, element.unwrap(), value));
                }
                element = None;
                index += len + 1;
            }
            b'-' | b'0'..=b'9' if element.is_some() => {
                let (num, len) = std::str::from_utf8(
                    &bytes
                        .iter()
                        .skip(index)
                        .take_while(|c| c.is_ascii_digit() || **c == b'-')
                        .copied()
                        .collect::<Vec<_>>(),
                )
                .map(|v| (v.parse::<i16>().map_err(|e| e.to_string()), v.len()))
                .map_err(|e| e.to_string())?;
                let num = num?;
                if num != 0 {
                    if !result.add((element.unwrap(), None, num)) {
                        return Err(format!("An element without a defined mass ({}) was used in a pro forma molecular formula ({})", element.unwrap(), value));
                    }
                }
                element = None;
                index += len;
            }
            b' ' => index += 1,
            _ => {
                let mut found = false;
                for possible in ELEMENT_PARSE_LIST {
                    if value[index..].starts_with(possible.0) {
                        element = Some(possible.1);
                        index += possible.0.len();
                        found = true;
                        break;
                    }
                }
                if !found {
                    return Err(format!(
                        "Could not parse ProForma elemental formula, broke down at index: {index}"
                    ));
                }
            }
        }
    }
    if let Some(element) = element {
        if !result.add((element, None, 1)) {
            return Err(format!( "An element without a defined mass ({element}) was used in a pro forma molecular formula ({value})",));
        }
    }
    Ok(result)
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
pub fn end_of_enclosure(chars: &[u8], start: usize, open: u8, close: u8) -> Option<usize> {
    let mut state = 1;
    for (i, ch) in chars[start..].iter().enumerate() {
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
/// Get the next number, returns length in bytes and the number
pub fn next_num(chars: &[u8], mut start: usize, allow_only_sign: bool) -> Option<(usize, isize)> {
    let mut sign = 1;
    let mut sign_set = false;
    if chars[start] == b'-' {
        sign = -1;
        start += 1;
        sign_set = true;
    } else if chars[start] == b'+' {
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
            .unwrap();
        Some((usize::from(sign_set) + len, sign * num))
    }
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
