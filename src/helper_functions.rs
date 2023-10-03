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
        match self {
            Some(o) => o.map(|v| Some(v)),
            None => Ok(None),
        }
    }
}

#[allow(dead_code)]
pub fn parse_named_counter<T: Copy>(
    value: &str,
    names: &[(&str, T)],
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
                if value[index..].starts_with(name.0) {
                    index += name.0.len();
                    let num = &value[index..]
                        .chars()
                        .skip_while(char::is_ascii_whitespace)
                        .take_while(|c| c.is_ascii_digit() || (allow_negative && *c == '-'))
                        .collect::<String>()
                        .trim()
                        .to_string();
                    if num.is_empty() {
                        output.push((name.1, 1));
                    } else {
                        output.push((name.1, num.parse().unwrap()));
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

                result.add((element.unwrap(), isotope, num));
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
                    result.add((element.unwrap(), 0, num));
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
        result.add((element, 0, 1));
    }
    Ok(result)
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
