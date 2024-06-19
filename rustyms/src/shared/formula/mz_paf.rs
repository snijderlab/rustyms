use crate::{
    error::{Context, CustomError},
    helper_functions::{explain_number_error, next_number},
    MolecularFormula, ELEMENT_PARSE_LIST,
};
use std::ops::RangeBounds;

impl MolecularFormula {
    /// Parse mzPAF formulas: `C13H9N1O`.
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It can panic if the string contains non UTF8 symbols.
    #[allow(dead_code)]
    pub fn from_mz_paf(value: &str) -> Result<Self, CustomError> {
        Self::from_mz_paf_inner(value, ..)
    }

    /// See [`Self::from_mz_paf`]. This is a variant to help in parsing a part of a larger line.
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It can panic if the string contains non UTF8 symbols.
    pub fn from_mz_paf_inner(
        value: &str,
        range: impl RangeBounds<usize>,
    ) -> Result<Self, CustomError> {
        let mut index = match range.start_bound() {
            std::ops::Bound::Unbounded => 0,
            std::ops::Bound::Included(s) => *s,
            std::ops::Bound::Excluded(s) => s + 1,
        };
        let end = match range.end_bound() {
            std::ops::Bound::Unbounded => value.len().saturating_sub(1),
            std::ops::Bound::Included(s) => *s,
            std::ops::Bound::Excluded(s) => s.saturating_sub(1),
        };
        let mut result = Self::default();
        while index <= end {
            let start_index = index;
            let element_text: String = value[index..]
                .chars()
                .take(2)
                .collect::<String>()
                .to_ascii_lowercase();
            let mut element = None;
            for possible in ELEMENT_PARSE_LIST {
                if element_text.starts_with(possible.0) {
                    element = Some(possible.1);
                    index += possible.0.len();
                    break;
                }
            }
            let end_index = index;
            if let Some(element) = element.take() {
                let mut amount = 1;
                if let Some(number) = next_number::<false, false, u16>(value, index..) {
                    index += number.0;
                    amount = number.2.map_err(|err| {
                        CustomError::error(
                            "Invalid formula",
                            format!("The element amount number {}", explain_number_error(&err)),
                            Context::line_range(0, value, start_index..=end_index),
                        )
                    })?;
                }
                if !result.add((element, None, i32::from(amount))) {
                    return Err(CustomError::error(
                        "Invalid mzPAF molecular formula",
                        "An element with undefined mass was used",
                        Context::line(
                            0,
                            value,
                            index,
                            value[index..]
                                .chars()
                                .next()
                                .map(char::len_utf8)
                                .unwrap_or_default(),
                        ),
                    ));
                }
            } else {
                return Err(CustomError::error(
                    "Invalid mzPAF molecular formula",
                    "Not a valid character in formula",
                    Context::line(
                        0,
                        value,
                        index,
                        value[index..]
                            .chars()
                            .next()
                            .map(char::len_utf8)
                            .unwrap_or_default(),
                    ),
                ));
            }
        }
        Ok(result)
    }
}
