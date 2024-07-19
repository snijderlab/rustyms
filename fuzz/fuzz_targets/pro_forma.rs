use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data) {
            let _ = rustyms::CompoundPeptidoform::pro_forma(s, None);
            let _ = rustyms::LinearPeptide::sloppy_pro_forma(
                s,
                0..s.len(),
                None,
                rustyms::SloppyParsingParameters::default(),
            );
        }
    });
}
