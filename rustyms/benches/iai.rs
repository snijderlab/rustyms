use std::hint::black_box;

use rustyms::align::*;
use rustyms::system::dalton;
use rustyms::system::Mass;
use rustyms::Simple;
use rustyms::*;

use iai_callgrind::{
    library_benchmark, library_benchmark_group, main, LibraryBenchmarkConfig, Tool, ValgrindTool,
};

#[inline(never)]
fn setup(a: &str, b: &str) -> (LinearPeptide<Simple>, LinearPeptide<Simple>) {
    let _force_elements_init = black_box(AminoAcid::A.formulas(SequencePosition::default(), 0));
    (
        LinearPeptide::pro_forma(a, None).unwrap().simple().unwrap(),
        LinearPeptide::pro_forma(b, None).unwrap().simple().unwrap(),
    )
}

#[inline(never)]
fn setup_simple() -> (LinearPeptide<Simple>, LinearPeptide<Simple>) {
    setup("ANAGRS", "AGGQRS")
}

#[inline(never)]
fn setup_igha() -> (LinearPeptide<Simple>, LinearPeptide<Simple>) {
    setup("ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNSSQDVTVPCRVPPPPPCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGATFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAQPWNHGETFTCTAAHPELKTPLTANITKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTYAVTSILRVAAEDWKKGETFSCMVGHEALPLAFTQKTIDRMAGSCCVADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGKREGPQY", "ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVTVPCPVPSTPPTPSPSTPPTPSPSCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGVTFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAEPWNHGKTFTCTAAYPESKTPLTATLSKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTFAVTSILRVAAEDWKKGDTFSCMVGHEALPLAFTQKTIDRLADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGNREGPQY")
}

#[library_benchmark]
#[bench::simple_1(setup_simple())]
#[bench::igha_1(setup_igha())]
pub fn align_1(setup: (LinearPeptide<Simple>, LinearPeptide<Simple>)) {
    align::<1, Simple, Simple>(
        &setup.0,
        &setup.1,
        matrix::BLOSUM62,
        Tolerance::new_absolute(Mass::new::<dalton>(0.01)),
        AlignType::GLOBAL,
    );
}

#[library_benchmark]
#[bench::simple_4(setup_simple())]
#[bench::igha_4(setup_igha())]
#[bench::ambiguous_not(setup("ANQRS", "ANQRS"))]
#[bench::ambiguous_a(setup("ANZRS", "ANQRS"))]
#[bench::ambiguous_b(setup("ANQRS", "ABQRS"))]
#[bench::ambiguous_ab(setup("ANZRS", "ABQRS"))]
// #[bench::igha_8(setup_igha(Some(8)))]
pub fn align_4(setup: (LinearPeptide<Simple>, LinearPeptide<Simple>)) {
    align::<4, Simple, Simple>(
        &setup.0,
        &setup.1,
        matrix::BLOSUM62,
        Tolerance::new_absolute(Mass::new::<dalton>(0.01)),
        AlignType::GLOBAL,
    );
}

#[library_benchmark]
#[bench::simple_unbounded(setup_simple())]
// #[bench::igha_8(setup_igha(Some(8)))]
pub fn align_unbounded(setup: (LinearPeptide<Simple>, LinearPeptide<Simple>)) {
    align::<{ u16::MAX }, Simple, Simple>(
        &setup.0,
        &setup.1,
        matrix::BLOSUM62,
        Tolerance::new_absolute(Mass::new::<dalton>(0.01)),
        AlignType::GLOBAL,
    );
}

library_benchmark_group!(name = alignment; benchmarks = align_1, align_4, align_unbounded);

main!(config = LibraryBenchmarkConfig::default()
.tool(Tool::new(ValgrindTool::DHAT)).tool(Tool::new(ValgrindTool::Massif)); library_benchmark_groups = alignment);
