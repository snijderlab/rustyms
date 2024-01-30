use rustyms::align::*;
use rustyms::system::dalton;
use rustyms::system::Mass;
use rustyms::*;

use iai_callgrind::{library_benchmark, library_benchmark_group, main};

fn setup(
    a: &str,
    b: &str,
    max_steps: Option<usize>,
) -> (LinearPeptide, LinearPeptide, Option<usize>) {
    (
        ComplexPeptide::pro_forma(a).unwrap().singular().unwrap(),
        ComplexPeptide::pro_forma(b).unwrap().singular().unwrap(),
        max_steps,
    )
}

fn setup_simple(max_steps: Option<usize>) -> (LinearPeptide, LinearPeptide, Option<usize>) {
    setup("ANAGRS", "AGGQRS", max_steps)
}

fn setup_igha(max_steps: Option<usize>) -> (LinearPeptide, LinearPeptide, Option<usize>) {
    setup("ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNSSQDVTVPCRVPPPPPCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGATFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAQPWNHGETFTCTAAHPELKTPLTANITKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTYAVTSILRVAAEDWKKGETFSCMVGHEALPLAFTQKTIDRMAGSCCVADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGKREGPQY", "ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVTVPCPVPSTPPTPSPSTPPTPSPSCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGVTFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAEPWNHGKTFTCTAAYPESKTPLTATLSKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTFAVTSILRVAAEDWKKGDTFSCMVGHEALPLAFTQKTIDRLADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGNREGPQY", max_steps)
}

#[library_benchmark]
#[bench::simple_1(setup_simple(Some(1)))]
#[bench::simple_4(setup_simple(Some(4)))]
#[bench::simple_unbounded(setup_simple(None))]
#[bench::igha_1(setup_igha(Some(1)))]
#[bench::igha_4(setup_igha(Some(4)))]
#[bench::ambiguous_not(setup("ANQRS", "ANQRS", Some(4)))]
#[bench::ambiguous_a(setup("ANZRS", "ANQRS", Some(4)))]
#[bench::ambiguous_b(setup("ANQRS", "ABQRS", Some(4)))]
#[bench::ambiguous_ab(setup("ANZRS", "ABQRS", Some(4)))]
// #[bench::igha_8(setup_igha(Some(8)))]
pub fn align(setup: (LinearPeptide, LinearPeptide, Option<usize>)) {
    rustyms::align::align(
        setup.0,
        setup.1,
        align::BLOSUM62,
        Tolerance::new_absolute(Mass::new::<dalton>(0.01)),
        Type::GLOBAL,
        setup.2,
    );
}

library_benchmark_group!(name = alignment; benchmarks = align);

main!(library_benchmark_groups = alignment);
