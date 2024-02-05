use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use rustyms::align::*;
use rustyms::system::dalton;
use rustyms::system::Mass;
use rustyms::*;

pub fn simple_1(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ANAGRS")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("AGGQRS")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("simple 1", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                Some(1),
            )
        })
    });
}

pub fn simple_4(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ANAGRS")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("AGGQRS")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("simple 4", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                Some(4),
            )
        })
    });
}

pub fn simple_unbounded(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ANAGRS")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("AGGQRS")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("simple ∞", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                None,
            )
        })
    });
}

pub fn igha_1(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNSSQDVTVPCRVPPPPPCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGATFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAQPWNHGETFTCTAAHPELKTPLTANITKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTYAVTSILRVAAEDWKKGETFSCMVGHEALPLAFTQKTIDRMAGSCCVADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGKREGPQY")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVTVPCPVPSTPPTPSPSTPPTPSPSCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGVTFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAEPWNHGKTFTCTAAYPESKTPLTATLSKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTFAVTSILRVAAEDWKKGDTFSCMVGHEALPLAFTQKTIDRLADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGNREGPQY")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("igha 1", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                Some(1),
            )
        })
    });
}

pub fn igha_4(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNSSQDVTVPCRVPPPPPCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGATFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAQPWNHGETFTCTAAHPELKTPLTANITKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTYAVTSILRVAAEDWKKGETFSCMVGHEALPLAFTQKTIDRMAGSCCVADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGKREGPQY")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVTVPCPVPSTPPTPSPSTPPTPSPSCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGVTFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAEPWNHGKTFTCTAAYPESKTPLTATLSKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTFAVTSILRVAAEDWKKGDTFSCMVGHEALPLAFTQKTIDRLADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGNREGPQY")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("igha 4", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                Some(4),
            )
        })
    });
}

pub fn igha_8(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNSSQDVTVPCRVPPPPPCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGATFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAQPWNHGETFTCTAAHPELKTPLTANITKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTYAVTSILRVAAEDWKKGETFSCMVGHEALPLAFTQKTIDRMAGSCCVADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGKREGPQY")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVTVPCPVPSTPPTPSPSTPPTPSPSCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGVTFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAEPWNHGKTFTCTAAYPESKTPLTATLSKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTFAVTSILRVAAEDWKKGDTFSCMVGHEALPLAFTQKTIDRLADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGNREGPQY")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("igha 8", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                Some(8),
            )
        })
    });
}

pub fn igha_32(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNSSQDVTVPCRVPPPPPCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGATFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAQPWNHGETFTCTAAHPELKTPLTANITKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTYAVTSILRVAAEDWKKGETFSCMVGHEALPLAFTQKTIDRMAGSCCVADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGKREGPQY")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVTVPCPVPSTPPTPSPSTPPTPSPSCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGVTFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAEPWNHGKTFTCTAAYPESKTPLTATLSKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTFAVTSILRVAAEDWKKGDTFSCMVGHEALPLAFTQKTIDRLADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGNREGPQY")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("igha 32", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                Some(32),
            )
        })
    });
}

pub fn igha_unbounded(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNSSQDVTVPCRVPPPPPCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGATFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAQPWNHGETFTCTAAHPELKTPLTANITKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTYAVTSILRVAAEDWKKGETFSCMVGHEALPLAFTQKTIDRMAGSCCVADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGKREGPQY")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVTVPCPVPSTPPTPSPSTPPTPSPSCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGVTFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAEPWNHGKTFTCTAAYPESKTPLTATLSKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTFAVTSILRVAAEDWKKGDTFSCMVGHEALPLAFTQKTIDRLADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGNREGPQY")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("igha ∞", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                None,
            )
        })
    });
}

pub fn simple_not_ambiguous(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ANQRS")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("ANQRS")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("not ambiguous", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                Some(4),
            )
        })
    });
}

pub fn simple_ambiguous_a(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ANZRS")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("ANQRS")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("ambiguous A", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                Some(4),
            )
        })
    });
}

pub fn simple_ambiguous_b(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ANQRS")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("ABQRS")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("ambiguous B", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                Some(4),
            )
        })
    });
}

pub fn simple_ambiguous_ab(c: &mut Criterion) {
    let a = ComplexPeptide::pro_forma("ANZRS")
        .unwrap()
        .singular()
        .unwrap();
    let b = ComplexPeptide::pro_forma("ABQRS")
        .unwrap()
        .singular()
        .unwrap();
    c.bench_function("ambiguous AB", |bencher| {
        bencher.iter(|| {
            align(
                a.clone(),
                b.clone(),
                align::BLOSUM62,
                Tolerance::new_absolute(Mass::new::<dalton>(0.1)),
                AlignType::GLOBAL,
                Some(4),
            )
        })
    });
}

criterion_group!(
    name = alignments;
    config = Criterion::default()
        .significance_level(0.1)
        .sample_size(200)
        .measurement_time(Duration::from_secs(20));
    targets = simple_1,
    simple_4,
    simple_unbounded,
    igha_1,
    igha_4,
    igha_8,
    // igha_32,
    // igha_unbounded,
    simple_not_ambiguous,
    simple_ambiguous_a,
    simple_ambiguous_b,
    simple_ambiguous_ab,
);
criterion_main!(alignments);
