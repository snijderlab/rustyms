use crate::*;

#[test]
fn triple_a() {
    // Compare rustyms with https://proteomicsresource.washington.edu/cgi-bin/fragment.cgi
    #[allow(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (44.049476, "a+1"),
        (72.044390, "b+1"),
        (89.070939, "c+1"),
        (115.086589, "a+2"),
        (143.081504, "b+2"),
        (160.108053, "c+2"),
        (187.071333, "x+2"),
        (161.092069, "y+2"),
        (144.065520, "z+2"),
        (145.073345, "z·+2"),
        (116.034220, "x+1"),
        (90.054955, "y+1"),
        (73.028406, "z+1"),
        (74.036231, "z·+1"),
        (232.129183, "precursor"),
    ];
    let model = Model {
        a: (Location::SkipC(1), vec![]),
        b: (Location::SkipC(1), vec![]),
        c: (Location::SkipC(1), vec![]),
        d: (Location::None, vec![]),
        v: (Location::None, vec![]),
        w: (Location::None, vec![]),
        x: (Location::SkipN(1), vec![]),
        y: (Location::SkipN(1), vec![]),
        z: (Location::SkipN(1), vec![]),
        precursor: vec![],
        ppm: MassOverCharge::new::<mz>(20.0),
    };
    test(
        theoretical_fragments,
        &Peptide::pro_forma("AAA").unwrap(),
        &model,
    );
}

#[test]
fn with_modifications() {
    // Compare rustyms with https://proteomicsresource.washington.edu/cgi-bin/fragment.cgi mods: -17.02655@[ 15.99491@
    #[allow(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (84.044389, "a1"),
        (155.081503, "a2"),
        (226.118617, "a3"),
        (112.039304, "b1"),
        (183.076418, "b2"),
        (254.113532, "b3"),
        (129.065853, "c1"),
        (200.102967, "c2"),
        (271.140081, "c3"),
        (334.106728, "x3"),
        (263.069614, "x2"),
        (192.0325, "x1"),
        (308.127463, "y3"),
        (237.09035, "y2"),
        (166.053236, "y1"),
        (291.100914, "z3"),
        (220.063801, "z2"),
        (149.026687, "z1"),
        (292.108739, "z·3"),
        (221.071626, "z·2"),
        (150.034512, "z·1"),
        (419.159491, "precursor"),
    ];
    let model = Model {
        a: (Location::SkipC(1), vec![]),
        b: (Location::SkipC(1), vec![]),
        c: (Location::SkipC(1), vec![]),
        d: (Location::None, vec![]),
        v: (Location::None, vec![]),
        w: (Location::None, vec![]),
        x: (Location::SkipN(1), vec![]),
        y: (Location::SkipN(1), vec![]),
        z: (Location::SkipN(1), vec![]),
        precursor: vec![],
        ppm: MassOverCharge::new::<mz>(20.0),
    };
    test(
        theoretical_fragments,
        &Peptide::pro_forma("[Gln->pyro-Glu]-QAAM[Oxidation]").unwrap(),
        &model,
    );
}

#[test]
fn with_possible_modifications() {
    // Compare rustyms with https://proteomicsresource.washington.edu/cgi-bin/fragment.cgi mods: -17.02655@[ 15.99491@
    #[allow(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (148.042671, "M1-b1"),
        (150.058326, "M1-y1"),
        (132.047761, "M2-b1"),
        (166.053236, "M2-y2"),
        (297.09372, "precursor"),
    ];
    let model = Model {
        a: (Location::None, vec![]),
        b: (Location::SkipC(1), vec![]),
        c: (Location::None, vec![]),
        d: (Location::None, vec![]),
        v: (Location::None, vec![]),
        w: (Location::None, vec![]),
        x: (Location::None, vec![]),
        y: (Location::SkipN(1), vec![]),
        z: (Location::None, vec![]),
        precursor: vec![],
        ppm: MassOverCharge::new::<mz>(20.0),
    };
    test(
        theoretical_fragments,
        &Peptide::pro_forma("M[Oxidation#a1]M[#a1]").unwrap(),
        &model,
    );
}

#[test]
fn all_aminoacids() {
    // Compare rustyms with https://proteomicsresource.washington.edu/cgi-bin/fragment.cgi
    #[allow(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (44.049476, "a1"),
        (200.150587, "a2"),
        (314.193514, "a3"),
        (429.220457, "a4"),
        (532.229642, "a5"),
        (660.288219, "a6"),
        (789.330812, "a7"),
        (846.352276, "a8"),
        (983.411188, "a9"),
        (1096.495252, "a10"),
        (1209.579316, "a11"),
        (1337.674279, "a12"),
        (1468.714764, "a13"),
        (1615.783178, "a14"),
        (1712.835942, "a15"),
        (1799.86797, "a16"),
        (1900.915649, "a17"),
        (2086.994962, "a18"),
        (2250.05829, "a19"),
        (72.04439, "b1"),
        (228.145501, "b2"),
        (342.188429, "b3"),
        (457.215372, "b4"),
        (560.224556, "b5"),
        (688.283134, "b6"),
        (817.325727, "b7"),
        (874.347191, "b8"),
        (1011.406103, "b9"),
        (1124.490167, "b10"),
        (1237.574231, "b11"),
        (1365.669194, "b12"),
        (1496.709678, "b13"),
        (1643.778092, "b14"),
        (1740.830856, "b15"),
        (1827.862885, "b16"),
        (1928.910563, "b17"),
        (2114.989876, "b18"),
        (2278.053205, "b19"),
        (89.070939, "c1"),
        (245.17205, "c2"),
        (359.214978, "c3"),
        (474.241921, "c4"),
        (577.251105, "c5"),
        (705.309683, "c6"),
        (834.352276, "c7"),
        (891.37374, "c8"),
        (1028.432652, "c9"),
        (1141.516716, "c10"),
        (1254.60078, "c11"),
        (1382.695743, "c12"),
        (1513.736227, "c13"),
        (1660.804641, "c14"),
        (1757.857405, "c15"),
        (1844.889434, "c16"),
        (1945.937112, "c17"),
        (2132.016425, "c18"),
        (2295.079754, "c19"),
        (2350.074334, "x19"),
        (2193.973223, "x18"),
        (2079.930296, "x17"),
        (1964.903353, "x16"),
        (1861.894168, "x15"),
        (1733.83559, "x14"),
        (1604.792997, "x13"),
        (1547.771534, "x12"),
        (1410.712622, "x11"),
        (1297.628558, "x10"),
        (1184.544494, "x9"),
        (1056.449531, "x8"),
        (925.409046, "x7"),
        (778.340632, "x6"),
        (681.287868, "x5"),
        (594.25584, "x4"),
        (493.208161, "x3"),
        (307.128848, "x2"),
        (144.06552, "x1"),
        (2324.09507, "y19"),
        (2167.993959, "y18"),
        (2053.951031, "y17"),
        (1938.924088, "y16"),
        (1835.914903, "y15"),
        (1707.856326, "y14"),
        (1578.813733, "y13"),
        (1521.792269, "y12"),
        (1384.733357, "y11"),
        (1271.649293, "y10"),
        (1158.565229, "y9"),
        (1030.470266, "y8"),
        (899.429781, "y7"),
        (752.361367, "y6"),
        (655.308604, "y5"),
        (568.276575, "y4"),
        (467.228897, "y3"),
        (281.149584, "y2"),
        (118.086255, "y1"),
        (2307.06852, "z19"),
        (2150.967409, "z18"),
        (2036.924482, "z17"),
        (1921.897539, "z16"),
        (1818.888354, "z15"),
        (1690.829777, "z14"),
        (1561.787184, "z13"),
        (1504.76572, "z12"),
        (1367.706808, "z11"),
        (1254.622744, "z10"),
        (1141.53868, "z9"),
        (1013.443717, "z8"),
        (882.403232, "z7"),
        (735.334818, "z6"),
        (638.282055, "z5"),
        (551.250026, "z4"),
        (450.202348, "z3"),
        (264.123035, "z2"),
        (101.059706, "z1"),
        (2308.076345, "z·19"),
        (2151.975234, "z·18"),
        (2037.932307, "z·17"),
        (1922.905364, "z·16"),
        (1819.896179, "z·15"),
        (1691.837602, "z·14"),
        (1562.795009, "z·13"),
        (1505.773545, "z·12"),
        (1368.714633, "z·11"),
        (1255.630569, "z·10"),
        (1142.546505, "z·9"),
        (1014.451542, "z·8"),
        (883.411057, "z·7"),
        (736.342643, "z·6"),
        (639.28988, "z·5"),
        (552.257851, "z·4"),
        (451.210173, "z·3"),
        (265.13086, "z·2"),
        (102.067531, "z·1"),
        (2395.132183, "precursor"),
    ];
    let model = Model {
        a: (Location::SkipC(1), vec![]),
        b: (Location::SkipC(1), vec![]),
        c: (Location::SkipC(1), vec![]),
        d: (Location::None, vec![]),
        v: (Location::None, vec![]),
        w: (Location::None, vec![]),
        x: (Location::SkipN(1), vec![]),
        y: (Location::SkipN(1), vec![]),
        z: (Location::SkipN(1), vec![]),
        precursor: vec![],
        ppm: MassOverCharge::new::<mz>(20.0),
    };
    test(
        theoretical_fragments,
        &Peptide::pro_forma("ARNDCQEGHILKMFPSTWYV").unwrap(),
        &model,
    );
}

fn test(theoretical_fragments: &[(f64, &str)], peptide: &Peptide, model: &Model) {
    let mut calculated_fragments = peptide
        .generate_theoretical_fragments::<MonoIsotopic>(Charge::new::<e>(1.0), model)
        .unwrap();
    dbg!(&calculated_fragments);
    let mut combined = Vec::new();
    let mut found = Vec::new();
    let mut this_found;
    for goal in theoretical_fragments {
        this_found = false;
        for fragment in 0..calculated_fragments.len() {
            if calculated_fragments[fragment]
                .theoretical_mass
                .ppm(da(goal.0))
                < 20.0
            {
                combined.push((goal.0, goal.1, calculated_fragments.remove(fragment)));
                this_found = true;
                break;
            }
        }
        found.push(this_found);
    }
    dbg!(&combined);

    for left in &calculated_fragments {
        println!("Excess fragments: {left}");
    }
    let not_found: Vec<usize> = found
        .iter()
        .enumerate()
        .filter_map(|(index, is_found)| (!is_found).then_some(index))
        .collect();

    for index in &not_found {
        let (mass, name) = theoretical_fragments[*index];
        println!("Not found: {mass} {name}");
    }
    assert_eq!(not_found.len(), 0, "Not all needed fragments are found");
    assert_eq!(
        calculated_fragments.len(),
        0,
        "Not all generated fragments are accounted for"
    );
}
