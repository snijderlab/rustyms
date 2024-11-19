#![allow(clippy::missing_panics_doc)]
use std::io::BufReader;

use crate::identification::{test_format, DeepNovoFamilyData, DeepNovoFamilyVersion};

#[test]
fn deepnovo() {
    match test_format::<DeepNovoFamilyData>(
        BufReader::new(DEEPNOVO_V0_0_1.as_bytes()),
        None,
        false,
        true,
        Some(DeepNovoFamilyVersion::DeepNovoV0_0_1),
    ) {
        Ok(n) => assert_eq!(n, 19),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

#[test]
fn pointnovo() {
    match test_format::<DeepNovoFamilyData>(
        BufReader::new(POINTNOVO_V0_0_1.as_bytes()),
        None,
        false,
        true,
        Some(DeepNovoFamilyVersion::PointNovoFamily),
    ) {
        Ok(n) => assert_eq!(n, 20),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

#[test]
fn biatnovo() {
    match test_format::<DeepNovoFamilyData>(
        BufReader::new(BIATNOVO_V0_1.as_bytes()),
        None,
        false,
        true,
        Some(DeepNovoFamilyVersion::PointNovoFamily),
    ) {
        Ok(n) => assert_eq!(n, 20),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

#[test]
fn pgpointnovo() {
    match test_format::<DeepNovoFamilyData>(
        BufReader::new(PGPOINTNOVO_V1_0_6.as_bytes()),
        None,
        false,
        true,
        Some(DeepNovoFamilyVersion::PointNovoFamily),
    ) {
        Ok(n) => assert_eq!(n, 20),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

const DEEPNOVO_V0_0_1: &str = "scan	predicted_sequence	predicted_score	predicted_position_score
14	V,A,P,E,I,P,V,Y,I,N,E,V,A,I,V,I,I,P,A,I,A,A,R	-1.53	-2.72,-0.55,-1.77,-0.85,-1.36,-0.46,-1.95,-0.89,-2.23,-1.98,-0.38,-1.17,-2.06,-1.45,-1.43,-2.26,-1.98,-1.47,-2.14,-1.91,-2.26,-1.89
15	I,I,P,K,P,I,V,I,I,G,A,Nmod,G,T,T,V,I,V,G,M,P,A,G,A,K	-1.05	-3.45,-4.20,-2.48,-0.37,-3.40,-1.80,-1.82,-1.67,-0.16,-0.39,-1.86,-0.00,-0.01,-0.00,-0.00,-0.03,-0.12,-0.12,-0.31,-1.24,-0.55,-0.68,-0.45,-1.16
2535	E,E,E,T,H,I,K,P,S,N,T,K,V,D,K	-0.42	-0.96,-0.10,-2.23,-0.04,-0.31,-0.06,-0.00,-0.03,-0.32,-1.52,-0.09,-0.30,-0.30,-0.06
2536	N,F,P,I,N,P,P,P,E,N,S,I,Cmod,K,K	-0.78	-1.86,-0.29,-1.46,-1.12,-0.04,-0.61,-0.40,-0.56,-1.87,-0.79,-0.47,-1.51,-0.46,-0.07
2537	A,P,Q,S,G,S,V,S,S,D,D,T,Q,N,Q,H,E,N,Q,H,T,N,S,T,Y,R	-0.61	-0.55,-1.66,-0.75,-0.50,-1.33,-1.52,-0.94,-0.13,-0.59,-0.63,-0.04,-0.23,-0.93,-0.90,-0.00,-0.51,-0.24,-0.67,-0.58,-0.15,-1.00,-0.48,-0.32,-1.13,-0.14
2538	H,G,G,N,S,S,S,N,N,G,T,G,M,P,P,Q,D,A,G,G,G,V,G,S,A,F,Y,D,I,H	-0.55	-0.73,-0.69,-0.82,-0.11,-0.86,-0.73,-0.03,-0.21,-0.45,-1.39,-0.05,-0.01,-0.40,-1.24,-0.07,-1.46,-0.18,-0.42,-0.01,-0.00,-0.96,-0.22,-0.13,-0.52,-0.17,-0.41,-0.23,-1.93,-1.91
2539	H,I,D,Cmod,A,A,Q,D,Q,Y,T,Q,A,Q,Y,N,Q,H,T,N,S,T,Y,R	-0.57	-0.04,-2.49,-0.26,-0.66,-0.68,-0.08,-1.62,-0.64,-0.08,-0.96,-0.43,-1.19,-1.03,-0.05,-0.18,-0.63,-0.10,-0.04,-0.28,-0.07,-0.13,-1.84,-0.15
2540	S,G,S,D,Y,Y,Y,G,F,Y,A,G,P,V,S,P,Y,T,H,A,T,T,S,G,D,I,H	-0.76	-2.37,-0.13,-1.19,-0.87,-0.27,-0.23,-1.02,-0.15,-0.26,-0.41,-1.13,-0.92,-1.05,-0.65,-0.35,-0.56,-0.17,-0.07,-0.29,-0.71,-1.07,-0.71,-0.95,-1.30,-2.76,-0.83
2541	H,V,E,G,E,S,I,D,A,Y,D,S,S,S,S,S,T,P,S,R,P,A,S,S,S,N,M,G,K	-0.59	-0.56,-1.03,-0.07,-1.77,-0.68,-0.17,-0.68,-0.62,-0.01,-0.12,-0.43,-1.43,-1.30,-0.12,-0.20,-0.71,-0.92,-0.32,-0.40,-0.35,-0.11,-1.63,-0.05,-0.73,-0.81,-1.31,-0.28,-0.22
2542	H,S,G,S,S,S,S,T,S,S,N,S,S,S,S,S,N,T,V,Y,T,H,A,S,S,S,S,I,I,H	-0.65	-0.71,-1.20,-1.20,-0.85,-0.43,-0.27,-0.49,-0.11,-0.58,-0.35,-0.48,-1.00,-1.16,-0.30,-0.10,-0.19,-1.26,-0.03,-0.03,-0.07,-0.09,-0.15,-1.08,-0.19,-1.13,-1.28,-1.28,-1.91,-1.43
2546	D,N,M,Q,H,I,K,P,S,N,T,K,V,D,K	-0.55	-0.89,-2.42,-1.65,-0.47,-0.72,-0.03,-0.00,-0.07,-0.54,-0.69,-0.43,-0.18,-0.09,-0.08
2547	E,T,P,G,A,N,T,P,V,P,S,N,T,K,V,D,K	-0.65	-1.01,-0.38,-1.43,-0.46,-1.41,-1.73,-0.00,-0.40,-0.02,-0.07,-0.69,-1.44,-0.76,-0.31,-0.77,-0.08
2550	S,F,Y,H,I,K,G,E,I,T,K,V,D,K	-0.53	-0.58,-2.46,-0.26,-0.74,-2.08,-0.00,-0.04,-0.12,-0.85,-0.12,-0.14,-0.04,-0.00
2551	I,A,N,V,N,E,I,S,I,N,H,V,V,S,R	-0.57	-0.84,-0.10,-0.27,-0.24,-0.52,-0.52,-0.41,-0.30,-0.27,-1.09,-1.42,-0.12,-0.22,-2.20
2552	H,V,E,Cmod,P,A,A,I,Y,Cmod,P,S,Y,V,T,S,I,A,G,Y,G,S,S,Y,D	-0.67	-0.86,-0.81,-0.00,-0.54,-1.68,-0.50,-0.86,-0.41,-0.02,-0.11,-1.76,-0.18,-0.08,-0.01,-0.50,-0.83,-0.51,-0.05,-0.12,-0.07,-1.89,-1.09,-0.86,-3.08
2553	N,A,E,E,G,Q,F,D,G,T,D,P,E,A,W,T,T,V,A,G,Y,P,G,V,A,A,N	-0.70	-2.09,-0.01,-0.14,-0.26,-1.12,-0.44,-0.54,-0.25,-0.15,-1.36,-0.28,-0.79,-0.03,-0.05,-1.10,-0.85,-0.38,-0.91,-0.02,-0.33,-0.03,-2.22,-1.48,-0.32,-1.46,-2.26
2554	H,V,E,Cmod,P,A,A,I,Y,Cmod,G,A,G,P,Y,Y,T,Y,Y,S,E,S,W,G,K	-0.59	-0.95,-1.26,-0.01,-0.37,-1.77,-0.64,-1.03,-0.34,-0.02,-0.13,-0.80,-1.25,-1.40,-0.78,-0.03,-0.08,-0.56,-0.22,-0.25,-0.04,-0.53,-1.27,-0.96,-0.00
2559	I,S,H,V,S,S,G,G,G,A,S,E,V,I,I,V,G,A,D	-0.39	-0.03,-1.84,-0.93,-1.22,-0.48,-0.84,-0.00,-0.00,-0.01,-0.00,-0.00,-1.15,-0.40,-0.11,-0.00,-0.14,-0.14,-0.03
2560	I,I,G,G,I,M,T,Y,G,D,F,N,Q,T,P,K	-0.59	-0.26,-1.40,-1.23,-0.65,-0.37,-0.11,-0.04,-0.49,-0.22,-1.55,-0.72,-0.66,-0.97,-0.57,-0.27";

const PGPOINTNOVO_V1_0_6: &str = "feature_id	feature_area	predicted_sequence	predicted_score	predicted_position_score	precursor_mz	precursor_charge	protein_access_id	scan_list_middle	scan_list_original	predicted_score_max
F1:15	1.0	S,A,Q,E,S,S,S,S,S,S,S,L,L,L,V,L,L,L,L,L,K,K,K	-2.65	-2.29,-2.19,-2.25,-2.31,-2.31,-2.46,-2.45,-2.45,-2.44,-2.45,-2.45,-3.09,-3.13,-3.09,-2.20,-3.06,-3.09,-3.09,-3.09,-3.11,-3.93,-4.07,-0.07	1216.224487304688	2.0	DENOVO	F1:15	F1:15	-2.65
F1:2538	1.0	T,K,P,E,S,A,G,P,G,N,T,P,Y,G,E,E,G,E,E,S,E,D,T,S,T,Y,R	-1.47	-2.49,-0.85,-0.44,-2.59,-1.56,-1.82,-2.12,-1.74,-1.68,-1.79,-0.75,-1.07,-1.94,-2.38,-2.27,-0.74,-2.58,-1.07,-0.83,-1.86,-1.83,-1.09,-0.23,-1.49,-0.40,-2.00,-0.07	722.81494140625	4.0	DENOVO	F1:2538	F1:2538	-1.47
F1:2542	1.0	S,G,G,S,D,A,S,E,S,D,G,Y,T,E,P,V,K,T,Y,T,H,G,P,S,G,T,S,Q,P,A	-1.39	-2.43,-1.00,-1.54,-2.41,-0.66,-1.09,-2.54,-1.19,-2.21,-0.43,-0.32,-1.73,-0.88,-0.85,-2.35,-2.54,-2.14,-1.07,-0.43,-0.39,-0.16,-0.89,-1.83,-2.40,-1.98,-2.38,-1.40,-1.70,-0.51,-0.16	743.328918457031	4.0	DENOVO	F1:2542	F1:2542	-1.39
F1:2551	1.0	L,A,N,V,G,G,H,K,F,G,Q,P,E,V,V,L	-1.02	-1.89,-0.11,-0.60,-0.24,-1.42,-0.21,-1.24,-0.21,-0.64,-0.64,-0.26,-0.87,-2.82,-1.72,-3.19,-0.30	416.983551025391	4.0	DENOVO	F1:2551	F1:2551	-1.02
F1:2559	1.0	L,C(Carbamidomethylation),N,V,D,H,K,P,N,S,T,K,V,D,K	-0.93	-0.67,-2.27,-1.34,-1.17,-0.85,-0.75,-0.29,-0.00,-1.49,-2.77,-0.78,-0.79,-0.13,-0.44,-0.14	585.632873535156	3.0	DENOVO	F1:2559	F1:2559	-0.93
F1:2563	1.0	H,L,V,Q,G,G,G,S,S,D,T,E,T,D,P,G,E,T,Y,T,K,G,G,Q,S,D,P,S,Y	-1.85	-0.04,-1.24,-2.50,-2.99,-2.38,-2.52,-2.59,-2.37,-2.53,-1.99,-0.62,-1.55,-2.65,-2.11,-0.82,-1.96,-2.66,-0.05,-0.67,-0.81,-2.08,-2.51,-2.51,-2.01,-2.12,-2.54,-0.45,-2.29,-2.24	990.76806640625	3.0	DENOVO	F1:2563	F1:2563	-1.85
F1:2570	1.0	S,S,Q,A,H,G,N,G,K,P,S,N,K,T,V,D,K,K	-0.94	-0.00,-0.05,-2.98,-1.33,-1.74,-1.25,-1.65,-2.62,-0.78,-0.01,-0.36,-0.92,-0.44,-0.89,-0.33,-0.17,-1.35,-0.08	377.402893066406	5.0	DENOVO	F1:2570	F1:2570	-0.94
F1:2581	1.0	K,H,G,S,L,V,D,P,G,S,N,P,D,Y,R	-1.46	-0.23,-2.53,-1.29,-2.12,-1.65,-1.41,-2.03,-1.56,-2.01,-2.66,-0.06,-0.01,-2.85,-0.12,-1.37	411.207000732422	4.0	DENOVO	F1:2581	F1:2581	-1.46
F1:2585	1.0	E,H,M,H,E,A,L,A,E,D,H,Y,T,Q,K	-0.84	-2.38,-1.29,-1.16,-0.38,-0.06,-0.22,-0.57,-1.43,-2.17,-0.79,-0.03,-1.93,-0.06,-0.14,-0.02	460.457855224609	4.0	DENOVO	F1:2585	F1:2585	-0.84
F1:2590	1.0	C(Carbamidomethylation),D,I,V,D,H,K,P,S,N,T,K	-0.72	-0.04,-3.48,-2.24,-0.62,-0.21,-1.05,-0.06,-0.01,-0.11,-0.42,-0.20,-0.20	354.176544189453	4.0	DENOVO	F1:2590	F1:2590	-0.72
F1:2594	1.0	A,P,P,A,L,A,T,Y,L,P,A,A,S,S,P,P,A,P,A,G,L	-1.51	-0.01,-0.80,-2.14,-2.12,-1.96,-2.77,-1.35,-1.74,-1.36,-1.92,-2.33,-2.38,-2.15,-1.63,-0.07,-0.78,-1.27,-0.81,-0.64,-0.40,-3.00	484.016082763672	4.0	DENOVO	F1:2594	F1:2594	-1.51
F1:2603	1.0	S,A,D,G,S,Q,S,A,V,S,E,P,S,E,E,C(Carbamidomethylation),S,D,S,A,E,T,P,P	-1.78	-1.63,-0.53,-2.68,-1.32,-2.32,-1.76,-2.42,-1.15,-0.47,-2.59,-1.71,-1.38,-2.34,-2.34,-2.45,-0.66,-2.75,-0.28,-2.28,-1.87,-2.52,-1.41,-2.84,-0.90	606.746704101563	4.0	DENOVO	F1:2603	F1:2603	-1.78
F1:2607	1.0	H,L,D,Y,H,G,R,P,G,D,Y,P,S,G,D,T,Y,P,A,E,S,E,S,T,Y,R	-1.58	-0.04,-1.30,-2.86,-1.92,-1.93,-2.12,-2.35,-1.15,-2.51,-2.15,-0.62,-2.08,-1.17,-2.55,-2.47,-1.15,-1.11,-0.44,-1.74,-2.32,-1.01,-2.83,-0.29,-0.32,-0.78,-1.76	743.328735351563	4.0	DENOVO	F1:2607	F1:2607	-1.58
F1:2614	1.0	A,N,I,V,N,H,K,P,S,N,T,K,V,D,K	-0.49	-0.01,-0.38,-2.11,-0.20,-1.51,-0.46,-0.28,-0.01,-0.03,-0.18,-0.42,-1.33,-0.11,-0.11,-0.29	416.982208251953	4.0	DENOVO	F1:2614	F1:2614	-0.49
F1:13999	1.0	S,L,A,A,P,S,V,F,L,F,P,Q,S,G,L,M(Oxidation),H,K	-1.11	-1.94,-1.17,-0.19,-0.02,-0.29,-2.37,-0.90,-0.04,-0.96,-0.03,-0.95,-2.41,-2.35,-0.73,-0.83,-1.65,-3.04,-0.01	973.518920898438	2.0	DENOVO	F1:13999	F1:13999	-1.11
F1:14017	1.0	V,P,P,V,S,E,S,S,V,E,L,D,P,E,A,A,N,P,D,S,S,A,A,S,V,R	-1.92	-0.05,-1.81,-0.45,-3.02,-2.47,-2.24,-2.27,-2.21,-2.20,-1.36,-1.70,-2.59,-2.18,-2.38,-2.06,-2.78,-0.31,-1.25,-2.63,-2.17,-2.25,-1.73,-1.25,-2.30,-2.58,-1.74	870.756530761719	3.0	DENOVO	F1:14017	F1:14017	-1.92
F1:14027	1.0				376.314300537109	2.0		F1:14027	F1:14027	
F1:14039	1.0				419.316162109375	2.0		F1:14039	F1:14039	
F1:14051	1.0	P,K,K,P,L,L,T,Q,L,K	-2.2	-0.00,-3.37,-4.03,-0.94,-2.76,-2.46,-2.31,-2.42,-2.55,-1.16	583.386840820313	2.0	DENOVO	F1:14051	F1:14051	-2.2
F1:14056	1.0	L,T,P,A,S,G,G,A,T,T,Y,L,S,Q,E,Y,S,K	-1.76	-0.50,-0.94,-1.32,-1.77,-2.15,-2.37,-2.56,-1.95,-2.29,-2.81,-1.06,-1.10,-2.52,-2.12,-2.35,-1.03,-0.77,-2.10	937.463256835938	2.0	DENOVO	F1:14056	F1:14056	-1.76";

const POINTNOVO_V0_0_1: &str = "feature_id	feature_area	predicted_sequence	predicted_score	predicted_position_score	precursor_mz	precursor_charge	protein_access_id	scan_list_middle	scan_list_original	predicted_score_max
F1:14	1.0	K,K,K,K,A,L,L,G,L,P,A,L,S,G,L,G,G,E,E,E,H,K,K	-2.16	-0.00,-0.32,-2.52,-3.19,-3.02,-2.87,-2.69,-2.87,-2.25,-2.30,-3.06,-2.20,-2.60,-1.99,-2.69,-0.31,-2.15,-2.69,-2.62,-1.17,-1.68,-2.85,-1.69	1216.224487304688	2.0	DENOVO	F1:14	F1:14	-2.16
F1:15	1.0	K,K,K,K,K,K,K,R,L,L,L,S,G,G,S,E,D,G,E,S,S,K	-2.54	-0.00,-3.10,-3.20,-3.22,-3.22,-3.22,-3.22,-2.96,-2.84,-2.84,-2.84,-2.68,-2.73,-2.74,-2.68,-2.44,-1.39,-2.75,-2.48,-0.93,-2.74,-1.69	1216.224487304688	2.0	DENOVO	F1:15	F1:15	-2.54
F1:2535	1.0	C(Carbamidomethylation),D,L,V,N,H,K,P,S,N,T,K,V,D,K	-0.28	-1.53,-0.30,-0.95,-0.08,-0.10,-0.05,-0.01,-0.01,-0.16,-0.40,-0.09,-0.09,-0.11,-0.39,-0.00	439.477081298828	4.0	DENOVO	F1:2535	F1:2535	-0.28
F1:2536	1.0	C(Carbamidomethylation),D,L,V,N,H,K,P,S,N,T,K,V,D,K	-0.70	-2.64,-0.87,-0.97,-0.07,-0.67,-0.41,-0.10,-0.70,-1.15,-0.84,-0.43,-0.78,-0.59,-0.30,-0.00	439.477081298828	4.0	DENOVO	F1:2536	F1:2536	-0.70
F1:2537	1.0	S,S,M(Oxidation),Y,Y,Y,S,D,T,K,N,S,Q,T,R,E,Q,Y,N,S,T,Y,R	-1.07	-0.00,-0.06,-3.51,-1.88,-0.15,-0.18,-1.67,-2.24,-0.33,-0.26,-1.25,-2.26,-2.17,-0.55,-0.75,-2.08,-1.79,-0.89,-0.67,-0.79,-0.44,-0.24,-0.48	722.81494140625	4.0	DENOVO	F1:2537	F1:2537	-1.07
F1:2538	1.0	K,T,C(Carbamidomethylation),L,E,D,P,D,D,E,N,R,M,E,D,L,E,G,G,G,T,S,T,Y,R	-1.51	-0.00,-1.00,-2.37,-1.88,-1.33,-1.33,-1.15,-2.13,-1.33,-2.44,-1.87,-1.03,-1.00,-1.85,-2.36,-0.82,-1.71,-2.45,-1.24,-2.25,-0.76,-1.36,-0.35,-1.54,-2.23	722.81494140625	4.0	DENOVO	F1:2538	F1:2538	-1.51
F1:2539	1.0	M,T,M,F,Y,C(Carbamidomethylation),G,F,H,K,T,V,Q,C(Carbamidomethylation),Y,R,Y,N,S,T,Y,R	-1.05	-0.01,-1.07,-3.40,-1.65,-0.62,-0.29,-1.07,-0.57,-0.50,-0.72,-1.58,-0.75,-0.47,-1.22,-1.71,-2.60,-1.30,-0.74,-0.47,-1.17,-0.30,-0.81	729.072875976563	4.0	DENOVO	F1:2539	F1:2539	-1.05
F1:2540	1.0	Y,S,D,C(Carbamidomethylation),F,Y,A,R,G,F,G,Q,W,Y,T,A,Q,Q,G,T,S,T,Y,R	-1.37	-0.00,-0.90,-2.04,-0.19,-1.94,-1.16,-0.37,-2.18,-1.22,-2.14,-0.92,-1.43,-3.03,-1.37,-1.81,-0.66,-2.23,-2.25,-0.87,-1.01,-1.71,-0.19,-0.81,-2.47	729.072875976563	4.0	DENOVO	F1:2540	F1:2540	-1.37
F1:2541	1.0	G,Y,A,D,C(Carbamidomethylation),Y,C(Carbamidomethylation),L,E,N,A,R,V,D,Y,D,Y,Q,F,T,K,Y,R	-1.21	-0.02,-0.92,-1.60,-0.53,-0.16,-0.18,-0.59,-1.38,-2.80,-0.57,-1.55,-0.62,-2.67,-1.47,-1.64,-1.65,-2.69,-0.86,-0.98,-1.71,-1.46,-0.65,-1.04	743.328918457031	4.0	DENOVO	F1:2541	F1:2541	-1.21
F1:2542	1.0	S,D,Y,G,F,C(Carbamidomethylation),Y,G,Y,T,E,R,K,L,S,S,A,E,N,S,E,S,T,Y,R	-1.49	-2.65,-1.36,-0.67,-2.68,-1.22,-1.67,-0.36,-0.50,-2.38,-0.90,-0.27,-2.81,-2.61,-1.93,-0.73,-2.66,-1.05,-1.74,-1.17,-2.83,-0.82,-1.73,-1.10,-1.32,-0.07	743.328918457031	4.0	DENOVO	F1:2542	F1:2542	-1.49
F1:2546	1.0	C(Carbamidomethylation),D,V,L,N,H,K,P,S,N,T,K,V,D,K	-0.38	-1.54,-0.84,-0.78,-0.72,-0.64,-0.18,-0.07,-0.05,-0.06,-0.36,-0.03,-0.07,-0.03,-0.36,-0.00	351.783843994141	5.0	DENOVO	F1:2546	F1:2546	-0.38
F1:2547	1.0	C(Carbamidomethylation),D,L,V,N,H,K,P,S,N,T,K,V,D,K	-0.32	-2.11,-0.76,-0.81,-0.01,-0.31,-0.05,-0.10,-0.03,-0.07,-0.20,-0.01,-0.03,-0.03,-0.21,-0.00	351.783843994141	5.0	DENOVO	F1:2547	F1:2547	-0.32
F1:2550	1.0	L,A,N,V,N,H,K,P,S,N,T,K,V,D,K	-0.26	-1.31,-0.13,-0.60,-0.76,-0.12,-0.10,-0.01,-0.01,-0.24,-0.28,-0.05,-0.07,-0.05,-0.18,-0.00	416.983551025391	4.0	DENOVO	F1:2550	F1:2550	-0.26
F1:2551	1.0	L,A,N,V,N,H,K,P,S,N,T,K,V,D,K	-0.62	-1.96,-0.33,-0.56,-0.20,-1.55,-0.20,-0.04,-0.28,-1.70,-1.03,-0.29,-0.35,-0.64,-0.22,-0.00	416.983551025391	4.0	DENOVO	F1:2551	F1:2551	-0.62
F1:2552	1.0	G,D,D,G,E,A,Y,G,Q,R,L,C(Carbamidomethylation),D,K,Y,A,G,R,D,N,S,T,Y,R	-1.24	-0.01,-0.97,-1.18,-2.16,-2.16,-2.30,-1.91,-1.97,-1.46,-0.37,-1.00,-1.89,-2.22,-1.12,-1.32,-0.34,-1.85,-0.91,-2.38,-0.41,-0.46,-0.35,-0.27,-0.86	692.5576171875	4.0	DENOVO	F1:2552	F1:2552	-1.24
F1:2553	1.0	G,V,D,F,S,S,E,F,T,S,M,G,H,P,S,T,Y,T,L,E,S,T,N,V,S,S	-1.59	-1.93,-1.98,-1.43,-2.04,-1.00,-1.17,-1.10,-1.06,-1.16,-2.62,-1.84,-2.04,-0.62,-1.87,-2.16,-1.85,-0.83,-2.04,-1.71,-0.94,-2.05,-1.21,-2.49,-2.38,-1.76,-0.13	692.5576171875	4.0	DENOVO	F1:2553	F1:2553	-1.59
F1:2554	1.0	H,G,G,N,C(Carbamidomethylation),Y,Y,H,T,K,S,E,G,R,Y,T,S,G,A,Y,N,S,T,Y,R	-1.32	-0.00,-1.38,-2.28,-3.39,-0.10,-0.63,-0.17,-0.66,-1.20,-1.36,-0.84,-2.53,-1.43,-2.81,-1.42,-0.71,-2.34,-0.64,-2.47,-2.29,-0.42,-0.22,-2.63,-0.29,-0.75	733.071838378906	4.0	DENOVO	F1:2554	F1:2554	-1.32
F1:2559	1.0	L,C(Carbamidomethylation),D,V,N,H,K,P,S,N,T,K,V,D,K	-0.59	-0.71,-1.25,-0.14,-0.20,-1.84,-0.05,-0.01,-0.34,-2.51,-1.04,-0.15,-0.22,-0.14,-0.19,-0.00	585.632873535156	3.0	DENOVO	F1:2559	F1:2559	-0.59
F1:2560	1.0	G,T,D,N,V,N,H,K,R,G,K,S,A,Q,D,K	-1.11	-0.01,-0.75,-2.93,-1.00,-1.65,-1.23,-0.54,-1.03,-1.09,-0.81,-0.60,-0.89,-1.61,-2.78,-0.69,-0.22	585.632873535156	3.0	DENOVO	F1:2560	F1:2560	-1.11
F1:2659	1.0	H,G,G,G,G,Y,S,D,P,D,E,E,T,T,G,M(Oxidation),C(Carbamidomethylation),D,K,N,P,A,C(Carbamidomethylation),Y,W	-1.79	-0.00,-1.47,-2.56,-2.58,-2.81,-1.59,-2.01,-1.74,-1.34,-2.20,-2.61,-2.40,-0.89,-2.46,-0.40,-1.85,-1.51,-1.83,-2.75,-0.88,-2.13,-2.09,-0.60,-1.31,-2.61	705.77099609375	4.0	DENOVO	F1:2659	F1:2659	-1.79";

const BIATNOVO_V0_1: &str = "feature_id	feature_area	predicted_sequence	predicted_score	predicted_position_score	precursor_mz	precursor_charge	protein_access_id	scan_list_middle	scan_list_original	predicted_score_max
F1:14	1.0	Q(Deamidation),W,W,W,Q,L,G,L,V,G,L,V,L,V,L,V,L,L,Y,E	-5.75	-6.29,-5.10,-5.09,-5.82,-5.88,-4.68,-6.96,-5.82,-4.32,-7.10,-6.12,-4.90,-7.30,-5.04,-7.19,-5.08,-5.17,-5.51,-6.77,-4.75	1216.224487304688	2.0	DENOVO	F1:14	F1:14	-5.75
F1:15	1.0	K,W,W,W,W,W,R,E,K,N,P,E,N,P,E,K,E	-5.15	-5.84,-5.48,-4.69,-4.79,-4.85,-5.54,-5.78,-2.96,-6.02,-5.10,-5.64,-4.75,-5.65,-5.56,-3.40,-6.13,-5.45	1216.224487304688	2.0	DENOVO	F1:15	F1:15	-5.15
F1:2535	1.0	V,G,N,P,E,K,E,K,N,P,E,K,E,K,E	-5.34	-5.97,-5.65,-4.93,-6.26,-3.66,-5.94,-4.91,-6.10,-5.47,-5.29,-3.16,-6.08,-4.76,-6.27,-5.72	439.477081298828	4.0	DENOVO	F1:2535	F1:2535	-5.34
F1:2536	1.0	K,S,L,L,N,N,P,E,K,D,N,P,E,K,E	-5.40	-4.54,-4.42,-6.67,-6.83,-5.07,-5.49,-5.59,-3.77,-5.43,-5.92,-5.66,-5.40,-3.83,-6.18,-6.22	439.477081298828	4.0	DENOVO	F1:2536	F1:2536	-5.40
F1:2537	1.0	K,W,W,W,W,W,W,W,W,N,V,E,N,P,E,N,P,E,P,E	-5.40	-5.68,-5.49,-4.77,-4.92,-4.99,-5.03,-5.16,-5.31,-7.11,-6.07,-4.35,-6.05,-5.94,-5.63,-4.79,-5.69,-5.25,-5.15,-5.94,-4.56	722.81494140625	4.0	DENOVO	F1:2537	F1:2537	-5.40
F1:2538	1.0	K,S,C(Carbamidomethylation),Q,M(Oxidation),E,L,Q,M(Oxidation),E,L,Q,M(Oxidation),E,L,Q,M(Oxidation),E,L,F,F,M	-5.63	-4.27,-5.17,-7.02,-6.03,-6.01,-4.82,-4.44,-6.22,-5.84,-4.90,-4.56,-6.28,-5.78,-4.82,-4.88,-6.37,-5.90,-4.72,-5.62,-6.91,-6.31,-6.94	722.81494140625	4.0	DENOVO	F1:2538	F1:2538	-5.63
F1:2539	1.0	Q(Deamidation),W,W,W,W,Q,M(Oxidation),E,L,Q,L,Q,L,Q,M(Oxidation),E,L,Q,M(Oxidation),N(Deamidation),N	-5.54	-6.03,-4.87,-4.89,-5.13,-5.95,-5.86,-6.05,-5.26,-4.49,-6.02,-4.77,-6.02,-4.66,-5.78,-5.62,-5.20,-4.37,-6.20,-5.03,-5.82,-8.29	729.072875976563	4.0	DENOVO	F1:2539	F1:2539	-5.54
F1:2540	1.0	K,W,W,W,W,W,P,E,K,D,N,N,V,E,N,P,E,N,P,E,P,E	-5.40	-5.67,-5.34,-4.68,-4.87,-5.00,-7.01,-6.05,-4.32,-4.97,-6.33,-5.60,-6.29,-4.27,-5.82,-5.41,-5.89,-4.69,-5.54,-5.48,-5.17,-5.75,-4.54	729.072875976563	4.0	DENOVO	F1:2540	F1:2540	-5.40
F1:2541	1.0	S,V,R,V,E,N,N,P,E,N,P,E,N,P,E,N,P,E,N,P,E,N,P,E,P,E	-5.29	-8.11,-5.22,-5.80,-6.10,-6.05,-5.92,-5.91,-5.69,-4.66,-5.19,-5.45,-4.59,-5.19,-5.35,-4.46,-5.16,-5.28,-4.47,-5.12,-5.22,-4.31,-5.26,-4.82,-4.91,-5.15,-4.13	743.328918457031	4.0	DENOVO	F1:2541	F1:2541	-5.29
F1:2542	1.0	P,E,V,P,E,N,P,E,K,D,N,N,N,V,E,N,P,E,N,P,E,N,P,E,P,E	-5.25	-5.95,-5.49,-6.04,-5.99,-4.98,-4.93,-5.21,-3.82,-4.69,-6.77,-6.12,-6.17,-5.94,-4.68,-5.72,-4.88,-5.34,-4.48,-4.84,-5.31,-4.30,-5.01,-4.92,-4.94,-5.67,-4.18	743.328918457031	4.0	DENOVO	F1:2542	F1:2542	-5.25
F1:2546	1.0	A,R,V,D,N,V,R,E,K,D,N,N,V,P,E	-5.52	-7.99,-5.70,-6.14,-5.38,-5.64,-4.07,-5.71,-3.50,-5.49,-6.11,-6.33,-5.94,-4.40,-6.04,-4.30	351.783843994141	5.0	DENOVO	F1:2546	F1:2546	-5.52
F1:2547	1.0	K,S,C(Carbamidomethylation),Q,M(Oxidation),E,L,Q,L,G,S,V,F,L,V	-5.65	-4.40,-4.94,-6.56,-6.03,-5.96,-5.45,-4.54,-6.17,-4.97,-6.20,-5.89,-4.23,-7.04,-6.85,-5.58	351.783843994141	5.0	DENOVO	F1:2547	F1:2547	-5.65
F1:2550	1.0	K,H,R,I,N,P,E,K,E,K,E,K,E	-5.48	-5.03,-5.81,-6.21,-7.72,-5.95,-5.44,-3.08,-6.34,-3.90,-6.24,-4.03,-6.35,-5.13	416.983551025391	4.0	DENOVO	F1:2550	F1:2550	-5.48
F1:2551	1.0	Q(Deamidation),W,W,W,G,L,G,L,V,L,L,Y,E	-5.70	-6.58,-5.01,-5.00,-5.23,-7.30,-4.69,-6.78,-5.97,-4.85,-5.31,-5.56,-6.95,-4.85	416.983551025391	4.0	DENOVO	F1:2551	F1:2551	-5.70
F1:2552	1.0	K,W,W,W,P,E,K,D,N,V,E,N,P,E,N,P,E,N,P,E,P,E	-5.41	-5.67,-5.64,-4.98,-7.11,-6.12,-4.31,-4.93,-6.35,-5.33,-4.62,-5.78,-5.56,-5.85,-4.62,-5.46,-5.82,-4.63,-5.61,-5.39,-5.09,-5.67,-4.50	692.5576171875	4.0	DENOVO	F1:2552	F1:2552	-5.41
F1:2553	1.0	R,V,P,E,N,P,E,N,P,E,N,P,E,N,P,E,N,P,E,N,P,E,P,E	-5.19	-6.05,-6.42,-6.79,-4.80,-5.33,-5.87,-4.81,-5.04,-5.66,-4.71,-4.97,-5.51,-4.66,-5.03,-5.46,-4.49,-4.98,-5.37,-4.37,-5.01,-4.92,-4.83,-5.18,-4.20	692.5576171875	4.0	DENOVO	F1:2553	F1:2553	-5.19
F1:2554	1.0	H,E,K,N,N,N,N,N,N,P,E,K,D,N,N,V,E,N,P,E,N,P,E,P,E	-5.40	-10.59,-3.62,-6.28,-5.28,-6.29,-6.10,-5.92,-5.78,-5.45,-5.73,-3.61,-5.23,-5.73,-5.48,-5.38,-4.29,-5.45,-4.43,-5.74,-4.34,-4.55,-5.34,-4.78,-5.44,-4.23	733.071838378906	4.0	DENOVO	F1:2554	F1:2554	-5.40
F1:2559	1.0	K,L,G,K,E,K,D,N,P,E,N,P,E,K,E	-5.44	-4.61,-5.79,-6.11,-5.95,-4.30,-5.23,-5.88,-6.02,-5.62,-5.16,-5.39,-5.38,-3.87,-6.12,-6.14	585.632873535156	3.0	DENOVO	F1:2559	F1:2559	-5.44
F1:2560	1.0	K,L,G,K,E,K,D,N,P,E,N,P,E,K,E	-5.44	-4.64,-5.81,-6.09,-5.96,-4.29,-5.21,-5.89,-6.04,-5.61,-5.13,-5.37,-5.37,-3.92,-6.15,-6.14	585.632873535156	3.0	DENOVO	F1:2560	F1:2560	-5.44
F1:2660	1.0	K,S,C(Carbamidomethylation),Q,M(Oxidation),E,L,Q,M(Oxidation),E,L,Q,M(Oxidation),N(Deamidation),M,S,C(Carbamidomethylation),Q,M(Oxidation),N(Deamidation),T,M(Oxidation)	-5.79	-4.59,-4.85,-6.96,-6.41,-6.07,-5.14,-4.30,-6.58,-5.87,-5.37,-4.26,-6.37,-4.43,-5.89,-6.37,-4.95,-6.45,-6.08,-4.51,-5.20,-9.10,-7.58	709.518188476563	4.0	DENOVO	F1:2660	F1:2660	-5.79";