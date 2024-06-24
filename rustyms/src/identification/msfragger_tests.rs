#![allow(clippy::missing_panics_doc)]
use std::io::BufReader;

use super::IdentifiedPeptideSource;

use super::{csv::parse_csv_raw, msfragger, IdentifiedPeptide, MSFraggerData};

#[test]
fn msfragger_v21() {
    let reader = BufReader::new(DATA_V21.as_bytes());
    let lines = parse_csv_raw(reader, b'\t', None).unwrap();
    for line in lines.map(std::result::Result::unwrap) {
        println!("{line}");
        let _read: IdentifiedPeptide = MSFraggerData::parse_specific(&line, &msfragger::V21)
            .unwrap()
            .into();
    }
}

#[test]
fn msfragger_detect() {
    let reader = BufReader::new(DATA_V21.as_bytes());
    let lines = parse_csv_raw(reader, b'\t', None).unwrap();
    for line in lines.map(std::result::Result::unwrap) {
        println!("{line}");
        let result = MSFraggerData::parse(&line).unwrap();
        let _read: IdentifiedPeptide = result.0.into();
        assert_eq!(result.1, &msfragger::V21);
    }
}

const DATA_V21: &str = r"Spectrum                                                                 	Spectrum.File                                                      	Peptide                               	Modified sequence                              	Extended.Peptide                                        	Prev.AA 	Next.AA 	Peptide.Length 	Charge 	Retention 	Observed.Mass 	Calibrated.Observed.Mass 	Observed.M.Z 	Calibrated.Observed.M.Z 	Calculated.Peptide.Mass 	Calculated.M.Z 	Delta.Mass 	Expectation  	Hyperscore 	Nextscore 	PeptideProphet.Probability 	Number.of.Enzymatic.Termini 	Number.of.Missed.Cleavages 	Protein.Start 	Protein.End 	Intensity     	Assigned.Modifications                                            	Observed.Modifications 	Purity 	Is.Unique 	Protein                   	Protein.ID 	Entry.Name  	Gene     	Protein.Description                                                    	Mapped.Genes                                  	Mapped.Proteins                                                                                                                         	condition  	group
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.24699.24699.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	ACVFGNEPK                             	AC[751]VFGNEPK                                 	GICFPIVK.ACVFGNEPK.ASDEVPIA                             	K       	A       	             9 	     3 	2113.766  	    1611.9152 	               1611.9067 	    538.3123 	               538.3095 	              1611.9058 	      538.3092 	9e-04      	0.02054707   	    17.249 	   14.328 	                    0.9938 	                          2 	                         0 	           79 	         87 	 126014272    	2C(DB16 (C))                                                      	                       	     0 	true      	sp|Q8TB61|S35B2_HUMAN     	Q8TB61     	S35B2_HUMAN 	SLC35B2  	Adenosine 3'-phospho 5'-phosphosulfate transporter 1                   	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.24859.24859.2 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	ACVFGNEPK                             	AC[751]VFGNEPK                                 	GICFPIVK.ACVFGNEPK.ASDEVPIA                             	K       	A       	             9 	     2 	2126.352  	    1612.9125 	               1612.9115 	    807.4635 	               807.463  	              1611.9058 	      806.9602 	1.0057     	0.03336185   	    15.407 	   11.99  	                    0.9511 	                          2 	                         0 	           79 	         87 	  14438628    	2C(DB16 (C))                                                      	                       	     0 	true      	sp|Q8TB61|S35B2_HUMAN     	Q8TB61     	S35B2_HUMAN 	SLC35B2  	Adenosine 3'-phospho 5'-phosphosulfate transporter 1                   	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.24885.24885.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	ACVFGNEPK                             	AC[751]VFGNEPK                                 	GICFPIVK.ACVFGNEPK.ASDEVPIA                             	K       	A       	             9 	     3 	2128.4214 	    1611.9141 	               1611.9067 	    538.312  	               538.3095 	              1611.9058 	      538.3092 	9e-04      	0.01162117   	    19.239 	   14.355 	                    0.9958 	                          2 	                         0 	           79 	         87 	  48494552    	2C(DB16 (C))                                                      	                       	     0 	true      	sp|Q8TB61|S35B2_HUMAN     	Q8TB61     	S35B2_HUMAN 	SLC35B2  	Adenosine 3'-phospho 5'-phosphosulfate transporter 1                   	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.29547.29547.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	AILPCIK                               	AILPC[751]IK                                   	KPSAIQQR.AILPCIK.GYDVIAQA                               	R       	G       	             7 	     3 	2498.9243 	    1405.9248 	               1405.9152 	    469.6489 	               469.6457 	              1404.9141 	      469.312  	1.0011     	0.118168     	    11.21  	    0     	                    0.8232 	                          2 	                         0 	           62 	         68 	  30299220    	5C(DB16 (C))                                                      	                       	     0 	false     	sp|P60842|IF4A1_HUMAN     	P60842     	IF4A1_HUMAN 	EIF4A1   	Eukaryotic initiation factor 4A-I                                      	EIF4A2, USP35                                 	sp|Q14240|IF4A2_HUMAN, sp|Q9P2H5|UBP35_HUMAN                                                                                            	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.29648.29648.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	AILPCIK                               	AILPC[751]IK                                   	KPSAIQQR.AILPCIK.GYDVIAQA                               	R       	G       	             7 	     3 	2507.1746 	    1404.9186 	               1404.9143 	    469.3135 	               469.312  	              1404.9141 	      469.312  	2e-04      	0.07998345   	    13.888 	    8.34  	                    0.9762 	                          2 	                         0 	           62 	         68 	  27065336    	5C(DB16 (C))                                                      	                       	     0 	false     	sp|P60842|IF4A1_HUMAN     	P60842     	IF4A1_HUMAN 	EIF4A1   	Eukaryotic initiation factor 4A-I                                      	EIF4A2, USP35                                 	sp|Q14240|IF4A2_HUMAN, sp|Q9P2H5|UBP35_HUMAN                                                                                            	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.33383.33383.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	AILPCIK                               	AILPC[779]IK                                   	KPSAIQQR.AILPCIK.GYDVIAQA                               	R       	G       	             7 	     3 	2813.4368 	    1432.9509 	               1432.9463 	    478.6576 	               478.656  	              1432.9453 	      478.6557 	9e-04      	0.083946     	    10.804 	    0     	                    0.9751 	                          2 	                         0 	           62 	         68 	   1130153.5  	5C(DB18 (C))                                                      	                       	     0 	false     	sp|P60842|IF4A1_HUMAN     	P60842     	IF4A1_HUMAN 	EIF4A1   	Eukaryotic initiation factor 4A-I                                      	EIF4A2, USP35                                 	sp|Q14240|IF4A2_HUMAN, sp|Q9P2H5|UBP35_HUMAN                                                                                            	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.33985.33985.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	ALFLLPCISR                            	ALFLLPC[751]ISR                                	PCVVIIAK.ALFLLPCISR.RIARIRRG                            	K       	R       	            10 	     3 	2864.4353 	    1782.1174 	               1782.1111 	    595.0464 	               595.0443 	              1780.1047 	      594.3755 	2.0063     	0.01829153   	    18.451 	   10.272 	                    0.9753 	                          2 	                         0 	          447 	        456 	  11851491    	7C(DB16 (C))                                                      	                       	     0 	true      	sp|Q658P3|STEA3_HUMAN     	Q658P3     	STEA3_HUMAN 	STEAP3   	Metalloreductase STEAP3                                                	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.26237.26237.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	APTFSYGPDGNGFSLGCSK                   	APTFSYGPDGNGFSLGC[751]SK                       	RTTIESFR.APTFSYGPDGNGFSLGCSK.NWRQVFGD                   	R       	N       	            19 	     3 	2234.9429 	    2553.3066 	               2553.3008 	    852.1095 	               852.1075 	              2552.2983 	      851.7734 	1.0024     	0            	    51.206 	   18.487 	                    1      	                          2 	                         0 	          247 	        265 	         0    	17C(DB16 (C))                                                     	                       	     0 	true      	sp|Q5W0Z9|ZDH20_HUMAN     	Q5W0Z9     	ZDH20_HUMAN 	ZDHHC20  	Palmitoyltransferase ZDHHC20                                           	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.24082.24082.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	AQCPIVER                              	AQC[751]PIVER                                  	YAAKRFRK.AQCPIVER.LTNSMMMH                              	K       	L       	             8 	     3 	2065.576  	    1562.9282 	               1562.9207 	    521.9833 	               521.9808 	              1562.9218 	      521.9812 	-0.0011    	0.005384842  	    17.493 	   16.374 	                    0.9959 	                          2 	                         0 	           64 	         71 	   8443916    	3C(DB16 (C))                                                      	                       	     0 	true      	sp|P46782|RS5_HUMAN       	P46782     	RS5_HUMAN   	RPS5     	Small ribosomal subunit protein uS7                                    	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.21551.21551.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	ARPAAAVLCR                            	ARPAAAVLC[751]R                                	ASPPEKPR.ARPAAAVLCR.GPVEPIVF                            	R       	G       	            10 	     3 	1870.3008 	    1675.041  	               1675.0337 	    559.3543 	               559.3518 	              1675.033  	      559.3516 	7e-04      	0.2447654    	    13.945 	   11.075 	                    0.5583 	                          2 	                         1 	           13 	         22 	 191776512    	9C(DB16 (C))                                                      	                       	     0 	true      	sp|Q96NT5|PCFT_HUMAN      	Q96NT5     	PCFT_HUMAN  	SLC46A1  	Proton-coupled folate transporter                                      	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.25432.25432.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	ASLCISTK                              	ASLC[751]ISTK                                  	KITGTANK.ASLCISTK.KEVEKMNK                              	K       	K       	             8 	     3 	2172.4705 	    1470.8951 	               1470.886  	    491.3056 	               491.3026 	              1469.8892 	      490.9703 	0.9968     	0.04118603   	    15.041 	   14.661 	                    0.5564 	                          2 	                         0 	          426 	        433 	   5650528    	4C(DB16 (C))                                                      	                       	     0 	true      	sp|P09874|PARP1_HUMAN     	P09874     	PARP1_HUMAN 	PARP1    	Poly [ADP-ribose] polymerase 1                                         	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.29715.29715.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	AVLFCLSEDK                            	AVLFC[751]LSEDK                                	EEVKKRKK.AVLFCLSEDK.KNIIIEEG                            	K       	K       	            10 	     3 	2512.806  	    1772.0212 	               1772.0122 	    591.681  	               591.678  	              1772.0156 	      591.6791 	-0.0034    	0.2732122    	    13.754 	   10.448 	                    0.734  	                          2 	                         0 	           35 	         44 	   2499257.5  	5C(DB16 (C))                                                      	                       	     0 	true      	sp|P23528|COF1_HUMAN      	P23528     	COF1_HUMAN  	CFL1     	Cofilin-1                                                              	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.25665.25665.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	AVLFCLSEDKK                           	AVLFC[751]LSEDKK                               	EEVKKRKK.AVLFCLSEDKK.NIIIEEGK                           	K       	N       	            11 	     3 	2190.4666 	    1901.1077 	               1901.1052 	    634.7098 	               634.709  	              1900.1106 	      634.3775 	0.9946     	0.001829507  	    21.61  	   19.1   	                    0.8597 	                          2 	                         0 	           35 	         45 	  13520978    	5C(DB16 (C))                                                      	                       	     0 	true      	sp|P23528|COF1_HUMAN      	P23528     	COF1_HUMAN  	CFL1     	Cofilin-1                                                              	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.30771.30771.2 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	CDPDYLR                               	C[751]DPDYLR                                   	RRGIAGIR.CDPDYLR.GAIGRIKV                               	R       	G       	             7 	     2 	2599.2131 	    1528.8337 	               1528.8315 	    765.4241 	               765.423  	              1528.8323 	      765.4234 	-7e-04     	0.03641469   	    13.606 	    9.168 	                    0.9889 	                          2 	                         0 	           45 	         51 	   1738923.9  	1C(DB16 (C))                                                      	                       	     0 	true      	sp|Q8IZR5|CKLF4_HUMAN     	Q8IZR5     	CKLF4_HUMAN 	CMTM4    	CKLF-like MARVEL transmembrane domain-containing protein 4             	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.26533.26533.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	CTTEAEQDIEEEK                         	C[751]TTEAEQDIEEEK                             	IIKKKGYR.CTTEAEQDIEEEK.VEKIEIND                         	R       	V       	            13 	     3 	2258.5857 	    2173.0984 	               2173.089  	    725.3734 	               725.3703 	              2172.0872 	      725.0363 	1.0019     	1.00208e-06  	    31.691 	   21.442 	                    1      	                          2 	                         0 	           88 	        100 	   6883843.5  	1C(DB16 (C))                                                      	                       	     0 	true      	sp|Q8IUW5|RELL1_HUMAN     	Q8IUW5     	RELL1_HUMAN 	RELL1    	RELT-like protein 1                                                    	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.29517.29517.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	DLELTPNSGTLCGSLSGK                    	DLELTPNSGTLC[751]GSLSGK                        	KRKAFVKR.DLELTPNSGTLCGSLSGK.KKKRMMYI                    	R       	K       	            18 	     3 	2496.6162 	    2440.3413 	               2440.3306 	    814.4544 	               814.4508 	              2439.3293 	      814.117  	1.0012     	1.521751e-08 	    37.105 	   20.013 	                    1      	                          2 	                         0 	          323 	        340 	  15977042    	12C(DB16 (C))                                                     	                       	     0 	true      	sp|Q14168|MPP2_HUMAN      	Q14168     	MPP2_HUMAN  	MPP2     	MAGUK p55 subfamily member 2                                           	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.29527.29527.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	DLELTPNSGTLCGSLSGK                    	DLELTPNSGTLC[751]GSLSGK                        	KRKAFVKR.DLELTPNSGTLCGSLSGK.KKKRMMYI                    	R       	K       	            18 	     3 	2497.3723 	    2439.3384 	               2439.3296 	    814.1201 	               814.1171 	              2439.3293 	      814.117  	2e-04      	1.40093e-09  	    36.959 	   17.891 	                    1      	                          2 	                         0 	          323 	        340 	  16202080    	12C(DB16 (C))                                                     	                       	     0 	true      	sp|Q14168|MPP2_HUMAN      	Q14168     	MPP2_HUMAN  	MPP2     	MAGUK p55 subfamily member 2                                           	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.33219.33219.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	DLELTPNSGTLCGSLSGK                    	DLELTPNSGTLC[779]GSLSGK                        	KRKAFVKR.DLELTPNSGTLCGSLSGK.KKKRMMYI                    	R       	K       	            18 	     3 	2800.523  	    2468.3755 	               2468.3638 	    823.7991 	               823.7952 	              2467.3606 	      823.4608 	1.0031     	1.56772e-09  	    39.677 	   17.93  	                    1      	                          2 	                         0 	          323 	        340 	   3411319.8  	12C(DB18 (C))                                                     	                       	     0 	true      	sp|Q14168|MPP2_HUMAN      	Q14168     	MPP2_HUMAN  	MPP2     	MAGUK p55 subfamily member 2                                           	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
20240206_EX1_UM3_6579035_SA_EXT00_SS16_pH7p5-single_res30k.27348.27348.3 	D:\02-February\6579035\SS16\MSFragger3\4_res30k_1\interact.pep.xml 	DMEDSVGLDCK                           	DMEDSVGLDC[751]K                               	EDIARISK.DMEDSVGLDCK.RYIGITVA                           	K       	R       	            11 	     3 	2321.7913 	    1858.9524 	               1858.9421 	    620.6581 	               620.6546 	              1858.9419 	      620.6546 	2e-04      	0.001836538  	    21.444 	    0     	                    0.9976 	                          2 	                         0 	           98 	        108 	   7916283.5  	10C(DB16 (C))                                                     	                       	     0 	true      	sp|Q8TAA9|VANG1_HUMAN     	Q8TAA9     	VANG1_HUMAN 	VANGL1   	Vang-like protein 1                                                    	                                              	                                                                                                                                        	4_res30k_1 	4_res30k
";
