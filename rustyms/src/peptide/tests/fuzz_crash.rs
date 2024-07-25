use crate::parse_test;

parse_test!(ne "ESNCe/", fuzz_1);
parse_test!(ne "/", fuzz_2);
parse_test!(ne "[XLMOd:02001#XL1]", fuzz_3);
parse_test!(ne "[XLOOD:02001#X1]", fuzz_4);
parse_test!(ne "<[TMT6plex]@K,N-term>AtPEILTCNSIGCL/", fuzz_5);
parse_test!(ne "-", fuzz_6);
parse_test!(ne "Y(", fuzz_7);
parse_test!(ne "((", fuzz_8);
parse_test!(ne "()", fuzz_9);
parse_test!(ne "///", fuzz_10);
parse_test!(ne "S///////", fuzz_11);
parse_test!(ne "S//////////(", fuzz_12);
parse_test!(ne "[XLMO01#XL1]", fuzz_13);
parse_test!(ne "SrS/////////", fuzz_14);
parse_test!(ne "S/////////0?", fuzz_15);
parse_test!(ne "[XLMOD:07001#XL1]", fuzz_16);
parse_test!(ne "S//////////ETKvXLSES/0/", fuzz_17);
parse_test!(ne "S///////////////7777777777777777777777777777////E", fuzz_18);
parse_test!(ne "S//Q+++++++/", fuzz_19);
parse_test!(ne "-[Phospho]ddddd(", fuzz_20);
parse_test!(ne "KSPEK/3[1]SE1#XX+2Na", fuzz_21);
parse_test!(ne
    "KgPEKME/3[+222222222222222222222222222222Na1#XX1]SE",
    fuzz_22
);
parse_test!(ne "S//QES/0/0+K/3EK/3EK/3[1]f5[+", fuzz_23);
parse_test!(ne "S//+++/", fuzz_24);
parse_test!(ne "SEK[XLMOD:001#XQ1]SESPE/", fuzz_25);
parse_test!(ne "[[ZPX[[[[[[1]k]]]]]]]", fuzz_26);
parse_test!(ne "[Formula:[09F6]^", fuzz_27);
parse_test!(ne "[Formulaaaaaaaaaaaaaaaaaaaaaa:[13C2][#XL1]SESho]", fuzz_28);
parse_test!(ne "[Formula:[13C2][/EMEVTK[#XL1]", fuzz_29);
parse_test!(ne "[Formul[:\\13C8]^", fuzz_30);
parse_test!(ne "[Fula:[13C23C2]", fuzz_31);
parse_test!(ne "[FoFormula:[13C[EMEVEESPEK/3[+2Na+,+H+]", fuzz_32);
parse_test!(ne "[Formula0m//13CC2]", fuzz_33);
parse_test!(ne "[F][XLMOD:02001#XL1]", fuzz_34);
parse_test!(ne "K[1#Phospho]ddKphoK/", fuzz_35);
parse_test!(ne "K[1#Phospho]tdKSPEK/", fuzz_36);
parse_test!(ne "SEQUEN[Formula:[13C2][12C-2]H2N]CENCE//EME//E/", fuzz_37);
parse_test!(ne "<[TMT6plex]@K,N-term>ddGddd//E/////", fuzz_38);
parse_test!(ne
    "<[TMT6plex]@K,N-term>dSEK[XLMOD:020#XXXXXXLLLLLLLLLL1]UENCE///////////////////3///////",
    fuzz_39
);
parse_test!(ne
    "[NCE//EMEVTK[XLMO[NCE//EMEVTK[XLMOD:02001UUUUUUJ[phlspho]",
    fuzz_40
);
parse_test!(ne "[+2EMEVEENa1#Xa+, 
2SPEK/3[+2N+H+]", fuzz_41);
parse_test!(ne "[PhosK[phlspho]d2NaE#Xa+, 
+H+]", fuzz_42);
parse_test!(ne
    "UE+[Formula:}13Cd [12C-2UE+[Formula:}13Cd [1}13Cd [12C-2BE2NH2N]",
    fuzz_43
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]@K,N-terMT6plex]@K,N-term]@K,N-termS//m>S/",
    fuzz_44
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]@K,N-terMT6plex]@<,N-term]@K,N-term>S//m>S/",
    fuzz_45
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]@K,N-terMT6@K,N-term]@K,N-terMT6plexS//m>S/",
    fuzz_46
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]  KK[Phos -tWrMT6]lex]@K,N-term]@K,N-termY//m>S/",
    fuzz_47
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]@Kpho]erMT6plex<@K,N-term]@K,N-term>S//m>S/",
    fuzz_48
);
parse_test!(ne
    "<[TMT6plex]@K,N-term]@K,N-terMT6plex]ex]@K,C-term>S//Q+++++++//?",
    fuzz_49
);
parse_test!(ne "EMEVEESPEK/3[+5555555555555555555555555Na1#Xa+, ]", fuzz_50);
parse_test!(ne "EMEVEESPEK/3[+26255555555555555555[-hospho]dAddd]", fuzz_51);
parse_test!(ne "DMEVEESPEK/3[+2,2555555555555555555555dNa1#Xa+, ]", fuzz_52);
parse_test!(ne "EMEV+[X://Q++[XS//Q+++++++////////R22222/5NXa+, ]", fuzz_53);
parse_test!(ne "[Formula:   CE[Formula:   C2]", fuzz_54);
parse_test!(ne "EK/3[+2Na+,+H+]PAK/3K/", fuzz_55);
parse_test!(ne "SEK[XLMOD:02001#XL1]UENCE//////////////MEV//E///", fuzz_56);
parse_test!(ne
    "EMEVESEQU//Q++[X://Q++[XS//Q+++++++////////R+/////0//R:/201#X2]",
    fuzz_57
);
parse_test!(ne
    "EMEVESEQU//Q++[j://^++[XS//Q+++++++/3[+2NSRa//////////////R+/R:0201#X2]",
    fuzz_58
);
parse_test!(ne "[TMT6plex]", fuzz_59);
parse_test!(ne "[pphlpho]", fuzz_60);
parse_test!(ne "SSSHSSSSSSSM+S+NPS+S+NPSS+?", fuzz_61);
parse_test!(ne "SESPEK/3[+2N1111111111N1111111111										11111a+,2N]CCE", fuzz_62);
parse_test!(ne "?", fuzz_63);
parse_test!(ne "[TMT6plex]?", fuzz_64);
parse_test!(ne "[]?", fuzz_65);
parse_test!(ne "{}{}", fuzz_66);
parse_test!(ne "EMEVEESH+?", fuzz_67);
parse_test!(ne "EM/3A+H+?", fuzz_68);
parse_test!(ne "{}{}{6}", fuzz_69);
parse_test!(ne "EM/3[+SPEKEVVEE2Na+,+VW( EE3Na+,+HPFKAUBE2Nb+,+HPFKAUBE2`{b,++]", fuzz_70);
parse_test!(ne "EM/3[+SPEKEVEE2Na+,+VEE2Na+,+HPFKAUBE2Nb+,+HPFKAUBW( E2Nb+,++]", fuzz_71);
parse_test!(ne "EM/3[+SPEKEVEE2Na+,+VEE2Na+,+HPFKAUBW( E2Nb+,+HPFKAMOD:0200L11[X]2|||||UBE2Nb++++]", fuzz_72);
parse_test!(ne "{}", fuzz_73);
parse_test!(ne "SEGUEN[Formula:[13C2][12C-2]H W$Ґq]", fuzz_74);
parse_test!(ne "SEGUEN[Formula:[13C2][12C-2]H 555555555W$Ґ555W]", fuzz_75);
parse_test!(ne "SEGUEN[Formula:[13C2][12C-2]Hǂ> 555555555555q]", fuzz_76);
parse_test!(ne "EMEVEESPEK/3[+2W$à+H+]", fuzz_77);
parse_test!(ne "EEMEVEESPEK/3[+2NW$àa+,H+]", fuzz_78);
parse_test!(ne "EMA+HBSPEK/3[+2N<ǂ>A5HAa+,+H+]", fuzz_79);
parse_test!(ne "E/3[+2NE+,+sssss<ǂ>A5HAsssssssssssssssssssssssssssssssssssssssssssssssssssssH+][U:Phossssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssspho]?", fuzz_80);
parse_test!(ne "EMEVEESPEK/3[+2N  皚+]", fuzz_81);
parse_test!(ne "EMEVEVEESPEKrmEESPEJ/3[+ᚚ2/52Na000N00000000000000Na000000Na+,+H+]", fuzz_82);
parse_test!(ne "Z+[13B2]^12777777777777", fuzz_83);
parse_test!(ne "EMEVEESPEK/3[+2Nᛚ-teoedMEVEESPEK/3[+2Nᛚ-teoedAa+,+H+]Aa+,+H+]", fuzz_84);
parse_test!(ne "FGCICNDBXTDOALBAHJFT/3[+2NaSPEKEVE00W$Ґ00000000000000110Na+,+H+]", fuzz_85);
parse_test!(ne "EMEVEESPEKEVEESPEK/3[+[13C2][12C-2]HW$à2Neeedeeeeeeeeeea+,+H+]/", fuzz_86);
parse_test!(ne "HMEVEESPEKEVEESPEK/3Z+[13C2]^-2E-2]H2Na/,+H+]/3[", fuzz_87);
parse_test!(ne "EVEESEVEESPEK/3[+2NaE33W$àEKKB+,+)+]PEK/3[+2NaE33W$àEKKB+,+)+]", fuzz_88);
parse_test!(ne "FJAVEXBDVGDGEYAH/3[+2NaSPSPEKEVE0EKEVE00Na+,+H+01,010010Na+,+H+11,000000100EKEVE00Na+,+H+00,11000001E00Na+,+H+10,010001N<ǂ>a+,+H+11,100000111EKEVE00Na+,*H+00,11000D10]", fuzz_89);
parse_test!(ne "BYNDUCDBCEQUEEeEVEESPEK/3[+2Na+,+H+][i[Oxidatinn^on]^-2]SK[#XL1]PEKCE//M[XLMOR:02001LT////QBT", fuzz_90);
parse_test!(ne "EMECIDWDDMJLOAIAPCFQRBO/3[+2EEEEE0000000000000000+,+HNaSPEKEVBFEEEEEEEEEEEEW$à000000+,+HNaSPEKEVEEEEEEEEEEEEEEEEE0000000000+,+HNaSPEKEVBFEEEEEE8EEEEEEE000000/00+,+H+]0", fuzz_91);
parse_test!(ne "EMEVESPEK/3[+2eW$Ґa+,+HH+]", fuzz_92);
parse_test!(ne "<d><d>NEKEVEESPEMEVEESPEK/3[8P  皚H+]/333333[+2N", fuzz_93);
parse_test!(ne "EMEgEgEE//V//HMEgMEVRAL[Formula:[13C2][12C-2]N4EEMEVEESPEKEEMEVENSPSEM[:Oyidatiom3<cEVENSPSEM[:Pxidatiom3<cvU:Oxidatiom3<cn]W$à{+2Na+,+B+N+K,1,c-term>A+H+MEVSPnU:Oxidatiom3<cn]{+2Na+,+B+N+g,K,c-ter5]WEN[Formula:W1SSSSSSSSSSSSSSS-2EHEVRAL[Formula:[13C2[12C-2]N45]WEN[Formula:W1SSSSSSSSSSSSSSS", fuzz_94);
parse_test!(ne "jTRA[Formula:[1SB2]ymEMK[X:[07C     ---,  򆼡,-----   ]hZyvLMOD:0203[13AU]EVTKSESPgK/31 2Na+    ---,  򆼡,-----   ]hZyvLMeeeee~eeeeeeeeeeeeMeeee eeILTClSIGDE+K+K[U:iTidatioOxidation]nQEVNhxiL1]EVTKSESPrK/", fuzz_95);
parse_test!(ne "EEiTRA[formula:[13C2]KKKKKKKKKKKKKKKKtion]hZA[Formula:o]iTRA[formula:[13333333333333333333333333333333333333333333333333333333333333iTRAQ4plex]-EM[U:Oatdxi33333333333333EMM[U:QdatioɒLdation]nQEVNhS[xxidati]n]EEEM33333333333333333333333333333C2]KKKn]hZ", fuzz_96);
parse_test!(ne "<d>", fuzz_97);
parse_test!(ne "<[TMT2plex]@K,N-Term<Rf>_ZSU\\`R>", fuzz_98);
parse_test!(ne "<[TMT6plex]@K,N-term>", fuzz_99);
parse_test!(ne "<[TMT6plex]@K,N-term>", fuzz_100);
parse_test!(ne "<[TMT6plex]@K,N-termKATPEILTCNSIGoooooooooooooooooooooooCL>", fuzz_101);
parse_test!(ne "<[TMt6plex]@K,N-termK,N-term)A)AqaqIRORXO6plATPEILqq,N-termK,N-term>", fuzz_102);
parse_test!(ne "<111TC>", fuzz_103);
parse_test!(ne "<d>", fuzz_104);
parse_test!(ne "<[TMT6plex]@K,N>", fuzz_105);
parse_test!(ne "<d>", fuzz_106);
parse_test!(ne "<111TC><100Y>", fuzz_107);
parse_test!(ne "<[TMt6plex]@K,N-termK,c-termK,K,N-termK,c-termKNr,c-term<111TC>TEEESSSSSSSSSSSSSSSSSSSSSSSSSSSSSSatLoln->", fuzz_108);
parse_test!(ne "<111TC><111TE><111TE><111TC>", fuzz_109);
parse_test!(ne "<11c>", fuzz_110);
parse_test!(ne "<111nB>", fuzz_111);
parse_test!(ne "<111SB><111SB><111SB><111SB><111SB>", fuzz_112);
parse_test!(ne "<111TC><111TC><111TC><11c><111TC><111TC><111TC><111TC><11c><111TC><111TC><111TC><111TC><111TC><111TC><11c><111TC><111TC><111TC><111TC><10b><111TC><111TC><111TC>", fuzz_113);
parse_test!(ne "<[TMT6plex]@L,N-term<lI><lOlex>TC>", fuzz_114);
parse_test!(ne "<[TMT6plex]@K,k>", fuzz_115);
parse_test!(ne "<[]@K,k>", fuzz_116);
parse_test!(ne "<[TMt6plex]@m>", fuzz_117);
parse_test!(ne "<[TMt6plex]@K,N-termK,N-term,N-term4F:K,r,N-term>", fuzz_118);
parse_test!(ne "[Glycan:SSSFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDEC5FFFFFFF7777777777777777777SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS767777FFFFFDEC5HexS0FFFFFFFFFFFFFFFFSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDEC5FFFFDEC5FFFFFFF7777777777777777777SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS767777FFFFFDEC5HexS0FFFFFFFFFFFFFFFFSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDEC5FFFFFFF7777777777777777777SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSnSSSSSSSSSSSSSSSSSSSSSSSSSSSS767777FFFFFDECFFFFFFF7777777777777777777SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSnSSSSSSSSSSSSSSSSSSSSSSSSSSSS767777FFFFFDEC5FFFFDEC5FFFFFF0]?EMHHHHHHHHFFFFFFFFFNE", fuzz_119);
parse_test!(ne "()[Glycan:H0]AA/LAA)[Glycan:H0](LAA)AA(LAA)[Glycan:H0]AA(LAA)[Glycan:H0](LAA)-[Me[For:[12C-P]4?Q[[Glycan:HexSHHSHgHHHEPP-[Me[", fuzz_120);
parse_test!(ne "<[TMt6plex]@K,N-termK,c-termK,K,N-termK,c-termKNrmG,c-termK,K,N-termK,c-termoK,K,K,N-termKKeq,N-termK,c-termK,K,N-termK,c,K,N-term-tertK,N-term>", fuzz_121);
parse_test!(ne "<[TMt6plex]@K,N-termerm,N-term4F:K,r,N-termcan:K,m,N-term4pxrm,C-term4@111nB><111nB><111nB><111nB><111nB><111nB><111nB><111nB><111nB><111nB>", fuzz_122);
parse_test!(ne "<[TMT2plex]@K,m>", fuzz_123);
parse_test!(ne "<[TMT6plex]@K,C-termn:H>", fuzz_124);
parse_test!(ne "<d><d>", fuzz_125);
parse_test!(ne "<[TMT6plex]@K,N-term>", fuzz_126);
parse_test!(ne "()[3][3]", fuzz_127);
