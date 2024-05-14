use rand::distributions::{Distribution, Standard};

use crate::{
    glycan::{BaseSugar, GlycanStructure, GlycanSubstituent, MonoSaccharide},
    modification::SimpleModification,
    system::{dalton, Mass, OrderedMass},
    Element, MolecularFormula,
};

impl Distribution<SimpleModification> for Standard {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> SimpleModification {
        match rng.gen_range(0..=3) {
            0 => SimpleModification::Mass(rng.gen()),
            1 => SimpleModification::Formula(rng.gen()),
            2 => {
                let mut glycans = Vec::new();
                for _ in 0..rng.gen_range(0..32) {
                    glycans.push(rng.gen());
                }
                SimpleModification::Glycan(glycans)
            }
            3 => SimpleModification::GlycanStructure(rng.gen()),
            _ => todo!(),
        }
    }
}

impl Distribution<GlycanStructure> for Standard {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> GlycanStructure {
        let mut branches = Vec::new();
        for _ in 0..rng.gen_range(0..=2) {
            branches.push(rng.gen());
        }
        GlycanStructure::new(rng.gen(), branches)
    }
}

impl Distribution<MonoSaccharide> for Standard {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> MonoSaccharide {
        let mut substituents = Vec::new();
        for _ in 0..rng.gen_range(0..5) {
            substituents.push(rng.gen());
        }
        MonoSaccharide::new(rng.gen(), &substituents)
    }
}

impl Distribution<MolecularFormula> for Standard {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> MolecularFormula {
        let mut formula = MolecularFormula::default();
        for _ in 0..rng.gen_range(0..32) {
            let element: Element = rng.gen();
            let isotope = rng.gen();
            if element.is_valid(isotope) {
                let _ = formula.add((element, isotope, rng.gen()));
            }
        }
        formula
    }
}

impl Distribution<GlycanSubstituent> for Standard {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> GlycanSubstituent {
        match rng.gen_range(1..=44) {
            1 => GlycanSubstituent::Acetimidoyl,
            2 => GlycanSubstituent::Acetyl,
            3 => GlycanSubstituent::AcetylAlanyl,
            4 => GlycanSubstituent::AcetylGlutaminyl,
            5 => GlycanSubstituent::Acid,
            6 => GlycanSubstituent::Alanyl,
            7 => GlycanSubstituent::Alcohol,
            8 => GlycanSubstituent::Amino,
            9 => GlycanSubstituent::Aric,
            10 => GlycanSubstituent::CargoxyEthylidene,
            11 => GlycanSubstituent::Deoxy,
            12 => GlycanSubstituent::Didehydro,
            13 => GlycanSubstituent::DiHydroxyButyryl,
            14 => GlycanSubstituent::DiMethyl,
            15 => GlycanSubstituent::DiMethylAcetimidoyl,
            16 => GlycanSubstituent::DiMethylGlyceryl,
            17 => GlycanSubstituent::Element(rng.gen()),
            18 => GlycanSubstituent::Ethanolamine,
            19 => GlycanSubstituent::EtOH,
            20 => GlycanSubstituent::Formyl,
            21 => GlycanSubstituent::Glyceryl,
            22 => GlycanSubstituent::Glycolyl,
            23 => GlycanSubstituent::Glycyl,
            24 => GlycanSubstituent::HydroxyButyryl,
            25 => GlycanSubstituent::HydroxyMethyl,
            26 => GlycanSubstituent::Lac,
            27 => GlycanSubstituent::Lactyl,
            28 => GlycanSubstituent::Methyl,
            29 => GlycanSubstituent::MethylAcetimidoyl,
            30 => GlycanSubstituent::MethylGlutamyl,
            31 => GlycanSubstituent::NAcetyl,
            32 => GlycanSubstituent::NDiMe,
            33 => GlycanSubstituent::NFo,
            34 => GlycanSubstituent::NGlycolyl,
            35 => GlycanSubstituent::OCarboxyEthyl,
            36 => GlycanSubstituent::PCholine,
            37 => GlycanSubstituent::Phosphate,
            38 => GlycanSubstituent::Pyruvyl,
            39 => GlycanSubstituent::Suc,
            40 => GlycanSubstituent::Sulfate,
            41 => GlycanSubstituent::Tauryl,
            42 => GlycanSubstituent::Ulo,
            43 => GlycanSubstituent::Ulof,
            _ => GlycanSubstituent::Water,
        }
    }
}

impl Distribution<BaseSugar> for Standard {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> BaseSugar {
        match rng.gen_range(0..=9) {
            0 => BaseSugar::None,
            1 => BaseSugar::Sugar,
            2 => BaseSugar::Triose,
            3 => BaseSugar::Tetrose(None),
            4 => BaseSugar::Pentose(None),
            5 => BaseSugar::Hexose(None),
            6 => BaseSugar::Heptose(None),
            7 => BaseSugar::Octose,
            8 => BaseSugar::Nonose,
            _ => BaseSugar::Decose,
        }
    }
}

impl Distribution<Element> for Standard {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> Element {
        Element::try_from(rng.gen_range(Element::Electron as usize..Element::Og as usize)).unwrap()
    }
}

impl Distribution<OrderedMass> for Standard {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> OrderedMass {
        Mass::new::<dalton>(rng.gen_range(f64::MIN..f64::MAX)).into()
    }
}
