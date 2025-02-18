use super::{
    BaseSugar, GlycanStructure, GlycanSubstituent, HeptoseIsomer, HexoseIsomer, MonoSaccharide,
    PentoseIsomer,
};

impl MonoSaccharide {
    fn get_shape(&self) -> (Shape, Colour, String) {
        let modifications = if self.furanose {
            "f".to_string()
        } else {
            String::new()
        };
        match (&self.base_sugar, self.substituents.as_slice()) {
            (BaseSugar::Hexose(isomer), []) => {
                let (s, c) = match isomer {
                    None => (Shape::Circle, Colour::White),
                    Some(HexoseIsomer::Glucose) => (Shape::Circle, Colour::Blue),
                    Some(HexoseIsomer::Mannose) => (Shape::Circle, Colour::Green),
                    Some(HexoseIsomer::Galactose) => (Shape::Circle, Colour::Yellow),
                    Some(HexoseIsomer::Gulose) => (Shape::Circle, Colour::Orange),
                    Some(HexoseIsomer::Altrose) => (Shape::Circle, Colour::Pink),
                    Some(HexoseIsomer::Allose) => (Shape::Circle, Colour::Purple),
                    Some(HexoseIsomer::Talose) => (Shape::Circle, Colour::LightBlue),
                    Some(HexoseIsomer::Idose) => (Shape::Circle, Colour::Brown),
                    Some(HexoseIsomer::Psicose) => (Shape::Pentagon, Colour::Pink),
                    Some(HexoseIsomer::Fructose) => (Shape::Pentagon, Colour::Green),
                    Some(HexoseIsomer::Sorbose) => (Shape::Pentagon, Colour::Orange),
                    Some(HexoseIsomer::Tagatose) => (Shape::Pentagon, Colour::Yellow),
                };
                (s, c, modifications)
            }
            (BaseSugar::Hexose(isomer), [GlycanSubstituent::NAcetyl]) => {
                let c = match isomer {
                    Some(HexoseIsomer::Glucose) => Colour::Blue,
                    Some(HexoseIsomer::Mannose) => Colour::Green,
                    Some(HexoseIsomer::Galactose) => Colour::Yellow,
                    Some(HexoseIsomer::Gulose) => Colour::Orange,
                    Some(HexoseIsomer::Altrose) => Colour::Pink,
                    Some(HexoseIsomer::Allose) => Colour::Purple,
                    Some(HexoseIsomer::Talose) => Colour::LightBlue,
                    Some(HexoseIsomer::Idose) => Colour::Brown,
                    Some(_) | None => Colour::White,
                };
                (Shape::Square, c, modifications)
            }
            (BaseSugar::Hexose(isomer), [GlycanSubstituent::Amino]) => {
                let c = match isomer {
                    Some(HexoseIsomer::Glucose) => Colour::Blue,
                    Some(HexoseIsomer::Mannose) => Colour::Green,
                    Some(HexoseIsomer::Galactose) => Colour::Yellow,
                    Some(HexoseIsomer::Gulose) => Colour::Orange,
                    Some(HexoseIsomer::Altrose) => Colour::Pink,
                    Some(HexoseIsomer::Allose) => Colour::Purple,
                    Some(HexoseIsomer::Talose) => Colour::LightBlue,
                    Some(HexoseIsomer::Idose) => Colour::Brown,
                    Some(_) | None => Colour::White,
                };
                (Shape::CrossedSquare, c, modifications)
            }
            (BaseSugar::Hexose(isomer), [GlycanSubstituent::Acid]) => {
                let c = match isomer {
                    Some(HexoseIsomer::Glucose) => Colour::Blue,
                    Some(HexoseIsomer::Mannose) => Colour::Green,
                    Some(HexoseIsomer::Galactose) => Colour::Yellow,
                    Some(HexoseIsomer::Gulose) => Colour::Orange,
                    Some(HexoseIsomer::Altrose) => Colour::Pink,
                    Some(HexoseIsomer::Allose) => Colour::Purple,
                    Some(HexoseIsomer::Talose) => Colour::LightBlue,
                    Some(HexoseIsomer::Idose) => Colour::Brown,
                    Some(_) | None => Colour::White,
                };
                (Shape::DividedDiamond, c, modifications)
            }
            (BaseSugar::Hexose(isomer), [GlycanSubstituent::Deoxy]) => {
                let c = match isomer {
                    Some(HexoseIsomer::Glucose) => Colour::Blue,
                    Some(HexoseIsomer::Mannose) => Colour::Green,
                    Some(HexoseIsomer::Galactose) => Colour::Red,
                    Some(HexoseIsomer::Gulose) => Colour::Orange,
                    Some(HexoseIsomer::Altrose) => Colour::Pink,
                    Some(HexoseIsomer::Talose) => Colour::LightBlue,
                    Some(_) | None => Colour::White,
                };
                (Shape::Triangle, c, modifications)
            }
            (BaseSugar::Hexose(isomer), [GlycanSubstituent::Deoxy, GlycanSubstituent::NAcetyl]) => {
                let c = match isomer {
                    Some(HexoseIsomer::Glucose) => Colour::Blue,
                    Some(HexoseIsomer::Mannose) => Colour::Green,
                    Some(HexoseIsomer::Galactose) => Colour::Red,
                    Some(HexoseIsomer::Altrose) => Colour::Pink,
                    Some(HexoseIsomer::Talose) => Colour::LightBlue,
                    Some(_) | None => Colour::White,
                };
                (Shape::DividedTriangle, c, modifications)
            }
            (BaseSugar::Hexose(isomer), [GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy]) => {
                let c = match isomer {
                    Some(HexoseIsomer::Glucose) => Colour::Blue,
                    Some(HexoseIsomer::Mannose) => Colour::Green,
                    Some(HexoseIsomer::Galactose) => Colour::Orange,
                    Some(HexoseIsomer::Altrose) => Colour::Pink,
                    Some(HexoseIsomer::Allose) => Colour::Purple,
                    Some(HexoseIsomer::Talose) => Colour::LightBlue,
                    Some(_) | None => Colour::White,
                };
                (Shape::Rectangle, c, modifications)
            }
            (BaseSugar::Pentose(isomer), []) => (
                Shape::Star,
                match isomer {
                    None | Some(PentoseIsomer::Xylulose) => Colour::White,
                    Some(PentoseIsomer::Arabinose) => Colour::Green,
                    Some(PentoseIsomer::Lyxose) => Colour::Yellow,
                    Some(PentoseIsomer::Xylose) => Colour::Orange,
                    Some(PentoseIsomer::Ribose) => Colour::Pink,
                },
                modifications,
            ),
            (
                BaseSugar::Nonose,
                [GlycanSubstituent::Acid, GlycanSubstituent::Amino, GlycanSubstituent::Amino, GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
            ) => (Shape::FlatDiamond, Colour::White, modifications), // This does not detect which exact isomer it is, as that information is currently not tracked
            (BaseSugar::Nonose, mods)
                if mods.contains(&GlycanSubstituent::Amino)
                    && mods.contains(&GlycanSubstituent::Acid) =>
            {
                (
                    Shape::Diamond,
                    match mods {
                        [GlycanSubstituent::Amino, GlycanSubstituent::Deoxy, GlycanSubstituent::Acid] => {
                            Colour::Red
                        }
                        [GlycanSubstituent::Amino, GlycanSubstituent::Acid] => Colour::Brown,
                        [GlycanSubstituent::Amino, GlycanSubstituent::Acetyl, GlycanSubstituent::Acid] => {
                            Colour::LightBlue
                        }
                        [GlycanSubstituent::Amino, GlycanSubstituent::Glycolyl, GlycanSubstituent::Acid] => {
                            Colour::Purple
                        }
                        _ => Colour::White,
                    },
                    modifications,
                )
            }
            (
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                [GlycanSubstituent::Deoxy, GlycanSubstituent::Amino, GlycanSubstituent::Amino],
            ) => (Shape::Hexagon, Colour::Blue, modifications),
            (BaseSugar::Heptose(Some(HeptoseIsomer::GlyceroMannoHeptopyranose)), []) => {
                (Shape::Hexagon, Colour::Green, modifications)
            }
            (BaseSugar::Octose, [GlycanSubstituent::Deoxy, GlycanSubstituent::Acid]) => {
                (Shape::Hexagon, Colour::Yellow, modifications)
            }
            (
                BaseSugar::Heptose(None),
                [GlycanSubstituent::Deoxy, GlycanSubstituent::Acid, GlycanSubstituent::Acid],
            ) => (Shape::Hexagon, Colour::Orange, modifications),
            (
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                [GlycanSubstituent::NAcetyl, GlycanSubstituent::OCarboxyEthyl],
            ) => (Shape::Hexagon, Colour::Purple, modifications),
            (
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                [GlycanSubstituent::NGlycolyl, GlycanSubstituent::OCarboxyEthyl],
            ) => (Shape::Hexagon, Colour::LightBlue, modifications),
            (
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                &[GlycanSubstituent::Amino, GlycanSubstituent::OCarboxyEthyl],
            ) => (Shape::Hexagon, Colour::Brown, modifications),
            _ => (Shape::Hexagon, Colour::White, modifications),
        }
    }
}

enum Colour {
    White,
    Blue,
    Green,
    Yellow,
    Orange,
    Pink,
    Purple,
    LightBlue,
    Brown,
    Red,
}

impl Colour {
    /// Represented as percentages 0..=100
    fn cmyk(&self) -> [u8; 4] {
        match self {
            Self::White => [0, 0, 0, 0],
            Self::Blue => [100, 50, 0, 0],
            Self::Green => [100, 0, 100, 0],
            Self::Yellow => [0, 15, 100, 0],
            Self::Orange => [0, 65, 100, 0],
            Self::Pink => [0, 47, 24, 0],
            Self::Purple => [38, 88, 0, 0],
            Self::LightBlue => [41, 5, 3, 0],
            Self::Brown => [32, 48, 76, 13],
            Self::Red => [0, 100, 100, 0],
        }
    }

    /// Represented as bytes 0..=255
    fn rgb(&self) -> [u8; 3] {
        match self {
            Self::White => [255, 255, 255],
            Self::Blue => [0, 144, 188],
            Self::Green => [0, 166, 81],
            Self::Yellow => [255, 212, 0],
            Self::Orange => [244, 121, 32],
            Self::Pink => [246, 158, 161],
            Self::Purple => [165, 67, 153],
            Self::LightBlue => [143, 204, 233],
            Self::Brown => [161, 122, 77],
            Self::Red => [237, 28, 36],
        }
    }
}

enum Shape {
    Circle,
    Square,
    CrossedSquare,
    DividedDiamond,
    Triangle,
    DividedTriangle,
    Rectangle,
    Star,
    Diamond,
    FlatDiamond,
    Hexagon,
    Pentagon,
}

impl GlycanStructure {
    fn render(&self) -> RenderedMonosaccharide {
        todo!()
    }
}

struct RenderedMonosaccharide {
    /// The depth along the main axis of the glycan, starting at 0 at the top (in the leaves)
    y: usize,
    /// The sideways placement starting at 0 at the leftmost uppermost monosaccharide, 1.0 is the width of one monosaccharide
    x: f32,
    /// The total width of the (sub)tree with all of its branches and sides
    width: f32,
    shape: Shape,
    colour: Colour,
    modifications: String,
    branches: Vec<RenderedMonosaccharide>,
    sides: Vec<RenderedMonosaccharide>,
}

impl RenderedMonosaccharide {
    fn transpose(&mut self, y: usize, x: f32) {
        self.y += y;
        self.x += x;
        for branch in &mut self.branches {
            branch.transpose(y, x);
        }
        for side in &mut self.sides {
            side.transpose(y, x);
        }
    }
}
