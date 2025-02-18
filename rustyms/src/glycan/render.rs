use std::{f32::consts::PI, fmt::Write};

use super::{
    BaseSugar, GlycanStructure, GlycanSubstituent, HeptoseIsomer, HexoseIsomer, MonoSaccharide,
    PentoseIsomer,
};

impl MonoSaccharide {
    fn get_shape(&self) -> (Shape, Colour, String, String) {
        // Common substitutions
        let mut nacetyl = 0;
        let mut acid = 0;
        let mut amino = 0;
        let mut deoxy = 0;
        // Additional needed substitutions
        let mut acetyl = 0;
        let mut glycolyl = 0;
        let mut nglycolyl = 0;
        let mut o_carboxy_ethyl = 0;
        let mut inner_modifications = if self.furanose {
            "f".to_string()
        } else {
            String::new()
        };
        let mut outer_modifications = String::new();
        for m in &self.substituents {
            match m {
                GlycanSubstituent::NAcetyl => nacetyl += 1,
                GlycanSubstituent::Acid => acid += 1,
                GlycanSubstituent::Amino => amino += 1,
                GlycanSubstituent::Deoxy => deoxy += 1,
                GlycanSubstituent::Acetyl => acetyl += 1,
                GlycanSubstituent::Glycolyl => glycolyl += 1,
                GlycanSubstituent::OCarboxyEthyl => o_carboxy_ethyl += 1,
                GlycanSubstituent::NGlycolyl => nglycolyl += 1,
                GlycanSubstituent::Didehydro => inner_modifications.push_str("en"), // Missing symbols: o for alditols, an for anhydro, on for lactone, am for lactam
                _ => outer_modifications.push_str(m.notation()), // This does overlap with mods that are sometimes expected for things
            }
        }
        let mut outer_mods = |acetyl: usize, glycolyl: usize, nglycolyl: usize, o_carboxy_ethyl: usize| {
            [GlycanSubstituent::Acetyl.notation().repeat(acetyl),GlycanSubstituent::Glycolyl.notation().repeat(glycolyl),GlycanSubstituent::NGlycolyl.notation().repeat(nglycolyl),GlycanSubstituent::OCarboxyEthyl.notation().repeat(o_carboxy_ethyl),outer_modifications].join("")
        };
        match (&self.base_sugar, nacetyl, acid, amino, deoxy) {
            (BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), 1, 0, 0, 0)
                if o_carboxy_ethyl > 0 =>
            {
                (
                    Shape::Hexagon,
                    Colour::Purple,
                    inner_modifications,
                    outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl - 1),
                )
            }
            (BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), 0, 0, 0, 0)
                if nglycolyl > 0 && o_carboxy_ethyl > 0 =>
            {
                (
                    Shape::Hexagon,
                    Colour::LightBlue,
                    inner_modifications,
                    outer_mods(acetyl, glycolyl, nglycolyl - 1, o_carboxy_ethyl - 1),
                )
            }
            (BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), 0, 0, 1, 0)
                if o_carboxy_ethyl > 0 =>
            {
                (
                    Shape::Hexagon,
                    Colour::Brown,
                    inner_modifications,
                    outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl - 1),
                )
            }
            (BaseSugar::Hexose(isomer), 0, 0, 0, 0) => {
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
                (s, c, inner_modifications, outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl))
            }
            (BaseSugar::Hexose(isomer), 1, 0, 0, 0) => {
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
                (Shape::Square, c, inner_modifications, outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl))
            }
            (BaseSugar::Hexose(isomer), 0, 0, 1, 0) => {
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
                (
                    Shape::CrossedSquare,
                    c,
                    inner_modifications,
                    outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
                )
            }
            (BaseSugar::Hexose(isomer), 0, 1, 0, 0) => {
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
                (
                    Shape::DividedDiamond,
                    c,
                    inner_modifications,
                    outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
                )
            }
            (BaseSugar::Hexose(isomer), 0, 0, 0, 1) => {
                let c = match isomer {
                    Some(HexoseIsomer::Glucose) => Colour::Blue,
                    Some(HexoseIsomer::Mannose) => Colour::Green,
                    Some(HexoseIsomer::Galactose) => Colour::Red,
                    Some(HexoseIsomer::Gulose) => Colour::Orange,
                    Some(HexoseIsomer::Altrose) => Colour::Pink,
                    Some(HexoseIsomer::Talose) => Colour::LightBlue,
                    Some(_) | None => Colour::White,
                };
                (Shape::Triangle, c, inner_modifications, outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl))
            }
            (BaseSugar::Hexose(isomer), 1, 0, 0, 1) => {
                let c = match isomer {
                    Some(HexoseIsomer::Glucose) => Colour::Blue,
                    Some(HexoseIsomer::Mannose) => Colour::Green,
                    Some(HexoseIsomer::Galactose) => Colour::Red,
                    Some(HexoseIsomer::Altrose) => Colour::Pink,
                    Some(HexoseIsomer::Talose) => Colour::LightBlue,
                    Some(_) | None => Colour::White,
                };
                (
                    Shape::DividedTriangle,
                    c,
                    inner_modifications,
                    outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
                )
            }
            (BaseSugar::Hexose(isomer), 0, 0, 0, 2) => {
                let c = match isomer {
                    Some(HexoseIsomer::Glucose) => Colour::Blue,
                    Some(HexoseIsomer::Mannose) => Colour::Green,
                    Some(HexoseIsomer::Galactose) => Colour::Orange,
                    Some(HexoseIsomer::Altrose) => Colour::Pink,
                    Some(HexoseIsomer::Allose) => Colour::Purple,
                    Some(HexoseIsomer::Talose) => Colour::LightBlue,
                    Some(_) | None => Colour::White,
                };
                (
                    Shape::Rectangle,
                    c,
                    inner_modifications,
                    outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
                )
            }
            (BaseSugar::Pentose(isomer), 0, 0, 0, 0) => (
                Shape::Star,
                match isomer {
                    None | Some(PentoseIsomer::Xylulose) => Colour::White,
                    Some(PentoseIsomer::Arabinose) => Colour::Green,
                    Some(PentoseIsomer::Lyxose) => Colour::Yellow,
                    Some(PentoseIsomer::Xylose) => Colour::Orange,
                    Some(PentoseIsomer::Ribose) => Colour::Pink,
                },
                inner_modifications,
                outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
            ),
            (BaseSugar::Nonose, 0, 1, 2, 2) => (
                Shape::FlatDiamond,
                Colour::White,
                inner_modifications,
                outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
            ), // This does not detect which exact isomer it is, as that information is currently not tracked
            (BaseSugar::Nonose, 0, a, b, _) if a > 0 && b > 0 => (
                Shape::Diamond,
                match (acid, amino, deoxy) {
                    (1, 1, 1) => Colour::Red, // This could either be Sia (Red) or Kdo (Green) for now defaults to red until the isomeric state is tracked
                    (1, 1, 0) if acetyl > 0 => {
                        Colour::Purple
                    }
                    (1, 1, 0) if glycolyl > 0 => {
                        Colour::LightBlue
                    }
                    (1, 1, 0) => Colour::Brown,
                    _ => Colour::White,
                },
                inner_modifications,
                outer_mods(acetyl.saturating_sub(1), glycolyl.saturating_sub(1), nglycolyl, o_carboxy_ethyl),
            ),
            (BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), 0, 0, 2, 1) => (
                Shape::Hexagon,
                Colour::Blue,
                inner_modifications,
                outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
            ),
            (BaseSugar::Heptose(Some(HeptoseIsomer::GlyceroMannoHeptopyranose)), 0, 0, 0, 0) => (
                Shape::Hexagon,
                Colour::Green,
                inner_modifications,
                outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
            ),
            (BaseSugar::Octose, 0, 1, 0, 1) => (
                Shape::Hexagon,
                Colour::Yellow,
                inner_modifications,
                outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
            ),
            (BaseSugar::Heptose(None), 0, 2, 0, 1) => (
                Shape::Hexagon,
                Colour::Orange,
                inner_modifications,
                outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
            ),
            _ => (
                Shape::Hexagon,
                Colour::White,
                inner_modifications,
                outer_mods(acetyl, glycolyl, nglycolyl, o_carboxy_ethyl),
            ),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
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

#[derive(Debug, PartialEq, Eq)]
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

impl Shape {
    fn height(&self) -> f32 {
        match self {
            Self::Rectangle | Self::FlatDiamond | Self::Hexagon => 0.5,
            _ => 1.0,
        }
    }
}

impl GlycanStructure {
    fn render(&self) -> RenderedMonosaccharide {
        let (shape, colour, inner_modifications, outer_modifications) = self.sugar.get_shape();
        if self.branches.is_empty() {
            RenderedMonosaccharide {
                y: 0,
                x: 0.0,
                mid_point: 0.5,
                width: 1.0,
                shape,
                colour,
                inner_modifications,
                outer_modifications,
                branches: Vec::new(),
                sides: Vec::new(),
            }
        } else {
            let mut depth = 0;
            let mut branches = Vec::new();
            let mut sides = Vec::new();
            for branch in &self.branches {
                let rendered = branch.render();
                if rendered.is_sideways() && sides.len() < 2 {
                    sides.push(rendered);
                } else {
                    depth = depth.max(rendered.y);
                    branches.push(rendered);
                }
            }
            // Update all branch placements
            let mut displacement = 0.0;
            for branch in &mut branches {
                branch.transpose(depth - branch.y, displacement);
                displacement += branch.width;
            }
            depth += 1;
            // Determine the center point for this sugar
            let mut center = match branches.len() {
                0 => 0.5,
                1 => branches[0].mid_point,
                _ => {
                    let left_midpoint = branches[0].mid_point;
                    let right_midpoint =
                        branches.last().unwrap().x + branches.last().unwrap().mid_point;
                    (left_midpoint + right_midpoint) / 2.0
                }
            };
            let mut width = branches.last().map_or(1.0, |b| b.x + b.width);
            if !sides.is_empty() {
                sides[0].transpose(depth, center + 0.5);
                width = width.max(center + 0.5 + sides[0].width);
            }
            if sides.len() == 2 {
                let mut x = center - 0.5 - sides[1].width;
                if x < 0.0 {
                    let shift = -x;
                    center += shift;
                    for branch in &mut branches {
                        branch.transpose(0, shift);
                    }
                    sides[0].transpose(0, shift);
                    width += shift;
                    x = 0.0;
                }
                sides[1].transpose(depth, x);
            }
            RenderedMonosaccharide {
                y: depth,
                x: 0.0,
                mid_point: center,
                width,
                shape,
                colour,
                inner_modifications,
                outer_modifications,
                branches,
                sides,
            }
        }
    }
}

struct RenderedMonosaccharide {
    /// The depth of this sugar along the main axis of the glycan, starting at 0 at the top (in the leaves)
    y: usize,
    /// The sideways placement of this whole tree starting at 0 at the leftmost monosaccharide, 1.0 is the width of one monosaccharide
    x: f32,
    /// The sideways placement of this sugar within this tree, for the absolute sideways placement of this sugar add this to `x`
    mid_point: f32,
    /// The total width of the (sub)tree with all of its branches and sides
    width: f32,
    shape: Shape,
    colour: Colour,
    inner_modifications: String,
    outer_modifications: String,
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

    /// Check if this sugar should be rendered to the side of the parent sugar
    fn is_sideways(&self) -> bool {
        self.colour == Colour::Red
            && self.shape == Shape::Triangle
            && self.branches.is_empty()
            && self.sides.is_empty()
    }

    fn to_svg(&self) -> String {
        const COLUMN_SIZE: f32 = 10.0;
        const SUGAR_SIZE: f32 = 5.0;
        const STROKE_SIZE: f32 = 0.5;

        fn render_elements(buffer: &mut String, element: &RenderedMonosaccharide) {
            // First all lines to get good stacking behaviour
            for branch in &element.branches {
                write!(
                    buffer,
                    "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                    (element.x + element.mid_point) * COLUMN_SIZE,
                    (element.y as f32 + 0.5) * COLUMN_SIZE,
                    (branch.x + branch.mid_point) * COLUMN_SIZE,
                    // (branch.y as f32).mul_add(COLUMN_SIZE, branch.shape.height().mul_add(SUGAR_SIZE, -branch.shape.height().mul_add(SUGAR_SIZE, -COLUMN_SIZE) / 2.0))
                     (branch.y as f32 + 0.5) * COLUMN_SIZE
                )
                .unwrap();
            }
            for branch in &element.sides {
                write!(
                    buffer,
                    "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                    (element.x + element.mid_point) * COLUMN_SIZE,
                    (element.y as f32 + 0.5) * COLUMN_SIZE,
                    (branch.x + branch.mid_point) * COLUMN_SIZE,
                    (branch.y as f32 + 0.5) * COLUMN_SIZE
                )
                .unwrap();
            }
            // Render the sugar
            let colour = element.colour.rgb();
            let fill = format!("rgb({},{},{})", colour[0], colour[1], colour[2]);
            match element.shape {
                Shape::Circle => write!(
                    buffer,
                    "<circle r=\"{}\" cx=\"{}\" cy=\"{}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                    SUGAR_SIZE / 2.0,
                    (element.x + element.mid_point) * COLUMN_SIZE,
                    (element.y as f32 + 0.5) * COLUMN_SIZE
                )
                .unwrap(),
                Shape::Square => write!(
                    buffer,
                    "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                    (element.x + element.mid_point).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 2.0),
                    (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 2.0),
                    SUGAR_SIZE,
                    SUGAR_SIZE
                )
                .unwrap(),
                Shape::Rectangle => write!(
                    buffer,
                    "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                    (element.x + element.mid_point).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 2.0),
                    (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 4.0),
                    SUGAR_SIZE,
                    SUGAR_SIZE / 2.0
                )
                .unwrap(),
                Shape::CrossedSquare => {
                    let x1 = (element.x + element.mid_point).mul_add(COLUMN_SIZE,- SUGAR_SIZE / 2.0);
                    let y1 = (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 2.0);
                    let x2 = x1 + SUGAR_SIZE;
                    let y2 = y1 + SUGAR_SIZE;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y1} {x2} {y1} {x2} {y2}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\" stroke-linejoin=\"bevel\"/><polygon points=\"{x1} {y1} {x1} {y2} {x2} {y2}\" fill=\"white\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\" stroke-linejoin=\"bevel\"/>",
                )
                .unwrap();},
                Shape::DividedDiamond => {
                    let x1 = (element.x + element.mid_point).mul_add(COLUMN_SIZE,- SUGAR_SIZE / 2.0);
                    let x2 = x1 + SUGAR_SIZE / 2.0;
                    let x3 = x1 + SUGAR_SIZE;
                    let y1 = (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 2.0);
                    let y2 = y1 + SUGAR_SIZE / 2.0;
                    let y3 = y1 + SUGAR_SIZE;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y2} {x2} {y1} {x3} {y2}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\" stroke-linejoin=\"bevel\"/><polygon points=\"{x1} {y2} {x2} {y3} {x3} {y2}\" fill=\"white\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\" stroke-linejoin=\"bevel\"/>",
                )
                .unwrap();},
                Shape::Triangle => {
                    let x1 = (element.x + element.mid_point).mul_add(COLUMN_SIZE,- SUGAR_SIZE / 2.0);
                    let x2 = x1 + SUGAR_SIZE / 2.0;
                    let x3 = x1 + SUGAR_SIZE;
                    let y1 = (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 2.0);
                    let y2 = y1 + SUGAR_SIZE;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y2} {x2} {y1} {x3} {y2}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                )
                .unwrap();},
                Shape::DividedTriangle => {
                    let x1 = (element.x + element.mid_point).mul_add(COLUMN_SIZE,- SUGAR_SIZE / 2.0);
                    let x2 = x1 + SUGAR_SIZE / 2.0;
                    let x3 = x1 + SUGAR_SIZE;
                    let y1 = (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 2.0);
                    let y2 = y1 + SUGAR_SIZE;

                    write!(
                    buffer,
                    "<polygon points=\"{x2} {y1} {x3} {y2} {x2} {y1}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\" stroke-linejoin=\"bevel\"/><polygon points=\"{x2} {y1} {x1} {y2} {x2} {y1}\" fill=\"white\"  stroke=\"black\" stroke-width=\"{STROKE_SIZE}\" stroke-linejoin=\"bevel\"/>",
                )
                .unwrap();},
                Shape::Diamond => {
                    let x1 = (element.x + element.mid_point).mul_add(COLUMN_SIZE,- SUGAR_SIZE / 2.0);
                    let x2 = x1 + SUGAR_SIZE / 2.0;
                    let x3 = x1 + SUGAR_SIZE;
                    let y1 = (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 2.0);
                    let y2 = y1 + SUGAR_SIZE / 2.0;
                    let y3 = y1 + SUGAR_SIZE;

                    write!(
                    buffer,
                    "<polygon points=\"{x2} {y1} {x3} {y2} {x2} {y3} {x1} {y2}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                )
                .unwrap();},
                Shape::FlatDiamond => {
                    let x1 = (element.x + element.mid_point).mul_add(COLUMN_SIZE,- SUGAR_SIZE / 2.0);
                    let x2 = x1 + SUGAR_SIZE / 2.0;
                    let x3 = x1 + SUGAR_SIZE;
                    let y1 = (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 4.0);
                    let y2 = y1 + SUGAR_SIZE / 4.0;
                    let y3 = y1 + SUGAR_SIZE / 2.0;

                    write!(
                    buffer,
                    "<polygon points=\"{x2} {y1} {x3} {y2} {x2} {y3} {x1} {y2}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                )
                .unwrap();},
                Shape::Hexagon => {
                    let x1 = (element.x + element.mid_point).mul_add(COLUMN_SIZE,- SUGAR_SIZE / 2.0);
                    let x2 = x1 + SUGAR_SIZE / 3.0;
                    let x3 = x1 + SUGAR_SIZE / 3.0 * 2.0;
                    let x4 = x1 + SUGAR_SIZE;
                    let y1 = (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 4.0);
                    let y2 = y1 + SUGAR_SIZE / 4.0;
                    let y3 = y1 + SUGAR_SIZE / 2.0;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y2} {x2} {y1} {x3} {y1} {x4} {y2} {x3} {y3} {x2} {y2}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                )
                .unwrap();},
                Shape::Pentagon => {
                    let a = (18.0/360.0 * 2.0 * PI).cos() * SUGAR_SIZE / 2.0;
                    let b = (18.0/360.0 * 2.0 * PI).sin() * SUGAR_SIZE / 2.0;
                    let c = (36.0/360.0*2.0*PI).cos() * SUGAR_SIZE / 2.0;
                    let d = (36.0/360.0*2.0*PI).sin() * SUGAR_SIZE / 2.0;
                    let base_x = (element.x + element.mid_point) * COLUMN_SIZE;
                    let x1 = base_x - a;
                    let x2 = base_x - d;
                    let x3 = base_x;
                    let x4 = base_x + d;
                    let x5 = base_x + a;
                    let y1 = (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 2.0);
                    let y2 = y1 + SUGAR_SIZE / 2.0 - b;
                    let y3 = y1 + SUGAR_SIZE / 2.0 + c;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y2} {x3} {y1} {x5} {y2} {x4} {y3} {x2} {y3}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                )
                .unwrap();},
                Shape::Star => {
                    const PHI: f32 = 1.618033988749894848204586834365638118_f32;
                    let a = (18.0/360.0 * 2.0 * PI).cos() * SUGAR_SIZE / 2.0;
                    let b = (18.0/360.0 * 2.0 * PI).sin() * SUGAR_SIZE / 2.0;
                    let c = (36.0/360.0*2.0*PI).cos() * SUGAR_SIZE / 2.0;
                    let d = (36.0/360.0*2.0*PI).sin() * SUGAR_SIZE / 2.0;
                    let e = d / (36.0/360.0*2.0*PI).tan();
                    let i = 2.0 * a / (54.0/360.0*2.0*PI).sin() / (1.0 + 1.0 / PHI);
                    let f = (18.0/360.0 * 2.0 * PI).cos() * i - SUGAR_SIZE / 2.0; 
                    let g = (18.0/360.0 * 2.0 * PI).sin() * i; 
                    let h = (SUGAR_SIZE / 2.0 - b) * (18.0/360.0*2.0*PI).tan();
                    let j = (18.0/360.0*2.0*PI).tan() * g;
                    // dbg!((a,b,c,d,e,i,f,g,h,j));
                    let base_x = (element.x + element.mid_point) * COLUMN_SIZE;
                    let x1 = base_x - a;
                    let x2 = base_x - d;
                    let x3 = base_x - g;
                    let x4 = base_x - h;
                    let x5 = base_x;
                    let x6 = base_x + h;
                    let x7 = base_x + g;
                    let x8 = base_x + d;
                    let x9 = base_x + a;
                    let y1 = (element.y as f32 + 0.5).mul_add(COLUMN_SIZE, - SUGAR_SIZE / 2.0);
                    let y2 = y1 + SUGAR_SIZE / 2.0 - b;
                    let y3 = y1 + SUGAR_SIZE / 2.0 + j;
                    let y4 = y1 + SUGAR_SIZE / 2.0 + f;
                    let y5 = y1 + SUGAR_SIZE / 2.0 + c;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y2} {x4} {y2} {x5} {y1} {x6} {y2} {x9} {y2} {x7} {y3} {x8} {y5} {x5} {y4} {x2} {y5} {x3} {y3}\" fill=\"{fill}\" stroke=\"black\" stroke-width=\"{STROKE_SIZE}\"/>",
                )
                .unwrap();},
            }
            if !element.inner_modifications.is_empty() {
                write!(buffer, "<text x=\"{}\" y=\"{}\" fill=\"black\" text-anchor=\"middle\" font-style=\"italic\" font-size=\"{}px\" dominant-baseline=\"middle\">{}</text>",
                    (element.x + element.mid_point) * COLUMN_SIZE,
                    (element.y as f32 + 0.5) * COLUMN_SIZE, 
                    SUGAR_SIZE / 2.0,
                    element.inner_modifications).unwrap();
            }
            if !element.outer_modifications.is_empty() {
                write!(buffer, "<text x=\"{}\" y=\"{}\" fill=\"black\" text-anchor=\"middle\" font-size=\"{}px\" dominant-baseline=\"ideographic\">{}</text>",
                    (element.x + element.mid_point) * COLUMN_SIZE,
                    (element.y as f32).mul_add(COLUMN_SIZE, -element.shape.height().mul_add(SUGAR_SIZE, -COLUMN_SIZE) / 2.0), 
                    SUGAR_SIZE / 2.0,
                    element.outer_modifications).unwrap();
            }
            // Render all connected sugars
            for branch in &element.branches {
                render_elements(buffer, branch);
            }
            for branch in &element.sides {
                render_elements(buffer, branch);
            }
        }

        let mut picture = String::new();
        write!(
            &mut picture,
            "<svg width=\"{}\" height=\"{}\">",
            self.width * COLUMN_SIZE,
            (self.y as f32 + 1.0) * COLUMN_SIZE
        )
        .unwrap();

        render_elements(&mut picture, self);

        write!(&mut picture, "</svg>").unwrap();

        picture
    }
}

#[test]
fn test_rendering() {
    let structure = GlycanStructure::from_short_iupac("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", 0..105, 0).unwrap(); // G01670UQ
    let rendered = structure.render();
    let svg = rendered.to_svg();
    println!("{svg}");
    let structure =
        GlycanStructure::from_short_iupac("Fuc(?1-?)Gal(?1-?)GalNAc(?1-", 0..28, 0).unwrap(); // G13523IF
    let rendered = structure.render();
    let svg = rendered.to_svg();
    println!("{svg}");
    let structure = GlycanStructure::from_short_iupac(
        "GlcN(b1-4)GlcNAc(b1-4)GlcNAc(b1-4)GlcNAc6S(?1-",
        0..46,
        0,
    )
    .unwrap(); // G00613DO
    let rendered = structure.render();
    let svg = rendered.to_svg();
    println!("{svg}");
    let structure = GlycanStructure::from_short_iupac("Neu5Gc(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Gal(a1-3)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Neu5Ac(a2-8)Neu5Ac(a2-3/6)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-", 0..256, 0).unwrap(); // G00621IU
    let rendered = structure.render();
    let svg = rendered.to_svg();
    println!("{svg}");
    let structure = GlycanStructure::from_short_iupac("Rha2,3,4Ac3(a1-2)[Xyl(b1-3)]Ara(a1-", 0..35, 0).unwrap(); // G01464QV
    let rendered = structure.render();
    let svg = rendered.to_svg();
    println!("{svg}");
    let structure = GlycanStructure::from_short_iupac("Fruf(b2-1a)[Glc(a1-2)Glc(a1-2)Glc(a1-2)Glc(a1-2)Glc(a1-2)Glc(a1-2)]Glc", 0..70, 0).unwrap(); // G04421VO
    let rendered = structure.render();
    let svg = rendered.to_svg();
    println!("{svg}");
    let structure = GlycanStructure::from_short_iupac("Kdn(a2-3)Gal(b1-4)ManNAc(b1-2)[Kdn(a2-3)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[GlcNAc(b1-4)][Kdn(a2-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Gc(a2-3)Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-", 0..203, 0).unwrap(); // G04458LN
    let rendered = structure.render();
    let svg = rendered.to_svg();
    println!("{svg}");
    todo!()
}
