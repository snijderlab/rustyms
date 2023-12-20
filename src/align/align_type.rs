use std::str::FromStr;

use serde::{Deserialize, Serialize};

/// The type of alignment to perform
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Type {
    state: u8,
}

impl Type {
    /// Global alignment, which tries to find the best alignment to link both sequences fully to each other, like the Needleman Wunsch algorithm
    pub const GLOBAL: Self = Self::new(true, true, true, true);
    /// Local alignment, which tries to find the best patch of both sequences to align to each other, this could lead to trailing ends on both sides of both sequences, like the Smith Waterman
    pub const LOCAL: Self = Self::new(false, false, false, false);
    /// Hybrid alignment, the first sequence will be fully aligned to the second sequence, this could lead to trailing ends on the second sequence but not on the first.
    pub const GLOBAL_A: Self = Self::new(true, false, true, false);
    /// Hybrid alignment, the second sequence will be fully aligned to the first sequence, this could lead to trailing ends on the first sequence but not on the second.
    pub const GLOBAL_B: Self = Self::new(false, true, false, true);
    /// Hybrid alignment, globally align the left (start) side of the alignment
    pub const GLOBAL_LEFT: Self = Self::new(true, true, false, false);
    /// Hybrid alignment, globally align the right (end) side of the alignment
    pub const GLOBAL_RIGHT: Self = Self::new(false, false, true, true);
    /// Hybrid alignment, extend sequence a with sequence b (▀▀██▄▄)
    pub const EXTEND_A: Self = Self::new(false, true, true, false);
    /// Hybrid alignment, extend sequence b with sequence a (▄▄██▀▀)
    pub const EXTEND_B: Self = Self::new(true, false, false, true);

    const LEFT_A: u8 = 0;
    const LEFT_B: u8 = 1;
    const RIGHT_A: u8 = 2;
    const RIGHT_B: u8 = 3;

    /// Create a new alignment type
    #[allow(clippy::fn_params_excessive_bools)]
    pub const fn new(left_a: bool, left_b: bool, right_a: bool, right_b: bool) -> Self {
        const fn fu8(v: bool) -> u8 {
            if v {
                1
            } else {
                0
            }
        }
        Self {
            state: (fu8(left_a) << Self::LEFT_A)
                | (fu8(left_b) << Self::LEFT_B)
                | (fu8(right_a) << Self::RIGHT_A)
                | (fu8(right_b) << Self::RIGHT_B),
        }
    }

    /// Check if left a is global
    pub const fn left_a(&self) -> bool {
        self.state >> Self::LEFT_A & 1 == 1
    }
    /// Check if left b is global
    pub const fn left_b(&self) -> bool {
        self.state >> Self::LEFT_B & 1 == 1
    }
    /// Check if right a is global
    pub const fn right_a(&self) -> bool {
        self.state >> Self::RIGHT_A & 1 == 1
    }
    /// Check if right b is global
    pub const fn right_b(&self) -> bool {
        self.state >> Self::RIGHT_B & 1 == 1
    }

    /// Get a descriptive name for this alignment type
    pub const fn description(&self) -> &'static str {
        match (self.left_a(), self.left_b(), self.right_a(), self.right_b()) {
            (true, true, true, true) => "global",
            (false, false, false, false) => "local",
            (true, false, true, false) => "global A",
            (false, true, false, true) => "global B",
            (true, true, false, false) => "global left",
            (false, false, true, true) => "global right",
            (false, true, true, false) => "extend A",
            (true, false, false, true) => "extend B",
            _ => "special",
        }
    }

    /// Get a concise symbol for this alignment type (four quadrants, each filled if that spot is global)
    pub const fn symbol(&self) -> char {
        match (self.left_a(), self.left_b(), self.right_a(), self.right_b()) {
            (false, false, false, false) => ' ',
            (true, false, false, false) => '▘',
            (false, true, false, false) => '▖',
            (false, false, true, false) => '▝',
            (false, false, false, true) => '▗',
            (true, false, true, false) => '▀',
            (false, true, false, true) => '▄',
            (true, true, false, false) => '▌',
            (false, false, true, true) => '▐',
            (false, true, true, false) => '▞',
            (true, false, false, true) => '▚',
            (false, true, true, true) => '▟',
            (true, false, true, true) => '▜',
            (true, true, false, true) => '▙',
            (true, true, true, false) => '▛',
            (true, true, true, true) => '█',
        }
    }
}

impl std::fmt::Display for Type {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {}", self.symbol(), self.description())
    }
}

impl FromStr for Type {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "0" | "0000" | " " | "l" | "local" => Ok(Self::new(false, false, false, false)),
            "8" | "1000" | "▘" => Ok(Self::new(true, false, false, false)),
            "4" | "0100" | "▖" => Ok(Self::new(false, true, false, false)),
            "2" | "0010" | "▝" => Ok(Self::new(false, false, true, false)),
            "1" | "0001" | "▗" => Ok(Self::new(false, false, false, true)),
            "10" | "1010" | "▀" | "a" | "global a" | "global_a" => {
                Ok(Self::new(true, false, true, false))
            }
            "5" | "0101" | "▄" | "b" | "global b" | "global_b" => {
                Ok(Self::new(false, true, false, true))
            }
            "12" | "1100" | "▌" | "left" | "global left" | "global_left" => {
                Ok(Self::new(true, true, false, false))
            }
            "3" | "0011" | "▐" | "right" | "global right" | "global_right" => {
                Ok(Self::new(false, false, true, true))
            }
            "6" | "0110" | "▞" | "ea" | "extend a" | "extend_a" => {
                Ok(Self::new(false, true, true, false))
            }
            "9" | "1001" | "▚" | "eb" | "extend b" | "extend_b" => {
                Ok(Self::new(true, false, false, true))
            }
            "7" | "0111" | "▟" => Ok(Self::new(false, true, true, true)),
            "11" | "1011" | "▜" => Ok(Self::new(true, false, true, true)),
            "13" | "1101" | "▙" => Ok(Self::new(true, true, false, true)),
            "14" | "1110" | "▛" => Ok(Self::new(true, true, true, false)),
            "15" | "1111" | "█" | "g" | "global" => Ok(Self::new(true, true, true, true)),
            _ => Err(()),
        }
    }
}

impl Default for Type {
    /// Defaults to global alignment
    fn default() -> Self {
        Self::GLOBAL
    }
}

impl Type {
    pub(super) const fn global(self) -> bool {
        !matches!(self, Self::LOCAL)
    }
}
