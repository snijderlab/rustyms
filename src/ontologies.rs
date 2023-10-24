use crate::aminoacids::*;
use crate::element::*;
use crate::glycan::*;
use crate::modification::*;
use crate::placement_rules::*;

include!(concat!(env!("OUT_DIR"), "/unimod.rs"));
include!(concat!(env!("OUT_DIR"), "/psi-mod.rs"));
//include!(concat!(env!("OUT_DIR"), "/gnome.rs"));
pub const GNOME_ONTOLOGY: &[(usize, &str, Modification)] = unsafe {
    std::mem::transmute::<&[u8], &[(usize, String, GNOmeModification)]>(include_bytes!(concat!(
        env!("OUT_DIR"),
        "/gnome.dat"
    )))
};
