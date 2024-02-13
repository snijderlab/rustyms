#![allow(non_snake_case,non_upper_case_globals)]
use std::sync::OnceLock;
use crate::shared::{Germlines, Species};
/// Get the germlines for any of the available species. See the main documentation for which species have which data available.
pub fn germlines(species: Species) -> Option<&'static Germlines> {match species {
Species::BosTaurus => Some(lock_BosTaurus()),
Species::CamelusDromedarius => Some(lock_CamelusDromedarius()),
Species::CanisLupusFamiliaris => Some(lock_CanisLupusFamiliaris()),
Species::CapraHircus => Some(lock_CapraHircus()),
Species::DanioRerio => Some(lock_DanioRerio()),
Species::EquusCaballus => Some(lock_EquusCaballus()),
Species::FelisCatus => Some(lock_FelisCatus()),
Species::GallusGallus => Some(lock_GallusGallus()),
Species::GorillaGorilla => Some(lock_GorillaGorilla()),
Species::GorillaGorillaGorilla => Some(lock_GorillaGorillaGorilla()),
Species::HomoSapiens => Some(lock_HomoSapiens()),
Species::IctalurusPunctatus => Some(lock_IctalurusPunctatus()),
Species::LemurCatta => Some(lock_LemurCatta()),
Species::MacacaFascicularis => Some(lock_MacacaFascicularis()),
Species::MacacaMulatta => Some(lock_MacacaMulatta()),
Species::MusCookii => Some(lock_MusCookii()),
Species::MusMinutoides => Some(lock_MusMinutoides()),
Species::MusMusculus => Some(lock_MusMusculus()),
Species::MusMusculusDomesticus => Some(lock_MusMusculusDomesticus()),
Species::MusPahari => Some(lock_MusPahari()),
Species::MusSaxicola => Some(lock_MusSaxicola()),
Species::MusSpretus => Some(lock_MusSpretus()),
Species::MustelaPutoriusFuro => Some(lock_MustelaPutoriusFuro()),
Species::OncorhynchusMykiss => Some(lock_OncorhynchusMykiss()),
Species::OrnithorhynchusAnatinus => Some(lock_OrnithorhynchusAnatinus()),
Species::OryctolagusCuniculus => Some(lock_OryctolagusCuniculus()),
Species::OvisAries => Some(lock_OvisAries()),
Species::PongoAbelii => Some(lock_PongoAbelii()),
Species::PongoPygmaeus => Some(lock_PongoPygmaeus()),
Species::RattusNorvegicus => Some(lock_RattusNorvegicus()),
Species::SalmoSalar => Some(lock_SalmoSalar()),
Species::SusScrofa => Some(lock_SusScrofa()),
Species::VicugnaPacos => Some(lock_VicugnaPacos()),
_=>None}}
/// Get all germlines in one iterator, see the main documentation for more information about the available germlines
pub fn all_germlines() -> impl std::iter::Iterator<Item = &'static Germlines> {
std::iter::once(lock_BosTaurus())
.chain(std::iter::once(lock_CamelusDromedarius()))
.chain(std::iter::once(lock_CanisLupusFamiliaris()))
.chain(std::iter::once(lock_CapraHircus()))
.chain(std::iter::once(lock_DanioRerio()))
.chain(std::iter::once(lock_EquusCaballus()))
.chain(std::iter::once(lock_FelisCatus()))
.chain(std::iter::once(lock_GallusGallus()))
.chain(std::iter::once(lock_GorillaGorilla()))
.chain(std::iter::once(lock_GorillaGorillaGorilla()))
.chain(std::iter::once(lock_HomoSapiens()))
.chain(std::iter::once(lock_IctalurusPunctatus()))
.chain(std::iter::once(lock_LemurCatta()))
.chain(std::iter::once(lock_MacacaFascicularis()))
.chain(std::iter::once(lock_MacacaMulatta()))
.chain(std::iter::once(lock_MusCookii()))
.chain(std::iter::once(lock_MusMinutoides()))
.chain(std::iter::once(lock_MusMusculus()))
.chain(std::iter::once(lock_MusMusculusDomesticus()))
.chain(std::iter::once(lock_MusPahari()))
.chain(std::iter::once(lock_MusSaxicola()))
.chain(std::iter::once(lock_MusSpretus()))
.chain(std::iter::once(lock_MustelaPutoriusFuro()))
.chain(std::iter::once(lock_OncorhynchusMykiss()))
.chain(std::iter::once(lock_OrnithorhynchusAnatinus()))
.chain(std::iter::once(lock_OryctolagusCuniculus()))
.chain(std::iter::once(lock_OvisAries()))
.chain(std::iter::once(lock_PongoAbelii()))
.chain(std::iter::once(lock_PongoPygmaeus()))
.chain(std::iter::once(lock_RattusNorvegicus()))
.chain(std::iter::once(lock_SalmoSalar()))
.chain(std::iter::once(lock_SusScrofa()))
.chain(std::iter::once(lock_VicugnaPacos()))
}
/// Get all germlines in one parallel iterator, see the main documentation for more information about the available germlines
#[cfg(feature = "rayon")]
use rayon::prelude::*;
#[cfg(feature = "rayon")]
pub fn par_germlines() -> impl rayon::prelude::ParallelIterator<Item = &'static Germlines> {
rayon::iter::once(lock_BosTaurus())
.chain(rayon::iter::once(lock_CamelusDromedarius()))
.chain(rayon::iter::once(lock_CanisLupusFamiliaris()))
.chain(rayon::iter::once(lock_CapraHircus()))
.chain(rayon::iter::once(lock_DanioRerio()))
.chain(rayon::iter::once(lock_EquusCaballus()))
.chain(rayon::iter::once(lock_FelisCatus()))
.chain(rayon::iter::once(lock_GallusGallus()))
.chain(rayon::iter::once(lock_GorillaGorilla()))
.chain(rayon::iter::once(lock_GorillaGorillaGorilla()))
.chain(rayon::iter::once(lock_HomoSapiens()))
.chain(rayon::iter::once(lock_IctalurusPunctatus()))
.chain(rayon::iter::once(lock_LemurCatta()))
.chain(rayon::iter::once(lock_MacacaFascicularis()))
.chain(rayon::iter::once(lock_MacacaMulatta()))
.chain(rayon::iter::once(lock_MusCookii()))
.chain(rayon::iter::once(lock_MusMinutoides()))
.chain(rayon::iter::once(lock_MusMusculus()))
.chain(rayon::iter::once(lock_MusMusculusDomesticus()))
.chain(rayon::iter::once(lock_MusPahari()))
.chain(rayon::iter::once(lock_MusSaxicola()))
.chain(rayon::iter::once(lock_MusSpretus()))
.chain(rayon::iter::once(lock_MustelaPutoriusFuro()))
.chain(rayon::iter::once(lock_OncorhynchusMykiss()))
.chain(rayon::iter::once(lock_OrnithorhynchusAnatinus()))
.chain(rayon::iter::once(lock_OryctolagusCuniculus()))
.chain(rayon::iter::once(lock_OvisAries()))
.chain(rayon::iter::once(lock_PongoAbelii()))
.chain(rayon::iter::once(lock_PongoPygmaeus()))
.chain(rayon::iter::once(lock_RattusNorvegicus()))
.chain(rayon::iter::once(lock_SalmoSalar()))
.chain(rayon::iter::once(lock_SusScrofa()))
.chain(rayon::iter::once(lock_VicugnaPacos()))
}
static LOCK_BosTaurus: OnceLock<Germlines> = OnceLock::new();
fn lock_BosTaurus()->&'static Germlines{LOCK_BosTaurus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic bovine.bin")).unwrap()})}
static LOCK_CamelusDromedarius: OnceLock<Germlines> = OnceLock::new();
fn lock_CamelusDromedarius()->&'static Germlines{LOCK_CamelusDromedarius.get_or_init(|| {bincode::deserialize(include_bytes!("Arabian camel.bin")).unwrap()})}
static LOCK_CanisLupusFamiliaris: OnceLock<Germlines> = OnceLock::new();
fn lock_CanisLupusFamiliaris()->&'static Germlines{LOCK_CanisLupusFamiliaris.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic dog.bin")).unwrap()})}
static LOCK_CapraHircus: OnceLock<Germlines> = OnceLock::new();
fn lock_CapraHircus()->&'static Germlines{LOCK_CapraHircus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic goat.bin")).unwrap()})}
static LOCK_DanioRerio: OnceLock<Germlines> = OnceLock::new();
fn lock_DanioRerio()->&'static Germlines{LOCK_DanioRerio.get_or_init(|| {bincode::deserialize(include_bytes!("Zebrafish.bin")).unwrap()})}
static LOCK_EquusCaballus: OnceLock<Germlines> = OnceLock::new();
fn lock_EquusCaballus()->&'static Germlines{LOCK_EquusCaballus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic horse.bin")).unwrap()})}
static LOCK_FelisCatus: OnceLock<Germlines> = OnceLock::new();
fn lock_FelisCatus()->&'static Germlines{LOCK_FelisCatus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic cat.bin")).unwrap()})}
static LOCK_GallusGallus: OnceLock<Germlines> = OnceLock::new();
fn lock_GallusGallus()->&'static Germlines{LOCK_GallusGallus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic chicken.bin")).unwrap()})}
static LOCK_GorillaGorilla: OnceLock<Germlines> = OnceLock::new();
fn lock_GorillaGorilla()->&'static Germlines{LOCK_GorillaGorilla.get_or_init(|| {bincode::deserialize(include_bytes!("Western gorilla.bin")).unwrap()})}
static LOCK_GorillaGorillaGorilla: OnceLock<Germlines> = OnceLock::new();
fn lock_GorillaGorillaGorilla()->&'static Germlines{LOCK_GorillaGorillaGorilla.get_or_init(|| {bincode::deserialize(include_bytes!("Western lowland gorilla.bin")).unwrap()})}
static LOCK_HomoSapiens: OnceLock<Germlines> = OnceLock::new();
fn lock_HomoSapiens()->&'static Germlines{LOCK_HomoSapiens.get_or_init(|| {bincode::deserialize(include_bytes!("Human.bin")).unwrap()})}
static LOCK_IctalurusPunctatus: OnceLock<Germlines> = OnceLock::new();
fn lock_IctalurusPunctatus()->&'static Germlines{LOCK_IctalurusPunctatus.get_or_init(|| {bincode::deserialize(include_bytes!("Channel catfish.bin")).unwrap()})}
static LOCK_LemurCatta: OnceLock<Germlines> = OnceLock::new();
fn lock_LemurCatta()->&'static Germlines{LOCK_LemurCatta.get_or_init(|| {bincode::deserialize(include_bytes!("Ring-tailed lemur.bin")).unwrap()})}
static LOCK_MacacaFascicularis: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaFascicularis()->&'static Germlines{LOCK_MacacaFascicularis.get_or_init(|| {bincode::deserialize(include_bytes!("Crab-eating macaque.bin")).unwrap()})}
static LOCK_MacacaMulatta: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaMulatta()->&'static Germlines{LOCK_MacacaMulatta.get_or_init(|| {bincode::deserialize(include_bytes!("Rhesus monkey.bin")).unwrap()})}
static LOCK_MusCookii: OnceLock<Germlines> = OnceLock::new();
fn lock_MusCookii()->&'static Germlines{LOCK_MusCookii.get_or_init(|| {bincode::deserialize(include_bytes!("Cook's mouse.bin")).unwrap()})}
static LOCK_MusMinutoides: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMinutoides()->&'static Germlines{LOCK_MusMinutoides.get_or_init(|| {bincode::deserialize(include_bytes!("Southern African pygmy mouse.bin")).unwrap()})}
static LOCK_MusMusculus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMusculus()->&'static Germlines{LOCK_MusMusculus.get_or_init(|| {bincode::deserialize(include_bytes!("House mouse.bin")).unwrap()})}
static LOCK_MusMusculusDomesticus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMusculusDomesticus()->&'static Germlines{LOCK_MusMusculusDomesticus.get_or_init(|| {bincode::deserialize(include_bytes!("Western European house mouse.bin")).unwrap()})}
static LOCK_MusPahari: OnceLock<Germlines> = OnceLock::new();
fn lock_MusPahari()->&'static Germlines{LOCK_MusPahari.get_or_init(|| {bincode::deserialize(include_bytes!("Shrew mouse.bin")).unwrap()})}
static LOCK_MusSaxicola: OnceLock<Germlines> = OnceLock::new();
fn lock_MusSaxicola()->&'static Germlines{LOCK_MusSaxicola.get_or_init(|| {bincode::deserialize(include_bytes!("Spiny mouse.bin")).unwrap()})}
static LOCK_MusSpretus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusSpretus()->&'static Germlines{LOCK_MusSpretus.get_or_init(|| {bincode::deserialize(include_bytes!("Western wild mouse.bin")).unwrap()})}
static LOCK_MustelaPutoriusFuro: OnceLock<Germlines> = OnceLock::new();
fn lock_MustelaPutoriusFuro()->&'static Germlines{LOCK_MustelaPutoriusFuro.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic ferret.bin")).unwrap()})}
static LOCK_OncorhynchusMykiss: OnceLock<Germlines> = OnceLock::new();
fn lock_OncorhynchusMykiss()->&'static Germlines{LOCK_OncorhynchusMykiss.get_or_init(|| {bincode::deserialize(include_bytes!("Rainbow trout.bin")).unwrap()})}
static LOCK_OrnithorhynchusAnatinus: OnceLock<Germlines> = OnceLock::new();
fn lock_OrnithorhynchusAnatinus()->&'static Germlines{LOCK_OrnithorhynchusAnatinus.get_or_init(|| {bincode::deserialize(include_bytes!("Platypus.bin")).unwrap()})}
static LOCK_OryctolagusCuniculus: OnceLock<Germlines> = OnceLock::new();
fn lock_OryctolagusCuniculus()->&'static Germlines{LOCK_OryctolagusCuniculus.get_or_init(|| {bincode::deserialize(include_bytes!("Rabbit.bin")).unwrap()})}
static LOCK_OvisAries: OnceLock<Germlines> = OnceLock::new();
fn lock_OvisAries()->&'static Germlines{LOCK_OvisAries.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic sheep.bin")).unwrap()})}
static LOCK_PongoAbelii: OnceLock<Germlines> = OnceLock::new();
fn lock_PongoAbelii()->&'static Germlines{LOCK_PongoAbelii.get_or_init(|| {bincode::deserialize(include_bytes!("Sumatran orangutan.bin")).unwrap()})}
static LOCK_PongoPygmaeus: OnceLock<Germlines> = OnceLock::new();
fn lock_PongoPygmaeus()->&'static Germlines{LOCK_PongoPygmaeus.get_or_init(|| {bincode::deserialize(include_bytes!("Bornean orangutan.bin")).unwrap()})}
static LOCK_RattusNorvegicus: OnceLock<Germlines> = OnceLock::new();
fn lock_RattusNorvegicus()->&'static Germlines{LOCK_RattusNorvegicus.get_or_init(|| {bincode::deserialize(include_bytes!("Norway rat.bin")).unwrap()})}
static LOCK_SalmoSalar: OnceLock<Germlines> = OnceLock::new();
fn lock_SalmoSalar()->&'static Germlines{LOCK_SalmoSalar.get_or_init(|| {bincode::deserialize(include_bytes!("Atlantic salmon.bin")).unwrap()})}
static LOCK_SusScrofa: OnceLock<Germlines> = OnceLock::new();
fn lock_SusScrofa()->&'static Germlines{LOCK_SusScrofa.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic pig.bin")).unwrap()})}
static LOCK_VicugnaPacos: OnceLock<Germlines> = OnceLock::new();
fn lock_VicugnaPacos()->&'static Germlines{LOCK_VicugnaPacos.get_or_init(|| {bincode::deserialize(include_bytes!("Alpaca.bin")).unwrap()})}
