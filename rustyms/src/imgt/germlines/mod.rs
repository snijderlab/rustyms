// @generated
#![allow(non_snake_case,non_upper_case_globals)]
use std::sync::OnceLock;
use super::shared::{Germlines, Species};
/// Get the germlines for any of the available species. See the main documentation for which species have which data available.
pub fn germlines(species: Species) -> Option<&'static Germlines> {match species {
Species::AnarhichasMinor => Some(lock_AnarhichasMinor()),
Species::BosTaurus => Some(lock_BosTaurus()),
Species::CamelusDromedarius => Some(lock_CamelusDromedarius()),
Species::CanisLupusFamiliaris => Some(lock_CanisLupusFamiliaris()),
Species::CapraHircus => Some(lock_CapraHircus()),
Species::CarcharhinusPlumbeus => Some(lock_CarcharhinusPlumbeus()),
Species::CercocebusAtys => Some(lock_CercocebusAtys()),
Species::ChaenocephalusAceratus => Some(lock_ChaenocephalusAceratus()),
Species::CyprinusCarpio => Some(lock_CyprinusCarpio()),
Species::DanioRerio => Some(lock_DanioRerio()),
Species::DicentrarchusLabrax => Some(lock_DicentrarchusLabrax()),
Species::EquusCaballus => Some(lock_EquusCaballus()),
Species::FelisCatus => Some(lock_FelisCatus()),
Species::GadusMorhua => Some(lock_GadusMorhua()),
Species::GallusGallus => Some(lock_GallusGallus()),
Species::GasterosteusAculeatus => Some(lock_GasterosteusAculeatus()),
Species::GinglymostomaCirratum => Some(lock_GinglymostomaCirratum()),
Species::GorillaGorilla => Some(lock_GorillaGorilla()),
Species::GorillaGorillaGorilla => Some(lock_GorillaGorillaGorilla()),
Species::HeterodontusFrancisci => Some(lock_HeterodontusFrancisci()),
Species::HomoSapiens => Some(lock_HomoSapiens()),
Species::HydrolagusColliei => Some(lock_HydrolagusColliei()),
Species::HylobatesLar => Some(lock_HylobatesLar()),
Species::IctalurusPunctatus => Some(lock_IctalurusPunctatus()),
Species::LemurCatta => Some(lock_LemurCatta()),
Species::LeucorajaErinacea => Some(lock_LeucorajaErinacea()),
Species::MacacaArctoides => Some(lock_MacacaArctoides()),
Species::MacacaCyclopis => Some(lock_MacacaCyclopis()),
Species::MacacaFascicularis => Some(lock_MacacaFascicularis()),
Species::MacacaMulatta => Some(lock_MacacaMulatta()),
Species::MacacaNemestrina => Some(lock_MacacaNemestrina()),
Species::MacacaSilenus => Some(lock_MacacaSilenus()),
Species::MacacaThibetana => Some(lock_MacacaThibetana()),
Species::MesocricetusAuratus => Some(lock_MesocricetusAuratus()),
Species::MonodelphisDomestica => Some(lock_MonodelphisDomestica()),
Species::MusCookii => Some(lock_MusCookii()),
Species::MusMinutoides => Some(lock_MusMinutoides()),
Species::MusMusculus => Some(lock_MusMusculus()),
Species::MusMusculusCastaneus => Some(lock_MusMusculusCastaneus()),
Species::MusMusculusDomesticus => Some(lock_MusMusculusDomesticus()),
Species::MusMusculusMolossinus => Some(lock_MusMusculusMolossinus()),
Species::MusMusculusMusculus => Some(lock_MusMusculusMusculus()),
Species::MusPahari => Some(lock_MusPahari()),
Species::MusSaxicola => Some(lock_MusSaxicola()),
Species::MusSp => Some(lock_MusSp()),
Species::MusSpretus => Some(lock_MusSpretus()),
Species::MustelaPutoriusFuro => Some(lock_MustelaPutoriusFuro()),
Species::NeogaleVison => Some(lock_NeogaleVison()),
Species::NototheniaCoriiceps => Some(lock_NototheniaCoriiceps()),
Species::OncorhynchusMykiss => Some(lock_OncorhynchusMykiss()),
Species::OrnithorhynchusAnatinus => Some(lock_OrnithorhynchusAnatinus()),
Species::OryctolagusCuniculus => Some(lock_OryctolagusCuniculus()),
Species::OryctolagusCuniculusAlgirus => Some(lock_OryctolagusCuniculusAlgirus()),
Species::OryctolagusCuniculusCuniculus => Some(lock_OryctolagusCuniculusCuniculus()),
Species::OvisAries => Some(lock_OvisAries()),
Species::PanTroglodytes => Some(lock_PanTroglodytes()),
Species::PapioAnubisAnubis => Some(lock_PapioAnubisAnubis()),
Species::PongoAbelii => Some(lock_PongoAbelii()),
Species::PongoPygmaeus => Some(lock_PongoPygmaeus()),
Species::ProtopterusAethiopicus => Some(lock_ProtopterusAethiopicus()),
Species::RajaEglanteria => Some(lock_RajaEglanteria()),
Species::RattusNorvegicus => Some(lock_RattusNorvegicus()),
Species::RattusRattus => Some(lock_RattusRattus()),
Species::SalmoSalar => Some(lock_SalmoSalar()),
Species::SalmoTrutta => Some(lock_SalmoTrutta()),
Species::SeriolaQuinqueradiata => Some(lock_SeriolaQuinqueradiata()),
Species::SinipercaChuatsi => Some(lock_SinipercaChuatsi()),
Species::SusScrofa => Some(lock_SusScrofa()),
Species::TrematomusBernacchii => Some(lock_TrematomusBernacchii()),
Species::VicugnaPacos => Some(lock_VicugnaPacos()),
Species::XenopusLaevisOrGilli => Some(lock_XenopusLaevisOrGilli()),
_=>None}}
/// Get all germlines in one iterator, see the main documentation for more information about the available germlines
pub fn all_germlines() -> impl std::iter::Iterator<Item = &'static Germlines> {
[
lock_AnarhichasMinor(),
lock_BosTaurus(),
lock_CamelusDromedarius(),
lock_CanisLupusFamiliaris(),
lock_CapraHircus(),
lock_CarcharhinusPlumbeus(),
lock_CercocebusAtys(),
lock_ChaenocephalusAceratus(),
lock_CyprinusCarpio(),
lock_DanioRerio(),
lock_DicentrarchusLabrax(),
lock_EquusCaballus(),
lock_FelisCatus(),
lock_GadusMorhua(),
lock_GallusGallus(),
lock_GasterosteusAculeatus(),
lock_GinglymostomaCirratum(),
lock_GorillaGorilla(),
lock_GorillaGorillaGorilla(),
lock_HeterodontusFrancisci(),
lock_HomoSapiens(),
lock_HydrolagusColliei(),
lock_HylobatesLar(),
lock_IctalurusPunctatus(),
lock_LemurCatta(),
lock_LeucorajaErinacea(),
lock_MacacaArctoides(),
lock_MacacaCyclopis(),
lock_MacacaFascicularis(),
lock_MacacaMulatta(),
lock_MacacaNemestrina(),
lock_MacacaSilenus(),
lock_MacacaThibetana(),
lock_MesocricetusAuratus(),
lock_MonodelphisDomestica(),
lock_MusCookii(),
lock_MusMinutoides(),
lock_MusMusculus(),
lock_MusMusculusCastaneus(),
lock_MusMusculusDomesticus(),
lock_MusMusculusMolossinus(),
lock_MusMusculusMusculus(),
lock_MusPahari(),
lock_MusSaxicola(),
lock_MusSp(),
lock_MusSpretus(),
lock_MustelaPutoriusFuro(),
lock_NeogaleVison(),
lock_NototheniaCoriiceps(),
lock_OncorhynchusMykiss(),
lock_OrnithorhynchusAnatinus(),
lock_OryctolagusCuniculus(),
lock_OryctolagusCuniculusAlgirus(),
lock_OryctolagusCuniculusCuniculus(),
lock_OvisAries(),
lock_PanTroglodytes(),
lock_PapioAnubisAnubis(),
lock_PongoAbelii(),
lock_PongoPygmaeus(),
lock_ProtopterusAethiopicus(),
lock_RajaEglanteria(),
lock_RattusNorvegicus(),
lock_RattusRattus(),
lock_SalmoSalar(),
lock_SalmoTrutta(),
lock_SeriolaQuinqueradiata(),
lock_SinipercaChuatsi(),
lock_SusScrofa(),
lock_TrematomusBernacchii(),
lock_VicugnaPacos(),
lock_XenopusLaevisOrGilli(),
].into_iter()
}
/// Get all germlines in one parallel iterator, see the main documentation for more information about the available germlines
#[cfg(feature = "rayon")]
use rayon::prelude::*;
#[cfg(feature = "rayon")]
pub fn par_germlines() -> impl rayon::prelude::ParallelIterator<Item = &'static Germlines> {
[
lock_AnarhichasMinor(),
lock_BosTaurus(),
lock_CamelusDromedarius(),
lock_CanisLupusFamiliaris(),
lock_CapraHircus(),
lock_CarcharhinusPlumbeus(),
lock_CercocebusAtys(),
lock_ChaenocephalusAceratus(),
lock_CyprinusCarpio(),
lock_DanioRerio(),
lock_DicentrarchusLabrax(),
lock_EquusCaballus(),
lock_FelisCatus(),
lock_GadusMorhua(),
lock_GallusGallus(),
lock_GasterosteusAculeatus(),
lock_GinglymostomaCirratum(),
lock_GorillaGorilla(),
lock_GorillaGorillaGorilla(),
lock_HeterodontusFrancisci(),
lock_HomoSapiens(),
lock_HydrolagusColliei(),
lock_HylobatesLar(),
lock_IctalurusPunctatus(),
lock_LemurCatta(),
lock_LeucorajaErinacea(),
lock_MacacaArctoides(),
lock_MacacaCyclopis(),
lock_MacacaFascicularis(),
lock_MacacaMulatta(),
lock_MacacaNemestrina(),
lock_MacacaSilenus(),
lock_MacacaThibetana(),
lock_MesocricetusAuratus(),
lock_MonodelphisDomestica(),
lock_MusCookii(),
lock_MusMinutoides(),
lock_MusMusculus(),
lock_MusMusculusCastaneus(),
lock_MusMusculusDomesticus(),
lock_MusMusculusMolossinus(),
lock_MusMusculusMusculus(),
lock_MusPahari(),
lock_MusSaxicola(),
lock_MusSp(),
lock_MusSpretus(),
lock_MustelaPutoriusFuro(),
lock_NeogaleVison(),
lock_NototheniaCoriiceps(),
lock_OncorhynchusMykiss(),
lock_OrnithorhynchusAnatinus(),
lock_OryctolagusCuniculus(),
lock_OryctolagusCuniculusAlgirus(),
lock_OryctolagusCuniculusCuniculus(),
lock_OvisAries(),
lock_PanTroglodytes(),
lock_PapioAnubisAnubis(),
lock_PongoAbelii(),
lock_PongoPygmaeus(),
lock_ProtopterusAethiopicus(),
lock_RajaEglanteria(),
lock_RattusNorvegicus(),
lock_RattusRattus(),
lock_SalmoSalar(),
lock_SalmoTrutta(),
lock_SeriolaQuinqueradiata(),
lock_SinipercaChuatsi(),
lock_SusScrofa(),
lock_TrematomusBernacchii(),
lock_VicugnaPacos(),
lock_XenopusLaevisOrGilli(),
].into_par_iter()
}
static LOCK_AnarhichasMinor: OnceLock<Germlines> = OnceLock::new();
fn lock_AnarhichasMinor()->&'static Germlines{LOCK_AnarhichasMinor.get_or_init(|| {bincode::deserialize(include_bytes!("Spotted wolffish.bin")).unwrap()})}
static LOCK_BosTaurus: OnceLock<Germlines> = OnceLock::new();
fn lock_BosTaurus()->&'static Germlines{LOCK_BosTaurus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic bovine.bin")).unwrap()})}
static LOCK_CamelusDromedarius: OnceLock<Germlines> = OnceLock::new();
fn lock_CamelusDromedarius()->&'static Germlines{LOCK_CamelusDromedarius.get_or_init(|| {bincode::deserialize(include_bytes!("Arabian camel.bin")).unwrap()})}
static LOCK_CanisLupusFamiliaris: OnceLock<Germlines> = OnceLock::new();
fn lock_CanisLupusFamiliaris()->&'static Germlines{LOCK_CanisLupusFamiliaris.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic dog.bin")).unwrap()})}
static LOCK_CapraHircus: OnceLock<Germlines> = OnceLock::new();
fn lock_CapraHircus()->&'static Germlines{LOCK_CapraHircus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic goat.bin")).unwrap()})}
static LOCK_CarcharhinusPlumbeus: OnceLock<Germlines> = OnceLock::new();
fn lock_CarcharhinusPlumbeus()->&'static Germlines{LOCK_CarcharhinusPlumbeus.get_or_init(|| {bincode::deserialize(include_bytes!("Sandbar shark.bin")).unwrap()})}
static LOCK_CercocebusAtys: OnceLock<Germlines> = OnceLock::new();
fn lock_CercocebusAtys()->&'static Germlines{LOCK_CercocebusAtys.get_or_init(|| {bincode::deserialize(include_bytes!("Sooty mangabey.bin")).unwrap()})}
static LOCK_ChaenocephalusAceratus: OnceLock<Germlines> = OnceLock::new();
fn lock_ChaenocephalusAceratus()->&'static Germlines{LOCK_ChaenocephalusAceratus.get_or_init(|| {bincode::deserialize(include_bytes!("Blackfin icefish.bin")).unwrap()})}
static LOCK_CyprinusCarpio: OnceLock<Germlines> = OnceLock::new();
fn lock_CyprinusCarpio()->&'static Germlines{LOCK_CyprinusCarpio.get_or_init(|| {bincode::deserialize(include_bytes!("Common carp.bin")).unwrap()})}
static LOCK_DanioRerio: OnceLock<Germlines> = OnceLock::new();
fn lock_DanioRerio()->&'static Germlines{LOCK_DanioRerio.get_or_init(|| {bincode::deserialize(include_bytes!("Zebrafish.bin")).unwrap()})}
static LOCK_DicentrarchusLabrax: OnceLock<Germlines> = OnceLock::new();
fn lock_DicentrarchusLabrax()->&'static Germlines{LOCK_DicentrarchusLabrax.get_or_init(|| {bincode::deserialize(include_bytes!("European seabass.bin")).unwrap()})}
static LOCK_EquusCaballus: OnceLock<Germlines> = OnceLock::new();
fn lock_EquusCaballus()->&'static Germlines{LOCK_EquusCaballus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic horse.bin")).unwrap()})}
static LOCK_FelisCatus: OnceLock<Germlines> = OnceLock::new();
fn lock_FelisCatus()->&'static Germlines{LOCK_FelisCatus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic cat.bin")).unwrap()})}
static LOCK_GadusMorhua: OnceLock<Germlines> = OnceLock::new();
fn lock_GadusMorhua()->&'static Germlines{LOCK_GadusMorhua.get_or_init(|| {bincode::deserialize(include_bytes!("Atlantic cod.bin")).unwrap()})}
static LOCK_GallusGallus: OnceLock<Germlines> = OnceLock::new();
fn lock_GallusGallus()->&'static Germlines{LOCK_GallusGallus.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic chicken.bin")).unwrap()})}
static LOCK_GasterosteusAculeatus: OnceLock<Germlines> = OnceLock::new();
fn lock_GasterosteusAculeatus()->&'static Germlines{LOCK_GasterosteusAculeatus.get_or_init(|| {bincode::deserialize(include_bytes!("Three-spined stickleback.bin")).unwrap()})}
static LOCK_GinglymostomaCirratum: OnceLock<Germlines> = OnceLock::new();
fn lock_GinglymostomaCirratum()->&'static Germlines{LOCK_GinglymostomaCirratum.get_or_init(|| {bincode::deserialize(include_bytes!("Nurse shark.bin")).unwrap()})}
static LOCK_GorillaGorilla: OnceLock<Germlines> = OnceLock::new();
fn lock_GorillaGorilla()->&'static Germlines{LOCK_GorillaGorilla.get_or_init(|| {bincode::deserialize(include_bytes!("Western gorilla.bin")).unwrap()})}
static LOCK_GorillaGorillaGorilla: OnceLock<Germlines> = OnceLock::new();
fn lock_GorillaGorillaGorilla()->&'static Germlines{LOCK_GorillaGorillaGorilla.get_or_init(|| {bincode::deserialize(include_bytes!("Western lowland gorilla.bin")).unwrap()})}
static LOCK_HeterodontusFrancisci: OnceLock<Germlines> = OnceLock::new();
fn lock_HeterodontusFrancisci()->&'static Germlines{LOCK_HeterodontusFrancisci.get_or_init(|| {bincode::deserialize(include_bytes!("Horn shark.bin")).unwrap()})}
static LOCK_HomoSapiens: OnceLock<Germlines> = OnceLock::new();
fn lock_HomoSapiens()->&'static Germlines{LOCK_HomoSapiens.get_or_init(|| {bincode::deserialize(include_bytes!("Human.bin")).unwrap()})}
static LOCK_HydrolagusColliei: OnceLock<Germlines> = OnceLock::new();
fn lock_HydrolagusColliei()->&'static Germlines{LOCK_HydrolagusColliei.get_or_init(|| {bincode::deserialize(include_bytes!("Spotted ratfish.bin")).unwrap()})}
static LOCK_HylobatesLar: OnceLock<Germlines> = OnceLock::new();
fn lock_HylobatesLar()->&'static Germlines{LOCK_HylobatesLar.get_or_init(|| {bincode::deserialize(include_bytes!("Common gibbon.bin")).unwrap()})}
static LOCK_IctalurusPunctatus: OnceLock<Germlines> = OnceLock::new();
fn lock_IctalurusPunctatus()->&'static Germlines{LOCK_IctalurusPunctatus.get_or_init(|| {bincode::deserialize(include_bytes!("Channel catfish.bin")).unwrap()})}
static LOCK_LemurCatta: OnceLock<Germlines> = OnceLock::new();
fn lock_LemurCatta()->&'static Germlines{LOCK_LemurCatta.get_or_init(|| {bincode::deserialize(include_bytes!("Ring-tailed lemur.bin")).unwrap()})}
static LOCK_LeucorajaErinacea: OnceLock<Germlines> = OnceLock::new();
fn lock_LeucorajaErinacea()->&'static Germlines{LOCK_LeucorajaErinacea.get_or_init(|| {bincode::deserialize(include_bytes!("Little skate.bin")).unwrap()})}
static LOCK_MacacaArctoides: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaArctoides()->&'static Germlines{LOCK_MacacaArctoides.get_or_init(|| {bincode::deserialize(include_bytes!("Stump-tailed macaque.bin")).unwrap()})}
static LOCK_MacacaCyclopis: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaCyclopis()->&'static Germlines{LOCK_MacacaCyclopis.get_or_init(|| {bincode::deserialize(include_bytes!("Taiwan macaque.bin")).unwrap()})}
static LOCK_MacacaFascicularis: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaFascicularis()->&'static Germlines{LOCK_MacacaFascicularis.get_or_init(|| {bincode::deserialize(include_bytes!("Crab-eating macaque.bin")).unwrap()})}
static LOCK_MacacaMulatta: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaMulatta()->&'static Germlines{LOCK_MacacaMulatta.get_or_init(|| {bincode::deserialize(include_bytes!("Rhesus monkey.bin")).unwrap()})}
static LOCK_MacacaNemestrina: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaNemestrina()->&'static Germlines{LOCK_MacacaNemestrina.get_or_init(|| {bincode::deserialize(include_bytes!("Pig-tailed macaque.bin")).unwrap()})}
static LOCK_MacacaSilenus: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaSilenus()->&'static Germlines{LOCK_MacacaSilenus.get_or_init(|| {bincode::deserialize(include_bytes!("Liontail macaque.bin")).unwrap()})}
static LOCK_MacacaThibetana: OnceLock<Germlines> = OnceLock::new();
fn lock_MacacaThibetana()->&'static Germlines{LOCK_MacacaThibetana.get_or_init(|| {bincode::deserialize(include_bytes!("Pere David's macaque.bin")).unwrap()})}
static LOCK_MesocricetusAuratus: OnceLock<Germlines> = OnceLock::new();
fn lock_MesocricetusAuratus()->&'static Germlines{LOCK_MesocricetusAuratus.get_or_init(|| {bincode::deserialize(include_bytes!("Golden hamster.bin")).unwrap()})}
static LOCK_MonodelphisDomestica: OnceLock<Germlines> = OnceLock::new();
fn lock_MonodelphisDomestica()->&'static Germlines{LOCK_MonodelphisDomestica.get_or_init(|| {bincode::deserialize(include_bytes!("Gray short-tailed opossum.bin")).unwrap()})}
static LOCK_MusCookii: OnceLock<Germlines> = OnceLock::new();
fn lock_MusCookii()->&'static Germlines{LOCK_MusCookii.get_or_init(|| {bincode::deserialize(include_bytes!("Cook's mouse.bin")).unwrap()})}
static LOCK_MusMinutoides: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMinutoides()->&'static Germlines{LOCK_MusMinutoides.get_or_init(|| {bincode::deserialize(include_bytes!("Southern African pygmy mouse.bin")).unwrap()})}
static LOCK_MusMusculus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMusculus()->&'static Germlines{LOCK_MusMusculus.get_or_init(|| {bincode::deserialize(include_bytes!("House mouse.bin")).unwrap()})}
static LOCK_MusMusculusCastaneus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMusculusCastaneus()->&'static Germlines{LOCK_MusMusculusCastaneus.get_or_init(|| {bincode::deserialize(include_bytes!("Southeastern Asian house mouse.bin")).unwrap()})}
static LOCK_MusMusculusDomesticus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMusculusDomesticus()->&'static Germlines{LOCK_MusMusculusDomesticus.get_or_init(|| {bincode::deserialize(include_bytes!("Western European house mouse.bin")).unwrap()})}
static LOCK_MusMusculusMolossinus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMusculusMolossinus()->&'static Germlines{LOCK_MusMusculusMolossinus.get_or_init(|| {bincode::deserialize(include_bytes!("Japanese wild mouse.bin")).unwrap()})}
static LOCK_MusMusculusMusculus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusMusculusMusculus()->&'static Germlines{LOCK_MusMusculusMusculus.get_or_init(|| {bincode::deserialize(include_bytes!("Eastern European house mouse.bin")).unwrap()})}
static LOCK_MusPahari: OnceLock<Germlines> = OnceLock::new();
fn lock_MusPahari()->&'static Germlines{LOCK_MusPahari.get_or_init(|| {bincode::deserialize(include_bytes!("Shrew mouse.bin")).unwrap()})}
static LOCK_MusSaxicola: OnceLock<Germlines> = OnceLock::new();
fn lock_MusSaxicola()->&'static Germlines{LOCK_MusSaxicola.get_or_init(|| {bincode::deserialize(include_bytes!("Spiny mouse.bin")).unwrap()})}
static LOCK_MusSp: OnceLock<Germlines> = OnceLock::new();
fn lock_MusSp()->&'static Germlines{LOCK_MusSp.get_or_init(|| {bincode::deserialize(include_bytes!("Mice.bin")).unwrap()})}
static LOCK_MusSpretus: OnceLock<Germlines> = OnceLock::new();
fn lock_MusSpretus()->&'static Germlines{LOCK_MusSpretus.get_or_init(|| {bincode::deserialize(include_bytes!("Western wild mouse.bin")).unwrap()})}
static LOCK_MustelaPutoriusFuro: OnceLock<Germlines> = OnceLock::new();
fn lock_MustelaPutoriusFuro()->&'static Germlines{LOCK_MustelaPutoriusFuro.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic ferret.bin")).unwrap()})}
static LOCK_NeogaleVison: OnceLock<Germlines> = OnceLock::new();
fn lock_NeogaleVison()->&'static Germlines{LOCK_NeogaleVison.get_or_init(|| {bincode::deserialize(include_bytes!("American mink.bin")).unwrap()})}
static LOCK_NototheniaCoriiceps: OnceLock<Germlines> = OnceLock::new();
fn lock_NototheniaCoriiceps()->&'static Germlines{LOCK_NototheniaCoriiceps.get_or_init(|| {bincode::deserialize(include_bytes!("Black rockcod.bin")).unwrap()})}
static LOCK_OncorhynchusMykiss: OnceLock<Germlines> = OnceLock::new();
fn lock_OncorhynchusMykiss()->&'static Germlines{LOCK_OncorhynchusMykiss.get_or_init(|| {bincode::deserialize(include_bytes!("Rainbow trout.bin")).unwrap()})}
static LOCK_OrnithorhynchusAnatinus: OnceLock<Germlines> = OnceLock::new();
fn lock_OrnithorhynchusAnatinus()->&'static Germlines{LOCK_OrnithorhynchusAnatinus.get_or_init(|| {bincode::deserialize(include_bytes!("Platypus.bin")).unwrap()})}
static LOCK_OryctolagusCuniculus: OnceLock<Germlines> = OnceLock::new();
fn lock_OryctolagusCuniculus()->&'static Germlines{LOCK_OryctolagusCuniculus.get_or_init(|| {bincode::deserialize(include_bytes!("Rabbit.bin")).unwrap()})}
static LOCK_OryctolagusCuniculusAlgirus: OnceLock<Germlines> = OnceLock::new();
fn lock_OryctolagusCuniculusAlgirus()->&'static Germlines{LOCK_OryctolagusCuniculusAlgirus.get_or_init(|| {bincode::deserialize(include_bytes!("European rabbit.bin")).unwrap()})}
static LOCK_OryctolagusCuniculusCuniculus: OnceLock<Germlines> = OnceLock::new();
fn lock_OryctolagusCuniculusCuniculus()->&'static Germlines{LOCK_OryctolagusCuniculusCuniculus.get_or_init(|| {bincode::deserialize(include_bytes!("Rabbit.bin")).unwrap()})}
static LOCK_OvisAries: OnceLock<Germlines> = OnceLock::new();
fn lock_OvisAries()->&'static Germlines{LOCK_OvisAries.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic sheep.bin")).unwrap()})}
static LOCK_PanTroglodytes: OnceLock<Germlines> = OnceLock::new();
fn lock_PanTroglodytes()->&'static Germlines{LOCK_PanTroglodytes.get_or_init(|| {bincode::deserialize(include_bytes!("Chimpanzee.bin")).unwrap()})}
static LOCK_PapioAnubisAnubis: OnceLock<Germlines> = OnceLock::new();
fn lock_PapioAnubisAnubis()->&'static Germlines{LOCK_PapioAnubisAnubis.get_or_init(|| {bincode::deserialize(include_bytes!("Olive baboon anubis.bin")).unwrap()})}
static LOCK_PongoAbelii: OnceLock<Germlines> = OnceLock::new();
fn lock_PongoAbelii()->&'static Germlines{LOCK_PongoAbelii.get_or_init(|| {bincode::deserialize(include_bytes!("Sumatran orangutan.bin")).unwrap()})}
static LOCK_PongoPygmaeus: OnceLock<Germlines> = OnceLock::new();
fn lock_PongoPygmaeus()->&'static Germlines{LOCK_PongoPygmaeus.get_or_init(|| {bincode::deserialize(include_bytes!("Bornean orangutan.bin")).unwrap()})}
static LOCK_ProtopterusAethiopicus: OnceLock<Germlines> = OnceLock::new();
fn lock_ProtopterusAethiopicus()->&'static Germlines{LOCK_ProtopterusAethiopicus.get_or_init(|| {bincode::deserialize(include_bytes!("Marbled lungfish.bin")).unwrap()})}
static LOCK_RajaEglanteria: OnceLock<Germlines> = OnceLock::new();
fn lock_RajaEglanteria()->&'static Germlines{LOCK_RajaEglanteria.get_or_init(|| {bincode::deserialize(include_bytes!("Clearnose skate.bin")).unwrap()})}
static LOCK_RattusNorvegicus: OnceLock<Germlines> = OnceLock::new();
fn lock_RattusNorvegicus()->&'static Germlines{LOCK_RattusNorvegicus.get_or_init(|| {bincode::deserialize(include_bytes!("Norway rat.bin")).unwrap()})}
static LOCK_RattusRattus: OnceLock<Germlines> = OnceLock::new();
fn lock_RattusRattus()->&'static Germlines{LOCK_RattusRattus.get_or_init(|| {bincode::deserialize(include_bytes!("Black rat.bin")).unwrap()})}
static LOCK_SalmoSalar: OnceLock<Germlines> = OnceLock::new();
fn lock_SalmoSalar()->&'static Germlines{LOCK_SalmoSalar.get_or_init(|| {bincode::deserialize(include_bytes!("Atlantic salmon.bin")).unwrap()})}
static LOCK_SalmoTrutta: OnceLock<Germlines> = OnceLock::new();
fn lock_SalmoTrutta()->&'static Germlines{LOCK_SalmoTrutta.get_or_init(|| {bincode::deserialize(include_bytes!("River trout.bin")).unwrap()})}
static LOCK_SeriolaQuinqueradiata: OnceLock<Germlines> = OnceLock::new();
fn lock_SeriolaQuinqueradiata()->&'static Germlines{LOCK_SeriolaQuinqueradiata.get_or_init(|| {bincode::deserialize(include_bytes!("Japanese amberjack.bin")).unwrap()})}
static LOCK_SinipercaChuatsi: OnceLock<Germlines> = OnceLock::new();
fn lock_SinipercaChuatsi()->&'static Germlines{LOCK_SinipercaChuatsi.get_or_init(|| {bincode::deserialize(include_bytes!("Mandarin fish.bin")).unwrap()})}
static LOCK_SusScrofa: OnceLock<Germlines> = OnceLock::new();
fn lock_SusScrofa()->&'static Germlines{LOCK_SusScrofa.get_or_init(|| {bincode::deserialize(include_bytes!("Domestic pig.bin")).unwrap()})}
static LOCK_TrematomusBernacchii: OnceLock<Germlines> = OnceLock::new();
fn lock_TrematomusBernacchii()->&'static Germlines{LOCK_TrematomusBernacchii.get_or_init(|| {bincode::deserialize(include_bytes!("Emerald rockcod.bin")).unwrap()})}
static LOCK_VicugnaPacos: OnceLock<Germlines> = OnceLock::new();
fn lock_VicugnaPacos()->&'static Germlines{LOCK_VicugnaPacos.get_or_init(|| {bincode::deserialize(include_bytes!("Alpaca.bin")).unwrap()})}
static LOCK_XenopusLaevisOrGilli: OnceLock<Germlines> = OnceLock::new();
fn lock_XenopusLaevisOrGilli()->&'static Germlines{LOCK_XenopusLaevisOrGilli.get_or_init(|| {bincode::deserialize(include_bytes!("African or Cape clawed frog.bin")).unwrap()})}
