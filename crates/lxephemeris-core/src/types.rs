//! Domain types, `BodyId` namespaces, `CalcOptions`, `HouseSystem` (23 codes),
//! `AyanamsaModel` (12 codes), and the 68 minor bodies registry (`REVIEW_AND_SPEC_ASTRA.md`).

use crate::error::EphemerisError;
use serde::{Deserialize, Serialize};

pub use crate::events::body_display_name;

/// Exact IAU 2012 astronomical unit in kilometers.
pub const AU_KM: f64 = 149_597_870.700;

/// Exact SI seconds per Julian day.
pub const SECONDS_PER_DAY: f64 = 86_400.0;

/// Julian Date of J2000.0 epoch (`2000-01-01T12:00:00 TT/TDB`).
pub const JD_J2000: f64 = 2_451_545.0;

/// Julian Date of `1850-01-01T00:00:00 TT` (`2433282.5` is 1950; `1850-01-01T00:00:00` is `2396758.5`).
pub const JD_1850_01_01: f64 = 2_396_758.5;

/// Julian Date of `2100-01-01T00:00:00 TT` (`2488069.5`).
pub const JD_2100_01_01: f64 = 2_488_069.5;

/// Separate typed identifier spaces for celestial objects (`REVIEW_AND_SPEC_ASTRA.md` §5).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BodyId {
    Planet(PlanetId),
    MinorPlanet(MpcNumber),
    LunarPoint(LunarPointId),
    Hypothetical(HypotheticalId),
    FixedStar(HipId),
    /// Ecliptic obliquity and IAU 2000A nutation angles (`SE_ECL_NUT` / `R2.5`).
    EclNut,
}

/// Major solar system bodies and reference barycenters in `DE440s`.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PlanetId {
    Sun,
    Moon,
    Mercury,
    Venus,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
    EarthMoonBarycenter,
    SolarSystemBarycenter,
}

/// Physical representation of the target in the underlying ephemeris kernel.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TargetKind {
    /// Physical center of mass of the body (`10`, `199`, `299`, `399`, `301`).
    BodyCenter,
    /// Barycenter of the planetary system (`0`, `3`, `4..=9` in `de440s.bsp`).
    SystemBarycenter,
    /// Minor body center of mass (`2000000 + MPC` or Horizons SPK ID).
    MinorBodyCenter,
    /// Mathematical / osculating / mean / interpolated point.
    SyntheticPoint,
    /// Astrometric catalog star (`VizieR I/311` / `I/239` + `XHIP V/137D`).
    CatalogStar,
}

impl PlanetId {
    /// Returns the NAIF SPK target ID used in `de440s.bsp`.
    ///
    /// Per `ARCHITECTURE.md` §2.2.3:
    /// - `Mercury (199)`, `Venus (299)`, `Earth (399)`, `Moon (301)`, `Sun (10)` use `BodyCenter`
    /// - `Mars (4)`, `Jupiter (5)`, `Saturn (6)`, `Uranus (7)`, `Neptune (8)`, `Pluto (9)`,
    ///   `EMB (3)`, `SSB (0)` use `SystemBarycenter`.
    #[must_use]
    pub const fn naif_target_id(self) -> i32 {
        match self {
            Self::SolarSystemBarycenter => 0,
            Self::EarthMoonBarycenter => 3,
            Self::Mars => 4,
            Self::Jupiter => 5,
            Self::Saturn => 6,
            Self::Uranus => 7,
            Self::Neptune => 8,
            Self::Pluto => 9,
            Self::Sun => 10,
            Self::Mercury => 199,
            Self::Venus => 299,
            Self::Moon => 301,
            Self::Earth => 399,
        }
    }

    /// Returns the native NAIF center ID of this segment in `de440s.bsp`.
    #[must_use]
    pub const fn native_de440s_center_id(self) -> i32 {
        match self {
            Self::Mercury => 1,
            Self::Venus => 2,
            Self::Earth | Self::Moon => 3,
            Self::SolarSystemBarycenter
            | Self::EarthMoonBarycenter
            | Self::Mars
            | Self::Jupiter
            | Self::Saturn
            | Self::Uranus
            | Self::Neptune
            | Self::Pluto
            | Self::Sun => 0,
        }
    }

    /// Returns the [`TargetKind`] in `de440s.bsp`.
    #[must_use]
    pub const fn target_kind(self) -> TargetKind {
        match self {
            Self::Sun | Self::Moon | Self::Mercury | Self::Venus | Self::Earth => {
                TargetKind::BodyCenter
            }
            Self::Mars
            | Self::Jupiter
            | Self::Saturn
            | Self::Uranus
            | Self::Neptune
            | Self::Pluto
            | Self::EarthMoonBarycenter
            | Self::SolarSystemBarycenter => TargetKind::SystemBarycenter,
        }
    }

    /// Canonical name of the body.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Sun => "Sun",
            Self::Moon => "Moon",
            Self::Mercury => "Mercury",
            Self::Venus => "Venus",
            Self::Earth => "Earth",
            Self::Mars => "Mars",
            Self::Jupiter => "Jupiter",
            Self::Saturn => "Saturn",
            Self::Uranus => "Uranus",
            Self::Neptune => "Neptune",
            Self::Pluto => "Pluto",
            Self::EarthMoonBarycenter => "EarthMoonBarycenter",
            Self::SolarSystemBarycenter => "SolarSystemBarycenter",
        }
    }
}

/// Minor Planet Center (MPC) catalog number (e.g. `1` for Ceres, `2060` for Chiron, `90377` for Sedna).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct MpcNumber(pub u32);

impl MpcNumber {
    /// Standard 8-digit NAIF SPK ID for numbered minor planets in JPL Horizons (`20_000_000 + MPC`).
    /// Note: Actual runtime lookup also validates against `manifest.json` / DAF segment descriptors.
    #[must_use]
    pub const fn standard_naif_spk_id(self) -> i32 {
        20_000_000_i32.saturating_add(self.0 as i32)
    }
}

/// Category in the `astra` minor bodies registry (`REVIEW_AND_SPEC_ASTRA.md` §6).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MinorBodyCategory {
    Asteroid,
    DwarfPlanet,
}

impl MinorBodyCategory {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Asteroid => "asteroid",
            Self::DwarfPlanet => "dwarf_planet",
        }
    }
}

/// Static descriptor for each of the 68 minor bodies in `REVIEW_AND_SPEC_ASTRA.md` §6.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct MinorBodySpec {
    pub mpc: MpcNumber,
    pub name: &'static str,
    pub category: MinorBodyCategory,
}

/// Complete registry of all 68 minor bodies required by `astra` (`REVIEW_AND_SPEC_ASTRA.md` §6).
pub const MINOR_BODIES_68: [MinorBodySpec; 68] = [
    MinorBodySpec { mpc: MpcNumber(1), name: "Ceres", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(2), name: "Pallas", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(3), name: "Juno", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(4), name: "Vesta", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(5), name: "Astraea", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(6), name: "Hebe", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(7), name: "Iris", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(8), name: "Flora", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(9), name: "Metis", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(10), name: "Hygiea", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(11), name: "Parthenope", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(12), name: "Victoria", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(13), name: "Egeria", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(14), name: "Irene", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(15), name: "Eunomia", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(16), name: "Psyche", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(18), name: "Melpomene", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(19), name: "Fortuna", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(26), name: "Proserpina", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(30), name: "Urania", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(37), name: "Fides", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(40), name: "Harmonia", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(42), name: "Isis", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(55), name: "Pandora", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(78), name: "Diana", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(80), name: "Sappho", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(100), name: "Hekate", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(105), name: "Artemis", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(114), name: "Kassandra", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(128), name: "Nemesis", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(227), name: "Philosophia", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(251), name: "Sophia", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(259), name: "Aletheia", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(344), name: "Desiderata", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(375), name: "Ursula", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(399), name: "Persephone", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(407), name: "Arachne", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(433), name: "Eros", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(447), name: "Valentine", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(471), name: "Papagena", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(490), name: "Veritas", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(504), name: "Cora", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(763), name: "Cupido", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(944), name: "Hidalgo", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(1181), name: "Lilith", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(1221), name: "Amor", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(1387), name: "Kama", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(1388), name: "Aphrodite", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(1489), name: "Attila", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(1566), name: "Icarus", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(1580), name: "Betulia", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(1585), name: "Union", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(1862), name: "Apollo", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(1930), name: "Lucifer", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(2060), name: "Chiron", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(2101), name: "Adonis", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(3671), name: "Dionysus", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(3782), name: "Celle", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(3811), name: "Karma", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(4386), name: "Lust", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(5129), name: "Groom", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(5145), name: "Pholus", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(6583), name: "Destinn", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(7066), name: "Nessus", category: MinorBodyCategory::Asteroid },
    MinorBodySpec { mpc: MpcNumber(90377), name: "Sedna", category: MinorBodyCategory::DwarfPlanet },
    MinorBodySpec { mpc: MpcNumber(136108), name: "Haumea", category: MinorBodyCategory::DwarfPlanet },
    MinorBodySpec { mpc: MpcNumber(136199), name: "Eris", category: MinorBodyCategory::DwarfPlanet },
    MinorBodySpec { mpc: MpcNumber(136472), name: "Makemake", category: MinorBodyCategory::DwarfPlanet },
];

/// Calculated lunar node and apogee points (`REVIEW_AND_SPEC_ASTRA.md` §9).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum LunarPointId {
    MeanNode,
    TrueNode,
    MeanSouthNode,
    TrueSouthNode,
    MeanLilith,
    TrueLilith,
    InterpolatedLilith,
    InterpolatedPerigee,
}

impl LunarPointId {
    #[must_use]
    pub const fn astra_id(self) -> &'static str {
        match self {
            Self::MeanNode => "Mean_Node",
            Self::TrueNode => "True_Node",
            Self::MeanSouthNode => "Mean_South_Node",
            Self::TrueSouthNode => "True_South_Node",
            Self::MeanLilith => "Mean_Lilith",
            Self::TrueLilith => "True_Lilith",
            Self::InterpolatedLilith => "Interpolated_Lilith",
            Self::InterpolatedPerigee => "Interpolated_Perigee",
        }
    }
}

/// 8 Uranian (Hamburg School) planets + `Transpluto` (`Isis`) (`REVIEW_AND_SPEC_ASTRA.md` §10).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum HypotheticalId {
    Cupido,
    Hades,
    Zeus,
    Kronos,
    Apollon,
    Admetos,
    Vulkanus,
    Poseidon,
    Transpluto,
}

impl HypotheticalId {
    #[must_use]
    pub const fn astra_id(self) -> &'static str {
        match self {
            Self::Cupido => "Cupido",
            Self::Hades => "Hades",
            Self::Zeus => "Zeus",
            Self::Kronos => "Kronos",
            Self::Apollon => "Apollon",
            Self::Admetos => "Admetos",
            Self::Vulkanus => "Vulkanus",
            Self::Poseidon => "Poseidon",
            Self::Transpluto => "Transpluto",
        }
    }
}

/// Hipparcos catalog identifier (`VizieR I/311` / `I/239`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct HipId(pub u32);

/// Geodetic observer on WGS84 ellipsoid (`REVIEW_AND_SPEC_ASTRA.md` §5).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Observer {
    /// Geodetic latitude in degrees (`-90.0..=90.0`, north positive).
    pub latitude_deg: f64,
    /// Geodetic longitude in degrees (`-180.0..=180.0`, east positive).
    pub longitude_deg: f64,
    /// Height above WGS84 reference ellipsoid in meters.
    pub altitude_m: f64,
}

impl Observer {
    /// Creates and validates a geodetic WGS84 [`Observer`].
    pub fn new(
        latitude_deg: f64,
        longitude_deg: f64,
        altitude_m: f64,
    ) -> Result<Self, EphemerisError> {
        if !latitude_deg.is_finite()
            || !(-90.0..=90.0).contains(&latitude_deg)
            || !longitude_deg.is_finite()
            || !(-360.0..=360.0).contains(&longitude_deg)
            || !altitude_m.is_finite()
            || !(-15_000.0..=100_000_000.0).contains(&altitude_m)
        {
            return Err(EphemerisError::InvalidInput {
                message: format!(
                    "Invalid observer coordinates: lat={latitude_deg}, lon={longitude_deg}, alt={altitude_m}"
                ),
            });
        }
        Ok(Self {
            latitude_deg,
            longitude_deg,
            altitude_m,
        })
    }
}

/// Coordinate origin for astrometric calculations.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub enum Origin {
    #[default]
    Geocentric,
    Topocentric(Observer),
    Heliocentric,
    Barycentric,
}

/// Fundamental reference plane.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum Plane {
    #[default]
    Ecliptic,
    Equatorial,
    Horizontal,
}

/// Reference equinox for coordinate axes.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum Equinox {
    #[default]
    OfDate,
    J2000,
}

/// Nutation mode (`True` includes IAU 2000A nutation $\Delta\psi, \Delta\varepsilon$; `Mean` omits nutation).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum NutationMode {
    #[default]
    True,
    Mean,
}

/// Aberration and light-time correction mode.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum AbMode {
    /// Light-time + relativistic gravitational deflection + annual/diurnal aberration.
    #[default]
    Apparent,
    /// Light-time + relativistic gravitational deflection, without stellar aberration (`SEFLG_NOABERR`).
    NoAberration,
    /// Light-time + annual/diurnal aberration, without gravitational deflection (`SEFLG_NOGDEFL`).
    NoDeflection,
    /// Light-time only (no gravitational deflection or stellar aberration, `SEFLG_NOABERR | SEFLG_NOGDEFL`).
    Astrometric,
    /// Instantaneous geometric vector at epoch $t$ ($\tau = 0$, `SEFLG_TRUEPOS`).
    Geometric,
}

/// Zodiacal reference mode (`Tropical` or `Sidereal(AyanamsaModel)`).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum Zodiac {
    #[default]
    Tropical,
    Sidereal(AyanamsaModel),
}

/// Complete calculation options for a single or batch astrometric query (`ARCHITECTURE.md` §2.2.3).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CalcOptions {
    pub origin: Origin,
    pub plane: Plane,
    pub equinox: Equinox,
    pub nutation: NutationMode,
    pub aberration: AbMode,
    pub zodiac: Zodiac,
    pub include_speed: bool,
}

impl Default for CalcOptions {
    fn default() -> Self {
        Self {
            origin: Origin::Geocentric,
            plane: Plane::Ecliptic,
            equinox: Equinox::OfDate,
            nutation: NutationMode::True,
            aberration: AbMode::Apparent,
            zodiac: Zodiac::Tropical,
            include_speed: true,
        }
    }
}

/// All 12 Vedic / Sidereal ayanamsa models required by `astra` (`REVIEW_AND_SPEC_ASTRA.md` §8).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum AyanamsaModel {
    /// `0` — Fagan-Bradley
    FaganBradley,
    /// `1` — Lahiri (Chitra Paksha, IAE 1989)
    #[default]
    Lahiri,
    /// `2` — De Luce
    DeLuce,
    /// `3` — B. V. Raman
    Raman,
    /// `4` — Usha-Shashi
    Ushashashi,
    /// `5` — Krishnamurti (KP)
    Krishnamurti,
    /// `6` — Djwhal Khul
    DjwhalKhul,
    /// `7` — Sri Yukteshwar
    Yukteshwar,
    /// `8` — J. N. Bhasin
    JnBhasin,
    /// `9` — Babylonian (Kugler 1)
    BabylonKugler1,
    /// `10` — Babylonian (Kugler 2)
    BabylonKugler2,
    /// `11` — Babylonian (Kugler 3)
    BabylonKugler3,
    /// `12` — Babylonian (Huber)
    BabylonHuber,
    /// `13` — Babylonian (Eta Piscium)
    BabylonEtpsc,
    /// `14` — Babylonian (Aldebaran at 15° Taurus)
    Aldebaran15Tau,
    /// `15` — Hipparchos
    Hipparchos,
    /// `16` — Sassanian
    Sassanian,
    /// `17` — Galactic Center at 0° Sagittarius
    GalCent0Sag,
    /// `18` — J2000.0
    J2000,
    /// `19` — J1900.0
    J1900,
    /// `20` — B1950.0
    B1950,
    /// `21` — Surya Siddhanta
    Suryasiddhanta,
    /// `22` — Surya Siddhanta (Mean Sun)
    SuryasiddhantaMsun,
    /// `23` — Aryabhata
    Aryabhata,
    /// `24` — Aryabhata (Mean Sun)
    AryabhataMsun,
    /// `25` — Surya Siddhanta Revati
    SsRevati,
    /// `26` — Surya Siddhanta Citra
    SsCitra,
    /// `27` — True Citra (Spica at 180°)
    TrueCitra,
    /// `28` — True Revati ($\zeta$ Psc at 359°50')
    TrueRevati,
    /// `29` — True Pushya ($\delta$ Cnc at 106°)
    TruePushya,
    /// `30` — Galactic Center (Gil Brand)
    GalCentRGilbrand,
    /// `31` — Galactic Equator (IAU 1958)
    GalEqIau1958,
    /// `32` — Galactic Equator (True)
    GalEqTrue,
    /// `33` — Galactic Equator (Mula)
    GalEqMula,
    /// `34` — Galactic Alignment (Mardyks)
    GalAlignMardyks,
    /// `35` — True Mula (Chandra Hari)
    TrueMula,
    /// `36` — Galactic Center Mula (Wilhelm)
    GalCentMulaWilhelm,
    /// `37` — Aryabhata 522 CE
    Aryabhata522,
    /// `38` — Babylonian (Britton)
    BabylonBritton,
    /// `39` — True Sheratan
    TrueSheratan,
    /// `40` — True Aldebaran
    TrueAldebaran,
    /// `41` — Galactic Equator (Fiorenza)
    GalEquFiorenza,
    /// `42` — Vettius Valens Moon
    ValensMoon,
    /// `255` — Custom user-defined epoch anchor `(t0_jd_tt, ayan_t0_deg)` stored as IEEE-754 bits.
    User {
        t0_jd_tt_bits: u64,
        ayan_t0_deg_bits: u64,
    },
}

impl AyanamsaModel {
    /// Creates a custom user-defined ayanamsa (`SE_SIDM_USER = 255`) anchored at `(t0_jd_tt, ayan_t0_deg)`.
    pub fn user(t0_jd_tt: f64, ayan_t0_deg: f64) -> Result<Self, EphemerisError> {
        if !t0_jd_tt.is_finite() || !ayan_t0_deg.is_finite() {
            return Err(EphemerisError::InvalidInput {
                message: format!(
                    "Custom ayanamsa requires finite (t0_jd_tt, ayan_t0_deg), got ({t0_jd_tt}, {ayan_t0_deg})"
                ),
            });
        }
        let norm_bits = |x: f64| if x == 0.0 { 0.0_f64.to_bits() } else { x.to_bits() };
        Ok(Self::User {
            t0_jd_tt_bits: norm_bits(t0_jd_tt),
            ayan_t0_deg_bits: norm_bits(ayan_t0_deg),
        })
    }

    /// All 12 canonical `astra` ayanamsa models in registry order (`REVIEW_AND_SPEC_ASTRA.md` §8).
    pub const ALL: [Self; 12] = [
        Self::Lahiri,
        Self::Krishnamurti,
        Self::Raman,
        Self::Yukteshwar,
        Self::TrueCitra,
        Self::JnBhasin,
        Self::Ushashashi,
        Self::TrueRevati,
        Self::TruePushya,
        Self::SsCitra,
        Self::Suryasiddhanta,
        Self::FaganBradley,
    ];

    /// All 43 named Swiss Ephemeris `SE_SIDM_*` ayanamsa models (`0..=42`).
    pub const ALL_EXTENDED: [Self; 43] = [
        Self::FaganBradley,
        Self::Lahiri,
        Self::DeLuce,
        Self::Raman,
        Self::Ushashashi,
        Self::Krishnamurti,
        Self::DjwhalKhul,
        Self::Yukteshwar,
        Self::JnBhasin,
        Self::BabylonKugler1,
        Self::BabylonKugler2,
        Self::BabylonKugler3,
        Self::BabylonHuber,
        Self::BabylonEtpsc,
        Self::Aldebaran15Tau,
        Self::Hipparchos,
        Self::Sassanian,
        Self::GalCent0Sag,
        Self::J2000,
        Self::J1900,
        Self::B1950,
        Self::Suryasiddhanta,
        Self::SuryasiddhantaMsun,
        Self::Aryabhata,
        Self::AryabhataMsun,
        Self::SsRevati,
        Self::SsCitra,
        Self::TrueCitra,
        Self::TrueRevati,
        Self::TruePushya,
        Self::GalCentRGilbrand,
        Self::GalEqIau1958,
        Self::GalEqTrue,
        Self::GalEqMula,
        Self::GalAlignMardyks,
        Self::TrueMula,
        Self::GalCentMulaWilhelm,
        Self::Aryabhata522,
        Self::BabylonBritton,
        Self::TrueSheratan,
        Self::TrueAldebaran,
        Self::GalEquFiorenza,
        Self::ValensMoon,
    ];

    /// Parses an `astra` or `SE_SIDM_*` ayanamsa string identifier strictly.
    /// Unknown names return [`EphemerisError::UnsupportedModel`] without falling back to Lahiri.
    pub fn from_code(code: &str) -> Result<Self, EphemerisError> {
        match code {
            "fagan_bradley" => Ok(Self::FaganBradley),
            "lahiri" => Ok(Self::Lahiri),
            "deluce" => Ok(Self::DeLuce),
            "raman" => Ok(Self::Raman),
            "ushashashi" => Ok(Self::Ushashashi),
            "krishnamurti" => Ok(Self::Krishnamurti),
            "djwhal_khul" => Ok(Self::DjwhalKhul),
            "yukteshwar" => Ok(Self::Yukteshwar),
            "jn_bhasin" => Ok(Self::JnBhasin),
            "babyl_kugler1" => Ok(Self::BabylonKugler1),
            "babyl_kugler2" => Ok(Self::BabylonKugler2),
            "babyl_kugler3" => Ok(Self::BabylonKugler3),
            "babyl_huber" => Ok(Self::BabylonHuber),
            "babyl_etpsc" => Ok(Self::BabylonEtpsc),
            "aldebaran_15tau" => Ok(Self::Aldebaran15Tau),
            "hipparchos" => Ok(Self::Hipparchos),
            "sassanian" => Ok(Self::Sassanian),
            "galcent_0sag" => Ok(Self::GalCent0Sag),
            "j2000" => Ok(Self::J2000),
            "j1900" => Ok(Self::J1900),
            "b1950" => Ok(Self::B1950),
            "suryasiddhanta" => Ok(Self::Suryasiddhanta),
            "suryasiddhanta_msun" => Ok(Self::SuryasiddhantaMsun),
            "aryabhata" => Ok(Self::Aryabhata),
            "aryabhata_msun" => Ok(Self::AryabhataMsun),
            "ss_revati" => Ok(Self::SsRevati),
            "ss_citra" => Ok(Self::SsCitra),
            "true_citra" => Ok(Self::TrueCitra),
            "true_revati" => Ok(Self::TrueRevati),
            "true_pushya" => Ok(Self::TruePushya),
            "galcent_rgilbrand" => Ok(Self::GalCentRGilbrand),
            "galeq_iau1958" => Ok(Self::GalEqIau1958),
            "galeq_true" => Ok(Self::GalEqTrue),
            "galeq_mula" => Ok(Self::GalEqMula),
            "galalign_mardyks" => Ok(Self::GalAlignMardyks),
            "true_mula" => Ok(Self::TrueMula),
            "galcent_mula_wilhelm" => Ok(Self::GalCentMulaWilhelm),
            "aryabhata_522" => Ok(Self::Aryabhata522),
            "babyl_britton" => Ok(Self::BabylonBritton),
            "true_sheratan" => Ok(Self::TrueSheratan),
            "true_aldebaran" => Ok(Self::TrueAldebaran),
            "galequ_fiorenza" => Ok(Self::GalEquFiorenza),
            "valens_moon" => Ok(Self::ValensMoon),
            other => Err(EphemerisError::UnsupportedModel {
                model: other.to_owned(),
            }),
        }
    }

    /// Returns the canonical string identifier.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::FaganBradley => "fagan_bradley",
            Self::Lahiri => "lahiri",
            Self::DeLuce => "deluce",
            Self::Raman => "raman",
            Self::Ushashashi => "ushashashi",
            Self::Krishnamurti => "krishnamurti",
            Self::DjwhalKhul => "djwhal_khul",
            Self::Yukteshwar => "yukteshwar",
            Self::JnBhasin => "jn_bhasin",
            Self::BabylonKugler1 => "babyl_kugler1",
            Self::BabylonKugler2 => "babyl_kugler2",
            Self::BabylonKugler3 => "babyl_kugler3",
            Self::BabylonHuber => "babyl_huber",
            Self::BabylonEtpsc => "babyl_etpsc",
            Self::Aldebaran15Tau => "aldebaran_15tau",
            Self::Hipparchos => "hipparchos",
            Self::Sassanian => "sassanian",
            Self::GalCent0Sag => "galcent_0sag",
            Self::J2000 => "j2000",
            Self::J1900 => "j1900",
            Self::B1950 => "b1950",
            Self::Suryasiddhanta => "suryasiddhanta",
            Self::SuryasiddhantaMsun => "suryasiddhanta_msun",
            Self::Aryabhata => "aryabhata",
            Self::AryabhataMsun => "aryabhata_msun",
            Self::SsRevati => "ss_revati",
            Self::SsCitra => "ss_citra",
            Self::TrueCitra => "true_citra",
            Self::TrueRevati => "true_revati",
            Self::TruePushya => "true_pushya",
            Self::GalCentRGilbrand => "galcent_rgilbrand",
            Self::GalEqIau1958 => "galeq_iau1958",
            Self::GalEqTrue => "galeq_true",
            Self::GalEqMula => "galeq_mula",
            Self::GalAlignMardyks => "galalign_mardyks",
            Self::TrueMula => "true_mula",
            Self::GalCentMulaWilhelm => "galcent_mula_wilhelm",
            Self::Aryabhata522 => "aryabhata_522",
            Self::BabylonBritton => "babyl_britton",
            Self::TrueSheratan => "true_sheratan",
            Self::TrueAldebaran => "true_aldebaran",
            Self::GalEquFiorenza => "galequ_fiorenza",
            Self::ValensMoon => "valens_moon",
            Self::User { .. } => "user",
        }
    }

    /// Returns the legacy numerical crosswalk ID (`SE_SIDM_*`).
    #[must_use]
    pub const fn legacy_id(self) -> u32 {
        match self {
            Self::FaganBradley => 0,
            Self::Lahiri => 1,
            Self::DeLuce => 2,
            Self::Raman => 3,
            Self::Ushashashi => 4,
            Self::Krishnamurti => 5,
            Self::DjwhalKhul => 6,
            Self::Yukteshwar => 7,
            Self::JnBhasin => 8,
            Self::BabylonKugler1 => 9,
            Self::BabylonKugler2 => 10,
            Self::BabylonKugler3 => 11,
            Self::BabylonHuber => 12,
            Self::BabylonEtpsc => 13,
            Self::Aldebaran15Tau => 14,
            Self::Hipparchos => 15,
            Self::Sassanian => 16,
            Self::GalCent0Sag => 17,
            Self::J2000 => 18,
            Self::J1900 => 19,
            Self::B1950 => 20,
            Self::Suryasiddhanta => 21,
            Self::SuryasiddhantaMsun => 22,
            Self::Aryabhata => 23,
            Self::AryabhataMsun => 24,
            Self::SsRevati => 25,
            Self::SsCitra => 26,
            Self::TrueCitra => 27,
            Self::TrueRevati => 28,
            Self::TruePushya => 29,
            Self::GalCentRGilbrand => 30,
            Self::GalEqIau1958 => 31,
            Self::GalEqTrue => 32,
            Self::GalEqMula => 33,
            Self::GalAlignMardyks => 34,
            Self::TrueMula => 35,
            Self::GalCentMulaWilhelm => 36,
            Self::Aryabhata522 => 37,
            Self::BabylonBritton => 38,
            Self::TrueSheratan => 39,
            Self::TrueAldebaran => 40,
            Self::GalEquFiorenza => 41,
            Self::ValensMoon => 42,
            Self::User { .. } => 255,
        }
    }

    /// Resolves a legacy numerical `SE_SIDM_*` ID (`0..=42` or `255`) into an [`AyanamsaModel`].
    pub fn from_legacy_id(id: u32) -> Result<Self, EphemerisError> {
        match id {
            0 => Ok(Self::FaganBradley),
            1 => Ok(Self::Lahiri),
            2 => Ok(Self::DeLuce),
            3 => Ok(Self::Raman),
            4 => Ok(Self::Ushashashi),
            5 => Ok(Self::Krishnamurti),
            6 => Ok(Self::DjwhalKhul),
            7 => Ok(Self::Yukteshwar),
            8 => Ok(Self::JnBhasin),
            9 => Ok(Self::BabylonKugler1),
            10 => Ok(Self::BabylonKugler2),
            11 => Ok(Self::BabylonKugler3),
            12 => Ok(Self::BabylonHuber),
            13 => Ok(Self::BabylonEtpsc),
            14 => Ok(Self::Aldebaran15Tau),
            15 => Ok(Self::Hipparchos),
            16 => Ok(Self::Sassanian),
            17 => Ok(Self::GalCent0Sag),
            18 => Ok(Self::J2000),
            19 => Ok(Self::J1900),
            20 => Ok(Self::B1950),
            21 => Ok(Self::Suryasiddhanta),
            22 => Ok(Self::SuryasiddhantaMsun),
            23 => Ok(Self::Aryabhata),
            24 => Ok(Self::AryabhataMsun),
            25 => Ok(Self::SsRevati),
            26 => Ok(Self::SsCitra),
            27 => Ok(Self::TrueCitra),
            28 => Ok(Self::TrueRevati),
            29 => Ok(Self::TruePushya),
            30 => Ok(Self::GalCentRGilbrand),
            31 => Ok(Self::GalEqIau1958),
            32 => Ok(Self::GalEqTrue),
            33 => Ok(Self::GalEqMula),
            34 => Ok(Self::GalAlignMardyks),
            35 => Ok(Self::TrueMula),
            36 => Ok(Self::GalCentMulaWilhelm),
            37 => Ok(Self::Aryabhata522),
            38 => Ok(Self::BabylonBritton),
            39 => Ok(Self::TrueSheratan),
            40 => Ok(Self::TrueAldebaran),
            41 => Ok(Self::GalEquFiorenza),
            42 => Ok(Self::ValensMoon),
            43 => Self::user(JD_J2000, 23.842_323_260_326_51),
            44 => Self::user(JD_J2000, 23.863_481_230_643_23),
            45 => Self::user(JD_J2000, 23.780_364_984_916_844),
            46 => Self::user(JD_J2000, 23.856_789_016_285_802),
            255 => Self::user(JD_J2000, 0.0),
            other => Err(EphemerisError::UnsupportedModel {
                model: format!("SE_SIDM({other})"),
            }),
        }
    }
}

/// Polar behavior category for house systems (`REVIEW_AND_SPEC_ASTRA.md` §7).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PolarCategory {
    HardFail,
    SoftFail,
    PolarSafe,
}

/// Policy for handling circumpolar / high-latitude house calculation failures (`ARCHITECTURE.md` §2.4 & `R4.4`).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum PolarPolicy {
    /// Default policy: return `Err(EphemerisError::UndefinedHouseSystem)` for `HardFail` systems
    /// above the polar limit or when semi-arc equations do not converge.
    #[default]
    Strict,
    /// Explicitly requested fallback: substitute Porphyry (`O`) when a `HardFail` system is undefined,
    /// recording `requested_house_system`, `effective_house_system`, and `HouseWarning::PolarFallbackToPorphyry`.
    FallbackToPorphyry,
    /// Astra production polar policy (`R4.4`): substitute Whole Sign (`W`) at the true latitude (without
    /// latitude clamping) when a `HardFail` system (`P, K, i`) is undefined inside the polar circle OR
    /// when a `SoftFail` system (`C, R, B, T, L, I, Y`) produces a non-monotonic (inverted) 12-cusp ring.
    FallbackToWholeSign,
}

/// Diagnostic warning emitted during house cusp computation.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum HouseWarning {
    PolarFallbackToPorphyry,
    PolarFallbackToWholeSign,
    PolarUnreliable,
    NearPoleBoundary,
}

/// All 23 case-sensitive house systems in `astra` (`REVIEW_AND_SPEC_ASTRA.md` §7).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum HouseSystem {
    /// `P` — Placidus (`hard_fail`)
    #[default]
    Placidus,
    /// `K` — Koch (`hard_fail`)
    Koch,
    /// `O` — Porphyry (`polar_safe`)
    Porphyry,
    /// `R` — Regiomontanus (`soft_fail`)
    Regiomontanus,
    /// `C` — Campanus (`soft_fail`)
    Campanus,
    /// `A` — Equal (from ASC) (`polar_safe`)
    EqualFromAsc,
    /// `W` — Whole Sign (`polar_safe`)
    WholeSign,
    /// `B` — Alcabitius (`soft_fail`)
    Alcabitius,
    /// `M` — Morinus (`polar_safe`)
    Morinus,
    /// `T` — Polich-Page (Topocentric) (`soft_fail`)
    PolichPage,
    /// `U` — Krusinski-Pisa-Goelzer (`polar_safe`)
    KrusinskiPisaGoelzer,
    /// `X` — Meridian / Axial Rotation (`polar_safe`)
    Meridian,
    /// `V` — Vehlow Equal (`polar_safe`)
    VehlowEqual,
    /// `L` — Pullen SD (Sinusoidal Delta) (`soft_fail`)
    PullenSinusoidalDelta,
    /// `Q` — Pullen SR (Sinusoidal Ratio) (`polar_safe`)
    PullenSinusoidalRatio,
    /// `I` — Sunshine (Treindl) (`soft_fail`)
    SunshineTreindl,
    /// `i` — Sunshine (Makransky alt) (`hard_fail`)
    SunshineMakransky,
    /// `F` — Carter Poli-Equatorial (`polar_safe`)
    CarterPoliEquatorial,
    /// `N` — Whole Sign (0° Aries) (`polar_safe`)
    WholeSignAries,
    /// `S` — Sripati (`polar_safe`)
    Sripati,
    /// `D` — Equal (from MC) (`polar_safe`)
    EqualFromMc,
    /// `Y` — APC (Apparent Pole) (`soft_fail`)
    Apc,
    /// `H` — Azimuthal / Horizontal (`polar_safe`)
    Azimuthal,
}

impl HouseSystem {
    /// All 23 case-sensitive house systems in `astra` registry order.
    pub const ALL: [Self; 23] = [
        Self::Placidus,
        Self::Koch,
        Self::Porphyry,
        Self::Regiomontanus,
        Self::Campanus,
        Self::EqualFromAsc,
        Self::WholeSign,
        Self::Alcabitius,
        Self::Morinus,
        Self::PolichPage,
        Self::KrusinskiPisaGoelzer,
        Self::Meridian,
        Self::VehlowEqual,
        Self::PullenSinusoidalDelta,
        Self::PullenSinusoidalRatio,
        Self::SunshineTreindl,
        Self::SunshineMakransky,
        Self::CarterPoliEquatorial,
        Self::WholeSignAries,
        Self::Sripati,
        Self::EqualFromMc,
        Self::Apc,
        Self::Azimuthal,
    ];

    /// Parses a single case-sensitive character code into a [`HouseSystem`].
    pub fn from_code(code: char) -> Result<Self, EphemerisError> {
        match code {
            'P' => Ok(Self::Placidus),
            'K' => Ok(Self::Koch),
            'O' => Ok(Self::Porphyry),
            'R' => Ok(Self::Regiomontanus),
            'C' => Ok(Self::Campanus),
            'A' => Ok(Self::EqualFromAsc),
            'W' => Ok(Self::WholeSign),
            'B' => Ok(Self::Alcabitius),
            'M' => Ok(Self::Morinus),
            'T' => Ok(Self::PolichPage),
            'U' => Ok(Self::KrusinskiPisaGoelzer),
            'X' => Ok(Self::Meridian),
            'V' => Ok(Self::VehlowEqual),
            'L' => Ok(Self::PullenSinusoidalDelta),
            'Q' => Ok(Self::PullenSinusoidalRatio),
            'I' => Ok(Self::SunshineTreindl),
            'i' => Ok(Self::SunshineMakransky),
            'F' => Ok(Self::CarterPoliEquatorial),
            'N' => Ok(Self::WholeSignAries),
            'S' => Ok(Self::Sripati),
            'D' => Ok(Self::EqualFromMc),
            'Y' => Ok(Self::Apc),
            'H' => Ok(Self::Azimuthal),
            other => Err(EphemerisError::UnsupportedModel {
                model: format!("HouseSystem({other})"),
            }),
        }
    }

    /// Returns the single case-sensitive character code.
    #[must_use]
    pub const fn code(self) -> char {
        match self {
            Self::Placidus => 'P',
            Self::Koch => 'K',
            Self::Porphyry => 'O',
            Self::Regiomontanus => 'R',
            Self::Campanus => 'C',
            Self::EqualFromAsc => 'A',
            Self::WholeSign => 'W',
            Self::Alcabitius => 'B',
            Self::Morinus => 'M',
            Self::PolichPage => 'T',
            Self::KrusinskiPisaGoelzer => 'U',
            Self::Meridian => 'X',
            Self::VehlowEqual => 'V',
            Self::PullenSinusoidalDelta => 'L',
            Self::PullenSinusoidalRatio => 'Q',
            Self::SunshineTreindl => 'I',
            Self::SunshineMakransky => 'i',
            Self::CarterPoliEquatorial => 'F',
            Self::WholeSignAries => 'N',
            Self::Sripati => 'S',
            Self::EqualFromMc => 'D',
            Self::Apc => 'Y',
            Self::Azimuthal => 'H',
        }
    }

    /// Exact human-readable name from `astra`.
    #[must_use]
    pub const fn astra_name(self) -> &'static str {
        match self {
            Self::Placidus => "Placidus",
            Self::Koch => "Koch",
            Self::Porphyry => "Porphyry",
            Self::Regiomontanus => "Regiomontanus",
            Self::Campanus => "Campanus",
            Self::EqualFromAsc => "Equal (from ASC)",
            Self::WholeSign => "Whole Sign",
            Self::Alcabitius => "Alcabitius",
            Self::Morinus => "Morinus",
            Self::PolichPage => "Polich-Page (Topocentric)",
            Self::KrusinskiPisaGoelzer => "Krusinski-Pisa-Goelzer",
            Self::Meridian => "Meridian / Axial Rotation",
            Self::VehlowEqual => "Vehlow Equal",
            Self::PullenSinusoidalDelta => "Pullen SD (Sinusoidal Delta)",
            Self::PullenSinusoidalRatio => "Pullen SR (Sinusoidal Ratio)",
            Self::SunshineTreindl => "Sunshine (Treindl)",
            Self::SunshineMakransky => "Sunshine (Makransky alt)",
            Self::CarterPoliEquatorial => "Carter Poli-Equatorial",
            Self::WholeSignAries => "Whole Sign (0° Aries)",
            Self::Sripati => "Sripati",
            Self::EqualFromMc => "Equal (from MC)",
            Self::Apc => "APC (Apparent Pole)",
            Self::Azimuthal => "Azimuthal / Horizontal",
        }
    }

    /// Returns the polar behavior classification from `REVIEW_AND_SPEC_ASTRA.md` §7.
    #[must_use]
    pub const fn polar_category(self) -> PolarCategory {
        match self {
            Self::Placidus | Self::Koch | Self::SunshineMakransky => PolarCategory::HardFail,
            Self::Regiomontanus
            | Self::Campanus
            | Self::Alcabitius
            | Self::PolichPage
            | Self::PullenSinusoidalDelta
            | Self::SunshineTreindl
            | Self::Apc => PolarCategory::SoftFail,
            Self::Porphyry
            | Self::EqualFromAsc
            | Self::WholeSign
            | Self::Morinus
            | Self::KrusinskiPisaGoelzer
            | Self::Meridian
            | Self::VehlowEqual
            | Self::PullenSinusoidalRatio
            | Self::CarterPoliEquatorial
            | Self::WholeSignAries
            | Self::Sripati
            | Self::EqualFromMc
            | Self::Azimuthal => PolarCategory::PolarSafe,
        }
    }
}

impl std::fmt::Display for PlanetId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

impl std::str::FromStr for PlanetId {
    type Err = EphemerisError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "sun" => Ok(Self::Sun),
            "moon" => Ok(Self::Moon),
            "mercury" => Ok(Self::Mercury),
            "venus" => Ok(Self::Venus),
            "earth" => Ok(Self::Earth),
            "mars" => Ok(Self::Mars),
            "jupiter" => Ok(Self::Jupiter),
            "saturn" => Ok(Self::Saturn),
            "uranus" => Ok(Self::Uranus),
            "neptune" => Ok(Self::Neptune),
            "pluto" => Ok(Self::Pluto),
            "emb" | "earthmoonbarycenter" | "earth_moon_barycenter" => {
                Ok(Self::EarthMoonBarycenter)
            }
            "ssb" | "solarsystembarycenter" | "solar_system_barycenter" => {
                Ok(Self::SolarSystemBarycenter)
            }
            _ => Err(EphemerisError::UnsupportedBody {
                body: s.to_owned(),
            }),
        }
    }
}

impl std::fmt::Display for LunarPointId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.astra_id())
    }
}

impl std::str::FromStr for LunarPointId {
    type Err = EphemerisError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "mean_node" | "meannode" => Ok(Self::MeanNode),
            "true_node" | "truenode" => Ok(Self::TrueNode),
            "mean_south_node" | "meansouthnode" => Ok(Self::MeanSouthNode),
            "true_south_node" | "truesouthnode" => Ok(Self::TrueSouthNode),
            "mean_lilith" | "meanlilith" => Ok(Self::MeanLilith),
            "true_lilith" | "truelilith" => Ok(Self::TrueLilith),
            "interpolated_lilith" | "interpolatedlilith" | "intp_apog" => {
                Ok(Self::InterpolatedLilith)
            }
            "interpolated_perigee" | "interpolatedperigee" | "intp_perg" => {
                Ok(Self::InterpolatedPerigee)
            }
            _ => Err(EphemerisError::UnsupportedBody {
                body: s.to_owned(),
            }),
        }
    }
}

impl std::fmt::Display for HypotheticalId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.astra_id())
    }
}

impl std::str::FromStr for HypotheticalId {
    type Err = EphemerisError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "cupido" => Ok(Self::Cupido),
            "hades" => Ok(Self::Hades),
            "zeus" => Ok(Self::Zeus),
            "kronos" => Ok(Self::Kronos),
            "apollon" => Ok(Self::Apollon),
            "admetos" => Ok(Self::Admetos),
            "vulkanus" => Ok(Self::Vulkanus),
            "poseidon" => Ok(Self::Poseidon),
            "transpluto" | "isis" => Ok(Self::Transpluto),
            _ => Err(EphemerisError::UnsupportedBody {
                body: s.to_owned(),
            }),
        }
    }
}

impl std::fmt::Display for BodyId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Planet(p) => write!(f, "{p}"),
            Self::MinorPlanet(mpc) => {
                if let Some(spec) = MINOR_BODIES_68.iter().find(|m| m.mpc == *mpc) {
                    f.write_str(spec.name)
                } else {
                    write!(f, "mpc:{}", mpc.0)
                }
            }
            Self::LunarPoint(lp) => write!(f, "{lp}"),
            Self::Hypothetical(h) => write!(f, "{h}"),
            Self::FixedStar(hip) => write!(f, "hip:{}", hip.0),
            Self::EclNut => f.write_str("EclNut"),
        }
    }
}

impl std::str::FromStr for BodyId {
    type Err = EphemerisError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim();
        let lower = trimmed.to_ascii_lowercase();
        if lower == "eclnut" || lower == "ecl_nut" {
            return Ok(Self::EclNut);
        }
        if let Some(mpc_s) = lower.strip_prefix("mpc:") {
            let num = mpc_s
                .trim()
                .parse::<u32>()
                .map_err(|_| EphemerisError::UnsupportedBody {
                    body: trimmed.to_owned(),
                })?;
            return Ok(Self::MinorPlanet(MpcNumber(num)));
        }
        if let Some(hip_s) = lower.strip_prefix("hip:") {
            let num = hip_s
                .trim()
                .parse::<u32>()
                .map_err(|_| EphemerisError::UnsupportedBody {
                    body: trimmed.to_owned(),
                })?;
            return Ok(Self::FixedStar(HipId(num)));
        }
        if let Ok(p) = trimmed.parse::<PlanetId>() {
            return Ok(Self::Planet(p));
        }
        if let Ok(lp) = trimmed.parse::<LunarPointId>() {
            return Ok(Self::LunarPoint(lp));
        }
        if let Ok(h) = trimmed.parse::<HypotheticalId>() {
            return Ok(Self::Hypothetical(h));
        }
        if let Some(spec) = MINOR_BODIES_68
            .iter()
            .find(|m| m.name.eq_ignore_ascii_case(trimmed))
        {
            return Ok(Self::MinorPlanet(spec.mpc));
        }
        if let Some(star) = crate::stars::find_star_by_name_or_bayer(trimmed) {
            return Ok(Self::FixedStar(star.hip_id()));
        }
        Err(EphemerisError::UnsupportedBody {
            body: trimmed.to_owned(),
        })
    }
}

impl std::fmt::Display for AyanamsaModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.code())
    }
}

impl std::str::FromStr for AyanamsaModel {
    type Err = EphemerisError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_code(s.trim())
    }
}

impl std::fmt::Display for HouseSystem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.code())
    }
}

impl std::str::FromStr for HouseSystem {
    type Err = EphemerisError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut chars = s.chars();
        if let (Some(c), None) = (chars.next(), chars.next()) {
            Self::from_code(c)
        } else {
            Err(EphemerisError::UnsupportedModel {
                model: format!("HouseSystem({s})"),
            })
        }
    }
}

impl std::fmt::Display for PolarPolicy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Self::Strict => "strict",
            Self::FallbackToPorphyry => "fallback_to_porphyry",
            Self::FallbackToWholeSign => "fallback_to_whole_sign",
        };
        f.write_str(s)
    }
}

impl std::str::FromStr for PolarPolicy {
    type Err = EphemerisError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "strict" => Ok(Self::Strict),
            "porphyry" | "fallback_to_porphyry" | "o" => Ok(Self::FallbackToPorphyry),
            "whole_sign" | "wholesign" | "fallback_to_whole_sign" | "w" => {
                Ok(Self::FallbackToWholeSign)
            }
            other => Err(EphemerisError::InvalidInput {
                message: format!("Unknown polar policy: {other}"),
            }),
        }
    }
}
