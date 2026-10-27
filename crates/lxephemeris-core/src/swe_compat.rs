//! Generalized Swiss Ephemeris (`pyswisseph` / `swisseph`) migration compatibility layer (`SweAdapter`).
//!
//! Provides a clean-room, instance-scoped adapter (zero process-global state) mapping legacy
//! numerical body codes (`SE_SUN..SE_PLUTO`, `SE_MEAN_NODE..SE_INTP_APOG`, `SE_CHIRON..SE_VESTA`,
//! `SE_CUPIDO..SE_ISIS`, `SE_AST_OFFSET + mpc`, `SE_ECL_NUT`), bitflags (`SEFLG_*`), and sidereal
//! mode codes (`SE_SIDM_*`) directly onto [`Engine`].

use crate::astrometry::{
    AtmosphereParams, CalcResult, DiscMode, HorizontalCoords, RiseTransEvent, RiseTransOptions,
};
use crate::eclipses::{EclipseOptions, LunarEclipseKind, SolarEclipseKind};
use crate::engine::Engine;
use crate::error::EphemerisError;
use crate::events::body_display_name;
use crate::heliacal::HeliacalOptions;
use crate::houses::HouseOptions;
use crate::time::{julday, revjul, CalendarDate, CalendarKind, TimeInput};
use crate::types::{
    AbMode, AyanamsaModel, BodyId, CalcOptions, Equinox, HouseSystem, HypotheticalId, LunarPointId,
    MpcNumber, NutationMode, Observer, Origin, Plane, PlanetId, PolarPolicy, Zodiac,
    MINOR_BODIES_68, SECONDS_PER_DAY,
};
use std::sync::Arc;

// === Legacy numerical body constants ===
pub const SE_ECL_NUT: i32 = -1;
pub const SE_SUN: i32 = 0;
pub const SE_MOON: i32 = 1;
pub const SE_MERCURY: i32 = 2;
pub const SE_VENUS: i32 = 3;
pub const SE_MARS: i32 = 4;
pub const SE_JUPITER: i32 = 5;
pub const SE_SATURN: i32 = 6;
pub const SE_URANUS: i32 = 7;
pub const SE_NEPTUNE: i32 = 8;
pub const SE_PLUTO: i32 = 9;
pub const SE_MEAN_NODE: i32 = 10;
pub const SE_TRUE_NODE: i32 = 11;
pub const SE_MEAN_APOG: i32 = 12;
pub const SE_OSCU_APOG: i32 = 13;
pub const SE_EARTH: i32 = 14;
pub const SE_CHIRON: i32 = 15;
pub const SE_PHOLUS: i32 = 16;
pub const SE_CERES: i32 = 17;
pub const SE_PALLAS: i32 = 18;
pub const SE_JUNO: i32 = 19;
pub const SE_VESTA: i32 = 20;
pub const SE_INTP_APOG: i32 = 21;
pub const SE_INTP_PERG: i32 = 22;
pub const SE_FICT_OFFSET: i32 = 40;
pub const SE_CUPIDO: i32 = 40;
pub const SE_HADES: i32 = 41;
pub const SE_ZEUS: i32 = 42;
pub const SE_KRONOS: i32 = 43;
pub const SE_APOLLON: i32 = 44;
pub const SE_ADMETOS: i32 = 45;
pub const SE_VULKANUS: i32 = 46;
pub const SE_POSEIDON: i32 = 47;
pub const SE_ISIS: i32 = 48;
pub const SE_AST_OFFSET: i32 = 10_000;

// === Legacy calculation bitflags ===
pub const SEFLG_JPLEPH: u32 = 1;
pub const SEFLG_SWIEPH: u32 = 2;
pub const SEFLG_MOSEPH: u32 = 4;
pub const SEFLG_HELCTR: u32 = 8;
pub const SEFLG_TRUEPOS: u32 = 16;
pub const SEFLG_J2000: u32 = 32;
pub const SEFLG_NONUT: u32 = 64;
pub const SEFLG_SPEED3: u32 = 128;
pub const SEFLG_SPEED: u32 = 256;
pub const SEFLG_NOGDEFL: u32 = 512;
pub const SEFLG_NOABERR: u32 = 1024;
pub const SEFLG_EQUATORIAL: u32 = 2048;
pub const SEFLG_XYZ: u32 = 4096;
pub const SEFLG_RADIANS: u32 = 8192;
pub const SEFLG_BARYCTR: u32 = 16384;
pub const SEFLG_TOPOCTR: u32 = 32768;
pub const SEFLG_SIDEREAL: u32 = 65536;
pub const SEFLG_ICRS: u32 = 131072;

// === Legacy sidereal mode constants (all 43 named models + SE_SIDM_USER = 255) ===
pub const SE_SIDM_FAGAN_BRADLEY: u32 = 0;
pub const SE_SIDM_LAHIRI: u32 = 1;
pub const SE_SIDM_DELUCE: u32 = 2;
pub const SE_SIDM_RAMAN: u32 = 3;
pub const SE_SIDM_USHASHASHI: u32 = 4;
pub const SE_SIDM_KRISHNAMURTI: u32 = 5;
pub const SE_SIDM_DJWHAL_KHUL: u32 = 6;
pub const SE_SIDM_YUKTESHWAR: u32 = 7;
pub const SE_SIDM_JN_BHASIN: u32 = 8;
pub const SE_SIDM_BABYL_KUGLER1: u32 = 9;
pub const SE_SIDM_BABYL_KUGLER2: u32 = 10;
pub const SE_SIDM_BABYL_KUGLER3: u32 = 11;
pub const SE_SIDM_BABYL_HUBER: u32 = 12;
pub const SE_SIDM_BABYL_ETPSC: u32 = 13;
pub const SE_SIDM_ALDEBARAN_15TAU: u32 = 14;
pub const SE_SIDM_HIPPARCHOS: u32 = 15;
pub const SE_SIDM_SASSANIAN: u32 = 16;
pub const SE_SIDM_GALCENT_0SAG: u32 = 17;
pub const SE_SIDM_J2000: u32 = 18;
pub const SE_SIDM_J1900: u32 = 19;
pub const SE_SIDM_B1950: u32 = 20;
pub const SE_SIDM_SURYASIDDHANTA: u32 = 21;
pub const SE_SIDM_SURYASIDDHANTA_MSUN: u32 = 22;
pub const SE_SIDM_ARYABHATA: u32 = 23;
pub const SE_SIDM_ARYABHATA_MSUN: u32 = 24;
pub const SE_SIDM_SS_REVATI: u32 = 25;
pub const SE_SIDM_SS_CITRA: u32 = 26;
pub const SE_SIDM_TRUE_CITRA: u32 = 27;
pub const SE_SIDM_TRUE_REVATI: u32 = 28;
pub const SE_SIDM_TRUE_PUSHYA: u32 = 29;
pub const SE_SIDM_GALCENT_RGILBRAND: u32 = 30;
pub const SE_SIDM_GALEQ_IAU1958: u32 = 31;
pub const SE_SIDM_GALEQ_TRUE: u32 = 32;
pub const SE_SIDM_GALEQ_MULA: u32 = 33;
pub const SE_SIDM_GALALIGN_MARDYKS: u32 = 34;
pub const SE_SIDM_TRUE_MULA: u32 = 35;
pub const SE_SIDM_GALCENT_MULA_WILHELM: u32 = 36;
pub const SE_SIDM_ARYABHATA_522: u32 = 37;
pub const SE_SIDM_BABYL_BRITTON: u32 = 38;
pub const SE_SIDM_TRUE_SHERATAN: u32 = 39;
pub const SE_SIDM_TRUE_ALDEBARAN: u32 = 40;
pub const SE_SIDM_GALEQU_FIORENZA: u32 = 41;
pub const SE_SIDM_VALENS_MOON: u32 = 42;
pub const SE_SIDM_LAHIRI_1940: u32 = 43;
pub const SE_SIDM_LAHIRI_VP285: u32 = 44;
pub const SE_SIDM_KRISHNAMURTI_VP291: u32 = 45;
pub const SE_SIDM_KP_SENTHILATHIBAN: u32 = 45;
pub const SE_SIDM_LAHIRI_ICRC: u32 = 46;
pub const SE_SIDM_GALCENT_COCHRANE: u32 = 47;
pub const SE_SIDM_VEDIC_SHEORAN: u32 = 48;
pub const SE_SIDM_USER: u32 = 255;
pub const SE_SIDBIT_ECL_T0: u32 = 256;
pub const SE_SIDBIT_SSY_PLANE: u32 = 512;
pub const SE_SIDBIT_USER_UT: u32 = 1024;
pub const SE_SIDBIT_ECL_DATE: u32 = 2048;
pub const SE_SIDBIT_NO_PREC_OFFSET: u32 = 4096;
pub const SE_SIDBIT_PREC_ORIG: u32 = 8192;

// === Legacy rise/transit and horizontal coordinate constants ===
pub const SE_CALC_RISE: i32 = 1;
pub const SE_CALC_SET: i32 = 2;
pub const SE_CALC_MTRANSIT: i32 = 4;
pub const SE_CALC_ITRANSIT: i32 = 8;
pub const SE_BIT_GEOCTR_NO_ECL_LAT: i32 = 128;
pub const SE_BIT_DISC_CENTER: i32 = 256;
pub const SE_BIT_DISC_BOTTOM: i32 = 8192;
pub const SE_BIT_NO_REFRACTION: i32 = 512;
pub const SE_BIT_CIVIL_TWILIGHT: i32 = 1024;
pub const SE_BIT_NAUTIC_TWILIGHT: i32 = 2048;
pub const SE_BIT_ASTRO_TWILIGHT: i32 = 4096;
pub const SE_BIT_FIXED_DISC_SIZE: i32 = 16384;
pub const SE_BIT_HINDU_RISING: i32 = SE_BIT_DISC_CENTER | SE_BIT_NO_REFRACTION | 128;
pub const SE_ECL2HOR: i32 = 0;
pub const SE_EQU2HOR: i32 = 1;
pub const SE_HOR2ECL: i32 = 0;
pub const SE_HOR2EQU: i32 = 1;

// === Legacy calendar flags ===
pub const SE_JUL_CAL: i32 = 0;
pub const SE_GREG_CAL: i32 = 1;

// === Legacy eclipse & heliacal constants ===
pub const SE_ECL_CENTRAL: i32 = 1;
pub const SE_ECL_NONCENTRAL: i32 = 2;
pub const SE_ECL_TOTAL: i32 = 4;
pub const SE_ECL_ANNULAR: i32 = 8;
pub const SE_ECL_PARTIAL: i32 = 16;
pub const SE_ECL_ANNULAR_TOTAL: i32 = 32;
pub const SE_ECL_HYBRID: i32 = 32;
pub const SE_ECL_PENUMBRAL: i32 = 64;
pub const SE_ECL_ALLTYPES_SOLAR: i32 =
    SE_ECL_CENTRAL | SE_ECL_NONCENTRAL | SE_ECL_TOTAL | SE_ECL_ANNULAR | SE_ECL_PARTIAL | SE_ECL_ANNULAR_TOTAL;
pub const SE_ECL_ALLTYPES_LUNAR: i32 = SE_ECL_TOTAL | SE_ECL_PARTIAL | SE_ECL_PENUMBRAL;
pub const SE_ECL_VISIBLE: i32 = 128;
pub const SE_ECL_MAX_VISIBLE: i32 = 256;
pub const SE_ECL_1ST_VISIBLE: i32 = 512;
pub const SE_ECL_2ND_VISIBLE: i32 = 1024;
pub const SE_ECL_3RD_VISIBLE: i32 = 2048;
pub const SE_ECL_4TH_VISIBLE: i32 = 4096;
pub const SE_ECL_OCC_BEG_DAYLIGHT: i32 = 8192;
pub const SE_ECL_OCC_END_DAYLIGHT: i32 = 16384;
pub const SE_ECL_ONE_TRY: i32 = 32768;

pub const SE_HELIACAL_RISING: i32 = 1;
pub const SE_HELIACAL_SETTING: i32 = 2;
pub const SE_EVENING_FIRST: i32 = 3;
pub const SE_MORNING_LAST: i32 = 4;
pub const SE_HELFLAG_HIGH_PRECISION: i32 = 256;

/// Maps a legacy integer body code `ipl` to a typed [`BodyId`].
pub fn body_id_from_swe_ipl(ipl: i32) -> Result<BodyId, EphemerisError> {
    match ipl {
        SE_ECL_NUT => Ok(BodyId::EclNut),
        SE_SUN => Ok(BodyId::Planet(PlanetId::Sun)),
        SE_MOON => Ok(BodyId::Planet(PlanetId::Moon)),
        SE_MERCURY => Ok(BodyId::Planet(PlanetId::Mercury)),
        SE_VENUS => Ok(BodyId::Planet(PlanetId::Venus)),
        SE_MARS => Ok(BodyId::Planet(PlanetId::Mars)),
        SE_JUPITER => Ok(BodyId::Planet(PlanetId::Jupiter)),
        SE_SATURN => Ok(BodyId::Planet(PlanetId::Saturn)),
        SE_URANUS => Ok(BodyId::Planet(PlanetId::Uranus)),
        SE_NEPTUNE => Ok(BodyId::Planet(PlanetId::Neptune)),
        SE_PLUTO => Ok(BodyId::Planet(PlanetId::Pluto)),
        SE_MEAN_NODE => Ok(BodyId::LunarPoint(LunarPointId::MeanNode)),
        SE_TRUE_NODE => Ok(BodyId::LunarPoint(LunarPointId::TrueNode)),
        SE_MEAN_APOG => Ok(BodyId::LunarPoint(LunarPointId::MeanLilith)),
        SE_OSCU_APOG => Ok(BodyId::LunarPoint(LunarPointId::TrueLilith)),
        SE_EARTH => Ok(BodyId::Planet(PlanetId::Earth)),
        SE_CHIRON => Ok(BodyId::MinorPlanet(MpcNumber(2060))),
        SE_PHOLUS => Ok(BodyId::MinorPlanet(MpcNumber(5145))),
        SE_CERES => Ok(BodyId::MinorPlanet(MpcNumber(1))),
        SE_PALLAS => Ok(BodyId::MinorPlanet(MpcNumber(2))),
        SE_JUNO => Ok(BodyId::MinorPlanet(MpcNumber(3))),
        SE_VESTA => Ok(BodyId::MinorPlanet(MpcNumber(4))),
        SE_INTP_APOG => Ok(BodyId::LunarPoint(LunarPointId::InterpolatedLilith)),
        SE_INTP_PERG => Ok(BodyId::LunarPoint(LunarPointId::InterpolatedPerigee)),
        SE_CUPIDO => Ok(BodyId::Hypothetical(HypotheticalId::Cupido)),
        SE_HADES => Ok(BodyId::Hypothetical(HypotheticalId::Hades)),
        SE_ZEUS => Ok(BodyId::Hypothetical(HypotheticalId::Zeus)),
        SE_KRONOS => Ok(BodyId::Hypothetical(HypotheticalId::Kronos)),
        SE_APOLLON => Ok(BodyId::Hypothetical(HypotheticalId::Apollon)),
        SE_ADMETOS => Ok(BodyId::Hypothetical(HypotheticalId::Admetos)),
        SE_VULKANUS => Ok(BodyId::Hypothetical(HypotheticalId::Vulkanus)),
        SE_POSEIDON => Ok(BodyId::Hypothetical(HypotheticalId::Poseidon)),
        SE_ISIS => Ok(BodyId::Hypothetical(HypotheticalId::Transpluto)),
        n if n > SE_AST_OFFSET => Ok(BodyId::MinorPlanet(MpcNumber((n - SE_AST_OFFSET) as u32))),
        other => Err(EphemerisError::UnsupportedBody {
            body: format!("SE_IPL({other})"),
        }),
    }
}

/// Converts legacy `SEFLG_*` bitflags into [`CalcOptions`].
pub fn calc_options_from_swe_flags(
    flags: u32,
    topo_observer: Option<Observer>,
    sid_mode: AyanamsaModel,
) -> Result<CalcOptions, EphemerisError> {
    let origin = if flags & SEFLG_TOPOCTR != 0 {
        let obs = topo_observer.ok_or_else(|| EphemerisError::InvalidInput {
            message: "SEFLG_TOPOCTR set without prior set_topo(lon, lat, alt)".to_owned(),
        })?;
        Origin::Topocentric(obs)
    } else if flags & SEFLG_HELCTR != 0 {
        Origin::Heliocentric
    } else if flags & SEFLG_BARYCTR != 0 {
        Origin::Barycentric
    } else {
        Origin::Geocentric
    };

    let plane = if flags & SEFLG_EQUATORIAL != 0 {
        Plane::Equatorial
    } else {
        Plane::Ecliptic
    };

    let equinox = if flags & SEFLG_J2000 != 0 {
        Equinox::J2000
    } else {
        Equinox::OfDate
    };

    let nutation = if (flags & SEFLG_NONUT != 0) || (flags & SEFLG_J2000 != 0) {
        NutationMode::Mean
    } else {
        NutationMode::True
    };

    let aberration = if flags & SEFLG_TRUEPOS != 0 {
        AbMode::Geometric
    } else if (flags & SEFLG_NOABERR != 0) && (flags & SEFLG_NOGDEFL != 0) {
        AbMode::Astrometric
    } else if flags & SEFLG_NOABERR != 0 {
        AbMode::NoAberration
    } else if flags & SEFLG_NOGDEFL != 0 {
        AbMode::NoDeflection
    } else {
        AbMode::Apparent
    };

    let zodiac = if (flags & SEFLG_SIDEREAL != 0) && plane == Plane::Ecliptic {
        Zodiac::Sidereal(sid_mode)
    } else {
        Zodiac::Tropical
    };

    let include_speed = (flags & (SEFLG_SPEED | 128)) != 0;

    Ok(CalcOptions {
        origin,
        plane,
        equinox,
        nutation,
        aberration,
        zodiac,
        include_speed,
    })
}

/// Formats a [`CalcResult`] into the 6-element `[v0, v1, v2, v3, v4, v5]` array expected by `swe.calc_ut`.
#[must_use]
pub fn format_swe_state6(res: &CalcResult, flags: u32) -> [f64; 6] {
    if res.body == BodyId::EclNut {
        if flags & SEFLG_RADIANS != 0 {
            return [
                res.lon_deg.to_radians(),
                res.lat_deg.to_radians(),
                res.dist_au.to_radians(),
                res.lon_speed_deg_day.to_radians(),
                0.0,
                0.0,
            ];
        }
        return [
            res.lon_deg,
            res.lat_deg,
            res.dist_au,
            res.lon_speed_deg_day,
            0.0,
            0.0,
        ];
    }

    if flags & SEFLG_XYZ != 0 {
        [
            res.cartesian_au[0],
            res.cartesian_au[1],
            res.cartesian_au[2],
            res.cartesian_vel_au_day[0],
            res.cartesian_vel_au_day[1],
            res.cartesian_vel_au_day[2],
        ]
    } else if flags & SEFLG_RADIANS != 0 {
        [
            res.lon_deg.to_radians(),
            res.lat_deg.to_radians(),
            res.dist_au,
            res.lon_speed_deg_day.to_radians(),
            res.lat_speed_deg_day.to_radians(),
            res.dist_speed_au_day,
        ]
    } else {
        [
            res.lon_deg,
            res.lat_deg,
            res.dist_au,
            res.lon_speed_deg_day,
            res.lat_speed_deg_day,
            res.dist_speed_au_day,
        ]
    }
}

/// Instance-scoped adapter for drop-in migration from `pyswisseph` / `swisseph` without global state.
#[derive(Clone)]
pub struct SweAdapter {
    engine: Arc<Engine>,
    topo_observer: Option<Observer>,
    sid_mode: AyanamsaModel,
    polar_policy: PolarPolicy,
}

impl SweAdapter {
    /// Creates a new [`SweAdapter`] wrapping `engine` (defaults to `Lahiri` sidereal mode and `FallbackToPorphyry` polar policy, matching Swiss Ephemeris).
    #[must_use]
    pub fn new(engine: Arc<Engine>) -> Self {
        Self {
            engine,
            topo_observer: None,
            sid_mode: AyanamsaModel::Lahiri,
            polar_policy: PolarPolicy::FallbackToPorphyry,
        }
    }

    /// Returns a reference to the underlying [`Engine`].
    #[must_use]
    pub fn engine(&self) -> &Engine {
        &self.engine
    }

    /// Sets the instance topocentric observer (`lon_deg` East positive, `lat_deg` North positive, `alt_m` meters),
    /// matching `swe.set_topo(geolon, geolat, geoalt)`.
    pub fn set_topo(
        &mut self,
        lon_deg: f64,
        lat_deg: f64,
        alt_m: f64,
    ) -> Result<(), EphemerisError> {
        self.topo_observer = Some(Observer::new(lat_deg, lon_deg, alt_m)?);
        Ok(())
    }

    /// Sets the instance sidereal mode from a `SE_SIDM_*` integer constant.
    pub fn set_sid_mode(&mut self, sid_mode: u32) -> Result<(), EphemerisError> {
        self.sid_mode = AyanamsaModel::from_legacy_id(sid_mode & 0xFF)?;
        Ok(())
    }

    /// Sets the instance sidereal mode with optional custom epoch `t0` (`JD_TT`) and initial ayanamsa
    /// `ayan_t0` (degrees) when `sid_mode == SE_SIDM_USER` (`255`), matching `swe_set_sid_mode(sid_mode, t0, ayan_t0)`.
    pub fn set_sid_mode_custom(
        &mut self,
        sid_mode: u32,
        t0: f64,
        ayan_t0: f64,
    ) -> Result<(), EphemerisError> {
        let base_mode = sid_mode & 0xFF;
        if base_mode == SE_SIDM_USER {
            self.sid_mode = AyanamsaModel::user(t0, ayan_t0)?;
        } else {
            self.sid_mode = AyanamsaModel::from_legacy_id(base_mode)?;
        }
        Ok(())
    }

    /// Sets the instance sidereal mode from an `astra` ayanamsa string code (e.g. `"lahiri"`).
    pub fn set_sid_mode_by_name(&mut self, code: &str) -> Result<(), EphemerisError> {
        self.sid_mode = AyanamsaModel::from_code(code)?;
        Ok(())
    }

    /// Sets the instance [`PolarPolicy`] for `houses_ex` / `houses_armc`.
    pub fn set_polar_policy(&mut self, policy: PolarPolicy) {
        self.polar_policy = policy;
    }

    /// Returns the currently configured [`AyanamsaModel`].
    #[must_use]
    pub const fn sid_mode(&self) -> AyanamsaModel {
        self.sid_mode
    }

    /// Returns the currently configured topocentric [`Observer`], if any.
    #[must_use]
    pub const fn topo_observer(&self) -> Option<Observer> {
        self.topo_observer
    }

    const fn swe_star_anchor(model: AyanamsaModel) -> Option<(crate::types::HipId, f64)> {
        match model {
            AyanamsaModel::TrueCitra => Some((crate::types::HipId(65474), 180.0)),
            AyanamsaModel::TrueRevati => Some((crate::types::HipId(5737), 359.0 + 50.0 / 60.0)),
            AyanamsaModel::TruePushya => Some((crate::types::HipId(42911), 106.0)),
            _ => None,
        }
    }

    fn swe_ayanamsa_val(
        &self,
        t_in: TimeInput,
        flags: u32,
    ) -> Result<crate::ayanamsa::AyanamsaValue, EphemerisError> {
        let nut = if (flags & SEFLG_NONUT) != 0 {
            NutationMode::Mean
        } else {
            NutationMode::True
        };
        let mut val = self.engine.ayanamsa(t_in, self.sid_mode, nut)?;
        if let Some((hip, target_lon)) = Self::swe_star_anchor(self.sid_mode) {
            if (flags & (SEFLG_NOABERR | SEFLG_TRUEPOS)) == 0 {
                let star_opts = CalcOptions {
                    equinox: Equinox::OfDate,
                    zodiac: Zodiac::Tropical,
                    aberration: AbMode::Apparent,
                    nutation: NutationMode::Mean,
                    origin: Origin::Geocentric,
                    plane: Plane::Ecliptic,
                    include_speed: false,
                };
                let body = BodyId::FixedStar(hip);
                let use_moshier = (flags & SEFLG_MOSEPH != 0) || self.engine.is_moshier_mode();
                let star_res = if use_moshier {
                    self.engine.calc_moshier(body, t_in, star_opts)?
                } else {
                    match self.engine.calc(body, t_in, star_opts) {
                        Ok(r) => r,
                        Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => self
                            .engine
                            .calc_moshier(body, t_in, star_opts)
                            .map_err(|_| orig_err)?,
                        Err(e) => return Err(e),
                    }
                };
                let mean_deg = (star_res.lon_deg - target_lon).rem_euclid(360.0);
                let dpsi = (val.true_deg - val.mean_deg + 180.0).rem_euclid(360.0) - 180.0;
                let true_deg = (mean_deg + dpsi).rem_euclid(360.0);
                val.mean_deg = mean_deg;
                val.true_deg = true_deg;
                val.effective_deg = if nut == NutationMode::Mean {
                    mean_deg
                } else {
                    true_deg
                };
            }
        }
        Ok(val)
    }

    fn apply_swe_star_sidereal_correction(
        &self,
        t_in: TimeInput,
        flags: u32,
        res: &mut CalcResult,
    ) -> Result<(), EphemerisError> {
        if res.body != BodyId::EclNut
            && (flags & SEFLG_SIDEREAL) != 0
            && (flags & SEFLG_EQUATORIAL) == 0
            && (flags & (SEFLG_NOABERR | SEFLG_TRUEPOS)) == 0
            && Self::swe_star_anchor(self.sid_mode).is_some()
        {
            let nut = if (flags & SEFLG_NONUT != 0) || (flags & SEFLG_J2000 != 0) {
                NutationMode::Mean
            } else {
                NutationMode::True
            };
            let base = self.engine.ayanamsa(t_in, self.sid_mode, nut)?;
            let swe_val = self.swe_ayanamsa_val(t_in, flags)?;
            let delta_deg =
                (swe_val.effective_deg - base.effective_deg + 180.0).rem_euclid(360.0) - 180.0;
            res.lon_deg = (res.lon_deg - delta_deg).rem_euclid(360.0);
            if (flags & SEFLG_XYZ) != 0 {
                let (s, c) = (-delta_deg.to_radians()).sin_cos();
                let (x, y) = (res.cartesian_au[0], res.cartesian_au[1]);
                res.cartesian_au[0] = x * c - y * s;
                res.cartesian_au[1] = x * s + y * c;
                let (vx, vy) = (res.cartesian_vel_au_day[0], res.cartesian_vel_au_day[1]);
                res.cartesian_vel_au_day[0] = vx * c - vy * s;
                res.cartesian_vel_au_day[1] = vx * s + vy * c;
            }
        }
        Ok(())
    }

    fn calc_time_input_inner(
        &self,
        t_in: TimeInput,
        epoch_jd: f64,
        ipl: i32,
        flags: u32,
    ) -> Result<([f64; 6], i32), EphemerisError> {
        let body = body_id_from_swe_ipl(ipl)?;
        let opts = calc_options_from_swe_flags(flags, self.topo_observer, self.sid_mode)?;
        let use_moshier = (flags & SEFLG_MOSEPH != 0) || self.engine.is_moshier_mode();
        if ipl == SE_EARTH && opts.origin == Origin::Geocentric {
            let retflag = if use_moshier {
                ((flags & !(SEFLG_JPLEPH | SEFLG_SWIEPH)) | SEFLG_MOSEPH) as i32
            } else {
                (flags | SEFLG_SWIEPH) as i32
            };
            return Ok(([0.0; 6], retflag));
        }
        if matches!(body, BodyId::LunarPoint(_))
            && matches!(opts.origin, Origin::Heliocentric | Origin::Barycentric)
        {
            let retflag = if use_moshier {
                ((flags & !(SEFLG_JPLEPH | SEFLG_SWIEPH)) | SEFLG_MOSEPH) as i32
            } else {
                (flags | SEFLG_SWIEPH) as i32
            };
            return Ok(([0.0; 6], retflag));
        }
        let disallow_moshier_fallback = matches!(body, BodyId::MinorPlanet(_))
            || ipl == 15
            || (ipl >= 17 && ipl != 21 && ipl != 22);
        let (mut res, used_moshier) = if use_moshier {
            (self.engine.calc_moshier(body, t_in, opts)?, true)
        } else {
            match self.engine.calc(body, t_in, opts) {
                Ok(r) => (r, false),
                Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                    if disallow_moshier_fallback {
                        return Err(orig_err);
                    }
                    match self.engine.calc_moshier(body, t_in, opts) {
                        Ok(r) => (r, true),
                        Err(_) => return Err(orig_err),
                    }
                }
                Err(e) => return Err(e),
            }
        };
        if !use_moshier && used_moshier && disallow_moshier_fallback {
            return Err(EphemerisError::OutOfCoverage {
                body: format!("ipl={ipl}"),
                epoch_jd_tdb: epoch_jd,
                start_jd_tdb: 2_396_758.5,
                end_jd_tdb: 2_488_069.5,
            });
        }
        self.apply_swe_star_sidereal_correction(t_in, flags, &mut res)?;
        let retflag = if used_moshier {
            ((flags & !(SEFLG_JPLEPH | SEFLG_SWIEPH)) | SEFLG_MOSEPH) as i32
        } else {
            (flags | SEFLG_SWIEPH) as i32
        };
        Ok((format_swe_state6(&res, flags), retflag))
    }

    fn mean_apog_ecl_time_input(
        &self,
        t_in: TimeInput,
        epoch_jd: f64,
        nonut: bool,
    ) -> Result<(f64, f64, f64, f64, f64), EphemerisError> {
        let bfl = SEFLG_SWIEPH | SEFLG_SPEED | if nonut { SEFLG_NONUT } else { 0 };
        let (lx_apog, _) = self.calc_time_input_inner(t_in, epoch_jd, SE_MEAN_APOG, bfl)?;
        let (node, _) = self.calc_time_input_inner(t_in, epoch_jd, SE_MEAN_NODE, bfl)?;
        let inc = 5.145_396_4_f64.to_radians();
        let u = (lx_apog[0] - node[0]).to_radians();
        let du = (lx_apog[3] - node[3]).to_radians();
        let (sin_inc, cos_inc) = (inc.sin(), inc.cos());
        let (sin_u, cos_u) = (u.sin(), u.cos());
        let lat_rad = (sin_inc * sin_u).clamp(-1.0, 1.0).asin();
        let lat = lat_rad.to_degrees();
        let lon = (node[0] + (cos_inc * sin_u).atan2(cos_u).to_degrees()).rem_euclid(360.0);
        let dlon = node[3] + (du * cos_inc / (1.0 - (sin_inc * sin_u).powi(2))).to_degrees();
        let dlat = (du * sin_inc * cos_u / lat_rad.cos()).to_degrees();
        Ok((lon, lat, 0.002_710_625_131_722_546, dlon, dlat))
    }

    fn calc_with_time_input(
        &self,
        t_in: TimeInput,
        epoch_jd: f64,
        is_ut: bool,
        ipl: i32,
        flags: u32,
    ) -> Result<([f64; 6], i32), EphemerisError> {
        let is_hel_or_bary =
            (flags & SEFLG_TOPOCTR) == 0 && (flags & (SEFLG_HELCTR | SEFLG_BARYCTR)) != 0;
        if ipl == SE_MEAN_APOG && !is_hel_or_bary {
            let (_, retflag) = self.calc_time_input_inner(t_in, epoch_jd, ipl, flags)?;
            let nonut = (flags & (SEFLG_NONUT | SEFLG_J2000)) != 0;
            let (mut lon, mut lat, mut dist, mut dlon, mut dlat) =
                self.mean_apog_ecl_time_input(t_in, epoch_jd, nonut)?;
            let (j2000_dlon_off, j2000_dlat_off) = if (flags & SEFLG_J2000) != 0 {
                let (lx_j, _) = self.calc_time_input_inner(
                    t_in,
                    epoch_jd,
                    SE_MEAN_APOG,
                    SEFLG_J2000 | SEFLG_SWIEPH | SEFLG_SPEED,
                )?;
                let (lx_d, _) = self.calc_time_input_inner(
                    t_in,
                    epoch_jd,
                    SE_MEAN_APOG,
                    SEFLG_NONUT | SEFLG_SWIEPH | SEFLG_SPEED,
                )?;
                let dl = (lx_j[0] - lx_d[0] + 180.0).rem_euclid(360.0) - 180.0;
                let db = lx_j[1] - lx_d[1];
                lon = (lon + dl).rem_euclid(360.0);
                lat += db;
                dlon += lx_j[3] - lx_d[3];
                dlat += lx_j[4] - lx_d[4];
                (dl, db)
            } else {
                (0.0, 0.0)
            };
            if (flags & (SEFLG_SPEED | 128)) == 0 {
                dlon = 0.0;
                dlat = 0.0;
            }
            if (flags & SEFLG_EQUATORIAL) != 0 {
                let eps = if (flags & SEFLG_J2000) != 0 {
                    23.439_291_111_111_11
                } else {
                    let (ecl_nut, _) = self.calc_time_input_inner(t_in, epoch_jd, SE_ECL_NUT, 0)?;
                    if (flags & SEFLG_NONUT) != 0 {
                        ecl_nut[1]
                    } else {
                        ecl_nut[0]
                    }
                };
                let (ra, dec, r_dist) = Self::cotrans(lon, lat, dist, -eps);
                let (dra, ddec) = if (flags & (SEFLG_SPEED | 128)) != 0 {
                    let dt = 1e-4;
                    let t_p = if is_ut {
                        TimeInput::from_ut1_jd(epoch_jd + dt)
                    } else {
                        TimeInput::from_tt_jd(epoch_jd + dt)
                    };
                    let t_m = if is_ut {
                        TimeInput::from_ut1_jd(epoch_jd - dt)
                    } else {
                        TimeInput::from_tt_jd(epoch_jd - dt)
                    };
                    let (lon_p, lat_p, d_p, _, _) =
                        self.mean_apog_ecl_time_input(t_p, epoch_jd + dt, nonut)?;
                    let (lon_m, lat_m, d_m, _, _) =
                        self.mean_apog_ecl_time_input(t_m, epoch_jd - dt, nonut)?;
                    let (ra_p, dec_p, _) = Self::cotrans(
                        (lon_p + j2000_dlon_off).rem_euclid(360.0),
                        lat_p + j2000_dlat_off,
                        d_p,
                        -eps,
                    );
                    let (ra_m, dec_m, _) = Self::cotrans(
                        (lon_m + j2000_dlon_off).rem_euclid(360.0),
                        lat_m + j2000_dlat_off,
                        d_m,
                        -eps,
                    );
                    (
                        (((ra_p - ra_m + 180.0).rem_euclid(360.0)) - 180.0) / (2.0 * dt),
                        (dec_p - dec_m) / (2.0 * dt),
                    )
                } else {
                    (0.0, 0.0)
                };
                lon = ra;
                lat = dec;
                dist = r_dist;
                dlon = dra;
                dlat = ddec;
            } else if (flags & SEFLG_SIDEREAL) != 0 {
                let aya_res = self.swe_ayanamsa_val(t_in, flags)?;
                lon = (lon - aya_res.effective_deg).rem_euclid(360.0);
                if (flags & (SEFLG_SPEED | 128)) != 0 {
                    dlon -= aya_res.rate_deg_day;
                }
            }
            if (flags & SEFLG_XYZ) != 0 {
                let (lr, br) = (lon.to_radians(), lat.to_radians());
                let (vlr, vbr) = (dlon.to_radians(), dlat.to_radians());
                let (cl, sl) = (lr.cos(), lr.sin());
                let (cb, sb) = (br.cos(), br.sin());
                let x = dist * cb * cl;
                let y = dist * cb * sl;
                let z = dist * sb;
                let vx = -dist * sb * cl * vbr - dist * cb * sl * vlr;
                let vy = -dist * sb * sl * vbr + dist * cb * cl * vlr;
                let vz = dist * cb * vbr;
                return Ok(([x, y, z, vx, vy, vz], retflag));
            }
            if (flags & SEFLG_RADIANS) != 0 {
                return Ok((
                    [
                        lon.to_radians(),
                        lat.to_radians(),
                        dist,
                        dlon.to_radians(),
                        dlat.to_radians(),
                        0.0,
                    ],
                    retflag,
                ));
            }
            return Ok(([lon, lat, dist, dlon, dlat, 0.0], retflag));
        }
        self.calc_time_input_inner(t_in, epoch_jd, ipl, flags)
    }

    /// Evaluates `ipl` at `jd_ut` (`UT1`), returning `([v0, v1, v2, v3, v4, v5], retflag)`.
    pub fn calc_ut(
        &self,
        jd_ut: f64,
        ipl: i32,
        flags: u32,
    ) -> Result<([f64; 6], i32), EphemerisError> {
        self.calc_with_time_input(TimeInput::from_ut1_jd(jd_ut), jd_ut, true, ipl, flags)
    }

    /// Evaluates `ipl` at `jd_tt` (`TT` / Ephemeris Time), returning `([v0, v1, v2, v3, v4, v5], retflag)`.
    pub fn calc(
        &self,
        jd_tt: f64,
        ipl: i32,
        flags: u32,
    ) -> Result<([f64; 6], i32), EphemerisError> {
        self.calc_with_time_input(TimeInput::from_tt_jd(jd_tt), jd_tt, false, ipl, flags)
    }

    /// Evaluates a fixed star by name or Bayer designation at `jd_ut` (`UT1`), returning
    /// `([v0, v1, v2, v3, v4, v5], star_name_bayer, retflag)`.
    pub fn fixstar2_ut(
        &self,
        star: &str,
        jd_ut: f64,
        flags: u32,
    ) -> Result<([f64; 6], String, i32), EphemerisError> {
        let meta = self.engine.star_info(star)?;
        let opts = calc_options_from_swe_flags(flags, self.topo_observer, self.sid_mode)?;
        let body = BodyId::FixedStar(crate::types::HipId(meta.hip));
        let t_in = TimeInput::from_ut1_jd(jd_ut);
        let use_moshier = (flags & SEFLG_MOSEPH != 0) || self.engine.is_moshier_mode();
        let (mut res, used_moshier) = if use_moshier {
            (self.engine.calc_moshier(body, t_in, opts)?, true)
        } else {
            match self.engine.calc(body, t_in, opts) {
                Ok(r) => (r, false),
                Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                    match self.engine.calc_moshier(body, t_in, opts) {
                        Ok(r) => (r, true),
                        Err(_) => return Err(orig_err),
                    }
                }
                Err(e) => return Err(e),
            }
        };
        self.apply_swe_star_sidereal_correction(t_in, flags, &mut res)?;
        let retflag = if used_moshier {
            ((flags & !(SEFLG_JPLEPH | SEFLG_SWIEPH)) | SEFLG_MOSEPH) as i32
        } else {
            (flags | SEFLG_SWIEPH) as i32
        };
        Ok((format_swe_state6(&res, flags), meta.canonical_name, retflag))
    }

    /// Returns the catalog visual magnitude `vmag` of a fixed star (`swe.fixstar2_mag(star)`).
    pub fn fixstar2_mag(&self, star: &str) -> Result<f64, EphemerisError> {
        Ok(self.engine.star_info(star)?.vmag)
    }

    /// Evaluates house cusps (`[f64; 12]`) and angles (`[f64; 8]`) at `jd_ut` (`UT1`), matching
    /// `swe.houses_ex(jd_ut, lat, lon, hsys, flags)`.
    pub fn houses_ex(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
        hsys: char,
        flags: u32,
    ) -> Result<([f64; 12], [f64; 8]), EphemerisError> {
        let hsys = if hsys == 'E' { 'A' } else { hsys };
        let alt_m = self.topo_observer.map_or(0.0, |o| o.altitude_m);
        let obs = Observer::new(lat_deg, lon_deg, alt_m)?;
        let system = HouseSystem::from_code(hsys)?;
        let zodiac = if flags & SEFLG_SIDEREAL != 0 {
            Zodiac::Sidereal(self.sid_mode)
        } else {
            Zodiac::Tropical
        };
        let nutation = if flags & SEFLG_NONUT != 0 {
            NutationMode::Mean
        } else {
            NutationMode::True
        };
        let opts = HouseOptions {
            system,
            polar_policy: self.polar_policy,
            zodiac,
            nutation,
            include_speed: false,
        };
        let t_in = TimeInput::from_ut1_jd(jd_ut);
        let res = self.engine.houses(t_in, obs, opts)?;
        let mut cusps = res.cusps_deg;
        let mut ascmc = [
            res.angles.asc_deg,
            res.angles.mc_deg,
            res.angles.armc_deg,
            res.angles.vertex_deg,
            res.angles.equatorial_asc_deg,
            res.angles.co_asc_koch_deg,
            res.angles.co_asc_munkasey_deg,
            res.angles.polar_asc_deg,
        ];
        if (flags & SEFLG_SIDEREAL) != 0 && Self::swe_star_anchor(self.sid_mode).is_some() {
            let base = self.engine.ayanamsa(t_in, self.sid_mode, nutation)?;
            let swe_val = self.swe_ayanamsa_val(t_in, flags)?;
            let delta =
                (swe_val.effective_deg - base.effective_deg + 180.0).rem_euclid(360.0) - 180.0;
            for idx in [0usize, 1, 3, 4, 5, 6, 7] {
                if let Some(a) = ascmc.get_mut(idx) {
                    *a = (*a - delta).rem_euclid(360.0);
                }
            }
            if system == HouseSystem::WholeSign {
                let sign_base = (ascmc[0] / 30.0).floor() * 30.0;
                for (i, c) in cusps.iter_mut().enumerate() {
                    *c = (sign_base + (i as f64) * 30.0).rem_euclid(360.0);
                }
            } else if system != HouseSystem::WholeSignAries {
                for c in &mut cusps {
                    *c = (*c - delta).rem_euclid(360.0);
                }
            }
        }
        if (flags & SEFLG_RADIANS) != 0 {
            for c in &mut cusps {
                *c = c.to_radians();
            }
            for a in &mut ascmc {
                *a = a.to_radians();
            }
        }
        Ok((cusps, ascmc))
    }

    /// Evaluates house cusps (`[f64; 12]`), angles (`[f64; 8]`), cusp speeds (`[f64; 12]`), and
    /// angle speeds (`[f64; 8]`) at `jd_ut` (`UT1`), matching `swe.houses_ex2(jd_ut, lat, lon, hsys, flags)`.
    #[allow(clippy::type_complexity)]
    pub fn houses_ex2(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
        hsys: char,
        flags: u32,
    ) -> Result<([f64; 12], [f64; 8], [f64; 12], [f64; 8]), EphemerisError> {
        let hsys = if hsys == 'E' { 'A' } else { hsys };
        let alt_m = self.topo_observer.map_or(0.0, |o| o.altitude_m);
        let obs = Observer::new(lat_deg, lon_deg, alt_m)?;
        let system = HouseSystem::from_code(hsys)?;
        let zodiac = if flags & SEFLG_SIDEREAL != 0 {
            Zodiac::Sidereal(self.sid_mode)
        } else {
            Zodiac::Tropical
        };
        let nutation = if flags & SEFLG_NONUT != 0 {
            NutationMode::Mean
        } else {
            NutationMode::True
        };
        let opts = HouseOptions {
            system,
            polar_policy: self.polar_policy,
            zodiac,
            nutation,
            include_speed: true,
        };
        let t_in = TimeInput::from_ut1_jd(jd_ut);
        let res = self.engine.houses(t_in, obs, opts)?;
        let mut cusps = res.cusps_deg;
        let mut ascmc = [
            res.angles.asc_deg,
            res.angles.mc_deg,
            res.angles.armc_deg,
            res.angles.vertex_deg,
            res.angles.equatorial_asc_deg,
            res.angles.co_asc_koch_deg,
            res.angles.co_asc_munkasey_deg,
            res.angles.polar_asc_deg,
        ];
        if (flags & SEFLG_SIDEREAL) != 0 && Self::swe_star_anchor(self.sid_mode).is_some() {
            let base = self.engine.ayanamsa(t_in, self.sid_mode, nutation)?;
            let swe_val = self.swe_ayanamsa_val(t_in, flags)?;
            let delta =
                (swe_val.effective_deg - base.effective_deg + 180.0).rem_euclid(360.0) - 180.0;
            for idx in [0usize, 1, 3, 4, 5, 6, 7] {
                if let Some(a) = ascmc.get_mut(idx) {
                    *a = (*a - delta).rem_euclid(360.0);
                }
            }
            if system == HouseSystem::WholeSign {
                let sign_base = (ascmc[0] / 30.0).floor() * 30.0;
                for (i, c) in cusps.iter_mut().enumerate() {
                    *c = (sign_base + (i as f64) * 30.0).rem_euclid(360.0);
                }
            } else if system != HouseSystem::WholeSignAries {
                for c in &mut cusps {
                    *c = (*c - delta).rem_euclid(360.0);
                }
            }
        }
        let mut cusp_speed = res.cusp_speeds_deg_day;
        let ascmc_speed = [
            res.angle_speeds.asc_deg,
            res.angle_speeds.mc_deg,
            res.angle_speeds.armc_deg,
            res.angle_speeds.vertex_deg,
            res.angle_speeds.equatorial_asc_deg,
            res.angle_speeds.co_asc_koch_deg,
            res.angle_speeds.co_asc_munkasey_deg,
            res.angle_speeds.polar_asc_deg,
        ];
        if system == HouseSystem::WholeSign && (flags & SEFLG_SIDEREAL) != 0 {
            cusp_speed = [ascmc_speed[0]; 12];
        } else if matches!(
            system,
            HouseSystem::WholeSign
                | HouseSystem::WholeSignAries
                | HouseSystem::KrusinskiPisaGoelzer
        ) {
            cusp_speed = [0.0; 12];
            cusp_speed[0] = ascmc_speed[0];
            cusp_speed[3] = ascmc_speed[1];
            cusp_speed[6] = ascmc_speed[0];
            cusp_speed[9] = ascmc_speed[1];
        }
        if (flags & SEFLG_RADIANS) != 0 {
            for c in &mut cusps {
                *c = c.to_radians();
            }
            for a in &mut ascmc {
                *a = a.to_radians();
            }
        }
        Ok((cusps, ascmc, cusp_speed, ascmc_speed))
    }

    /// Evaluates house cusps (`[f64; 12]`) and angles (`[f64; 8]`) directly from `armc_deg`, `lat_deg`,
    /// and `eps_deg`, matching `swe.houses_armc(armc, lat, eps, hsys)`.
    pub fn houses_armc(
        &self,
        armc_deg: f64,
        lat_deg: f64,
        eps_deg: f64,
        hsys: char,
    ) -> Result<([f64; 12], [f64; 8]), EphemerisError> {
        self.houses_armc_ex(armc_deg, lat_deg, eps_deg, hsys, None)
    }

    /// Evaluates house cusps (`[f64; 12]`) and angles (`[f64; 8]`) from `armc_deg`, `lat_deg`,
    /// `eps_deg`, and optional solar declination `sun_decl_deg` (`ascmc[9]` in `swe_houses_armc`
    /// for Sunshine house systems `'I'` and `'i'`).
    pub fn houses_armc_ex(
        &self,
        armc_deg: f64,
        lat_deg: f64,
        eps_deg: f64,
        hsys: char,
        sun_decl_deg: Option<f64>,
    ) -> Result<([f64; 12], [f64; 8]), EphemerisError> {
        let hsys = if hsys == 'E' { 'A' } else { hsys };
        let req_sys = HouseSystem::from_code(hsys)?;
        let (system, decl_deg) = if let Some(d) = sun_decl_deg {
            (req_sys, d)
        } else if req_sys == HouseSystem::SunshineTreindl {
            (HouseSystem::Porphyry, 0.0)
        } else if req_sys == HouseSystem::SunshineMakransky {
            match self.polar_policy {
                PolarPolicy::Strict => {
                    return Err(EphemerisError::UndefinedHouseSystem {
                        house_system: hsys,
                        latitude_deg: lat_deg,
                        reason: "Sunshine house systems ('I', 'i') require solar declination (use houses_ex with JD_UT or pass ascmc9 to houses_armc)".to_owned(),
                    });
                }
                PolarPolicy::FallbackToWholeSign => (HouseSystem::WholeSign, 0.0),
                PolarPolicy::FallbackToPorphyry => (HouseSystem::Porphyry, 0.0),
            }
        } else {
            (req_sys, 0.0)
        };
        let opts = HouseOptions {
            system,
            polar_policy: self.polar_policy,
            zodiac: Zodiac::Tropical,
            nutation: NutationMode::True,
            include_speed: false,
        };
        let res = self
            .engine
            .houses_armc(armc_deg, lat_deg, eps_deg, decl_deg, opts)?;
        let ascmc = [
            res.angles.asc_deg,
            res.angles.mc_deg,
            res.angles.armc_deg,
            res.angles.vertex_deg,
            res.angles.equatorial_asc_deg,
            res.angles.co_asc_koch_deg,
            res.angles.co_asc_munkasey_deg,
            res.angles.polar_asc_deg,
        ];
        Ok((res.cusps_deg, ascmc))
    }

    /// Evaluates the configured sidereal ayanamsa at `jd_ut` (`UT1`) with `flags`, returning `(retflag, ayanamsa_deg)`
    /// matching `swe.get_ayanamsa_ex_ut(jd_ut, flags)`.
    pub fn get_ayanamsa_ex_ut(&self, jd_ut: f64, flags: u32) -> Result<(i32, f64), EphemerisError> {
        let val = self.swe_ayanamsa_val(TimeInput::from_ut1_jd(jd_ut), flags)?;
        let retflag = if (flags & SEFLG_MOSEPH != 0) || self.engine.is_moshier_mode() {
            ((flags & !(SEFLG_JPLEPH | SEFLG_SWIEPH)) | SEFLG_MOSEPH) as i32
        } else {
            (flags | SEFLG_SWIEPH) as i32
        };
        Ok((retflag, val.effective_deg))
    }

    /// Evaluates the mean sidereal ayanamsa in degrees at `jd_ut` (`UT1`), matching `swe.get_ayanamsa_ut(jd_ut)`.
    pub fn get_ayanamsa_ut(&self, jd_ut: f64) -> Result<f64, EphemerisError> {
        let val = self.swe_ayanamsa_val(TimeInput::from_ut1_jd(jd_ut), SEFLG_NONUT)?;
        Ok(val.mean_deg)
    }

    /// Evaluates the mean sidereal ayanamsa in degrees at `jd_et` (`TT`), matching `swe.get_ayanamsa(jd_et)`.
    pub fn get_ayanamsa(&self, jd_et: f64) -> Result<f64, EphemerisError> {
        let val = self.swe_ayanamsa_val(TimeInput::from_tt_jd(jd_et), SEFLG_NONUT)?;
        Ok(val.mean_deg)
    }

    /// Evaluates the sidereal ayanamsa in degrees at `jd_et` (`TT`) with flags (`SEFLG_NONUT` -> mean,
    /// default -> true nutated), returning `(retflag, ayanamsa_deg)` matching `swe.get_ayanamsa_ex(jd_et, flags)`.
    pub fn get_ayanamsa_ex(&self, jd_et: f64, flags: u32) -> Result<(i32, f64), EphemerisError> {
        let val = self.swe_ayanamsa_val(TimeInput::from_tt_jd(jd_et), flags)?;
        let retflag = if (flags & SEFLG_MOSEPH != 0) || self.engine.is_moshier_mode() {
            ((flags & !(SEFLG_JPLEPH | SEFLG_SWIEPH)) | SEFLG_MOSEPH) as i32
        } else {
            (flags | SEFLG_SWIEPH) as i32
        };
        Ok((retflag, val.effective_deg))
    }

    /// Returns $\Delta T = \text{TT} - \text{UT1}$ in **days** at `jd_ut`, matching `swe.deltat(jd_ut)`.
    #[must_use]
    pub fn deltat(&self, jd_ut: f64) -> f64 {
        self.engine
            .delta_t(TimeInput::from_ut1_jd(jd_ut))
            .delta_t_sec
            / SECONDS_PER_DAY
    }

    /// Returns $\Delta T = \text{TT} - \text{UT1}$ in **seconds** at `jd_ut`.
    #[must_use]
    pub fn deltat_sec(&self, jd_ut: f64) -> f64 {
        self.engine
            .delta_t(TimeInput::from_ut1_jd(jd_ut))
            .delta_t_sec
    }

    /// Returns Greenwich Apparent Sidereal Time (`GAST`) in **hours** `[0, 24)` at `jd_ut`, matching `swe.sidtime(jd_ut)`.
    #[must_use]
    pub fn sidtime(&self, jd_ut: f64) -> f64 {
        self.engine
            .sidtime(TimeInput::from_ut1_jd(jd_ut), NutationMode::True, 0.0)
            .gast_hours
    }

    /// Returns Greenwich Sidereal Time in **hours** `[0, 24)` at `jd_ut` with explicit obliquity `eps_deg` and
    /// nutation in longitude `nut_lon_deg` (both in degrees), matching `swe.sidtime0(jd_ut, eps, nut)`.
    #[must_use]
    pub fn sidtime0(&self, jd_ut: f64, eps_deg: f64, nut_lon_deg: f64) -> f64 {
        let gmst_hours = self
            .engine
            .sidtime(TimeInput::from_ut1_jd(jd_ut), NutationMode::Mean, 0.0)
            .gmst_hours;
        let eq_eq_hours = (nut_lon_deg * eps_deg.to_radians().cos()) / 15.0;
        (gmst_hours + eq_eq_hours).rem_euclid(24.0)
    }

    /// Evaluates horizontal coordinates (`azalt`) at `jd_ut` (`UT1`) for `geopos = [lon_deg, lat_deg, alt_m]`,
    /// `calc_flag` (`SE_ECL2HOR` or `SE_EQU2HOR`), and `xin = [lon_or_ra_deg, lat_or_dec_deg]`.
    ///
    /// Note: Following Swiss Ephemeris `swe_azalt` convention (`R2.6`), `azimuth_deg` in `SweAdapter::azalt`
    /// is measured **from South towards West (`0° = S, 90° = W, 180° = N, 270° = E`)**, and `atpress <= 0.0`
    /// estimates standard barometric pressure from elevation `geopos[2]`.
    pub fn azalt(
        &self,
        jd_ut: f64,
        calc_flag: i32,
        geopos: [f64; 3],
        atpress: f64,
        attemp: f64,
        xin: [f64; 2],
    ) -> Result<HorizontalCoords, EphemerisError> {
        let obs = Observer::new(geopos[1], geopos[0], geopos[2])?;
        let pressure_hpa = if atpress <= 0.0 {
            let h = geopos[2].clamp(-500.0, 11_000.0);
            1013.25 * (1.0 - 0.0065 * h / 288.15).powf(5.25588)
        } else {
            atpress
        };
        let atmo = AtmosphereParams {
            pressure_hpa,
            temperature_c: attemp,
            ..AtmosphereParams::default()
        };
        let t_in = TimeInput::from_ut1_jd(jd_ut);
        let mut hz = if calc_flag == SE_ECL2HOR {
            self.engine
                .azalt_ecliptic(obs, t_in, xin[0], xin[1], atmo)
        } else {
            self.engine
                .azalt_equatorial(obs, t_in, xin[0], xin[1], atmo)
        };
        // Convert North->East azimuth to Swiss Ephemeris South->West azimuth
        hz.azimuth_deg = (hz.azimuth_deg + 180.0).rem_euclid(360.0);
        Ok(hz)
    }

    /// Searches for the next rise, set, or meridian transit (`JD_UT1`), matching `swe.rise_trans`.
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans(
        &self,
        jd_ut: f64,
        ipl: i32,
        starname: Option<&str>,
        epheflag: u32,
        rsmi: i32,
        geopos: [f64; 3],
        atpress: f64,
        attemp: f64,
    ) -> Result<Option<f64>, EphemerisError> {
        let obs = Observer::new(geopos[1], geopos[0], geopos[2])?;
        let body = if let Some(s) = starname.filter(|s| !s.trim().is_empty()) {
            let entry = self.engine.resolve_star_entry(s)?;
            BodyId::FixedStar(entry.hip_id())
        } else {
            body_id_from_swe_ipl(ipl)?
        };
        let ev = if (rsmi & SE_CALC_SET) != 0 {
            RiseTransEvent::Set
        } else if (rsmi & SE_CALC_MTRANSIT) != 0 {
            RiseTransEvent::UpperTransit
        } else if (rsmi & SE_CALC_ITRANSIT) != 0 {
            RiseTransEvent::LowerTransit
        } else {
            RiseTransEvent::Rise
        };
        let custom_horizon_alt_deg = if (rsmi & SE_BIT_CIVIL_TWILIGHT) != 0 {
            Some(-6.0)
        } else if (rsmi & SE_BIT_NAUTIC_TWILIGHT) != 0 {
            Some(-12.0)
        } else if (rsmi & SE_BIT_ASTRO_TWILIGHT) != 0 {
            Some(-18.0)
        } else {
            None
        };
        let disc_mode = if custom_horizon_alt_deg.is_some() || (rsmi & SE_BIT_DISC_CENTER) != 0 {
            DiscMode::Center
        } else if (rsmi & SE_BIT_DISC_BOTTOM) != 0 {
            DiscMode::BottomLimb
        } else {
            DiscMode::UpperLimb
        };
        let pressure_hpa = if custom_horizon_alt_deg.is_some() || (rsmi & SE_BIT_NO_REFRACTION) != 0
        {
            0.0
        } else if atpress <= 0.0 {
            let h = geopos[2].clamp(-500.0, 11_000.0);
            1013.25 * (1.0 - 0.0065 * h / 288.15).powf(5.25588)
        } else {
            atpress
        };
        let opts = RiseTransOptions {
            max_search_days: 2.0,
            custom_horizon_alt_deg,
            disc_mode,
            pressure_hpa,
            temperature_c: attemp,
        };
        if (epheflag & SEFLG_MOSEPH) != 0 || self.engine.is_moshier_mode() {
            return self.engine.rise_trans_moshier(body, obs, jd_ut, ev, opts);
        }
        match self.engine.rise_trans(body, obs, jd_ut, ev, opts) {
            Ok(res) => Ok(res),
            Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                match self.engine.rise_trans_moshier(body, obs, jd_ut, ev, opts) {
                    Ok(res) => Ok(res),
                    Err(_) => Err(orig_err),
                }
            }
            Err(e) => Err(e),
        }
    }

    /// Converts calendar date `(year, month, day, hour)` to Julian Date (`gregflag == 1` Gregorian, `0` Julian).
    pub fn julday(
        year: i32,
        month: u8,
        day: u8,
        hour: f64,
        gregflag: i32,
    ) -> Result<f64, EphemerisError> {
        let cal = if gregflag == SE_JUL_CAL {
            CalendarKind::Julian
        } else {
            CalendarKind::Gregorian
        };
        julday(year, month, day, hour, cal)
    }

    /// Inverse calendar conversion from `jd` to [`CalendarDate`] (`gregflag == 1` Gregorian, `0` Julian).
    pub fn revjul(jd: f64, gregflag: i32) -> Result<CalendarDate, EphemerisError> {
        let cal = if gregflag == SE_JUL_CAL {
            CalendarKind::Julian
        } else {
            CalendarKind::Gregorian
        };
        revjul(jd, cal)
    }

    /// Evaluates planetary, solar, lunar, or minor-body phenomena at `jd_ut` (`UT1`), returning
    /// `(attr[0..20], retflag)` matching `swe_pheno_ut(jd_ut, ipl, flags)`:
    /// - `attr[0]`: phase angle $\alpha$ (deg)
    /// - `attr[1]`: phase (illuminated fraction $k \in [0, 1]$)
    /// - `attr[2]`: elongation $\psi$ (deg)
    /// - `attr[3]`: apparent diameter of disc $d$ (deg)
    /// - `attr[4]`: apparent visual magnitude $V$
    /// - `attr[5]`: equatorial horizontal parallax (deg)
    pub fn pheno_ut(
        &self,
        jd_ut: f64,
        ipl: i32,
        flags: u32,
    ) -> Result<([f64; 20], i32), EphemerisError> {
        let body = body_id_from_swe_ipl(ipl)?;
        let opts = calc_options_from_swe_flags(flags, self.topo_observer, self.sid_mode)?;
        let t_in = TimeInput::from_ut1_jd(jd_ut);
        let use_moshier = (flags & SEFLG_MOSEPH != 0) || self.engine.is_moshier_mode();
        let (ph, used_moshier) = if use_moshier {
            (self.engine.phenomena_moshier(body, t_in, opts)?, true)
        } else {
            match self.engine.phenomena(body, t_in, opts) {
                Ok(p) => (p, false),
                Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                    match self.engine.phenomena_moshier(body, t_in, opts) {
                        Ok(p) => (p, true),
                        Err(_) => return Err(orig_err),
                    }
                }
                Err(e) => return Err(e),
            }
        };
        let retflag = if used_moshier {
            ((flags & !(SEFLG_JPLEPH | SEFLG_SWIEPH)) | SEFLG_MOSEPH) as i32
        } else {
            (flags | SEFLG_SWIEPH) as i32
        };
        Ok((ph.to_swe_attr20(), retflag))
    }

    /// Evaluates planetary, solar, lunar, or minor-body phenomena at `jd_et` (`TT`), returning
    /// `(attr[0..20], retflag)` matching `swe_pheno(jd_et, ipl, flags)`.
    pub fn pheno(
        &self,
        jd_et: f64,
        ipl: i32,
        flags: u32,
    ) -> Result<([f64; 20], i32), EphemerisError> {
        let body = body_id_from_swe_ipl(ipl)?;
        let opts = calc_options_from_swe_flags(flags, self.topo_observer, self.sid_mode)?;
        let t_in = TimeInput::from_tt_jd(jd_et);
        let use_moshier = (flags & SEFLG_MOSEPH != 0) || self.engine.is_moshier_mode();
        let (ph, used_moshier) = if use_moshier {
            (self.engine.phenomena_moshier(body, t_in, opts)?, true)
        } else {
            match self.engine.phenomena(body, t_in, opts) {
                Ok(p) => (p, false),
                Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                    match self.engine.phenomena_moshier(body, t_in, opts) {
                        Ok(p) => (p, true),
                        Err(_) => return Err(orig_err),
                    }
                }
                Err(e) => return Err(e),
            }
        };
        let retflag = if used_moshier {
            ((flags & !(SEFLG_JPLEPH | SEFLG_SWIEPH)) | SEFLG_MOSEPH) as i32
        } else {
            (flags | SEFLG_SWIEPH) as i32
        };
        Ok((ph.to_swe_attr20(), retflag))
    }

    /// Evaluates planetary or lunar nodes and apsides at `jd_ut` (`UT1`), returning
    /// `(xnasc[0..6], xndsc[0..6], xperi[0..6], xaphe[0..6], retflag)` matching `swe_nod_aps_ut(jd_ut, ipl, flags, method)`.
    #[allow(clippy::type_complexity)]
    pub fn nod_aps_ut(
        &self,
        jd_ut: f64,
        ipl: i32,
        flags: u32,
        method: u32,
    ) -> Result<([f64; 6], [f64; 6], [f64; 6], [f64; 6], i32), EphemerisError> {
        let body = body_id_from_swe_ipl(ipl)?;
        let opts = calc_options_from_swe_flags(flags, self.topo_observer, self.sid_mode)?;
        let t_in = TimeInput::from_ut1_jd(jd_ut);
        let use_moshier = (flags & SEFLG_MOSEPH != 0) || self.engine.is_moshier_mode();
        let (raw_res, used_moshier) = if use_moshier {
            (
                self.engine.nod_aps_moshier(body, t_in, opts, method)?,
                true,
            )
        } else {
            match self.engine.nod_aps(body, t_in, opts, method) {
                Ok(r) => (r, false),
                Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                    match self.engine.nod_aps_moshier(body, t_in, opts, method) {
                        Ok(r) => (r, true),
                        Err(_) => return Err(orig_err),
                    }
                }
                Err(e) => return Err(e),
            }
        };
        let res = raw_res.format_with_flags(flags);
        let retflag = if used_moshier {
            ((flags & !(SEFLG_JPLEPH | SEFLG_SWIEPH)) | SEFLG_MOSEPH) as i32
        } else {
            (flags | SEFLG_SWIEPH) as i32
        };
        Ok((
            res.asc_node,
            res.desc_node,
            res.periapsis,
            res.apoapsis,
            retflag,
        ))
    }

    /// Evaluates planetary or lunar nodes and apsides at `jd_et` (`TT`), returning
    /// `(xnasc[0..6], xndsc[0..6], xperi[0..6], xaphe[0..6], retflag)` matching `swe_nod_aps(jd_et, ipl, flags, method)`.
    #[allow(clippy::type_complexity)]
    pub fn nod_aps(
        &self,
        jd_et: f64,
        ipl: i32,
        flags: u32,
        method: u32,
    ) -> Result<([f64; 6], [f64; 6], [f64; 6], [f64; 6], i32), EphemerisError> {
        let body = body_id_from_swe_ipl(ipl)?;
        let opts = calc_options_from_swe_flags(flags, self.topo_observer, self.sid_mode)?;
        let t_in = TimeInput::from_tt_jd(jd_et);
        let use_moshier = (flags & SEFLG_MOSEPH != 0) || self.engine.is_moshier_mode();
        let (raw_res, used_moshier) = if use_moshier {
            (
                self.engine.nod_aps_moshier(body, t_in, opts, method)?,
                true,
            )
        } else {
            match self.engine.nod_aps(body, t_in, opts, method) {
                Ok(r) => (r, false),
                Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                    match self.engine.nod_aps_moshier(body, t_in, opts, method) {
                        Ok(r) => (r, true),
                        Err(_) => return Err(orig_err),
                    }
                }
                Err(e) => return Err(e),
            }
        };
        let res = raw_res.format_with_flags(flags);
        let retflag = if used_moshier {
            ((flags & !(SEFLG_JPLEPH | SEFLG_SWIEPH)) | SEFLG_MOSEPH) as i32
        } else {
            (flags | SEFLG_SWIEPH) as i32
        };
        Ok((
            res.asc_node,
            res.desc_node,
            res.periapsis,
            res.apoapsis,
            retflag,
        ))
    }

    /// Searches for the next global Lunar Occultation of a planet or fixed star starting from `jd_ut_start`,
    /// returning `(retflag, tret[0..10])` matching `swe_lun_occult_when_glob`.
    pub fn lun_occult_when_glob(
        &self,
        jd_ut_start: f64,
        ipl: i32,
        starname: Option<&str>,
        flags: u32,
        backward: bool,
    ) -> Result<(i32, [f64; 10]), EphemerisError> {
        self.lun_occult_when_glob_ex(jd_ut_start, ipl, starname, flags, 0, backward)
    }

    /// Searches for the next global Lunar Occultation of a planet or fixed star matching `ifltype` (`0` = any)
    /// starting from `jd_ut_start`, returning `(retflag, tret[0..10])` matching `swe_lun_occult_when_glob`.
    pub fn lun_occult_when_glob_ex(
        &self,
        jd_ut_start: f64,
        ipl: i32,
        starname: Option<&str>,
        _flags: u32,
        ifltype: i32,
        backward: bool,
    ) -> Result<(i32, [f64; 10]), EphemerisError> {
        let target = if let Some(s) = starname.filter(|s| !s.trim().is_empty()) {
            let meta = self.engine.star_info(s)?;
            BodyId::FixedStar(crate::types::HipId(meta.hip))
        } else {
            body_id_from_swe_ipl(ipl)?
        };
        let kind_mask =
            ifltype & (SE_ECL_TOTAL | SE_ECL_ANNULAR | SE_ECL_PARTIAL | SE_ECL_ANNULAR_TOTAL);
        let cent_mask = ifltype & (SE_ECL_CENTRAL | SE_ECL_NONCENTRAL);
        let mut cursor = jd_ut_start;
        for _ in 0..600 {
            let info = self.engine.next_lunar_occultation(
                target,
                TimeInput::from_ut1_jd(cursor),
                backward,
                3652.5,
            )?;
            let retflag = info.swe_retflag();
            let matches_kind = kind_mask == 0 || (retflag & kind_mask) != 0;
            let matches_cent = cent_mask == 0 || (retflag & cent_mask) != 0;
            if matches_kind && matches_cent {
                return Ok((retflag, info.to_swe_glob_tret10()));
            }
            cursor = info.jd_ut1_max + if backward { -5.0 } else { 5.0 };
        }
        Err(EphemerisError::NonConvergence {
            solver: "lun_occult_when_glob".to_owned(),
            iterations: 600,
            residual: 0.0,
        })
    }

    /// Searches for the next locally visible Lunar Occultation of a planet or fixed star at `geopos = [lon, lat, alt]`
    /// starting from `jd_ut_start`, returning `(retflag, tret[0..10], attr[0..20])` matching `swe_lun_occult_when_loc`.
    pub fn lun_occult_when_loc(
        &self,
        jd_ut_start: f64,
        ipl: i32,
        starname: Option<&str>,
        flags: u32,
        geopos: [f64; 3],
        backward: bool,
    ) -> Result<(i32, [f64; 10], [f64; 20]), EphemerisError> {
        let obs = Observer::new(geopos[1], geopos[0], geopos[2])?;
        let target = if let Some(s) = starname.filter(|s| !s.trim().is_empty()) {
            let meta = self.engine.star_info(s)?;
            BodyId::FixedStar(crate::types::HipId(meta.hip))
        } else {
            body_id_from_swe_ipl(ipl)?
        };
        let info = self.engine.next_local_lunar_occultation(
            target,
            obs,
            TimeInput::from_ut1_jd(jd_ut_start),
            backward,
            3652.5,
        )?;
        let mut tret = info.to_swe_tret10();
        let mut attr = [0.0; 20];
        let mut rf = match info.kind {
            crate::occultations::OccultationKind::Total => SE_ECL_TOTAL,
            crate::occultations::OccultationKind::Annular => SE_ECL_ANNULAR,
            crate::occultations::OccultationKind::Partial
            | crate::occultations::OccultationKind::Grazing => SE_ECL_PARTIAL,
        };
        let is_vis_moon = |t_ut: f64| -> bool {
            if t_ut <= 0.0 {
                return false;
            }
            self.engine
                .calc_azalt(
                    BodyId::Planet(PlanetId::Moon),
                    TimeInput::from_ut1_jd(t_ut),
                    obs,
                    AtmosphereParams::default(),
                )
                .is_ok_and(|hz| hz.altitude_apparent_deg >= -0.266)
        };
        let vis_max = is_vis_moon(tret[0]);
        let vis_c1 = is_vis_moon(tret[1]);
        let vis_c2 = is_vis_moon(tret[2]);
        let vis_c3 = is_vis_moon(tret[3]);
        let vis_c4 = is_vis_moon(tret[4]);
        if vis_max {
            rf |= SE_ECL_VISIBLE | SE_ECL_MAX_VISIBLE;
        }
        if vis_c1 {
            rf |= SE_ECL_VISIBLE | SE_ECL_1ST_VISIBLE;
        }
        if vis_c2 {
            rf |= SE_ECL_VISIBLE | SE_ECL_2ND_VISIBLE;
        }
        if vis_c3 {
            rf |= SE_ECL_VISIBLE | SE_ECL_3RD_VISIBLE;
        }
        if vis_c4 {
            rf |= SE_ECL_VISIBLE | SE_ECL_4TH_VISIBLE;
        }
        if !vis_c1 && vis_c4 && tret[1] > 0.0 {
            if let Ok(Some(t_rise)) = self.rise_trans(
                tret[1] - 0.05,
                SE_MOON,
                None,
                flags,
                SE_CALC_RISE,
                geopos,
                0.0,
                15.0,
            ) {
                if t_rise >= tret[1] - 0.01 && t_rise <= tret[4] + 0.01 {
                    tret[5] = t_rise;
                }
            }
        }
        if vis_c1 && !vis_c4 && tret[1] > 0.0 {
            if let Ok(Some(t_set)) = self.rise_trans(
                tret[1] - 0.05,
                SE_MOON,
                None,
                flags,
                SE_CALC_SET,
                geopos,
                0.0,
                15.0,
            ) {
                if t_set >= tret[1] - 0.01 && t_set <= tret[4] + 0.01 {
                    tret[6] = t_set;
                }
            }
        }
        if tret[1] > 0.0 {
            if let Ok(hz_sun) = self.engine.calc_azalt(
                BodyId::Planet(PlanetId::Sun),
                TimeInput::from_ut1_jd(tret[1]),
                obs,
                AtmosphereParams::default(),
            ) {
                if hz_sun.altitude_apparent_deg >= -0.266 {
                    rf |= SE_ECL_OCC_BEG_DAYLIGHT;
                }
            }
        }
        if tret[4] > 0.0 {
            if let Ok(hz_sun) = self.engine.calc_azalt(
                BodyId::Planet(PlanetId::Sun),
                TimeInput::from_ut1_jd(tret[4]),
                obs,
                AtmosphereParams::default(),
            ) {
                if hz_sun.altitude_apparent_deg >= -0.266 {
                    rf |= SE_ECL_OCC_END_DAYLIGHT;
                }
            }
        }
        attr[0] = info.magnitude;
        attr[2] = if info.kind == crate::occultations::OccultationKind::Total {
            1.0
        } else {
            info.magnitude.clamp(0.0, 1.0)
        };
        attr[3] = -3474.8;
        if let Ok(hz_moon) = self.engine.calc_azalt(
            BodyId::Planet(PlanetId::Moon),
            TimeInput::from_ut1_jd(tret[0]),
            obs,
            AtmosphereParams::default(),
        ) {
            attr[4] = (hz_moon.azimuth_deg + 180.0).rem_euclid(360.0);
            attr[5] = hz_moon.altitude_true_deg;
            attr[6] = hz_moon.altitude_apparent_deg;
        }
        let topo_opts = CalcOptions {
            origin: Origin::Topocentric(obs),
            plane: Plane::Equatorial,
            equinox: Equinox::OfDate,
            nutation: NutationMode::True,
            aberration: AbMode::Apparent,
            zodiac: Zodiac::Tropical,
            include_speed: false,
        };
        let t_max_in = TimeInput::from_ut1_jd(tret[0]);
        let res_m = self
            .engine
            .calc(BodyId::Planet(PlanetId::Moon), t_max_in, topo_opts)
            .or_else(|_| {
                self.engine
                    .calc_moshier(BodyId::Planet(PlanetId::Moon), t_max_in, topo_opts)
            });
        let res_t = self
            .engine
            .calc(target, t_max_in, topo_opts)
            .or_else(|_| self.engine.calc_moshier(target, t_max_in, topo_opts));
        if let (Ok(rm), Ok(rt)) = (res_m, res_t) {
            let dra =
                ((rm.lon_deg - rt.lon_deg + 180.0).rem_euclid(360.0) - 180.0) * rm.lat_deg.to_radians().cos();
            let ddec = rm.lat_deg - rt.lat_deg;
            attr[7] = dra.hypot(ddec);
            let r_moon_deg = (1737.4 / (rm.dist_au * crate::types::AU_KM).max(1000.0))
                .clamp(-1.0, 1.0)
                .asin()
                .to_degrees();
            if let BodyId::Planet(p) = target {
                let r_pl_km = match p {
                    PlanetId::Mercury => 2439.7,
                    PlanetId::Venus => 6051.8,
                    PlanetId::Mars => 3389.5,
                    PlanetId::Jupiter => 69911.0,
                    PlanetId::Saturn => 58232.0,
                    PlanetId::Uranus => 25362.0,
                    PlanetId::Neptune => 24622.0,
                    PlanetId::Pluto => 1188.3,
                    _ => 0.0,
                };
                if r_pl_km > 0.0 {
                    let r_t_deg = (r_pl_km / (rt.dist_au * crate::types::AU_KM).max(1000.0))
                        .clamp(-1.0, 1.0)
                        .asin()
                        .to_degrees();
                    if r_t_deg > 0.0 {
                        attr[1] = r_moon_deg / r_t_deg;
                    }
                }
            }
        }
        Ok((rf, tret, attr))
    }

    /// Computes the 3D fractional mundane house position $H \in [1.0, 13.0)$ of `xpin = [ecl_lon_deg, ecl_lat_deg]`,
    /// matching `swe_house_pos(armc, geolat, eps, hsys, xpin)`.
    pub fn house_pos(
        &self,
        armc_deg: f64,
        geolat_deg: f64,
        eps_deg: f64,
        hsys: char,
        xpin: [f64; 2],
    ) -> Result<f64, EphemerisError> {
        let hsys = if hsys == 'E' { 'A' } else { hsys };
        let system = HouseSystem::from_code(hsys)?;
        self.engine
            .house_pos(armc_deg, geolat_deg, eps_deg, system, xpin[0], xpin[1])
    }

    /// Evaluates the Equation of Time ($E = \text{Apparent Solar Time} - \text{Mean Solar Time}$)
    /// in **days** at `jd_ut` (`UT1`), matching `swe_time_equ(jd_ut)`.
    pub fn time_equ(&self, jd_ut: f64) -> Result<f64, EphemerisError> {
        self.engine.time_equation(TimeInput::from_ut1_jd(jd_ut))
    }

    /// Inverse horizontal coordinate transformation (`swe_azalt_rev`): converts `xin = [azimuth_south_west_deg, true_alt_deg]`
    /// at `jd_ut` (`UT1`) and `geopos = [lon_deg, lat_deg, alt_m]` into `[ecl_lon_deg, ecl_lat_deg]` (`calc_flag == SE_HOR2ECL`)
    /// or `[ra_deg, dec_deg]` (`calc_flag == SE_HOR2EQU`).
    pub fn azalt_rev(
        &self,
        jd_ut: f64,
        calc_flag: i32,
        geopos: [f64; 3],
        xin: [f64; 2],
    ) -> Result<[f64; 2], EphemerisError> {
        let obs = Observer::new(geopos[1], geopos[0], geopos[2])?;
        let az_north_east_deg = (xin[0] + 180.0).rem_euclid(360.0);
        let plane = if calc_flag == SE_HOR2ECL {
            Plane::Ecliptic
        } else {
            Plane::Equatorial
        };
        Ok(self.engine.azalt_rev(
            obs,
            TimeInput::from_ut1_jd(jd_ut),
            az_north_east_deg,
            xin[1],
            plane,
        ))
    }

    /// Legacy no-op compatibility method for `swe.set_ephe_path(path)`.
    pub fn set_ephe_path(&mut self, _path: &str) {}

    /// Legacy no-op compatibility method for `swe.close()`.
    pub fn close(&mut self) {}

    /// Returns the human-readable planet/body name for legacy `ipl` (`swe.get_planet_name(ipl)`),
    /// or `String::new()` if `ipl` is unknown or an unbundled MPC number.
    #[must_use]
    pub fn get_planet_name(ipl: i32) -> String {
        match body_id_from_swe_ipl(ipl) {
            Ok(BodyId::MinorPlanet(mpc)) => MINOR_BODIES_68
                .binary_search_by_key(&mpc, |s| s.mpc)
                .ok()
                .and_then(|idx| MINOR_BODIES_68.get(idx))
                .map_or_else(String::new, |s| s.name.to_owned()),
            Ok(body) => body_display_name(body),
            Err(_) => String::new(),
        }
    }

    /// Returns the canonical Swiss Ephemeris ayanamsa name for `sid_mode` (`0..=46`),
    /// matching `swe_get_ayanamsa_name(isidmode)` (or `""` if unknown).
    #[must_use]
    pub const fn get_ayanamsa_name(sid_mode: i32) -> &'static str {
        match sid_mode {
            0 => "Fagan/Bradley",
            1 => "Lahiri",
            2 => "De Luce",
            3 => "Raman",
            4 => "Usha/Shashi",
            5 => "Krishnamurti",
            6 => "Djwhal Khul",
            7 => "Yukteshwar",
            8 => "J.N. Bhasin",
            9 => "Babylonian/Kugler 1",
            10 => "Babylonian/Kugler 2",
            11 => "Babylonian/Kugler 3",
            12 => "Babylonian/Huber",
            13 => "Babylonian/Eta Piscium",
            14 => "Babylonian/Aldebaran = 15 Tau",
            15 => "Hipparchos",
            16 => "Sassanian",
            17 => "Galact. Center = 0 Sag",
            18 => "J2000",
            19 => "J1900",
            20 => "B1950",
            21 => "Suryasiddhanta",
            22 => "Suryasiddhanta, mean Sun",
            23 => "Aryabhata",
            24 => "Aryabhata, mean Sun",
            25 => "SS Revati",
            26 => "SS Citra",
            27 => "True Citra",
            28 => "True Revati",
            29 => "True Pushya (PVRN Rao)",
            30 => "Galactic Center (Gil Brand)",
            31 => "Galactic Equator (IAU1958)",
            32 => "Galactic Equator",
            33 => "Galactic Equator mid-Mula",
            34 => "Skydram (Mardyks)",
            35 => "True Mula (Chandra Hari)",
            36 => "Dhruva/Gal.Center/Mula (Wilhelm)",
            37 => "Aryabhata 522",
            38 => "Babylonian/Britton",
            39 => "\"Vedic\"/Sheoran",
            40 => "Cochrane (Gal.Center = 0 Cap)",
            41 => "Galactic Equator (Fiorenza)",
            42 => "Vettius Valens",
            43 => "Lahiri 1940",
            44 => "Lahiri VP285",
            45 => "Krishnamurti-Senthilathiban",
            46 => "Lahiri ICRC",
            _ => "",
        }
    }

    /// Returns the human-readable house system name for `hsys` (`swe.house_name(hsys)`), or `"Unknown"`.
    #[must_use]
    pub const fn house_name(hsys: char) -> &'static str {
        match hsys {
            'P' => "Placidus",
            'W' => "equal/ whole sign",
            'K' => "Koch",
            'A' | 'E' => "equal",
            'R' => "Regiomontanus",
            'C' => "Campanus",
            'V' => "equal/Vehlow",
            'D' => "equal (MC)",
            'N' => "equal/1=Aries",
            'M' => "Morinus",
            'B' => "Alcabitius",
            'O' => "Porphyry",
            'S' => "Sripati",
            'T' => "Polich/Page",
            'X' => "axial rotation system/Meridian houses",
            'H' => "horizon/azimut",
            'I' => "Sunshine",
            'i' => "Sunshine/alt.",
            'F' => "Carter poli-equ.",
            'L' => "Pullen SD",
            'Q' => "Pullen SR",
            'U' => "Krusinski-Pisa-Goelzer",
            'Y' => "APC houses",
            'G' => "Gauquelin sectors",
            _ => "Unknown",
        }
    }

    /// Signed angular difference `p1 - p2` in degrees normalized to `[-180, 180)`,
    /// matching Swiss Ephemeris `swe_difdeg2n(p1, p2)`.
    #[must_use]
    pub fn difdeg2n(p1: f64, p2: f64) -> f64 {
        (p1 - p2 + 180.0).rem_euclid(360.0) - 180.0
    }

    /// Normalizes an angle in degrees to `[0, 360)`, matching Swiss Ephemeris `swe_degnorm(x)`.
    #[must_use]
    pub fn degnorm(x: f64) -> f64 {
        x.rem_euclid(360.0)
    }

    /// Coordinate transformation between ecliptic and equatorial spherical coordinates
    /// `(lon_or_ra_deg, lat_or_dec_deg, dist)` with obliquity `eps` (degrees),
    /// matching Swiss Ephemeris `swe_cotrans(xpo, xpn, eps)`:
    /// - `eps > 0`: equatorial -> ecliptic
    /// - `eps < 0`: ecliptic -> equatorial
    #[must_use]
    pub fn cotrans(lon: f64, lat: f64, dist: f64, eps: f64) -> (f64, f64, f64) {
        let lon_rad = lon.to_radians();
        let lat_rad = lat.to_radians();
        let eps_rad = eps.to_radians();
        let x = lat_rad.cos() * lon_rad.cos();
        let y = lat_rad.cos() * lon_rad.sin();
        let z = lat_rad.sin();
        let y_new = y * eps_rad.cos() + z * eps_rad.sin();
        let z_new = -y * eps_rad.sin() + z * eps_rad.cos();
        let out_lon = y_new.atan2(x).to_degrees().rem_euclid(360.0);
        let out_lat = z_new.clamp(-1.0, 1.0).asin().to_degrees();
        (out_lon, out_lat, dist)
    }

    /// Searches forward (`backward == false`) or backward (`backward == true`) from `jd_ut_start`
    /// for the next global Solar Eclipse matching `ifltype` (`0` = any), returning `(retflag, tret[0..10])`
    /// matching `swe.sol_eclipse_when_glob`.
    pub fn sol_eclipse_when_glob(
        &self,
        jd_ut_start: f64,
        _flags: u32,
        ifltype: i32,
        backward: bool,
    ) -> Result<(i32, [f64; 10]), EphemerisError> {
        let kind_mask =
            ifltype & (SE_ECL_TOTAL | SE_ECL_ANNULAR | SE_ECL_PARTIAL | SE_ECL_ANNULAR_TOTAL);
        let cent_mask = ifltype & (SE_ECL_CENTRAL | SE_ECL_NONCENTRAL);
        let mut cursor = jd_ut_start;
        for _ in 0..500 {
            let info = self.engine.next_solar_eclipse(
                TimeInput::from_ut1_jd(cursor),
                backward,
                None,
                EclipseOptions::default(),
            )?;
            let cent_bit = if info.central_latitude_deg.is_some() {
                SE_ECL_CENTRAL
            } else {
                SE_ECL_NONCENTRAL
            };
            let retflag = match info.kind {
                SolarEclipseKind::Total => SE_ECL_TOTAL | cent_bit,
                SolarEclipseKind::Annular => SE_ECL_ANNULAR | cent_bit,
                SolarEclipseKind::Hybrid => SE_ECL_ANNULAR_TOTAL | cent_bit,
                SolarEclipseKind::Partial => SE_ECL_PARTIAL | SE_ECL_NONCENTRAL,
            };
            let matches_kind = kind_mask == 0 || (retflag & kind_mask) != 0;
            let matches_cent = cent_mask == 0 || (retflag & cent_mask) != 0;
            if matches_kind && matches_cent {
                let dt_days = info.delta_t_sec / SECONDS_PER_DAY;
                let mut tret = [0.0; 10];
                tret[0] = info.max_eclipse_jd_ut1;
                tret[1] = info.t_noon_jd_tt.map_or(0.0, |t| t - dt_days);
                tret[2] = info.p1_jd_tt.map_or(0.0, |t| t - dt_days);
                tret[3] = info.p4_jd_tt.map_or(0.0, |t| t - dt_days);
                if info.kind != SolarEclipseKind::Partial {
                    tret[4] = info
                        .u1_jd_tt
                        .or(info.u2_jd_tt)
                        .map_or(0.0, |t| t - dt_days);
                    tret[5] = info
                        .u4_jd_tt
                        .or(info.u3_jd_tt)
                        .map_or(0.0, |t| t - dt_days);
                    tret[6] = info
                        .c_begin_jd_tt
                        .or(info.u2_jd_tt)
                        .or(info.u1_jd_tt)
                        .map_or(0.0, |t| t - dt_days);
                    tret[7] = info
                        .c_end_jd_tt
                        .or(info.u3_jd_tt)
                        .or(info.u4_jd_tt)
                        .map_or(0.0, |t| t - dt_days);
                }
                return Ok((retflag, tret));
            }
            cursor = info.max_eclipse_jd_ut1 + if backward { -10.0 } else { 10.0 };
        }
        Err(EphemerisError::NonConvergence {
            solver: "sol_eclipse_when_glob".to_owned(),
            iterations: 500,
            residual: 0.0,
        })
    }

    /// Evaluates instantaneous geographic central-line position and attributes for a Solar Eclipse at `jd_ut`,
    /// returning `(retflag, geopos[0..10], attr[0..20])` matching `swe.sol_eclipse_where` (`retflag == 0` if no eclipse at `jd_ut`).
    pub fn sol_eclipse_where(
        &self,
        jd_ut: f64,
        flags: u32,
    ) -> Result<(i32, [f64; 10], [f64; 20]), EphemerisError> {
        let use_moshier = (flags & SEFLG_MOSEPH) != 0;
        let Some((kind, is_central, geopos, attr)) = self.engine.solar_eclipse_where_at_ut1(
            jd_ut,
            EclipseOptions::default(),
            use_moshier,
        )?
        else {
            return Ok((0, [0.0; 10], [0.0; 20]));
        };
        let cent_bit = if is_central {
            SE_ECL_CENTRAL
        } else {
            SE_ECL_NONCENTRAL
        };
        let retflag = match kind {
            SolarEclipseKind::Total => SE_ECL_TOTAL | cent_bit,
            SolarEclipseKind::Annular => SE_ECL_ANNULAR | cent_bit,
            SolarEclipseKind::Hybrid => SE_ECL_ANNULAR_TOTAL | cent_bit,
            SolarEclipseKind::Partial => SE_ECL_PARTIAL | SE_ECL_NONCENTRAL,
        };
        Ok((retflag, geopos, attr))
    }

    /// Searches forward (`backward == false`) or backward (`backward == true`) from `jd_ut_start`
    /// for the next locally visible Solar Eclipse at `geopos = [lon_deg, lat_deg, alt_m]`,
    /// returning `(retflag, tret[0..10], attr[0..20])` matching `swe.sol_eclipse_when_loc`.
    pub fn sol_eclipse_when_loc(
        &self,
        jd_ut_start: f64,
        geopos: [f64; 3],
        flags: u32,
        backward: bool,
    ) -> Result<(i32, [f64; 10], [f64; 20]), EphemerisError> {
        let obs = Observer::new(geopos[1], geopos[0], geopos[2])?;
        let info = match self.engine.next_local_solar_eclipse(
            TimeInput::from_ut1_jd(jd_ut_start),
            obs,
            backward,
            EclipseOptions::default(),
        ) {
            Ok(inf) => inf,
            Err(_) => self.engine.next_solar_eclipse(
                TimeInput::from_ut1_jd(jd_ut_start),
                backward,
                Some(obs),
                EclipseOptions::default(),
            )?,
        };
        let dt_days = info.delta_t_sec / SECONDS_PER_DAY;
        let mut tret = [0.0; 10];
        let mut attr = [0.0; 20];
        let retflag = if let Some(ref loc) = info.local {
            tret[0] = loc.max_eclipse_jd_ut1;
            tret[1] = loc.c1_jd_tt.map_or(0.0, |t| t - dt_days);
            tret[2] = loc.c2_jd_tt.map_or(0.0, |t| t - dt_days);
            tret[3] = loc.c3_jd_tt.map_or(0.0, |t| t - dt_days);
            tret[4] = loc.c4_jd_tt.map_or(0.0, |t| t - dt_days);
            let mut rf = match loc.local_kind {
                SolarEclipseKind::Total => SE_ECL_TOTAL,
                SolarEclipseKind::Annular => SE_ECL_ANNULAR,
                SolarEclipseKind::Hybrid => SE_ECL_ANNULAR_TOTAL,
                SolarEclipseKind::Partial => SE_ECL_PARTIAL,
            };
            let is_vis_sun = |t_ut: f64| -> bool {
                if t_ut <= 0.0 {
                    return false;
                }
                self.engine
                    .calc_azalt(
                        BodyId::Planet(PlanetId::Sun),
                        TimeInput::from_ut1_jd(t_ut),
                        obs,
                        AtmosphereParams::default(),
                    )
                    .is_ok_and(|hz| hz.altitude_apparent_deg >= -0.266)
            };
            let vis_max = is_vis_sun(tret[0]);
            let vis_c1 = is_vis_sun(tret[1]);
            let vis_c2 = is_vis_sun(tret[2]);
            let vis_c3 = is_vis_sun(tret[3]);
            let vis_c4 = is_vis_sun(tret[4]);
            if vis_max {
                rf |= SE_ECL_VISIBLE | SE_ECL_MAX_VISIBLE;
            }
            if vis_c1 {
                rf |= SE_ECL_VISIBLE | SE_ECL_1ST_VISIBLE;
            }
            if vis_c2 {
                rf |= SE_ECL_VISIBLE | SE_ECL_2ND_VISIBLE;
            }
            if vis_c3 {
                rf |= SE_ECL_VISIBLE | SE_ECL_3RD_VISIBLE;
            }
            if vis_c4 {
                rf |= SE_ECL_VISIBLE | SE_ECL_4TH_VISIBLE;
            }
            if !vis_c1 && vis_c4 && tret[1] > 0.0 {
                if let Ok(Some(t_rise)) = self.rise_trans(
                    tret[1] - 0.05,
                    SE_SUN,
                    None,
                    flags,
                    SE_CALC_RISE,
                    geopos,
                    0.0,
                    15.0,
                ) {
                    if t_rise >= tret[1] - 0.01 && t_rise <= tret[4] + 0.01 {
                        tret[5] = t_rise;
                        if !vis_max {
                            tret[0] = t_rise;
                        }
                    }
                }
            }
            if vis_c1 && !vis_c4 && tret[1] > 0.0 {
                if let Ok(Some(t_set)) = self.rise_trans(
                    tret[1] - 0.05,
                    SE_SUN,
                    None,
                    flags,
                    SE_CALC_SET,
                    geopos,
                    0.0,
                    15.0,
                ) {
                    if t_set >= tret[1] - 0.01 && t_set <= tret[4] + 0.01 {
                        tret[6] = t_set;
                        if !vis_max {
                            tret[0] = t_set;
                        }
                    }
                }
            }
            if tret[5] > 0.0 || tret[6] > 0.0 || info.central_latitude_deg.is_none() {
                rf |= SE_ECL_NONCENTRAL;
            }
            attr[0] = loc.magnitude;
            attr[2] = loc.obscuration;
            attr[4] = (loc.sun_azimuth_deg + 180.0).rem_euclid(360.0);
            attr[5] = loc.sun_altitude_deg;
            attr[6] = loc.sun_altitude_deg;
            attr[8] = loc.magnitude;
            attr[9] = f64::from(info.saros_series);
            attr[10] = f64::from(info.saros_member);
            let t_eval = TimeInput::from_ut1_jd(tret[0]);
            if let Ok(hz_sun) = self.engine.calc_azalt(
                BodyId::Planet(PlanetId::Sun),
                t_eval,
                obs,
                AtmosphereParams::default(),
            ) {
                attr[4] = (hz_sun.azimuth_deg + 180.0).rem_euclid(360.0);
                attr[5] = hz_sun.altitude_true_deg;
                attr[6] = hz_sun.altitude_apparent_deg;
            }
            let topo_opts = CalcOptions {
                origin: Origin::Topocentric(obs),
                plane: Plane::Equatorial,
                equinox: Equinox::OfDate,
                nutation: NutationMode::True,
                aberration: AbMode::Apparent,
                zodiac: Zodiac::Tropical,
                include_speed: false,
            };
            let res_s = self
                .engine
                .calc(BodyId::Planet(PlanetId::Sun), t_eval, topo_opts)
                .or_else(|_| {
                    self.engine
                        .calc_moshier(BodyId::Planet(PlanetId::Sun), t_eval, topo_opts)
                });
            let res_m = self
                .engine
                .calc(BodyId::Planet(PlanetId::Moon), t_eval, topo_opts)
                .or_else(|_| {
                    self.engine
                        .calc_moshier(BodyId::Planet(PlanetId::Moon), t_eval, topo_opts)
                });
            if let (Ok(rs), Ok(rm)) = (res_s, res_m) {
                let d_sun_km = (rs.dist_au * crate::types::AU_KM).max(1000.0);
                let d_moon_km = (rm.dist_au * crate::types::AU_KM).max(1000.0);
                let s_sun = (696_000.0 / d_sun_km).clamp(-1.0, 1.0).asin().to_degrees();
                let s_moon = (1737.4 / d_moon_km).clamp(-1.0, 1.0).asin().to_degrees();
                let dra = ((rm.lon_deg - rs.lon_deg + 180.0).rem_euclid(360.0) - 180.0)
                    * rs.lat_deg.to_radians().cos();
                let ddec = rm.lat_deg - rs.lat_deg;
                let sep_deg = dra.hypot(ddec);
                attr[7] = sep_deg;
                if s_sun > 0.0 {
                    let ratio = s_moon / s_sun;
                    attr[1] = ratio;
                    attr[3] = -2.0 * d_moon_km * (s_moon - s_sun).to_radians().sin();
                    if !vis_max {
                        let mag = ((s_sun + s_moon - sep_deg) / (2.0 * s_sun)).max(0.0);
                        attr[0] = mag;
                        attr[8] = mag;
                        let r1 = s_sun;
                        let r2 = s_moon;
                        let d = sep_deg;
                        attr[2] = if d >= r1 + r2 {
                            0.0
                        } else if d <= (r1 - r2).abs() {
                            (r2 / r1).powi(2)
                        } else {
                            let a1 = ((d * d + r1 * r1 - r2 * r2) / (2.0 * d * r1))
                                .clamp(-1.0, 1.0)
                                .acos();
                            let a2 = ((d * d + r2 * r2 - r1 * r1) / (2.0 * d * r2))
                                .clamp(-1.0, 1.0)
                                .acos();
                            let term = (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2);
                            let area = r1 * r1 * a1 + r2 * r2 * a2 - 0.5 * term.max(0.0).sqrt();
                            (area / (std::f64::consts::PI * r1 * r1)).clamp(0.0, 1.0)
                        };
                    } else if loc.local_kind != SolarEclipseKind::Partial {
                        attr[2] = ratio * ratio;
                        attr[8] = ratio;
                    }
                }
            }
            rf
        } else {
            tret[0] = info.max_eclipse_jd_ut1;
            attr[9] = f64::from(info.saros_series);
            attr[10] = f64::from(info.saros_member);
            0
        };
        Ok((retflag, tret, attr))
    }

    /// Searches forward (`backward == false`) or backward (`backward == true`) from `jd_ut_start`
    /// for the next Lunar Eclipse matching `ifltype` (`0` = any), returning `(retflag, tret[0..10])`
    /// matching `swe.lun_eclipse_when`.
    pub fn lun_eclipse_when(
        &self,
        jd_ut_start: f64,
        _flags: u32,
        ifltype: i32,
        backward: bool,
    ) -> Result<(i32, [f64; 10]), EphemerisError> {
        let kind_mask = ifltype & (SE_ECL_TOTAL | SE_ECL_PARTIAL | SE_ECL_PENUMBRAL);
        let mut cursor = jd_ut_start;
        for _ in 0..500 {
            let info = self.engine.next_lunar_eclipse(
                TimeInput::from_ut1_jd(cursor),
                backward,
                EclipseOptions::default(),
            )?;
            let retflag = match info.kind {
                LunarEclipseKind::Total => SE_ECL_TOTAL,
                LunarEclipseKind::Partial => SE_ECL_PARTIAL,
                LunarEclipseKind::Penumbral => SE_ECL_PENUMBRAL,
            };
            if kind_mask == 0 || (retflag & kind_mask) != 0 {
                let dt_days = info.delta_t_sec / SECONDS_PER_DAY;
                let mut tret = [0.0; 10];
                tret[0] = info.max_eclipse_jd_ut1;
                tret[2] = info.u1_jd_tt.map_or(0.0, |t| t - dt_days);
                tret[3] = info.u4_jd_tt.map_or(0.0, |t| t - dt_days);
                tret[4] = info.u2_jd_tt.map_or(0.0, |t| t - dt_days);
                tret[5] = info.u3_jd_tt.map_or(0.0, |t| t - dt_days);
                tret[6] = info.p1_jd_tt.map_or(0.0, |t| t - dt_days);
                tret[7] = info.p4_jd_tt.map_or(0.0, |t| t - dt_days);
                return Ok((retflag, tret));
            }
            cursor = info.max_eclipse_jd_ut1 + if backward { -10.0 } else { 10.0 };
        }
        Err(EphemerisError::NonConvergence {
            solver: "lun_eclipse_when".to_owned(),
            iterations: 500,
            residual: 0.0,
        })
    }

    /// Evaluates instantaneous magnitudes and Saros metadata for a Lunar Eclipse at `jd_ut`,
    /// returning `(retflag, attr[0..20])` matching `swe.lun_eclipse_how` (`retflag == 0` if no eclipse at `jd_ut`).
    pub fn lun_eclipse_how(
        &self,
        jd_ut: f64,
        geopos: [f64; 3],
        flags: u32,
    ) -> Result<(i32, [f64; 20]), EphemerisError> {
        let use_moshier = (flags & SEFLG_MOSEPH) != 0;
        let Some((kind, mut attr)) = self.engine.lunar_eclipse_how_at_ut1(
            jd_ut,
            EclipseOptions::default(),
            use_moshier,
        )?
        else {
            return Ok((0, [0.0; 20]));
        };
        let mut retflag = match kind {
            LunarEclipseKind::Total => SE_ECL_TOTAL,
            LunarEclipseKind::Partial => SE_ECL_PARTIAL,
            LunarEclipseKind::Penumbral => SE_ECL_PENUMBRAL,
        };
        if let Ok(obs) = Observer::new(geopos[1], geopos[0], geopos[2]) {
            if let Ok(hz_moon) = self.engine.calc_azalt(
                BodyId::Planet(PlanetId::Moon),
                TimeInput::from_ut1_jd(jd_ut),
                obs,
                AtmosphereParams::default(),
            ) {
                attr[4] = (hz_moon.azimuth_deg + 180.0).rem_euclid(360.0);
                attr[5] = hz_moon.altitude_true_deg;
                attr[6] = hz_moon.altitude_apparent_deg;
                if hz_moon.altitude_apparent_deg < -0.5 {
                    retflag = 0;
                }
            }
        }
        Ok((retflag, attr))
    }

    /// Searches for the next heliacal visibility event (`event_type`: `1 = SE_HELIACAL_RISING`,
    /// `2 = SE_HELIACAL_SETTING`, `3 = SE_EVENING_FIRST`, `4 = SE_MORNING_LAST`) of `object_name`
    /// starting at `jd_ut_start` for `geopos = [lon_deg, lat_deg, alt_m]`, returning
    /// `([trise_opt, trise_best, trise_end], retflag)` matching `swe.heliacal_ut`.
    #[allow(clippy::too_many_arguments)]
    pub fn heliacal_ut(
        &self,
        jd_ut_start: f64,
        geopos: [f64; 3],
        atmo: [f64; 4],
        obs_arr: [f64; 6],
        object_name: &str,
        event_type: i32,
        helflag: i32,
    ) -> Result<([f64; 3], i32), EphemerisError> {
        let body = if let Ok(meta) = self.engine.star_info(object_name) {
            BodyId::FixedStar(crate::types::HipId(meta.hip))
        } else {
            object_name.trim().parse::<BodyId>()?
        };
        let obs = Observer::new(geopos[1], geopos[0], geopos[2])?;
        let pressure_hpa = if atmo[0] > 0.0 { atmo[0] } else { 1013.25 };
        let temperature_c = if atmo[0] > 0.0 || atmo[1] != 0.0 {
            atmo[1]
        } else {
            15.0
        };
        let relative_humidity = if atmo[2] > 1.0 {
            (atmo[2] / 100.0).clamp(0.0, 1.0)
        } else if atmo[2] > 0.0 {
            atmo[2]
        } else {
            0.5
        };
        // In Swiss Ephemeris swe_heliacal_ut, dattm[3] >= 1.0 is meteorological visibility range V_m in km
        // (Koschmieder relation k_V ≈ 0.14 + 3.0 / V_m), whereas 0 < dattm[3] < 1.0 is direct k_V (mag/airmass).
        let k_v = if atmo[3] >= 1.0 {
            (0.14 + 3.0 / atmo[3]).clamp(0.10, 1.0)
        } else if atmo[3] > 0.0 {
            atmo[3]
        } else {
            0.20
        };
        let observer_age_years = if obs_arr[0] > 0.0 { obs_arr[0] } else { 36.0 };
        let snellen_ratio = if obs_arr[1] > 0.0 { obs_arr[1] } else { 1.0 };
        let hel_opts = HeliacalOptions {
            k_v,
            observer_age_years,
            snellen_ratio,
            relative_humidity,
            temperature_c,
            pressure_hpa,
            max_search_days: 1095.75,
        };
        if body == BodyId::Planet(PlanetId::Moon)
            && matches!(event_type, SE_HELIACAL_RISING | SE_HELIACAL_SETTING)
        {
            return Err(EphemerisError::InvalidInput {
                message: format!(
                    "heliacal_ut event_type {event_type} is not valid for Moon (use 3=SE_EVENING_FIRST or 4=SE_MORNING_LAST)"
                ),
            });
        }
        if matches!(event_type, SE_EVENING_FIRST | SE_MORNING_LAST)
            && matches!(
                body,
                BodyId::Planet(
                    PlanetId::Mars
                        | PlanetId::Jupiter
                        | PlanetId::Saturn
                        | PlanetId::Uranus
                        | PlanetId::Neptune
                        | PlanetId::Pluto
                )
            )
        {
            return Err(EphemerisError::InvalidInput {
                message: format!(
                    "heliacal_ut event_type {event_type} is only defined for inner bodies (Moon, Mercury, Venus), not '{object_name}'"
                ),
            });
        }
        let flags_u32 = helflag as u32;
        let use_moshier = (flags_u32 & SEFLG_MOSEPH != 0) || self.engine.is_moshier_mode();
        let (hel, used_moshier) = if use_moshier {
            (
                self.engine
                    .heliacal_phenomena_moshier(body, jd_ut_start, obs, hel_opts)?,
                true,
            )
        } else {
            match self
                .engine
                .heliacal_phenomena(body, jd_ut_start, obs, hel_opts)
            {
                Ok(h) => (h, false),
                Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                    match self
                        .engine
                        .heliacal_phenomena_moshier(body, jd_ut_start, obs, hel_opts)
                    {
                        Ok(h) => (h, true),
                        Err(_) => return Err(orig_err),
                    }
                }
                Err(e) => return Err(e),
            }
        };
        let jd_opt = match event_type {
            SE_HELIACAL_RISING => hel.heliacal_rise_jd_ut1,
            SE_HELIACAL_SETTING => hel.heliacal_set_jd_ut1,
            SE_EVENING_FIRST => hel.evening_rise_jd_ut1,
            SE_MORNING_LAST => hel.morning_set_jd_ut1,
            other => {
                return Err(EphemerisError::InvalidInput {
                    message: format!("Unknown heliacal event_type {other} (expected 1..=4)"),
                });
            }
        };
        let jd_ev = jd_opt.ok_or_else(|| EphemerisError::NonConvergence {
            solver: "heliacal_ut".to_owned(),
            iterations: 548,
            residual: 0.0,
        })?;
        // Physical twilight visibility window [t_start, t_opt, t_end] bounded by object horizon crossing
        // `jd_ev` and the actual topocentric sunrise/sunset on that date.
        let is_morning = event_type == SE_HELIACAL_RISING || event_type == SE_MORNING_LAST;
        let fallback_dt = (hel.arcus_visionis_deg * 0.5 / 360.0).clamp(0.004, 0.025);
        let dt_twilight = if is_morning {
            self.rise_trans(
                jd_ev - 0.1,
                SE_SUN,
                None,
                flags_u32,
                SE_CALC_RISE | SE_BIT_DISC_CENTER,
                geopos,
                pressure_hpa,
                temperature_c,
            )
            .ok()
            .flatten()
            .filter(|&t_sr| t_sr > jd_ev && t_sr - jd_ev < 0.25)
            .map_or(fallback_dt, |t_sr| (0.65 * (t_sr - jd_ev)).max(0.002))
        } else {
            self.rise_trans(
                jd_ev - 0.25,
                SE_SUN,
                None,
                flags_u32,
                SE_CALC_SET | SE_BIT_DISC_CENTER,
                geopos,
                pressure_hpa,
                temperature_c,
            )
            .ok()
            .flatten()
            .filter(|&t_ss| t_ss < jd_ev && jd_ev - t_ss < 0.25)
            .map_or(fallback_dt, |t_ss| (0.65 * (jd_ev - t_ss)).max(0.002))
        };
        let dret = if is_morning {
            [jd_ev, jd_ev + 0.5 * dt_twilight, jd_ev + dt_twilight]
        } else {
            [jd_ev - dt_twilight, jd_ev - 0.5 * dt_twilight, jd_ev]
        };
        let retflag = if used_moshier {
            SEFLG_MOSEPH as i32
        } else {
            SEFLG_SWIEPH as i32
        };
        Ok((dret, retflag))
    }
}
