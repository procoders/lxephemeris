//! Complete Node.js / TypeScript bindings (`napi-rs`) for `lxephemeris`.
//!
//! Provides `JsEngine`, `JsSweAdapter`, `JsSplitJd`, `JsObserver`, `JsCalcOptions`, `JsCalcResult`,
//! `JsBatchItemResult`, `JsDeltaTInfo`, `JsSiderealTimeInfo`, `JsEclNutInfo`, `JsStarMetadata`,
//! `JsHeliacalResult`, `JsAyanamsaValue`, `JsHouseResult`, `JsHorizontalCoords`,
//! `JsSolarEclipseInfo`, `JsLunarEclipseInfo`, `JsCalendarDate`, `JsLongitudeCrossingEvent`,
//! `JsSolarTermEvent`, `JsSeasonsOfYear`, `JsLunarPhaseEvent`, `JsLunarPhaseSummary`,
//! `JsStationEvent`, `JsDeclinationExtremumEvent`, `JsParanEvent`, `julday`, `revjul`,
//! `housesArmc`, `houseOfLon`, `starInfo`, `listStars`, `gauquelinSectors`, flat `Float64Array`-compatible
//! batch buffers, and structured `[ErrorCode]` errors.

#![cfg_attr(not(feature = "napi-addon"), forbid(unsafe_code))]

use lxephemeris_core::{
    azalt as core_azalt, body_display_name, find_star_by_name_or_bayer,
    gauquelin_sectors as core_gauquelin_sectors, house_of_lon as core_house_of_lon,
    houses_armc_with_sidereal, julday as core_julday, load_custom_stars_file,
    map_zodiac_placement as core_map_zodiac_placement, revjul as core_revjul, AbMode,
    AspectPreset, AtmosphereParams, AyanamsaModel, AyanamsaValue, BodyId, CalcOptions, CalcResult,
    CalendarKind, DeclinationExtremumEvent, DeltaTInfo, DiscMode, EclNutInfo, EclipseOptions,
    Engine, EngineConfig, EphemerisError, Equinox, HeliacalOptions, HeliacalResult, HipId,
    HorizontalCoords, HouseOptions, HouseSystem, HousesResult, HypotheticalId,
    LongitudeCrossingEvent, LunarEclipseInfo, LunarOccultationInfo, LunarPhaseEvent,
    LunarPhaseSummary, LunarPointId, MpcNumber, NodesApsidesResult, NutationMode, Observer, Origin,
    ParanEvent, ParanOptions, PhenomenaResult, Plane, PlanetId, PolarPolicy, RiseTransEvent,
    RiseTransOptions, SeasonsOfYear, ShadowEnlargementRule, SiderealTimeInfo, SolarEclipseInfo,
    SolarTermEvent, SplitJd, StarMetadata, StationEvent, TimeInput, Zodiac, MINOR_BODIES_68,
};
#[cfg(feature = "swe-compat")]
use lxephemeris_core::SweAdapter;
#[cfg(feature = "napi-addon")]
use napi_derive::napi;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::sync::Arc;
#[cfg(feature = "swe-compat")]
use std::sync::Mutex;

/// Parses a string body identifier into a [`BodyId`].
pub fn parse_body_id(s: &str) -> Result<BodyId, EphemerisError> {
    let trimmed = s.trim();
    match trimmed.to_ascii_lowercase().as_str() {
        "sun" => Ok(BodyId::Planet(PlanetId::Sun)),
        "moon" => Ok(BodyId::Planet(PlanetId::Moon)),
        "mercury" => Ok(BodyId::Planet(PlanetId::Mercury)),
        "venus" => Ok(BodyId::Planet(PlanetId::Venus)),
        "earth" => Ok(BodyId::Planet(PlanetId::Earth)),
        "mars" => Ok(BodyId::Planet(PlanetId::Mars)),
        "jupiter" => Ok(BodyId::Planet(PlanetId::Jupiter)),
        "saturn" => Ok(BodyId::Planet(PlanetId::Saturn)),
        "uranus" => Ok(BodyId::Planet(PlanetId::Uranus)),
        "neptune" => Ok(BodyId::Planet(PlanetId::Neptune)),
        "pluto" => Ok(BodyId::Planet(PlanetId::Pluto)),
        "ecl_nut" | "eclnut" => Ok(BodyId::EclNut),
        "mean_node" | "meannode" => Ok(BodyId::LunarPoint(LunarPointId::MeanNode)),
        "true_node" | "truenode" => Ok(BodyId::LunarPoint(LunarPointId::TrueNode)),
        "mean_south_node" | "meansouthnode" => Ok(BodyId::LunarPoint(LunarPointId::MeanSouthNode)),
        "true_south_node" | "truesouthnode" => Ok(BodyId::LunarPoint(LunarPointId::TrueSouthNode)),
        "mean_lilith" | "meanlilith" | "mean_apogee" => {
            Ok(BodyId::LunarPoint(LunarPointId::MeanLilith))
        }
        "true_lilith" | "truelilith" | "osc_apogee" => {
            Ok(BodyId::LunarPoint(LunarPointId::TrueLilith))
        }
        "interpolated_lilith" | "interpolatedlilith" | "intp_apogee" | "intp_apog" => {
            Ok(BodyId::LunarPoint(LunarPointId::InterpolatedLilith))
        }
        "interpolated_perigee" | "interpolatedperigee" | "intp_perigee" | "intp_perg" => {
            Ok(BodyId::LunarPoint(LunarPointId::InterpolatedPerigee))
        }
        "cupido" => Ok(BodyId::Hypothetical(HypotheticalId::Cupido)),
        "hades" => Ok(BodyId::Hypothetical(HypotheticalId::Hades)),
        "zeus" => Ok(BodyId::Hypothetical(HypotheticalId::Zeus)),
        "kronos" => Ok(BodyId::Hypothetical(HypotheticalId::Kronos)),
        "apollon" => Ok(BodyId::Hypothetical(HypotheticalId::Apollon)),
        "admetos" => Ok(BodyId::Hypothetical(HypotheticalId::Admetos)),
        "vulkanus" => Ok(BodyId::Hypothetical(HypotheticalId::Vulkanus)),
        "poseidon" => Ok(BodyId::Hypothetical(HypotheticalId::Poseidon)),
        "transpluto" | "isis" => Ok(BodyId::Hypothetical(HypotheticalId::Transpluto)),
        other => {
            if let Some(mpc_str) = other.strip_prefix("mpc:") {
                let mpc: u32 = mpc_str
                    .parse()
                    .map_err(|_| EphemerisError::UnsupportedBody {
                        body: trimmed.to_owned(),
                    })?;
                return Ok(BodyId::MinorPlanet(MpcNumber(mpc)));
            }
            if let Some(hip_str) = other.strip_prefix("hip:") {
                let hip: u32 = hip_str
                    .parse()
                    .map_err(|_| EphemerisError::UnsupportedBody {
                        body: trimmed.to_owned(),
                    })?;
                return Ok(BodyId::FixedStar(HipId(hip)));
            }
            if let Some(spec) = MINOR_BODIES_68
                .iter()
                .find(|m| m.name.eq_ignore_ascii_case(trimmed))
            {
                return Ok(BodyId::MinorPlanet(spec.mpc));
            }
            if let Some(star) = find_star_by_name_or_bayer(trimmed) {
                return Ok(BodyId::FixedStar(star.hip_id()));
            }
            Err(EphemerisError::UnsupportedBody {
                body: trimmed.to_owned(),
            })
        }
    }
}

/// Parses a time scale string (`"TT"`, `"UT1"`, `"UTC"`, `"TDB"`) and [`JsSplitJd`] into [`TimeInput`].
pub fn parse_time_input(scale: &str, jd: JsSplitJd) -> Result<TimeInput, EphemerisError> {
    let jd2 = jd.jd2.unwrap_or(0.0);
    if !jd.jd1.is_finite() || !jd2.is_finite() {
        return Err(EphemerisError::InvalidInput {
            message: format!("Non-finite Julian Date: jd1={}, jd2={jd2}", jd.jd1),
        });
    }
    let split = SplitJd::new(jd.jd1, jd2);
    match scale.trim().to_ascii_uppercase().as_str() {
        "TT" => Ok(TimeInput::Tt(split)),
        "UT1" | "UT" => Ok(TimeInput::Ut1(split)),
        "UTC" => Ok(TimeInput::Utc(split)),
        "TDB" | "ET" => Ok(TimeInput::Tdb(split)),
        other => Err(EphemerisError::InvalidInput {
            message: format!("Unknown time scale '{other}' (expected TT, UT1, UTC, or TDB)"),
        }),
    }
}

fn parse_calendar_kind(calendar: &str) -> Result<CalendarKind, EphemerisError> {
    match calendar.trim().to_ascii_lowercase().as_str() {
        "gregorian" | "greg" | "g" => Ok(CalendarKind::Gregorian),
        "julian" | "jul" | "j" => Ok(CalendarKind::Julian),
        "auto" | "papal" | "astro_auto" => Ok(CalendarKind::AstroAuto),
        other => Err(EphemerisError::InvalidInput {
            message: format!("Unknown calendar '{other}'"),
        }),
    }
}

fn parse_polar_policy(s: &str) -> Result<PolarPolicy, EphemerisError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "strict" => Ok(PolarPolicy::Strict),
        "fallback" | "fallback_to_porphyry" | "porphyry" => Ok(PolarPolicy::FallbackToPorphyry),
        "fallback_whole_sign" | "fallback_to_whole_sign" | "whole_sign" | "wholesign" => {
            Ok(PolarPolicy::FallbackToWholeSign)
        }
        other => Err(EphemerisError::InvalidInput {
            message: format!("Unknown polar_policy '{other}'"),
        }),
    }
}

fn parse_disc_mode(s: &str) -> Result<DiscMode, EphemerisError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "upper_limb" | "upper" | "limb" => Ok(DiscMode::UpperLimb),
        "center" | "disc_center" => Ok(DiscMode::Center),
        "bottom_limb" | "bottom" | "disc_bottom" => Ok(DiscMode::BottomLimb),
        other => Err(EphemerisError::InvalidInput {
            message: format!("Unknown disc_mode '{other}'"),
        }),
    }
}

fn parse_zodiac(zodiac: &str, ayanamsa: &str) -> Result<Zodiac, EphemerisError> {
    match zodiac.trim().to_ascii_lowercase().as_str() {
        "tropical" => Ok(Zodiac::Tropical),
        "sidereal" => {
            let m = AyanamsaModel::from_code(ayanamsa)?;
            Ok(Zodiac::Sidereal(m))
        }
        other => {
            let m = AyanamsaModel::from_code(other)?;
            Ok(Zodiac::Sidereal(m))
        }
    }
}

fn parse_nutation_mode(s: &str) -> Result<NutationMode, EphemerisError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "true" => Ok(NutationMode::True),
        "mean" | "none" | "false" => Ok(NutationMode::Mean),
        other => Err(EphemerisError::InvalidInput {
            message: format!("Unknown nutation mode '{other}' (expected 'true' or 'mean')"),
        }),
    }
}

fn parse_shadow_rule(s: &str) -> Result<ShadowEnlargementRule, EphemerisError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "danjon" | "danjon_1_over_85" | "1/85" => Ok(ShadowEnlargementRule::Danjon1Over85),
        "traditional" | "traditional_2_percent" | "2%" => {
            Ok(ShadowEnlargementRule::Traditional2Percent)
        }
        other => Err(EphemerisError::InvalidInput {
            message: format!("Unknown shadow_rule '{other}' (expected 'danjon' or 'traditional')"),
        }),
    }
}

fn parse_aspect_preset(s: Option<&str>) -> Result<AspectPreset, EphemerisError> {
    match s.unwrap_or("major").trim().to_ascii_lowercase().as_str() {
        "major" | "ptolemaic" => Ok(AspectPreset::Major),
        "standard" | "minor" => Ok(AspectPreset::Standard),
        "kepler" | "harmonic" => Ok(AspectPreset::Kepler),
        "all" => Ok(AspectPreset::All),
        other => Err(EphemerisError::InvalidInput {
            message: format!(
                "Unknown aspect preset '{other}' (expected 'major', 'standard', 'kepler', or 'all')"
            ),
        }),
    }
}

fn default_chart_bodies() -> Vec<BodyId> {
    vec![
        BodyId::Planet(PlanetId::Sun),
        BodyId::Planet(PlanetId::Moon),
        BodyId::Planet(PlanetId::Mercury),
        BodyId::Planet(PlanetId::Venus),
        BodyId::Planet(PlanetId::Mars),
        BodyId::Planet(PlanetId::Jupiter),
        BodyId::Planet(PlanetId::Saturn),
        BodyId::Planet(PlanetId::Uranus),
        BodyId::Planet(PlanetId::Neptune),
        BodyId::Planet(PlanetId::Pluto),
        BodyId::LunarPoint(LunarPointId::MeanNode),
        BodyId::LunarPoint(LunarPointId::TrueNode),
        BodyId::LunarPoint(LunarPointId::MeanLilith),
    ]
}

fn parse_body_list_or_default(bodies: Option<&[String]>) -> Result<Vec<BodyId>, EphemerisError> {
    match bodies {
        Some(list) => list.iter().map(|b| parse_body_id(b)).collect(),
        None => Ok(default_chart_bodies()),
    }
}

fn snake_to_camel(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut capitalize_next = false;
    for ch in s.chars() {
        if ch == '_' {
            capitalize_next = true;
        } else if capitalize_next {
            out.extend(ch.to_uppercase());
            capitalize_next = false;
        } else {
            out.push(ch);
        }
    }
    out
}

fn enrich_json_value(val: serde_json::Value) -> serde_json::Value {
    match val {
        serde_json::Value::Array(arr) => {
            serde_json::Value::Array(arr.into_iter().map(enrich_json_value).collect())
        }
        serde_json::Value::Object(mut map) => {
            for (k, v) in map.iter_mut() {
                if matches!(k.as_str(), "body" | "body1" | "body2" | "last_aspect_body") {
                    if let Ok(bid) = serde_json::from_value::<BodyId>(v.clone()) {
                        *v = serde_json::Value::String(body_display_name(bid));
                        continue;
                    }
                }
                let taken = std::mem::replace(v, serde_json::Value::Null);
                *v = enrich_json_value(taken);
            }
            let mut camel_additions = Vec::new();
            for (k, v) in &map {
                if k.contains('_') {
                    let camel = snake_to_camel(k);
                    if !map.contains_key(&camel) {
                        camel_additions.push((camel, v.clone()));
                    }
                }
            }
            for (k, v) in camel_additions {
                map.insert(k, v);
            }
            serde_json::Value::Object(map)
        }
        other => other,
    }
}

fn to_js_value<T: Serialize>(item: &T) -> Result<serde_json::Value, String> {
    let v = serde_json::to_value(item).map_err(|e| e.to_string())?;
    Ok(enrich_json_value(v))
}

fn to_js_value_vec<T: Serialize>(items: &[T]) -> Result<Vec<serde_json::Value>, String> {
    items.iter().map(to_js_value).collect()
}

/// Formats an [`EphemerisError`] as a structured error string `"[ErrorCode] message"`.
#[must_use]
pub fn format_js_error(err: &EphemerisError) -> String {
    format!("[{}] {err}", err.code().as_str())
}

/// Two-part Julian Date input for Node.js / TypeScript.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsSplitJd {
    pub jd1: f64,
    pub jd2: Option<f64>,
}

/// Geodetic WGS84 observer for Node.js / TypeScript.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsObserver {
    pub latitude_deg: f64,
    pub longitude_deg: f64,
    pub altitude_m: Option<f64>,
}

/// Calculation options object for Node.js / TypeScript.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct JsCalcOptions {
    pub origin: Option<String>,
    pub plane: Option<String>,
    pub equinox: Option<String>,
    pub nutation: Option<String>,
    pub aberration: Option<String>,
    pub zodiac: Option<String>,
    pub ayanamsa: Option<String>,
    pub include_speed: Option<bool>,
    pub observer: Option<JsObserver>,
}

impl JsCalcOptions {
    pub fn to_core(&self) -> Result<CalcOptions, EphemerisError> {
        let orig = match self
            .origin
            .as_deref()
            .unwrap_or("geocentric")
            .to_ascii_lowercase()
            .as_str()
        {
            "geocentric" => Origin::Geocentric,
            "topocentric" => {
                let obs = self.observer.ok_or_else(|| EphemerisError::InvalidInput {
                    message: "origin='topocentric' requires an observer".to_owned(),
                })?;
                Origin::Topocentric(Observer::new(
                    obs.latitude_deg,
                    obs.longitude_deg,
                    obs.altitude_m.unwrap_or(0.0),
                )?)
            }
            "heliocentric" => Origin::Heliocentric,
            "barycentric" => Origin::Barycentric,
            other => {
                return Err(EphemerisError::InvalidInput {
                    message: format!("Unknown origin '{other}'"),
                })
            }
        };

        let pl = match self
            .plane
            .as_deref()
            .unwrap_or("ecliptic")
            .to_ascii_lowercase()
            .as_str()
        {
            "ecliptic" => Plane::Ecliptic,
            "equatorial" => Plane::Equatorial,
            "horizontal" => Plane::Horizontal,
            other => {
                return Err(EphemerisError::InvalidInput {
                    message: format!("Unknown plane '{other}'"),
                })
            }
        };

        let eq = match self
            .equinox
            .as_deref()
            .unwrap_or("of_date")
            .to_ascii_lowercase()
            .as_str()
        {
            "of_date" | "ofdate" => Equinox::OfDate,
            "j2000" | "icrf" => Equinox::J2000,
            other => {
                return Err(EphemerisError::InvalidInput {
                    message: format!("Unknown equinox '{other}'"),
                })
            }
        };

        let nut = parse_nutation_mode(self.nutation.as_deref().unwrap_or("true"))?;

        let ab_mode = match self
            .aberration
            .as_deref()
            .unwrap_or("apparent")
            .to_ascii_lowercase()
            .as_str()
        {
            "apparent" => AbMode::Apparent,
            "no_aberration" | "noaberr" => AbMode::NoAberration,
            "no_deflection" | "nogdefl" => AbMode::NoDeflection,
            "astrometric" => AbMode::Astrometric,
            "geometric" | "true" | "truepos" => AbMode::Geometric,
            other => {
                return Err(EphemerisError::InvalidInput {
                    message: format!("Unknown aberration '{other}'"),
                })
            }
        };

        let zod = parse_zodiac(
            self.zodiac.as_deref().unwrap_or("tropical"),
            self.ayanamsa.as_deref().unwrap_or("lahiri"),
        )?;

        Ok(CalcOptions {
            origin: orig,
            plane: pl,
            equinox: eq,
            nutation: nut,
            aberration: ab_mode,
            zodiac: zod,
            include_speed: self.include_speed.unwrap_or(true),
        })
    }
}

/// House calculation options object for Node.js / TypeScript.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct JsHouseOptions {
    pub system: Option<String>,
    pub polar_policy: Option<String>,
    pub zodiac: Option<String>,
    pub ayanamsa: Option<String>,
    pub nutation: Option<String>,
    pub include_speed: Option<bool>,
}

impl JsHouseOptions {
    pub fn to_core(&self) -> Result<HouseOptions, EphemerisError> {
        let sys_str = self.system.as_deref().unwrap_or("P");
        let ch = sys_str.chars().next().unwrap_or(' ');
        let sys = HouseSystem::from_code(ch)?;
        let pol = parse_polar_policy(self.polar_policy.as_deref().unwrap_or("strict"))?;
        let zod = parse_zodiac(
            self.zodiac.as_deref().unwrap_or("tropical"),
            self.ayanamsa.as_deref().unwrap_or("lahiri"),
        )?;
        let nut = parse_nutation_mode(self.nutation.as_deref().unwrap_or("true"))?;
        Ok(HouseOptions {
            system: sys,
            polar_policy: pol,
            zodiac: zod,
            nutation: nut,
            include_speed: self.include_speed.unwrap_or(true),
        })
    }
}

/// Single-body astrometric calculation result exposed to Node.js / TypeScript.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsCalcResult {
    pub target_id: i32,
    pub center_id: i32,
    pub jd_tt: f64,
    pub jd_ut1: f64,
    pub jd_tdb: f64,
    pub delta_t_sec: f64,
    pub time_extrapolated: bool,
    pub time_table_version: String,
    pub lon_deg: f64,
    pub lat_deg: f64,
    pub dist_au: f64,
    pub dist_km: f64,
    pub lon_speed_deg_day: f64,
    pub lat_speed_deg_day: f64,
    pub dist_speed_au_day: f64,
    pub dist_speed_km_day: f64,
    pub magnitude_v: Option<f64>,
    pub cartesian_au: Vec<f64>,
    pub cartesian_vel_au_day: Vec<f64>,
    pub warnings: Vec<String>,
}

impl From<CalcResult> for JsCalcResult {
    fn from(r: CalcResult) -> Self {
        Self {
            target_id: r.target_id,
            center_id: r.center_id,
            jd_tt: r.jd_tt,
            jd_ut1: r.jd_ut1,
            jd_tdb: r.jd_tdb,
            delta_t_sec: r.delta_t_sec,
            time_extrapolated: r.time_extrapolated,
            time_table_version: r.time_table_version,
            lon_deg: r.lon_deg,
            lat_deg: r.lat_deg,
            dist_au: r.dist_au,
            dist_km: r.dist_km,
            lon_speed_deg_day: r.lon_speed_deg_day,
            lat_speed_deg_day: r.lat_speed_deg_day,
            dist_speed_au_day: r.dist_speed_au_day,
            dist_speed_km_day: r.dist_speed_km_day,
            magnitude_v: r.magnitude_v,
            cartesian_au: r.cartesian_au.to_vec(),
            cartesian_vel_au_day: r.cartesian_vel_au_day.to_vec(),
            warnings: r.warnings,
        }
    }
}

/// Per-item batch result entry for non-aborting batch queries in Node.js / TypeScript.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsBatchItemResult {
    pub ok: bool,
    pub result: Option<JsCalcResult>,
    pub error_code: Option<String>,
    pub error_message: Option<String>,
}

/// $\Delta T = \text{TT} - \text{UT1}$ and $\text{TT} - \text{TDB}$ result (`R1.2`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsDeltaTInfo {
    pub delta_t_sec: f64,
    pub tt_minus_tdb_sec: f64,
    pub extrapolated: bool,
    pub table_version: String,
    pub warning: Option<String>,
}

impl From<DeltaTInfo> for JsDeltaTInfo {
    fn from(d: DeltaTInfo) -> Self {
        let warning = if d.time_extrapolated {
            Some("TimeExtrapolated".to_owned())
        } else {
            None
        };
        Self {
            delta_t_sec: d.delta_t_sec,
            tt_minus_tdb_sec: d.tt_minus_tdb_sec,
            extrapolated: d.time_extrapolated,
            table_version: d.time_table_version,
            warning,
        }
    }
}

/// Greenwich and local Sidereal Time result (`R1.4`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsSiderealTimeInfo {
    pub gmst_hours: f64,
    pub gast_hours: f64,
    pub local_sidereal_hours: f64,
    pub armc_deg: f64,
}

impl From<SiderealTimeInfo> for JsSiderealTimeInfo {
    fn from(s: SiderealTimeInfo) -> Self {
        Self {
            gmst_hours: s.gmst_hours,
            gast_hours: s.gast_hours,
            local_sidereal_hours: s.local_sidereal_hours,
            armc_deg: s.armc_deg,
        }
    }
}

/// Ecliptic obliquity and IAU 2006/2000A nutation angles (`R2.2`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsEclNutInfo {
    pub true_obliquity_deg: f64,
    pub mean_obliquity_deg: f64,
    pub nutation_lon_deg: f64,
    pub nutation_obl_deg: f64,
}

impl From<EclNutInfo> for JsEclNutInfo {
    fn from(e: EclNutInfo) -> Self {
        Self {
            true_obliquity_deg: e.true_obliquity_deg,
            mean_obliquity_deg: e.mean_obliquity_deg,
            nutation_lon_deg: e.nutation_lon_deg,
            nutation_obl_deg: e.nutation_obl_deg,
        }
    }
}

/// Fixed star catalog metadata (`R6.1`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsStarMetadata {
    pub hip: u32,
    pub canonical_name: String,
    pub bayer: String,
    pub vmag: f64,
    pub ra_rad: f64,
    pub dec_rad: f64,
    pub pm_ra_mas_yr: f64,
    pub pm_dec_mas_yr: f64,
    pub plx_mas: f64,
    pub rv_kms: f64,
    pub rv_missing: bool,
}

impl From<StarMetadata> for JsStarMetadata {
    fn from(s: StarMetadata) -> Self {
        Self {
            hip: s.hip,
            canonical_name: s.canonical_name,
            bayer: s.bayer,
            vmag: s.vmag,
            ra_rad: s.ra_rad,
            dec_rad: s.dec_rad,
            pm_ra_mas_yr: s.pm_ra_mas_yr,
            pm_dec_mas_yr: s.pm_dec_mas_yr,
            plx_mas: s.plx_mas,
            rv_kms: s.rv_kms,
            rv_missing: s.rv_missing,
        }
    }
}

/// Heliacal phenomena result (`R6.3`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsHeliacalResult {
    pub target_name: String,
    pub visual_magnitude: f64,
    pub arcus_visionis_deg: f64,
    pub heliacal_rise_jd_ut1: Option<f64>,
    pub morning_set_jd_ut1: Option<f64>,
    pub evening_rise_jd_ut1: Option<f64>,
    pub heliacal_set_jd_ut1: Option<f64>,
}

impl From<HeliacalResult> for JsHeliacalResult {
    fn from(h: HeliacalResult) -> Self {
        Self {
            target_name: h.target_name,
            visual_magnitude: h.visual_magnitude,
            arcus_visionis_deg: h.arcus_visionis_deg,
            heliacal_rise_jd_ut1: h.heliacal_rise_jd_ut1,
            morning_set_jd_ut1: h.morning_set_jd_ut1,
            evening_rise_jd_ut1: h.evening_rise_jd_ut1,
            heliacal_set_jd_ut1: h.heliacal_set_jd_ut1,
        }
    }
}

/// Vedic / Sidereal ayanamsa evaluation result.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsAyanamsaValue {
    pub model_code: String,
    pub jd_tt: f64,
    pub mean_deg: f64,
    pub true_deg: f64,
    pub effective_deg: f64,
    pub rate_deg_day: f64,
    pub dpsi_arcsec: f64,
}

impl From<AyanamsaValue> for JsAyanamsaValue {
    fn from(v: AyanamsaValue) -> Self {
        Self {
            model_code: v.model.code().to_owned(),
            jd_tt: v.jd_tt,
            mean_deg: v.mean_deg,
            true_deg: v.true_deg,
            effective_deg: v.effective_deg,
            rate_deg_day: v.rate_deg_day,
            dpsi_arcsec: v.dpsi_arcsec,
        }
    }
}

/// 12 house cusps + 8 chart angles + angular velocities exposed to Node.js / TypeScript.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsHouseResult {
    pub system_code: String,
    pub effective_system_code: String,
    pub polar_house_fallback: bool,
    pub cusps_deg: Vec<f64>,
    pub cusp_speeds_deg_day: Vec<f64>,
    pub asc_deg: f64,
    pub mc_deg: f64,
    pub armc_deg: f64,
    pub vertex_deg: f64,
    pub equatorial_asc_deg: f64,
    pub co_asc_koch_deg: f64,
    pub co_asc_munkasey_deg: f64,
    pub polar_asc_deg: f64,
    pub asc_speed_deg_day: f64,
    pub mc_speed_deg_day: f64,
    pub warnings: Vec<String>,
}

impl From<HousesResult> for JsHouseResult {
    fn from(r: HousesResult) -> Self {
        Self {
            system_code: r.requested_house_system.code().to_string(),
            effective_system_code: r.effective_house_system.code().to_string(),
            polar_house_fallback: r.polar_house_fallback,
            cusps_deg: r.cusps_deg.to_vec(),
            cusp_speeds_deg_day: r.cusp_speeds_deg_day.to_vec(),
            asc_deg: r.angles.asc_deg,
            mc_deg: r.angles.mc_deg,
            armc_deg: r.angles.armc_deg,
            vertex_deg: r.angles.vertex_deg,
            equatorial_asc_deg: r.angles.equatorial_asc_deg,
            co_asc_koch_deg: r.angles.co_asc_koch_deg,
            co_asc_munkasey_deg: r.angles.co_asc_munkasey_deg,
            polar_asc_deg: r.angles.polar_asc_deg,
            asc_speed_deg_day: r.angle_speeds.asc_deg,
            mc_speed_deg_day: r.angle_speeds.mc_deg,
            warnings: r.warnings.into_iter().map(|w| format!("{w:?}")).collect(),
        }
    }
}

/// Topocentric horizontal coordinates (`azalt`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsHorizontalCoords {
    pub azimuth_deg: f64,
    pub altitude_true_deg: f64,
    pub altitude_apparent_deg: f64,
    pub hour_angle_deg: f64,
    pub near_horizon: bool,
}

impl From<HorizontalCoords> for JsHorizontalCoords {
    fn from(h: HorizontalCoords) -> Self {
        Self {
            azimuth_deg: h.azimuth_deg,
            altitude_true_deg: h.altitude_true_deg,
            altitude_apparent_deg: h.altitude_apparent_deg,
            hour_angle_deg: h.hour_angle_deg,
            near_horizon: h.near_horizon,
        }
    }
}

/// Global and local solar eclipse circumstances exposed to Node.js / TypeScript.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsSolarEclipseInfo {
    pub kind: String,
    pub saros_series: i32,
    pub saros_member: u32,
    pub max_eclipse_jd_tt: f64,
    pub max_eclipse_jd_ut1: f64,
    pub delta_t_sec: f64,
    pub gamma: f64,
    pub u_radius_earth: f64,
    pub penumbra_radius_earth: f64,
    pub magnitude: f64,
    pub central_duration_sec: Option<f64>,
    pub path_width_km: Option<f64>,
    pub central_latitude_deg: Option<f64>,
    pub central_longitude_deg: Option<f64>,
    pub p1_jd_tt: Option<f64>,
    pub u1_jd_tt: Option<f64>,
    pub u2_jd_tt: Option<f64>,
    pub u3_jd_tt: Option<f64>,
    pub u4_jd_tt: Option<f64>,
    pub p4_jd_tt: Option<f64>,
    pub local_visible_at_max: Option<bool>,
    pub local_is_central: Option<bool>,
    pub local_kind: Option<String>,
    pub local_c1_jd_tt: Option<f64>,
    pub local_c2_jd_tt: Option<f64>,
    pub local_max_jd_tt: Option<f64>,
    pub local_c3_jd_tt: Option<f64>,
    pub local_c4_jd_tt: Option<f64>,
    pub local_magnitude: Option<f64>,
    pub local_obscuration: Option<f64>,
    pub local_central_duration_sec: Option<f64>,
    pub local_sun_altitude_deg: Option<f64>,
    pub local_sun_azimuth_deg: Option<f64>,
}

impl From<SolarEclipseInfo> for JsSolarEclipseInfo {
    fn from(s: SolarEclipseInfo) -> Self {
        let loc = s.local;
        Self {
            kind: format!("{:?}", s.kind),
            saros_series: s.saros_series,
            saros_member: s.saros_member,
            max_eclipse_jd_tt: s.max_eclipse_jd_tt,
            max_eclipse_jd_ut1: s.max_eclipse_jd_ut1,
            delta_t_sec: s.delta_t_sec,
            gamma: s.gamma,
            u_radius_earth: s.u_radius_earth,
            penumbra_radius_earth: s.penumbra_radius_earth,
            magnitude: s.magnitude,
            central_duration_sec: s.central_duration_sec,
            path_width_km: s.path_width_km,
            central_latitude_deg: s.central_latitude_deg,
            central_longitude_deg: s.central_longitude_deg,
            p1_jd_tt: s.p1_jd_tt,
            u1_jd_tt: s.u1_jd_tt,
            u2_jd_tt: s.u2_jd_tt,
            u3_jd_tt: s.u3_jd_tt,
            u4_jd_tt: s.u4_jd_tt,
            p4_jd_tt: s.p4_jd_tt,
            local_visible_at_max: loc.as_ref().map(|l| l.visible_at_max),
            local_is_central: loc.as_ref().map(|l| l.is_central_locally),
            local_kind: loc.as_ref().map(|l| format!("{:?}", l.local_kind)),
            local_c1_jd_tt: loc.as_ref().and_then(|l| l.c1_jd_tt),
            local_c2_jd_tt: loc.as_ref().and_then(|l| l.c2_jd_tt),
            local_max_jd_tt: loc.as_ref().map(|l| l.max_eclipse_jd_tt),
            local_c3_jd_tt: loc.as_ref().and_then(|l| l.c3_jd_tt),
            local_c4_jd_tt: loc.as_ref().and_then(|l| l.c4_jd_tt),
            local_magnitude: loc.as_ref().map(|l| l.magnitude),
            local_obscuration: loc.as_ref().map(|l| l.obscuration),
            local_central_duration_sec: loc.as_ref().map(|l| l.central_duration_sec),
            local_sun_altitude_deg: loc.as_ref().map(|l| l.sun_altitude_deg),
            local_sun_azimuth_deg: loc.as_ref().map(|l| l.sun_azimuth_deg),
        }
    }
}

/// Lunar eclipse circumstances exposed to Node.js / TypeScript.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsLunarEclipseInfo {
    pub kind: String,
    pub saros_series: i32,
    pub saros_member: u32,
    pub shadow_rule: String,
    pub max_eclipse_jd_tt: f64,
    pub max_eclipse_jd_ut1: f64,
    pub delta_t_sec: f64,
    pub gamma: f64,
    pub penumbral_magnitude: f64,
    pub umbral_magnitude: f64,
    pub penumbral_semi_duration_min: f64,
    pub partial_semi_duration_min: f64,
    pub total_semi_duration_min: f64,
    pub p1_jd_tt: Option<f64>,
    pub u1_jd_tt: Option<f64>,
    pub u2_jd_tt: Option<f64>,
    pub u3_jd_tt: Option<f64>,
    pub u4_jd_tt: Option<f64>,
    pub p4_jd_tt: Option<f64>,
}

impl From<LunarEclipseInfo> for JsLunarEclipseInfo {
    fn from(l: LunarEclipseInfo) -> Self {
        Self {
            kind: format!("{:?}", l.kind),
            saros_series: l.saros_series,
            saros_member: l.saros_member,
            shadow_rule: format!("{:?}", l.shadow_rule),
            max_eclipse_jd_tt: l.max_eclipse_jd_tt,
            max_eclipse_jd_ut1: l.max_eclipse_jd_ut1,
            delta_t_sec: l.delta_t_sec,
            gamma: l.gamma,
            penumbral_magnitude: l.penumbral_magnitude,
            umbral_magnitude: l.umbral_magnitude,
            penumbral_semi_duration_min: l.penumbral_semi_duration_min,
            partial_semi_duration_min: l.partial_semi_duration_min,
            total_semi_duration_min: l.total_semi_duration_min,
            p1_jd_tt: l.p1_jd_tt,
            u1_jd_tt: l.u1_jd_tt,
            u2_jd_tt: l.u2_jd_tt,
            u3_jd_tt: l.u3_jd_tt,
            u4_jd_tt: l.u4_jd_tt,
            p4_jd_tt: l.p4_jd_tt,
        }
    }
}

/// Calendar date returned by `revjul`.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsCalendarDate {
    pub year: i32,
    pub month: u8,
    pub day: u8,
    pub hour_ut: f64,
}

/// Result returned by `JsSweAdapter.calcUt` / `JsSweAdapter.calc`.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsSweCalcOut {
    pub xx: Vec<f64>,
    pub values: Vec<f64>,
    pub retflag: i32,
}

/// Result returned by `JsSweAdapter.fixstar2Ut`.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsSweFixstarOut {
    pub xx: Vec<f64>,
    pub values: Vec<f64>,
    pub star_name: String,
    pub retflag: i32,
}

/// Result returned by `JsSweAdapter.getAyanamsaExUt`.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsSweAyanamsaOut {
    pub retflag: i32,
    pub ayanamsa: f64,
}

/// Result returned by `JsSweAdapter.housesEx` / `houses` / `housesArmc`.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsSweHousesOut {
    pub cusps: Vec<f64>,
    pub ascmc: Vec<f64>,
}

/// Planetary, lunar, or minor-body photometric and illumination phenomena result.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsPhenomenaResult {
    pub phase_angle_deg: f64,
    pub phase_illuminated: f64,
    pub elongation_deg: f64,
    pub apparent_diameter_deg: f64,
    pub apparent_magnitude_v: f64,
    pub horizontal_parallax_deg: f64,
    pub attr20: Vec<f64>,
}

impl From<PhenomenaResult> for JsPhenomenaResult {
    fn from(p: PhenomenaResult) -> Self {
        let attr20 = p.to_swe_attr20().to_vec();
        Self {
            phase_angle_deg: p.phase_angle_deg,
            phase_illuminated: p.phase_illuminated,
            elongation_deg: p.elongation_deg,
            apparent_diameter_deg: p.apparent_diameter_deg,
            apparent_magnitude_v: p.apparent_magnitude_v,
            horizontal_parallax_deg: p.horizontal_parallax_deg,
            attr20,
        }
    }
}

/// Planetary or lunar orbital nodes and apsides calculation result (`swe_nod_aps`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsNodApsResult {
    pub ascending_node: Vec<f64>,
    pub descending_node: Vec<f64>,
    pub periapsis: Vec<f64>,
    pub apoapsis: Vec<f64>,
    pub retflag: i32,
}

impl From<NodesApsidesResult> for JsNodApsResult {
    fn from(n: NodesApsidesResult) -> Self {
        Self {
            ascending_node: n.asc_node.to_vec(),
            descending_node: n.desc_node.to_vec(),
            periapsis: n.periapsis.to_vec(),
            apoapsis: n.apoapsis.to_vec(),
            retflag: 0,
        }
    }
}

/// Global or local lunar occultation circumstances and contact times.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsLunarOccultationInfo {
    pub kind: String,
    pub target: String,
    pub jd_ut1_max: f64,
    pub jd_tt_max: f64,
    pub t1_first_contact_jd_ut1: Option<f64>,
    pub t2_inner_immersion_jd_ut1: Option<f64>,
    pub t3_inner_emersion_jd_ut1: Option<f64>,
    pub t4_last_contact_jd_ut1: Option<f64>,
    pub gamma: f64,
    pub magnitude: f64,
    pub moon_altitude_deg: Option<f64>,
    pub tret: Vec<f64>,
}

impl From<LunarOccultationInfo> for JsLunarOccultationInfo {
    fn from(l: LunarOccultationInfo) -> Self {
        let tret = l.to_swe_tret10().to_vec();
        Self {
            kind: format!("{:?}", l.kind),
            target: format!("{:?}", l.target),
            jd_ut1_max: l.jd_ut1_max,
            jd_tt_max: l.jd_tt_max,
            t1_first_contact_jd_ut1: l.t1_first_contact_jd_ut1,
            t2_inner_immersion_jd_ut1: l.t2_inner_immersion_jd_ut1,
            t3_inner_emersion_jd_ut1: l.t3_inner_emersion_jd_ut1,
            t4_last_contact_jd_ut1: l.t4_last_contact_jd_ut1,
            gamma: l.gamma,
            magnitude: l.magnitude,
            moon_altitude_deg: l.moon_altitude_deg,
            tret,
        }
    }
}

/// Result returned by `JsSweAdapter.phenoUt` / `JsSweAdapter.pheno`.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsSwePhenoOut {
    pub attr: Vec<f64>,
    pub retflag: i32,
}

/// Result returned by `JsSweAdapter.lunOccultWhenGlob`.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsSweOccultGlobOut {
    pub retflag: i32,
    pub tret: Vec<f64>,
}

/// Result returned by `JsSweAdapter.lunOccultWhenLoc`.
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsSweOccultLocOut {
    pub retflag: i32,
    pub tret: Vec<f64>,
    pub attr: Vec<f64>,
}

/// Result of solving for the exact epoch when a body reaches a target ecliptic longitude (`R9.1`, `R9.2`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsLongitudeCrossingEvent {
    pub body: String,
    pub target_lon_deg: f64,
    pub jd_tt: f64,
    pub jd_ut1: f64,
    pub delta_t_sec: f64,
    pub actual_lon_deg: f64,
    pub lat_deg: f64,
    pub dist_au: f64,
    pub dist_km: f64,
    pub lon_speed_deg_day: f64,
    pub is_retrograde: bool,
}

impl From<LongitudeCrossingEvent> for JsLongitudeCrossingEvent {
    fn from(e: LongitudeCrossingEvent) -> Self {
        Self {
            body: body_display_name(e.body),
            target_lon_deg: e.target_lon_deg,
            jd_tt: e.jd_tt,
            jd_ut1: e.jd_ut1,
            delta_t_sec: e.delta_t_sec,
            actual_lon_deg: e.actual_lon_deg,
            lat_deg: e.lat_deg,
            dist_au: e.dist_au,
            dist_km: e.dist_km,
            lon_speed_deg_day: e.lon_speed_deg_day,
            is_retrograde: e.is_retrograde,
        }
    }
}

/// Chinese 24 Solar Term (*Jieqi*) event (`R12`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsSolarTermEvent {
    pub term_index: u8,
    pub target_lon_deg: f64,
    pub pinyin_name: String,
    pub is_major_term: bool,
    pub jd_tt: f64,
    pub jd_ut1: f64,
    pub delta_t_sec: f64,
}

impl From<SolarTermEvent> for JsSolarTermEvent {
    fn from(e: SolarTermEvent) -> Self {
        Self {
            term_index: e.term_index,
            target_lon_deg: e.target_lon_deg,
            pinyin_name: e.pinyin_name.to_owned(),
            is_major_term: e.is_major_term,
            jd_tt: e.jd_tt,
            jd_ut1: e.jd_ut1,
            delta_t_sec: e.delta_t_sec,
        }
    }
}

/// The 4 cardinal astronomical seasons of a calendar year (`0°`, `90°`, `180°`, `270°` apparent solar longitude).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsSeasonsOfYear {
    pub year: i32,
    pub vernal_equinox_jd_ut1: f64,
    pub vernal_equinox_jd_tt: f64,
    pub summer_solstice_jd_ut1: f64,
    pub summer_solstice_jd_tt: f64,
    pub autumnal_equinox_jd_ut1: f64,
    pub autumnal_equinox_jd_tt: f64,
    pub winter_solstice_jd_ut1: f64,
    pub winter_solstice_jd_tt: f64,
}

impl From<SeasonsOfYear> for JsSeasonsOfYear {
    fn from(s: SeasonsOfYear) -> Self {
        Self {
            year: s.year,
            vernal_equinox_jd_ut1: s.vernal_equinox_jd_ut1,
            vernal_equinox_jd_tt: s.vernal_equinox_jd_tt,
            summer_solstice_jd_ut1: s.summer_solstice_jd_ut1,
            summer_solstice_jd_tt: s.summer_solstice_jd_tt,
            autumnal_equinox_jd_ut1: s.autumnal_equinox_jd_ut1,
            autumnal_equinox_jd_tt: s.autumnal_equinox_jd_tt,
            winter_solstice_jd_ut1: s.winter_solstice_jd_ut1,
            winter_solstice_jd_tt: s.winter_solstice_jd_tt,
        }
    }
}

/// Result of solving for a principal or custom lunar phase elongation (`R9.2`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsLunarPhaseEvent {
    pub kind: Option<String>,
    pub target_elongation_deg: f64,
    pub jd_tt: f64,
    pub jd_ut1: f64,
    pub delta_t_sec: f64,
    pub moon_lon_deg: f64,
    pub sun_lon_deg: f64,
    pub moon_dist_au: f64,
    pub moon_dist_km: f64,
}

impl From<LunarPhaseEvent> for JsLunarPhaseEvent {
    fn from(e: LunarPhaseEvent) -> Self {
        Self {
            kind: e.kind.map(|k| format!("{k:?}")),
            target_elongation_deg: e.target_elongation_deg,
            jd_tt: e.jd_tt,
            jd_ut1: e.jd_ut1,
            delta_t_sec: e.delta_t_sec,
            moon_lon_deg: e.moon_lon_deg,
            sun_lon_deg: e.sun_lon_deg,
            moon_dist_au: e.moon_dist_au,
            moon_dist_km: e.moon_dist_km,
        }
    }
}

/// Complete lunar phase, illumination, distance, age, and surrounding quarters summary (`R0`, `R2.2`, `R9.2`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsLunarPhaseSummary {
    pub jd_tt: f64,
    pub jd_ut1: f64,
    pub elongation_deg: f64,
    pub illumination_fraction: f64,
    pub phase_angle_deg: f64,
    pub moon_lon_deg: f64,
    pub moon_lat_deg: f64,
    pub moon_dist_au: f64,
    pub moon_dist_km: f64,
    pub apparent_diameter_deg: f64,
    pub magnitude_v: f64,
    pub age_days: f64,
    pub synodic_fraction: f64,
    pub phase_8_index: u8,
    pub phase_8_name: String,
    pub phase_28_index: u8,
    pub prev_new_moon_jd_ut1: f64,
    pub next_new_moon_jd_ut1: f64,
    pub next_first_quarter_jd_ut1: f64,
    pub next_full_moon_jd_ut1: f64,
    pub next_last_quarter_jd_ut1: f64,
}

impl From<LunarPhaseSummary> for JsLunarPhaseSummary {
    fn from(s: LunarPhaseSummary) -> Self {
        Self {
            jd_tt: s.jd_tt,
            jd_ut1: s.jd_ut1,
            elongation_deg: s.elongation_deg,
            illumination_fraction: s.illumination_fraction,
            phase_angle_deg: s.phase_angle_deg,
            moon_lon_deg: s.moon_lon_deg,
            moon_lat_deg: s.moon_lat_deg,
            moon_dist_au: s.moon_dist_au,
            moon_dist_km: s.moon_dist_km,
            apparent_diameter_deg: s.apparent_diameter_deg,
            magnitude_v: s.magnitude_v,
            age_days: s.age_days,
            synodic_fraction: s.synodic_fraction,
            phase_8_index: s.phase_8_index,
            phase_8_name: s.phase_8_name.to_owned(),
            phase_28_index: s.phase_28_index,
            prev_new_moon_jd_ut1: s.prev_new_moon_jd_ut1,
            next_new_moon_jd_ut1: s.next_new_moon_jd_ut1,
            next_first_quarter_jd_ut1: s.next_first_quarter_jd_ut1,
            next_full_moon_jd_ut1: s.next_full_moon_jd_ut1,
            next_last_quarter_jd_ut1: s.next_last_quarter_jd_ut1,
        }
    }
}

/// Planetary or lunar-point stationary event ($\dot{\lambda}(t) = 0$, `R9.2`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsStationEvent {
    pub body: String,
    pub kind: String,
    pub jd_tt: f64,
    pub jd_ut1: f64,
    pub delta_t_sec: f64,
    pub lon_deg: f64,
    pub lat_deg: f64,
    pub dist_au: f64,
    pub dist_km: f64,
}

impl From<StationEvent> for JsStationEvent {
    fn from(e: StationEvent) -> Self {
        Self {
            body: body_display_name(e.body),
            kind: format!("{:?}", e.kind),
            jd_tt: e.jd_tt,
            jd_ut1: e.jd_ut1,
            delta_t_sec: e.delta_t_sec,
            lon_deg: e.lon_deg,
            lat_deg: e.lat_deg,
            dist_au: e.dist_au,
            dist_km: e.dist_km,
        }
    }
}

/// Lunar or planetary declination extremum event ($\dot{\delta}(t) = 0$, `R9.2`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsDeclinationExtremumEvent {
    pub body: String,
    pub kind: String,
    pub jd_tt: f64,
    pub jd_ut1: f64,
    pub delta_t_sec: f64,
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub dist_au: f64,
    pub dist_km: f64,
}

impl From<DeclinationExtremumEvent> for JsDeclinationExtremumEvent {
    fn from(e: DeclinationExtremumEvent) -> Self {
        Self {
            body: body_display_name(e.body),
            kind: format!("{:?}", e.kind),
            jd_tt: e.jd_tt,
            jd_ut1: e.jd_ut1,
            delta_t_sec: e.delta_t_sec,
            ra_deg: e.ra_deg,
            dec_deg: e.dec_deg,
            dist_au: e.dist_au,
            dist_km: e.dist_km,
        }
    }
}

/// A single Bernadette Brady diurnal fixed-star paran aspect (`R6.2`).
#[cfg_attr(feature = "napi-addon", napi(object))]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JsParanEvent {
    pub star_name: String,
    pub star_bayer: String,
    pub star_hip: u32,
    pub planet_name: String,
    pub paran_type: String,
    pub star_angle: String,
    pub planet_angle: String,
    pub same_angle: bool,
    pub star_event_jd_ut1: f64,
    pub planet_event_jd_ut1: f64,
    pub time_orb_minutes: f64,
}

impl From<ParanEvent> for JsParanEvent {
    fn from(e: ParanEvent) -> Self {
        Self {
            star_name: e.star_name,
            star_bayer: e.star_bayer,
            star_hip: e.star_hip,
            planet_name: e.planet_name,
            paran_type: e.paran_type.as_str().to_owned(),
            star_angle: e.star_angle.as_str().to_owned(),
            planet_angle: e.planet_angle.as_str().to_owned(),
            same_angle: e.same_angle,
            star_event_jd_ut1: e.star_event_jd_ut1,
            planet_event_jd_ut1: e.planet_event_jd_ut1,
            time_orb_minutes: e.time_orb_minutes,
        }
    }
}

/// Named in-memory `.bsp` kernel buffer for `JsEngine.fromBytes`.
#[cfg(feature = "napi-addon")]
#[napi(object)]
pub struct JsKernelBytes {
    pub name: String,
    pub data: napi::bindgen_prelude::Buffer,
}

#[derive(Clone)]
enum EngineSource {
    DataDir {
        dir: PathBuf,
        config: EngineConfig,
    },
    BspPaths {
        paths: Vec<PathBuf>,
        config: EngineConfig,
    },
    Bytes {
        kernels: Vec<(String, Vec<u8>)>,
        config: EngineConfig,
    },
    Moshier {
        config: EngineConfig,
    },
    DataDirOrMoshier {
        dir: PathBuf,
        config: EngineConfig,
    },
    EnvOrMoshier {
        config: EngineConfig,
    },
}

/// Thread-safe Node.js `JsEngine` backed by `lxephemeris_core::Engine`.
#[cfg_attr(feature = "napi-addon", napi)]
#[derive(Clone)]
pub struct JsEngine {
    inner: Arc<std::sync::RwLock<Arc<Engine>>>,
    source: Arc<std::sync::Mutex<EngineSource>>,
    loaded_stars: Arc<std::sync::Mutex<Vec<PathBuf>>>,
}

impl JsEngine {
    #[must_use]
    pub fn engine(&self) -> Arc<Engine> {
        Arc::clone(&*self.inner.read().unwrap_or_else(|e| e.into_inner()))
    }

    /// Creates a `JsEngine`: if `data_dir` is provided, loads `from_data_dir`; otherwise uses `Engine::from_env_or_moshier`.
    pub fn new_core(
        data_dir: Option<String>,
        max_threads: Option<u32>,
        pin_time_table_version: Option<String>,
        verify_manifest_sha256: Option<bool>,
    ) -> Result<Self, String> {
        if let Some(dir_str) = data_dir {
            let trimmed = dir_str.trim();
            if !trimmed.is_empty() {
                return Self::from_data_dir_core(
                    trimmed.to_owned(),
                    max_threads,
                    pin_time_table_version,
                    verify_manifest_sha256,
                );
            }
        }
        let cfg = EngineConfig {
            max_threads: max_threads.map(|t| t as usize),
            pin_time_table_version,
            verify_manifest_sha256: verify_manifest_sha256.unwrap_or(false),
        };
        let eng = Engine::from_env_or_moshier(cfg.clone());
        Ok(Self {
            inner: Arc::new(std::sync::RwLock::new(Arc::new(eng))),
            source: Arc::new(std::sync::Mutex::new(EngineSource::EnvOrMoshier {
                config: cfg,
            })),
            loaded_stars: Arc::new(std::sync::Mutex::new(Vec::new())),
        })
    }

    /// Loads a `JsEngine` from `data_dir`.
    pub fn from_data_dir_core(
        data_dir: String,
        max_threads: Option<u32>,
        pin_time_table_version: Option<String>,
        verify_manifest_sha256: Option<bool>,
    ) -> Result<Self, String> {
        let cfg = EngineConfig {
            max_threads: max_threads.map(|t| t as usize),
            pin_time_table_version,
            verify_manifest_sha256: verify_manifest_sha256.unwrap_or(false),
        };
        let dir = PathBuf::from(data_dir);
        let eng = Engine::from_data_dir(&dir, cfg.clone()).map_err(|e| format_js_error(&e))?;
        Ok(Self {
            inner: Arc::new(std::sync::RwLock::new(Arc::new(eng))),
            source: Arc::new(std::sync::Mutex::new(EngineSource::DataDir {
                dir,
                config: cfg,
            })),
            loaded_stars: Arc::new(std::sync::Mutex::new(Vec::new())),
        })
    }

    /// Loads a `JsEngine` from explicit `.bsp` file paths.
    pub fn from_bsp_paths_core(
        paths: Vec<String>,
        max_threads: Option<u32>,
        pin_time_table_version: Option<String>,
    ) -> Result<Self, String> {
        let cfg = EngineConfig {
            max_threads: max_threads.map(|t| t as usize),
            pin_time_table_version,
            verify_manifest_sha256: false,
        };
        let pb_paths: Vec<PathBuf> = paths.into_iter().map(PathBuf::from).collect();
        let eng =
            Engine::from_bsp_paths(&pb_paths, cfg.clone()).map_err(|e| format_js_error(&e))?;
        Ok(Self {
            inner: Arc::new(std::sync::RwLock::new(Arc::new(eng))),
            source: Arc::new(std::sync::Mutex::new(EngineSource::BspPaths {
                paths: pb_paths,
                config: cfg,
            })),
            loaded_stars: Arc::new(std::sync::Mutex::new(Vec::new())),
        })
    }

    /// Loads a `JsEngine` from in-memory `(name, bytes)` kernel slices.
    pub fn from_bytes_core(
        kernels: &[(String, Vec<u8>)],
        max_threads: Option<u32>,
    ) -> Result<Self, String> {
        let cfg = EngineConfig {
            max_threads: max_threads.map(|t| t as usize),
            ..Default::default()
        };
        let eng = Engine::from_bytes(kernels, cfg.clone()).map_err(|e| format_js_error(&e))?;
        Ok(Self {
            inner: Arc::new(std::sync::RwLock::new(Arc::new(eng))),
            source: Arc::new(std::sync::Mutex::new(EngineSource::Bytes {
                kernels: kernels.to_vec(),
                config: cfg,
            })),
            loaded_stars: Arc::new(std::sync::Mutex::new(Vec::new())),
        })
    }

    /// Creates a fileless analytical ("Moshier") `JsEngine` (-3000..+3000 CE) requiring zero `.bsp` kernels.
    pub fn new_moshier_core(
        max_threads: Option<u32>,
        pin_time_table_version: Option<String>,
    ) -> Self {
        let cfg = EngineConfig {
            max_threads: max_threads.map(|t| t as usize),
            pin_time_table_version,
            verify_manifest_sha256: false,
        };
        let eng = Engine::new_moshier(cfg.clone());
        Self {
            inner: Arc::new(std::sync::RwLock::new(Arc::new(eng))),
            source: Arc::new(std::sync::Mutex::new(EngineSource::Moshier { config: cfg })),
            loaded_stars: Arc::new(std::sync::Mutex::new(Vec::new())),
        }
    }

    /// Loads a `JsEngine` from `data_dir` if kernels exist, or seamlessly falls back to Moshier analytical mode.
    pub fn from_data_dir_or_moshier_core(
        data_dir: String,
        max_threads: Option<u32>,
        pin_time_table_version: Option<String>,
        verify_manifest_sha256: Option<bool>,
    ) -> Result<Self, String> {
        let cfg = EngineConfig {
            max_threads: max_threads.map(|t| t as usize),
            pin_time_table_version,
            verify_manifest_sha256: verify_manifest_sha256.unwrap_or(false),
        };
        let dir = PathBuf::from(data_dir);
        let eng =
            Engine::from_data_dir_or_moshier(&dir, cfg.clone()).map_err(|e| format_js_error(&e))?;
        Ok(Self {
            inner: Arc::new(std::sync::RwLock::new(Arc::new(eng))),
            source: Arc::new(std::sync::Mutex::new(EngineSource::DataDirOrMoshier {
                dir,
                config: cfg,
            })),
            loaded_stars: Arc::new(std::sync::Mutex::new(Vec::new())),
        })
    }

    /// Returns capabilities metadata as a JSON string.
    pub fn capabilities_json_core(&self) -> Result<String, String> {
        serde_json::to_string(&self.engine().capabilities()).map_err(|e| e.to_string())
    }

    /// Returns coverage metadata for `body` as a JSON string.
    pub fn coverage_json_core(&self, body: &str) -> Result<String, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let cov = self
            .engine()
            .coverage(bid)
            .map_err(|e| format_js_error(&e))?;
        serde_json::to_string(&cov).map_err(|e| e.to_string())
    }

    /// Evaluates $\Delta T = \text{TT} - \text{UT1}$ and $\text{TT} - \text{TDB}$ (`R1.2`).
    pub fn delta_t_core(
        &self,
        jd: JsSplitJd,
        scale: Option<String>,
    ) -> Result<JsDeltaTInfo, String> {
        let t_in = parse_time_input(scale.as_deref().unwrap_or("UT1"), jd)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsDeltaTInfo::from(self.engine().delta_t(t_in)))
    }

    /// Evaluates GMST, GAST, Local Sidereal Time, and ARMC (`R1.4`).
    pub fn sidtime_core(
        &self,
        jd: JsSplitJd,
        scale: Option<String>,
        nutation: Option<String>,
        longitude_deg: Option<f64>,
    ) -> Result<JsSiderealTimeInfo, String> {
        let t_in = parse_time_input(scale.as_deref().unwrap_or("UT1"), jd)
            .map_err(|e| format_js_error(&e))?;
        let nut = parse_nutation_mode(nutation.as_deref().unwrap_or("true"))
            .map_err(|e| format_js_error(&e))?;
        let info = self
            .engine()
            .sidtime(t_in, nut, longitude_deg.unwrap_or(0.0));
        Ok(JsSiderealTimeInfo::from(info))
    }

    /// Evaluates mean/true obliquity and nutation in longitude/obliquity (`ECL_NUT`, `R2.2`).
    pub fn ecl_nut_core(
        &self,
        jd: JsSplitJd,
        scale: Option<String>,
    ) -> Result<JsEclNutInfo, String> {
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let info = self
            .engine()
            .ecl_nut(t_in)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsEclNutInfo::from(info))
    }

    /// Looks up catalog metadata (including `vmag`) for a fixed star (`R6.1`).
    pub fn star_info_core(&self, name_or_bayer: &str) -> Result<JsStarMetadata, String> {
        let info = self
            .engine()
            .star_info(name_or_bayer)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsStarMetadata::from(info))
    }

    /// Lists all 80 cataloged fixed stars (`R6.1`).
    pub fn list_stars_core(&self) -> Vec<JsStarMetadata> {
        self.engine()
            .list_stars()
            .into_iter()
            .map(JsStarMetadata::from)
            .collect()
    }

    /// Returns the active license metadata (`status`, `tier`, `isCommercialValid`, `licenseId`, `holder`, etc.).
    pub fn license_info_core(&self) -> Result<serde_json::Value, String> {
        let info = self.engine().license_info();
        Ok(serde_json::json!({
            "status": info.status.as_str(),
            "tier": info.tier,
            "isCommercialValid": info.is_commercial_valid,
            "licenseId": info.license_id,
            "holder": info.holder,
            "keyId": info.key_id,
            "issuedAt": info.issued_at,
            "updatesUntil": info.updates_until,
            "engineReleaseDate": info.engine_release_date,
            "message": info.message,
        }))
    }

    /// Evaluates a single body at `jd` (`JsSplitJd`) in `scale`.
    pub fn calc_core(
        &self,
        body: &str,
        jd: JsSplitJd,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> Result<JsCalcResult, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let res = self
            .engine()
            .calc(bid, t_in, opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsCalcResult::from(res))
    }

    /// Evaluates a fixed star by name, Bayer designation, or `"hip:<id>"`.
    pub fn calc_star_by_name_core(
        &self,
        name_or_bayer: &str,
        jd: JsSplitJd,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> Result<JsCalcResult, String> {
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let res = self
            .engine()
            .calc_star_by_name(name_or_bayer, t_in, opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsCalcResult::from(res))
    }

    /// Evaluates a batch of `bodies` across `epochs` (`Vec<JsSplitJd>`) in atomic mode.
    pub fn calc_batch_core(
        &self,
        bodies: &[String],
        epochs: &[JsSplitJd],
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<Vec<JsCalcResult>>, String> {
        let bids: Vec<BodyId> = bodies
            .iter()
            .map(|b| parse_body_id(b))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| format_js_error(&e))?;
        let sc = scale.as_deref().unwrap_or("TT");
        let t_inputs: Vec<TimeInput> = epochs
            .iter()
            .map(|&jd| parse_time_input(sc, jd))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let batch = self
            .engine()
            .calc_batch(&bids, &t_inputs, opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(batch
            .into_iter()
            .map(|row| row.into_iter().map(JsCalcResult::from).collect())
            .collect())
    }

    /// Evaluates a batch of `bodies` across `epochs` (`Vec<JsSplitJd>`) and packs the result into a
    /// contiguous `Vec<f64>` (`Float64Array` on Node.js) of shape `[epochs.len(), bodies.len(), 6]`.
    pub fn calc_batch_flat_f64_core(
        &self,
        bodies: &[String],
        epochs: &[JsSplitJd],
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<f64>, String> {
        let bids: Vec<BodyId> = bodies
            .iter()
            .map(|b| parse_body_id(b))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| format_js_error(&e))?;
        let sc = scale.as_deref().unwrap_or("TT");
        let t_inputs: Vec<TimeInput> = epochs
            .iter()
            .map(|&jd| parse_time_input(sc, jd))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let batch = self
            .engine()
            .calc_batch(&bids, &t_inputs, opts)
            .map_err(|e| format_js_error(&e))?;
        let mut out = Vec::with_capacity(t_inputs.len() * bids.len() * 6);
        for row in batch {
            for r in row {
                out.extend_from_slice(&[
                    r.lon_deg,
                    r.lat_deg,
                    r.dist_au,
                    r.lon_speed_deg_day,
                    r.lat_speed_deg_day,
                    r.dist_speed_au_day,
                ]);
            }
        }
        Ok(out)
    }

    /// Evaluates a batch of `bodies` across `epochs` in per-item mode without aborting on individual epoch errors.
    pub fn calc_batch_per_item_core(
        &self,
        bodies: &[String],
        epochs: &[JsSplitJd],
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<Vec<JsBatchItemResult>>, String> {
        let bids: Vec<BodyId> = bodies
            .iter()
            .map(|b| parse_body_id(b))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| format_js_error(&e))?;
        let sc = scale.as_deref().unwrap_or("TT");
        let t_inputs: Vec<TimeInput> = epochs
            .iter()
            .map(|&jd| parse_time_input(sc, jd))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let batch = self
            .engine()
            .calc_batch_per_item(&bids, &t_inputs, opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(batch
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|item| match item {
                        Ok(r) => JsBatchItemResult {
                            ok: true,
                            result: Some(JsCalcResult::from(r)),
                            error_code: None,
                            error_message: None,
                        },
                        Err(e) => JsBatchItemResult {
                            ok: false,
                            result: None,
                            error_code: Some(e.code().as_str().to_owned()),
                            error_message: Some(e.to_string()),
                        },
                    })
                    .collect()
            })
            .collect())
    }

    /// Evaluates any of the 12 Vedic / Sidereal ayanamsa models at `jd` in `scale`.
    pub fn ayanamsa_core(
        &self,
        model: &str,
        jd: JsSplitJd,
        scale: Option<String>,
        nutation: Option<String>,
    ) -> Result<JsAyanamsaValue, String> {
        let m = AyanamsaModel::from_code(model).map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let nut = parse_nutation_mode(nutation.as_deref().unwrap_or("true"))
            .map_err(|e| format_js_error(&e))?;
        let val = self
            .engine()
            .ayanamsa(t_in, m, nut)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsAyanamsaValue::from(val))
    }

    /// Evaluates any of the 23 house systems at `jd` in `scale` for `observer`.
    #[allow(clippy::too_many_arguments)]
    pub fn houses_core(
        &self,
        jd: JsSplitJd,
        observer: JsObserver,
        system: Option<String>,
        scale: Option<String>,
        polar_policy: Option<String>,
        zodiac: Option<String>,
        ayanamsa: Option<String>,
        nutation: Option<String>,
        include_speed: Option<bool>,
    ) -> Result<JsHouseResult, String> {
        let sys_str = system.as_deref().unwrap_or("P");
        let ch = sys_str.chars().next().unwrap_or(' ');
        let sys = HouseSystem::from_code(ch).map_err(|e| format_js_error(&e))?;
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let pol = parse_polar_policy(polar_policy.as_deref().unwrap_or("strict"))
            .map_err(|e| format_js_error(&e))?;
        let zod = parse_zodiac(
            zodiac.as_deref().unwrap_or("tropical"),
            ayanamsa.as_deref().unwrap_or("lahiri"),
        )
        .map_err(|e| format_js_error(&e))?;
        let nut = parse_nutation_mode(nutation.as_deref().unwrap_or("true"))
            .map_err(|e| format_js_error(&e))?;
        let opts = HouseOptions {
            system: sys,
            polar_policy: pol,
            zodiac: zod,
            nutation: nut,
            include_speed: include_speed.unwrap_or(true),
        };
        let res = self
            .engine()
            .houses(t_in, obs, opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsHouseResult::from(res))
    }

    /// Evaluates house cusps and chart angles across a batch of `epochs` (`R9.3`).
    #[allow(clippy::too_many_arguments)]
    pub fn houses_batch_core(
        &self,
        epochs: &[JsSplitJd],
        observer: JsObserver,
        system: Option<String>,
        scale: Option<String>,
        polar_policy: Option<String>,
        zodiac: Option<String>,
        ayanamsa: Option<String>,
        nutation: Option<String>,
        include_speed: Option<bool>,
    ) -> Result<Vec<JsHouseResult>, String> {
        let sys_str = system.as_deref().unwrap_or("P");
        let ch = sys_str.chars().next().unwrap_or(' ');
        let sys = HouseSystem::from_code(ch).map_err(|e| format_js_error(&e))?;
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let sc = scale.as_deref().unwrap_or("TT");
        let t_inputs: Vec<TimeInput> = epochs
            .iter()
            .map(|&jd| parse_time_input(sc, jd))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| format_js_error(&e))?;
        let pol = parse_polar_policy(polar_policy.as_deref().unwrap_or("strict"))
            .map_err(|e| format_js_error(&e))?;
        let zod = parse_zodiac(
            zodiac.as_deref().unwrap_or("tropical"),
            ayanamsa.as_deref().unwrap_or("lahiri"),
        )
        .map_err(|e| format_js_error(&e))?;
        let nut = parse_nutation_mode(nutation.as_deref().unwrap_or("true"))
            .map_err(|e| format_js_error(&e))?;
        let opts = HouseOptions {
            system: sys,
            polar_policy: pol,
            zodiac: zod,
            nutation: nut,
            include_speed: include_speed.unwrap_or(true),
        };
        let res = self
            .engine()
            .houses_batch(&t_inputs, obs, opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(res.into_iter().map(JsHouseResult::from).collect())
    }

    /// Evaluates topocentric horizontal coordinates (`azalt`) with optical refraction from `(ra_deg, dec_deg)`.
    #[allow(clippy::too_many_arguments)]
    pub fn azalt_core(
        &self,
        ra_deg: f64,
        dec_deg: f64,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        pressure_hpa: Option<f64>,
        temperature_c: Option<f64>,
    ) -> Result<JsHorizontalCoords, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let eng = self.engine();
        let moment = eng.at(t_in, Some(obs)).map_err(|e| format_js_error(&e))?;
        let atmo = AtmosphereParams {
            pressure_hpa: pressure_hpa.unwrap_or(1013.25),
            temperature_c: temperature_c.unwrap_or(15.0),
            ..Default::default()
        };
        let hz = core_azalt(obs, moment.center().gast_rad, ra_deg, dec_deg, atmo);
        Ok(JsHorizontalCoords::from(hz))
    }

    /// Evaluates topocentric horizontal coordinates (`azalt`) from ecliptic `(ecl_lon_deg, ecl_lat_deg)` (`R2.6`).
    #[allow(clippy::too_many_arguments)]
    pub fn azalt_ecliptic_core(
        &self,
        ecl_lon_deg: f64,
        ecl_lat_deg: f64,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        pressure_hpa: Option<f64>,
        temperature_c: Option<f64>,
    ) -> Result<JsHorizontalCoords, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let atmo = AtmosphereParams {
            pressure_hpa: pressure_hpa.unwrap_or(1013.25),
            temperature_c: temperature_c.unwrap_or(15.0),
            ..Default::default()
        };
        let hz = self
            .engine()
            .azalt_ecliptic(obs, t_in, ecl_lon_deg, ecl_lat_deg, atmo);
        Ok(JsHorizontalCoords::from(hz))
    }

    /// Evaluates topocentric horizontal coordinates (`azalt`) directly for `body` (`R2.6`).
    #[allow(clippy::too_many_arguments)]
    pub fn calc_azalt_core(
        &self,
        body: &str,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        pressure_hpa: Option<f64>,
        temperature_c: Option<f64>,
    ) -> Result<JsHorizontalCoords, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let atmo = AtmosphereParams {
            pressure_hpa: pressure_hpa.unwrap_or(1013.25),
            temperature_c: temperature_c.unwrap_or(15.0),
            ..Default::default()
        };
        let hz = self
            .engine()
            .calc_azalt(bid, t_in, obs, atmo)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsHorizontalCoords::from(hz))
    }

    /// Searches for the next rise, set, upper meridian transit, or lower meridian transit (`JD_UT1`).
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans_core(
        &self,
        body: &str,
        event: &str,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        custom_horizon_alt_deg: Option<f64>,
        max_search_days: Option<f64>,
        disc_mode: Option<String>,
        pressure_hpa: Option<f64>,
        temperature_c: Option<f64>,
    ) -> Result<Option<f64>, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("UT1"), jd)
            .map_err(|e| format_js_error(&e))?;
        let eng = self.engine();
        let moment = eng.at(t_in, Some(obs)).map_err(|e| format_js_error(&e))?;
        let jd_ut1_start = moment.center().ut1.jd_f64();
        let ev = match event.trim().to_ascii_lowercase().as_str() {
            "rise" => RiseTransEvent::Rise,
            "set" => RiseTransEvent::Set,
            "upper" | "upper_transit" | "culmination" => RiseTransEvent::UpperTransit,
            "lower" | "lower_transit" | "anti_culmination" => RiseTransEvent::LowerTransit,
            other => {
                return Err(format_js_error(&EphemerisError::InvalidInput {
                    message: format!("Unknown transit event '{other}'"),
                }))
            }
        };
        let dm = parse_disc_mode(disc_mode.as_deref().unwrap_or("upper_limb"))
            .map_err(|e| format_js_error(&e))?;
        let opts = RiseTransOptions {
            max_search_days: max_search_days.unwrap_or(2.0),
            custom_horizon_alt_deg,
            disc_mode: dm,
            pressure_hpa: pressure_hpa.unwrap_or(1013.25),
            temperature_c: temperature_c.unwrap_or(15.0),
        };
        self.engine()
            .rise_trans(bid, obs, jd_ut1_start, ev, opts)
            .map_err(|e| format_js_error(&e))
    }

    /// Computes heliacal risings and settings via Bradley Schaefer's visibility model (`R6.3`).
    #[allow(clippy::too_many_arguments)]
    pub fn heliacal_phenomena_core(
        &self,
        body: &str,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        extinction_k_v: Option<f64>,
        snellen_factor: Option<f64>,
        pressure_hpa: Option<f64>,
        temperature_c: Option<f64>,
        relative_humidity: Option<f64>,
        max_search_days: Option<f64>,
        observer_age_years: Option<f64>,
        snellen_ratio: Option<f64>,
    ) -> Result<JsHeliacalResult, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("UT1"), jd)
            .map_err(|e| format_js_error(&e))?;
        let eng = self.engine();
        let moment = eng.at(t_in, Some(obs)).map_err(|e| format_js_error(&e))?;
        let jd_ut1_start = moment.center().ut1.jd_f64();
        let opts = HeliacalOptions {
            k_v: extinction_k_v.unwrap_or(0.20),
            observer_age_years: observer_age_years.unwrap_or(36.0),
            snellen_ratio: snellen_ratio.or(snellen_factor).unwrap_or(1.0),
            pressure_hpa: pressure_hpa.unwrap_or(1013.25),
            temperature_c: temperature_c.unwrap_or(15.0),
            relative_humidity: relative_humidity.unwrap_or(0.50),
            max_search_days: max_search_days.unwrap_or(400.0),
        };
        let res = self
            .engine()
            .heliacal_phenomena(bid, jd_ut1_start, obs, opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsHeliacalResult::from(res))
    }

    /// Searches forward or backward for a solar eclipse from `jd` in `scale`.
    pub fn next_solar_eclipse_core(
        &self,
        jd: JsSplitJd,
        scale: Option<String>,
        backward: Option<bool>,
        observer: Option<JsObserver>,
        max_synodic_months: Option<u32>,
    ) -> Result<JsSolarEclipseInfo, String> {
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let obs = match observer {
            Some(o) => Some(
                Observer::new(o.latitude_deg, o.longitude_deg, o.altitude_m.unwrap_or(0.0))
                    .map_err(|e| format_js_error(&e))?,
            ),
            None => None,
        };
        let opts = EclipseOptions {
            max_synodic_months: max_synodic_months.unwrap_or(36),
            ..Default::default()
        };
        let info = self
            .engine()
            .next_solar_eclipse(t_in, backward.unwrap_or(false), obs, opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSolarEclipseInfo::from(info))
    }

    /// Searches forward or backward for the next solar eclipse that is locally visible (`R8.1`).
    pub fn next_local_solar_eclipse_core(
        &self,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        backward: Option<bool>,
        max_synodic_months: Option<u32>,
    ) -> Result<JsSolarEclipseInfo, String> {
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let opts = EclipseOptions {
            max_synodic_months: max_synodic_months.unwrap_or(120),
            ..Default::default()
        };
        let info = self
            .engine()
            .next_local_solar_eclipse(t_in, obs, backward.unwrap_or(false), opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSolarEclipseInfo::from(info))
    }

    /// Enumerates all solar eclipses in `[start_jd_tt, end_jd_tt]` (`R8.3`).
    pub fn solar_eclipses_in_range_core(
        &self,
        start_jd_tt: f64,
        end_jd_tt: f64,
        observer: Option<JsObserver>,
    ) -> Result<Vec<JsSolarEclipseInfo>, String> {
        let obs = match observer {
            Some(o) => Some(
                Observer::new(o.latitude_deg, o.longitude_deg, o.altitude_m.unwrap_or(0.0))
                    .map_err(|e| format_js_error(&e))?,
            ),
            None => None,
        };
        let list = self
            .engine()
            .solar_eclipses_in_range(
                TimeInput::from_tt_jd(start_jd_tt),
                TimeInput::from_tt_jd(end_jd_tt),
                obs,
                EclipseOptions::default(),
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(list.into_iter().map(JsSolarEclipseInfo::from).collect())
    }

    /// Searches forward or backward for a lunar eclipse from `jd` in `scale`.
    pub fn next_lunar_eclipse_core(
        &self,
        jd: JsSplitJd,
        scale: Option<String>,
        backward: Option<bool>,
        shadow_rule: Option<String>,
        max_synodic_months: Option<u32>,
    ) -> Result<JsLunarEclipseInfo, String> {
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let rule = parse_shadow_rule(shadow_rule.as_deref().unwrap_or("danjon"))
            .map_err(|e| format_js_error(&e))?;
        let opts = EclipseOptions {
            shadow_rule: rule,
            max_synodic_months: max_synodic_months.unwrap_or(36),
            ..Default::default()
        };
        let info = self
            .engine()
            .next_lunar_eclipse(t_in, backward.unwrap_or(false), opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsLunarEclipseInfo::from(info))
    }

    /// Enumerates all lunar eclipses in `[start_jd_tt, end_jd_tt]` (`R8.3`).
    pub fn lunar_eclipses_in_range_core(
        &self,
        start_jd_tt: f64,
        end_jd_tt: f64,
        shadow_rule: Option<String>,
    ) -> Result<Vec<JsLunarEclipseInfo>, String> {
        let rule = parse_shadow_rule(shadow_rule.as_deref().unwrap_or("danjon"))
            .map_err(|e| format_js_error(&e))?;
        let opts = EclipseOptions {
            shadow_rule: rule,
            ..Default::default()
        };
        let list = self
            .engine()
            .lunar_eclipses_in_range(
                TimeInput::from_tt_jd(start_jd_tt),
                TimeInput::from_tt_jd(end_jd_tt),
                opts,
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(list.into_iter().map(JsLunarEclipseInfo::from).collect())
    }

    /// Computes the 36 Gauquelin sectors (`Vec<f64>` of length 36) from `armc_deg`, `eps_true_deg`, `latitude_deg`.
    pub fn gauquelin_sectors_core(
        &self,
        armc_deg: f64,
        eps_true_deg: f64,
        latitude_deg: f64,
        polar_policy: Option<String>,
    ) -> Result<Vec<f64>, String> {
        gauquelin_sectors_core(armc_deg, eps_true_deg, latitude_deg, polar_policy)
    }

    /// Returns `true` if this engine is operating in pure fileless analytical ("Moshier") mode.
    pub fn is_moshier_mode_core(&self) -> bool {
        self.engine().is_moshier_mode()
    }

    /// Evaluates `body` at `jd` explicitly using the fileless analytical Moshier backend (`[-3000 CE, +3000 CE]`).
    pub fn calc_moshier_core(
        &self,
        body: &str,
        jd: JsSplitJd,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> Result<JsCalcResult, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let res = self
            .engine()
            .calc_moshier(bid, t_in, opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsCalcResult::from(res))
    }

    /// Loads additional custom stars from a `sefstars.txt` or JSON file.
    pub fn load_custom_stars_core(&self, path: String) -> Result<u32, String> {
        let p = PathBuf::from(&path);
        let stars = load_custom_stars_file(&p).map_err(|e| format_js_error(&e))?;
        let count = stars.len() as u32;

        let mut loaded = self.loaded_stars.lock().map_err(|e| e.to_string())?;
        loaded.push(p);

        let mut w = self.inner.write().map_err(|e| e.to_string())?;
        if let Some(eng_mut) = Arc::get_mut(&mut *w) {
            eng_mut.add_custom_stars(stars);
        } else {
            let src = self.source.lock().map_err(|e| e.to_string())?;
            let mut new_eng = match &*src {
                EngineSource::DataDir { dir, config } => {
                    Engine::from_data_dir(dir, config.clone()).map_err(|e| format_js_error(&e))?
                }
                EngineSource::BspPaths { paths, config } => {
                    Engine::from_bsp_paths(paths, config.clone())
                        .map_err(|e| format_js_error(&e))?
                }
                EngineSource::Bytes { kernels, config } => {
                    Engine::from_bytes(kernels, config.clone()).map_err(|e| format_js_error(&e))?
                }
                EngineSource::Moshier { config } => Engine::new_moshier(config.clone()),
                EngineSource::DataDirOrMoshier { dir, config } => {
                    Engine::from_data_dir_or_moshier(dir, config.clone())
                        .map_err(|e| format_js_error(&e))?
                }
                EngineSource::EnvOrMoshier { config } => {
                    Engine::from_env_or_moshier(config.clone())
                }
            };
            for star_path in loaded.iter() {
                let _ = new_eng.load_custom_stars(star_path);
            }
            *w = Arc::new(new_eng);
        }
        Ok(count)
    }

    /// Evaluates planetary, solar, lunar, or minor-body photometric and illumination phenomena.
    pub fn phenomena_core(
        &self,
        body: &str,
        jd: JsSplitJd,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> Result<JsPhenomenaResult, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let res = self
            .engine()
            .phenomena(bid, t_in, opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsPhenomenaResult::from(res))
    }

    /// Evaluates planetary or lunar orbital nodes and apsides.
    pub fn nod_aps_core(
        &self,
        body: &str,
        jd: JsSplitJd,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
        method_bits: Option<u32>,
    ) -> Result<JsNodApsResult, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let bits = method_bits.unwrap_or(0);
        let res = self
            .engine()
            .nod_aps(bid, t_in, opts, bits)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsNodApsResult::from(res))
    }

    /// Searches forward or backward for the next global lunar occultation.
    pub fn next_lunar_occultation_core(
        &self,
        target: &str,
        jd: JsSplitJd,
        scale: Option<String>,
        backward: Option<bool>,
        max_search_days: Option<f64>,
    ) -> Result<JsLunarOccultationInfo, String> {
        let eng = self.engine();
        let target_body = if let Ok(bid) = parse_body_id(target) {
            bid
        } else {
            let meta = eng.star_info(target).map_err(|e| format_js_error(&e))?;
            BodyId::FixedStar(HipId(meta.hip))
        };
        let t_in = parse_time_input(scale.as_deref().unwrap_or("UT1"), jd)
            .map_err(|e| format_js_error(&e))?;
        let info = eng
            .next_lunar_occultation(
                target_body,
                t_in,
                backward.unwrap_or(false),
                max_search_days.unwrap_or(3652.5),
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(JsLunarOccultationInfo::from(info))
    }

    /// Searches forward or backward for the next local lunar occultation.
    pub fn next_local_lunar_occultation_core(
        &self,
        target: &str,
        observer: JsObserver,
        jd: JsSplitJd,
        scale: Option<String>,
        backward: Option<bool>,
        max_search_days: Option<f64>,
    ) -> Result<JsLunarOccultationInfo, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let eng = self.engine();
        let target_body = if let Ok(bid) = parse_body_id(target) {
            bid
        } else {
            let meta = eng.star_info(target).map_err(|e| format_js_error(&e))?;
            BodyId::FixedStar(HipId(meta.hip))
        };
        let t_in = parse_time_input(scale.as_deref().unwrap_or("UT1"), jd)
            .map_err(|e| format_js_error(&e))?;
        let info = eng
            .next_local_lunar_occultation(
                target_body,
                obs,
                t_in,
                backward.unwrap_or(false),
                max_search_days.unwrap_or(3652.5),
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(JsLunarOccultationInfo::from(info))
    }

    /// Finds all global lunar occultations in range.
    pub fn lunar_occultations_in_range_core(
        &self,
        target: &str,
        start_jd: JsSplitJd,
        end_jd: JsSplitJd,
        scale: Option<String>,
    ) -> Result<Vec<JsLunarOccultationInfo>, String> {
        let sc = scale.as_deref().unwrap_or("UT1");
        let t_start = parse_time_input(sc, start_jd).map_err(|e| format_js_error(&e))?;
        let t_end = parse_time_input(sc, end_jd).map_err(|e| format_js_error(&e))?;
        let eng = self.engine();
        let target_body = if let Ok(bid) = parse_body_id(target) {
            bid
        } else {
            let meta = eng.star_info(target).map_err(|e| format_js_error(&e))?;
            BodyId::FixedStar(HipId(meta.hip))
        };
        let list = eng
            .lunar_occultations_in_range(target_body, t_start, t_end)
            .map_err(|e| format_js_error(&e))?;
        Ok(list.into_iter().map(JsLunarOccultationInfo::from).collect())
    }

    /// Finds all local lunar occultations in range.
    pub fn local_lunar_occultations_in_range_core(
        &self,
        target: &str,
        observer: JsObserver,
        start_jd: JsSplitJd,
        end_jd: JsSplitJd,
        scale: Option<String>,
    ) -> Result<Vec<JsLunarOccultationInfo>, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let sc = scale.as_deref().unwrap_or("UT1");
        let t_start = parse_time_input(sc, start_jd).map_err(|e| format_js_error(&e))?;
        let t_end = parse_time_input(sc, end_jd).map_err(|e| format_js_error(&e))?;
        let eng = self.engine();
        let target_body = if let Ok(bid) = parse_body_id(target) {
            bid
        } else {
            let meta = eng.star_info(target).map_err(|e| format_js_error(&e))?;
            BodyId::FixedStar(HipId(meta.hip))
        };
        let list = eng
            .local_lunar_occultations_in_range(target_body, obs, t_start, t_end)
            .map_err(|e| format_js_error(&e))?;
        Ok(list.into_iter().map(JsLunarOccultationInfo::from).collect())
    }

    /// Computes the 3D fractional mundane house position `[1.0, 13.0)`.
    pub fn house_pos_core(
        &self,
        armc_deg: f64,
        geolat_deg: f64,
        eps_deg: f64,
        system: &str,
        ecl_lon_deg: f64,
        ecl_lat_deg: f64,
    ) -> Result<f64, String> {
        let ch = system.chars().next().unwrap_or('P');
        let sys = HouseSystem::from_code(ch).map_err(|e| format_js_error(&e))?;
        self.engine()
            .house_pos(armc_deg, geolat_deg, eps_deg, sys, ecl_lon_deg, ecl_lat_deg)
            .map_err(|e| format_js_error(&e))
    }

    /// Evaluates the Equation of Time in days.
    pub fn time_equation_core(&self, jd: JsSplitJd, scale: Option<String>) -> Result<f64, String> {
        let t_in = parse_time_input(scale.as_deref().unwrap_or("TT"), jd)
            .map_err(|e| format_js_error(&e))?;
        self.engine()
            .time_equation(t_in)
            .map_err(|e| format_js_error(&e))
    }

    /// Inverse horizontal transformation (`swe_azalt_rev`).
    pub fn azalt_rev_core(
        &self,
        observer: JsObserver,
        jd: JsSplitJd,
        scale: Option<String>,
        azimuth_north_east_deg: f64,
        altitude_true_deg: f64,
        target_plane: Option<String>,
    ) -> Result<Vec<f64>, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let t_in = parse_time_input(scale.as_deref().unwrap_or("UT1"), jd)
            .map_err(|e| format_js_error(&e))?;
        let plane = match target_plane
            .as_deref()
            .unwrap_or("ecliptic")
            .to_ascii_lowercase()
            .as_str()
        {
            "equatorial" | "equ" => Plane::Equatorial,
            _ => Plane::Ecliptic,
        };
        let res =
            self.engine()
                .azalt_rev(obs, t_in, azimuth_north_east_deg, altitude_true_deg, plane);
        Ok(res.to_vec())
    }

    /// Finds the first instant starting from `jd_ut1` when `body`'s ecliptic longitude reaches `target_lon_deg` (`R9.1`, `R9.2`).
    pub fn search_body_longitude_core(
        &self,
        body: &str,
        jd_ut1: f64,
        target_lon_deg: f64,
        forward: Option<bool>,
        max_search_days: Option<f64>,
        options: Option<JsCalcOptions>,
    ) -> Result<JsLongitudeCrossingEvent, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let ev = self
            .engine()
            .search_body_longitude(
                bid,
                TimeInput::from_ut1_jd(jd_ut1),
                target_lon_deg,
                forward.unwrap_or(true),
                max_search_days.unwrap_or(400.0),
                &opts,
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(JsLongitudeCrossingEvent::from(ev))
    }

    /// Finds all crossings of `target_lon_deg` by `body` in `[start_jd_ut1, end_jd_ut1]`, including 3-pass retrograde returns (`R9.1`, `R9.2`).
    pub fn body_longitude_crossings_in_range_core(
        &self,
        body: &str,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        target_lon_deg: f64,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<JsLongitudeCrossingEvent>, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let list = self
            .engine()
            .body_longitude_crossings_in_range(
                bid,
                TimeInput::from_ut1_jd(start_jd_ut1),
                TimeInput::from_ut1_jd(end_jd_ut1),
                target_lon_deg,
                &opts,
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(list
            .into_iter()
            .map(JsLongitudeCrossingEvent::from)
            .collect())
    }

    /// Finds all Chinese 24 Solar Terms (*Jieqi*) in `[start_jd_ut1, end_jd_ut1]` (`R12`).
    pub fn solar_terms_in_range_core(
        &self,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
    ) -> Result<Vec<JsSolarTermEvent>, String> {
        let list = self
            .engine()
            .solar_terms_in_range(
                TimeInput::from_ut1_jd(start_jd_ut1),
                TimeInput::from_ut1_jd(end_jd_ut1),
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(list.into_iter().map(JsSolarTermEvent::from).collect())
    }

    /// Computes all 4 cardinal astronomical seasons (`0°`, `90°`, `180°`, `270°` apparent solar longitude) for `year`.
    pub fn seasons_in_year_core(&self, year: i32) -> Result<JsSeasonsOfYear, String> {
        let s = self
            .engine()
            .seasons_in_year(year)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSeasonsOfYear::from(s))
    }

    /// Computes the exact Human Design `-88°` prenatal solar arc epoch (`R9.2`).
    pub fn human_design_epoch_core(
        &self,
        birth_jd_ut1: f64,
    ) -> Result<JsLongitudeCrossingEvent, String> {
        let ev = self
            .engine()
            .human_design_epoch(TimeInput::from_ut1_jd(birth_jd_ut1))
            .map_err(|e| format_js_error(&e))?;
        Ok(JsLongitudeCrossingEvent::from(ev))
    }

    /// Finds the next or previous instant from `start_jd_ut1` when lunar-solar elongation reaches `target_elongation_deg` (`R9.2`).
    pub fn search_lunar_phase_core(
        &self,
        start_jd_ut1: f64,
        target_elongation_deg: f64,
        forward: Option<bool>,
    ) -> Result<JsLunarPhaseEvent, String> {
        let ev = self
            .engine()
            .search_lunar_phase(
                TimeInput::from_ut1_jd(start_jd_ut1),
                target_elongation_deg,
                forward.unwrap_or(true),
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(JsLunarPhaseEvent::from(ev))
    }

    /// Returns all 4 principal lunar phases occurring within `[start_jd_ut1, end_jd_ut1]` (`R9.2`).
    pub fn lunar_phases_in_range_core(
        &self,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
    ) -> Result<Vec<JsLunarPhaseEvent>, String> {
        let list = self
            .engine()
            .lunar_phases_in_range(
                TimeInput::from_ut1_jd(start_jd_ut1),
                TimeInput::from_ut1_jd(end_jd_ut1),
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(list.into_iter().map(JsLunarPhaseEvent::from).collect())
    }

    /// Computes the complete lunar phase summary (including 8-phase and Kerykeion 28-phase indices) at `jd_ut1` (`R0`, `R2.2`, `R9.2`).
    pub fn lunar_phase_summary_core(&self, jd_ut1: f64) -> Result<JsLunarPhaseSummary, String> {
        let s = self
            .engine()
            .lunar_phase_summary(TimeInput::from_ut1_jd(jd_ut1))
            .map_err(|e| format_js_error(&e))?;
        Ok(JsLunarPhaseSummary::from(s))
    }

    /// Finds all retrograde and direct stations ($\dot{\lambda}(t) = 0$) of `body` within `[start_jd_ut1, end_jd_ut1]` (`R9.2`).
    pub fn search_stations_in_range_core(
        &self,
        body: &str,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<JsStationEvent>, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let list = self
            .engine()
            .search_stations_in_range(
                bid,
                TimeInput::from_ut1_jd(start_jd_ut1),
                TimeInput::from_ut1_jd(end_jd_ut1),
                &opts,
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(list.into_iter().map(JsStationEvent::from).collect())
    }

    /// Finds all northern and southern declination extrema ($\dot{\delta}(t) = 0$) of `body` within `[start_jd_ut1, end_jd_ut1]` (`R9.2`).
    pub fn search_declination_extrema_in_range_core(
        &self,
        body: &str,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<JsDeclinationExtremumEvent>, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let list = self
            .engine()
            .search_declination_extrema_in_range(
                bid,
                TimeInput::from_ut1_jd(start_jd_ut1),
                TimeInput::from_ut1_jd(end_jd_ut1),
                &opts,
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(list
            .into_iter()
            .map(JsDeclinationExtremumEvent::from)
            .collect())
    }

    /// Computes Bernadette Brady diurnal fixed-star parans (`R6.2`) for `observer` starting at `jd_ut1`.
    pub fn find_parans_core(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
        orb_minutes: Option<f64>,
        bodies: Option<Vec<String>>,
        star_names: Option<Vec<String>>,
    ) -> Result<Vec<JsParanEvent>, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let mut paran_opts = ParanOptions::default();
        if let Some(orb) = orb_minutes {
            paran_opts.orb_minutes = orb;
        }
        if let Some(b_list) = bodies {
            paran_opts.bodies = b_list
                .iter()
                .map(|b| parse_body_id(b))
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| format_js_error(&e))?;
        }
        paran_opts.star_names = star_names;
        let list = self
            .engine()
            .find_parans(TimeInput::from_ut1_jd(jd_ut1), obs, &paran_opts)
            .map_err(|e| format_js_error(&e))?;
        Ok(list.into_iter().map(JsParanEvent::from).collect())
    }

    /// Finds the next planetary return (`Solar Return`, `Lunar Return`, `Saturn Return`, etc.) after `search_after_jd_ut1`.
    pub fn planetary_return_core(
        &self,
        body: &str,
        natal_jd_ut1: f64,
        search_after_jd_ut1: f64,
        options: Option<JsCalcOptions>,
    ) -> Result<serde_json::Value, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let ev = self
            .engine()
            .planetary_return(
                bid,
                TimeInput::from_ut1_jd(natal_jd_ut1),
                TimeInput::from_ut1_jd(search_after_jd_ut1),
                opts,
            )
            .map_err(|e| format_js_error(&e))?;
        to_js_value(&ev)
    }

    /// Finds all planetary returns of `body` in `[start_jd_ut1, end_jd_ut1]`, including 3-pass retrograde returns.
    pub fn planetary_returns_in_range_core(
        &self,
        body: &str,
        natal_jd_ut1: f64,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<serde_json::Value>, String> {
        let bid = parse_body_id(body).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let list = self
            .engine()
            .planetary_returns_in_range(
                bid,
                TimeInput::from_ut1_jd(natal_jd_ut1),
                TimeInput::from_ut1_jd(start_jd_ut1),
                TimeInput::from_ut1_jd(end_jd_ut1),
                opts,
            )
            .map_err(|e| format_js_error(&e))?;
        to_js_value_vec(&list)
    }

    /// Computes mutual natal/mundane ecliptic aspects among `bodies` at `jd_ut1`.
    pub fn calc_chart_aspects_core(
        &self,
        jd_ut1: f64,
        bodies: Option<Vec<String>>,
        preset: Option<String>,
        orb_scale: Option<f64>,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<serde_json::Value>, String> {
        let body_ids =
            parse_body_list_or_default(bodies.as_deref()).map_err(|e| format_js_error(&e))?;
        let pr = parse_aspect_preset(preset.as_deref()).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let list = self
            .engine()
            .calc_chart_aspects(
                TimeInput::from_ut1_jd(jd_ut1),
                &body_ids,
                pr,
                orb_scale.unwrap_or(1.0),
                opts,
            )
            .map_err(|e| format_js_error(&e))?;
        to_js_value_vec(&list)
    }

    /// Computes bipartite synastry or transit-to-natal ecliptic aspects between Chart A and Chart B.
    #[allow(clippy::too_many_arguments)]
    pub fn calc_synastry_aspects_core(
        &self,
        jd_ut1_a: f64,
        bodies_a: Vec<String>,
        jd_ut1_b: f64,
        bodies_b: Vec<String>,
        preset: Option<String>,
        orb_scale: Option<f64>,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<serde_json::Value>, String> {
        let bids_a =
            parse_body_list_or_default(Some(&bodies_a)).map_err(|e| format_js_error(&e))?;
        let bids_b =
            parse_body_list_or_default(Some(&bodies_b)).map_err(|e| format_js_error(&e))?;
        let pr = parse_aspect_preset(preset.as_deref()).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let list = self
            .engine()
            .calc_synastry_aspects(
                TimeInput::from_ut1_jd(jd_ut1_a),
                &bids_a,
                TimeInput::from_ut1_jd(jd_ut1_b),
                &bids_b,
                pr,
                orb_scale.unwrap_or(1.0),
                opts,
            )
            .map_err(|e| format_js_error(&e))?;
        to_js_value_vec(&list)
    }

    /// Computes Parallel and ContraParallel declination aspects among `bodies` at `jd_ut1`.
    pub fn calc_declination_aspects_core(
        &self,
        jd_ut1: f64,
        bodies: Option<Vec<String>>,
        max_orb_deg: Option<f64>,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<serde_json::Value>, String> {
        let body_ids =
            parse_body_list_or_default(bodies.as_deref()).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let list = self
            .engine()
            .calc_declination_aspects(
                TimeInput::from_ut1_jd(jd_ut1),
                &body_ids,
                max_orb_deg.unwrap_or(1.2),
                opts,
            )
            .map_err(|e| format_js_error(&e))?;
        to_js_value_vec(&list)
    }

    /// Searches forward or backward from `start_jd_ut1` for the exact aspect `target_angle_deg` between `body1` and `body2`.
    pub fn search_exact_aspect_core(
        &self,
        body1: &str,
        body2: &str,
        target_angle_deg: f64,
        start_jd_ut1: f64,
        forward: Option<bool>,
        options: Option<JsCalcOptions>,
    ) -> Result<serde_json::Value, String> {
        let b1 = parse_body_id(body1).map_err(|e| format_js_error(&e))?;
        let b2 = parse_body_id(body2).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let ev = self
            .engine()
            .search_exact_aspect(
                b1,
                b2,
                target_angle_deg,
                TimeInput::from_ut1_jd(start_jd_ut1),
                forward.unwrap_or(true),
                opts,
            )
            .map_err(|e| format_js_error(&e))?;
        to_js_value(&ev)
    }

    /// Finds all exact mutual sky-to-sky aspects among `bodies` in `[start_jd_ut1, end_jd_ut1]`.
    pub fn mutual_aspects_in_range_core(
        &self,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        bodies: Option<Vec<String>>,
        preset: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<serde_json::Value>, String> {
        let body_ids =
            parse_body_list_or_default(bodies.as_deref()).map_err(|e| format_js_error(&e))?;
        let pr = parse_aspect_preset(preset.as_deref()).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let list = self
            .engine()
            .mutual_aspects_in_range(
                &body_ids,
                pr,
                TimeInput::from_ut1_jd(start_jd_ut1),
                TimeInput::from_ut1_jd(end_jd_ut1),
                opts,
            )
            .map_err(|e| format_js_error(&e))?;
        to_js_value_vec(&list)
    }

    /// Finds all exact transit-to-natal aspect hits in `[start_jd_ut1, end_jd_ut1]`, including multi-pass retrograde cycles.
    #[allow(clippy::too_many_arguments)]
    pub fn natal_transits_in_range_core(
        &self,
        natal_jd_ut1: f64,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        transiting_bodies: Option<Vec<String>>,
        natal_bodies: Option<Vec<String>>,
        preset: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<serde_json::Value>, String> {
        let t_bodies = parse_body_list_or_default(transiting_bodies.as_deref())
            .map_err(|e| format_js_error(&e))?;
        let n_bodies =
            parse_body_list_or_default(natal_bodies.as_deref()).map_err(|e| format_js_error(&e))?;
        let pr = parse_aspect_preset(preset.as_deref()).map_err(|e| format_js_error(&e))?;
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let list = self
            .engine()
            .natal_transits_in_range(
                &t_bodies,
                TimeInput::from_ut1_jd(natal_jd_ut1),
                &n_bodies,
                pr,
                TimeInput::from_ut1_jd(start_jd_ut1),
                TimeInput::from_ut1_jd(end_jd_ut1),
                opts,
            )
            .map_err(|e| format_js_error(&e))?;
        to_js_value_vec(&list)
    }

    /// Finds all Void-of-Course Moon periods in `[start_jd_ut1, end_jd_ut1]`.
    pub fn void_of_course_moon_in_range_core(
        &self,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<JsCalcOptions>,
    ) -> Result<Vec<serde_json::Value>, String> {
        let opts = match options {
            Some(o) => o.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        let list = self
            .engine()
            .void_of_course_moon_in_range(
                TimeInput::from_ut1_jd(start_jd_ut1),
                TimeInput::from_ut1_jd(end_jd_ut1),
                None,
                opts,
            )
            .map_err(|e| format_js_error(&e))?;
        to_js_value_vec(&list)
    }

    /// Computes the complete Vedic Panchanga (Vara, Tithi, Nakshatra, Yoga, Karana) at `jd_ut1` for `observer`.
    pub fn panchanga_core(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
        ayanamsa: Option<String>,
    ) -> Result<serde_json::Value, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let ayan = AyanamsaModel::from_code(ayanamsa.as_deref().unwrap_or("lahiri"))
            .map_err(|e| format_js_error(&e))?;
        let p = self
            .engine()
            .panchanga(TimeInput::from_ut1_jd(jd_ut1), obs, ayan)
            .map_err(|e| format_js_error(&e))?;
        to_js_value(&p)
    }

    /// Computes the 120-year Vimshottari Dasha tree (9 Mahadashas × 9 Antardashas) from `birth_jd_ut1`.
    pub fn vimshottari_dasha_core(
        &self,
        birth_jd_ut1: f64,
        ayanamsa: Option<String>,
        year_length_days: Option<f64>,
    ) -> Result<serde_json::Value, String> {
        let ayan = AyanamsaModel::from_code(ayanamsa.as_deref().unwrap_or("lahiri"))
            .map_err(|e| format_js_error(&e))?;
        let tree = self
            .engine()
            .vimshottari_dasha(TimeInput::from_ut1_jd(birth_jd_ut1), ayan, year_length_days)
            .map_err(|e| format_js_error(&e))?;
        to_js_value(&tree)
    }

    /// Computes a single Parashara divisional chart (Varga) placement (`D1`, `D2`, `D3`, `D4`, `D7`, `D9`, `D10`, `D12`, `D16`, `D20`, `D24`, `D27`, `D30`, `D40`, `D45`, `D60`).
    pub fn calc_varga_core(
        &self,
        sidereal_lon_deg: f64,
        division: u32,
    ) -> Result<serde_json::Value, String> {
        let div_u16 = u16::try_from(division).map_err(|_| {
            format_js_error(&EphemerisError::InvalidInput {
                message: format!("Invalid varga division: {division}"),
            })
        })?;
        let v = self
            .engine()
            .calc_varga(sidereal_lon_deg, div_u16)
            .map_err(|e| format_js_error(&e))?;
        to_js_value(&v)
    }

    /// Computes all 16 classical Parashara Shodashavarga placements (`D1` through `D60`).
    pub fn calc_shodashavarga_core(
        &self,
        sidereal_lon_deg: f64,
    ) -> Result<Vec<serde_json::Value>, String> {
        let list = self
            .engine()
            .calc_shodashavarga(sidereal_lon_deg)
            .map_err(|e| format_js_error(&e))?;
        to_js_value_vec(&list)
    }

    /// Evaluates diurnal/nocturnal chart sect (`is_day_chart`) and true solar altitude at `jd_ut1` for `observer`.
    pub fn chart_sect_core(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
    ) -> Result<serde_json::Value, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let sect = self
            .engine()
            .chart_sect(TimeInput::from_ut1_jd(jd_ut1), obs)
            .map_err(|e| format_js_error(&e))?;
        to_js_value(&sect)
    }

    /// Computes the 24 unequal Chaldean planetary hours from preceding sunrise to next sunrise at `jd_ut1` for `observer`.
    pub fn planetary_hours_core(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
    ) -> Result<serde_json::Value, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let ph = self
            .engine()
            .planetary_hours(TimeInput::from_ut1_jd(jd_ut1), obs)
            .map_err(|e| format_js_error(&e))?;
        to_js_value(&ph)
    }

    /// Computes the 7 classical Hermetic / Arabic Lots (Fortune, Spirit, Eros, Victory, Necessity, Courage, Nemesis) with sect reversal and house placement.
    pub fn calc_lots_core(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
        house_options: Option<JsHouseOptions>,
    ) -> Result<Vec<serde_json::Value>, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let h_opts = match house_options {
            Some(ho) => ho.to_core().map_err(|e| format_js_error(&e))?,
            None => HouseOptions::default(),
        };
        let list = self
            .engine()
            .calc_lots(TimeInput::from_ut1_jd(jd_ut1), obs, h_opts)
            .map_err(|e| format_js_error(&e))?;
        to_js_value_vec(&list)
    }

    /// Computes a single-call complete Astrological & Human Design Chart Snapshot at `jd_ut1` for `observer`.
    #[allow(clippy::too_many_arguments)]
    pub fn chart_snapshot_core(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
        bodies: Option<Vec<String>>,
        house_options: Option<JsHouseOptions>,
        calc_options: Option<JsCalcOptions>,
        preset: Option<String>,
        orb_scale: Option<f64>,
    ) -> Result<serde_json::Value, String> {
        let obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m.unwrap_or(0.0),
        )
        .map_err(|e| format_js_error(&e))?;
        let body_ids =
            parse_body_list_or_default(bodies.as_deref()).map_err(|e| format_js_error(&e))?;
        let mut h_opts = match house_options {
            Some(ho) => ho.to_core().map_err(|e| format_js_error(&e))?,
            None => HouseOptions::default(),
        };
        let mut c_opts = match calc_options {
            Some(co) => co.to_core().map_err(|e| format_js_error(&e))?,
            None => CalcOptions::default(),
        };
        if matches!(h_opts.zodiac, Zodiac::Tropical) && !matches!(c_opts.zodiac, Zodiac::Tropical) {
            h_opts.zodiac = c_opts.zodiac;
        } else if matches!(c_opts.zodiac, Zodiac::Tropical)
            && !matches!(h_opts.zodiac, Zodiac::Tropical)
        {
            c_opts.zodiac = h_opts.zodiac;
        }
        let pr = parse_aspect_preset(preset.as_deref()).map_err(|e| format_js_error(&e))?;
        let snap = self
            .engine()
            .chart_snapshot(
                TimeInput::from_ut1_jd(jd_ut1),
                obs,
                &body_ids,
                h_opts,
                c_opts,
                pr,
                orb_scale.unwrap_or(1.0),
            )
            .map_err(|e| format_js_error(&e))?;
        to_js_value(&snap)
    }

    /// Maps an ecliptic longitude (`lon_deg`) into Zodiac Sign, Decan, Nakshatra (27 × 4 Pada), and Human Design 64 Gate / 6 Line / Color / Tone / Base.
    pub fn map_zodiac_placement_core(&self, lon_deg: f64) -> Result<serde_json::Value, String> {
        to_js_value(&core_map_zodiac_placement(lon_deg))
    }
}

#[cfg(feature = "napi-addon")]
#[napi]
impl JsEngine {
    /// Zero-threshold constructor: if `dataDir` is provided and non-empty, loads `fromDataDir`;
    /// otherwise checks `LXEPHEMERIS_DATA_DIR` / `SE_EPHE_PATH` or falls back to embedded Moshier (`new Engine()`).
    #[napi(constructor)]
    pub fn new(
        data_dir: Option<String>,
        max_threads: Option<u32>,
        pin_time_table_version: Option<String>,
        verify_manifest_sha256: Option<bool>,
    ) -> napi::Result<Self> {
        Self::new_core(
            data_dir,
            max_threads,
            pin_time_table_version,
            verify_manifest_sha256,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi(factory)]
    pub fn from_data_dir(
        data_dir: String,
        max_threads: Option<u32>,
        pin_time_table_version: Option<String>,
        verify_manifest_sha256: Option<bool>,
    ) -> napi::Result<Self> {
        Self::from_data_dir_core(
            data_dir,
            max_threads,
            pin_time_table_version,
            verify_manifest_sha256,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi(factory)]
    pub fn from_bsp_paths(
        paths: Vec<String>,
        max_threads: Option<u32>,
        pin_time_table_version: Option<String>,
    ) -> napi::Result<Self> {
        Self::from_bsp_paths_core(paths, max_threads, pin_time_table_version)
            .map_err(napi::Error::from_reason)
    }

    #[napi(factory)]
    pub fn from_bytes(kernels: Vec<JsKernelBytes>, max_threads: Option<u32>) -> napi::Result<Self> {
        let pairs: Vec<(String, Vec<u8>)> = kernels
            .into_iter()
            .map(|k| (k.name, k.data.to_vec()))
            .collect();
        Self::from_bytes_core(&pairs, max_threads).map_err(napi::Error::from_reason)
    }

    #[napi(factory)]
    pub fn new_moshier(
        max_threads: Option<u32>,
        pin_time_table_version: Option<String>,
    ) -> napi::Result<Self> {
        Ok(Self::new_moshier_core(max_threads, pin_time_table_version))
    }

    #[napi(factory)]
    pub fn from_data_dir_or_moshier(
        data_dir: String,
        max_threads: Option<u32>,
        pin_time_table_version: Option<String>,
        verify_manifest_sha256: Option<bool>,
    ) -> napi::Result<Self> {
        Self::from_data_dir_or_moshier_core(
            data_dir,
            max_threads,
            pin_time_table_version,
            verify_manifest_sha256,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn capabilities_json(&self) -> napi::Result<String> {
        self.capabilities_json_core()
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn license_info(&self) -> napi::Result<serde_json::Value> {
        self.license_info_core()
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn coverage_json(&self, body: String) -> napi::Result<String> {
        self.coverage_json_core(&body)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn delta_t(&self, jd: JsSplitJd, scale: Option<String>) -> napi::Result<JsDeltaTInfo> {
        self.delta_t_core(jd, scale)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn sidtime(
        &self,
        jd: JsSplitJd,
        scale: Option<String>,
        nutation: Option<String>,
        longitude_deg: Option<f64>,
    ) -> napi::Result<JsSiderealTimeInfo> {
        self.sidtime_core(jd, scale, nutation, longitude_deg)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn ecl_nut(&self, jd: JsSplitJd, scale: Option<String>) -> napi::Result<JsEclNutInfo> {
        self.ecl_nut_core(jd, scale)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn star_info(&self, name_or_bayer: String) -> napi::Result<JsStarMetadata> {
        self.star_info_core(&name_or_bayer)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn list_stars(&self) -> Vec<JsStarMetadata> {
        self.list_stars_core()
    }

    #[napi]
    pub fn calc(
        &self,
        body: String,
        jd: JsSplitJd,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<JsCalcResult> {
        self.calc_core(&body, jd, scale, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn calc_star_by_name(
        &self,
        name_or_bayer: String,
        jd: JsSplitJd,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<JsCalcResult> {
        self.calc_star_by_name_core(&name_or_bayer, jd, scale, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn calc_batch(
        &self,
        bodies: Vec<String>,
        epochs: Vec<JsSplitJd>,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<Vec<JsCalcResult>>> {
        self.calc_batch_core(&bodies, &epochs, scale, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn calc_batch_flat_f64(
        &self,
        bodies: Vec<String>,
        epochs: Vec<JsSplitJd>,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<napi::bindgen_prelude::Float64Array> {
        let vec = self
            .calc_batch_flat_f64_core(&bodies, &epochs, scale, options)
            .map_err(napi::Error::from_reason)?;
        Ok(napi::bindgen_prelude::Float64Array::new(vec))
    }

    #[napi]
    pub fn calc_batch_per_item(
        &self,
        bodies: Vec<String>,
        epochs: Vec<JsSplitJd>,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<Vec<JsBatchItemResult>>> {
        self.calc_batch_per_item_core(&bodies, &epochs, scale, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn ayanamsa(
        &self,
        model: String,
        jd: JsSplitJd,
        scale: Option<String>,
        nutation: Option<String>,
    ) -> napi::Result<JsAyanamsaValue> {
        self.ayanamsa_core(&model, jd, scale, nutation)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn houses(
        &self,
        jd: JsSplitJd,
        observer: JsObserver,
        system: Option<String>,
        scale: Option<String>,
        polar_policy: Option<String>,
        zodiac: Option<String>,
        ayanamsa: Option<String>,
        nutation: Option<String>,
        include_speed: Option<bool>,
    ) -> napi::Result<JsHouseResult> {
        self.houses_core(
            jd,
            observer,
            system,
            scale,
            polar_policy,
            zodiac,
            ayanamsa,
            nutation,
            include_speed,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn houses_batch(
        &self,
        epochs: Vec<JsSplitJd>,
        observer: JsObserver,
        system: Option<String>,
        scale: Option<String>,
        polar_policy: Option<String>,
        zodiac: Option<String>,
        ayanamsa: Option<String>,
        nutation: Option<String>,
        include_speed: Option<bool>,
    ) -> napi::Result<Vec<JsHouseResult>> {
        self.houses_batch_core(
            &epochs,
            observer,
            system,
            scale,
            polar_policy,
            zodiac,
            ayanamsa,
            nutation,
            include_speed,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn azalt(
        &self,
        ra_deg: f64,
        dec_deg: f64,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        pressure_hpa: Option<f64>,
        temperature_c: Option<f64>,
    ) -> napi::Result<JsHorizontalCoords> {
        self.azalt_core(
            ra_deg,
            dec_deg,
            jd,
            observer,
            scale,
            pressure_hpa,
            temperature_c,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn azalt_ecliptic(
        &self,
        ecl_lon_deg: f64,
        ecl_lat_deg: f64,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        pressure_hpa: Option<f64>,
        temperature_c: Option<f64>,
    ) -> napi::Result<JsHorizontalCoords> {
        self.azalt_ecliptic_core(
            ecl_lon_deg,
            ecl_lat_deg,
            jd,
            observer,
            scale,
            pressure_hpa,
            temperature_c,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn calc_azalt(
        &self,
        body: String,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        pressure_hpa: Option<f64>,
        temperature_c: Option<f64>,
    ) -> napi::Result<JsHorizontalCoords> {
        self.calc_azalt_core(&body, jd, observer, scale, pressure_hpa, temperature_c)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans(
        &self,
        body: String,
        event: String,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        custom_horizon_alt_deg: Option<f64>,
        max_search_days: Option<f64>,
        disc_mode: Option<String>,
        pressure_hpa: Option<f64>,
        temperature_c: Option<f64>,
    ) -> napi::Result<Option<f64>> {
        self.rise_trans_core(
            &body,
            &event,
            jd,
            observer,
            scale,
            custom_horizon_alt_deg,
            max_search_days,
            disc_mode,
            pressure_hpa,
            temperature_c,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn heliacal_phenomena(
        &self,
        body: String,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        extinction_k_v: Option<f64>,
        snellen_factor: Option<f64>,
        pressure_hpa: Option<f64>,
        temperature_c: Option<f64>,
        relative_humidity: Option<f64>,
        max_search_days: Option<f64>,
        observer_age_years: Option<f64>,
        snellen_ratio: Option<f64>,
    ) -> napi::Result<JsHeliacalResult> {
        self.heliacal_phenomena_core(
            &body,
            jd,
            observer,
            scale,
            extinction_k_v,
            snellen_factor,
            pressure_hpa,
            temperature_c,
            relative_humidity,
            max_search_days,
            observer_age_years,
            snellen_ratio,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn next_solar_eclipse(
        &self,
        jd: JsSplitJd,
        scale: Option<String>,
        backward: Option<bool>,
        observer: Option<JsObserver>,
        max_synodic_months: Option<u32>,
    ) -> napi::Result<JsSolarEclipseInfo> {
        self.next_solar_eclipse_core(jd, scale, backward, observer, max_synodic_months)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn prev_solar_eclipse(
        &self,
        jd: JsSplitJd,
        scale: Option<String>,
        observer: Option<JsObserver>,
        max_synodic_months: Option<u32>,
    ) -> napi::Result<JsSolarEclipseInfo> {
        self.next_solar_eclipse_core(jd, scale, Some(true), observer, max_synodic_months)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn next_local_solar_eclipse(
        &self,
        jd: JsSplitJd,
        observer: JsObserver,
        scale: Option<String>,
        backward: Option<bool>,
        max_synodic_months: Option<u32>,
    ) -> napi::Result<JsSolarEclipseInfo> {
        self.next_local_solar_eclipse_core(jd, observer, scale, backward, max_synodic_months)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn solar_eclipses_in_range(
        &self,
        start_jd_tt: f64,
        end_jd_tt: f64,
        observer: Option<JsObserver>,
    ) -> napi::Result<Vec<JsSolarEclipseInfo>> {
        self.solar_eclipses_in_range_core(start_jd_tt, end_jd_tt, observer)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn next_lunar_eclipse(
        &self,
        jd: JsSplitJd,
        scale: Option<String>,
        backward: Option<bool>,
        shadow_rule: Option<String>,
        max_synodic_months: Option<u32>,
    ) -> napi::Result<JsLunarEclipseInfo> {
        self.next_lunar_eclipse_core(jd, scale, backward, shadow_rule, max_synodic_months)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn prev_lunar_eclipse(
        &self,
        jd: JsSplitJd,
        scale: Option<String>,
        shadow_rule: Option<String>,
        max_synodic_months: Option<u32>,
    ) -> napi::Result<JsLunarEclipseInfo> {
        self.next_lunar_eclipse_core(jd, scale, Some(true), shadow_rule, max_synodic_months)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn lunar_eclipses_in_range(
        &self,
        start_jd_tt: f64,
        end_jd_tt: f64,
        shadow_rule: Option<String>,
    ) -> napi::Result<Vec<JsLunarEclipseInfo>> {
        self.lunar_eclipses_in_range_core(start_jd_tt, end_jd_tt, shadow_rule)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn gauquelin_sectors(
        &self,
        armc_deg: f64,
        eps_true_deg: f64,
        latitude_deg: f64,
        polar_policy: Option<String>,
    ) -> napi::Result<Vec<f64>> {
        self.gauquelin_sectors_core(armc_deg, eps_true_deg, latitude_deg, polar_policy)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn is_moshier_mode(&self) -> bool {
        self.is_moshier_mode_core()
    }

    #[napi]
    pub fn calc_moshier(
        &self,
        body: String,
        jd: JsSplitJd,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<JsCalcResult> {
        self.calc_moshier_core(&body, jd, scale, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn load_custom_stars(&self, path: String) -> napi::Result<u32> {
        self.load_custom_stars_core(path)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn phenomena(
        &self,
        body: String,
        jd: JsSplitJd,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<JsPhenomenaResult> {
        self.phenomena_core(&body, jd, scale, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn nod_aps(
        &self,
        body: String,
        jd: JsSplitJd,
        scale: Option<String>,
        options: Option<JsCalcOptions>,
        method_bits: Option<u32>,
    ) -> napi::Result<JsNodApsResult> {
        self.nod_aps_core(&body, jd, scale, options, method_bits)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn next_lunar_occultation(
        &self,
        target: String,
        jd: JsSplitJd,
        scale: Option<String>,
        backward: Option<bool>,
        max_search_days: Option<f64>,
    ) -> napi::Result<JsLunarOccultationInfo> {
        self.next_lunar_occultation_core(&target, jd, scale, backward, max_search_days)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn next_local_lunar_occultation(
        &self,
        target: String,
        observer: JsObserver,
        jd: JsSplitJd,
        scale: Option<String>,
        backward: Option<bool>,
        max_search_days: Option<f64>,
    ) -> napi::Result<JsLunarOccultationInfo> {
        self.next_local_lunar_occultation_core(
            &target,
            observer,
            jd,
            scale,
            backward,
            max_search_days,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn lunar_occultations_in_range(
        &self,
        target: String,
        start_jd: JsSplitJd,
        end_jd: JsSplitJd,
        scale: Option<String>,
    ) -> napi::Result<Vec<JsLunarOccultationInfo>> {
        self.lunar_occultations_in_range_core(&target, start_jd, end_jd, scale)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn local_lunar_occultations_in_range(
        &self,
        target: String,
        observer: JsObserver,
        start_jd: JsSplitJd,
        end_jd: JsSplitJd,
        scale: Option<String>,
    ) -> napi::Result<Vec<JsLunarOccultationInfo>> {
        self.local_lunar_occultations_in_range_core(&target, observer, start_jd, end_jd, scale)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn house_pos(
        &self,
        armc_deg: f64,
        geolat_deg: f64,
        eps_deg: f64,
        system: String,
        ecl_lon_deg: f64,
        ecl_lat_deg: f64,
    ) -> napi::Result<f64> {
        self.house_pos_core(
            armc_deg,
            geolat_deg,
            eps_deg,
            &system,
            ecl_lon_deg,
            ecl_lat_deg,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn time_equation(&self, jd: JsSplitJd, scale: Option<String>) -> napi::Result<f64> {
        self.time_equation_core(jd, scale)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn azalt_rev(
        &self,
        observer: JsObserver,
        jd: JsSplitJd,
        scale: Option<String>,
        azimuth_north_east_deg: f64,
        altitude_true_deg: f64,
        target_plane: Option<String>,
    ) -> napi::Result<Vec<f64>> {
        self.azalt_rev_core(
            observer,
            jd,
            scale,
            azimuth_north_east_deg,
            altitude_true_deg,
            target_plane,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn search_body_longitude(
        &self,
        body: String,
        jd_ut1: f64,
        target_lon_deg: f64,
        forward: Option<bool>,
        max_search_days: Option<f64>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<JsLongitudeCrossingEvent> {
        self.search_body_longitude_core(
            &body,
            jd_ut1,
            target_lon_deg,
            forward,
            max_search_days,
            options,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn body_longitude_crossings_in_range(
        &self,
        body: String,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        target_lon_deg: f64,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<JsLongitudeCrossingEvent>> {
        self.body_longitude_crossings_in_range_core(
            &body,
            start_jd_ut1,
            end_jd_ut1,
            target_lon_deg,
            options,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn solar_terms_in_range(
        &self,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
    ) -> napi::Result<Vec<JsSolarTermEvent>> {
        self.solar_terms_in_range_core(start_jd_ut1, end_jd_ut1)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn seasons_in_year(&self, year: i32) -> napi::Result<JsSeasonsOfYear> {
        self.seasons_in_year_core(year)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn human_design_epoch(&self, birth_jd_ut1: f64) -> napi::Result<JsLongitudeCrossingEvent> {
        self.human_design_epoch_core(birth_jd_ut1)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn search_lunar_phase(
        &self,
        start_jd_ut1: f64,
        target_elongation_deg: f64,
        forward: Option<bool>,
    ) -> napi::Result<JsLunarPhaseEvent> {
        self.search_lunar_phase_core(start_jd_ut1, target_elongation_deg, forward)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn lunar_phases_in_range(
        &self,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
    ) -> napi::Result<Vec<JsLunarPhaseEvent>> {
        self.lunar_phases_in_range_core(start_jd_ut1, end_jd_ut1)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn lunar_phase_summary(&self, jd_ut1: f64) -> napi::Result<JsLunarPhaseSummary> {
        self.lunar_phase_summary_core(jd_ut1)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn search_stations_in_range(
        &self,
        body: String,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<JsStationEvent>> {
        self.search_stations_in_range_core(&body, start_jd_ut1, end_jd_ut1, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn search_declination_extrema_in_range(
        &self,
        body: String,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<JsDeclinationExtremumEvent>> {
        self.search_declination_extrema_in_range_core(&body, start_jd_ut1, end_jd_ut1, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn find_parans(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
        orb_minutes: Option<f64>,
        bodies: Option<Vec<String>>,
        star_names: Option<Vec<String>>,
    ) -> napi::Result<Vec<JsParanEvent>> {
        self.find_parans_core(jd_ut1, observer, orb_minutes, bodies, star_names)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn planetary_return(
        &self,
        body: String,
        natal_jd_ut1: f64,
        search_after_jd_ut1: f64,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<serde_json::Value> {
        self.planetary_return_core(&body, natal_jd_ut1, search_after_jd_ut1, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn planetary_returns_in_range(
        &self,
        body: String,
        natal_jd_ut1: f64,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<serde_json::Value>> {
        self.planetary_returns_in_range_core(
            &body,
            natal_jd_ut1,
            start_jd_ut1,
            end_jd_ut1,
            options,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn calc_chart_aspects(
        &self,
        jd_ut1: f64,
        bodies: Option<Vec<String>>,
        preset: Option<String>,
        orb_scale: Option<f64>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<serde_json::Value>> {
        self.calc_chart_aspects_core(jd_ut1, bodies, preset, orb_scale, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn calc_synastry_aspects(
        &self,
        jd_ut1_a: f64,
        bodies_a: Vec<String>,
        jd_ut1_b: f64,
        bodies_b: Vec<String>,
        preset: Option<String>,
        orb_scale: Option<f64>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<serde_json::Value>> {
        self.calc_synastry_aspects_core(
            jd_ut1_a, bodies_a, jd_ut1_b, bodies_b, preset, orb_scale, options,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn calc_declination_aspects(
        &self,
        jd_ut1: f64,
        bodies: Option<Vec<String>>,
        max_orb_deg: Option<f64>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<serde_json::Value>> {
        self.calc_declination_aspects_core(jd_ut1, bodies, max_orb_deg, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn search_exact_aspect(
        &self,
        body1: String,
        body2: String,
        target_angle_deg: f64,
        start_jd_ut1: f64,
        forward: Option<bool>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<serde_json::Value> {
        self.search_exact_aspect_core(
            &body1,
            &body2,
            target_angle_deg,
            start_jd_ut1,
            forward,
            options,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn mutual_aspects_in_range(
        &self,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        bodies: Option<Vec<String>>,
        preset: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<serde_json::Value>> {
        self.mutual_aspects_in_range_core(start_jd_ut1, end_jd_ut1, bodies, preset, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn natal_transits_in_range(
        &self,
        natal_jd_ut1: f64,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        transiting_bodies: Option<Vec<String>>,
        natal_bodies: Option<Vec<String>>,
        preset: Option<String>,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<serde_json::Value>> {
        self.natal_transits_in_range_core(
            natal_jd_ut1,
            start_jd_ut1,
            end_jd_ut1,
            transiting_bodies,
            natal_bodies,
            preset,
            options,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn void_of_course_moon_in_range(
        &self,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<JsCalcOptions>,
    ) -> napi::Result<Vec<serde_json::Value>> {
        self.void_of_course_moon_in_range_core(start_jd_ut1, end_jd_ut1, options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn panchanga(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
        ayanamsa: Option<String>,
    ) -> napi::Result<serde_json::Value> {
        self.panchanga_core(jd_ut1, observer, ayanamsa)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn vimshottari_dasha(
        &self,
        birth_jd_ut1: f64,
        ayanamsa: Option<String>,
        year_length_days: Option<f64>,
    ) -> napi::Result<serde_json::Value> {
        self.vimshottari_dasha_core(birth_jd_ut1, ayanamsa, year_length_days)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn calc_varga(
        &self,
        sidereal_lon_deg: f64,
        division: u32,
    ) -> napi::Result<serde_json::Value> {
        self.calc_varga_core(sidereal_lon_deg, division)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn calc_shodashavarga(
        &self,
        sidereal_lon_deg: f64,
    ) -> napi::Result<Vec<serde_json::Value>> {
        self.calc_shodashavarga_core(sidereal_lon_deg)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn chart_sect(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
    ) -> napi::Result<serde_json::Value> {
        self.chart_sect_core(jd_ut1, observer)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn planetary_hours(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
    ) -> napi::Result<serde_json::Value> {
        self.planetary_hours_core(jd_ut1, observer)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn calc_lots(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
        house_options: Option<JsHouseOptions>,
    ) -> napi::Result<Vec<serde_json::Value>> {
        self.calc_lots_core(jd_ut1, observer, house_options)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn chart_snapshot(
        &self,
        jd_ut1: f64,
        observer: JsObserver,
        bodies: Option<Vec<String>>,
        house_options: Option<JsHouseOptions>,
        calc_options: Option<JsCalcOptions>,
        preset: Option<String>,
        orb_scale: Option<f64>,
    ) -> napi::Result<serde_json::Value> {
        self.chart_snapshot_core(
            jd_ut1,
            observer,
            bodies,
            house_options,
            calc_options,
            preset,
            orb_scale,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn map_zodiac_placement(&self, lon_deg: f64) -> napi::Result<serde_json::Value> {
        self.map_zodiac_placement_core(lon_deg)
            .map_err(napi::Error::from_reason)
    }
}

#[cfg(feature = "swe-compat")]
mod swe_node {
use super::*;

/// Swiss Ephemeris Migration Adapter (`JsSweAdapter`) for Node.js / TypeScript.
#[cfg_attr(feature = "napi-addon", napi)]
pub struct JsSweAdapter {
    inner: Mutex<SweAdapter>,
}

impl JsSweAdapter {
    pub fn from_engine_core(engine: &JsEngine) -> Self {
        Self {
            inner: Mutex::new(SweAdapter::new(engine.engine())),
        }
    }

    pub fn new_moshier_core() -> Self {
        let eng = JsEngine::new_moshier_core(None, None);
        Self::from_engine_core(&eng)
    }

    pub fn from_data_dir_core(data_dir: String) -> Result<Self, String> {
        let eng = JsEngine::from_data_dir_core(data_dir, None, None, None)?;
        Ok(Self::from_engine_core(&eng))
    }

    pub fn from_data_dir_or_moshier_core(data_dir: String) -> Result<Self, String> {
        let eng = JsEngine::from_data_dir_or_moshier_core(data_dir, None, None, None)?;
        Ok(Self::from_engine_core(&eng))
    }

    pub fn set_sid_mode_core(
        &self,
        sid_mode: u32,
        t0: Option<f64>,
        ayan_t0: Option<f64>,
    ) -> Result<(), String> {
        let mut g = self.inner.lock().map_err(|e| e.to_string())?;
        if sid_mode == SE_SIDM_USER || t0.is_some() || ayan_t0.is_some() {
            g.set_sid_mode_custom(sid_mode, t0.unwrap_or(0.0), ayan_t0.unwrap_or(0.0))
                .map_err(|e| format_js_error(&e))
        } else {
            g.set_sid_mode(sid_mode).map_err(|e| format_js_error(&e))
        }
    }

    pub fn set_topo_core(
        &self,
        lon_deg: f64,
        lat_deg: f64,
        alt_m: Option<f64>,
    ) -> Result<(), String> {
        let mut g = self.inner.lock().map_err(|e| e.to_string())?;
        g.set_topo(lon_deg, lat_deg, alt_m.unwrap_or(0.0))
            .map_err(|e| format_js_error(&e))
    }

    pub fn set_polar_policy_core(&self, policy: &str) -> Result<(), String> {
        let pol = parse_polar_policy(policy).map_err(|e| format_js_error(&e))?;
        let mut g = self.inner.lock().map_err(|e| e.to_string())?;
        g.set_polar_policy(pol);
        Ok(())
    }

    pub fn deltat_core(&self, jd_ut: f64) -> Result<f64, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        Ok(g.deltat(jd_ut))
    }

    pub fn deltat_sec_core(&self, jd_ut: f64) -> Result<f64, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        Ok(g.deltat_sec(jd_ut))
    }

    pub fn sidtime_core(&self, jd_ut: f64) -> Result<f64, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        Ok(g.sidtime(jd_ut))
    }

    pub fn sidtime0_core(&self, jd_ut: f64, eps_deg: f64, nut_lon_deg: f64) -> Result<f64, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        Ok(g.sidtime0(jd_ut, eps_deg, nut_lon_deg))
    }

    pub fn get_ayanamsa_ut_core(&self, jd_ut: f64) -> Result<f64, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        g.get_ayanamsa_ut(jd_ut).map_err(|e| format_js_error(&e))
    }

    pub fn get_ayanamsa_core(&self, jd_et: f64) -> Result<f64, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        g.get_ayanamsa(jd_et).map_err(|e| format_js_error(&e))
    }

    pub fn get_ayanamsa_ex_ut_core(
        &self,
        jd_ut: f64,
        iflag: Option<u32>,
    ) -> Result<JsSweAyanamsaOut, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (retflag, ayanamsa) = g
            .get_ayanamsa_ex_ut(jd_ut, iflag.unwrap_or(2))
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSweAyanamsaOut { retflag, ayanamsa })
    }

    pub fn get_ayanamsa_ex_core(
        &self,
        jd_et: f64,
        iflag: Option<u32>,
    ) -> Result<JsSweAyanamsaOut, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (retflag, ayanamsa) = g
            .get_ayanamsa_ex(jd_et, iflag.unwrap_or(2))
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSweAyanamsaOut { retflag, ayanamsa })
    }

    pub fn calc_ut_core(
        &self,
        jd_ut: f64,
        ipl: i32,
        iflag: Option<u32>,
    ) -> Result<JsSweCalcOut, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (xx, retflag) = g
            .calc_ut(jd_ut, ipl, iflag.unwrap_or(258))
            .map_err(|e| format_js_error(&e))?;
        let vec = xx.to_vec();
        Ok(JsSweCalcOut {
            xx: vec.clone(),
            values: vec,
            retflag,
        })
    }

    pub fn calc_core(
        &self,
        jd_et: f64,
        ipl: i32,
        iflag: Option<u32>,
    ) -> Result<JsSweCalcOut, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (xx, retflag) = g
            .calc(jd_et, ipl, iflag.unwrap_or(258))
            .map_err(|e| format_js_error(&e))?;
        let vec = xx.to_vec();
        Ok(JsSweCalcOut {
            xx: vec.clone(),
            values: vec,
            retflag,
        })
    }

    pub fn fixstar2_ut_core(
        &self,
        star: &str,
        jd_ut: f64,
        iflag: Option<u32>,
    ) -> Result<JsSweFixstarOut, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (xx, name, retflag) = g
            .fixstar2_ut(star, jd_ut, iflag.unwrap_or(258))
            .map_err(|e| format_js_error(&e))?;
        let vec = xx.to_vec();
        Ok(JsSweFixstarOut {
            xx: vec.clone(),
            values: vec,
            star_name: name,
            retflag,
        })
    }

    pub fn fixstar2_mag_core(&self, star: &str) -> Result<f64, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        g.fixstar2_mag(star).map_err(|e| format_js_error(&e))
    }

    pub fn houses_ex_core(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
        hsys: Option<String>,
        iflag: Option<u32>,
    ) -> Result<JsSweHousesOut, String> {
        let ch = hsys.as_deref().unwrap_or("P").chars().next().unwrap_or('P');
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (cusps, ascmc) = g
            .houses_ex(jd_ut, lat_deg, lon_deg, ch, iflag.unwrap_or(0))
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSweHousesOut {
            cusps: cusps.to_vec(),
            ascmc: ascmc.to_vec(),
        })
    }

    pub fn houses_core(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
        hsys: Option<String>,
    ) -> Result<JsSweHousesOut, String> {
        self.houses_ex_core(jd_ut, lat_deg, lon_deg, hsys, Some(0))
    }

    pub fn houses_armc_core(
        &self,
        armc_deg: f64,
        lat_deg: f64,
        eps_deg: f64,
        hsys: Option<String>,
        sun_decl_deg: Option<f64>,
    ) -> Result<JsSweHousesOut, String> {
        let ch = hsys.as_deref().unwrap_or("P").chars().next().unwrap_or('P');
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (cusps, ascmc) = g
            .houses_armc_ex(armc_deg, lat_deg, eps_deg, ch, sun_decl_deg)
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSweHousesOut {
            cusps: cusps.to_vec(),
            ascmc: ascmc.to_vec(),
        })
    }

    pub fn azalt_core(
        &self,
        jd_ut: f64,
        calc_flag: i32,
        geopos: &[f64],
        atpress: f64,
        attemp: f64,
        xin: &[f64],
    ) -> Result<Vec<f64>, String> {
        let (Some(&glon), Some(&glat), Some(&x0), Some(&x1)) =
            (geopos.first(), geopos.get(1), xin.first(), xin.get(1))
        else {
            return Err(format_js_error(&EphemerisError::InvalidInput {
                message: "geopos and xin must have at least 2 elements".to_owned(),
            }));
        };
        let gpos = [glon, glat, geopos.get(2).copied().unwrap_or(0.0)];
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let out = g
            .azalt(jd_ut, calc_flag, gpos, atpress, attemp, [x0, x1])
            .map_err(|e| format_js_error(&e))?;
        Ok(vec![
            out.azimuth_deg,
            out.altitude_true_deg,
            out.altitude_apparent_deg,
        ])
    }

    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans_core(
        &self,
        jd_ut: f64,
        ipl: i32,
        rsmi: i32,
        geopos: &[f64],
        atpress: Option<f64>,
        attemp: Option<f64>,
        starname: Option<String>,
        epheflag: Option<u32>,
    ) -> Result<Option<f64>, String> {
        let (Some(&glon), Some(&glat)) = (geopos.first(), geopos.get(1)) else {
            return Err(format_js_error(&EphemerisError::InvalidInput {
                message: "geopos must have at least 2 elements [lon, lat, alt]".to_owned(),
            }));
        };
        let galt = geopos.get(2).copied().unwrap_or(0.0);
        let gpos = [glon, glat, galt];
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let press = atpress.unwrap_or(1013.25);
        let temp = attemp.unwrap_or(15.0);
        g.rise_trans(
            jd_ut,
            ipl,
            starname.as_deref(),
            epheflag.unwrap_or(258),
            rsmi,
            gpos,
            press,
            temp,
        )
        .map_err(|e| format_js_error(&e))
    }

    pub fn pheno_ut_core(
        &self,
        jd_ut: f64,
        ipl: i32,
        iflag: Option<u32>,
    ) -> Result<JsSwePhenoOut, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (attr, retflag) = g
            .pheno_ut(jd_ut, ipl, iflag.unwrap_or(0))
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSwePhenoOut {
            attr: attr.to_vec(),
            retflag,
        })
    }

    pub fn pheno_core(
        &self,
        jd_et: f64,
        ipl: i32,
        iflag: Option<u32>,
    ) -> Result<JsSwePhenoOut, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (attr, retflag) = g
            .pheno(jd_et, ipl, iflag.unwrap_or(0))
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSwePhenoOut {
            attr: attr.to_vec(),
            retflag,
        })
    }

    pub fn nod_aps_ut_core(
        &self,
        jd_ut: f64,
        ipl: i32,
        iflag: Option<u32>,
        method: Option<u32>,
    ) -> Result<JsNodApsResult, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (xnasc, xndsc, xperi, xaphe, retflag) = g
            .nod_aps_ut(jd_ut, ipl, iflag.unwrap_or(0), method.unwrap_or(0))
            .map_err(|e| format_js_error(&e))?;
        Ok(JsNodApsResult {
            ascending_node: xnasc.to_vec(),
            descending_node: xndsc.to_vec(),
            periapsis: xperi.to_vec(),
            apoapsis: xaphe.to_vec(),
            retflag,
        })
    }

    pub fn nod_aps_core(
        &self,
        jd_et: f64,
        ipl: i32,
        iflag: Option<u32>,
        method: Option<u32>,
    ) -> Result<JsNodApsResult, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (xnasc, xndsc, xperi, xaphe, retflag) = g
            .nod_aps(jd_et, ipl, iflag.unwrap_or(0), method.unwrap_or(0))
            .map_err(|e| format_js_error(&e))?;
        Ok(JsNodApsResult {
            ascending_node: xnasc.to_vec(),
            descending_node: xndsc.to_vec(),
            periapsis: xperi.to_vec(),
            apoapsis: xaphe.to_vec(),
            retflag,
        })
    }

    pub fn lun_occult_when_glob_core(
        &self,
        jd_ut_start: f64,
        ipl: i32,
        starname: Option<String>,
        iflag: Option<u32>,
        ifltype: Option<i32>,
        backward: Option<bool>,
    ) -> Result<JsSweOccultGlobOut, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (retflag, tret) = g
            .lun_occult_when_glob_ex(
                jd_ut_start,
                ipl,
                starname.as_deref(),
                iflag.unwrap_or(0),
                ifltype.unwrap_or(0),
                backward.unwrap_or(false),
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSweOccultGlobOut {
            retflag,
            tret: tret.to_vec(),
        })
    }

    pub fn lun_occult_when_loc_core(
        &self,
        jd_ut_start: f64,
        ipl: i32,
        starname: Option<String>,
        iflag: Option<u32>,
        geopos: &[f64],
        backward: Option<bool>,
    ) -> Result<JsSweOccultLocOut, String> {
        let (Some(&p0), Some(&p1), Some(&p2)) = (geopos.first(), geopos.get(1), geopos.get(2))
        else {
            return Err("geopos must have at least 3 elements [lon, lat, alt]".to_string());
        };
        let pos = [p0, p1, p2];
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let (retflag, tret, attr) = g
            .lun_occult_when_loc(
                jd_ut_start,
                ipl,
                starname.as_deref(),
                iflag.unwrap_or(0),
                pos,
                backward.unwrap_or(false),
            )
            .map_err(|e| format_js_error(&e))?;
        Ok(JsSweOccultLocOut {
            retflag,
            tret: tret.to_vec(),
            attr: attr.to_vec(),
        })
    }

    pub fn house_pos_core(
        &self,
        armc_deg: f64,
        geolat_deg: f64,
        eps_deg: f64,
        xpin: &[f64],
        hsys: Option<String>,
    ) -> Result<f64, String> {
        let (Some(&x0), Some(&x1)) = (xpin.first(), xpin.get(1)) else {
            return Err("xpin must have at least 2 elements [lon, lat]".to_string());
        };
        let hsys_ch = hsys.as_deref().unwrap_or("P").chars().next().unwrap_or('P');
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        g.house_pos(armc_deg, geolat_deg, eps_deg, hsys_ch, [x0, x1])
            .map_err(|e| format_js_error(&e))
    }

    pub fn time_equ_core(&self, jd_et: f64) -> Result<f64, String> {
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        g.time_equ(jd_et).map_err(|e| format_js_error(&e))
    }

    pub fn azalt_rev_core(
        &self,
        jd_ut: f64,
        calc_flag: i32,
        geopos: &[f64],
        xin: &[f64],
    ) -> Result<Vec<f64>, String> {
        let (Some(&g0), Some(&g1), Some(&g2)) = (geopos.first(), geopos.get(1), geopos.get(2))
        else {
            return Err("geopos must have at least 3 elements [lon, lat, alt]".to_string());
        };
        let (Some(&x0), Some(&x1)) = (xin.first(), xin.get(1)) else {
            return Err("xin must have at least 2 elements [azimuth, altitude]".to_string());
        };
        let g = self.inner.lock().map_err(|e| e.to_string())?;
        let res = g
            .azalt_rev(jd_ut, calc_flag, [g0, g1, g2], [x0, x1])
            .map_err(|e| format_js_error(&e))?;
        Ok(res.to_vec())
    }
}

#[cfg(feature = "napi-addon")]
#[napi]
impl JsSweAdapter {
    #[napi(constructor)]
    pub fn new(engine: &JsEngine) -> Self {
        Self::from_engine_core(engine)
    }

    #[napi(factory)]
    pub fn from_engine(engine: &JsEngine) -> Self {
        Self::from_engine_core(engine)
    }

    #[napi(factory)]
    pub fn new_moshier() -> napi::Result<Self> {
        Ok(Self::new_moshier_core())
    }

    #[napi(factory)]
    pub fn from_data_dir(data_dir: String) -> napi::Result<Self> {
        Self::from_data_dir_core(data_dir).map_err(napi::Error::from_reason)
    }

    #[napi(factory)]
    pub fn from_data_dir_or_moshier(data_dir: String) -> napi::Result<Self> {
        Self::from_data_dir_or_moshier_core(data_dir).map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn set_sid_mode(
        &self,
        sid_mode: u32,
        t0: Option<f64>,
        ayan_t0: Option<f64>,
    ) -> napi::Result<()> {
        self.set_sid_mode_core(sid_mode, t0, ayan_t0)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn set_topo(&self, lon_deg: f64, lat_deg: f64, alt_m: Option<f64>) -> napi::Result<()> {
        self.set_topo_core(lon_deg, lat_deg, alt_m)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn set_polar_policy(&self, policy: String) -> napi::Result<()> {
        self.set_polar_policy_core(&policy)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn deltat(&self, jd_ut: f64) -> napi::Result<f64> {
        self.deltat_core(jd_ut).map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn deltat_sec(&self, jd_ut: f64) -> napi::Result<f64> {
        self.deltat_sec_core(jd_ut)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn sidtime(&self, jd_ut: f64) -> napi::Result<f64> {
        self.sidtime_core(jd_ut).map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn sidtime0(&self, jd_ut: f64, eps_deg: f64, nut_lon_deg: f64) -> napi::Result<f64> {
        self.sidtime0_core(jd_ut, eps_deg, nut_lon_deg)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn get_ayanamsa_ut(&self, jd_ut: f64) -> napi::Result<f64> {
        self.get_ayanamsa_ut_core(jd_ut)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn get_ayanamsa(&self, jd_et: f64) -> napi::Result<f64> {
        self.get_ayanamsa_core(jd_et)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn get_ayanamsa_ex_ut(
        &self,
        jd_ut: f64,
        iflag: Option<u32>,
    ) -> napi::Result<JsSweAyanamsaOut> {
        self.get_ayanamsa_ex_ut_core(jd_ut, iflag)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn get_ayanamsa_ex(
        &self,
        jd_et: f64,
        iflag: Option<u32>,
    ) -> napi::Result<JsSweAyanamsaOut> {
        self.get_ayanamsa_ex_core(jd_et, iflag)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn calc_ut(&self, jd_ut: f64, ipl: i32, iflag: Option<u32>) -> napi::Result<JsSweCalcOut> {
        self.calc_ut_core(jd_ut, ipl, iflag)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn calc(&self, jd_et: f64, ipl: i32, iflag: Option<u32>) -> napi::Result<JsSweCalcOut> {
        self.calc_core(jd_et, ipl, iflag)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn fixstar2_ut(
        &self,
        star: String,
        jd_ut: f64,
        iflag: Option<u32>,
    ) -> napi::Result<JsSweFixstarOut> {
        self.fixstar2_ut_core(&star, jd_ut, iflag)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn fixstar2_mag(&self, star: String) -> napi::Result<f64> {
        self.fixstar2_mag_core(&star)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn houses_ex(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
        hsys: Option<String>,
        iflag: Option<u32>,
    ) -> napi::Result<JsSweHousesOut> {
        self.houses_ex_core(jd_ut, lat_deg, lon_deg, hsys, iflag)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn houses(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
        hsys: Option<String>,
    ) -> napi::Result<JsSweHousesOut> {
        self.houses_core(jd_ut, lat_deg, lon_deg, hsys)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn houses_armc(
        &self,
        armc_deg: f64,
        lat_deg: f64,
        eps_deg: f64,
        hsys: Option<String>,
        sun_decl_deg: Option<f64>,
    ) -> napi::Result<JsSweHousesOut> {
        self.houses_armc_core(armc_deg, lat_deg, eps_deg, hsys, sun_decl_deg)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn azalt(
        &self,
        jd_ut: f64,
        calc_flag: i32,
        geopos: Vec<f64>,
        atpress: f64,
        attemp: f64,
        xin: Vec<f64>,
    ) -> napi::Result<Vec<f64>> {
        self.azalt_core(jd_ut, calc_flag, &geopos, atpress, attemp, &xin)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans(
        &self,
        jd_ut: f64,
        ipl: i32,
        rsmi: i32,
        geopos: Vec<f64>,
        atpress: Option<f64>,
        attemp: Option<f64>,
        starname: Option<String>,
        epheflag: Option<u32>,
    ) -> napi::Result<Option<f64>> {
        self.rise_trans_core(
            jd_ut, ipl, rsmi, &geopos, atpress, attemp, starname, epheflag,
        )
        .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn pheno_ut(
        &self,
        jd_ut: f64,
        ipl: i32,
        iflag: Option<u32>,
    ) -> napi::Result<JsSwePhenoOut> {
        self.pheno_ut_core(jd_ut, ipl, iflag)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn pheno(&self, jd_et: f64, ipl: i32, iflag: Option<u32>) -> napi::Result<JsSwePhenoOut> {
        self.pheno_core(jd_et, ipl, iflag)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn nod_aps_ut(
        &self,
        jd_ut: f64,
        ipl: i32,
        iflag: Option<u32>,
        method: Option<u32>,
    ) -> napi::Result<JsNodApsResult> {
        self.nod_aps_ut_core(jd_ut, ipl, iflag, method)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn nod_aps(
        &self,
        jd_et: f64,
        ipl: i32,
        iflag: Option<u32>,
        method: Option<u32>,
    ) -> napi::Result<JsNodApsResult> {
        self.nod_aps_core(jd_et, ipl, iflag, method)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn lun_occult_when_glob(
        &self,
        jd_ut_start: f64,
        ipl: i32,
        starname: Option<String>,
        iflag: Option<u32>,
        ifltype_or_backward: Option<napi::bindgen_prelude::Either<i32, bool>>,
        backward: Option<bool>,
    ) -> napi::Result<JsSweOccultGlobOut> {
        let (ifltype, bw) = match ifltype_or_backward {
            Some(napi::bindgen_prelude::Either::A(t)) => (Some(t), backward),
            Some(napi::bindgen_prelude::Either::B(b)) => (None, Some(b)),
            None => (None, backward),
        };
        self.lun_occult_when_glob_core(jd_ut_start, ipl, starname, iflag, ifltype, bw)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn lun_occult_when_loc(
        &self,
        jd_ut_start: f64,
        ipl: i32,
        starname: Option<String>,
        iflag: Option<u32>,
        geopos: Vec<f64>,
        backward: Option<bool>,
    ) -> napi::Result<JsSweOccultLocOut> {
        self.lun_occult_when_loc_core(jd_ut_start, ipl, starname, iflag, &geopos, backward)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn house_pos(
        &self,
        armc_deg: f64,
        geolat_deg: f64,
        eps_deg: f64,
        xpin: Vec<f64>,
        hsys: Option<String>,
    ) -> napi::Result<f64> {
        self.house_pos_core(armc_deg, geolat_deg, eps_deg, &xpin, hsys)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn time_equ(&self, jd_et: f64) -> napi::Result<f64> {
        self.time_equ_core(jd_et).map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn azalt_rev(
        &self,
        jd_ut: f64,
        calc_flag: i32,
        geopos: Vec<f64>,
        xin: Vec<f64>,
    ) -> napi::Result<Vec<f64>> {
        self.azalt_rev_core(jd_ut, calc_flag, &geopos, &xin)
            .map_err(napi::Error::from_reason)
    }

    #[napi]
    pub fn julday(
        year: i32,
        month: u8,
        day: u8,
        hour: Option<f64>,
        gregflag: Option<i32>,
    ) -> napi::Result<f64> {
        SweAdapter::julday(year, month, day, hour.unwrap_or(0.0), gregflag.unwrap_or(1))
            .map_err(|e| napi::Error::from_reason(format_js_error(&e)))
    }

    #[napi]
    pub fn revjul(jd: f64, gregflag: Option<i32>) -> napi::Result<JsCalendarDate> {
        let d = SweAdapter::revjul(jd, gregflag.unwrap_or(1))
            .map_err(|e| napi::Error::from_reason(format_js_error(&e)))?;
        Ok(JsCalendarDate {
            year: d.year,
            month: d.month,
            day: d.day,
            hour_ut: d.hour,
        })
    }
}

// Export Swiss Ephemeris Migration Constants for Node.js / TypeScript
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_ECL_NUT: i32 = lxephemeris_core::SE_ECL_NUT;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SUN: i32 = lxephemeris_core::SE_SUN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_MOON: i32 = lxephemeris_core::SE_MOON;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_MERCURY: i32 = lxephemeris_core::SE_MERCURY;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_VENUS: i32 = lxephemeris_core::SE_VENUS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_MARS: i32 = lxephemeris_core::SE_MARS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_JUPITER: i32 = lxephemeris_core::SE_JUPITER;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SATURN: i32 = lxephemeris_core::SE_SATURN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_URANUS: i32 = lxephemeris_core::SE_URANUS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_NEPTUNE: i32 = lxephemeris_core::SE_NEPTUNE;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_PLUTO: i32 = lxephemeris_core::SE_PLUTO;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_MEAN_NODE: i32 = lxephemeris_core::SE_MEAN_NODE;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_TRUE_NODE: i32 = lxephemeris_core::SE_TRUE_NODE;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_MEAN_APOG: i32 = lxephemeris_core::SE_MEAN_APOG;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_OSCU_APOG: i32 = lxephemeris_core::SE_OSCU_APOG;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_EARTH: i32 = lxephemeris_core::SE_EARTH;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_CHIRON: i32 = lxephemeris_core::SE_CHIRON;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_PHOLUS: i32 = lxephemeris_core::SE_PHOLUS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_CERES: i32 = lxephemeris_core::SE_CERES;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_PALLAS: i32 = lxephemeris_core::SE_PALLAS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_JUNO: i32 = lxephemeris_core::SE_JUNO;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_VESTA: i32 = lxephemeris_core::SE_VESTA;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_INTP_APOG: i32 = lxephemeris_core::SE_INTP_APOG;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_INTP_PERG: i32 = lxephemeris_core::SE_INTP_PERG;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_FICT_OFFSET: i32 = lxephemeris_core::SE_FICT_OFFSET;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_CUPIDO: i32 = lxephemeris_core::SE_CUPIDO;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_HADES: i32 = lxephemeris_core::SE_HADES;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_ZEUS: i32 = lxephemeris_core::SE_ZEUS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_KRONOS: i32 = lxephemeris_core::SE_KRONOS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_APOLLON: i32 = lxephemeris_core::SE_APOLLON;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_ADMETOS: i32 = lxephemeris_core::SE_ADMETOS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_VULKANUS: i32 = lxephemeris_core::SE_VULKANUS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_POSEIDON: i32 = lxephemeris_core::SE_POSEIDON;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_ISIS: i32 = lxephemeris_core::SE_ISIS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_AST_OFFSET: i32 = lxephemeris_core::SE_AST_OFFSET;

#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_SWIEPH: u32 = lxephemeris_core::SEFLG_SWIEPH;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_MOSEPH: u32 = lxephemeris_core::SEFLG_MOSEPH;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_HELCTR: u32 = lxephemeris_core::SEFLG_HELCTR;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_TRUEPOS: u32 = lxephemeris_core::SEFLG_TRUEPOS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_J2000: u32 = lxephemeris_core::SEFLG_J2000;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_NONUT: u32 = lxephemeris_core::SEFLG_NONUT;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_SPEED: u32 = lxephemeris_core::SEFLG_SPEED;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_NOGDEFL: u32 = lxephemeris_core::SEFLG_NOGDEFL;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_NOABERR: u32 = lxephemeris_core::SEFLG_NOABERR;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_BARYCTR: u32 = lxephemeris_core::SEFLG_BARYCTR;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_TOPOCTR: u32 = lxephemeris_core::SEFLG_TOPOCTR;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_EQUATORIAL: u32 = lxephemeris_core::SEFLG_EQUATORIAL;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_XYZ: u32 = lxephemeris_core::SEFLG_XYZ;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_RADIANS: u32 = lxephemeris_core::SEFLG_RADIANS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SEFLG_SIDEREAL: u32 = lxephemeris_core::SEFLG_SIDEREAL;

#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_FAGAN_BRADLEY: u32 = lxephemeris_core::SE_SIDM_FAGAN_BRADLEY;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_LAHIRI: u32 = lxephemeris_core::SE_SIDM_LAHIRI;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_RAMAN: u32 = lxephemeris_core::SE_SIDM_RAMAN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_USHASHASHI: u32 = lxephemeris_core::SE_SIDM_USHASHASHI;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_KRISHNAMURTI: u32 = lxephemeris_core::SE_SIDM_KRISHNAMURTI;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_YUKTESHWAR: u32 = lxephemeris_core::SE_SIDM_YUKTESHWAR;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_JN_BHASIN: u32 = lxephemeris_core::SE_SIDM_JN_BHASIN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_SURYASIDDHANTA: u32 = lxephemeris_core::SE_SIDM_SURYASIDDHANTA;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_SS_CITRA: u32 = lxephemeris_core::SE_SIDM_SS_CITRA;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_TRUE_CITRA: u32 = lxephemeris_core::SE_SIDM_TRUE_CITRA;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_TRUE_REVATI: u32 = lxephemeris_core::SE_SIDM_TRUE_REVATI;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_TRUE_PUSHYA: u32 = lxephemeris_core::SE_SIDM_TRUE_PUSHYA;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_USER: u32 = lxephemeris_core::SE_SIDM_USER;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_DELUCE: u32 = lxephemeris_core::SE_SIDM_DELUCE;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_ARYABHATA: u32 = lxephemeris_core::SE_SIDM_ARYABHATA;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_ARYABHATA_MSUN: u32 = lxephemeris_core::SE_SIDM_ARYABHATA_MSUN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_SS_REVATI: u32 = lxephemeris_core::SE_SIDM_SS_REVATI;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_SURYASIDDHANTA_MSUN: u32 = lxephemeris_core::SE_SIDM_SURYASIDDHANTA_MSUN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_HIPPARCHOS: u32 = lxephemeris_core::SE_SIDM_HIPPARCHOS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_SASSANIAN: u32 = lxephemeris_core::SE_SIDM_SASSANIAN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_J1900: u32 = lxephemeris_core::SE_SIDM_J1900;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_B1950: u32 = lxephemeris_core::SE_SIDM_B1950;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_J2000: u32 = lxephemeris_core::SE_SIDM_J2000;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_BABYL_KUGLER1: u32 = lxephemeris_core::SE_SIDM_BABYL_KUGLER1;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_BABYL_KUGLER2: u32 = lxephemeris_core::SE_SIDM_BABYL_KUGLER2;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_BABYL_KUGLER3: u32 = lxephemeris_core::SE_SIDM_BABYL_KUGLER3;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_BABYL_HUBER: u32 = lxephemeris_core::SE_SIDM_BABYL_HUBER;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_BABYL_ETPSC: u32 = lxephemeris_core::SE_SIDM_BABYL_ETPSC;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_ALDEBARAN_15TAU: u32 = lxephemeris_core::SE_SIDM_ALDEBARAN_15TAU;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_GALCENT_RGILBRAND: u32 = lxephemeris_core::SE_SIDM_GALCENT_RGILBRAND;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_GALCENT_0SAG: u32 = lxephemeris_core::SE_SIDM_GALCENT_0SAG;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_GALALIGN_MARDYKS: u32 = lxephemeris_core::SE_SIDM_GALALIGN_MARDYKS;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_TRUE_MULA: u32 = lxephemeris_core::SE_SIDM_TRUE_MULA;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_GALCENT_MULA_WILHELM: u32 = lxephemeris_core::SE_SIDM_GALCENT_MULA_WILHELM;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_ARYABHATA_522: u32 = lxephemeris_core::SE_SIDM_ARYABHATA_522;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_BABYL_BRITTON: u32 = lxephemeris_core::SE_SIDM_BABYL_BRITTON;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_TRUE_SHERATAN: u32 = lxephemeris_core::SE_SIDM_TRUE_SHERATAN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_TRUE_ALDEBARAN: u32 = lxephemeris_core::SE_SIDM_TRUE_ALDEBARAN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_DJWHAL_KHUL: u32 = lxephemeris_core::SE_SIDM_DJWHAL_KHUL;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_VALENS_MOON: u32 = lxephemeris_core::SE_SIDM_VALENS_MOON;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_GALEQ_IAU1958: u32 = lxephemeris_core::SE_SIDM_GALEQ_IAU1958;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_GALEQ_TRUE: u32 = lxephemeris_core::SE_SIDM_GALEQ_TRUE;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_GALEQU_FIORENZA: u32 = lxephemeris_core::SE_SIDM_GALEQU_FIORENZA;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_GALEQ_MULA: u32 = lxephemeris_core::SE_SIDM_GALEQ_MULA;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_LAHIRI_1940: u32 = lxephemeris_core::SE_SIDM_LAHIRI_1940;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_LAHIRI_VP285: u32 = lxephemeris_core::SE_SIDM_LAHIRI_VP285;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_KP_SENTHILATHIBAN: u32 = lxephemeris_core::SE_SIDM_KP_SENTHILATHIBAN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDM_LAHIRI_ICRC: u32 = lxephemeris_core::SE_SIDM_LAHIRI_ICRC;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDBIT_ECL_T0: u32 = lxephemeris_core::SE_SIDBIT_ECL_T0;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDBIT_SSY_PLANE: u32 = lxephemeris_core::SE_SIDBIT_SSY_PLANE;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDBIT_USER_UT: u32 = lxephemeris_core::SE_SIDBIT_USER_UT;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDBIT_ECL_DATE: u32 = lxephemeris_core::SE_SIDBIT_ECL_DATE;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDBIT_NO_PREC_OFFSET: u32 = lxephemeris_core::SE_SIDBIT_NO_PREC_OFFSET;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_SIDBIT_PREC_ORIG: u32 = lxephemeris_core::SE_SIDBIT_PREC_ORIG;

#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_CALC_RISE: i32 = lxephemeris_core::SE_CALC_RISE;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_CALC_SET: i32 = lxephemeris_core::SE_CALC_SET;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_CALC_MTRANSIT: i32 = lxephemeris_core::SE_CALC_MTRANSIT;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_CALC_ITRANSIT: i32 = lxephemeris_core::SE_CALC_ITRANSIT;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_BIT_DISC_CENTER: i32 = lxephemeris_core::SE_BIT_DISC_CENTER;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_BIT_DISC_BOTTOM: i32 = lxephemeris_core::SE_BIT_DISC_BOTTOM;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_BIT_GEOCTR_NO_ECL_LAT: i32 = lxephemeris_core::SE_BIT_GEOCTR_NO_ECL_LAT;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_BIT_NO_REFRACTION: i32 = lxephemeris_core::SE_BIT_NO_REFRACTION;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_BIT_CIVIL_TWILIGHT: i32 = lxephemeris_core::SE_BIT_CIVIL_TWILIGHT;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_BIT_NAUTIC_TWILIGHT: i32 = lxephemeris_core::SE_BIT_NAUTIC_TWILIGHT;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_BIT_ASTRO_TWILIGHT: i32 = lxephemeris_core::SE_BIT_ASTRO_TWILIGHT;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_BIT_FIXED_DISC_SIZE: i32 = lxephemeris_core::SE_BIT_FIXED_DISC_SIZE;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_BIT_HINDU_RISING: i32 = lxephemeris_core::SE_BIT_HINDU_RISING;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_ECL2HOR: i32 = lxephemeris_core::SE_ECL2HOR;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_EQU2HOR: i32 = lxephemeris_core::SE_EQU2HOR;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_HOR2ECL: i32 = lxephemeris_core::SE_HOR2ECL;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_HOR2EQU: i32 = lxephemeris_core::SE_HOR2EQU;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_NODBIT_MEAN: u32 = lxephemeris_core::SE_NODBIT_MEAN;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_NODBIT_OSCU: u32 = lxephemeris_core::SE_NODBIT_OSCU;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_NODBIT_OSCU_BAR: u32 = lxephemeris_core::SE_NODBIT_OSCU_BAR;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_NODBIT_FOPOINT: u32 = lxephemeris_core::SE_NODBIT_FOPOINT;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_JUL_CAL: i32 = lxephemeris_core::SE_JUL_CAL;
#[cfg_attr(feature = "napi-addon", napi)]
pub const SE_GREG_CAL: i32 = lxephemeris_core::SE_GREG_CAL;
}
#[cfg(feature = "swe-compat")]
pub use swe_node::*;

/// Converts a calendar date `(year, month, day, hour_ut)` into a Julian Date (`f64`).
pub fn julday_core(
    year: i32,
    month: u8,
    day: u8,
    hour_ut: Option<f64>,
    calendar: Option<String>,
) -> Result<f64, String> {
    let cal = parse_calendar_kind(calendar.as_deref().unwrap_or("gregorian"))
        .map_err(|e| format_js_error(&e))?;
    core_julday(year, month, day, hour_ut.unwrap_or(0.0), cal).map_err(|e| format_js_error(&e))
}

/// Converts a Julian Date (`jd`) into a calendar date `JsCalendarDate`.
pub fn revjul_core(jd: f64, calendar: Option<String>) -> Result<JsCalendarDate, String> {
    let cal = parse_calendar_kind(calendar.as_deref().unwrap_or("gregorian"))
        .map_err(|e| format_js_error(&e))?;
    let d = core_revjul(jd, cal).map_err(|e| format_js_error(&e))?;
    Ok(JsCalendarDate {
        year: d.year,
        month: d.month,
        day: d.day,
        hour_ut: d.hour,
    })
}

/// Evaluates house cusps and chart angles directly from `armc_deg`, `eps_true_deg`, and `latitude_deg`.
#[allow(clippy::too_many_arguments)]
pub fn houses_armc_core(
    armc_deg: f64,
    eps_true_deg: f64,
    latitude_deg: f64,
    system: Option<String>,
    polar_policy: Option<String>,
    ayanamsa_deg: Option<f64>,
    sun_decl_deg: Option<f64>,
    include_speed: Option<bool>,
) -> Result<JsHouseResult, String> {
    let sys_str = system.as_deref().unwrap_or("P");
    let ch = sys_str.chars().next().unwrap_or(' ');
    let sys = HouseSystem::from_code(ch).map_err(|e| format_js_error(&e))?;
    let pol = parse_polar_policy(polar_policy.as_deref().unwrap_or("strict"))
        .map_err(|e| format_js_error(&e))?;
    let opts = HouseOptions {
        system: sys,
        polar_policy: pol,
        zodiac: Zodiac::Tropical,
        nutation: NutationMode::True,
        include_speed: include_speed.unwrap_or(true),
    };
    let res = houses_armc_with_sidereal(
        armc_deg,
        latitude_deg,
        eps_true_deg,
        sun_decl_deg.unwrap_or(0.0),
        ayanamsa_deg.unwrap_or(0.0),
        opts,
    )
    .map_err(|e| format_js_error(&e))?;
    Ok(JsHouseResult::from(res))
}

/// Determines the house number (`1..=12`) for `lon_deg` given 12 cusps (`R4.5`).
pub fn house_of_lon_core(
    lon_deg: f64,
    cusps_deg: &[f64],
    eps_deg: Option<f64>,
) -> Result<u8, String> {
    let arr: [f64; 12] = cusps_deg.try_into().map_err(|_| {
        format_js_error(&EphemerisError::InvalidInput {
            message: format!("Expected 12 house cusps, got {}", cusps_deg.len()),
        })
    })?;
    Ok(core_house_of_lon(lon_deg, &arr, eps_deg.unwrap_or(1e-9)))
}

/// Computes the 36 Gauquelin sectors (`Vec<f64>` of length 36) from `armc_deg`, `eps_true_deg`, `latitude_deg`.
pub fn gauquelin_sectors_core(
    armc_deg: f64,
    eps_true_deg: f64,
    latitude_deg: f64,
    polar_policy: Option<String>,
) -> Result<Vec<f64>, String> {
    let pol = parse_polar_policy(polar_policy.as_deref().unwrap_or("strict"))
        .map_err(|e| format_js_error(&e))?;
    let res = core_gauquelin_sectors(armc_deg, latitude_deg, eps_true_deg, pol)
        .map_err(|e| format_js_error(&e))?;
    Ok(res.sectors_deg.to_vec())
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn julday(
    year: i32,
    month: u8,
    day: u8,
    hour_ut: Option<f64>,
    calendar: Option<String>,
) -> napi::Result<f64> {
    julday_core(year, month, day, hour_ut, calendar).map_err(napi::Error::from_reason)
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn revjul(jd: f64, calendar: Option<String>) -> napi::Result<JsCalendarDate> {
    revjul_core(jd, calendar).map_err(napi::Error::from_reason)
}

#[cfg(feature = "napi-addon")]
#[napi]
#[allow(clippy::too_many_arguments)]
pub fn houses_armc(
    armc_deg: f64,
    eps_true_deg: f64,
    latitude_deg: f64,
    system: Option<String>,
    polar_policy: Option<String>,
    ayanamsa_deg: Option<f64>,
    sun_decl_deg: Option<f64>,
    include_speed: Option<bool>,
) -> napi::Result<JsHouseResult> {
    houses_armc_core(
        armc_deg,
        eps_true_deg,
        latitude_deg,
        system,
        polar_policy,
        ayanamsa_deg,
        sun_decl_deg,
        include_speed,
    )
    .map_err(napi::Error::from_reason)
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn house_of_lon(lon_deg: f64, cusps_deg: Vec<f64>, eps_deg: Option<f64>) -> napi::Result<u8> {
    house_of_lon_core(lon_deg, &cusps_deg, eps_deg).map_err(napi::Error::from_reason)
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn star_info(name_or_bayer: String) -> napi::Result<JsStarMetadata> {
    lxephemeris_core::star_info(&name_or_bayer)
        .map(JsStarMetadata::from)
        .map_err(|e| napi::Error::from_reason(format_js_error(&e)))
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn list_stars() -> Vec<JsStarMetadata> {
    lxephemeris_core::list_stars()
        .into_iter()
        .map(JsStarMetadata::from)
        .collect()
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn gauquelin_sectors(
    armc_deg: f64,
    eps_true_deg: f64,
    latitude_deg: f64,
    polar_policy: Option<String>,
) -> napi::Result<Vec<f64>> {
    gauquelin_sectors_core(armc_deg, eps_true_deg, latitude_deg, polar_policy)
        .map_err(napi::Error::from_reason)
}

/// Converts civil wall-clock time `(year, month, day, hour, minute, second)` with `tz_offset_hours` into a UTC Julian Date (`f64`).
#[allow(clippy::too_many_arguments)]
pub fn utc_from_local_core(
    year: i32,
    month: u32,
    day: u32,
    hour: Option<u32>,
    minute: Option<u32>,
    second: Option<f64>,
    tz_offset_hours: Option<f64>,
    gregorian: Option<bool>,
) -> Result<f64, String> {
    let m = u8::try_from(month).map_err(|_| {
        format_js_error(&EphemerisError::InvalidInput {
            message: format!("Month out of range: {month}"),
        })
    })?;
    let d = u8::try_from(day).map_err(|_| {
        format_js_error(&EphemerisError::InvalidInput {
            message: format!("Day out of range: {day}"),
        })
    })?;
    let h = u8::try_from(hour.unwrap_or(0)).map_err(|_| {
        format_js_error(&EphemerisError::InvalidInput {
            message: format!("Hour out of range: {:?}", hour),
        })
    })?;
    let min = u8::try_from(minute.unwrap_or(0)).map_err(|_| {
        format_js_error(&EphemerisError::InvalidInput {
            message: format!("Minute out of range: {:?}", minute),
        })
    })?;
    let sec = second.unwrap_or(0.0);
    let tz = tz_offset_hours.unwrap_or(0.0);
    if h > 23 || min > 59 || !sec.is_finite() || !(0.0..=61.0).contains(&sec) || !tz.is_finite() {
        return Err(format_js_error(&EphemerisError::InvalidInput {
            message: format!(
                "Invalid local time components: hour={h}, minute={min}, second={sec}, tz_offset_hours={tz}"
            ),
        }));
    }
    let cal = if gregorian.unwrap_or(true) {
        CalendarKind::Gregorian
    } else {
        CalendarKind::Julian
    };
    let hour_ut = f64::from(h) + f64::from(min) / 60.0 + sec / 3600.0 - tz;
    core_julday(year, m, d, hour_ut, cal).map_err(|e| format_js_error(&e))
}

/// Maps an ecliptic longitude (`lon_deg`) into Zodiac Sign, Decan, Nakshatra (27 × 4 Pada), and Human Design 64 Gate / 6 Line / Color / Tone / Base.
pub fn map_zodiac_placement_core(lon_deg: f64) -> Result<serde_json::Value, String> {
    to_js_value(&core_map_zodiac_placement(lon_deg))
}

#[cfg(feature = "napi-addon")]
#[napi]
#[allow(clippy::too_many_arguments)]
pub fn utc_from_local(
    year: i32,
    month: u32,
    day: u32,
    hour: Option<u32>,
    minute: Option<u32>,
    second: Option<f64>,
    tz_offset_hours: Option<f64>,
    gregorian: Option<bool>,
) -> napi::Result<f64> {
    utc_from_local_core(
        year,
        month,
        day,
        hour,
        minute,
        second,
        tz_offset_hours,
        gregorian,
    )
    .map_err(napi::Error::from_reason)
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn map_zodiac_placement(lon_deg: f64) -> napi::Result<serde_json::Value> {
    map_zodiac_placement_core(lon_deg).map_err(napi::Error::from_reason)
}

pub fn verify_license_key_core(token: &str) -> serde_json::Value {
    let info = lxephemeris_core::verify_license_key(token);
    serde_json::json!({
        "status": info.status.as_str(),
        "tier": info.tier,
        "isCommercialValid": info.is_commercial_valid,
        "licenseId": info.license_id,
        "holder": info.holder,
        "keyId": info.key_id,
        "issuedAt": info.issued_at,
        "updatesUntil": info.updates_until,
        "engineReleaseDate": info.engine_release_date,
        "message": info.message,
    })
}

pub fn register_commercial_license_core(token: &str) -> serde_json::Value {
    let info = lxephemeris_core::register_commercial_license(token);
    serde_json::json!({
        "status": info.status.as_str(),
        "tier": info.tier,
        "isCommercialValid": info.is_commercial_valid,
        "licenseId": info.license_id,
        "holder": info.holder,
        "keyId": info.key_id,
        "issuedAt": info.issued_at,
        "updatesUntil": info.updates_until,
        "engineReleaseDate": info.engine_release_date,
        "message": info.message,
    })
}

pub fn license_info_global_core() -> serde_json::Value {
    let info = lxephemeris_core::resolve_license_info(None);
    serde_json::json!({
        "status": info.status.as_str(),
        "tier": info.tier,
        "isCommercialValid": info.is_commercial_valid,
        "licenseId": info.license_id,
        "holder": info.holder,
        "keyId": info.key_id,
        "issuedAt": info.issued_at,
        "updatesUntil": info.updates_until,
        "engineReleaseDate": info.engine_release_date,
        "message": info.message,
    })
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn verify_license_key(token: String) -> serde_json::Value {
    verify_license_key_core(&token)
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn register_commercial_license(token: String) -> serde_json::Value {
    register_commercial_license_core(&token)
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn clear_commercial_license() {
    lxephemeris_core::clear_commercial_license();
}

#[cfg(feature = "napi-addon")]
#[napi]
pub fn license_info() -> serde_json::Value {
    license_info_global_core()
}

#[cfg(test)]
#[allow(clippy::indexing_slicing)]
mod tests {
    use super::*;
    use lxephemeris_core::JD_J2000;
    #[cfg(feature = "swe-compat")]
    use lxephemeris_core::{SEFLG_SPEED, SE_SUN};
    use std::path::Path;

    fn workspace_data_dir() -> String {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(|p| p.parent())
            .unwrap_or_else(|| Path::new("."))
            .join("data")
            .display()
            .to_string()
    }

    #[test]
    fn test_nodejs_vertical_slice_single_batch_and_errors() {
        let data_dir = workspace_data_dir();
        if !Path::new(&data_dir).join("de440s.bsp").exists() {
            return;
        }
        let engine = JsEngine::from_data_dir_core(
            data_dir,
            Some(2),
            Some("embedded".to_owned()),
            Some(false),
        )
        .expect("load JsEngine");

        let jd_j2000 = JsSplitJd {
            jd1: JD_J2000,
            jd2: Some(0.0),
        };

        // 1. Single calculation (Sun apparent ecliptic longitude at J2000.0 = 280.36816528 deg)
        let sun = engine
            .calc_core("Sun", jd_j2000, Some("TT".to_owned()), None)
            .expect("calc Sun");
        assert!((sun.lon_deg - 280.368_165_28).abs() < 1e-6);
        assert!(sun.lon_speed_deg_day > 1.01 && sun.lon_speed_deg_day < 1.03);
        assert!(sun.dist_km > 147_000_000.0 && sun.dist_km < 147_200_000.0);
        assert!(sun.warnings.is_empty());

        // 2. Named minor planet & fixed star
        let chiron = engine
            .calc_core("Chiron", jd_j2000, Some("TT".to_owned()), None)
            .expect("calc Chiron");
        assert!(chiron.lon_deg >= 0.0 && chiron.lon_deg < 360.0);

        let regulus = engine
            .calc_star_by_name_core("Regulus", jd_j2000, Some("TT".to_owned()), None)
            .expect("calc Regulus");
        assert!((regulus.lon_deg - 149.829_040_6).abs() < 1e-6);

        // 3. Ayanamsa & houses & eclipses & calendar
        let lahiri = engine
            .ayanamsa_core("lahiri", jd_j2000, Some("TT".to_owned()), None)
            .expect("ayanamsa lahiri");
        assert!((lahiri.mean_deg - 23.857_092_3).abs() < 1e-4);

        let obs = JsObserver {
            latitude_deg: 51.5074,
            longitude_deg: -0.1278,
            altitude_m: Some(0.0),
        };
        let h_p = engine
            .houses_core(
                jd_j2000,
                obs,
                Some("P".to_owned()),
                Some("TT".to_owned()),
                None,
                None,
                None,
                None,
                None,
            )
            .expect("houses Placidus");
        assert_eq!(h_p.cusps_deg.len(), 12);

        let sol = engine
            .next_solar_eclipse_core(
                JsSplitJd {
                    jd1: 2_460_400.0,
                    jd2: Some(0.0),
                },
                Some("TT".to_owned()),
                Some(false),
                None,
                None,
            )
            .expect("next solar eclipse");
        assert_eq!(sol.kind, "Total");
        assert_eq!(sol.saros_series, 139);
        assert_eq!(sol.saros_member, 30);

        let jd = julday_core(2000, 1, 1, Some(12.0), None).expect("julday");
        assert!((jd - JD_J2000).abs() < 1e-12);
        let d = revjul_core(jd, None).expect("revjul");
        assert_eq!((d.year, d.month, d.day), (2000, 1, 1));

        // 4. Flat Float64Array batch calculation
        let bodies = vec!["Sun".to_owned(), "Moon".to_owned(), "mpc:2060".to_owned()];
        let epochs = vec![
            jd_j2000,
            JsSplitJd {
                jd1: JD_J2000,
                jd2: Some(1.0),
            },
        ];
        let flat = engine
            .calc_batch_flat_f64_core(&bodies, &epochs, Some("TT".to_owned()), None)
            .expect("calc_batch_flat_f64");
        assert_eq!(flat.len(), 2 * 3 * 6);
        assert!((flat[0] - sun.lon_deg).abs() < 1e-12);

        // 5. Per-item batch with out-of-coverage epoch
        let mixed_epochs = vec![
            jd_j2000,
            JsSplitJd {
                jd1: 2_200_000.0,
                jd2: Some(0.0),
            },
        ];
        let per_item = engine
            .calc_batch_per_item_core(&bodies, &mixed_epochs, Some("TT".to_owned()), None)
            .expect("calc_batch_per_item");
        assert!(per_item[0][0].ok);
        assert!(!per_item[1][0].ok);
        assert_eq!(per_item[1][0].error_code.as_deref(), Some("OutOfCoverage"));

        // 6. Structured error codes
        let err_body = engine
            .calc_core("NotABody", jd_j2000, None, None)
            .expect_err("must fail on unknown body");
        assert!(err_body.contains("[UnsupportedBody]"));

        // 7. SweAdapter check
        #[cfg(feature = "swe-compat")]
        {
            let swe = JsSweAdapter::from_engine_core(&engine);
            let out = swe
                .calc_core(JD_J2000, SE_SUN, Some(SEFLG_SPEED))
                .expect("swe.calc");
            assert!((out.xx[0] - sun.lon_deg).abs() < 1e-12);
        }
    }

    #[test]
    fn test_wave_b_and_c_features() {
        let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let repo_root = manifest_dir
            .parent()
            .and_then(|p| p.parent())
            .expect("repo root");
        let data_dir = repo_root.join("data");
        if !data_dir.join("de440s.bsp").exists() {
            return;
        }

        let engine = JsEngine::from_data_dir_core(
            data_dir.display().to_string(),
            Some(2),
            None,
            None,
        )
        .expect("from_data_dir_core");
        assert!(!engine.is_moshier_mode_core());

        // 1. Moshier mode engine
        let mosh_eng = JsEngine::new_moshier_core(None, None);
        assert!(mosh_eng.is_moshier_mode_core());
        let jd_j2000 = JsSplitJd {
            jd1: JD_J2000,
            jd2: Some(0.0),
        };
        let mosh_sun = mosh_eng
            .calc_core("Sun", jd_j2000, Some("TT".to_owned()), None)
            .expect("calc Sun in Moshier");
        assert!((mosh_sun.lon_deg - 280.368_165).abs() < 0.05);

        // calc_moshier on SPK engine
        let mosh_via_spk = engine
            .calc_moshier_core("Sun", jd_j2000, Some("TT".to_owned()), None)
            .expect("calc_moshier");
        assert!((mosh_via_spk.lon_deg - mosh_sun.lon_deg).abs() < 1e-10);

        // 2. Phenomena
        let pheno = engine
            .phenomena_core("Venus", jd_j2000, Some("TT".to_owned()), None)
            .expect("phenomena Venus");
        assert!((0.0..=180.0).contains(&pheno.phase_angle_deg));
        assert!((0.0..=1.0).contains(&pheno.phase_illuminated));
        assert_eq!(pheno.attr20.len(), 20);

        // 3. Nodes and Apsides
        let nod = engine
            .nod_aps_core("Mars", jd_j2000, Some("TT".to_owned()), None, None)
            .expect("nod_aps Mars");
        assert_eq!(nod.ascending_node.len(), 6);
        assert_eq!(nod.descending_node.len(), 6);
        assert_eq!(nod.periapsis.len(), 6);
        assert_eq!(nod.apoapsis.len(), 6);

        // 4. Lunar Occultation of Aldebaran in 2015
        let jd_2015 = 2_457_023.5;
        let jd_2015_split = JsSplitJd {
            jd1: jd_2015,
            jd2: Some(0.0),
        };
        let occ = engine
            .next_lunar_occultation_core(
                "Aldebaran",
                jd_2015_split,
                Some("UT1".to_owned()),
                Some(false),
                Some(60.0),
            )
            .expect("next_lunar_occultation");
        assert!(occ.jd_ut1_max > jd_2015);
        assert_eq!(occ.tret.len(), 10);

        // 5. House position
        let hpos = engine
            .house_pos_core(100.0, 51.5, 23.44, "P", 120.0, 0.0)
            .expect("house_pos");
        assert!((1.0..13.0).contains(&hpos));

        // 6. Time equation
        let teq = engine
            .time_equation_core(jd_j2000, Some("TT".to_owned()))
            .expect("time_equation");
        assert!(teq.abs() < 0.05);

        // 7. azalt_rev
        let obs = JsObserver {
            latitude_deg: 51.5074,
            longitude_deg: -0.1278,
            altitude_m: Some(0.0),
        };
        let rev = engine
            .azalt_rev_core(
                obs,
                jd_j2000,
                Some("UT1".to_owned()),
                180.0,
                45.0,
                Some("ecliptic".to_owned()),
            )
            .expect("azalt_rev");
        assert_eq!(rev.len(), 2);

        // 8. SweAdapter Wave B & C methods
        #[cfg(feature = "swe-compat")]
        {
            let swe = JsSweAdapter::from_engine_core(&engine);

            // Custom sidereal mode
            swe.set_sid_mode_core(SE_SIDM_USER, Some(JD_J2000), Some(24.0))
                .expect("set_sid_mode_core custom");

            let p_ut = swe
                .pheno_ut_core(JD_J2000, SE_VENUS, None)
                .expect("swe pheno_ut");
            assert_eq!(p_ut.attr.len(), 20);

            let p_et = swe.pheno_core(JD_J2000, SE_VENUS, None).expect("swe pheno");
            assert_eq!(p_et.attr.len(), 20);

            let na_ut = swe
                .nod_aps_ut_core(JD_J2000, SE_MARS, None, None)
                .expect("swe nod_aps_ut");
            assert_eq!(na_ut.ascending_node.len(), 6);

            let na_et = swe
                .nod_aps_core(JD_J2000, SE_MARS, None, None)
                .expect("swe nod_aps");
            assert_eq!(na_et.ascending_node.len(), 6);

            let swe_occ_glob = swe
                .lun_occult_when_glob_core(jd_2015, 0, Some("Aldebaran".to_owned()), None, None, None)
                .expect("swe lun_occult_when_glob");
            assert_eq!(swe_occ_glob.tret.len(), 10);

            let geopos = vec![0.0, 51.4769, 0.0];
            let swe_occ_loc = swe
                .lun_occult_when_loc_core(
                    jd_2015,
                    0,
                    Some("Aldebaran".to_owned()),
                    None,
                    &geopos,
                    None,
                )
                .expect("swe lun_occult_when_loc");
            assert_eq!(swe_occ_loc.tret.len(), 10);
            assert_eq!(swe_occ_loc.attr.len(), 20);

            let swe_hpos = swe
                .house_pos_core(100.0, 51.5, 23.44, &[120.0, 0.0], Some("P".to_owned()))
                .expect("swe house_pos");
            assert!((swe_hpos - hpos).abs() < 1e-6);

            let swe_teq = swe.time_equ_core(JD_J2000).expect("swe time_equ");
            assert!((swe_teq - teq).abs() < 1e-6);

            let xin = vec![0.0, 45.0];
            let swe_rev = swe
                .azalt_rev_core(JD_J2000, SE_HOR2ECL, &geopos, &xin)
                .expect("swe azalt_rev");
            assert_eq!(swe_rev.len(), 2);

            // SweAdapter new_moshier
            let mosh_swe = JsSweAdapter::new_moshier_core();
            let mosh_calc = mosh_swe
                .calc_ut_core(JD_J2000, SE_SUN, None)
                .expect("mosh_swe calc_ut");
            assert_eq!(mosh_calc.xx.len(), 6);

            // Constant checks
            assert_eq!(SE_SIDM_USER, 255);
            assert_eq!(SE_NODBIT_MEAN, 1);
            assert_eq!(SE_HOR2ECL, 0);
            assert_eq!(SE_HOR2EQU, 1);
        }
    }

    #[test]
    fn test_wave_f_astronomical_events_and_parans() {
        let data_dir = workspace_data_dir();
        if !Path::new(&data_dir).join("de440s.bsp").exists() {
            return;
        }
        let engine = JsEngine::from_data_dir_core(
            data_dir,
            Some(2),
            Some("embedded".to_owned()),
            Some(false),
        )
        .expect("load JsEngine");

        // 1. search_body_longitude (Vernal Equinox 2024 ~ 2024-03-20 03:06 UT)
        let jd_2024_jan_1 = julday_core(2024, 1, 1, Some(0.0), None).expect("jd 2024-01-01");
        let vernal = engine
            .search_body_longitude_core("Sun", jd_2024_jan_1, 0.0, Some(true), Some(120.0), None)
            .expect("search_body_longitude Sun 0 deg");
        assert_eq!(vernal.body, "Sun");
        assert!((vernal.actual_lon_deg - 0.0).abs() < 1e-6 || (vernal.actual_lon_deg - 360.0).abs() < 1e-6);
        assert!(!vernal.is_retrograde);

        // 2. body_longitude_crossings_in_range (Mercury triple pass in 2024)
        let jd_2024_end = jd_2024_jan_1 + 366.0;
        let crossings = engine
            .body_longitude_crossings_in_range_core(
                "Mercury",
                jd_2024_jan_1,
                jd_2024_end,
                20.0,
                None,
            )
            .expect("body_longitude_crossings_in_range Mercury");
        assert!(crossings.len() >= 3);
        assert!(crossings.iter().any(|c| c.is_retrograde));

        // 3. solar_terms_in_range (24 Chinese Solar Terms in 2024)
        let terms = engine
            .solar_terms_in_range_core(jd_2024_jan_1, jd_2024_end)
            .expect("solar_terms_in_range 2024");
        assert_eq!(terms.len(), 24);
        assert!(terms.iter().any(|t| t.pinyin_name == "Chunfen" && t.is_major_term));
        assert!(terms.iter().any(|t| t.pinyin_name == "Dongzhi" && t.is_major_term));

        // 4. seasons_in_year
        let seasons = engine.seasons_in_year_core(2024).expect("seasons_in_year 2024");
        assert_eq!(seasons.year, 2024);
        assert!((seasons.vernal_equinox_jd_ut1 - vernal.jd_ut1).abs() < 1e-6);
        assert!(seasons.summer_solstice_jd_ut1 > seasons.vernal_equinox_jd_ut1);
        assert!(seasons.autumnal_equinox_jd_ut1 > seasons.summer_solstice_jd_ut1);
        assert!(seasons.winter_solstice_jd_ut1 > seasons.autumnal_equinox_jd_ut1);

        // 5. human_design_epoch (-88 deg solar arc)
        let hd = engine
            .human_design_epoch_core(JD_J2000)
            .expect("human_design_epoch");
        assert_eq!(hd.body, "Sun");
        assert!((JD_J2000 - hd.jd_ut1) > 85.0 && (JD_J2000 - hd.jd_ut1) < 92.0);

        // 6. search_lunar_phase
        let full_moon = engine
            .search_lunar_phase_core(jd_2024_jan_1, 180.0, Some(true))
            .expect("search_lunar_phase FullMoon");
        assert_eq!(full_moon.kind.as_deref(), Some("FullMoon"));
        assert!((full_moon.target_elongation_deg - 180.0).abs() < 1e-9);
        assert!(full_moon.moon_dist_km > 350_000.0 && full_moon.moon_dist_km < 410_000.0);

        // 7. lunar_phases_in_range
        let phases = engine
            .lunar_phases_in_range_core(jd_2024_jan_1, jd_2024_jan_1 + 32.0)
            .expect("lunar_phases_in_range");
        assert!(phases.len() >= 4);

        // 8. lunar_phase_summary
        let summary = engine
            .lunar_phase_summary_core(JD_J2000)
            .expect("lunar_phase_summary");
        assert!((0.0..360.0).contains(&summary.elongation_deg));
        assert!((0..=7).contains(&summary.phase_8_index));
        assert!((1..=28).contains(&summary.phase_28_index));
        assert!(summary.prev_new_moon_jd_ut1 <= JD_J2000);
        assert!(summary.next_new_moon_jd_ut1 >= JD_J2000);

        // 9. search_stations_in_range
        let stations = engine
            .search_stations_in_range_core("Mercury", jd_2024_jan_1, jd_2024_end, None)
            .expect("search_stations_in_range Mercury");
        assert!(stations.len() >= 6);
        assert!(stations.iter().any(|s| s.kind == "StationRetrograde"));
        assert!(stations.iter().any(|s| s.kind == "StationDirect"));

        // 10. search_declination_extrema_in_range
        let dec_extrema = engine
            .search_declination_extrema_in_range_core(
                "Moon",
                jd_2024_jan_1,
                jd_2024_jan_1 + 30.0,
                None,
            )
            .expect("search_declination_extrema_in_range Moon");
        assert!(dec_extrema.len() >= 2);
        assert!(dec_extrema.iter().any(|e| e.kind == "NorthernMaximum" && e.dec_deg > 18.0));
        assert!(dec_extrema.iter().any(|e| e.kind == "SouthernMinimum" && e.dec_deg < -18.0));

        // 11. find_parans
        let obs_london = JsObserver {
            latitude_deg: 51.5074,
            longitude_deg: -0.1278,
            altitude_m: Some(0.0),
        };
        let parans = engine
            .find_parans_core(
                JD_J2000,
                obs_london,
                Some(10.0),
                Some(vec!["Sun".to_owned(), "Moon".to_owned(), "Mars".to_owned()]),
                Some(vec![
                    "Regulus".to_owned(),
                    "Sirius".to_owned(),
                    "Spica".to_owned(),
                    "Aldebaran".to_owned(),
                ]),
            )
            .expect("find_parans");
        for p in &parans {
            assert!(p.time_orb_minutes <= 10.0);
            assert_eq!(p.paran_type, p.star_angle);
        }
    }

    #[test]
    fn test_nodejs_wave_h_almanac_aspects_and_snapshot() {
        // Zero-arg constructor (from_env_or_moshier)
        let engine = JsEngine::new_core(None, None, None, None).expect("new_core zero-arg");
        let jd_utc = utc_from_local_core(2000, 1, 1, Some(13), Some(0), Some(0.0), Some(1.0), Some(true))
            .expect("utc_from_local_core");
        assert!((jd_utc - 2_451_545.0).abs() < 1e-9);

        let sr = engine
            .planetary_return_core("Sun", JD_J2000, JD_J2000 + 350.0, None)
            .expect("planetary_return_core");
        assert_eq!(sr["body"], "Sun");
        assert!(sr["jdUt1"].as_f64().expect("jdUt1") > JD_J2000 + 360.0);

        let chart_aspects = engine
            .calc_chart_aspects_core(JD_J2000, None, Some("ptolemaic".to_owned()), Some(1.0), None)
            .expect("calc_chart_aspects_core");
        assert!(!chart_aspects.is_empty());

        let obs = JsObserver {
            latitude_deg: 28.6139,
            longitude_deg: 77.2090,
            altitude_m: Some(216.0),
        };
        let panch = engine
            .panchanga_core(JD_J2000, obs, Some("lahiri".to_owned()))
            .expect("panchanga_core");
        assert!(panch["tithi"]["tithiIndex"].as_u64().expect("tithiIndex") >= 1);

        let ph = engine
            .planetary_hours_core(JD_J2000, obs)
            .expect("planetary_hours_core");
        assert_eq!(ph["hours"].as_array().expect("hours").len(), 24);

        let lots = engine
            .calc_lots_core(JD_J2000, obs, None)
            .expect("calc_lots_core");
        assert_eq!(lots.len(), 7);

        let snap = engine
            .chart_snapshot_core(JD_J2000, obs, None, None, None, None, None)
            .expect("chart_snapshot_core");
        assert_eq!(snap["lots"].as_array().expect("lots").len(), 7);
        assert_eq!(snap["bodies"].as_array().expect("bodies").len(), 13);

        let snap_empty = engine
            .chart_snapshot_core(JD_J2000, obs, Some(vec![]), None, None, None, None)
            .expect("chart_snapshot_core empty bodies");
        assert_eq!(snap_empty["bodies"].as_array().expect("bodies").len(), 0);
        assert_eq!(snap_empty["aspects"].as_array().expect("aspects").len(), 0);

        let zp = map_zodiac_placement_core(302.0).expect("map_zodiac_placement_core");
        assert_eq!(zp["hdGate"], 41);
        assert_eq!(zp["hdLine"], 1);
    }
}
