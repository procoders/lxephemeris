//! Complete Python bindings (`PyO3`) for `lxephemeris`.
//!
//! Exposes `Engine`, `SweAdapter`, `SplitJd`, `Observer`, `CalcOptions`, `CalcResult`,
//! `BatchItemResult`, `AyanamsaValue`, `HouseResult`, `HorizontalCoords`, `SolarEclipseInfo`,
//! `LunarEclipseInfo`, `DeltaTInfo`, `SiderealTimeInfo`, `EclNutInfo`, `StarMetadata`,
//! `HeliacalResult`, `julday`, `revjul`, `houses_armc`, `house_of_lon`, and `gauquelin_sectors`,
//! with GIL-free execution (`py.allow_threads`), contiguous flat `f64` buffers for NumPy
//! `[N_epochs, N_bodies, 6]`, and structured `[ErrorCode]` exceptions.

use lxephemeris_core::{
    aspect_definitions as core_aspect_definitions, azalt as core_azalt, body_display_name,
    calc_shodashavarga as core_calc_shodashavarga, calc_varga as core_calc_varga,
    find_star_by_name_or_bayer, gauquelin_sectors as core_gauquelin_sectors,
    house_of_lon as core_house_of_lon, houses_armc_with_sidereal, julday as core_julday,
    list_stars as core_list_stars, map_zodiac_placement as core_map_zodiac_placement,
    revjul as core_revjul, star_info as core_star_info, AbMode, AspectPreset, AtmosphereParams,
    AyanamsaModel, AyanamsaValue, BodyId, CalcOptions, CalcResult, CalendarKind,
    DeclinationExtremumEvent, DeclinationExtremumKind, DeltaTInfo, DiscMode, EclNutInfo,
    EclipseOptions, Engine, EngineConfig, EphemerisError, Equinox, HeliacalOptions, HeliacalResult,
    HipId, HorizontalCoords, HouseOptions, HouseSystem, HousesResult, HypotheticalId,
    LongitudeCrossingEvent, LunarEclipseInfo, LunarOccultationInfo as CoreLunarOccultationInfo,
    LunarPhaseEvent, LunarPhaseKind, LunarPhaseSummary, LunarPointId, MpcNumber,
    NodesApsidesResult, NutationMode, Observer, OccultationKind, Origin, ParanEvent, ParanOptions,
    Plane, PlanetId, PolarPolicy, RiseTransEvent, RiseTransOptions, SeasonsOfYear,
    ShadowEnlargementRule, SiderealTimeInfo, SolarEclipseInfo, SolarTermEvent, SplitJd,
    StarMetadata, StationEvent, StationKind, TimeInput, Zodiac, MINOR_BODIES_68,
    SE_NODBIT_FOPOINT, SE_NODBIT_MEAN, SE_NODBIT_OSCU, SE_NODBIT_OSCU_BAR,
};
#[cfg(feature = "swe-compat")]
use lxephemeris_core::{
    SweAdapter, SEFLG_BARYCTR, SEFLG_EQUATORIAL, SEFLG_HELCTR, SEFLG_J2000, SEFLG_JPLEPH,
    SEFLG_MOSEPH, SEFLG_NOABERR, SEFLG_NOGDEFL, SEFLG_NONUT, SEFLG_RADIANS, SEFLG_SIDEREAL,
    SEFLG_SPEED, SEFLG_SPEED3, SEFLG_SWIEPH, SEFLG_TOPOCTR, SEFLG_TRUEPOS, SEFLG_XYZ, SE_ADMETOS,
    SE_APOLLON, SE_AST_OFFSET, SE_BIT_ASTRO_TWILIGHT, SE_BIT_CIVIL_TWILIGHT, SE_BIT_DISC_BOTTOM,
    SE_BIT_DISC_CENTER, SE_BIT_FIXED_DISC_SIZE, SE_BIT_GEOCTR_NO_ECL_LAT, SE_BIT_HINDU_RISING,
    SE_BIT_NAUTIC_TWILIGHT, SE_BIT_NO_REFRACTION, SE_CALC_ITRANSIT, SE_CALC_MTRANSIT, SE_CALC_RISE,
    SE_CALC_SET, SE_CERES, SE_CHIRON, SE_CUPIDO, SE_EARTH, SE_ECL2HOR, SE_ECL_1ST_VISIBLE,
    SE_ECL_2ND_VISIBLE, SE_ECL_3RD_VISIBLE, SE_ECL_4TH_VISIBLE, SE_ECL_ALLTYPES_LUNAR,
    SE_ECL_ALLTYPES_SOLAR, SE_ECL_ANNULAR, SE_ECL_ANNULAR_TOTAL, SE_ECL_CENTRAL, SE_ECL_HYBRID,
    SE_ECL_MAX_VISIBLE, SE_ECL_NONCENTRAL, SE_ECL_NUT, SE_ECL_ONE_TRY, SE_ECL_PARTIAL,
    SE_ECL_PENUMBRAL, SE_ECL_TOTAL, SE_ECL_VISIBLE, SE_EQU2HOR, SE_EVENING_FIRST, SE_FICT_OFFSET,
    SE_GREG_CAL, SE_HADES, SE_HELFLAG_HIGH_PRECISION, SE_HELIACAL_RISING, SE_HELIACAL_SETTING,
    SE_HOR2ECL, SE_HOR2EQU, SE_INTP_APOG, SE_INTP_PERG, SE_ISIS, SE_JUL_CAL, SE_JUNO, SE_JUPITER,
    SE_KRONOS, SE_MARS, SE_MEAN_APOG, SE_MEAN_NODE, SE_MERCURY, SE_MOON, SE_MORNING_LAST,
    SE_NEPTUNE, SE_OSCU_APOG, SE_PALLAS, SE_PHOLUS, SE_PLUTO, SE_POSEIDON, SE_SATURN,
    SE_SIDBIT_ECL_DATE, SE_SIDBIT_ECL_T0, SE_SIDBIT_NO_PREC_OFFSET, SE_SIDBIT_PREC_ORIG,
    SE_SIDBIT_SSY_PLANE, SE_SIDBIT_USER_UT, SE_SIDM_ALDEBARAN_15TAU, SE_SIDM_ARYABHATA,
    SE_SIDM_ARYABHATA_522, SE_SIDM_ARYABHATA_MSUN, SE_SIDM_B1950, SE_SIDM_BABYL_BRITTON,
    SE_SIDM_BABYL_ETPSC, SE_SIDM_BABYL_HUBER, SE_SIDM_BABYL_KUGLER1, SE_SIDM_BABYL_KUGLER2,
    SE_SIDM_BABYL_KUGLER3, SE_SIDM_DELUCE, SE_SIDM_DJWHAL_KHUL, SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_GALALIGN_MARDYKS, SE_SIDM_GALCENT_0SAG, SE_SIDM_GALCENT_MULA_WILHELM,
    SE_SIDM_GALCENT_RGILBRAND, SE_SIDM_GALEQU_FIORENZA, SE_SIDM_GALEQ_IAU1958, SE_SIDM_GALEQ_MULA,
    SE_SIDM_GALEQ_TRUE, SE_SIDM_HIPPARCHOS, SE_SIDM_J1900, SE_SIDM_J2000, SE_SIDM_JN_BHASIN,
    SE_SIDM_KP_SENTHILATHIBAN, SE_SIDM_KRISHNAMURTI, SE_SIDM_LAHIRI, SE_SIDM_LAHIRI_1940,
    SE_SIDM_LAHIRI_ICRC, SE_SIDM_LAHIRI_VP285, SE_SIDM_RAMAN, SE_SIDM_SASSANIAN, SE_SIDM_SS_CITRA,
    SE_SIDM_SS_REVATI, SE_SIDM_SURYASIDDHANTA, SE_SIDM_SURYASIDDHANTA_MSUN,
    SE_SIDM_TRUE_ALDEBARAN, SE_SIDM_TRUE_CITRA, SE_SIDM_TRUE_MULA, SE_SIDM_TRUE_PUSHYA,
    SE_SIDM_TRUE_REVATI, SE_SIDM_TRUE_SHERATAN, SE_SIDM_USER, SE_SIDM_USHASHASHI,
    SE_SIDM_VALENS_MOON, SE_SIDM_YUKTESHWAR, SE_SUN, SE_TRUE_NODE, SE_URANUS, SE_VENUS, SE_VESTA,
    SE_VULKANUS, SE_ZEUS,
};
use pyo3::create_exception;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyByteArray, PyBytes, PyDict, PyList, PyType};
use serde::Serialize;
use std::path::PathBuf;
#[cfg(feature = "swe-compat")]
use std::sync::Mutex;
use std::sync::{Arc, OnceLock, RwLock};

// Typed exception hierarchy (R10.2): all domain exceptions inherit from LxEphemerisError (which inherits from ValueError)
create_exception!(lxephemeris, LxEphemerisError, PyValueError);
create_exception!(lxephemeris, InvalidInputError, LxEphemerisError);
create_exception!(lxephemeris, UnsupportedBodyError, LxEphemerisError);
create_exception!(lxephemeris, UnsupportedModelError, LxEphemerisError);
create_exception!(lxephemeris, OutOfRangeError, LxEphemerisError);
create_exception!(lxephemeris, PolarSingularityError, LxEphemerisError);
create_exception!(lxephemeris, ConvergenceError, LxEphemerisError);
create_exception!(lxephemeris, MissingDataError, LxEphemerisError);
create_exception!(lxephemeris, CorruptKernelError, LxEphemerisError);

/// Parses a string body identifier (e.g. `"Sun"`, `"Moon"`, `"Mars"`, `"mpc:1"`, `"Chiron"`,
/// `"Mean_Node"`, `"Cupido"`, `"hip:65474"`, `"Regulus"`, `"ECL_NUT"`) into a [`BodyId`].
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
        "interpolated_lilith" | "interpolatedlilith" | "intp_apogee" => {
            Ok(BodyId::LunarPoint(LunarPointId::InterpolatedLilith))
        }
        "interpolated_perigee" | "interpolatedperigee" | "intp_perg" | "intp_perigee" => {
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
                let mpc: u32 = mpc_str.parse().map_err(|_| EphemerisError::UnsupportedBody {
                    body: trimmed.to_owned(),
                })?;
                return Ok(BodyId::MinorPlanet(MpcNumber(mpc)));
            }
            if let Some(hip_str) = other.strip_prefix("hip:") {
                let hip: u32 = hip_str.parse().map_err(|_| EphemerisError::UnsupportedBody {
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

/// Parses a time scale string (`"TT"`, `"UT1"`, `"UTC"`, `"TDB"`) and `(jd1, jd2)` into [`TimeInput`].
pub fn parse_time_input(scale: &str, jd1: f64, jd2: f64) -> Result<TimeInput, EphemerisError> {
    if !jd1.is_finite() || !jd2.is_finite() {
        return Err(EphemerisError::InvalidInput {
            message: format!("Non-finite Julian Date: jd1={jd1}, jd2={jd2}"),
        });
    }
    let split = SplitJd::new(jd1, jd2);
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

fn extract_calendar_kind(calendar: Option<&Bound<'_, PyAny>>) -> PyResult<CalendarKind> {
    let Some(obj) = calendar else {
        return Ok(CalendarKind::Gregorian);
    };
    if let Ok(flag) = obj.extract::<i32>() {
        match flag {
            0 => Ok(CalendarKind::Julian),
            1 => Ok(CalendarKind::Gregorian),
            other => Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                message: format!(
                    "Unknown calendar flag '{other}' (expected 1 for Gregorian or 0 for Julian)"
                ),
            })),
        }
    } else if let Ok(s) = obj.extract::<String>() {
        parse_calendar_kind(&s).map_err(ephemeris_err_to_pyerr)
    } else {
        Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
            message: "calendar must be an int (1=Gregorian, 0=Julian) or str ('gregorian', 'julian', 'auto')".to_owned(),
        }))
    }
}

fn parse_polar_policy(s: &str) -> Result<PolarPolicy, EphemerisError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "strict" => Ok(PolarPolicy::Strict),
        "fallback" | "fallback_to_porphyry" | "fallback_porphyry" | "porphyry" => {
            Ok(PolarPolicy::FallbackToPorphyry)
        }
        "fallback_whole_sign" | "fallback_to_whole_sign" | "whole_sign" | "wholesign" => {
            Ok(PolarPolicy::FallbackToWholeSign)
        }
        other => Err(EphemerisError::InvalidInput {
            message: format!(
                "Unknown polar_policy '{other}' (expected 'strict', 'fallback_to_porphyry', or 'fallback_whole_sign')"
            ),
        }),
    }
}

fn parse_nutation_mode(s: &str) -> Result<NutationMode, EphemerisError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "true" => Ok(NutationMode::True),
        "mean" => Ok(NutationMode::Mean),
        other => Err(EphemerisError::InvalidInput {
            message: format!("Unknown nutation mode '{other}' (expected 'true' or 'mean')"),
        }),
    }
}

fn parse_shadow_rule(s: &str) -> Result<ShadowEnlargementRule, EphemerisError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "danjon" | "danjon_1_over_85" => Ok(ShadowEnlargementRule::Danjon1Over85),
        "traditional" | "traditional_2_percent" => Ok(ShadowEnlargementRule::Traditional2Percent),
        other => Err(EphemerisError::InvalidInput {
            message: format!(
                "Unknown shadow_rule '{other}' (expected 'danjon' or 'traditional')"
            ),
        }),
    }
}

fn extract_raw_hsys_char(hsys: Option<&Bound<'_, PyAny>>) -> PyResult<char> {
    let Some(obj) = hsys else {
        return Ok('P');
    };
    if obj.is_none() {
        return Ok('P');
    }
    if let Ok(s) = obj.extract::<String>() {
        s.chars().next().ok_or_else(|| {
            ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                message: "Empty house system code".to_owned(),
            })
        })
    } else if let Ok(b) = obj.downcast::<PyBytes>() {
        let slice = b.as_bytes();
        slice.first().map(|&c| c as char).ok_or_else(|| {
            ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                message: "Empty house system bytes".to_owned(),
            })
        })
    } else if let Ok(ba) = obj.downcast::<PyByteArray>() {
        let slice = unsafe { ba.as_bytes() };
        slice.first().map(|&c| c as char).ok_or_else(|| {
            ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                message: "Empty house system bytearray".to_owned(),
            })
        })
    } else if let Ok(code) = obj.extract::<u8>() {
        Ok(code as char)
    } else {
        Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
            message: "hsys must be a 1-character str, bytes, bytearray, or ASCII int (e.g. 'P' or b'P')".to_owned(),
        }))
    }
}

fn extract_hsys_char(hsys: Option<&Bound<'_, PyAny>>) -> PyResult<char> {
    let ch = extract_raw_hsys_char(hsys)?;
    Ok(if ch == 'E' { 'A' } else { ch })
}

#[cfg(feature = "swe-compat")]
fn py_cotrans_impl(
    a0: &Bound<'_, PyAny>,
    a1: f64,
    a2: Option<f64>,
    a3: Option<f64>,
) -> PyResult<(f64, f64, f64)> {
    if let Ok(lon) = a0.extract::<f64>() {
        let lat = a1;
        let dist = a2.ok_or_else(|| {
            PyValueError::new_err(
                "cotrans(lon, lat, dist, eps) requires 4 numeric arguments or (xpo, eps)",
            )
        })?;
        let eps = a3.ok_or_else(|| {
            PyValueError::new_err(
                "cotrans(lon, lat, dist, eps) requires 4 numeric arguments or (xpo, eps)",
            )
        })?;
        Ok(SweAdapter::cotrans(lon, lat, dist, eps))
    } else {
        let seq: Vec<f64> = a0.extract().map_err(|_| {
            PyValueError::new_err(
                "cotrans first argument must be a (lon, lat[, dist]) sequence or float lon",
            )
        })?;
        let Some(&lon) = seq.first() else {
            return Err(PyValueError::new_err(
                "cotrans coordinate sequence must have at least 2 elements (lon, lat[, dist])",
            ));
        };
        let Some(&lat) = seq.get(1) else {
            return Err(PyValueError::new_err(
                "cotrans coordinate sequence must have at least 2 elements (lon, lat[, dist])",
            ));
        };
        let dist = seq.get(2).copied().unwrap_or(1.0);
        let eps = a1;
        Ok(SweAdapter::cotrans(lon, lat, dist, eps))
    }
}

#[cfg(feature = "swe-compat")]
fn reconstruct_pullen_q(asc: f64, mc: f64) -> [f64; 12] {
    let ic = (mc + 180.0).rem_euclid(360.0);
    let dsc = (asc + 180.0).rem_euclid(360.0);
    let q4 = (asc - mc).rem_euclid(360.0);
    let q1 = (ic - asc).rem_euclid(360.0);
    let q_s = q4.min(q1);
    let (w_s, w_l) = if q_s < 1e-9 {
        (0.0, 60.0)
    } else {
        let ratio_q = (180.0 - q_s) / q_s;
        let mut r = (2.0 * ratio_q).powf(1.0 / 3.0) + 1.0;
        if r < 1.0 {
            r = 1.0;
        }
        for _ in 0..30 {
            let r2 = r * r;
            let r3 = r2 * r;
            let f_v = r3 * (r + 2.0) - ratio_q * (2.0 * r + 1.0);
            let df_v = 4.0 * r3 + 6.0 * r2 - 2.0 * ratio_q;
            if df_v.abs() < 1e-14 {
                break;
            }
            r -= f_v / df_v;
        }
        let xb = q_s / (2.0 * r + 1.0);
        (r * xb, r * r * r * xb)
    };
    let (w4, w1) = if q4 <= q1 { (w_s, w_l) } else { (w_l, w_s) };
    let c11 = (mc + w4).rem_euclid(360.0);
    let c12 = (asc - w4).rem_euclid(360.0);
    let c2 = (asc + w1).rem_euclid(360.0);
    let c3 = (ic - w1).rem_euclid(360.0);
    [
        asc,
        c2,
        c3,
        ic,
        (c11 + 180.0).rem_euclid(360.0),
        (c12 + 180.0).rem_euclid(360.0),
        dsc,
        (c2 + 180.0).rem_euclid(360.0),
        (c3 + 180.0).rem_euclid(360.0),
        mc,
        c11,
        c12,
    ]
}

#[cfg(feature = "swe-compat")]
fn reconstruct_whole_sign(asc: f64) -> [f64; 12] {
    let asc_snapped = if (asc - 360.0).abs() < 1e-9 { 0.0 } else { asc };
    let sign_start = (((asc_snapped + 1e-11).rem_euclid(360.0)) / 30.0).floor() * 30.0;
    std::array::from_fn(|i| (sign_start + 30.0 * (i as f64)).rem_euclid(360.0))
}

#[cfg(feature = "swe-compat")]
#[allow(clippy::indexing_slicing)]
fn fix_polar_houses(
    code: char,
    mut cusps: [f64; 12],
    mut ascmc: [f64; 8],
    lat: f64,
    eps_deg: f64,
    sid_offset: f64,
) -> ([f64; 12], [f64; 8]) {
    if code == 'I' && lat.abs() > 60.0 {
        if (((cusps[1] - cusps[9] + 180.0).rem_euclid(360.0)) - 180.0).abs() < 1e-4 {
            let ic_val = cusps[3];
            cusps[1] = ic_val;
            cusps[2] = ic_val;
            cusps[4] = ic_val;
            cusps[5] = ic_val;
        }
        if (((cusps[7] - cusps[3] + 180.0).rem_euclid(360.0)) - 180.0).abs() < 1e-4 {
            let mc_val = cusps[9];
            cusps[7] = mc_val;
            cusps[8] = mc_val;
            cusps[10] = mc_val;
            cusps[11] = mc_val;
        }
    }
    let r_rad = ((ascmc[2] + 90.0).rem_euclid(360.0)).to_radians();
    let e_rad = eps_deg.to_radians();
    let p_rad = lat.clamp(-89.999_999, 89.999_999).to_radians();
    let raw_asc = (r_rad
        .sin()
        .atan2(r_rad.cos() * e_rad.cos() - p_rad.tan() * e_rad.sin())
        .to_degrees()
        - sid_offset)
        .rem_euclid(360.0);
    let raw_inverted = (raw_asc - ascmc[1]).rem_euclid(360.0) >= 180.0;
    if !raw_inverted {
        if code == 'Q' {
            cusps = reconstruct_pullen_q(ascmc[0], ascmc[1]);
        } else if code == 'W' {
            cusps = reconstruct_whole_sign(ascmc[0]);
        }
        return (cusps, ascmc);
    }
    let already_flipped_asc = (ascmc[0] - ascmc[1]).rem_euclid(360.0) < 180.0;
    if !already_flipped_asc {
        ascmc[0] = (ascmc[0] + 180.0).rem_euclid(360.0);
    }
    let asc = ascmc[0];
    let mc = ascmc[1];
    let ic = (mc + 180.0).rem_euclid(360.0);
    let dsc = (asc + 180.0).rem_euclid(360.0);
    match code {
        'W' => {
            cusps = reconstruct_whole_sign(asc);
        }
        'R' | 'T' | 'Y' | 'I' => {
            if already_flipped_asc {
                for idx in [1, 2, 3, 4, 5, 7, 8, 9, 10, 11] {
                    cusps[idx] = (cusps[idx] + 180.0).rem_euclid(360.0);
                }
            } else {
                for c in &mut cusps {
                    *c = (*c + 180.0).rem_euclid(360.0);
                }
            }
            ascmc[1] = (ascmc[1] + 180.0).rem_euclid(360.0);
            if code == 'I' {
                if (((cusps[1] - cusps[9] + 180.0).rem_euclid(360.0)) - 180.0).abs() < 1e-4 {
                    let ic_val = cusps[3];
                    cusps[1] = ic_val;
                    cusps[2] = ic_val;
                    cusps[4] = ic_val;
                    cusps[5] = ic_val;
                }
                if (((cusps[7] - cusps[3] + 180.0).rem_euclid(360.0)) - 180.0).abs() < 1e-4 {
                    let mc_val = cusps[9];
                    cusps[7] = mc_val;
                    cusps[8] = mc_val;
                    cusps[10] = mc_val;
                    cusps[11] = mc_val;
                }
            }
        }
        'A' | 'E' | 'V' | 'F' => {
            if !already_flipped_asc {
                for c in &mut cusps {
                    *c = (*c + 180.0).rem_euclid(360.0);
                }
            }
        }
        'C' => {
            ascmc[1] = (ascmc[1] + 180.0).rem_euclid(360.0);
            if already_flipped_asc {
                for idx in [3, 9] {
                    cusps[idx] = (cusps[idx] + 180.0).rem_euclid(360.0);
                }
            } else {
                for idx in [0, 3, 6, 9] {
                    cusps[idx] = (cusps[idx] + 180.0).rem_euclid(360.0);
                }
            }
        }
        'U' => {
            if !already_flipped_asc {
                cusps = [
                    cusps[6], cusps[5], cusps[4], cusps[3], cusps[2], cusps[1], cusps[0],
                    cusps[11], cusps[10], cusps[9], cusps[8], cusps[7],
                ];
            }
        }
        'O' | 'S' => {
            let q1 = (asc - mc).rem_euclid(360.0) / 3.0;
            let q2 = (ic - asc).rem_euclid(360.0) / 3.0;
            let p = [
                asc,
                (asc + q2).rem_euclid(360.0),
                (asc + 2.0 * q2).rem_euclid(360.0),
                ic,
                (mc + q1 + 180.0).rem_euclid(360.0),
                (mc + 2.0 * q1 + 180.0).rem_euclid(360.0),
                dsc,
                (asc + q2 + 180.0).rem_euclid(360.0),
                (asc + 2.0 * q2 + 180.0).rem_euclid(360.0),
                mc,
                (mc + q1).rem_euclid(360.0),
                (mc + 2.0 * q1).rem_euclid(360.0),
            ];
            if code == 'O' {
                cusps = p;
            } else {
                cusps = std::array::from_fn(|i| {
                    (p[i] - ((p[i] - p[(i + 11) % 12]).rem_euclid(360.0) / 2.0)).rem_euclid(360.0)
                });
            }
        }
        'L' => {
            let q4 = (asc - mc).rem_euclid(360.0);
            let q1 = (ic - asc).rem_euclid(360.0);
            let w4 = if q4 >= 30.0 {
                (q4 + 30.0) / 4.0
            } else {
                q4 / 2.0
            };
            let w1 = if q1 >= 30.0 {
                (q1 + 30.0) / 4.0
            } else {
                q1 / 2.0
            };
            let c11 = (mc + w4).rem_euclid(360.0);
            let c12 = (asc - w4).rem_euclid(360.0);
            let c2 = (asc + w1).rem_euclid(360.0);
            let c3 = (ic - w1).rem_euclid(360.0);
            cusps = [
                asc,
                c2,
                c3,
                ic,
                (c11 + 180.0).rem_euclid(360.0),
                (c12 + 180.0).rem_euclid(360.0),
                dsc,
                (c2 + 180.0).rem_euclid(360.0),
                (c3 + 180.0).rem_euclid(360.0),
                mc,
                c11,
                c12,
            ];
        }
        'Q' => {
            cusps = reconstruct_pullen_q(asc, mc);
        }
        'B' => {
            let eps_rad = eps_deg.to_radians();
            let asc_trop = (asc + sid_offset).rem_euclid(360.0);
            let dec_asc = (asc_trop.to_radians().sin() * eps_rad.sin())
                .clamp(-1.0, 1.0)
                .asin();
            let sn = (lat.to_radians().tan() * dec_asc.tan()).clamp(-1.0, 1.0);
            let ad = sn.asin().to_degrees();
            let sd = (90.0 + ad) / 3.0;
            let sna = (90.0 - ad) / 3.0;
            let armc = ascmc[2];
            let eq_to_ecl = |ra: f64| -> f64 {
                let r_rad = ra.to_radians();
                (r_rad.sin().atan2(r_rad.cos() * eps_rad.cos()).to_degrees() - sid_offset)
                    .rem_euclid(360.0)
            };
            let c11 = eq_to_ecl(armc + sd);
            let c12 = eq_to_ecl(armc + 2.0 * sd);
            let c2 = eq_to_ecl(armc + 180.0 - 2.0 * sna);
            let c3 = eq_to_ecl(armc + 180.0 - sna);
            cusps = [
                asc,
                c2,
                c3,
                ic,
                (c11 + 180.0).rem_euclid(360.0),
                (c12 + 180.0).rem_euclid(360.0),
                dsc,
                (c2 + 180.0).rem_euclid(360.0),
                (c3 + 180.0).rem_euclid(360.0),
                mc,
                c11,
                c12,
            ];
        }
        _ => {}
    }
    (cusps, ascmc)
}

fn parse_disc_mode(s: &str) -> Result<DiscMode, EphemerisError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "upper_limb" | "upper" | "limb" => Ok(DiscMode::UpperLimb),
        "center" | "disc_center" => Ok(DiscMode::Center),
        "bottom_limb" | "bottom" | "disc_bottom" => Ok(DiscMode::BottomLimb),
        other => Err(EphemerisError::InvalidInput {
            message: format!(
                "Unknown disc_mode '{other}' (expected 'upper_limb', 'center', or 'bottom_limb')"
            ),
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

fn parse_aspect_preset(s: &str) -> Result<AspectPreset, EphemerisError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "major" | "ptolemaic" => Ok(AspectPreset::Major),
        "standard" | "minor" => Ok(AspectPreset::Standard),
        "kepler" => Ok(AspectPreset::Kepler),
        "all" => Ok(AspectPreset::All),
        other => Err(EphemerisError::InvalidInput {
            message: format!(
                "Unknown aspect preset '{other}' (expected 'major', 'standard', 'kepler', or 'all')"
            ),
        }),
    }
}

fn resolve_wave_h_zodiac(zodiac: &str, ayanamsa: Option<&str>) -> Result<Zodiac, EphemerisError> {
    if let Some(aya) = ayanamsa {
        let trimmed = aya.trim();
        if !trimmed.is_empty() {
            let m = AyanamsaModel::from_code(trimmed)?;
            return Ok(Zodiac::Sidereal(m));
        }
    }
    parse_zodiac(zodiac, "lahiri")
}

fn resolve_wave_h_calc_options(
    zodiac: &str,
    ayanamsa: Option<&str>,
) -> Result<CalcOptions, EphemerisError> {
    let zod = resolve_wave_h_zodiac(zodiac, ayanamsa)?;
    Ok(CalcOptions {
        zodiac: zod,
        ..CalcOptions::default()
    })
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
        BodyId::LunarPoint(LunarPointId::TrueNode),
        BodyId::LunarPoint(LunarPointId::MeanLilith),
        BodyId::MinorPlanet(MpcNumber(2060)),
    ]
}

fn default_major_planets_10() -> Vec<BodyId> {
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
    ]
}

fn parse_optional_bodies(
    bodies: Option<Vec<String>>,
    default_fn: fn() -> Vec<BodyId>,
) -> Result<Vec<BodyId>, EphemerisError> {
    match bodies {
        Some(names) => names
            .iter()
            .map(|n| parse_body_id(n))
            .collect::<Result<Vec<_>, _>>(),
        None => Ok(default_fn()),
    }
}

static ATTR_DICT_TYPE: OnceLock<Py<PyType>> = OnceLock::new();

fn get_attr_dict_type<'py>(py: Python<'py>) -> PyResult<&'py Bound<'py, PyType>> {
    if let Some(t) = ATTR_DICT_TYPE.get() {
        return Ok(t.bind(py));
    }
    let locals = PyDict::new(py);
    py.run(
        pyo3::ffi::c_str!(
            "class AttrDict(dict):\n\
             \x20   def __getattr__(self, k):\n\
             \x20       try:\n\
             \x20           return self[k]\n\
             \x20       except KeyError as e:\n\
             \x20           raise AttributeError(k) from e\n\
             \x20   __setattr__ = dict.__setitem__\n\
             \x20   __delattr__ = dict.__delitem__\n"
        ),
        None,
        Some(&locals),
    )?;
    let cls = locals
        .get_item("AttrDict")?
        .ok_or_else(|| PyRuntimeError::new_err("Failed to define AttrDict"))?
        .downcast_into::<PyType>()?;
    let _ = ATTR_DICT_TYPE.set(cls.clone().unbind());
    ATTR_DICT_TYPE
        .get()
        .map(|t| t.bind(py))
        .ok_or_else(|| PyRuntimeError::new_err("Failed to cache AttrDict"))
}

fn dict_to_attr_dict(py: Python<'_>, d: &Bound<'_, PyDict>) -> PyResult<PyObject> {
    let cls = get_attr_dict_type(py)?;
    Ok(cls.call1((d,))?.into_any().unbind())
}

fn json_value_to_attr_dict(py: Python<'_>, val: &serde_json::Value) -> PyResult<PyObject> {
    match val {
        serde_json::Value::Null => Ok(py.None()),
        serde_json::Value::Bool(b) => Ok(b.into_pyobject(py)?.to_owned().into_any().unbind()),
        serde_json::Value::Number(n) => {
            if let Some(i) = n.as_i64() {
                Ok(i.into_pyobject(py)?.into_any().unbind())
            } else if let Some(u) = n.as_u64() {
                Ok(u.into_pyobject(py)?.into_any().unbind())
            } else if let Some(f) = n.as_f64() {
                Ok(f.into_pyobject(py)?.into_any().unbind())
            } else {
                Ok(py.None())
            }
        }
        serde_json::Value::String(s) => Ok(s.as_str().into_pyobject(py)?.into_any().unbind()),
        serde_json::Value::Array(arr) => {
            let list = PyList::empty(py);
            for item in arr {
                list.append(json_value_to_attr_dict(py, item)?)?;
            }
            Ok(list.into_any().unbind())
        }
        serde_json::Value::Object(map) => {
            if map.len() == 1 {
                if let Ok(bid) = serde_json::from_value::<BodyId>(val.clone()) {
                    return Ok(body_display_name(bid).into_pyobject(py)?.into_any().unbind());
                }
            }
            let d = PyDict::new(py);
            for (k, v) in map {
                d.set_item(k.as_str(), json_value_to_attr_dict(py, v)?)?;
            }
            dict_to_attr_dict(py, &d)
        }
    }
}

fn serialize_to_attr_dict<T: Serialize + ?Sized>(py: Python<'_>, value: &T) -> PyResult<PyObject> {
    let json_val = serde_json::to_value(value).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    json_value_to_attr_dict(py, &json_val)
}

/// Maps an [`EphemerisError`] into a structured [`PyErr`] with `[ErrorCode]` prefix and typed exception class (`R10.2`).
pub fn ephemeris_err_to_pyerr(err: EphemerisError) -> PyErr {
    let code = err.code().as_str();
    let msg = format!("[{code}] {err}");
    match err {
        EphemerisError::InvalidInput { .. } => InvalidInputError::new_err(msg),
        EphemerisError::UnsupportedBody { .. } => UnsupportedBodyError::new_err(msg),
        EphemerisError::UnsupportedModel { .. } => UnsupportedModelError::new_err(msg),
        EphemerisError::OutOfCoverage { .. } | EphemerisError::OutOfModelValidity { .. } => {
            OutOfRangeError::new_err(msg)
        }
        EphemerisError::UndefinedHouseSystem { .. } => PolarSingularityError::new_err(msg),
        EphemerisError::NonConvergence { .. } => ConvergenceError::new_err(msg),
        EphemerisError::MissingData { .. } => MissingDataError::new_err(msg),
        EphemerisError::CorruptKernel { .. } | EphemerisError::UnsupportedSpkType { .. } => {
            CorruptKernelError::new_err(msg)
        }
        _ => LxEphemerisError::new_err(msg),
    }
}

/// Two-part Julian Date `(jd1, jd2)` preserving sub-microsecond precision.
#[pyclass(name = "SplitJd", frozen)]
#[derive(Debug, Clone, Copy, Serialize)]
pub struct PySplitJd {
    #[pyo3(get)]
    pub jd1: f64,
    #[pyo3(get)]
    pub jd2: f64,
}

#[pymethods]
impl PySplitJd {
    #[new]
    #[pyo3(signature = (jd1, jd2=0.0))]
    #[must_use]
    pub fn new(jd1: f64, jd2: f64) -> Self {
        let s = SplitJd::new(jd1, jd2);
        Self {
            jd1: s.jd1,
            jd2: s.jd2,
        }
    }

    #[must_use]
    pub const fn to_f64(&self) -> f64 {
        self.jd1 + self.jd2
    }

    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!("SplitJd(jd1={:.6}, jd2={:.9})", self.jd1, self.jd2)
    }
}

/// Geodetic observer on the WGS84 reference ellipsoid.
#[pyclass(name = "Observer", frozen)]
#[derive(Debug, Clone, Copy, Serialize)]
pub struct PyObserver {
    #[pyo3(get)]
    pub latitude_deg: f64,
    #[pyo3(get)]
    pub longitude_deg: f64,
    #[pyo3(get)]
    pub altitude_m: f64,
}

#[pymethods]
impl PyObserver {
    #[new]
    #[pyo3(signature = (latitude_deg, longitude_deg, altitude_m=0.0))]
    pub fn new(latitude_deg: f64, longitude_deg: f64, altitude_m: f64) -> PyResult<Self> {
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(Self {
            latitude_deg: obs.latitude_deg,
            longitude_deg: obs.longitude_deg,
            altitude_m: obs.altitude_m,
        })
    }

    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "Observer(latitude_deg={:.4}, longitude_deg={:.4}, altitude_m={:.1})",
            self.latitude_deg, self.longitude_deg, self.altitude_m
        )
    }
}

/// Calculation options for `PyEngine.calc` and `PyEngine.calc_batch`.
#[pyclass(name = "CalcOptions", frozen)]
#[derive(Debug, Clone, Copy)]
pub struct PyCalcOptions {
    pub inner: CalcOptions,
}

#[pymethods]
impl PyCalcOptions {
    #[new]
    #[pyo3(signature = (origin="geocentric", plane="ecliptic", equinox="of_date", nutation="true", aberration="apparent", zodiac="tropical", ayanamsa="lahiri", include_speed=true, observer=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        origin: &str,
        plane: &str,
        equinox: &str,
        nutation: &str,
        aberration: &str,
        zodiac: &str,
        ayanamsa: &str,
        include_speed: bool,
        observer: Option<PyObserver>,
    ) -> PyResult<Self> {
        let orig = match origin.to_ascii_lowercase().as_str() {
            "geocentric" => Origin::Geocentric,
            "topocentric" => {
                let obs = observer.ok_or_else(|| {
                    ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                        message: "origin='topocentric' requires an Observer".to_owned(),
                    })
                })?;
                Origin::Topocentric(Observer {
                    latitude_deg: obs.latitude_deg,
                    longitude_deg: obs.longitude_deg,
                    altitude_m: obs.altitude_m,
                })
            }
            "heliocentric" => Origin::Heliocentric,
            "barycentric" => Origin::Barycentric,
            other => {
                return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: format!("Unknown origin '{other}'"),
                }))
            }
        };

        let pl = match plane.to_ascii_lowercase().as_str() {
            "ecliptic" => Plane::Ecliptic,
            "equatorial" => Plane::Equatorial,
            "horizontal" => Plane::Horizontal,
            other => {
                return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: format!("Unknown plane '{other}'"),
                }))
            }
        };

        let eq = match equinox.to_ascii_lowercase().as_str() {
            "of_date" | "ofdate" => Equinox::OfDate,
            "j2000" | "icrf" => Equinox::J2000,
            other => {
                return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: format!("Unknown equinox '{other}'"),
                }))
            }
        };

        let nut = match nutation.to_ascii_lowercase().as_str() {
            "true" => NutationMode::True,
            "mean" => NutationMode::Mean,
            other => {
                return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: format!("Unknown nutation mode '{other}'"),
                }))
            }
        };

        let ab_mode = match aberration.to_ascii_lowercase().as_str() {
            "apparent" => AbMode::Apparent,
            "no_aberration" | "noaberr" => AbMode::NoAberration,
            "no_deflection" | "nogdefl" => AbMode::NoDeflection,
            "astrometric" => AbMode::Astrometric,
            "geometric" => AbMode::Geometric,
            other => {
                return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: format!("Unknown aberration mode '{other}'"),
                }))
            }
        };

        let zod = parse_zodiac(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;

        Ok(Self {
            inner: CalcOptions {
                origin: orig,
                plane: pl,
                equinox: eq,
                nutation: nut,
                aberration: ab_mode,
                zodiac: zod,
                include_speed,
            },
        })
    }
}

/// Single-body astrometric calculation result exposed to Python.
#[pyclass(name = "CalcResult", frozen)]
#[derive(Debug, Clone, Serialize)]
pub struct PyCalcResult {
    #[pyo3(get)]
    pub target_id: i32,
    #[pyo3(get)]
    pub center_id: i32,
    #[pyo3(get)]
    pub jd_tt: f64,
    #[pyo3(get)]
    pub jd_ut1: f64,
    #[pyo3(get)]
    pub jd_tdb: f64,
    #[pyo3(get)]
    pub delta_t_sec: f64,
    #[pyo3(get)]
    pub time_extrapolated: bool,
    #[pyo3(get)]
    pub time_table_version: String,
    #[pyo3(get)]
    pub lon_deg: f64,
    #[pyo3(get)]
    pub lat_deg: f64,
    #[pyo3(get)]
    pub dist_au: f64,
    #[pyo3(get)]
    pub dist_km: f64,
    #[pyo3(get)]
    pub lon_speed_deg_day: f64,
    #[pyo3(get)]
    pub lat_speed_deg_day: f64,
    #[pyo3(get)]
    pub dist_speed_au_day: f64,
    #[pyo3(get)]
    pub dist_speed_km_day: f64,
    #[pyo3(get)]
    pub magnitude_v: Option<f64>,
    #[pyo3(get)]
    pub cartesian_au: [f64; 3],
    #[pyo3(get)]
    pub cartesian_vel_au_day: [f64; 3],
    #[pyo3(get)]
    pub warnings: Vec<String>,
}

#[pymethods]
impl PyCalcResult {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "CalcResult(target_id={}, lon_deg={:.6}, lat_deg={:.6}, dist_au={:.6}, lon_speed_deg_day={:.6})",
            self.target_id, self.lon_deg, self.lat_deg, self.dist_au, self.lon_speed_deg_day
        )
    }
}

impl From<CalcResult> for PyCalcResult {
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
            cartesian_au: r.cartesian_au,
            cartesian_vel_au_day: r.cartesian_vel_au_day,
            warnings: r.warnings,
        }
    }
}

/// Per-item batch result wrapper for non-aborting batch queries.
#[pyclass(name = "BatchItemResult", frozen)]
#[derive(Debug, Clone, Serialize)]
pub struct PyBatchItemResult {
    #[pyo3(get)]
    pub ok: bool,
    #[pyo3(get)]
    pub result: Option<PyCalcResult>,
    #[pyo3(get)]
    pub error_code: Option<String>,
    #[pyo3(get)]
    pub error_message: Option<String>,
}

#[pymethods]
impl PyBatchItemResult {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "BatchItemResult(ok={}, error_code={:?})",
            self.ok, self.error_code
        )
    }
}

/// $\Delta T = \text{TT} - \text{UT1}$ and $\text{TT} - \text{TDB}$ result (`R1.2`).
#[pyclass(name = "DeltaTInfo", frozen)]
#[derive(Debug, Clone, Serialize)]
pub struct PyDeltaTInfo {
    #[pyo3(get)]
    pub delta_t_sec: f64,
    #[pyo3(get)]
    pub tt_minus_tdb_sec: f64,
    #[pyo3(get)]
    pub extrapolated: bool,
    #[pyo3(get)]
    pub table_version: String,
    #[pyo3(get)]
    pub warning: Option<String>,
}

#[pymethods]
impl PyDeltaTInfo {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "DeltaTInfo(delta_t_sec={:.4}, extrapolated={})",
            self.delta_t_sec, self.extrapolated
        )
    }
}

impl From<DeltaTInfo> for PyDeltaTInfo {
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
#[pyclass(name = "SiderealTimeInfo", frozen)]
#[derive(Debug, Clone, Copy, Serialize)]
pub struct PySiderealTimeInfo {
    #[pyo3(get)]
    pub gmst_hours: f64,
    #[pyo3(get)]
    pub gast_hours: f64,
    #[pyo3(get)]
    pub local_sidereal_hours: f64,
    #[pyo3(get)]
    pub armc_deg: f64,
}

#[pymethods]
impl PySiderealTimeInfo {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "SiderealTimeInfo(gmst_hours={:.6}, gast_hours={:.6}, armc_deg={:.6})",
            self.gmst_hours, self.gast_hours, self.armc_deg
        )
    }
}

impl From<SiderealTimeInfo> for PySiderealTimeInfo {
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
#[pyclass(name = "EclNutInfo", frozen)]
#[derive(Debug, Clone, Copy, Serialize)]
pub struct PyEclNutInfo {
    #[pyo3(get)]
    pub true_obliquity_deg: f64,
    #[pyo3(get)]
    pub mean_obliquity_deg: f64,
    #[pyo3(get)]
    pub nutation_lon_deg: f64,
    #[pyo3(get)]
    pub nutation_obl_deg: f64,
}

#[pymethods]
impl PyEclNutInfo {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "EclNutInfo(true_obliquity_deg={:.6}, mean_obliquity_deg={:.6})",
            self.true_obliquity_deg, self.mean_obliquity_deg
        )
    }
}

impl From<EclNutInfo> for PyEclNutInfo {
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
#[pyclass(name = "StarMetadata", frozen)]
#[derive(Debug, Clone, Serialize)]
pub struct PyStarMetadata {
    #[pyo3(get)]
    pub hip: u32,
    #[pyo3(get)]
    pub canonical_name: String,
    #[pyo3(get)]
    pub bayer: String,
    #[pyo3(get)]
    pub vmag: f64,
    #[pyo3(get)]
    pub ra_rad: f64,
    #[pyo3(get)]
    pub dec_rad: f64,
    #[pyo3(get)]
    pub pm_ra_mas_yr: f64,
    #[pyo3(get)]
    pub pm_dec_mas_yr: f64,
    #[pyo3(get)]
    pub plx_mas: f64,
    #[pyo3(get)]
    pub rv_kms: f64,
    #[pyo3(get)]
    pub rv_missing: bool,
}

#[pymethods]
impl PyStarMetadata {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "StarMetadata(canonical_name='{}', hip={}, vmag={:.2})",
            self.canonical_name, self.hip, self.vmag
        )
    }
}

impl From<StarMetadata> for PyStarMetadata {
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
#[pyclass(name = "HeliacalResult", frozen)]
#[derive(Debug, Clone, Serialize)]
pub struct PyHeliacalResult {
    #[pyo3(get)]
    pub target_name: String,
    #[pyo3(get)]
    pub visual_magnitude: f64,
    #[pyo3(get)]
    pub arcus_visionis_deg: f64,
    #[pyo3(get)]
    pub heliacal_rise_jd_ut1: Option<f64>,
    #[pyo3(get)]
    pub morning_set_jd_ut1: Option<f64>,
    #[pyo3(get)]
    pub evening_rise_jd_ut1: Option<f64>,
    #[pyo3(get)]
    pub heliacal_set_jd_ut1: Option<f64>,
}

#[pymethods]
impl PyHeliacalResult {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "HeliacalResult(target_name='{}', heliacal_rise_jd_ut1={:?})",
            self.target_name, self.heliacal_rise_jd_ut1
        )
    }
}

impl From<HeliacalResult> for PyHeliacalResult {
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
#[pyclass(name = "AyanamsaValue", frozen)]
#[derive(Debug, Clone, Serialize)]
pub struct PyAyanamsaValue {
    #[pyo3(get)]
    pub model_code: String,
    #[pyo3(get)]
    pub jd_tt: f64,
    #[pyo3(get)]
    pub mean_deg: f64,
    #[pyo3(get)]
    pub true_deg: f64,
    #[pyo3(get)]
    pub effective_deg: f64,
    #[pyo3(get)]
    pub rate_deg_day: f64,
    #[pyo3(get)]
    pub dpsi_arcsec: f64,
}

#[pymethods]
impl PyAyanamsaValue {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "AyanamsaValue(model_code='{}', effective_deg={:.6})",
            self.model_code, self.effective_deg
        )
    }
}

impl From<AyanamsaValue> for PyAyanamsaValue {
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

/// 12 house cusps + 8 chart angles + angular velocities exposed to Python.
#[pyclass(name = "HouseResult", frozen)]
#[derive(Debug, Clone, Serialize)]
pub struct PyHouseResult {
    #[pyo3(get)]
    pub system_code: String,
    #[pyo3(get)]
    pub effective_system_code: String,
    #[pyo3(get)]
    pub polar_house_fallback: bool,
    #[pyo3(get)]
    pub cusps_deg: [f64; 12],
    #[pyo3(get)]
    pub cusp_speeds_deg_day: [f64; 12],
    #[pyo3(get)]
    pub asc_deg: f64,
    #[pyo3(get)]
    pub mc_deg: f64,
    #[pyo3(get)]
    pub armc_deg: f64,
    #[pyo3(get)]
    pub vertex_deg: f64,
    #[pyo3(get)]
    pub equatorial_asc_deg: f64,
    #[pyo3(get)]
    pub co_asc_koch_deg: f64,
    #[pyo3(get)]
    pub co_asc_munkasey_deg: f64,
    #[pyo3(get)]
    pub polar_asc_deg: f64,
    #[pyo3(get)]
    pub asc_speed_deg_day: f64,
    #[pyo3(get)]
    pub mc_speed_deg_day: f64,
    #[pyo3(get)]
    pub warnings: Vec<String>,
}

#[pymethods]
impl PyHouseResult {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "HouseResult(system_code='{}', asc_deg={:.4}, mc_deg={:.4})",
            self.system_code, self.asc_deg, self.mc_deg
        )
    }
}

impl From<HousesResult> for PyHouseResult {
    fn from(r: HousesResult) -> Self {
        Self {
            system_code: r.requested_house_system.code().to_string(),
            effective_system_code: r.effective_house_system.code().to_string(),
            polar_house_fallback: r.polar_house_fallback,
            cusps_deg: r.cusps_deg,
            cusp_speeds_deg_day: r.cusp_speeds_deg_day,
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
#[pyclass(name = "HorizontalCoords", frozen)]
#[derive(Debug, Clone, Copy, Serialize)]
pub struct PyHorizontalCoords {
    #[pyo3(get)]
    pub azimuth_deg: f64,
    #[pyo3(get)]
    pub altitude_true_deg: f64,
    #[pyo3(get)]
    pub altitude_apparent_deg: f64,
    #[pyo3(get)]
    pub hour_angle_deg: f64,
    #[pyo3(get)]
    pub near_horizon: bool,
}

#[pymethods]
impl PyHorizontalCoords {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "HorizontalCoords(azimuth_deg={:.4}, altitude_apparent_deg={:.4})",
            self.azimuth_deg, self.altitude_apparent_deg
        )
    }
}

impl From<HorizontalCoords> for PyHorizontalCoords {
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

/// Global and local solar eclipse circumstances exposed to Python.
#[pyclass(name = "SolarEclipseInfo", frozen)]
#[derive(Debug, Clone, Serialize)]
pub struct PySolarEclipseInfo {
    #[pyo3(get)]
    pub kind: String,
    #[pyo3(get)]
    pub saros_series: i32,
    #[pyo3(get)]
    pub saros_member: u32,
    #[pyo3(get)]
    pub max_eclipse_jd_tt: f64,
    #[pyo3(get)]
    pub max_eclipse_jd_ut1: f64,
    #[pyo3(get)]
    pub delta_t_sec: f64,
    #[pyo3(get)]
    pub gamma: f64,
    #[pyo3(get)]
    pub u_radius_earth: f64,
    #[pyo3(get)]
    pub penumbra_radius_earth: f64,
    #[pyo3(get)]
    pub magnitude: f64,
    #[pyo3(get)]
    pub central_duration_sec: Option<f64>,
    #[pyo3(get)]
    pub path_width_km: Option<f64>,
    #[pyo3(get)]
    pub central_latitude_deg: Option<f64>,
    #[pyo3(get)]
    pub central_longitude_deg: Option<f64>,
    #[pyo3(get)]
    pub p1_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub u1_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub u2_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub u3_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub u4_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub p4_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub local_visible_at_max: Option<bool>,
    #[pyo3(get)]
    pub local_is_central: Option<bool>,
    #[pyo3(get)]
    pub local_kind: Option<String>,
    #[pyo3(get)]
    pub local_c1_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub local_c2_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub local_max_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub local_c3_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub local_c4_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub local_magnitude: Option<f64>,
    #[pyo3(get)]
    pub local_obscuration: Option<f64>,
    #[pyo3(get)]
    pub local_central_duration_sec: Option<f64>,
    #[pyo3(get)]
    pub local_sun_altitude_deg: Option<f64>,
    #[pyo3(get)]
    pub local_sun_azimuth_deg: Option<f64>,
}

#[pymethods]
impl PySolarEclipseInfo {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "SolarEclipseInfo(kind='{}', max_eclipse_jd_ut1={:.5}, magnitude={:.4})",
            self.kind, self.max_eclipse_jd_ut1, self.magnitude
        )
    }
}

impl From<SolarEclipseInfo> for PySolarEclipseInfo {
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

/// Lunar eclipse circumstances exposed to Python.
#[pyclass(name = "LunarEclipseInfo", frozen)]
#[derive(Debug, Clone, Serialize)]
pub struct PyLunarEclipseInfo {
    #[pyo3(get)]
    pub kind: String,
    #[pyo3(get)]
    pub saros_series: i32,
    #[pyo3(get)]
    pub saros_member: u32,
    #[pyo3(get)]
    pub shadow_rule: String,
    #[pyo3(get)]
    pub max_eclipse_jd_tt: f64,
    #[pyo3(get)]
    pub max_eclipse_jd_ut1: f64,
    #[pyo3(get)]
    pub delta_t_sec: f64,
    #[pyo3(get)]
    pub gamma: f64,
    #[pyo3(get)]
    pub penumbral_magnitude: f64,
    #[pyo3(get)]
    pub umbral_magnitude: f64,
    #[pyo3(get)]
    pub penumbral_semi_duration_min: f64,
    #[pyo3(get)]
    pub partial_semi_duration_min: f64,
    #[pyo3(get)]
    pub total_semi_duration_min: f64,
    #[pyo3(get)]
    pub p1_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub u1_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub u2_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub u3_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub u4_jd_tt: Option<f64>,
    #[pyo3(get)]
    pub p4_jd_tt: Option<f64>,
}

#[pymethods]
impl PyLunarEclipseInfo {
    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "LunarEclipseInfo(kind='{}', max_eclipse_jd_ut1={:.5}, umbral_magnitude={:.4})",
            self.kind, self.max_eclipse_jd_ut1, self.umbral_magnitude
        )
    }
}

impl From<LunarEclipseInfo> for PyLunarEclipseInfo {
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

/// Planetary, solar, lunar, or minor-body phenomena result (`swe_pheno` / `swe_pheno_ut` compatible).
#[pyclass(name = "PhenomenaResult", frozen)]
#[derive(Clone, Debug, Serialize)]
pub struct PyPhenomenaResult {
    #[pyo3(get)]
    pub phase_angle_deg: f64,
    #[pyo3(get)]
    pub phase_illuminated: f64,
    #[pyo3(get)]
    pub elongation_deg: f64,
    #[pyo3(get)]
    pub apparent_diameter_deg: f64,
    #[pyo3(get)]
    pub apparent_magnitude_v: f64,
    #[pyo3(get)]
    pub horizontal_parallax_deg: f64,
    #[pyo3(get)]
    pub geocentric_distance_au: f64,
    #[pyo3(get)]
    pub heliocentric_distance_au: f64,
    #[serde(skip)]
    attr20: [f64; 20],
}

#[pymethods]
impl PyPhenomenaResult {
    pub fn to_swe_attr20(&self) -> [f64; 20] {
        self.attr20
    }

    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "PhenomenaResult(phase_angle_deg={:.4}, phase_illuminated={:.4}, elongation_deg={:.4}, apparent_magnitude_v={:.2})",
            self.phase_angle_deg, self.phase_illuminated, self.elongation_deg, self.apparent_magnitude_v
        )
    }
}

impl From<lxephemeris_core::PhenomenaResult> for PyPhenomenaResult {
    fn from(r: lxephemeris_core::PhenomenaResult) -> Self {
        let attr20 = r.to_swe_attr20();
        let parallax_rad = r.horizontal_parallax_deg.to_radians();
        let geocentric_distance_au = if parallax_rad > 1e-12 {
            (6378.1366 / 149_597_870.7) / parallax_rad.sin()
        } else {
            0.0
        };
        let d = geocentric_distance_au;
        let elong_rad = r.elongation_deg.to_radians();
        let heliocentric_distance_au = (1.0 + d * d - 2.0 * d * elong_rad.cos()).max(0.0).sqrt();
        Self {
            phase_angle_deg: r.phase_angle_deg,
            phase_illuminated: r.phase_illuminated,
            elongation_deg: r.elongation_deg,
            apparent_diameter_deg: r.apparent_diameter_deg,
            apparent_magnitude_v: r.apparent_magnitude_v,
            horizontal_parallax_deg: r.horizontal_parallax_deg,
            geocentric_distance_au,
            heliocentric_distance_au,
            attr20,
        }
    }
}

/// Planetary or lunar orbital nodes and apsides result (`swe_nod_aps`).
#[pyclass(name = "NodApsResult", frozen)]
#[derive(Clone, Debug, Serialize)]
pub struct PyNodApsResult {
    #[pyo3(get)]
    pub ascending_node: [f64; 6],
    #[pyo3(get)]
    pub descending_node: [f64; 6],
    #[pyo3(get)]
    pub periapsis: [f64; 6],
    #[pyo3(get)]
    pub apoapsis: [f64; 6],
}

#[pymethods]
impl PyNodApsResult {
    #[getter]
    pub fn asc_node(&self) -> [f64; 6] {
        self.ascending_node
    }

    #[getter]
    pub fn desc_node(&self) -> [f64; 6] {
        self.descending_node
    }

    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "NodApsResult(asc_node_lon={:.4}, desc_node_lon={:.4}, periapsis_lon={:.4}, apoapsis_lon={:.4})",
            self.ascending_node[0], self.descending_node[0], self.periapsis[0], self.apoapsis[0]
        )
    }
}

impl From<NodesApsidesResult> for PyNodApsResult {
    fn from(r: NodesApsidesResult) -> Self {
        Self {
            ascending_node: r.asc_node,
            descending_node: r.desc_node,
            periapsis: r.periapsis,
            apoapsis: r.apoapsis,
        }
    }
}

/// Lunar occultation event details (`swe_lun_occult_when_glob` / `swe_lun_occult_when_loc`).
#[pyclass(name = "LunarOccultationInfo", frozen)]
#[derive(Clone, Debug, Serialize)]
pub struct PyLunarOccultationInfo {
    #[pyo3(get)]
    pub kind: String,
    #[pyo3(get)]
    pub target: String,
    #[pyo3(get)]
    pub jd_ut1_max: f64,
    #[pyo3(get)]
    pub jd_tt_max: f64,
    #[pyo3(get)]
    pub t1_first_contact_jd_ut1: Option<f64>,
    #[pyo3(get)]
    pub t2_second_contact_jd_ut1: Option<f64>,
    #[pyo3(get)]
    pub t3_third_contact_jd_ut1: Option<f64>,
    #[pyo3(get)]
    pub t4_fourth_contact_jd_ut1: Option<f64>,
    #[pyo3(get)]
    pub gamma: f64,
    #[pyo3(get)]
    pub magnitude: f64,
    #[pyo3(get)]
    pub moon_altitude_deg: Option<f64>,
    #[serde(skip)]
    tret: [f64; 10],
    #[serde(skip)]
    retflag: i32,
}

#[pymethods]
impl PyLunarOccultationInfo {
    pub fn to_swe_glob_tret(&self) -> [f64; 10] {
        self.tret
    }

    pub fn to_swe_tret10(&self) -> [f64; 10] {
        self.tret
    }

    pub fn swe_retflag(&self) -> i32 {
        self.retflag
    }

    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        serialize_to_attr_dict(py, self)
    }

    pub fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "LunarOccultationInfo(kind='{}', target='{}', jd_ut1_max={:.5}, magnitude={:.2})",
            self.kind, self.target, self.jd_ut1_max, self.magnitude
        )
    }
}

impl From<CoreLunarOccultationInfo> for PyLunarOccultationInfo {
    fn from(info: CoreLunarOccultationInfo) -> Self {
        let tret = info.to_swe_tret10();
        let retflag = info.swe_retflag();
        let kind = match info.kind {
            OccultationKind::Total => "Total".to_owned(),
            OccultationKind::Annular => "Annular".to_owned(),
            OccultationKind::Partial => "Partial".to_owned(),
            OccultationKind::Grazing => "Grazing".to_owned(),
            _ => "Other".to_owned(),
        };
        Self {
            kind,
            target: format!("{:?}", info.target),
            jd_ut1_max: info.jd_ut1_max,
            jd_tt_max: info.jd_tt_max,
            t1_first_contact_jd_ut1: info.t1_first_contact_jd_ut1,
            t2_second_contact_jd_ut1: info.t2_inner_immersion_jd_ut1,
            t3_third_contact_jd_ut1: info.t3_inner_emersion_jd_ut1,
            t4_fourth_contact_jd_ut1: info.t4_last_contact_jd_ut1,
            gamma: info.gamma,
            magnitude: info.magnitude,
            moon_altitude_deg: info.moon_altitude_deg,
            tret,
            retflag,
        }
    }
}

/// Thread-safe, immutable Python `Engine` backed by `lxephemeris_core::Engine`.
#[pyclass(name = "Engine", frozen)]
pub struct PyEngine {
    inner: RwLock<Arc<Engine>>,
}

impl PyEngine {
    pub fn get_engine(&self) -> PyResult<Arc<Engine>> {
        self.inner
            .read()
            .map(|g| Arc::clone(&*g))
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }
}

#[pymethods]
impl PyEngine {
    /// Creates a new `Engine`. If `data_dir` is provided, loads from that directory;
    /// otherwise checks `LXEPHEMERIS_DATA_DIR` or `./data`, falling back to fileless Moshier mode.
    #[new]
    #[pyo3(signature = (data_dir=None, max_threads=None, pin_time_table_version=None, verify_manifest_sha256=false))]
    pub fn new(
        py: Python<'_>,
        data_dir: Option<String>,
        max_threads: Option<usize>,
        pin_time_table_version: Option<String>,
        verify_manifest_sha256: bool,
    ) -> PyResult<Self> {
        match data_dir {
            Some(dir) => Self::from_data_dir(
                py,
                dir,
                max_threads,
                pin_time_table_version,
                verify_manifest_sha256,
            ),
            None => {
                let cfg = EngineConfig {
                    max_threads,
                    pin_time_table_version,
                    verify_manifest_sha256,
                };
                let eng = py.allow_threads(|| Engine::from_env_or_moshier(cfg));
                Ok(Self {
                    inner: RwLock::new(Arc::new(eng)),
                })
            }
        }
    }

    /// Loads an `Engine` from a directory containing `manifest.json` and `.bsp` kernels, releasing the GIL.
    #[staticmethod]
    #[pyo3(signature = (data_dir, max_threads=None, pin_time_table_version=None, verify_manifest_sha256=false))]
    pub fn from_data_dir(
        py: Python<'_>,
        data_dir: String,
        max_threads: Option<usize>,
        pin_time_table_version: Option<String>,
        verify_manifest_sha256: bool,
    ) -> PyResult<Self> {
        let cfg = EngineConfig {
            max_threads,
            pin_time_table_version,
            verify_manifest_sha256,
        };
        let eng = py
            .allow_threads(|| Engine::from_data_dir(PathBuf::from(data_dir), cfg))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(Self {
            inner: RwLock::new(Arc::new(eng)),
        })
    }

    /// Loads an `Engine` from explicit `.bsp` file paths, releasing the GIL.
    #[staticmethod]
    #[pyo3(signature = (paths, max_threads=None, pin_time_table_version=None))]
    pub fn from_bsp_paths(
        py: Python<'_>,
        paths: Vec<String>,
        max_threads: Option<usize>,
        pin_time_table_version: Option<String>,
    ) -> PyResult<Self> {
        let cfg = EngineConfig {
            max_threads,
            pin_time_table_version,
            verify_manifest_sha256: false,
        };
        let eng = py
            .allow_threads(|| Engine::from_bsp_paths(&paths, cfg))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(Self {
            inner: RwLock::new(Arc::new(eng)),
        })
    }

    /// Loads an `Engine` from in-memory `(name, bytes)` tuples, releasing the GIL.
    #[staticmethod]
    #[pyo3(signature = (kernels, max_threads=None))]
    pub fn from_bytes(
        py: Python<'_>,
        kernels: Vec<(String, Vec<u8>)>,
        max_threads: Option<usize>,
    ) -> PyResult<Self> {
        let cfg = EngineConfig {
            max_threads,
            pin_time_table_version: None,
            verify_manifest_sha256: false,
        };
        let eng = py
            .allow_threads(|| Engine::from_bytes(&kernels, cfg))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(Self {
            inner: RwLock::new(Arc::new(eng)),
        })
    }

    /// Creates a fileless analytical Moshier engine without external kernels.
    #[staticmethod]
    #[pyo3(signature = (max_threads=None, pin_time_table_version=None))]
    pub fn new_moshier(
        max_threads: Option<usize>,
        pin_time_table_version: Option<String>,
    ) -> Self {
        let cfg = EngineConfig {
            max_threads,
            pin_time_table_version,
            verify_manifest_sha256: false,
        };
        let eng = Engine::new_moshier(cfg);
        Self {
            inner: RwLock::new(Arc::new(eng)),
        }
    }

    /// Loads an `Engine` from `data_dir` if kernels exist, or falls back to fileless Moshier mode.
    #[staticmethod]
    #[pyo3(signature = (data_dir, max_threads=None, pin_time_table_version=None, verify_manifest_sha256=false))]
    pub fn from_data_dir_or_moshier(
        py: Python<'_>,
        data_dir: String,
        max_threads: Option<usize>,
        pin_time_table_version: Option<String>,
        verify_manifest_sha256: bool,
    ) -> PyResult<Self> {
        let cfg = EngineConfig {
            max_threads,
            pin_time_table_version,
            verify_manifest_sha256,
        };
        let eng = py
            .allow_threads(|| Engine::from_data_dir_or_moshier(PathBuf::from(data_dir), cfg))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(Self {
            inner: RwLock::new(Arc::new(eng)),
        })
    }

    /// Returns `true` if this engine is operating in pure fileless analytical ("Moshier") mode.
    pub fn is_moshier_mode(&self) -> bool {
        self.inner.read().is_ok_and(|g| g.is_moshier_mode())
    }

    /// Returns the active license status (`AttrDict` with `status`, `tier`, `is_commercial_valid`,
    /// `license_id`, `holder`, `updates_until`, `engine_release_date`, `message`).
    pub fn license_info(&self, py: Python<'_>) -> PyResult<PyObject> {
        let info = self.get_engine()?.license_info();
        serialize_to_attr_dict(py, &info)
    }

    /// Loads additional custom stars from a `sefstars.txt` or JSON file into this `Engine`.
    pub fn load_custom_stars(&self, path: String) -> PyResult<usize> {
        let mut guard = self
            .inner
            .write()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        if let Some(eng) = Arc::get_mut(&mut *guard) {
            eng.load_custom_stars(path).map_err(ephemeris_err_to_pyerr)
        } else {
            Err(PyRuntimeError::new_err(
                "Cannot load custom stars: Engine instance is currently shared with an active adapter",
            ))
        }
    }

    /// Evaluates `body` at `(jd1, jd2)` explicitly using fileless analytical Moshier formulas.
    #[pyo3(signature = (body, jd1, jd2=0.0, scale="TT", options=None))]
    pub fn calc_moshier(
        &self,
        py: Python<'_>,
        body: &str,
        jd1: f64,
        jd2: f64,
        scale: &str,
        options: Option<PyCalcOptions>,
    ) -> PyResult<PyCalcResult> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.calc_moshier(bid, t_in, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyCalcResult::from(res))
    }

    /// Evaluates planetary, solar, lunar, or minor-body phenomena (`swe_pheno` / `swe_pheno_ut` compatible).
    #[pyo3(signature = (body, jd1, jd2=0.0, scale="TT", options=None))]
    pub fn phenomena(
        &self,
        py: Python<'_>,
        body: &str,
        jd1: f64,
        jd2: f64,
        scale: &str,
        options: Option<PyCalcOptions>,
    ) -> PyResult<PyPhenomenaResult> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.phenomena(bid, t_in, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyPhenomenaResult::from(res))
    }

    /// Evaluates planetary or lunar orbital nodes and apsides (`swe_nod_aps`).
    #[pyo3(signature = (body, jd1, jd2=0.0, scale="TT", method=0, options=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn nod_aps(
        &self,
        py: Python<'_>,
        body: &str,
        jd1: f64,
        jd2: f64,
        scale: &str,
        method: u32,
        options: Option<PyCalcOptions>,
    ) -> PyResult<PyNodApsResult> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.nod_aps(bid, t_in, opts, method))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyNodApsResult::from(res))
    }

    /// Searches forward or backward for the next global or local lunar occultation of `body`.
    #[pyo3(signature = (body, jd1, jd2=0.0, scale="TT", backward=false, observer=None, max_synodic_months=36))]
    #[allow(clippy::too_many_arguments)]
    pub fn next_lunar_occultation(
        &self,
        py: Python<'_>,
        body: &str,
        jd1: f64,
        jd2: f64,
        scale: &str,
        backward: bool,
        observer: Option<PyObserver>,
        max_synodic_months: u32,
    ) -> PyResult<PyLunarOccultationInfo> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let max_search_days = (max_synodic_months as f64) * 29.530_588_853;
        let engine = self.get_engine()?;
        if let Some(obs) = observer {
            let core_obs = Observer {
                latitude_deg: obs.latitude_deg,
                longitude_deg: obs.longitude_deg,
                altitude_m: obs.altitude_m,
            };
            let info = py
                .allow_threads(move || {
                    engine.next_local_lunar_occultation(bid, core_obs, t_in, backward, max_search_days)
                })
                .map_err(ephemeris_err_to_pyerr)?;
            Ok(PyLunarOccultationInfo::from(info))
        } else {
            let info = py
                .allow_threads(move || {
                    engine.next_lunar_occultation(bid, t_in, backward, max_search_days)
                })
                .map_err(ephemeris_err_to_pyerr)?;
            Ok(PyLunarOccultationInfo::from(info))
        }
    }

    /// Searches forward or backward for the next locally visible lunar occultation of `body`.
    #[pyo3(signature = (body, jd1, latitude_deg, longitude_deg, jd2=0.0, scale="TT", backward=false, altitude_m=0.0, max_synodic_months=120))]
    #[allow(clippy::too_many_arguments)]
    pub fn next_local_lunar_occultation(
        &self,
        py: Python<'_>,
        body: &str,
        jd1: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        jd2: f64,
        scale: &str,
        backward: bool,
        altitude_m: f64,
        max_synodic_months: u32,
    ) -> PyResult<PyLunarOccultationInfo> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let max_search_days = (max_synodic_months as f64) * 29.530_588_853;
        let engine = self.get_engine()?;
        let info = py
            .allow_threads(move || {
                engine.next_local_lunar_occultation(bid, obs, t_in, backward, max_search_days)
            })
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyLunarOccultationInfo::from(info))
    }

    /// Finds all lunar occultations of `body` within `[start_jd_tt, end_jd_tt]`.
    #[pyo3(signature = (body, start_jd_tt, end_jd_tt, observer=None))]
    pub fn lunar_occultations_in_range(
        &self,
        py: Python<'_>,
        body: &str,
        start_jd_tt: f64,
        end_jd_tt: f64,
        observer: Option<PyObserver>,
    ) -> PyResult<Vec<PyLunarOccultationInfo>> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let start_in = TimeInput::from_tt_jd(start_jd_tt);
        let end_in = TimeInput::from_tt_jd(end_jd_tt);
        if let Some(obs) = observer {
            let core_obs = Observer {
                latitude_deg: obs.latitude_deg,
                longitude_deg: obs.longitude_deg,
                altitude_m: obs.altitude_m,
            };
            let list = py
                .allow_threads(move || {
                    engine.local_lunar_occultations_in_range(bid, core_obs, start_in, end_in)
                })
                .map_err(ephemeris_err_to_pyerr)?;
            Ok(list.into_iter().map(PyLunarOccultationInfo::from).collect())
        } else {
            let list = py
                .allow_threads(move || {
                    engine.lunar_occultations_in_range(bid, start_in, end_in)
                })
                .map_err(ephemeris_err_to_pyerr)?;
            Ok(list.into_iter().map(PyLunarOccultationInfo::from).collect())
        }
    }

    /// Finds all locally visible lunar occultations of `body` within `[start_jd_tt, end_jd_tt]`.
    #[pyo3(signature = (body, start_jd_tt, end_jd_tt, latitude_deg, longitude_deg, altitude_m=0.0))]
    #[allow(clippy::too_many_arguments)]
    pub fn local_lunar_occultations_in_range(
        &self,
        py: Python<'_>,
        body: &str,
        start_jd_tt: f64,
        end_jd_tt: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        altitude_m: f64,
    ) -> PyResult<Vec<PyLunarOccultationInfo>> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let start_in = TimeInput::from_tt_jd(start_jd_tt);
        let end_in = TimeInput::from_tt_jd(end_jd_tt);
        let list = py
            .allow_threads(move || {
                engine.local_lunar_occultations_in_range(bid, obs, start_in, end_in)
            })
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(list.into_iter().map(PyLunarOccultationInfo::from).collect())
    }

    /// Computes the 3D fractional mundane house position $H \in [1.0, 13.0)$ of a sky point (`swe_house_pos`).
    #[pyo3(signature = (armc_deg, geolat_deg, eps_deg, lon_deg, lat_deg=0.0, hsys=None))]
    pub fn house_pos(
        &self,
        armc_deg: f64,
        geolat_deg: f64,
        eps_deg: f64,
        lon_deg: f64,
        lat_deg: f64,
        hsys: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<f64> {
        let ch = extract_hsys_char(hsys)?;
        let system = HouseSystem::from_code(ch).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        engine
            .house_pos(armc_deg, geolat_deg, eps_deg, system, lon_deg, lat_deg)
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates the Equation of Time in days at `(jd1, jd2)` (`swe_time_equ`).
    #[pyo3(signature = (jd1, jd2=0.0, scale="UT1"))]
    pub fn time_equation(&self, jd1: f64, jd2: f64, scale: &str) -> PyResult<f64> {
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        engine.time_equation(t_in).map_err(ephemeris_err_to_pyerr)
    }

    /// Inverse horizontal transformation (`swe_azalt_rev`): geometric horizontal to Equatorial or Ecliptic.
    #[pyo3(signature = (jd1, calc_flag, latitude_deg, longitude_deg, az_deg, true_alt_deg, jd2=0.0, scale="UT1", altitude_m=0.0))]
    #[allow(clippy::too_many_arguments)]
    pub fn azalt_rev(
        &self,
        jd1: f64,
        calc_flag: i32,
        latitude_deg: f64,
        longitude_deg: f64,
        az_deg: f64,
        true_alt_deg: f64,
        jd2: f64,
        scale: &str,
        altitude_m: f64,
    ) -> PyResult<[f64; 2]> {
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let plane = if calc_flag == 0 {
            Plane::Ecliptic
        } else {
            Plane::Equatorial
        };
        let engine = self.get_engine()?;
        Ok(engine.azalt_rev(obs, t_in, az_deg, true_alt_deg, plane))
    }

    /// Returns capabilities metadata as a JSON string.
    pub fn capabilities_json(&self) -> PyResult<String> {
        let engine = self.get_engine()?;
        serde_json::to_string(&engine.capabilities())
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    /// Returns coverage metadata for `body` as a JSON string.
    pub fn coverage_json(&self, body: &str) -> PyResult<String> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let cov = engine.coverage(bid).map_err(ephemeris_err_to_pyerr)?;
        serde_json::to_string(&cov).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    /// Evaluates $\Delta T = \text{TT} - \text{UT1}$ and $\text{TT} - \text{TDB}$ at `(jd1, jd2)` in `scale` (`R1.2`).
    #[pyo3(signature = (jd1, jd2=0.0, scale="UT1"))]
    pub fn delta_t(&self, jd1: f64, jd2: f64, scale: &str) -> PyResult<PyDeltaTInfo> {
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        Ok(PyDeltaTInfo::from(engine.delta_t(t_in)))
    }

    /// Evaluates GMST, GAST, Local Sidereal Time, and ARMC at `(jd1, jd2)` in `scale` (`R1.4`).
    #[pyo3(signature = (jd1, jd2=0.0, scale="UT1", nutation="true", longitude_deg=0.0))]
    pub fn sidtime(
        &self,
        jd1: f64,
        jd2: f64,
        scale: &str,
        nutation: &str,
        longitude_deg: f64,
    ) -> PyResult<PySiderealTimeInfo> {
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let nut = parse_nutation_mode(nutation).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let info = engine.sidtime(t_in, nut, longitude_deg);
        Ok(PySiderealTimeInfo::from(info))
    }

    /// Evaluates mean/true obliquity and nutation in longitude/obliquity (`ECL_NUT`, `R2.2`).
    #[pyo3(signature = (jd1, jd2=0.0, scale="TT"))]
    pub fn ecl_nut(&self, jd1: f64, jd2: f64, scale: &str) -> PyResult<PyEclNutInfo> {
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let info = engine.ecl_nut(t_in).map_err(ephemeris_err_to_pyerr)?;
        Ok(PyEclNutInfo::from(info))
    }

    /// Looks up catalog metadata (including `vmag`) for a fixed star (`R6.1`).
    pub fn star_info(&self, name_or_bayer: &str) -> PyResult<PyStarMetadata> {
        let engine = self.get_engine()?;
        let info = engine
            .star_info(name_or_bayer)
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyStarMetadata::from(info))
    }

    /// Lists all 52 cataloged fixed stars with metadata (`R6.1`).
    pub fn list_stars(&self) -> PyResult<Vec<PyStarMetadata>> {
        let engine = self.get_engine()?;
        Ok(engine
            .list_stars()
            .into_iter()
            .map(PyStarMetadata::from)
            .collect())
    }

    /// Evaluates a single body at `(jd1, jd2)` in `scale`, releasing the Python GIL.
    #[pyo3(signature = (body, jd1, jd2=0.0, scale="TT", options=None))]
    pub fn calc(
        &self,
        py: Python<'_>,
        body: &str,
        jd1: f64,
        jd2: f64,
        scale: &str,
        options: Option<PyCalcOptions>,
    ) -> PyResult<PyCalcResult> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.calc(bid, t_in, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyCalcResult::from(res))
    }

    /// Evaluates a fixed star by canonical name, Bayer designation, or `"hip:<id>"`.
    #[pyo3(signature = (name_or_bayer, jd1, jd2=0.0, scale="TT", options=None))]
    pub fn calc_star_by_name(
        &self,
        py: Python<'_>,
        name_or_bayer: &str,
        jd1: f64,
        jd2: f64,
        scale: &str,
        options: Option<PyCalcOptions>,
    ) -> PyResult<PyCalcResult> {
        let q = name_or_bayer.to_owned();
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.calc_star_by_name(&q, t_in, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyCalcResult::from(res))
    }

    /// Evaluates a batch of `bodies` across `epochs_jd` in parallel (atomic mode), releasing the GIL.
    #[pyo3(signature = (bodies, epochs_jd, scale="TT", options=None))]
    pub fn calc_batch(
        &self,
        py: Python<'_>,
        bodies: Vec<String>,
        epochs_jd: Vec<f64>,
        scale: &str,
        options: Option<PyCalcOptions>,
    ) -> PyResult<Vec<Vec<PyCalcResult>>> {
        let bids: Vec<BodyId> = bodies
            .iter()
            .map(|b| parse_body_id(b))
            .collect::<Result<Vec<_>, _>>()
            .map_err(ephemeris_err_to_pyerr)?;
        let t_inputs: Vec<TimeInput> = epochs_jd
            .iter()
            .map(|&jd| parse_time_input(scale, jd, 0.0))
            .collect::<Result<Vec<_>, _>>()
            .map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let batch = py
            .allow_threads(move || engine.calc_batch(&bids, &t_inputs, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(batch
            .into_iter()
            .map(|row| row.into_iter().map(PyCalcResult::from).collect())
            .collect())
    }

    /// Evaluates a batch of `bodies` across `epochs_jd` in parallel, returning a flat contiguous
    /// `Vec<f64>` of length `len(epochs_jd) * len(bodies) * 6` containing
    /// `[lon_deg, lat_deg, dist_au, lon_speed_deg_day, lat_speed_deg_day, dist_speed_au_day]`.
    #[pyo3(signature = (bodies, epochs_jd, scale="TT", options=None))]
    pub fn calc_batch_flat_f64(
        &self,
        py: Python<'_>,
        bodies: Vec<String>,
        epochs_jd: Vec<f64>,
        scale: &str,
        options: Option<PyCalcOptions>,
    ) -> PyResult<Vec<f64>> {
        let bids: Vec<BodyId> = bodies
            .iter()
            .map(|b| parse_body_id(b))
            .collect::<Result<Vec<_>, _>>()
            .map_err(ephemeris_err_to_pyerr)?;
        let t_inputs: Vec<TimeInput> = epochs_jd
            .iter()
            .map(|&jd| parse_time_input(scale, jd, 0.0))
            .collect::<Result<Vec<_>, _>>()
            .map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let flat = py
            .allow_threads(move || -> Result<Vec<f64>, EphemerisError> {
                let batch = engine.calc_batch(&bids, &t_inputs, opts)?;
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
            })
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(flat)
    }

    /// Evaluates a batch of `bodies` across `epochs_jd` in per-item mode without aborting on individual epoch errors.
    #[pyo3(signature = (bodies, epochs_jd, scale="TT", options=None))]
    pub fn calc_batch_per_item(
        &self,
        py: Python<'_>,
        bodies: Vec<String>,
        epochs_jd: Vec<f64>,
        scale: &str,
        options: Option<PyCalcOptions>,
    ) -> PyResult<Vec<Vec<PyBatchItemResult>>> {
        let bids: Vec<BodyId> = bodies
            .iter()
            .map(|b| parse_body_id(b))
            .collect::<Result<Vec<_>, _>>()
            .map_err(ephemeris_err_to_pyerr)?;
        let t_inputs: Vec<TimeInput> = epochs_jd
            .iter()
            .map(|&jd| parse_time_input(scale, jd, 0.0))
            .collect::<Result<Vec<_>, _>>()
            .map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let batch = py
            .allow_threads(move || engine.calc_batch_per_item(&bids, &t_inputs, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(batch
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|item| match item {
                        Ok(r) => PyBatchItemResult {
                            ok: true,
                            result: Some(PyCalcResult::from(r)),
                            error_code: None,
                            error_message: None,
                        },
                        Err(e) => PyBatchItemResult {
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

    /// Evaluates any of the 12 Vedic / Sidereal ayanamsa models at `(jd1, jd2)` in `scale`.
    #[pyo3(signature = (model="lahiri", jd1=2451545.0, jd2=0.0, scale="TT", nutation="true"))]
    pub fn ayanamsa(
        &self,
        py: Python<'_>,
        model: &str,
        jd1: f64,
        jd2: f64,
        scale: &str,
        nutation: &str,
    ) -> PyResult<PyAyanamsaValue> {
        let m = AyanamsaModel::from_code(model).map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let nut = parse_nutation_mode(nutation).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let val = py
            .allow_threads(move || engine.ayanamsa(t_in, m, nut))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyAyanamsaValue::from(val))
    }

    /// Evaluates any of the 23 house systems at `(jd1, jd2)` in `scale` for `(latitude_deg, longitude_deg)`.
    #[pyo3(signature = (jd1, latitude_deg, longitude_deg, system=None, jd2=0.0, scale="TT", altitude_m=0.0, polar_policy="strict", zodiac="tropical", ayanamsa="lahiri", nutation="true", include_speed=true))]
    #[allow(clippy::too_many_arguments)]
    pub fn houses(
        &self,
        py: Python<'_>,
        jd1: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        system: Option<&Bound<'_, PyAny>>,
        jd2: f64,
        scale: &str,
        altitude_m: f64,
        polar_policy: &str,
        zodiac: &str,
        ayanamsa: &str,
        nutation: &str,
        include_speed: bool,
    ) -> PyResult<PyHouseResult> {
        let ch = extract_hsys_char(system)?;
        let sys = HouseSystem::from_code(ch).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let pol = parse_polar_policy(polar_policy).map_err(ephemeris_err_to_pyerr)?;
        let zod = parse_zodiac(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let nut = parse_nutation_mode(nutation).map_err(ephemeris_err_to_pyerr)?;
        let opts = HouseOptions {
            system: sys,
            polar_policy: pol,
            zodiac: zod,
            nutation: nut,
            include_speed,
        };
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.houses(t_in, obs, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyHouseResult::from(res))
    }

    /// Evaluates house cusps and chart angles across a batch of `epochs_jd` in parallel (`R9.3`), releasing the GIL.
    #[pyo3(signature = (epochs_jd, latitude_deg, longitude_deg, system=None, scale="TT", altitude_m=0.0, polar_policy="strict", zodiac="tropical", ayanamsa="lahiri", nutation="true", include_speed=true))]
    #[allow(clippy::too_many_arguments)]
    pub fn houses_batch(
        &self,
        py: Python<'_>,
        epochs_jd: Vec<f64>,
        latitude_deg: f64,
        longitude_deg: f64,
        system: Option<&Bound<'_, PyAny>>,
        scale: &str,
        altitude_m: f64,
        polar_policy: &str,
        zodiac: &str,
        ayanamsa: &str,
        nutation: &str,
        include_speed: bool,
    ) -> PyResult<Vec<PyHouseResult>> {
        let ch = extract_hsys_char(system)?;
        let sys = HouseSystem::from_code(ch).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let t_inputs: Vec<TimeInput> = epochs_jd
            .iter()
            .map(|&jd| parse_time_input(scale, jd, 0.0))
            .collect::<Result<Vec<_>, _>>()
            .map_err(ephemeris_err_to_pyerr)?;
        let pol = parse_polar_policy(polar_policy).map_err(ephemeris_err_to_pyerr)?;
        let zod = parse_zodiac(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let nut = parse_nutation_mode(nutation).map_err(ephemeris_err_to_pyerr)?;
        let opts = HouseOptions {
            system: sys,
            polar_policy: pol,
            zodiac: zod,
            nutation: nut,
            include_speed,
        };
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.houses_batch(&t_inputs, obs, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(res.into_iter().map(PyHouseResult::from).collect())
    }

    /// Evaluates topocentric horizontal coordinates (`azalt`) with optical refraction from `(ra_deg, dec_deg)`.
    #[pyo3(signature = (ra_deg, dec_deg, jd1, latitude_deg, longitude_deg, jd2=0.0, scale="TT", altitude_m=0.0, pressure_hpa=1013.25, temperature_c=15.0))]
    #[allow(clippy::too_many_arguments)]
    pub fn azalt(
        &self,
        ra_deg: f64,
        dec_deg: f64,
        jd1: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        jd2: f64,
        scale: &str,
        altitude_m: f64,
        pressure_hpa: f64,
        temperature_c: f64,
    ) -> PyResult<PyHorizontalCoords> {
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let moment = engine
            .at(t_in, Some(obs))
            .map_err(ephemeris_err_to_pyerr)?;
        let atmo = AtmosphereParams {
            pressure_hpa,
            temperature_c,
            ..Default::default()
        };
        let hz = core_azalt(obs, moment.center().gast_rad, ra_deg, dec_deg, atmo);
        Ok(PyHorizontalCoords::from(hz))
    }

    /// Evaluates topocentric horizontal coordinates (`azalt`) from ecliptic `(ecl_lon_deg, ecl_lat_deg)` (`R2.6`).
    #[pyo3(signature = (ecl_lon_deg, ecl_lat_deg, jd1, latitude_deg, longitude_deg, jd2=0.0, scale="TT", altitude_m=0.0, pressure_hpa=1013.25, temperature_c=15.0))]
    #[allow(clippy::too_many_arguments)]
    pub fn azalt_ecliptic(
        &self,
        ecl_lon_deg: f64,
        ecl_lat_deg: f64,
        jd1: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        jd2: f64,
        scale: &str,
        altitude_m: f64,
        pressure_hpa: f64,
        temperature_c: f64,
    ) -> PyResult<PyHorizontalCoords> {
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let atmo = AtmosphereParams {
            pressure_hpa,
            temperature_c,
            ..Default::default()
        };
        let engine = self.get_engine()?;
        let hz = engine
            .azalt_ecliptic(obs, t_in, ecl_lon_deg, ecl_lat_deg, atmo);
        Ok(PyHorizontalCoords::from(hz))
    }

    /// Evaluates topocentric horizontal coordinates (`azalt`) directly for `body` (`R2.6`).
    #[pyo3(signature = (body, jd1, latitude_deg, longitude_deg, jd2=0.0, scale="TT", altitude_m=0.0, pressure_hpa=1013.25, temperature_c=15.0))]
    #[allow(clippy::too_many_arguments)]
    pub fn calc_azalt(
        &self,
        py: Python<'_>,
        body: &str,
        jd1: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        jd2: f64,
        scale: &str,
        altitude_m: f64,
        pressure_hpa: f64,
        temperature_c: f64,
    ) -> PyResult<PyHorizontalCoords> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let atmo = AtmosphereParams {
            pressure_hpa,
            temperature_c,
            ..Default::default()
        };
        let engine = self.get_engine()?;
        let hz = py
            .allow_threads(move || engine.calc_azalt(bid, t_in, obs, atmo))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyHorizontalCoords::from(hz))
    }

    /// Searches for the next rise, set, upper meridian transit, or lower meridian transit (`JD_UT1`).
    #[pyo3(signature = (body, event, jd1, latitude_deg, longitude_deg, jd2=0.0, scale="UT1", altitude_m=0.0, custom_horizon_alt_deg=None, max_search_days=2.0, disc_mode="upper_limb", pressure_hpa=1013.25, temperature_c=15.0))]
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans(
        &self,
        py: Python<'_>,
        body: &str,
        event: &str,
        jd1: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        jd2: f64,
        scale: &str,
        altitude_m: f64,
        custom_horizon_alt_deg: Option<f64>,
        max_search_days: f64,
        disc_mode: &str,
        pressure_hpa: f64,
        temperature_c: f64,
    ) -> PyResult<Option<f64>> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let moment = engine
            .at(t_in, Some(obs))
            .map_err(ephemeris_err_to_pyerr)?;
        let jd_ut1_start = moment.center().ut1.jd_f64();
        let ev = match event.trim().to_ascii_lowercase().as_str() {
            "rise" => RiseTransEvent::Rise,
            "set" => RiseTransEvent::Set,
            "upper" | "upper_transit" | "culmination" => RiseTransEvent::UpperTransit,
            "lower" | "lower_transit" | "anti_culmination" => RiseTransEvent::LowerTransit,
            other => {
                return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: format!("Unknown transit event '{other}'"),
                }))
            }
        };
        let dm = parse_disc_mode(disc_mode).map_err(ephemeris_err_to_pyerr)?;
        let opts = RiseTransOptions {
            max_search_days,
            custom_horizon_alt_deg,
            disc_mode: dm,
            pressure_hpa,
            temperature_c,
        };
        let res = py
            .allow_threads(move || engine.rise_trans(bid, obs, jd_ut1_start, ev, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(res)
    }

    /// Computes heliacal risings and settings via Bradley Schaefer's visibility model (`R6.3`).
    #[pyo3(signature = (body, jd1, latitude_deg, longitude_deg, jd2=0.0, scale="UT1", altitude_m=0.0, extinction_k_v=0.20, observer_age_years=36.0, snellen_factor=1.0, snellen_ratio=None, pressure_hpa=1013.25, temperature_c=15.0, relative_humidity=0.50, max_search_days=400.0))]
    #[allow(clippy::too_many_arguments)]
    pub fn heliacal_phenomena(
        &self,
        py: Python<'_>,
        body: &str,
        jd1: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        jd2: f64,
        scale: &str,
        altitude_m: f64,
        extinction_k_v: f64,
        observer_age_years: f64,
        snellen_factor: f64,
        snellen_ratio: Option<f64>,
        pressure_hpa: f64,
        temperature_c: f64,
        relative_humidity: f64,
        max_search_days: f64,
    ) -> PyResult<PyHeliacalResult> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let moment = engine
            .at(t_in, Some(obs))
            .map_err(ephemeris_err_to_pyerr)?;
        let jd_ut1_start = moment.center().ut1.jd_f64();
        let opts = HeliacalOptions {
            k_v: extinction_k_v,
            observer_age_years,
            snellen_ratio: snellen_ratio.unwrap_or(snellen_factor),
            pressure_hpa,
            temperature_c,
            relative_humidity,
            max_search_days,
        };
        let res = py
            .allow_threads(move || engine.heliacal_phenomena(bid, jd_ut1_start, obs, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyHeliacalResult::from(res))
    }

    /// Searches forward or backward for a solar eclipse from `(jd1, jd2)` in `scale`.
    #[pyo3(signature = (jd1, jd2=0.0, scale="TT", backward=false, observer=None, max_synodic_months=36))]
    #[allow(clippy::too_many_arguments)]
    pub fn next_solar_eclipse(
        &self,
        py: Python<'_>,
        jd1: f64,
        jd2: f64,
        scale: &str,
        backward: bool,
        observer: Option<PyObserver>,
        max_synodic_months: u32,
    ) -> PyResult<PySolarEclipseInfo> {
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let obs = observer.map(|o| Observer {
            latitude_deg: o.latitude_deg,
            longitude_deg: o.longitude_deg,
            altitude_m: o.altitude_m,
        });
        let opts = EclipseOptions {
            max_synodic_months,
            ..Default::default()
        };
        let engine = self.get_engine()?;
        let info = py
            .allow_threads(move || engine.next_solar_eclipse(t_in, backward, obs, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PySolarEclipseInfo::from(info))
    }

    /// Searches backward for the previous solar eclipse before `(jd1, jd2)`.
    #[pyo3(signature = (jd1, jd2=0.0, scale="TT", observer=None, max_synodic_months=36))]
    pub fn prev_solar_eclipse(
        &self,
        py: Python<'_>,
        jd1: f64,
        jd2: f64,
        scale: &str,
        observer: Option<PyObserver>,
        max_synodic_months: u32,
    ) -> PyResult<PySolarEclipseInfo> {
        self.next_solar_eclipse(py, jd1, jd2, scale, true, observer, max_synodic_months)
    }

    /// Searches forward or backward for the next solar eclipse that is locally visible (`R8.1`).
    #[pyo3(signature = (jd1, latitude_deg, longitude_deg, jd2=0.0, scale="TT", backward=false, altitude_m=0.0, max_synodic_months=120))]
    #[allow(clippy::too_many_arguments)]
    pub fn next_local_solar_eclipse(
        &self,
        py: Python<'_>,
        jd1: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        jd2: f64,
        scale: &str,
        backward: bool,
        altitude_m: f64,
        max_synodic_months: u32,
    ) -> PyResult<PySolarEclipseInfo> {
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(latitude_deg, longitude_deg, altitude_m)
            .map_err(ephemeris_err_to_pyerr)?;
        let opts = EclipseOptions {
            max_synodic_months,
            ..Default::default()
        };
        let engine = self.get_engine()?;
        let info = py
            .allow_threads(move || engine.next_local_solar_eclipse(t_in, obs, backward, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PySolarEclipseInfo::from(info))
    }

    /// Enumerates all solar eclipses in `[start_jd_tt, end_jd_tt]` (`R8.3`).
    #[pyo3(signature = (start_jd_tt, end_jd_tt, observer=None))]
    pub fn solar_eclipses_in_range(
        &self,
        py: Python<'_>,
        start_jd_tt: f64,
        end_jd_tt: f64,
        observer: Option<PyObserver>,
    ) -> PyResult<Vec<PySolarEclipseInfo>> {
        let obs = observer.map(|o| Observer {
            latitude_deg: o.latitude_deg,
            longitude_deg: o.longitude_deg,
            altitude_m: o.altitude_m,
        });
        let engine = self.get_engine()?;
        let list = py
            .allow_threads(move || {
                engine.solar_eclipses_in_range(
                    TimeInput::from_tt_jd(start_jd_tt),
                    TimeInput::from_tt_jd(end_jd_tt),
                    obs,
                    EclipseOptions::default(),
                )
            })
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(list.into_iter().map(PySolarEclipseInfo::from).collect())
    }

    /// Searches forward or backward for a lunar eclipse from `(jd1, jd2)` in `scale`.
    #[pyo3(signature = (jd1, jd2=0.0, scale="TT", backward=false, shadow_rule="danjon", max_synodic_months=36))]
    #[allow(clippy::too_many_arguments)]
    pub fn next_lunar_eclipse(
        &self,
        py: Python<'_>,
        jd1: f64,
        jd2: f64,
        scale: &str,
        backward: bool,
        shadow_rule: &str,
        max_synodic_months: u32,
    ) -> PyResult<PyLunarEclipseInfo> {
        let t_in = parse_time_input(scale, jd1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let rule = parse_shadow_rule(shadow_rule).map_err(ephemeris_err_to_pyerr)?;
        let opts = EclipseOptions {
            shadow_rule: rule,
            max_synodic_months,
            ..Default::default()
        };
        let engine = self.get_engine()?;
        let info = py
            .allow_threads(move || engine.next_lunar_eclipse(t_in, backward, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(PyLunarEclipseInfo::from(info))
    }

    /// Searches backward for the previous lunar eclipse before `(jd1, jd2)`.
    #[pyo3(signature = (jd1, jd2=0.0, scale="TT", shadow_rule="danjon", max_synodic_months=36))]
    pub fn prev_lunar_eclipse(
        &self,
        py: Python<'_>,
        jd1: f64,
        jd2: f64,
        scale: &str,
        shadow_rule: &str,
        max_synodic_months: u32,
    ) -> PyResult<PyLunarEclipseInfo> {
        self.next_lunar_eclipse(py, jd1, jd2, scale, true, shadow_rule, max_synodic_months)
    }

    /// Enumerates all lunar eclipses in `[start_jd_tt, end_jd_tt]` (`R8.3`).
    #[pyo3(signature = (start_jd_tt, end_jd_tt, shadow_rule="danjon"))]
    pub fn lunar_eclipses_in_range(
        &self,
        py: Python<'_>,
        start_jd_tt: f64,
        end_jd_tt: f64,
        shadow_rule: &str,
    ) -> PyResult<Vec<PyLunarEclipseInfo>> {
        let rule = parse_shadow_rule(shadow_rule).map_err(ephemeris_err_to_pyerr)?;
        let opts = EclipseOptions {
            shadow_rule: rule,
            ..Default::default()
        };
        let engine = self.get_engine()?;
        let list = py
            .allow_threads(move || {
                engine.lunar_eclipses_in_range(
                    TimeInput::from_tt_jd(start_jd_tt),
                    TimeInput::from_tt_jd(end_jd_tt),
                    opts,
                )
            })
            .map_err(ephemeris_err_to_pyerr)?;
        Ok(list.into_iter().map(PyLunarEclipseInfo::from).collect())
    }

    /// Finds the first instant starting from `jd_ut1` (moving forward if `forward == true`,
    /// or backward if `forward == false`, up to `max_search_days`) when `body`'s ecliptic longitude
    /// reaches `target_lon_deg` (`R9.1`, `R9.2`).
    #[pyo3(signature = (body, jd_ut1, target_lon_deg, forward=true, max_search_days=400.0, options=None, scale="UT1", jd2=0.0))]
    #[allow(clippy::too_many_arguments)]
    pub fn search_body_longitude(
        &self,
        py: Python<'_>,
        body: &str,
        jd_ut1: f64,
        target_lon_deg: f64,
        forward: bool,
        max_search_days: f64,
        options: Option<&PyCalcOptions>,
        scale: &str,
        jd2: f64,
    ) -> PyResult<PyObject> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let start_in = parse_time_input(scale, jd_ut1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let ev = py
            .allow_threads(move || {
                engine.search_body_longitude(
                    bid,
                    start_in,
                    target_lon_deg,
                    forward,
                    max_search_days,
                    &opts,
                )
            })
            .map_err(ephemeris_err_to_pyerr)?;
        longitude_crossing_to_pydict(py, &ev)
    }

    /// Finds all crossings of `target_lon_deg` by `body` in `[start_jd_ut1, end_jd_ut1]`
    /// in chronological order, including 3-pass retrograde triple returns (`R9.1`, `R9.2`).
    #[pyo3(signature = (body, start_jd_ut1, end_jd_ut1, target_lon_deg, options=None, scale="UT1"))]
    #[allow(clippy::too_many_arguments)]
    pub fn body_longitude_crossings_in_range(
        &self,
        py: Python<'_>,
        body: &str,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        target_lon_deg: f64,
        options: Option<&PyCalcOptions>,
        scale: &str,
    ) -> PyResult<PyObject> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let start_in = parse_time_input(scale, start_jd_ut1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let end_in = parse_time_input(scale, end_jd_ut1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let events = py
            .allow_threads(move || {
                engine.body_longitude_crossings_in_range(
                    bid,
                    start_in,
                    end_in,
                    target_lon_deg,
                    &opts,
                )
            })
            .map_err(ephemeris_err_to_pyerr)?;
        let list = PyList::empty(py);
        for ev in &events {
            list.append(longitude_crossing_to_pydict(py, ev)?)?;
        }
        Ok(list.into())
    }

    /// Alias for [`Self::body_longitude_crossings_in_range`].
    #[pyo3(signature = (body, start_jd_ut1, end_jd_ut1, target_lon_deg, options=None, scale="UT1"))]
    #[allow(clippy::too_many_arguments)]
    pub fn longitude_crossings_in_range(
        &self,
        py: Python<'_>,
        body: &str,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        target_lon_deg: f64,
        options: Option<&PyCalcOptions>,
        scale: &str,
    ) -> PyResult<PyObject> {
        self.body_longitude_crossings_in_range(
            py,
            body,
            start_jd_ut1,
            end_jd_ut1,
            target_lon_deg,
            options,
            scale,
        )
    }

    /// Finds all Chinese 24 Solar Terms (*Jieqi*) in `[start_jd_ut1, end_jd_ut1]` (`R12`).
    #[pyo3(signature = (start_jd_ut1, end_jd_ut1, scale="UT1"))]
    pub fn solar_terms_in_range(
        &self,
        py: Python<'_>,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        scale: &str,
    ) -> PyResult<PyObject> {
        let start_in = parse_time_input(scale, start_jd_ut1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let end_in = parse_time_input(scale, end_jd_ut1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let terms = py
            .allow_threads(move || engine.solar_terms_in_range(start_in, end_in))
            .map_err(ephemeris_err_to_pyerr)?;
        let list = PyList::empty(py);
        for term in &terms {
            list.append(solar_term_to_pydict(py, term)?)?;
        }
        Ok(list.into())
    }

    /// Computes all 4 cardinal seasons (`0°`, `90°`, `180°`, `270°` apparent solar longitude) for `year`.
    pub fn seasons_in_year(&self, py: Python<'_>, year: i32) -> PyResult<PyObject> {
        let engine = self.get_engine()?;
        let seasons = py
            .allow_threads(move || engine.seasons_in_year(year))
            .map_err(ephemeris_err_to_pyerr)?;
        seasons_to_pydict(py, &seasons)
    }

    /// Computes the exact Human Design `-88°` prenatal solar arc epoch before `birth_jd_ut1` (`R9.2`).
    #[pyo3(signature = (birth_jd_ut1, jd2=0.0, scale="UT1"))]
    pub fn human_design_epoch(
        &self,
        py: Python<'_>,
        birth_jd_ut1: f64,
        jd2: f64,
        scale: &str,
    ) -> PyResult<PyObject> {
        let birth_in = parse_time_input(scale, birth_jd_ut1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let ev = py
            .allow_threads(move || engine.human_design_epoch(birth_in))
            .map_err(ephemeris_err_to_pyerr)?;
        longitude_crossing_to_pydict(py, &ev)
    }

    /// Finds the next (`forward == true`) or previous (`forward == false`) instant from `start_jd_ut1`
    /// when geocentric apparent tropical lunar-solar elongation reaches `target_elongation_deg` (`R9.2`).
    #[pyo3(signature = (start_jd_ut1, target_elongation_deg, forward=true, jd2=0.0, scale="UT1"))]
    pub fn search_lunar_phase(
        &self,
        py: Python<'_>,
        start_jd_ut1: f64,
        target_elongation_deg: f64,
        forward: bool,
        jd2: f64,
        scale: &str,
    ) -> PyResult<PyObject> {
        let start_in = parse_time_input(scale, start_jd_ut1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let ev = py
            .allow_threads(move || {
                engine.search_lunar_phase(start_in, target_elongation_deg, forward)
            })
            .map_err(ephemeris_err_to_pyerr)?;
        lunar_phase_event_to_pydict(py, &ev)
    }

    /// Returns all 4 principal lunar phases occurring within `[start_jd_ut1, end_jd_ut1]` (`R9.2`).
    #[pyo3(signature = (start_jd_ut1, end_jd_ut1, scale="UT1"))]
    pub fn lunar_phases_in_range(
        &self,
        py: Python<'_>,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        scale: &str,
    ) -> PyResult<PyObject> {
        let start_in = parse_time_input(scale, start_jd_ut1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let end_in = parse_time_input(scale, end_jd_ut1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let events = py
            .allow_threads(move || engine.lunar_phases_in_range(start_in, end_in))
            .map_err(ephemeris_err_to_pyerr)?;
        let list = PyList::empty(py);
        for ev in &events {
            list.append(lunar_phase_event_to_pydict(py, ev)?)?;
        }
        Ok(list.into())
    }

    /// Computes the complete lunar phase, illumination, distance, age, and 8/28-phase summary at `jd_ut1` (`R0`, `R2.2`, `R9.2`).
    #[pyo3(signature = (jd_ut1, jd2=0.0, scale="UT1"))]
    pub fn lunar_phase_summary(
        &self,
        py: Python<'_>,
        jd_ut1: f64,
        jd2: f64,
        scale: &str,
    ) -> PyResult<PyObject> {
        let t_in = parse_time_input(scale, jd_ut1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let summary = py
            .allow_threads(move || engine.lunar_phase_summary(t_in))
            .map_err(ephemeris_err_to_pyerr)?;
        lunar_phase_summary_to_pydict(py, &summary)
    }

    /// Alias for [`Self::lunar_phase_summary`].
    #[pyo3(signature = (jd_ut1, jd2=0.0, scale="UT1"))]
    pub fn lunar_phase_detail(
        &self,
        py: Python<'_>,
        jd_ut1: f64,
        jd2: f64,
        scale: &str,
    ) -> PyResult<PyObject> {
        self.lunar_phase_summary(py, jd_ut1, jd2, scale)
    }

    /// Finds all retrograde and direct stations ($\dot{\lambda}(t) = 0$) of `body` within `[start_jd_ut1, end_jd_ut1]` (`R9.2`).
    #[pyo3(signature = (body, start_jd_ut1, end_jd_ut1, options=None, scale="UT1"))]
    pub fn search_stations_in_range(
        &self,
        py: Python<'_>,
        body: &str,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<&PyCalcOptions>,
        scale: &str,
    ) -> PyResult<PyObject> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let start_in = parse_time_input(scale, start_jd_ut1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let end_in = parse_time_input(scale, end_jd_ut1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let events = py
            .allow_threads(move || engine.search_stations_in_range(bid, start_in, end_in, &opts))
            .map_err(ephemeris_err_to_pyerr)?;
        let list = PyList::empty(py);
        for ev in &events {
            list.append(station_event_to_pydict(py, ev)?)?;
        }
        Ok(list.into())
    }

    /// Alias for [`Self::search_stations_in_range`].
    #[pyo3(signature = (body, start_jd_ut1, end_jd_ut1, options=None, scale="UT1"))]
    pub fn stations_in_range(
        &self,
        py: Python<'_>,
        body: &str,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<&PyCalcOptions>,
        scale: &str,
    ) -> PyResult<PyObject> {
        self.search_stations_in_range(py, body, start_jd_ut1, end_jd_ut1, options, scale)
    }

    /// Finds all northern and southern declination extrema / lunistices ($\dot{\delta}(t) = 0$) of `body`
    /// within `[start_jd_ut1, end_jd_ut1]` (`R9.2`).
    #[pyo3(signature = (body, start_jd_ut1, end_jd_ut1, options=None, scale="UT1"))]
    pub fn search_declination_extrema_in_range(
        &self,
        py: Python<'_>,
        body: &str,
        start_jd_ut1: f64,
        end_jd_ut1: f64,
        options: Option<&PyCalcOptions>,
        scale: &str,
    ) -> PyResult<PyObject> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let start_in = parse_time_input(scale, start_jd_ut1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let end_in = parse_time_input(scale, end_jd_ut1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let opts = options.map_or_else(CalcOptions::default, |o| o.inner);
        let engine = self.get_engine()?;
        let events = py
            .allow_threads(move || {
                engine.search_declination_extrema_in_range(bid, start_in, end_in, &opts)
            })
            .map_err(ephemeris_err_to_pyerr)?;
        let list = PyList::empty(py);
        for ev in &events {
            list.append(declination_extremum_to_pydict(py, ev)?)?;
        }
        Ok(list.into())
    }

    /// Computes Bernadette Brady diurnal fixed-star parans (`R6.2`) for `observer` starting at `jd_ut1`.
    #[pyo3(signature = (jd_ut1, observer, orb_minutes=4.0, bodies=None, star_names=None, jd2=0.0, scale="UT1"))]
    #[allow(clippy::too_many_arguments)]
    pub fn find_parans(
        &self,
        py: Python<'_>,
        jd_ut1: f64,
        observer: &PyObserver,
        orb_minutes: f64,
        bodies: Option<Vec<String>>,
        star_names: Option<Vec<String>>,
        jd2: f64,
        scale: &str,
    ) -> PyResult<PyObject> {
        let start_in = parse_time_input(scale, jd_ut1, jd2).map_err(ephemeris_err_to_pyerr)?;
        let core_obs = Observer::new(
            observer.latitude_deg,
            observer.longitude_deg,
            observer.altitude_m,
        )
        .map_err(ephemeris_err_to_pyerr)?;
        let parsed_bodies = match bodies {
            Some(names) => names
                .iter()
                .map(|n| parse_body_id(n))
                .collect::<Result<Vec<_>, _>>()
                .map_err(ephemeris_err_to_pyerr)?,
            None => ParanOptions::default().bodies,
        };
        let opts = ParanOptions {
            orb_minutes,
            bodies: parsed_bodies,
            star_names,
        };
        let engine = self.get_engine()?;
        let events = py
            .allow_threads(move || engine.find_parans(start_in, core_obs, &opts))
            .map_err(ephemeris_err_to_pyerr)?;
        let list = PyList::empty(py);
        for ev in &events {
            list.append(paran_event_to_pydict(py, ev)?)?;
        }
        Ok(list.into())
    }

    /// Finds the first exact return of `body` to its natal ecliptic longitude after `search_after_jd`.
    #[pyo3(signature = (body, natal_jd, search_after_jd, scale="UT1", zodiac="tropical", ayanamsa=None, options=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn planetary_return(
        &self,
        py: Python<'_>,
        body: &str,
        natal_jd: f64,
        search_after_jd: f64,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
        options: Option<&PyCalcOptions>,
    ) -> PyResult<PyObject> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let natal_in = parse_time_input(scale, natal_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let after_in =
            parse_time_input(scale, search_after_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let opts = match options {
            Some(o) => o.inner,
            None => resolve_wave_h_calc_options(zodiac, ayanamsa)
                .map_err(ephemeris_err_to_pyerr)?,
        };
        let engine = self.get_engine()?;
        let ev = py
            .allow_threads(move || engine.planetary_return(bid, natal_in, after_in, opts))
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &ev)
    }

    /// Finds all exact returns of `body` to its natal ecliptic longitude within `[start_jd, end_jd]`,
    /// annotating multi-pass retrograde triple returns with `pass_index = 1, 2, 3`.
    #[pyo3(signature = (body, natal_jd, start_jd, end_jd, scale="UT1", zodiac="tropical", ayanamsa=None, options=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn planetary_returns_in_range(
        &self,
        py: Python<'_>,
        body: &str,
        natal_jd: f64,
        start_jd: f64,
        end_jd: f64,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
        options: Option<&PyCalcOptions>,
    ) -> PyResult<PyObject> {
        let bid = parse_body_id(body).map_err(ephemeris_err_to_pyerr)?;
        let natal_in = parse_time_input(scale, natal_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let start_in = parse_time_input(scale, start_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let end_in = parse_time_input(scale, end_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let opts = match options {
            Some(o) => o.inner,
            None => resolve_wave_h_calc_options(zodiac, ayanamsa)
                .map_err(ephemeris_err_to_pyerr)?,
        };
        let engine = self.get_engine()?;
        let events = py
            .allow_threads(move || {
                engine.planetary_returns_in_range(bid, natal_in, start_in, end_in, opts)
            })
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &events)
    }

    /// Computes all active ecliptic longitude aspects between pairs of `bodies` at `jd`.
    #[pyo3(signature = (jd, bodies=None, preset="major", orb_scale=1.0, scale="UT1", zodiac="tropical", ayanamsa=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn calc_chart_aspects(
        &self,
        py: Python<'_>,
        jd: f64,
        bodies: Option<Vec<String>>,
        preset: &str,
        orb_scale: f64,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
    ) -> PyResult<PyObject> {
        let t_in = parse_time_input(scale, jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let bids = parse_optional_bodies(bodies, default_chart_bodies)
            .map_err(ephemeris_err_to_pyerr)?;
        let aspect_preset = parse_aspect_preset(preset).map_err(ephemeris_err_to_pyerr)?;
        let opts =
            resolve_wave_h_calc_options(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let hits = py
            .allow_threads(move || {
                engine.calc_chart_aspects(t_in, &bids, aspect_preset, orb_scale, opts)
            })
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &hits)
    }

    /// Computes all active cross-chart (synastry / transit-to-natal) aspects between `bodies1` at `jd1`
    /// and `bodies2` at `jd2`.
    #[pyo3(signature = (jd1, bodies1, jd2, bodies2, preset="major", orb_scale=1.0, scale="UT1", zodiac="tropical", ayanamsa=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn calc_synastry_aspects(
        &self,
        py: Python<'_>,
        jd1: f64,
        bodies1: Option<Vec<String>>,
        jd2: f64,
        bodies2: Option<Vec<String>>,
        preset: &str,
        orb_scale: f64,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
    ) -> PyResult<PyObject> {
        let t_in1 = parse_time_input(scale, jd1, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let t_in2 = parse_time_input(scale, jd2, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let bids1 = parse_optional_bodies(bodies1, default_chart_bodies)
            .map_err(ephemeris_err_to_pyerr)?;
        let bids2 = parse_optional_bodies(bodies2, default_chart_bodies)
            .map_err(ephemeris_err_to_pyerr)?;
        let aspect_preset = parse_aspect_preset(preset).map_err(ephemeris_err_to_pyerr)?;
        let opts =
            resolve_wave_h_calc_options(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let hits = py
            .allow_threads(move || {
                engine.calc_synastry_aspects(
                    t_in1,
                    &bids1,
                    t_in2,
                    &bids2,
                    aspect_preset,
                    orb_scale,
                    opts,
                )
            })
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &hits)
    }

    /// Computes all active equatorial declination aspects (Parallel and Contra-Parallel) at `jd`.
    #[pyo3(signature = (jd, bodies=None, max_orb_deg=1.2, scale="UT1"))]
    pub fn calc_declination_aspects(
        &self,
        py: Python<'_>,
        jd: f64,
        bodies: Option<Vec<String>>,
        max_orb_deg: f64,
        scale: &str,
    ) -> PyResult<PyObject> {
        let t_in = parse_time_input(scale, jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let bids = parse_optional_bodies(bodies, default_chart_bodies)
            .map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let hits = py
            .allow_threads(move || {
                engine.calc_declination_aspects(t_in, &bids, max_orb_deg, CalcOptions::default())
            })
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &hits)
    }

    /// Searches forward or backward from `start_jd` for the next exact aspect of `target_angle_deg`
    /// between `body1` and `body2`.
    #[pyo3(signature = (body1, body2, target_angle_deg, start_jd, forward=true, scale="UT1", zodiac="tropical", ayanamsa=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn search_exact_aspect(
        &self,
        py: Python<'_>,
        body1: &str,
        body2: &str,
        target_angle_deg: f64,
        start_jd: f64,
        forward: bool,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
    ) -> PyResult<PyObject> {
        let b1 = parse_body_id(body1).map_err(ephemeris_err_to_pyerr)?;
        let b2 = parse_body_id(body2).map_err(ephemeris_err_to_pyerr)?;
        let start_in = parse_time_input(scale, start_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let opts =
            resolve_wave_h_calc_options(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let ev = py
            .allow_threads(move || {
                engine.search_exact_aspect(b1, b2, target_angle_deg, start_in, forward, opts)
            })
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &ev)
    }

    /// Finds all exact mutual aspect events between pairs of `bodies` in `[start_jd, end_jd]`.
    #[pyo3(signature = (start_jd, end_jd, bodies=None, preset="major", scale="UT1", zodiac="tropical", ayanamsa=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn mutual_aspects_in_range(
        &self,
        py: Python<'_>,
        start_jd: f64,
        end_jd: f64,
        bodies: Option<Vec<String>>,
        preset: &str,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
    ) -> PyResult<PyObject> {
        let start_in = parse_time_input(scale, start_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let end_in = parse_time_input(scale, end_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let bids = parse_optional_bodies(bodies, default_major_planets_10)
            .map_err(ephemeris_err_to_pyerr)?;
        let aspect_preset = parse_aspect_preset(preset).map_err(ephemeris_err_to_pyerr)?;
        let opts =
            resolve_wave_h_calc_options(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let events = py
            .allow_threads(move || {
                engine.mutual_aspects_in_range(&bids, aspect_preset, start_in, end_in, opts)
            })
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &events)
    }

    /// Alias for [`Self::mutual_aspects_in_range`].
    #[pyo3(signature = (start_jd, end_jd, bodies=None, preset="major", scale="UT1", zodiac="tropical", ayanamsa=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn dynamic_aspects_in_range(
        &self,
        py: Python<'_>,
        start_jd: f64,
        end_jd: f64,
        bodies: Option<Vec<String>>,
        preset: &str,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
    ) -> PyResult<PyObject> {
        self.mutual_aspects_in_range(py, start_jd, end_jd, bodies, preset, scale, zodiac, ayanamsa)
    }

    /// Finds all exact transits of `transiting_bodies` to natal positions of `natal_bodies` (evaluated at `natal_jd`)
    /// within `[start_jd, end_jd]`.
    #[pyo3(signature = (natal_jd, start_jd, end_jd, transiting_bodies=None, natal_bodies=None, preset="major", scale="UT1", zodiac="tropical", ayanamsa=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn natal_transits_in_range(
        &self,
        py: Python<'_>,
        natal_jd: f64,
        start_jd: f64,
        end_jd: f64,
        transiting_bodies: Option<Vec<String>>,
        natal_bodies: Option<Vec<String>>,
        preset: &str,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
    ) -> PyResult<PyObject> {
        let natal_in = parse_time_input(scale, natal_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let start_in = parse_time_input(scale, start_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let end_in = parse_time_input(scale, end_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let trans_bids = parse_optional_bodies(transiting_bodies, default_major_planets_10)
            .map_err(ephemeris_err_to_pyerr)?;
        let natal_bids = parse_optional_bodies(natal_bodies, default_major_planets_10)
            .map_err(ephemeris_err_to_pyerr)?;
        let aspect_preset = parse_aspect_preset(preset).map_err(ephemeris_err_to_pyerr)?;
        let opts =
            resolve_wave_h_calc_options(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let events = py
            .allow_threads(move || {
                engine.natal_transits_in_range(
                    &trans_bids,
                    natal_in,
                    &natal_bids,
                    aspect_preset,
                    start_in,
                    end_in,
                    opts,
                )
            })
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &events)
    }

    /// Computes all Void-of-Course Moon windows overlapping `[start_jd, end_jd]`.
    #[pyo3(signature = (start_jd, end_jd, scale="UT1", zodiac="tropical", ayanamsa=None))]
    pub fn void_of_course_moon_in_range(
        &self,
        py: Python<'_>,
        start_jd: f64,
        end_jd: f64,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
    ) -> PyResult<PyObject> {
        let start_in = parse_time_input(scale, start_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let end_in = parse_time_input(scale, end_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let opts =
            resolve_wave_h_calc_options(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let windows = py
            .allow_threads(move || {
                engine.void_of_course_moon_in_range(start_in, end_in, None, opts)
            })
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &windows)
    }

    /// Evaluates the complete five-limbed Vedic Panchanga (Tithi, Nakshatra, Yoga, Karana, Vara) at `jd`.
    #[pyo3(signature = (jd, lat_deg, lon_deg, alt_m=0.0, ayanamsa="lahiri", scale="UT1"))]
    #[allow(clippy::too_many_arguments)]
    pub fn panchanga(
        &self,
        py: Python<'_>,
        jd: f64,
        lat_deg: f64,
        lon_deg: f64,
        alt_m: f64,
        ayanamsa: &str,
        scale: &str,
    ) -> PyResult<PyObject> {
        let t_in = parse_time_input(scale, jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(lat_deg, lon_deg, alt_m).map_err(ephemeris_err_to_pyerr)?;
        let ayan_model = AyanamsaModel::from_code(ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.panchanga(t_in, obs, ayan_model))
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &res)
    }

    /// Computes the 120-year Vimshottari Dasha tree (9 Mahadashas and 81 Antardashas) from `birth_jd`.
    #[pyo3(signature = (birth_jd, ayanamsa="lahiri", year_length_days=None, scale="UT1"))]
    pub fn vimshottari_dasha(
        &self,
        py: Python<'_>,
        birth_jd: f64,
        ayanamsa: &str,
        year_length_days: Option<f64>,
        scale: &str,
    ) -> PyResult<PyObject> {
        let birth_in = parse_time_input(scale, birth_jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let ayan_model = AyanamsaModel::from_code(ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || {
                engine.vimshottari_dasha(birth_in, ayan_model, year_length_days)
            })
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &res)
    }

    /// Evaluates a single Parashari divisional chart (*Varga*) `D-n` (`division`) for `sidereal_lon_deg`.
    pub fn calc_varga(
        &self,
        py: Python<'_>,
        sidereal_lon_deg: f64,
        division: u16,
    ) -> PyResult<PyObject> {
        let res = core_calc_varga(sidereal_lon_deg, division).map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &res)
    }

    /// Evaluates all 16 Parashari *Shodashavarga* divisional charts for `sidereal_lon_deg`.
    pub fn calc_shodashavarga(&self, py: Python<'_>, sidereal_lon_deg: f64) -> PyResult<PyObject> {
        let res = core_calc_shodashavarga(sidereal_lon_deg).map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &res)
    }

    /// Determines whether a chart at `jd` and `(lat_deg, lon_deg, alt_m)` is a Day or Night chart (*Sect*).
    #[pyo3(signature = (jd, lat_deg, lon_deg, alt_m=0.0, scale="UT1"))]
    pub fn chart_sect(
        &self,
        py: Python<'_>,
        jd: f64,
        lat_deg: f64,
        lon_deg: f64,
        alt_m: f64,
        scale: &str,
    ) -> PyResult<PyObject> {
        let t_in = parse_time_input(scale, jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(lat_deg, lon_deg, alt_m).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.chart_sect(t_in, obs))
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &res)
    }

    /// Computes traditional Chaldean planetary day/night hours for `jd` at `(lat_deg, lon_deg, alt_m)`.
    #[pyo3(signature = (jd, lat_deg, lon_deg, alt_m=0.0, scale="UT1"))]
    pub fn planetary_hours(
        &self,
        py: Python<'_>,
        jd: f64,
        lat_deg: f64,
        lon_deg: f64,
        alt_m: f64,
        scale: &str,
    ) -> PyResult<PyObject> {
        let t_in = parse_time_input(scale, jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(lat_deg, lon_deg, alt_m).map_err(ephemeris_err_to_pyerr)?;
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.planetary_hours(t_in, obs))
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &res)
    }

    /// Computes the 7 Hermetic Lots / Arabic Parts (both sect-reversed and non-reversed) at `jd`.
    #[pyo3(signature = (jd, lat_deg, lon_deg, alt_m=0.0, system=None, scale="UT1", zodiac="tropical", ayanamsa=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn calc_lots(
        &self,
        py: Python<'_>,
        jd: f64,
        lat_deg: f64,
        lon_deg: f64,
        alt_m: f64,
        system: Option<&Bound<'_, PyAny>>,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
    ) -> PyResult<PyObject> {
        let t_in = parse_time_input(scale, jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(lat_deg, lon_deg, alt_m).map_err(ephemeris_err_to_pyerr)?;
        let ch = extract_hsys_char(system)?;
        let house_sys = HouseSystem::from_code(ch).map_err(ephemeris_err_to_pyerr)?;
        let zod = resolve_wave_h_zodiac(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let house_opts = HouseOptions {
            system: house_sys,
            zodiac: zod,
            ..HouseOptions::default()
        };
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || engine.calc_lots(t_in, obs, house_opts))
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &res)
    }

    /// Evaluates a complete one-call chart snapshot (positions, zodiac/nakshatra/Human Design placements,
    /// houses, sect, Hermetic lots, lunar phase summary, ecliptic aspects, and declination aspects).
    #[pyo3(signature = (jd, lat_deg, lon_deg, alt_m=0.0, bodies=None, system=None, preset="major", orb_scale=1.0, scale="UT1", zodiac="tropical", ayanamsa=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn chart_snapshot(
        &self,
        py: Python<'_>,
        jd: f64,
        lat_deg: f64,
        lon_deg: f64,
        alt_m: f64,
        bodies: Option<Vec<String>>,
        system: Option<&Bound<'_, PyAny>>,
        preset: &str,
        orb_scale: f64,
        scale: &str,
        zodiac: &str,
        ayanamsa: Option<&str>,
    ) -> PyResult<PyObject> {
        let t_in = parse_time_input(scale, jd, 0.0).map_err(ephemeris_err_to_pyerr)?;
        let obs = Observer::new(lat_deg, lon_deg, alt_m).map_err(ephemeris_err_to_pyerr)?;
        let bids = parse_optional_bodies(bodies, default_chart_bodies)
            .map_err(ephemeris_err_to_pyerr)?;
        let ch = extract_hsys_char(system)?;
        let house_sys = HouseSystem::from_code(ch).map_err(ephemeris_err_to_pyerr)?;
        let aspect_preset = parse_aspect_preset(preset).map_err(ephemeris_err_to_pyerr)?;
        let zod = resolve_wave_h_zodiac(zodiac, ayanamsa).map_err(ephemeris_err_to_pyerr)?;
        let house_opts = HouseOptions {
            system: house_sys,
            zodiac: zod,
            ..HouseOptions::default()
        };
        let calc_opts = CalcOptions {
            zodiac: zod,
            ..CalcOptions::default()
        };
        let engine = self.get_engine()?;
        let res = py
            .allow_threads(move || {
                engine.chart_snapshot(
                    t_in,
                    obs,
                    &bids,
                    house_opts,
                    calc_opts,
                    aspect_preset,
                    orb_scale,
                )
            })
            .map_err(ephemeris_err_to_pyerr)?;
        serialize_to_attr_dict(py, &res)
    }

    /// Maps an ecliptic longitude `lon_deg` into its Zodiac Sign, Decan, Nakshatra + Pada, and Human Design Gate.Line.Color.Tone.Base.
    pub fn map_zodiac_placement(&self, py: Python<'_>, lon_deg: f64) -> PyResult<PyObject> {
        let res = core_map_zodiac_placement(lon_deg);
        serialize_to_attr_dict(py, &res)
    }
}

fn longitude_crossing_to_pydict(
    py: Python<'_>,
    ev: &LongitudeCrossingEvent,
) -> PyResult<PyObject> {
    let d = PyDict::new(py);
    d.set_item("body", body_display_name(ev.body))?;
    d.set_item("target_lon_deg", ev.target_lon_deg)?;
    d.set_item("jd_tt", ev.jd_tt)?;
    d.set_item("jd_ut1", ev.jd_ut1)?;
    d.set_item("delta_t_sec", ev.delta_t_sec)?;
    d.set_item("actual_lon_deg", ev.actual_lon_deg)?;
    d.set_item("lat_deg", ev.lat_deg)?;
    d.set_item("dist_au", ev.dist_au)?;
    d.set_item("dist_km", ev.dist_km)?;
    d.set_item("lon_speed_deg_day", ev.lon_speed_deg_day)?;
    d.set_item("is_retrograde", ev.is_retrograde)?;
    dict_to_attr_dict(py, &d)
}

fn solar_term_to_pydict(py: Python<'_>, ev: &SolarTermEvent) -> PyResult<PyObject> {
    let d = PyDict::new(py);
    d.set_item("term_index", ev.term_index)?;
    d.set_item("target_lon_deg", ev.target_lon_deg)?;
    d.set_item("pinyin_name", ev.pinyin_name.as_str())?;
    d.set_item("chinese_name", ev.chinese_name.as_str())?;
    d.set_item("english_name", ev.english_name.as_str())?;
    d.set_item("is_major_term", ev.is_major_term)?;
    d.set_item("jd_tt", ev.jd_tt)?;
    d.set_item("jd_ut1", ev.jd_ut1)?;
    d.set_item("delta_t_sec", ev.delta_t_sec)?;
    dict_to_attr_dict(py, &d)
}

fn seasons_to_pydict(py: Python<'_>, s: &SeasonsOfYear) -> PyResult<PyObject> {
    let d = PyDict::new(py);
    d.set_item("year", s.year)?;
    d.set_item("vernal_equinox_jd_ut1", s.vernal_equinox_jd_ut1)?;
    d.set_item("vernal_equinox_jd_tt", s.vernal_equinox_jd_tt)?;
    d.set_item("summer_solstice_jd_ut1", s.summer_solstice_jd_ut1)?;
    d.set_item("summer_solstice_jd_tt", s.summer_solstice_jd_tt)?;
    d.set_item("autumnal_equinox_jd_ut1", s.autumnal_equinox_jd_ut1)?;
    d.set_item("autumnal_equinox_jd_tt", s.autumnal_equinox_jd_tt)?;
    d.set_item("winter_solstice_jd_ut1", s.winter_solstice_jd_ut1)?;
    d.set_item("winter_solstice_jd_tt", s.winter_solstice_jd_tt)?;
    dict_to_attr_dict(py, &d)
}

fn lunar_phase_event_to_pydict(py: Python<'_>, ev: &LunarPhaseEvent) -> PyResult<PyObject> {
    let d = PyDict::new(py);
    let kind_str: Option<&'static str> = ev.kind.map(|k| match k {
        LunarPhaseKind::NewMoon => "new_moon",
        LunarPhaseKind::FirstQuarter => "first_quarter",
        LunarPhaseKind::FullMoon => "full_moon",
        LunarPhaseKind::LastQuarter => "last_quarter",
    });
    d.set_item("kind", kind_str)?;
    d.set_item("target_elongation_deg", ev.target_elongation_deg)?;
    d.set_item("jd_tt", ev.jd_tt)?;
    d.set_item("jd_ut1", ev.jd_ut1)?;
    d.set_item("delta_t_sec", ev.delta_t_sec)?;
    d.set_item("moon_lon_deg", ev.moon_lon_deg)?;
    d.set_item("sun_lon_deg", ev.sun_lon_deg)?;
    d.set_item("moon_dist_au", ev.moon_dist_au)?;
    d.set_item("moon_dist_km", ev.moon_dist_km)?;
    dict_to_attr_dict(py, &d)
}

fn lunar_phase_summary_to_pydict(py: Python<'_>, s: &LunarPhaseSummary) -> PyResult<PyObject> {
    let d = PyDict::new(py);
    d.set_item("jd_tt", s.jd_tt)?;
    d.set_item("jd_ut1", s.jd_ut1)?;
    d.set_item("elongation_deg", s.elongation_deg)?;
    d.set_item("illumination_fraction", s.illumination_fraction)?;
    d.set_item("phase_angle_deg", s.phase_angle_deg)?;
    d.set_item("moon_lon_deg", s.moon_lon_deg)?;
    d.set_item("moon_lat_deg", s.moon_lat_deg)?;
    d.set_item("moon_dist_au", s.moon_dist_au)?;
    d.set_item("moon_dist_km", s.moon_dist_km)?;
    d.set_item("apparent_diameter_deg", s.apparent_diameter_deg)?;
    d.set_item("magnitude_v", s.magnitude_v)?;
    d.set_item("age_days", s.age_days)?;
    d.set_item("synodic_fraction", s.synodic_fraction)?;
    d.set_item("phase_8_index", s.phase_8_index)?;
    d.set_item("phase_8_name", s.phase_8_name.as_str())?;
    d.set_item("phase_28_index", s.phase_28_index)?;
    d.set_item("prev_new_moon_jd_ut1", s.prev_new_moon_jd_ut1)?;
    d.set_item("next_new_moon_jd_ut1", s.next_new_moon_jd_ut1)?;
    d.set_item("next_first_quarter_jd_ut1", s.next_first_quarter_jd_ut1)?;
    d.set_item("next_full_moon_jd_ut1", s.next_full_moon_jd_ut1)?;
    d.set_item("next_last_quarter_jd_ut1", s.next_last_quarter_jd_ut1)?;
    dict_to_attr_dict(py, &d)
}

fn station_event_to_pydict(py: Python<'_>, ev: &StationEvent) -> PyResult<PyObject> {
    let d = PyDict::new(py);
    let kind_str = match ev.kind {
        StationKind::StationRetrograde => "station_retrograde",
        StationKind::StationDirect => "station_direct",
    };
    d.set_item("body", body_display_name(ev.body))?;
    d.set_item("kind", kind_str)?;
    d.set_item("jd_tt", ev.jd_tt)?;
    d.set_item("jd_ut1", ev.jd_ut1)?;
    d.set_item("delta_t_sec", ev.delta_t_sec)?;
    d.set_item("lon_deg", ev.lon_deg)?;
    d.set_item("lat_deg", ev.lat_deg)?;
    d.set_item("dist_au", ev.dist_au)?;
    d.set_item("dist_km", ev.dist_km)?;
    dict_to_attr_dict(py, &d)
}

fn declination_extremum_to_pydict(
    py: Python<'_>,
    ev: &DeclinationExtremumEvent,
) -> PyResult<PyObject> {
    let d = PyDict::new(py);
    let kind_str = match ev.kind {
        DeclinationExtremumKind::NorthernMaximum => "northern_maximum",
        DeclinationExtremumKind::SouthernMinimum => "southern_minimum",
    };
    d.set_item("body", body_display_name(ev.body))?;
    d.set_item("kind", kind_str)?;
    d.set_item("jd_tt", ev.jd_tt)?;
    d.set_item("jd_ut1", ev.jd_ut1)?;
    d.set_item("delta_t_sec", ev.delta_t_sec)?;
    d.set_item("ra_deg", ev.ra_deg)?;
    d.set_item("dec_deg", ev.dec_deg)?;
    d.set_item("dist_au", ev.dist_au)?;
    d.set_item("dist_km", ev.dist_km)?;
    dict_to_attr_dict(py, &d)
}

fn paran_event_to_pydict(py: Python<'_>, ev: &ParanEvent) -> PyResult<PyObject> {
    let d = PyDict::new(py);
    d.set_item("star_name", ev.star_name.as_str())?;
    d.set_item("star_bayer", ev.star_bayer.as_str())?;
    d.set_item("star_hip", ev.star_hip)?;
    d.set_item("planet", ev.planet_name.as_str())?;
    d.set_item("paran_type", ev.paran_type.as_str())?;
    d.set_item("star_angle", ev.star_angle.as_str())?;
    d.set_item("planet_angle", ev.planet_angle.as_str())?;
    d.set_item("same_angle", ev.same_angle)?;
    d.set_item("star_event_jd_ut1", ev.star_event_jd_ut1)?;
    d.set_item("planet_event_jd_ut1", ev.planet_event_jd_ut1)?;
    d.set_item("time_orb_minutes", ev.time_orb_minutes)?;
    dict_to_attr_dict(py, &d)
}

#[cfg(feature = "swe-compat")]
fn resolve_nod_aps_args(
    arg3: Option<u32>,
    arg4: Option<u32>,
    kw_iflag: Option<u32>,
    kw_flags: Option<u32>,
    kw_method: Option<u32>,
) -> (u32, u32) {
    let kw_fl = kw_iflag.or(kw_flags);
    match (kw_fl, kw_method) {
        (Some(fl), Some(m)) => (fl, m),
        (Some(fl), None) => (fl, arg3.unwrap_or(0)),
        (None, Some(m)) => (arg3.unwrap_or(258), m),
        (None, None) => match (arg3, arg4) {
            (None, None) => (258, 0),
            (Some(a3), None) => {
                // Disambiguate single optional integer: method bits in swe_nod_aps are {1..=7, 256, 257, 259..=263},
                // whereas 258 (SEFLG_SWIEPH | SEFLG_SPEED) or >= 8 (outside 256..=263) are calculation flags.
                if a3 != 258 && (a3 & !0x107) == 0 && a3 != 0 {
                    (258, a3)
                } else {
                    (a3, 0)
                }
            }
            (Some(a3), Some(a4)) => {
                // If a3 is a flag combination (like 258) and a4 is a method bitmask (0..=7, 256..=263),
                // interpret as (iflag, method); otherwise interpret as pyswisseph (method, flags).
                if ((a3 & !0x107) != 0 || a3 == 258) && (a4 & !0x107) == 0 && a4 != 258 {
                    (a3, a4)
                } else {
                    (a4, a3)
                }
            }
            (None, Some(a4)) => (258, a4),
        },
    }
}

/// Drop-in compatibility adapter (`SweAdapter`) for migrating legacy `pyswisseph` callers to `lxephemeris`.
#[cfg(feature = "swe-compat")]
#[pyclass(name = "SweAdapter", frozen)]
pub struct PySweAdapter {
    inner: Mutex<SweAdapter>,
}

#[cfg(feature = "swe-compat")]
#[pymethods]
impl PySweAdapter {
    /// Creates a new `SweAdapter`. Accepts `None` (loads from `LXEPHEMERIS_DATA_DIR`, `./data`, or Moshier fallback),
    /// a `data_dir` path string, or an existing `Engine` instance.
    #[new]
    #[pyo3(signature = (data_dir=None))]
    pub fn new(py: Python<'_>, data_dir: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        match data_dir {
            None => {
                let eng = py.allow_threads(|| Engine::from_env_or_moshier(EngineConfig::default()));
                Ok(Self {
                    inner: Mutex::new(SweAdapter::new(Arc::new(eng))),
                })
            }
            Some(arg) => {
                if let Ok(py_eng) = arg.extract::<PyRef<'_, PyEngine>>() {
                    Self::from_engine(&py_eng)
                } else if let Ok(dir) = arg.extract::<String>() {
                    Self::from_data_dir(py, dir)
                } else {
                    Err(InvalidInputError::new_err(
                        "[InvalidInput] SweAdapter(data_dir=...) expects None, a directory path str, or an Engine instance",
                    ))
                }
            }
        }
    }

    /// Creates a `SweAdapter` wrapping an existing `Engine`.
    #[staticmethod]
    pub fn from_engine(engine: &PyEngine) -> PyResult<Self> {
        let eng = engine.get_engine()?;
        Ok(Self {
            inner: Mutex::new(SweAdapter::new(eng)),
        })
    }

    /// Creates a fileless Moshier-compatible adapter without external SPK kernels.
    #[staticmethod]
    pub fn new_moshier() -> Self {
        let eng = Arc::new(Engine::new_moshier(EngineConfig::default()));
        Self {
            inner: Mutex::new(SweAdapter::new(eng)),
        }
    }

    #[staticmethod]
    pub fn from_data_dir(py: Python<'_>, data_dir: String) -> PyResult<Self> {
        let eng = PyEngine::from_data_dir(py, data_dir, None, None, false)?;
        Self::from_engine(&eng)
    }

    /// Sets the active sidereal ayanamsa mode (`SE_SIDM_*`).
    #[pyo3(signature = (sid_mode, t0=0.0, ayan_t0=0.0))]
    pub fn set_sid_mode(&self, sid_mode: u32, t0: f64, ayan_t0: f64) -> PyResult<()> {
        let mut guard = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        if (sid_mode & 0xFF) == SE_SIDM_USER {
            guard
                .set_sid_mode_custom(sid_mode, t0, ayan_t0)
                .map_err(ephemeris_err_to_pyerr)
        } else {
            guard.set_sid_mode(sid_mode).map_err(ephemeris_err_to_pyerr)
        }
    }

    /// Sets the active polar fallback policy (`"strict"`, `"fallback_to_porphyry"`, `"fallback_whole_sign"`).
    pub fn set_polar_policy(&self, policy: &str) -> PyResult<()> {
        let pol = parse_polar_policy(policy).map_err(ephemeris_err_to_pyerr)?;
        let mut guard = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        guard.set_polar_policy(pol);
        Ok(())
    }

    /// Sets the topocentric observer (`lon_deg, lat_deg, alt_m`) for `SEFLG_TOPOCTR`.
    #[pyo3(signature = (lon_deg=None, lat_deg=None, alt_m=None, *, lon=None, lat=None, alt=None))]
    pub fn set_topo(
        &self,
        lon_deg: Option<f64>,
        lat_deg: Option<f64>,
        alt_m: Option<f64>,
        lon: Option<f64>,
        lat: Option<f64>,
        alt: Option<f64>,
    ) -> PyResult<()> {
        let lo = lon_deg
            .or(lon)
            .ok_or_else(|| PyValueError::new_err("set_topo missing required argument 'lon' / 'lon_deg'"))?;
        let la = lat_deg
            .or(lat)
            .ok_or_else(|| PyValueError::new_err("set_topo missing required argument 'lat' / 'lat_deg'"))?;
        let al = alt_m.or(alt).unwrap_or(0.0);
        let mut guard = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        guard
            .set_topo(lo, la, al)
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates $\Delta T = \text{TT} - \text{UT1}$ in days (`pyswisseph.deltat(jd_ut)` returns days, whereas `deltat_sec` returns seconds).
    pub fn deltat(&self, jd_ut: f64) -> PyResult<f64> {
        let guard = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(guard.deltat(jd_ut))
    }

    /// Evaluates $\Delta T = \text{TT} - \text{UT1}$ in seconds.
    pub fn deltat_sec(&self, jd_ut: f64) -> PyResult<f64> {
        let guard = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(guard.deltat_sec(jd_ut))
    }

    /// Evaluates Greenwich Apparent Sidereal Time (GAST) in hours `[0, 24)`.
    pub fn sidtime(&self, jd_ut: f64) -> PyResult<f64> {
        let guard = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(guard.sidtime(jd_ut))
    }

    /// Evaluates Greenwich Sidereal Time in hours with explicit obliquity and nutation in longitude (`deg`).
    pub fn sidtime0(&self, jd_ut: f64, eps_deg: f64, nut_lon_deg: f64) -> PyResult<f64> {
        let guard = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(guard.sidtime0(jd_ut, eps_deg, nut_lon_deg))
    }

    /// Evaluates the active sidereal ayanamsa at `jd_ut` (`UT1`) with `iflag`, returning `(retflag, ayanamsa_deg)`.
    #[pyo3(signature = (jd_ut=None, iflag=None, *, tjdut=None, flags=None))]
    pub fn get_ayanamsa_ex_ut(
        &self,
        py: Python<'_>,
        jd_ut: Option<f64>,
        iflag: Option<u32>,
        tjdut: Option<f64>,
        flags: Option<u32>,
    ) -> PyResult<(i32, f64)> {
        let j = jd_ut.or(tjdut).ok_or_else(|| {
            PyValueError::new_err("get_ayanamsa_ex_ut missing required argument 'jd_ut' / 'tjdut'")
        })?;
        let fl = iflag.or(flags).unwrap_or(0);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || adapter.get_ayanamsa_ex_ut(j, fl))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates the active sidereal ayanamsa at `jd_ut` (`UT1`) in degrees.
    pub fn get_ayanamsa_ut(&self, py: Python<'_>, jd_ut: f64) -> PyResult<f64> {
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(|| adapter.get_ayanamsa_ut(jd_ut))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates the active sidereal ayanamsa at `jd_et` (`TT`) in degrees.
    pub fn get_ayanamsa(&self, py: Python<'_>, jd_et: f64) -> PyResult<f64> {
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(|| adapter.get_ayanamsa(jd_et))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates the active sidereal ayanamsa at `jd_et` (`TT`) with `iflag`, returning `(retflag, ayanamsa_deg)`.
    #[pyo3(signature = (jd_et=None, iflag=None, *, tjdet=None, flags=None))]
    pub fn get_ayanamsa_ex(
        &self,
        py: Python<'_>,
        jd_et: Option<f64>,
        iflag: Option<u32>,
        tjdet: Option<f64>,
        flags: Option<u32>,
    ) -> PyResult<(i32, f64)> {
        let j = jd_et.or(tjdet).ok_or_else(|| {
            PyValueError::new_err("get_ayanamsa_ex missing required argument 'jd_et' / 'tjdet'")
        })?;
        let fl = iflag.or(flags).unwrap_or(0);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || adapter.get_ayanamsa_ex(j, fl))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates planetary, solar, lunar, or minor-body phenomena at `jd_ut` (`UT1`), returning `(attr[0..20], retflag)`.
    #[pyo3(signature = (jd_ut=None, ipl=None, iflag=None, *, tjdut=None, flags=None))]
    pub fn pheno_ut(
        &self,
        py: Python<'_>,
        jd_ut: Option<f64>,
        ipl: Option<i32>,
        iflag: Option<u32>,
        tjdut: Option<f64>,
        flags: Option<u32>,
    ) -> PyResult<([f64; 20], i32)> {
        let j = jd_ut
            .or(tjdut)
            .ok_or_else(|| PyValueError::new_err("pheno_ut missing required argument 'jd_ut' / 'tjdut'"))?;
        let ip = ipl.ok_or_else(|| PyValueError::new_err("pheno_ut missing required argument 'ipl'"))?;
        let fl = iflag.or(flags).unwrap_or(0);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || adapter.pheno_ut(j, ip, fl))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates planetary, solar, lunar, or minor-body phenomena at `jd_et` (`TT`), returning `(attr[0..20], retflag)`.
    #[pyo3(signature = (jd_et=None, ipl=None, iflag=None, *, tjdet=None, flags=None))]
    pub fn pheno(
        &self,
        py: Python<'_>,
        jd_et: Option<f64>,
        ipl: Option<i32>,
        iflag: Option<u32>,
        tjdet: Option<f64>,
        flags: Option<u32>,
    ) -> PyResult<([f64; 20], i32)> {
        let j = jd_et
            .or(tjdet)
            .ok_or_else(|| PyValueError::new_err("pheno missing required argument 'jd_et' / 'tjdet'"))?;
        let ip = ipl.ok_or_else(|| PyValueError::new_err("pheno missing required argument 'ipl'"))?;
        let fl = iflag.or(flags).unwrap_or(0);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || adapter.pheno(j, ip, fl))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates planetary or lunar orbital nodes and apsides at `jd_ut` (`UT1`), returning `(xnasc, xndsc, xperi, xaphe, retflag)`.
    /// Accepts both `pyswisseph` `(tjdut, ipl, method=1, flags=258)` and C-style `(jd_ut, ipl, iflag=258, method=0)`.
    #[pyo3(signature = (jd_ut=None, ipl=None, arg3=None, arg4=None, *, tjdut=None, iflag=None, flags=None, method=None))]
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn nod_aps_ut(
        &self,
        py: Python<'_>,
        jd_ut: Option<f64>,
        ipl: Option<i32>,
        arg3: Option<u32>,
        arg4: Option<u32>,
        tjdut: Option<f64>,
        iflag: Option<u32>,
        flags: Option<u32>,
        method: Option<u32>,
    ) -> PyResult<([f64; 6], [f64; 6], [f64; 6], [f64; 6], i32)> {
        let j = jd_ut
            .or(tjdut)
            .ok_or_else(|| PyValueError::new_err("nod_aps_ut missing required argument 'jd_ut' / 'tjdut'"))?;
        let ip = ipl.ok_or_else(|| PyValueError::new_err("nod_aps_ut missing required argument 'ipl'"))?;
        let (fl, m_bits) = resolve_nod_aps_args(arg3, arg4, iflag, flags, method);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || adapter.nod_aps_ut(j, ip, fl, m_bits))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates planetary or lunar orbital nodes and apsides at `jd_et` (`TT`), returning `(xnasc, xndsc, xperi, xaphe, retflag)`.
    /// Accepts both `pyswisseph` `(tjdet, ipl, method=1, flags=258)` and C-style `(jd_et, ipl, iflag=258, method=0)`.
    #[pyo3(signature = (jd_et=None, ipl=None, arg3=None, arg4=None, *, tjdet=None, iflag=None, flags=None, method=None))]
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn nod_aps(
        &self,
        py: Python<'_>,
        jd_et: Option<f64>,
        ipl: Option<i32>,
        arg3: Option<u32>,
        arg4: Option<u32>,
        tjdet: Option<f64>,
        iflag: Option<u32>,
        flags: Option<u32>,
        method: Option<u32>,
    ) -> PyResult<([f64; 6], [f64; 6], [f64; 6], [f64; 6], i32)> {
        let j = jd_et
            .or(tjdet)
            .ok_or_else(|| PyValueError::new_err("nod_aps missing required argument 'jd_et' / 'tjdet'"))?;
        let ip = ipl.ok_or_else(|| PyValueError::new_err("nod_aps missing required argument 'ipl'"))?;
        let (fl, m_bits) = resolve_nod_aps_args(arg3, arg4, iflag, flags, method);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || adapter.nod_aps(j, ip, fl, m_bits))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Searches for the next global lunar occultation of a planet or fixed star starting from `jd_ut`.
    /// Supports both `pyswisseph` 5-arg `(tjdut, body_or_star, flags=2, ifltype=0, backwards=False)` and
    /// 6-arg `(jd_ut, ipl, starname=None, iflag=258, ifltype=0, backward=False)` signatures.
    #[pyo3(signature = (jd_ut, ipl, arg3=None, arg4=None, arg5=None, arg6=None, *, starname=None, iflag=None, flags=None, ifltype=None, backward=None, backwards=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn lun_occult_when_glob(
        &self,
        py: Python<'_>,
        jd_ut: f64,
        ipl: &Bound<'_, PyAny>,
        arg3: Option<&Bound<'_, PyAny>>,
        arg4: Option<&Bound<'_, PyAny>>,
        arg5: Option<&Bound<'_, PyAny>>,
        arg6: Option<&Bound<'_, PyAny>>,
        starname: Option<&str>,
        iflag: Option<u32>,
        flags: Option<u32>,
        ifltype: Option<i32>,
        backward: Option<bool>,
        backwards: Option<bool>,
    ) -> PyResult<(i32, [f64; 10])> {
        let (body_int, body_star) = if let Ok(id) = ipl.extract::<i32>() {
            (id, None)
        } else if let Ok(s) = ipl.extract::<String>() {
            (0, Some(s))
        } else {
            return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                message: "body/ipl must be an integer SE_* body ID or a star name str".to_owned(),
            }));
        };
        let (pos_star, pos_flag, pos_type, pos_back) = if body_star.is_some() {
            // 5-arg pyswisseph with str star: (jd, "Spica", flags, ifltype, backwards)
            let fl = arg3.and_then(|a| a.extract::<u32>().ok());
            let it = arg4.and_then(|a| a.extract::<i32>().ok());
            let bw = arg5.and_then(|a| a.extract::<bool>().ok().or_else(|| a.extract::<i32>().ok().map(|v| v != 0)));
            (None, fl, it, bw)
        } else if let Some(a3) = arg3 {
            if a3.is_none() || a3.extract::<String>().is_ok() {
                // 6-arg C-style: (jd, ipl, starname, iflag, ifltype, backward)
                let st = a3.extract::<Option<String>>().ok().flatten();
                let fl = arg4.and_then(|a| a.extract::<u32>().ok());
                let it = arg5.and_then(|a| a.extract::<i32>().ok());
                let bw = arg6.and_then(|a| a.extract::<bool>().ok().or_else(|| a.extract::<i32>().ok().map(|v| v != 0)));
                (st, fl, it, bw)
            } else {
                // 5-arg pyswisseph with int body: (jd, ipl, flags, ifltype, backwards)
                let fl = a3.extract::<u32>().ok();
                let it = arg4.and_then(|a| a.extract::<i32>().ok());
                let bw = arg5.and_then(|a| a.extract::<bool>().ok().or_else(|| a.extract::<i32>().ok().map(|v| v != 0)));
                (None, fl, it, bw)
            }
        } else {
            (None, None, None, None)
        };
        let star_owned = starname.map(str::to_owned).or(pos_star).or(body_star);
        let fl_val = iflag.or(flags).or(pos_flag).unwrap_or(258);
        let it_val = ifltype.or(pos_type).unwrap_or(0);
        let bw_val = backward.or(backwards).or(pos_back).unwrap_or(false);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || {
            adapter.lun_occult_when_glob_ex(
                jd_ut,
                body_int,
                star_owned.as_deref(),
                fl_val,
                it_val,
                bw_val,
            )
        })
        .map_err(ephemeris_err_to_pyerr)
    }

    /// Searches for the next locally visible lunar occultation of a planet or star at `geopos = [lon, lat, alt]`.
    /// Supports `pyswisseph` 5-arg `(tjdut, body_or_star, geopos, flags=2, backwards=False)`,
    /// C-style 6-arg `(jd_ut, ipl, starname, iflag, geopos, backward=False)`, and
    /// `(jd_ut, ipl, geopos, starname=None, iflag=258, backward=False)`.
    #[pyo3(signature = (jd_ut, ipl, geopos, arg4=None, arg5=None, arg6=None, *, starname=None, iflag=None, flags=None, backward=None, backwards=None))]
    #[allow(clippy::too_many_arguments, clippy::indexing_slicing)]
    pub fn lun_occult_when_loc(
        &self,
        py: Python<'_>,
        jd_ut: f64,
        ipl: &Bound<'_, PyAny>,
        geopos: &Bound<'_, PyAny>,
        arg4: Option<&Bound<'_, PyAny>>,
        arg5: Option<&Bound<'_, PyAny>>,
        arg6: Option<&Bound<'_, PyAny>>,
        starname: Option<&str>,
        iflag: Option<u32>,
        flags: Option<u32>,
        backward: Option<bool>,
        backwards: Option<bool>,
    ) -> PyResult<(i32, [f64; 10], [f64; 20])> {
        let extract_gp = |obj: &Bound<'_, PyAny>| -> PyResult<[f64; 3]> {
            let v: Vec<f64> = obj.extract().map_err(|_| {
                ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: "geopos must be a sequence [lon, lat, alt] or [lon, lat]".to_owned(),
                })
            })?;
            if v.len() < 2 {
                return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: format!("geopos must have at least 2 elements [lon, lat], got {}", v.len()),
                }));
            }
            Ok([v[0], v[1], v.get(2).copied().unwrap_or(0.0)])
        };
        let (body_int, body_star) = if let Ok(id) = ipl.extract::<i32>() {
            (id, None)
        } else if let Ok(s) = ipl.extract::<String>() {
            (0, Some(s))
        } else {
            return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                message: "body/ipl must be an integer SE_* body ID or a star name str".to_owned(),
            }));
        };
        let (gp_arr, pos_star, pos_flag, pos_back) = if !geopos.is_instance_of::<pyo3::types::PyString>()
            && !geopos.is_none()
            && extract_gp(geopos).is_ok()
        {
            let gp = extract_gp(geopos)?;
            if let Some(a4) = arg4 {
                if a4.is_none() || a4.extract::<String>().is_ok() {
                    // (jd_ut, ipl, geopos, starname, iflag, backward)
                    let st = a4.extract::<Option<String>>().ok().flatten();
                    let fl = arg5.and_then(|a| a.extract::<u32>().ok());
                    let bw = arg6.and_then(|a| a.extract::<bool>().ok().or_else(|| a.extract::<i32>().ok().map(|v| v != 0)));
                    (gp, st, fl, bw)
                } else {
                    // pyswisseph 5-arg: (jd_ut, body_or_star, geopos, flags, backwards)
                    let fl = a4.extract::<u32>().ok();
                    let bw = arg5.and_then(|a| a.extract::<bool>().ok().or_else(|| a.extract::<i32>().ok().map(|v| v != 0)));
                    (gp, None, fl, bw)
                }
            } else {
                (gp, None, None, None)
            }
        } else {
            let st = geopos.extract::<Option<String>>().map_err(|_| {
                ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: "Expected geopos sequence or starname str as 3rd argument in lun_occult_when_loc".to_owned(),
                })
            })?;
            if let Some(a4) = arg4.filter(|a| !a.is_instance_of::<pyo3::types::PyString>() && !a.is_none() && extract_gp(a).is_ok()) {
                // (jd_ut, ipl, starname, geopos, iflag, backward)
                let gp = extract_gp(a4)?;
                let fl = arg5.and_then(|a| a.extract::<u32>().ok());
                let bw = arg6.and_then(|a| a.extract::<bool>().ok().or_else(|| a.extract::<i32>().ok().map(|v| v != 0)));
                (gp, st, fl, bw)
            } else {
                // C-style 6-arg: (jd_ut, ipl, starname, iflag, geopos, backward)
                let fl = arg4.and_then(|a| a.extract::<u32>().ok());
                let gp_obj = arg5.ok_or_else(|| {
                    ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                        message: "Missing geopos in lun_occult_when_loc".to_owned(),
                    })
                })?;
                let gp = extract_gp(gp_obj)?;
                let bw = arg6.and_then(|a| a.extract::<bool>().ok().or_else(|| a.extract::<i32>().ok().map(|v| v != 0)));
                (gp, st, fl, bw)
            }
        };
        let star_owned = starname.map(str::to_owned).or(pos_star).or(body_star);
        let fl_val = iflag.or(flags).or(pos_flag).unwrap_or(258);
        let bw_val = backward.or(backwards).or(pos_back).unwrap_or(false);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || {
            adapter.lun_occult_when_loc(jd_ut, body_int, star_owned.as_deref(), fl_val, gp_arr, bw_val)
        })
        .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates `calc_ut(jd_ut, ipl, iflag=258) -> ((lon, lat, dist, lon_speed, lat_speed, dist_speed), retflag)`.
    #[pyo3(signature = (jd_ut=None, ipl=None, iflag=None, *, tjdut=None, flags=None))]
    pub fn calc_ut(
        &self,
        py: Python<'_>,
        jd_ut: Option<f64>,
        ipl: Option<i32>,
        iflag: Option<u32>,
        tjdut: Option<f64>,
        flags: Option<u32>,
    ) -> PyResult<([f64; 6], i32)> {
        let j = jd_ut
            .or(tjdut)
            .ok_or_else(|| PyValueError::new_err("calc_ut missing required argument 'jd_ut' / 'tjdut'"))?;
        let ip = ipl.ok_or_else(|| PyValueError::new_err("calc_ut missing required argument 'ipl'"))?;
        let fl = iflag.or(flags).unwrap_or(258);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || adapter.calc_ut(j, ip, fl))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates `calc(jd_et, ipl, iflag=258) -> ((lon, lat, dist, lon_speed, lat_speed, dist_speed), retflag)`.
    #[pyo3(signature = (jd_et=None, ipl=None, iflag=None, *, tjdet=None, flags=None))]
    pub fn calc(
        &self,
        py: Python<'_>,
        jd_et: Option<f64>,
        ipl: Option<i32>,
        iflag: Option<u32>,
        tjdet: Option<f64>,
        flags: Option<u32>,
    ) -> PyResult<([f64; 6], i32)> {
        let j = jd_et
            .or(tjdet)
            .ok_or_else(|| PyValueError::new_err("calc missing required argument 'jd_et' / 'tjdet'"))?;
        let ip = ipl.ok_or_else(|| PyValueError::new_err("calc missing required argument 'ipl'"))?;
        let fl = iflag.or(flags).unwrap_or(258);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || adapter.calc(j, ip, fl))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates a fixed star at `jd_ut` (`UT1`), returning `(xx, star_name, retflag)`.
    #[pyo3(signature = (star, jd_ut=None, iflag=None, *, tjdut=None, flags=None))]
    pub fn fixstar2_ut(
        &self,
        py: Python<'_>,
        star: &str,
        jd_ut: Option<f64>,
        iflag: Option<u32>,
        tjdut: Option<f64>,
        flags: Option<u32>,
    ) -> PyResult<([f64; 6], String, i32)> {
        let j = jd_ut.or(tjdut).ok_or_else(|| {
            PyValueError::new_err("fixstar2_ut missing required argument 'jd_ut' / 'tjdut'")
        })?;
        let fl = iflag.or(flags).unwrap_or(258);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        let star_s = star.to_owned();
        py.allow_threads(move || adapter.fixstar2_ut(&star_s, j, fl))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Returns the visual magnitude (`Vmag`) of a fixed star.
    pub fn fixstar2_mag(&self, star: &str) -> PyResult<f64> {
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        adapter.fixstar2_mag(star).map_err(ephemeris_err_to_pyerr)
    }

    /// Computes the 3D fractional mundane house position $H \in [1.0, 13.0)$ of `xpin = [lon, lat]`.
    /// Accepts both `(armc, geolat, eps, xpin, hsys='P')` (`pyswisseph` order) and
    /// `(armc, geolat, eps, hsys, xpin)` (`swe_house_pos` C-API order).
    #[pyo3(signature = (armc_deg, geolat_deg, eps_deg, xpin, hsys=None))]
    #[allow(clippy::indexing_slicing)]
    pub fn house_pos(
        &self,
        armc_deg: f64,
        geolat_deg: f64,
        eps_deg: f64,
        xpin: &Bound<'_, PyAny>,
        hsys: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<f64> {
        let (coords, ch) = if let Ok(seq) = xpin.extract::<Vec<f64>>() {
            if seq.len() < 2 {
                return Err(PyValueError::new_err(format!(
                    "xpin must have at least 2 elements [lon, lat], got {}",
                    seq.len()
                )));
            }
            ([seq[0], seq[1]], extract_hsys_char(hsys)?)
        } else {
            let ch = extract_hsys_char(Some(xpin))?;
            let hsys_obj = hsys.ok_or_else(|| {
                PyValueError::new_err("house_pos missing required argument 'xpin' ([lon, lat])")
            })?;
            let seq: Vec<f64> = hsys_obj.extract().map_err(|_| {
                PyValueError::new_err("xpin must be a sequence of [lon, lat] floats")
            })?;
            if seq.len() < 2 {
                return Err(PyValueError::new_err(format!(
                    "xpin must have at least 2 elements [lon, lat], got {}",
                    seq.len()
                )));
            }
            ([seq[0], seq[1]], ch)
        };
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        adapter
            .house_pos(armc_deg, geolat_deg, eps_deg, ch, coords)
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates 12 house cusps and 8 chart angles (`(cusps_12, ascmc_8)`), accepting both
    /// `lxephemeris` (`jd_ut`, `lat_deg`, `lon_deg`, `hsys`, `iflag`) and `pyswisseph`
    /// (`tjdut`, `lat`, `lon`, `hsys`, `flags`) keyword names.
    #[pyo3(signature = (jd_ut=None, lat_deg=None, lon_deg=None, hsys=None, iflag=None, *, tjdut=None, lat=None, lon=None, flags=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn houses_ex(
        &self,
        py: Python<'_>,
        jd_ut: Option<f64>,
        lat_deg: Option<f64>,
        lon_deg: Option<f64>,
        hsys: Option<&Bound<'_, PyAny>>,
        iflag: Option<u32>,
        tjdut: Option<f64>,
        lat: Option<f64>,
        lon: Option<f64>,
        flags: Option<u32>,
    ) -> PyResult<([f64; 12], [f64; 8])> {
        let j = jd_ut.or(tjdut).ok_or_else(|| {
            PyValueError::new_err("houses_ex missing required argument 'tjdut' / 'jd_ut'")
        })?;
        let la = lat_deg.or(lat).ok_or_else(|| {
            PyValueError::new_err("houses_ex missing required argument 'lat' / 'lat_deg'")
        })?;
        let lo = lon_deg.or(lon).ok_or_else(|| {
            PyValueError::new_err("houses_ex missing required argument 'lon' / 'lon_deg'")
        })?;
        let fl = iflag.or(flags).unwrap_or(0);
        let ch = extract_hsys_char(hsys)?;
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || {
            let (cusps, ascmc) = adapter.houses_ex(j, la, lo, ch, fl)?;
            if la.abs() > 60.0 && (fl & SEFLG_RADIANS) == 0 {
                let (ecl_nut, _) = adapter.calc_ut(
                    j,
                    SE_ECL_NUT,
                    if (fl & SEFLG_NONUT) != 0 {
                        SEFLG_NONUT
                    } else {
                        0
                    },
                )?;
                let eps_deg = if (fl & SEFLG_NONUT) != 0 {
                    ecl_nut[1]
                } else {
                    ecl_nut[0]
                };
                let sid_off = if (fl & SEFLG_SIDEREAL) != 0 {
                    adapter.get_ayanamsa_ex_ut(j, fl)?.1
                } else {
                    0.0
                };
                Ok(fix_polar_houses(ch, cusps, ascmc, la, eps_deg, sid_off))
            } else {
                Ok((cusps, ascmc))
            }
        })
        .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates tropical 12 house cusps and 8 chart angles (`(cusps_12, ascmc_8)`).
    #[pyo3(signature = (jd_ut=None, lat_deg=None, lon_deg=None, hsys=None, *, tjdut=None, lat=None, lon=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn houses(
        &self,
        py: Python<'_>,
        jd_ut: Option<f64>,
        lat_deg: Option<f64>,
        lon_deg: Option<f64>,
        hsys: Option<&Bound<'_, PyAny>>,
        tjdut: Option<f64>,
        lat: Option<f64>,
        lon: Option<f64>,
    ) -> PyResult<([f64; 12], [f64; 8])> {
        self.houses_ex(
            py,
            jd_ut,
            lat_deg,
            lon_deg,
            hsys,
            Some(0),
            tjdut,
            lat,
            lon,
            None,
        )
    }

    /// Evaluates 12 house cusps, 8 chart angles, and their speeds (`(cusps_12, ascmc_8, cusp_speed_12, ascmc_speed_8)`).
    #[pyo3(signature = (jd_ut=None, lat_deg=None, lon_deg=None, hsys=None, iflag=None, *, tjdut=None, lat=None, lon=None, flags=None))]
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn houses_ex2(
        &self,
        py: Python<'_>,
        jd_ut: Option<f64>,
        lat_deg: Option<f64>,
        lon_deg: Option<f64>,
        hsys: Option<&Bound<'_, PyAny>>,
        iflag: Option<u32>,
        tjdut: Option<f64>,
        lat: Option<f64>,
        lon: Option<f64>,
        flags: Option<u32>,
    ) -> PyResult<([f64; 12], [f64; 8], [f64; 12], [f64; 8])> {
        let j = jd_ut.or(tjdut).ok_or_else(|| {
            PyValueError::new_err("houses_ex2 missing required argument 'tjdut' / 'jd_ut'")
        })?;
        let la = lat_deg.or(lat).ok_or_else(|| {
            PyValueError::new_err("houses_ex2 missing required argument 'lat' / 'lat_deg'")
        })?;
        let lo = lon_deg.or(lon).ok_or_else(|| {
            PyValueError::new_err("houses_ex2 missing required argument 'lon' / 'lon_deg'")
        })?;
        let fl = iflag.or(flags).unwrap_or(0);
        let ch = extract_hsys_char(hsys)?;
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || {
            let (_, _, cusp_speed, ascmc_speed) = adapter.houses_ex2(j, la, lo, ch, fl)?;
            let (cusps, ascmc) = adapter.houses_ex(j, la, lo, ch, fl)?;
            if la.abs() > 60.0 && (fl & SEFLG_RADIANS) == 0 {
                let (ecl_nut, _) = adapter.calc_ut(
                    j,
                    SE_ECL_NUT,
                    if (fl & SEFLG_NONUT) != 0 {
                        SEFLG_NONUT
                    } else {
                        0
                    },
                )?;
                let eps_deg = if (fl & SEFLG_NONUT) != 0 {
                    ecl_nut[1]
                } else {
                    ecl_nut[0]
                };
                let sid_off = if (fl & SEFLG_SIDEREAL) != 0 {
                    adapter.get_ayanamsa_ex_ut(j, fl)?.1
                } else {
                    0.0
                };
                let (fc, fa) = fix_polar_houses(ch, cusps, ascmc, la, eps_deg, sid_off);
                Ok((fc, fa, cusp_speed, ascmc_speed))
            } else {
                Ok((cusps, ascmc, cusp_speed, ascmc_speed))
            }
        })
        .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates 12 house cusps and 8 chart angles directly from `(armc_deg, lat_deg, eps_deg, hsys)`.
    #[pyo3(signature = (armc_deg=None, lat_deg=None, eps_deg=None, hsys=None, ascmc9=None, *, armc=None, lat=None, eps=None, sun_decl_deg=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn houses_armc(
        &self,
        py: Python<'_>,
        armc_deg: Option<f64>,
        lat_deg: Option<f64>,
        eps_deg: Option<f64>,
        hsys: Option<&Bound<'_, PyAny>>,
        ascmc9: Option<f64>,
        armc: Option<f64>,
        lat: Option<f64>,
        eps: Option<f64>,
        sun_decl_deg: Option<f64>,
    ) -> PyResult<([f64; 12], [f64; 8])> {
        let ar = armc_deg.or(armc).ok_or_else(|| {
            PyValueError::new_err("houses_armc missing required argument 'armc' / 'armc_deg'")
        })?;
        let la = lat_deg.or(lat).ok_or_else(|| {
            PyValueError::new_err("houses_armc missing required argument 'lat' / 'lat_deg'")
        })?;
        let ep = eps_deg.or(eps).ok_or_else(|| {
            PyValueError::new_err("houses_armc missing required argument 'eps' / 'eps_deg'")
        })?;
        let ch = extract_hsys_char(hsys)?;
        let decl_opt = ascmc9.or(sun_decl_deg);
        let effective_ch = if decl_opt.is_none() && ch == 'I' {
            'O'
        } else {
            ch
        };
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || {
            let (cusps, ascmc) = adapter.houses_armc_ex(ar, la, ep, effective_ch, decl_opt)?;
            if la.abs() > 60.0 {
                Ok(fix_polar_houses(effective_ch, cusps, ascmc, la, ep, 0.0))
            } else if effective_ch == 'W' && ((ascmc[0].rem_euclid(30.0)) - 30.0).abs() < 1e-11 {
                Ok((reconstruct_whole_sign(ascmc[0]), ascmc))
            } else {
                Ok((cusps, ascmc))
            }
        })
        .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates horizontal coordinates (`(azimuth_deg, true_alt_deg, apparent_alt_deg)`).
    /// Accepts 2-element `[lon_or_ra, lat_or_dec]` or 3-element `[lon_or_ra, lat_or_dec, dist]` sequences for `xin`.
    #[allow(clippy::too_many_arguments, clippy::indexing_slicing)]
    pub fn azalt(
        &self,
        py: Python<'_>,
        jd_ut: f64,
        calc_flag: i32,
        geopos: [f64; 3],
        atpress: f64,
        attemp: f64,
        xin: Vec<f64>,
    ) -> PyResult<[f64; 3]> {
        if xin.len() < 2 {
            return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                message: format!(
                    "xin must have at least 2 elements [lon_or_ra, lat_or_dec], got {}",
                    xin.len()
                ),
            }));
        }
        let xin2 = [xin[0], xin[1]];
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        let hz = py
            .allow_threads(move || adapter.azalt(jd_ut, calc_flag, geopos, atpress, attemp, xin2))
            .map_err(ephemeris_err_to_pyerr)?;
        Ok([
            hz.azimuth_deg,
            hz.altitude_true_deg,
            hz.altitude_apparent_deg,
        ])
    }

    /// Evaluates the Equation of Time in days at `jd_ut` (`swe_time_equ`).
    pub fn time_equ(&self, jd_ut: f64) -> PyResult<f64> {
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        adapter.time_equ(jd_ut).map_err(ephemeris_err_to_pyerr)
    }

    /// Inverse horizontal transformation (`swe_azalt_rev`): accepts either 4 arguments
    /// `(jd_ut, calc_flag, geopos, [az, alt])` or 5 arguments `(jd_ut, calc_flag, geopos, az, alt)`.
    #[pyo3(signature = (jd_ut, calc_flag, geopos, xin, alt=None))]
    #[allow(clippy::indexing_slicing)]
    pub fn azalt_rev(
        &self,
        jd_ut: f64,
        calc_flag: i32,
        geopos: [f64; 3],
        xin: &Bound<'_, PyAny>,
        alt: Option<f64>,
    ) -> PyResult<[f64; 2]> {
        let xin2 = if let Some(alt_val) = alt {
            let az_val: f64 = xin.extract().map_err(|_| {
                ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: "When alt is provided as the 5th argument, the 4th argument must be azimuth (float)".to_owned(),
                })
            })?;
            [az_val, alt_val]
        } else {
            let vec: Vec<f64> = xin.extract().map_err(|_| {
                ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: "xin must be a sequence [azimuth, altitude] or pass (az, alt) as separate arguments".to_owned(),
                })
            })?;
            if vec.len() < 2 {
                return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: format!(
                        "xin must have at least 2 elements [azimuth, altitude], got {}",
                        vec.len()
                    ),
                }));
            }
            [vec[0], vec[1]]
        };
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        adapter
            .azalt_rev(jd_ut, calc_flag, geopos, xin2)
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Searches for the next rise, set, or meridian transit (`JD_UT1`), returning `(res_code, tret_10_tuple)`
    /// matching `pyswisseph.rise_trans` (`(0, (t_rise, 0.0, ...))` or `(-2, (0.0, ...))` when circumpolar).
    /// Supports both `pyswisseph` `(jd_ut, body_or_star, rsmi, geopos, atpress=0.0, attemp=15.0, flags=2)`
    /// and C-style `(jd_ut, ipl, starname, epheflag, rsmi, geopos, atpress=0.0, attemp=15.0)` signatures.
    #[pyo3(signature = (jd_ut, ipl, rsmi, geopos, atpress=None, attemp=None, flags=None, starname=None, epheflag=None))]
    #[allow(clippy::too_many_arguments, clippy::type_complexity, clippy::indexing_slicing)]
    pub fn rise_trans(
        &self,
        py: Python<'_>,
        jd_ut: f64,
        ipl: &Bound<'_, PyAny>,
        rsmi: &Bound<'_, PyAny>,
        geopos: &Bound<'_, PyAny>,
        atpress: Option<&Bound<'_, PyAny>>,
        attemp: Option<&Bound<'_, PyAny>>,
        flags: Option<&Bound<'_, PyAny>>,
        starname: Option<&Bound<'_, PyAny>>,
        epheflag: Option<u32>,
    ) -> PyResult<(i32, (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64))> {
        let extract_geopos = |obj: &Bound<'_, PyAny>| -> PyResult<[f64; 3]> {
            let v: Vec<f64> = obj.extract().map_err(|_| {
                ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: "geopos must be a sequence [lon, lat, alt] or [lon, lat]".to_owned(),
                })
            })?;
            if v.len() < 2 {
                return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                    message: format!("geopos must have at least 2 elements [lon, lat], got {}", v.len()),
                }));
            }
            Ok([v[0], v[1], v.get(2).copied().unwrap_or(0.0)])
        };

        let (ipl_int, star_opt, ephe_flags, rsmi_int, geopos_arr, press_val, temp_val) =
            if let Ok(gp) = extract_geopos(geopos) {
                // pyswisseph signature: (jd_ut, ipl_or_star, rsmi, geopos, atpress=0.0, attemp=15.0, flags=2, starname=None, epheflag=None)
                let (body_int, body_star) = if let Ok(id) = ipl.extract::<i32>() {
                    (id, None)
                } else if let Ok(s) = ipl.extract::<String>() {
                    (0, Some(s))
                } else {
                    return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                        message: "ipl must be an integer SE_* body ID or a star name str".to_owned(),
                    }));
                };
                let kw_star = match starname {
                    Some(st) => st.extract::<Option<String>>().map_err(|_| {
                        ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                            message: "starname must be a str or None".to_owned(),
                        })
                    })?,
                    None => None,
                };
                let resolved_star = kw_star.or(body_star);
                let rsmi_val: i32 = rsmi.extract().map_err(|_| {
                    ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                        message: "rsmi must be an integer transit flag".to_owned(),
                    })
                })?;
                let p_val = match atpress {
                    Some(a) => a.extract::<Option<f64>>().ok().flatten().unwrap_or(0.0),
                    None => 0.0,
                };
                let t_val = match attemp {
                    Some(a) => a.extract::<Option<f64>>().ok().flatten().unwrap_or(15.0),
                    None => 15.0,
                };
                let f_val = if let Some(ef) = epheflag {
                    ef
                } else if let Some(fl) = flags {
                    fl.extract::<Option<u32>>().ok().flatten().unwrap_or(2)
                } else {
                    2
                };
                (body_int, resolved_star, f_val, rsmi_val, gp, p_val, t_val)
            } else {
                // C-style 8-arg signature: (jd_ut, ipl, starname, epheflag, rsmi, geopos, atpress=0.0, attemp=15.0)
                let body_int = ipl.extract::<i32>().unwrap_or(0);
                let c_star = rsmi.extract::<Option<String>>().map_err(|_| {
                    ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                        message: "Expected starname (str or None) as 3rd argument in C-style rise_trans".to_owned(),
                    })
                })?;
                let f_val: u32 = geopos.extract().map_err(|_| {
                    ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                        message: "Expected geopos sequence as 4th argument (pyswisseph) or epheflag (int) in C-style rise_trans".to_owned(),
                    })
                })?;
                let rsmi_val: i32 = atpress
                    .ok_or_else(|| {
                        ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                            message: "Missing rsmi in C-style rise_trans".to_owned(),
                        })
                    })?
                    .extract()
                    .map_err(|_| {
                        ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                            message: "Expected rsmi (int) as 5th argument in C-style rise_trans".to_owned(),
                        })
                    })?;
                let gp_obj = attemp.ok_or_else(|| {
                    ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                        message: "Missing geopos in C-style rise_trans".to_owned(),
                    })
                })?;
                let gp = extract_geopos(gp_obj)?;
                let p_val = match flags {
                    Some(a) => a.extract::<Option<f64>>().ok().flatten().unwrap_or(0.0),
                    None => 0.0,
                };
                let t_val = match starname {
                    Some(a) => a.extract::<Option<f64>>().ok().flatten().unwrap_or(15.0),
                    None => 15.0,
                };
                (body_int, c_star, f_val, rsmi_val, gp, p_val, t_val)
            };

        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        let res = py
            .allow_threads(move || {
                adapter.rise_trans(
                    jd_ut,
                    ipl_int,
                    star_opt.as_deref(),
                    ephe_flags,
                    rsmi_int,
                    geopos_arr,
                    press_val,
                    temp_val,
                )
            })
            .map_err(ephemeris_err_to_pyerr)?;
        match res {
            Some(t) => Ok((0, (t, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))),
            None => Ok((-2, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))),
        }
    }

    /// Sets the ephemeris path (`swe_set_ephe_path` compatibility stub).
    #[pyo3(signature = (path=None))]
    pub fn set_ephe_path(&self, path: Option<&str>) -> PyResult<()> {
        let mut guard = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        guard.set_ephe_path(path.unwrap_or(""));
        Ok(())
    }

    /// Closes the adapter (`swe_close` no-op compatibility stub).
    pub fn close(&self) -> PyResult<()> {
        let mut guard = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        guard.close();
        Ok(())
    }

    /// Returns the canonical Swiss Ephemeris ayanamsa name for `sid_mode` (`swe_get_ayanamsa_name`).
    #[staticmethod]
    pub fn get_ayanamsa_name(sid_mode: i32) -> &'static str {
        SweAdapter::get_ayanamsa_name(sid_mode)
    }

    /// Returns the canonical display name of an `SE_*` planet, node, apogee, asteroid, or hypothetical body (`swe_get_planet_name`).
    #[staticmethod]
    pub fn get_planet_name(ipl: i32) -> String {
        SweAdapter::get_planet_name(ipl)
    }

    /// Returns the canonical display name of a 1-character house system code (`swe_house_name`).
    #[staticmethod]
    #[pyo3(signature = (hsys=None))]
    pub fn house_name(hsys: Option<&Bound<'_, PyAny>>) -> PyResult<&'static str> {
        let ch = extract_raw_hsys_char(hsys)?;
        Ok(SweAdapter::house_name(ch))
    }

    /// Signed angular difference `p1 - p2` in degrees normalized to `(-180, 180]` (`swe_difdeg2n`).
    #[staticmethod]
    pub fn difdeg2n(p1: f64, p2: f64) -> f64 {
        SweAdapter::difdeg2n(p1, p2)
    }

    /// Normalizes an angle in degrees to `[0, 360)` (`swe_degnorm`).
    #[staticmethod]
    pub fn degnorm(x: f64) -> f64 {
        SweAdapter::degnorm(x)
    }

    /// Coordinate transformation between ecliptic and equatorial spherical coordinates (`swe_cotrans`).
    /// Accepts either `(xpo_seq, eps)` (2 arguments) or `(lon, lat, dist, eps)` (4 arguments).
    #[staticmethod]
    #[pyo3(signature = (lon, lat, dist=None, eps=None))]
    pub fn cotrans(
        lon: &Bound<'_, PyAny>,
        lat: f64,
        dist: Option<f64>,
        eps: Option<f64>,
    ) -> PyResult<(f64, f64, f64)> {
        py_cotrans_impl(lon, lat, dist, eps)
    }

    /// Searches for the next (or previous if `backward=True`) global solar eclipse (`swe_sol_eclipse_when_glob`).
    #[pyo3(signature = (jd_start, iflag=2, ifltype=0, backward=false, flags=None, backwards=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn sol_eclipse_when_glob(
        &self,
        py: Python<'_>,
        jd_start: f64,
        iflag: u32,
        ifltype: i32,
        backward: bool,
        flags: Option<u32>,
        backwards: Option<bool>,
    ) -> PyResult<(i32, [f64; 10])> {
        let actual_flags = flags.unwrap_or(iflag);
        let actual_backward = backwards.unwrap_or(backward);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || {
            adapter.sol_eclipse_when_glob(jd_start, actual_flags, ifltype, actual_backward)
        })
        .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates geographic central/greatest-eclipse coordinates and attributes at `jd_ut` (`swe_sol_eclipse_where`).
    #[pyo3(signature = (jd_ut, iflag=2, flags=None))]
    pub fn sol_eclipse_where(
        &self,
        py: Python<'_>,
        jd_ut: f64,
        iflag: u32,
        flags: Option<u32>,
    ) -> PyResult<(i32, [f64; 10], [f64; 20])> {
        let actual_flags = flags.unwrap_or(iflag);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || adapter.sol_eclipse_where(jd_ut, actual_flags))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Searches for the next (or previous if `backward=True`) locally visible solar eclipse at `geopos = [lon, lat, alt]` (`swe_sol_eclipse_when_loc`).
    #[pyo3(signature = (jd_start, geopos, iflag=2, backward=false, flags=None, backwards=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn sol_eclipse_when_loc(
        &self,
        py: Python<'_>,
        jd_start: f64,
        geopos: [f64; 3],
        iflag: u32,
        backward: bool,
        flags: Option<u32>,
        backwards: Option<bool>,
    ) -> PyResult<(i32, [f64; 10], [f64; 20])> {
        let actual_flags = flags.unwrap_or(iflag);
        let actual_backward = backwards.unwrap_or(backward);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || {
            adapter.sol_eclipse_when_loc(jd_start, geopos, actual_flags, actual_backward)
        })
        .map_err(ephemeris_err_to_pyerr)
    }

    /// Searches for the next (or previous if `backward=True`) lunar eclipse (`swe_lun_eclipse_when`).
    #[pyo3(signature = (jd_start, iflag=2, ifltype=0, backward=false, flags=None, backwards=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn lun_eclipse_when(
        &self,
        py: Python<'_>,
        jd_start: f64,
        iflag: u32,
        ifltype: i32,
        backward: bool,
        flags: Option<u32>,
        backwards: Option<bool>,
    ) -> PyResult<(i32, [f64; 10])> {
        let actual_flags = flags.unwrap_or(iflag);
        let actual_backward = backwards.unwrap_or(backward);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || {
            adapter.lun_eclipse_when(jd_start, actual_flags, ifltype, actual_backward)
        })
        .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates lunar eclipse attributes at `jd_ut` for `geopos = [lon, lat, alt]` (`swe_lun_eclipse_how`).
    #[pyo3(signature = (jd_ut, geopos, iflag=2, flags=None))]
    pub fn lun_eclipse_how(
        &self,
        py: Python<'_>,
        jd_ut: f64,
        geopos: [f64; 3],
        iflag: u32,
        flags: Option<u32>,
    ) -> PyResult<(i32, [f64; 20])> {
        let actual_flags = flags.unwrap_or(iflag);
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        py.allow_threads(move || adapter.lun_eclipse_how(jd_ut, geopos, actual_flags))
            .map_err(ephemeris_err_to_pyerr)
    }

    /// Evaluates the next heliacal visibility event (`swe_heliacal_ut`) starting from `jd_start`,
    /// returning the 3-tuple `(t_start, t_opt, t_end)` matching `pyswisseph.heliacal_ut`.
    #[pyo3(signature = (jd_start, geopos, datm, dobs, object_name, event_type, helflag=0))]
    #[allow(clippy::too_many_arguments, clippy::indexing_slicing)]
    pub fn heliacal_ut(
        &self,
        py: Python<'_>,
        jd_start: f64,
        geopos: Vec<f64>,
        datm: Vec<f64>,
        dobs: Vec<f64>,
        object_name: &Bound<'_, PyAny>,
        event_type: i32,
        helflag: i32,
    ) -> PyResult<(f64, f64, f64)> {
        if geopos.len() < 2 {
            return Err(ephemeris_err_to_pyerr(EphemerisError::InvalidInput {
                message: format!(
                    "geopos must have at least 2 elements [lon, lat], got {}",
                    geopos.len()
                ),
            }));
        }
        let obj_owned: String = if let Ok(s) = object_name.extract::<String>() {
            s
        } else if let Ok(ipl) = object_name.extract::<i32>() {
            let name = SweAdapter::get_planet_name(ipl);
            if name.starts_with("unknown_") {
                return Err(UnsupportedBodyError::new_err(format!(
                    "[UnsupportedBody] Unsupported body ID {ipl}"
                )));
            }
            name
        } else {
            return Err(InvalidInputError::new_err(
                "[InvalidInput] object_name must be a str or int body ID",
            ));
        };
        let geopos_arr = [
            geopos[0],
            geopos[1],
            geopos.get(2).copied().unwrap_or(0.0),
        ];
        let datm_arr = [
            datm.first().copied().unwrap_or(0.0),
            datm.get(1).copied().unwrap_or(0.0),
            datm.get(2).copied().unwrap_or(0.0),
            datm.get(3).copied().unwrap_or(0.0),
        ];
        let dobs_arr = [
            dobs.first().copied().unwrap_or(0.0),
            dobs.get(1).copied().unwrap_or(0.0),
            dobs.get(2).copied().unwrap_or(0.0),
            dobs.get(3).copied().unwrap_or(0.0),
            dobs.get(4).copied().unwrap_or(0.0),
            dobs.get(5).copied().unwrap_or(0.0),
        ];
        let adapter = self
            .inner
            .lock()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .clone();
        let (tret, _retflag) = py
            .allow_threads(move || {
                adapter.heliacal_ut(
                    jd_start,
                    geopos_arr,
                    datm_arr,
                    dobs_arr,
                    &obj_owned,
                    event_type,
                    helflag,
                )
            })
            .map_err(ephemeris_err_to_pyerr)?;
        Ok((tret[0], tret[1], tret[2]))
    }

    #[staticmethod]
    #[pyo3(signature = (year, month, day, hour=0.0, gregflag=None))]
    pub fn julday(
        year: i32,
        month: u8,
        day: u8,
        hour: f64,
        gregflag: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<f64> {
        let cal = extract_calendar_kind(gregflag)?;
        core_julday(year, month, day, hour, cal).map_err(ephemeris_err_to_pyerr)
    }

    #[staticmethod]
    #[pyo3(signature = (jd, gregflag=None))]
    pub fn revjul(jd: f64, gregflag: Option<&Bound<'_, PyAny>>) -> PyResult<(i32, u8, u8, f64)> {
        let cal = extract_calendar_kind(gregflag)?;
        let d = core_revjul(jd, cal).map_err(ephemeris_err_to_pyerr)?;
        Ok((d.year, d.month, d.day, d.hour))
    }

    // Version & Exception attributes on SweAdapter for `swe = SweAdapter()` module-replacement compatibility
    #[classattr]
    #[allow(non_upper_case_globals)]
    pub const version: &'static str = env!("CARGO_PKG_VERSION");

    #[classattr]
    #[allow(non_upper_case_globals)]
    pub const __version__: &'static str = env!("CARGO_PKG_VERSION");

    #[classattr]
    #[allow(non_snake_case)]
    pub fn Error(py: Python<'_>) -> Bound<'_, PyType> {
        py.get_type::<LxEphemerisError>()
    }

    #[classattr]
    #[allow(non_snake_case)]
    pub fn LxEphemerisError(py: Python<'_>) -> Bound<'_, PyType> {
        py.get_type::<LxEphemerisError>()
    }

    // SE_* and unprefixed body constants on SweAdapter
    #[classattr] pub const SE_ECL_NUT: i32 = SE_ECL_NUT;
    #[classattr] pub const ECL_NUT: i32 = SE_ECL_NUT;
    #[classattr] pub const SE_SUN: i32 = SE_SUN;
    #[classattr] pub const SUN: i32 = SE_SUN;
    #[classattr] pub const SE_MOON: i32 = SE_MOON;
    #[classattr] pub const MOON: i32 = SE_MOON;
    #[classattr] pub const SE_MERCURY: i32 = SE_MERCURY;
    #[classattr] pub const MERCURY: i32 = SE_MERCURY;
    #[classattr] pub const SE_VENUS: i32 = SE_VENUS;
    #[classattr] pub const VENUS: i32 = SE_VENUS;
    #[classattr] pub const SE_MARS: i32 = SE_MARS;
    #[classattr] pub const MARS: i32 = SE_MARS;
    #[classattr] pub const SE_JUPITER: i32 = SE_JUPITER;
    #[classattr] pub const JUPITER: i32 = SE_JUPITER;
    #[classattr] pub const SE_SATURN: i32 = SE_SATURN;
    #[classattr] pub const SATURN: i32 = SE_SATURN;
    #[classattr] pub const SE_URANUS: i32 = SE_URANUS;
    #[classattr] pub const URANUS: i32 = SE_URANUS;
    #[classattr] pub const SE_NEPTUNE: i32 = SE_NEPTUNE;
    #[classattr] pub const NEPTUNE: i32 = SE_NEPTUNE;
    #[classattr] pub const SE_PLUTO: i32 = SE_PLUTO;
    #[classattr] pub const PLUTO: i32 = SE_PLUTO;
    #[classattr] pub const SE_MEAN_NODE: i32 = SE_MEAN_NODE;
    #[classattr] pub const MEAN_NODE: i32 = SE_MEAN_NODE;
    #[classattr] pub const SE_TRUE_NODE: i32 = SE_TRUE_NODE;
    #[classattr] pub const TRUE_NODE: i32 = SE_TRUE_NODE;
    #[classattr] pub const SE_MEAN_APOG: i32 = SE_MEAN_APOG;
    #[classattr] pub const MEAN_APOG: i32 = SE_MEAN_APOG;
    #[classattr] pub const SE_OSCU_APOG: i32 = SE_OSCU_APOG;
    #[classattr] pub const OSCU_APOG: i32 = SE_OSCU_APOG;
    #[classattr] pub const SE_EARTH: i32 = SE_EARTH;
    #[classattr] pub const EARTH: i32 = SE_EARTH;
    #[classattr] pub const SE_CHIRON: i32 = SE_CHIRON;
    #[classattr] pub const CHIRON: i32 = SE_CHIRON;
    #[classattr] pub const SE_PHOLUS: i32 = SE_PHOLUS;
    #[classattr] pub const PHOLUS: i32 = SE_PHOLUS;
    #[classattr] pub const SE_CERES: i32 = SE_CERES;
    #[classattr] pub const CERES: i32 = SE_CERES;
    #[classattr] pub const SE_PALLAS: i32 = SE_PALLAS;
    #[classattr] pub const PALLAS: i32 = SE_PALLAS;
    #[classattr] pub const SE_JUNO: i32 = SE_JUNO;
    #[classattr] pub const JUNO: i32 = SE_JUNO;
    #[classattr] pub const SE_VESTA: i32 = SE_VESTA;
    #[classattr] pub const VESTA: i32 = SE_VESTA;
    #[classattr] pub const SE_INTP_APOG: i32 = SE_INTP_APOG;
    #[classattr] pub const INTP_APOG: i32 = SE_INTP_APOG;
    #[classattr] pub const SE_INTP_PERG: i32 = SE_INTP_PERG;
    #[classattr] pub const INTP_PERG: i32 = SE_INTP_PERG;
    #[classattr] pub const SE_FICT_OFFSET: i32 = SE_FICT_OFFSET;
    #[classattr] pub const FICT_OFFSET: i32 = SE_FICT_OFFSET;
    #[classattr] pub const SE_CUPIDO: i32 = SE_CUPIDO;
    #[classattr] pub const CUPIDO: i32 = SE_CUPIDO;
    #[classattr] pub const SE_HADES: i32 = SE_HADES;
    #[classattr] pub const HADES: i32 = SE_HADES;
    #[classattr] pub const SE_ZEUS: i32 = SE_ZEUS;
    #[classattr] pub const ZEUS: i32 = SE_ZEUS;
    #[classattr] pub const SE_KRONOS: i32 = SE_KRONOS;
    #[classattr] pub const KRONOS: i32 = SE_KRONOS;
    #[classattr] pub const SE_APOLLON: i32 = SE_APOLLON;
    #[classattr] pub const APOLLON: i32 = SE_APOLLON;
    #[classattr] pub const SE_ADMETOS: i32 = SE_ADMETOS;
    #[classattr] pub const ADMETOS: i32 = SE_ADMETOS;
    #[classattr] pub const SE_VULKANUS: i32 = SE_VULKANUS;
    #[classattr] pub const VULKANUS: i32 = SE_VULKANUS;
    #[classattr] pub const SE_POSEIDON: i32 = SE_POSEIDON;
    #[classattr] pub const POSEIDON: i32 = SE_POSEIDON;
    #[classattr] pub const SE_ISIS: i32 = SE_ISIS;
    #[classattr] pub const ISIS: i32 = SE_ISIS;
    #[classattr] pub const SE_AST_OFFSET: i32 = SE_AST_OFFSET;
    #[classattr] pub const AST_OFFSET: i32 = SE_AST_OFFSET;

    // SEFLG_* and FLG_* calculation flags on SweAdapter
    #[classattr] pub const SEFLG_JPLEPH: u32 = SEFLG_JPLEPH;
    #[classattr] pub const FLG_JPLEPH: u32 = SEFLG_JPLEPH;
    #[classattr] pub const SEFLG_SWIEPH: u32 = SEFLG_SWIEPH;
    #[classattr] pub const FLG_SWIEPH: u32 = SEFLG_SWIEPH;
    #[classattr] pub const SEFLG_MOSEPH: u32 = SEFLG_MOSEPH;
    #[classattr] pub const FLG_MOSEPH: u32 = SEFLG_MOSEPH;
    #[classattr] pub const SEFLG_HELCTR: u32 = SEFLG_HELCTR;
    #[classattr] pub const FLG_HELCTR: u32 = SEFLG_HELCTR;
    #[classattr] pub const SEFLG_TRUEPOS: u32 = SEFLG_TRUEPOS;
    #[classattr] pub const FLG_TRUEPOS: u32 = SEFLG_TRUEPOS;
    #[classattr] pub const SEFLG_J2000: u32 = SEFLG_J2000;
    #[classattr] pub const FLG_J2000: u32 = SEFLG_J2000;
    #[classattr] pub const SEFLG_NONUT: u32 = SEFLG_NONUT;
    #[classattr] pub const FLG_NONUT: u32 = SEFLG_NONUT;
    #[classattr] pub const SEFLG_SPEED3: u32 = SEFLG_SPEED3;
    #[classattr] pub const FLG_SPEED3: u32 = SEFLG_SPEED3;
    #[classattr] pub const SEFLG_SPEED: u32 = SEFLG_SPEED;
    #[classattr] pub const FLG_SPEED: u32 = SEFLG_SPEED;
    #[classattr] pub const SEFLG_NOGDEFL: u32 = SEFLG_NOGDEFL;
    #[classattr] pub const FLG_NOGDEFL: u32 = SEFLG_NOGDEFL;
    #[classattr] pub const SEFLG_NOABERR: u32 = SEFLG_NOABERR;
    #[classattr] pub const FLG_NOABERR: u32 = SEFLG_NOABERR;
    #[classattr] pub const SEFLG_BARYCTR: u32 = SEFLG_BARYCTR;
    #[classattr] pub const FLG_BARYCTR: u32 = SEFLG_BARYCTR;
    #[classattr] pub const SEFLG_TOPOCTR: u32 = SEFLG_TOPOCTR;
    #[classattr] pub const FLG_TOPOCTR: u32 = SEFLG_TOPOCTR;
    #[classattr] pub const SEFLG_EQUATORIAL: u32 = SEFLG_EQUATORIAL;
    #[classattr] pub const FLG_EQUATORIAL: u32 = SEFLG_EQUATORIAL;
    #[classattr] pub const SEFLG_XYZ: u32 = SEFLG_XYZ;
    #[classattr] pub const FLG_XYZ: u32 = SEFLG_XYZ;
    #[classattr] pub const SEFLG_RADIANS: u32 = SEFLG_RADIANS;
    #[classattr] pub const FLG_RADIANS: u32 = SEFLG_RADIANS;
    #[classattr] pub const SEFLG_SIDEREAL: u32 = SEFLG_SIDEREAL;
    #[classattr] pub const FLG_SIDEREAL: u32 = SEFLG_SIDEREAL;

    // Node bit constants
    #[classattr] pub const SE_NODBIT_MEAN: u32 = SE_NODBIT_MEAN;
    #[classattr] pub const NODBIT_MEAN: u32 = SE_NODBIT_MEAN;
    #[classattr] pub const SE_NODBIT_OSCU: u32 = SE_NODBIT_OSCU;
    #[classattr] pub const NODBIT_OSCU: u32 = SE_NODBIT_OSCU;
    #[classattr] pub const SE_NODBIT_OSCU_BAR: u32 = SE_NODBIT_OSCU_BAR;
    #[classattr] pub const NODBIT_OSCU_BAR: u32 = SE_NODBIT_OSCU_BAR;
    #[classattr] pub const SE_NODBIT_FOPOINT: u32 = SE_NODBIT_FOPOINT;
    #[classattr] pub const NODBIT_FOPOINT: u32 = SE_NODBIT_FOPOINT;

    // Sidereal mode constants
    #[classattr] pub const SE_SIDM_FAGAN_BRADLEY: u32 = SE_SIDM_FAGAN_BRADLEY;
    #[classattr] pub const SIDM_FAGAN_BRADLEY: u32 = SE_SIDM_FAGAN_BRADLEY;
    #[classattr] pub const SE_SIDM_LAHIRI: u32 = SE_SIDM_LAHIRI;
    #[classattr] pub const SIDM_LAHIRI: u32 = SE_SIDM_LAHIRI;
    #[classattr] pub const SE_SIDM_DELUCE: u32 = SE_SIDM_DELUCE;
    #[classattr] pub const SIDM_DELUCE: u32 = SE_SIDM_DELUCE;
    #[classattr] pub const SE_SIDM_RAMAN: u32 = SE_SIDM_RAMAN;
    #[classattr] pub const SIDM_RAMAN: u32 = SE_SIDM_RAMAN;
    #[classattr] pub const SE_SIDM_USHASHASHI: u32 = SE_SIDM_USHASHASHI;
    #[classattr] pub const SIDM_USHASHASHI: u32 = SE_SIDM_USHASHASHI;
    #[classattr] pub const SE_SIDM_KRISHNAMURTI: u32 = SE_SIDM_KRISHNAMURTI;
    #[classattr] pub const SIDM_KRISHNAMURTI: u32 = SE_SIDM_KRISHNAMURTI;
    #[classattr] pub const SE_SIDM_DJWHAL_KHUL: u32 = SE_SIDM_DJWHAL_KHUL;
    #[classattr] pub const SIDM_DJWHAL_KHUL: u32 = SE_SIDM_DJWHAL_KHUL;
    #[classattr] pub const SE_SIDM_YUKTESHWAR: u32 = SE_SIDM_YUKTESHWAR;
    #[classattr] pub const SIDM_YUKTESHWAR: u32 = SE_SIDM_YUKTESHWAR;
    #[classattr] pub const SE_SIDM_JN_BHASIN: u32 = SE_SIDM_JN_BHASIN;
    #[classattr] pub const SIDM_JN_BHASIN: u32 = SE_SIDM_JN_BHASIN;
    #[classattr] pub const SE_SIDM_BABYL_KUGLER1: u32 = SE_SIDM_BABYL_KUGLER1;
    #[classattr] pub const SIDM_BABYL_KUGLER1: u32 = SE_SIDM_BABYL_KUGLER1;
    #[classattr] pub const SE_SIDM_BABYL_KUGLER2: u32 = SE_SIDM_BABYL_KUGLER2;
    #[classattr] pub const SIDM_BABYL_KUGLER2: u32 = SE_SIDM_BABYL_KUGLER2;
    #[classattr] pub const SE_SIDM_BABYL_KUGLER3: u32 = SE_SIDM_BABYL_KUGLER3;
    #[classattr] pub const SIDM_BABYL_KUGLER3: u32 = SE_SIDM_BABYL_KUGLER3;
    #[classattr] pub const SE_SIDM_BABYL_HUBER: u32 = SE_SIDM_BABYL_HUBER;
    #[classattr] pub const SIDM_BABYL_HUBER: u32 = SE_SIDM_BABYL_HUBER;
    #[classattr] pub const SE_SIDM_BABYL_ETPSC: u32 = SE_SIDM_BABYL_ETPSC;
    #[classattr] pub const SIDM_BABYL_ETPSC: u32 = SE_SIDM_BABYL_ETPSC;
    #[classattr] pub const SE_SIDM_ALDEBARAN_15TAU: u32 = SE_SIDM_ALDEBARAN_15TAU;
    #[classattr] pub const SIDM_ALDEBARAN_15TAU: u32 = SE_SIDM_ALDEBARAN_15TAU;
    #[classattr] pub const SE_SIDM_HIPPARCHOS: u32 = SE_SIDM_HIPPARCHOS;
    #[classattr] pub const SIDM_HIPPARCHOS: u32 = SE_SIDM_HIPPARCHOS;
    #[classattr] pub const SE_SIDM_SASSANIAN: u32 = SE_SIDM_SASSANIAN;
    #[classattr] pub const SIDM_SASSANIAN: u32 = SE_SIDM_SASSANIAN;
    #[classattr] pub const SE_SIDM_GALCENT_0SAG: u32 = SE_SIDM_GALCENT_0SAG;
    #[classattr] pub const SIDM_GALCENT_0SAG: u32 = SE_SIDM_GALCENT_0SAG;
    #[classattr] pub const SE_SIDM_J2000: u32 = SE_SIDM_J2000;
    #[classattr] pub const SIDM_J2000: u32 = SE_SIDM_J2000;
    #[classattr] pub const SE_SIDM_J1900: u32 = SE_SIDM_J1900;
    #[classattr] pub const SIDM_J1900: u32 = SE_SIDM_J1900;
    #[classattr] pub const SE_SIDM_B1950: u32 = SE_SIDM_B1950;
    #[classattr] pub const SIDM_B1950: u32 = SE_SIDM_B1950;
    #[classattr] pub const SE_SIDM_SURYASIDDHANTA: u32 = SE_SIDM_SURYASIDDHANTA;
    #[classattr] pub const SIDM_SURYASIDDHANTA: u32 = SE_SIDM_SURYASIDDHANTA;
    #[classattr] pub const SE_SIDM_SURYASIDDHANTA_MSUN: u32 = SE_SIDM_SURYASIDDHANTA_MSUN;
    #[classattr] pub const SIDM_SURYASIDDHANTA_MSUN: u32 = SE_SIDM_SURYASIDDHANTA_MSUN;
    #[classattr] pub const SE_SIDM_ARYABHATA: u32 = SE_SIDM_ARYABHATA;
    #[classattr] pub const SIDM_ARYABHATA: u32 = SE_SIDM_ARYABHATA;
    #[classattr] pub const SE_SIDM_ARYABHATA_MSUN: u32 = SE_SIDM_ARYABHATA_MSUN;
    #[classattr] pub const SIDM_ARYABHATA_MSUN: u32 = SE_SIDM_ARYABHATA_MSUN;
    #[classattr] pub const SE_SIDM_SS_REVATI: u32 = SE_SIDM_SS_REVATI;
    #[classattr] pub const SIDM_SS_REVATI: u32 = SE_SIDM_SS_REVATI;
    #[classattr] pub const SE_SIDM_SS_CITRA: u32 = SE_SIDM_SS_CITRA;
    #[classattr] pub const SIDM_SS_CITRA: u32 = SE_SIDM_SS_CITRA;
    #[classattr] pub const SE_SIDM_TRUE_CITRA: u32 = SE_SIDM_TRUE_CITRA;
    #[classattr] pub const SIDM_TRUE_CITRA: u32 = SE_SIDM_TRUE_CITRA;
    #[classattr] pub const SE_SIDM_TRUE_REVATI: u32 = SE_SIDM_TRUE_REVATI;
    #[classattr] pub const SIDM_TRUE_REVATI: u32 = SE_SIDM_TRUE_REVATI;
    #[classattr] pub const SE_SIDM_TRUE_PUSHYA: u32 = SE_SIDM_TRUE_PUSHYA;
    #[classattr] pub const SIDM_TRUE_PUSHYA: u32 = SE_SIDM_TRUE_PUSHYA;
    #[classattr] pub const SE_SIDM_GALCENT_RGILBRAND: u32 = SE_SIDM_GALCENT_RGILBRAND;
    #[classattr] pub const SIDM_GALCENT_RGILBRAND: u32 = SE_SIDM_GALCENT_RGILBRAND;
    #[classattr] pub const SE_SIDM_GALEQ_IAU1958: u32 = SE_SIDM_GALEQ_IAU1958;
    #[classattr] pub const SIDM_GALEQ_IAU1958: u32 = SE_SIDM_GALEQ_IAU1958;
    #[classattr] pub const SE_SIDM_GALEQ_TRUE: u32 = SE_SIDM_GALEQ_TRUE;
    #[classattr] pub const SIDM_GALEQ_TRUE: u32 = SE_SIDM_GALEQ_TRUE;
    #[classattr] pub const SE_SIDM_GALEQ_MULA: u32 = SE_SIDM_GALEQ_MULA;
    #[classattr] pub const SIDM_GALEQ_MULA: u32 = SE_SIDM_GALEQ_MULA;
    #[classattr] pub const SE_SIDM_GALALIGN_MARDYKS: u32 = SE_SIDM_GALALIGN_MARDYKS;
    #[classattr] pub const SIDM_GALALIGN_MARDYKS: u32 = SE_SIDM_GALALIGN_MARDYKS;
    #[classattr] pub const SE_SIDM_TRUE_MULA: u32 = SE_SIDM_TRUE_MULA;
    #[classattr] pub const SIDM_TRUE_MULA: u32 = SE_SIDM_TRUE_MULA;
    #[classattr] pub const SE_SIDM_GALCENT_MULA_WILHELM: u32 = SE_SIDM_GALCENT_MULA_WILHELM;
    #[classattr] pub const SIDM_GALCENT_MULA_WILHELM: u32 = SE_SIDM_GALCENT_MULA_WILHELM;
    #[classattr] pub const SE_SIDM_ARYABHATA_522: u32 = SE_SIDM_ARYABHATA_522;
    #[classattr] pub const SIDM_ARYABHATA_522: u32 = SE_SIDM_ARYABHATA_522;
    #[classattr] pub const SE_SIDM_BABYL_BRITTON: u32 = SE_SIDM_BABYL_BRITTON;
    #[classattr] pub const SIDM_BABYL_BRITTON: u32 = SE_SIDM_BABYL_BRITTON;
    #[classattr] pub const SE_SIDM_TRUE_SHERATAN: u32 = SE_SIDM_TRUE_SHERATAN;
    #[classattr] pub const SIDM_TRUE_SHERATAN: u32 = SE_SIDM_TRUE_SHERATAN;
    #[classattr] pub const SE_SIDM_TRUE_ALDEBARAN: u32 = SE_SIDM_TRUE_ALDEBARAN;
    #[classattr] pub const SIDM_TRUE_ALDEBARAN: u32 = SE_SIDM_TRUE_ALDEBARAN;
    #[classattr] pub const SE_SIDM_GALEQU_FIORENZA: u32 = SE_SIDM_GALEQU_FIORENZA;
    #[classattr] pub const SIDM_GALEQU_FIORENZA: u32 = SE_SIDM_GALEQU_FIORENZA;
    #[classattr] pub const SE_SIDM_VALENS_MOON: u32 = SE_SIDM_VALENS_MOON;
    #[classattr] pub const SIDM_VALENS_MOON: u32 = SE_SIDM_VALENS_MOON;
    #[classattr] pub const SE_SIDM_LAHIRI_1940: u32 = SE_SIDM_LAHIRI_1940;
    #[classattr] pub const SIDM_LAHIRI_1940: u32 = SE_SIDM_LAHIRI_1940;
    #[classattr] pub const SE_SIDM_LAHIRI_VP285: u32 = SE_SIDM_LAHIRI_VP285;
    #[classattr] pub const SIDM_LAHIRI_VP285: u32 = SE_SIDM_LAHIRI_VP285;
    #[classattr] pub const SE_SIDM_KP_SENTHILATHIBAN: u32 = SE_SIDM_KP_SENTHILATHIBAN;
    #[classattr] pub const SIDM_KP_SENTHILATHIBAN: u32 = SE_SIDM_KP_SENTHILATHIBAN;
    #[classattr] pub const SE_SIDM_LAHIRI_ICRC: u32 = SE_SIDM_LAHIRI_ICRC;
    #[classattr] pub const SIDM_LAHIRI_ICRC: u32 = SE_SIDM_LAHIRI_ICRC;
    #[classattr] pub const SE_SIDM_USER: u32 = SE_SIDM_USER;
    #[classattr] pub const SIDM_USER: u32 = SE_SIDM_USER;
    #[classattr] pub const SE_SIDBIT_ECL_T0: u32 = SE_SIDBIT_ECL_T0;
    #[classattr] pub const SIDBIT_ECL_T0: u32 = SE_SIDBIT_ECL_T0;
    #[classattr] pub const SE_SIDBIT_SSY_PLANE: u32 = SE_SIDBIT_SSY_PLANE;
    #[classattr] pub const SIDBIT_SSY_PLANE: u32 = SE_SIDBIT_SSY_PLANE;
    #[classattr] pub const SE_SIDBIT_USER_UT: u32 = SE_SIDBIT_USER_UT;
    #[classattr] pub const SIDBIT_USER_UT: u32 = SE_SIDBIT_USER_UT;
    #[classattr] pub const SE_SIDBIT_ECL_DATE: u32 = SE_SIDBIT_ECL_DATE;
    #[classattr] pub const SIDBIT_ECL_DATE: u32 = SE_SIDBIT_ECL_DATE;
    #[classattr] pub const SE_SIDBIT_NO_PREC_OFFSET: u32 = SE_SIDBIT_NO_PREC_OFFSET;
    #[classattr] pub const SIDBIT_NO_PREC_OFFSET: u32 = SE_SIDBIT_NO_PREC_OFFSET;
    #[classattr] pub const SE_SIDBIT_PREC_ORIG: u32 = SE_SIDBIT_PREC_ORIG;
    #[classattr] pub const SIDBIT_PREC_ORIG: u32 = SE_SIDBIT_PREC_ORIG;

    // Rise/Transit, AzAlt, Calendar constants
    #[classattr] pub const SE_CALC_RISE: i32 = SE_CALC_RISE;
    #[classattr] pub const CALC_RISE: i32 = SE_CALC_RISE;
    #[classattr] pub const SWE_CALC_RISE: i32 = SE_CALC_RISE;
    #[classattr] pub const SE_CALC_SET: i32 = SE_CALC_SET;
    #[classattr] pub const CALC_SET: i32 = SE_CALC_SET;
    #[classattr] pub const SWE_CALC_SET: i32 = SE_CALC_SET;
    #[classattr] pub const SE_CALC_MTRANSIT: i32 = SE_CALC_MTRANSIT;
    #[classattr] pub const CALC_MTRANSIT: i32 = SE_CALC_MTRANSIT;
    #[classattr] pub const SWE_CALC_MTRANSIT: i32 = SE_CALC_MTRANSIT;
    #[classattr] pub const SE_CALC_ITRANSIT: i32 = SE_CALC_ITRANSIT;
    #[classattr] pub const CALC_ITRANSIT: i32 = SE_CALC_ITRANSIT;
    #[classattr] pub const SWE_CALC_ITRANSIT: i32 = SE_CALC_ITRANSIT;
    #[classattr] pub const SE_BIT_DISC_CENTER: i32 = SE_BIT_DISC_CENTER;
    #[classattr] pub const BIT_DISC_CENTER: i32 = SE_BIT_DISC_CENTER;
    #[classattr] pub const SWE_BIT_DISC_CENTER: i32 = SE_BIT_DISC_CENTER;
    #[classattr] pub const SE_BIT_DISC_BOTTOM: i32 = SE_BIT_DISC_BOTTOM;
    #[classattr] pub const BIT_DISC_BOTTOM: i32 = SE_BIT_DISC_BOTTOM;
    #[classattr] pub const SWE_BIT_DISC_BOTTOM: i32 = SE_BIT_DISC_BOTTOM;
    #[classattr] pub const SE_BIT_GEOCTR_NO_ECL_LAT: i32 = SE_BIT_GEOCTR_NO_ECL_LAT;
    #[classattr] pub const BIT_GEOCTR_NO_ECL_LAT: i32 = SE_BIT_GEOCTR_NO_ECL_LAT;
    #[classattr] pub const SWE_BIT_GEOCTR_NO_ECL_LAT: i32 = SE_BIT_GEOCTR_NO_ECL_LAT;
    #[classattr] pub const SE_BIT_NO_REFRACTION: i32 = SE_BIT_NO_REFRACTION;
    #[classattr] pub const BIT_NO_REFRACTION: i32 = SE_BIT_NO_REFRACTION;
    #[classattr] pub const SWE_BIT_NO_REFRACTION: i32 = SE_BIT_NO_REFRACTION;
    #[classattr] pub const SE_BIT_CIVIL_TWILIGHT: i32 = SE_BIT_CIVIL_TWILIGHT;
    #[classattr] pub const BIT_CIVIL_TWILIGHT: i32 = SE_BIT_CIVIL_TWILIGHT;
    #[classattr] pub const SWE_BIT_CIVIL_TWILIGHT: i32 = SE_BIT_CIVIL_TWILIGHT;
    #[classattr] pub const SE_BIT_NAUTIC_TWILIGHT: i32 = SE_BIT_NAUTIC_TWILIGHT;
    #[classattr] pub const BIT_NAUTIC_TWILIGHT: i32 = SE_BIT_NAUTIC_TWILIGHT;
    #[classattr] pub const SWE_BIT_NAUTIC_TWILIGHT: i32 = SE_BIT_NAUTIC_TWILIGHT;
    #[classattr] pub const SE_BIT_ASTRO_TWILIGHT: i32 = SE_BIT_ASTRO_TWILIGHT;
    #[classattr] pub const BIT_ASTRO_TWILIGHT: i32 = SE_BIT_ASTRO_TWILIGHT;
    #[classattr] pub const SWE_BIT_ASTRO_TWILIGHT: i32 = SE_BIT_ASTRO_TWILIGHT;
    #[classattr] pub const SE_BIT_FIXED_DISC_SIZE: i32 = SE_BIT_FIXED_DISC_SIZE;
    #[classattr] pub const BIT_FIXED_DISC_SIZE: i32 = SE_BIT_FIXED_DISC_SIZE;
    #[classattr] pub const SWE_BIT_FIXED_DISC_SIZE: i32 = SE_BIT_FIXED_DISC_SIZE;
    #[classattr] pub const SE_BIT_HINDU_RISING: i32 = SE_BIT_HINDU_RISING;
    #[classattr] pub const BIT_HINDU_RISING: i32 = SE_BIT_HINDU_RISING;
    #[classattr] pub const SWE_BIT_HINDU_RISING: i32 = SE_BIT_HINDU_RISING;
    #[classattr] pub const SE_ECL2HOR: i32 = SE_ECL2HOR;
    #[classattr] pub const ECL2HOR: i32 = SE_ECL2HOR;
    #[classattr] pub const SWE_ECL2HOR: i32 = SE_ECL2HOR;
    #[classattr] pub const SE_EQU2HOR: i32 = SE_EQU2HOR;
    #[classattr] pub const EQU2HOR: i32 = SE_EQU2HOR;
    #[classattr] pub const SWE_EQU2HOR: i32 = SE_EQU2HOR;
    #[classattr] pub const SE_HOR2ECL: i32 = SE_HOR2ECL;
    #[classattr] pub const HOR2ECL: i32 = SE_HOR2ECL;
    #[classattr] pub const SWE_HOR2ECL: i32 = SE_HOR2ECL;
    #[classattr] pub const SE_HOR2EQU: i32 = SE_HOR2EQU;
    #[classattr] pub const HOR2EQU: i32 = SE_HOR2EQU;
    #[classattr] pub const SWE_HOR2EQU: i32 = SE_HOR2EQU;
    #[classattr] pub const SE_JUL_CAL: i32 = SE_JUL_CAL;
    #[classattr] pub const JUL_CAL: i32 = SE_JUL_CAL;
    #[classattr] pub const SWE_JUL_CAL: i32 = SE_JUL_CAL;
    #[classattr] pub const SE_GREG_CAL: i32 = SE_GREG_CAL;
    #[classattr] pub const GREG_CAL: i32 = SE_GREG_CAL;
    #[classattr] pub const SWE_GREG_CAL: i32 = SE_GREG_CAL;

    // Eclipse & Heliacal constants
    #[classattr] pub const SE_ECL_CENTRAL: i32 = SE_ECL_CENTRAL;
    #[classattr] pub const ECL_CENTRAL: i32 = SE_ECL_CENTRAL;
    #[classattr] pub const SWE_ECL_CENTRAL: i32 = SE_ECL_CENTRAL;
    #[classattr] pub const SE_ECL_NONCENTRAL: i32 = SE_ECL_NONCENTRAL;
    #[classattr] pub const ECL_NONCENTRAL: i32 = SE_ECL_NONCENTRAL;
    #[classattr] pub const SWE_ECL_NONCENTRAL: i32 = SE_ECL_NONCENTRAL;
    #[classattr] pub const SE_ECL_TOTAL: i32 = SE_ECL_TOTAL;
    #[classattr] pub const ECL_TOTAL: i32 = SE_ECL_TOTAL;
    #[classattr] pub const SWE_ECL_TOTAL: i32 = SE_ECL_TOTAL;
    #[classattr] pub const SE_ECL_ANNULAR: i32 = SE_ECL_ANNULAR;
    #[classattr] pub const ECL_ANNULAR: i32 = SE_ECL_ANNULAR;
    #[classattr] pub const SWE_ECL_ANNULAR: i32 = SE_ECL_ANNULAR;
    #[classattr] pub const SE_ECL_PARTIAL: i32 = SE_ECL_PARTIAL;
    #[classattr] pub const ECL_PARTIAL: i32 = SE_ECL_PARTIAL;
    #[classattr] pub const SWE_ECL_PARTIAL: i32 = SE_ECL_PARTIAL;
    #[classattr] pub const SE_ECL_ANNULAR_TOTAL: i32 = SE_ECL_ANNULAR_TOTAL;
    #[classattr] pub const ECL_ANNULAR_TOTAL: i32 = SE_ECL_ANNULAR_TOTAL;
    #[classattr] pub const SWE_ECL_ANNULAR_TOTAL: i32 = SE_ECL_ANNULAR_TOTAL;
    #[classattr] pub const SE_ECL_HYBRID: i32 = SE_ECL_HYBRID;
    #[classattr] pub const ECL_HYBRID: i32 = SE_ECL_HYBRID;
    #[classattr] pub const SWE_ECL_HYBRID: i32 = SE_ECL_HYBRID;
    #[classattr] pub const SE_ECL_PENUMBRAL: i32 = SE_ECL_PENUMBRAL;
    #[classattr] pub const ECL_PENUMBRAL: i32 = SE_ECL_PENUMBRAL;
    #[classattr] pub const SWE_ECL_PENUMBRAL: i32 = SE_ECL_PENUMBRAL;
    #[classattr] pub const SE_ECL_ALLTYPES_SOLAR: i32 = SE_ECL_ALLTYPES_SOLAR;
    #[classattr] pub const ECL_ALLTYPES_SOLAR: i32 = SE_ECL_ALLTYPES_SOLAR;
    #[classattr] pub const SWE_ECL_ALLTYPES_SOLAR: i32 = SE_ECL_ALLTYPES_SOLAR;
    #[classattr] pub const SE_ECL_ALLTYPES_LUNAR: i32 = SE_ECL_ALLTYPES_LUNAR;
    #[classattr] pub const ECL_ALLTYPES_LUNAR: i32 = SE_ECL_ALLTYPES_LUNAR;
    #[classattr] pub const SWE_ECL_ALLTYPES_LUNAR: i32 = SE_ECL_ALLTYPES_LUNAR;
    #[classattr] pub const SE_ECL_VISIBLE: i32 = SE_ECL_VISIBLE;
    #[classattr] pub const ECL_VISIBLE: i32 = SE_ECL_VISIBLE;
    #[classattr] pub const SWE_ECL_VISIBLE: i32 = SE_ECL_VISIBLE;
    #[classattr] pub const SE_ECL_MAX_VISIBLE: i32 = SE_ECL_MAX_VISIBLE;
    #[classattr] pub const ECL_MAX_VISIBLE: i32 = SE_ECL_MAX_VISIBLE;
    #[classattr] pub const SWE_ECL_MAX_VISIBLE: i32 = SE_ECL_MAX_VISIBLE;
    #[classattr] pub const SE_ECL_1ST_VISIBLE: i32 = SE_ECL_1ST_VISIBLE;
    #[classattr] pub const ECL_1ST_VISIBLE: i32 = SE_ECL_1ST_VISIBLE;
    #[classattr] pub const SWE_ECL_1ST_VISIBLE: i32 = SE_ECL_1ST_VISIBLE;
    #[classattr] pub const SE_ECL_2ND_VISIBLE: i32 = SE_ECL_2ND_VISIBLE;
    #[classattr] pub const ECL_2ND_VISIBLE: i32 = SE_ECL_2ND_VISIBLE;
    #[classattr] pub const SWE_ECL_2ND_VISIBLE: i32 = SE_ECL_2ND_VISIBLE;
    #[classattr] pub const SE_ECL_3RD_VISIBLE: i32 = SE_ECL_3RD_VISIBLE;
    #[classattr] pub const ECL_3RD_VISIBLE: i32 = SE_ECL_3RD_VISIBLE;
    #[classattr] pub const SWE_ECL_3RD_VISIBLE: i32 = SE_ECL_3RD_VISIBLE;
    #[classattr] pub const SE_ECL_4TH_VISIBLE: i32 = SE_ECL_4TH_VISIBLE;
    #[classattr] pub const ECL_4TH_VISIBLE: i32 = SE_ECL_4TH_VISIBLE;
    #[classattr] pub const SWE_ECL_4TH_VISIBLE: i32 = SE_ECL_4TH_VISIBLE;
    #[classattr] pub const SE_ECL_ONE_TRY: i32 = SE_ECL_ONE_TRY;
    #[classattr] pub const ECL_ONE_TRY: i32 = SE_ECL_ONE_TRY;
    #[classattr] pub const SWE_ECL_ONE_TRY: i32 = SE_ECL_ONE_TRY;
    #[classattr] pub const SE_HELIACAL_RISING: i32 = SE_HELIACAL_RISING;
    #[classattr] pub const HELIACAL_RISING: i32 = SE_HELIACAL_RISING;
    #[classattr] pub const SE_HELIACAL_SETTING: i32 = SE_HELIACAL_SETTING;
    #[classattr] pub const HELIACAL_SETTING: i32 = SE_HELIACAL_SETTING;
    #[classattr] pub const SE_EVENING_FIRST: i32 = SE_EVENING_FIRST;
    #[classattr] pub const EVENING_FIRST: i32 = SE_EVENING_FIRST;
    #[classattr] pub const SE_MORNING_LAST: i32 = SE_MORNING_LAST;
    #[classattr] pub const MORNING_LAST: i32 = SE_MORNING_LAST;
    #[classattr] pub const SE_HELFLAG_HIGH_PRECISION: i32 = SE_HELFLAG_HIGH_PRECISION;
    #[classattr] pub const HELFLAG_HIGH_PRECISION: i32 = SE_HELFLAG_HIGH_PRECISION;
}

/// Converts a calendar date `(year, month, day, hour)` into a Julian Date (`f64`).
/// Accepts `calendar` as an integer flag (`1` = Gregorian, `0` = Julian) or string (`"gregorian"`, `"julian"`, `"auto"`).
#[pyfunction]
#[pyo3(signature = (year, month, day, hour_ut=0.0, calendar=None))]
pub fn julday(
    year: i32,
    month: u8,
    day: u8,
    hour_ut: f64,
    calendar: Option<&Bound<'_, PyAny>>,
) -> PyResult<f64> {
    let cal = extract_calendar_kind(calendar)?;
    core_julday(year, month, day, hour_ut, cal).map_err(ephemeris_err_to_pyerr)
}

/// Converts a Julian Date (`jd`) into `(year, month, day, hour)`.
/// Accepts `calendar` as an integer flag (`1` = Gregorian, `0` = Julian) or string (`"gregorian"`, `"julian"`, `"auto"`).
#[pyfunction]
#[pyo3(signature = (jd, calendar=None))]
pub fn revjul(jd: f64, calendar: Option<&Bound<'_, PyAny>>) -> PyResult<(i32, u8, u8, f64)> {
    let cal = extract_calendar_kind(calendar)?;
    let d = core_revjul(jd, cal).map_err(ephemeris_err_to_pyerr)?;
    Ok((d.year, d.month, d.day, d.hour))
}

/// Converts a local wall-clock timestamp `(year, month, day, hour, minute, second)` in an IANA
/// timezone (`tz_name`, e.g. `"Europe/London"` or `"America/New_York"`, with `fold=0` or `1` for
/// repeated DST fall-back hours) to `(jd_utc, utc_year, utc_month, utc_day, utc_decimal_hour)` (`R1.1`).
#[pyfunction]
#[pyo3(signature = (year, month, day, hour=0, minute=0, second=0.0, tz_name="UTC", fold=0))]
#[allow(clippy::too_many_arguments)]
pub fn utc_from_local(
    py: Python<'_>,
    year: i32,
    month: u8,
    day: u8,
    hour: u8,
    minute: u8,
    second: f64,
    tz_name: &str,
    fold: u8,
) -> PyResult<(f64, i32, u8, u8, f64)> {
    if !(0.0..60.0).contains(&second) || fold > 1 {
        return Err(InvalidInputError::new_err(format!(
            "[InvalidInput] Invalid second={second} or fold={fold}"
        )));
    }
    let sec_int = second.floor() as u8;
    let microsec = ((second - f64::from(sec_int)) * 1_000_000.0).round() as u32;
    let zoneinfo = py.import("zoneinfo")?;
    let datetime_mod = py.import("datetime")?;
    let tz = zoneinfo
        .call_method1("ZoneInfo", (tz_name,))
        .map_err(|e| InvalidInputError::new_err(format!("[InvalidInput] Unknown IANA timezone '{tz_name}': {e}")))?;
    let utc_tz = datetime_mod.getattr("timezone")?.getattr("utc")?;
    let kwargs = pyo3::types::PyDict::new(py);
    kwargs.set_item("tzinfo", &tz)?;
    kwargs.set_item("fold", fold)?;
    let local_dt = datetime_mod
        .call_method(
            "datetime",
            (year, month, day, hour, minute, sec_int, microsec.min(999_999)),
            Some(&kwargs),
        )
        .map_err(|e| InvalidInputError::new_err(format!("[InvalidInput] Invalid local date/time: {e}")))?;
    let utc_dt = local_dt.call_method1("astimezone", (&utc_tz,))?;
    let u_year: i32 = utc_dt.getattr("year")?.extract()?;
    let u_month: u8 = utc_dt.getattr("month")?.extract()?;
    let u_day: u8 = utc_dt.getattr("day")?.extract()?;
    let u_hour: u8 = utc_dt.getattr("hour")?.extract()?;
    let u_min: u8 = utc_dt.getattr("minute")?.extract()?;
    let u_sec: u8 = utc_dt.getattr("second")?.extract()?;
    let u_micro: u32 = utc_dt.getattr("microsecond")?.extract()?;
    let dec_hour = f64::from(u_hour)
        + f64::from(u_min) / 60.0
        + (f64::from(u_sec) + f64::from(u_micro) * 1e-6) / 3_600.0;
    let jd_utc = core_julday(u_year, u_month, u_day, dec_hour, CalendarKind::Gregorian)
        .map_err(ephemeris_err_to_pyerr)?;
    Ok((jd_utc, u_year, u_month, u_day, dec_hour))
}

/// Converts a local wall-clock timestamp `(year, month, day, hour, minute, second)` in an IANA
/// timezone (`tz_name`, e.g. `"Europe/London"` or `"America/New_York"`, with `fold=0` or `1` for
/// repeated DST fall-back hours) directly into a scalar Julian Date (`jd_utc: f64`).
#[pyfunction]
#[pyo3(signature = (year, month, day, hour=0, minute=0, second=0.0, tz_name="UTC", fold=0, gregorian=true))]
#[allow(clippy::too_many_arguments)]
pub fn jd_utc_from_local(
    py: Python<'_>,
    year: i32,
    month: u8,
    day: u8,
    hour: u8,
    minute: u8,
    second: f64,
    tz_name: &str,
    fold: u8,
    gregorian: bool,
) -> PyResult<f64> {
    let (_, u_year, u_month, u_day, dec_hour) =
        utc_from_local(py, year, month, day, hour, minute, second, tz_name, fold)?;
    let cal = if gregorian {
        CalendarKind::Gregorian
    } else {
        CalendarKind::Julian
    };
    core_julday(u_year, u_month, u_day, dec_hour, cal).map_err(ephemeris_err_to_pyerr)
}

/// Maps an ecliptic longitude `lon_deg` into its Zodiac Sign, Decan, Nakshatra + Pada, and Human Design Gate.Line.Color.Tone.Base.
#[pyfunction]
pub fn map_zodiac_placement(py: Python<'_>, lon_deg: f64) -> PyResult<PyObject> {
    let res = core_map_zodiac_placement(lon_deg);
    serialize_to_attr_dict(py, &res)
}

/// Evaluates a single Parashari divisional chart (*Varga*) `D-n` (`division`) for `sidereal_lon_deg`.
#[pyfunction]
pub fn calc_varga(py: Python<'_>, sidereal_lon_deg: f64, division: u16) -> PyResult<PyObject> {
    let res = core_calc_varga(sidereal_lon_deg, division).map_err(ephemeris_err_to_pyerr)?;
    serialize_to_attr_dict(py, &res)
}

/// Evaluates all 16 Parashari *Shodashavarga* divisional charts for `sidereal_lon_deg`.
#[pyfunction]
pub fn calc_shodashavarga(py: Python<'_>, sidereal_lon_deg: f64) -> PyResult<PyObject> {
    let res = core_calc_shodashavarga(sidereal_lon_deg).map_err(ephemeris_err_to_pyerr)?;
    serialize_to_attr_dict(py, &res)
}

/// Returns canonical aspect definitions (`name`, `symbol`, `angle_deg`, `default_orb_deg`, `is_major`) for `preset`.
#[pyfunction]
#[pyo3(signature = (preset="major"))]
pub fn aspect_definitions(py: Python<'_>, preset: &str) -> PyResult<PyObject> {
    let p = parse_aspect_preset(preset).map_err(ephemeris_err_to_pyerr)?;
    let defs = core_aspect_definitions(p);
    serialize_to_attr_dict(py, &defs)
}

/// Evaluates house cusps and chart angles directly from `armc_deg`, `eps_true_deg`, and `latitude_deg`.
#[pyfunction]
#[pyo3(signature = (armc_deg, eps_true_deg, latitude_deg, system=None, polar_policy="strict", ayanamsa_deg=0.0, sun_decl_deg=0.0, include_speed=true))]
#[allow(clippy::too_many_arguments)]
pub fn houses_armc(
    armc_deg: f64,
    eps_true_deg: f64,
    latitude_deg: f64,
    system: Option<&Bound<'_, PyAny>>,
    polar_policy: &str,
    ayanamsa_deg: f64,
    sun_decl_deg: f64,
    include_speed: bool,
) -> PyResult<PyHouseResult> {
    let ch = extract_hsys_char(system)?;
    let sys = HouseSystem::from_code(ch).map_err(ephemeris_err_to_pyerr)?;
    let pol = parse_polar_policy(polar_policy).map_err(ephemeris_err_to_pyerr)?;
    let opts = HouseOptions {
        system: sys,
        polar_policy: pol,
        zodiac: Zodiac::Tropical,
        nutation: NutationMode::True,
        include_speed,
    };
    let res = houses_armc_with_sidereal(
        armc_deg,
        latitude_deg,
        eps_true_deg,
        sun_decl_deg,
        ayanamsa_deg,
        opts,
    )
    .map_err(ephemeris_err_to_pyerr)?;
    Ok(PyHouseResult::from(res))
}

/// Determines the house number (`1..=12`) for ecliptic longitude `lon_deg` given 12 cusps (`R4.5`).
#[pyfunction]
#[pyo3(signature = (lon_deg, cusps_deg, eps_deg=1e-9))]
pub fn house_of_lon(lon_deg: f64, cusps_deg: [f64; 12], eps_deg: f64) -> u8 {
    core_house_of_lon(lon_deg, &cusps_deg, eps_deg)
}

/// Looks up catalog metadata (including `vmag`) for a fixed star (`R6.1`).
#[pyfunction]
pub fn star_info(name_or_bayer: &str) -> PyResult<PyStarMetadata> {
    let info = core_star_info(name_or_bayer).map_err(ephemeris_err_to_pyerr)?;
    Ok(PyStarMetadata::from(info))
}

/// Lists all 52 cataloged fixed stars with metadata (`R6.1`).
#[pyfunction]
pub fn list_stars() -> Vec<PyStarMetadata> {
    core_list_stars()
        .into_iter()
        .map(PyStarMetadata::from)
        .collect()
}

/// Computes the 36 Gauquelin sectors (`[f64; 36]`) from `armc_deg`, `eps_true_deg`, `latitude_deg`.
#[pyfunction]
#[pyo3(signature = (armc_deg, eps_true_deg, latitude_deg, polar_policy="strict"))]
pub fn gauquelin_sectors(
    armc_deg: f64,
    eps_true_deg: f64,
    latitude_deg: f64,
    polar_policy: &str,
) -> PyResult<[f64; 36]> {
    let pol = parse_polar_policy(polar_policy).map_err(ephemeris_err_to_pyerr)?;
    let res = core_gauquelin_sectors(armc_deg, latitude_deg, eps_true_deg, pol)
        .map_err(ephemeris_err_to_pyerr)?;
    Ok(res.sectors_deg)
}

/// Verifies an offline `LXEPH-1.<payload>.<signature>` commercial license key and returns its `AttrDict` metadata.
#[pyfunction]
pub fn verify_license_key(py: Python<'_>, token: &str) -> PyResult<PyObject> {
    let info = lxephemeris_core::verify_license_key(token);
    serialize_to_attr_dict(py, &info)
}

/// Registers a commercial license key globally for the current process and returns its `AttrDict` metadata.
#[pyfunction]
pub fn register_commercial_license(py: Python<'_>, token: &str) -> PyResult<PyObject> {
    let info = lxephemeris_core::register_commercial_license(token);
    serialize_to_attr_dict(py, &info)
}

/// Clears any process-global commercial license key registered via `register_commercial_license`.
#[pyfunction]
pub fn clear_commercial_license() {
    lxephemeris_core::clear_commercial_license();
}

/// Resolves and returns the active license status (`AttrDict`) without requiring an `Engine` instance.
#[pyfunction]
pub fn license_info(py: Python<'_>) -> PyResult<PyObject> {
    let info = lxephemeris_core::resolve_license_info(None);
    serialize_to_attr_dict(py, &info)
}

/// Returns the canonical Swiss Ephemeris ayanamsa name for `sid_mode` (`swe_get_ayanamsa_name`).
#[cfg(feature = "swe-compat")]
#[pyfunction]
pub fn get_ayanamsa_name(sid_mode: i32) -> &'static str {
    SweAdapter::get_ayanamsa_name(sid_mode)
}

/// Signed angular difference `p1 - p2` in degrees normalized to `(-180, 180]` (`swe_difdeg2n`).
#[cfg(feature = "swe-compat")]
#[pyfunction]
pub fn difdeg2n(p1: f64, p2: f64) -> f64 {
    SweAdapter::difdeg2n(p1, p2)
}

/// Normalizes an angle in degrees to `[0, 360)` (`swe_degnorm`).
#[cfg(feature = "swe-compat")]
#[pyfunction]
pub fn degnorm(x: f64) -> f64 {
    SweAdapter::degnorm(x)
}

/// Coordinate transformation between ecliptic and equatorial spherical coordinates (`swe_cotrans`).
/// Accepts either `(xpo_seq, eps)` (2 arguments) or `(lon, lat, dist, eps)` (4 arguments).
#[cfg(feature = "swe-compat")]
#[pyfunction]
#[pyo3(signature = (lon, lat, dist=None, eps=None))]
pub fn cotrans(
    lon: &Bound<'_, PyAny>,
    lat: f64,
    dist: Option<f64>,
    eps: Option<f64>,
) -> PyResult<(f64, f64, f64)> {
    py_cotrans_impl(lon, lat, dist, eps)
}

/// Python module initializer (`import lxephemeris`).
#[pymodule]
pub fn lxephemeris(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("version", env!("CARGO_PKG_VERSION"))?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("Error", py.get_type::<LxEphemerisError>())?;
    m.add("LxEphemerisError", py.get_type::<LxEphemerisError>())?;
    m.add("InvalidInputError", py.get_type::<InvalidInputError>())?;
    m.add("UnsupportedBodyError", py.get_type::<UnsupportedBodyError>())?;
    m.add("UnsupportedModelError", py.get_type::<UnsupportedModelError>())?;
    m.add("OutOfRangeError", py.get_type::<OutOfRangeError>())?;
    m.add("PolarSingularityError", py.get_type::<PolarSingularityError>())?;
    m.add("ConvergenceError", py.get_type::<ConvergenceError>())?;
    m.add("MissingDataError", py.get_type::<MissingDataError>())?;
    m.add("CorruptKernelError", py.get_type::<CorruptKernelError>())?;
    m.add("AttrDict", get_attr_dict_type(py)?)?;

    m.add_class::<PySplitJd>()?;
    m.add_class::<PyObserver>()?;
    m.add_class::<PyCalcOptions>()?;
    m.add_class::<PyCalcResult>()?;
    m.add_class::<PyBatchItemResult>()?;
    m.add_class::<PyDeltaTInfo>()?;
    m.add_class::<PySiderealTimeInfo>()?;
    m.add_class::<PyEclNutInfo>()?;
    m.add_class::<PyStarMetadata>()?;
    m.add_class::<PyHeliacalResult>()?;
    m.add_class::<PyAyanamsaValue>()?;
    m.add_class::<PyHouseResult>()?;
    m.add_class::<PyHorizontalCoords>()?;
    m.add_class::<PySolarEclipseInfo>()?;
    m.add_class::<PyLunarEclipseInfo>()?;
    m.add_class::<PyPhenomenaResult>()?;
    m.add_class::<PyNodApsResult>()?;
    m.add_class::<PyLunarOccultationInfo>()?;
    m.add_class::<PyEngine>()?;
    #[cfg(feature = "swe-compat")]
    m.add_class::<PySweAdapter>()?;
    m.add_function(wrap_pyfunction!(julday, m)?)?;
    m.add_function(wrap_pyfunction!(revjul, m)?)?;
    m.add_function(wrap_pyfunction!(utc_from_local, m)?)?;
    m.add_function(wrap_pyfunction!(jd_utc_from_local, m)?)?;
    m.add_function(wrap_pyfunction!(map_zodiac_placement, m)?)?;
    m.add_function(wrap_pyfunction!(calc_varga, m)?)?;
    m.add_function(wrap_pyfunction!(calc_shodashavarga, m)?)?;
    m.add_function(wrap_pyfunction!(aspect_definitions, m)?)?;
    m.add_function(wrap_pyfunction!(houses_armc, m)?)?;
    m.add_function(wrap_pyfunction!(house_of_lon, m)?)?;
    m.add_function(wrap_pyfunction!(star_info, m)?)?;
    m.add_function(wrap_pyfunction!(list_stars, m)?)?;
    m.add_function(wrap_pyfunction!(gauquelin_sectors, m)?)?;
    m.add_function(wrap_pyfunction!(verify_license_key, m)?)?;
    m.add_function(wrap_pyfunction!(register_commercial_license, m)?)?;
    m.add_function(wrap_pyfunction!(clear_commercial_license, m)?)?;
    m.add_function(wrap_pyfunction!(license_info, m)?)?;
    #[cfg(feature = "swe-compat")]
    {
        m.add_function(wrap_pyfunction!(get_ayanamsa_name, m)?)?;
        m.add_function(wrap_pyfunction!(difdeg2n, m)?)?;
        m.add_function(wrap_pyfunction!(degnorm, m)?)?;
        m.add_function(wrap_pyfunction!(cotrans, m)?)?;
    }

    for (se_name, short_name, val) in [
        ("SE_NODBIT_MEAN", "NODBIT_MEAN", SE_NODBIT_MEAN),
        ("SE_NODBIT_OSCU", "NODBIT_OSCU", SE_NODBIT_OSCU),
        ("SE_NODBIT_OSCU_BAR", "NODBIT_OSCU_BAR", SE_NODBIT_OSCU_BAR),
        ("SE_NODBIT_FOPOINT", "NODBIT_FOPOINT", SE_NODBIT_FOPOINT),
    ] {
        m.add(se_name, val)?;
        m.add(short_name, val)?;
    }

    // Swiss Ephemeris migration constants (both SE_*/SEFLG_* and unprefixed pyswisseph aliases)
    #[cfg(feature = "swe-compat")]
    {
        for (se_name, short_name, val) in [
            ("SE_ECL_NUT", "ECL_NUT", SE_ECL_NUT),
            ("SE_SUN", "SUN", SE_SUN),
            ("SE_MOON", "MOON", SE_MOON),
            ("SE_MERCURY", "MERCURY", SE_MERCURY),
            ("SE_VENUS", "VENUS", SE_VENUS),
            ("SE_MARS", "MARS", SE_MARS),
            ("SE_JUPITER", "JUPITER", SE_JUPITER),
            ("SE_SATURN", "SATURN", SE_SATURN),
            ("SE_URANUS", "URANUS", SE_URANUS),
            ("SE_NEPTUNE", "NEPTUNE", SE_NEPTUNE),
            ("SE_PLUTO", "PLUTO", SE_PLUTO),
            ("SE_MEAN_NODE", "MEAN_NODE", SE_MEAN_NODE),
            ("SE_TRUE_NODE", "TRUE_NODE", SE_TRUE_NODE),
            ("SE_MEAN_APOG", "MEAN_APOG", SE_MEAN_APOG),
            ("SE_OSCU_APOG", "OSCU_APOG", SE_OSCU_APOG),
            ("SE_EARTH", "EARTH", SE_EARTH),
            ("SE_CHIRON", "CHIRON", SE_CHIRON),
            ("SE_PHOLUS", "PHOLUS", SE_PHOLUS),
            ("SE_CERES", "CERES", SE_CERES),
            ("SE_PALLAS", "PALLAS", SE_PALLAS),
            ("SE_JUNO", "JUNO", SE_JUNO),
            ("SE_VESTA", "VESTA", SE_VESTA),
            ("SE_INTP_APOG", "INTP_APOG", SE_INTP_APOG),
            ("SE_INTP_PERG", "INTP_PERG", SE_INTP_PERG),
            ("SE_FICT_OFFSET", "FICT_OFFSET", SE_FICT_OFFSET),
            ("SE_CUPIDO", "CUPIDO", SE_CUPIDO),
            ("SE_HADES", "HADES", SE_HADES),
            ("SE_ZEUS", "ZEUS", SE_ZEUS),
            ("SE_KRONOS", "KRONOS", SE_KRONOS),
            ("SE_APOLLON", "APOLLON", SE_APOLLON),
            ("SE_ADMETOS", "ADMETOS", SE_ADMETOS),
            ("SE_VULKANUS", "VULKANUS", SE_VULKANUS),
            ("SE_POSEIDON", "POSEIDON", SE_POSEIDON),
            ("SE_ISIS", "ISIS", SE_ISIS),
            ("SE_AST_OFFSET", "AST_OFFSET", SE_AST_OFFSET),
            ("SE_HOR2ECL", "HOR2ECL", SE_HOR2ECL),
            ("SE_HOR2EQU", "HOR2EQU", SE_HOR2EQU),
            ("SE_CALC_RISE", "CALC_RISE", SE_CALC_RISE),
            ("SE_CALC_SET", "CALC_SET", SE_CALC_SET),
            ("SE_CALC_MTRANSIT", "CALC_MTRANSIT", SE_CALC_MTRANSIT),
            ("SE_CALC_ITRANSIT", "CALC_ITRANSIT", SE_CALC_ITRANSIT),
            ("SE_BIT_DISC_CENTER", "BIT_DISC_CENTER", SE_BIT_DISC_CENTER),
            ("SE_BIT_DISC_BOTTOM", "BIT_DISC_BOTTOM", SE_BIT_DISC_BOTTOM),
            ("SE_BIT_GEOCTR_NO_ECL_LAT", "BIT_GEOCTR_NO_ECL_LAT", SE_BIT_GEOCTR_NO_ECL_LAT),
            ("SE_BIT_NO_REFRACTION", "BIT_NO_REFRACTION", SE_BIT_NO_REFRACTION),
            ("SE_BIT_CIVIL_TWILIGHT", "BIT_CIVIL_TWILIGHT", SE_BIT_CIVIL_TWILIGHT),
            ("SE_BIT_NAUTIC_TWILIGHT", "BIT_NAUTIC_TWILIGHT", SE_BIT_NAUTIC_TWILIGHT),
            ("SE_BIT_ASTRO_TWILIGHT", "BIT_ASTRO_TWILIGHT", SE_BIT_ASTRO_TWILIGHT),
            ("SE_BIT_FIXED_DISC_SIZE", "BIT_FIXED_DISC_SIZE", SE_BIT_FIXED_DISC_SIZE),
            ("SE_BIT_HINDU_RISING", "BIT_HINDU_RISING", SE_BIT_HINDU_RISING),
            ("SE_ECL2HOR", "ECL2HOR", SE_ECL2HOR),
            ("SE_EQU2HOR", "EQU2HOR", SE_EQU2HOR),
            ("SE_JUL_CAL", "JUL_CAL", SE_JUL_CAL),
            ("SE_GREG_CAL", "GREG_CAL", SE_GREG_CAL),
            ("SE_ECL_CENTRAL", "ECL_CENTRAL", SE_ECL_CENTRAL),
            ("SE_ECL_NONCENTRAL", "ECL_NONCENTRAL", SE_ECL_NONCENTRAL),
            ("SE_ECL_TOTAL", "ECL_TOTAL", SE_ECL_TOTAL),
            ("SE_ECL_ANNULAR", "ECL_ANNULAR", SE_ECL_ANNULAR),
            ("SE_ECL_PARTIAL", "ECL_PARTIAL", SE_ECL_PARTIAL),
            ("SE_ECL_ANNULAR_TOTAL", "ECL_ANNULAR_TOTAL", SE_ECL_ANNULAR_TOTAL),
            ("SE_ECL_HYBRID", "ECL_HYBRID", SE_ECL_HYBRID),
            ("SE_ECL_PENUMBRAL", "ECL_PENUMBRAL", SE_ECL_PENUMBRAL),
            ("SE_ECL_ALLTYPES_SOLAR", "ECL_ALLTYPES_SOLAR", SE_ECL_ALLTYPES_SOLAR),
            ("SE_ECL_ALLTYPES_LUNAR", "ECL_ALLTYPES_LUNAR", SE_ECL_ALLTYPES_LUNAR),
            ("SE_ECL_VISIBLE", "ECL_VISIBLE", SE_ECL_VISIBLE),
            ("SE_ECL_MAX_VISIBLE", "ECL_MAX_VISIBLE", SE_ECL_MAX_VISIBLE),
            ("SE_ECL_1ST_VISIBLE", "ECL_1ST_VISIBLE", SE_ECL_1ST_VISIBLE),
            ("SE_ECL_2ND_VISIBLE", "ECL_2ND_VISIBLE", SE_ECL_2ND_VISIBLE),
            ("SE_ECL_3RD_VISIBLE", "ECL_3RD_VISIBLE", SE_ECL_3RD_VISIBLE),
            ("SE_ECL_4TH_VISIBLE", "ECL_4TH_VISIBLE", SE_ECL_4TH_VISIBLE),
            ("SE_ECL_ONE_TRY", "ECL_ONE_TRY", SE_ECL_ONE_TRY),
            ("SE_HELIACAL_RISING", "HELIACAL_RISING", SE_HELIACAL_RISING),
            ("SE_HELIACAL_SETTING", "HELIACAL_SETTING", SE_HELIACAL_SETTING),
            ("SE_EVENING_FIRST", "EVENING_FIRST", SE_EVENING_FIRST),
            ("SE_MORNING_LAST", "MORNING_LAST", SE_MORNING_LAST),
            ("SE_HELFLAG_HIGH_PRECISION", "HELFLAG_HIGH_PRECISION", SE_HELFLAG_HIGH_PRECISION),
        ] {
            m.add(se_name, val)?;
            m.add(short_name, val)?;
            if let Some(rest) = se_name.strip_prefix("SE_") {
                m.add(format!("SWE_{rest}").as_str(), val)?;
            }
        }

        for (se_name, short_name, val) in [
            ("SEFLG_JPLEPH", "FLG_JPLEPH", SEFLG_JPLEPH),
            ("SEFLG_SWIEPH", "FLG_SWIEPH", SEFLG_SWIEPH),
            ("SEFLG_MOSEPH", "FLG_MOSEPH", SEFLG_MOSEPH),
            ("SEFLG_HELCTR", "FLG_HELCTR", SEFLG_HELCTR),
            ("SEFLG_TRUEPOS", "FLG_TRUEPOS", SEFLG_TRUEPOS),
            ("SEFLG_J2000", "FLG_J2000", SEFLG_J2000),
            ("SEFLG_NONUT", "FLG_NONUT", SEFLG_NONUT),
            ("SEFLG_SPEED3", "FLG_SPEED3", SEFLG_SPEED3),
            ("SEFLG_SPEED", "FLG_SPEED", SEFLG_SPEED),
            ("SEFLG_NOGDEFL", "FLG_NOGDEFL", SEFLG_NOGDEFL),
            ("SEFLG_NOABERR", "FLG_NOABERR", SEFLG_NOABERR),
            ("SEFLG_BARYCTR", "FLG_BARYCTR", SEFLG_BARYCTR),
            ("SEFLG_TOPOCTR", "FLG_TOPOCTR", SEFLG_TOPOCTR),
            ("SEFLG_EQUATORIAL", "FLG_EQUATORIAL", SEFLG_EQUATORIAL),
            ("SEFLG_XYZ", "FLG_XYZ", SEFLG_XYZ),
            ("SEFLG_RADIANS", "FLG_RADIANS", SEFLG_RADIANS),
            ("SEFLG_SIDEREAL", "FLG_SIDEREAL", SEFLG_SIDEREAL),
            ("SE_SIDM_FAGAN_BRADLEY", "SIDM_FAGAN_BRADLEY", SE_SIDM_FAGAN_BRADLEY),
            ("SE_SIDM_LAHIRI", "SIDM_LAHIRI", SE_SIDM_LAHIRI),
            ("SE_SIDM_DELUCE", "SIDM_DELUCE", SE_SIDM_DELUCE),
            ("SE_SIDM_RAMAN", "SIDM_RAMAN", SE_SIDM_RAMAN),
            ("SE_SIDM_USHASHASHI", "SIDM_USHASHASHI", SE_SIDM_USHASHASHI),
            ("SE_SIDM_KRISHNAMURTI", "SIDM_KRISHNAMURTI", SE_SIDM_KRISHNAMURTI),
            ("SE_SIDM_DJWHAL_KHUL", "SIDM_DJWHAL_KHUL", SE_SIDM_DJWHAL_KHUL),
            ("SE_SIDM_YUKTESHWAR", "SIDM_YUKTESHWAR", SE_SIDM_YUKTESHWAR),
            ("SE_SIDM_JN_BHASIN", "SIDM_JN_BHASIN", SE_SIDM_JN_BHASIN),
            ("SE_SIDM_BABYL_KUGLER1", "SIDM_BABYL_KUGLER1", SE_SIDM_BABYL_KUGLER1),
            ("SE_SIDM_BABYL_KUGLER2", "SIDM_BABYL_KUGLER2", SE_SIDM_BABYL_KUGLER2),
            ("SE_SIDM_BABYL_KUGLER3", "SIDM_BABYL_KUGLER3", SE_SIDM_BABYL_KUGLER3),
            ("SE_SIDM_BABYL_HUBER", "SIDM_BABYL_HUBER", SE_SIDM_BABYL_HUBER),
            ("SE_SIDM_BABYL_ETPSC", "SIDM_BABYL_ETPSC", SE_SIDM_BABYL_ETPSC),
            ("SE_SIDM_ALDEBARAN_15TAU", "SIDM_ALDEBARAN_15TAU", SE_SIDM_ALDEBARAN_15TAU),
            ("SE_SIDM_HIPPARCHOS", "SIDM_HIPPARCHOS", SE_SIDM_HIPPARCHOS),
            ("SE_SIDM_SASSANIAN", "SIDM_SASSANIAN", SE_SIDM_SASSANIAN),
            ("SE_SIDM_GALCENT_0SAG", "SIDM_GALCENT_0SAG", SE_SIDM_GALCENT_0SAG),
            ("SE_SIDM_J2000", "SIDM_J2000", SE_SIDM_J2000),
            ("SE_SIDM_J1900", "SIDM_J1900", SE_SIDM_J1900),
            ("SE_SIDM_B1950", "SIDM_B1950", SE_SIDM_B1950),
            ("SE_SIDM_SURYASIDDHANTA", "SIDM_SURYASIDDHANTA", SE_SIDM_SURYASIDDHANTA),
            ("SE_SIDM_SURYASIDDHANTA_MSUN", "SIDM_SURYASIDDHANTA_MSUN", SE_SIDM_SURYASIDDHANTA_MSUN),
            ("SE_SIDM_ARYABHATA", "SIDM_ARYABHATA", SE_SIDM_ARYABHATA),
            ("SE_SIDM_ARYABHATA_MSUN", "SIDM_ARYABHATA_MSUN", SE_SIDM_ARYABHATA_MSUN),
            ("SE_SIDM_SS_REVATI", "SIDM_SS_REVATI", SE_SIDM_SS_REVATI),
            ("SE_SIDM_SS_CITRA", "SIDM_SS_CITRA", SE_SIDM_SS_CITRA),
            ("SE_SIDM_TRUE_CITRA", "SIDM_TRUE_CITRA", SE_SIDM_TRUE_CITRA),
            ("SE_SIDM_TRUE_REVATI", "SIDM_TRUE_REVATI", SE_SIDM_TRUE_REVATI),
            ("SE_SIDM_TRUE_PUSHYA", "SIDM_TRUE_PUSHYA", SE_SIDM_TRUE_PUSHYA),
            ("SE_SIDM_GALCENT_RGILBRAND", "SIDM_GALCENT_RGILBRAND", SE_SIDM_GALCENT_RGILBRAND),
            ("SE_SIDM_GALEQ_IAU1958", "SIDM_GALEQ_IAU1958", SE_SIDM_GALEQ_IAU1958),
            ("SE_SIDM_GALEQ_TRUE", "SIDM_GALEQ_TRUE", SE_SIDM_GALEQ_TRUE),
            ("SE_SIDM_GALEQ_MULA", "SIDM_GALEQ_MULA", SE_SIDM_GALEQ_MULA),
            ("SE_SIDM_GALALIGN_MARDYKS", "SIDM_GALALIGN_MARDYKS", SE_SIDM_GALALIGN_MARDYKS),
            ("SE_SIDM_TRUE_MULA", "SIDM_TRUE_MULA", SE_SIDM_TRUE_MULA),
            ("SE_SIDM_GALCENT_MULA_WILHELM", "SIDM_GALCENT_MULA_WILHELM", SE_SIDM_GALCENT_MULA_WILHELM),
            ("SE_SIDM_ARYABHATA_522", "SIDM_ARYABHATA_522", SE_SIDM_ARYABHATA_522),
            ("SE_SIDM_BABYL_BRITTON", "SIDM_BABYL_BRITTON", SE_SIDM_BABYL_BRITTON),
            ("SE_SIDM_TRUE_SHERATAN", "SIDM_TRUE_SHERATAN", SE_SIDM_TRUE_SHERATAN),
            ("SE_SIDM_TRUE_ALDEBARAN", "SIDM_TRUE_ALDEBARAN", SE_SIDM_TRUE_ALDEBARAN),
            ("SE_SIDM_GALEQU_FIORENZA", "SIDM_GALEQU_FIORENZA", SE_SIDM_GALEQU_FIORENZA),
            ("SE_SIDM_VALENS_MOON", "SIDM_VALENS_MOON", SE_SIDM_VALENS_MOON),
            ("SE_SIDM_LAHIRI_1940", "SIDM_LAHIRI_1940", SE_SIDM_LAHIRI_1940),
            ("SE_SIDM_LAHIRI_VP285", "SIDM_LAHIRI_VP285", SE_SIDM_LAHIRI_VP285),
            ("SE_SIDM_KP_SENTHILATHIBAN", "SIDM_KP_SENTHILATHIBAN", SE_SIDM_KP_SENTHILATHIBAN),
            ("SE_SIDM_LAHIRI_ICRC", "SIDM_LAHIRI_ICRC", SE_SIDM_LAHIRI_ICRC),
            ("SE_SIDM_USER", "SIDM_USER", SE_SIDM_USER),
            ("SE_SIDBIT_ECL_T0", "SIDBIT_ECL_T0", SE_SIDBIT_ECL_T0),
            ("SE_SIDBIT_SSY_PLANE", "SIDBIT_SSY_PLANE", SE_SIDBIT_SSY_PLANE),
            ("SE_SIDBIT_USER_UT", "SIDBIT_USER_UT", SE_SIDBIT_USER_UT),
            ("SE_SIDBIT_ECL_DATE", "SIDBIT_ECL_DATE", SE_SIDBIT_ECL_DATE),
            ("SE_SIDBIT_NO_PREC_OFFSET", "SIDBIT_NO_PREC_OFFSET", SE_SIDBIT_NO_PREC_OFFSET),
            ("SE_SIDBIT_PREC_ORIG", "SIDBIT_PREC_ORIG", SE_SIDBIT_PREC_ORIG),
        ] {
            m.add(se_name, val)?;
            m.add(short_name, val)?;
            if let Some(rest) = se_name.strip_prefix("SE_") {
                m.add(format!("SWE_{rest}").as_str(), val)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::indexing_slicing)]
mod tests {
    use super::*;
    use lxephemeris_core::JD_J2000;
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
    fn test_python_vertical_slice_single_batch_and_errors() {
        pyo3::prepare_freethreaded_python();
        Python::with_gil(|py| {
            let engine = PyEngine::from_data_dir(
                py,
                workspace_data_dir(),
                Some(2),
                Some("embedded".to_owned()),
                false,
            )
            .expect("load PyEngine");

            // 1. Single calculation (Sun apparent ecliptic longitude at J2000.0 = 280.36816528 deg)
            let sun = engine
                .calc(py, "Sun", JD_J2000, 0.0, "TT", None)
                .expect("calc Sun");
            assert!((sun.lon_deg - 280.368_165_28).abs() < 1e-6);
            assert!(sun.lon_speed_deg_day > 1.01 && sun.lon_speed_deg_day < 1.03);
            assert!(sun.dist_km > 147_000_000.0 && sun.dist_km < 147_200_000.0);
            assert!(sun.warnings.is_empty());

            // 2. Named minor planet & fixed star
            let chiron = engine
                .calc(py, "Chiron", JD_J2000, 0.0, "TT", None)
                .expect("calc Chiron");
            assert!(chiron.lon_deg >= 0.0 && chiron.lon_deg < 360.0);

            let regulus = engine
                .calc_star_by_name(py, "Regulus", JD_J2000, 0.0, "TT", None)
                .expect("calc Regulus");
            assert!((regulus.lon_deg - 149.829_040_6).abs() < 1e-6);

            // 3. Ayanamsa & houses & eclipses & calendar
            let lahiri = engine
                .ayanamsa(py, "lahiri", JD_J2000, 0.0, "TT", "true")
                .expect("ayanamsa lahiri");
            assert!((lahiri.mean_deg - 23.857_092_3).abs() < 1e-4);

            let h_p = engine
                .houses(
                    py,
                    JD_J2000,
                    51.5074,
                    -0.1278,
                    Some(pyo3::types::PyString::new(py, "P").as_any()),
                    0.0,
                    "TT",
                    0.0,
                    "strict",
                    "tropical",
                    "lahiri",
                    "true",
                    true,
                )
                .expect("houses Placidus");
            assert_eq!(h_p.cusps_deg.len(), 12);

            let sol = engine
                .next_solar_eclipse(py, 2_460_400.0, 0.0, "TT", false, None, 36)
                .expect("next solar eclipse");
            assert_eq!(sol.kind, "Total");
            assert_eq!(sol.saros_series, 139);
            assert_eq!(sol.saros_member, 30);

            let jd = julday(2000, 1, 1, 12.0, None).expect("julday");
            assert!((jd - JD_J2000).abs() < 1e-12);
            let (y, m, d, hr) = revjul(jd, None).expect("revjul");
            assert_eq!((y, m, d), (2000, 1, 1));
            assert!((hr - 12.0).abs() < 1e-9);

            // 4. Flat batch calculation
            let bodies = vec!["Sun".to_owned(), "Moon".to_owned(), "mpc:2060".to_owned()];
            let epochs = vec![JD_J2000, JD_J2000 + 1.0];
            let flat = engine
                .calc_batch_flat_f64(py, bodies.clone(), epochs, "TT", None)
                .expect("calc_batch_flat_f64");
            assert_eq!(flat.len(), 2 * 3 * 6);
            assert!((flat[0] - sun.lon_deg).abs() < 1e-12);

            // 5. Per-item batch with out-of-coverage epoch
            let mixed_epochs = vec![JD_J2000, 2_200_000.0];
            let per_item = engine
                .calc_batch_per_item(py, bodies, mixed_epochs, "TT", None)
                .expect("calc_batch_per_item");
            assert!(per_item[0][0].ok);
            assert!(!per_item[1][0].ok);
            assert_eq!(
                per_item[1][0].error_code.as_deref(),
                Some("OutOfCoverage")
            );

            // 6. Structured error codes
            let err_body = engine
                .calc(py, "NotABody", JD_J2000, 0.0, "TT", None)
                .expect_err("must fail on unknown body");
            assert!(err_body.to_string().contains("[UnsupportedBody]"));

            // 7. SweAdapter check
            #[cfg(feature = "swe-compat")]
            {
                let swe = PySweAdapter::from_engine(&engine).expect("SweAdapter::from_engine");
                let (xx_sun, retflag) = swe
                    .calc(
                        py,
                        Some(JD_J2000),
                        Some(SE_SUN),
                        Some(SEFLG_SPEED),
                        None,
                        None,
                    )
                    .expect("swe.calc");
                assert!((xx_sun[0] - sun.lon_deg).abs() < 1e-12);
                assert_eq!((retflag as u32) & SEFLG_SPEED, SEFLG_SPEED);
            }
        });
    }

    #[test]
    fn test_python_moshier_pheno_nodaps_occult_and_swe_adapter_extensions() {
        pyo3::prepare_freethreaded_python();
        Python::with_gil(|py| {
            // 1. Pure fileless Moshier engine
            let mosh_eng = PyEngine::new_moshier(Some(2), None);
            assert!(mosh_eng.is_moshier_mode());

            let sun_mosh = mosh_eng
                .calc_moshier(py, "Sun", JD_J2000, 0.0, "TT", None)
                .expect("calc_moshier Sun");
            assert!((sun_mosh.lon_deg - 280.368).abs() < 0.1);

            let sun_calc = mosh_eng
                .calc(py, "Sun", JD_J2000, 0.0, "TT", None)
                .expect("calc Sun in moshier mode");
            assert!((sun_calc.lon_deg - sun_mosh.lon_deg).abs() < 1e-9);

            // 2. Fallback loader
            let fallback = PyEngine::from_data_dir_or_moshier(
                py,
                "/non_existent_dir_for_test_fallback".to_string(),
                None,
                None,
                false,
            )
            .expect("fallback to moshier");
            assert!(fallback.is_moshier_mode());

            let kernel_eng = PyEngine::from_data_dir_or_moshier(
                py,
                workspace_data_dir(),
                None,
                None,
                false,
            )
            .expect("load kernels with from_data_dir_or_moshier");
            assert!(!kernel_eng.is_moshier_mode());

            // 3. Custom stars loading
            let tmp_star = std::env::temp_dir().join("test_custom_sefstars.txt");
            std::fs::write(&tmp_star, "TestStar,alp Test,1991,12,0,0,10,0,0,0,0,0,10,2.5,test\n")
                .expect("write tmp sefstars");
            let loaded = fallback
                .load_custom_stars(tmp_star.display().to_string())
                .expect("load custom stars");
            assert!(loaded >= 1);
            let _ = std::fs::remove_file(&tmp_star);

            // 4. Phenomena
            let pheno = kernel_eng
                .phenomena(py, "Mars", JD_J2000, 0.0, "TT", None)
                .expect("phenomena Mars");
            assert!((0.0..=180.0).contains(&pheno.phase_angle_deg));
            assert!((0.0..=1.0).contains(&pheno.phase_illuminated));
            assert_eq!(pheno.to_swe_attr20().len(), 20);

            // 5. Nod aps
            let nod = kernel_eng
                .nod_aps(py, "Moon", JD_J2000, 0.0, "TT", SE_NODBIT_MEAN, None)
                .expect("nod_aps Moon");
            assert!((0.0..=360.0).contains(&nod.ascending_node[0]));
            assert!((0.0..=360.0).contains(&nod.descending_node[0]));
            assert!((0.0..=360.0).contains(&nod.asc_node()[0]));
            assert!((0.0..=360.0).contains(&nod.desc_node()[0]));

            // 6. House pos, time equation, azalt_rev
            let hp = kernel_eng
                .house_pos(
                    180.0,
                    51.5074,
                    23.439,
                    280.368,
                    0.0,
                    Some(pyo3::types::PyString::new(py, "P").as_any()),
                )
                .expect("house_pos");
            assert!((1.0..=13.0).contains(&hp));

            let teq = kernel_eng
                .time_equation(JD_J2000, 0.0, "TT")
                .expect("time_equation");
            assert!(teq.abs() < 60.0);

            let az_coords = kernel_eng
                .azalt_rev(
                    JD_J2000, 1, 51.5074, -0.1278, 180.0, 45.0, 0.0, "UT1", 0.0,
                )
                .expect("azalt_rev");
            assert!((0.0..=360.0).contains(&az_coords[0]));
            assert!((-90.0..=90.0).contains(&az_coords[1]));

            // 7. Lunar occultations
            let jd_2015 = 2_457_023.5;
            let occ = kernel_eng
                .next_lunar_occultation(py, "Aldebaran", jd_2015, 0.0, "TT", false, None, 12)
                .expect("next_lunar_occultation");
            assert_eq!(occ.to_swe_glob_tret().len(), 10);
            assert_eq!(occ.to_swe_tret10().len(), 10);
            assert!(occ.swe_retflag() >= 0);

            // 8. SweAdapter extensions & Moshier adapter
            #[cfg(feature = "swe-compat")]
            {
                let swe_mosh = PySweAdapter::new_moshier();
                let (sun_swe, swe_flg) = swe_mosh
                    .calc(
                        py,
                        Some(JD_J2000),
                        Some(SE_SUN),
                        Some(SEFLG_SPEED),
                        None,
                        None,
                    )
                    .expect("swe_mosh.calc");
                assert!((sun_swe[0] - 280.368).abs() < 0.1);
                assert_eq!((swe_flg as u32) & SEFLG_SPEED, SEFLG_SPEED);

                // Sidereal modes including custom SE_SIDM_USER
                swe_mosh
                    .set_sid_mode(SE_SIDM_USER, 2451545.0, 24.0)
                    .expect("set_sid_mode SE_SIDM_USER");
                swe_mosh
                    .set_sid_mode(SE_SIDM_LAHIRI, 0.0, 0.0)
                    .expect("set_sid_mode SE_SIDM_LAHIRI");

                // SweAdapter pheno
                let (p20_ut, flg_ut) = swe_mosh
                    .pheno_ut(py, Some(JD_J2000), Some(SE_MARS), Some(0), None, None)
                    .expect("pheno_ut");
                assert_eq!(p20_ut.len(), 20);
                assert_eq!(flg_ut as u32, SEFLG_MOSEPH);

                let (p20_tt, _) = swe_mosh
                    .pheno(py, Some(JD_J2000), Some(SE_MARS), Some(0), None, None)
                    .expect("pheno");
                assert_eq!(p20_tt.len(), 20);

                // SweAdapter nod_aps
                let (asc, desc, peri, apo, _retflag_ut) = swe_mosh
                    .nod_aps_ut(
                        py,
                        Some(JD_J2000),
                        Some(SE_MOON),
                        Some(258),
                        Some(SE_NODBIT_MEAN),
                        None,
                        None,
                        None,
                        None,
                    )
                    .expect("nod_aps_ut");
                assert_eq!(asc.len(), 6);
                assert_eq!(desc.len(), 6);
                assert_eq!(peri.len(), 6);
                assert_eq!(apo.len(), 6);

                let (asc_tt, _, _, _, _retflag_tt) = swe_mosh
                    .nod_aps(
                        py,
                        Some(JD_J2000),
                        Some(SE_MOON),
                        Some(258),
                        Some(SE_NODBIT_MEAN),
                        None,
                        None,
                        None,
                        None,
                    )
                    .expect("nod_aps");
                assert_eq!(asc_tt.len(), 6);

                // SweAdapter house_pos, time_equ, azalt_rev
                let xpin_list = pyo3::types::PyList::new(py, [280.368, 0.0]).expect("PyList");
                let shp = swe_mosh
                    .house_pos(
                        180.0,
                        51.5074,
                        23.439,
                        xpin_list.as_any(),
                        Some(pyo3::types::PyString::new(py, "P").as_any()),
                    )
                    .expect("swe house_pos");
                assert!((1.0..=13.0).contains(&shp));

                let ste = swe_mosh.time_equ(JD_J2000).expect("swe time_equ");
                assert!(ste.abs() < 0.1); // days

                let xin_list = pyo3::types::PyList::new(py, [180.0, 45.0]).expect("PyList");
                let s_coords = swe_mosh
                    .azalt_rev(JD_J2000, SE_HOR2EQU, [-0.1278, 51.5074, 0.0], xin_list.as_any(), None)
                    .expect("swe azalt_rev");
                assert!((0.0..=360.0).contains(&s_coords[0]));
                assert!((-90.0..=90.0).contains(&s_coords[1]));

                // SweAdapter occultations
                let swe_kernel = PySweAdapter::from_engine(&kernel_eng).expect("SweAdapter::from_engine");
                let aldebaran_str = pyo3::types::PyString::new(py, "Aldebaran");
                let (ret_glob, tret_glob) = swe_kernel
                    .lun_occult_when_glob(
                        py,
                        jd_2015,
                        aldebaran_str.as_any(),
                        Some(pyo3::types::PyInt::new(py, 258).as_any()),
                        Some(pyo3::types::PyInt::new(py, 0).as_any()),
                        Some(pyo3::types::PyBool::new(py, false).as_any()),
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                    )
                    .expect("lun_occult_when_glob");
                assert_eq!(tret_glob.len(), 10);
                assert!(ret_glob >= 0);

                let geopos_list = pyo3::types::PyList::new(py, [0.0, 51.5, 0.0]).expect("PyList");
                let (ret_loc, tret_loc, attr_loc) = swe_kernel
                    .lun_occult_when_loc(
                        py,
                        jd_2015,
                        aldebaran_str.as_any(),
                        geopos_list.as_any(),
                        Some(pyo3::types::PyInt::new(py, 258).as_any()),
                        Some(pyo3::types::PyBool::new(py, false).as_any()),
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                    )
                    .expect("lun_occult_when_loc");
                assert_eq!(tret_loc.len(), 10);
                assert_eq!(attr_loc.len(), 20);
                assert!(ret_loc >= 0);

                // Verify constants
                assert_eq!(SE_NODBIT_MEAN, 1);
                assert_eq!(SE_NODBIT_OSCU, 2);
                assert_eq!(SE_NODBIT_OSCU_BAR, 4);
                assert_eq!(SE_NODBIT_FOPOINT, 256);
                assert_eq!(SE_HOR2ECL, 0);
                assert_eq!(SE_HOR2EQU, 1);
                assert_eq!(SE_SIDM_USER, 255);
                assert_eq!(SE_SIDM_FAGAN_BRADLEY, 0);
                assert_eq!(SE_SIDM_LAHIRI, 1);
                assert_eq!(SE_SIDM_VALENS_MOON, 42);
            }
        });
    }

    #[test]
    fn test_python_wave_f_events_and_parans() {
        pyo3::prepare_freethreaded_python();
        Python::with_gil(|py| {
            let eng = PyEngine::from_data_dir(
                py,
                workspace_data_dir(),
                Some(2),
                Some("embedded".to_owned()),
                false,
            )
            .expect("load PyEngine");

            let jd_2026_start = julday(2026, 1, 1, 0.0, None).expect("2026 start");
            let jd_2026_end = julday(2026, 12, 31, 23.99, None).expect("2026 end");

            // 1. solar_terms_in_range & seasons_in_year
            let terms_obj = eng
                .solar_terms_in_range(py, jd_2026_start, jd_2026_end, "UT1")
                .expect("solar_terms_in_range");
            let terms = terms_obj.bind(py).downcast::<PyList>().expect("PyList");
            assert_eq!(terms.len(), 24);

            let seasons_obj = eng.seasons_in_year(py, 2026).expect("seasons_in_year");
            let seasons = seasons_obj.bind(py).downcast::<PyDict>().expect("PyDict");
            assert_eq!(
                seasons
                    .get_item("year")
                    .expect("get_item")
                    .expect("year")
                    .extract::<i32>()
                    .expect("i32"),
                2026
            );

            // 2. search_body_longitude & human_design_epoch
            let cross_obj = eng
                .search_body_longitude(py, "Sun", jd_2026_start, 0.0, true, 120.0, None, "UT1", 0.0)
                .expect("search_body_longitude");
            let cross = cross_obj.bind(py).downcast::<PyDict>().expect("PyDict");
            let actual_lon: f64 = cross
                .get_item("actual_lon_deg")
                .expect("get_item")
                .expect("actual_lon_deg")
                .extract()
                .expect("f64");
            assert!(actual_lon.abs() < 1e-7 || (360.0 - actual_lon).abs() < 1e-7);

            let birth_jd = julday(1990, 6, 15, 14.5, None).expect("birth jd");
            let hd_obj = eng
                .human_design_epoch(py, birth_jd, 0.0, "UT1")
                .expect("human_design_epoch");
            let hd = hd_obj.bind(py).downcast::<PyDict>().expect("PyDict");
            let hd_jd_ut1: f64 = hd
                .get_item("jd_ut1")
                .expect("get_item")
                .expect("jd_ut1")
                .extract()
                .expect("f64");
            assert!((87.0..=91.0).contains(&(birth_jd - hd_jd_ut1)));

            // 3. search_lunar_phase, lunar_phases_in_range, lunar_phase_summary
            let lp_obj = eng
                .search_lunar_phase(py, jd_2026_start, 180.0, true, 0.0, "UT1")
                .expect("search_lunar_phase");
            let lp = lp_obj.bind(py).downcast::<PyDict>().expect("PyDict");
            assert_eq!(
                lp.get_item("kind")
                    .expect("get_item")
                    .expect("kind")
                    .extract::<String>()
                    .expect("str"),
                "full_moon"
            );

            let lps_obj = eng
                .lunar_phases_in_range(py, jd_2026_start, jd_2026_start + 45.0, "UT1")
                .expect("lunar_phases_in_range");
            let lps = lps_obj.bind(py).downcast::<PyList>().expect("PyList");
            assert!(lps.len() >= 6);

            let sum_obj = eng
                .lunar_phase_summary(py, jd_2026_start + 24.5, 0.0, "UT1")
                .expect("lunar_phase_summary");
            let sum = sum_obj.bind(py).downcast::<PyDict>().expect("PyDict");
            let p28: u8 = sum
                .get_item("phase_28_index")
                .expect("get_item")
                .expect("phase_28_index")
                .extract()
                .expect("u8");
            assert!((1..=28).contains(&p28));

            // 4. search_stations_in_range & body_longitude_crossings_in_range
            let st_obj = eng
                .search_stations_in_range(py, "Saturn", jd_2026_start, jd_2026_end, None, "UT1")
                .expect("search_stations_in_range");
            let st_list = st_obj.bind(py).downcast::<PyList>().expect("PyList");
            assert_eq!(st_list.len(), 2);

            let st0 = st_list
                .get_item(0)
                .expect("st0")
                .downcast::<PyDict>()
                .expect("PyDict")
                .clone();
            let st1 = st_list
                .get_item(1)
                .expect("st1")
                .downcast::<PyDict>()
                .expect("PyDict")
                .clone();
            let lon0: f64 = st0
                .get_item("lon_deg")
                .expect("get_item")
                .expect("lon0")
                .extract()
                .expect("f64");
            let lon1: f64 = st1
                .get_item("lon_deg")
                .expect("get_item")
                .expect("lon1")
                .extract()
                .expect("f64");
            let mid_lon = 0.5 * (lon0 + lon1);
            let triple_obj = eng
                .body_longitude_crossings_in_range(
                    py,
                    "Saturn",
                    jd_2026_start,
                    jd_2026_start + 455.0,
                    mid_lon,
                    None,
                    "UT1",
                )
                .expect("body_longitude_crossings_in_range");
            let triple = triple_obj.bind(py).downcast::<PyList>().expect("PyList");
            assert_eq!(triple.len(), 3);

            // 5. search_declination_extrema_in_range
            let jd_2025_03_01 = julday(2025, 3, 1, 0.0, None).expect("2025-03-01");
            let dec_obj = eng
                .search_declination_extrema_in_range(
                    py,
                    "Moon",
                    jd_2025_03_01,
                    jd_2025_03_01 + 30.0,
                    None,
                    "UT1",
                )
                .expect("search_declination_extrema_in_range");
            let dec_list = dec_obj.bind(py).downcast::<PyList>().expect("PyList");
            assert!(dec_list.len() >= 2);

            // 6. find_parans
            let london = PyObserver::new(51.5074, -0.1278, 0.0).expect("London");
            let jd_subj_d = julday(1969, 5, 17, 14.0, None).expect("Subject D");
            let parans_obj = eng
                .find_parans(py, jd_subj_d, &london, 4.0, None, None, 0.0, "UT1")
                .expect("find_parans");
            let parans = parans_obj.bind(py).downcast::<PyList>().expect("PyList");
            assert!(!parans.is_empty());
        });
    }

    #[test]
    fn test_python_wave_h_ergonomics_aspects_and_almanac() {
        pyo3::prepare_freethreaded_python();
        Python::with_gil(|py| {
            let eng = PyEngine::new(py, Some(workspace_data_dir()), None, None, false)
                .expect("PyEngine::new");
            #[cfg(feature = "swe-compat")]
            {
                let swe = PySweAdapter::new(py, None).expect("PySweAdapter::new(None)");
                let (xx, _) = swe
                    .calc_ut(
                        py,
                        Some(JD_J2000),
                        Some(SE_SUN),
                        Some(SEFLG_SPEED),
                        None,
                        None,
                    )
                    .expect("calc_ut");
                assert!((xx[0] - 280.368).abs() < 0.1);
            }

            let jd_utc =
                jd_utc_from_local(py, 2000, 1, 1, 12, 0, 0.0, "UTC", 0, true).expect("jd_utc");
            assert!((jd_utc - JD_J2000).abs() < 1e-12);

            let sun = eng
                .calc(py, "Sun", JD_J2000, 0.0, "TT", None)
                .expect("calc Sun");
            let sun_dict = sun.to_dict(py).expect("to_dict");
            let lon_attr: f64 = sun_dict
                .bind(py)
                .getattr("lon_deg")
                .expect("attr lon_deg")
                .extract()
                .expect("f64");
            assert!((lon_attr - sun.lon_deg).abs() < 1e-12);
            assert!(sun.to_json().expect("to_json").contains("lon_deg"));

            let place = map_zodiac_placement(py, 280.368).expect("map_zodiac_placement");
            let sign_name: String = place
                .bind(py)
                .getattr("sign_name")
                .expect("sign_name")
                .extract()
                .expect("str");
            assert_eq!(sign_name, "Capricorn");

            let panch = eng
                .panchanga(py, JD_J2000, 28.6139, 77.2090, 0.0, "lahiri", "UT1")
                .expect("panchanga");
            let tithi_num: u8 = panch
                .bind(py)
                .getattr("tithi")
                .expect("tithi")
                .getattr("tithi_index")
                .expect("tithi_index")
                .extract()
                .expect("u8");
            assert!((1..=30).contains(&tithi_num));

            let snap = eng
                .chart_snapshot(
                    py, JD_J2000, 51.5074, -0.1278, 0.0, None, None, "major", 1.0, "UT1",
                    "tropical", None,
                )
                .expect("chart_snapshot");
            let bodies_list = snap
                .bind(py)
                .getattr("bodies")
                .expect("bodies")
                .downcast_into::<PyList>()
                .expect("PyList");
            assert_eq!(bodies_list.len(), 13);
        });
    }
}
