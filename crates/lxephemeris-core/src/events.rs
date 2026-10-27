//! High-precision astronomical and astrological event solvers (`Wave F`):
//! - Ecliptic longitude crossings (`search_body_longitude`, `body_longitude_crossings_in_range`),
//!   supporting Tropical and Sidereal zodiacs and retrograde triple-pass returns (`R9.1`, `R9.2`).
//! - Chinese 24 Solar Terms (`solar_terms_in_range`, `SolarTermEvent`) and 4 cardinal seasons
//!   (`seasons_in_year`, `SeasonsOfYear`), permanently eliminating the `ephem.Sun().hlong` $180^\circ$
//!   inversion bug (`R12`).
//! - Human Design `-88°` solar arc epoch solver (`human_design_epoch`, `R9.2`).
//! - Principal and custom lunar phases (`search_lunar_phase`, `lunar_phases_in_range`) and
//!   comprehensive lunar phase summary (`lunar_phase_summary`, `LunarPhaseSummary`) with
//!   8-phase and 28-phase (`kerykeion` replacement) indices (`R0`, `R2.2`, `R9.2`).
//! - Planetary retrograde and direct stations $\dot{\lambda}(t) = 0$ (`search_stations_in_range`)
//!   and lunar/planetary declination extrema / standstills $\dot{\delta}(t) = 0$
//!   (`search_declination_extrema_in_range`, `R9.2`).
//! - Bernadette Brady diurnal fixed-star parans (`find_parans`, `ParanEvent`, `ParanAngle`, `ParanOptions`, `R6.2`).

use crate::astrometry::{
    azalt, calc_body_in_moment, AtmosphereParams, CalcResult, RiseTransEvent, C_KM_S,
};
use crate::engine::Engine;
use crate::error::EphemerisError;
use crate::frame::{cart_to_spherical, vec_add, vec_unit};
use crate::moment::Moment;
use crate::stars::{eval_star_unit_vector, find_star_by_hip, StarEntry};
use crate::time::{gast_rad, julday_split, CalendarKind, TimeInput, Tt, Ut1};
use crate::types::{
    AbMode, BodyId, CalcOptions, Equinox, LunarPointId, NutationMode, Observer, Origin, Plane,
    PlanetId, Zodiac, MINOR_BODIES_68, SECONDS_PER_DAY,
};
use serde::{Deserialize, Serialize};
use sofars::astro::pvtob;

/// Canonical Pinyin names of the 24 Chinese Solar Terms (*Jieqi*), indexed `0..=23` starting
/// from `0°` geocentric apparent tropical solar longitude (`"Chunfen"` / Vernal Equinox).
pub const SOLAR_TERM_PINYIN_NAMES: [&str; 24] = [
    "Chunfen",     // 0:   0° (Vernal Equinox, Zhongqi major term)
    "Qingming",    // 1:  15°
    "Guyu",        // 2:  30° (Zhongqi major term)
    "Lixia",       // 3:  45°
    "Xiaoman",     // 4:  60° (Zhongqi major term)
    "Mangzhong",   // 5:  75°
    "Xiazhi",      // 6:  90° (Summer Solstice, Zhongqi major term)
    "Xiaoshu",     // 7: 105°
    "Dashu",       // 8: 120° (Zhongqi major term)
    "Liqiu",       // 9: 135°
    "Chushu",      // 10: 150° (Zhongqi major term)
    "Bailu",       // 11: 165°
    "Qiufen",      // 12: 180° (Autumnal Equinox, Zhongqi major term)
    "Hanlu",       // 13: 195°
    "Shuangjiang", // 14: 210° (Zhongqi major term)
    "Lidong",      // 15: 225°
    "Xiaoxue",     // 16: 240° (Zhongqi major term)
    "Daxue",       // 17: 255°
    "Dongzhi",     // 18: 270° (Winter Solstice, Zhongqi major term)
    "Xiaohan",     // 19: 285°
    "Dahan",       // 20: 300° (Zhongqi major term)
    "Lichun",      // 21: 315°
    "Yushui",      // 22: 330° (Zhongqi major term)
    "Jingzhe",     // 23: 345°
];

/// Canonical Chinese character names of the 24 Solar Terms (*Jieqi*), indexed `0..=23` starting
/// from `0°` (`"春分"` / Vernal Equinox).
pub const SOLAR_TERM_CHINESE_NAMES: [&str; 24] = [
    "春分", // 0:   0°
    "清明", // 1:  15°
    "谷雨", // 2:  30°
    "立夏", // 3:  45°
    "小满", // 4:  60°
    "芒种", // 5:  75°
    "夏至", // 6:  90°
    "小暑", // 7: 105°
    "大暑", // 8: 120°
    "立秋", // 9: 135°
    "处暑", // 10: 150°
    "白露", // 11: 165°
    "秋分", // 12: 180°
    "寒露", // 13: 195°
    "霜降", // 14: 210°
    "立冬", // 15: 225°
    "小雪", // 16: 240°
    "大雪", // 17: 255°
    "冬至", // 18: 270°
    "小寒", // 19: 285°
    "大寒", // 20: 300°
    "立春", // 21: 315°
    "雨水", // 22: 330°
    "惊蛰", // 23: 345°
];

/// Canonical English names of the 24 Solar Terms (*Jieqi*), indexed `0..=23` starting
/// from `0°` (`"Vernal Equinox"`).
pub const SOLAR_TERM_ENGLISH_NAMES: [&str; 24] = [
    "Vernal Equinox",       // 0:   0°
    "Clear and Bright",     // 1:  15°
    "Grain Rain",           // 2:  30°
    "Start of Summer",      // 3:  45°
    "Grain Buds",           // 4:  60°
    "Grain in Ear",         // 5:  75°
    "Summer Solstice",      // 6:  90°
    "Minor Heat",           // 7: 105°
    "Major Heat",           // 8: 120°
    "Start of Autumn",      // 9: 135°
    "End of Heat",          // 10: 150°
    "White Dew",            // 11: 165°
    "Autumnal Equinox",     // 12: 180°
    "Cold Dew",             // 13: 195°
    "Frost's Descent",      // 14: 210°
    "Start of Winter",      // 15: 225°
    "Minor Snow",           // 16: 240°
    "Major Snow",           // 17: 255°
    "Winter Solstice",      // 18: 270°
    "Minor Cold",           // 19: 285°
    "Major Cold",           // 20: 300°
    "Start of Spring",      // 21: 315°
    "Rain Water",           // 22: 330°
    "Awakening of Insects", // 23: 345°
];

/// Canonical names of the 8 lunar phase octants (`0..=7`, $45^\circ$ each).
pub const LUNAR_PHASE_8_NAMES: [&str; 8] = [
    "New Moon",
    "Waxing Crescent",
    "First Quarter",
    "Waxing Gibbous",
    "Full Moon",
    "Waning Gibbous",
    "Last Quarter",
    "Waning Crescent",
];

/// Wraps an angle in degrees into `[-180.0, +180.0)`.
#[inline]
#[must_use]
pub fn wrap_pm180(deg: f64) -> f64 {
    (deg + 180.0).rem_euclid(360.0) - 180.0
}

/// Wraps an angle in radians into `[-pi, +pi)`.
#[inline]
fn wrap_pm_pi(rad: f64) -> f64 {
    (rad + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI
}

/// Human-readable canonical name for a [`BodyId`].
#[must_use]
pub fn body_display_name(body: BodyId) -> String {
    match body {
        BodyId::Planet(p) => p.name().to_owned(),
        BodyId::MinorPlanet(mpc) => MINOR_BODIES_68
            .binary_search_by_key(&mpc, |s| s.mpc)
            .ok()
            .and_then(|idx| MINOR_BODIES_68.get(idx))
            .map_or_else(|| format!("MPC {}", mpc.0), |s| s.name.to_owned()),
        BodyId::LunarPoint(lp) => lp.astra_id().to_owned(),
        BodyId::Hypothetical(h) => h.astra_id().to_owned(),
        BodyId::FixedStar(hip) => {
            find_star_by_hip(hip).map_or_else(|| format!("HIP {}", hip.0), |s| s.name.clone())
        }
        BodyId::EclNut => "EclNut".to_owned(),
    }
}

// ============================================================================
// A. Body Longitude Crossings, Solar Terms, Seasons, Returns, Human Design
// ============================================================================

/// Result of solving for the exact epoch when a body reaches a target ecliptic longitude (`R9.1`, `R9.2`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LongitudeCrossingEvent {
    /// Target body identifier.
    pub body: BodyId,
    /// Requested target ecliptic longitude in degrees `[0, 360)`.
    pub target_lon_deg: f64,
    /// Terrestrial Time (`TT`) Julian Date of the crossing.
    pub jd_tt: f64,
    /// Universal Time (`UT1`) Julian Date of the crossing.
    pub jd_ut1: f64,
    /// $\Delta T = \text{TT} - \text{UT1}$ in seconds at the crossing.
    pub delta_t_sec: f64,
    /// Evaluated ecliptic longitude at the crossing in degrees `[0, 360)`.
    pub actual_lon_deg: f64,
    /// Evaluated ecliptic latitude at the crossing in degrees `[-90, +90]`.
    pub lat_deg: f64,
    /// Radial distance in AU at the crossing.
    pub dist_au: f64,
    /// Radial distance in kilometers at the crossing (`dist_au * AU_KM`).
    pub dist_km: f64,
    /// Longitudinal angular velocity $\dot{\lambda}$ in degrees/day at the crossing.
    pub lon_speed_deg_day: f64,
    /// `true` if the body is in apparent retrograde motion ($\dot{\lambda} < 0$) at the crossing.
    pub is_retrograde: bool,
}

impl LongitudeCrossingEvent {
    fn from_calc_result(target_lon_deg: f64, res: &CalcResult) -> Self {
        Self {
            body: res.body,
            target_lon_deg,
            jd_tt: res.jd_tt,
            jd_ut1: res.jd_ut1,
            delta_t_sec: res.delta_t_sec,
            actual_lon_deg: res.lon_deg,
            lat_deg: res.lat_deg,
            dist_au: res.dist_au,
            dist_km: res.dist_km,
            lon_speed_deg_day: res.lon_speed_deg_day,
            is_retrograde: res.lon_speed_deg_day < 0.0,
        }
    }
}

/// Characteristic coarse search step (in days) for ecliptic longitude root-finding.
fn characteristic_longitude_step_days(body: BodyId) -> f64 {
    match body {
        BodyId::Planet(PlanetId::Moon) => 0.5,
        BodyId::Planet(
            PlanetId::Sun | PlanetId::Mercury | PlanetId::Venus | PlanetId::Earth,
        ) => 1.0,
        BodyId::LunarPoint(_) => 0.5,
        BodyId::Planet(_) | BodyId::MinorPlanet(_) | BodyId::Hypothetical(_) => 2.0,
        BodyId::FixedStar(_) => 10.0,
        BodyId::EclNut => 1.0,
    }
}

/// Refines a bracketed longitude crossing on `[t_a, t_b]` (where `fa * fb <= 0` and `|fa|, |fb| < 90°`)
/// using a safeguarded hybrid **Newton-Illinois** solver to `< 1e-9` degrees (`< 0.001 s` for Moon, `< 0.01 s` for Sun).
#[allow(clippy::too_many_arguments)]
fn refine_longitude_crossing(
    engine: &Engine,
    body: BodyId,
    target_lon_deg: f64,
    t_a: f64,
    res_a: CalcResult,
    t_b: f64,
    res_b: CalcResult,
    eval_opts: CalcOptions,
) -> Result<LongitudeCrossingEvent, EphemerisError> {
    let mut a = t_a;
    let mut b = t_b;
    let mut fa = wrap_pm180(res_a.lon_deg - target_lon_deg);
    let mut fb = wrap_pm180(res_b.lon_deg - target_lon_deg);
    let mut va = res_a.lon_speed_deg_day;
    let mut vb = res_b.lon_speed_deg_day;

    if fa.abs() < 1.0e-9 {
        return Ok(LongitudeCrossingEvent::from_calc_result(target_lon_deg, &res_a));
    }
    if fb.abs() < 1.0e-9 {
        return Ok(LongitudeCrossingEvent::from_calc_result(target_lon_deg, &res_b));
    }

    let mut best_res = if fa.abs() <= fb.abs() { res_a } else { res_b };
    let mut side = 0i8;

    for _ in 0..50 {
        let span = (b - a).abs();
        if span < 1.0e-10 {
            return Ok(LongitudeCrossingEvent::from_calc_result(
                target_lon_deg,
                &best_res,
            ));
        }

        // Choose safeguarded Newton step from closer endpoint if derivative is well-conditioned
        // and lands strictly inside the active bracket; otherwise use Illinois false position.
        let (x_curr, f_curr, v_curr) = if fa.abs() <= fb.abs() {
            (a, fa, va)
        } else {
            (b, fb, vb)
        };
        let t_min = a.min(b);
        let t_max = a.max(b);
        let pad = 0.02 * span;

        let c_newton = if v_curr.abs() > 1.0e-5 {
            Some(x_curr - f_curr / v_curr)
        } else {
            None
        };

        let c = match c_newton {
            Some(cn) if cn > t_min + pad && cn < t_max - pad && (cn - x_curr).abs() < 0.75 * span => {
                cn
            }
            _ => {
                let denom = fb - fa;
                let cf = if denom.abs() > 1.0e-15 {
                    (a * fb - b * fa) / denom
                } else {
                    0.5 * (a + b)
                };
                cf.clamp(t_min + 0.01 * span, t_max - 0.01 * span)
            }
        };

        let res_c = engine.calc(body, TimeInput::from_tt_jd(c), eval_opts)?;
        let fc = wrap_pm180(res_c.lon_deg - target_lon_deg);
        let vc = res_c.lon_speed_deg_day;

        if fc.abs() < fa.abs().min(fb.abs()) {
            best_res = res_c.clone();
        }
        if fc.abs() < 1.0e-9 {
            return Ok(LongitudeCrossingEvent::from_calc_result(target_lon_deg, &res_c));
        }

        if (fc > 0.0) == (fb > 0.0) {
            b = c;
            fb = fc;
            vb = vc;
            if side == -1 {
                fa *= 0.5;
            }
            side = -1;
        } else {
            a = c;
            fa = fc;
            va = vc;
            if side == 1 {
                fb *= 0.5;
            }
            side = 1;
        }
    }

    if wrap_pm180(best_res.lon_deg - target_lon_deg).abs() < 1.0e-6 {
        return Ok(LongitudeCrossingEvent::from_calc_result(
            target_lon_deg,
            &best_res,
        ));
    }

    Err(EphemerisError::NonConvergence {
        solver: "refine_longitude_crossing".to_owned(),
        iterations: 50,
        residual: (b - a).abs(),
    })
}

/// Inspects a single step `[t_prev, t_curr]` (which may be forward `t_curr > t_prev` or backward `t_curr < t_prev`)
/// and returns all crossings of `target` in directional order. If a retrograde station (`speed_prev * speed_curr <= 0`)
/// occurs inside `[t_prev, t_curr]` near `target`, subdivides into 8 sub-steps so grazing double crossings are never missed.
#[allow(clippy::too_many_arguments)]
fn collect_crossings_in_step(
    engine: &Engine,
    body: BodyId,
    target: f64,
    t_prev: f64,
    res_prev: &CalcResult,
    t_curr: f64,
    res_curr: &CalcResult,
    eval_opts: CalcOptions,
    out: &mut Vec<LongitudeCrossingEvent>,
) -> Result<(), EphemerisError> {
    let f_prev = wrap_pm180(res_prev.lon_deg - target);
    let f_curr = wrap_pm180(res_curr.lon_deg - target);
    let v_prev = res_prev.lon_speed_deg_day;
    let v_curr = res_curr.lon_speed_deg_day;

    // If velocity changes sign inside this step and we are within 20° of target, subdivide into 8 sub-steps
    // and split the exact stationary sub-interval at v(t_station) = 0 so that even grazing double crossings
    // within 1e-8° of the station extremum are bracketed on two strictly monotonic halves.
    if v_prev * v_curr <= 0.0
        && (f_prev.abs() < 20.0 || f_curr.abs() < 20.0)
        && (t_curr - t_prev).abs() > 1.0e-7
    {
        let sub_n = 8usize;
        let mut sub_t_prev = t_prev;
        let mut sub_res_prev = res_prev.clone();
        let mut sub_f_prev = f_prev;

        for k in 1..=sub_n {
            let sub_t_curr = if k == sub_n {
                t_curr
            } else {
                t_prev + (t_curr - t_prev) * (k as f64 / sub_n as f64)
            };
            let sub_res_curr = if k == sub_n {
                res_curr.clone()
            } else {
                engine.calc(body, TimeInput::from_tt_jd(sub_t_curr), eval_opts)?
            };
            let sub_f_curr = wrap_pm180(sub_res_curr.lon_deg - target);

            let sub_v_prev = sub_res_prev.lon_speed_deg_day;
            let sub_v_curr = sub_res_curr.lon_speed_deg_day;
            if sub_v_prev * sub_v_curr <= 0.0 && (sub_t_curr - sub_t_prev).abs() > 1.0e-8 {
                // Locate the exact station extremum t_stat where lon_speed_deg_day == 0
                let mut a = sub_t_prev;
                let mut b = sub_t_curr;
                let mut va = sub_v_prev;
                let mut vb = sub_v_curr;
                let mut stat_res = sub_res_prev.clone();
                let mut side = 0i8;
                for _ in 0..35 {
                    let width = (b - a).abs();
                    if width <= 1.0e-9 || vb.abs() < 1.0e-11 {
                        break;
                    }
                    let denom = vb - va;
                    let mut c = if denom.abs() > 1.0e-14 {
                        (a * vb - b * va) / denom
                    } else {
                        0.5 * (a + b)
                    };
                    let lo = a.min(b);
                    let hi = a.max(b);
                    let margin = 0.05 * width;
                    if c <= lo + margin || c >= hi - margin {
                        c = 0.5 * (a + b);
                    }
                    let rc = engine.calc(body, TimeInput::from_tt_jd(c), eval_opts)?;
                    let vc = rc.lon_speed_deg_day;
                    stat_res = rc;
                    if vc.abs() < 1.0e-11 {
                        a = c;
                        b = c;
                        break;
                    }
                    if (vc > 0.0) == (vb > 0.0) {
                        b = c;
                        vb = vc;
                        if side == -1 {
                            va *= 0.5;
                        }
                        side = -1;
                    } else {
                        a = c;
                        va = vc;
                        if side == 1 {
                            vb *= 0.5;
                        }
                        side = 1;
                    }
                }
                let t_stat = 0.5 * (a + b);
                if (t_stat - stat_res.jd_tt).abs() > 1.0e-9 {
                    stat_res = engine.calc(body, TimeInput::from_tt_jd(t_stat), eval_opts)?;
                }
                let f_stat = wrap_pm180(stat_res.lon_deg - target);

                // First monotonic half: [sub_t_prev, t_stat]
                if sub_f_prev.abs() < 90.0 && f_stat.abs() < 90.0 && sub_f_prev * f_stat <= 0.0 {
                    let ev = refine_longitude_crossing(
                        engine,
                        body,
                        target,
                        sub_t_prev,
                        sub_res_prev.clone(),
                        t_stat,
                        stat_res.clone(),
                        eval_opts,
                    )?;
                    if out
                        .last()
                        .map_or(true, |last| (ev.jd_tt - last.jd_tt).abs() > 1.0e-6)
                    {
                        out.push(ev);
                    }
                }
                // Second monotonic half: [t_stat, sub_t_curr]
                if f_stat.abs() < 90.0 && sub_f_curr.abs() < 90.0 && f_stat * sub_f_curr <= 0.0 {
                    let ev = refine_longitude_crossing(
                        engine,
                        body,
                        target,
                        t_stat,
                        stat_res,
                        sub_t_curr,
                        sub_res_curr.clone(),
                        eval_opts,
                    )?;
                    if out
                        .last()
                        .map_or(true, |last| (ev.jd_tt - last.jd_tt).abs() > 1.0e-6)
                    {
                        out.push(ev);
                    }
                }
            } else if sub_f_prev.abs() < 90.0
                && sub_f_curr.abs() < 90.0
                && sub_f_prev * sub_f_curr <= 0.0
            {
                let ev = refine_longitude_crossing(
                    engine,
                    body,
                    target,
                    sub_t_prev,
                    sub_res_prev.clone(),
                    sub_t_curr,
                    sub_res_curr.clone(),
                    eval_opts,
                )?;
                if out
                    .last()
                    .map_or(true, |last| (ev.jd_tt - last.jd_tt).abs() > 1.0e-6)
                {
                    out.push(ev);
                }
            }

            sub_t_prev = sub_t_curr;
            sub_res_prev = sub_res_curr;
            sub_f_prev = sub_f_curr;
        }
        return Ok(());
    }

    if f_prev.abs() < 90.0 && f_curr.abs() < 90.0 && f_prev * f_curr <= 0.0 {
        let ev = refine_longitude_crossing(
            engine,
            body,
            target,
            t_prev,
            res_prev.clone(),
            t_curr,
            res_curr.clone(),
            eval_opts,
        )?;
        if out
            .last()
            .map_or(true, |last| (ev.jd_tt - last.jd_tt).abs() > 1.0e-5)
        {
            out.push(ev);
        }
    }

    Ok(())
}

/// Finds the first instant starting from `start_input` (moving forward if `forward == true`,
/// or backward if `forward == false`, up to `max_search_days`) when `body`'s ecliptic longitude
/// reaches `target_lon_deg.rem_euclid(360.0)` under `options` (`R9.1`, `R9.2`).
pub fn search_body_longitude(
    engine: &Engine,
    body: BodyId,
    start_input: TimeInput,
    target_lon_deg: f64,
    forward: bool,
    max_search_days: f64,
    options: &CalcOptions,
) -> Result<LongitudeCrossingEvent, EphemerisError> {
    if body == BodyId::EclNut {
        return Err(EphemerisError::UnsupportedBody {
            body: "EclNut".to_owned(),
        });
    }
    if !target_lon_deg.is_finite()
        || !max_search_days.is_finite()
        || !(1.0e-6..=365_250.0).contains(&max_search_days)
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid search_body_longitude parameters: target_lon_deg={target_lon_deg}, max_search_days={max_search_days}"
            ),
        });
    }

    let target = target_lon_deg.rem_euclid(360.0);
    let mut eval_opts = *options;
    eval_opts.plane = Plane::Ecliptic;
    eval_opts.include_speed = true;

    let obs = match eval_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let (_, tt_start, _, _) = engine.time_table().resolve_scales(start_input, obs);
    let jd_tt_start = tt_start.jd_f64();

    let dir = if forward { 1.0 } else { -1.0 };
    let base_step = characteristic_longitude_step_days(body)
        .min(max_search_days * 0.25)
        .max(1.0e-3);
    let total_steps = ((max_search_days / base_step).ceil() as usize).max(2);

    let mut t_prev = jd_tt_start;
    let mut res_prev = engine.calc(body, TimeInput::from_tt_jd(t_prev), eval_opts)?;
    let mut found = Vec::new();

    for step in 1..=total_steps {
        let dt = ((step as f64) * base_step).min(max_search_days);
        let t_curr = jd_tt_start + dir * dt;
        let res_curr = engine.calc(body, TimeInput::from_tt_jd(t_curr), eval_opts)?;

        collect_crossings_in_step(
            engine,
            body,
            target,
            t_prev,
            &res_prev,
            t_curr,
            &res_curr,
            eval_opts,
            &mut found,
        )?;

        if let Some(first) = found.first().cloned() {
            return Ok(first);
        }

        t_prev = t_curr;
        res_prev = res_curr;
    }

    Err(EphemerisError::NonConvergence {
        solver: "search_body_longitude".to_owned(),
        iterations: total_steps as u32,
        residual: wrap_pm180(res_prev.lon_deg - target).abs(),
    })
}

/// Finds **all** crossings of `target_lon_deg.rem_euclid(360.0)` by `body` in `[start_input, end_input]`
/// in chronological order, including 3-pass retrograde triple returns (`R9.1`, `R9.2`).
pub fn body_longitude_crossings_in_range(
    engine: &Engine,
    body: BodyId,
    start_input: TimeInput,
    end_input: TimeInput,
    target_lon_deg: f64,
    options: &CalcOptions,
) -> Result<Vec<LongitudeCrossingEvent>, EphemerisError> {
    if body == BodyId::EclNut {
        return Err(EphemerisError::UnsupportedBody {
            body: "EclNut".to_owned(),
        });
    }
    if !target_lon_deg.is_finite() {
        return Err(EphemerisError::InvalidInput {
            message: format!("Invalid target_lon_deg: {target_lon_deg}"),
        });
    }

    let mut eval_opts = *options;
    eval_opts.plane = Plane::Ecliptic;
    eval_opts.include_speed = true;

    let obs = match eval_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let (_, tt_start, _, _) = engine.time_table().resolve_scales(start_input, obs);
    let (_, tt_end, _, _) = engine.time_table().resolve_scales(end_input, obs);
    let jd_tt_start = tt_start.jd_f64();
    let jd_tt_end = tt_end.jd_f64();

    if !jd_tt_start.is_finite() || !jd_tt_end.is_finite() || jd_tt_end < jd_tt_start {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid time range in body_longitude_crossings_in_range: start={jd_tt_start}, end={jd_tt_end}"
            ),
        });
    }
    if (jd_tt_end - jd_tt_start) < 1.0e-9 {
        return Ok(Vec::new());
    }

    let target = target_lon_deg.rem_euclid(360.0);
    let span = jd_tt_end - jd_tt_start;
    let step_days = characteristic_longitude_step_days(body)
        .min(span * 0.5)
        .max(1.0e-3);
    let total_steps = ((span / step_days).ceil() as usize).max(1);

    let mut t_prev = jd_tt_start;
    let mut res_prev = engine.calc(body, TimeInput::from_tt_jd(t_prev), eval_opts)?;
    let mut events = Vec::new();

    for step in 1..=total_steps {
        let t_curr = if step == total_steps {
            jd_tt_end
        } else {
            (jd_tt_start + (step as f64) * step_days).min(jd_tt_end)
        };
        let res_curr = engine.calc(body, TimeInput::from_tt_jd(t_curr), eval_opts)?;

        collect_crossings_in_step(
            engine,
            body,
            target,
            t_prev,
            &res_prev,
            t_curr,
            &res_curr,
            eval_opts,
            &mut events,
        )?;

        t_prev = t_curr;
        res_prev = res_curr;
    }

    Ok(events)
}

/// Chinese 24 Solar Term (*Jieqi*) event (`R12`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolarTermEvent {
    /// Solar term index `0..=23` (`0` = `0°` Chunfen / Vernal Equinox, `6` = `90°` Xiazhi / Summer Solstice,
    /// `12` = `180°` Qiufen / Autumnal Equinox, `18` = `270°` Dongzhi / Winter Solstice).
    pub term_index: u8,
    /// Exact target geocentric apparent tropical solar longitude in degrees (`0.0, 15.0, ..., 345.0`).
    pub target_lon_deg: f64,
    /// Standardized Pinyin name (`"Chunfen"`, `"Qingming"`, ..., `"Dongzhi"`, ..., `"Jingzhe"`).
    pub pinyin_name: String,
    /// Canonical Chinese name (`"春分"`, `"清明"`, ..., `"冬至"`, ..., `"惊蛰"`).
    pub chinese_name: String,
    /// Canonical English name (`"Vernal Equinox"`, ..., `"Minor Cold"`, `"Major Cold"`, `"Start of Spring"`, ..., `"Winter Solstice"`).
    pub english_name: String,
    /// `true` for even indices (`term_index % 2 == 0`), i.e., *Zhongqi* (principal terms at multiples of $30^\circ$).
    pub is_major_term: bool,
    /// Terrestrial Time (`TT`) Julian Date.
    pub jd_tt: f64,
    /// Universal Time (`UT1`) Julian Date.
    pub jd_ut1: f64,
    /// $\Delta T = \text{TT} - \text{UT1}$ in seconds.
    pub delta_t_sec: f64,
}

/// Finds all Chinese 24 Solar Terms (where **geocentric apparent tropical** solar longitude
/// $\lambda_\odot = k \cdot 15^\circ$) in `[start_input, end_input]` in chronological order (`R12`).
pub fn solar_terms_in_range(
    engine: &Engine,
    start_input: TimeInput,
    end_input: TimeInput,
) -> Result<Vec<SolarTermEvent>, EphemerisError> {
    let (_, tt_start, _, _) = engine.time_table().resolve_scales(start_input, None);
    let (_, tt_end, _, _) = engine.time_table().resolve_scales(end_input, None);
    let jd_tt_start = tt_start.jd_f64();
    let jd_tt_end = tt_end.jd_f64();

    if !jd_tt_start.is_finite() || !jd_tt_end.is_finite() || jd_tt_end < jd_tt_start {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid time range in solar_terms_in_range: start={jd_tt_start}, end={jd_tt_end}"
            ),
        });
    }
    if (jd_tt_end - jd_tt_start) < 1.0e-9 {
        return Ok(Vec::new());
    }

    let sun_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Ecliptic,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: true,
    };
    let sun = BodyId::Planet(PlanetId::Sun);

    // Coarse step of 5.0 days (~4.9° solar motion per step, strictly < 15° spacing between terms)
    let span = jd_tt_end - jd_tt_start;
    let step_days = 5.0_f64.min(span * 0.5).max(1.0e-3);
    let total_steps = ((span / step_days).ceil() as usize).max(1);

    let mut t_prev = jd_tt_start;
    let mut res_prev = engine.calc(sun, TimeInput::from_tt_jd(t_prev), sun_opts)?;
    let mut terms: Vec<SolarTermEvent> = Vec::new();

    for step in 1..=total_steps {
        let t_curr = if step == total_steps {
            jd_tt_end
        } else {
            (jd_tt_start + (step as f64) * step_days).min(jd_tt_end)
        };
        let res_curr = engine.calc(sun, TimeInput::from_tt_jd(t_curr), sun_opts)?;

        for (idx, &pinyin_name) in SOLAR_TERM_PINYIN_NAMES.iter().enumerate() {
            let term_index = idx as u8;
            let target_lon_deg = f64::from(term_index) * 15.0;
            let f_prev = wrap_pm180(res_prev.lon_deg - target_lon_deg);
            let f_curr = wrap_pm180(res_curr.lon_deg - target_lon_deg);

            if f_prev <= 0.0 && f_curr >= 0.0 && f_prev.abs() < 30.0 && f_curr.abs() < 30.0 {
                let crossing = refine_longitude_crossing(
                    engine,
                    sun,
                    target_lon_deg,
                    t_prev,
                    res_prev.clone(),
                    t_curr,
                    res_curr.clone(),
                    sun_opts,
                )?;
                if crossing.jd_tt >= jd_tt_start - 1.0e-9
                    && crossing.jd_tt <= jd_tt_end + 1.0e-9
                    && terms
                        .last()
                        .map_or(true, |last| (crossing.jd_tt - last.jd_tt).abs() > 1.0)
                {
                    let chinese_name = SOLAR_TERM_CHINESE_NAMES
                        .get(idx)
                        .copied()
                        .unwrap_or("")
                        .to_owned();
                    let english_name = SOLAR_TERM_ENGLISH_NAMES
                        .get(idx)
                        .copied()
                        .unwrap_or("")
                        .to_owned();
                    terms.push(SolarTermEvent {
                        term_index,
                        target_lon_deg,
                        pinyin_name: pinyin_name.to_owned(),
                        chinese_name,
                        english_name,
                        is_major_term: term_index % 2 == 0,
                        jd_tt: crossing.jd_tt,
                        jd_ut1: crossing.jd_ut1,
                        delta_t_sec: crossing.delta_t_sec,
                    });
                }
            }
        }

        t_prev = t_curr;
        res_prev = res_curr;
    }

    terms.sort_by(|a, b| a.jd_tt.total_cmp(&b.jd_tt));
    Ok(terms)
}

/// The 4 cardinal astronomical seasons of a calendar year (`0°`, `90°`, `180°`, `270°` apparent solar longitude).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SeasonsOfYear {
    /// Calendar year (astronomical numbering).
    pub year: i32,
    /// Vernal Equinox (`0°` Aries ingress) in `JD(UT1)`.
    pub vernal_equinox_jd_ut1: f64,
    /// Vernal Equinox (`0°` Aries ingress) in `JD(TT)`.
    pub vernal_equinox_jd_tt: f64,
    /// Summer Solstice (`90°` Cancer ingress) in `JD(UT1)`.
    pub summer_solstice_jd_ut1: f64,
    /// Summer Solstice (`90°` Cancer ingress) in `JD(TT)`.
    pub summer_solstice_jd_tt: f64,
    /// Autumnal Equinox (`180°` Libra ingress) in `JD(UT1)`.
    pub autumnal_equinox_jd_ut1: f64,
    /// Autumnal Equinox (`180°` Libra ingress) in `JD(TT)`.
    pub autumnal_equinox_jd_tt: f64,
    /// Winter Solstice (`270°` Capricorn ingress) in `JD(UT1)`.
    pub winter_solstice_jd_ut1: f64,
    /// Winter Solstice (`270°` Capricorn ingress) in `JD(TT)`.
    pub winter_solstice_jd_tt: f64,
}

/// Computes all 4 cardinal seasons (`0°`, `90°`, `180°`, `270°` geocentric apparent tropical solar longitude)
/// for `year`.
pub fn seasons_in_year(engine: &Engine, year: i32) -> Result<SeasonsOfYear, EphemerisError> {
    let sun_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Ecliptic,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: true,
    };
    let sun = BodyId::Planet(PlanetId::Sun);

    // Anchor search in early February to accommodate both Gregorian and historical Julian drift across [-3000, +3000]
    let feb_1 = julday_split(year, 2, 1, 0.0, CalendarKind::AstroAuto)?;
    let vernal = search_body_longitude(
        engine,
        sun,
        TimeInput::Tt(feb_1),
        0.0,
        true,
        90.0,
        &sun_opts,
    )?;
    let summer = search_body_longitude(
        engine,
        sun,
        TimeInput::from_tt_jd(vernal.jd_tt + 75.0),
        90.0,
        true,
        35.0,
        &sun_opts,
    )?;
    let autumn = search_body_longitude(
        engine,
        sun,
        TimeInput::from_tt_jd(summer.jd_tt + 75.0),
        180.0,
        true,
        35.0,
        &sun_opts,
    )?;
    let winter = search_body_longitude(
        engine,
        sun,
        TimeInput::from_tt_jd(autumn.jd_tt + 75.0),
        270.0,
        true,
        35.0,
        &sun_opts,
    )?;

    Ok(SeasonsOfYear {
        year,
        vernal_equinox_jd_ut1: vernal.jd_ut1,
        vernal_equinox_jd_tt: vernal.jd_tt,
        summer_solstice_jd_ut1: summer.jd_ut1,
        summer_solstice_jd_tt: summer.jd_tt,
        autumnal_equinox_jd_ut1: autumn.jd_ut1,
        autumnal_equinox_jd_tt: autumn.jd_tt,
        winter_solstice_jd_ut1: winter.jd_ut1,
        winter_solstice_jd_tt: winter.jd_tt,
    })
}

/// Computes the exact Human Design "Design" prenatal epoch (`R9.2`): evaluates the Sun's
/// geocentric apparent tropical longitude $\lambda_{\text{birth}}$ at `birth_input` and solves
/// backward (~88–89 days) for the exact instant when $\lambda_\odot(t_{\text{design}}) = (\lambda_{\text{birth}} - 88^\circ) \pmod{360^\circ}$.
pub fn human_design_epoch(
    engine: &Engine,
    birth_input: TimeInput,
) -> Result<LongitudeCrossingEvent, EphemerisError> {
    let sun_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Ecliptic,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: true,
    };
    let sun = BodyId::Planet(PlanetId::Sun);
    let sun_birth = engine.calc(sun, birth_input, sun_opts)?;
    let target_lon_deg = (sun_birth.lon_deg - 88.0).rem_euclid(360.0);

    // Jump back 84 days first (~83°–85° of solar motion) and search backward up to 15 days
    search_body_longitude(
        engine,
        sun,
        TimeInput::from_tt_jd(sun_birth.jd_tt - 84.0),
        target_lon_deg,
        false,
        15.0,
        &sun_opts,
    )
}

/// Exact planetary or lunar return event (Tropical or Sidereal across all 44 ayanamsa models).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlanetaryReturnEvent {
    /// Returning body identifier.
    pub body: BodyId,
    /// Natal ecliptic longitude $\lambda_{\text{natal}} \in [0^\circ, 360^\circ)$ evaluated under `options`.
    pub natal_lon_deg: f64,
    /// Universal Time (`UT1`) Julian Date of the return pass.
    pub jd_ut1: f64,
    /// Terrestrial Time (`TT`) Julian Date of the return pass.
    pub jd_tt: f64,
    /// Evaluated ecliptic longitude at the return pass in degrees `[0, 360)`.
    pub actual_lon_deg: f64,
    /// Longitudinal angular velocity $\dot{\lambda}$ in degrees/day at the return pass.
    pub lon_speed_deg_day: f64,
    /// `true` if the body is in apparent retrograde motion ($\dot{\lambda} < 0$) at the return pass.
    pub is_retrograde: bool,
    /// Pass index within a retrograde loop return series (`1` for single/first pass, `2` for retrograde
    /// second pass, `3` for direct third pass; resets to `1` on subsequent orbital returns).
    pub pass_index: u32,
}

/// Characteristic geocentric ecliptic return period in days for grouping multi-pass retrograde returns
/// and sizing search windows.
#[must_use]
pub fn characteristic_orbital_period_days(body: BodyId) -> f64 {
    match body {
        BodyId::Planet(PlanetId::Moon) => 27.321_661,
        BodyId::Planet(
            PlanetId::Sun
            | PlanetId::Mercury
            | PlanetId::Venus
            | PlanetId::Earth
            | PlanetId::EarthMoonBarycenter,
        ) => 365.256_363,
        BodyId::Planet(PlanetId::Mars) => 686.979_59,
        BodyId::Planet(PlanetId::Jupiter) => 4_332.589,
        BodyId::Planet(PlanetId::Saturn) => 10_759.22,
        BodyId::Planet(PlanetId::Uranus) => 30_688.5,
        BodyId::Planet(PlanetId::Neptune) => 60_195.0,
        BodyId::Planet(PlanetId::Pluto) => 90_560.0,
        BodyId::Planet(PlanetId::SolarSystemBarycenter) => 365.256_363,
        BodyId::LunarPoint(
            LunarPointId::MeanNode
            | LunarPointId::TrueNode
            | LunarPointId::MeanSouthNode
            | LunarPointId::TrueSouthNode,
        ) => 6_798.383,
        BodyId::LunarPoint(_) => 3_232.605,
        BodyId::MinorPlanet(mpc) => match mpc.0 {
            2060 => 18_520.0,   // Chiron (~50.7 yr)
            5145 => 33_600.0,   // Pholus (~92 yr)
            10199 => 23_010.0,  // Chariklo (~63 yr)
            136199 => 204_175.0, // Eris (~559 yr)
            90377 => 365_250.0, // Sedna (clamped)
            _ => 1_680.0,       // Main-belt asteroid typical ~4.6 yr
        },
        BodyId::Hypothetical(_) => 36_525.0,
        BodyId::FixedStar(_) => 365_250.0,
        BodyId::EclNut => 365.256_363,
    }
}

/// Finds the first exact return of `body` to its natal ecliptic longitude after `search_after_input`
/// under `options` (supporting both Tropical and Sidereal returns in all ayanamsa models).
pub fn planetary_return(
    engine: &Engine,
    body: BodyId,
    natal_input: TimeInput,
    search_after_input: TimeInput,
    options: CalcOptions,
) -> Result<PlanetaryReturnEvent, EphemerisError> {
    if body == BodyId::EclNut {
        return Err(EphemerisError::UnsupportedBody {
            body: "EclNut".to_owned(),
        });
    }
    let mut eval_opts = options;
    eval_opts.plane = Plane::Ecliptic;
    eval_opts.include_speed = true;

    let natal_res = engine.calc(body, natal_input, eval_opts)?;
    let max_search_days = (characteristic_orbital_period_days(body) * 1.5).clamp(400.0, 365_250.0);
    let crossing = search_body_longitude(
        engine,
        body,
        search_after_input,
        natal_res.lon_deg,
        true,
        max_search_days,
        &eval_opts,
    )?;

    Ok(PlanetaryReturnEvent {
        body,
        natal_lon_deg: natal_res.lon_deg,
        jd_ut1: crossing.jd_ut1,
        jd_tt: crossing.jd_tt,
        actual_lon_deg: crossing.actual_lon_deg,
        lon_speed_deg_day: crossing.lon_speed_deg_day,
        is_retrograde: crossing.is_retrograde,
        pass_index: 1,
    })
}

/// Finds all returns of `body` to its natal ecliptic longitude within `[start_input, end_input]`
/// under `options`, numbering retrograde loop passes sequentially (`pass_index: 1, 2, 3`) and
/// resetting `pass_index` to `1` whenever more than half of `body`'s orbital period has elapsed
/// since the previous pass.
pub fn planetary_returns_in_range(
    engine: &Engine,
    body: BodyId,
    natal_input: TimeInput,
    start_input: TimeInput,
    end_input: TimeInput,
    options: CalcOptions,
) -> Result<Vec<PlanetaryReturnEvent>, EphemerisError> {
    if body == BodyId::EclNut {
        return Err(EphemerisError::UnsupportedBody {
            body: "EclNut".to_owned(),
        });
    }
    let mut eval_opts = options;
    eval_opts.plane = Plane::Ecliptic;
    eval_opts.include_speed = true;

    let natal_res = engine.calc(body, natal_input, eval_opts)?;
    let crossings = body_longitude_crossings_in_range(
        engine,
        body,
        start_input,
        end_input,
        natal_res.lon_deg,
        &eval_opts,
    )?;

    let half_period_days = 0.5 * characteristic_orbital_period_days(body);
    let mut out = Vec::with_capacity(crossings.len());
    let mut prev_jd_tt: Option<f64> = None;
    let mut current_pass: u32 = 0;

    for c in crossings {
        if let Some(prev) = prev_jd_tt {
            if (c.jd_tt - prev) > half_period_days {
                current_pass = 1;
            } else {
                current_pass += 1;
            }
        } else {
            current_pass = 1;
        }
        prev_jd_tt = Some(c.jd_tt);

        out.push(PlanetaryReturnEvent {
            body,
            natal_lon_deg: natal_res.lon_deg,
            jd_ut1: c.jd_ut1,
            jd_tt: c.jd_tt,
            actual_lon_deg: c.actual_lon_deg,
            lon_speed_deg_day: c.lon_speed_deg_day,
            is_retrograde: c.is_retrograde,
            pass_index: current_pass,
        });
    }

    Ok(out)
}

// ============================================================================
// B. Lunar Phases, Moon Age, and Kerykeion 8 / 28 Phase Summary
// ============================================================================

/// The 4 principal lunar phases (`0°`, `90°`, `180°`, `270°` geocentric apparent lunar-solar elongation).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum LunarPhaseKind {
    /// New Moon ($\lambda_{\text{Moon}} - \lambda_\odot = 0^\circ$).
    NewMoon,
    /// First Quarter ($\lambda_{\text{Moon}} - \lambda_\odot = 90^\circ$).
    FirstQuarter,
    /// Full Moon ($\lambda_{\text{Moon}} - \lambda_\odot = 180^\circ$).
    FullMoon,
    /// Last (Third) Quarter ($\lambda_{\text{Moon}} - \lambda_\odot = 270^\circ$).
    LastQuarter,
}

impl LunarPhaseKind {
    /// Target lunar-solar elongation in degrees (`0.0`, `90.0`, `180.0`, `270.0`).
    #[must_use]
    pub const fn target_elongation_deg(self) -> f64 {
        match self {
            Self::NewMoon => 0.0,
            Self::FirstQuarter => 90.0,
            Self::FullMoon => 180.0,
            Self::LastQuarter => 270.0,
        }
    }

    /// Maps an elongation in degrees to [`LunarPhaseKind`] if it matches one of the 4 principal quarters.
    #[must_use]
    pub fn from_elongation_deg(deg: f64) -> Option<Self> {
        let norm = deg.rem_euclid(360.0);
        if norm.abs() < 1.0e-7 || (norm - 360.0).abs() < 1.0e-7 {
            Some(Self::NewMoon)
        } else if (norm - 90.0).abs() < 1.0e-7 {
            Some(Self::FirstQuarter)
        } else if (norm - 180.0).abs() < 1.0e-7 {
            Some(Self::FullMoon)
        } else if (norm - 270.0).abs() < 1.0e-7 {
            Some(Self::LastQuarter)
        } else {
            None
        }
    }
}

/// Result of solving for a principal or custom lunar phase elongation (`R9.2`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LunarPhaseEvent {
    /// Principal lunar phase kind (`Some` for $0^\circ, 90^\circ, 180^\circ, 270^\circ$, `None` for custom elongation).
    pub kind: Option<LunarPhaseKind>,
    /// Target geocentric apparent tropical elongation $(\lambda_{\text{Moon}} - \lambda_\odot) \pmod{360^\circ}$.
    pub target_elongation_deg: f64,
    /// Terrestrial Time (`TT`) Julian Date.
    pub jd_tt: f64,
    /// Universal Time (`UT1`) Julian Date.
    pub jd_ut1: f64,
    /// $\Delta T = \text{TT} - \text{UT1}$ in seconds.
    pub delta_t_sec: f64,
    /// Geocentric apparent tropical lunar longitude $\lambda_{\text{Moon}}$ in degrees `[0, 360)`.
    pub moon_lon_deg: f64,
    /// Geocentric apparent tropical solar longitude $\lambda_\odot$ in degrees `[0, 360)`.
    pub sun_lon_deg: f64,
    /// Lunar geocentric distance in AU.
    pub moon_dist_au: f64,
    /// Lunar geocentric distance in kilometers (`moon_dist_au * AU_KM`, `R2.2`).
    pub moon_dist_km: f64,
}

/// Evaluates `(elongation_deg, rel_speed_deg_day, moon_calc, sun_calc)` at `jd_tt`.
fn eval_lunar_elongation_at_tt(
    engine: &Engine,
    jd_tt: f64,
) -> Result<(f64, f64, CalcResult, CalcResult), EphemerisError> {
    let opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Ecliptic,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: true,
    };
    let moment = engine.at(TimeInput::from_tt_jd(jd_tt), None)?;
    let moon = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Moon), opts)?;
    let sun = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Sun), opts)?;
    let elong = (moon.lon_deg - sun.lon_deg).rem_euclid(360.0);
    let rel_speed = moon.lon_speed_deg_day - sun.lon_speed_deg_day;
    Ok((elong, rel_speed, moon, sun))
}

/// Refines a bracketed lunar elongation crossing on `[t_a, t_b]` via hybrid Newton-Illinois.
#[allow(clippy::too_many_arguments)]
fn refine_lunar_elongation(
    engine: &Engine,
    target_elong_deg: f64,
    t_a: f64,
    elong_a: f64,
    speed_a: f64,
    moon_a: CalcResult,
    sun_a: CalcResult,
    t_b: f64,
    elong_b: f64,
    speed_b: f64,
    moon_b: CalcResult,
    sun_b: CalcResult,
) -> Result<LunarPhaseEvent, EphemerisError> {
    let mut a = t_a;
    let mut b = t_b;
    let mut fa = wrap_pm180(elong_a - target_elong_deg);
    let mut fb = wrap_pm180(elong_b - target_elong_deg);
    let mut va = speed_a;
    let mut vb = speed_b;

    let make_event = |moon: &CalcResult, sun: &CalcResult| LunarPhaseEvent {
        kind: LunarPhaseKind::from_elongation_deg(target_elong_deg),
        target_elongation_deg: target_elong_deg,
        jd_tt: moon.jd_tt,
        jd_ut1: moon.jd_ut1,
        delta_t_sec: moon.delta_t_sec,
        moon_lon_deg: moon.lon_deg,
        sun_lon_deg: sun.lon_deg,
        moon_dist_au: moon.dist_au,
        moon_dist_km: moon.dist_km,
    };

    if fa.abs() < 1.0e-9 {
        return Ok(make_event(&moon_a, &sun_a));
    }
    if fb.abs() < 1.0e-9 {
        return Ok(make_event(&moon_b, &sun_b));
    }

    let mut best_pair = if fa.abs() <= fb.abs() {
        (moon_a, sun_a)
    } else {
        (moon_b, sun_b)
    };
    let mut side = 0i8;

    for _ in 0..40 {
        let span = (b - a).abs();
        if span < 1.0e-10 {
            return Ok(make_event(&best_pair.0, &best_pair.1));
        }

        let (x_curr, f_curr, v_curr) = if fa.abs() <= fb.abs() {
            (a, fa, va)
        } else {
            (b, fb, vb)
        };
        let t_min = a.min(b);
        let t_max = a.max(b);
        let pad = 0.02 * span;

        let c_newton = if v_curr.abs() > 1.0 {
            Some(x_curr - f_curr / v_curr)
        } else {
            None
        };
        let c = match c_newton {
            Some(cn) if cn > t_min + pad && cn < t_max - pad => cn,
            _ => {
                let denom = fb - fa;
                if denom.abs() > 1.0e-15 {
                    ((a * fb - b * fa) / denom).clamp(t_min + pad, t_max - pad)
                } else {
                    0.5 * (a + b)
                }
            }
        };

        let (elong_c, vc, moon_c, sun_c) = eval_lunar_elongation_at_tt(engine, c)?;
        let fc = wrap_pm180(elong_c - target_elong_deg);

        if fc.abs() < fa.abs().min(fb.abs()) {
            best_pair = (moon_c.clone(), sun_c.clone());
        }
        if fc.abs() < 1.0e-9 {
            return Ok(make_event(&moon_c, &sun_c));
        }

        if (fc > 0.0) == (fb > 0.0) {
            b = c;
            fb = fc;
            vb = vc;
            if side == -1 {
                fa *= 0.5;
            }
            side = -1;
        } else {
            a = c;
            fa = fc;
            va = vc;
            if side == 1 {
                fb *= 0.5;
            }
            side = 1;
        }
    }

    Ok(make_event(&best_pair.0, &best_pair.1))
}

/// Finds the next (`forward == true`) or previous (`forward == false`) instant from `start_input`
/// when geocentric apparent tropical lunar-solar elongation $(\lambda_{\text{Moon}} - \lambda_\odot) \pmod{360^\circ}$
/// reaches `target_elongation_deg.rem_euclid(360.0)` to $< 10^{-8^\circ}$ (`< 0.002 s`).
pub fn search_lunar_phase(
    engine: &Engine,
    start_input: TimeInput,
    target_elongation_deg: f64,
    forward: bool,
) -> Result<LunarPhaseEvent, EphemerisError> {
    if !target_elongation_deg.is_finite() {
        return Err(EphemerisError::InvalidInput {
            message: format!("Invalid target_elongation_deg: {target_elongation_deg}"),
        });
    }
    let target = target_elongation_deg.rem_euclid(360.0);
    let (_, tt_start, _, _) = engine.time_table().resolve_scales(start_input, None);
    let jd_tt_start = tt_start.jd_f64();

    let dir = if forward { 1.0 } else { -1.0 };
    let step_days = 1.0_f64;
    let max_steps = 35usize; // synodic month is ~29.53 days

    let mut t_prev = jd_tt_start;
    let (mut elong_prev, mut speed_prev, mut moon_prev, mut sun_prev) =
        eval_lunar_elongation_at_tt(engine, t_prev)?;

    let delta_deg = if forward {
        (target - elong_prev).rem_euclid(360.0)
    } else {
        (elong_prev - target).rem_euclid(360.0)
    };
    if delta_deg > 35.0 {
        let jump_days = (delta_deg / 12.190_749_1 - 2.2).max(0.0);
        if jump_days > 0.0 {
            t_prev = jd_tt_start + dir * jump_days;
            (elong_prev, speed_prev, moon_prev, sun_prev) =
                eval_lunar_elongation_at_tt(engine, t_prev)?;
        }
    }
    let scan_origin = t_prev;

    for step in 1..=max_steps {
        let t_curr = scan_origin + dir * (step as f64) * step_days;
        let (elong_curr, speed_curr, moon_curr, sun_curr) =
            eval_lunar_elongation_at_tt(engine, t_curr)?;

        let f_prev = wrap_pm180(elong_prev - target);
        let f_curr = wrap_pm180(elong_curr - target);

        if f_prev.abs() < 60.0 && f_curr.abs() < 60.0 && f_prev * f_curr <= 0.0 {
            return refine_lunar_elongation(
                engine,
                target,
                t_prev,
                elong_prev,
                speed_prev,
                moon_prev,
                sun_prev,
                t_curr,
                elong_curr,
                speed_curr,
                moon_curr,
                sun_curr,
            );
        }

        t_prev = t_curr;
        elong_prev = elong_curr;
        speed_prev = speed_curr;
        moon_prev = moon_curr;
        sun_prev = sun_curr;
    }

    Err(EphemerisError::NonConvergence {
        solver: "search_lunar_phase".to_owned(),
        iterations: max_steps as u32,
        residual: wrap_pm180(elong_prev - target).abs(),
    })
}

/// Returns all 4 principal lunar phases (`NewMoon`, `FirstQuarter`, `FullMoon`, `LastQuarter`)
/// occurring within `[start_input, end_input]` in chronological order (`R9.2`).
pub fn lunar_phases_in_range(
    engine: &Engine,
    start_input: TimeInput,
    end_input: TimeInput,
) -> Result<Vec<LunarPhaseEvent>, EphemerisError> {
    let (_, tt_start, _, _) = engine.time_table().resolve_scales(start_input, None);
    let (_, tt_end, _, _) = engine.time_table().resolve_scales(end_input, None);
    let jd_tt_start = tt_start.jd_f64();
    let jd_tt_end = tt_end.jd_f64();

    if !jd_tt_start.is_finite() || !jd_tt_end.is_finite() || jd_tt_end < jd_tt_start {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid time range in lunar_phases_in_range: start={jd_tt_start}, end={jd_tt_end}"
            ),
        });
    }
    if (jd_tt_end - jd_tt_start) < 1.0e-9 {
        return Ok(Vec::new());
    }

    let span = jd_tt_end - jd_tt_start;
    let step_days = 3.0_f64.min(span * 0.5).max(1.0e-3);
    let total_steps = ((span / step_days).ceil() as usize).max(1);

    let mut t_prev = jd_tt_start;
    let (mut elong_prev, mut speed_prev, mut moon_prev, mut sun_prev) =
        eval_lunar_elongation_at_tt(engine, t_prev)?;
    let mut events: Vec<LunarPhaseEvent> = Vec::new();

    for step in 1..=total_steps {
        let t_curr = if step == total_steps {
            jd_tt_end
        } else {
            (jd_tt_start + (step as f64) * step_days).min(jd_tt_end)
        };
        let (elong_curr, speed_curr, moon_curr, sun_curr) =
            eval_lunar_elongation_at_tt(engine, t_curr)?;

        for kind in [
            LunarPhaseKind::NewMoon,
            LunarPhaseKind::FirstQuarter,
            LunarPhaseKind::FullMoon,
            LunarPhaseKind::LastQuarter,
        ] {
            let target = kind.target_elongation_deg();
            let f_prev = wrap_pm180(elong_prev - target);
            let f_curr = wrap_pm180(elong_curr - target);

            if f_prev <= 0.0 && f_curr >= 0.0 && f_prev.abs() < 45.0 && f_curr.abs() < 45.0 {
                let ev = refine_lunar_elongation(
                    engine,
                    target,
                    t_prev,
                    elong_prev,
                    speed_prev,
                    moon_prev.clone(),
                    sun_prev.clone(),
                    t_curr,
                    elong_curr,
                    speed_curr,
                    moon_curr.clone(),
                    sun_curr.clone(),
                )?;
                if ev.jd_tt >= jd_tt_start - 1.0e-9
                    && ev.jd_tt <= jd_tt_end + 1.0e-9
                    && events
                        .last()
                        .map_or(true, |last| (ev.jd_tt - last.jd_tt).abs() > 1.0)
                {
                    events.push(ev);
                }
            }
        }

        t_prev = t_curr;
        elong_prev = elong_curr;
        speed_prev = speed_curr;
        moon_prev = moon_curr;
        sun_prev = sun_curr;
    }

    events.sort_by(|a, b| a.jd_tt.total_cmp(&b.jd_tt));
    Ok(events)
}

/// Complete lunar phase, illumination, distance, age, and surrounding quarters summary (`R0`, `R2.2`, `R9.2`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LunarPhaseSummary {
    /// Terrestrial Time (`TT`) Julian Date.
    pub jd_tt: f64,
    /// Universal Time (`UT1`) Julian Date.
    pub jd_ut1: f64,
    /// Geocentric apparent tropical lunar-solar elongation $(\lambda_{\text{Moon}} - \lambda_\odot) \pmod{360^\circ}$ in `[0, 360)`.
    pub elongation_deg: f64,
    /// Illuminated fraction of the lunar disc $k = (1 + \cos\alpha)/2 \in [0, 1]$.
    pub illumination_fraction: f64,
    /// 3D Sun–Moon–Observer phase angle $\alpha$ in degrees `[0, 180]`.
    pub phase_angle_deg: f64,
    /// Geocentric apparent tropical lunar longitude in degrees `[0, 360)`.
    pub moon_lon_deg: f64,
    /// Geocentric apparent tropical lunar latitude in degrees `[-90, +90]`.
    pub moon_lat_deg: f64,
    /// Lunar geocentric distance in AU.
    pub moon_dist_au: f64,
    /// Lunar geocentric distance in kilometers (`moon_dist_au * AU_KM`, `R2.2`).
    pub moon_dist_km: f64,
    /// Apparent lunar equatorial disc diameter in degrees.
    pub apparent_diameter_deg: f64,
    /// Apparent visual magnitude $V$ of the Moon.
    pub magnitude_v: f64,
    /// Exact lunar age in days since the preceding New Moon (`jd_ut1 - prev_new_moon_jd_ut1`).
    pub age_days: f64,
    /// Fractional progress through the current synodic month `(jd_ut1 - prev_new_moon_jd_ut1) / (next_new_moon_jd_ut1 - prev_new_moon_jd_ut1)`.
    pub synodic_fraction: f64,
    /// 8-phase octant index `0..=7` (`((elongation_deg / 45.0).floor() as u8) % 8`).
    pub phase_8_index: u8,
    /// Canonical name of the 8-phase octant (`"New Moon"`, `"Waxing Crescent"`, ..., `"Waning Crescent"`).
    pub phase_8_name: String,
    /// 28-phase index `1..=28` (`(((elongation_deg / (360.0 / 28.0)).floor() as u8) % 28) + 1`),
    /// replacing `kerykeion.calculate_moon_phase`.
    pub phase_28_index: u8,
    /// Preceding New Moon timestamp in `JD(UT1)`.
    pub prev_new_moon_jd_ut1: f64,
    /// Following New Moon timestamp in `JD(UT1)`.
    pub next_new_moon_jd_ut1: f64,
    /// Next First Quarter timestamp in `JD(UT1)`.
    pub next_first_quarter_jd_ut1: f64,
    /// Next Full Moon timestamp in `JD(UT1)`.
    pub next_full_moon_jd_ut1: f64,
    /// Next Last Quarter timestamp in `JD(UT1)`.
    pub next_last_quarter_jd_ut1: f64,
}

/// Computes the complete [`LunarPhaseSummary`] at `input` (`R0`, `R2.2`, `R9.2`).
pub fn lunar_phase_summary(
    engine: &Engine,
    input: TimeInput,
) -> Result<LunarPhaseSummary, EphemerisError> {
    let opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Ecliptic,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: true,
    };
    let moment = engine.at(input, None)?;
    let moon = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Moon), opts)?;
    let sun = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Sun), opts)?;
    let pheno = engine.phenomena_in_moment(&moment, BodyId::Planet(PlanetId::Moon), opts)?;

    let elongation_deg = (moon.lon_deg - sun.lon_deg).rem_euclid(360.0);
    let phase_8_index = ((elongation_deg / 45.0).floor() as u8) % 8;
    let phase_8_name = LUNAR_PHASE_8_NAMES
        .get(usize::from(phase_8_index))
        .copied()
        .unwrap_or("New Moon")
        .to_owned();
    let phase_28_index = (((elongation_deg / (360.0 / 28.0)).floor() as u8) % 28) + 1;

    let prev_nm = search_lunar_phase(engine, input, 0.0, false)?;
    let next_nm = search_lunar_phase(engine, input, 0.0, true)?;
    let next_fq = search_lunar_phase(engine, input, 90.0, true)?;
    let next_fm = search_lunar_phase(engine, input, 180.0, true)?;
    let next_lq = search_lunar_phase(engine, input, 270.0, true)?;

    let age_days = (moon.jd_ut1 - prev_nm.jd_ut1).max(0.0);
    let synodic_span = (next_nm.jd_ut1 - prev_nm.jd_ut1).max(1.0e-9);
    let synodic_fraction = (age_days / synodic_span).clamp(0.0, 1.0);

    Ok(LunarPhaseSummary {
        jd_tt: moon.jd_tt,
        jd_ut1: moon.jd_ut1,
        elongation_deg,
        illumination_fraction: pheno.phase_illuminated,
        phase_angle_deg: pheno.phase_angle_deg,
        moon_lon_deg: moon.lon_deg,
        moon_lat_deg: moon.lat_deg,
        moon_dist_au: moon.dist_au,
        moon_dist_km: moon.dist_km,
        apparent_diameter_deg: pheno.apparent_diameter_deg,
        magnitude_v: pheno.apparent_magnitude_v,
        age_days,
        synodic_fraction,
        phase_8_index,
        phase_8_name,
        phase_28_index,
        prev_new_moon_jd_ut1: prev_nm.jd_ut1,
        next_new_moon_jd_ut1: next_nm.jd_ut1,
        next_first_quarter_jd_ut1: next_fq.jd_ut1,
        next_full_moon_jd_ut1: next_fm.jd_ut1,
        next_last_quarter_jd_ut1: next_lq.jd_ut1,
    })
}

// ============================================================================
// C. Planetary Stations (\dot{\lambda} = 0) and Declination Extrema (\dot{\delta} = 0)
// ============================================================================

/// Kind of planetary stationary point ($\dot{\lambda}(t) = 0$).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum StationKind {
    /// Station Retrograde (direct $\to$ retrograde: $\dot{\lambda}$ crosses $0$ from $+$ to $-$).
    StationRetrograde,
    /// Station Direct (retrograde $\to$ direct: $\dot{\lambda}$ crosses $0$ from $-$ to $+$).
    StationDirect,
}

/// Planetary or lunar-point stationary event ($\dot{\lambda}(t) = 0$, `R9.2`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StationEvent {
    /// Body identifier.
    pub body: BodyId,
    /// Station transition kind (`StationRetrograde` or `StationDirect`).
    pub kind: StationKind,
    /// Terrestrial Time (`TT`) Julian Date.
    pub jd_tt: f64,
    /// Universal Time (`UT1`) Julian Date.
    pub jd_ut1: f64,
    /// $\Delta T = \text{TT} - \text{UT1}$ in seconds.
    pub delta_t_sec: f64,
    /// Ecliptic longitude in degrees `[0, 360)` at the station.
    pub lon_deg: f64,
    /// Ecliptic latitude in degrees `[-90, +90]` at the station.
    pub lat_deg: f64,
    /// Radial distance in AU at the station.
    pub dist_au: f64,
    /// Radial distance in kilometers at the station.
    pub dist_km: f64,
}

/// Finds all retrograde and direct stations ($\dot{\lambda}(t) = 0$) of `body` within `[start_input, end_input]`
/// under `options` (`R9.2`). Works seamlessly in both `Zodiac::Tropical` and `Zodiac::Sidereal(model)`.
pub fn search_stations_in_range(
    engine: &Engine,
    body: BodyId,
    start_input: TimeInput,
    end_input: TimeInput,
    options: &CalcOptions,
) -> Result<Vec<StationEvent>, EphemerisError> {
    let mut eval_opts = *options;
    eval_opts.plane = Plane::Ecliptic;
    eval_opts.include_speed = true;

    let obs = match eval_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let (_, tt_start, _, _) = engine.time_table().resolve_scales(start_input, obs);
    let (_, tt_end, _, _) = engine.time_table().resolve_scales(end_input, obs);
    let jd_tt_start = tt_start.jd_f64();
    let jd_tt_end = tt_end.jd_f64();

    if !jd_tt_start.is_finite() || !jd_tt_end.is_finite() || jd_tt_end < jd_tt_start {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid time range in search_stations_in_range: start={jd_tt_start}, end={jd_tt_end}"
            ),
        });
    }

    // Bodies with strictly monotonic ecliptic motion (Sun, Moon, Mean Node, Mean Lilith, Fixed Stars)
    // do not have retrograde stations.
    if matches!(
        body,
        BodyId::Planet(
            PlanetId::Sun
                | PlanetId::Moon
                | PlanetId::Earth
                | PlanetId::EarthMoonBarycenter
                | PlanetId::SolarSystemBarycenter
        ) | BodyId::LunarPoint(
            LunarPointId::MeanNode | LunarPointId::MeanSouthNode | LunarPointId::MeanLilith
        ) | BodyId::EclNut
    ) || (jd_tt_end - jd_tt_start) < 1.0e-9
    {
        return Ok(Vec::new());
    }

    let span = jd_tt_end - jd_tt_start;
    let base_step: f64 = if matches!(eval_opts.origin, Origin::Topocentric(_)) {
        0.1
    } else {
        match body {
            BodyId::FixedStar(_) => 0.25,
            BodyId::LunarPoint(_) => 0.5,
            BodyId::Planet(PlanetId::Mercury) => 2.0,
            BodyId::Planet(PlanetId::Venus) => 3.0,
            _ => 5.0,
        }
    };
    let step_days = base_step.min(span * 0.5).max(1.0e-3);
    let total_steps = ((span / step_days).ceil() as usize).max(1);
    let v_tol = if matches!(body, BodyId::FixedStar(_)) {
        1.0e-12
    } else {
        1.0e-9
    };

    let mut t_prev = jd_tt_start;
    let mut res_prev = engine.calc(body, TimeInput::from_tt_jd(t_prev), eval_opts)?;
    let mut events: Vec<StationEvent> = Vec::new();

    for step in 1..=total_steps {
        let t_curr = if step == total_steps {
            jd_tt_end
        } else {
            (jd_tt_start + (step as f64) * step_days).min(jd_tt_end)
        };
        let res_curr = engine.calc(body, TimeInput::from_tt_jd(t_curr), eval_opts)?;

        let v_prev = res_prev.lon_speed_deg_day;
        let v_curr = res_curr.lon_speed_deg_day;

        let kind_opt = if v_prev > 0.0 && v_curr <= 0.0 {
            Some(StationKind::StationRetrograde)
        } else if v_prev < 0.0 && v_curr >= 0.0 {
            Some(StationKind::StationDirect)
        } else {
            None
        };

        if let Some(kind) = kind_opt {
            // Refine root of lon_speed_deg_day(t) == 0 via the Illinois method
            let mut a = t_prev;
            let mut b = t_curr;
            let mut fa = v_prev;
            let mut fb = v_curr;
            let mut best_res = if fa.abs() <= fb.abs() {
                res_prev.clone()
            } else {
                res_curr.clone()
            };
            let mut side = 0i8;

            for _ in 0..30 {
                if (b - a).abs() < 1.0e-7 || fa.abs().min(fb.abs()) < v_tol {
                    break;
                }
                let denom = fb - fa;
                let c = if denom.abs() > 1.0e-15 {
                    ((a * fb - b * fa) / denom).clamp(a + 0.01 * (b - a), b - 0.01 * (b - a))
                } else {
                    0.5 * (a + b)
                };
                let res_c = engine.calc(body, TimeInput::from_tt_jd(c), eval_opts)?;
                let fc = res_c.lon_speed_deg_day;
                if fc.abs() < fa.abs().min(fb.abs()) {
                    best_res = res_c.clone();
                }
                if fc.abs() < v_tol {
                    best_res = res_c;
                    break;
                }
                if (fc > 0.0) == (fb > 0.0) {
                    b = c;
                    fb = fc;
                    if side == -1 {
                        fa *= 0.5;
                    }
                    side = -1;
                } else {
                    a = c;
                    fa = fc;
                    if side == 1 {
                        fb *= 0.5;
                    }
                    side = 1;
                }
            }

            if events
                .last()
                .map_or(true, |last| (best_res.jd_tt - last.jd_tt).abs() > 1.0e-3)
            {
                events.push(StationEvent {
                    body,
                    kind,
                    jd_tt: best_res.jd_tt,
                    jd_ut1: best_res.jd_ut1,
                    delta_t_sec: best_res.delta_t_sec,
                    lon_deg: best_res.lon_deg,
                    lat_deg: best_res.lat_deg,
                    dist_au: best_res.dist_au,
                    dist_km: best_res.dist_km,
                });
            }
        }

        t_prev = t_curr;
        res_prev = res_curr;
    }

    Ok(events)
}

/// Kind of declination extremum ($\dot{\delta}(t) = 0$).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DeclinationExtremumKind {
    /// Northern declination maximum ($\dot{\delta}$ crosses $0$ from $+$ to $-$, northern lunistice / standstill).
    NorthernMaximum,
    /// Southern declination minimum ($\dot{\delta}$ crosses $0$ from $-$ to $+$, southern lunistice / standstill).
    SouthernMinimum,
}

/// Lunar or planetary declination extremum event ($\dot{\delta}(t) = 0$, `R9.2`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DeclinationExtremumEvent {
    /// Body identifier.
    pub body: BodyId,
    /// Extremum kind (`NorthernMaximum` or `SouthernMinimum`).
    pub kind: DeclinationExtremumKind,
    /// Terrestrial Time (`TT`) Julian Date.
    pub jd_tt: f64,
    /// Universal Time (`UT1`) Julian Date.
    pub jd_ut1: f64,
    /// $\Delta T = \text{TT} - \text{UT1}$ in seconds.
    pub delta_t_sec: f64,
    /// True Right Ascension of Date $\alpha$ in degrees `[0, 360)`.
    pub ra_deg: f64,
    /// True Declination of Date $\delta$ in degrees `[-90, +90]`.
    pub dec_deg: f64,
    /// Radial distance in AU.
    pub dist_au: f64,
    /// Radial distance in kilometers.
    pub dist_km: f64,
}

/// Finds all northern and southern declination extrema ($\dot{\delta}(t) = 0$) of `body`
/// within `[start_input, end_input]` (`R9.2`).
pub fn search_declination_extrema_in_range(
    engine: &Engine,
    body: BodyId,
    start_input: TimeInput,
    end_input: TimeInput,
    options: &CalcOptions,
) -> Result<Vec<DeclinationExtremumEvent>, EphemerisError> {
    if body == BodyId::EclNut {
        return Err(EphemerisError::UnsupportedBody {
            body: "EclNut".to_owned(),
        });
    }

    let mut eval_opts = *options;
    eval_opts.plane = Plane::Equatorial;
    eval_opts.equinox = Equinox::OfDate;
    eval_opts.zodiac = Zodiac::Tropical;
    eval_opts.include_speed = true;

    let obs = match eval_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let (_, tt_start, _, _) = engine.time_table().resolve_scales(start_input, obs);
    let (_, tt_end, _, _) = engine.time_table().resolve_scales(end_input, obs);
    let jd_tt_start = tt_start.jd_f64();
    let jd_tt_end = tt_end.jd_f64();

    if !jd_tt_start.is_finite() || !jd_tt_end.is_finite() || jd_tt_end < jd_tt_start {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid time range in search_declination_extrema_in_range: start={jd_tt_start}, end={jd_tt_end}"
            ),
        });
    }
    if (jd_tt_end - jd_tt_start) < 1.0e-9 {
        return Ok(Vec::new());
    }

    let span = jd_tt_end - jd_tt_start;
    let base_step: f64 = if matches!(eval_opts.origin, Origin::Topocentric(_))
        && matches!(body, BodyId::Planet(PlanetId::Moon) | BodyId::LunarPoint(_))
    {
        0.125
    } else {
        match body {
            BodyId::Planet(PlanetId::Moon) | BodyId::LunarPoint(_) => 1.5,
            BodyId::Planet(PlanetId::Mercury | PlanetId::Venus) => 2.0,
            _ => 5.0,
        }
    };
    let step_days = base_step.min(span * 0.5).max(1.0e-3);
    let total_steps = ((span / step_days).ceil() as usize).max(1);

    let mut t_prev = jd_tt_start;
    let mut res_prev = engine.calc(body, TimeInput::from_tt_jd(t_prev), eval_opts)?;
    let mut events: Vec<DeclinationExtremumEvent> = Vec::new();

    for step in 1..=total_steps {
        let t_curr = if step == total_steps {
            jd_tt_end
        } else {
            (jd_tt_start + (step as f64) * step_days).min(jd_tt_end)
        };
        let res_curr = engine.calc(body, TimeInput::from_tt_jd(t_curr), eval_opts)?;

        let v_prev = res_prev.lat_speed_deg_day;
        let v_curr = res_curr.lat_speed_deg_day;

        let kind_opt = if v_prev > 0.0 && v_curr <= 0.0 {
            Some(DeclinationExtremumKind::NorthernMaximum)
        } else if v_prev < 0.0 && v_curr >= 0.0 {
            Some(DeclinationExtremumKind::SouthernMinimum)
        } else {
            None
        };

        if let Some(kind) = kind_opt {
            let mut a = t_prev;
            let mut b = t_curr;
            let mut fa = v_prev;
            let mut fb = v_curr;
            let mut best_res = if fa.abs() <= fb.abs() {
                res_prev.clone()
            } else {
                res_curr.clone()
            };
            let mut side = 0i8;

            for _ in 0..30 {
                if (b - a).abs() < 1.0e-7 || fa.abs().min(fb.abs()) < 1.0e-7 {
                    break;
                }
                let denom = fb - fa;
                let c = if denom.abs() > 1.0e-15 {
                    ((a * fb - b * fa) / denom).clamp(a + 0.01 * (b - a), b - 0.01 * (b - a))
                } else {
                    0.5 * (a + b)
                };
                let res_c = engine.calc(body, TimeInput::from_tt_jd(c), eval_opts)?;
                let fc = res_c.lat_speed_deg_day;
                if fc.abs() < fa.abs().min(fb.abs()) {
                    best_res = res_c.clone();
                }
                if fc.abs() < 1.0e-7 {
                    best_res = res_c;
                    break;
                }
                if (fc > 0.0) == (fb > 0.0) {
                    b = c;
                    fb = fc;
                    if side == -1 {
                        fa *= 0.5;
                    }
                    side = -1;
                } else {
                    a = c;
                    fa = fc;
                    if side == 1 {
                        fb *= 0.5;
                    }
                    side = 1;
                }
            }

            if events
                .last()
                .map_or(true, |last| (best_res.jd_tt - last.jd_tt).abs() > 1.0e-3)
            {
                events.push(DeclinationExtremumEvent {
                    body,
                    kind,
                    jd_tt: best_res.jd_tt,
                    jd_ut1: best_res.jd_ut1,
                    delta_t_sec: best_res.delta_t_sec,
                    ra_deg: best_res.lon_deg,
                    dec_deg: best_res.lat_deg,
                    dist_au: best_res.dist_au,
                    dist_km: best_res.dist_km,
                });
            }
        }

        t_prev = t_curr;
        res_prev = res_curr;
    }

    Ok(events)
}

// ============================================================================
// D. Brady Diurnal Fixed-Star Parans (R6.2)
// ============================================================================

/// Diurnal cardinal angle for Bernadette Brady fixed-star parans (`R6.2`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ParanAngle {
    /// Geometric center rising across the eastern horizon ($a = 0^\circ, \dot{a} > 0$).
    #[serde(rename = "rising")]
    Rising,
    /// Upper meridian transit / culmination ($H = 0^\circ$).
    #[serde(rename = "culminating")]
    Culminating,
    /// Geometric center setting across the western horizon ($a = 0^\circ, \dot{a} < 0$).
    #[serde(rename = "setting")]
    Setting,
    /// Lower meridian transit / anti-culmination ($H = \pm 180^\circ$).
    #[serde(rename = "anti-culminating")]
    AntiCulminating,
}

impl ParanAngle {
    /// All 4 diurnal paran angles in canonical order.
    pub const ALL: [Self; 4] = [
        Self::Rising,
        Self::Culminating,
        Self::Setting,
        Self::AntiCulminating,
    ];

    /// Canonical string identifier per `2026-06-08-fixed-stars-parans-endpoint-design.md`
    /// (`"rising"`, `"culminating"`, `"setting"`, `"anti-culminating"`).
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Rising => "rising",
            Self::Culminating => "culminating",
            Self::Setting => "setting",
            Self::AntiCulminating => "anti-culminating",
        }
    }

    /// Maps this [`ParanAngle`] to the corresponding [`RiseTransEvent`].
    #[must_use]
    pub const fn to_rise_trans_event(self) -> RiseTransEvent {
        match self {
            Self::Rising => RiseTransEvent::Rise,
            Self::Culminating => RiseTransEvent::UpperTransit,
            Self::Setting => RiseTransEvent::Set,
            Self::AntiCulminating => RiseTransEvent::LowerTransit,
        }
    }
}

/// A single Bernadette Brady diurnal fixed-star paran aspect (`R6.2`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ParanEvent {
    /// Canonical star name (e.g. `"Regulus"`, `"Spica"`, `"Sirius"`).
    pub star_name: String,
    /// Bayer/Flamsteed identifier (e.g. `"alLeo"`, `"alVir"`, `"alCMa"`).
    pub star_bayer: String,
    /// Hipparcos catalog identifier (`HIP`).
    pub star_hip: u32,
    /// Planet / body identifier.
    pub planet: BodyId,
    /// Canonical planet name (e.g. `"Sun"`, `"Moon"`, `"Mars"`).
    pub planet_name: String,
    /// Primary paran classification angle (always equal to `star_angle` per FE contract).
    pub paran_type: ParanAngle,
    /// Diurnal angle of the fixed star (`Rising`, `Culminating`, `Setting`, `AntiCulminating`).
    pub star_angle: ParanAngle,
    /// Simultaneous diurnal angle of the planet (`Rising`, `Culminating`, `Setting`, `AntiCulminating`).
    pub planet_angle: ParanAngle,
    /// `true` if `star_angle == planet_angle` (co-rising, co-culminating, co-setting, or co-anti-culminating).
    pub same_angle: bool,
    /// Exact `JD(UT1)` when the star crosses `star_angle`.
    pub star_event_jd_ut1: f64,
    /// Exact `JD(UT1)` when the planet crosses `planet_angle`.
    pub planet_event_jd_ut1: f64,
    /// Circular diurnal time difference in minutes $\in [0, 720]$.
    pub time_orb_minutes: f64,
}

/// Configuration options for [`find_parans`] (`R6.2`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ParanOptions {
    /// Maximum circular diurnal time orb in minutes (default `4.0` minutes $\approx 1^\circ$ of Earth rotation).
    pub orb_minutes: f64,
    /// Natal bodies to include in the paran scan (default: the 10 major bodies `Sun..=Pluto`).
    pub bodies: Vec<BodyId>,
    /// Optional filter of star names or Bayer identifiers (`None` checks all catalog stars).
    pub star_names: Option<Vec<String>>,
}

impl Default for ParanOptions {
    fn default() -> Self {
        Self {
            orb_minutes: 4.0,
            bodies: vec![
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
            ],
            star_names: None,
        }
    }
}

/// Solves all 4 diurnal geometric-center crossings `(Rising, Culminating, Setting, AntiCulminating)`
/// in `[jd_ut1_start, jd_ut1_start + 1.05]` for a fixed star using 2-epoch (`t0`, `t0 + 1d`)
/// full 3D Hipparcos-2 apparent equatorial coordinates (`pmpx + ldsun + ab + r_bpn`), exact
/// instantaneous diurnal aberration (`pvtob`), and `GAST` Newton root-polishing (`< 0.005 s`
/// agreement with `rise_trans` across all declinations including `Polaris`).
fn solve_star_diurnal_events(
    star: &StarEntry,
    observer: Observer,
    jd_ut1_start: f64,
    m0: &Moment<'_>,
    m1: &Moment<'_>,
) -> [Option<(ParanAngle, [f64; 3])>; 4] {
    // Evaluate geocentric apparent True Equatorial of Date unit vectors at m0 and m1.
    // Over 24 hours, geocentric star position (proper motion + annual parallax + solar deflection
    // + annual aberration + precession + nutation) varies smoothly and linearly by < 0.15".
    // Topocentric diurnal aberration rotates 360° per sidereal day with GAST and is applied
    // analytically at each candidate (u, gast) below.
    let geo_eq_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };
    let u0_geo = eval_star_unit_vector(star, m0.center(), &geo_eq_opts);
    let u1_geo = eval_star_unit_vector(star, m1.center(), &geo_eq_opts);

    let phi = observer.latitude_deg.to_radians();
    let lon_rad = observer.longitude_deg.to_radians();
    let (sin_phi, cos_phi) = phi.sin_cos();

    let gast0 = m0.center().gast_rad;
    let gast1 = m1.center().gast_rad;
    // Mean sidereal rate over the 24h day (~2*pi * 1.00273790935 rad/day)
    let omega_sid = wrap_pm_pi(gast1 - gast0) + std::f64::consts::TAU;
    let sidereal_day = std::f64::consts::TAU / omega_sid;

    let dt0_days = m0.center().delta_t_info.delta_t_sec / SECONDS_PER_DAY;
    let ddt_days =
        (m1.center().delta_t_info.delta_t_sec - m0.center().delta_t_info.delta_t_sec) / SECONDS_PER_DAY;

    // Evaluates topocentric apparent (ra_rad, dec_rad) at fractional day `u_days` and sidereal angle `gast`.
    let topo_ra_dec_at = |u_days: f64, gast: f64| -> (f64, f64) {
        let u_interp = [
            u0_geo[0] + u_days * (u1_geo[0] - u0_geo[0]),
            u0_geo[1] + u_days * (u1_geo[1] - u0_geo[1]),
            u0_geo[2] + u_days * (u1_geo[2] - u0_geo[2]),
        ];
        let mut pv_m = [[0.0_f64; 3]; 2];
        pvtob(
            lon_rad,
            phi,
            observer.altitude_m,
            0.0,
            0.0,
            0.0,
            gast,
            &mut pv_m,
        );
        let v_diurnal_over_c = [
            pv_m[1][0] / (1000.0 * C_KM_S),
            pv_m[1][1] / (1000.0 * C_KM_S),
            pv_m[1][2] / (1000.0 * C_KM_S),
        ];
        let u_topo = vec_unit(&vec_add(&u_interp, &v_diurnal_over_c));
        let sph = cart_to_spherical(&u_topo);
        (sph.lon_deg.to_radians(), sph.lat_deg.to_radians())
    };

    let (ra0_rad, dec0_rad) = topo_ra_dec_at(0.0, gast0);

    let target_ha_for_dec = |angle: ParanAngle, dec_rad: f64| -> Option<f64> {
        match angle {
            ParanAngle::Culminating => Some(0.0),
            ParanAngle::AntiCulminating => Some(std::f64::consts::PI),
            ParanAngle::Rising | ParanAngle::Setting => {
                let (sin_dec, cos_dec) = dec_rad.sin_cos();
                let denom = cos_phi * cos_dec;
                if denom.abs() <= 1.0e-12 {
                    return None;
                }
                let cos_h0 = (-sin_phi * sin_dec) / denom;
                if !(-1.0..1.0).contains(&cos_h0) && cos_h0.abs() >= 1.0 {
                    return None;
                }
                let h0 = cos_h0.clamp(-1.0, 1.0).acos();
                Some(if angle == ParanAngle::Rising { -h0 } else { h0 })
            }
        }
    };

    let solve_angle = |angle: ParanAngle| -> Option<(ParanAngle, [f64; 3])> {
        let h_targ_0 = target_ha_for_dec(angle, dec0_rad)?;
        let h_start = gast0 + lon_rad - ra0_rad;
        let dh0 = (h_targ_0 - h_start).rem_euclid(std::f64::consts::TAU);
        let u0 = dh0 / omega_sid;

        let refine_u = |mut u: f64| -> Option<f64> {
            for _ in 0..3 {
                let t_ut1 = jd_ut1_start + u;
                let t_tt = t_ut1 + dt0_days + u * ddt_days;
                let gast = gast_rad(Ut1::from_jd(t_ut1), Tt::from_jd(t_tt));
                let (ra, dec) = topo_ra_dec_at(u, gast);
                let h_curr = gast + lon_rad - ra;
                let h_targ = target_ha_for_dec(angle, dec)?;
                let err = wrap_pm_pi(h_targ - h_curr);
                u += err / omega_sid;
            }
            Some(jd_ut1_start + u)
        };

        let t_mid = refine_u(u0)?;
        let t_prev = refine_u(u0 - sidereal_day).unwrap_or(t_mid - sidereal_day);
        let t_next = refine_u(u0 + sidereal_day).unwrap_or(t_mid + sidereal_day);
        Some((angle, [t_prev, t_mid, t_next]))
    };

    [
        solve_angle(ParanAngle::Rising),
        solve_angle(ParanAngle::Culminating),
        solve_angle(ParanAngle::Setting),
        solve_angle(ParanAngle::AntiCulminating),
    ]
}

type BodyDiurnalRow = (BodyId, String, [Option<(ParanAngle, f64)>; 4]);

/// Solves all 4 diurnal geometric-center crossings `(Rising, Culminating, Setting, AntiCulminating)`
/// in `[jd_ut1_start, jd_ut1_start + 1.05]` for each body in `bodies`, sharing the precomputed
/// 30-minute coarse [`Moment`] grid across all bodies and refining each bracket with the exact
/// Illinois method of [`crate::astrometry::rise_trans`].
fn solve_bodies_diurnal_events(
    engine: &Engine,
    bodies: &[BodyId],
    observer: Observer,
    jd_ut1_start: f64,
) -> Result<Vec<BodyDiurnalRow>, EphemerisError> {
    let max_search_days = 1.05_f64;
    let step_days = 1.0 / 48.0;
    let total_steps = ((max_search_days / step_days).ceil() as usize).max(2);

    // Precompute the coarse Moment grid once across the 25.2-hour window
    let mut coarse_grid: Vec<(f64, Moment<'_>)> = Vec::with_capacity(total_steps + 1);
    for step in 0..=total_steps {
        let t = (jd_ut1_start + (step as f64) * step_days).min(jd_ut1_start + max_search_days);
        let m = engine.at(TimeInput::from_ut1_jd(t), Some(observer))?;
        coarse_grid.push((t, m));
    }

    let calc_opts = CalcOptions {
        origin: Origin::Topocentric(observer),
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };

    let eval_alt_ha_in_moment =
        |m: &Moment<'_>, body: BodyId| -> Result<(f64, f64), EphemerisError> {
            let res = calc_body_in_moment(m, body, &calc_opts)?;
            let hz = azalt(
                observer,
                m.center().gast_rad,
                res.lon_deg,
                res.lat_deg,
                AtmosphereParams::airless(),
            );
            Ok((hz.altitude_true_deg, hz.hour_angle_deg))
        };

    let eval_alt_ha_at_ut1 = |jd_ut1: f64, body: BodyId| -> Result<(f64, f64), EphemerisError> {
        let m = engine.at(TimeInput::from_ut1_jd(jd_ut1), Some(observer))?;
        eval_alt_ha_in_moment(&m, body)
    };

    let mut out = Vec::with_capacity(bodies.len());
    let Some(((t_first, m_first), rest_grid)) = coarse_grid.split_first() else {
        return Ok(out);
    };

    for &body in bodies {
        let mut brackets: [Option<(f64, f64, f64, f64)>; 4] = [None; 4];
        let (mut alt_prev, mut ha_prev) = eval_alt_ha_in_moment(m_first, body)?;
        let mut t_prev = *t_first;

        for (t_curr_ref, m_curr) in rest_grid {
            let t_curr = *t_curr_ref;
            let (alt_curr, ha_curr) = eval_alt_ha_in_moment(m_curr, body)?;

            if brackets[0].is_none() && alt_prev < 0.0 && alt_curr >= 0.0 {
                brackets[0] = Some((t_prev, alt_prev, t_curr, alt_curr));
            }
            if brackets[1].is_none()
                && ha_prev < 0.0
                && ha_curr >= 0.0
                && (ha_curr - ha_prev).abs() < 90.0
            {
                brackets[1] = Some((t_prev, ha_prev, t_curr, ha_curr));
            }
            if brackets[2].is_none() && alt_prev > 0.0 && alt_curr <= 0.0 {
                brackets[2] = Some((t_prev, alt_prev, t_curr, alt_curr));
            }
            if brackets[3].is_none() {
                let l_prev = (ha_prev + 360.0).rem_euclid(360.0) - 180.0;
                let l_curr = (ha_curr + 360.0).rem_euclid(360.0) - 180.0;
                if l_prev < 0.0 && l_curr >= 0.0 && (l_curr - l_prev).abs() < 90.0 {
                    brackets[3] = Some((t_prev, l_prev, t_curr, l_curr));
                }
            }

            if brackets.iter().all(Option::is_some) {
                break;
            }

            t_prev = t_curr;
            alt_prev = alt_curr;
            ha_prev = ha_curr;
        }

        let mut solved_angles: [Option<(ParanAngle, f64)>; 4] = [None; 4];
        for ((angle, bracket_opt), slot) in ParanAngle::ALL
            .into_iter()
            .zip(brackets)
            .zip(solved_angles.iter_mut())
        {
            if let Some((t_a, f_a, t_b, f_b)) = bracket_opt {
                let f_scalar = |t: f64| -> Result<f64, EphemerisError> {
                    let (alt, ha) = eval_alt_ha_at_ut1(t, body)?;
                    Ok(match angle {
                        ParanAngle::Rising | ParanAngle::Setting => alt,
                        ParanAngle::Culminating => ha,
                        ParanAngle::AntiCulminating => (ha + 360.0).rem_euclid(360.0) - 180.0,
                    })
                };

                let mut a = t_a;
                let mut b = t_b;
                let mut fa = f_a;
                let mut fb = f_b;
                let mut side = 0i8;
                let mut root = if fa.abs() <= fb.abs() { a } else { b };

                for _ in 0..20 {
                    if (b - a).abs() < 1.0e-7 || fa.abs().min(fb.abs()) < 1.0e-6 {
                        root = if fa.abs() <= fb.abs() { a } else { b };
                        break;
                    }
                    let denom = fb - fa;
                    let c = if denom.abs() > 1.0e-15 {
                        (a * fb - b * fa) / denom
                    } else {
                        0.5 * (a + b)
                    };
                    let fc = f_scalar(c)?;
                    if fc.abs() < 1.0e-6 {
                        root = c;
                        break;
                    }
                    if (fc > 0.0) == (fb > 0.0) {
                        b = c;
                        fb = fc;
                        if side == -1 {
                            fa *= 0.5;
                        }
                        side = -1;
                    } else {
                        a = c;
                        fa = fc;
                        if side == 1 {
                            fb *= 0.5;
                        }
                        side = 1;
                    }
                    root = if fa.abs() <= fb.abs() { a } else { b };
                }
                *slot = Some((angle, root));
            }
        }

        out.push((body, body_display_name(body), solved_angles));
    }

    Ok(out)
}

/// Computes Bernadette Brady diurnal fixed-star parans (`R6.2`) for `observer` starting at `start_input`:
/// - Uses airless (`pressure_hpa = 0.0`) geometric center (`DiscMode::Center`) horizon crossings.
/// - Propagates all fixed stars with full 3D Hipparcos-2 proper motion, parallax, radial velocity,
///   gravitational deflection, and aberration to the epoch of calculation.
/// - Handles circumpolar stars and planets gracefully (emitting `Culminating` and `AntiCulminating`
///   without failing when `Rising`/`Setting` do not occur).
/// - Sorts results ascending by `time_orb_minutes` (with deterministic tie-breaking by `star_name`, `planet_name`).
pub fn find_parans(
    engine: &Engine,
    start_input: TimeInput,
    observer: Observer,
    options: &ParanOptions,
) -> Result<Vec<ParanEvent>, EphemerisError> {
    let _ = Observer::new(
        observer.latitude_deg,
        observer.longitude_deg,
        observer.altitude_m,
    )?;
    if !options.orb_minutes.is_finite() || !(0.0..=720.0).contains(&options.orb_minutes) {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid ParanOptions::orb_minutes: {} (must be in [0.0, 720.0])",
                options.orb_minutes
            ),
        });
    }

    let stars: Vec<StarEntry> = match &options.star_names {
        Some(names) => names
            .iter()
            .map(|name| engine.resolve_star_entry(name))
            .collect::<Result<Vec<_>, _>>()?,
        None => engine.all_star_entries(),
    };

    if stars.is_empty() || options.bodies.is_empty() {
        return Ok(Vec::new());
    }

    let (ut1_start, _, _, _) = engine
        .time_table()
        .resolve_scales(start_input, Some(observer));
    let jd_ut1_start = ut1_start.jd_f64();

    // 1. Solve diurnal crossings for all requested planets/bodies
    let body_events = solve_bodies_diurnal_events(engine, &options.bodies, observer, jd_ut1_start)?;

    // 2. Solve diurnal crossings for all requested fixed stars
    let m0 = engine.at(TimeInput::from_ut1_jd(jd_ut1_start), Some(observer))?;
    let m1 = engine.at(TimeInput::from_ut1_jd(jd_ut1_start + 1.0), Some(observer))?;

    let mut parans = Vec::new();

    for star in &stars {
        let star_angles = solve_star_diurnal_events(star, observer, jd_ut1_start, &m0, &m1);

        for (star_angle, t_star_candidates) in star_angles.into_iter().flatten() {
            for (planet_id, planet_name, planet_angles) in &body_events {
                for (planet_angle, t_planet) in planet_angles.iter().copied().flatten() {
                    let mut best_t_star = t_star_candidates[1];
                    let mut best_diff_days = (best_t_star - t_planet).abs();
                    for &cand in &t_star_candidates {
                        let d = (cand - t_planet).abs();
                        if d < best_diff_days {
                            best_diff_days = d;
                            best_t_star = cand;
                        }
                    }
                    let circ_min = best_diff_days * 1440.0;

                    if circ_min <= options.orb_minutes {
                        parans.push(ParanEvent {
                            star_name: star.name.clone(),
                            star_bayer: star.bayer.clone(),
                            star_hip: star.hip,
                            planet: *planet_id,
                            planet_name: planet_name.clone(),
                            paran_type: star_angle,
                            star_angle,
                            planet_angle,
                            same_angle: star_angle == planet_angle,
                            star_event_jd_ut1: best_t_star,
                            planet_event_jd_ut1: t_planet,
                            time_orb_minutes: circ_min,
                        });
                    }
                }
            }
        }
    }

    parans.sort_by(|a, b| {
        a.time_orb_minutes
            .total_cmp(&b.time_orb_minutes)
            .then_with(|| a.star_name.cmp(&b.star_name))
            .then_with(|| a.planet_name.cmp(&b.planet_name))
            .then_with(|| a.star_angle.as_str().cmp(b.star_angle.as_str()))
            .then_with(|| a.planet_angle.as_str().cmp(b.planet_angle.as_str()))
    });

    Ok(parans)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astrometry::RiseTransOptions;
    use crate::engine::EngineConfig;
    use crate::stars::find_star_by_name_or_bayer;

    #[test]
    fn test_fast_star_diurnal_solver_matches_rise_trans_to_sub_tenth_second() {
        let engine = Engine::new_moshier(EngineConfig::default());
        let london = Observer::new(51.5074, -0.1278, 0.0).expect("London observer");
        let jd_start = 2_440_359.083_333_333_5; // 1969-05-17 14:00 UT
        let m0 = engine
            .at(TimeInput::from_ut1_jd(jd_start), Some(london))
            .expect("m0");
        let m1 = engine
            .at(TimeInput::from_ut1_jd(jd_start + 1.0), Some(london))
            .expect("m1");
        let mut rt_opts = RiseTransOptions::geometric_center();
        rt_opts.max_search_days = 1.05;

        for name in ["Regulus", "Sirius", "Spica", "Aldebaran", "Dubhe", "Polaris"] {
            let star = find_star_by_name_or_bayer(name).expect("star in catalog");
            let fast = solve_star_diurnal_events(star, london, jd_start, &m0, &m1);

            for (angle, fast_opt) in ParanAngle::ALL.into_iter().zip(fast) {
                let exact = engine
                    .rise_trans(
                        BodyId::FixedStar(star.hip_id()),
                        london,
                        jd_start,
                        angle.to_rise_trans_event(),
                        rt_opts,
                    )
                    .expect("rise_trans");
                assert_eq!(
                    fast_opt.is_some(),
                    exact.is_some(),
                    "Star {name} {angle:?} existence mismatch: fast={fast_opt:?}, exact={exact:?}"
                );
                if let (Some((_, t_fast)), Some(t_exact)) = (fast_opt, exact) {
                    let diff_sec = (t_fast[1] - t_exact).abs() * SECONDS_PER_DAY;
                    assert!(
                        diff_sec < 0.05,
                        "Star {name} {angle:?} mismatch: fast={:?}, exact={t_exact}, diff={diff_sec:.4} s",
                        t_fast
                    );
                }
            }
        }
    }

    #[test]
    fn test_solar_term_and_lunar_phase_summary_deserialize_owned() {
        let engine = Engine::new_moshier(EngineConfig::default());
        let start = TimeInput::from_ut1_jd(2_461_041.5); // 2026-01-01
        let end = TimeInput::from_ut1_jd(2_461_406.5);   // 2027-01-01

        let terms = solar_terms_in_range(&engine, start, end).expect("solar_terms_in_range");
        assert_eq!(terms.len(), 24);

        let xiaohan = terms.first().expect("Xiaohan");
        assert_eq!(xiaohan.term_index, 19);
        assert_eq!(xiaohan.pinyin_name, "Xiaohan");
        assert_eq!(xiaohan.chinese_name, "小寒");
        assert_eq!(xiaohan.english_name, "Minor Cold");

        let dahan = terms.get(1).expect("Dahan");
        assert_eq!(dahan.term_index, 20);
        assert_eq!(dahan.pinyin_name, "Dahan");
        assert_eq!(dahan.chinese_name, "大寒");
        assert_eq!(dahan.english_name, "Major Cold");

        let lichun = terms.get(2).expect("Lichun");
        assert_eq!(lichun.term_index, 21);
        assert_eq!(lichun.pinyin_name, "Lichun");
        assert_eq!(lichun.chinese_name, "立春");
        assert_eq!(lichun.english_name, "Start of Spring");

        let chunfen = terms.iter().find(|t| t.term_index == 0).expect("Chunfen");
        assert_eq!(chunfen.english_name, "Vernal Equinox");
        assert_eq!(chunfen.chinese_name, "春分");

        let dongzhi = terms.iter().find(|t| t.term_index == 18).expect("Dongzhi");
        assert_eq!(dongzhi.english_name, "Winter Solstice");
        assert_eq!(dongzhi.chinese_name, "冬至");

        // Verify DeserializeOwned round-trip for SolarTermEvent
        let json_terms = serde_json::to_string(&terms).expect("serialize SolarTermEvent");
        let terms_back: Vec<SolarTermEvent> =
            serde_json::from_str(&json_terms).expect("deserialize SolarTermEvent");
        assert_eq!(terms_back.len(), terms.len());
        for (a, b) in terms_back.iter().zip(terms.iter()) {
            assert_eq!(a.term_index, b.term_index);
            assert_eq!(a.pinyin_name, b.pinyin_name);
            assert_eq!(a.chinese_name, b.chinese_name);
            assert_eq!(a.english_name, b.english_name);
            assert_eq!(a.is_major_term, b.is_major_term);
            assert!((a.jd_tt - b.jd_tt).abs() < 1e-8);
            assert!((a.jd_ut1 - b.jd_ut1).abs() < 1e-8);
        }

        // Verify DeserializeOwned round-trip for LunarPhaseSummary
        let summary = lunar_phase_summary(&engine, start).expect("lunar_phase_summary");
        assert!(!summary.phase_8_name.is_empty());
        let json_summary = serde_json::to_string(&summary).expect("serialize LunarPhaseSummary");
        let summary_back: LunarPhaseSummary =
            serde_json::from_str(&json_summary).expect("deserialize LunarPhaseSummary");
        assert_eq!(summary_back.phase_8_index, summary.phase_8_index);
        assert_eq!(summary_back.phase_8_name, summary.phase_8_name);
        assert_eq!(summary_back.phase_28_index, summary.phase_28_index);
        assert!((summary_back.jd_ut1 - summary.jd_ut1).abs() < 1e-8);
    }

    #[test]
    fn test_planetary_return_and_returns_in_range_with_retrograde_passes() {
        use crate::types::AyanamsaModel;

        let engine = Engine::new_moshier(EngineConfig::default());
        let natal = TimeInput::from_ut1_jd(2_447_892.5); // 1990-01-01 00:00 UT
        let search_after = TimeInput::from_ut1_jd(2_461_000.0); // late 2025

        // 1. Tropical Solar Return in Jan 2026
        let sr_trop = planetary_return(
            &engine,
            BodyId::Planet(PlanetId::Sun),
            natal,
            search_after,
            CalcOptions::default(),
        )
        .expect("Tropical Solar Return");
        assert_eq!(sr_trop.pass_index, 1);
        assert!(!sr_trop.is_retrograde);
        assert!(wrap_pm180(sr_trop.actual_lon_deg - sr_trop.natal_lon_deg).abs() < 1.0e-7);

        // 2. Sidereal (Lahiri) Precession-Corrected Solar Return in Jan 2026
        let sid_opts = CalcOptions {
            zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
            ..Default::default()
        };
        let sr_sid = planetary_return(
            &engine,
            BodyId::Planet(PlanetId::Sun),
            natal,
            search_after,
            sid_opts,
        )
        .expect("Sidereal Lahiri Solar Return");
        assert_eq!(sr_sid.pass_index, 1);
        assert!(wrap_pm180(sr_sid.actual_lon_deg - sr_sid.natal_lon_deg).abs() < 1.0e-7);
        // Over 36 years (1990 -> 2026), ~0.5° of general precession shifts the sidereal solar return
        // by ~0.5 days later than the tropical solar return!
        let precession_shift_days = sr_sid.jd_ut1 - sr_trop.jd_ut1;
        assert!(
            (0.45..=0.58).contains(&precession_shift_days),
            "Expected ~0.51d sidereal solar return shift after 36 yr, got {precession_shift_days:.4}d"
        );

        // 3. Multi-year Solar Returns in range [2025-12-01, 2028-02-01] -> 3 returns, all pass_index == 1
        let returns_3yr = planetary_returns_in_range(
            &engine,
            BodyId::Planet(PlanetId::Sun),
            natal,
            TimeInput::from_ut1_jd(2_461_010.5),
            TimeInput::from_ut1_jd(2_461_802.5),
            CalcOptions::default(),
        )
        .expect("3-year solar returns");
        assert_eq!(returns_3yr.len(), 3);
        assert!(returns_3yr.iter().all(|r| r.pass_index == 1));

        // 4. Retrograde 3-pass Saturn Return series
        // Find Saturn's 2026 retrograde loop midpoint and set natal epoch ~29.45 years earlier
        // or at the midpoint epoch itself, then scan across 2026-2027.
        let start_2026 = TimeInput::from_ut1_jd(2_461_041.5);
        let end_2027 = TimeInput::from_ut1_jd(2_461_467.5);
        let sat_stations = search_stations_in_range(
            &engine,
            BodyId::Planet(PlanetId::Saturn),
            start_2026,
            end_2027,
            &CalcOptions::default(),
        )
        .expect("Saturn stations");
        let st0 = sat_stations.first().expect("st0");
        let st1 = sat_stations.get(1).expect("st1");
        let mid_sat_jd = 0.5 * (st0.jd_tt + st1.jd_tt);
        let sat_returns = planetary_returns_in_range(
            &engine,
            BodyId::Planet(PlanetId::Saturn),
            TimeInput::from_tt_jd(mid_sat_jd),
            start_2026,
            end_2027,
            CalcOptions::default(),
        )
        .expect("Saturn 3-pass return");
        assert_eq!(sat_returns.len(), 3);
        let r0 = sat_returns.first().expect("r0");
        let r1 = sat_returns.get(1).expect("r1");
        let r2 = sat_returns.get(2).expect("r2");
        assert_eq!(r0.pass_index, 1);
        assert!(!r0.is_retrograde);
        assert_eq!(r1.pass_index, 2);
        assert!(r1.is_retrograde);
        assert_eq!(r2.pass_index, 3);
        assert!(!r2.is_retrograde);
    }

    #[test]
    fn test_grazing_retrograde_station_crossings() {
        let engine = Engine::new_moshier(EngineConfig::default());
        let opts = CalcOptions::default();
        let start = TimeInput::from_ut1_jd(2_451_545.0);
        let end = TimeInput::from_ut1_jd(2_451_665.0);
        let stations = search_stations_in_range(
            &engine,
            BodyId::Planet(PlanetId::Mercury),
            start,
            end,
            &opts,
        )
        .expect("Mercury stations");
        let st = stations.first().expect("first Mercury station");
        // Target 1e-4° inside the retrograde station turnaround extremum must produce all 3 passes
        // (direct before station, retrograde immediately after station, and direct after direct station)
        let target = (st.lon_deg - 1.0e-4).rem_euclid(360.0);
        let hits = body_longitude_crossings_in_range(
            &engine,
            BodyId::Planet(PlanetId::Mercury),
            start,
            end,
            target,
            &opts,
        )
        .expect("Mercury grazing crossings");
        assert_eq!(
            hits.len(),
            3,
            "Expected 3 crossings for grazing retrograde station target={target}, got {}",
            hits.len()
        );
        let first_hit = search_body_longitude(
            &engine,
            BodyId::Planet(PlanetId::Mercury),
            start,
            target,
            true,
            120.0,
            &opts,
        )
        .expect("first crossing via search_body_longitude");
        let h0 = hits.first().expect("h0");
        assert!(
            (first_hit.jd_tt - h0.jd_tt).abs() < 1.0e-6,
            "search_body_longitude missed first grazing crossing: got {}, expected {}",
            first_hit.jd_tt,
            h0.jd_tt
        );
    }
}
