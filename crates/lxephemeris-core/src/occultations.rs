//! Global and local lunar occultations of planets, minor bodies, and fixed stars
//! (`swe_lun_occult_when_glob` / `swe_lun_occult_when_loc` compatibility; `ARCHITECTURE.md` §2.6).
//!
//! - **Global Lunar Occultations (`next_lunar_occultation`)**:
//!   Generalizes the Besselian Fundamental Plane solver where the occulted background
//!   source is an arbitrary target (`Planet`, `MinorPlanet`, or `FixedStar`).
//!   For fixed stars ($d \to \infty$), the shadow cylinder has parallel rays ($\tan f_1 = \tan f_2 = 0$)
//!   and cylinder radius equal to the Moon's radius $k_{\text{Moon}} = R_{\text{Moon}} / R_\oplus = 0.27240088$.
//!   An occultation occurs anywhere on Earth whenever the fundamental plane minimum axis distance
//!   $|\gamma| \le 1.0 + L_1$.
//! - **Local Topocentric Lunar Occultations (`next_local_lunar_occultation`)**:
//!   Evaluates topocentric apparent directions $\hat{\mathbf{u}}_{\text{Moon}}(t)$ and $\hat{\mathbf{u}}_{\text{target}}(t)$
//!   and angular separation $\theta(t)$ from a WGS84 observer, detecting contact times $C_1, C_2, C_3, C_4$
//!   when $\theta(t) \le s_{\text{Moon}}(t) + s_{\text{target}}(t)$ and the Moon is above the local horizon
//!   ($h_{\text{Moon}} \ge -0.5^\circ$).

use crate::astrometry::{
    eval_spk_target_vector_au, resolve_spk_target, AtmosphereParams,
};
use crate::error::EphemerisError;
use crate::frame::{cart_to_spherical, vec_dot, vec_norm, vec_scale, vec_sub, vec_unit, Vec3};
use crate::moment::EpochContext;
use crate::spk::SpkSource;
use crate::stars::{eval_star_unit_vector, find_star_by_hip};
use crate::time::{TimeInput, TimeTable};
use crate::types::{
    AbMode, BodyId, CalcOptions, Equinox, NutationMode, Observer, Origin, Plane, PlanetId, Zodiac,
    AU_KM,
};
use serde::{Deserialize, Serialize};

/// WGS84 Earth equatorial radius in km.
pub const R_EARTH_EQ_KM: f64 = 6_378.137;

/// Mean lunar volumetric radius in km (IAU).
pub const R_MOON_KM: f64 = 1_737.40;

/// Lunar radius in units of Earth equatorial radii (`R_MOON_KM / R_EARTH_EQ_KM`).
pub const K_MOON: f64 = R_MOON_KM / R_EARTH_EQ_KM;

/// WGS84 Earth flattening $f = 1 / 298.257223563$.
pub const WGS84_FLATTENING: f64 = 1.0 / 298.257_223_563;

/// WGS84 first eccentricity squared $e^2 = 2f - f^2$.
pub const WGS84_E2: f64 = 2.0 * WGS84_FLATTENING - WGS84_FLATTENING * WGS84_FLATTENING;

// === Legacy Swiss Ephemeris return flag constants ===
/// Central occultation (shadow axis intersects Earth or observer is in umbra).
pub const SE_ECL_CENTRAL: i32 = 1;
/// Non-central occultation (shadow axis misses Earth center).
pub const SE_ECL_NONCENTRAL: i32 = 2;
/// Total occultation phase.
pub const SE_ECL_TOTAL: i32 = 4;
/// Annular occultation phase.
pub const SE_ECL_ANNULAR: i32 = 8;
/// Partial or grazing occultation phase.
pub const SE_ECL_PARTIAL: i32 = 16;
/// Occultation or eclipse is visible above local horizon ($h \ge -0.5^\circ$).
pub const SE_ECL_VISIBLE: i32 = 128;
/// Alias for `SE_ECL_VISIBLE` in lunar occultations.
pub const SE_ECL_OCCULTATION: i32 = 128;

/// Classification of a lunar occultation.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OccultationKind {
    /// Target disk is completely obscured by the Moon.
    Total,
    /// Target disk is larger than the Moon (e.g. Sun during annular eclipse).
    Annular,
    /// Moon's limb covers part of the target disk.
    Partial,
    /// Moon's limb grazes the target within a fraction of a percent of disk radius.
    Grazing,
}

/// Circumstances and contact times of a lunar occultation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LunarOccultationInfo {
    /// Geometric classification of the occultation.
    pub kind: OccultationKind,
    /// Occulted background body.
    pub target: BodyId,
    /// Instant of maximum occultation (`JD_UT1`).
    pub jd_ut1_max: f64,
    /// Instant of maximum occultation (`JD_TT`).
    pub jd_tt_max: f64,
    /// First contact / immersion begin (`JD_UT1`, `None` if grazing/miss).
    pub t1_first_contact_jd_ut1: Option<f64>,
    /// Second contact / total immersion begin (`JD_UT1`, `None` if partial/grazing).
    pub t2_inner_immersion_jd_ut1: Option<f64>,
    /// Third contact / total emersion begin (`JD_UT1`, `None` if partial/grazing).
    pub t3_inner_emersion_jd_ut1: Option<f64>,
    /// Fourth contact / emersion end (`JD_UT1`, `None` if grazing/miss).
    pub t4_last_contact_jd_ut1: Option<f64>,
    /// Minimum shadow-axis distance in Earth equatorial radii (global) or apparent angular separation in degrees (local).
    pub gamma: f64,
    /// Fraction of target diameter covered by Moon ($\ge 1.0$ for total occultation of a point star).
    pub magnitude: f64,
    /// Apparent topocentric altitude of the Moon in degrees at maximum occultation (`None` for global).
    pub moon_altitude_deg: Option<f64>,
    /// Instant of meridian noon (`x(t) = 0`) in `JD_UT1` for global occultations.
    #[serde(default)]
    pub t_noon_jd_ut1: Option<f64>,
    /// Instant when central shadow axis enters Earth (`d = R_earth`) in `JD_UT1` for global occultations.
    #[serde(default)]
    pub t_central_begin_jd_ut1: Option<f64>,
    /// Instant when central shadow axis leaves Earth (`d = R_earth`) in `JD_UT1` for global occultations.
    #[serde(default)]
    pub t_central_end_jd_ut1: Option<f64>,
}

impl LunarOccultationInfo {
    /// Formats local `tret[0..10]` with `[jd_ut1_max, t1, t2, t3, t4, 0, 0, 0, 0, 0]` (`swe_lun_occult_when_loc` compatibility).
    #[must_use]
    pub fn to_swe_tret10(&self) -> [f64; 10] {
        let t1 = self.t1_first_contact_jd_ut1.unwrap_or(0.0);
        let t2 = self.t2_inner_immersion_jd_ut1.unwrap_or(0.0);
        let t3 = self.t3_inner_emersion_jd_ut1.unwrap_or(0.0);
        let t4 = self.t4_last_contact_jd_ut1.unwrap_or(0.0);
        [
            self.jd_ut1_max,
            t1,
            t2,
            t3,
            t4,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    }

    /// Formats global `tret[0..10]` matching `swe_lun_occult_when_glob`:
    /// `[jd_ut1_max, t_noon, t1_begin, t4_end, t2_tot_begin, t3_tot_end, t_cent_begin, t_cent_end, 0, 0]`.
    #[must_use]
    pub fn to_swe_glob_tret10(&self) -> [f64; 10] {
        let t1 = self.t1_first_contact_jd_ut1.unwrap_or(0.0);
        let t4 = self.t4_last_contact_jd_ut1.unwrap_or(0.0);
        let is_central = self.gamma.abs() <= 1.0;
        let t_noon = self
            .t_noon_jd_ut1
            .unwrap_or(if is_central { self.jd_ut1_max } else { 0.0 });
        if matches!(self.kind, OccultationKind::Partial | OccultationKind::Grazing) {
            return [
                self.jd_ut1_max,
                t_noon,
                t1,
                t4,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ];
        }
        let t2 = self.t2_inner_immersion_jd_ut1.unwrap_or(t1);
        let t3 = self.t3_inner_emersion_jd_ut1.unwrap_or(t4);
        let tc1 = self
            .t_central_begin_jd_ut1
            .unwrap_or(if is_central { t2 } else { 0.0 });
        let tc2 = self
            .t_central_end_jd_ut1
            .unwrap_or(if is_central { t3 } else { 0.0 });
        [
            self.jd_ut1_max,
            t_noon,
            t1,
            t4,
            t2,
            t3,
            tc1,
            tc2,
            0.0,
            0.0,
        ]
    }

    /// Returns the Swiss Ephemeris return flag bitmask:
    /// `SE_ECL_TOTAL = 4`, `SE_ECL_ANNULAR = 8`, `SE_ECL_PARTIAL = 16`,
    /// `SE_ECL_CENTRAL = 1`, `SE_ECL_NONCENTRAL = 2`, `SE_ECL_VISIBLE = 128`.
    #[must_use]
    pub fn swe_retflag(&self) -> i32 {
        let mut flag = 0;
        match self.kind {
            OccultationKind::Total => flag |= SE_ECL_TOTAL,
            OccultationKind::Annular => flag |= SE_ECL_ANNULAR,
            OccultationKind::Partial | OccultationKind::Grazing => flag |= SE_ECL_PARTIAL,
        }
        if self.gamma.abs() <= 1.0 {
            flag |= SE_ECL_CENTRAL;
        } else {
            flag |= SE_ECL_NONCENTRAL;
        }
        if let Some(alt) = self.moon_altitude_deg {
            if alt >= -0.5 {
                flag |= SE_ECL_VISIBLE;
            }
        }
        flag
    }
}

/// Returns the nominal equatorial radius in km for a solar system body or fixed star.
#[must_use]
pub fn target_equatorial_radius_km(target: BodyId) -> f64 {
    match target {
        BodyId::Planet(p) => match p {
            PlanetId::Sun => 696_000.0,
            PlanetId::Mercury => 2_439.7,
            PlanetId::Venus => 6_051.8,
            PlanetId::Earth => R_EARTH_EQ_KM,
            PlanetId::Moon => R_MOON_KM,
            PlanetId::Mars => 3_396.2,
            PlanetId::Jupiter => 71_492.0,
            PlanetId::Saturn => 60_268.0,
            PlanetId::Uranus => 25_559.0,
            PlanetId::Neptune => 24_764.0,
            PlanetId::Pluto => 1_188.3,
            PlanetId::EarthMoonBarycenter | PlanetId::SolarSystemBarycenter => 0.0,
        },
        BodyId::MinorPlanet(mpc) => match mpc.0 {
            1 => 469.7,      // Ceres
            2 => 256.0,      // Pallas
            3 => 123.0,      // Juno
            4 => 262.7,      // Vesta
            10 => 215.0,     // Hygiea
            2060 => 109.0,   // Chiron
            90377 => 498.0,  // Sedna
            136108 => 816.0, // Haumea
            136199 => 1163.0,// Eris
            136472 => 715.0, // Makemake
            _ => 50.0,       // generic asteroid fallback radius
        },
        BodyId::FixedStar(_) => 0.0,
        BodyId::Hypothetical(_) | BodyId::LunarPoint(_) | BodyId::EclNut => 0.0,
    }
}

/// Validates that `target` is an occultable physical body.
fn validate_target(target: BodyId) -> Result<(), EphemerisError> {
    match target {
        BodyId::Planet(p) => match p {
            PlanetId::Moon => Err(EphemerisError::InvalidInput {
                message: "Target body cannot be the Moon in a lunar occultation calculation".to_owned(),
            }),
            PlanetId::Earth => Err(EphemerisError::InvalidInput {
                message: "Target body cannot be the Earth in a lunar occultation calculation".to_owned(),
            }),
            PlanetId::EarthMoonBarycenter | PlanetId::SolarSystemBarycenter => {
                Err(EphemerisError::InvalidInput {
                    message: format!("Target body cannot be a barycenter ({p:?})"),
                })
            }
            _ => Ok(()),
        },
        BodyId::LunarPoint(lp) => Err(EphemerisError::InvalidInput {
            message: format!("Target cannot be a synthetic lunar point ({lp:?})"),
        }),
        BodyId::EclNut => Err(EphemerisError::InvalidInput {
            message: "Target cannot be EclNut in an occultation calculation".to_owned(),
        }),
        _ => Ok(()),
    }
}

/// Validates search time range bounds.
fn validate_search_bounds(start_jd_ut1: f64, max_search_days: f64) -> Result<(), EphemerisError> {
    if !start_jd_ut1.is_finite() || start_jd_ut1 < -1_000_000.0 || start_jd_ut1 > 5_000_000.0 {
        return Err(EphemerisError::InvalidInput {
            message: format!("Invalid start_jd_ut1: {start_jd_ut1}"),
        });
    }
    if !max_search_days.is_finite() || max_search_days <= 0.0 || max_search_days > 36_525.0 {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid max_search_days: {max_search_days} (must be in 0.0..=36525.0)"
            ),
        });
    }
    Ok(())
}

fn tt_to_ut1(time_table: &TimeTable, jd_tt: f64, observer: Option<Observer>) -> f64 {
    let (ut1, _, _, _) = time_table.resolve_scales(TimeInput::from_tt_jd(jd_tt), observer);
    ut1.jd_f64()
}

fn ut1_to_tt(time_table: &TimeTable, jd_ut1: f64, observer: Option<Observer>) -> f64 {
    let (_, tt, _, _) = time_table.resolve_scales(TimeInput::from_ut1_jd(jd_ut1), observer);
    tt.jd_f64()
}

/// Instantaneous Besselian fundamental plane state for a lunar occultation.
#[derive(Debug, Clone, Copy)]
struct BesselianInstant {
    x: f64,
    y: f64,
    #[allow(dead_code)]
    z: f64,
    #[allow(dead_code)]
    sin_d: f64,
    #[allow(dead_code)]
    cos_d: f64,
    #[allow(dead_code)]
    mu_rad: f64,
    #[allow(dead_code)]
    tan_f1: f64,
    #[allow(dead_code)]
    tan_f2: f64,
    l1: f64,
    l2: f64,
    l2_signed: f64,
    rho1: f64,
}

impl BesselianInstant {
    /// Scaled shadow axis distance $\gamma_{\text{scaled}} = \sqrt{x^2 + (y/\rho_1)^2}$.
    fn gamma_scaled(&self) -> f64 {
        let y1 = self.y / self.rho1;
        (self.x * self.x + y1 * y1).sqrt()
    }

    /// Unscaled distance $\Delta = \sqrt{x^2 + y^2}$ and Earth limb radius $R_{\text{earth}}$ in direction $(x, y)$.
    fn delta_and_earth_radius(&self) -> (f64, f64) {
        let delta = (self.x * self.x + self.y * self.y).sqrt();
        let g = self.gamma_scaled();
        let r_earth = if g > 1e-12 { delta / g } else { 1.0 };
        (delta, r_earth)
    }
}

/// Evaluates the Besselian fundamental plane state at `jd_tt` for `target`.
fn eval_occultation_besselian(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    jd_tt: f64,
    target: BodyId,
) -> Result<BesselianInstant, EphemerisError> {
    let ctx = EpochContext::new(spk, time_table, TimeInput::from_tt_jd(jd_tt), None)?;
    let opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };

    let r_moon_au = eval_spk_target_vector_au(spk, &ctx, 301, &opts)?;
    let r_moon_km = vec_scale(&r_moon_au, AU_KM);
    let r_moon_er = vec_scale(&r_moon_km, 1.0 / R_EARTH_EQ_KM);

    let (k_hat, g_norm, k_targ) = match target {
        BodyId::FixedStar(hip) => {
            let star = find_star_by_hip(hip).ok_or_else(|| EphemerisError::UnsupportedBody {
                body: format!("FixedStar(HIP {})", hip.0),
            })?;
            let u_star = eval_star_unit_vector(star, &ctx, &opts);
            (u_star, f64::INFINITY, 0.0)
        }
        BodyId::Planet(p) => {
            let target_id = p.naif_target_id();
            let r_targ_au = eval_spk_target_vector_au(spk, &ctx, target_id, &opts)?;
            let r_targ_km = vec_scale(&r_targ_au, AU_KM);
            let r_targ_er = vec_scale(&r_targ_km, 1.0 / R_EARTH_EQ_KM);
            let g_vec = vec_sub(&r_targ_er, &r_moon_er);
            let g_len = vec_norm(&g_vec);
            let k = if g_len > 1e-12 {
                vec_scale(&g_vec, 1.0 / g_len)
            } else {
                [1.0, 0.0, 0.0]
            };
            let r_targ = target_equatorial_radius_km(target) / R_EARTH_EQ_KM;
            (k, g_len, r_targ)
        }
        BodyId::MinorPlanet(mpc) => {
            let (target_id, _, _) = resolve_spk_target(BodyId::MinorPlanet(mpc))?;
            let r_targ_au = eval_spk_target_vector_au(spk, &ctx, target_id, &opts)?;
            let r_targ_km = vec_scale(&r_targ_au, AU_KM);
            let r_targ_er = vec_scale(&r_targ_km, 1.0 / R_EARTH_EQ_KM);
            let g_vec = vec_sub(&r_targ_er, &r_moon_er);
            let g_len = vec_norm(&g_vec);
            let k = if g_len > 1e-12 {
                vec_scale(&g_vec, 1.0 / g_len)
            } else {
                [1.0, 0.0, 0.0]
            };
            let r_targ = target_equatorial_radius_km(target) / R_EARTH_EQ_KM;
            (k, g_len, r_targ)
        }
        other => {
            return Err(EphemerisError::InvalidInput {
                message: format!("Unsupported occultation target {other:?}"),
            });
        }
    };

    let sin_d = k_hat[2].clamp(-1.0, 1.0);
    let d_rad = sin_d.asin();
    let cos_d = d_rad.cos();
    let a_rad = k_hat[1].atan2(k_hat[0]).rem_euclid(std::f64::consts::TAU);
    let (sin_a, cos_a) = a_rad.sin_cos();

    let i_hat: Vec3 = [-sin_a, cos_a, 0.0];
    let j_hat: Vec3 = [-cos_a * sin_d, -sin_a * sin_d, cos_d];

    let x = vec_dot(&r_moon_er, &i_hat);
    let y = vec_dot(&r_moon_er, &j_hat);
    let z = vec_dot(&r_moon_er, &k_hat);

    let mu_rad = (ctx.gast_rad - a_rad).rem_euclid(std::f64::consts::TAU);

    let (tan_f1, tan_f2, l1, l2, l2_signed) = if !g_norm.is_finite()
        || g_norm > 1e12
        || (k_targ + K_MOON) / g_norm < 1e-12
    {
        (0.0, 0.0, K_MOON, K_MOON, -K_MOON)
    } else {
        let sin_f1 = ((k_targ + K_MOON) / g_norm).clamp(-0.999_999, 0.999_999);
        let sin_f2 = ((k_targ - K_MOON) / g_norm).clamp(-0.999_999, 0.999_999);
        let cos_f1 = (1.0 - sin_f1 * sin_f1).max(1e-12).sqrt();
        let cos_f2 = (1.0 - sin_f2 * sin_f2).max(1e-12).sqrt();
        let tf1 = sin_f1 / cos_f1;
        let tf2 = sin_f2 / cos_f2;
        let l1_val = z * tf1 + K_MOON / cos_f1;
        let l2_val = z * tf2 - K_MOON / cos_f2;
        (tf1, tf2, l1_val, l2_val.abs(), l2_val)
    };

    let rho1 = (1.0 - WGS84_E2 * cos_d * cos_d).max(1e-12).sqrt();

    Ok(BesselianInstant {
        x,
        y,
        z,
        sin_d,
        cos_d,
        mu_rad,
        tan_f1,
        tan_f2,
        l1,
        l2,
        l2_signed,
        rho1,
    })
}

/// Evaluates the signed geocentric Right Ascension difference $\Delta\alpha(t) = \alpha_{\text{Moon}} - \alpha_{\text{target}}$ in degrees.
fn geocentric_ra_diff_deg(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    jd_tt: f64,
    target: BodyId,
) -> Result<f64, EphemerisError> {
    let ctx = EpochContext::new(spk, time_table, TimeInput::from_tt_jd(jd_tt), None)?;
    let opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };
    let r_moon_au = eval_spk_target_vector_au(spk, &ctx, 301, &opts)?;
    let ra_moon_deg = r_moon_au[1].atan2(r_moon_au[0]).to_degrees().rem_euclid(360.0);

    let ra_targ_deg = match target {
        BodyId::FixedStar(hip) => {
            let star = find_star_by_hip(hip).ok_or_else(|| EphemerisError::UnsupportedBody {
                body: format!("FixedStar(HIP {})", hip.0),
            })?;
            let u_star = eval_star_unit_vector(star, &ctx, &opts);
            u_star[1].atan2(u_star[0]).to_degrees().rem_euclid(360.0)
        }
        BodyId::Planet(p) => {
            let r = eval_spk_target_vector_au(spk, &ctx, p.naif_target_id(), &opts)?;
            r[1].atan2(r[0]).to_degrees().rem_euclid(360.0)
        }
        BodyId::MinorPlanet(mpc) => {
            let (target_id, _, _) = resolve_spk_target(BodyId::MinorPlanet(mpc))?;
            let r = eval_spk_target_vector_au(spk, &ctx, target_id, &opts)?;
            r[1].atan2(r[0]).to_degrees().rem_euclid(360.0)
        }
        other => {
            return Err(EphemerisError::InvalidInput {
                message: format!("Unsupported occultation target {other:?}"),
            });
        }
    };

    let diff = (ra_moon_deg - ra_targ_deg + 180.0).rem_euclid(360.0) - 180.0;
    Ok(diff)
}

/// Refines a bracketed root $f(t) = 0$ using bisection to `tol_days`.
fn bisection_root<F>(
    mut f: F,
    mut a: f64,
    mut b: f64,
    tol_days: f64,
    max_iter: u32,
) -> Result<Option<f64>, EphemerisError>
where
    F: FnMut(f64) -> Result<f64, EphemerisError>,
{
    let mut fa = f(a)?;
    let fb = f(b)?;
    if fa.abs() < 1e-13 {
        return Ok(Some(a));
    }
    if fb.abs() < 1e-13 {
        return Ok(Some(b));
    }
    if fa * fb > 0.0 {
        return Ok(None);
    }
    for _ in 0..max_iter {
        let mid = 0.5 * (a + b);
        if (b - a).abs() <= tol_days {
            return Ok(Some(mid));
        }
        let fmid = f(mid)?;
        if fmid.abs() < 1e-13 {
            return Ok(Some(mid));
        }
        if fa * fmid <= 0.0 {
            b = mid;
        } else {
            a = mid;
            fa = fmid;
        }
    }
    Ok(Some(0.5 * (a + b)))
}

/// Finds the minimum of a scalar function $f(t)$ on $[a, b]$ via Golden Section search to `tol_days`.
fn golden_section_min<F>(
    mut f: F,
    mut a: f64,
    mut b: f64,
    tol_days: f64,
    max_iter: u32,
) -> Result<f64, EphemerisError>
where
    F: FnMut(f64) -> Result<f64, EphemerisError>,
{
    const R: f64 = 0.618_033_988_749_895; // (sqrt(5) - 1) / 2
    const C: f64 = 1.0 - R;
    let mut x1 = a + C * (b - a);
    let mut x2 = a + R * (b - a);
    let mut f1 = f(x1)?;
    let mut f2 = f(x2)?;
    for _ in 0..max_iter {
        if (b - a).abs() <= tol_days {
            return Ok(0.5 * (a + b));
        }
        if f1 < f2 {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + C * (b - a);
            f1 = f(x1)?;
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + R * (b - a);
            f2 = f(x2)?;
        }
    }
    Ok(0.5 * (a + b))
}

/// Searches for the next global lunar occultation of `target` on or after `start_jd_ut1` within `max_search_days`.
pub fn next_lunar_occultation(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    target: BodyId,
    start_jd_ut1: f64,
    max_search_days: f64,
) -> Result<Option<LunarOccultationInfo>, EphemerisError> {
    validate_target(target)?;
    validate_search_bounds(start_jd_ut1, max_search_days)?;

    let t_start_tt = ut1_to_tt(time_table, start_jd_ut1, None);
    let t_end_tt = t_start_tt + max_search_days;

    let pad_days = 0.25_f64;
    let t_scan_start = t_start_tt - pad_days;
    let t_scan_end = t_end_tt + pad_days;
    let scan_span = t_scan_end - t_scan_start;

    let step_days = 0.25_f64;
    let total_steps = ((scan_span / step_days).ceil() as usize).max(1);

    let mut t_prev = t_scan_start;
    let mut diff_prev = geocentric_ra_diff_deg(spk, time_table, t_prev, target)?;

    for step in 1..=total_steps {
        let t_curr = (t_scan_start + (step as f64) * step_days).min(t_scan_end);
        let diff_curr = geocentric_ra_diff_deg(spk, time_table, t_curr, target)?;

        // Conjunction occurred if diff crosses 0 from negative to positive
        if diff_prev <= 0.0 && diff_curr > 0.0 && (diff_curr - diff_prev) < 45.0 {
            // Refine exact RA conjunction time via bisection
            let mut f_ra = |t: f64| geocentric_ra_diff_deg(spk, time_table, t, target);
            if let Some(t_conj) = bisection_root(&mut f_ra, t_prev, t_curr, 1.0e-5, 50)? {
                // Find minimum separation t_max around t_conj on [t_conj - 0.12, t_conj + 0.12]
                let mut f_delta2 = |t: f64| -> Result<f64, EphemerisError> {
                    let b = eval_occultation_besselian(spk, time_table, t, target)?;
                    Ok(b.x * b.x + b.y * b.y)
                };
                let t_max =
                    golden_section_min(&mut f_delta2, t_conj - 0.12, t_conj + 0.12, 1.0e-8, 45)?;

                let jd_ut1_max = tt_to_ut1(time_table, t_max, None);
                // Only consider occultations where jd_ut1_max > start_jd_ut1 and within max_search_days
                if jd_ut1_max > start_jd_ut1 && jd_ut1_max <= start_jd_ut1 + max_search_days {
                    let b_max = eval_occultation_besselian(spk, time_table, t_max, target)?;
                    let (delta_max, r_earth_max) = b_max.delta_and_earth_radius();

                    // Check if shadow cone/cylinder intersects Earth anywhere
                    if delta_max <= r_earth_max + b_max.l1 {
                        let gamma_abs = delta_max;
                        let gamma = if b_max.y >= 0.0 {
                            gamma_abs
                        } else {
                            -gamma_abs
                        };

                        // Solve global external contact times t1, t4
                        let mut f_contact_ext = |t: f64| -> Result<f64, EphemerisError> {
                            let b = eval_occultation_besselian(spk, time_table, t, target)?;
                            let (d, re) = b.delta_and_earth_radius();
                            Ok(d - (re + b.l1))
                        };
                        let t1_tt = bisection_root(
                            &mut f_contact_ext,
                            t_max - 0.25,
                            t_max,
                            1.0e-8,
                            50,
                        )?;
                        let t4_tt = bisection_root(
                            &mut f_contact_ext,
                            t_max,
                            t_max + 0.25,
                            1.0e-8,
                            50,
                        )?;

                        // Solve meridian noon (b.x == 0.0) and central axis contacts (d == re)
                        let mut f_noon = |t: f64| -> Result<f64, EphemerisError> {
                            let b = eval_occultation_besselian(spk, time_table, t, target)?;
                            Ok(b.x)
                        };
                        let t_noon_tt =
                            bisection_root(&mut f_noon, t_max - 0.25, t_max + 0.25, 1.0e-8, 50)?;
                        let (tc1_tt, tc2_tt) = if delta_max < r_earth_max {
                            let mut f_contact_axis = |t: f64| -> Result<f64, EphemerisError> {
                                let b = eval_occultation_besselian(spk, time_table, t, target)?;
                                let (d, re) = b.delta_and_earth_radius();
                                Ok(d - re)
                            };
                            let tc1 = bisection_root(
                                &mut f_contact_axis,
                                t_max - 0.20,
                                t_max,
                                1.0e-8,
                                50,
                            )?;
                            let tc2 = bisection_root(
                                &mut f_contact_axis,
                                t_max,
                                t_max + 0.20,
                                1.0e-8,
                                50,
                            )?;
                            (tc1, tc2)
                        } else {
                            (None, None)
                        };

                        // Solve global internal contact times t2, t3 if total/annular
                        let is_star = matches!(target, BodyId::FixedStar(_));
                        let (kind, t2_tt, t3_tt, magnitude) = if is_star {
                            // Point star: R_star = 0, l1 = l2 = R_moon / R_earth ~ 0.2725
                            let grazing_margin = 0.005;
                            let k = if delta_max > r_earth_max + b_max.l2 - grazing_margin {
                                OccultationKind::Grazing
                            } else {
                                OccultationKind::Total
                            };
                            (k, t1_tt, t4_tt, 1.0)
                        } else if delta_max < r_earth_max + b_max.l2 {
                            let mut f_contact_int = |t: f64| -> Result<f64, EphemerisError> {
                                let b = eval_occultation_besselian(spk, time_table, t, target)?;
                                let (d, re) = b.delta_and_earth_radius();
                                Ok(d - (re + b.l2))
                            };
                            let t2 = bisection_root(
                                &mut f_contact_int,
                                t_max - 0.20,
                                t_max,
                                1.0e-8,
                                50,
                            )?;
                            let t3 = bisection_root(
                                &mut f_contact_int,
                                t_max,
                                t_max + 0.20,
                                1.0e-8,
                                50,
                            )?;
                            let k = if b_max.l2_signed <= 0.0 {
                                OccultationKind::Total
                            } else {
                                OccultationKind::Annular
                            };
                            let mag = (r_earth_max + b_max.l1 - delta_max) / (b_max.l1 + b_max.l2);
                            (k, t2, t3, mag)
                        } else if delta_max <= r_earth_max + b_max.l1 {
                            let k = if delta_max > r_earth_max + b_max.l1 - 0.02 {
                                OccultationKind::Grazing
                            } else {
                                OccultationKind::Partial
                            };
                            let mag = (r_earth_max + b_max.l1 - delta_max) / (b_max.l1 + b_max.l2);
                            (k, None, None, mag)
                        } else {
                            (OccultationKind::Partial, None, None, 0.0)
                        };

                        let t1_ut1 = t1_tt.map(|t| tt_to_ut1(time_table, t, None));
                        let t2_ut1 = t2_tt.map(|t| tt_to_ut1(time_table, t, None));
                        let t3_ut1 = t3_tt.map(|t| tt_to_ut1(time_table, t, None));
                        let t4_ut1 = t4_tt.map(|t| tt_to_ut1(time_table, t, None));
                        let t_noon_ut1 = t_noon_tt.map(|t| tt_to_ut1(time_table, t, None));
                        let tc1_ut1 = tc1_tt.map(|t| tt_to_ut1(time_table, t, None));
                        let tc2_ut1 = tc2_tt.map(|t| tt_to_ut1(time_table, t, None));

                        return Ok(Some(LunarOccultationInfo {
                            kind,
                            target,
                            jd_ut1_max,
                            jd_tt_max: t_max,
                            t1_first_contact_jd_ut1: t1_ut1,
                            t2_inner_immersion_jd_ut1: t2_ut1,
                            t3_inner_emersion_jd_ut1: t3_ut1,
                            t4_last_contact_jd_ut1: t4_ut1,
                            gamma,
                            magnitude,
                            moon_altitude_deg: None,
                            t_noon_jd_ut1: t_noon_ut1,
                            t_central_begin_jd_ut1: tc1_ut1,
                            t_central_end_jd_ut1: tc2_ut1,
                        }));
                    }
                }
            }
        }
        t_prev = t_curr;
        diff_prev = diff_curr;
    }

    Ok(None)
}

/// Topocentric angular circumstances between Moon and target.
#[derive(Debug, Clone, Copy)]
struct TopocentricCircumstance {
    sep_rad: f64,
    s_moon_rad: f64,
    s_targ_rad: f64,
    moon_alt_deg: f64,
    ra_diff_deg: f64,
}

/// Evaluates topocentric apparent separation and radii between Moon and target at `observer`.
fn eval_local_circumstance(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    jd_tt: f64,
    target: BodyId,
    observer: Observer,
) -> Result<TopocentricCircumstance, EphemerisError> {
    let ctx = EpochContext::new(
        spk,
        time_table,
        TimeInput::from_tt_jd(jd_tt),
        Some(observer),
    )?;
    let opts = CalcOptions {
        origin: Origin::Topocentric(observer),
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };

    let r_moon_au = eval_spk_target_vector_au(spk, &ctx, 301, &opts)?;
    let dist_moon_km = vec_norm(&r_moon_au) * AU_KM;
    let u_moon = vec_unit(&r_moon_au);
    let s_moon_rad = (R_MOON_KM / dist_moon_km).clamp(-1.0, 1.0).asin();

    // Topocentric apparent altitude of the Moon
    let r_moon_hz = crate::astrometry::equatorial_to_horizontal_vector_au(
        observer,
        ctx.gast_rad,
        &r_moon_au,
        AtmosphereParams::default(),
    );
    let sph_moon_hz = cart_to_spherical(&r_moon_hz);
    let moon_alt_deg = sph_moon_hz.lat_deg;

    let (u_targ, s_targ_rad) = match target {
        BodyId::FixedStar(hip) => {
            let star = find_star_by_hip(hip).ok_or_else(|| EphemerisError::UnsupportedBody {
                body: format!("FixedStar(HIP {})", hip.0),
            })?;
            let u_star = eval_star_unit_vector(star, &ctx, &opts);
            (u_star, 0.0)
        }
        BodyId::Planet(p) => {
            let r = eval_spk_target_vector_au(spk, &ctx, p.naif_target_id(), &opts)?;
            let d_km = vec_norm(&r) * AU_KM;
            let u = vec_unit(&r);
            let r_eq = target_equatorial_radius_km(target);
            let s = if d_km > 0.0 {
                (r_eq / d_km).clamp(-1.0, 1.0).asin()
            } else {
                0.0
            };
            (u, s)
        }
        BodyId::MinorPlanet(mpc) => {
            let (target_id, _, _) = resolve_spk_target(BodyId::MinorPlanet(mpc))?;
            let r = eval_spk_target_vector_au(spk, &ctx, target_id, &opts)?;
            let d_km = vec_norm(&r) * AU_KM;
            let u = vec_unit(&r);
            let r_eq = target_equatorial_radius_km(target);
            let s = if d_km > 0.0 {
                (r_eq / d_km).clamp(-1.0, 1.0).asin()
            } else {
                0.0
            };
            (u, s)
        }
        other => {
            return Err(EphemerisError::InvalidInput {
                message: format!("Unsupported occultation target {other:?}"),
            });
        }
    };

    let cos_sep = vec_dot(&u_moon, &u_targ).clamp(-1.0, 1.0);
    let sep_rad = cos_sep.acos();

    let ra_moon_deg = u_moon[1].atan2(u_moon[0]).to_degrees().rem_euclid(360.0);
    let ra_targ_deg = u_targ[1].atan2(u_targ[0]).to_degrees().rem_euclid(360.0);
    let ra_diff_deg = (ra_moon_deg - ra_targ_deg + 180.0).rem_euclid(360.0) - 180.0;

    Ok(TopocentricCircumstance {
        sep_rad,
        s_moon_rad,
        s_targ_rad,
        moon_alt_deg,
        ra_diff_deg,
    })
}

/// Searches for the next local topocentric lunar occultation of `target` on or after `start_jd_ut1`
/// within `max_search_days` at `observer`.
pub fn next_local_lunar_occultation(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    target: BodyId,
    start_jd_ut1: f64,
    max_search_days: f64,
    observer: Observer,
) -> Result<Option<LunarOccultationInfo>, EphemerisError> {
    validate_target(target)?;
    validate_search_bounds(start_jd_ut1, max_search_days)?;

    let t_start_tt = ut1_to_tt(time_table, start_jd_ut1, Some(observer));
    let t_end_tt = t_start_tt + max_search_days;

    let pad_days = 0.25_f64;
    let t_scan_start = t_start_tt - pad_days;
    let t_scan_end = t_end_tt + pad_days;
    let scan_span = t_scan_end - t_scan_start;

    let step_days = 0.25_f64;
    let total_steps = ((scan_span / step_days).ceil() as usize).max(1);

    let mut t_prev = t_scan_start;
    let mut diff_prev =
        eval_local_circumstance(spk, time_table, t_prev, target, observer)?.ra_diff_deg;

    for step in 1..=total_steps {
        let t_curr = (t_scan_start + (step as f64) * step_days).min(t_scan_end);
        let diff_curr =
            eval_local_circumstance(spk, time_table, t_curr, target, observer)?.ra_diff_deg;

        if diff_prev <= 0.0 && diff_curr > 0.0 && (diff_curr - diff_prev) < 45.0 {
            // Refine topocentric RA conjunction
            let mut f_ra = |t: f64| {
                eval_local_circumstance(spk, time_table, t, target, observer)
                    .map(|c| c.ra_diff_deg)
            };
            if let Some(t_conj) = bisection_root(&mut f_ra, t_prev, t_curr, 1.0e-5, 50)? {
                // Find minimum separation t_max around t_conj on [t_conj - 0.12, t_conj + 0.12]
                let mut f_sep = |t: f64| {
                    eval_local_circumstance(spk, time_table, t, target, observer)
                        .map(|c| c.sep_rad)
                };
                let t_max =
                    golden_section_min(&mut f_sep, t_conj - 0.12, t_conj + 0.12, 1.0e-8, 45)?;

                let jd_ut1_max = tt_to_ut1(time_table, t_max, Some(observer));
                // Only consider occultations where jd_ut1_max > start_jd_ut1 and within max_search_days
                if jd_ut1_max > start_jd_ut1 && jd_ut1_max <= start_jd_ut1 + max_search_days {
                    let c_max = eval_local_circumstance(spk, time_table, t_max, target, observer)?;
                    let sum_radii = c_max.s_moon_rad + c_max.s_targ_rad;

                    // Condition 1: centers approach closer than sum of semi-diameters
                    if c_max.sep_rad <= sum_radii {
                        let diff_radii = (c_max.s_moon_rad - c_max.s_targ_rad).abs();
                        let is_total = c_max.sep_rad < diff_radii;

                        let kind = if is_total {
                            if c_max.s_moon_rad >= c_max.s_targ_rad {
                                OccultationKind::Total
                            } else {
                                OccultationKind::Annular
                            }
                        } else {
                            let grazing_margin = 0.05 * c_max.s_targ_rad.max(0.002_f64.to_radians());
                            if c_max.sep_rad >= sum_radii - grazing_margin {
                                OccultationKind::Grazing
                            } else {
                                OccultationKind::Partial
                            }
                        };

                        // Solve for external contacts t1, t4
                        let mut f_c14 = |t: f64| {
                            let c = eval_local_circumstance(spk, time_table, t, target, observer)?;
                            Ok(c.sep_rad - (c.s_moon_rad + c.s_targ_rad))
                        };
                        let t1_tt = bisection_root(
                            &mut f_c14,
                            t_max - 0.08,
                            t_max,
                            1.0e-8,
                            50,
                        )?;
                        let t4_tt = bisection_root(
                            &mut f_c14,
                            t_max,
                            t_max + 0.08,
                            1.0e-8,
                            50,
                        )?;

                        let (t2_tt, t3_tt) = if is_total {
                            let mut f_c23 = |t: f64| {
                                let c = eval_local_circumstance(spk, time_table, t, target, observer)?;
                                Ok(c.sep_rad - (c.s_moon_rad - c.s_targ_rad).abs())
                            };
                            let t2 = bisection_root(
                                &mut f_c23,
                                t_max - 0.08,
                                t_max,
                                1.0e-8,
                                50,
                            )?;
                            let t3 = bisection_root(
                                &mut f_c23,
                                t_max,
                                t_max + 0.08,
                                1.0e-8,
                                50,
                            )?;
                            (t2, t3)
                        } else {
                            (None, None)
                        };

                        // Condition 2: Moon is above the observer's horizon (altitude >= -0.5°) at max or any contact
                        let mut visible_any = c_max.moon_alt_deg >= -0.5;
                        if !visible_any {
                            for t_c in [t1_tt, t2_tt, t3_tt, t4_tt].into_iter().flatten() {
                                if eval_local_circumstance(spk, time_table, t_c, target, observer)?
                                    .moon_alt_deg
                                    >= -0.5
                                {
                                    visible_any = true;
                                    break;
                                }
                            }
                        }

                        if visible_any {
                            let magnitude = if is_total && c_max.s_moon_rad >= c_max.s_targ_rad {
                                1.0
                            } else if c_max.s_targ_rad > 1e-9 {
                                (sum_radii - c_max.sep_rad) / (2.0 * c_max.s_targ_rad)
                            } else {
                                1.0
                            };

                            let jd_ut1_max = tt_to_ut1(time_table, t_max, Some(observer));
                            let t1_ut1 = t1_tt.map(|t| tt_to_ut1(time_table, t, Some(observer)));
                            let t2_ut1 = t2_tt.map(|t| tt_to_ut1(time_table, t, Some(observer)));
                            let t3_ut1 = t3_tt.map(|t| tt_to_ut1(time_table, t, Some(observer)));
                            let t4_ut1 = t4_tt.map(|t| tt_to_ut1(time_table, t, Some(observer)));

                            return Ok(Some(LunarOccultationInfo {
                                kind,
                                target,
                                jd_ut1_max,
                                jd_tt_max: t_max,
                                t1_first_contact_jd_ut1: t1_ut1,
                                t2_inner_immersion_jd_ut1: t2_ut1,
                                t3_inner_emersion_jd_ut1: t3_ut1,
                                t4_last_contact_jd_ut1: t4_ut1,
                                gamma: c_max.sep_rad.to_degrees(),
                                magnitude,
                                moon_altitude_deg: Some(c_max.moon_alt_deg),
                                t_noon_jd_ut1: None,
                                t_central_begin_jd_ut1: None,
                                t_central_end_jd_ut1: None,
                            }));
                        }
                    }
                }
            }
        }
        t_prev = t_curr;
        diff_prev = diff_curr;
    }

    Ok(None)
}

/// Finds all global lunar occultations of `target` within `[start_jd_ut1, end_jd_ut1]`.
pub fn lunar_occultations_in_range(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    target: BodyId,
    start_jd_ut1: f64,
    end_jd_ut1: f64,
) -> Result<Vec<LunarOccultationInfo>, EphemerisError> {
    if !start_jd_ut1.is_finite() || !end_jd_ut1.is_finite() || end_jd_ut1 < start_jd_ut1 {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid occultation search range: start_jd_ut1={start_jd_ut1}, end_jd_ut1={end_jd_ut1}"
            ),
        });
    }

    let start_tt = ut1_to_tt(time_table, start_jd_ut1, None);
    let end_tt = ut1_to_tt(time_table, end_jd_ut1, None);

    let mut results = Vec::new();
    let mut current_start = start_jd_ut1;

    while current_start <= end_jd_ut1 {
        let remaining_days = end_jd_ut1 - current_start;
        if remaining_days <= 0.0 {
            break;
        }
        let search_days = remaining_days.min(365.25);
        match next_lunar_occultation(spk, time_table, target, current_start, search_days)? {
            Some(info) => {
                if info.jd_tt_max >= start_tt && info.jd_tt_max <= end_tt {
                    results.push(info.clone());
                }
                if info.jd_tt_max > end_tt {
                    break;
                }
                current_start = (info.jd_ut1_max + 0.1).max(current_start + 0.25);
            }
            None => {
                current_start += search_days;
            }
        }
    }

    Ok(results)
}

/// Finds all local topocentric lunar occultations of `target` within `[start_jd_ut1, end_jd_ut1]` at `observer`.
pub fn local_lunar_occultations_in_range(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    target: BodyId,
    start_jd_ut1: f64,
    end_jd_ut1: f64,
    observer: Observer,
) -> Result<Vec<LunarOccultationInfo>, EphemerisError> {
    if !start_jd_ut1.is_finite() || !end_jd_ut1.is_finite() || end_jd_ut1 < start_jd_ut1 {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid occultation search range: start_jd_ut1={start_jd_ut1}, end_jd_ut1={end_jd_ut1}"
            ),
        });
    }

    let start_tt = ut1_to_tt(time_table, start_jd_ut1, Some(observer));
    let end_tt = ut1_to_tt(time_table, end_jd_ut1, Some(observer));

    let mut results = Vec::new();
    let mut current_start = start_jd_ut1;

    while current_start <= end_jd_ut1 {
        let remaining_days = end_jd_ut1 - current_start;
        if remaining_days <= 0.0 {
            break;
        }
        let search_days = remaining_days.min(365.25);
        match next_local_lunar_occultation(
            spk,
            time_table,
            target,
            current_start,
            search_days,
            observer,
        )? {
            Some(info) => {
                if info.jd_tt_max >= start_tt && info.jd_tt_max <= end_tt {
                    results.push(info.clone());
                }
                if info.jd_tt_max > end_tt {
                    break;
                }
                current_start = (info.jd_ut1_max + 0.1).max(current_start + 0.25);
            }
            None => {
                current_start += search_days;
            }
        }
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{HipId, LunarPointId, MpcNumber};

    #[test]
    fn test_target_validation() {
        assert!(validate_target(BodyId::Planet(PlanetId::Moon)).is_err());
        assert!(validate_target(BodyId::Planet(PlanetId::Earth)).is_err());
        assert!(validate_target(BodyId::Planet(PlanetId::EarthMoonBarycenter)).is_err());
        assert!(validate_target(BodyId::Planet(PlanetId::SolarSystemBarycenter)).is_err());
        assert!(validate_target(BodyId::LunarPoint(LunarPointId::MeanNode)).is_err());
        assert!(validate_target(BodyId::EclNut).is_err());

        assert!(validate_target(BodyId::Planet(PlanetId::Mars)).is_ok());
        assert!(validate_target(BodyId::Planet(PlanetId::Jupiter)).is_ok());
        assert!(validate_target(BodyId::MinorPlanet(MpcNumber(1))).is_ok());
        assert!(validate_target(BodyId::FixedStar(HipId(21421))).is_ok());
    }

    #[test]
    fn test_search_bounds_validation() {
        assert!(validate_search_bounds(f64::NAN, 30.0).is_err());
        assert!(validate_search_bounds(2451545.0, 0.0).is_err());
        assert!(validate_search_bounds(2451545.0, -10.0).is_err());
        assert!(validate_search_bounds(2451545.0, 40_000.0).is_err());
        assert!(validate_search_bounds(2451545.0, 365.25).is_ok());
    }

    #[test]
    fn test_target_radii() {
        assert!((target_equatorial_radius_km(BodyId::Planet(PlanetId::Sun)) - 696_000.0).abs() < 1e-6);
        assert!((target_equatorial_radius_km(BodyId::Planet(PlanetId::Mars)) - 3396.2).abs() < 1e-6);
        assert!((target_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(1))) - 469.7).abs() < 1e-6);
        assert!((target_equatorial_radius_km(BodyId::FixedStar(HipId(21421))) - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_lunar_occultation_info_helpers() {
        let info = LunarOccultationInfo {
            kind: OccultationKind::Total,
            target: BodyId::FixedStar(HipId(21421)),
            jd_ut1_max: 2451545.5,
            jd_tt_max: 2451545.5007,
            t1_first_contact_jd_ut1: Some(2451545.48),
            t2_inner_immersion_jd_ut1: Some(2451545.49),
            t3_inner_emersion_jd_ut1: Some(2451545.51),
            t4_last_contact_jd_ut1: Some(2451545.52),
            gamma: 0.45,
            magnitude: 1.15,
            moon_altitude_deg: Some(35.5),
            t_noon_jd_ut1: None,
            t_central_begin_jd_ut1: None,
            t_central_end_jd_ut1: None,
        };

        let tret = info.to_swe_tret10();
        assert!((tret[0] - 2451545.5).abs() < 1e-8);
        assert!((tret[1] - 2451545.48).abs() < 1e-8);
        assert!((tret[2] - 2451545.49).abs() < 1e-8);
        assert!((tret[3] - 2451545.51).abs() < 1e-8);
        assert!((tret[4] - 2451545.52).abs() < 1e-8);
        assert!((tret[5] - 0.0).abs() < 1e-8);

        let retflag = info.swe_retflag();
        // Total (4) | Central (1) | Visible (128) = 133
        assert_eq!(retflag & SE_ECL_TOTAL, SE_ECL_TOTAL);
        assert_eq!(retflag & SE_ECL_CENTRAL, SE_ECL_CENTRAL);
        assert_eq!(retflag & SE_ECL_VISIBLE, SE_ECL_VISIBLE);
    }

    #[test]
    fn test_noncentral_stellar_occultation_gamma_1_05() {
        let info = LunarOccultationInfo {
            kind: OccultationKind::Total,
            target: BodyId::FixedStar(HipId(21421)), // Aldebaran
            jd_ut1_max: 2451545.5,
            jd_tt_max: 2451545.5007,
            t1_first_contact_jd_ut1: Some(2451545.48),
            t2_inner_immersion_jd_ut1: None,
            t3_inner_emersion_jd_ut1: None,
            t4_last_contact_jd_ut1: Some(2451545.52),
            gamma: 1.05,
            magnitude: 1.0,
            moon_altitude_deg: None,
            t_noon_jd_ut1: None,
            t_central_begin_jd_ut1: None,
            t_central_end_jd_ut1: None,
        };

        let retflag = info.swe_retflag();
        // Non-central total stellar occultation: SE_ECL_TOTAL | SE_ECL_NONCENTRAL (4 | 2 = 6)
        assert_eq!(retflag & SE_ECL_TOTAL, SE_ECL_TOTAL);
        assert_eq!(retflag & SE_ECL_NONCENTRAL, SE_ECL_NONCENTRAL);
        assert_eq!(retflag & SE_ECL_CENTRAL, 0);
        assert_eq!(retflag & SE_ECL_PARTIAL, 0);
    }

    #[test]
    fn test_golden_section_and_bisection() -> Result<(), EphemerisError> {
        // Test parabolic minimum of (x - 3)^2 + 5 on [0, 5]
        let f_parab = |x: f64| -> Result<f64, EphemerisError> { Ok((x - 3.0) * (x - 3.0) + 5.0) };
        let min_x = golden_section_min(f_parab, 0.0, 5.0, 1e-8, 50)?;
        assert!((min_x - 3.0).abs() < 1e-6);

        // Test root of x^3 - 8 on [0, 4] -> root = 2
        let f_cube = |x: f64| -> Result<f64, EphemerisError> { Ok(x * x * x - 8.0) };
        let root = bisection_root(f_cube, 0.0, 4.0, 1e-8, 50)?;
        match root {
            Some(r) => assert!((r - 2.0).abs() < 1e-6),
            None => {
                return Err(EphemerisError::InvalidInput {
                    message: "bisection_root expected Some, got None".to_owned(),
                })
            }
        }
        Ok(())
    }
}
