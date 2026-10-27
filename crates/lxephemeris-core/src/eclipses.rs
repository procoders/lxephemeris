//! Solar and Lunar Eclipses (`ARCHITECTURE.md` §2.6):
//! - Besselian fundamental plane geometry from astrometric geocentric (`AbMode::Astrometric`)
//!   `DE440s` Sun (`10`) and Moon (`301`) states.
//! - Physical constants: $k_{\text{mean}} = 0.2725076$, $k_{\text{umbra}} = 0.272281$,
//!   WGS84 flattening $f = 1/298.257223563$, and atmospheric shadow enlargement via
//!   **Danjon ($1/85$)** (default) or traditional Chauvenet ($2\%$).
//! - Illinois bracketed root-finding to $10^{-10}\text{ d}$ ($< 10\text{ µs}$) for maximum eclipse,
//!   global contacts ($P_1, U_1, U_2, U_3, U_4, P_4$), central line coordinates, and local
//!   topocentric circumstances ($C_1, C_2, \text{Max}, C_3, C_4$).

use crate::astrometry::{azalt, eval_spk_target_vector_au, AtmosphereParams};
use crate::error::EphemerisError;
use crate::frame::{vec_dot, vec_norm, vec_scale, vec_sub, Vec3};
use crate::moment::EpochContext;
use crate::spk::SpkSource;
use crate::time::{TimeInput, TimeTable};
use crate::types::{AbMode, CalcOptions, Equinox, NutationMode, Observer, Origin, Plane, Zodiac, AU_KM};
use serde::{Deserialize, Serialize};

/// WGS84 Earth equatorial radius in km.
pub const R_EARTH_EQ_KM: f64 = 6_378.137;

/// WGS84 Earth flattening $f = 1 / 298.257223563$.
pub const WGS84_FLATTENING: f64 = 1.0 / 298.257_223_563;

/// WGS84 first eccentricity squared $e^2 = 2f - f^2$.
pub const WGS84_E2: f64 = 2.0 * WGS84_FLATTENING - WGS84_FLATTENING * WGS84_FLATTENING;

/// Solar nominal photospheric radius in km (IAU).
pub const R_SUN_KM: f64 = 696_000.0;

/// Solar radius in units of Earth equatorial radii (`R_SUN_KM / R_EARTH_EQ_KM`).
pub const K_SUN: f64 = R_SUN_KM / R_EARTH_EQ_KM;

/// Mean lunar radius in units of Earth equatorial radii (penumbral/partial solar & lunar contacts).
pub const K_MOON_MEAN: f64 = 0.272_507_6;

/// Lunar valley-floor umbral radius in units of Earth equatorial radii (IAU 1982 total/annular solar contacts).
pub const K_MOON_UMBRA: f64 = 0.272_281_0;

/// Earth atmospheric shadow enlargement rule for lunar eclipses (`ARCHITECTURE.md` §2.6).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum ShadowEnlargementRule {
    /// Danjon (1951) $1/85$ atmospheric enlargement of Earth's shadow radius (default).
    #[default]
    Danjon1Over85,
    /// Traditional Chauvenet $2\%$ ($1.02$) shadow enlargement rule.
    Traditional2Percent,
}

/// Options for solar and lunar eclipse calculations.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EclipseOptions {
    /// Shadow enlargement rule for lunar eclipses (`Danjon1Over85` default).
    pub shadow_rule: ShadowEnlargementRule,
    /// Aberration mode for Besselian fundamental plane evaluation (`AbMode::Apparent` default, or `AbMode::Astrometric`).
    pub aberration: AbMode,
    /// Maximum number of synodic months to search in `next_solar_eclipse` / `next_lunar_eclipse` (default `36`).
    pub max_synodic_months: u32,
}

impl Default for EclipseOptions {
    fn default() -> Self {
        Self {
            shadow_rule: ShadowEnlargementRule::Danjon1Over85,
            aberration: AbMode::Apparent,
            max_synodic_months: 36,
        }
    }
}

/// Classification of a solar eclipse.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SolarEclipseKind {
    Total,
    Annular,
    Hybrid,
    Partial,
}

/// Classification of a lunar eclipse.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum LunarEclipseKind {
    Total,
    Partial,
    Penumbral,
}

/// Local topocentric circumstances of a solar eclipse at a specific [`Observer`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LocalSolarEclipse {
    pub observer: Observer,
    /// `true` if the eclipse is visible above the observer's horizon at maximum local eclipse.
    pub visible_at_max: bool,
    /// `true` if the observer is inside the umbral/antumbral path (`Total` or `Annular` locally).
    pub is_central_locally: bool,
    /// Local eclipse type (`Total`, `Annular`, or `Partial`).
    pub local_kind: SolarEclipseKind,
    /// Instant of maximum local eclipse (`JD_TT`).
    pub max_eclipse_jd_tt: f64,
    /// Instant of maximum local eclipse (`JD_UT1`).
    pub max_eclipse_jd_ut1: f64,
    /// Maximum fraction of the Sun's diameter covered locally.
    pub magnitude: f64,
    /// Maximum fraction of the Sun's disk area covered locally `[0, 1]`.
    pub obscuration: f64,
    /// True altitude of the Sun in degrees at `max_eclipse_jd_tt`.
    pub sun_altitude_deg: f64,
    /// Azimuth of the Sun in degrees `[0, 360)` at `max_eclipse_jd_tt`.
    pub sun_azimuth_deg: f64,
    /// First contact $C_1$ (`JD_TT`, start of partial phase).
    pub c1_jd_tt: Option<f64>,
    /// Second contact $C_2$ (`JD_TT`, start of totality/annularity, if central locally).
    pub c2_jd_tt: Option<f64>,
    /// Third contact $C_3$ (`JD_TT`, end of totality/annularity, if central locally).
    pub c3_jd_tt: Option<f64>,
    /// Fourth contact $C_4$ (`JD_TT`, end of partial phase).
    pub c4_jd_tt: Option<f64>,
    /// Duration of local totality or annularity in seconds (`0.0` if partial).
    pub central_duration_sec: f64,
}

/// Global and optional local circumstances of a Solar Eclipse.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolarEclipseInfo {
    pub kind: SolarEclipseKind,
    /// Saros series number (e.g., `145` for 1999-08-11 and 2017-08-21, `139` for 2024-04-08, `R8.1`).
    pub saros_series: i32,
    /// Member ordinal within the Saros series (`1..=85`, `R8.1`).
    pub saros_member: u32,
    /// Instant of greatest global eclipse (`JD_TT`).
    pub max_eclipse_jd_tt: f64,
    /// Instant of greatest global eclipse (`JD_UT1`).
    pub max_eclipse_jd_ut1: f64,
    /// $\Delta T = \text{TT} - \text{UT1}$ in seconds at maximum eclipse.
    pub delta_t_sec: f64,
    /// Minimum distance $\gamma$ of the shadow axis from the Earth's center in Earth equatorial radii
    /// (positive North, negative South).
    pub gamma: f64,
    /// Umbral radius $u = L_2(t_{\text{max}})$ in Earth equatorial radii (negative for Total).
    pub u_radius_earth: f64,
    /// Penumbral radius $L_1(t_{\text{max}})$ in Earth equatorial radii.
    pub penumbra_radius_earth: f64,
    /// Maximum eclipse magnitude.
    pub magnitude: f64,
    /// Geographic latitude of the central line at `max_eclipse_jd_tt` (`None` if non-central `Partial`).
    pub central_latitude_deg: Option<f64>,
    /// Geographic longitude `[-180, +180]` of the central line at `max_eclipse_jd_tt` (`None` if `Partial`).
    pub central_longitude_deg: Option<f64>,
    /// Duration of totality or annularity in seconds on the central line at `max_eclipse_jd_tt` (`None` if `Partial`, `R8.1`).
    pub central_duration_sec: Option<f64>,
    /// Width of the central umbral/antumbral shadow path in kilometers at `max_eclipse_jd_tt` (`None` if `Partial`, `R8.1`).
    pub path_width_km: Option<f64>,
    /// First external penumbral contact $P_1$ (`JD_TT`).
    pub p1_jd_tt: Option<f64>,
    /// First external umbral contact $U_1$ (`JD_TT`).
    pub u1_jd_tt: Option<f64>,
    /// First internal umbral contact $U_2$ (`JD_TT`).
    pub u2_jd_tt: Option<f64>,
    /// Last internal umbral contact $U_3$ (`JD_TT`).
    pub u3_jd_tt: Option<f64>,
    /// Last external umbral contact $U_4$ (`JD_TT`).
    pub u4_jd_tt: Option<f64>,
    /// Last external penumbral contact $P_4$ (`JD_TT`).
    pub p4_jd_tt: Option<f64>,
    /// Meridian noon contact (`x(t) = 0`, `JD_TT`).
    #[serde(default)]
    pub t_noon_jd_tt: Option<f64>,
    /// Central axis first contact (`d = R_earth`, `JD_TT`).
    #[serde(default)]
    pub c_begin_jd_tt: Option<f64>,
    /// Central axis last contact (`d = R_earth`, `JD_TT`).
    #[serde(default)]
    pub c_end_jd_tt: Option<f64>,
    /// Optional local topocentric circumstances if an [`Observer`] was supplied.
    pub local: Option<LocalSolarEclipse>,
}

/// Complete circumstances and contact times of a Lunar Eclipse.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LunarEclipseInfo {
    pub kind: LunarEclipseKind,
    pub shadow_rule: ShadowEnlargementRule,
    /// Saros series number (e.g., `126` for 2003-11-09, `123` for 2025-03-14, `R8.1`).
    pub saros_series: i32,
    /// Member ordinal within the Saros series (`1..=85`, `R8.1`).
    pub saros_member: u32,
    /// Instant of greatest eclipse (`JD_TT`).
    pub max_eclipse_jd_tt: f64,
    /// Instant of greatest eclipse (`JD_UT1`).
    pub max_eclipse_jd_ut1: f64,
    /// $\Delta T = \text{TT} - \text{UT1}$ in seconds at maximum eclipse.
    pub delta_t_sec: f64,
    /// Minimum distance $\gamma$ of the Moon's center from the Earth's shadow axis in Earth equatorial radii.
    pub gamma: f64,
    /// Penumbral eclipse magnitude at greatest eclipse.
    pub penumbral_magnitude: f64,
    /// Umbral eclipse magnitude at greatest eclipse (negative for purely penumbral eclipses).
    pub umbral_magnitude: f64,
    /// Penumbral semi-duration in minutes (`0.5 * (P4 - P1)`).
    pub penumbral_semi_duration_min: f64,
    /// Partial (umbral) semi-duration in minutes (`0.5 * (U4 - U1)`, `0.0` if penumbral).
    pub partial_semi_duration_min: f64,
    /// Total semi-duration in minutes (`0.5 * (U3 - U2)`, `0.0` if not total).
    pub total_semi_duration_min: f64,
    /// First penumbral contact $P_1$ (`JD_TT`).
    pub p1_jd_tt: Option<f64>,
    /// First umbral contact $U_1$ (`JD_TT`, `None` if penumbral).
    pub u1_jd_tt: Option<f64>,
    /// Start of totality $U_2$ (`JD_TT`, `None` if not total).
    pub u2_jd_tt: Option<f64>,
    /// End of totality $U_3$ (`JD_TT`, `None` if not total).
    pub u3_jd_tt: Option<f64>,
    /// Last umbral contact $U_4$ (`JD_TT`, `None` if penumbral).
    pub u4_jd_tt: Option<f64>,
    /// Last penumbral contact $P_4$ (`JD_TT`).
    pub p4_jd_tt: Option<f64>,
}

/// Instantaneous Besselian fundamental plane state for a Solar Eclipse at `jd_tt`.
#[derive(Debug, Clone, Copy)]
struct SolarBesselianInstant {
    x: f64,
    y: f64,
    #[allow(dead_code)]
    z: f64,
    sin_d: f64,
    cos_d: f64,
    mu_rad: f64,
    tan_f1: f64,
    tan_f2: f64,
    l1: f64,
    l2: f64,
    rho1: f64,
    rho2: f64,
}

impl SolarBesselianInstant {
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

    /// Evaluates central line `(lat_deg, lon_deg, zeta, L1, L2)` if `gamma_scaled() <= 1.0`.
    fn central_line(&self) -> Option<(f64, f64, f64, f64, f64)> {
        let y1 = self.y / self.rho1;
        let r2 = self.x * self.x + y1 * y1;
        if r2 > 1.0 {
            return None;
        }
        let zeta1 = (1.0 - r2).max(0.0).sqrt();
        let sin_d1 = self.sin_d / self.rho1;
        let cos_d1 = ((1.0 - WGS84_E2).sqrt() * self.cos_d) / self.rho1;

        let sin_phi1 = y1 * cos_d1 + zeta1 * sin_d1;
        let cos_phi1_sin_h = self.x;
        let cos_phi1_cos_h = -y1 * sin_d1 + zeta1 * cos_d1;

        let h_rad = cos_phi1_sin_h.atan2(cos_phi1_cos_h);
        let cos_phi1 = (cos_phi1_sin_h * cos_phi1_sin_h + cos_phi1_cos_h * cos_phi1_cos_h)
            .sqrt()
            .max(1e-15);
        let tan_phi_geod = (sin_phi1 / cos_phi1) / (1.0 - WGS84_E2).sqrt();
        let lat_deg = tan_phi_geod.atan().to_degrees();
        let lon_east_deg = (h_rad - self.mu_rad).to_degrees();
        let lon_norm_deg = (lon_east_deg + 180.0).rem_euclid(360.0) - 180.0;

        let cos_dd = cos_d1 * self.cos_d + sin_d1 * self.sin_d;
        let sin_dd = sin_d1 * self.cos_d - cos_d1 * self.sin_d;
        let zeta = self.rho2 * (zeta1 * cos_dd - y1 * sin_dd);
        let big_l1 = self.l1 - zeta * self.tan_f1;
        let big_l2 = self.l2 - zeta * self.tan_f2;

        Some((lat_deg, lon_norm_deg, zeta, big_l1, big_l2))
    }

    /// Evaluates topocentric observer coordinates $(\xi, \eta, \zeta)$ and shadow distance $m_{\text{loc}}$.
    fn local_observer(&self, obs: Observer) -> (f64, f64, f64, f64, f64, f64) {
        let phi_rad = obs.latitude_deg.to_radians();
        let lon_rad = obs.longitude_deg.to_radians();
        let (sin_phi, cos_phi) = phi_rad.sin_cos();
        let u = ((1.0 - WGS84_E2).sqrt() * phi_rad.tan()).atan();
        let (sin_u, cos_u) = u.sin_cos();
        let h_er = (obs.altitude_m / 1000.0) / R_EARTH_EQ_KM;
        let rho_sin_phip = (1.0 - WGS84_E2).sqrt() * sin_u + h_er * sin_phi;
        let rho_cos_phip = cos_u + h_er * cos_phi;

        let theta = self.mu_rad + lon_rad;
        let (sin_t, cos_t) = theta.sin_cos();

        let xi = rho_cos_phip * sin_t;
        let eta = rho_sin_phip * self.cos_d - rho_cos_phip * cos_t * self.sin_d;
        let zeta = rho_sin_phip * self.sin_d + rho_cos_phip * cos_t * self.cos_d;

        let dx = self.x - xi;
        let dy = self.y - eta;
        let m_loc = (dx * dx + dy * dy).sqrt();
        let big_l1 = self.l1 - zeta * self.tan_f1;
        let big_l2 = self.l2 - zeta * self.tan_f2;
        (xi, eta, zeta, m_loc, big_l1, big_l2)
    }
}

/// Evaluates the [`SolarBesselianInstant`] at `jd_tt` using geocentric Sun and Moon vectors with `aberration`.
fn eval_solar_besselian(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    jd_tt: f64,
    aberration: AbMode,
) -> Result<SolarBesselianInstant, EphemerisError> {
    let ctx = EpochContext::new(spk, time_table, TimeInput::from_tt_jd(jd_tt), None)?;
    let opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };
    let r_sun_au = eval_spk_target_vector_au(spk, &ctx, 10, &opts)?;
    let r_moon_au = eval_spk_target_vector_au(spk, &ctx, 301, &opts)?;

    // Convert from AU to Earth equatorial radii R_EARTH_EQ_KM
    let scale_er = AU_KM / R_EARTH_EQ_KM;
    let r_sun = vec_scale(&r_sun_au, scale_er);
    let r_moon = vec_scale(&r_moon_au, scale_er);

    let g_vec = vec_sub(&r_sun, &r_moon);
    let g_norm = vec_norm(&g_vec);
    let k_hat = vec_scale(&g_vec, 1.0 / g_norm);

    let sin_d = k_hat[2].clamp(-1.0, 1.0);
    let d_rad = sin_d.asin();
    let cos_d = d_rad.cos();
    let a_rad = k_hat[1].atan2(k_hat[0]).rem_euclid(std::f64::consts::TAU);
    let (sin_a, cos_a) = a_rad.sin_cos();

    let i_hat: Vec3 = [-sin_a, cos_a, 0.0];
    let j_hat: Vec3 = [-cos_a * sin_d, -sin_a * sin_d, cos_d];

    let x = vec_dot(&r_moon, &i_hat);
    let y = vec_dot(&r_moon, &j_hat);
    let z = vec_dot(&r_moon, &k_hat);

    let mu_rad = (ctx.gast_rad - a_rad).rem_euclid(std::f64::consts::TAU);

    let sin_f1 = ((K_SUN + K_MOON_MEAN) / g_norm).clamp(-0.999_999, 0.999_999);
    let sin_f2 = ((K_SUN - K_MOON_UMBRA) / g_norm).clamp(-0.999_999, 0.999_999);
    let tan_f1 = sin_f1 / (1.0 - sin_f1 * sin_f1).sqrt();
    let tan_f2 = sin_f2 / (1.0 - sin_f2 * sin_f2).sqrt();

    let c1 = z + K_MOON_MEAN / sin_f1;
    let c2 = z - K_MOON_UMBRA / sin_f2;
    let l1 = c1 * tan_f1;
    let l2 = c2 * tan_f2;

    let rho1 = (1.0 - WGS84_E2 * cos_d * cos_d).sqrt();
    let rho2 = (1.0 - WGS84_E2 * sin_d * sin_d).sqrt();

    Ok(SolarBesselianInstant {
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
        rho1,
        rho2,
    })
}

/// Refines a scalar root $f(t) = 0$ on `[a, b]` using the Illinois method with bisection safeguard to `tol_days`.
fn illinois_root<F>(
    mut f: F,
    mut a: f64,
    mut b: f64,
    tol_days: f64,
    max_iter: u32,
) -> Result<f64, EphemerisError>
where
    F: FnMut(f64) -> Result<f64, EphemerisError>,
{
    let mut fa = f(a)?;
    let mut fb = f(b)?;
    if fa.abs() < 1e-13 {
        return Ok(a);
    }
    if fb.abs() < 1e-13 {
        return Ok(b);
    }

    let eff_tol = tol_days.max(4.0 * f64::EPSILON * a.abs().max(b.abs()));
    let mut side = 0i8;
    for _ in 0..max_iter {
        let width = (b - a).abs();
        if width <= eff_tol || fb.abs() < 1e-13 {
            return Ok(0.5 * (a + b));
        }
        let denom = fb - fa;
        let mut c = if denom.abs() > 1e-15 {
            (a * fb - b * fa) / denom
        } else {
            0.5 * (a + b)
        };
        let lo = a.min(b);
        let hi = a.max(b);
        let margin = 0.05 * width;
        if c <= lo + margin || c >= hi - margin {
            c = 0.5 * (a + b);
        }
        let fc = f(c)?;
        if fc.abs() < 1e-13 {
            return Ok(c);
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
    Err(EphemerisError::NonConvergence {
        solver: "illinois_root".to_owned(),
        iterations: max_iter,
        residual: fb.abs(),
    })
}

/// Finds a root of `f(t) = 0` on `[t_start, t_end]` by scanning `steps` sub-intervals and refining via Illinois.
fn scan_and_refine_root<F>(
    mut f: F,
    t_start: f64,
    t_end: f64,
    steps: usize,
    tol_days: f64,
) -> Result<Option<f64>, EphemerisError>
where
    F: FnMut(f64) -> Result<f64, EphemerisError>,
{
    let dt = (t_end - t_start) / (steps as f64);
    let mut t_prev = t_start;
    let mut f_prev = f(t_prev)?;
    for i in 1..=steps {
        let t_curr = t_start + (i as f64) * dt;
        let f_curr = f(t_curr)?;
        if f_prev * f_curr <= 0.0 {
            let root = illinois_root(&mut f, t_prev, t_curr, tol_days, 80)?;
            return Ok(Some(root));
        }
        t_prev = t_curr;
        f_prev = f_curr;
    }
    Ok(None)
}

/// Exact NASA GSFC Five Millennium Canon reference member offsets $M_{\text{solar},0}(s) = 1 - \Delta i_{\text{first}}(s)$
/// for all Solar Saros series $s \in [0, 180]$ (`-2000 CE .. +3000 CE`, anchored at $\Delta i = 0 \iff k = -5 + 358(s - 145)$).
const SOLAR_SAROS_M0_0_TO_180: [i64; 181] = [
    43, 40, 41, 40, 37, 38, 37, 34, 35, 36, 32, 35, 46, 41, 43, 45, 41, 41, 42, 39, 38, 39, 35, 35,
    36, 32, 32, 33, 30, 30, 41, 29, 39, 42, 40, 39, 40, 38, 36, 37, 35, 33, 34, 32, 30, 31, 29, 27,
    30, 27, 26, 39, 39, 35, 37, 37, 34, 35, 34, 31, 32, 31, 28, 29, 28, 25, 27, 26, 23, 30, 37, 31,
    35, 35, 32, 33, 33, 29, 30, 30, 26, 27, 27, 24, 24, 25, 21, 23, 34, 21, 31, 34, 31, 30, 31, 29,
    28, 28, 26, 25, 25, 23, 22, 23, 20, 20, 24, 20, 22, 31, 30, 28, 29, 28, 26, 27, 25, 23, 24, 23,
    20, 21, 20, 17, 20, 21, 16, 28, 30, 25, 27, 27, 24, 25, 25, 22, 22, 22, 19, 19, 20, 16, 17, 19,
    14, 21, 28, 25, 25, 26, 24, 23, 23, 21, 20, 21, 18, 17, 18, 16, 15, 17, 14, 14, 26, 25, 22, 25,
    24, 21, 22, 21, 18, 19, 18, 16, 17, 16, 13, 15, 15,
];

/// Exact NASA GSFC Five Millennium Canon reference member offsets $M_{\text{lunar},0}(s) = 1 - \Delta i_{\text{first}}(s)$
/// for all Lunar Saros series $s \in [1, 180]$ (`-2000 CE .. +3000 CE`, anchored at $\Delta i = 0 \iff k = 47 + 358(s - 126)$).
const LUNAR_SAROS_M0_1_TO_180: [i64; 180] = [
    54, 53, 57, 63, 54, 65, 65, 61, 63, 62, 59, 60, 59, 56, 57, 56, 53, 54, 53, 50, 52, 52, 48, 61,
    63, 58, 60, 60, 57, 58, 58, 54, 55, 54, 51, 52, 52, 48, 49, 50, 46, 48, 60, 47, 57, 59, 56, 55,
    56, 53, 53, 53, 50, 49, 50, 47, 47, 48, 44, 45, 51, 44, 51, 56, 55, 53, 54, 52, 50, 51, 50, 47,
    48, 47, 45, 46, 44, 42, 44, 44, 41, 53, 54, 50, 52, 52, 48, 49, 49, 46, 46, 46, 43, 43, 43, 40,
    41, 43, 38, 46, 52, 48, 49, 50, 47, 47, 48, 45, 44, 45, 42, 42, 42, 39, 39, 41, 37, 38, 49, 47,
    46, 49, 47, 45, 46, 44, 43, 43, 42, 40, 41, 39, 37, 39, 37, 35, 43, 47, 41, 46, 47, 43, 44, 44,
    41, 42, 41, 38, 39, 39, 36, 37, 37, 33, 36, 46, 34, 44, 46, 42, 43, 43, 40, 40, 41, 38, 37, 38,
    35, 35, 36, 33, 33, 43, 36, 40, 44, 42, 41, 42,
];

/// Computes the exact `(saros_series, saros_member)` for a Solar Eclipse via the van den Bergh / Meeus
/// Saros–Inex congruence $\Delta k = 223\,\Delta i + 358\,\Delta s$ ($38 \times 358 - 61 \times 223 = 1$,
/// anchored to `1999-08-11` = lunation $k_0 = -5$, Solar Saros `145`, member `21`).
#[must_use]
pub fn compute_solar_saros(max_eclipse_jd_tt: f64, gamma: f64) -> (i32, u32) {
    if !max_eclipse_jd_tt.is_finite() || !gamma.is_finite() || max_eclipse_jd_tt.abs() > 1.0e9 {
        return (0, 1);
    }
    const SYNODIC_MONTH_DAYS: f64 = 29.530_588_853;
    const JD_NM_2000_01_06: f64 = 2_451_550.26;
    let k = ((max_eclipse_jd_tt - JD_NM_2000_01_06) / SYNODIC_MONTH_DAYS).round() as i64;
    let dk = k.saturating_sub(-5); // Relative to 1999-08-11 (k0 = -5, Saros 145, member 21)
    let ds_mod = dk.saturating_mul(38).rem_euclid(223);
    let m = (((dk as f64 / 358.0) - (ds_mod as f64)) / 223.0).round() as i64;
    let ds = ds_mod.saturating_add(m.saturating_mul(223));
    let saros_series = 145_i64
        .saturating_add(ds)
        .clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
    let di = dk.saturating_sub(ds.saturating_mul(358)) / 223;
    let saros_member = if let Some(&m0) = usize::try_from(saros_series)
        .ok()
        .and_then(|idx| SOLAR_SAROS_M0_0_TO_180.get(idx))
    {
        m0.saturating_add(di).clamp(1, 90) as u32
    } else {
        let sign = if saros_series.rem_euclid(2) != 0 {
            1.0
        } else {
            -1.0
        };
        let gamma_0 = gamma.clamp(-10.0, 10.0) + sign * 0.0422 * (di as f64);
        let member_0 = (((1.56 - sign * gamma_0) / 0.0422).round() as i64).saturating_add(1);
        member_0.saturating_add(di).clamp(1, 90) as u32
    };
    (saros_series, saros_member)
}

/// Computes the exact `(saros_series, saros_member)` for a Lunar Eclipse via the van den Bergh / Meeus
/// Saros–Inex congruence $\Delta k = 223\,\Delta i + 358\,\Delta s$
/// (anchored to `2003-11-09` = Full Moon lunation $k_0 = 47$, Lunar Saros `126`, member `44`).
#[must_use]
pub fn compute_lunar_saros(max_eclipse_jd_tt: f64, gamma: f64) -> (i32, u32) {
    if !max_eclipse_jd_tt.is_finite() || !gamma.is_finite() || max_eclipse_jd_tt.abs() > 1.0e9 {
        return (0, 1);
    }
    const SYNODIC_MONTH_DAYS: f64 = 29.530_588_853;
    const JD_FM_2000_01_21: f64 = 2_451_564.70;
    let k = ((max_eclipse_jd_tt - JD_FM_2000_01_21) / SYNODIC_MONTH_DAYS).round() as i64;
    let dk = k.saturating_sub(47); // Relative to 2003-11-09 (k0 = 47, Saros 126, member 44)
    let ds_mod = dk.saturating_mul(38).rem_euclid(223);
    let m = (((dk as f64 / 358.0) - (ds_mod as f64)) / 223.0).round() as i64;
    let ds = ds_mod.saturating_add(m.saturating_mul(223));
    let saros_series = 126_i64
        .saturating_add(ds)
        .clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
    let di = dk.saturating_sub(ds.saturating_mul(358)) / 223;
    let saros_member = if let Some(&m0) = usize::try_from(saros_series - 1)
        .ok()
        .and_then(|idx| LUNAR_SAROS_M0_1_TO_180.get(idx))
    {
        m0.saturating_add(di).clamp(1, 90) as u32
    } else {
        // Odd Lunar Saros series move South-to-North (d(gamma)/d(di) > 0), opposite to Solar Saros series
        let sign = if saros_series.rem_euclid(2) != 0 {
            -1.0
        } else {
            1.0
        };
        let gamma_0 = gamma.clamp(-10.0, 10.0) + sign * 0.0422 * (di as f64);
        let member_0 = (((1.56 - sign * gamma_0) / 0.0422).round() as i64).saturating_add(1);
        member_0.saturating_add(di).clamp(1, 90) as u32
    };
    (saros_series, saros_member)
}

/// Evaluates whether a Solar Eclipse occurs near `approx_syzygy_jd_tt`, returning `Ok(Some(SolarEclipseInfo))`
/// if an eclipse takes place, or `Ok(None)` if the New Moon misses the Earth's shadow cone.
pub fn solar_eclipse_at(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    approx_syzygy_jd_tt: f64,
    observer: Option<Observer>,
    options: EclipseOptions,
) -> Result<Option<SolarEclipseInfo>, EphemerisError> {
    let ab = options.aberration;
    // 1. Parabolic 3-point pre-estimate of the minimum of Delta^2(t) = x(t)^2 + y(t)^2
    let h_pre = 0.10;
    let b_m = eval_solar_besselian(spk, time_table, approx_syzygy_jd_tt - h_pre, ab)?;
    let b_0 = eval_solar_besselian(spk, time_table, approx_syzygy_jd_tt, ab)?;
    let b_p = eval_solar_besselian(spk, time_table, approx_syzygy_jd_tt + h_pre, ab)?;
    let fm = b_m.x * b_m.x + b_m.y * b_m.y;
    let f0 = b_0.x * b_0.x + b_0.y * b_0.y;
    let fp = b_p.x * b_p.x + b_p.y * b_p.y;
    let denom = fp - 2.0 * f0 + fm;
    let dt_est = if denom > 1e-12 {
        (-0.5 * h_pre * (fp - fm) / denom).clamp(-0.55, 0.55)
    } else {
        0.0
    };
    let t_est = approx_syzygy_jd_tt + dt_est;
    let b_est = eval_solar_besselian(spk, time_table, t_est, ab)?;
    if b_est.x.hypot(b_est.y) > 1.75 {
        return Ok(None);
    }

    // Refine exact minimum of Delta^2(t) = x(t)^2 + y(t)^2 around t_est
    let h_deriv = 1.0e-4;
    let d_delta2 = |t: f64| -> Result<f64, EphemerisError> {
        let bm = eval_solar_besselian(spk, time_table, t - h_deriv, ab)?;
        let bp = eval_solar_besselian(spk, time_table, t + h_deriv, ab)?;
        let dm2 = bm.x * bm.x + bm.y * bm.y;
        let dp2 = bp.x * bp.x + bp.y * bp.y;
        Ok(dp2 - dm2)
    };

    let Some(t_max) = scan_and_refine_root(
        d_delta2,
        t_est - 0.03,
        t_est + 0.03,
        4,
        1.0e-10,
    )?
    else {
        return Ok(None);
    };

    let b_max = eval_solar_besselian(spk, time_table, t_max, ab)?;
    let (delta_max, r_earth_max) = b_max.delta_and_earth_radius();

    // Check if penumbra touches Earth at maximum eclipse
    if delta_max >= r_earth_max + b_max.l1 {
        return Ok(None);
    }

    let gamma_abs = delta_max;
    let gamma = if b_max.y >= 0.0 { gamma_abs } else { -gamma_abs };

    // Determine global contacts P1, P4, U1, U4, U2, U3
    let f_p = |t: f64| -> Result<f64, EphemerisError> {
        let b = eval_solar_besselian(spk, time_table, t, ab)?;
        let (d, re) = b.delta_and_earth_radius();
        Ok(d - (re + b.l1))
    };
    let p1_jd_tt = scan_and_refine_root(f_p, t_max - 0.25, t_max, 8, 1.0e-10)?;
    let p4_jd_tt = scan_and_refine_root(f_p, t_max, t_max + 0.25, 8, 1.0e-10)?;

    let f_u_ext = |t: f64| -> Result<f64, EphemerisError> {
        let b = eval_solar_besselian(spk, time_table, t, ab)?;
        let (d, re) = b.delta_and_earth_radius();
        Ok(d - (re + b.l2.abs()))
    };
    let u1_jd_tt = if delta_max < r_earth_max + b_max.l2.abs() {
        scan_and_refine_root(f_u_ext, t_max - 0.20, t_max, 8, 1.0e-10)?
    } else {
        None
    };
    let u4_jd_tt = if delta_max < r_earth_max + b_max.l2.abs() {
        scan_and_refine_root(f_u_ext, t_max, t_max + 0.20, 8, 1.0e-10)?
    } else {
        None
    };

    let f_u_int = |t: f64| -> Result<f64, EphemerisError> {
        let b = eval_solar_besselian(spk, time_table, t, ab)?;
        let (d, re) = b.delta_and_earth_radius();
        Ok(d - (re - b.l2.abs()))
    };
    let u2_jd_tt = if delta_max < r_earth_max - b_max.l2.abs() {
        scan_and_refine_root(f_u_int, t_max - 0.18, t_max, 8, 1.0e-10)?
    } else {
        None
    };
    let u3_jd_tt = if delta_max < r_earth_max - b_max.l2.abs() {
        scan_and_refine_root(f_u_int, t_max, t_max + 0.18, 8, 1.0e-10)?
    } else {
        None
    };

    let f_noon = |t: f64| -> Result<f64, EphemerisError> {
        let b = eval_solar_besselian(spk, time_table, t, ab)?;
        Ok(b.x)
    };
    let t_noon_jd_tt = scan_and_refine_root(f_noon, t_max - 0.25, t_max + 0.25, 8, 1.0e-10)?;

    let f_c_axis = |t: f64| -> Result<f64, EphemerisError> {
        let b = eval_solar_besselian(spk, time_table, t, ab)?;
        let (d, re) = b.delta_and_earth_radius();
        Ok(d - re)
    };
    let (c_begin_jd_tt, c_end_jd_tt) = if delta_max <= r_earth_max {
        (
            scan_and_refine_root(f_c_axis, t_max - 0.20, t_max, 8, 1.0e-10)?,
            scan_and_refine_root(f_c_axis, t_max, t_max + 0.20, 8, 1.0e-10)?,
        )
    } else {
        (None, None)
    };

    // Determine central line, magnitude, u_radius, central_duration_sec, path_width_km, and SolarEclipseKind
    let (
        kind,
        central_latitude_deg,
        central_longitude_deg,
        central_duration_sec,
        path_width_km,
        u_radius_earth,
        penumbra_radius_earth,
        magnitude,
    ) = if let Some((lat_deg, lon_deg, zeta, big_l1, big_l2)) = b_max.central_line() {
        // Central eclipse: check if Hybrid (total at max under mean or valley-floor lunar radius, annular at ends where zeta ~ 0)
        let big_l2_mean = big_l2 - (K_MOON_MEAN - K_MOON_UMBRA);
        let kind = if (big_l2 < 0.0 || big_l2_mean < 0.0) && b_max.l2 > 0.0 {
            SolarEclipseKind::Hybrid
        } else if big_l2 < 0.0 {
            SolarEclipseKind::Total
        } else {
            SolarEclipseKind::Annular
        };
        let mag = (big_l1 - big_l2) / (big_l1 + big_l2);

        // Compute topocentric shadow velocity on the central line at t_max for central_duration_sec and path_width_km
        let central_obs = Observer {
            latitude_deg: lat_deg,
            longitude_deg: lon_deg,
            altitude_m: 0.0,
        };
        let dt_v = 1.0e-4;
        let b_vm = eval_solar_besselian(spk, time_table, t_max - dt_v, ab)?;
        let b_vp = eval_solar_besselian(spk, time_table, t_max + dt_v, ab)?;
        let (xim, etam, _, _, _, _) = b_vm.local_observer(central_obs);
        let (xip, etap, _, _, _, _) = b_vp.local_observer(central_obs);
        let vx = ((b_vp.x - xip) - (b_vm.x - xim)) / (2.0 * dt_v);
        let vy = ((b_vp.y - etap) - (b_vm.y - etam)) / (2.0 * dt_v);
        let v_rel = (vx * vx + vy * vy).sqrt().max(1e-9);
        let dur_sec = (2.0 * big_l2.abs() / v_rel) * 86_400.0;

        let nx = -vy / v_rel;
        let ny = vx / v_rel;
        let proj = b_max.x * nx + b_max.y * ny;
        let stretch = (1.0 + (proj * proj) / (zeta * zeta).max(0.01)).sqrt();
        let width_km = 2.0 * big_l2.abs() * R_EARTH_EQ_KM * stretch;

        (
            kind,
            Some(lat_deg),
            Some(lon_deg),
            Some(dur_sec),
            Some(width_km),
            big_l2,
            big_l1,
            mag,
        )
    } else if delta_max < r_earth_max + b_max.l2.abs() {
        // Non-central Total or Annular (umbra grazes polar region without axis hitting Earth)
        let kind = if b_max.l2 < 0.0 {
            SolarEclipseKind::Total
        } else {
            SolarEclipseKind::Annular
        };
        let mag = (r_earth_max + b_max.l1 - delta_max) / (b_max.l1 + b_max.l2);
        (
            kind,
            None,
            None,
            None,
            None,
            b_max.l2,
            b_max.l1,
            mag,
        )
    } else {
        // Partial eclipse
        let mag = (r_earth_max + b_max.l1 - delta_max) / (b_max.l1 + b_max.l2);
        (
            SolarEclipseKind::Partial,
            None,
            None,
            None,
            None,
            b_max.l2,
            b_max.l1,
            mag,
        )
    };

    let (saros_series, saros_member) = compute_solar_saros(t_max, gamma);

    let (ut1_max, _, _, dt_info) =
        time_table.resolve_scales(TimeInput::from_tt_jd(t_max), observer);

    // Optional local topocentric circumstances
    let local = if let Some(obs) = observer {
        eval_local_solar_eclipse(spk, time_table, t_max, obs, ab)?
    } else {
        None
    };

    Ok(Some(SolarEclipseInfo {
        kind,
        saros_series,
        saros_member,
        max_eclipse_jd_tt: t_max,
        max_eclipse_jd_ut1: ut1_max.jd_f64(),
        delta_t_sec: dt_info.delta_t_sec,
        gamma,
        u_radius_earth,
        penumbra_radius_earth,
        magnitude,
        central_latitude_deg,
        central_longitude_deg,
        central_duration_sec,
        path_width_km,
        p1_jd_tt,
        u1_jd_tt,
        u2_jd_tt,
        u3_jd_tt,
        u4_jd_tt,
        p4_jd_tt,
        t_noon_jd_tt,
        c_begin_jd_tt,
        c_end_jd_tt,
        local,
    }))
}

/// Evaluates local topocentric solar eclipse circumstances for `obs` around global `t_max_global`.
fn eval_local_solar_eclipse(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    t_max_global: f64,
    obs: Observer,
    ab: AbMode,
) -> Result<Option<LocalSolarEclipse>, EphemerisError> {
    let h_deriv = 1.0e-4;
    let d_m2 = |t: f64| -> Result<f64, EphemerisError> {
        let bm = eval_solar_besselian(spk, time_table, t - h_deriv, ab)?.local_observer(obs);
        let bp = eval_solar_besselian(spk, time_table, t + h_deriv, ab)?.local_observer(obs);
        Ok(bp.3 * bp.3 - bm.3 * bm.3)
    };

    let Some(t_loc_max) = scan_and_refine_root(
        d_m2,
        t_max_global - 0.20,
        t_max_global + 0.20,
        10,
        1.0e-10,
    )?
    else {
        return Ok(None);
    };

    let b_loc = eval_solar_besselian(spk, time_table, t_loc_max, ab)?;
    let (_, _, _, m_max, big_l1, big_l2) = b_loc.local_observer(obs);
    if m_max >= big_l1 {
        return Ok(None);
    }

    let is_central_locally = m_max < big_l2.abs();
    let local_kind = if is_central_locally {
        if big_l2 < 0.0 {
            SolarEclipseKind::Total
        } else {
            SolarEclipseKind::Annular
        }
    } else {
        SolarEclipseKind::Partial
    };

    let magnitude = (big_l1 - m_max) / (big_l1 + big_l2);
    // Exact Euclidean circle-circle intersection area divided by solar disk area pi * r_sun^2
    let r_sun = 0.5 * (big_l1 + big_l2);
    let r_moon = 0.5 * (big_l1 - big_l2);
    let k_ratio = (big_l1 - big_l2) / (big_l1 + big_l2);
    let obscuration = if is_central_locally {
        if local_kind == SolarEclipseKind::Total {
            1.0
        } else {
            (k_ratio * k_ratio).clamp(0.0, 1.0)
        }
    } else if m_max >= r_sun + r_moon || r_sun <= 0.0 || r_moon <= 0.0 {
        0.0
    } else if m_max <= (r_sun - r_moon).abs() {
        (k_ratio * k_ratio).clamp(0.0, 1.0)
    } else {
        let cos_alpha = ((r_sun * r_sun + m_max * m_max - r_moon * r_moon)
            / (2.0 * r_sun * m_max))
            .clamp(-1.0, 1.0);
        let cos_beta = ((r_moon * r_moon + m_max * m_max - r_sun * r_sun)
            / (2.0 * r_moon * m_max))
            .clamp(-1.0, 1.0);
        let alpha = cos_alpha.acos();
        let beta = cos_beta.acos();
        let area_overlap = r_sun * r_sun * (alpha - alpha.sin() * cos_alpha)
            + r_moon * r_moon * (beta - beta.sin() * cos_beta);
        (area_overlap / (std::f64::consts::PI * r_sun * r_sun)).clamp(0.0, 1.0)
    };

    let f_c14 = |t: f64| -> Result<f64, EphemerisError> {
        let (_, _, _, m, l1, _) = eval_solar_besselian(spk, time_table, t, ab)?.local_observer(obs);
        Ok(m - l1)
    };
    let c1_jd_tt = scan_and_refine_root(f_c14, t_loc_max - 0.18, t_loc_max, 8, 1.0e-10)?;
    let c4_jd_tt = scan_and_refine_root(f_c14, t_loc_max, t_loc_max + 0.18, 8, 1.0e-10)?;

    let (c2_jd_tt, c3_jd_tt, central_duration_sec) = if is_central_locally {
        let f_c23 = |t: f64| -> Result<f64, EphemerisError> {
            let (_, _, _, m, _, l2) = eval_solar_besselian(spk, time_table, t, ab)?.local_observer(obs);
            Ok(m - l2.abs())
        };
        let c2 = scan_and_refine_root(f_c23, t_loc_max - 0.02, t_loc_max, 6, 1.0e-10)?;
        let c3 = scan_and_refine_root(f_c23, t_loc_max, t_loc_max + 0.02, 6, 1.0e-10)?;
        let dur = match (c2, c3) {
            (Some(t2), Some(t3)) => (t3 - t2) * 86_400.0,
            _ => 0.0,
        };
        (c2, c3, dur)
    } else {
        (None, None, 0.0)
    };

    // Sun altitude and azimuth at local maximum
    let ctx_loc = EpochContext::new(spk, time_table, TimeInput::from_tt_jd(t_loc_max), Some(obs))?;
    let sun_eq_opts = CalcOptions {
        origin: Origin::Topocentric(obs),
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };
    let r_sun_eq = eval_spk_target_vector_au(spk, &ctx_loc, 10, &sun_eq_opts)?;
    let sph_sun = crate::frame::cart_to_spherical(&r_sun_eq);
    let hz = azalt(
        obs,
        ctx_loc.gast_rad,
        sph_sun.lon_deg,
        sph_sun.lat_deg,
        AtmosphereParams::default(),
    );

    Ok(Some(LocalSolarEclipse {
        observer: obs,
        visible_at_max: hz.altitude_apparent_deg > 0.0,
        is_central_locally,
        local_kind,
        max_eclipse_jd_tt: t_loc_max,
        max_eclipse_jd_ut1: ctx_loc.ut1.jd_f64(),
        magnitude,
        obscuration,
        sun_altitude_deg: hz.altitude_true_deg,
        sun_azimuth_deg: hz.azimuth_deg,
        c1_jd_tt,
        c2_jd_tt,
        c3_jd_tt,
        c4_jd_tt,
        central_duration_sec,
    }))
}

/// Instantaneous Lunar Eclipse shadow state at `jd_tt`.
#[derive(Debug, Clone, Copy)]
struct LunarShadowInstant {
    /// Angular distance $\rho_{\text{rad}}$ between Moon's center and Earth's shadow axis in radians.
    rho_rad: f64,
    /// North-South component $y_{\text{rad}}$ (positive when Moon is North of shadow center).
    y_rad: f64,
    /// Moon's equatorial horizontal parallax $\pi_\mathbb{C}$ in radians.
    pi_moon_rad: f64,
    /// Moon's apparent semi-diameter $s_\mathbb{C}$ in radians.
    s_moon_rad: f64,
    /// Penumbral contact radius $L_1 = r_p + s_\mathbb{C}$ in radians.
    l1_rad: f64,
    /// Outer umbral contact radius $L_2 = r_u + s_\mathbb{C}$ in radians.
    l2_rad: f64,
    /// Inner (total) umbral contact radius $L_3 = r_u - s_\mathbb{C}$ in radians.
    l3_rad: f64,
}

/// Evaluates the [`LunarShadowInstant`] at `jd_tt` using geocentric Sun and Moon vectors with `aberration`.
fn eval_lunar_shadow(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    jd_tt: f64,
    shadow_rule: ShadowEnlargementRule,
    aberration: AbMode,
) -> Result<LunarShadowInstant, EphemerisError> {
    let ctx = EpochContext::new(spk, time_table, TimeInput::from_tt_jd(jd_tt), None)?;
    let opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };
    let r_sun_au = eval_spk_target_vector_au(spk, &ctx, 10, &opts)?;
    let r_moon_au = eval_spk_target_vector_au(spk, &ctx, 301, &opts)?;

    let r_sun_km = vec_scale(&r_sun_au, AU_KM);
    let r_moon_km = vec_scale(&r_moon_au, AU_KM);
    let dist_sun_km = vec_norm(&r_sun_km);
    let dist_moon_km = vec_norm(&r_moon_km);

    // Shadow axis unit vector k_hat = -r_sun / |r_sun|
    let k_hat = vec_scale(&r_sun_km, -1.0 / dist_sun_km);
    let sin_d = k_hat[2].clamp(-1.0, 1.0);
    let cos_d = (1.0 - sin_d * sin_d).max(1e-15).sqrt();
    let a_rad = k_hat[1].atan2(k_hat[0]);
    let (sin_a, cos_a) = a_rad.sin_cos();

    let _i_hat: Vec3 = [-sin_a, cos_a, 0.0];
    let j_hat: Vec3 = [-cos_a * sin_d, -sin_a * sin_d, cos_d];

    let u_moon = vec_scale(&r_moon_km, 1.0 / dist_moon_km);
    let y_rad = vec_dot(&u_moon, &j_hat);
    let cos_rho = vec_dot(&u_moon, &k_hat).clamp(-1.0, 1.0);
    let rho_rad = cos_rho.acos();

    let pi_moon_rad = (R_EARTH_EQ_KM / dist_moon_km).clamp(-1.0, 1.0).asin();
    let pi_sun_rad = (R_EARTH_EQ_KM / dist_sun_km).clamp(-1.0, 1.0).asin();
    let s_sun_rad = (R_SUN_KM / dist_sun_km).clamp(-1.0, 1.0).asin();
    let s_moon_rad = ((K_MOON_MEAN * R_EARTH_EQ_KM) / dist_moon_km)
        .clamp(-1.0, 1.0)
        .asin();

    // Mean Earth oblateness factor at lunar latitude (0.99833 * pi_moon)
    let pi_1 = 0.998_33 * pi_moon_rad;

    let (r_umbra, r_penumbra) = match shadow_rule {
        ShadowEnlargementRule::Danjon1Over85 => {
            let enlarged_pi = (86.0 / 85.0) * pi_1;
            (
                enlarged_pi + pi_sun_rad - s_sun_rad,
                enlarged_pi + pi_sun_rad + s_sun_rad,
            )
        }
        ShadowEnlargementRule::Traditional2Percent => (
            1.02 * (pi_1 + pi_sun_rad - s_sun_rad),
            1.02 * (pi_1 + pi_sun_rad + s_sun_rad),
        ),
    };

    Ok(LunarShadowInstant {
        rho_rad,
        y_rad,
        pi_moon_rad,
        s_moon_rad,
        l1_rad: r_penumbra + s_moon_rad,
        l2_rad: r_umbra + s_moon_rad,
        l3_rad: r_umbra - s_moon_rad,
    })
}

/// Evaluates whether a Lunar Eclipse occurs near `approx_syzygy_jd_tt`, returning `Ok(Some(LunarEclipseInfo))`
/// if an eclipse takes place, or `Ok(None)` if the Full Moon misses the Earth's penumbra.
pub fn lunar_eclipse_at(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    approx_syzygy_jd_tt: f64,
    options: EclipseOptions,
) -> Result<Option<LunarEclipseInfo>, EphemerisError> {
    let ab = options.aberration;
    // 1. Parabolic 3-point pre-estimate of the minimum of rho^2(t)
    let h_pre = 0.10;
    let s_m = eval_lunar_shadow(spk, time_table, approx_syzygy_jd_tt - h_pre, options.shadow_rule, ab)?;
    let s_0 = eval_lunar_shadow(spk, time_table, approx_syzygy_jd_tt, options.shadow_rule, ab)?;
    let s_p = eval_lunar_shadow(spk, time_table, approx_syzygy_jd_tt + h_pre, options.shadow_rule, ab)?;
    let fm = s_m.rho_rad * s_m.rho_rad;
    let f0 = s_0.rho_rad * s_0.rho_rad;
    let fp = s_p.rho_rad * s_p.rho_rad;
    let denom = fp - 2.0 * f0 + fm;
    let dt_est = if denom > 1e-15 {
        (-0.5 * h_pre * (fp - fm) / denom).clamp(-0.55, 0.55)
    } else {
        0.0
    };
    let t_est = approx_syzygy_jd_tt + dt_est;
    let s_est = eval_lunar_shadow(spk, time_table, t_est, options.shadow_rule, ab)?;
    if s_est.rho_rad > s_est.l1_rad + 0.005 {
        return Ok(None);
    }

    let h_deriv = 1.0e-4;
    let d_rho2 = |t: f64| -> Result<f64, EphemerisError> {
        let sm = eval_lunar_shadow(spk, time_table, t - h_deriv, options.shadow_rule, ab)?;
        let sp = eval_lunar_shadow(spk, time_table, t + h_deriv, options.shadow_rule, ab)?;
        Ok(sp.rho_rad * sp.rho_rad - sm.rho_rad * sm.rho_rad)
    };

    let Some(t_max) = scan_and_refine_root(
        d_rho2,
        t_est - 0.03,
        t_est + 0.03,
        4,
        1.0e-10,
    )?
    else {
        return Ok(None);
    };

    let s_max = eval_lunar_shadow(spk, time_table, t_max, options.shadow_rule, ab)?;
    if s_max.rho_rad >= s_max.l1_rad {
        return Ok(None);
    }

    let gamma_abs = s_max.rho_rad / s_max.pi_moon_rad;
    let gamma = if s_max.y_rad >= 0.0 {
        gamma_abs
    } else {
        -gamma_abs
    };

    let penumbral_magnitude = (s_max.l1_rad - s_max.rho_rad) / (2.0 * s_max.s_moon_rad);
    let umbral_magnitude = (s_max.l2_rad - s_max.rho_rad) / (2.0 * s_max.s_moon_rad);

    let kind = if s_max.rho_rad <= s_max.l3_rad {
        LunarEclipseKind::Total
    } else if s_max.rho_rad <= s_max.l2_rad {
        LunarEclipseKind::Partial
    } else {
        LunarEclipseKind::Penumbral
    };

    let f_p = |t: f64| -> Result<f64, EphemerisError> {
        let s = eval_lunar_shadow(spk, time_table, t, options.shadow_rule, ab)?;
        Ok(s.rho_rad - s.l1_rad)
    };
    let p1_jd_tt = scan_and_refine_root(f_p, t_max - 0.30, t_max, 8, 1.0e-10)?;
    let p4_jd_tt = scan_and_refine_root(f_p, t_max, t_max + 0.30, 8, 1.0e-10)?;

    let (u1_jd_tt, u4_jd_tt) = if kind != LunarEclipseKind::Penumbral {
        let f_u = |t: f64| -> Result<f64, EphemerisError> {
            let s = eval_lunar_shadow(spk, time_table, t, options.shadow_rule, ab)?;
            Ok(s.rho_rad - s.l2_rad)
        };
        (
            scan_and_refine_root(f_u, t_max - 0.20, t_max, 8, 1.0e-10)?,
            scan_and_refine_root(f_u, t_max, t_max + 0.20, 8, 1.0e-10)?,
        )
    } else {
        (None, None)
    };

    let (u2_jd_tt, u3_jd_tt) = if kind == LunarEclipseKind::Total {
        let f_t = |t: f64| -> Result<f64, EphemerisError> {
            let s = eval_lunar_shadow(spk, time_table, t, options.shadow_rule, ab)?;
            Ok(s.rho_rad - s.l3_rad)
        };
        (
            scan_and_refine_root(f_t, t_max - 0.10, t_max, 8, 1.0e-10)?,
            scan_and_refine_root(f_t, t_max, t_max + 0.10, 8, 1.0e-10)?,
        )
    } else {
        (None, None)
    };

    let semi_dur_min = |t1: Option<f64>, t2: Option<f64>| -> f64 {
        match (t1, t2) {
            (Some(a), Some(b)) => 0.5 * (b - a) * 1_440.0,
            _ => 0.0,
        }
    };

    let (saros_series, saros_member) = compute_lunar_saros(t_max, gamma);
    let (ut1_max, _, _, dt_info) = time_table.resolve_scales(TimeInput::from_tt_jd(t_max), None);

    Ok(Some(LunarEclipseInfo {
        kind,
        shadow_rule: options.shadow_rule,
        saros_series,
        saros_member,
        max_eclipse_jd_tt: t_max,
        max_eclipse_jd_ut1: ut1_max.jd_f64(),
        delta_t_sec: dt_info.delta_t_sec,
        gamma,
        penumbral_magnitude,
        umbral_magnitude,
        penumbral_semi_duration_min: semi_dur_min(p1_jd_tt, p4_jd_tt),
        partial_semi_duration_min: semi_dur_min(u1_jd_tt, u4_jd_tt),
        total_semi_duration_min: semi_dur_min(u2_jd_tt, u3_jd_tt),
        p1_jd_tt,
        u1_jd_tt,
        u2_jd_tt,
        u3_jd_tt,
        u4_jd_tt,
        p4_jd_tt,
    }))
}

/// Evaluates the wrapped ecliptic elongation $\Delta\lambda(t) = (\lambda_\mathbb{C}(t) - \lambda_\odot(t) - \text{target\_deg})$ in `[-180, +180)`.
///
/// Uses direct J2000 ecliptic projection from SPK states to avoid expensive IAU 2000A nutation during coarse syzygy bracketing.
fn elongation_diff_deg(
    spk: &dyn SpkSource,
    _time_table: &TimeTable,
    jd_tt: f64,
    target_deg: f64,
) -> Result<f64, EphemerisError> {
    let st_earth = spk.state_icrf_km_split(399, 0, jd_tt, 0.0)?;
    let st_sun = spk.state_icrf_km_split(10, 0, jd_tt, 0.0)?;
    let st_moon = spk.state_icrf_km_split(301, 0, jd_tt, 0.0)?;
    let r_sun = vec_sub(&st_sun.pos_km, &st_earth.pos_km);
    let r_moon = vec_sub(&st_moon.pos_km, &st_earth.pos_km);

    // Rotate from ICRF equatorial to J2000 ecliptic (eps0 = 23.43929111 deg)
    let (sin_e, cos_e) = 23.439_291_11_f64.to_radians().sin_cos();
    let sun_lon = (r_sun[1] * cos_e + r_sun[2] * sin_e)
        .atan2(r_sun[0])
        .to_degrees();
    let moon_lon = (r_moon[1] * cos_e + r_moon[2] * sin_e)
        .atan2(r_moon[0])
        .to_degrees();
    Ok((moon_lon - sun_lon - target_deg + 180.0).rem_euclid(360.0) - 180.0)
}

/// Searches forward (`backward == false`) or backward (`backward == true`) from `start_jd_tt`
/// for the next Solar Eclipse, stepping through every New Moon without skipping grazing events.
pub fn next_solar_eclipse(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    start_jd_tt: f64,
    backward: bool,
    observer: Option<Observer>,
    options: EclipseOptions,
) -> Result<SolarEclipseInfo, EphemerisError> {
    if !start_jd_tt.is_finite()
        || !(1..=240_000).contains(&options.max_synodic_months)
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid next_solar_eclipse input: start_jd_tt={start_jd_tt}, max_synodic_months={}",
                options.max_synodic_months
            ),
        });
    }
    let step_days = if backward { -5.0 } else { 5.0 };
    let max_steps = (options.max_synodic_months as usize) * 6;
    // Start 1 day before start_jd_tt (in search direction) so an eclipse whose mean syzygy
    // slightly precedes start_jd_tt while max_eclipse_jd_tt >= start_jd_tt is never missed.
    let mut t_prev = if backward {
        start_jd_tt + 1.0
    } else {
        start_jd_tt - 1.0
    };
    let mut e_prev = elongation_diff_deg(spk, time_table, t_prev, 0.0)?;

    for _ in 0..max_steps {
        let t_curr = t_prev + step_days;
        let e_curr = elongation_diff_deg(spk, time_table, t_curr, 0.0)?;

        let crossed = if backward {
            e_prev >= 0.0 && e_curr < 0.0 && (e_prev - e_curr) < 120.0
        } else {
            e_prev <= 0.0 && e_curr > 0.0 && (e_curr - e_prev) < 120.0
        };

        if crossed {
            let t_syzygy = illinois_root(
                |t| elongation_diff_deg(spk, time_table, t, 0.0),
                t_prev.min(t_curr),
                t_prev.max(t_curr),
                1.0e-5,
                40,
            )?;
            if let Some(info) = solar_eclipse_at(spk, time_table, t_syzygy, observer, options)? {
                if (!backward && info.max_eclipse_jd_tt >= start_jd_tt - 1e-4)
                    || (backward && info.max_eclipse_jd_tt <= start_jd_tt + 1e-4)
                {
                    return Ok(info);
                }
            }
        }

        t_prev = t_curr;
        e_prev = e_curr;
    }

    Err(EphemerisError::NonConvergence {
        solver: "next_solar_eclipse".to_owned(),
        iterations: max_steps as u32,
        residual: 0.0,
    })
}

/// Searches forward (`backward == false`) or backward (`backward == true`) from `start_jd_tt`
/// for the next Lunar Eclipse, stepping through every Full Moon without skipping penumbral events.
pub fn next_lunar_eclipse(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    start_jd_tt: f64,
    backward: bool,
    options: EclipseOptions,
) -> Result<LunarEclipseInfo, EphemerisError> {
    if !start_jd_tt.is_finite()
        || !(1..=240_000).contains(&options.max_synodic_months)
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid next_lunar_eclipse input: start_jd_tt={start_jd_tt}, max_synodic_months={}",
                options.max_synodic_months
            ),
        });
    }
    let step_days = if backward { -5.0 } else { 5.0 };
    let max_steps = (options.max_synodic_months as usize) * 6;
    let mut t_prev = if backward {
        start_jd_tt + 1.0
    } else {
        start_jd_tt - 1.0
    };
    let mut e_prev = elongation_diff_deg(spk, time_table, t_prev, 180.0)?;

    for _ in 0..max_steps {
        let t_curr = t_prev + step_days;
        let e_curr = elongation_diff_deg(spk, time_table, t_curr, 180.0)?;

        let crossed = if backward {
            e_prev >= 0.0 && e_curr < 0.0 && (e_prev - e_curr) < 120.0
        } else {
            e_prev <= 0.0 && e_curr > 0.0 && (e_curr - e_prev) < 120.0
        };

        if crossed {
            let t_syzygy = illinois_root(
                |t| elongation_diff_deg(spk, time_table, t, 180.0),
                t_prev.min(t_curr),
                t_prev.max(t_curr),
                1.0e-5,
                40,
            )?;
            if let Some(info) = lunar_eclipse_at(spk, time_table, t_syzygy, options)? {
                if (!backward && info.max_eclipse_jd_tt >= start_jd_tt - 1e-4)
                    || (backward && info.max_eclipse_jd_tt <= start_jd_tt + 1e-4)
                {
                    return Ok(info);
                }
            }
        }

        t_prev = t_curr;
        e_prev = e_curr;
    }

    Err(EphemerisError::NonConvergence {
        solver: "next_lunar_eclipse".to_owned(),
        iterations: max_steps as u32,
        residual: 0.0,
    })
}

/// Searches forward (`backward == false`) or backward (`backward == true`) from `start_jd_tt`
/// for the next Solar Eclipse that is **locally visible above the horizon** at `observer`
/// (`local.visible_at_max && local.magnitude > 0.0`, `R8.2`).
pub fn next_local_solar_eclipse(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    start_jd_tt: f64,
    observer: Observer,
    backward: bool,
    options: EclipseOptions,
) -> Result<SolarEclipseInfo, EphemerisError> {
    if !start_jd_tt.is_finite()
        || !(1..=240_000).contains(&options.max_synodic_months)
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid next_local_solar_eclipse input: start_jd_tt={start_jd_tt}, max_synodic_months={}",
                options.max_synodic_months
            ),
        });
    }
    let step_days = if backward { -5.0 } else { 5.0 };
    let search_months = options.max_synodic_months.max(1200) as usize;
    let max_steps = search_months * 6;
    let mut t_prev = if backward {
        start_jd_tt + 1.0
    } else {
        start_jd_tt - 1.0
    };
    let mut e_prev = elongation_diff_deg(spk, time_table, t_prev, 0.0)?;

    let sun_visible_at = |jd_tt: f64| -> Result<bool, EphemerisError> {
        let ctx = EpochContext::new(spk, time_table, TimeInput::from_tt_jd(jd_tt), Some(observer))?;
        let sun_eq_opts = CalcOptions {
            origin: Origin::Topocentric(observer),
            plane: Plane::Equatorial,
            equinox: Equinox::OfDate,
            nutation: NutationMode::True,
            aberration: AbMode::Apparent,
            zodiac: Zodiac::Tropical,
            include_speed: false,
        };
        let r_sun_eq = eval_spk_target_vector_au(spk, &ctx, 10, &sun_eq_opts)?;
        let sph_sun = crate::frame::cart_to_spherical(&r_sun_eq);
        let hz = azalt(
            observer,
            ctx.gast_rad,
            sph_sun.lon_deg,
            sph_sun.lat_deg,
            AtmosphereParams::default(),
        );
        Ok(hz.altitude_apparent_deg > 0.0)
    };

    for _ in 0..max_steps {
        let t_curr = t_prev + step_days;
        let e_curr = elongation_diff_deg(spk, time_table, t_curr, 0.0)?;

        let crossed = if backward {
            e_prev >= 0.0 && e_curr < 0.0 && (e_prev - e_curr) < 120.0
        } else {
            e_prev <= 0.0 && e_curr > 0.0 && (e_curr - e_prev) < 120.0
        };

        if crossed {
            let t_syzygy = illinois_root(
                |t| elongation_diff_deg(spk, time_table, t, 0.0),
                t_prev.min(t_curr),
                t_prev.max(t_curr),
                1.0e-5,
                40,
            )?;
            if let Some(info) =
                solar_eclipse_at(spk, time_table, t_syzygy, Some(observer), options)?
            {
                let in_dir = (!backward && info.max_eclipse_jd_tt >= start_jd_tt - 1e-4)
                    || (backward && info.max_eclipse_jd_tt <= start_jd_tt + 1e-4);
                if in_dir {
                    if let Some(ref loc) = info.local {
                        if loc.magnitude > 0.0 {
                            let visible_any = loc.visible_at_max
                                || match loc.c1_jd_tt {
                                    Some(t_c1) => sun_visible_at(t_c1)?,
                                    None => false,
                                }
                                || match loc.c4_jd_tt {
                                    Some(t_c4) => sun_visible_at(t_c4)?,
                                    None => false,
                                };
                            if visible_any {
                                return Ok(info);
                            }
                        }
                    }
                }
            }
        }

        t_prev = t_curr;
        e_prev = e_curr;
    }

    Err(EphemerisError::NonConvergence {
        solver: "next_local_solar_eclipse".to_owned(),
        iterations: max_steps as u32,
        residual: 0.0,
    })
}

/// Finds all Solar Eclipses occurring within `[start_jd_tt, end_jd_tt]` (`R8.3`).
pub fn solar_eclipses_in_range(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    start_jd_tt: f64,
    end_jd_tt: f64,
    observer: Option<Observer>,
    options: EclipseOptions,
) -> Result<Vec<SolarEclipseInfo>, EphemerisError> {
    if !start_jd_tt.is_finite()
        || !end_jd_tt.is_finite()
        || end_jd_tt < start_jd_tt
        || (end_jd_tt - start_jd_tt) > 3_652_500.0
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid eclipse range: start_jd_tt={start_jd_tt}, end_jd_tt={end_jd_tt}"
            ),
        });
    }
    let mut results = Vec::new();
    let step_days = 5.0;
    let total_steps = (((end_jd_tt - start_jd_tt + 10.0) / step_days).ceil() as usize).max(1);
    let mut t_prev = start_jd_tt - 5.0;
    let mut e_prev = elongation_diff_deg(spk, time_table, t_prev, 0.0)?;

    for _ in 0..=total_steps {
        let t_curr = t_prev + step_days;
        let e_curr = elongation_diff_deg(spk, time_table, t_curr, 0.0)?;
        if e_prev <= 0.0 && e_curr > 0.0 && (e_curr - e_prev) < 120.0 {
            let t_syzygy = illinois_root(
                |t| elongation_diff_deg(spk, time_table, t, 0.0),
                t_prev,
                t_curr,
                1.0e-5,
                40,
            )?;
            if let Some(info) = solar_eclipse_at(spk, time_table, t_syzygy, observer, options)? {
                if info.max_eclipse_jd_tt >= start_jd_tt && info.max_eclipse_jd_tt <= end_jd_tt {
                    results.push(info);
                }
            }
        }
        t_prev = t_curr;
        e_prev = e_curr;
    }

    Ok(results)
}

/// Finds all Lunar Eclipses occurring within `[start_jd_tt, end_jd_tt]` (`R8.3`).
pub fn lunar_eclipses_in_range(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    start_jd_tt: f64,
    end_jd_tt: f64,
    options: EclipseOptions,
) -> Result<Vec<LunarEclipseInfo>, EphemerisError> {
    if !start_jd_tt.is_finite()
        || !end_jd_tt.is_finite()
        || end_jd_tt < start_jd_tt
        || (end_jd_tt - start_jd_tt) > 3_652_500.0
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid eclipse range: start_jd_tt={start_jd_tt}, end_jd_tt={end_jd_tt}"
            ),
        });
    }
    let mut results = Vec::new();
    let step_days = 5.0;
    let total_steps = (((end_jd_tt - start_jd_tt + 10.0) / step_days).ceil() as usize).max(1);
    let mut t_prev = start_jd_tt - 5.0;
    let mut e_prev = elongation_diff_deg(spk, time_table, t_prev, 180.0)?;

    for _ in 0..=total_steps {
        let t_curr = t_prev + step_days;
        let e_curr = elongation_diff_deg(spk, time_table, t_curr, 180.0)?;
        if e_prev <= 0.0 && e_curr > 0.0 && (e_curr - e_prev) < 120.0 {
            let t_syzygy = illinois_root(
                |t| elongation_diff_deg(spk, time_table, t, 180.0),
                t_prev,
                t_curr,
                1.0e-5,
                40,
            )?;
            if let Some(info) = lunar_eclipse_at(spk, time_table, t_syzygy, options)? {
                if info.max_eclipse_jd_tt >= start_jd_tt && info.max_eclipse_jd_tt <= end_jd_tt {
                    results.push(info);
                }
            }
        }
        t_prev = t_curr;
        e_prev = e_curr;
    }

    Ok(results)
}

/// Evaluates the instantaneous geographic position and attributes of a Solar Eclipse at `jd_ut1` (`swe_sol_eclipse_where`),
/// returning `Ok(None)` if no solar eclipse is active on Earth at `jd_ut1`.
#[allow(clippy::type_complexity)]
pub fn solar_eclipse_where_at_ut1(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    jd_ut1: f64,
    options: EclipseOptions,
) -> Result<Option<(SolarEclipseKind, bool, [f64; 10], [f64; 20])>, EphemerisError> {
    if !jd_ut1.is_finite() {
        return Err(EphemerisError::InvalidInput {
            message: format!("Non-finite jd_ut1={jd_ut1} in solar_eclipse_where_at_ut1"),
        });
    }
    let (_, tt, _, _) = time_table.resolve_scales(TimeInput::from_ut1_jd(jd_ut1), None);
    let jd_tt = tt.jd_f64();
    let b = eval_solar_besselian(spk, time_table, jd_tt, options.aberration)?;
    let (delta, r_earth) = b.delta_and_earth_radius();
    if delta >= r_earth + b.l1 {
        return Ok(None);
    }

    let mut geopos = [0.0; 10];
    let mut attr = [0.0; 20];
    let (kind, is_central, mag, ratio, u_diam_km) =
        if let Some((lat_deg, lon_deg, _zeta, big_l1, big_l2)) = b.central_line() {
            let big_l2_mean = big_l2 - (K_MOON_MEAN - K_MOON_UMBRA);
            let kind = if (big_l2 < 0.0 || big_l2_mean < 0.0) && b.l2 > 0.0 {
                SolarEclipseKind::Hybrid
            } else if big_l2 < 0.0 {
                SolarEclipseKind::Total
            } else {
                SolarEclipseKind::Annular
            };
            geopos[0] = lon_deg;
            geopos[1] = lat_deg;
            let mag = (big_l1 - big_l2) / (big_l1 + big_l2);
            let ratio = (big_l1 - big_l2) / (big_l1 + big_l2);
            let u_diam_km = 2.0 * big_l2.abs() * R_EARTH_EQ_KM;
            (kind, true, mag, ratio, u_diam_km)
        } else {
            let y1 = b.y / b.rho1;
            let g_scaled = (b.x * b.x + y1 * y1).sqrt().max(1e-12);
            let ux = b.x / g_scaled;
            let uy1 = y1 / g_scaled;
            let sin_d1 = b.sin_d / b.rho1;
            let cos_d1 = ((1.0 - WGS84_E2).sqrt() * b.cos_d) / b.rho1;
            let sin_phi1 = uy1 * cos_d1;
            let cos_phi1_sin_h = ux;
            let cos_phi1_cos_h = -uy1 * sin_d1;
            let h_rad = cos_phi1_sin_h.atan2(cos_phi1_cos_h);
            let cos_phi1 = (cos_phi1_sin_h * cos_phi1_sin_h + cos_phi1_cos_h * cos_phi1_cos_h)
                .sqrt()
                .max(1e-15);
            let tan_phi_geod = (sin_phi1 / cos_phi1) / (1.0 - WGS84_E2).sqrt();
            let lat_deg = tan_phi_geod.atan().to_degrees();
            let lon_east_deg = (h_rad - b.mu_rad).to_degrees();
            let lon_deg = (lon_east_deg + 180.0).rem_euclid(360.0) - 180.0;
            geopos[0] = lon_deg;
            geopos[1] = lat_deg;
            let kind = if delta < r_earth + b.l2.abs() {
                if b.l2 < 0.0 {
                    SolarEclipseKind::Total
                } else {
                    SolarEclipseKind::Annular
                }
            } else {
                SolarEclipseKind::Partial
            };
            let mag = (r_earth + b.l1 - delta) / (b.l1 + b.l2);
            let ratio = (b.l1 - b.l2) / (b.l1 + b.l2);
            let u_diam_km = 2.0 * b.l2.abs() * R_EARTH_EQ_KM;
            (kind, false, mag, ratio, u_diam_km)
        };

    let (saros_series, saros_member) =
        compute_solar_saros(jd_tt, if b.y >= 0.0 { delta } else { -delta });
    attr[0] = mag;
    attr[1] = ratio;
    attr[2] = if kind == SolarEclipseKind::Total {
        1.0
    } else {
        (ratio * ratio).clamp(0.0, 1.0) * mag.clamp(0.0, 1.0)
    };
    attr[3] = u_diam_km;
    attr[8] = mag;
    attr[9] = f64::from(saros_series);
    attr[10] = f64::from(saros_member);

    Ok(Some((kind, is_central, geopos, attr)))
}

/// Evaluates the instantaneous attributes of a Lunar Eclipse at `jd_ut1` (`swe_lun_eclipse_how`),
/// returning `Ok(None)` if the Moon is outside the Earth's penumbra at `jd_ut1`.
pub fn lunar_eclipse_how_at_ut1(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    jd_ut1: f64,
    options: EclipseOptions,
) -> Result<Option<(LunarEclipseKind, [f64; 20])>, EphemerisError> {
    if !jd_ut1.is_finite() {
        return Err(EphemerisError::InvalidInput {
            message: format!("Non-finite jd_ut1={jd_ut1} in lunar_eclipse_how_at_ut1"),
        });
    }
    let (_, tt, _, _) = time_table.resolve_scales(TimeInput::from_ut1_jd(jd_ut1), None);
    let jd_tt = tt.jd_f64();
    let s = eval_lunar_shadow(spk, time_table, jd_tt, options.shadow_rule, options.aberration)?;
    if s.rho_rad >= s.l1_rad {
        return Ok(None);
    }

    let umbral_mag = (s.l2_rad - s.rho_rad) / (2.0 * s.s_moon_rad);
    let penumbral_mag = (s.l1_rad - s.rho_rad) / (2.0 * s.s_moon_rad);
    let kind = if s.rho_rad <= s.l3_rad {
        LunarEclipseKind::Total
    } else if s.rho_rad < s.l2_rad {
        LunarEclipseKind::Partial
    } else {
        LunarEclipseKind::Penumbral
    };

    let mut attr = [0.0; 20];
    attr[0] = umbral_mag.max(0.0);
    attr[1] = penumbral_mag;
    attr[7] = s.rho_rad.to_degrees();
    attr[8] = umbral_mag.max(0.0);
    if let Some(info) = lunar_eclipse_at(spk, time_table, jd_tt, options)? {
        attr[9] = f64::from(info.saros_series);
        attr[10] = f64::from(info.saros_member);
    }

    Ok(Some((kind, attr)))
}
