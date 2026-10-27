//! All 23 case-sensitive Astrological House Systems (`P, K, O, R, C, A, W, B, M, T, U, X, V, L, Q, I, i, F, N, S, D, Y, H`),
//! `houses_armc`, 8 chart angles (`ASC`, `MC`, `ARMC`, `Vertex`, `Equatorial ASC`, `Co-ASC Koch`,
//! `Co-ASC Munkasey`, `Polar ASC`), cusp & angle speeds (`deg/day`), 36 Gauquelin sectors,
//! and strict polar policy enforcement (`PolarPolicy::Strict` default vs `PolarPolicy::FallbackToPorphyry`).
//!
//! Clean-room mathematical implementation based on:
//! - Holden, R. W. (1977), *The Elements of House Division*
//! - Munkasey, M. P. (1991), *An Astrological House Formulary*
//! - Permissive references: `jayeshmepani/jpl-ephemeris` (MIT) and `TheDaniel166/moira` (MIT).

use crate::ayanamsa::eval_ayanamsa;
use crate::error::EphemerisError;
use crate::frame::FrameTransform;
use crate::moment::Moment;
use crate::spk::SpkSource;
use crate::time::{TimeInput, TimeTable};
use crate::types::{
    AbMode, CalcOptions, Equinox, HouseSystem, HouseWarning, NutationMode, Observer, Origin, Plane,
    PlanetId, PolarCategory, PolarPolicy, Zodiac,
};
use serde::{Deserialize, Serialize};
use sofars::erst::gst06;

/// Options for house cusp and angle computation (`Engine::houses` and `houses_armc`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct HouseOptions {
    /// Requested house system (`P, K, O, R, C, A, W, B, M, T, U, X, V, L, Q, I, i, F, N, S, D, Y, H`).
    pub system: HouseSystem,
    /// Behaviour inside the polar circle ($|\phi| \ge 90^\circ - \varepsilon$); defaults to `PolarPolicy::Strict`.
    pub polar_policy: PolarPolicy,
    /// Tropical or Sidereal (`AyanamsaModel`) zodiac.
    pub zodiac: Zodiac,
    /// True (`ε0 + Δε`) or Mean (`ε0`) obliquity of the ecliptic.
    pub nutation: NutationMode,
    /// Whether to compute angular velocities (`deg/day`) for cusps and angles.
    pub include_speed: bool,
}

impl Default for HouseOptions {
    fn default() -> Self {
        Self {
            system: HouseSystem::Placidus,
            polar_policy: PolarPolicy::Strict,
            zodiac: Zodiac::Tropical,
            nutation: NutationMode::True,
            include_speed: true,
        }
    }
}

/// The 8 canonical astrological chart angles in degrees `[0, 360)`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HouseAngles {
    /// Ascendant (`ASC`, Eastern horizon intersection with the ecliptic).
    pub asc_deg: f64,
    /// Midheaven (`MC`, Upper meridian intersection with the ecliptic).
    pub mc_deg: f64,
    /// Right Ascension of the Midheaven (`RAMC` / `ARMC`, local apparent sidereal time in degrees `[0, 360)`).
    pub armc_deg: f64,
    /// Vertex (Western Prime Vertical intersection with the ecliptic).
    pub vertex_deg: f64,
    /// Equatorial Ascendant (`East Point`, intersection of celestial equator East point with ecliptic).
    pub equatorial_asc_deg: f64,
    /// Co-Ascendant according to Walter Koch.
    pub co_asc_koch_deg: f64,
    /// Co-Ascendant according to Michael Munkasey.
    pub co_asc_munkasey_deg: f64,
    /// Polar Ascendant according to Michael Munkasey.
    pub polar_asc_deg: f64,
}

/// Complete house division result with provenance, polar fallback tracking, warnings, and speeds.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HousesResult {
    /// The house system originally requested by the caller.
    pub requested_house_system: HouseSystem,
    /// The house system actually evaluated (`Porphyry` or `WholeSign` if polar fallback triggered).
    pub effective_house_system: HouseSystem,
    /// `true` when `effective_house_system != requested_house_system` due to polar or cusp-inversion fallback (`R4.4`).
    pub polar_house_fallback: bool,
    /// Optional human-readable explanation when `effective_house_system != requested_house_system`.
    pub fallback_reason: Option<String>,
    /// Polar or geometric warnings (`PolarFallbackToPorphyry`, `PolarFallbackToWholeSign`, `PolarUnreliable`, `NearPoleBoundary`).
    pub warnings: Vec<HouseWarning>,
    /// 12 house cusps in degrees `[0, 360)` (index `0` = House 1, ..., index `11` = House 12).
    pub cusps_deg: [f64; 12],
    /// Angular velocities of the 12 house cusps in degrees per day (`0.0` if `include_speed == false`).
    pub cusp_speeds_deg_day: [f64; 12],
    /// The 8 chart angles in degrees `[0, 360)`.
    pub angles: HouseAngles,
    /// Angular velocities of the 8 chart angles in degrees per day.
    pub angle_speeds: HouseAngles,
}

/// Returns `true` if the 12 house cusps form a strictly counter-clockwise zodiacal ring
/// ($0^\circ < (c_{(i+1)\bmod 12} - c_i) \bmod 360^\circ < 180^\circ$ and $\sum_{i=0}^{11} \Delta_i = 360^\circ$, `R4.4`).
#[must_use]
pub fn is_cusp_ring_monotonic(cusps_deg: &[f64; 12]) -> bool {
    let mut total_arc = 0.0;
    for i in 0..12 {
        let c_curr = cusps_deg.get(i).copied().unwrap_or(0.0);
        let c_next = cusps_deg.get((i + 1) % 12).copied().unwrap_or(0.0);
        let step = (c_next - c_curr).rem_euclid(360.0);
        if !(1e-6..180.0).contains(&step) {
            return false;
        }
        total_arc += step;
    }
    (total_arc - 360.0).abs() < 1e-4
}

/// Determines which house (`1..=12`) contains ecliptic longitude `lon_deg` using half-open intervals
/// $[c_i - \varepsilon, c_{i+1} - \varepsilon)$ with boundary tolerance `eps_deg` (default `1e-9°`, `R4.5`),
/// ensuring that a body within `eps_deg` of cusp $i+1$ is placed into house $i+1$ rather than falling
/// into house $i$ due to floating-point rounding.
#[must_use]
pub fn house_of_lon(lon_deg: f64, cusps_deg: &[f64; 12], eps_deg: f64) -> u8 {
    if !lon_deg.is_finite() {
        return 1;
    }
    let lon = lon_deg.rem_euclid(360.0);
    let tol = if eps_deg.is_finite() {
        eps_deg.abs()
    } else {
        1.0e-9
    };

    // 1. Check if `lon` is within `tol` of any cusp `c_i`; if so, assign to House `i + 1`.
    for (i, &c) in cusps_deg.iter().enumerate() {
        let diff = (lon - c.rem_euclid(360.0) + 180.0).rem_euclid(360.0) - 180.0;
        if diff.abs() <= tol {
            return (i + 1) as u8;
        }
    }

    // 2. Otherwise, check which half-open arc [c_i, c_{(i+1)%12}) contains `lon`.
    for i in 0..12 {
        let c_start = cusps_deg.get(i).copied().unwrap_or(0.0).rem_euclid(360.0);
        let c_end = cusps_deg
            .get((i + 1) % 12)
            .copied()
            .unwrap_or(0.0)
            .rem_euclid(360.0);
        let span = (c_end - c_start).rem_euclid(360.0);
        let offset = (lon - c_start).rem_euclid(360.0);
        if offset < span {
            return (i + 1) as u8;
        }
    }

    1
}

mod serde_array36 {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    pub fn serialize<S>(arr: &[f64; 36], serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        arr.as_slice().serialize(serializer)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<[f64; 36], D::Error>
    where
        D: Deserializer<'de>,
    {
        let v: Vec<f64> = Vec::deserialize(deserializer)?;
        v.try_into()
            .map_err(|_| serde::de::Error::custom("expected exactly 36 elements"))
    }
}

/// 36 Gauquelin sectors (`1..=36`, clockwise from Ascendant to MC in Placidus semi-arc space).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GauquelinResult {
    /// 36 sector cusps in degrees `[0, 360)`.
    #[serde(with = "serde_array36")]
    pub sectors_deg: [f64; 36],
    /// Warnings emitted during semi-arc evaluation.
    pub warnings: Vec<HouseWarning>,
}

/// Fundamental ecliptic projection `asc1(x_deg, f_lat_deg, eps_rad)` (Holden 1977 / Munkasey 1991):
///
/// $$\lambda = \operatorname{atan2}\bigl(\sin x,\; \cos x \cos\varepsilon - \tan f \sin\varepsilon\bigr) \bmod 360^\circ$$
#[must_use]
pub fn asc1(x_deg: f64, f_lat_deg: f64, eps_rad: f64) -> f64 {
    let (sin_x, cos_x) = x_deg.to_radians().sin_cos();
    let (sin_e, cos_e) = eps_rad.sin_cos();
    let f_clamped = f_lat_deg.clamp(-89.999_999, 89.999_999);
    let tan_f = f_clamped.to_radians().tan();
    sin_x
        .atan2(cos_x * cos_e - tan_f * sin_e)
        .to_degrees()
        .rem_euclid(360.0)
}

/// Evaluates the 8 [`HouseAngles`] from `armc_deg`, geographic `lat_deg`, and obliquity `eps_rad` (`MODEL_PASSPORTS.md` §5.2).
#[must_use]
pub fn compute_angles(armc_deg: f64, lat_deg: f64, eps_rad: f64) -> HouseAngles {
    let armc = armc_deg.rem_euclid(360.0);
    let mc_deg = asc1(armc, 0.0, eps_rad);
    let mut asc_deg = asc1(armc + 90.0, lat_deg, eps_rad);
    let equatorial_asc_deg = asc1(armc + 90.0, 0.0, eps_rad);
    let polar_asc_deg = asc1((armc + 180.0).rem_euclid(360.0) + 90.0, lat_deg, eps_rad);
    let co_asc_koch_deg = (polar_asc_deg + 180.0).rem_euclid(360.0);

    let (mut vertex_deg, co_asc_munkasey_deg) = if lat_deg.abs() > 1e-12 {
        let cot_lat = 1.0 / lat_deg.clamp(-89.999_999, 89.999_999).to_radians().tan();
        let (sin_r, cos_r) = armc.to_radians().sin_cos();
        let (sin_e, cos_e) = eps_rad.sin_cos();
        let v = (-cos_r)
            .atan2(sin_r * cos_e - cot_lat * sin_e)
            .to_degrees()
            .rem_euclid(360.0);
        let cm = cos_r
            .atan2(-sin_r * cos_e - cot_lat * sin_e)
            .to_degrees()
            .rem_euclid(360.0);
        (v, cm)
    } else {
        (180.0, 0.0)
    };

    // Inside the polar circles (|lat| >= 90 - eps), when the ecliptic pole is below the horizon,
    // asc1(armc + 90, lat, eps) returns the western (setting, dh/dt < 0) horizon intersection where
    // (asc - mc) mod 360 >= 180. Per Swiss Ephemeris (`swehouse.c`) and canonical definition
    // (Lilly, Christian Astrology I, p. 51: "he ascends, or then ariseth, and is visible in our Horizon"),
    // the Ascendant is always the eastern (physically rising, dh/dt > 0) horizon intersection where
    // (asc - mc) mod 360 < 180.
    if (asc_deg - mc_deg).rem_euclid(360.0) >= 180.0 {
        asc_deg = (asc_deg + 180.0).rem_euclid(360.0);
    }
    let v_rel = (vertex_deg - mc_deg).rem_euclid(360.0);
    if v_rel > 1.0e-12 && v_rel < 180.0 {
        vertex_deg = (vertex_deg + 180.0).rem_euclid(360.0);
    }

    HouseAngles {
        asc_deg,
        mc_deg,
        armc_deg: armc,
        vertex_deg,
        equatorial_asc_deg,
        co_asc_koch_deg,
        co_asc_munkasey_deg,
        polar_asc_deg,
    }
}

/// Computes Porphyry (`O`) 12 house cusps from `mc_deg` and `asc_deg`.
#[must_use]
pub fn porphyry_cusps(mc_deg: f64, asc_deg: f64) -> [f64; 12] {
    let ic_deg = (mc_deg + 180.0).rem_euclid(360.0);
    let dsc_deg = (asc_deg + 180.0).rem_euclid(360.0);
    let q1 = (asc_deg - mc_deg).rem_euclid(360.0) / 3.0;
    let q2 = (ic_deg - asc_deg).rem_euclid(360.0) / 3.0;

    let c10 = mc_deg;
    let c11 = (mc_deg + q1).rem_euclid(360.0);
    let c12 = (mc_deg + 2.0 * q1).rem_euclid(360.0);
    let c1 = asc_deg;
    let c2 = (asc_deg + q2).rem_euclid(360.0);
    let c3 = (asc_deg + 2.0 * q2).rem_euclid(360.0);
    let c4 = ic_deg;
    let c5 = (c11 + 180.0).rem_euclid(360.0);
    let c6 = (c12 + 180.0).rem_euclid(360.0);
    let c7 = dsc_deg;
    let c8 = (c2 + 180.0).rem_euclid(360.0);
    let c9 = (c3 + 180.0).rem_euclid(360.0);

    [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12]
}

/// Solves $g(\Delta\alpha) = \Delta\alpha - f \arccos(-\sin(\text{RAMC} + \Delta\alpha)\tan\phi\tan\varepsilon) = 0$
/// (`MODEL_PASSPORTS.md` §5.3.1), returning $\Delta\alpha$ in degrees via the Illinois method ($10^{-13}\text{ rad}$ tolerance, max 50 iterations).
fn solve_placidus_offset_deg(
    armc_deg: f64,
    lat_deg: f64,
    eps_rad: f64,
    f: f64,
) -> Result<f64, EphemerisError> {
    let tt = lat_deg.to_radians().tan() * eps_rad.tan();
    if tt.abs() >= 1.0 {
        return Err(EphemerisError::UndefinedHouseSystem {
            house_system: 'P',
            latitude_deg: lat_deg,
            reason: "Placidus semi-arc undefined inside polar circle (|tan(lat)*tan(eps)| >= 1)"
                .to_owned(),
        });
    }
    let r = armc_deg.to_radians();
    let eval_g = |da: f64| -> Result<f64, EphemerisError> {
        let arg = -(r + da).sin() * tt;
        if arg.abs() >= 1.0 {
            return Err(EphemerisError::UndefinedHouseSystem {
                house_system: 'P',
                latitude_deg: lat_deg,
                reason: "Circumpolar ecliptic point in Placidus".to_owned(),
            });
        }
        Ok(da - f * arg.acos())
    };
    let (mut a, mut b) = if f > 0.0 {
        (0.0, std::f64::consts::PI)
    } else {
        (-std::f64::consts::PI, 0.0)
    };
    let mut fa = eval_g(a)?;
    let mut fb = eval_g(b)?;
    let tol_rad = 1.0e-13;
    let mut side = 0i8;
    for _ in 0..50 {
        if (b - a).abs() <= tol_rad || fb.abs() <= tol_rad {
            return Ok(b.to_degrees());
        }
        let denom = fb - fa;
        let c = if denom.abs() > 1e-15 {
            (a * fb - b * fa) / denom
        } else {
            0.5 * (a + b)
        };
        let fc = eval_g(c)?;
        if fc.abs() <= tol_rad {
            return Ok(c.to_degrees());
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
        solver: "solve_placidus_offset_deg".to_owned(),
        iterations: 50,
        residual: fb.abs(),
    })
}

/// Spherical linear interpolation (`slerp`) between unit vectors `u` and `v` at `t` in `[0, 1]`.
fn slerp_vec3(u: [f64; 3], v: [f64; 3], t: f64) -> [f64; 3] {
    let dot = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]).clamp(-1.0, 1.0);
    let omega = dot.acos();
    if omega.abs() < 1e-12 {
        return u;
    }
    let sin_om = omega.sin();
    let w1 = ((1.0 - t) * omega).sin() / sin_om;
    let w2 = (t * omega).sin() / sin_om;
    [
        w1 * u[0] + w2 * v[0],
        w1 * u[1] + w2 * v[1],
        w1 * u[2] + w2 * v[2],
    ]
}

/// Computes the ecliptic longitude `[0, 360)` where the great circle with `normal` intersects the ecliptic
/// in the half-space of `target_dir` (`MODEL_PASSPORTS.md` §5.3.5 & §5.3.23).
fn great_circle_ecl_lon(normal: [f64; 3], target_dir: [f64; 3], eps_rad: f64) -> f64 {
    let (sin_e, cos_e) = eps_rad.sin_cos();
    let mut lam = (-normal[0]).atan2(normal[1] * cos_e + normal[2] * sin_e);
    let (sin_l, cos_l) = lam.sin_cos();
    let v = [cos_l, sin_l * cos_e, sin_l * sin_e];
    let dot = v[0] * target_dir[0] + v[1] * target_dir[1] + v[2] * target_dir[2];
    if dot < 0.0 {
        lam += std::f64::consts::PI;
    }
    lam.to_degrees().rem_euclid(360.0)
}

/// Evaluates the 12 house cusps for `system` at `(armc_deg, lat_deg, eps_rad, sun_dec_deg)`
/// according to `MODEL_PASSPORTS.md` §5.3.
fn eval_raw_house_cusps(
    system: HouseSystem,
    _armc_deg: f64,
    lat_deg: f64,
    eps_rad: f64,
    sun_dec_deg: f64,
    angles: &HouseAngles,
) -> Result<[f64; 12], EphemerisError> {
    let mc = angles.mc_deg;
    let asc = angles.asc_deg;
    let ic = (mc + 180.0).rem_euclid(360.0);
    let dsc = (asc + 180.0).rem_euclid(360.0);
    let armc = angles.armc_deg;

    let wrap360 = |x: f64| -> f64 {
        let r = x.rem_euclid(360.0);
        if r >= 360.0 { 0.0 } else { r }
    };

    let sym_cusps = |c11: f64, c12: f64, c2: f64, c3: f64| -> [f64; 12] {
        [
            wrap360(asc),
            wrap360(c2),
            wrap360(c3),
            wrap360(ic),
            wrap360(c11 + 180.0),
            wrap360(c12 + 180.0),
            wrap360(dsc),
            wrap360(c2 + 180.0),
            wrap360(c3 + 180.0),
            wrap360(mc),
            wrap360(c11),
            wrap360(c12),
        ]
    };

    let r_rad = armc.to_radians();
    let p_rad = lat_deg.clamp(-89.999_999, 89.999_999).to_radians();
    let east = [-r_rad.sin(), r_rad.cos(), 0.0];
    let zenith = [
        p_rad.cos() * r_rad.cos(),
        p_rad.cos() * r_rad.sin(),
        p_rad.sin(),
    ];
    let north = crate::frame::vec_cross(&zenith, &east);

    match system {
        HouseSystem::Placidus => {
            let d11 = solve_placidus_offset_deg(armc, lat_deg, eps_rad, 1.0 / 3.0)?;
            let d12 = solve_placidus_offset_deg(armc, lat_deg, eps_rad, 2.0 / 3.0)?;
            let d8 = solve_placidus_offset_deg(armc, lat_deg, eps_rad, -2.0 / 3.0)?;
            let d9 = solve_placidus_offset_deg(armc, lat_deg, eps_rad, -1.0 / 3.0)?;
            let c11 = asc1(armc + d11, 0.0, eps_rad);
            let c12 = asc1(armc + d12, 0.0, eps_rad);
            let c8 = asc1(armc + d8, 0.0, eps_rad);
            let c9 = asc1(armc + d9, 0.0, eps_rad);
            let c2 = (c8 + 180.0).rem_euclid(360.0);
            let c3 = (c9 + 180.0).rem_euclid(360.0);
            Ok(sym_cusps(c11, c12, c2, c3))
        }
        HouseSystem::Koch => {
            let sin_dec_mc = mc.to_radians().sin() * eps_rad.sin();
            let cos_dec_mc = (1.0 - sin_dec_mc * sin_dec_mc).max(1e-15).sqrt();
            let tan_dec_mc = sin_dec_mc / cos_dec_mc;
            let x_ad = lat_deg.to_radians().tan() * tan_dec_mc;
            if x_ad.abs() >= 1.0 || lat_deg.abs() >= 90.0 - eps_rad.to_degrees() {
                return Err(EphemerisError::UndefinedHouseSystem {
                    house_system: 'K',
                    latitude_deg: lat_deg,
                    reason: "Koch MC ascensional difference undefined (|tan(lat)*tan(dec_MC)| >= 1)"
                        .to_owned(),
                });
            }
            let ad_mc = x_ad.asin().to_degrees();
            let c11 = asc1(armc + 30.0 - 2.0 * ad_mc / 3.0, lat_deg, eps_rad);
            let c12 = asc1(armc + 60.0 - ad_mc / 3.0, lat_deg, eps_rad);
            let c2 = asc1(armc + 120.0 + ad_mc / 3.0, lat_deg, eps_rad);
            let c3 = asc1(armc + 150.0 + 2.0 * ad_mc / 3.0, lat_deg, eps_rad);
            Ok(sym_cusps(c11, c12, c2, c3))
        }
        HouseSystem::Porphyry => Ok(porphyry_cusps(mc, asc)),
        HouseSystem::Regiomontanus => {
            let tan_lat = p_rad.tan();
            let regio_cusp = |h_deg: f64| -> f64 {
                let rh_rad = (armc + h_deg).to_radians();
                let h_rad = h_deg.to_radians();
                let (sin_rh, cos_rh) = rh_rad.sin_cos();
                let (sin_e, cos_e) = eps_rad.sin_cos();
                sin_rh
                    .atan2(cos_rh * cos_e - h_rad.sin() * tan_lat * sin_e)
                    .to_degrees()
                    .rem_euclid(360.0)
            };
            Ok(sym_cusps(
                regio_cusp(30.0),
                regio_cusp(60.0),
                regio_cusp(120.0),
                regio_cusp(150.0),
            ))
        }
        HouseSystem::Campanus => {
            let camp_cusp = |h_deg: f64| -> f64 {
                let elev = (90.0 - h_deg).to_radians();
                let (sin_el, cos_el) = elev.sin_cos();
                let pt = [
                    cos_el * east[0] + sin_el * zenith[0],
                    cos_el * east[1] + sin_el * zenith[1],
                    cos_el * east[2] + sin_el * zenith[2],
                ];
                let normal = crate::frame::vec_cross(&north, &pt);
                great_circle_ecl_lon(normal, pt, eps_rad)
            };
            Ok(sym_cusps(
                camp_cusp(30.0),
                camp_cusp(60.0),
                camp_cusp(120.0),
                camp_cusp(150.0),
            ))
        }
        HouseSystem::EqualFromAsc => {
            let mut out = [0.0; 12];
            for (i, slot) in out.iter_mut().enumerate() {
                *slot = (asc + (i as f64) * 30.0).rem_euclid(360.0);
            }
            Ok(out)
        }
        HouseSystem::WholeSign => {
            let base = (asc / 30.0).floor() * 30.0;
            let mut out = [0.0; 12];
            for (i, slot) in out.iter_mut().enumerate() {
                *slot = (base + (i as f64) * 30.0).rem_euclid(360.0);
            }
            Ok(out)
        }
        HouseSystem::Alcabitius => {
            let sin_dec_asc = asc.to_radians().sin() * eps_rad.sin();
            let cos_dec_asc = (1.0 - sin_dec_asc * sin_dec_asc).max(1e-15).sqrt();
            let tan_dec_asc = sin_dec_asc / cos_dec_asc;
            let x_ad = (p_rad.tan() * tan_dec_asc).clamp(-1.0, 1.0);
            let sday = 90.0 + x_ad.asin().to_degrees();
            let snight = 180.0 - sday;
            let c11 = asc1(armc + sday / 3.0, 0.0, eps_rad);
            let c12 = asc1(armc + 2.0 * sday / 3.0, 0.0, eps_rad);
            let c2 = asc1(armc + sday + snight / 3.0, 0.0, eps_rad);
            let c3 = asc1(armc + sday + 2.0 * snight / 3.0, 0.0, eps_rad);
            Ok(sym_cusps(c11, c12, c2, c3))
        }
        HouseSystem::Morinus => {
            let cos_e = eps_rad.cos();
            let mut out = [0.0; 12];
            for (i, slot) in out.iter_mut().enumerate() {
                let ra_deg = armc + 90.0 + (i as f64) * 30.0;
                let (sin_a, cos_a) = ra_deg.to_radians().sin_cos();
                *slot = (sin_a * cos_e).atan2(cos_a).to_degrees().rem_euclid(360.0);
            }
            Ok(out)
        }
        HouseSystem::PolichPage => {
            let tan_lat = p_rad.tan();
            let f1 = (tan_lat / 3.0).atan().to_degrees();
            let f2 = (2.0 * tan_lat / 3.0).atan().to_degrees();
            let c11 = asc1(armc + 30.0, f1, eps_rad);
            let c12 = asc1(armc + 60.0, f2, eps_rad);
            let c2 = asc1(armc + 120.0, f2, eps_rad);
            let c3 = asc1(armc + 150.0, f1, eps_rad);
            Ok(sym_cusps(c11, c12, c2, c3))
        }
        HouseSystem::KrusinskiPisaGoelzer => {
            let (sin_l, cos_l) = asc.to_radians().sin_cos();
            let (sin_e, cos_e) = eps_rad.sin_cos();
            let v_asc = [cos_l, sin_l * cos_e, sin_l * sin_e];
            let nadir = [-zenith[0], -zenith[1], -zenith[2]];
            let krus_pt = |target: [f64; 3], frac: f64| -> f64 {
                let p = slerp_vec3(v_asc, target, frac);
                let ra_deg = p[1].atan2(p[0]).to_degrees();
                asc1(ra_deg, 0.0, eps_rad)
            };
            let c11 = krus_pt(zenith, 2.0 / 3.0);
            let c12 = krus_pt(zenith, 1.0 / 3.0);
            let c2 = krus_pt(nadir, 1.0 / 3.0);
            let c3 = krus_pt(nadir, 2.0 / 3.0);
            Ok(sym_cusps(c11, c12, c2, c3))
        }
        HouseSystem::Meridian => {
            let mut out = [0.0; 12];
            for (i, slot) in out.iter_mut().enumerate() {
                *slot = asc1(armc + 90.0 + (i as f64) * 30.0, 0.0, eps_rad);
            }
            Ok(out)
        }
        HouseSystem::VehlowEqual => {
            let base = (asc - 15.0).rem_euclid(360.0);
            let mut out = [0.0; 12];
            for (i, slot) in out.iter_mut().enumerate() {
                *slot = (base + (i as f64) * 30.0).rem_euclid(360.0);
            }
            Ok(out)
        }
        HouseSystem::PullenSinusoidalDelta => {
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
            Ok(sym_cusps(mc + w4, asc - w4, asc + w1, ic - w1))
        }
        HouseSystem::PullenSinusoidalRatio => {
            let q4 = (asc - mc).rem_euclid(360.0);
            let q1 = (ic - asc).rem_euclid(360.0);
            let q_small = q4.min(q1);
            let (w_small, w_large) = if q_small < 1.0e-9 {
                (0.0, 60.0)
            } else {
                let ratio_q = (180.0 - q_small) / q_small;
                let mut r = ((2.0 * ratio_q).cbrt() + 1.0).max(1.0);
                for _ in 0..30 {
                    let r2 = r * r;
                    let r3 = r2 * r;
                    let f_val = r3 * (r + 2.0) - ratio_q * (2.0 * r + 1.0);
                    let df_val = 4.0 * r3 + 6.0 * r2 - 2.0 * ratio_q;
                    if df_val.abs() < 1.0e-14 {
                        break;
                    }
                    r -= f_val / df_val;
                }
                let x_base = q_small / (2.0 * r + 1.0);
                (r * x_base, r * r * r * x_base)
            };
            let (w4_sr, w1_sr) = if q4 <= q1 {
                (w_small, w_large)
            } else {
                (w_large, w_small)
            };
            Ok(sym_cusps(mc + w4_sr, asc - w4_sr, asc + w1_sr, ic - w1_sr))
        }
        HouseSystem::SunshineTreindl | HouseSystem::SunshineMakransky => {
            let sin_kv = sun_dec_deg.to_radians().tan() * p_rad.tan();
            if system == HouseSystem::SunshineMakransky
                && (sin_kv.abs() >= 1.0 || lat_deg.abs() >= 90.0 - eps_rad.to_degrees())
            {
                return Err(EphemerisError::UndefinedHouseSystem {
                    house_system: 'i',
                    latitude_deg: lat_deg,
                    reason: "Sunshine (Makransky) solar semi-arc undefined inside polar circle or |tan(lat)*tan(dec_sun)| >= 1"
                        .to_owned(),
                });
            }
            let kv_clamped = sin_kv.clamp(-1.0, 1.0).asin();
            let (sin_e, cos_e) = eps_rad.sin_cos();
            let (sin_r, cos_r) = r_rad.sin_cos();
            let tan_p = p_rad.tan();
            let mut cusps = [0.0; 12];
            cusps[0] = asc;
            cusps[3] = ic;
            cusps[6] = dsc;
            cusps[9] = mc;
            for h in [2usize, 3, 5, 6, 8, 9, 11, 12] {
                if sin_kv >= 1.0 && h < 8 {
                    if let Some(slot) = cusps.get_mut(h - 1) {
                        *slot = ic;
                    }
                    continue;
                }
                if sin_kv <= -1.0 && h >= 8 {
                    if let Some(slot) = cusps.get_mut(h - 1) {
                        *slot = mc;
                    }
                    continue;
                }
                let h_f = h as f64;
                let a = if h < 8 {
                    kv_clamped
                        + r_rad
                        + std::f64::consts::FRAC_PI_2
                        + ((h_f - 1.0) / 3.0) * (std::f64::consts::FRAC_PI_2 - kv_clamped)
                } else {
                    kv_clamped
                        + r_rad
                        + std::f64::consts::FRAC_PI_2
                        + ((h_f - 13.0) / 3.0) * (std::f64::consts::FRAC_PI_2 + kv_clamped)
                };
                let y_val = sin_kv * sin_r + a.sin();
                let x_val = cos_e * (sin_kv * cos_r + a.cos()) + sin_e * tan_p * (r_rad - a).sin();
                if let Some(slot) = cusps.get_mut(h - 1) {
                    *slot = y_val.atan2(x_val).to_degrees().rem_euclid(360.0);
                }
            }
            Ok(cusps)
        }
        HouseSystem::CarterPoliEquatorial => {
            let (sin_l, cos_l) = asc.to_radians().sin_cos();
            let ra_asc_deg = (sin_l * eps_rad.cos())
                .atan2(cos_l)
                .to_degrees()
                .rem_euclid(360.0);
            let mut out = [0.0; 12];
            for (i, slot) in out.iter_mut().enumerate() {
                *slot = asc1(ra_asc_deg + (i as f64) * 30.0, 0.0, eps_rad);
            }
            Ok(out)
        }
        HouseSystem::WholeSignAries => {
            let mut out = [0.0; 12];
            for (i, slot) in out.iter_mut().enumerate() {
                *slot = ((i as f64) * 30.0).rem_euclid(360.0);
            }
            Ok(out)
        }
        HouseSystem::Sripati => {
            let p = porphyry_cusps(mc, asc);
            let mut out = [0.0; 12];
            for i in 0..12 {
                let prev = p.get((i + 11) % 12).copied().unwrap_or(0.0);
                let curr = p.get(i).copied().unwrap_or(0.0);
                let span = (curr - prev).rem_euclid(360.0);
                if let Some(slot) = out.get_mut(i) {
                    *slot = (prev + 0.5 * span).rem_euclid(360.0);
                }
            }
            Ok(out)
        }
        HouseSystem::EqualFromMc => {
            let c1 = (mc + 90.0).rem_euclid(360.0);
            let mut out = [0.0; 12];
            for (i, slot) in out.iter_mut().enumerate() {
                *slot = (c1 + (i as f64) * 30.0).rem_euclid(360.0);
            }
            Ok(out)
        }
        HouseSystem::Apc => {
            let tan_p = p_rad.tan();
            let tan_e = eps_rad.tan();
            let (sin_r, cos_r) = r_rad.sin_cos();
            let (sin_e, cos_e) = eps_rad.sin_cos();
            let num = tan_p * tan_e * cos_r;
            let denom = 1.0 + tan_p * tan_e * sin_r;
            let kv_apc = if denom.abs() < 1.0e-12 {
                if num.abs() < 1.0e-12 {
                    0.0
                } else if num >= 0.0 {
                    std::f64::consts::FRAC_PI_2
                } else {
                    -std::f64::consts::FRAC_PI_2
                }
            } else {
                (num / denom).atan()
            };
            let sin_kv = kv_apc.sin();
            let mut cusps = [0.0; 12];
            cusps[0] = asc;
            cusps[3] = ic;
            cusps[6] = dsc;
            cusps[9] = mc;
            for h in [2usize, 3, 5, 6, 8, 9, 11, 12] {
                let h_f = h as f64;
                let a = if h < 8 {
                    kv_apc
                        + r_rad
                        + std::f64::consts::FRAC_PI_2
                        + ((h_f - 1.0) / 3.0) * (std::f64::consts::FRAC_PI_2 - kv_apc)
                } else {
                    kv_apc
                        + r_rad
                        + std::f64::consts::FRAC_PI_2
                        + ((h_f - 13.0) / 3.0) * (std::f64::consts::FRAC_PI_2 + kv_apc)
                };
                let y_val = sin_kv * sin_r + a.sin();
                let x_val = cos_e * (sin_kv * cos_r + a.cos()) + sin_e * tan_p * (r_rad - a).sin();
                if let Some(slot) = cusps.get_mut(h - 1) {
                    *slot = y_val.atan2(x_val).to_degrees().rem_euclid(360.0);
                }
            }
            Ok(cusps)
        }
        HouseSystem::Azimuthal => {
            let sgn = if lat_deg > 0.0 { 1.0 } else { -1.0 };
            let (sin_e, cos_e) = eps_rad.sin_cos();
            let horiz_cusp = |h_deg: f64| -> f64 {
                let (sin_t, cos_t) = (h_deg - 90.0).to_radians().sin_cos();
                let nx = sgn * cos_t * north[0] - sin_t * east[0];
                let ny = sgn * cos_t * north[1] - sin_t * east[1];
                let nz = sgn * cos_t * north[2] - sin_t * east[2];
                let mut lam = (-nx)
                    .atan2(ny * cos_e + nz * sin_e)
                    .to_degrees()
                    .rem_euclid(360.0);
                if lat_deg != 0.0 && (lam - mc).rem_euclid(360.0) >= 180.0 {
                    lam = (lam + 180.0).rem_euclid(360.0);
                }
                lam
            };
            let c11_h = horiz_cusp(30.0);
            let c12_h = horiz_cusp(60.0);
            let c1_h = horiz_cusp(90.0);
            let c2_h = horiz_cusp(120.0);
            let c3_h = horiz_cusp(150.0);
            Ok([
                c1_h,
                c2_h,
                c3_h,
                ic,
                (c11_h + 180.0).rem_euclid(360.0),
                (c12_h + 180.0).rem_euclid(360.0),
                (c1_h + 180.0).rem_euclid(360.0),
                (c2_h + 180.0).rem_euclid(360.0),
                (c3_h + 180.0).rem_euclid(360.0),
                mc,
                c11_h,
                c12_h,
            ])
        }
    }
}

/// Evaluates house cusps and angles directly from `armc_deg`, geographic `latitude_deg`,
/// ecliptic `obliquity_deg`, `sun_declination_deg`, and `sidereal_offset_deg` (`ARCHITECTURE.md` §2.4).
pub fn houses_armc_with_sidereal(
    armc_deg: f64,
    latitude_deg: f64,
    obliquity_deg: f64,
    sun_declination_deg: f64,
    sidereal_offset_deg: f64,
    options: HouseOptions,
) -> Result<HousesResult, EphemerisError> {
    if !armc_deg.is_finite()
        || !latitude_deg.is_finite()
        || !(-90.0..=90.0).contains(&latitude_deg)
        || !obliquity_deg.is_finite()
        || obliquity_deg <= 0.0
        || obliquity_deg >= 90.0
        || !sun_declination_deg.is_finite()
        || !(-90.0..=90.0).contains(&sun_declination_deg)
        || !sidereal_offset_deg.is_finite()
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid houses_armc inputs: armc={armc_deg}, lat={latitude_deg}, obl={obliquity_deg}, sun_dec={sun_declination_deg}, sid={sidereal_offset_deg}"
            ),
        });
    }

    let req_sys = options.system;
    let abs_lat = latitude_deg.abs();
    let polar_limit_deg = 90.0 - obliquity_deg;

    // Exact geographic poles |lat| >= 89.999999°: local horizon and meridian are undefined;
    // only 'N' (Whole Sign 0° Aries) is defined at exact poles (`ARCHITECTURE.md` §2.4).
    if abs_lat >= 89.999_999 && req_sys != HouseSystem::WholeSignAries {
        return Err(EphemerisError::UndefinedHouseSystem {
            house_system: req_sys.code(),
            latitude_deg,
            reason: "Geographic pole (|latitude| == 90°): local horizon and meridian are undefined"
                .to_owned(),
        });
    }

    let eps_rad = obliquity_deg.to_radians();
    let mut angles = compute_angles(armc_deg, latitude_deg, eps_rad);
    let mut warnings = Vec::new();
    if abs_lat >= 89.0 {
        warnings.push(HouseWarning::NearPoleBoundary);
    }

    let in_polar_circle = abs_lat >= polar_limit_deg;
    let mut effective_system = req_sys;
    let mut fallback_reason = None;

    let mut cusps_deg = match req_sys.polar_category() {
        PolarCategory::HardFail => {
            match eval_raw_house_cusps(
                req_sys,
                armc_deg,
                latitude_deg,
                eps_rad,
                sun_declination_deg,
                &angles,
            ) {
                Ok(c) => {
                    if !is_cusp_ring_monotonic(&c)
                        && !(options.polar_policy == PolarPolicy::FallbackToPorphyry
                            && req_sys == HouseSystem::PolichPage)
                    {
                        match options.polar_policy {
                            PolarPolicy::FallbackToWholeSign => {
                                effective_system = HouseSystem::WholeSign;
                                fallback_reason = Some(format!(
                                    "Cusp-inversion fallback from '{}' to 'W' (Whole Sign) at latitude {:.4}°",
                                    req_sys.code(),
                                    latitude_deg
                                ));
                                warnings.push(HouseWarning::PolarFallbackToWholeSign);
                                eval_raw_house_cusps(
                                    HouseSystem::WholeSign,
                                    armc_deg,
                                    latitude_deg,
                                    eps_rad,
                                    sun_declination_deg,
                                    &angles,
                                )?
                            }
                            PolarPolicy::FallbackToPorphyry => {
                                effective_system = HouseSystem::Porphyry;
                                fallback_reason = Some(format!(
                                    "Cusp-inversion fallback from '{}' to 'O' (Porphyry) at latitude {:.4}°",
                                    req_sys.code(),
                                    latitude_deg
                                ));
                                warnings.push(HouseWarning::PolarFallbackToPorphyry);
                                porphyry_cusps(angles.mc_deg, angles.asc_deg)
                            }
                            PolarPolicy::Strict => c,
                        }
                    } else {
                        c
                    }
                }
                Err(err) => match options.polar_policy {
                    PolarPolicy::Strict => return Err(err),
                    PolarPolicy::FallbackToPorphyry => {
                        effective_system = HouseSystem::Porphyry;
                        fallback_reason = Some(format!(
                            "Polar fallback from '{}' to 'O' (Porphyry) at latitude {:.4}°",
                            req_sys.code(),
                            latitude_deg
                        ));
                        warnings.push(HouseWarning::PolarFallbackToPorphyry);
                        porphyry_cusps(angles.mc_deg, angles.asc_deg)
                    }
                    PolarPolicy::FallbackToWholeSign => {
                        effective_system = HouseSystem::WholeSign;
                        fallback_reason = Some(format!(
                            "Polar fallback from '{}' to 'W' (Whole Sign) at latitude {:.4}°",
                            req_sys.code(),
                            latitude_deg
                        ));
                        warnings.push(HouseWarning::PolarFallbackToWholeSign);
                        eval_raw_house_cusps(
                            HouseSystem::WholeSign,
                            armc_deg,
                            latitude_deg,
                            eps_rad,
                            sun_declination_deg,
                            &angles,
                        )?
                    }
                },
            }
        }
        PolarCategory::SoftFail => {
            if in_polar_circle {
                warnings.push(HouseWarning::PolarUnreliable);
            }
            let c = eval_raw_house_cusps(
                req_sys,
                armc_deg,
                latitude_deg,
                eps_rad,
                sun_declination_deg,
                &angles,
            )?;
            let raw_asc_inverted = (asc1(armc_deg + 90.0, latitude_deg, eps_rad) - angles.mc_deg)
                .rem_euclid(360.0)
                >= 180.0;
            let ring_inverted = (!is_cusp_ring_monotonic(&c)
                && !(options.polar_policy == PolarPolicy::FallbackToPorphyry
                    && req_sys == HouseSystem::PolichPage
                    && !in_polar_circle))
                || (in_polar_circle && raw_asc_inverted);
            if ring_inverted {
                match options.polar_policy {
                    PolarPolicy::FallbackToWholeSign => {
                        effective_system = HouseSystem::WholeSign;
                        fallback_reason = Some(format!(
                            "Polar/cusp-inversion fallback from '{}' to 'W' (Whole Sign) at latitude {:.4}°",
                            req_sys.code(),
                            latitude_deg
                        ));
                        warnings.push(HouseWarning::PolarFallbackToWholeSign);
                        eval_raw_house_cusps(
                            HouseSystem::WholeSign,
                            armc_deg,
                            latitude_deg,
                            eps_rad,
                            sun_declination_deg,
                            &angles,
                        )?
                    }
                    PolarPolicy::FallbackToPorphyry | PolarPolicy::Strict => c,
                }
            } else {
                c
            }
        }
        PolarCategory::PolarSafe => {
            let c = eval_raw_house_cusps(
                req_sys,
                armc_deg,
                latitude_deg,
                eps_rad,
                sun_declination_deg,
                &angles,
            )?;
            if options.polar_policy == PolarPolicy::FallbackToWholeSign && !is_cusp_ring_monotonic(&c)
            {
                effective_system = HouseSystem::WholeSign;
                fallback_reason = Some(format!(
                    "Polar/cusp-inversion fallback from '{}' to 'W' (Whole Sign) at latitude {:.4}°",
                    req_sys.code(),
                    latitude_deg
                ));
                warnings.push(HouseWarning::PolarFallbackToWholeSign);
                eval_raw_house_cusps(
                    HouseSystem::WholeSign,
                    armc_deg,
                    latitude_deg,
                    eps_rad,
                    sun_declination_deg,
                    &angles,
                )?
            } else {
                c
            }
        }
    };

    // Apply sidereal offset if requested
    if sidereal_offset_deg.abs() > 0.0 {
        if effective_system == HouseSystem::WholeSign {
            let sid_asc = (angles.asc_deg - sidereal_offset_deg).rem_euclid(360.0);
            let base = (sid_asc / 30.0).floor() * 30.0;
            for (i, slot) in cusps_deg.iter_mut().enumerate() {
                *slot = (base + (i as f64) * 30.0).rem_euclid(360.0);
            }
        } else if effective_system != HouseSystem::WholeSignAries {
            for slot in &mut cusps_deg {
                *slot = (*slot - sidereal_offset_deg).rem_euclid(360.0);
            }
        }
        angles.asc_deg = (angles.asc_deg - sidereal_offset_deg).rem_euclid(360.0);
        angles.mc_deg = (angles.mc_deg - sidereal_offset_deg).rem_euclid(360.0);
        angles.vertex_deg = (angles.vertex_deg - sidereal_offset_deg).rem_euclid(360.0);
        angles.equatorial_asc_deg =
            (angles.equatorial_asc_deg - sidereal_offset_deg).rem_euclid(360.0);
        angles.co_asc_koch_deg = (angles.co_asc_koch_deg - sidereal_offset_deg).rem_euclid(360.0);
        angles.co_asc_munkasey_deg =
            (angles.co_asc_munkasey_deg - sidereal_offset_deg).rem_euclid(360.0);
        angles.polar_asc_deg = (angles.polar_asc_deg - sidereal_offset_deg).rem_euclid(360.0);
    }

    // Compute cusp and angle speeds (deg/day) via finite difference in ARMC
    // (d(ARMC)/dt = 360.98564736629 deg/day_UT1)
    let (cusp_speeds_deg_day, angle_speeds) = if options.include_speed {
        let darmc_dt = 360.985_647_366_29;
        let h_days = 1.0e-4;
        let d_armc = darmc_dt * h_days;
        let mut opts_no_spd = options;
        opts_no_spd.include_speed = false;

        if let (Ok(r_m), Ok(r_p)) = (
            houses_armc_with_sidereal(
                (armc_deg - d_armc).rem_euclid(360.0),
                latitude_deg,
                obliquity_deg,
                sun_declination_deg,
                sidereal_offset_deg,
                opts_no_spd,
            ),
            houses_armc_with_sidereal(
                (armc_deg + d_armc).rem_euclid(360.0),
                latitude_deg,
                obliquity_deg,
                sun_declination_deg,
                sidereal_offset_deg,
                opts_no_spd,
            ),
        ) {
            let diff_rate = |a_m: f64, a_p: f64| -> f64 {
                let d = (a_p - a_m + 180.0).rem_euclid(360.0) - 180.0;
                d / (2.0 * h_days)
            };
            let mut c_spd = [0.0; 12];
            if effective_system == HouseSystem::WholeSign
                || effective_system == HouseSystem::WholeSignAries
            {
                c_spd = [0.0; 12];
            } else {
                for i in 0..12 {
                    let cm = r_m.cusps_deg.get(i).copied().unwrap_or(0.0);
                    let cp = r_p.cusps_deg.get(i).copied().unwrap_or(0.0);
                    if let Some(slot) = c_spd.get_mut(i) {
                        *slot = diff_rate(cm, cp);
                    }
                }
            }
            let a_spd = HouseAngles {
                asc_deg: diff_rate(r_m.angles.asc_deg, r_p.angles.asc_deg),
                mc_deg: diff_rate(r_m.angles.mc_deg, r_p.angles.mc_deg),
                armc_deg: darmc_dt,
                vertex_deg: diff_rate(r_m.angles.vertex_deg, r_p.angles.vertex_deg),
                equatorial_asc_deg: diff_rate(
                    r_m.angles.equatorial_asc_deg,
                    r_p.angles.equatorial_asc_deg,
                ),
                co_asc_koch_deg: diff_rate(
                    r_m.angles.co_asc_koch_deg,
                    r_p.angles.co_asc_koch_deg,
                ),
                co_asc_munkasey_deg: diff_rate(
                    r_m.angles.co_asc_munkasey_deg,
                    r_p.angles.co_asc_munkasey_deg,
                ),
                polar_asc_deg: diff_rate(r_m.angles.polar_asc_deg, r_p.angles.polar_asc_deg),
            };
            (c_spd, a_spd)
        } else {
            (
                [0.0; 12],
                HouseAngles {
                    asc_deg: 0.0,
                    mc_deg: 0.0,
                    armc_deg: darmc_dt,
                    vertex_deg: 0.0,
                    equatorial_asc_deg: 0.0,
                    co_asc_koch_deg: 0.0,
                    co_asc_munkasey_deg: 0.0,
                    polar_asc_deg: 0.0,
                },
            )
        }
    } else {
        (
            [0.0; 12],
            HouseAngles {
                asc_deg: 0.0,
                mc_deg: 0.0,
                armc_deg: 0.0,
                vertex_deg: 0.0,
                equatorial_asc_deg: 0.0,
                co_asc_koch_deg: 0.0,
                co_asc_munkasey_deg: 0.0,
                polar_asc_deg: 0.0,
            },
        )
    };

    Ok(HousesResult {
        requested_house_system: req_sys,
        effective_house_system: effective_system,
        polar_house_fallback: effective_system != req_sys,
        fallback_reason,
        warnings,
        cusps_deg,
        cusp_speeds_deg_day,
        angles,
        angle_speeds,
    })
}

/// Evaluates house cusps and angles directly from `armc_deg`, `latitude_deg`, `obliquity_deg`,
/// and `sun_declination_deg` (Tropical).
pub fn houses_armc(
    armc_deg: f64,
    latitude_deg: f64,
    obliquity_deg: f64,
    sun_declination_deg: f64,
    options: HouseOptions,
) -> Result<HousesResult, EphemerisError> {
    houses_armc_with_sidereal(
        armc_deg,
        latitude_deg,
        obliquity_deg,
        sun_declination_deg,
        0.0,
        options,
    )
}

/// Evaluates house cusps and angles inside an existing [`Moment`] for `observer`, reusing the
/// precomputed `GAST`, obliquity, and frame matrices at zero `nut00a` cost.
pub fn calc_houses_in_moment(
    moment: &Moment<'_>,
    observer: Observer,
    options: HouseOptions,
) -> Result<HousesResult, EphemerisError> {
    let center = moment.center();
    let (eps_rad, gst_rad) = match options.nutation {
        NutationMode::True => (center.frame.true_eps_rad, center.gast_rad),
        NutationMode::Mean => (center.frame.eps0_rad, center.gmst_rad),
    };
    let armc_deg = (gst_rad.to_degrees() + observer.longitude_deg).rem_euclid(360.0);

    let sun_dec_deg = if matches!(
        options.system,
        HouseSystem::SunshineTreindl | HouseSystem::SunshineMakransky | HouseSystem::Apc
    ) {
        let sun_opts = CalcOptions {
            origin: Origin::Geocentric,
            plane: Plane::Equatorial,
            equinox: Equinox::OfDate,
            nutation: options.nutation,
            aberration: AbMode::Apparent,
            zodiac: Zodiac::Tropical,
            include_speed: false,
        };
        let sun_res = crate::astrometry::calc_spk_body_in_moment(
            moment,
            crate::types::BodyId::Planet(PlanetId::Sun),
            &sun_opts,
        )?;
        sun_res.lat_deg
    } else {
        0.0
    };

    let sidereal_offset_deg = match options.zodiac {
        Zodiac::Tropical => 0.0,
        Zodiac::Sidereal(model) => moment.ayanamsa(model, options.nutation)?.effective_deg,
    };

    houses_armc_with_sidereal(
        armc_deg,
        observer.latitude_deg,
        eps_rad.to_degrees(),
        sun_dec_deg,
        sidereal_offset_deg,
        options,
    )
}

/// Evaluates house cusps and angles at `input` for `observer` using `spk` and `time_table`.
pub fn calc_houses(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    input: TimeInput,
    observer: Observer,
    options: HouseOptions,
) -> Result<HousesResult, EphemerisError> {
    let (ut1, tt, _, _) = time_table.resolve_scales(input, Some(observer));
    let jd_tt = tt.jd_f64();

    let ft = FrameTransform::at_tt(tt);
    let (eps_rad, gst_rad) = match options.nutation {
        NutationMode::True => (
            ft.true_eps_rad,
            gst06(ut1.0.jd1, ut1.0.jd2, tt.0.jd1, tt.0.jd2, &ft.r_bpn),
        ),
        NutationMode::Mean => (ft.eps0_rad, crate::time::gmst_rad(ut1, tt)),
    };
    let armc_deg = (gst_rad.to_degrees() + observer.longitude_deg).rem_euclid(360.0);

    // Solar declination (only needed for Sunshine I/i and APC Y, or cheap from Moment)
    let sun_dec_deg = if matches!(
        options.system,
        HouseSystem::SunshineTreindl | HouseSystem::SunshineMakransky | HouseSystem::Apc
    ) {
        let m = crate::moment::Moment::new(spk, time_table, input, Some(observer))?;
        let sun_opts = CalcOptions {
            origin: Origin::Geocentric,
            plane: Plane::Equatorial,
            equinox: Equinox::OfDate,
            nutation: options.nutation,
            aberration: AbMode::Apparent,
            zodiac: Zodiac::Tropical,
            include_speed: false,
        };
        let sun_res = crate::astrometry::calc_spk_body_in_moment(
            &m,
            crate::types::BodyId::Planet(PlanetId::Sun),
            &sun_opts,
        )?;
        sun_res.lat_deg
    } else {
        0.0
    };

    let sidereal_offset_deg = match options.zodiac {
        Zodiac::Tropical => 0.0,
        Zodiac::Sidereal(model) => {
            eval_ayanamsa(spk, model, jd_tt, options.nutation)?.effective_deg
        }
    };

    houses_armc_with_sidereal(
        armc_deg,
        observer.latitude_deg,
        eps_rad.to_degrees(),
        sun_dec_deg,
        sidereal_offset_deg,
        options,
    )
}

/// Evaluates the 36 Gauquelin sectors (`1..=36`, clockwise from Ascendant through MC, DSC, IC)
/// at `(armc_deg, latitude_deg, obliquity_deg)` with `polar_policy`.
pub fn gauquelin_sectors(
    armc_deg: f64,
    latitude_deg: f64,
    obliquity_deg: f64,
    polar_policy: PolarPolicy,
) -> Result<GauquelinResult, EphemerisError> {
    if !armc_deg.is_finite()
        || !latitude_deg.is_finite()
        || !(-90.0..=90.0).contains(&latitude_deg)
        || !obliquity_deg.is_finite()
        || obliquity_deg <= 0.0
        || obliquity_deg >= 90.0
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid gauquelin_sectors inputs: armc={armc_deg}, lat={latitude_deg}, obl={obliquity_deg}"
            ),
        });
    }
    let abs_lat = latitude_deg.abs();
    if abs_lat >= 89.999_999 {
        return Err(EphemerisError::UndefinedHouseSystem {
            house_system: 'G',
            latitude_deg,
            reason: "Geographic pole (|latitude| == 90°): Gauquelin sectors are undefined"
                .to_owned(),
        });
    }
    let eps_rad = obliquity_deg.to_radians();
    let angles = compute_angles(armc_deg, latitude_deg, eps_rad);
    let mc = angles.mc_deg;
    let asc = angles.asc_deg;
    let ic = (mc + 180.0).rem_euclid(360.0);
    let dsc = (asc + 180.0).rem_euclid(360.0);

    let mut warnings = Vec::new();
    let use_porphyry = if abs_lat >= 90.0 - obliquity_deg {
        match polar_policy {
            PolarPolicy::Strict => {
                return Err(EphemerisError::UndefinedHouseSystem {
                    house_system: 'G',
                    latitude_deg,
                    reason: "Gauquelin Placidus semi-arc undefined inside polar circle".to_owned(),
                });
            }
            PolarPolicy::FallbackToPorphyry => {
                warnings.push(HouseWarning::PolarFallbackToPorphyry);
                true
            }
            PolarPolicy::FallbackToWholeSign => {
                warnings.push(HouseWarning::PolarFallbackToWholeSign);
                true
            }
        }
    } else {
        false
    };

    let mut sectors = [0.0; 36];
    let q_diurnal = (asc - mc).rem_euclid(360.0);
    let q_nocturnal = (ic - asc).rem_euclid(360.0);

    sectors[0] = asc;
    sectors[9] = mc;
    sectors[18] = dsc;
    sectors[27] = ic;

    for idx in 1..9 {
        let f_pos = ((9 - idx) as f64) / 9.0;
        let c_rising = if use_porphyry {
            (mc + f_pos * q_diurnal).rem_euclid(360.0)
        } else {
            let d = solve_placidus_offset_deg(armc_deg, latitude_deg, eps_rad, f_pos)
                .unwrap_or(f_pos * q_diurnal);
            asc1(armc_deg + d, 0.0, eps_rad)
        };
        if let Some(slot) = sectors.get_mut(idx) {
            *slot = c_rising;
        }
        if let Some(slot) = sectors.get_mut(idx + 18) {
            *slot = (c_rising + 180.0).rem_euclid(360.0);
        }

        let f_neg = -(idx as f64) / 9.0;
        let c_setting = if use_porphyry {
            let frac_from_asc = ((9 - idx) as f64) / 9.0;
            (asc + frac_from_asc * q_nocturnal + 180.0).rem_euclid(360.0)
        } else {
            let d = solve_placidus_offset_deg(armc_deg, latitude_deg, eps_rad, f_neg)
                .unwrap_or(0.0);
            asc1(armc_deg + d, 0.0, eps_rad)
        };
        if let Some(slot) = sectors.get_mut(9 + idx) {
            *slot = c_setting;
        }
        if let Some(slot) = sectors.get_mut(27 + idx) {
            *slot = (c_setting + 180.0).rem_euclid(360.0);
        }
    }

    Ok(GauquelinResult {
        sectors_deg: sectors,
        warnings,
    })
}

/// Computes the 3D fractional mundane house position $H \in [1.0, 13.0)$ of a sky point
/// with ecliptic coordinates `(ecl_lon_deg, ecl_lat_deg)` at `(armc_deg, latitude_deg, obliquity_deg)`
/// for `system`, matching `swe_house_pos` (`Wave C`).
///
/// Accounts for non-zero ecliptic latitude $\beta$ in 3D mundane house systems (`Placidus`,
/// `Regiomontanus`, `Campanus`, `Horizontal`, `Meridian`, `Morinus`, `Koch`, `Alcabitius`),
/// and interpolates across the 12 house cusps for ecliptic-partitioned systems.
pub fn house_pos(
    armc_deg: f64,
    latitude_deg: f64,
    obliquity_deg: f64,
    system: HouseSystem,
    ecl_lon_deg: f64,
    ecl_lat_deg: f64,
) -> Result<f64, EphemerisError> {
    if !armc_deg.is_finite()
        || !latitude_deg.is_finite()
        || !(-90.0..=90.0).contains(&latitude_deg)
        || !obliquity_deg.is_finite()
        || obliquity_deg <= 0.0
        || obliquity_deg >= 90.0
        || !ecl_lon_deg.is_finite()
        || !ecl_lat_deg.is_finite()
        || !(-90.0..=90.0).contains(&ecl_lat_deg)
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid house_pos inputs: armc={armc_deg}, lat={latitude_deg}, obl={obliquity_deg}, lon={ecl_lon_deg}, lat_ecl={ecl_lat_deg}"
            ),
        });
    }

    let armc = armc_deg.rem_euclid(360.0);
    let lon = ecl_lon_deg.rem_euclid(360.0);
    let lat_rad = ecl_lat_deg.to_radians();
    let lon_rad = lon.to_radians();
    let eps_rad = obliquity_deg.to_radians();
    let phi_rad = latitude_deg.clamp(-89.999_999, 89.999_999).to_radians();

    // Ecliptic unit vector and equatorial unit vector
    let (sin_l, cos_l) = lon_rad.sin_cos();
    let (sin_b, cos_b) = lat_rad.sin_cos();
    let (sin_e, cos_e) = eps_rad.sin_cos();

    let ux = cos_b * cos_l;
    let uy = cos_b * sin_l * cos_e - sin_b * sin_e;
    let uz = cos_b * sin_l * sin_e + sin_b * cos_e;

    let alpha_deg = uy.atan2(ux).to_degrees().rem_euclid(360.0);
    let sin_delta = uz.clamp(-1.0, 1.0);
    let delta_rad = sin_delta.asin();

    // Hour angle H in [-180, 180): positive West of meridian, negative East of meridian
    let h_deg = (armc - alpha_deg + 180.0).rem_euclid(360.0) - 180.0;
    let h_rad = h_deg.to_radians();

    let wrap_house = |h: f64| -> f64 {
        let r = (h - 1.0).rem_euclid(12.0) + 1.0;
        if r >= 13.0 { 1.0 } else { r }
    };

    match system {
        HouseSystem::Placidus => {
            let tan_prod = delta_rad.tan() * phi_rad.tan();
            if tan_prod.abs() >= 1.0 {
                // Fall back to Porphyry interpolation when body is circumpolar in Placidus
                let res = houses_armc(
                    armc,
                    latitude_deg,
                    obliquity_deg,
                    0.0,
                    HouseOptions {
                        system: HouseSystem::Porphyry,
                        polar_policy: PolarPolicy::FallbackToPorphyry,
                        ..HouseOptions::default()
                    },
                )?;
                return Ok(interp_cusps_lon(lon, &res.cusps_deg));
            }
            let ad_deg = tan_prod.clamp(-1.0, 1.0).asin().to_degrees();
            let sda = (90.0 + ad_deg).max(1.0e-9);
            let sna = (90.0 - ad_deg).max(1.0e-9);
            let pos = if h_deg.abs() <= sda {
                if h_deg <= 0.0 {
                    10.0 + 3.0 * (-h_deg / sda)
                } else {
                    10.0 - 3.0 * (h_deg / sda)
                }
            } else if h_deg <= 0.0 {
                1.0 + 3.0 * ((-h_deg - sda) / sna)
            } else {
                7.0 - 3.0 * ((h_deg - sda) / sna)
            };
            Ok(wrap_house(pos))
        }
        HouseSystem::PolichPage => {
            let tan_prod = delta_rad.tan() * phi_rad.tan();
            let hm_deg = if h_deg.abs() < 1.0e-13 {
                0.0
            } else if (h_deg.abs() - 180.0).abs() < 1.0e-13 {
                180.0
            } else {
                let above_horizon = (h_rad.cos() + tan_prod) >= 0.0;
                let (mut lo, mut hi): (f64, f64) = if h_deg >= 0.0 {
                    if above_horizon {
                        (0.0, 90.0)
                    } else {
                        (90.0, 180.0)
                    }
                } else if above_horizon {
                    (-90.0, 0.0)
                } else {
                    (-180.0, -90.0)
                };
                let mut hm: f64 = 0.5 * (lo + hi);
                let dt = if above_horizon {
                    -1.0 / 90.0
                } else {
                    1.0 / 90.0
                };
                for _ in 0..60 {
                    let tau = if hm.abs() <= 90.0 {
                        -hm / 90.0
                    } else if hm > 90.0 {
                        (hm - 180.0) / 90.0
                    } else {
                        (hm + 180.0) / 90.0
                    };
                    let diff_rad = (hm - h_deg).to_radians();
                    let g = diff_rad.sin() - tau * tan_prod;
                    if g.abs() < 1.0e-15 || (hi - lo) < 1.0e-13 {
                        break;
                    }
                    if g < 0.0 {
                        lo = hm;
                    } else {
                        hi = hm;
                    }
                    let dg = diff_rad.cos() * (std::f64::consts::PI / 180.0) - dt * tan_prod;
                    let cand = if dg.abs() > 1.0e-9 {
                        hm - g / dg
                    } else {
                        0.5 * (lo + hi)
                    };
                    hm = if cand > lo && cand < hi {
                        cand
                    } else {
                        0.5 * (lo + hi)
                    };
                }
                hm
            };
            Ok(wrap_house(10.0 - hm_deg / 30.0))
        }
        HouseSystem::Regiomontanus => {
            let qr_deg = h_rad
                .sin()
                .atan2(h_rad.cos() + delta_rad.tan() * phi_rad.tan())
                .to_degrees();
            let pos = 10.0 - qr_deg / 30.0;
            Ok(wrap_house(pos))
        }
        HouseSystem::Campanus => {
            let r_rad = armc.to_radians();
            let east = [-r_rad.sin(), r_rad.cos(), 0.0];
            let zenith = [
                phi_rad.cos() * r_rad.cos(),
                phi_rad.cos() * r_rad.sin(),
                phi_rad.sin(),
            ];
            let u_eq = [ux, uy, uz];
            let ue = crate::frame::vec_dot(&u_eq, &east);
            let uz_hz = crate::frame::vec_dot(&u_eq, &zenith);
            let theta_deg = (-uz_hz).atan2(ue).to_degrees().rem_euclid(360.0);
            Ok(wrap_house(1.0 + theta_deg / 30.0))
        }
        HouseSystem::Azimuthal => {
            let r_rad = armc.to_radians();
            let east = [-r_rad.sin(), r_rad.cos(), 0.0];
            let zenith = [
                phi_rad.cos() * r_rad.cos(),
                phi_rad.cos() * r_rad.sin(),
                phi_rad.sin(),
            ];
            let north = crate::frame::vec_cross(&zenith, &east);
            let u_eq = [ux, uy, uz];
            let ue = crate::frame::vec_dot(&u_eq, &east);
            let un = crate::frame::vec_dot(&u_eq, &north);
            let theta_deg = un.atan2(ue).to_degrees().rem_euclid(360.0);
            Ok(wrap_house(1.0 + theta_deg / 30.0))
        }
        HouseSystem::Meridian => {
            let pos = 10.0 - h_deg / 30.0;
            Ok(wrap_house(pos))
        }
        HouseSystem::Morinus => {
            let alpha_mor_deg = (uy * cos_e + uz * sin_e)
                .atan2(ux * cos_e)
                .to_degrees()
                .rem_euclid(360.0);
            let diff = (alpha_mor_deg - armc - 90.0).rem_euclid(360.0);
            Ok(wrap_house(1.0 + diff / 30.0))
        }
        HouseSystem::Alcabitius => {
            let angles = compute_angles(armc, latitude_deg, eps_rad);
            let sin_dec_asc = angles.asc_deg.to_radians().sin() * sin_e;
            let cos_dec_asc = (1.0 - sin_dec_asc * sin_dec_asc).max(1.0e-15).sqrt();
            let tan_dec_asc = sin_dec_asc / cos_dec_asc;
            let ad_asc = (phi_rad.tan() * tan_dec_asc)
                .clamp(-1.0, 1.0)
                .asin()
                .to_degrees();
            let sday = (90.0 + ad_asc).max(1.0e-9);
            let snight = (90.0 - ad_asc).max(1.0e-9);
            let a_east = (alpha_deg - armc).rem_euclid(360.0);
            let pos = if a_east < sday {
                10.0 + 3.0 * (a_east / sday)
            } else if a_east < 180.0 {
                1.0 + 3.0 * ((a_east - sday) / snight)
            } else if a_east < 180.0 + sday {
                4.0 + 3.0 * ((a_east - 180.0) / sday)
            } else {
                7.0 + 3.0 * ((a_east - 180.0 - sday) / snight)
            };
            Ok(wrap_house(pos))
        }
        HouseSystem::Koch => {
            let angles = compute_angles(armc, latitude_deg, eps_rad);
            let sin_dec_mc = angles.mc_deg.to_radians().sin() * sin_e;
            let cos_dec_mc = (1.0 - sin_dec_mc * sin_dec_mc).max(1.0e-15).sqrt();
            let tan_dec_mc = sin_dec_mc / cos_dec_mc;
            let x_ad_mc = phi_rad.tan() * tan_dec_mc;
            let x_ad_body = delta_rad.tan() * phi_rad.tan();
            if x_ad_mc.abs() >= 1.0
                || x_ad_body.abs() >= 1.0
                || latitude_deg.abs() >= 90.0 - obliquity_deg
            {
                let res = houses_armc(
                    armc,
                    latitude_deg,
                    obliquity_deg,
                    0.0,
                    HouseOptions {
                        system: HouseSystem::Porphyry,
                        polar_policy: PolarPolicy::FallbackToPorphyry,
                        ..HouseOptions::default()
                    },
                )?;
                return Ok(interp_cusps_lon(lon, &res.cusps_deg));
            }
            let ad_mc = x_ad_mc.asin().to_degrees();
            let ad_body = x_ad_body.asin().to_degrees();
            let sda_mc = (90.0 + ad_mc).max(1.0e-9);
            let pos = if h_deg <= 0.0 {
                10.0 + 3.0 * ((-h_deg - ad_body + ad_mc) / sda_mc)
            } else {
                4.0 + 3.0 * (((180.0 - h_deg) + ad_body + ad_mc) / sda_mc)
            };
            Ok(wrap_house(pos))
        }
        _ => {
            let res = houses_armc(
                armc,
                latitude_deg,
                obliquity_deg,
                0.0,
                HouseOptions {
                    system,
                    polar_policy: PolarPolicy::FallbackToPorphyry,
                    ..HouseOptions::default()
                },
            )?;
            Ok(interp_cusps_lon(lon, &res.cusps_deg))
        }
    }
}

/// Interpolates fractional house position `[1.0, 13.0)` of `lon_deg` across a 12-cusp ring.
pub(crate) fn interp_cusps_lon(lon_deg: f64, cusps_deg: &[f64; 12]) -> f64 {
    let lon = lon_deg.rem_euclid(360.0);
    let idx = (house_of_lon(lon, cusps_deg, 1.0e-9) as usize)
        .saturating_sub(1)
        .min(11);
    let c_start = cusps_deg.get(idx).copied().unwrap_or(0.0);
    let c_end = cusps_deg.get((idx + 1) % 12).copied().unwrap_or(30.0);
    let span = (c_end - c_start).rem_euclid(360.0).max(1.0e-9);
    let offset = (lon - c_start).rem_euclid(360.0);
    let frac = (offset / span).clamp(0.0, 0.999_999_999_999);
    (idx as f64) + 1.0 + frac
}
