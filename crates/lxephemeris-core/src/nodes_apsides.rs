//! Planetary and lunar nodes and apsides (`swe_nod_aps` / `swe_nod_aps_ut`).
//!
//! Provides clean-room celestial mechanics for calculating:
//! - Ascending and descending orbital nodes in the reference ecliptic plane
//! - Periapsis (perihelion / perigee) and apoapsis (aphelion / apogee)
//! - Second focal point of elliptical orbits (`SE_NODBIT_FOPOINT`)
//! - Mean lunar and planetary orbits via Simon et al. (1994) / Standish (1995) polynomials (`SE_NODBIT_MEAN`)
//! - Instantaneous osculating Keplerian orbits around the primary (`SE_NODBIT_OSCU`) or SSB (`SE_NODBIT_OSCU_BAR`)
//! - Observer shifts (geocentric, topocentric, heliocentric, barycentric)
//! - 4th-order 4-point symmetric stencil velocities ($h = 0.25\text{ d}$ for planets, $h = 0.0416667\text{ d}$ for Moon)
//! - Equinox of date / J2000, true nutation, and sidereal ayanamsa frame rotations.

use crate::astrometry::{equatorial_to_horizontal_vector_au, AtmosphereParams};
use crate::ayanamsa::eval_ayanamsa;
use crate::error::EphemerisError;
use crate::frame::{
    cart_state_to_spherical, cart_to_spherical, mat_transpose, mat_vec, spherical_state_to_cart,
    stencil_4pt_derivative, vec_cross, vec_dot, vec_norm, SphericalState, Vec3,
};
use crate::moment::EpochContext;
use crate::spk::SpkSource;
use crate::time::{TimeInput, TimeTable};
use crate::types::{
    BodyId, CalcOptions, Equinox, Observer, Origin, Plane, PlanetId, Zodiac, AU_KM, JD_J2000,
};
use serde::{Deserialize, Serialize};
use sofars::fundargs::{faf03, fal03, faom03};

// === Swiss Ephemeris `swe_nod_aps` bitmask constants ===

/// Mean node and apsis elements (Simon et al. 1994 / Standish 1995).
pub const SE_NODBIT_MEAN: u32 = 1;

/// Osculating node and apsis elements around the central body.
pub const SE_NODBIT_OSCU: u32 = 2;

/// Osculating node and apsis elements around the Solar System Barycenter (SSB).
pub const SE_NODBIT_OSCU_BAR: u32 = 4;

/// Returns the second focal point of the ellipse instead of the apoapsis point.
pub const SE_NODBIT_FOPOINT: u32 = 256;

// === Stencil step sizes for velocities ===

/// Stencil step size $h$ for planetary and asteroid nodes/apsides ($0.25\text{ d} = 6\text{ h}$).
pub const H_NODES_APSIDES_PLANET_DAYS: f64 = 0.25;

/// Stencil step size $h$ for lunar nodes and apsides ($0.0416667\text{ d} = 1\text{ h}$).
pub const H_NODES_APSIDES_MOON_DAYS: f64 = 0.0416667;

/// Stencil step size $h$ for horizontal plane nodes/apsides ($5\times 10^{-4}\text{ d} = 43.2\text{ s}$).
/// Reduces diurnal $\sin(\omega h)/(\omega h)$ rotation error to $< 1.7\text{ ppm}$.
pub const H_NODES_APSIDES_HORIZONTAL_DAYS: f64 = 5.0e-4;

// === Standard Gravitational Parameters (DE440 in km^3 / s^2) ===

/// Sun gravitational parameter $GM_\odot$ from DE440 (Park et al. 2021).
pub const GM_SUN: f64 = 132_712_440_041.939_38;

/// Mercury gravitational parameter $GM$ from DE440.
pub const GM_MERCURY: f64 = 22_031.780_000;

/// Venus gravitational parameter $GM$ from DE440.
pub const GM_VENUS: f64 = 324_858.592_000;

/// Earth gravitational parameter $GM_\oplus$ from DE440.
pub const GM_EARTH: f64 = 398_600.435_436;

/// Moon gravitational parameter $GM_{\text{Moon}}$ from DE440.
pub const GM_MOON: f64 = 4_902.800_066;

/// Earth-Moon combined gravitational parameter $\mu_{\text{EM}} = GM_\oplus + GM_{\text{Moon}}$.
pub const GM_EM: f64 = GM_EARTH + GM_MOON; // 403_503.235_502

/// Mars system gravitational parameter $GM$ from DE440.
pub const GM_MARS: f64 = 42_828.375_214;

/// Jupiter system gravitational parameter $GM$ from DE440.
pub const GM_JUPITER: f64 = 126_712_764.800;

/// Saturn system gravitational parameter $GM$ from DE440.
pub const GM_SATURN: f64 = 37_940_585.200;

/// Uranus system gravitational parameter $GM$ from DE440.
pub const GM_URANUS: f64 = 5_794_548.600;

/// Neptune system gravitational parameter $GM$ from DE440.
pub const GM_NEPTUNE: f64 = 6_836_527.100;

/// Pluto system gravitational parameter $GM$ from DE440.
pub const GM_PLUTO: f64 = 977.000;

/// Sum of all major planetary gravitational parameters in DE440.
pub const GM_PLANETS_SUM: f64 = GM_MERCURY
    + GM_VENUS
    + GM_EM
    + GM_MARS
    + GM_JUPITER
    + GM_SATURN
    + GM_URANUS
    + GM_NEPTUNE
    + GM_PLUTO;

/// Total Solar System gravitational parameter $\mu_{\text{SSB}} = GM_\odot + \sum GM_{\text{planets}}$.
pub const GM_SSB: f64 = GM_SUN + GM_PLANETS_SUM; // 132_890_518_666.629_38

// === Output Types ===

/// Result of a planetary or lunar nodes and apsides calculation (`swe_nod_aps`).
///
/// Each 6-element array contains:
/// - `[lon_deg, lat_deg, dist_au, lon_speed_deg_day, lat_speed_deg_day, dist_speed_au_day]`
///   in spherical coordinates (or Right Ascension / Declination when `Plane::Equatorial`),
/// - or Cartesian `[x, y, z, vx, vy, vz]` (in AU and AU/day) when formatted via `format_swe_state6`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NodesApsidesResult {
    /// Ascending node vector: `[lon, lat, dist, lon_speed, lat_speed, dist_speed]`.
    pub asc_node: [f64; 6],
    /// Descending node vector: `[lon, lat, dist, lon_speed, lat_speed, dist_speed]`.
    pub desc_node: [f64; 6],
    /// Periapsis (perihelion / perigee) vector: `[lon, lat, dist, lon_speed, lat_speed, dist_speed]`.
    pub periapsis: [f64; 6],
    /// Apoapsis (aphelion / apogee) vector (or second focal point if `SE_NODBIT_FOPOINT` is set).
    pub apoapsis: [f64; 6],
}

impl NodesApsidesResult {
    /// Formats each of the 4 points with Swiss Ephemeris `SEFLG_*` bitflags:
    /// - `flags & SEFLG_XYZ != 0`: outputs Cartesian `[x, y, z, vx, vy, vz]` in AU and AU/day.
    /// - `flags & SEFLG_RADIANS != 0`: outputs `[lon_rad, lat_rad, dist_au, lon_speed_rad_day, lat_speed_rad_day, dist_speed_au_day]`.
    #[must_use]
    pub fn format_with_flags(&self, flags: u32) -> Self {
        const SEFLG_XYZ: u32 = 4096;
        const SEFLG_RADIANS: u32 = 8192;

        if flags & SEFLG_XYZ != 0 {
            let include_speed = (flags & 256) != 0;
            let fmt_xyz = |s: [f64; 6]| -> [f64; 6] {
                let sph = SphericalState {
                    lon_deg: s[0],
                    lat_deg: s[1],
                    dist_au: s[2],
                    lon_speed_deg_day: s[3],
                    lat_speed_deg_day: s[4],
                    dist_speed_au_day: s[5],
                };
                let (pos, vel) = spherical_state_to_cart(&sph, include_speed);
                [pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]]
            };
            Self {
                asc_node: fmt_xyz(self.asc_node),
                desc_node: fmt_xyz(self.desc_node),
                periapsis: fmt_xyz(self.periapsis),
                apoapsis: fmt_xyz(self.apoapsis),
            }
        } else if flags & SEFLG_RADIANS != 0 {
            let fmt_rad = |s: [f64; 6]| -> [f64; 6] {
                [
                    s[0].to_radians(),
                    s[1].to_radians(),
                    s[2],
                    s[3].to_radians(),
                    s[4].to_radians(),
                    s[5],
                ]
            };
            Self {
                asc_node: fmt_rad(self.asc_node),
                desc_node: fmt_rad(self.desc_node),
                periapsis: fmt_rad(self.periapsis),
                apoapsis: fmt_rad(self.apoapsis),
            }
        } else {
            self.clone()
        }
    }
}

// === Standish (1995) / Simon et al. (1994) Mean Orbital Elements ===

#[derive(Debug, Clone, Copy)]
struct StandishPlanetElements {
    a0: f64,
    a_rate: f64,
    e0: f64,
    e_rate: f64,
    i0_deg: f64,
    i_rate_deg: f64,
    om0_deg: f64,
    om_rate_deg: f64,
    w0_deg: f64,
    w_rate_deg: f64,
}

fn standish_planet_elements(p: PlanetId) -> Option<StandishPlanetElements> {
    match p {
        PlanetId::Mercury => Some(StandishPlanetElements {
            a0: 0.38709927,
            a_rate: 0.00000037,
            e0: 0.20563593,
            e_rate: 0.00002527,
            i0_deg: 7.00497902,
            i_rate_deg: -0.00594749,
            om0_deg: 48.33076593,
            om_rate_deg: -0.12534081,
            w0_deg: 77.45779628,
            w_rate_deg: 0.16047689,
        }),
        PlanetId::Venus => Some(StandishPlanetElements {
            a0: 0.72333566,
            a_rate: 0.00000026,
            e0: 0.00677672,
            e_rate: -0.00004107,
            i0_deg: 3.39467605,
            i_rate_deg: -0.00078890,
            om0_deg: 76.67984255,
            om_rate_deg: -0.27769418,
            w0_deg: 131.60246718,
            w_rate_deg: 0.00268329,
        }),
        PlanetId::Earth | PlanetId::EarthMoonBarycenter | PlanetId::Sun => {
            Some(StandishPlanetElements {
                a0: 1.00000261,
                a_rate: 0.00000562,
                e0: 0.01671123,
                e_rate: -0.00004392,
                i0_deg: 0.00000000,
                i_rate_deg: 0.00000000,
                om0_deg: 0.00000000,
                om_rate_deg: 0.00000000,
                w0_deg: 102.93768193,
                w_rate_deg: 0.32327364,
            })
        }
        PlanetId::Mars => Some(StandishPlanetElements {
            a0: 1.52371034,
            a_rate: 0.00001847,
            e0: 0.09339410,
            e_rate: 0.00007882,
            i0_deg: 1.84969142,
            i_rate_deg: -0.00813131,
            om0_deg: 49.55953891,
            om_rate_deg: -0.29257343,
            w0_deg: 336.05637041,
            w_rate_deg: 0.44441088,
        }),
        PlanetId::Jupiter => Some(StandishPlanetElements {
            a0: 5.20288700,
            a_rate: -0.00011607,
            e0: 0.04838624,
            e_rate: -0.00013253,
            i0_deg: 1.30439695,
            i_rate_deg: -0.00183714,
            om0_deg: 100.47390909,
            om_rate_deg: 0.20469106,
            w0_deg: 14.72847983,
            w_rate_deg: 0.21252668,
        }),
        PlanetId::Saturn => Some(StandishPlanetElements {
            a0: 9.53667594,
            a_rate: -0.00125060,
            e0: 0.05386179,
            e_rate: -0.00050991,
            i0_deg: 2.48599187,
            i_rate_deg: 0.00193609,
            om0_deg: 113.66242448,
            om_rate_deg: -0.28867794,
            w0_deg: 92.59887831,
            w_rate_deg: -0.41897216,
        }),
        PlanetId::Uranus => Some(StandishPlanetElements {
            a0: 19.18916464,
            a_rate: -0.00196176,
            e0: 0.04725744,
            e_rate: -0.00004397,
            i0_deg: 0.77263783,
            i_rate_deg: -0.00242939,
            om0_deg: 74.01692503,
            om_rate_deg: 0.04240589,
            w0_deg: 170.95427630,
            w_rate_deg: 0.40805281,
        }),
        PlanetId::Neptune => Some(StandishPlanetElements {
            a0: 30.06992276,
            a_rate: 0.00026291,
            e0: 0.00859048,
            e_rate: 0.00005105,
            i0_deg: 1.77004347,
            i_rate_deg: 0.00035372,
            om0_deg: 131.78422574,
            om_rate_deg: -0.00598750,
            w0_deg: 44.96476227,
            w_rate_deg: -0.32241464,
        }),
        PlanetId::Pluto | PlanetId::Moon | PlanetId::SolarSystemBarycenter => None,
    }
}

fn planet_gm(p: PlanetId) -> f64 {
    match p {
        PlanetId::Mercury => GM_MERCURY,
        PlanetId::Venus => GM_VENUS,
        PlanetId::Earth => GM_EM,
        PlanetId::Mars => GM_MARS,
        PlanetId::Jupiter => GM_JUPITER,
        PlanetId::Saturn => GM_SATURN,
        PlanetId::Uranus => GM_URANUS,
        PlanetId::Neptune => GM_NEPTUNE,
        PlanetId::Pluto => GM_PLUTO,
        PlanetId::Moon => GM_MOON,
        PlanetId::Sun => GM_SUN,
        PlanetId::EarthMoonBarycenter => GM_EM,
        PlanetId::SolarSystemBarycenter => GM_SSB,
    }
}

fn resolve_target_and_center_id(
    body: BodyId,
    method_bits: u32,
) -> Result<(i32, i32, f64), EphemerisError> {
    match body {
        BodyId::Planet(p) => match p {
            PlanetId::Moon => Ok((301, 399, GM_EM)),
            PlanetId::Earth | PlanetId::EarthMoonBarycenter | PlanetId::Sun => {
                Ok((3, 10, GM_SUN + GM_EM))
            }
            PlanetId::SolarSystemBarycenter => Err(EphemerisError::UnsupportedBody {
                body: "PlanetId::SolarSystemBarycenter has no orbital nodes".to_owned(),
            }),
            other_p => {
                let is_outer_bary = matches!(
                    other_p,
                    PlanetId::Saturn | PlanetId::Uranus | PlanetId::Neptune | PlanetId::Pluto
                );
                let use_bary = (method_bits & SE_NODBIT_OSCU_BAR) != 0 && is_outer_bary;
                let center_id = if use_bary { 0 } else { 10 };
                let target_id = other_p.naif_target_id();
                let gm = if use_bary {
                    GM_SSB
                } else {
                    GM_SUN + planet_gm(other_p)
                };
                Ok((target_id, center_id, gm))
            }
        },
        BodyId::MinorPlanet(mpc) => {
            let center_id = if method_bits & SE_NODBIT_OSCU_BAR != 0 {
                0
            } else {
                10
            };
            let target_id = mpc.standard_naif_spk_id();
            let gm = if method_bits & SE_NODBIT_OSCU_BAR != 0 {
                GM_SSB
            } else {
                GM_SUN
            };
            Ok((target_id, center_id, gm))
        }
        other => Err(EphemerisError::UnsupportedBody {
            body: format!("{other:?}"),
        }),
    }
}

fn shift_time_input(input: TimeInput, delta_days: f64) -> TimeInput {
    match input {
        TimeInput::Tt(s) => TimeInput::Tt(s.add_days(delta_days)),
        TimeInput::Ut1(s) => TimeInput::Ut1(s.add_days(delta_days)),
        TimeInput::Utc(s) => TimeInput::Utc(s.add_days(delta_days)),
        TimeInput::Tdb(s) => TimeInput::Tdb(s.add_days(delta_days)),
    }
}

/// Evaluates the 4 orbital points `[asc_node, desc_node, periapsis, apoapsis]`
/// (in AU, in the requested output frame) at a single instant `t_input`.
fn eval_points_at_time(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    body: BodyId,
    t_input: TimeInput,
    options: &CalcOptions,
    method_bits: u32,
    observer: Option<Observer>,
) -> Result<([Vec3; 4], EpochContext), EphemerisError> {
    let ctx = EpochContext::new(spk, time_table, t_input, observer)?;
    let is_mean = (method_bits & SE_NODBIT_MEAN) != 0
        || (method_bits & (SE_NODBIT_OSCU | SE_NODBIT_OSCU_BAR)) == 0;
    let eff_method_bits = if is_mean {
        method_bits & !SE_NODBIT_OSCU_BAR
    } else {
        method_bits
    };
    let (target_id, center_id, mu) = resolve_target_and_center_id(body, eff_method_bits)?;

    // Reference ecliptic matrix: J2000 Ecliptic if requested, else Mean Ecliptic of Date.
    let r_ecl = if options.equinox == Equinox::J2000 {
        ctx.frame.r_ecl_j2000
    } else {
        ctx.frame.r_ecl_mean
    };

    // Compute the 4 points relative to the orbital center (in km, in the r_ecl frame).
    let (mut pts_center_km, used_ecl_mat) = if is_mean && body == BodyId::Planet(PlanetId::Moon) {
        // Simon et al. (1994) mean lunar orbit
        let big_t = (ctx.tt.jd_f64() - JD_J2000) / 36_525.0;
        let om_rad = faom03(big_t);
        let f_rad = faf03(big_t);
        let l_rad = fal03(big_t);
        let varpi_rad = f_rad + om_rad - l_rad;
        let omega_rad = varpi_rad - om_rad; // f_rad - l_rad

        let a = 384_400.0_f64;
        let e = 0.0549006_f64;
        let inc = 5.145396_f64.to_radians();
        let p = a * (1.0 - e * e);

        let (sin_om, cos_om) = om_rad.sin_cos();
        let (sin_w, cos_w) = omega_rad.sin_cos();
        let (sin_i, cos_i) = inc.sin_cos();

        let r_asc = p / (1.0 + e * cos_w);
        let r_desc = p / (1.0 - e * cos_w);
        let r_peri = a * (1.0 - e);

        let pos_asc = [r_asc * cos_om, r_asc * sin_om, 0.0];
        let pos_desc = [-r_desc * cos_om, -r_desc * sin_om, 0.0];

        let xp = r_peri * cos_w;
        let yp = r_peri * sin_w;
        let pos_peri = [
            cos_om * xp - sin_om * yp * cos_i,
            sin_om * xp + cos_om * yp * cos_i,
            yp * sin_i,
        ];
        let p_hat = [
            pos_peri[0] / r_peri,
            pos_peri[1] / r_peri,
            pos_peri[2] / r_peri,
        ];

        let pos_apo = if method_bits & SE_NODBIT_FOPOINT != 0 {
            let d_foc = -2.0 * a * e;
            [d_foc * p_hat[0], d_foc * p_hat[1], d_foc * p_hat[2]]
        } else {
            let r_apo = -a * (1.0 + e);
            [r_apo * p_hat[0], r_apo * p_hat[1], r_apo * p_hat[2]]
        };

        // Lunar mean elements are in Mean Ecliptic of Date.
        ([pos_asc, pos_desc, pos_peri, pos_apo], ctx.frame.r_ecl_mean)
    } else if is_mean
        && matches!(body, BodyId::Planet(p) if standish_planet_elements(p).is_some())
    {
        // Standish (1995) / Simon et al. (1994) mean planetary orbit
        let p_id = match body {
            BodyId::Planet(p) => p,
            _ => PlanetId::Mercury,
        };
        let el = standish_planet_elements(p_id).ok_or_else(|| EphemerisError::UnsupportedBody {
            body: format!("{body:?}"),
        })?;
        let big_t = (ctx.tt.jd_f64() - JD_J2000) / 36_525.0;

        let a_au = el.a0 + el.a_rate * big_t;
        let e = (el.e0 + el.e_rate * big_t).clamp(0.0, 0.999_999);
        let inc = (el.i0_deg + el.i_rate_deg * big_t).to_radians();
        let om_rad = (el.om0_deg + el.om_rate_deg * big_t).to_radians();
        let varpi_rad = (el.w0_deg + el.w_rate_deg * big_t).to_radians();
        let omega_rad = varpi_rad - om_rad;

        let a_km = a_au * AU_KM;
        let p_km = a_km * (1.0 - e * e);

        let (sin_om, cos_om) = om_rad.sin_cos();
        let (sin_w, cos_w) = omega_rad.sin_cos();
        let (sin_i, cos_i) = inc.sin_cos();

        let r_asc = p_km / (1.0 + e * cos_w);
        let r_desc = p_km / (1.0 - e * cos_w);
        let r_peri = a_km * (1.0 - e);

        let pos_asc = [r_asc * cos_om, r_asc * sin_om, 0.0];
        let pos_desc = [-r_desc * cos_om, -r_desc * sin_om, 0.0];

        let xp = r_peri * cos_w;
        let yp = r_peri * sin_w;
        let pos_peri = [
            cos_om * xp - sin_om * yp * cos_i,
            sin_om * xp + cos_om * yp * cos_i,
            yp * sin_i,
        ];
        let p_hat = [
            pos_peri[0] / r_peri,
            pos_peri[1] / r_peri,
            pos_peri[2] / r_peri,
        ];

        let pos_apo = if method_bits & SE_NODBIT_FOPOINT != 0 {
            let d_foc = -2.0 * a_km * e;
            [d_foc * p_hat[0], d_foc * p_hat[1], d_foc * p_hat[2]]
        } else {
            let r_apo = -a_km * (1.0 + e);
            [r_apo * p_hat[0], r_apo * p_hat[1], r_apo * p_hat[2]]
        };

        // Standish elements are in J2000 Ecliptic.
        ([pos_asc, pos_desc, pos_peri, pos_apo], ctx.frame.r_ecl_j2000)
    } else {
        // Instantaneous osculating Keplerian orbit from SPK state
        let st = spk.state_icrf_km_split(target_id, center_id, ctx.tdb.0.jd1, ctx.tdb.0.jd2)?;
        let r_vec = mat_vec(&r_ecl, &st.pos_km);
        let v_vec = mat_vec(&r_ecl, &st.vel_km_s);

        let r = vec_norm(&r_vec);
        let v2 = vec_dot(&v_vec, &v_vec);
        let h_vec = vec_cross(&r_vec, &v_vec);
        let h = vec_norm(&h_vec);

        if r < 1.0e-12 || h < 1.0e-12 {
            return Err(EphemerisError::InvalidInput {
                message: "Degenerate orbit (zero radius or angular momentum)".to_owned(),
            });
        }

        let v_cross_h = vec_cross(&v_vec, &h_vec);
        let inv_mu = 1.0 / mu;
        let inv_r = 1.0 / r;
        let e_vec = [
            v_cross_h[0] * inv_mu - r_vec[0] * inv_r,
            v_cross_h[1] * inv_mu - r_vec[1] * inv_r,
            v_cross_h[2] * inv_mu - r_vec[2] * inv_r,
        ];
        let e = vec_norm(&e_vec);

        let inv_a = (2.0 / r) - (v2 / mu);
        if inv_a <= 1.0e-15 || e >= 1.0 {
            return Err(EphemerisError::InvalidInput {
                message: "Unbound parabolic or hyperbolic orbit is unsupported for nodes/apsides"
                    .to_owned(),
            });
        }
        let a = 1.0 / inv_a;
        let p = (h * h) / mu;

        // Line of nodes in the ecliptic plane: k x h = (-h_y, h_x, 0)
        let n_vec = [-h_vec[1], h_vec[0], 0.0];
        let n_norm = vec_norm(&n_vec);
        let n_hat = if n_norm < 1.0e-12 * h {
            [1.0, 0.0, 0.0]
        } else {
            [n_vec[0] / n_norm, n_vec[1] / n_norm, 0.0]
        };

        let e_dot_n = vec_dot(&e_vec, &n_hat);
        let denom_asc = (1.0 + e_dot_n).max(1.0e-12);
        let denom_desc = (1.0 - e_dot_n).max(1.0e-12);

        let r_asc = p / denom_asc;
        let r_desc = p / denom_desc;

        let pos_asc = [r_asc * n_hat[0], r_asc * n_hat[1], 0.0];
        let pos_desc = [-r_desc * n_hat[0], -r_desc * n_hat[1], 0.0];

        let p_hat = if e < 1.0e-12 {
            n_hat
        } else {
            [e_vec[0] / e, e_vec[1] / e, e_vec[2] / e]
        };

        let r_peri = a * (1.0 - e);
        let pos_peri = [r_peri * p_hat[0], r_peri * p_hat[1], r_peri * p_hat[2]];

        let pos_apo = if method_bits & SE_NODBIT_FOPOINT != 0 {
            let d_foc = -2.0 * a * e;
            [d_foc * p_hat[0], d_foc * p_hat[1], d_foc * p_hat[2]]
        } else {
            let r_apo = -a * (1.0 + e);
            [r_apo * p_hat[0], r_apo * p_hat[1], r_apo * p_hat[2]]
        };

        ([pos_asc, pos_desc, pos_peri, pos_apo], r_ecl)
    };

    if body == BodyId::Planet(PlanetId::Sun)
        && matches!(options.origin, Origin::Geocentric | Origin::Topocentric(_))
    {
        for pt in &mut pts_center_km {
            pt[0] = -pt[0];
            pt[1] = -pt[1];
            pt[2] = -pt[2];
        }
    }

    // Vector from Observer to Orbital Center in ICRF (km): r_{Obs -> C} = r_{C, bary} - r_{Obs, bary}
    let obs_to_c_icrf_km = if body == BodyId::Planet(PlanetId::Sun) {
        [0.0, 0.0, 0.0]
    } else {
        let p_center_bary = match center_id {
            0 => [0.0, 0.0, 0.0],
            10 => ctx.sun_bary_km.pos_km,
            399 => ctx.earth_bary_km.pos_km,
            _ => [0.0, 0.0, 0.0],
        };
        let p_obs_bary = match options.origin {
            Origin::Geocentric => ctx.earth_bary_km.pos_km,
            Origin::Topocentric(_) => ctx.observer_bary_km.pos_km,
            Origin::Heliocentric => ctx.sun_bary_km.pos_km,
            Origin::Barycentric => [0.0, 0.0, 0.0],
        };
        [
            p_center_bary[0] - p_obs_bary[0],
            p_center_bary[1] - p_obs_bary[1],
            p_center_bary[2] - p_obs_bary[2],
        ]
    };

    // Rotate observer shift vector into the used ecliptic frame
    let obs_to_c_ecl_km = mat_vec(&used_ecl_mat, &obs_to_c_icrf_km);
    let used_ecl_mat_t = mat_transpose(&used_ecl_mat);
    let r_out = ctx
        .frame
        .icrf_to_frame_matrix(options.plane, options.equinox, options.nutation);

    let mut out_points = [[0.0_f64; 3]; 4];
    for (i, pt_km) in pts_center_km.iter().enumerate() {
        // Observer-relative position in the ecliptic frame (in AU)
        let pt_obs_ecl_au = [
            (pt_km[0] + obs_to_c_ecl_km[0]) / AU_KM,
            (pt_km[1] + obs_to_c_ecl_km[1]) / AU_KM,
            (pt_km[2] + obs_to_c_ecl_km[2]) / AU_KM,
        ];
        // Rotate back to ICRF via R_ecl^T, then into output frame via R_out
        let pt_icrf_au = mat_vec(&used_ecl_mat_t, &pt_obs_ecl_au);
        let mut pt_out_au = mat_vec(&r_out, &pt_icrf_au);

        if options.plane == Plane::Horizontal {
            let obs = match options.origin {
                Origin::Topocentric(o) => o,
                _ => ctx.observer.ok_or_else(|| EphemerisError::InvalidInput {
                    message: "Plane::Horizontal requires a topocentric Observer".to_owned(),
                })?,
            };
            pt_out_au = equatorial_to_horizontal_vector_au(
                obs,
                ctx.gast_rad,
                &pt_out_au,
                AtmosphereParams::default(),
            );
        }

        if let Some(dest) = out_points.get_mut(i) {
            *dest = pt_out_au;
        }
    }

    Ok((out_points, ctx))
}

/// Computes the planetary or lunar orbital nodes and apsides (`swe_nod_aps` / `swe_nod_aps_ut`).
///
/// # Arguments
/// - `spk`: Abstract binary SPK source providing planetary ephemerides (`DE440s` / minor bodies).
/// - `time_table`: Source for $\Delta T$, `UT1-UTC`, and leap seconds.
/// - `body`: Celestial target body (`Planet` or `MinorPlanet`).
/// - `input`: Input time epoch in any supported time scale (`UT1`, `TT`, `UTC`, `TDB`).
/// - `options`: Calculation options (`Origin`, `Plane`, `Equinox`, `NutationMode`, `Zodiac`, `include_speed`).
/// - `method_bits`: Method selection bitmask:
///   - `SE_NODBIT_OSCU` (default, or when 0): instantaneous osculating Keplerian orbit.
///   - `SE_NODBIT_MEAN`: Simon et al. (1994) / Standish (1995) mean orbit.
///   - `SE_NODBIT_OSCU_BAR`: instantaneous osculating orbit around Solar System Barycenter (SSB).
///   - `SE_NODBIT_FOPOINT`: replaces apoapsis with the second focal point of the ellipse.
///
/// # Returns
/// A [`NodesApsidesResult`] with 4 six-element state arrays:
/// - `asc_node`: Ascending node `[lon, lat, dist, lon_speed, lat_speed, dist_speed]`.
/// - `desc_node`: Descending node `[lon, lat, dist, lon_speed, lat_speed, dist_speed]`.
/// - `periapsis`: Periapsis (perihelion / perigee).
/// - `apoapsis`: Apoapsis (aphelion / apogee) or second focal point.
pub fn calc_nodes_apsides(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    body: BodyId,
    input: TimeInput,
    options: CalcOptions,
    method_bits: u32,
) -> Result<NodesApsidesResult, EphemerisError> {
    let obs = match options.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };

    if options.plane == Plane::Horizontal && obs.is_none() {
        return Err(EphemerisError::InvalidInput {
            message: "Plane::Horizontal requires a topocentric Observer".to_owned(),
        });
    }

    // Step size for stencil velocities
    let h_days = if options.plane == Plane::Horizontal {
        H_NODES_APSIDES_HORIZONTAL_DAYS
    } else if body == BodyId::Planet(PlanetId::Moon) {
        H_NODES_APSIDES_MOON_DAYS
    } else {
        H_NODES_APSIDES_PLANET_DAYS
    };

    let (pts_0, ctx_0) =
        eval_points_at_time(spk, time_table, body, input, &options, method_bits, obs)?;

    // Sidereal ayanamsa offset and secular rate if applicable
    let sidereal_info = if options.plane == Plane::Ecliptic {
        if let Zodiac::Sidereal(model) = options.zodiac {
            if options.equinox == Equinox::J2000 {
                return Err(EphemerisError::InvalidInput {
                    message: "Zodiac::Sidereal cannot be combined with Equinox::J2000 (ayanamsa is defined relative to Equinox::OfDate and would double-count precession)".to_owned(),
                });
            }
            let aya = eval_ayanamsa(spk, model, ctx_0.tt.jd_f64(), options.nutation)?;
            Some((aya.effective_deg, aya.rate_deg_day))
        } else {
            None
        }
    } else {
        None
    };

    if options.include_speed {
        let t_m2 = shift_time_input(input, -2.0 * h_days);
        let t_m1 = shift_time_input(input, -h_days);
        let t_p1 = shift_time_input(input, h_days);
        let t_p2 = shift_time_input(input, 2.0 * h_days);

        let (pts_m2, _) =
            eval_points_at_time(spk, time_table, body, t_m2, &options, method_bits, obs)?;
        let (pts_m1, _) =
            eval_points_at_time(spk, time_table, body, t_m1, &options, method_bits, obs)?;
        let (pts_p1, _) =
            eval_points_at_time(spk, time_table, body, t_p1, &options, method_bits, obs)?;
        let (pts_p2, _) =
            eval_points_at_time(spk, time_table, body, t_p2, &options, method_bits, obs)?;

        let mut res_arr = [[0.0_f64; 6]; 4];
        for i in 0..4 {
            let p0 = pts_0.get(i).copied().unwrap_or([0.0, 0.0, 0.0]);
            let pm2 = pts_m2.get(i).copied().unwrap_or([0.0, 0.0, 0.0]);
            let pm1 = pts_m1.get(i).copied().unwrap_or([0.0, 0.0, 0.0]);
            let pp1 = pts_p1.get(i).copied().unwrap_or([0.0, 0.0, 0.0]);
            let pp2 = pts_p2.get(i).copied().unwrap_or([0.0, 0.0, 0.0]);

            let v0 = stencil_4pt_derivative(&pm2, &pm1, &pp1, &pp2, h_days);
            let mut sph = cart_state_to_spherical(&p0, &v0);

            if let Some((aya_deg, aya_rate_deg_day)) = sidereal_info {
                sph.lon_deg = (sph.lon_deg - aya_deg).rem_euclid(360.0);
                sph.lon_speed_deg_day -= aya_rate_deg_day;
            }

            if let Some(dest) = res_arr.get_mut(i) {
                *dest = [
                    sph.lon_deg,
                    sph.lat_deg,
                    sph.dist_au,
                    sph.lon_speed_deg_day,
                    sph.lat_speed_deg_day,
                    sph.dist_speed_au_day,
                ];
            }
        }

        if matches!(
            body,
            BodyId::Planet(PlanetId::Sun | PlanetId::Earth | PlanetId::EarthMoonBarycenter)
        ) {
            res_arr[0] = [0.0; 6];
            res_arr[1] = [0.0; 6];
        }

        #[allow(clippy::indexing_slicing)]
        Ok(NodesApsidesResult {
            asc_node: res_arr[0],
            desc_node: res_arr[1],
            periapsis: res_arr[2],
            apoapsis: res_arr[3],
        })
    } else {
        let mut res_arr = [[0.0_f64; 6]; 4];
        for i in 0..4 {
            let p0 = pts_0.get(i).copied().unwrap_or([0.0, 0.0, 0.0]);
            let c = cart_to_spherical(&p0);
            let mut lon = c.lon_deg;
            if let Some((aya_deg, _)) = sidereal_info {
                lon = (lon - aya_deg).rem_euclid(360.0);
            }
            if let Some(dest) = res_arr.get_mut(i) {
                *dest = [lon, c.lat_deg, c.dist, 0.0, 0.0, 0.0];
            }
        }

        if matches!(
            body,
            BodyId::Planet(PlanetId::Sun | PlanetId::Earth | PlanetId::EarthMoonBarycenter)
        ) {
            res_arr[0] = [0.0; 6];
            res_arr[1] = [0.0; 6];
        }

        #[allow(clippy::indexing_slicing)]
        Ok(NodesApsidesResult {
            asc_node: res_arr[0],
            desc_node: res_arr[1],
            periapsis: res_arr[2],
            apoapsis: res_arr[3],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::moshier::MoshierAnalyticBackend;
    use crate::spk::AniseSpkBackend;
    use crate::types::{AyanamsaModel, NutationMode};

    fn get_test_spk_and_time() -> Result<(Box<dyn SpkSource>, TimeTable), EphemerisError> {
        let workspace_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(|p| p.parent())
            .unwrap_or_else(|| std::path::Path::new("."));
        let de440s_path = workspace_dir.join("data/de440s.bsp");
        let backend: Box<dyn SpkSource> = if de440s_path.exists() {
            Box::new(AniseSpkBackend::from_bsp_paths(&[de440s_path])?)
        } else {
            Box::new(MoshierAnalyticBackend::new())
        };
        let tt = TimeTable::embedded();
        Ok((backend, tt))
    }

    #[test]
    fn test_constants_match_swiss_ephemeris() {
        assert_eq!(SE_NODBIT_MEAN, 1);
        assert_eq!(SE_NODBIT_OSCU, 2);
        assert_eq!(SE_NODBIT_OSCU_BAR, 4);
        assert_eq!(SE_NODBIT_FOPOINT, 256);
    }

    #[test]
    fn test_lunar_osculating_nodes_and_apsides_j2000() -> Result<(), EphemerisError> {
        let (spk, time_table) = get_test_spk_and_time()?;
        let t_j2000 = TimeInput::from_tt_jd(JD_J2000);
        let opts = CalcOptions {
            origin: Origin::Geocentric,
            plane: Plane::Ecliptic,
            equinox: Equinox::OfDate,
            nutation: NutationMode::True,
            aberration: crate::types::AbMode::Geometric,
            zodiac: Zodiac::Tropical,
            include_speed: true,
        };

        let res = calc_nodes_apsides(
            &*spk,
            &time_table,
            BodyId::Planet(PlanetId::Moon),
            t_j2000,
            opts,
            SE_NODBIT_OSCU,
        )?;

        // 1. Ascending and descending node longitudes differ by exactly 180 degrees
        let lon_diff = (res.desc_node[0] - res.asc_node[0] + 360.0).rem_euclid(360.0);
        assert!(
            (lon_diff - 180.0).abs() < 1.0e-5,
            "Node longitudes must differ by 180 deg, got {lon_diff}"
        );

        // 2. Both nodes lie in the ecliptic plane (lat ~ 0 deg)
        assert!(
            res.asc_node[1].abs() < 1.0e-5,
            "Ascending node latitude must be ~0, got {}",
            res.asc_node[1]
        );
        assert!(
            res.desc_node[1].abs() < 1.0e-5,
            "Descending node latitude must be ~0, got {}",
            res.desc_node[1]
        );

        // 3. Node distance ~ 384,400 km (~0.00257 AU)
        let dist_km = res.asc_node[2] * AU_KM;
        assert!(
            (dist_km - 384_400.0).abs() < 30_000.0,
            "Ascending node distance unexpected: {dist_km} km"
        );

        // 4. Perigee distance < Apogee distance
        assert!(
            res.periapsis[2] < res.apoapsis[2],
            "Perigee distance ({}) must be less than apogee distance ({})",
            res.periapsis[2],
            res.apoapsis[2]
        );

        // 5. Speeds are computed and non-zero
        assert!(
            res.asc_node[3].abs() > 1.0e-6,
            "Node speed should be non-zero"
        );
        Ok(())
    }

    #[test]
    fn test_lunar_second_focal_point() -> Result<(), EphemerisError> {
        let (spk, time_table) = get_test_spk_and_time()?;
        let t_j2000 = TimeInput::from_tt_jd(JD_J2000);
        let opts = CalcOptions::default();

        let res_foc = calc_nodes_apsides(
            &*spk,
            &time_table,
            BodyId::Planet(PlanetId::Moon),
            t_j2000,
            opts,
            SE_NODBIT_OSCU | SE_NODBIT_FOPOINT,
        )?;

        // Second focal point distance from Earth is 2 * a * e ~ 2 * 384400 * 0.0549 ~ 42,200 km ~ 0.000282 AU
        let foc_dist_km = res_foc.apoapsis[2] * AU_KM;
        assert!(
            (foc_dist_km - 42_200.0).abs() < 15_000.0,
            "Lunar focal point distance should be ~42,200 km, got {foc_dist_km}"
        );
        Ok(())
    }

    #[test]
    fn test_lunar_mean_nodes_and_apsides() -> Result<(), EphemerisError> {
        let (spk, time_table) = get_test_spk_and_time()?;
        let t_j2000 = TimeInput::from_tt_jd(JD_J2000);
        let opts = CalcOptions {
            origin: Origin::Geocentric,
            plane: Plane::Ecliptic,
            equinox: Equinox::OfDate,
            nutation: NutationMode::Mean,
            aberration: crate::types::AbMode::Geometric,
            zodiac: Zodiac::Tropical,
            include_speed: true,
        };

        let res = calc_nodes_apsides(
            &*spk,
            &time_table,
            BodyId::Planet(PlanetId::Moon),
            t_j2000,
            opts,
            SE_NODBIT_MEAN,
        )?;

        // Simon et al. (1994) / faom03 at J2000 is 125.044555 deg
        let expected_node = faom03(0.0).to_degrees().rem_euclid(360.0);
        assert!(
            (res.asc_node[0] - expected_node).abs() < 1.0e-4,
            "Mean ascending node at J2000 expected {expected_node}, got {}",
            res.asc_node[0]
        );
        let lon_diff = (res.desc_node[0] - res.asc_node[0] + 360.0).rem_euclid(360.0);
        assert!((lon_diff - 180.0).abs() < 1.0e-5);
        Ok(())
    }

    #[test]
    fn test_mars_osculating_nodes_heliocentric() -> Result<(), EphemerisError> {
        let (spk, time_table) = get_test_spk_and_time()?;
        let t_j2000 = TimeInput::from_tt_jd(JD_J2000);
        let opts = CalcOptions {
            origin: Origin::Heliocentric,
            plane: Plane::Ecliptic,
            equinox: Equinox::J2000,
            nutation: NutationMode::Mean,
            aberration: crate::types::AbMode::Geometric,
            zodiac: Zodiac::Tropical,
            include_speed: true,
        };

        let res = calc_nodes_apsides(
            &*spk,
            &time_table,
            BodyId::Planet(PlanetId::Mars),
            t_j2000,
            opts,
            SE_NODBIT_OSCU,
        )?;

        // Mars heliocentric ascending node ~ 49.56 deg, latitude = 0
        assert!(
            (res.asc_node[0] - 49.56).abs() < 1.0,
            "Mars ascending node expected ~49.56, got {}",
            res.asc_node[0]
        );
        assert!(res.asc_node[1].abs() < 1.0e-5);

        // Mars perihelion ~ 1.38 AU, aphelion ~ 1.66 AU
        assert!(
            (res.periapsis[2] - 1.38).abs() < 0.05,
            "Mars perihelion distance expected ~1.38 AU, got {}",
            res.periapsis[2]
        );
        assert!(
            (res.apoapsis[2] - 1.66).abs() < 0.05,
            "Mars aphelion distance expected ~1.66 AU, got {}",
            res.apoapsis[2]
        );
        Ok(())
    }

    #[test]
    fn test_geocentric_observer_shift_for_planets() -> Result<(), EphemerisError> {
        let (spk, time_table) = get_test_spk_and_time()?;
        let t_j2000 = TimeInput::from_tt_jd(JD_J2000);

        let opts_helio = CalcOptions {
            origin: Origin::Heliocentric,
            plane: Plane::Ecliptic,
            equinox: Equinox::J2000,
            nutation: NutationMode::Mean,
            aberration: crate::types::AbMode::Geometric,
            zodiac: Zodiac::Tropical,
            include_speed: false,
        };
        let opts_geo = CalcOptions {
            origin: Origin::Geocentric,
            ..opts_helio
        };

        let res_h = calc_nodes_apsides(
            &*spk,
            &time_table,
            BodyId::Planet(PlanetId::Mars),
            t_j2000,
            opts_helio,
            SE_NODBIT_OSCU,
        )?;

        let res_g = calc_nodes_apsides(
            &*spk,
            &time_table,
            BodyId::Planet(PlanetId::Mars),
            t_j2000,
            opts_geo,
            SE_NODBIT_OSCU,
        )?;

        // Geocentric distance to Mars' node is distinct from heliocentric distance (shifted by Earth-Sun vector)
        assert!(
            (res_h.asc_node[2] - res_g.asc_node[2]).abs() > 0.1,
            "Geocentric shift must modify distance from observer"
        );
        Ok(())
    }

    #[test]
    fn test_sidereal_zodiac_and_equatorial_plane() -> Result<(), EphemerisError> {
        let (spk, time_table) = get_test_spk_and_time()?;
        let t_j2000 = TimeInput::from_tt_jd(JD_J2000);

        // Sidereal Lahiri on Moon mean nodes
        let opts_sid = CalcOptions {
            origin: Origin::Geocentric,
            plane: Plane::Ecliptic,
            equinox: Equinox::OfDate,
            nutation: NutationMode::Mean,
            aberration: crate::types::AbMode::Geometric,
            zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
            include_speed: true,
        };
        let res_sid = calc_nodes_apsides(
            &*spk,
            &time_table,
            BodyId::Planet(PlanetId::Moon),
            t_j2000,
            opts_sid,
            SE_NODBIT_MEAN,
        )?;

        // Lahiri ayanamsa at J2000 is ~23.85 deg, so Tropical 125.04 -> Sidereal ~101.19
        assert!(
            (res_sid.asc_node[0] - 101.19).abs() < 0.1,
            "Sidereal node longitude expected ~101.19, got {}",
            res_sid.asc_node[0]
        );

        // Equatorial plane: RA/Dec
        let opts_eq = CalcOptions {
            origin: Origin::Geocentric,
            plane: Plane::Equatorial,
            equinox: Equinox::OfDate,
            nutation: NutationMode::True,
            aberration: crate::types::AbMode::Geometric,
            zodiac: Zodiac::Tropical,
            include_speed: true,
        };
        let res_eq = calc_nodes_apsides(
            &*spk,
            &time_table,
            BodyId::Planet(PlanetId::Moon),
            t_j2000,
            opts_eq,
            SE_NODBIT_OSCU,
        )?;

        assert!((0.0..360.0).contains(&res_eq.asc_node[0]));
        assert!((-90.0..=90.0).contains(&res_eq.asc_node[1]));
        Ok(())
    }

    #[test]
    fn test_saturn_mean_perihelion_year_2100() -> Result<(), EphemerisError> {
        let (spk, time_table) = get_test_spk_and_time()?;
        let t_2100 = TimeInput::from_tt_jd(JD_J2000 + 36525.0);

        let opts = CalcOptions {
            origin: Origin::Heliocentric,
            plane: Plane::Ecliptic,
            equinox: Equinox::J2000,
            nutation: NutationMode::Mean,
            aberration: crate::types::AbMode::Geometric,
            zodiac: Zodiac::Tropical,
            include_speed: false,
        };
        let res = calc_nodes_apsides(
            &*spk,
            &time_table,
            BodyId::Planet(PlanetId::Saturn),
            t_2100,
            opts,
            SE_NODBIT_MEAN,
        )?;

        // Periapsis longitude: w0 = 92.59887831, rate = -0.41897216 deg/century
        // Expected orbital varpi at J2000 + 36525.0: 92.59887831 - 0.41897216 = 92.17990615 deg.
        // Reduced to ecliptic (inc = 2.486 deg): lambda_peri = 92.1981 deg.
        let perihelion_lon = res.periapsis[0];
        assert!(
            (perihelion_lon - 92.1981).abs() < 0.01,
            "Expected Saturn mean perihelion ~92.20 deg in year 2100, got {perihelion_lon}"
        );
        Ok(())
    }

    #[test]
    fn test_horizontal_plane_node_speed_stencil() -> Result<(), EphemerisError> {
        let (spk, time_table) = get_test_spk_and_time()?;
        let t_j2000 = TimeInput::from_tt_jd(JD_J2000);
        let obs = Observer::new(10.0, 50.0, 100.0).expect("valid observer");

        let opts_hz = CalcOptions {
            origin: Origin::Topocentric(obs),
            plane: Plane::Horizontal,
            equinox: Equinox::OfDate,
            nutation: NutationMode::True,
            aberration: crate::types::AbMode::Geometric,
            zodiac: Zodiac::Tropical,
            include_speed: true,
        };
        let res_hz = calc_nodes_apsides(
            &*spk,
            &time_table,
            BodyId::Planet(PlanetId::Mars),
            t_j2000,
            opts_hz,
            SE_NODBIT_OSCU,
        )?;

        // Topocentric horizontal azimuth speed should reflect diurnal rotation (~360 deg/day)
        let az_speed = res_hz.asc_node[3];
        assert!(
            az_speed.is_finite() && az_speed.abs() > 100.0,
            "Azimuth speed in horizontal plane should reflect diurnal rate, got {az_speed}"
        );
        Ok(())
    }
}
