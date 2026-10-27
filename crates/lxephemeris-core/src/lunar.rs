//! 4 Lunar Nodes (`Mean_Node`, `True_Node`, `Mean_South_Node`, `True_South_Node`) and
//! 3 Lilith / Black Moon Apogee variants (`Mean_Lilith`, `True_Lilith`, `Interpolated_Lilith`),
//! strictly implementing `MODEL_PASSPORTS.md` §2:
//! - IERS Conventions 2003 / Simon et al. (1994) fundamental argument polynomials (`faom03`, `faf03`, `fal03`).
//! - Osculating elements from geometric geocentric `DE440s` Moon state $(\mathbf{r}, \mathbf{v})$ with
//!   $\mu_{\text{EM}} = GM_\oplus + GM_{\text{Moon}} = 403503.235625\text{ km}^3/\text{s}^2$.
//! - `Interpolated_Lilith` (`natural_apogee_c2_spline_v1`): Brent root-finding on radial velocity
//!   $\dot{r}(\tau) = 0$ over $[t - 60\text{ d}, t + 60\text{ d}]$ + Thomas $C^2$ natural cubic spline.

use crate::astrometry::{equatorial_to_horizontal_vector_au, AtmosphereParams, CalcResult};
use crate::error::EphemerisError;
use crate::frame::{
    cart_state_to_spherical, cart_to_spherical, is_within_iau2006_domain, mat_add_scaled, mat_mul,
    mat_transpose, mat_vec, rot_x, stencil_4pt_derivative, vec_add, vec_cross, vec_dot, vec_norm,
    vec_scale, vec_sub, Vec3,
};
use crate::moment::{Moment, H_PLANET_DAYS};
use crate::spk::SpkSource;
use crate::time::{tt_to_tdb, Tt};
use crate::types::{
    BodyId, CalcOptions, Equinox, LunarPointId, NutationMode, Origin, Plane, TargetKind, AU_KM,
    JD_J2000,
};
use sofars::fundargs::{faf03, fal03, faom03};
use sofars::pnp::{bp06, ltp, ltpecl, numat, nut00a, obl06};

/// Combined Earth + Moon gravitational parameter $\mu_{\text{EM}} = GM_\oplus + GM_{\text{Moon}}$
/// in $\text{km}^3/\text{s}^2$ from the JPL DE440 header (Park et al. 2021, Table 2).
pub const MU_EM_KM3_S2: f64 = 403_503.235_625;

/// Nominal mean Earth-Moon distance in AU (`385000.56 km / AU_KM`).
pub const MEAN_MOON_DIST_AU: f64 = 385_000.56 / AU_KM;

/// Safe lower bound (`JD_TT`) for `Interpolated_Lilith` on `de440s.bsp` (`1850-03-01T00:00:00 TT`).
pub const JD_INTERP_LILITH_SAFE_MIN: f64 = 2_396_817.5;

/// Safe upper bound (`JD_TT`) for `Interpolated_Lilith` on `de440s.bsp` (`2149-11-01T00:00:00 TT`).
pub const JD_INTERP_LILITH_SAFE_MAX: f64 = 2_506_270.5;

/// Computes the 3x3 rotation matrix from ICRF to Mean or True Ecliptic of Date at `jd_tt`.
fn icrf_to_ecl_matrix_at_tt(
    jd_tt: f64,
    nutation: bool,
    cached_nutation: Option<(f64, f64)>,
) -> [[f64; 3]; 3] {
    let tt = Tt::from_jd(jd_tt);
    let (eps0, r_bp) = if is_within_iau2006_domain(jd_tt) {
        let e0 = obl06(tt.0.jd1, tt.0.jd2);
        let (_, _, rbp) = bp06(tt.0.jd1, tt.0.jd2);
        (e0, rbp)
    } else {
        let epj = 2000.0 + (jd_tt - JD_J2000) / 365.25;
        let rbp = ltp(epj);
        let eq_pole = [rbp[2][0], rbp[2][1], rbp[2][2]];
        let ecl_pole = ltpecl(epj);
        let e0 = vec_dot(&eq_pole, &ecl_pole).clamp(-1.0, 1.0).acos();
        (e0, rbp)
    };

    if nutation {
        let (dpsi, deps) = cached_nutation.unwrap_or_else(|| nut00a(tt.0.jd1, tt.0.jd2));
        let r_nut = numat(eps0, dpsi, deps);
        let r_eq = mat_mul(&r_nut, &r_bp);
        mat_mul(&rot_x(eps0 + deps), &r_eq)
    } else {
        mat_mul(&rot_x(eps0), &r_bp)
    }
}

/// Evaluates the geometric geocentric Moon state `(r_km, v_km_s)` in ICRF at `jd_tt`,
/// falling back to the analytical Moshier source if `jd_tt` is slightly outside the SPK kernel boundary.
fn moon_geo_icrf_km_s(spk: &dyn SpkSource, jd_tt: f64) -> Result<(Vec3, Vec3), EphemerisError> {
    let tt = Tt::from_jd(jd_tt);
    let tdb = tt_to_tdb(tt, None);
    let (st_moon, st_earth) = match (
        spk.state_icrf_km_split(301, 0, tdb.0.jd1, tdb.0.jd2),
        spk.state_icrf_km_split(399, 0, tdb.0.jd1, tdb.0.jd2),
    ) {
        (Ok(m), Ok(e)) => (m, e),
        (Err(EphemerisError::OutOfCoverage { .. }), _)
        | (_, Err(EphemerisError::OutOfCoverage { .. })) => {
            let mosh = crate::moshier::MoshierAnalyticBackend::new();
            (
                mosh.state_icrf_km_split(301, 0, tdb.0.jd1, tdb.0.jd2)?,
                mosh.state_icrf_km_split(399, 0, tdb.0.jd1, tdb.0.jd2)?,
            )
        }
        (Err(e), _) | (_, Err(e)) => return Err(e),
    };
    Ok((
        vec_sub(&st_moon.pos_km, &st_earth.pos_km),
        vec_sub(&st_moon.vel_km_s, &st_earth.vel_km_s),
    ))
}

/// Evaluates osculating `(true_node_lon_deg, true_lilith_lon_deg, true_lilith_lat_deg, apogee_dist_au)`
/// from a geometric geocentric Moon state `(r_icrf, v_icrf)` and `r_ecl` matrix.
fn osculating_from_state_and_matrix(
    r_icrf: &Vec3,
    v_icrf: &Vec3,
    r_ecl: &[[f64; 3]; 3],
) -> (f64, f64, f64, f64) {
    let r = mat_vec(r_ecl, r_icrf);
    let v = mat_vec(r_ecl, v_icrf);

    let h_vec = vec_cross(&r, &v);
    let node_lon_deg = h_vec[0].atan2(-h_vec[1]).to_degrees().rem_euclid(360.0);

    let r_norm = vec_norm(&r);
    let v_norm2 = vec_dot(&v, &v);
    let v_cross_h = vec_cross(&v, &h_vec);
    let e_vec = vec_sub(
        &vec_scale(&v_cross_h, 1.0 / MU_EM_KM3_S2),
        &vec_scale(&r, 1.0 / r_norm),
    );
    let a_vec = vec_scale(&e_vec, -1.0);
    let sph = cart_to_spherical(&a_vec);

    let inv_a = 2.0 / r_norm - v_norm2 / MU_EM_KM3_S2;
    let apo_dist_au = if inv_a > 1e-12 {
        ((1.0 / inv_a) * (1.0 + sph.dist)) / AU_KM
    } else {
        MEAN_MOON_DIST_AU
    };

    (node_lon_deg, sph.lon_deg, sph.lat_deg, apo_dist_au)
}

/// Deterministic Brent's root-finding algorithm on `[a, b]` matching `MODEL_PASSPORTS.md` §2.3.
fn brent_root<F>(
    mut f: F,
    mut a: f64,
    mut b: f64,
    tol: f64,
    max_iter: u32,
) -> Result<f64, EphemerisError>
where
    F: FnMut(f64) -> Result<f64, EphemerisError>,
{
    let mut fa = f(a)?;
    let mut fb = f(b)?;
    if fa == 0.0 {
        return Ok(a);
    }
    if fb == 0.0 {
        return Ok(b);
    }
    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;

    for _ in 0..max_iter {
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol1 = 2.0 * f64::EPSILON * b.abs() + 0.5 * tol;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 || fb == 0.0 {
            return Ok(b);
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (mut p, mut q) = if (a - c).abs() == 0.0 {
                (2.0 * xm * s, 1.0 - s)
            } else {
                let q_val = fa / fc;
                let r_val = fb / fc;
                (
                    s * (2.0 * xm * q_val * (q_val - r_val) - (b - a) * (r_val - 1.0)),
                    (q_val - 1.0) * (r_val - 1.0) * (s - 1.0),
                )
            };
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            if 2.0 * p < (3.0 * xm * q - (tol1 * q).abs()).min((e * q).abs()) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if d.abs() > tol1 {
            b += d;
        } else {
            b += tol1.copysign(xm);
        }
        fb = f(b)?;
    }

    Err(EphemerisError::NonConvergence {
        solver: "interpolated_lilith_brent".to_owned(),
        iterations: max_iter,
        residual: fb.abs(),
    })
}

/// Natural $C^2$ cubic spline ($S''(x_0) = S''(x_{n-1}) = 0$) solved via the Thomas algorithm
/// using stack-allocated buffers (`n <= 32`), returning `(value, derivative)` at `x_eval`.
#[allow(clippy::indexing_slicing)]
fn natural_cubic_spline_eval(xs: &[f64], ys: &[f64], x_eval: f64) -> (f64, f64) {
    let n = xs.len().min(32);
    if n < 2 {
        return (ys.first().copied().unwrap_or(0.0), 0.0);
    }
    let mut h = [0.0_f64; 32];
    for i in 0..(n - 1) {
        h[i] = xs[i + 1] - xs[i];
    }
    let mut alpha = [0.0_f64; 32];
    for i in 1..(n - 1) {
        alpha[i] =
            (3.0 / h[i]) * (ys[i + 1] - ys[i]) - (3.0 / h[i - 1]) * (ys[i] - ys[i - 1]);
    }
    let mut l = [0.0_f64; 32];
    let mut mu = [0.0_f64; 32];
    let mut z = [0.0_f64; 32];
    l[0] = 1.0;
    for i in 1..(n - 1) {
        l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    let mut c = [0.0_f64; 32];
    let mut b = [0.0_f64; 32];
    let mut d = [0.0_f64; 32];
    for j in (0..(n - 1)).rev() {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (ys[j + 1] - ys[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }

    let mut idx = if x_eval >= xs[n - 1] {
        n - 2
    } else {
        0
    };
    for i in 0..(n - 1) {
        if xs[i] <= x_eval && x_eval <= xs[i + 1] {
            idx = i;
            break;
        }
    }
    let dx = x_eval - xs[idx];
    let val = ys[idx] + b[idx] * dx + c[idx] * dx * dx + d[idx] * dx * dx * dx;
    let deriv = b[idx] + 2.0 * c[idx] * dx + 3.0 * d[idx] * dx * dx;
    (val, deriv)
}

/// Evaluates `Interpolated_Lilith` (`natural_apogee_c2_spline_v1`) at `jd_tt` (`MODEL_PASSPORTS.md` §2.3)
/// using a knot-bracketed $C^2$ quintic-Hermite blend of adjacent 5-knot natural cubic splines so that
/// position, velocity, and acceleration are strictly $C^2$ continuous across all apogee knots and free of
/// moving-window boundary jumps.
#[allow(clippy::indexing_slicing)]
fn eval_interpolated_apsis_c2(
    spk: &dyn SpkSource,
    jd_tt: f64,
    nutation: bool,
    want_apogee: bool,
) -> Result<(f64, f64, f64, f64, f64, f64), EphemerisError> {
    let radial_rate = |t_tt: f64| -> Result<f64, EphemerisError> {
        let (r, v) = moon_geo_icrf_km_s(spk, t_tt)?;
        Ok(vec_dot(&r, &v) / vec_norm(&r))
    };

    let t_start = jd_tt - 96.0;
    let step = 2.0;
    let mut grid_t = [0.0_f64; 97];
    let mut grid_g = [0.0_f64; 97];
    for i in 0..=96 {
        let t = t_start + (i as f64) * step;
        grid_t[i] = t;
        grid_g[i] = radial_rate(t)?;
    }

    let mut knots_t = [0.0_f64; 16];
    let mut knots_lon = [0.0_f64; 16];
    let mut knots_lat = [0.0_f64; 16];
    let mut knots_dist = [0.0_f64; 16];
    let mut n_knots = 0usize;

    for i in 0..96 {
        if grid_g[i] * grid_g[i + 1] <= 0.0 && n_knots < 16 {
            let is_target_apsis = if want_apogee {
                grid_g[i] > grid_g[i + 1]
            } else {
                grid_g[i] < grid_g[i + 1]
            };
            if is_target_apsis {
                let t_root = brent_root(&radial_rate, grid_t[i], grid_t[i + 1], 1e-11, 100)?;
                let (r_icrf, _) = moon_geo_icrf_km_s(spk, t_root)?;
                let r_ecl = icrf_to_ecl_matrix_at_tt(t_root, false, None);
                let sph = cart_to_spherical(&mat_vec(&r_ecl, &r_icrf));
                knots_lon[n_knots] = sph.lon_deg;
                knots_lat[n_knots] = sph.lat_deg;
                knots_dist[n_knots] = sph.dist / AU_KM;
                knots_t[n_knots] = t_root;
                n_knots += 1;
            }
        }
    }

    if n_knots < 3 {
        return Err(EphemerisError::NonConvergence {
            solver: "interpolated_lilith_apsis_knots".to_owned(),
            iterations: 96,
            residual: n_knots as f64,
        });
    }

    let mut unwrapped_lon = [0.0_f64; 16];
    unwrapped_lon[0] = knots_lon[0];
    for k in 1..n_knots {
        let prev = unwrapped_lon[k - 1];
        let diff = (knots_lon[k] - prev + 180.0).rem_euclid(360.0) - 180.0;
        unwrapped_lon[k] = prev + diff;
    }

    // Locate the active apogee interval [t_k, t_{k+1}] containing jd_tt
    let mut k_iv = 0usize;
    for k in 0..(n_knots - 1) {
        if knots_t[k] <= jd_tt && jd_tt < knots_t[k + 1] {
            k_iv = k;
            break;
        }
    }

    let (lon_mean, dlon, lat_mean, dlat, dist_au, ddist) = if k_iv >= 2 && k_iv + 3 < n_knots {
        let kt_l = &knots_t[(k_iv - 2)..=(k_iv + 2)];
        let kt_r = &knots_t[(k_iv - 1)..=(k_iv + 3)];
        let h_iv = (knots_t[k_iv + 1] - knots_t[k_iv]).max(1e-12);
        let u = ((jd_tt - knots_t[k_iv]) / h_iv).clamp(0.0, 1.0);
        // C^2 quintic Hermite step: sigma(0)=0, sigma(1)=1, sigma'(0)=sigma'(1)=sigma''(0)=sigma''(1)=0
        let sigma = u * u * u * (10.0 + u * (-15.0 + 6.0 * u));
        let dsigma_dt = (30.0 * u * u * (1.0 - u) * (1.0 - u)) / h_iv;

        let blend = |ys_all: &[f64]| -> (f64, f64) {
            let (v_l, dv_l) = natural_cubic_spline_eval(kt_l, &ys_all[(k_iv - 2)..=(k_iv + 2)], jd_tt);
            let (v_r, dv_r) = natural_cubic_spline_eval(kt_r, &ys_all[(k_iv - 1)..=(k_iv + 3)], jd_tt);
            let val = (1.0 - sigma) * v_l + sigma * v_r;
            let deriv = (1.0 - sigma) * dv_l + sigma * dv_r + dsigma_dt * (v_r - v_l);
            (val, deriv)
        };

        let (lm, dlm) = blend(&unwrapped_lon[..n_knots]);
        let (bm, dbm) = blend(&knots_lat[..n_knots]);
        let (rm, drm) = blend(&knots_dist[..n_knots]);
        (lm, dlm, bm, dbm, rm, drm)
    } else {
        let kt = &knots_t[..n_knots];
        let (lm, dlm) = natural_cubic_spline_eval(kt, &unwrapped_lon[..n_knots], jd_tt);
        let (bm, dbm) = natural_cubic_spline_eval(kt, &knots_lat[..n_knots], jd_tt);
        let (rm, drm) = natural_cubic_spline_eval(kt, &knots_dist[..n_knots], jd_tt);
        (lm, dlm, bm, dbm, rm, drm)
    };

    let mut lon_out = lon_mean.rem_euclid(360.0);
    if nutation {
        let tt = Tt::from_jd(jd_tt);
        let (dpsi_rad, _) = nut00a(tt.0.jd1, tt.0.jd2);
        lon_out = (lon_out + dpsi_rad.to_degrees()).rem_euclid(360.0);
    }

    Ok((lon_out, lat_mean, dist_au, dlon, dlat, ddist))
}

/// Evaluates any of the 7 [`LunarPointId`] variants inside a [`Moment`].
pub fn calc_lunar_point_in_moment(
    moment: &Moment<'_>,
    point: LunarPointId,
    options: &CalcOptions,
    sidereal_offset_and_rate: Option<(f64, f64)>,
) -> Result<CalcResult, EphemerisError> {
    let center = moment.center();
    let jd_tt = center.tt.jd_f64();
    let nut_on = options.nutation == NutationMode::True;
    let dpsi_deg = center.frame.dpsi_rad.to_degrees();

    let (mut lon_deg, mut lat_deg, mut dist_au, mut lon_speed_deg_day, mut lat_speed_deg_day, mut dist_speed_au_day) =
        match point {
            LunarPointId::MeanNode | LunarPointId::MeanSouthNode => {
                let big_t = (jd_tt - JD_J2000) / 36_525.0;
                let om_mean_deg = faom03(big_t).to_degrees().rem_euclid(360.0);
                let mut lon = if nut_on {
                    (om_mean_deg + dpsi_deg).rem_euclid(360.0)
                } else {
                    om_mean_deg
                };
                if point == LunarPointId::MeanSouthNode {
                    lon = (lon + 180.0).rem_euclid(360.0);
                }
                let rate = (-6962890.5431
                    + 2.0 * 7.4722 * big_t
                    + 3.0 * 0.007_702 * big_t * big_t
                    - 4.0 * 0.000_059_39 * big_t * big_t * big_t)
                    / (3_600.0 * 36_525.0);
                (lon, 0.0, MEAN_MOON_DIST_AU, rate, 0.0, 0.0)
            }
            LunarPointId::MeanLilith => {
                let big_t = (jd_tt - JD_J2000) / 36_525.0;
                let f_rad = faf03(big_t);
                let om_rad = faom03(big_t);
                let l_rad = fal03(big_t);
                let mlil_mean_deg = ((f_rad + om_rad - l_rad).to_degrees() + 180.0).rem_euclid(360.0);
                let lon = if nut_on {
                    (mlil_mean_deg + dpsi_deg).rem_euclid(360.0)
                } else {
                    mlil_mean_deg
                };
                let rate = (14648449.0869
                    - 2.0 * 37.1582 * big_t
                    - 3.0 * 0.044_970 * big_t * big_t
                    + 4.0 * 0.000_189_48 * big_t * big_t * big_t)
                    / (3_600.0 * 36_525.0);
                (lon, 0.0, MEAN_MOON_DIST_AU, rate, 0.0, 0.0)
            }
            LunarPointId::TrueNode | LunarPointId::TrueSouthNode | LunarPointId::TrueLilith => {
                // At t_0, reuse center.earth_bary_km and center.frame (r_ecl_true / r_ecl_mean)
                let st_moon_0 =
                    moment
                        .spk()
                        .state_icrf_km_split(301, 0, center.tdb.0.jd1, center.tdb.0.jd2)?;
                let r_icrf_0 = vec_sub(&st_moon_0.pos_km, &center.earth_bary_km.pos_km);
                let v_icrf_0 = vec_sub(&st_moon_0.vel_km_s, &center.earth_bary_km.vel_km_s);
                let r_ecl_0 = if nut_on {
                    &center.frame.r_ecl_true
                } else {
                    &center.frame.r_ecl_mean
                };
                let (tn_0, tl_0, tlat_0, apo_dist_0) =
                    osculating_from_state_and_matrix(&r_icrf_0, &v_icrf_0, r_ecl_0);

                let (tn_speed, tl_speed, tlat_speed, apo_speed) = if options.include_speed {
                    let h = 1.0 / 24.0;
                    let nut_for_speed = if sidereal_offset_and_rate.is_some() {
                        false
                    } else {
                        nut_on
                    };
                    let dpsi_0 = center.frame.dpsi_rad;
                    let deps_0 = center.frame.deps_rad;
                    let eps0_0 = center.frame.eps0_rad;
                    let (ddpsi_rate, ddeps_rate) = if nut_for_speed {
                        moment.nutation_rate_rad_day()
                    } else {
                        (0.0, 0.0)
                    };
                    let (deps0_rate, dr_bp_rate) = moment.precession_rate_per_day();

                    let eval_shifted_osculating =
                        |dt_days: f64| -> Result<(f64, f64, f64, f64), EphemerisError> {
                            let tdb_k = center.tdb.0.add_days(dt_days);
                            let st_moon = moment.spk().state_icrf_km_split(
                                301,
                                0,
                                tdb_k.jd1,
                                tdb_k.jd2,
                            )?;
                            let st_earth = moment.spk().state_icrf_km_split(
                                399,
                                0,
                                tdb_k.jd1,
                                tdb_k.jd2,
                            )?;
                            let r_icrf = vec_sub(&st_moon.pos_km, &st_earth.pos_km);
                            let v_icrf = vec_sub(&st_moon.vel_km_s, &st_earth.vel_km_s);

                            let eps0_k = eps0_0 + dt_days * deps0_rate;
                            let r_bp_k = mat_add_scaled(&center.frame.r_bp, dr_bp_rate, dt_days);
                            let r_ecl_k = if nut_for_speed {
                                let dpsi_k = dpsi_0 + dt_days * ddpsi_rate;
                                let deps_k = deps_0 + dt_days * ddeps_rate;
                                let r_nut = numat(eps0_k, dpsi_k, deps_k);
                                let r_eq = mat_mul(&r_nut, &r_bp_k);
                                mat_mul(&rot_x(eps0_k + deps_k), &r_eq)
                            } else {
                                mat_mul(&rot_x(eps0_k), &r_bp_k)
                            };
                            Ok(osculating_from_state_and_matrix(
                                &r_icrf, &v_icrf, &r_ecl_k,
                            ))
                        };

                    let (tn_m2, tl_m2, lat_m2, apo_m2) = eval_shifted_osculating(-2.0 * h)?;
                    let (tn_m1, tl_m1, lat_m1, apo_m1) = eval_shifted_osculating(-h)?;
                    let (tn_p1, tl_p1, lat_p1, apo_p1) = eval_shifted_osculating(h)?;
                    let (tn_p2, tl_p2, lat_p2, apo_p2) = eval_shifted_osculating(2.0 * h)?;

                    let stencil_angle = |m2: f64, m1: f64, p1: f64, p2: f64| -> f64 {
                        let d1 = (p1 - m1 + 180.0).rem_euclid(360.0) - 180.0;
                        let d2 = (p2 - m2 + 180.0).rem_euclid(360.0) - 180.0;
                        (8.0 * d1 - d2) / (12.0 * h)
                    };
                    (
                        stencil_angle(tn_m2, tn_m1, tn_p1, tn_p2),
                        stencil_angle(tl_m2, tl_m1, tl_p1, tl_p2),
                        (lat_m2 - 8.0 * lat_m1 + 8.0 * lat_p1 - lat_p2) / (12.0 * h),
                        (apo_m2 - 8.0 * apo_m1 + 8.0 * apo_p1 - apo_p2) / (12.0 * h),
                    )
                } else {
                    (0.0, 0.0, 0.0, 0.0)
                };

                match point {
                    LunarPointId::TrueNode => (tn_0, 0.0, MEAN_MOON_DIST_AU, tn_speed, 0.0, 0.0),
                    LunarPointId::TrueSouthNode => (
                        (tn_0 + 180.0).rem_euclid(360.0),
                        0.0,
                        MEAN_MOON_DIST_AU,
                        tn_speed,
                        0.0,
                        0.0,
                    ),
                    LunarPointId::TrueLilith => {
                        (tl_0, tlat_0, apo_dist_0, tl_speed, tlat_speed, apo_speed)
                    }
                    _ => (0.0, 0.0, MEAN_MOON_DIST_AU, 0.0, 0.0, 0.0),
                }
            }
            LunarPointId::InterpolatedLilith => {
                eval_interpolated_apsis_c2(moment.spk(), jd_tt, nut_on, true)?
            }
            LunarPointId::InterpolatedPerigee => {
                eval_interpolated_apsis_c2(moment.spk(), jd_tt, nut_on, false)?
            }
        };

    if !options.include_speed {
        lon_speed_deg_day = 0.0;
        lat_speed_deg_day = 0.0;
        dist_speed_au_day = 0.0;
    } else if nut_on
        && sidereal_offset_and_rate.is_none()
        && matches!(
            point,
            LunarPointId::MeanNode
                | LunarPointId::MeanSouthNode
                | LunarPointId::MeanLilith
                | LunarPointId::InterpolatedLilith
                | LunarPointId::InterpolatedPerigee
        )
    {
        lon_speed_deg_day += moment.nutation_rate_rad_day().0.to_degrees();
    }

    let sph_to_cart_pos_vel =
        |l_deg: f64, b_deg: f64, r_au: f64, dl_deg: f64, db_deg: f64, dr_au: f64| -> (Vec3, Vec3) {
            let l_r = l_deg.to_radians();
            let b_r = b_deg.to_radians();
            let dl_r = dl_deg.to_radians();
            let db_r = db_deg.to_radians();
            let (sin_l, cos_l) = l_r.sin_cos();
            let (sin_b, cos_b) = b_r.sin_cos();
            let p = [r_au * cos_b * cos_l, r_au * cos_b * sin_l, r_au * sin_b];
            let v = if options.include_speed {
                [
                    dr_au * cos_b * cos_l
                        - r_au * db_r * sin_b * cos_l
                        - r_au * dl_r * cos_b * sin_l,
                    dr_au * cos_b * sin_l
                        - r_au * db_r * sin_b * sin_l
                        + r_au * dl_r * cos_b * cos_l,
                    dr_au * sin_b + r_au * db_r * cos_b,
                ]
            } else {
                [0.0, 0.0, 0.0]
            };
            (p, v)
        };

    let (pos_out_au, vel_out_au_day) =
        if options.plane == Plane::Ecliptic && options.equinox == Equinox::OfDate {
            if let Some((aya_deg, aya_rate_deg_day)) = sidereal_offset_and_rate {
                lon_deg = (lon_deg - aya_deg).rem_euclid(360.0);
                if options.include_speed {
                    lon_speed_deg_day -= aya_rate_deg_day;
                }
            }
            sph_to_cart_pos_vel(
                lon_deg,
                lat_deg,
                dist_au,
                lon_speed_deg_day,
                lat_speed_deg_day,
                dist_speed_au_day,
            )
        } else {
            let (pos_ecl_au, vel_ecl_au_day) = sph_to_cart_pos_vel(
                lon_deg,
                lat_deg,
                dist_au,
                lon_speed_deg_day,
                lat_speed_deg_day,
                dist_speed_au_day,
            );
            let r_ecl_from_icrf = if nut_on {
                center.frame.r_ecl_true
            } else {
                center.frame.r_ecl_mean
            };
            let r_icrf_from_ecl = mat_transpose(&r_ecl_from_icrf);
            let pos_icrf = mat_vec(&r_icrf_from_ecl, &pos_ecl_au);
            let vel_icrf = mat_vec(&r_icrf_from_ecl, &vel_ecl_au_day);
            let r_target =
                center
                    .frame
                    .icrf_to_frame_matrix(options.plane, options.equinox, options.nutation);
            let p_frame = mat_vec(&r_target, &pos_icrf);
            let v_frame = mat_vec(&r_target, &vel_icrf);

            if options.plane == Plane::Horizontal {
                let obs = match options.origin {
                    Origin::Topocentric(o) => o,
                    _ => center.observer.ok_or_else(|| EphemerisError::InvalidInput {
                        message: "Plane::Horizontal requires a topocentric Observer".to_owned(),
                    })?,
                };
                let p_horiz = equatorial_to_horizontal_vector_au(
                    obs,
                    center.gast_rad,
                    &p_frame,
                    AtmosphereParams::default(),
                );
                if options.include_speed {
                    let stencil = moment.planet_stencil_contexts()?;
                    let h = H_PLANET_DAYS;
                    let eval_hz = |ctx_k: &crate::moment::EpochContext, dt_days: f64| -> Vec3 {
                        let p_k = vec_add(&p_frame, &vec_scale(&v_frame, dt_days));
                        equatorial_to_horizontal_vector_au(
                            obs,
                            ctx_k.gast_rad,
                            &p_k,
                            AtmosphereParams::default(),
                        )
                    };
                    let hz_m2 = eval_hz(&stencil[0], -2.0 * h);
                    let hz_m1 = eval_hz(&stencil[1], -h);
                    let hz_p1 = eval_hz(&stencil[2], h);
                    let hz_p2 = eval_hz(&stencil[3], 2.0 * h);
                    let v_horiz = stencil_4pt_derivative(&hz_m2, &hz_m1, &hz_p1, &hz_p2, h);
                    let sph_h = cart_state_to_spherical(&p_horiz, &v_horiz);
                    lon_deg = sph_h.lon_deg;
                    lat_deg = sph_h.lat_deg;
                    dist_au = sph_h.dist_au;
                    lon_speed_deg_day = sph_h.lon_speed_deg_day;
                    lat_speed_deg_day = sph_h.lat_speed_deg_day;
                    dist_speed_au_day = sph_h.dist_speed_au_day;
                    (p_horiz, v_horiz)
                } else {
                    let sph_h = cart_to_spherical(&p_horiz);
                    lon_deg = sph_h.lon_deg;
                    lat_deg = sph_h.lat_deg;
                    dist_au = sph_h.dist;
                    (p_horiz, [0.0, 0.0, 0.0])
                }
            } else {
                let mut sph_t = if options.include_speed {
                    cart_state_to_spherical(&p_frame, &v_frame)
                } else {
                    let c = cart_to_spherical(&p_frame);
                    crate::frame::SphericalState {
                        lon_deg: c.lon_deg,
                        lat_deg: c.lat_deg,
                        dist_au: c.dist,
                        lon_speed_deg_day: 0.0,
                        lat_speed_deg_day: 0.0,
                        dist_speed_au_day: 0.0,
                    }
                };
                let mut p_out = p_frame;
                let mut v_out = v_frame;
                if options.plane == Plane::Ecliptic {
                    if let Some((aya_deg, aya_rate_deg_day)) = sidereal_offset_and_rate {
                        sph_t.lon_deg = (sph_t.lon_deg - aya_deg).rem_euclid(360.0);
                        if options.include_speed {
                            sph_t.lon_speed_deg_day -= aya_rate_deg_day;
                        }
                        let (ps, vs) =
                            crate::frame::spherical_state_to_cart(&sph_t, options.include_speed);
                        p_out = ps;
                        v_out = vs;
                    }
                }
                lon_deg = sph_t.lon_deg;
                lat_deg = sph_t.lat_deg;
                dist_au = sph_t.dist_au;
                lon_speed_deg_day = sph_t.lon_speed_deg_day;
                lat_speed_deg_day = sph_t.lat_speed_deg_day;
                dist_speed_au_day = sph_t.dist_speed_au_day;
                (p_out, v_out)
            }
        };

    Ok(CalcResult {
        body: BodyId::LunarPoint(point),
        target_id: 301,
        center_id: 399,
        target_kind: TargetKind::SyntheticPoint,
        jd_tt,
        jd_ut1: center.ut1.jd_f64(),
        jd_tdb: center.tdb.jd_f64(),
        delta_t_sec: center.delta_t_info.delta_t_sec,
        time_extrapolated: center.delta_t_info.time_extrapolated,
        time_table_version: center.delta_t_info.time_table_version.clone(),
        lon_deg,
        lat_deg,
        dist_au,
        dist_km: dist_au * AU_KM,
        lon_speed_deg_day,
        lat_speed_deg_day,
        dist_speed_au_day,
        dist_speed_km_day: dist_speed_au_day * AU_KM,
        magnitude_v: None,
        cartesian_au: pos_out_au,
        cartesian_vel_au_day: vel_out_au_day,
        warnings: Vec::new(),
    })
}
