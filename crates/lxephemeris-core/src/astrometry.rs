//! Planetary astrometry (`DE440s` + minor bodies), light-time iteration (`MAX_LT_ITER = 10`),
//! multi-body relativistic gravitational light deflection (`sofars::astro::ld` for Sun, Jupiter,
//! Saturn, and Earth), stellar/diurnal aberration (`sofars::astro::ab`), 4-point 4th-order
//! Cartesian stencil velocities, horizontal coordinates (`azalt` with `sofars::astro::refco`),
//! and rise/set/meridian transit root-finding (`rise_trans` via the Illinois method).

use crate::error::EphemerisError;
use crate::frame::{
    cart_state_to_spherical, cart_to_spherical, mat_transpose, mat_vec, stencil_4pt_derivative,
    vec_add, vec_dot, vec_norm, vec_scale, vec_sub, vec_unit, Vec3,
};
use crate::hypothetical::calc_hypothetical_in_moment;
use crate::lunar::calc_lunar_point_in_moment;
use crate::moment::{EpochContext, Moment, H_MOON_DAYS, H_PLANET_DAYS};
use crate::spk::SpkSource;
use crate::stars::calc_star_in_moment;
use crate::time::{TimeInput, TimeTable};
use crate::types::{
    AbMode, BodyId, CalcOptions, Equinox, NutationMode, Observer, Origin, Plane, PlanetId,
    TargetKind, Zodiac, AU_KM, MINOR_BODIES_68, SECONDS_PER_DAY,
};
use serde::{Deserialize, Serialize};
use sofars::astro::{ab, ld, pvtob, refco};

/// Speed of light in vacuum in km/s (IAU 2009 / CODATA).
pub const C_KM_S: f64 = 299_792.458;

/// Speed of light in vacuum in AU/day (`~173.14463267424034 AU/day`).
pub const C_AU_DAY: f64 = (C_KM_S * SECONDS_PER_DAY) / AU_KM;

/// Hard iteration limit for light-time fixed-point iteration (`ARCHITECTURE.md` §2.2.1).
pub const MAX_LT_ITER: u32 = 10;

/// Light-time convergence tolerance in days ($10^{-12}\text{ d} \approx 86.4\text{ ns}$).
pub const LT_CONVERGENCE_DAYS: f64 = 1.0e-12;

/// Solar mass in solar mass units.
pub const BM_SUN: f64 = 1.0;

/// Jupiter system mass in solar mass units (IAU / DE440).
pub const BM_JUPITER: f64 = 1.0 / 1_047.348_644;

/// Saturn system mass in solar mass units (IAU / DE440).
pub const BM_SATURN: f64 = 1.0 / 3_497.898;

/// Earth mass in solar mass units (IAU / DE440).
pub const BM_EARTH: f64 = 1.0 / 332_946.048_7;

/// Complete astrometric calculation result with provenance metadata (`ARCHITECTURE.md` §2.2.3).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CalcResult {
    /// Requested body identifier.
    pub body: BodyId,
    /// NAIF SPK target ID (`10`, `301`, `199`, `299`, `399`, `4..9`, `2000000 + mpc`, etc.).
    pub target_id: i32,
    /// NAIF SPK origin/center ID (`399` for Geocentric/Topocentric, `10` for Heliocentric, `0` for Barycentric).
    pub center_id: i32,
    /// Physical target kind (`BodyCenter`, `SystemBarycenter`, `AnalyticPoint`, `CatalogStar`).
    pub target_kind: TargetKind,
    /// Terrestrial Time (`TT`) Julian Date.
    pub jd_tt: f64,
    /// Universal Time (`UT1`) Julian Date.
    pub jd_ut1: f64,
    /// Barycentric Dynamical Time (`TDB`) Julian Date.
    pub jd_tdb: f64,
    /// $\Delta T = \text{TT} - \text{UT1}$ in seconds.
    pub delta_t_sec: f64,
    /// `true` if $\Delta T$ was evaluated outside empirical IERS/spline coverage.
    pub time_extrapolated: bool,
    /// Active time table version identifier.
    pub time_table_version: String,
    /// Ecliptic longitude, Right Ascension, or Azimuth in degrees `[0, 360)`.
    pub lon_deg: f64,
    /// Ecliptic latitude, Declination, or Altitude in degrees `[-90, +90]`.
    pub lat_deg: f64,
    /// Radial distance in AU.
    pub dist_au: f64,
    /// Radial distance in kilometers (`dist_au * AU_KM`, `R2.2`).
    pub dist_km: f64,
    /// Longitudinal angular velocity in degrees/day (`0.0` if `include_speed == false`).
    pub lon_speed_deg_day: f64,
    /// Latitudinal angular velocity in degrees/day (`0.0` if `include_speed == false`).
    pub lat_speed_deg_day: f64,
    /// Radial velocity in AU/day (`0.0` if `include_speed == false`).
    pub dist_speed_au_day: f64,
    /// Radial velocity in kilometers/day (`dist_speed_au_day * AU_KM`, `R2.2`).
    pub dist_speed_km_day: f64,
    /// Visual magnitude $V$ (populated for fixed stars `R6.1`, `None` for non-stellar points).
    pub magnitude_v: Option<f64>,
    /// Output frame 3D Cartesian position vector `[x, y, z]` in AU.
    pub cartesian_au: Vec3,
    /// Output frame 3D Cartesian velocity vector `[vx, vy, vz]` in AU/day.
    pub cartesian_vel_au_day: Vec3,
    /// Quality warning codes (e.g., `"NonPositiveParallax"` or `"MissingRadialVelocity"` for catalog stars).
    pub warnings: Vec<String>,
}

/// Resolves `(target_id, target_kind, default_stencil_h_days)` for an SPK-backed [`BodyId`].
pub fn resolve_spk_target(body: BodyId) -> Result<(i32, TargetKind, f64), EphemerisError> {
    match body {
        BodyId::Planet(p) => {
            let h = if p == PlanetId::Moon {
                H_MOON_DAYS
            } else {
                H_PLANET_DAYS
            };
            Ok((p.naif_target_id(), p.target_kind(), h))
        }
        BodyId::MinorPlanet(mpc) => {
            let spec = MINOR_BODIES_68
                .binary_search_by_key(&mpc, |s| s.mpc)
                .map(|idx| MINOR_BODIES_68.get(idx))
                .ok()
                .flatten()
                .ok_or_else(|| EphemerisError::UnsupportedBody {
                    body: format!("{body:?}"),
                })?;
            Ok((spec.mpc.standard_naif_spk_id(), TargetKind::MinorBodyCenter, H_PLANET_DAYS))
        }
        other => Err(EphemerisError::UnsupportedBody {
            body: format!("{other:?}"),
        }),
    }
}

/// Resolves the observer's barycentric ICRF state `(pos_km, vel_km_s)` and `center_id` for `origin`.
fn resolve_origin_bary_state_km(ctx: &EpochContext, origin: Origin) -> (Vec3, Vec3, i32) {
    match origin {
        Origin::Geocentric => (
            ctx.earth_bary_km.pos_km,
            ctx.earth_bary_km.vel_km_s,
            399,
        ),
        Origin::Topocentric(obs) => {
            if ctx.observer == Some(obs) {
                (
                    ctx.observer_bary_km.pos_km,
                    ctx.observer_bary_km.vel_km_s,
                    399,
                )
            } else {
                let lon_rad = obs.longitude_deg.to_radians();
                let lat_rad = obs.latitude_deg.to_radians();
                let mut pv_m = [[0.0_f64; 3]; 2];
                pvtob(
                    lon_rad,
                    lat_rad,
                    obs.altitude_m,
                    0.0,
                    0.0,
                    0.0,
                    ctx.gast_rad,
                    &mut pv_m,
                );
                let pos_true_eq_km = [
                    pv_m[0][0] / 1000.0,
                    pv_m[0][1] / 1000.0,
                    pv_m[0][2] / 1000.0,
                ];
                let vel_true_eq_km_s = [
                    pv_m[1][0] / 1000.0,
                    pv_m[1][1] / 1000.0,
                    pv_m[1][2] / 1000.0,
                ];
                let r_bpn_t = mat_transpose(&ctx.frame.r_bpn);
                let geo_pos_km = mat_vec(&r_bpn_t, &pos_true_eq_km);
                let geo_vel_km_s = mat_vec(&r_bpn_t, &vel_true_eq_km_s);
                (
                    vec_add(&ctx.earth_bary_km.pos_km, &geo_pos_km),
                    vec_add(&ctx.earth_bary_km.vel_km_s, &geo_vel_km_s),
                    399,
                )
            }
        }
        Origin::Heliocentric => (ctx.sun_bary_km.pos_km, ctx.sun_bary_km.vel_km_s, 10),
        Origin::Barycentric => ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0),
    }
}

/// Applies relativistic gravitational light deflection (`sofars::astro::ld`) from a single
/// deflecting body of mass `bm` (in solar units) at position `p_deflector_km`.
fn apply_single_deflector(
    u_unit: Vec3,
    bm: f64,
    p_targ_km: &Vec3,
    p_obs_km: &Vec3,
    p_deflector_km: &Vec3,
) -> Vec3 {
    let q_vec = vec_scale(&vec_sub(p_targ_km, p_deflector_km), 1.0 / AU_KM);
    let e_vec = vec_scale(&vec_sub(p_obs_km, p_deflector_km), 1.0 / AU_KM);
    let em = vec_norm(&e_vec);
    let qm = vec_norm(&q_vec);
    if em <= 1e-12 || qm <= 1e-12 {
        return u_unit;
    }
    let q_unit = vec_scale(&q_vec, 1.0 / qm);
    let e_unit = vec_scale(&e_vec, 1.0 / em);
    let deflected = ld(bm, u_unit, q_unit, e_unit, em, 0.0);
    vec_unit(&deflected)
}

/// Applies relativistic gravitational light deflection and stellar aberration to an ICRF astrometric vector `rho_au`.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn apply_aberration_and_deflection_icrf(
    rho_au: &Vec3,
    p_targ_km: &Vec3,
    p_obs_km: &Vec3,
    v_obs_kms: &Vec3,
    ctx: &EpochContext,
    target_id: i32,
    effective_origin: Origin,
    aberration: AbMode,
    tau_days: f64,
) -> Vec3 {
    let dist_au = vec_norm(rho_au);
    if aberration == AbMode::Astrometric || dist_au <= 1e-15 {
        return *rho_au;
    }
    let mut u_unit = vec_scale(rho_au, 1.0 / dist_au);
    let is_earth_based = matches!(effective_origin, Origin::Geocentric | Origin::Topocentric(_));

    if is_earth_based {
        let apply_deflection = matches!(aberration, AbMode::Apparent | AbMode::NoAberration);
        let apply_aberration = matches!(aberration, AbMode::Apparent | AbMode::NoDeflection);

        // 1. Multi-body relativistic gravitational light deflection (Sun, Jupiter, Saturn, +Earth)
        if apply_deflection && target_id != 10 && target_id != 301 {
            let u_initial = u_unit;
            let closest_approach_tau = |p_def_0_km: &Vec3| -> f64 {
                let d_vec_au = vec_scale(&vec_sub(p_def_0_km, p_obs_km), 1.0 / AU_KM);
                let proj_au = vec_dot(&d_vec_au, &u_initial);
                (proj_au / C_AU_DAY).clamp(0.0, tau_days)
            };

            // Sun (NAIF 10) at closest-approach retarded epoch
            let tau_sun = closest_approach_tau(&ctx.sun_bary_km.pos_km);
            let p_sun_ret = vec_sub(
                &ctx.sun_bary_km.pos_km,
                &vec_scale(&ctx.sun_bary_km.vel_km_s, tau_sun * SECONDS_PER_DAY),
            );
            u_unit = apply_single_deflector(u_unit, BM_SUN, p_targ_km, p_obs_km, &p_sun_ret);

            // Jupiter (NAIF 5) at closest-approach retarded epoch
            if target_id != 5 {
                let tau_jup = closest_approach_tau(&ctx.jup_bary_km.pos_km);
                let p_jup_ret = vec_sub(
                    &ctx.jup_bary_km.pos_km,
                    &vec_scale(&ctx.jup_bary_km.vel_km_s, tau_jup * SECONDS_PER_DAY),
                );
                u_unit = apply_single_deflector(u_unit, BM_JUPITER, p_targ_km, p_obs_km, &p_jup_ret);
            }

            // Saturn (NAIF 6) at closest-approach retarded epoch
            if target_id != 6 {
                let tau_sat = closest_approach_tau(&ctx.sat_bary_km.pos_km);
                let p_sat_ret = vec_sub(
                    &ctx.sat_bary_km.pos_km,
                    &vec_scale(&ctx.sat_bary_km.vel_km_s, tau_sat * SECONDS_PER_DAY),
                );
                u_unit = apply_single_deflector(u_unit, BM_SATURN, p_targ_km, p_obs_km, &p_sat_ret);
            }

            // Earth (NAIF 399) when topocentric
            if matches!(effective_origin, Origin::Topocentric(_)) && target_id != 399 {
                u_unit = apply_single_deflector(
                    u_unit,
                    BM_EARTH,
                    p_targ_km,
                    p_obs_km,
                    &ctx.earth_bary_km.pos_km,
                );
            }
        }

        // 2. Relativistic annual + diurnal aberration via sofars::astro::ab
        if apply_aberration {
            let v2_c2 = (vec_dot(v_obs_kms, v_obs_kms) / (C_KM_S * C_KM_S)).min(0.999_999);
            let bm1 = (1.0 - v2_c2).sqrt();
            let v_ab = vec_scale(v_obs_kms, 1.0 / C_KM_S);
            let s_dist_au = vec_norm(&vec_scale(
                &vec_sub(p_obs_km, &ctx.sun_bary_km.pos_km),
                1.0 / AU_KM,
            ));
            u_unit = ab(&u_unit, &v_ab, s_dist_au, bm1);
        }
    }

    vec_scale(&u_unit, dist_au)
}

/// Evaluates the 3D position vector `[x, y, z]` (in AU) of `target_id` in the output coordinate
/// frame specified by `options` at a single [`EpochContext`].
pub fn eval_spk_target_vector_au(
    spk: &dyn SpkSource,
    ctx: &EpochContext,
    target_id: i32,
    options: &CalcOptions,
) -> Result<Vec3, EphemerisError> {
    eval_spk_target_vector_au_seeded(spk, ctx, target_id, options, None).map(|(v, _)| v)
}

/// Evaluates the 3D position vector `[x, y, z]` (in AU) and converged light-time `tau_days` of
/// `target_id` at a single [`EpochContext`], optionally seeded with `seed_tau_days` from a nearby
/// stencil center to avoid a redundant `tau = 0` SPK lookup.
pub fn eval_spk_target_vector_au_seeded(
    spk: &dyn SpkSource,
    ctx: &EpochContext,
    target_id: i32,
    options: &CalcOptions,
    seed_tau_days: Option<f64>,
) -> Result<(Vec3, f64), EphemerisError> {
    let effective_origin = if options.plane == Plane::Horizontal
        && matches!(options.origin, Origin::Geocentric)
    {
        ctx.observer.map_or(Origin::Geocentric, Origin::Topocentric)
    } else {
        options.origin
    };
    let (p_obs_km, v_obs_kms, center_id) = resolve_origin_bary_state_km(ctx, effective_origin);
    if target_id == center_id && !matches!(effective_origin, Origin::Topocentric(_)) {
        return Ok(([0.0, 0.0, 0.0], 0.0));
    }

    let tdb_0 = ctx.tdb.0;
    let cached_st_0 = match target_id {
        10 => Some(ctx.sun_bary_km),
        5 => Some(ctx.jup_bary_km),
        6 => Some(ctx.sat_bary_km),
        399 => Some(ctx.earth_bary_km),
        _ => None,
    };

    let mut converged_tau_days = 0.0_f64;
    let vec_icrf_au = match options.aberration {
        AbMode::Geometric => {
            let st_targ = match cached_st_0 {
                Some(st) => st,
                None => spk.state_icrf_km_split(target_id, 0, tdb_0.jd1, tdb_0.jd2)?,
            };
            vec_scale(&vec_sub(&st_targ.pos_km, &p_obs_km), 1.0 / AU_KM)
        }
        AbMode::Astrometric | AbMode::NoAberration | AbMode::NoDeflection | AbMode::Apparent => {
            // Seed light-time from seed_tau_days (when available) or from tdb_0 state + velocity predictor
            let mut tau_days = if let Some(seed_tau) = seed_tau_days {
                seed_tau
            } else {
                let st_0 = match cached_st_0 {
                    Some(st) => st,
                    None => spk.state_icrf_km_split(target_id, 0, tdb_0.jd1, tdb_0.jd2)?,
                };
                let rho_0_au = vec_scale(&vec_sub(&st_0.pos_km, &p_obs_km), 1.0 / AU_KM);
                let tau_0_days = vec_norm(&rho_0_au) / C_AU_DAY;
                let p_pred_km = vec_sub(
                    &st_0.pos_km,
                    &vec_scale(&st_0.vel_km_s, tau_0_days * SECONDS_PER_DAY),
                );
                vec_norm(&vec_scale(&vec_sub(&p_pred_km, &p_obs_km), 1.0 / AU_KM)) / C_AU_DAY
            };

            // Fixed-point light-time iteration with hard MAX_LT_ITER = 10 limit
            let mut converged = false;
            let mut p_targ_km = [0.0_f64; 3];
            let mut rho_au = [0.0_f64; 3];

            let mut last_residual = 0.0_f64;
            for _ in 0..MAX_LT_ITER {
                let st_targ =
                    spk.state_icrf_km_split(target_id, 0, tdb_0.jd1, tdb_0.jd2 - tau_days)?;
                p_targ_km = st_targ.pos_km;
                rho_au = vec_scale(&vec_sub(&p_targ_km, &p_obs_km), 1.0 / AU_KM);
                let dist_au = vec_norm(&rho_au);
                let new_tau = dist_au / C_AU_DAY;
                let dtau = new_tau - tau_days;
                last_residual = dtau.abs();
                if last_residual < LT_CONVERGENCE_DAYS {
                    tau_days = new_tau;
                    converged = true;
                    break;
                }
                // Apply first-order velocity correction for remaining sub-millisecond dtau
                let p_corr_km = vec_sub(
                    &st_targ.pos_km,
                    &vec_scale(&st_targ.vel_km_s, dtau * SECONDS_PER_DAY),
                );
                let rho_corr_au = vec_scale(&vec_sub(&p_corr_km, &p_obs_km), 1.0 / AU_KM);
                let dist_corr_au = vec_norm(&rho_corr_au);
                let tau_corr = dist_corr_au / C_AU_DAY;
                if (tau_corr - new_tau).abs() < LT_CONVERGENCE_DAYS {
                    p_targ_km = p_corr_km;
                    rho_au = rho_corr_au;
                    tau_days = tau_corr;
                    converged = true;
                    break;
                }
                tau_days = tau_corr;
            }
            converged_tau_days = tau_days;

            if !converged {
                return Err(EphemerisError::NonConvergence {
                    solver: "light_time_iteration".to_owned(),
                    iterations: MAX_LT_ITER,
                    residual: last_residual,
                });
            }

            apply_aberration_and_deflection_icrf(
                &rho_au,
                &p_targ_km,
                &p_obs_km,
                &v_obs_kms,
                ctx,
                target_id,
                effective_origin,
                options.aberration,
                tau_days,
            )
        }
    };

    let r_mat = ctx
        .frame
        .icrf_to_frame_matrix(options.plane, options.equinox, options.nutation);
    let vec_frame_au = mat_vec(&r_mat, &vec_icrf_au);
    if options.plane == Plane::Horizontal {
        let obs = match options.origin {
            Origin::Topocentric(o) => o,
            _ => ctx.observer.ok_or_else(|| EphemerisError::InvalidInput {
                message: "Plane::Horizontal requires a topocentric Observer".to_owned(),
            })?,
        };
        Ok((
            equatorial_to_horizontal_vector_au(
                obs,
                ctx.gast_rad,
                &vec_frame_au,
                AtmosphereParams::default(),
            ),
            converged_tau_days,
        ))
    } else {
        Ok((vec_frame_au, converged_tau_days))
    }
}

/// Transforms a 3D vector in True Equatorial of Date (`pos_eq_au`) into a 3D Horizontal vector
/// `[x_north, y_east, z_up]` (in AU) whose spherical coordinates are `(azimuth_deg, altitude_apparent_deg, dist_au)`.
#[must_use]
pub fn equatorial_to_horizontal_vector_au(
    observer: Observer,
    gast_rad: f64,
    pos_eq_au: &Vec3,
    atmosphere: AtmosphereParams,
) -> Vec3 {
    let sph = cart_to_spherical(pos_eq_au);
    let hz = azalt(
        observer,
        gast_rad,
        sph.lon_deg,
        sph.lat_deg,
        atmosphere,
    );
    let az_rad = hz.azimuth_deg.to_radians();
    let alt_rad = hz.altitude_apparent_deg.to_radians();
    let (sin_az, cos_az) = az_rad.sin_cos();
    let (sin_alt, cos_alt) = alt_rad.sin_cos();
    [
        sph.dist * cos_alt * cos_az,
        sph.dist * cos_alt * sin_az,
        sph.dist * sin_alt,
    ]
}

/// Evaluates an SPK-backed body (`Planet` or `MinorPlanet`) inside a [`Moment`], using the
/// 4-point 4th-order Cartesian stencil when `options.include_speed` is `true`.
pub fn calc_spk_body_in_moment(
    moment: &Moment<'_>,
    body: BodyId,
    options: &CalcOptions,
) -> Result<CalcResult, EphemerisError> {
    let (target_id, target_kind, h_days) = resolve_spk_target(body)?;
    let center = moment.center();
    let (p_obs_km, _, center_id) = resolve_origin_bary_state_km(center, options.origin);

    let (pos_0, tau_0) =
        eval_spk_target_vector_au_seeded(moment.spk(), center, target_id, options, None)?;

    let (vel_0, mut sph) = if options.include_speed {
        let stencil = if (h_days - H_MOON_DAYS).abs() < 1e-12 {
            moment.moon_stencil_contexts()?
        } else {
            moment.planet_stencil_contexts()?
        };
        let mut stencil_opts = *options;
        if options.plane == Plane::Ecliptic && matches!(options.zodiac, Zodiac::Sidereal(_)) {
            stencil_opts.nutation = NutationMode::Mean;
        }
        let seed = (tau_0 > 0.0).then_some(tau_0);
        let (pos_m2, _) = eval_spk_target_vector_au_seeded(
            moment.spk(),
            &stencil[0],
            target_id,
            &stencil_opts,
            seed,
        )?;
        let (pos_m1, _) = eval_spk_target_vector_au_seeded(
            moment.spk(),
            &stencil[1],
            target_id,
            &stencil_opts,
            seed,
        )?;
        let (pos_p1, _) = eval_spk_target_vector_au_seeded(
            moment.spk(),
            &stencil[2],
            target_id,
            &stencil_opts,
            seed,
        )?;
        let (pos_p2, _) = eval_spk_target_vector_au_seeded(
            moment.spk(),
            &stencil[3],
            target_id,
            &stencil_opts,
            seed,
        )?;
        let v0 = stencil_4pt_derivative(&pos_m2, &pos_m1, &pos_p1, &pos_p2, h_days);
        let pos_for_rate = if stencil_opts.nutation != options.nutation {
            eval_spk_target_vector_au_seeded(moment.spk(), center, target_id, &stencil_opts, seed)?
                .0
        } else {
            pos_0
        };
        let mut s = cart_state_to_spherical(&pos_for_rate, &v0);
        if stencil_opts.nutation != options.nutation {
            let c0 = cart_to_spherical(&pos_0);
            s.lon_deg = c0.lon_deg;
            s.lat_deg = c0.lat_deg;
            s.dist_au = c0.dist;
        }
        (v0, s)
    } else {
        let c = cart_to_spherical(&pos_0);
        (
            [0.0, 0.0, 0.0],
            crate::frame::SphericalState {
                lon_deg: c.lon_deg,
                lat_deg: c.lat_deg,
                dist_au: c.dist,
                lon_speed_deg_day: 0.0,
                lat_speed_deg_day: 0.0,
                dist_speed_au_day: 0.0,
            },
        )
    };

    let mut cart_pos = pos_0;
    let mut cart_vel = vel_0;
    if let Some((aya_deg, aya_rate_deg_day)) = moment.sidereal_offset_and_rate(options)? {
        sph.lon_deg = (sph.lon_deg - aya_deg).rem_euclid(360.0);
        if options.include_speed {
            sph.lon_speed_deg_day -= aya_rate_deg_day;
        }
        let (p_sid, v_sid) = crate::frame::spherical_state_to_cart(&sph, options.include_speed);
        cart_pos = p_sid;
        cart_vel = v_sid;
    }

    let magnitude_v = if let Ok(st_t) = moment.spk().state_icrf_km_split(
        target_id,
        0,
        center.tdb.0.jd1,
        center.tdb.0.jd2 - tau_0,
    ) {
        let r_helio_au = vec_scale(
            &vec_sub(&st_t.pos_km, &center.sun_bary_km.pos_km),
            1.0 / AU_KM,
        );
        let delta_obs_au = vec_scale(&vec_sub(&st_t.pos_km, &p_obs_km), 1.0 / AU_KM);
        let r_obs_sun_au = vec_scale(&vec_sub(&center.sun_bary_km.pos_km, &p_obs_km), 1.0 / AU_KM);
        crate::phenomena::compute_apparent_magnitude_v(
            body,
            &r_helio_au,
            &delta_obs_au,
            &r_obs_sun_au,
        )
    } else {
        None
    };

    Ok(CalcResult {
        body,
        target_id,
        center_id,
        target_kind,
        jd_tt: center.tt.jd_f64(),
        jd_ut1: center.ut1.jd_f64(),
        jd_tdb: center.tdb.jd_f64(),
        delta_t_sec: center.delta_t_info.delta_t_sec,
        time_extrapolated: center.delta_t_info.time_extrapolated,
        time_table_version: center.delta_t_info.time_table_version.clone(),
        lon_deg: sph.lon_deg,
        lat_deg: sph.lat_deg,
        dist_au: sph.dist_au,
        dist_km: sph.dist_au * AU_KM,
        lon_speed_deg_day: sph.lon_speed_deg_day,
        lat_speed_deg_day: sph.lat_speed_deg_day,
        dist_speed_au_day: sph.dist_speed_au_day,
        dist_speed_km_day: sph.dist_speed_au_day * AU_KM,
        magnitude_v,
        cartesian_au: cart_pos,
        cartesian_vel_au_day: cart_vel,
        warnings: Vec::new(),
    })
}

/// Evaluates any [`BodyId`] (`Planet`, `MinorPlanet`, `LunarPoint`, `Hypothetical`, `FixedStar`, or `EclNut`)
/// inside a [`Moment`].
pub fn calc_body_in_moment(
    moment: &Moment<'_>,
    body: BodyId,
    options: &CalcOptions,
) -> Result<CalcResult, EphemerisError> {
    match body {
        BodyId::Planet(_) | BodyId::MinorPlanet(_) => calc_spk_body_in_moment(moment, body, options),
        BodyId::LunarPoint(lp) => {
            let sid = moment.sidereal_offset_and_rate(options)?;
            calc_lunar_point_in_moment(moment, lp, options, sid)
        }
        BodyId::Hypothetical(h) => {
            let sid = moment.sidereal_offset_and_rate(options)?;
            calc_hypothetical_in_moment(moment, h, options, sid)
        }
        BodyId::FixedStar(hip) => {
            let sid = moment.sidereal_offset_and_rate(options)?;
            calc_star_in_moment(moment, hip, options, sid)
        }
        BodyId::EclNut => {
            let center = moment.center();
            let info = center.frame.ecl_nut_info();
            Ok(CalcResult {
                body: BodyId::EclNut,
                target_id: -1,
                center_id: 399,
                target_kind: TargetKind::SyntheticPoint,
                jd_tt: center.tt.jd_f64(),
                jd_ut1: center.ut1.jd_f64(),
                jd_tdb: center.tdb.jd_f64(),
                delta_t_sec: center.delta_t_info.delta_t_sec,
                time_extrapolated: center.delta_t_info.time_extrapolated,
                time_table_version: center.delta_t_info.time_table_version.clone(),
                lon_deg: info.true_obliquity_deg,
                lat_deg: info.mean_obliquity_deg,
                dist_au: info.nutation_lon_deg,
                dist_km: 0.0,
                lon_speed_deg_day: info.nutation_obl_deg,
                lat_speed_deg_day: 0.0,
                dist_speed_au_day: 0.0,
                dist_speed_km_day: 0.0,
                magnitude_v: None,
                cartesian_au: [
                    info.true_obliquity_deg,
                    info.mean_obliquity_deg,
                    info.nutation_lon_deg,
                ],
                cartesian_vel_au_day: [info.nutation_obl_deg, 0.0, 0.0],
                warnings: Vec::new(),
            })
        }
    }
}

/// Atmospheric parameters for `sofars::astro::refco` optical refraction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AtmosphereParams {
    /// Ambient atmospheric pressure at observer in hPa / mbar (`0.0` disables refraction).
    pub pressure_hpa: f64,
    /// Ambient temperature at observer in degrees Celsius.
    pub temperature_c: f64,
    /// Relative humidity in `[0.0, 1.0]`.
    pub relative_humidity: f64,
    /// Effective optical wavelength in micrometers (`0.55` µm standard visual).
    pub wavelength_um: f64,
}

impl Default for AtmosphereParams {
    fn default() -> Self {
        Self {
            pressure_hpa: 1013.25,
            temperature_c: 15.0,
            relative_humidity: 0.5,
            wavelength_um: 0.55,
        }
    }
}

impl AtmosphereParams {
    /// Airless (unrefracted geometric horizon) atmospheric parameters.
    #[must_use]
    pub const fn airless() -> Self {
        Self {
            pressure_hpa: 0.0,
            temperature_c: 15.0,
            relative_humidity: 0.0,
            wavelength_um: 0.55,
        }
    }
}

/// Topocentric horizontal coordinates (Azimuth East of North `[0, 360)`, true and apparent Altitude).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HorizontalCoords {
    /// Azimuth in degrees `[0, 360)` measured East from North (`0° = N, 90° = E, 180° = S, 270° = W`).
    pub azimuth_deg: f64,
    /// Geometric (airless) altitude above the horizon in degrees `[-90, +90]`.
    pub altitude_true_deg: f64,
    /// Apparent (refracted) altitude above the horizon in degrees `[-90, +90]`.
    pub altitude_apparent_deg: f64,
    /// Local Hour Angle $H$ in degrees `[-180, +180)`.
    pub hour_angle_deg: f64,
    /// `true` if `|altitude_true_deg| < 5.0°` (sect boundary proximity flag `R7.3`).
    pub near_horizon: bool,
}

/// Transforms apparent topocentric True Equatorial of Date `(ra_deg, dec_deg)` into
/// [`HorizontalCoords`] at `gast_rad` for `observer`, applying IAU SOFA refraction (`refco`).
#[must_use]
pub fn azalt(
    observer: Observer,
    gast_rad: f64,
    ra_deg: f64,
    dec_deg: f64,
    atmosphere: AtmosphereParams,
) -> HorizontalCoords {
    let lon_rad = observer.longitude_deg.to_radians();
    let phi = observer.latitude_deg.to_radians();
    let ra_rad = ra_deg.to_radians();
    let dec_rad = dec_deg.to_radians();

    let last_rad = gast_rad + lon_rad;
    let h_rad = (last_rad - ra_rad + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
        - std::f64::consts::PI;

    let (sin_phi, cos_phi) = phi.sin_cos();
    let (sin_dec, cos_dec) = dec_rad.sin_cos();
    let (sin_h, cos_h) = h_rad.sin_cos();

    let sin_alt = (sin_phi * sin_dec + cos_phi * cos_dec * cos_h).clamp(-1.0, 1.0);
    let alt_true_rad = sin_alt.asin();

    let y = -cos_dec * sin_h;
    let x = cos_phi * sin_dec - sin_phi * cos_dec * cos_h;
    let az_deg = y.atan2(x).to_degrees().rem_euclid(360.0);

    let altitude_true_deg = alt_true_rad.to_degrees();
    let alt_app_rad = if atmosphere.pressure_hpa > 0.0 && altitude_true_deg > -5.0 {
        let (refa, refb) = refco(
            atmosphere.pressure_hpa,
            atmosphere.temperature_c,
            atmosphere.relative_humidity,
            atmosphere.wavelength_um,
        );
        // For alt >= 15° (z <= 75°), use IAU SOFA refco tan(z) + tan^3(z);
        // for low altitudes (-5° < alt < 15°), blend smoothly with Saemundsson (1986) / Bennett (1982)
        // low-altitude refraction so horizon refraction is ~29'–34' instead of breaking down at tan^3(z).
        let z_clamped = (std::f64::consts::FRAC_PI_2 - alt_true_rad)
            .clamp(0.0, 75.0_f64.to_radians());
        let tz = z_clamped.tan();
        let dz_sofa = refa * tz + refb * tz * tz * tz;

        let dz = if altitude_true_deg >= 15.0 {
            dz_sofa
        } else {
            let h_deg = altitude_true_deg.max(-4.5);
            let arg_rad = (h_deg + 10.3 / (h_deg + 5.11)).to_radians();
            let p_scale = atmosphere.pressure_hpa / 1010.0;
            let t_scale = 283.15 / (273.15 + atmosphere.temperature_c).max(100.0);
            let dz_low = ((1.02 / arg_rad.tan()) * p_scale * t_scale / 60.0).to_radians();
            if altitude_true_deg <= 10.0 {
                dz_low
            } else {
                // Smooth linear blend on [10°, 15°]
                let w = (altitude_true_deg - 10.0) / 5.0;
                (1.0 - w) * dz_low + w * dz_sofa
            }
        };
        (alt_true_rad + dz).min(std::f64::consts::FRAC_PI_2)
    } else {
        alt_true_rad
    };

    HorizontalCoords {
        azimuth_deg: az_deg,
        altitude_true_deg,
        altitude_apparent_deg: alt_app_rad.to_degrees(),
        hour_angle_deg: h_rad.to_degrees(),
        near_horizon: altitude_true_deg.abs() < 5.0,
    }
}

/// Transforms apparent topocentric True Ecliptic of Date `(ecl_lon_deg, ecl_lat_deg)` into
/// [`HorizontalCoords`] at `gast_rad` and `true_obliquity_rad` for `observer` (`R2.6 & R7.3`).
#[must_use]
pub fn azalt_from_ecliptic(
    observer: Observer,
    gast_rad: f64,
    true_obliquity_rad: f64,
    ecl_lon_deg: f64,
    ecl_lat_deg: f64,
    atmosphere: AtmosphereParams,
) -> HorizontalCoords {
    let (sin_l, cos_l) = ecl_lon_deg.to_radians().sin_cos();
    let (sin_b, cos_b) = ecl_lat_deg.to_radians().sin_cos();
    let (sin_e, cos_e) = true_obliquity_rad.sin_cos();
    let u_ecl = [cos_b * cos_l, cos_b * sin_l, sin_b];
    let u_eq = [
        u_ecl[0],
        u_ecl[1] * cos_e - u_ecl[2] * sin_e,
        u_ecl[1] * sin_e + u_ecl[2] * cos_e,
    ];
    let sph_eq = cart_to_spherical(&u_eq);
    azalt(
        observer,
        gast_rad,
        sph_eq.lon_deg,
        sph_eq.lat_deg,
        atmosphere,
    )
}

/// Inverse horizontal coordinate transformation (`swe_azalt_rev`): converts geometric
/// (`airless`) topocentric `(azimuth_north_east_deg, altitude_true_deg)` back into True Equatorial
/// `(ra_deg, dec_deg)` (`target_plane == Plane::Equatorial`) or True Ecliptic `(lon_deg, lat_deg)`
/// (`target_plane == Plane::Ecliptic`) at `gast_rad` and `true_obliquity_rad` for `observer`.
#[must_use]
pub fn azalt_rev(
    observer: Observer,
    gast_rad: f64,
    true_obliquity_rad: f64,
    azimuth_north_east_deg: f64,
    altitude_true_deg: f64,
    target_plane: Plane,
) -> [f64; 2] {
    let phi = observer.latitude_deg.to_radians();
    let lon_rad = observer.longitude_deg.to_radians();
    let az_rad = azimuth_north_east_deg.to_radians();
    let alt_rad = altitude_true_deg.clamp(-90.0, 90.0).to_radians();

    let (sin_phi, cos_phi) = phi.sin_cos();
    let (sin_az, cos_az) = az_rad.sin_cos();
    let (sin_alt, cos_alt) = alt_rad.sin_cos();

    let sin_dec = (sin_phi * sin_alt + cos_phi * cos_alt * cos_az).clamp(-1.0, 1.0);
    let dec_rad = sin_dec.asin();

    let sin_h = -cos_alt * sin_az;
    let cos_h = cos_phi * sin_alt - sin_phi * cos_alt * cos_az;
    let h_rad = sin_h.atan2(cos_h);

    let last_rad = gast_rad + lon_rad;
    let ra_rad = (last_rad - h_rad).rem_euclid(std::f64::consts::TAU);

    if target_plane == Plane::Equatorial {
        return [ra_rad.to_degrees().rem_euclid(360.0), dec_rad.to_degrees()];
    }

    let (sin_ra, cos_ra) = ra_rad.sin_cos();
    let cos_dec = dec_rad.cos();
    let u_eq = [cos_dec * cos_ra, cos_dec * sin_ra, sin_dec];
    let (sin_e, cos_e) = true_obliquity_rad.sin_cos();
    let u_ecl = [
        u_eq[0],
        u_eq[1] * cos_e + u_eq[2] * sin_e,
        -u_eq[1] * sin_e + u_eq[2] * cos_e,
    ];
    let sph_ecl = cart_to_spherical(&u_ecl);
    [sph_ecl.lon_deg, sph_ecl.lat_deg]
}

/// Evaluates the Equation of Time ($E = \text{Apparent Solar Time} - \text{Mean Solar Time}$)
/// in **days** (`[-0.02, +0.02]`), matching `swe_time_equ` (Meeus 1998 Ch. 28).
pub fn time_equation_days(moment: &Moment<'_>) -> Result<f64, EphemerisError> {
    let sun_eq_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };
    let sun_eq = calc_body_in_moment(moment, BodyId::Planet(PlanetId::Sun), &sun_eq_opts)?;
    let alpha_app_deg = sun_eq.lon_deg;

    let center = moment.center();
    let jd_tt = center.tt.jd_f64();
    let tau = (jd_tt - crate::types::JD_J2000) / 365_250.0;
    // Sun's mean longitude L0 (Simon et al. 1994 / Meeus Eq. 28.2, tau in Julian millennia)
    let l0_deg = (280.466_456_7
        + 360_007.698_277_9 * tau
        + 0.030_320_28 * tau * tau
        + (tau * tau * tau) / 49_931.0
        - (tau * tau * tau * tau) / 15_300.0
        - (tau * tau * tau * tau * tau) / 2_000_000.0)
        .rem_euclid(360.0);
    let dpsi_deg = center.frame.dpsi_rad.to_degrees();
    let cos_eps = center.frame.true_eps_rad.cos();
    let alpha_mean_deg = (l0_deg - 0.005_718_3 + dpsi_deg * cos_eps).rem_euclid(360.0);

    let diff_deg = (alpha_mean_deg - alpha_app_deg + 180.0).rem_euclid(360.0) - 180.0;
    Ok(diff_deg / 360.0)
}

/// Solar/Lunar/Planetary disc reference point for horizon crossings (`R6.2 & R7.1`).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum DiscMode {
    /// Upper limb of Sun/Moon (`-16'` semi-diameter offset), geometric center for point sources (default).
    #[default]
    UpperLimb,
    /// Geometric center of the body's disc (`0'` semi-diameter offset; used for astrological parans `R6.2`).
    Center,
    /// Bottom/lower limb of Sun/Moon (`+16'` semi-diameter offset).
    BottomLimb,
}

/// Horizon event type for [`rise_trans`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RiseTransEvent {
    /// Upper limb or center rises across the horizon ($\dot{a} > 0$).
    Rise,
    /// Upper limb or center sets across the horizon ($\dot{a} < 0$).
    Set,
    /// Upper meridian culmination ($H = 0$).
    UpperTransit,
    /// Lower meridian culmination ($H = \pm 180^\circ$).
    LowerTransit,
}

/// Options for [`rise_trans`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RiseTransOptions {
    /// Maximum search window in days (default `2.0` days).
    pub max_search_days: f64,
    /// Target true (airless) horizon altitude in degrees for `Rise`/`Set`
    /// (`None` computes horizon altitude from `disc_mode`, `pressure_hpa`, and `temperature_c`).
    pub custom_horizon_alt_deg: Option<f64>,
    /// Disc reference mode (`UpperLimb` default for Sun/Moon, `Center` for parans `R6.2`).
    pub disc_mode: DiscMode,
    /// Atmospheric pressure at observer in hPa (`1013.25` default for `-34'` refraction, `0.0` disables refraction).
    pub pressure_hpa: f64,
    /// Ambient temperature at observer in degrees Celsius (`15.0` default).
    pub temperature_c: f64,
}

impl Default for RiseTransOptions {
    fn default() -> Self {
        Self::almanac_upper_limb()
    }
}

impl RiseTransOptions {
    /// Standard almanac preset (`R7.1`): `UpperLimb` disc reference (true distance-dependent
    /// apparent semi-diameter for Sun/Moon) with standard atmospheric horizon refraction
    /// (`1013.25 hPa`, `15.0 °C` $\implies -0^\circ 34'$ refraction).
    #[must_use]
    pub const fn almanac_upper_limb() -> Self {
        Self {
            max_search_days: 2.0,
            custom_horizon_alt_deg: None,
            disc_mode: DiscMode::UpperLimb,
            pressure_hpa: 1013.25,
            temperature_c: 15.0,
        }
    }

    /// Geometric center preset (`R6.2` Brady diurnal parans): airless (`0.0 hPa`) geometric
    /// center crossing (`DiscMode::Center`) at $0^\circ$ true altitude.
    #[must_use]
    pub const fn geometric_center() -> Self {
        Self {
            max_search_days: 2.0,
            custom_horizon_alt_deg: None,
            disc_mode: DiscMode::Center,
            pressure_hpa: 0.0,
            temperature_c: 15.0,
        }
    }
}

/// Finds the next occurrence (`JD_UT1`) of `event` (`Rise`, `Set`, `UpperTransit`, `LowerTransit`)
/// after `jd_ut1_start` for `body` at `observer`, refining the root with the Illinois method to
/// $< 10^{-8}\text{ d}$ ($\approx 0.86\text{ ms}$). Returns `Ok(None)` if the body is circumpolar
/// (does not rise or set within `options.max_search_days`).
pub fn rise_trans(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    body: BodyId,
    observer: Observer,
    jd_ut1_start: f64,
    event: RiseTransEvent,
    options: RiseTransOptions,
) -> Result<Option<f64>, EphemerisError> {
    if !jd_ut1_start.is_finite()
        || !options.max_search_days.is_finite()
        || !(0.01..=36_525.0).contains(&options.max_search_days)
        || !options.pressure_hpa.is_finite()
        || !options.temperature_c.is_finite()
        || options
            .custom_horizon_alt_deg
            .is_some_and(|a| !a.is_finite() || !(-90.0..=90.0).contains(&a))
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid rise_trans options: jd_ut1_start={jd_ut1_start}, max_search_days={}",
                options.max_search_days
            ),
        });
    }
    let ref_deg = if options.pressure_hpa <= 0.0 {
        0.0
    } else {
        let p_factor = options.pressure_hpa / 1013.25;
        let t_factor = 288.15 / (273.15 + options.temperature_c).max(100.0);
        (34.0 / 60.0) * p_factor * t_factor
    };

    let calc_opts = CalcOptions {
        origin: Origin::Topocentric(observer),
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: crate::types::Zodiac::Tropical,
        include_speed: false,
    };

    let eval_at_ut1 = |jd_ut1: f64| -> Result<(f64, f64), EphemerisError> {
        let m = Moment::new(
            spk,
            time_table,
            TimeInput::from_ut1_jd(jd_ut1),
            Some(observer),
        )?;
        let res = calc_body_in_moment(&m, body, &calc_opts)?;
        let target_alt_deg = options.custom_horizon_alt_deg.unwrap_or_else(|| {
            // Compute true distance-dependent apparent semi-diameter from topocentric distance res.dist_au
            let dist_km = (res.dist_au * AU_KM).max(1_000.0);
            let sd_deg = match body {
                BodyId::Planet(PlanetId::Sun) => {
                    let s_deg = (695_700.0_f64 / dist_km).clamp(0.0, 1.0).asin().to_degrees();
                    match options.disc_mode {
                        DiscMode::UpperLimb => s_deg,
                        DiscMode::Center => 0.0,
                        DiscMode::BottomLimb => -s_deg,
                    }
                }
                BodyId::Planet(PlanetId::Moon) => {
                    let s_deg = (1_737.4_f64 / dist_km).clamp(0.0, 1.0).asin().to_degrees();
                    match options.disc_mode {
                        DiscMode::UpperLimb => s_deg,
                        DiscMode::Center => 0.0,
                        DiscMode::BottomLimb => -s_deg,
                    }
                }
                _ => 0.0,
            };
            -(ref_deg + sd_deg)
        });
        let hz = azalt(
            observer,
            m.center().gast_rad,
            res.lon_deg,
            res.lat_deg,
            AtmosphereParams::airless(),
        );
        Ok((hz.altitude_true_deg - target_alt_deg, hz.hour_angle_deg))
    };

    // Coarse step of 0.5 hours (1/48 day) across max_search_days
    let step_days = 1.0 / 48.0;
    let total_steps = ((options.max_search_days / step_days).ceil() as usize).max(2);

    let mut t_prev = jd_ut1_start;
    let (mut alt_prev, mut ha_prev) = eval_at_ut1(t_prev)?;

    for step in 1..=total_steps {
        let t_curr = (jd_ut1_start + (step as f64) * step_days)
            .min(jd_ut1_start + options.max_search_days);
        let (alt_curr, ha_curr) = eval_at_ut1(t_curr)?;

        let bracket_found = match event {
            RiseTransEvent::Rise => alt_prev < 0.0 && alt_curr >= 0.0,
            RiseTransEvent::Set => alt_prev > 0.0 && alt_curr <= 0.0,
            RiseTransEvent::UpperTransit => {
                ha_prev < 0.0 && ha_curr >= 0.0 && (ha_curr - ha_prev).abs() < 90.0
            }
            RiseTransEvent::LowerTransit => {
                // Shift hour angle by 180 deg into [-180, 180)
                let l_prev = (ha_prev + 360.0).rem_euclid(360.0) - 180.0;
                let l_curr = (ha_curr + 360.0).rem_euclid(360.0) - 180.0;
                l_prev < 0.0 && l_curr >= 0.0 && (l_curr - l_prev).abs() < 90.0
            }
        };

        if bracket_found {
            // Refine bracket [t_prev, t_curr] via the Illinois method
            let f_scalar = |t: f64| -> Result<f64, EphemerisError> {
                let (alt_diff, ha) = eval_at_ut1(t)?;
                Ok(match event {
                    RiseTransEvent::Rise | RiseTransEvent::Set => alt_diff,
                    RiseTransEvent::UpperTransit => ha,
                    RiseTransEvent::LowerTransit => (ha + 360.0).rem_euclid(360.0) - 180.0,
                })
            };

            let mut a = t_prev;
            let mut b = t_curr;
            let mut fa = f_scalar(a)?;
            let mut fb = f_scalar(b)?;
            let mut side = 0i8;

            for _ in 0..50 {
                if (b - a).abs() < 1e-8 || fb.abs() < 1e-9 {
                    return Ok(Some(b));
                }
                let denom = fb - fa;
                let c = if denom.abs() > 1e-15 {
                    (a * fb - b * fa) / denom
                } else {
                    0.5 * (a + b)
                };
                let fc = f_scalar(c)?;
                if fc.abs() < 1e-9 {
                    return Ok(Some(c));
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

            return Err(EphemerisError::NonConvergence {
                solver: "rise_trans_illinois".to_owned(),
                iterations: 50,
                residual: (b - a).abs(),
            });
        }

        t_prev = t_curr;
        alt_prev = alt_curr;
        ha_prev = ha_curr;
    }

    Ok(None)
}
