//! 8 Witte / Neely Uranian hypothetical bodies (`Cupido..Poseidon`) and `Transpluto` (`Isis`),
//! strictly implementing `MODEL_PASSPORTS.md` §3:
//! - James Neely (1980) osculating Keplerian elements at `J1900.0` (`JD 2415020.0`) for the 8
//!   Hamburg School (Alfred Witte / Friedrich Sieggrün) Uranians.
//! - M. E. Sevin / W. Strubell (1952) Keplerian elements for `Transpluto` (`Isis`).
//! - Unperturbed heliocentric two-body Keplerian motion rotated from source ecliptic (`equinox_jd`)
//!   to ICRF (`bp06` + `obl06`), combined with `DE440s` Earth/Sun state, 3-step light-time
//!   iteration, relativistic stellar aberration (`sofars::astro::ab`), and frame rotation.

use crate::astrometry::{
    equatorial_to_horizontal_vector_au, AtmosphereParams, CalcResult, C_AU_DAY, C_KM_S,
};
use crate::error::EphemerisError;
use crate::frame::{
    cart_state_to_spherical, cart_to_spherical, mat_mul, mat_transpose, mat_vec, rot_x,
    stencil_4pt_derivative, vec_dot, vec_norm, vec_scale, vec_sub, Vec3,
};
use crate::moment::{EpochContext, Moment, H_PLANET_DAYS};
use crate::types::{
    AbMode, BodyId, CalcOptions, HypotheticalId, Origin, Plane, TargetKind, AU_KM,
};
use serde::{Deserialize, Serialize};
use sofars::astro::ab;
use sofars::pnp::{bp06, obl06};

/// Daily Gaussian gravitational constant $k^\circ = 0.9856076686^\circ/\text{day}$ (`MODEL_PASSPORTS.md` §3.2).
pub const GAUSS_K_DEG_DAY: f64 = 0.985_607_668_6;

/// Canonical Keplerian orbital elements for a hypothetical body (`MODEL_PASSPORTS.md` §3.1).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HypotheticalElements {
    pub id: HypotheticalId,
    /// Element reference epoch $t_0$ (`JD_TT`).
    pub epoch_jd: f64,
    /// Reference ecliptic equinox `JD_TT` in which $(\omega, \Omega, i)$ are defined.
    pub equinox_jd: f64,
    /// Mean anomaly $M_0$ at $t_0$ in degrees.
    pub m0_deg: f64,
    /// Semi-major axis $a$ in AU.
    pub a_au: f64,
    /// Orbital eccentricity $e$.
    pub e: f64,
    /// Argument of perihelion $\omega$ in degrees.
    pub w_deg: f64,
    /// Longitude of ascending node $\Omega$ in degrees.
    pub om_deg: f64,
    /// Orbital inclination $i$ in degrees.
    pub i_deg: f64,
}

/// Returns the canonical [`HypotheticalElements`] for `id` (`MODEL_PASSPORTS.md` §3.1).
#[must_use]
pub const fn hypothetical_elements(id: HypotheticalId) -> HypotheticalElements {
    match id {
        HypotheticalId::Cupido => HypotheticalElements {
            id,
            epoch_jd: 2_415_020.0,
            equinox_jd: 2_415_020.0,
            m0_deg: 163.7409,
            a_au: 40.99837,
            e: 0.00460,
            w_deg: 171.4333,
            om_deg: 129.8325,
            i_deg: 1.0833,
        },
        HypotheticalId::Hades => HypotheticalElements {
            id,
            epoch_jd: 2_415_020.0,
            equinox_jd: 2_415_020.0,
            m0_deg: 27.6496,
            a_au: 50.66744,
            e: 0.00245,
            w_deg: 148.1796,
            om_deg: 161.3339,
            i_deg: 1.0500,
        },
        HypotheticalId::Zeus => HypotheticalElements {
            id,
            epoch_jd: 2_415_020.0,
            equinox_jd: 2_415_020.0,
            m0_deg: 165.1232,
            a_au: 59.21436,
            e: 0.00120,
            w_deg: 299.0440,
            om_deg: 0.0,
            i_deg: 0.0,
        },
        HypotheticalId::Kronos => HypotheticalElements {
            id,
            epoch_jd: 2_415_020.0,
            equinox_jd: 2_415_020.0,
            m0_deg: 169.0193,
            a_au: 64.81690,
            e: 0.00305,
            w_deg: 208.8801,
            om_deg: 0.0,
            i_deg: 0.0,
        },
        HypotheticalId::Apollon => HypotheticalElements {
            id,
            epoch_jd: 2_415_020.0,
            equinox_jd: 2_415_020.0,
            m0_deg: 138.0533,
            a_au: 70.29949,
            e: 0.0,
            w_deg: 0.0,
            om_deg: 0.0,
            i_deg: 0.0,
        },
        HypotheticalId::Admetos => HypotheticalElements {
            id,
            epoch_jd: 2_415_020.0,
            equinox_jd: 2_415_020.0,
            m0_deg: 351.3350,
            a_au: 73.62765,
            e: 0.0,
            w_deg: 0.0,
            om_deg: 0.0,
            i_deg: 0.0,
        },
        HypotheticalId::Vulkanus => HypotheticalElements {
            id,
            epoch_jd: 2_415_020.0,
            equinox_jd: 2_415_020.0,
            m0_deg: 55.8983,
            a_au: 77.25568,
            e: 0.0,
            w_deg: 0.0,
            om_deg: 0.0,
            i_deg: 0.0,
        },
        HypotheticalId::Poseidon => HypotheticalElements {
            id,
            epoch_jd: 2_415_020.0,
            equinox_jd: 2_415_020.0,
            m0_deg: 165.5163,
            a_au: 83.66907,
            e: 0.0,
            w_deg: 0.0,
            om_deg: 0.0,
            i_deg: 0.0,
        },
        HypotheticalId::Transpluto => HypotheticalElements {
            id,
            epoch_jd: 2_368_547.66,
            equinox_jd: 2_431_456.5,
            m0_deg: 0.0,
            a_au: 77.77500,
            e: 0.30000,
            w_deg: 0.7000,
            om_deg: 0.0,
            i_deg: 0.0,
        },
    }
}

/// Evaluates the heliocentric ICRF position vector (in AU) of `el` at `jd_tt` (`MODEL_PASSPORTS.md` §3.2).
#[must_use]
pub fn uranian_heliocentric_icrf_au(el: &HypotheticalElements, jd_tt: f64) -> Vec3 {
    let a = el.a_au;
    let e = el.e;
    let n = GAUSS_K_DEG_DAY / a.powf(1.5);
    let m_rad = (el.m0_deg + n * (jd_tt - el.epoch_jd))
        .rem_euclid(360.0)
        .to_radians();

    let mut ecc_anom = m_rad;
    for _ in 0..15 {
        let de = (ecc_anom - e * ecc_anom.sin() - m_rad) / (1.0 - e * ecc_anom.cos());
        ecc_anom -= de;
        if de.abs() < 1e-14 {
            break;
        }
    }

    let x_orb = a * (ecc_anom.cos() - e);
    let y_orb = a * (1.0 - e * e).sqrt() * ecc_anom.sin();

    let (sw, cw) = el.w_deg.to_radians().sin_cos();
    let (som, com) = el.om_deg.to_radians().sin_cos();
    let (si, ci) = el.i_deg.to_radians().sin_cos();

    let x_ecl = (com * cw - som * sw * ci) * x_orb + (-com * sw - som * cw * ci) * y_orb;
    let y_ecl = (som * cw + com * sw * ci) * x_orb + (-som * sw + com * cw * ci) * y_orb;
    let z_ecl = (sw * si) * x_orb + (cw * si) * y_orb;
    let v_ecl_eq = [x_ecl, y_ecl, z_ecl];

    // Rotate from mean ecliptic of equinox_jd to ICRF
    let eps0 = obl06(el.equinox_jd, 0.0);
    let (_, _, r_bp) = bp06(el.equinox_jd, 0.0);
    let r_ecl_eq = mat_mul(&rot_x(eps0), &r_bp);
    mat_vec(&mat_transpose(&r_ecl_eq), &v_ecl_eq)
}

/// Evaluates the 3D position vector `[x, y, z]` (in AU) of `el` in the output coordinate frame
/// specified by `options` at a single [`EpochContext`].
#[must_use]
pub fn eval_hypothetical_vector_au(
    el: &HypotheticalElements,
    ctx: &EpochContext,
    options: &CalcOptions,
) -> Vec3 {
    let jd_tt = ctx.tt.jd_f64();
    let (p_obs_km, v_obs_kms) = match options.origin {
        Origin::Geocentric => (ctx.earth_bary_km.pos_km, ctx.earth_bary_km.vel_km_s),
        Origin::Topocentric(_) => (ctx.observer_bary_km.pos_km, ctx.observer_bary_km.vel_km_s),
        Origin::Heliocentric => (ctx.sun_bary_km.pos_km, ctx.sun_bary_km.vel_km_s),
        Origin::Barycentric => ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
    };
    let r_obs_helio_au = vec_scale(&vec_sub(&p_obs_km, &ctx.sun_bary_km.pos_km), 1.0 / AU_KM);

    let vec_icrf_au = match options.aberration {
        AbMode::Geometric => {
            let r_obj_helio_au = uranian_heliocentric_icrf_au(el, jd_tt);
            vec_sub(&r_obj_helio_au, &r_obs_helio_au)
        }
        AbMode::Astrometric | AbMode::NoAberration | AbMode::NoDeflection | AbMode::Apparent => {
            let mut dt = 0.0_f64;
            let mut rho_au = [0.0_f64; 3];
            let mut dist_au = 0.0_f64;
            for _ in 0..3 {
                let r_obj_helio_au = uranian_heliocentric_icrf_au(el, jd_tt - dt);
                rho_au = vec_sub(&r_obj_helio_au, &r_obs_helio_au);
                dist_au = vec_norm(&rho_au);
                dt = dist_au / C_AU_DAY;
            }
            if matches!(options.aberration, AbMode::Apparent | AbMode::NoDeflection)
                && matches!(options.origin, Origin::Geocentric | Origin::Topocentric(_))
                && dist_au > 1e-12
            {
                let u_unit = vec_scale(&rho_au, 1.0 / dist_au);
                let v2_c2 = (vec_dot(&v_obs_kms, &v_obs_kms) / (C_KM_S * C_KM_S)).min(0.999_999);
                let bm1 = (1.0 - v2_c2).sqrt();
                let v_ab = vec_scale(&v_obs_kms, 1.0 / C_KM_S);
                let s_dist_au = vec_norm(&r_obs_helio_au);
                let u_ab = ab(&u_unit, &v_ab, s_dist_au, bm1);
                vec_scale(&u_ab, dist_au)
            } else {
                rho_au
            }
        }
    };

    let r_mat = ctx
        .frame
        .icrf_to_frame_matrix(options.plane, options.equinox, options.nutation);
    let vec_frame_au = mat_vec(&r_mat, &vec_icrf_au);
    if options.plane == Plane::Horizontal {
        let obs_opt = match options.origin {
            Origin::Topocentric(o) => Some(o),
            _ => ctx.observer,
        };
        if let Some(obs) = obs_opt {
            return equatorial_to_horizontal_vector_au(
                obs,
                ctx.gast_rad,
                &vec_frame_au,
                AtmosphereParams::default(),
            );
        }
    }
    vec_frame_au
}

/// Evaluates any of the 9 [`HypotheticalId`] bodies inside a [`Moment`].
pub fn calc_hypothetical_in_moment(
    moment: &Moment<'_>,
    id: HypotheticalId,
    options: &CalcOptions,
    sidereal_offset_and_rate: Option<(f64, f64)>,
) -> Result<CalcResult, EphemerisError> {
    let center = moment.center();
    if options.plane == Plane::Horizontal
        && !matches!(options.origin, Origin::Topocentric(_))
        && center.observer.is_none()
    {
        return Err(EphemerisError::InvalidInput {
            message: "Plane::Horizontal requires a topocentric Observer".to_owned(),
        });
    }
    let el = hypothetical_elements(id);
    let pos_0 = eval_hypothetical_vector_au(&el, center, options);

    let (vel_0, mut sph) = if options.include_speed {
        let stencil = moment.planet_stencil_contexts()?;
        let mut stencil_opts = *options;
        if options.plane == Plane::Ecliptic && sidereal_offset_and_rate.is_some() {
            stencil_opts.nutation = crate::types::NutationMode::Mean;
        }
        let pos_m2 = eval_hypothetical_vector_au(&el, &stencil[0], &stencil_opts);
        let pos_m1 = eval_hypothetical_vector_au(&el, &stencil[1], &stencil_opts);
        let pos_p1 = eval_hypothetical_vector_au(&el, &stencil[2], &stencil_opts);
        let pos_p2 = eval_hypothetical_vector_au(&el, &stencil[3], &stencil_opts);
        let v0 = stencil_4pt_derivative(&pos_m2, &pos_m1, &pos_p1, &pos_p2, H_PLANET_DAYS);
        let pos_for_rate = if stencil_opts.nutation != options.nutation {
            eval_hypothetical_vector_au(&el, center, &stencil_opts)
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
    if options.plane == Plane::Ecliptic {
        if let Some((aya_deg, aya_rate_deg_day)) = sidereal_offset_and_rate {
            sph.lon_deg = (sph.lon_deg - aya_deg).rem_euclid(360.0);
            if options.include_speed {
                sph.lon_speed_deg_day -= aya_rate_deg_day;
            }
            let (p_sid, v_sid) = crate::frame::spherical_state_to_cart(&sph, options.include_speed);
            cart_pos = p_sid;
            cart_vel = v_sid;
        }
    }

    let center_id = match options.origin {
        Origin::Geocentric | Origin::Topocentric(_) => 399,
        Origin::Heliocentric => 10,
        Origin::Barycentric => 0,
    };

    Ok(CalcResult {
        body: BodyId::Hypothetical(id),
        target_id: -(id as i32 + 1),
        center_id,
        target_kind: TargetKind::SyntheticPoint,
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
        magnitude_v: None,
        cartesian_au: cart_pos,
        cartesian_vel_au_day: cart_vel,
        warnings: Vec::new(),
    })
}
