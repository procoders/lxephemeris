//! `Moment` and `EpochContext`: per-epoch precomputed time scales, IAU 2006/2000A frame matrices,
//! sidereal time (`GAST`/`GMST`/`ERA`), Earth/Sun/Observer BCRS states, and lazy (`OnceLock`)
//! 4-point velocity stencil sub-contexts ($t_0 \pm h, t_0 \pm 2h$).

use crate::ayanamsa::{eval_ayanamsa_in_moment, AyanamsaValue};
use crate::error::EphemerisError;
use crate::frame::{
    mat_add_scaled, mat_transpose, mat_vec, nutation_rate_rad_day, precession_rate_per_day,
    vec_add, vec_scale, FrameTransform, Mat3,
};
use crate::spk::{SpkSource, SpkState};
use crate::time::{
    earth_rotation_angle_rad, gmst_rad, DeltaTInfo, Tdb, TimeInput, TimeTable, Tt, Ut1,
};
use crate::types::{
    AyanamsaModel, CalcOptions, NutationMode, Observer, Plane, Zodiac, SECONDS_PER_DAY,
};
use sofars::astro::pvtob;
use sofars::erst::gst06;
use std::sync::OnceLock;

/// Stencil step size $h$ in days for the Moon and lunar points ($10^{-4}\text{ d} = 8.64\text{ s}$).
pub const H_MOON_DAYS: f64 = 1.0e-4;

/// Stencil step size $h$ in days for the Sun, planets, 68 minor bodies, and Uranians ($10^{-3}\text{ d} = 86.4\text{ s}$).
pub const H_PLANET_DAYS: f64 = 1.0e-3;

/// Precomputed astrometric and frame context at a single instant.
#[derive(Debug, Clone)]
pub struct EpochContext {
    /// Universal Time (`UT1`).
    pub ut1: Ut1,
    /// Terrestrial Time (`TT`).
    pub tt: Tt,
    /// Barycentric Dynamical Time (`TDB`).
    pub tdb: Tdb,
    /// $\Delta T = \text{TT} - \text{UT1}$ metadata and extrapolation flag.
    pub delta_t_info: DeltaTInfo,
    /// Optional topocentric WGS84 observer.
    pub observer: Option<Observer>,
    /// Precomputed IAU 2006/2000A (or Vondrák 2011) rotation matrices and obliquity/nutation angles.
    pub frame: FrameTransform,
    /// Earth Rotation Angle (`ERA`, IAU 2000) in radians `[0, 2π)`.
    pub era_rad: f64,
    /// Greenwich Mean Sidereal Time (`GMST`, IAU 2006) in radians `[0, 2π)`.
    pub gmst_rad: f64,
    /// Greenwich Apparent Sidereal Time (`GAST`, IAU 2006/2000A) in radians `[0, 2π)`.
    pub gast_rad: f64,
    /// Barycentric (SSB=0) ICRF state of the Earth center (`399`) in km and km/s at `tdb`.
    pub earth_bary_km: SpkState,
    /// Barycentric (SSB=0) ICRF state of the Sun center (`10`) in km and km/s at `tdb`.
    pub sun_bary_km: SpkState,
    /// Barycentric (SSB=0) ICRF state of the Jupiter system barycenter (`5`) in km and km/s at `tdb`.
    pub jup_bary_km: SpkState,
    /// Barycentric (SSB=0) ICRF state of the Saturn system barycenter (`6`) in km and km/s at `tdb`.
    pub sat_bary_km: SpkState,
    /// Geocentric ICRF offset `(pos_km, vel_km_s)` of the topocentric observer (zero if geocentric).
    pub observer_geo_icrf_km: SpkState,
    /// Barycentric (SSB=0) ICRF state of the observer (`earth_bary_km + observer_geo_icrf_km`).
    pub observer_bary_km: SpkState,
}

impl EpochContext {
    /// Builds a complete [`EpochContext`] at `input` for `observer` using `spk` and `time_table`.
    pub fn new(
        spk: &dyn SpkSource,
        time_table: &TimeTable,
        input: TimeInput,
        observer: Option<Observer>,
    ) -> Result<Self, EphemerisError> {
        Self::new_inner(spk, time_table, input, observer, None, None)
    }

    /// Builds an [`EpochContext`] reusing pre-evaluated `(dpsi_rad, deps_rad)` and `(sun, jup, sat)` states when available.
    pub fn new_inner(
        spk: &dyn SpkSource,
        time_table: &TimeTable,
        input: TimeInput,
        observer: Option<Observer>,
        cached_nutation: Option<(f64, f64)>,
        cached_sun_giants: Option<(SpkState, SpkState, SpkState)>,
    ) -> Result<Self, EphemerisError> {
        let (ut1, tt, tdb, delta_t_info) = time_table.resolve_scales(input, observer);
        let frame = match cached_nutation {
            Some((dpsi, deps)) => FrameTransform::at_tt_with_nutation(tt, dpsi, deps),
            None => FrameTransform::at_tt(tt),
        };
        let era_rad = earth_rotation_angle_rad(ut1);
        let gmst_rad = gmst_rad(ut1, tt);
        let gast_rad = gst06(ut1.0.jd1, ut1.0.jd2, tt.0.jd1, tt.0.jd2, &frame.r_bpn);

        let tdb_split = tdb.0;
        let earth_bary_km = spk.state_icrf_km_split(399, 0, tdb_split.jd1, tdb_split.jd2)?;
        let (sun_bary_km, jup_bary_km, sat_bary_km) = match cached_sun_giants {
            Some(triple) => triple,
            None => (
                spk.state_icrf_km_split(10, 0, tdb_split.jd1, tdb_split.jd2)?,
                spk.state_icrf_km_split(5, 0, tdb_split.jd1, tdb_split.jd2)?,
                spk.state_icrf_km_split(6, 0, tdb_split.jd1, tdb_split.jd2)?,
            ),
        };

        let observer_geo_icrf_km = if let Some(obs) = observer {
            let lon_rad = obs.longitude_deg.to_radians();
            let lat_rad = obs.latitude_deg.to_radians();
            // pvtob with theta = gast_rad yields geocentric position [m] and velocity [m/s]
            // in the True Equator and Equinox of Date frame.
            let mut pv_m = [[0.0_f64; 3]; 2];
            pvtob(
                lon_rad,
                lat_rad,
                obs.altitude_m,
                0.0,
                0.0,
                0.0,
                gast_rad,
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
            // Rotate from True Equatorial of Date back to ICRF via R_BPN^T
            let r_bpn_t = mat_transpose(&frame.r_bpn);
            SpkState {
                pos_km: mat_vec(&r_bpn_t, &pos_true_eq_km),
                vel_km_s: mat_vec(&r_bpn_t, &vel_true_eq_km_s),
                target_id: 399,
                observer_id: 399,
            }
        } else {
            SpkState {
                pos_km: [0.0, 0.0, 0.0],
                vel_km_s: [0.0, 0.0, 0.0],
                target_id: 399,
                observer_id: 399,
            }
        };

        let observer_bary_km = SpkState {
            pos_km: vec_add(&earth_bary_km.pos_km, &observer_geo_icrf_km.pos_km),
            vel_km_s: vec_add(&earth_bary_km.vel_km_s, &observer_geo_icrf_km.vel_km_s),
            target_id: 399,
            observer_id: 0,
        };

        Ok(Self {
            ut1,
            tt,
            tdb,
            delta_t_info,
            observer,
            frame,
            era_rad,
            gmst_rad,
            gast_rad,
            earth_bary_km,
            sun_bary_km,
            jup_bary_km,
            sat_bary_km,
            observer_geo_icrf_km,
            observer_bary_km,
        })
    }
}

/// Thread-safe, per-epoch evaluation context (`engine.at(epoch, observer)`) that caches central
/// epoch $t_0$ quantities and lazily initializes shifted sub-contexts ($t_0 \pm h, t_0 \pm 2h$)
/// via `OnceLock` when velocities are requested.
pub struct Moment<'a> {
    spk: &'a dyn SpkSource,
    time_table: &'a TimeTable,
    center: EpochContext,
    nutation_rate: OnceLock<(f64, f64)>,
    precession_rate: OnceLock<(f64, Mat3)>,
    ayanamsa_cache: OnceLock<(AyanamsaModel, AyanamsaValue)>,
    moon_stencil: OnceLock<Result<[EpochContext; 4], EphemerisError>>,
    planet_stencil: OnceLock<Result<[EpochContext; 4], EphemerisError>>,
}

impl<'a> Moment<'a> {
    /// Creates a new [`Moment`] at `input` for an optional topocentric `observer`.
    pub fn new(
        spk: &'a dyn SpkSource,
        time_table: &'a TimeTable,
        input: TimeInput,
        observer: Option<Observer>,
    ) -> Result<Self, EphemerisError> {
        let center = EpochContext::new(spk, time_table, input, observer)?;
        Ok(Self {
            spk,
            time_table,
            center,
            nutation_rate: OnceLock::new(),
            precession_rate: OnceLock::new(),
            ayanamsa_cache: OnceLock::new(),
            moon_stencil: OnceLock::new(),
            planet_stencil: OnceLock::new(),
        })
    }

    /// Returns the cached nutation rate $(\dot{\Delta\psi}, \dot{\Delta\varepsilon})$ in radians/day at $t_0$.
    #[must_use]
    pub fn nutation_rate_rad_day(&self) -> (f64, f64) {
        *self
            .nutation_rate
            .get_or_init(|| nutation_rate_rad_day(self.center.tt.jd_f64()))
    }

    /// Returns the cached precession rate $(\dot{\varepsilon}_0, \dot{R}_{\text{bp}})$ per day at $t_0$.
    #[must_use]
    pub fn precession_rate_per_day(&self) -> &(f64, Mat3) {
        self.precession_rate.get_or_init(|| {
            precession_rate_per_day(
                self.center.tt,
                self.center.frame.eps0_rad,
                &self.center.frame.r_bp,
            )
        })
    }

    /// Returns the central [`EpochContext`] at $t_0$.
    #[must_use]
    pub const fn center(&self) -> &EpochContext {
        &self.center
    }

    /// Returns the underlying [`SpkSource`].
    #[must_use]
    pub fn spk(&self) -> &'a dyn SpkSource {
        self.spk
    }

    /// Returns the underlying [`TimeTable`].
    #[must_use]
    pub const fn time_table(&self) -> &'a TimeTable {
        self.time_table
    }

    /// Lazily initializes and returns the 4 shifted [`EpochContext`]s
    /// $[t_0 - 2h, t_0 - h, t_0 + h, t_0 + 2h]$ for the Moon stencil ($h = 10^{-4}\text{ d}$).
    pub fn moon_stencil_contexts(&self) -> Result<&[EpochContext; 4], EphemerisError> {
        self.moon_stencil
            .get_or_init(|| self.build_stencil(H_MOON_DAYS))
            .as_ref()
            .map_err(Clone::clone)
    }

    /// Lazily initializes and returns the 4 shifted [`EpochContext`]s
    /// $[t_0 - 2h, t_0 - h, t_0 + h, t_0 + 2h]$ for the planet/asteroid stencil ($h = 10^{-3}\text{ d}$).
    pub fn planet_stencil_contexts(&self) -> Result<&[EpochContext; 4], EphemerisError> {
        self.planet_stencil
            .get_or_init(|| self.build_stencil(H_PLANET_DAYS))
            .as_ref()
            .map_err(Clone::clone)
    }

    /// Builds a 4-point stencil $[t_0 - 2h, t_0 - h, t_0 + h, t_0 + 2h]$ for an arbitrary step `h_days`.
    pub fn build_stencil(&self, h_days: f64) -> Result<[EpochContext; 4], EphemerisError> {
        let c = &self.center;
        let obs = c.observer;
        let dpsi_0 = c.frame.dpsi_rad;
        let deps_0 = c.frame.deps_rad;
        let (ddpsi_rate, ddeps_rate) = self.nutation_rate_rad_day();
        let (deps0_rate, dr_bp_rate) = self.precession_rate_per_day();

        let shift_state = |st: &SpkState, dt_days: f64| -> SpkState {
            let dt_s = dt_days * SECONDS_PER_DAY;
            SpkState {
                pos_km: vec_add(&st.pos_km, &vec_scale(&st.vel_km_s, dt_s)),
                vel_km_s: st.vel_km_s,
                target_id: st.target_id,
                observer_id: st.observer_id,
            }
        };

        // Sidereal rotation rate in rad/day (ERA rate + general precession in RA)
        let omega_sid_rad_day = 1.002_737_811_911_354_6 * std::f64::consts::TAU
            + (4_612.156_534 / (3_600.0 * 36_525.0_f64)).to_radians();

        let build_shifted = |dt_days: f64| -> Result<EpochContext, EphemerisError> {
            let ut1 = Ut1(c.ut1.0.add_days(dt_days));
            let tt = Tt(c.tt.0.add_days(dt_days));
            let tdb = Tdb(c.tdb.0.add_days(dt_days));
            let dpsi = dpsi_0 + dt_days * ddpsi_rate;
            let deps = deps_0 + dt_days * ddeps_rate;
            let eps0 = c.frame.eps0_rad + dt_days * deps0_rate;
            let r_bp = mat_add_scaled(&c.frame.r_bp, dr_bp_rate, dt_days);
            let frame = FrameTransform::from_precomputed_bp(tt, eps0, r_bp, dpsi, deps);
            let era_rad = earth_rotation_angle_rad(ut1);
            let gmst_rad =
                (c.gmst_rad + dt_days * omega_sid_rad_day).rem_euclid(std::f64::consts::TAU);
            let gast_rad = (c.gast_rad
                + dt_days * omega_sid_rad_day
                + (dpsi - dpsi_0) * eps0.cos())
            .rem_euclid(std::f64::consts::TAU);

            let earth_bary_km =
                self.spk
                    .state_icrf_km_split(399, 0, tdb.0.jd1, tdb.0.jd2)?;
            let sun_bary_km = shift_state(&c.sun_bary_km, dt_days);
            let jup_bary_km = shift_state(&c.jup_bary_km, dt_days);
            let sat_bary_km = shift_state(&c.sat_bary_km, dt_days);

            let observer_geo_icrf_km = if let Some(o) = obs {
                let lon_rad = o.longitude_deg.to_radians();
                let lat_rad = o.latitude_deg.to_radians();
                let mut pv_m = [[0.0_f64; 3]; 2];
                pvtob(
                    lon_rad,
                    lat_rad,
                    o.altitude_m,
                    0.0,
                    0.0,
                    0.0,
                    gast_rad,
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
                let r_bpn_t = mat_transpose(&frame.r_bpn);
                SpkState {
                    pos_km: mat_vec(&r_bpn_t, &pos_true_eq_km),
                    vel_km_s: mat_vec(&r_bpn_t, &vel_true_eq_km_s),
                    target_id: 399,
                    observer_id: 399,
                }
            } else {
                SpkState {
                    pos_km: [0.0, 0.0, 0.0],
                    vel_km_s: [0.0, 0.0, 0.0],
                    target_id: 399,
                    observer_id: 399,
                }
            };

            let observer_bary_km = SpkState {
                pos_km: vec_add(&earth_bary_km.pos_km, &observer_geo_icrf_km.pos_km),
                vel_km_s: vec_add(&earth_bary_km.vel_km_s, &observer_geo_icrf_km.vel_km_s),
                target_id: 399,
                observer_id: 0,
            };

            Ok(EpochContext {
                ut1,
                tt,
                tdb,
                delta_t_info: c.delta_t_info.clone(),
                observer: obs,
                frame,
                era_rad,
                gmst_rad,
                gast_rad,
                earth_bary_km,
                sun_bary_km,
                jup_bary_km,
                sat_bary_km,
                observer_geo_icrf_km,
                observer_bary_km,
            })
        };

        let m2 = build_shifted(-2.0 * h_days)?;
        let m1 = build_shifted(-h_days)?;
        let p1 = build_shifted(h_days)?;
        let p2 = build_shifted(2.0 * h_days)?;
        Ok([m2, m1, p1, p2])
    }

    /// Evaluates the requested [`AyanamsaModel`] at this [`Moment`]'s central epoch, caching the
    /// first evaluated model on the [`Moment`].
    pub fn ayanamsa(
        &self,
        model: AyanamsaModel,
        nutation: NutationMode,
    ) -> Result<AyanamsaValue, EphemerisError> {
        if let Some((cached_model, cached_val)) = self.ayanamsa_cache.get() {
            if *cached_model == model {
                let mut out = cached_val.clone();
                out.effective_deg = match nutation {
                    NutationMode::True => out.true_deg,
                    NutationMode::Mean => out.mean_deg,
                };
                return Ok(out);
            }
        }
        let val = eval_ayanamsa_in_moment(self, model, nutation)?;
        let _ = self.ayanamsa_cache.set((model, val.clone()));
        Ok(val)
    }

    /// Returns `Some((effective_ayanamsa_deg, rate_deg_day))` if `options.plane == Plane::Ecliptic`
    /// and `options.zodiac == Zodiac::Sidereal(model)`, or `None` otherwise.
    ///
    /// Rejects `Equinox::J2000` combined with `Zodiac::Sidereal(_)` to prevent double-applying precession.
    pub fn sidereal_offset_and_rate(
        &self,
        options: &CalcOptions,
    ) -> Result<Option<(f64, f64)>, EphemerisError> {
        if options.plane == Plane::Ecliptic {
            if let Zodiac::Sidereal(model) = options.zodiac {
                if options.equinox == crate::types::Equinox::J2000 {
                    return Err(EphemerisError::InvalidInput {
                        message: "Zodiac::Sidereal cannot be combined with Equinox::J2000 (ayanamsa is defined relative to Equinox::OfDate and would double-count precession)".to_owned(),
                    });
                }
                let aya = self.ayanamsa(model, options.nutation)?;
                return Ok(Some((aya.effective_deg, aya.rate_deg_day)));
            }
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frame::vec_norm;
    use crate::spk::AniseSpkBackend;
    use crate::types::JD_J2000;

    #[test]
    fn test_moment_geocentric_and_topocentric_pvtob() {
        let workspace_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(|p| p.parent())
            .unwrap_or_else(|| std::path::Path::new("."));
        let de440s_path = workspace_dir.join("data/de440s.bsp");
        if !de440s_path.exists() {
            return;
        }
        let backend = AniseSpkBackend::from_bsp_paths(&[de440s_path]).expect("load de440s");
        let time_table = TimeTable::embedded();

        // Geocentric Moment
        let m_geo = Moment::new(
            &backend,
            &time_table,
            TimeInput::from_tt_jd(JD_J2000),
            None,
        )
        .expect("geo moment");
        assert_eq!(vec_norm(&m_geo.center().observer_geo_icrf_km.pos_km), 0.0);

        // Topocentric Moment at equator (lon=0, lat=0, alt=0)
        // WGS84 equatorial radius is 6378.137 km, rotational speed ~0.4651 km/s
        let obs = Observer::new(0.0, 0.0, 0.0).expect("equator observer");
        let m_topo = Moment::new(
            &backend,
            &time_table,
            TimeInput::from_tt_jd(JD_J2000),
            Some(obs),
        )
        .expect("topo moment");
        let r_obs_km = vec_norm(&m_topo.center().observer_geo_icrf_km.pos_km);
        let v_obs_km_s = vec_norm(&m_topo.center().observer_geo_icrf_km.vel_km_s);
        assert!(
            (r_obs_km - 6378.137).abs() < 0.01,
            "unexpected equatorial radius: {r_obs_km}"
        );
        assert!(
            (v_obs_km_s - 0.4651).abs() < 0.01,
            "unexpected equatorial velocity: {v_obs_km_s}"
        );

        // Lazy stencil initialization
        let stencil = m_topo.planet_stencil_contexts().expect("planet stencil");
        assert_eq!(stencil.len(), 4);
        assert!(
            (stencil[0].tt.jd_f64() - (JD_J2000 - 2.0 * H_PLANET_DAYS)).abs() < 1e-12
        );
    }
}
