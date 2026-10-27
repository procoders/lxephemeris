//! All 12 Vedic / Sidereal Ayanamsa models (`MODEL_PASSPORTS.md` §1):
//! - 9 epoch-anchored models propagated via Lieske et al. (1977, Table 5) / IAE (1989, p. 525)
//!   two-epoch general precession in longitude $p_A(t_0 \to t)$.
//! - 3 star-anchored models (`true_citra` = Spica HIP 65474 at $180^\circ$, `true_revati` = $\zeta$ Psc A
//!   HIP 5737 at $359^\circ 50'$, `true_pushya` = $\delta$ Cnc HIP 42911 at $106^\circ$) propagated via
//!   `sofars::astro::pmpx` with geocentric annual parallax and zero annual aberration to the mean
//!   ecliptic of date.

use crate::error::EphemerisError;
use crate::frame::{mat_vec, rot_x, vec_scale};
use crate::moment::Moment;
use crate::spk::SpkSource;
use crate::stars::{find_star_by_hip, DAS2R, JD_HIP_J1991_25};
use crate::time::{tt_to_tdb, Tt};
use crate::types::{AyanamsaModel, HipId, NutationMode, AU_KM, JD_J2000};
use serde::{Deserialize, Serialize};
use sofars::astro::pmpx;
use sofars::pnp::{bp06, ltpb, ltpecl, nut00a, nut80, obl06};

/// Complete evaluation of an [`AyanamsaModel`] at an epoch (`TT`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AyanamsaValue {
    pub model: AyanamsaModel,
    pub jd_tt: f64,
    /// Mean ayanamsa $A_{\text{mean}}(t)$ in degrees `[0, 360)` (without $\Delta\psi$).
    pub mean_deg: f64,
    /// True ayanamsa $A_{\text{true}}(t) = A_{\text{mean}}(t) + \Delta\psi(t)$ in degrees `[0, 360)`.
    pub true_deg: f64,
    /// Effective ayanamsa in degrees `[0, 360)` for the requested [`NutationMode`].
    pub effective_deg: f64,
    /// Secular rate $\dot{A}_{\text{mean}}(t)$ in degrees per day (`TT`).
    pub rate_deg_day: f64,
    /// Nutation in longitude $\Delta\psi(t)$ in arcseconds (IAU 2000A).
    pub dpsi_arcsec: f64,
}

/// Evaluates Lieske et al. (1977, Table 5) / Indian Astronomical Ephemeris (1989, p. 525)
/// two-epoch general precession in longitude $p_A(t_0 \to t)$ in arcseconds and secular rate in arcsec/day.
#[must_use]
pub fn iae_general_precession_arcsec(jd_t0: f64, jd_t: f64) -> (f64, f64) {
    let big_t = (jd_t0 - JD_J2000) / 36_525.0;
    let t = (jd_t - jd_t0) / 36_525.0;
    let c1 = 5029.0966 + 2.22226 * big_t - 0.000_042 * big_t * big_t;
    let c2 = 1.11161 - 0.000_127 * big_t;
    let c3 = -0.000_113;

    let p_arcsec = c1 * t + c2 * t * t + c3 * t * t * t;
    let rate_arcsec_day = (c1 + 2.0 * c2 * t + 3.0 * c3 * t * t) / 36_525.0;
    (p_arcsec, rate_arcsec_day)
}

/// Returns `(t0_jd_tt, a0_mean_deg)` for an epoch-defined [`AyanamsaModel`], or `None` if star-anchored.
#[must_use]
pub fn epoch_model_anchor(model: AyanamsaModel) -> Option<(f64, f64)> {
    match model {
        AyanamsaModel::FaganBradley => Some((
            2_433_282.423_459_05,
            24.0 + 2.0 / 60.0 + 31.36 / 3_600.0,
        )),
        AyanamsaModel::Lahiri => {
            // IAE 1989 p. 525: A_true(1956-03-21 00:00 TT) = 23° 15' 00.658" with IAU 1980 nutation subtracted
            let (dpsi_1956, _) = nut80(2_435_553.5, 0.0);
            let lahiri_mean_1956 = (23.0 * 3_600.0 + 15.0 * 60.0 + 0.658
                - dpsi_1956.to_degrees() * 3_600.0)
                / 3_600.0;
            Some((2_435_553.5, lahiri_mean_1956))
        }
        AyanamsaModel::DeLuce => Some((JD_J2000, 27.815_752_807)),
        AyanamsaModel::Raman => Some((2_415_020.0, 21.014_21)),
        AyanamsaModel::Ushashashi => Some((2_415_020.0, 18.660_96)),
        AyanamsaModel::Krishnamurti => Some((2_415_020.0, 22.363_659)),
        AyanamsaModel::DjwhalKhul => Some((2_494_459.5, 30.0)),
        AyanamsaModel::Yukteshwar => Some((2_415_020.0, 21.0 + 4.0 / 60.0 + 56.0 / 3_600.0)),
        AyanamsaModel::JnBhasin => Some((2_415_020.0, 21.0 + 21.0 / 60.0 + 56.0 / 3_600.0)),
        AyanamsaModel::BabylonKugler1 => Some((JD_J2000, 23.533_639_900)),
        AyanamsaModel::BabylonKugler2 => Some((JD_J2000, 24.933_639_900)),
        AyanamsaModel::BabylonKugler3 => Some((JD_J2000, 25.783_639_900)),
        AyanamsaModel::BabylonHuber => Some((JD_J2000, 24.733_639_900)),
        AyanamsaModel::BabylonEtpsc => Some((JD_J2000, 24.522_527_973)),
        AyanamsaModel::Hipparchos => Some((JD_J2000, 20.247_788_132)),
        AyanamsaModel::Sassanian => Some((JD_J2000, 19.992_959_400)),
        AyanamsaModel::GalCent0Sag => Some((JD_J2000, 26.846_045_806)),
        AyanamsaModel::J2000 => Some((JD_J2000, 0.0)),
        AyanamsaModel::J1900 => Some((2_415_020.0, 0.0)),
        AyanamsaModel::B1950 => Some((2_433_282.423_459_05, 0.0)),
        AyanamsaModel::Suryasiddhanta | AyanamsaModel::Aryabhata => {
            Some((1_903_396.877_069_8, 0.0))
        }
        AyanamsaModel::SuryasiddhantaMsun => Some((JD_J2000, 20.680_424_905)),
        AyanamsaModel::AryabhataMsun => Some((JD_J2000, 20.657_427_361)),
        AyanamsaModel::SsRevati => Some((JD_J2000, 20.103_388_395)),
        AyanamsaModel::SsCitra => Some((1_903_396.877_069_8, 2.0 + 6.0 / 60.0 + 38.536 / 3_600.0)),
        AyanamsaModel::GalCentRGilbrand => Some((JD_J2000, 22.469_104_789)),
        AyanamsaModel::GalEqIau1958 => Some((JD_J2000, 30.023_172_469)),
        AyanamsaModel::GalEqTrue => Some((JD_J2000, 30.076_092_180)),
        AyanamsaModel::GalEqMula => Some((JD_J2000, 23.409_425_513)),
        AyanamsaModel::GalAlignMardyks => Some((JD_J2000, 30.017_793_877)),
        AyanamsaModel::TrueMula => Some((JD_J2000, 24.579_972_884)),
        AyanamsaModel::GalCentMulaWilhelm => Some((JD_J2000, 20.039_233_112)),
        AyanamsaModel::Aryabhata522 => Some((1_911_797.5, 0.0)),
        AyanamsaModel::BabylonBritton => Some((JD_J2000, 24.615_752_807)),
        AyanamsaModel::TrueSheratan => Some((JD_J2000, 25.234_449_334)),
        AyanamsaModel::TrueAldebaran => Some((JD_J2000, 356.846_045_806)),
        AyanamsaModel::GalEquFiorenza => Some((JD_J2000, 25.000_019_094)),
        AyanamsaModel::ValensMoon => Some((JD_J2000, 22.795_609_328)),
        AyanamsaModel::Aldebaran15Tau => Some((JD_J2000, 24.758_923_900)),
        AyanamsaModel::User {
            t0_jd_tt_bits,
            ayan_t0_deg_bits,
        } => Some((
            f64::from_bits(t0_jd_tt_bits),
            f64::from_bits(ayan_t0_deg_bits),
        )),
        AyanamsaModel::TrueCitra
        | AyanamsaModel::TrueRevati
        | AyanamsaModel::TruePushya => None,
    }
}

/// Returns `(anchor_hip, target_sidereal_lon_deg)` for a star-anchored [`AyanamsaModel`], or `None` if epoch-defined.
#[must_use]
pub const fn star_model_anchor(model: AyanamsaModel) -> Option<(HipId, f64)> {
    match model {
        AyanamsaModel::TrueCitra => Some((HipId(65474), 180.0)),
        AyanamsaModel::TrueRevati => Some((HipId(5737), 359.0 + 50.0 / 60.0)),
        AyanamsaModel::TruePushya => Some((HipId(42911), 106.0)),
        _ => None,
    }
}

/// Evaluates the mean ecliptic longitude of date `[0, 360)` of `hip` at `jd_tt` using `pmpx`
/// with geocentric parallax and **zero annual aberration**.
fn star_mean_ecl_lon_no_ab(
    spk: &dyn SpkSource,
    hip: HipId,
    jd_tt: f64,
) -> Result<f64, EphemerisError> {
    let star = find_star_by_hip(hip).ok_or_else(|| EphemerisError::UnsupportedBody {
        body: format!("AnchorStar(HIP {})", hip.0),
    })?;
    let tt = Tt::from_jd(jd_tt);
    let tdb = tt_to_tdb(tt, None);
    let st_earth = spk.state_icrf_km_split(399, 0, tdb.0.jd1, tdb.0.jd2)?;
    let p_eb_au = vec_scale(&st_earth.pos_km, 1.0 / AU_KM);

    let ra1 = star.ra_1991_deg.to_radians();
    let de1 = star.de_1991_deg.to_radians();
    let cos_dec = de1.cos();
    let safe_cos_dec = if cos_dec.abs() < 1.0e-12 {
        1.0e-12_f64.copysign(cos_dec)
    } else {
        cos_dec
    };
    let pmr = (star.pm_ra_cosdec_mas_yr * 1e-3 * DAS2R) / safe_cos_dec;
    let pmd = star.pm_dec_mas_yr * 1e-3 * DAS2R;
    let px = (star.plx_mas * 1e-3).max(0.0);
    let pmt = (tdb.jd_f64() - JD_HIP_J1991_25) / 365.25;

    let u_icrs = pmpx(ra1, de1, pmr, pmd, px, star.rv_kms, pmt, p_eb_au);

    let r_ecl_mean = if crate::frame::is_within_iau2006_domain(jd_tt) {
        let eps0 = obl06(tt.0.jd1, tt.0.jd2);
        let (_, _, rbp) = bp06(tt.0.jd1, tt.0.jd2);
        crate::frame::mat_mul(&rot_x(eps0), &rbp)
    } else {
        let epj = 2000.0 + (jd_tt - JD_J2000) / 365.25;
        let rbp = ltpb(epj);
        let eq_pole = [rbp[2][0], rbp[2][1], rbp[2][2]];
        let ecl_pole = ltpecl(epj);
        let eps0 = crate::frame::vec_dot(&eq_pole, &ecl_pole)
            .clamp(-1.0, 1.0)
            .acos();
        crate::frame::mat_mul(&rot_x(eps0), &rbp)
    };

    let u_ecl = mat_vec(&r_ecl_mean, &u_icrs);
    Ok(u_ecl[1].atan2(u_ecl[0]).to_degrees().rem_euclid(360.0))
}

/// Evaluates the ecliptic-plane rotation correction `(offset_deg, rate_deg_day)` relative to `JD_J2000`
/// for galactic-equator node models (`GalEqIau1958`, `GalEqTrue`, `GalEqMula`, `GalCentMulaWilhelm`)
/// whose reference point is defined by the intersection of the Galactic plane with the moving ecliptic of date.
#[must_use]
fn model_secular_plane_correction_deg(model: AyanamsaModel, jd_tt: f64) -> (f64, f64) {
    let t_cy = (jd_tt - JD_J2000) / 36_525.0;
    match model {
        AyanamsaModel::GalEqIau1958 | AyanamsaModel::GalEqTrue | AyanamsaModel::GalEqMula => {
            let c1 = 0.007_323_5;
            let c2 = -0.000_014_56;
            (
                c1 * t_cy + c2 * t_cy * t_cy,
                (c1 + 2.0 * c2 * t_cy) / 36_525.0,
            )
        }
        AyanamsaModel::GalCentMulaWilhelm => {
            let c1 = 0.062_219_4;
            let c2 = 0.000_352_80;
            (
                c1 * t_cy + c2 * t_cy * t_cy,
                (c1 + 2.0 * c2 * t_cy) / 36_525.0,
            )
        }
        _ => (0.0, 0.0),
    }
}

/// Evaluates `model` at `jd_tt` with `nutation` mode, returning a full [`AyanamsaValue`].
pub fn eval_ayanamsa(
    spk: &dyn SpkSource,
    model: AyanamsaModel,
    jd_tt: f64,
    nutation: NutationMode,
) -> Result<AyanamsaValue, EphemerisError> {
    let tt = Tt::from_jd(jd_tt);
    let (dpsi_rad, _) = nut00a(tt.0.jd1, tt.0.jd2);
    let dpsi_deg = dpsi_rad.to_degrees();

    let (mean_deg, rate_deg_day) = if let Some((t0, a0_deg)) = epoch_model_anchor(model) {
        let (p_arcsec, rate_arcsec_day) = iae_general_precession_arcsec(t0, jd_tt);
        let (plane_off_deg, plane_rate_deg_day) = model_secular_plane_correction_deg(model, jd_tt);
        (
            (a0_deg + p_arcsec / 3_600.0 + plane_off_deg).rem_euclid(360.0),
            rate_arcsec_day / 3_600.0 + plane_rate_deg_day,
        )
    } else if let Some((hip, target_lon)) = star_model_anchor(model) {
        let lon_0 = star_mean_ecl_lon_no_ab(spk, hip, jd_tt)?;
        let lon_m = star_mean_ecl_lon_no_ab(spk, hip, jd_tt - 0.5)?;
        let lon_p = star_mean_ecl_lon_no_ab(spk, hip, jd_tt + 0.5)?;
        let mean = (lon_0 - target_lon).rem_euclid(360.0);
        let rate = (lon_p - lon_m + 180.0).rem_euclid(360.0) - 180.0;
        (mean, rate)
    } else {
        return Err(EphemerisError::UnsupportedModel {
            model: format!("{model:?}"),
        });
    };

    let true_deg = (mean_deg + dpsi_deg).rem_euclid(360.0);
    let effective_deg = match nutation {
        NutationMode::True => true_deg,
        NutationMode::Mean => mean_deg,
    };

    Ok(AyanamsaValue {
        model,
        jd_tt,
        mean_deg,
        true_deg,
        effective_deg,
        rate_deg_day,
        dpsi_arcsec: dpsi_deg * 3_600.0,
    })
}

/// Evaluates `model` inside a [`Moment`], reusing the `Moment`'s precomputed nutation $\Delta\psi$,
/// `FrameTransform`, and Earth barycentric state at zero `nut00a` or SPK overhead.
pub fn eval_ayanamsa_in_moment(
    moment: &Moment<'_>,
    model: AyanamsaModel,
    nutation: NutationMode,
) -> Result<AyanamsaValue, EphemerisError> {
    let center = moment.center();
    let jd_tt = center.tt.jd_f64();
    let dpsi_deg = center.frame.dpsi_rad.to_degrees();

    let (mean_deg, rate_deg_day) = if let Some((t0, a0_deg)) = epoch_model_anchor(model) {
        let (p_arcsec, rate_arcsec_day) = iae_general_precession_arcsec(t0, jd_tt);
        let (plane_off_deg, plane_rate_deg_day) = model_secular_plane_correction_deg(model, jd_tt);
        (
            (a0_deg + p_arcsec / 3_600.0 + plane_off_deg).rem_euclid(360.0),
            rate_arcsec_day / 3_600.0 + plane_rate_deg_day,
        )
    } else if let Some((hip, target_lon)) = star_model_anchor(model) {
        let star = find_star_by_hip(hip).ok_or_else(|| EphemerisError::UnsupportedBody {
            body: format!("AnchorStar(HIP {})", hip.0),
        })?;
        let lon_0 = crate::stars::star_mean_ecliptic_lon_no_ab_at_ctx(star, center);
        let lon_m = star_mean_ecl_lon_no_ab(moment.spk(), hip, jd_tt - 0.5)?;
        let lon_p = star_mean_ecl_lon_no_ab(moment.spk(), hip, jd_tt + 0.5)?;
        let mean = (lon_0 - target_lon).rem_euclid(360.0);
        let rate = (lon_p - lon_m + 180.0).rem_euclid(360.0) - 180.0;
        (mean, rate)
    } else {
        return Err(EphemerisError::UnsupportedModel {
            model: format!("{model:?}"),
        });
    };

    let true_deg = (mean_deg + dpsi_deg).rem_euclid(360.0);
    let effective_deg = match nutation {
        NutationMode::True => true_deg,
        NutationMode::Mean => mean_deg,
    };

    Ok(AyanamsaValue {
        model,
        jd_tt,
        mean_deg,
        true_deg,
        effective_deg,
        rate_deg_day,
        dpsi_arcsec: dpsi_deg * 3_600.0,
    })
}
