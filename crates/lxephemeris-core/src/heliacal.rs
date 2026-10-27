//! Heliacal risings and settings (`R6.3`) using Bradley E. Schaefer's astronomical visibility model
//! (*arcus visionis*, Schaefer 1987, 1993, 2000):
//! - Schaefer, B. E. (1987), "Heliacal Rise Phenomena", *JHA (Archaeoastronomy Suppl.)*, 18, S19–S33.
//! - Schaefer, B. E. (1993), "Astronomy and the Limits of Vision", *Vistas in Astronomy*, 36, 311–361.
//! - Schaefer, B. E. (2000), "New Methods and Techniques for Historical Astronomy and Archaeoastronomy",
//!   *Archaeoastronomy*, 15, 121–136.

use crate::astrometry::{azalt, calc_body_in_moment, AtmosphereParams};
use crate::error::EphemerisError;
use crate::frame::{vec_dot, vec_norm, vec_sub};
use crate::moment::Moment;
use crate::spk::SpkSource;
use crate::stars::find_star_by_hip;
use crate::time::{TimeInput, TimeTable};
use crate::types::{
    AbMode, BodyId, CalcOptions, Equinox, NutationMode, Observer, Origin, Plane, PlanetId, Zodiac,
    AU_KM,
};
use serde::{Deserialize, Serialize};

/// Atmospheric and physiological parameters for Bradley Schaefer's heliacal visibility model (`R6.3`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HeliacalOptions {
    /// Visual atmospheric extinction coefficient $k_V$ in magnitudes per airmass at sea level (`0.20` default).
    pub k_v: f64,
    /// Observer's age in years (`36.0` default; controls dark/twilight pupil diameter).
    pub observer_age_years: f64,
    /// Observer's Snellen visual acuity ratio (`1.0` default for standard 20/20 vision).
    pub snellen_ratio: f64,
    /// Relative humidity in `[0.0, 1.0]` (`0.5` default).
    pub relative_humidity: f64,
    /// Ambient temperature in degrees Celsius (`15.0` default).
    pub temperature_c: f64,
    /// Atmospheric pressure in hPa (`1013.25` default).
    pub pressure_hpa: f64,
    /// Search window in days from `start_jd_ut1` (`365.25` days default).
    pub max_search_days: f64,
}

impl Default for HeliacalOptions {
    fn default() -> Self {
        Self {
            k_v: 0.20,
            observer_age_years: 36.0,
            snellen_ratio: 1.0,
            relative_humidity: 0.5,
            temperature_c: 15.0,
            pressure_hpa: 1013.25,
            max_search_days: 365.25,
        }
    }
}

/// Result of a heliacal visibility cycle search over `[start_jd_ut1, start_jd_ut1 + max_search_days]` (`R6.3`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HeliacalResult {
    /// Target celestial body (`FixedStar`, `Planet`, or `MinorPlanet`).
    pub body: BodyId,
    /// Human-readable canonical name of the target body.
    pub target_name: String,
    /// Visual magnitude $m_V$ of the body at heliacal rising (or at `start_jd_ut1` if no rising occurs).
    pub visual_magnitude: f64,
    /// Critical arc of vision (*arcus visionis* $AV_{\text{crit}} = h_{\text{obj}} - h_\odot$) in degrees.
    pub arcus_visionis_deg: f64,
    /// Morning First (Heliacal Rising, `JD_UT1`): first morning appearance in the East before sunrise.
    pub heliacal_rise_jd_ut1: Option<f64>,
    /// Morning Last (Cosmical Setting, `JD_UT1`): last morning setting in the West before sunrise.
    pub morning_set_jd_ut1: Option<f64>,
    /// Evening First (Acronychal Rising, `JD_UT1`): first evening rising in the East after sunset.
    pub evening_rise_jd_ut1: Option<f64>,
    /// Evening Last (Heliacal Setting, `JD_UT1`): last evening disappearance in the West after sunset.
    pub heliacal_set_jd_ut1: Option<f64>,
}

/// Evaluates the apparent visual magnitude $m_V$ of `body` at `moment` using IAU / Mallama & Hilton (2018)
/// planetary phase curves or the Hipparcos visual magnitude for fixed stars.
pub fn estimate_visual_magnitude(moment: &Moment<'_>, body: BodyId) -> Result<f64, EphemerisError> {
    match body {
        BodyId::FixedStar(hip) => {
            let star = find_star_by_hip(hip).ok_or_else(|| EphemerisError::UnsupportedBody {
                body: format!("FixedStar(HIP {})", hip.0),
            })?;
            Ok(star.vmag)
        }
        BodyId::Planet(planet) => {
            if planet == PlanetId::Sun {
                return Ok(-26.74);
            }
            let ctx = moment.center();
            let tdb = ctx.tdb.0;
            let st_body =
                moment
                    .spk()
                    .state_icrf_km_split(planet.naif_target_id(), 0, tdb.jd1, tdb.jd2)?;
            let r_helio_km = vec_sub(&st_body.pos_km, &ctx.sun_bary_km.pos_km);
            let r_geo_km = vec_sub(&st_body.pos_km, &ctx.earth_bary_km.pos_km);
            let r_au = (vec_norm(&r_helio_km) / AU_KM).max(1e-6);
            let delta_au = (vec_norm(&r_geo_km) / AU_KM).max(1e-6);
            let cos_alpha = (vec_dot(&r_helio_km, &r_geo_km)
                / (vec_norm(&r_helio_km) * vec_norm(&r_geo_km)).max(1e-12))
            .clamp(-1.0, 1.0);
            let alpha_deg = cos_alpha.acos().to_degrees();
            let dist_term = 5.0 * (r_au * delta_au).log10();

            let vmag = match planet {
                PlanetId::Mercury => {
                    -0.60
                        + dist_term
                        + 0.0498 * alpha_deg
                        - 0.000_488 * alpha_deg * alpha_deg
                        + 0.000_003_02 * alpha_deg * alpha_deg * alpha_deg
                }
                PlanetId::Venus => {
                    -4.40
                        + dist_term
                        + 0.0009 * alpha_deg
                        + 0.000_239 * alpha_deg * alpha_deg
                        - 0.000_000_65 * alpha_deg * alpha_deg * alpha_deg
                }
                PlanetId::Mars => -1.52 + dist_term + 0.016 * alpha_deg,
                PlanetId::Jupiter => -9.40 + dist_term + 0.005 * alpha_deg,
                PlanetId::Saturn => -8.88 + dist_term + 0.044 * alpha_deg,
                PlanetId::Uranus => -7.19 + dist_term + 0.0028 * alpha_deg,
                PlanetId::Neptune => -6.87 + dist_term,
                PlanetId::Pluto => -1.00 + dist_term,
                PlanetId::Moon => {
                    0.23 + dist_term + 0.026 * alpha_deg + 4.0e-9 * alpha_deg.powi(4)
                }
                _ => 0.0,
            };
            Ok(vmag)
        }
        BodyId::MinorPlanet(mpc) => {
            let ctx = moment.center();
            let tdb = ctx.tdb.0;
            let naif_id = mpc.standard_naif_spk_id();
            let st_body = moment.spk().state_icrf_km_split(naif_id, 0, tdb.jd1, tdb.jd2)?;
            let r_helio_km = vec_sub(&st_body.pos_km, &ctx.sun_bary_km.pos_km);
            let r_geo_km = vec_sub(&st_body.pos_km, &ctx.earth_bary_km.pos_km);
            let r_au = (vec_norm(&r_helio_km) / AU_KM).max(1e-6);
            let delta_au = (vec_norm(&r_geo_km) / AU_KM).max(1e-6);
            let cos_alpha = (vec_dot(&r_helio_km, &r_geo_km)
                / (vec_norm(&r_helio_km) * vec_norm(&r_geo_km)).max(1e-12))
            .clamp(-1.0, 1.0);
            let alpha_rad = cos_alpha.acos();
            let tan_half = (0.5 * alpha_rad).tan().max(0.0);
            // Bowell et al. (1989) IAU (H, G) phase functions Phi_1 and Phi_2
            let phi1 = (-3.332 * tan_half.powf(0.631)).exp();
            let phi2 = (-1.862 * tan_half.powf(1.218)).exp();
            let (h_abs, g_slope) = match mpc.0 {
                1 => (3.34, 0.12),      // Ceres
                2 => (4.13, 0.11),      // Pallas
                3 => (5.33, 0.32),      // Juno
                4 => (3.20, 0.32),      // Vesta
                10 => (5.43, 0.15),     // Hygiea
                16 => (5.90, 0.20),     // Psyche
                433 => (10.40, 0.46),   // Eros
                2060 => (6.50, 0.15),   // Chiron
                5145 => (7.00, 0.15),   // Pholus
                7066 => (9.50, 0.15),   // Nessus
                90377 => (1.50, 0.15),  // Sedna
                136108 => (0.20, 0.15), // Haumea
                136199 => (-1.19, 0.15),// Eris
                136472 => (-0.20, 0.15),// Makemake
                _ => (6.50, 0.15),
            };
            let phase_flux = ((1.0 - g_slope) * phi1 + g_slope * phi2).max(1e-6);
            Ok(h_abs + 5.0 * (r_au * delta_au).log10() - 2.5 * phase_flux.log10())
        }
        other => Err(EphemerisError::UnsupportedBody {
            body: format!("Heliacal visibility not defined for {other:?}"),
        }),
    }
}

/// Computes the critical *arcus visionis* $AV_{\text{crit}}$ (in degrees of solar depression at the instant
/// the object crosses the geometric horizon $h_{\text{obj}} = 0^\circ$) according to Bradley Schaefer
/// (1987, 1993, 2000).
#[must_use]
pub fn schaefer_critical_arcus_visionis_deg(
    visual_magnitude: f64,
    azimuth_diff_deg: f64,
    observer_altitude_m: f64,
    options: &HeliacalOptions,
) -> f64 {
    // 1. Decompose sea-level extinction k_V into Rayleigh + Ozone (~0.14 mag/airmass, scaling with P/T)
    //    and Aerosol (~k_V - 0.14, scaling with elevation and Schaefer 1993 relative-humidity growth factor)
    let p_factor = (options.pressure_hpa / 1013.25).clamp(0.1, 1.5);
    let t_factor = (288.15 / (273.15 + options.temperature_c).max(180.0)).clamp(0.5, 1.5);
    let rh = options.relative_humidity.clamp(0.05, 0.95);
    // Normalized humidity aerosol growth factor (1.0 at standard RH = 0.50, per Schaefer 1993 eq. 18)
    let rh_factor = ((1.0 - 0.50) / (1.0 - rh)).powf(0.33).clamp(0.6, 2.2);
    let k_rayleigh_oz = 0.14 * p_factor * t_factor;
    let k_aerosol_base = (options.k_v - 0.14).max(0.02);
    let k_aerosol = k_aerosol_base * (-observer_altitude_m / 1_500.0).exp() * rh_factor;
    let k_eff = (k_rayleigh_oz + k_aerosol).clamp(0.05, 1.5);

    // 2. Observer pupil age factor (Schaefer 1993 eq. 6): reference age 36 yr, Snellen ratio S
    let age = options.observer_age_years.clamp(10.0, 90.0);
    let snellen = options.snellen_ratio.clamp(0.1, 4.0);
    let dm_age = 2.5 * (((age * age - 36.0 * 36.0) / 10_000.0) * std::f64::consts::LOG10_E);
    let dm_snellen = -5.0 * snellen.log10();

    // 3. Twilight azimuth glare term (Kastner 1976 / Schaefer 1987): brighter near solar azimuth (delta_az -> 0)
    let dm_az = 0.40 * azimuth_diff_deg.to_radians().cos();

    // 4. Effective magnitude normalized to standard conditions (k_V = 0.20, age = 36, S = 1.0, RH = 0.5)
    let m_eff = visual_magnitude + 14.5 * (k_eff - 0.20) + dm_age + dm_snellen + dm_az;

    // 5. Schaefer (1987, Fig. 2 & Table 1; 1993 Hecht threshold parameterization) critical solar depression at sighting
    (10.15 + 1.50 * m_eff + 0.045 * m_eff * m_eff).clamp(3.5, 22.0)
}

/// Computes the critical twilight sighting altitude $h_{\text{sight}}$ (in degrees of true altitude above the horizon)
/// at which a celestial body of `visual_magnitude` overcomes Rozenberg (1966) horizon airmass extinction $k_{\text{eff}} X(h)$.
#[must_use]
fn schaefer_critical_sighting_altitude_deg(
    visual_magnitude: f64,
    observer_altitude_m: f64,
    is_planet: bool,
    options: &HeliacalOptions,
) -> f64 {
    let p_factor = (options.pressure_hpa / 1013.25).clamp(0.1, 1.5);
    let t_factor = (288.15 / (273.15 + options.temperature_c).max(180.0)).clamp(0.5, 1.5);
    let rh = options.relative_humidity.clamp(0.05, 0.95);
    let rh_factor = ((1.0 - 0.50) / (1.0 - rh)).powf(0.33).clamp(0.6, 2.2);
    let k_rayleigh_oz = 0.14 * p_factor * t_factor;
    let k_aerosol_base = (options.k_v - 0.14).max(0.02);
    let k_aerosol = k_aerosol_base * (-observer_altitude_m / 1_500.0).exp() * rh_factor;
    let k_eff = (k_rayleigh_oz + k_aerosol).clamp(0.05, 1.5);

    let age = options.observer_age_years.clamp(10.0, 90.0);
    let snellen = options.snellen_ratio.clamp(0.1, 4.0);
    let dm_age = 2.5 * (((age * age - 36.0 * 36.0) / 10_000.0) * std::f64::consts::LOG10_E);
    let dm_snellen = -5.0 * snellen.log10();

    if is_planet {
        let m_eff = visual_magnitude + 8.0 * (k_eff - 0.20) + dm_age + dm_snellen;
        (2.65 + 0.55 * m_eff).clamp(1.8, 5.0)
    } else {
        let m_eff = visual_magnitude + 14.5 * (k_eff - 0.20) + dm_age + dm_snellen;
        (4.67 + 1.54 * m_eff + 0.04 * m_eff * m_eff).clamp(1.0, 15.0)
    }
}

#[derive(Debug, Clone, Copy)]
struct HorizonCrossingSample {
    t_cross_ut1: f64,
    margin_deg: f64,
    sun_ha_deg: f64,
    av_crit_deg: f64,
    vmag: f64,
}

fn format_body_target_name(body: BodyId) -> String {
    match body {
        BodyId::Planet(p) => p.name().to_owned(),
        BodyId::FixedStar(hip) => find_star_by_hip(hip).map_or_else(
            || format!("HIP {}", hip.0),
            |s| s.name.to_owned(),
        ),
        BodyId::MinorPlanet(mpc) => crate::types::MINOR_BODIES_68
            .iter()
            .find(|m| m.mpc == mpc)
            .map_or_else(|| format!("MPC {}", mpc.0), |m| m.name.to_owned()),
        other => format!("{other:?}"),
    }
}

/// Evaluates the daily horizon crossing (`is_rise == true` for Rising, `false` for Setting) near `jd_ut1_day`.
fn eval_daily_horizon_crossing(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    body: BodyId,
    observer: Observer,
    jd_ut1_day: f64,
    is_rise: bool,
    options: &HeliacalOptions,
) -> Result<Option<HorizonCrossingSample>, EphemerisError> {
    let eq_opts = CalcOptions {
        origin: Origin::Topocentric(observer),
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };

    // 1. Evaluate object RA/Dec at jd_ut1_day and solve diurnal sighting hour angle
    //    sin(h_sight) = sin(phi)*sin(dec) + cos(phi)*cos(dec)*cos(H)
    let m0 = Moment::new(
        spk,
        time_table,
        TimeInput::from_ut1_jd(jd_ut1_day),
        Some(observer),
    )?;
    let is_planet = matches!(body, BodyId::Planet(_));
    let obj0 = calc_body_in_moment(&m0, body, &eq_opts)?;
    let vmag0 = estimate_visual_magnitude(&m0, body)?;
    let h_sight_nom_deg =
        schaefer_critical_sighting_altitude_deg(vmag0, observer.altitude_m, is_planet, options);
    let phi_rad = observer.latitude_deg.to_radians();
    let dec0_rad = obj0.lat_deg.to_radians();
    let cos_h0_horizon = -phi_rad.tan() * dec0_rad.tan();
    if cos_h0_horizon.abs() >= 1.0 {
        return Ok(None);
    }
    let max_alt0_deg = (phi_rad.sin() * dec0_rad.sin() + phi_rad.cos() * dec0_rad.cos())
        .clamp(-1.0, 1.0)
        .asin()
        .to_degrees();
    let h_sight0_rad = h_sight_nom_deg.min(0.5 * max_alt0_deg.max(0.0)).to_radians();
    let denom0 = (phi_rad.cos() * dec0_rad.cos()).max(1.0e-12);
    let cos_h0 = ((h_sight0_rad.sin() - phi_rad.sin() * dec0_rad.sin()) / denom0).clamp(-0.999_999, 0.999_999);
    let h0_deg = cos_h0.acos().to_degrees();
    let target_ha_deg = if is_rise { -h0_deg } else { h0_deg };

    let last0_deg = (m0.center().gast_rad.to_degrees() + observer.longitude_deg).rem_euclid(360.0);
    let curr_ha0_deg = (last0_deg - obj0.lon_deg + 180.0).rem_euclid(360.0) - 180.0;
    let dha0_deg = (target_ha_deg - curr_ha0_deg + 180.0).rem_euclid(360.0) - 180.0;
    let t1 = jd_ut1_day + dha0_deg / 360.985_647_366_29;

    // 2. Second Newton refinement at t1 (handles planetary motion in RA/Dec)
    let m1 = Moment::new(spk, time_table, TimeInput::from_ut1_jd(t1), Some(observer))?;
    let obj1 = calc_body_in_moment(&m1, body, &eq_opts)?;
    let dec1_rad = obj1.lat_deg.to_radians();
    let cos_h1_horizon = -phi_rad.tan() * dec1_rad.tan();
    if cos_h1_horizon.abs() >= 1.0 {
        return Ok(None);
    }
    let max_alt1_deg = (phi_rad.sin() * dec1_rad.sin() + phi_rad.cos() * dec1_rad.cos())
        .clamp(-1.0, 1.0)
        .asin()
        .to_degrees();
    let h_sight1_rad = h_sight_nom_deg.min(0.5 * max_alt1_deg.max(0.0)).to_radians();
    let denom1 = (phi_rad.cos() * dec1_rad.cos()).max(1.0e-12);
    let cos_h1 = ((h_sight1_rad.sin() - phi_rad.sin() * dec1_rad.sin()) / denom1).clamp(-0.999_999, 0.999_999);
    let h1_deg = cos_h1.acos().to_degrees();
    let target_ha1_deg = if is_rise { -h1_deg } else { h1_deg };
    let last1_deg = (m1.center().gast_rad.to_degrees() + observer.longitude_deg).rem_euclid(360.0);
    let curr_ha1_deg = (last1_deg - obj1.lon_deg + 180.0).rem_euclid(360.0) - 180.0;
    let dha1_deg = (target_ha1_deg - curr_ha1_deg + 180.0).rem_euclid(360.0) - 180.0;
    let t_cross = t1 + dha1_deg / 360.985_647_366_29;

    // 3. Evaluate Sun and object at t_cross
    let m_cross = Moment::new(
        spk,
        time_table,
        TimeInput::from_ut1_jd(t_cross),
        Some(observer),
    )?;
    let obj_cross = calc_body_in_moment(&m_cross, body, &eq_opts)?;
    let sun_cross = calc_body_in_moment(&m_cross, BodyId::Planet(PlanetId::Sun), &eq_opts)?;

    let hz_obj = azalt(
        observer,
        m_cross.center().gast_rad,
        obj_cross.lon_deg,
        obj_cross.lat_deg,
        AtmosphereParams::airless(),
    );
    let hz_sun = azalt(
        observer,
        m_cross.center().gast_rad,
        sun_cross.lon_deg,
        sun_cross.lat_deg,
        AtmosphereParams::airless(),
    );

    let is_planet = matches!(body, BodyId::Planet(_));
    let vmag = estimate_visual_magnitude(&m_cross, body)?;
    let az_diff = (hz_obj.azimuth_deg - hz_sun.azimuth_deg + 180.0).rem_euclid(360.0) - 180.0;
    let av_crit = if is_planet {
        let p_factor = (options.pressure_hpa / 1013.25).clamp(0.1, 1.5);
        let t_factor = (288.15 / (273.15 + options.temperature_c).max(180.0)).clamp(0.5, 1.5);
        let rh = options.relative_humidity.clamp(0.05, 0.95);
        let rh_factor = ((1.0 - 0.50) / (1.0 - rh)).powf(0.33).clamp(0.6, 2.2);
        let k_rayleigh_oz = 0.14 * p_factor * t_factor;
        let k_aerosol_base = (options.k_v - 0.14).max(0.02);
        let k_aerosol = k_aerosol_base * (-observer.altitude_m / 1_500.0).exp() * rh_factor;
        let k_eff = (k_rayleigh_oz + k_aerosol).clamp(0.05, 1.5);
        let age = options.observer_age_years.clamp(10.0, 90.0);
        let snellen = options.snellen_ratio.clamp(0.1, 4.0);
        let dm_age = 2.5 * (((age * age - 36.0 * 36.0) / 10_000.0) * std::f64::consts::LOG10_E);
        let dm_snellen = -5.0 * snellen.log10();
        if body == BodyId::Planet(PlanetId::Moon) {
            (11.0 + 8.0 * (k_eff - 0.20) + dm_age + dm_snellen).clamp(8.5, 16.0)
        } else {
            let m_eff = vmag + 8.0 * (k_eff - 0.20) + dm_age + dm_snellen;
            (9.35 + 0.85 * m_eff).clamp(7.0, 15.0)
        }
    } else {
        schaefer_critical_arcus_visionis_deg(vmag, az_diff, observer.altitude_m, options)
    };
    let solar_depression_deg = -hz_sun.altitude_true_deg;
    let margin_deg = solar_depression_deg - av_crit;

    Ok(Some(HorizonCrossingSample {
        t_cross_ut1: t_cross,
        margin_deg,
        sun_ha_deg: hz_sun.hour_angle_deg,
        av_crit_deg: av_crit,
        vmag,
    }))
}

/// Computes the four annual heliacal phenomena (`heliacal_rise_jd_ut1`, `morning_set_jd_ut1`,
/// `evening_rise_jd_ut1`, `heliacal_set_jd_ut1`) for `body` (`FixedStar`, `Planet`, or `MinorPlanet`)
/// starting at `start_jd_ut1` for `observer` (`|latitude_deg| <= 60°`, `R6.3`).
pub fn heliacal_phenomena(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    body: BodyId,
    start_jd_ut1: f64,
    observer: Observer,
    options: HeliacalOptions,
) -> Result<HeliacalResult, EphemerisError> {
    if !start_jd_ut1.is_finite()
        || !options.max_search_days.is_finite()
        || !(1.0..=36_525.0).contains(&options.max_search_days)
        || !options.pressure_hpa.is_finite()
        || !options.temperature_c.is_finite()
        || !options.relative_humidity.is_finite()
        || !options.k_v.is_finite()
        || !options.observer_age_years.is_finite()
        || !options.snellen_ratio.is_finite()
        || options.snellen_ratio <= 0.0
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid heliacal options or start_jd_ut1={start_jd_ut1}, max_search_days={}",
                options.max_search_days
            ),
        });
    }
    if observer.latitude_deg.abs() > 60.0 {
        return Err(EphemerisError::OutOfModelValidity {
            model: "SchaeferHeliacalVisibility_MaxLat60".to_owned(),
            epoch_jd_tt: observer.latitude_deg,
            min_jd_tt: -60.0,
            max_jd_tt: 60.0,
        });
    }
    if body == BodyId::Planet(PlanetId::Sun) || body == BodyId::Planet(PlanetId::Earth) {
        return Err(EphemerisError::UnsupportedBody {
            body: format!("Heliacal phenomena not applicable to {body:?}"),
        });
    }

    let max_end_jd_ut1 = start_jd_ut1 + options.max_search_days;
    let is_moon = body == BodyId::Planet(PlanetId::Moon);
    let is_inner_planet = matches!(body, BodyId::Planet(PlanetId::Mercury | PlanetId::Venus));
    let is_inner_or_moon = is_inner_planet || is_moon;
    let step_days = if is_moon {
        0.5_f64.min(options.max_search_days)
    } else if is_inner_planet {
        1.0_f64.min(options.max_search_days)
    } else {
        2.0_f64.min(options.max_search_days)
    };
    let total_steps = ((options.max_search_days / step_days).ceil() as usize).max(1);

    let mut heliacal_rise_jd_ut1 = None;
    let mut morning_set_jd_ut1 = None;
    let mut evening_rise_jd_ut1 = None;
    let mut heliacal_set_jd_ut1 = None;

    let m_init = Moment::new(
        spk,
        time_table,
        TimeInput::from_ut1_jd(start_jd_ut1),
        Some(observer),
    )?;
    let mut ref_vmag = estimate_visual_magnitude(&m_init, body)?;
    let mut ref_av =
        schaefer_critical_arcus_visionis_deg(ref_vmag, 15.0, observer.altitude_m, &options);

    let refine_crossing = |d_a: f64, d_b: f64, is_rise: bool, appearing: bool| -> Result<Option<(f64, f64, f64)>, EphemerisError> {
        let mut samples: Vec<HorizonCrossingSample> = Vec::with_capacity(6);
        let span = d_b - d_a;
        for frac in [0.0, 0.25, 0.5, 0.75, 1.0, 1.25] {
            let d_probe = (d_a + frac * span).min(max_end_jd_ut1);
            if let Some(s) = eval_daily_horizon_crossing(
                spk,
                time_table,
                body,
                observer,
                d_probe,
                is_rise,
                &options,
            )? {
                if s.t_cross_ut1 >= start_jd_ut1
                    && s.t_cross_ut1 <= max_end_jd_ut1
                    && !samples
                        .iter()
                        .any(|prev| (prev.t_cross_ut1 - s.t_cross_ut1).abs() < 0.5)
                {
                    samples.push(s);
                }
            }
        }
        samples.sort_by(|x, y| x.t_cross_ut1.total_cmp(&y.t_cross_ut1));
        let chosen = if appearing {
            samples
                .iter()
                .copied()
                .find(|s| s.margin_deg >= 0.0)
                .or_else(|| {
                    if d_a >= start_jd_ut1 {
                        samples.last().copied()
                    } else {
                        None
                    }
                })
        } else {
            samples
                .iter()
                .copied()
                .rfind(|s| s.margin_deg >= 0.0)
                .or_else(|| {
                    if d_a >= start_jd_ut1 {
                        samples.first().copied()
                    } else {
                        None
                    }
                })
        };
        Ok(chosen.map(|s| (s.t_cross_ut1, s.av_crit_deg, s.vmag)))
    };

    let d_seed = start_jd_ut1 - 1.0;
    let mut prev_rise: Option<(f64, HorizonCrossingSample)> =
        eval_daily_horizon_crossing(spk, time_table, body, observer, d_seed, true, &options)?
            .filter(|s| s.t_cross_ut1 < start_jd_ut1)
            .map(|s| (d_seed, s));
    let mut prev_set: Option<(f64, HorizonCrossingSample)> =
        eval_daily_horizon_crossing(spk, time_table, body, observer, d_seed, false, &options)?
            .filter(|s| s.t_cross_ut1 < start_jd_ut1)
            .map(|s| (d_seed, s));

    for i in 0..=total_steps {
        let d_curr = (start_jd_ut1 + (i as f64) * step_days).min(max_end_jd_ut1);

        // Check Rise events in the East:
        // - Outer bodies & stars: Morning First (Heliacal Rise) and Evening First (Acronychal Rise)
        // - Inner planets (Mercury, Venus) & Moon: Morning First (Heliacal Rise) and Morning Last (waning crescent / disappearance)
        let need_rise = if is_moon {
            morning_set_jd_ut1.is_none()
        } else if is_inner_planet {
            heliacal_rise_jd_ut1.is_none() || morning_set_jd_ut1.is_none()
        } else {
            heliacal_rise_jd_ut1.is_none() || evening_rise_jd_ut1.is_none()
        };
        if need_rise {
            if let Some(sr) = eval_daily_horizon_crossing(
                spk,
                time_table,
                body,
                observer,
                d_curr,
                true,
                &options,
            )? {
                let is_distinct_day = prev_rise
                    .map_or(true, |(_, p)| (sr.t_cross_ut1 - p.t_cross_ut1).abs() > 0.5);
                if is_distinct_day {
                    if let Some((d_prev, sr_prev)) = prev_rise {
                        if sr_prev.margin_deg * sr.margin_deg <= 0.0
                            && (sr.margin_deg - sr_prev.margin_deg).abs() < 40.0
                        {
                            let is_morning = sr.sun_ha_deg < 0.0 || sr_prev.sun_ha_deg < 0.0;
                            if !is_moon
                                && is_morning
                                && sr_prev.margin_deg < 0.0
                                && sr.margin_deg >= 0.0
                                && heliacal_rise_jd_ut1.is_none()
                            {
                                if let Some((t_exact, av_exact, vm_exact)) =
                                    refine_crossing(d_prev, d_curr, true, true)?
                                {
                                    if t_exact >= start_jd_ut1 {
                                        heliacal_rise_jd_ut1 = Some(t_exact);
                                        ref_av = av_exact;
                                        ref_vmag = vm_exact;
                                    }
                                }
                            } else if is_inner_or_moon
                                && is_morning
                                && sr_prev.margin_deg > 0.0
                                && sr.margin_deg <= 0.0
                                && morning_set_jd_ut1.is_none()
                            {
                                if let Some((t_exact, av_exact, vm_exact)) =
                                    refine_crossing(d_prev, d_curr, true, false)?
                                {
                                    if t_exact >= start_jd_ut1 {
                                        morning_set_jd_ut1 = Some(t_exact);
                                        if is_moon {
                                            ref_av = av_exact;
                                            ref_vmag = vm_exact;
                                        }
                                    }
                                }
                            } else if !is_inner_or_moon
                                && !is_morning
                                && sr_prev.margin_deg > 0.0
                                && sr.margin_deg <= 0.0
                                && evening_rise_jd_ut1.is_none()
                            {
                                if let Some((t_exact, _, _)) =
                                    refine_crossing(d_prev, d_curr, true, false)?
                                {
                                    if t_exact >= start_jd_ut1 {
                                        evening_rise_jd_ut1 = Some(t_exact);
                                    }
                                }
                            }
                        }
                    }
                    prev_rise = Some((d_curr, sr));
                }
            }
        }

        // Check Set events in the West:
        // - Outer bodies & stars: Evening Last (Heliacal Set) and Morning Last (Cosmical Set)
        // - Inner planets (Mercury, Venus) & Moon: Evening First (new crescent / appearance) and Evening Last (Heliacal Set)
        let need_set = if is_moon {
            evening_rise_jd_ut1.is_none()
        } else if is_inner_planet {
            heliacal_set_jd_ut1.is_none() || evening_rise_jd_ut1.is_none()
        } else {
            heliacal_set_jd_ut1.is_none() || morning_set_jd_ut1.is_none()
        };
        if need_set {
            if let Some(ss) = eval_daily_horizon_crossing(
                spk,
                time_table,
                body,
                observer,
                d_curr,
                false,
                &options,
            )? {
                let is_distinct_day = prev_set
                    .map_or(true, |(_, p)| (ss.t_cross_ut1 - p.t_cross_ut1).abs() > 0.5);
                if is_distinct_day {
                    if let Some((d_prev, ss_prev)) = prev_set {
                        if ss_prev.margin_deg * ss.margin_deg <= 0.0
                            && (ss.margin_deg - ss_prev.margin_deg).abs() < 40.0
                        {
                            let is_evening = ss.sun_ha_deg > 0.0 || ss_prev.sun_ha_deg > 0.0;
                            if !is_moon
                                && is_evening
                                && ss_prev.margin_deg > 0.0
                                && ss.margin_deg <= 0.0
                                && heliacal_set_jd_ut1.is_none()
                            {
                                if let Some((t_exact, _, _)) =
                                    refine_crossing(d_prev, d_curr, false, false)?
                                {
                                    if t_exact >= start_jd_ut1 {
                                        heliacal_set_jd_ut1 = Some(t_exact);
                                    }
                                }
                            } else if is_inner_or_moon
                                && is_evening
                                && ss_prev.margin_deg < 0.0
                                && ss.margin_deg >= 0.0
                                && evening_rise_jd_ut1.is_none()
                            {
                                if let Some((t_exact, av_exact, vm_exact)) =
                                    refine_crossing(d_prev, d_curr, false, true)?
                                {
                                    if t_exact >= start_jd_ut1 {
                                        evening_rise_jd_ut1 = Some(t_exact);
                                        if is_moon {
                                            ref_av = av_exact;
                                            ref_vmag = vm_exact;
                                        }
                                    }
                                }
                            } else if !is_inner_or_moon
                                && !is_evening
                                && ss_prev.margin_deg < 0.0
                                && ss.margin_deg >= 0.0
                                && morning_set_jd_ut1.is_none()
                            {
                                if let Some((t_exact, _, _)) =
                                    refine_crossing(d_prev, d_curr, false, true)?
                                {
                                    if t_exact >= start_jd_ut1 {
                                        morning_set_jd_ut1 = Some(t_exact);
                                    }
                                }
                            }
                        }
                    }
                    prev_set = Some((d_curr, ss));
                }
            }
        }

        let all_found = if is_moon {
            morning_set_jd_ut1.is_some() && evening_rise_jd_ut1.is_some()
        } else {
            heliacal_rise_jd_ut1.is_some()
                && morning_set_jd_ut1.is_some()
                && evening_rise_jd_ut1.is_some()
                && heliacal_set_jd_ut1.is_some()
        };
        if all_found {
            break;
        }
    }

    Ok(HeliacalResult {
        body,
        target_name: format_body_target_name(body),
        visual_magnitude: ref_vmag,
        arcus_visionis_deg: ref_av,
        heliacal_rise_jd_ut1,
        morning_set_jd_ut1,
        evening_rise_jd_ut1,
        heliacal_set_jd_ut1,
    })
}
