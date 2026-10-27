//! Fixed Stars catalog (50 `astra` stars + 2 Vedic anchor stars from `VizieR I/311` Hipparcos-2
//! at epoch `J1991.25` and `XHIP V/137D` radial velocities) and astrometric propagation
//! (`sofars::astro::pmpx` $\to$ `ldsun` $\to$ `ab` $\to$ `bp06`/`num06a`).

use crate::astrometry::{
    equatorial_to_horizontal_vector_au, AtmosphereParams, CalcResult, C_KM_S,
};
use crate::error::EphemerisError;
use crate::frame::{
    cart_state_to_spherical, cart_to_spherical, mat_vec, stencil_4pt_derivative, vec_dot, vec_norm,
    vec_scale, vec_sub, Vec3,
};
use crate::moment::{EpochContext, Moment, H_PLANET_DAYS};
use crate::types::{
    AbMode, BodyId, CalcOptions, Equinox, HipId, Origin, Plane, TargetKind, AU_KM, JD_J2000,
};
use serde::{Deserialize, Serialize};
use sofars::astro::{ab, ldsun, pmpx};
use std::sync::OnceLock;

/// Hipparcos catalog reference epoch `J1991.25` (`JD 2448349.0625`).
pub const JD_HIP_J1991_25: f64 = 2448349.0625;

/// Arcseconds to radians conversion factor ($\pi / (180 \times 3600)$).
pub const DAS2R: f64 = 4.848_136_811_095_36e-6;

/// Parsec in Astronomical Units ($1\text{ pc} = 3600 \times 180 / \pi\text{ AU} \approx 206264.806\text{ AU}$).
pub const PARSEC_AU: f64 = 206_264.806_247_096_36;

/// Quality warning flag emitted when evaluating a fixed star.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum StarQualityWarning {
    /// Catalog parallax $\varpi \le 0$ was clamped to $0.0\text{ mas}$ (pure angular propagation).
    NonPositiveParallax,
    /// Catalog radial velocity was missing (`0.0 km/s` fallback used).
    MissingRadialVelocity,
}

impl StarQualityWarning {
    /// Canonical string code for [`CalcResult::warnings`].
    #[must_use]
    pub const fn as_code(self) -> &'static str {
        match self {
            Self::NonPositiveParallax => "NonPositiveParallax",
            Self::MissingRadialVelocity => "MissingRadialVelocity",
        }
    }
}

/// Canonical fixed star catalog entry (`VizieR I/311` Hipparcos-2 + `XHIP V/137D`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StarEntry {
    pub index: u32,
    pub is_astra_50: bool,
    pub name: String,
    pub bayer: String,
    pub hip: u32,
    pub hr: u32,
    pub hd: u32,
    pub comp: String,
    pub ra_1991_deg: f64,
    pub de_1991_deg: f64,
    pub ra_2000_deg: f64,
    pub de_2000_deg: f64,
    pub pm_ra_cosdec_mas_yr: f64,
    pub pm_dec_mas_yr: f64,
    pub plx_mas: f64,
    pub e_plx_mas: f64,
    pub rv_kms: f64,
    pub vmag: f64,
    pub flags: String,
}

/// Public star catalog metadata (`R6.1`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StarMetadata {
    pub hip: u32,
    pub canonical_name: String,
    pub bayer: String,
    pub vmag: f64,
    pub ra_rad: f64,
    pub dec_rad: f64,
    pub pm_ra_mas_yr: f64,
    pub pm_dec_mas_yr: f64,
    pub plx_mas: f64,
    pub rv_kms: f64,
    pub rv_missing: bool,
}

impl StarEntry {
    /// Returns the [`HipId`] of this star.
    #[must_use]
    pub const fn hip_id(&self) -> HipId {
        HipId(self.hip)
    }

    /// Converts this catalog entry into [`StarMetadata`] (`R6.1`).
    #[must_use]
    pub fn to_metadata(&self) -> StarMetadata {
        StarMetadata {
            hip: self.hip,
            canonical_name: self.name.clone(),
            bayer: self.bayer.clone(),
            vmag: self.vmag,
            ra_rad: self.ra_2000_deg.to_radians(),
            dec_rad: self.de_2000_deg.to_radians(),
            pm_ra_mas_yr: self.pm_ra_cosdec_mas_yr,
            pm_dec_mas_yr: self.pm_dec_mas_yr,
            plx_mas: self.plx_mas,
            rv_kms: self.rv_kms,
            rv_missing: self.flags.contains("MissingRadialVelocity"),
        }
    }

    /// Returns any quality warnings for this star (e.g., `NonPositiveParallax` if `plx_mas <= 0`).
    #[must_use]
    pub fn quality_warnings(&self) -> Vec<StarQualityWarning> {
        let mut w = Vec::new();
        if self.plx_mas <= 0.0 || self.flags.contains("NonPositiveParallax") {
            w.push(StarQualityWarning::NonPositiveParallax);
        }
        if self.flags.contains("MissingRadialVelocity") {
            w.push(StarQualityWarning::MissingRadialVelocity);
        }
        w
    }

    /// Nominal distance in AU derived from parallax (`1000 / plx_mas * PARSEC_AU`, or `1e9 AU` if `plx_mas <= 0.1`).
    #[must_use]
    pub fn nominal_distance_au(&self) -> f64 {
        if self.plx_mas > 0.1 {
            (1000.0 / self.plx_mas) * PARSEC_AU
        } else {
            1.0e9
        }
    }
}

/// Returns the embedded 52-star catalog (50 `astra` stars + $\zeta$ Psc A and $\delta$ Cnc).
#[must_use]
pub fn star_catalog_52() -> &'static [StarEntry] {
    static CATALOG: OnceLock<Vec<StarEntry>> = OnceLock::new();
    CATALOG.get_or_init(|| {
        let raw = include_str!("stars_catalog_52.json");
        serde_json::from_str(raw).unwrap_or_default()
    })
}

/// Supplementary embedded Hipparcos-2 (`VizieR I/311` + `XHIP V/137D`) catalog entry for
/// `Cor Caroli` ($\alpha^2\text{ CVn}$, `HIP 63125`, `HR 4915`, `HD 112413`).
fn cor_caroli_entry() -> &'static StarEntry {
    static COR_CAROLI: OnceLock<StarEntry> = OnceLock::new();
    COR_CAROLI.get_or_init(|| StarEntry {
        index: 53,
        is_astra_50: false,
        name: "Cor Caroli".to_owned(),
        bayer: "al2CVn".to_owned(),
        hip: 63125,
        hr: 4915,
        hd: 112413,
        comp: "A".to_owned(),
        ra_1991_deg: 194.007_618_1,
        de_1991_deg: 38.318_246_0,
        ra_2000_deg: 194.006_896_4,
        de_2000_deg: 38.318_379_7,
        pm_ra_cosdec_mas_yr: -232.91,
        pm_dec_mas_yr: 54.99,
        plx_mas: 28.39,
        e_plx_mas: 0.55,
        rv_kms: -3.3,
        vmag: 2.89,
        flags: String::new(),
    })
}

/// Looks up a [`StarEntry`] by its Hipparcos number (`HipId`).
#[must_use]
pub fn find_star_by_hip(hip: HipId) -> Option<&'static StarEntry> {
    if hip.0 == 63125 {
        return Some(cor_caroli_entry());
    }
    star_catalog_52().iter().find(|s| s.hip == hip.0)
}

/// Looks up a [`StarEntry`] by case-insensitive canonical `name`, `bayer` code, or `"hip:<id>"`.
#[must_use]
pub fn find_star_by_name_or_bayer(query: &str) -> Option<&'static StarEntry> {
    let trimmed = query.trim().trim_start_matches(',');
    let prefix = trimmed.split(',').next().unwrap_or(trimmed).trim();
    if let Some(hip_s) = prefix.to_ascii_lowercase().strip_prefix("hip:") {
        if let Ok(hip_num) = hip_s.trim().parse::<u32>() {
            return find_star_by_hip(HipId(hip_num));
        }
    }
    if let Ok(idx_num) = prefix.parse::<u32>() {
        if idx_num == 53 {
            return Some(cor_caroli_entry());
        }
        if let Some(s) = star_catalog_52().iter().find(|s| s.index == idx_num) {
            return Some(s);
        }
    }
    if let Some(s) = star_catalog_52().iter().find(|s| {
        s.name.eq_ignore_ascii_case(prefix)
            || s.bayer.eq_ignore_ascii_case(prefix)
            || s.name.eq_ignore_ascii_case(trimmed)
            || s.bayer.eq_ignore_ascii_case(trimmed)
    }) {
        return Some(s);
    }
    let cc = cor_caroli_entry();
    if cc.name.eq_ignore_ascii_case(prefix)
        || cc.bayer.eq_ignore_ascii_case(prefix)
        || cc.name.eq_ignore_ascii_case(trimmed)
        || cc.bayer.eq_ignore_ascii_case(trimmed)
        || prefix.eq_ignore_ascii_case("alCVn")
        || trimmed.eq_ignore_ascii_case("alCVn")
    {
        return Some(cc);
    }
    None
}

/// Returns [`StarMetadata`] for a star looked up by canonical name, Bayer designation, or `"hip:<id>"` (`R6.1`).
pub fn star_info(name_or_bayer: &str) -> Result<StarMetadata, EphemerisError> {
    let star = find_star_by_name_or_bayer(name_or_bayer).ok_or_else(|| {
        EphemerisError::UnsupportedBody {
            body: format!("FixedStar({name_or_bayer})"),
        }
    })?;
    Ok(star.to_metadata())
}

/// Returns [`StarMetadata`] for all 52 stars in the embedded catalog (`R6.1`).
#[must_use]
pub fn list_stars() -> Vec<StarMetadata> {
    star_catalog_52().iter().map(StarEntry::to_metadata).collect()
}

/// Evaluates the mean ecliptic longitude of date `[0, 360)` of `star` with proper motion and
/// geocentric annual parallax (`sofars::astro::pmpx`) and **zero annual aberration**.
///
/// Used by the 3 star-anchored Vedic ayanamsas (`true_citra`, `true_revati`, `true_pushya`;
/// `MODEL_PASSPORTS.md` §1.2).
#[must_use]
pub fn star_mean_ecliptic_lon_no_ab_at_ctx(star: &StarEntry, ctx: &EpochContext) -> f64 {
    let jd_tdb = ctx.tdb.jd_f64();
    let p_eb_au = vec_scale(&ctx.earth_bary_km.pos_km, 1.0 / AU_KM);
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
    let pmt = (jd_tdb - JD_HIP_J1991_25) / 365.25;

    let u_icrs = pmpx(ra1, de1, pmr, pmd, px, star.rv_kms, pmt, p_eb_au);
    let u_ecl = mat_vec(&ctx.frame.r_ecl_mean, &u_icrs);
    u_ecl[1].atan2(u_ecl[0]).to_degrees().rem_euclid(360.0)
}

/// Evaluates the unit direction vector of `star` in the output frame requested by `options` at `ctx`.
#[must_use]
pub fn eval_star_unit_vector(
    star: &StarEntry,
    ctx: &EpochContext,
    options: &CalcOptions,
) -> Vec3 {
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

    let u_icrs = if options.origin == Origin::Barycentric
        && options.equinox == Equinox::J2000
        && options.aberration == AbMode::Geometric
    {
        // Pure catalog ICRS at J2000.0 (space motion from J1991.25 to J2000.0, zero parallax)
        let pmt = (JD_J2000 - JD_HIP_J1991_25) / 365.25;
        pmpx(ra1, de1, pmr, pmd, 0.0, star.rv_kms, pmt, [0.0, 0.0, 0.0])
    } else {
        let jd_tdb = ctx.tdb.jd_f64();
        let (p_obs_km, v_obs_kms) = match options.origin {
            Origin::Geocentric => (ctx.earth_bary_km.pos_km, ctx.earth_bary_km.vel_km_s),
            Origin::Topocentric(_) => (ctx.observer_bary_km.pos_km, ctx.observer_bary_km.vel_km_s),
            Origin::Heliocentric => (ctx.sun_bary_km.pos_km, ctx.sun_bary_km.vel_km_s),
            Origin::Barycentric => ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
        };
        let p_eb_au = vec_scale(&p_obs_km, 1.0 / AU_KM);
        let px = if options.origin == Origin::Barycentric || options.aberration == AbMode::Geometric
        {
            0.0
        } else {
            (star.plx_mas * 1e-3).max(0.0)
        };
        let pmt = (jd_tdb - JD_HIP_J1991_25) / 365.25;
        let u_ppx = pmpx(ra1, de1, pmr, pmd, px, star.rv_kms, pmt, p_eb_au);

        if matches!(
            options.aberration,
            AbMode::Apparent | AbMode::NoAberration | AbMode::NoDeflection
        ) && matches!(options.origin, Origin::Geocentric | Origin::Topocentric(_))
        {
            let e_sun_vec =
                vec_scale(&vec_sub(&p_obs_km, &ctx.sun_bary_km.pos_km), 1.0 / AU_KM);
            let em = vec_norm(&e_sun_vec);
            let e_sun_unit = if em > 1e-12 {
                vec_scale(&e_sun_vec, 1.0 / em)
            } else {
                [1.0, 0.0, 0.0]
            };
            let u_ld = if matches!(options.aberration, AbMode::Apparent | AbMode::NoAberration) {
                ldsun(u_ppx, e_sun_unit, em)
            } else {
                u_ppx
            };
            if matches!(options.aberration, AbMode::Apparent | AbMode::NoDeflection) {
                let v2_c2 =
                    (vec_dot(&v_obs_kms, &v_obs_kms) / (C_KM_S * C_KM_S)).min(0.999_999);
                let bm1 = (1.0 - v2_c2).sqrt();
                let v_ab = vec_scale(&v_obs_kms, 1.0 / C_KM_S);
                ab(&u_ld, &v_ab, em, bm1)
            } else {
                u_ld
            }
        } else {
            u_ppx
        }
    };

    let r_mat = ctx
        .frame
        .icrf_to_frame_matrix(options.plane, options.equinox, options.nutation);
    let vec_frame = mat_vec(&r_mat, &u_icrs);
    if options.plane == Plane::Horizontal {
        let obs_opt = match options.origin {
            Origin::Topocentric(o) => Some(o),
            _ => ctx.observer,
        };
        if let Some(obs) = obs_opt {
            return equatorial_to_horizontal_vector_au(
                obs,
                ctx.gast_rad,
                &vec_frame,
                AtmosphereParams::default(),
            );
        }
    }
    vec_frame
}

/// Evaluates a fixed star (`HipId`) inside a [`Moment`].
pub fn calc_star_in_moment(
    moment: &Moment<'_>,
    hip: HipId,
    options: &CalcOptions,
    sidereal_offset_and_rate: Option<(f64, f64)>,
) -> Result<CalcResult, EphemerisError> {
    let star = find_star_by_hip(hip).ok_or_else(|| EphemerisError::UnsupportedBody {
        body: format!("FixedStar(HIP {})", hip.0),
    })?;
    calc_star_entry_in_moment(moment, star, options, sidereal_offset_and_rate)
}

/// Evaluates an arbitrary [`StarEntry`] (from the embedded catalog or a custom star catalog) inside a [`Moment`].
pub fn calc_star_entry_in_moment(
    moment: &Moment<'_>,
    star: &StarEntry,
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
    let hip = star.hip_id();
    let dist_au = star.nominal_distance_au();
    let u_0 = eval_star_unit_vector(star, center, options);
    let pos_0 = vec_scale(&u_0, dist_au);

    let (vel_0, mut sph) = if options.include_speed {
        let stencil = moment.planet_stencil_contexts()?;
        let mut stencil_opts = *options;
        if options.plane == Plane::Ecliptic && sidereal_offset_and_rate.is_some() {
            stencil_opts.nutation = crate::types::NutationMode::Mean;
        }
        let pos_m2 = vec_scale(&eval_star_unit_vector(star, &stencil[0], &stencil_opts), dist_au);
        let pos_m1 = vec_scale(&eval_star_unit_vector(star, &stencil[1], &stencil_opts), dist_au);
        let pos_p1 = vec_scale(&eval_star_unit_vector(star, &stencil[2], &stencil_opts), dist_au);
        let pos_p2 = vec_scale(&eval_star_unit_vector(star, &stencil[3], &stencil_opts), dist_au);
        let v0 = stencil_4pt_derivative(&pos_m2, &pos_m1, &pos_p1, &pos_p2, H_PLANET_DAYS);
        let pos_for_rate = if stencil_opts.nutation != options.nutation {
            vec_scale(&eval_star_unit_vector(star, center, &stencil_opts), dist_au)
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
        body: BodyId::FixedStar(hip),
        target_id: hip.0 as i32,
        center_id,
        target_kind: TargetKind::CatalogStar,
        jd_tt: center.tt.jd_f64(),
        jd_ut1: center.ut1.jd_f64(),
        jd_tdb: center.tdb.jd_f64(),
        delta_t_sec: center.delta_t_info.delta_t_sec,
        time_extrapolated: center.delta_t_info.time_extrapolated,
        time_table_version: center.delta_t_info.time_table_version.clone(),
        lon_deg: sph.lon_deg,
        lat_deg: sph.lat_deg,
        dist_au,
        dist_km: dist_au * AU_KM,
        lon_speed_deg_day: sph.lon_speed_deg_day,
        lat_speed_deg_day: sph.lat_speed_deg_day,
        dist_speed_au_day: sph.dist_speed_au_day,
        dist_speed_km_day: sph.dist_speed_au_day * AU_KM,
        magnitude_v: Some(star.vmag),
        cartesian_au: cart_pos,
        cartesian_vel_au_day: cart_vel,
        warnings: star
            .quality_warnings()
            .into_iter()
            .map(|w| w.as_code().to_owned())
            .collect(),
    })
}

/// Parses a custom star catalog from either JSON (`Vec<StarEntry>`) or `sefstars.txt`-compatible
/// comma-separated catalog text (`name, bayer, epoch, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s, pm_ra, pm_dec, rv, plx, vmag, ...`).
pub fn parse_sefstars_catalog(text: &str) -> Result<Vec<StarEntry>, EphemerisError> {
    let trimmed = text.trim();
    if trimmed.starts_with('[') {
        return serde_json::from_str::<Vec<StarEntry>>(trimmed).map_err(|e| {
            EphemerisError::InvalidInput {
                message: format!("Failed to parse JSON custom star catalog: {e}"),
            }
        });
    }

    let mut entries = Vec::new();
    for (line_idx, raw_line) in text.lines().enumerate() {
        let line = raw_line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split(',').map(str::trim).collect();
        if cols.len() < 14 {
            continue;
        }
        let name = cols.first().copied().unwrap_or("").to_owned();
        let bayer = cols.get(1).copied().unwrap_or("").to_owned();
        if name.is_empty() && bayer.is_empty() {
            continue;
        }
        let parse_f64 = |idx: usize| -> f64 {
            cols.get(idx)
                .and_then(|s| s.parse::<f64>().ok())
                .unwrap_or(0.0)
        };
        let ra_h = parse_f64(3);
        let ra_m = parse_f64(4);
        let ra_s = parse_f64(5);
        let ra_deg = (ra_h.abs() + ra_m.abs() / 60.0 + ra_s.abs() / 3600.0) * 15.0;

        let dec_str = cols.get(6).copied().unwrap_or("0");
        let dec_sign = if dec_str.starts_with('-') { -1.0 } else { 1.0 };
        let dec_d = parse_f64(6).abs();
        let dec_m = parse_f64(7).abs();
        let dec_s = parse_f64(8).abs();
        let de_deg = dec_sign * (dec_d + dec_m / 60.0 + dec_s / 3600.0);

        // In sefstars.txt: col 9 is pm_ra in seconds of RA/century or mas/yr, col 10 is pm_dec in arcsec/cy or mas/yr
        // If |pm| < 20.0 and epoch == 2000.0, sefstars stores pm_ra in 0.1s/cy or arcsec/cy; we support standard mas/yr or sefstars convention:
        let pm_ra_raw = parse_f64(9);
        let pm_dec_raw = parse_f64(10);
        let rv_kms = parse_f64(11);
        let plx_mas = parse_f64(12);
        let vmag = parse_f64(13);
        let hip = cols
            .get(16)
            .or_else(|| cols.get(15))
            .and_then(|s| s.parse::<u32>().ok())
            .unwrap_or(900_000 + line_idx as u32);

        // Convert J2000 position back by 8.75 yr to J1991.25 so pmpx propagates from J1991.25 consistently
        let dt_yr = (JD_J2000 - JD_HIP_J1991_25) / 365.25;
        let cos_dec = de_deg.to_radians().cos().abs().max(1.0e-6);
        let ra_1991_deg = (ra_deg - (pm_ra_raw / (3_600_000.0 * cos_dec)) * dt_yr).rem_euclid(360.0);
        let de_1991_deg = (de_deg - (pm_dec_raw / 3_600_000.0) * dt_yr).clamp(-90.0, 90.0);

        entries.push(StarEntry {
            index: (line_idx + 1) as u32,
            is_astra_50: false,
            name,
            bayer,
            hip,
            hr: 0,
            hd: 0,
            comp: String::new(),
            ra_1991_deg,
            de_1991_deg,
            ra_2000_deg: ra_deg,
            de_2000_deg: de_deg,
            pm_ra_cosdec_mas_yr: pm_ra_raw,
            pm_dec_mas_yr: pm_dec_raw,
            plx_mas,
            e_plx_mas: 0.0,
            rv_kms,
            vmag,
            flags: String::new(),
        });
    }

    Ok(entries)
}

/// Loads a custom star catalog from disk (`JSON` or `sefstars.txt`-compatible text).
pub fn load_custom_stars_file(
    path: impl AsRef<std::path::Path>,
) -> Result<Vec<StarEntry>, EphemerisError> {
    let p = path.as_ref();
    let text = std::fs::read_to_string(p).map_err(|e| EphemerisError::MissingData {
        detail: format!("Failed to read custom star catalog {}: {e}", p.display()),
    })?;
    parse_sefstars_catalog(&text)
}
