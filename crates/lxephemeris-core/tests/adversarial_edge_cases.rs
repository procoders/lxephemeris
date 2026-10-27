//! Adversarial & Pathological Numerical Edge-Case Torture Suite (Wave A).
//!
//! Verifies all 20 findings from the Codex (`gpt-6.1-sol`) and Opus (`pro`) Pre-v0.2
//! Open-Source & Computational Astronomy audits.
#![allow(clippy::indexing_slicing, clippy::field_reassign_with_default)]

use lxephemeris_core::{
    cart_state_to_spherical, cart_to_spherical, compute_lunar_saros, compute_solar_saros,
    eval_star_unit_vector, gauquelin_sectors, house_of_lon, houses_armc_with_sidereal, julday,
     revjul, star_mean_ecliptic_lon_no_ab_at_ctx, AbMode, AyanamsaModel, BodyId, CalcOptions,
    CalendarKind, EclipseOptions, Engine, EngineConfig, EphemerisError, Equinox, HeliacalOptions,
    HipId, HouseOptions, HouseSystem, HypotheticalId, LunarPointId, MpcNumber, NutationMode,
    Observer, Origin, Plane, PlanetId, PolarPolicy, RiseTransEvent, RiseTransOptions, SplitJd,
    StarEntry, TimeInput, TimeTable, Zodiac, JD_J2000,
};
use std::path::Path;

fn workspace_engine() -> Engine {
    let ws = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root");
    Engine::from_data_dir(ws.join("data"), EngineConfig::default()).expect("load engine")
}

#[test]
fn test_observer_altitude_bounds_and_non_finite_rejection() {
    // Pathological altitudes (1e308, NaN, Inf, -20 km) must be rejected before overflowing gd2gc/dtdb
    assert!(Observer::new(0.0, 0.0, 1.0e308).is_err());
    assert!(Observer::new(0.0, 0.0, f64::NAN).is_err());
    assert!(Observer::new(0.0, 0.0, f64::INFINITY).is_err());
    assert!(Observer::new(0.0, 0.0, -20_000.0).is_err());
    assert!(Observer::new(90.0001, 0.0, 0.0).is_err());
    assert!(Observer::new(0.0, 360.0001, 0.0).is_err());

    // Valid physical extremes (Mariana Trench -11 km and geostationary orbit 35,786 km) must succeed
    let trench = Observer::new(11.35, 142.2, -11_000.0).expect("Mariana trench observer");
    let geo_orbit = Observer::new(0.0, 0.0, 35_786_000.0).expect("GEO observer");

    let eng = workspace_engine();
    let m_trench = eng
        .at(TimeInput::from_tt_jd(JD_J2000), Some(trench))
        .expect("trench moment");
    let m_geo = eng
        .at(TimeInput::from_tt_jd(JD_J2000), Some(geo_orbit))
        .expect("geo moment");
    assert!(m_trench.center().observer_bary_km.pos_km[0].is_finite());
    assert!(m_geo.center().observer_bary_km.pos_km[0].is_finite());
}

#[test]
fn test_hostile_daf_type21_metadata_and_kq_bounds() {
    use lxephemeris_core::daf::eval_spk_type21;

    // 1. Fewer than 4 words
    assert!(matches!(
        eval_spk_type21(&[1.0, 2.0], 0.0, 20_000_001),
        Err(EphemerisError::CorruptKernel { .. })
    ));

    // 2. Non-finite et_s
    assert!(matches!(
        eval_spk_type21(&[0.0; 64], f64::NAN, 20_000_001),
        Err(EphemerisError::InvalidInput { .. })
    ));

    // 3. Fractional or huge N / MAXDIM
    let mut bad_meta = vec![0.0_f64; 80];
    let len = bad_meta.len();
    bad_meta[len - 2] = 15.5; // fractional MAXDIM
    bad_meta[len - 1] = 1.0;
    assert!(matches!(
        eval_spk_type21(&bad_meta, 0.0, 20_000_001),
        Err(EphemerisError::CorruptKernel { .. })
    ));

    // 4. Valid N=1, MAXDIM=15 (dlsize = 71), but hostile KQ[0] = 999999.0 or KQ[0] >= KQMAX1
    let maxdim = 15usize;
    let dlsize = 4 * maxdim + 11; // 71
    let mut hostile_rec = vec![0.0_f64; dlsize + 1 + 4];
    let total = hostile_rec.len();
    hostile_rec[total - 2] = maxdim as f64;
    hostile_rec[total - 1] = 1.0; // n_rec = 1
    hostile_rec[4 * maxdim + 7] = 10.0; // KQMAX1 = 10
    hostile_rec[4 * maxdim + 8] = 10.0; // KQ[0] = 10 >= KQMAX1 (invalid!)
    hostile_rec[4 * maxdim + 9] = 5.0;
    hostile_rec[4 * maxdim + 10] = 5.0;
    assert!(matches!(
        eval_spk_type21(&hostile_rec, 0.0, 20_000_001),
        Err(EphemerisError::CorruptKernel { .. })
    ));
}

#[test]
fn test_minor_planet_heliacal_naif_id_and_option_validation() {
    let eng = workspace_engine();
    let obs = Observer::new(30.0, 31.2, 50.0).expect("Cairo observer");

    // Ceres (MPC 1 -> NAIF 2000001) and Vesta (MPC 4 -> NAIF 2000004) must succeed
    let ceres_hel = eng
        .heliacal_phenomena(
            BodyId::MinorPlanet(MpcNumber(1)),
            JD_J2000,
            obs,
            HeliacalOptions::default(),
        )
        .expect("Ceres heliacal phenomena");
    assert_eq!(ceres_hel.target_name, "Ceres");
    assert!(ceres_hel.visual_magnitude.is_finite());
    assert!((6.0..10.5).contains(&ceres_hel.visual_magnitude));
    assert!(ceres_hel.arcus_visionis_deg.is_finite());

    // Invalid heliacal options (NaN start_jd_ut1, negative max_search_days, zero snellen_ratio) must fail cleanly
    assert!(eng
        .heliacal_phenomena(
            BodyId::Planet(PlanetId::Venus),
            f64::NAN,
            obs,
            HeliacalOptions::default()
        )
        .is_err());
    let mut bad_opts = HeliacalOptions::default();
    bad_opts.max_search_days = -10.0;
    assert!(eng
        .heliacal_phenomena(BodyId::Planet(PlanetId::Venus), JD_J2000, obs, bad_opts)
        .is_err());
    let mut bad_snellen = HeliacalOptions::default();
    bad_snellen.snellen_ratio = 0.0;
    assert!(eng
        .heliacal_phenomena(BodyId::Planet(PlanetId::Venus), JD_J2000, obs, bad_snellen)
        .is_err());
}

#[test]
fn test_split_jd_papal_gap_negative_jd_revjul_and_leapseconds_parser() {
    // 1. SplitJd overflow and NaN protection
    let overflow_split = SplitJd::new(1.0e308, 1.0e308);
    assert!(!overflow_split.is_valid_jd());
    let nan_split = SplitJd::new(f64::NAN, 0.0);
    assert!(!nan_split.is_valid_jd());

    // 2. Papal Gregorian reform gap (1582-10-05..14) rejected under AstroAuto
    for day in 5..=14 {
        assert!(
            julday(1582, 10, day, 12.0, CalendarKind::AstroAuto).is_err(),
            "1582-10-{day:02} must be rejected in AstroAuto"
        );
    }
    let jd_oct4 = julday(1582, 10, 4, 0.0, CalendarKind::AstroAuto).expect("1582-10-04");
    let jd_oct15 = julday(1582, 10, 15, 0.0, CalendarKind::AstroAuto).expect("1582-10-15");
    assert!((jd_oct15 - jd_oct4 - 1.0).abs() < 1e-12);

    // 3. Negative Julian Dates (JD < 0, Year <= -4713) round-trip in both Julian and Gregorian
    for &(year, month, day, hour) in &[
        (-4713, 1, 1, 12.0),
        (-4715, 6, 15, 6.0),
        (-5000, 3, 21, 0.0),
        (-10000, 12, 31, 18.0),
    ] {
        for cal in [CalendarKind::Julian, CalendarKind::Gregorian] {
            let jd = julday(year, month, day, hour, cal).expect("julday negative JD");
            assert!(jd < 0.0);
            let rev = revjul(jd, cal).expect("revjul negative JD");
            assert_eq!(
                (rev.year, rev.month, rev.day),
                (year, month, day),
                "Mismatch for {year}-{month}-{day} in {cal:?} (jd={jd})"
            );
            assert!((rev.hour - hour).abs() < 1e-8);
        }
    }

    // 4. Malformed leap-second lines where "TAI-UTC=" precedes "=JD" must not panic
    let malformed = "TAI-UTC= 37.0 =JD 2457754.5\n2457754.5 =JD 2457754.5 TAI-UTC= 37.0 S\n";
    let parsed = TimeTable::parse_leapseconds_table(malformed);
    assert_eq!(parsed.len(), 1);
    assert!((parsed[0].1 - 37.0).abs() < 1e-12);
}

#[test]
fn test_vector_norms_nan_propagation_and_celestial_pole_stars() {
    use lxephemeris_core::frame::{vec_norm, vec_unit};

    // Extreme finite vectors must not overflow x*x to Inf or underflow to 0
    let big = [1.0e200, 1.0e200, 1.0e200];
    let n_big = vec_norm(&big);
    assert!(n_big.is_finite());
    assert!((n_big - 3.0_f64.sqrt() * 1.0e200).abs() / n_big < 1e-14);
    let u_big = vec_unit(&big);
    assert!((u_big[0] - 1.0 / 3.0_f64.sqrt()).abs() < 1e-14);

    let tiny = [1.0e-200, 1.0e-200, 1.0e-200];
    let n_tiny = vec_norm(&tiny);
    assert!(n_tiny > 0.0 && n_tiny.is_finite());

    // NaN must propagate in vec_unit and cart_to_spherical instead of silently becoming 0.0
    let u_nan = vec_unit(&[f64::NAN, 1.0, 0.0]);
    assert!(u_nan[0].is_nan() && u_nan[1].is_nan() && u_nan[2].is_nan());
    let sph_nan = cart_to_spherical(&[f64::NAN, 1.0, 0.0]);
    assert!(sph_nan.lon_deg.is_nan() && sph_nan.lat_deg.is_nan() && sph_nan.dist.is_nan());

    // Exact celestial pole star (dec = +90.0 deg and -90.0 deg) must not divide by zero in pmr
    let eng = workspace_engine();
    let m = eng.at(TimeInput::from_tt_jd(JD_J2000), None).expect("moment");
    for dec_pole in [90.0, -90.0] {
        let pole_star = StarEntry {
            index: 999,
            is_astra_50: false,
            name: "PolarTest".to_owned(),
            bayer: "alpha Pol".to_owned(),
            hip: 999_999,
            hr: 0,
            hd: 0,
            comp: String::new(),
            ra_1991_deg: 45.0,
            de_1991_deg: dec_pole,
            ra_2000_deg: 45.0,
            de_2000_deg: dec_pole,
            pm_ra_cosdec_mas_yr: 25.0,
            pm_dec_mas_yr: -15.0,
            plx_mas: 10.0,
            e_plx_mas: 0.1,
            rv_kms: 5.0,
            vmag: 2.0,
            flags: String::new(),
        };
        let u = eval_star_unit_vector(&pole_star, m.center(), &CalcOptions::default());
        assert!(u[0].is_finite() && u[1].is_finite() && u[2].is_finite());
        let lon_no_ab = star_mean_ecliptic_lon_no_ab_at_ctx(&pole_star, m.center());
        assert!(lon_no_ab.is_finite() && (0.0..360.0).contains(&lon_no_ab));
    }
}

#[test]
fn test_houses_apc_polar_singularity_pullen_sr_and_mean_nutation_gmst() {
    // 1. APC ('Y') at the exact polar circle (lat = 90 - eps, RAMC = 270 deg) where 1 + tan(lat)*tan(eps)*sin(RAMC) == 0
    let eps_deg = 23.439_291_11;
    let lat_polar_circle = 90.0 - eps_deg;
    let apc_opts = HouseOptions {
        system: HouseSystem::Apc,
        polar_policy: PolarPolicy::FallbackToWholeSign,
        zodiac: Zodiac::Tropical,
        nutation: NutationMode::True,
        include_speed: true,
    };
    let apc_res = houses_armc_with_sidereal(270.0, lat_polar_circle, eps_deg, 0.0, 0.0, apc_opts)
        .expect("APC at polar circle singularity");
    for c in apc_res.cusps_deg {
        assert!(c.is_finite() && (0.0..360.0).contains(&c));
    }

    // 2. Pullen Sinusoidal Ratio ('Q') across extreme latitudes/angles
    let pullen_opts = HouseOptions {
        system: HouseSystem::PullenSinusoidalRatio,
        polar_policy: PolarPolicy::Strict,
        zodiac: Zodiac::Tropical,
        nutation: NutationMode::True,
        include_speed: true,
    };
    for lat in [0.0, 45.0, 66.56, 85.0] {
        for armc in [0.0, 90.0, 180.0, 270.0] {
            let q_res = houses_armc_with_sidereal(armc, lat, eps_deg, 0.0, 0.0, pullen_opts)
                .expect("Pullen SR");
            for c in q_res.cusps_deg {
                assert!(c.is_finite() && (0.0..360.0).contains(&c));
            }
        }
    }

    // 3. Mean houses (NutationMode::Mean / SEFLG_NONUT) must use GMST instead of GAST
    let eng = workspace_engine();
    let obs = Observer::new(51.5, 0.0, 0.0).expect("Greenwich");
    let h_true = eng
        .houses(
            TimeInput::from_ut1_jd(JD_J2000),
            obs,
            HouseOptions {
                system: HouseSystem::Placidus,
                polar_policy: PolarPolicy::Strict,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::True,
                include_speed: false,
            },
        )
        .expect("True houses");
    let h_mean = eng
        .houses(
            TimeInput::from_ut1_jd(JD_J2000),
            obs,
            HouseOptions {
                system: HouseSystem::Placidus,
                polar_policy: PolarPolicy::Strict,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::Mean,
                include_speed: false,
            },
        )
        .expect("Mean houses");
    // At J2000, |GAST - GMST| is ~12.78 arcseconds (~0.00355 deg)
    let armc_diff_arcsec = (h_true.angles.armc_deg - h_mean.angles.armc_deg).abs() * 3_600.0;
    assert!(
        (10.0..16.0).contains(&armc_diff_arcsec),
        "Expected equation of equinoxes difference (~12.78\"), got {armc_diff_arcsec}\""
    );

    // Also verify calc_houses_in_moment matches calc_houses for both True and Mean nutation
    let m = eng
        .at(TimeInput::from_ut1_jd(JD_J2000), Some(obs))
        .expect("moment");
    let hm_mean = eng
        .houses_in_moment(
            &m,
            obs,
            HouseOptions {
                system: HouseSystem::Placidus,
                polar_policy: PolarPolicy::Strict,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::Mean,
                include_speed: false,
            },
        )
        .expect("Mean houses in moment");
    assert!((hm_mean.angles.armc_deg - h_mean.angles.armc_deg).abs() < 1e-12);

    // 4. Non-finite inputs to houses_armc_with_sidereal, gauquelin_sectors, and house_of_lon
    assert!(
        houses_armc_with_sidereal(0.0, 45.0, 23.44, f64::NAN, 0.0, HouseOptions::default())
            .is_err()
    );
    assert!(
        houses_armc_with_sidereal(0.0, 45.0, 23.44, 0.0, f64::NAN, HouseOptions::default())
            .is_err()
    );
    assert!(gauquelin_sectors(f64::NAN, 45.0, 23.44, PolarPolicy::Strict).is_err());
    assert!(gauquelin_sectors(0.0, 90.0, 23.44, PolarPolicy::FallbackToPorphyry).is_err());
    assert_eq!(house_of_lon(f64::NAN, &h_true.cusps_deg, 1e-9), 1);
}

#[test]
fn test_sidereal_cartesian_vector_and_velocity_consistency() {
    let eng = workspace_engine();
    let m = eng
        .at(TimeInput::from_tt_jd(JD_J2000), None)
        .expect("moment at J2000");

    let sid_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Ecliptic,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
        include_speed: true,
    };

    for body in [
        BodyId::Planet(PlanetId::Sun),
        BodyId::Planet(PlanetId::Mars),
        BodyId::MinorPlanet(MpcNumber(1)),
        BodyId::LunarPoint(LunarPointId::MeanNode),
        BodyId::LunarPoint(LunarPointId::TrueLilith),
        BodyId::Hypothetical(HypotheticalId::Cupido),
        BodyId::FixedStar(HipId(65474)), // Spica
    ] {
        let res = eng
            .calc_in_moment(&m, body, sid_opts)
            .expect("sidereal calc");
        let sph_from_cart =
            cart_state_to_spherical(&res.cartesian_au, &res.cartesian_vel_au_day);
        let lon_diff =
            (sph_from_cart.lon_deg - res.lon_deg + 180.0).rem_euclid(360.0) - 180.0;
        assert!(
            lon_diff.abs() < 1e-11,
            "Sidereal Cartesian longitude mismatch for {body:?}: sph={} vs cart={}",
            res.lon_deg,
            sph_from_cart.lon_deg
        );
        assert!(
            (sph_from_cart.lat_deg - res.lat_deg).abs() < 1e-11,
            "Sidereal Cartesian latitude mismatch for {body:?}"
        );
        assert!(
            (sph_from_cart.lon_speed_deg_day - res.lon_speed_deg_day).abs() < 1e-9,
            "Sidereal Cartesian lon_speed mismatch for {body:?}: sph={} vs cart={}",
            res.lon_speed_deg_day,
            sph_from_cart.lon_speed_deg_day
        );
    }
}

#[test]
fn test_saros_overflow_eclipse_ranges_and_rise_trans_bounds() {
    // 1. Extreme dates/gammas in compute_solar_saros / compute_lunar_saros must never panic or overflow
    assert_eq!(compute_solar_saros(f64::NAN, 0.5), (0, 1));
    assert_eq!(compute_solar_saros(1.0e300, 0.5), (0, 1));
    assert_eq!(compute_lunar_saros(f64::NEG_INFINITY, 0.5), (0, 1));
    assert_eq!(compute_lunar_saros(1.0e300, -1.0e300), (0, 1));

    let eng = workspace_engine();
    let obs = Observer::new(51.5, 0.0, 0.0).expect("Greenwich");

    // 2. Eclipse range / search validation
    assert!(eng
        .solar_eclipses_in_range(
            TimeInput::from_tt_jd(f64::NAN),
            TimeInput::from_tt_jd(JD_J2000),
            None,
            EclipseOptions::default()
        )
        .is_err());
    let mut bad_ecl = EclipseOptions::default();
    bad_ecl.max_synodic_months = 0;
    assert!(eng
        .next_solar_eclipse(TimeInput::from_tt_jd(JD_J2000), false, None, bad_ecl)
        .is_err());

    // 3. rise_trans option validation
    let mut bad_rt = RiseTransOptions::default();
    bad_rt.max_search_days = -1.0;
    assert!(eng
        .rise_trans(
            BodyId::Planet(PlanetId::Sun),
            obs,
            JD_J2000,
            RiseTransEvent::Rise,
            bad_rt
        )
        .is_err());
}

#[test]
fn test_display_from_str_and_default_traits() {
    assert_eq!("Sun".parse::<BodyId>(), Ok(BodyId::Planet(PlanetId::Sun)));
    assert_eq!(
        "Ceres".parse::<BodyId>(),
        Ok(BodyId::MinorPlanet(MpcNumber(1)))
    );
    assert_eq!(
        "mpc:2060".parse::<BodyId>(),
        Ok(BodyId::MinorPlanet(MpcNumber(2060)))
    );
    assert_eq!(
        "True_Lilith".parse::<BodyId>(),
        Ok(BodyId::LunarPoint(LunarPointId::TrueLilith))
    );
    assert_eq!(
        "Cupido".parse::<BodyId>(),
        Ok(BodyId::Hypothetical(HypotheticalId::Cupido))
    );
    assert_eq!(
        "Spica".parse::<BodyId>(),
        Ok(BodyId::FixedStar(HipId(65474)))
    );
    assert_eq!("P".parse::<HouseSystem>(), Ok(HouseSystem::Placidus));
    assert_eq!("lahiri".parse::<AyanamsaModel>(), Ok(AyanamsaModel::Lahiri));
    assert_eq!(Plane::default(), Plane::Ecliptic);
    assert_eq!(Equinox::default(), Equinox::OfDate);
    assert_eq!(NutationMode::default(), NutationMode::True);
    assert_eq!(AbMode::default(), AbMode::Apparent);
    assert_eq!(Origin::default(), Origin::Geocentric);
    assert_eq!(Zodiac::default(), Zodiac::Tropical);
    assert_eq!(AyanamsaModel::default(), AyanamsaModel::Lahiri);
    assert_eq!(HouseSystem::default(), HouseSystem::Placidus);
}

#[test]
fn test_codex_round2_interpolated_lilith_c2_continuity_and_nutation_rate() {
    let eng = workspace_engine();

    // 1. Verify Interpolated_Lilith has zero moving-window jumps at t_apogee ± 60d
    for probe_jd in [2_451_555.025_848_6, 2_451_582.138_925_7] {
        let dt = 1.0e-6;
        let m_minus = eng
            .calc(
                BodyId::LunarPoint(LunarPointId::InterpolatedLilith),
                TimeInput::from_tt_jd(probe_jd - dt),
                CalcOptions::default(),
            )
            .expect("InterpolatedLilith -dt");
        let m_plus = eng
            .calc(
                BodyId::LunarPoint(LunarPointId::InterpolatedLilith),
                TimeInput::from_tt_jd(probe_jd + dt),
                CalcOptions::default(),
            )
            .expect("InterpolatedLilith +dt");
        let jump_arcsec = ((m_plus.lon_deg - m_minus.lon_deg + 180.0).rem_euclid(360.0) - 180.0)
            .abs()
            * 3600.0;
        assert!(
            jump_arcsec < 0.01,
            "InterpolatedLilith discontinuity at {probe_jd}: {jump_arcsec}\""
        );
        let fd_speed =
            ((m_plus.lon_deg - m_minus.lon_deg + 180.0).rem_euclid(360.0) - 180.0) / (2.0 * dt);
        assert!(
            (m_plus.lon_speed_deg_day - fd_speed).abs() < 1.0e-4,
            "InterpolatedLilith analytical vs FD speed mismatch at {probe_jd}: {} vs {fd_speed}",
            m_plus.lon_speed_deg_day
        );
    }

    // 2. Verify Tropical NutationMode::True includes nutation rate in MeanNode, MeanLilith, InterpolatedLilith
    let jd = 2_461_041.5;
    let dt = 1.0e-4;
    for lp in [
        LunarPointId::MeanNode,
        LunarPointId::MeanSouthNode,
        LunarPointId::MeanLilith,
        LunarPointId::InterpolatedLilith,
    ] {
        let r0 = eng
            .calc(
                BodyId::LunarPoint(lp),
                TimeInput::from_tt_jd(jd),
                CalcOptions::default(),
            )
            .expect("r0");
        let rm = eng
            .calc(
                BodyId::LunarPoint(lp),
                TimeInput::from_tt_jd(jd - dt),
                CalcOptions::default(),
            )
            .expect("rm");
        let rp = eng
            .calc(
                BodyId::LunarPoint(lp),
                TimeInput::from_tt_jd(jd + dt),
                CalcOptions::default(),
            )
            .expect("rp");
        let fd = ((rp.lon_deg - rm.lon_deg + 180.0).rem_euclid(360.0) - 180.0) / (2.0 * dt);
        assert!(
            (r0.lon_speed_deg_day - fd).abs() < 1.0e-6,
            "{lp:?} speed mismatch vs FD: got {}, fd {fd}",
            r0.lon_speed_deg_day
        );
    }
}

#[test]
fn test_codex_round2_saros_full_canon_heliacal_bounds_polich_page_and_swe_compat() {
    use lxephemeris_core::stars::parse_sefstars_catalog;
    use lxephemeris_core::swe_compat::*;
    use std::sync::Arc;

    let eng = Arc::new(workspace_engine());
    let mut swe = SweAdapter::new(Arc::clone(&eng));

    // 1. Full-canon Saros series & member numbers across -2000..+3000 (including Lunar 100 & 165, ancient/future series)
    for (jd_ut, exp_saros, exp_member) in [
        (1_355_906.454_362, 42, 33),
        (1_721_595.0, 66, 43),
        (2_461_265.24, 126, 48),
        (2_635_466.189_022, 173, 2),
    ] {
        let (saros, member) = compute_solar_saros(jd_ut, 0.0);
        assert_eq!(saros, exp_saros);
        assert_eq!(member, exp_member);
    }
    for (jd_ut, exp_saros, exp_member) in [
        (1_355_891.066_554, 16, 66),
        (2_270_276.738_716, 100, 60),
        (2_461_102.98, 133, 27),
        (2_634_772.492_923, 165, 6),
    ] {
        let (saros, member) = compute_lunar_saros(jd_ut, 0.0);
        assert_eq!(saros, exp_saros);
        assert_eq!(member, exp_member);
    }

    // 2. heliacal_phenomena strict max_search_days bound + outer planet event_type 3/4 rejection
    let obs = Observer::new(30.0, 31.2, 0.0).expect("Cairo");
    let mut hel_opts = HeliacalOptions::default();
    hel_opts.max_search_days = 1.0;
    let start_jd = 2_451_545.0;
    let hel_res = eng
        .heliacal_phenomena(BodyId::Planet(PlanetId::Mars), start_jd, obs, hel_opts)
        .expect("heliacal_phenomena");
    for t in [
        hel_res.heliacal_rise_jd_ut1,
        hel_res.heliacal_set_jd_ut1,
        hel_res.evening_rise_jd_ut1,
        hel_res.morning_set_jd_ut1,
    ]
    .into_iter()
    .flatten()
    {
        assert!(
            t >= start_jd && t <= start_jd + hel_opts.max_search_days,
            "heliacal event {t} exceeded [start, start + max_search_days]"
        );
    }
    assert!(swe
        .heliacal_ut(
            start_jd,
            [31.2, 30.0, 0.0],
            [1013.25, 15.0, 50.0, 0.2],
            [36.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            "Mars",
            SE_EVENING_FIRST,
            0
        )
        .is_err());
    assert!(swe
        .heliacal_ut(
            start_jd,
            [31.2, 30.0, 0.0],
            [1013.25, 15.0, 50.0, 0.2],
            [36.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            "Jupiter",
            SE_MORNING_LAST,
            0
        )
        .is_err());

    // 3. Polich-Page ('T') house_pos exactness on cusps and off-ecliptic points + Equator Vertex
    let (t_cusps, _) = swe.houses_armc(120.0, 51.5, 23.439_291_1, 'T').expect("T cusps");
    for (i, &c) in t_cusps.iter().enumerate() {
        let hp = swe
            .house_pos(120.0, 51.5, 23.439_291_1, 'T', [c, 0.0])
            .expect("T house_pos cusp");
        let diff = ((hp - (i as f64 + 1.0) + 6.0).rem_euclid(12.0) - 6.0).abs();
        assert!(diff < 1.0e-7, "Polich-Page cusp {} mapped to {hp}", i + 1);
    }
    let hp_off = swe
        .house_pos(120.0, 51.5, 23.439_291_1, 'T', [45.0, 5.0])
        .expect("T house_pos off-ecliptic");
    assert!((hp_off - 7.996_128_31).abs() < 1.0e-6);

    let (_, eq_ascmc_120) = swe.houses_armc(120.0, 0.0, 23.439_291_1, 'P').expect("eq 120");
    assert_eq!(eq_ascmc_120[3], 0.0);
    let (_, eq_ascmc_240) = swe.houses_armc(240.0, 0.0, 23.439_291_1, 'P').expect("eq 240");
    assert_eq!(eq_ascmc_240[3], 180.0);

    // 4. Custom star phenomena
    let custom_cat =
        "CustomTestStar,alTst,ICRS,10,08,22.311,11,58,01.95,-248.73,5.59,-5.9,41.13,1.35, 16, 999991\n";
    let custom_stars = parse_sefstars_catalog(custom_cat).expect("parse custom_cat");
    let mut eng_custom = workspace_engine();
    eng_custom.add_custom_stars(custom_stars);
    let pheno_custom = eng_custom
        .phenomena(
            BodyId::FixedStar(HipId(999_991)),
            TimeInput::from_tt_jd(JD_J2000),
            CalcOptions::default(),
        )
        .expect("custom star phenomena");
    assert!((pheno_custom.apparent_magnitude_v - 1.35).abs() < 1.0e-9);
    assert!(pheno_custom.elongation_deg > 0.0);

    // 5. SweAdapter::calc SE_MEAN_APOG 3D latitude, SEFLG_SPEED3, SEFLG_TOPOCTR precedence, heliocentric lunar points, houses_ex2 speeds
    let (apog_tt, _) = swe
        .calc(JD_J2000, SE_MEAN_APOG, SEFLG_SWIEPH | SEFLG_SPEED)
        .expect("calc SE_MEAN_APOG");
    let (apog_ut, _) = swe
        .calc_ut(JD_J2000, SE_MEAN_APOG, SEFLG_SWIEPH | SEFLG_SPEED)
        .expect("calc_ut SE_MEAN_APOG");
    assert!(apog_tt[1].abs() > 3.4);
    assert!((apog_tt[1] - apog_ut[1]).abs() < 0.001);

    let (sun_spd3, _) = swe
        .calc_ut(JD_J2000, SE_SUN, SEFLG_SWIEPH | 128)
        .expect("SEFLG_SPEED3");
    assert!(sun_spd3[3].abs() > 0.95);

    let (hel_node, _) = swe
        .calc_ut(JD_J2000, SE_MEAN_NODE, SEFLG_SWIEPH | SEFLG_HELCTR)
        .expect("hel node");
    assert_eq!(hel_node, [0.0; 6]);

    swe.set_topo(0.0, 51.5, 0.0).expect("set_topo");
    let (topo_only, _) = swe
        .calc_ut(JD_J2000, SE_MOON, SEFLG_SWIEPH | SEFLG_TOPOCTR)
        .expect("topo only");
    let (topo_hel, _) = swe
        .calc_ut(
            JD_J2000,
            SE_MOON,
            SEFLG_SWIEPH | SEFLG_TOPOCTR | SEFLG_HELCTR,
        )
        .expect("topo+hel");
    assert!((topo_only[0] - topo_hel[0]).abs() < 1.0e-12);

    let (cusps_rad, ascmc_rad, cspd_rad, aspd_rad) = swe
        .houses_ex2(JD_J2000, 51.5, 0.0, 'P', SEFLG_RADIANS)
        .expect("houses_ex2 radians");
    assert!(cusps_rad[0] <= std::f64::consts::TAU);
    assert!(ascmc_rad[0] <= std::f64::consts::TAU);
    // In pyswisseph, SEFLG_RADIANS converts cusps and ascmc to radians, but keeps speeds in deg/day
    assert!(cspd_rad[0] > 200.0 && cspd_rad[0] < 900.0);
    assert!((aspd_rad[2] - 360.985_647_366).abs() < 1.0e-4);
}

#[test]
fn test_codex_round3_engine_and_compat_regressions() {
    use lxephemeris_core::{
        SweAdapter, SEFLG_EQUATORIAL, SEFLG_J2000, SEFLG_SIDEREAL, SEFLG_SPEED, SEFLG_SWIEPH,
        SE_ECL_TOTAL, SE_EVENING_FIRST, SE_HELIACAL_RISING, SE_HELIACAL_SETTING, SE_MEAN_APOG,
        SE_MOON, SE_MORNING_LAST, SE_NODBIT_MEAN, SE_NODBIT_OSCU, SE_NODBIT_OSCU_BAR, SE_VENUS,
    };
    use std::sync::Arc;

    let eng = Arc::new(workspace_engine());
    let swe = SweAdapter::new(Arc::clone(&eng));

    // 1. Ancient Equation of Time at -2000-01-01 12:00 Gregorian (JD 990575.0): tau in Julian millennia, not centuries
    let jd_2000bce = julday(-2000, 1, 1, 12.0, CalendarKind::Gregorian).expect("julday -2000");
    let eot_sec = swe.time_equ(jd_2000bce).expect("time_equ -2000") * 86_400.0;
    assert!(
        (eot_sec - (-722.35)).abs() < 35.0,
        "Ancient equation of time at -2000 CE expected ~ -722s, got {eot_sec}s"
    );

    // 2. Short real polar daylight (12.97 min at lat = 67.39N, jd = 2461395.5) must not fall back to 12h polar day
    let obs_polar = Observer::new(67.39, 0.0, 0.0).expect("polar observer");
    let panch = eng
        .panchanga(
            TimeInput::from_ut1_jd(2_461_395.5),
            obs_polar,
            AyanamsaModel::Lahiri,
        )
        .expect("panchanga short polar daylight");
    assert!(
        !panch.polar_day_fallback,
        "Short polar daylight (12.97 min) must not trigger polar fallback"
    );
    let day_len_min = (panch.sunset_jd_ut1 - panch.sunrise_jd_ut1) * 1_440.0;
    assert!(
        (10.0..20.0).contains(&day_len_min),
        "Expected ~12.97 min daylight, got {day_len_min} min"
    );

    // 3. Apparent fixed-star stations (Regulus HIP 49669: 4 mean annual+EMB-reflex aberration stations, 14 true stations with fortnightly nutation)
    let mean_star_opts = CalcOptions {
        nutation: NutationMode::Mean,
        ..CalcOptions::default()
    };
    let star_stations_mean = eng
        .search_stations_in_range(
            BodyId::FixedStar(HipId(49669)),
            TimeInput::from_ut1_jd(JD_J2000),
            TimeInput::from_ut1_jd(JD_J2000 + 366.0),
            &mean_star_opts,
        )
        .expect("fixed star mean stations");
    assert_eq!(star_stations_mean.len(), 4);
    let star_stations_true = eng
        .search_stations_in_range(
            BodyId::FixedStar(HipId(49669)),
            TimeInput::from_ut1_jd(JD_J2000),
            TimeInput::from_ut1_jd(JD_J2000 + 366.0),
            &CalcOptions::default(),
        )
        .expect("fixed star true stations");
    assert_eq!(star_stations_true.len(), 14);

    // 4. Topocentric lunar declination extrema over 15 days must resolve diurnal parallax oscillations near lunistice
    let topo_opts = CalcOptions {
        origin: Origin::Topocentric(Observer::new(51.5, 0.0, 0.0).expect("Greenwich")),
        ..CalcOptions::default()
    };
    let topo_dec_ext = eng
        .search_declination_extrema_in_range(
            BodyId::Planet(PlanetId::Moon),
            TimeInput::from_ut1_jd(JD_J2000),
            TimeInput::from_ut1_jd(JD_J2000 + 15.0),
            &topo_opts,
        )
        .expect("topocentric lunar dec extrema");
    assert!(
        topo_dec_ext.len() >= 2,
        "Expected diurnal topocentric declination extrema near lunistice, got {}",
        topo_dec_ext.len()
    );

    // 5. Moon heliacal crescent visibility: event_type 3 & 4 converge, event_type 1 & 2 reject
    let gp = [31.2, 30.0, 0.0];
    let atm = [1013.25, 15.0, 50.0, 0.2];
    let obs_arr = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    assert!(swe
        .heliacal_ut(JD_J2000, gp, atm, obs_arr, "Moon", SE_HELIACAL_RISING, 0)
        .is_err());
    assert!(swe
        .heliacal_ut(JD_J2000, gp, atm, obs_arr, "Moon", SE_HELIACAL_SETTING, 0)
        .is_err());
    let (moon_ev3, _) = swe
        .heliacal_ut(JD_J2000, gp, atm, obs_arr, "Moon", SE_EVENING_FIRST, 0)
        .expect("Moon evening first crescent");
    let (moon_ev4, _) = swe
        .heliacal_ut(JD_J2000, gp, atm, obs_arr, "Moon", SE_MORNING_LAST, 0)
        .expect("Moon morning last crescent");
    assert!((moon_ev4[1] - 2_451_548.65).abs() < 1.2);
    assert!((moon_ev3[1] - 2_451_552.17).abs() < 1.2);

    // 6. nod_aps default method=0 matches SE_NODBIT_MEAN, and SE_NODBIT_OSCU_BAR matches SE_NODBIT_OSCU for Venus
    let (n_def, _, p_def, a_def, _) = swe
        .nod_aps(JD_J2000, SE_MOON, SEFLG_SWIEPH | SEFLG_SPEED, 0)
        .expect("nod_aps default 0");
    let (n_mean, _, p_mean, a_mean, _) = swe
        .nod_aps(JD_J2000, SE_MOON, SEFLG_SWIEPH | SEFLG_SPEED, SE_NODBIT_MEAN)
        .expect("nod_aps mean 1");
    assert!((n_def[0] - n_mean[0]).abs() < 1.0e-12);
    assert!((p_def[0] - p_mean[0]).abs() < 1.0e-12);
    assert!((a_def[0] - a_mean[0]).abs() < 1.0e-12);

    let (_, _, _, v_oscu, _) = swe
        .nod_aps(JD_J2000, SE_VENUS, SEFLG_SWIEPH, SE_NODBIT_OSCU)
        .expect("Venus oscu");
    let (_, _, _, v_bary, _) = swe
        .nod_aps(JD_J2000, SE_VENUS, SEFLG_SWIEPH, SE_NODBIT_OSCU_BAR)
        .expect("Venus oscu_bar");
    assert!((v_oscu[0] - v_bary[0]).abs() < 1.0e-12);

    // 7. Global lunar occultation tret layout: [max, noon, t1, t4, t2, t3, t_cent1, t_cent2, 0, 0]
    let (occ_rf, occ_tret) = swe
        .lun_occult_when_glob_ex(JD_J2000, SE_VENUS, None, SEFLG_SWIEPH, SE_ECL_TOTAL, false)
        .expect("Venus global occultation");
    assert_ne!(occ_rf & SE_ECL_TOTAL, 0);
    assert!(occ_tret[2] < occ_tret[4]);
    assert!(occ_tret[4] < occ_tret[0]);
    assert!(occ_tret[0] < occ_tret[5]);
    assert!(occ_tret[5] < occ_tret[3]);

    // 8. lun_eclipse_how horizon check & sol_eclipse_when_loc visibility flags/azimuth
    let (lhow_rf, _) = swe
        .lun_eclipse_how(2_460_394.800_615_769, [0.0, 51.5, 0.0], SEFLG_SWIEPH)
        .expect("lun_eclipse_how Greenwich below horizon");
    assert_eq!(lhow_rf, 0);
    let (sloc_rf, _, sloc_attr) = swe
        .sol_eclipse_when_loc(2_460_310.5, [-96.8, 32.8, 0.0], SEFLG_SWIEPH, false)
        .expect("Dallas 2024 total solar eclipse");
    assert_eq!(sloc_rf, 8068);
    assert!((sloc_attr[4] - 8.07).abs() < 0.5);

    // 9. SE_MEAN_APOG sidereal speed and J2000 equatorial speed
    let (apog_trop, _) = swe
        .calc_ut(JD_J2000, SE_MEAN_APOG, SEFLG_SWIEPH | SEFLG_SPEED)
        .expect("apog trop");
    let (apog_sid, _) = swe
        .calc_ut(
            JD_J2000,
            SE_MEAN_APOG,
            SEFLG_SWIEPH | SEFLG_SIDEREAL | SEFLG_SPEED,
        )
        .expect("apog sid");
    assert!(
        (apog_trop[3] - apog_sid[3]).abs() > 1.0e-5,
        "Sidereal SE_MEAN_APOG speed must include ayanamsa rate"
    );
    let (apog_eq_j2k, _) = swe
        .calc_ut(
            JD_J2000 - 36_525.0,
            SE_MEAN_APOG,
            SEFLG_SWIEPH | SEFLG_J2000 | SEFLG_EQUATORIAL | SEFLG_SPEED,
        )
        .expect("apog eq j2k");
    assert!(apog_eq_j2k[3].abs() > 0.05);

    // 10. Polich-Page ('T') cusps at -65°, 'W'/'U' cusp speeds, and difdeg2n(180, 0) == -180
    let (t_cusps_65, _, _, _) = swe
        .houses_ex2(2_461_041.5, -65.0, 13.0, 'T', 0)
        .expect("T houses at -65");
    assert!((t_cusps_65[5] - 342.1706).abs() < 0.01);
    let (_, w_ascmc, w_cspd, w_aspd) = swe
        .houses_ex2(JD_J2000, 51.5, 0.0, 'W', 0)
        .expect("W houses_ex2");
    let _ = w_ascmc;
    assert!((w_cspd[0] - w_aspd[0]).abs() < 1.0e-12);
    assert!((w_cspd[3] - w_aspd[1]).abs() < 1.0e-12);
    assert_eq!(w_cspd[1], 0.0);
    let (_, _, u_cspd, u_aspd) = swe
        .houses_ex2(JD_J2000, 51.5, 0.0, 'U', 0)
        .expect("U houses_ex2");
    assert!((u_cspd[0] - u_aspd[0]).abs() < 1.0e-12);
    assert!((u_cspd[9] - u_aspd[1]).abs() < 1.0e-12);
    assert_eq!(u_cspd[1], 0.0);
    assert_eq!(SweAdapter::difdeg2n(180.0, 0.0), -180.0);
}

#[test]
fn test_codex_round4_astra_regressions() {
    use lxephemeris_core::{
        AspectPreset, HouseOptions, SweAdapter, SEFLG_SIDEREAL, SEFLG_SPEED, SEFLG_SWIEPH,
        SE_BIT_ASTRO_TWILIGHT, SE_BIT_CIVIL_TWILIGHT, SE_BIT_NAUTIC_TWILIGHT, SE_CALC_RISE,
        SE_EARTH, SE_ECL_ANNULAR_TOTAL, SE_EVENING_FIRST, SE_INTP_PERG, SE_JUPITER, SE_MOON,
        SE_NODBIT_MEAN, SE_NODBIT_OSCU, SE_NODBIT_OSCU_BAR, SE_SIDBIT_ECL_T0, SE_SIDM_LAHIRI,
        SE_SIDM_LAHIRI_1940, SE_SIDM_LAHIRI_ICRC, SE_SIDM_LAHIRI_VP285, SE_SIDM_TRUE_CITRA, SE_SUN,
    };
    use std::sync::Arc;

    let eng = Arc::new(workspace_engine());
    let mut swe = SweAdapter::new(Arc::clone(&eng));

    // 1. SE_INTP_PERG = 22 (Interpolated Perigee)
    assert_eq!(SE_INTP_PERG, 22);
    assert_eq!(
        SweAdapter::get_planet_name(SE_INTP_PERG),
        "Interpolated_Perigee"
    );
    let (perg_xx, _) = swe
        .calc_ut(JD_J2000, SE_INTP_PERG, SEFLG_SWIEPH | SEFLG_SPEED)
        .expect("calc_ut SE_INTP_PERG");
    assert!(
        (90.0..94.0).contains(&perg_xx[0]),
        "Expected J2000 perigee lon ~91.9 deg, got {}",
        perg_xx[0]
    );
    assert!(
        (-3.5..-2.0).contains(&perg_xx[1]),
        "Expected J2000 perigee lat ~-2.8 deg, got {}",
        perg_xx[1]
    );
    assert!(
        (0.002_35..0.002_45).contains(&perg_xx[2]),
        "Expected J2000 perigee dist ~0.002388 AU, got {}",
        perg_xx[2]
    );

    // 2. Extended SE_SIDM_* 43..=46, SE_SIDBIT_* mask, and True Citra (27) Spica anchor at 180.0°
    for sid_id in [
        SE_SIDM_LAHIRI_1940,
        SE_SIDM_LAHIRI_VP285,
        45,
        SE_SIDM_LAHIRI_ICRC,
        SE_SIDM_LAHIRI | SE_SIDBIT_ECL_T0,
    ] {
        swe.set_sid_mode(sid_id).expect("set_sid_mode extended");
        let (_, aya) = swe
            .get_ayanamsa_ex_ut(JD_J2000, SEFLG_SWIEPH)
            .expect("get_ayanamsa_ex_ut");
        assert!((23.0..25.0).contains(&aya));
    }
    swe.set_sid_mode(SE_SIDM_TRUE_CITRA)
        .expect("set_sid_mode TRUE_CITRA");
    let (spica_xx, _, _) = swe
        .fixstar2_ut("Spica", JD_J2000, SEFLG_SWIEPH | SEFLG_SIDEREAL)
        .expect("Spica sidereal TRUE_CITRA");
    assert!(
        (spica_xx[0] - 180.0).abs() < 1.0e-4,
        "Spica in TRUE_CITRA must be anchored at 180.0 deg, got {}",
        spica_xx[0]
    );
    // Numeric star catalog lookup "1".."53"
    let (star1_xx, star1_name, _) = swe
        .fixstar2_ut("1", JD_J2000, SEFLG_SWIEPH)
        .expect("fixstar2_ut numeric '1'");
    assert!(!star1_name.is_empty());
    assert!((0.0..360.0).contains(&star1_xx[0]));
    swe.set_sid_mode(SE_SIDM_LAHIRI).expect("restore Lahiri");

    // 3. Geocentric SE_EARTH returns [0; 6]; nod_aps supports SE_SUN, SE_EARTH, and SE_NODBIT_OSCU_BAR
    let (earth_geo, _) = swe
        .calc_ut(JD_J2000, SE_EARTH, SEFLG_SWIEPH | SEFLG_SPEED)
        .expect("geocentric SE_EARTH");
    assert_eq!(earth_geo, [0.0; 6]);
    let (_, _, sun_peri, sun_aphe, _) = swe
        .nod_aps_ut(JD_J2000, SE_SUN, SEFLG_SWIEPH, SE_NODBIT_MEAN)
        .expect("nod_aps SE_SUN");
    let (_, _, earth_peri, earth_aphe, _) = swe
        .nod_aps_ut(JD_J2000, SE_EARTH, SEFLG_SWIEPH, SE_NODBIT_MEAN)
        .expect("nod_aps SE_EARTH geocentric");
    assert!((sun_peri[0] - 282.9124).abs() < 0.05);
    assert!((sun_aphe[0] - 102.9124).abs() < 0.05);
    assert!((earth_peri[0] - 191.6875).abs() < 0.05);
    assert!((earth_aphe[0] - 281.6590).abs() < 0.05);
    let (_, _, earth_hel_peri, earth_hel_aphe, _) = swe
        .nod_aps_ut(
            JD_J2000,
            SE_EARTH,
            SEFLG_SWIEPH | lxephemeris_core::SEFLG_HELCTR,
            SE_NODBIT_MEAN,
        )
        .expect("nod_aps SE_EARTH heliocentric");
    assert!(((sun_peri[0] - earth_hel_peri[0]).rem_euclid(360.0) - 180.0).abs() < 0.01);
    assert!(((sun_aphe[0] - earth_hel_aphe[0]).rem_euclid(360.0) - 180.0).abs() < 0.01);
    let (_, _, jup_oscu, _, _) = swe
        .nod_aps_ut(JD_J2000, SE_JUPITER, SEFLG_SWIEPH, SE_NODBIT_OSCU)
        .expect("Jupiter oscu");
    let (_, _, jup_bar, _, _) = swe
        .nod_aps_ut(JD_J2000, SE_JUPITER, SEFLG_SWIEPH, SE_NODBIT_OSCU_BAR)
        .expect("Jupiter oscu_bar");
    assert!((jup_oscu[0] - jup_bar[0]).abs() < 1.0e-9);
    let (_, _, moon_mean_combo, _, _) = swe
        .nod_aps_ut(
            JD_J2000,
            SE_MOON,
            SEFLG_SWIEPH,
            SE_NODBIT_MEAN | SE_NODBIT_OSCU,
        )
        .expect("Moon mean|oscu precedence");
    let (_, _, moon_mean_only, _, _) = swe
        .nod_aps_ut(JD_J2000, SE_MOON, SEFLG_SWIEPH, SE_NODBIT_MEAN)
        .expect("Moon mean only");
    assert!((moon_mean_combo[0] - moon_mean_only[0]).abs() < 1.0e-12);

    // 4. House systems: Azimuthal ('H') southern/equator, Sunshine ('I') polar, Whole Sign ('W') sidereal speeds, SoftFail polar
    let (h_south_cusps, h_south_ascmc) = swe
        .houses_ex(JD_J2000, -33.8688, 151.2093, 'H', 0)
        .expect("Azimuthal Sydney");
    assert!((h_south_cusps[0] - 167.8703).abs() < 0.01);
    assert!((h_south_cusps[9] - h_south_ascmc[1]).abs() < 1.0e-9);
    let (h_eq_cusps, h_eq_ascmc) = swe
        .houses_ex(JD_J2000, 0.0, 0.0, 'H', 0)
        .expect("Azimuthal equator");
    assert_eq!(h_eq_ascmc[6], 0.0); // co_asc_munkasey_deg at equator == 0
    assert!((h_eq_cusps[0] - 180.0).abs() < 1.0e-9);
    let (i_polar_cusps, _) = swe
        .houses_ex(2_460_482.0, 72.0, 25.0, 'I', 0)
        .expect("Sunshine polar day");
    assert!((i_polar_cusps[1] - i_polar_cusps[3]).abs() < 1.0e-6);
    assert!((i_polar_cusps[1] - 292.3989).abs() < 0.01);
    let w_sid_native = eng
        .houses(
            TimeInput::from_ut1_jd(JD_J2000),
            Observer::new(51.5, 0.0, 0.0).expect("Greenwich"),
            HouseOptions {
                system: HouseSystem::WholeSign,
                zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
                include_speed: true,
                ..HouseOptions::default()
            },
        )
        .expect("Whole Sign sidereal native speeds");
    assert_eq!(w_sid_native.cusp_speeds_deg_day, [0.0; 12]);
    let u_native = eng
        .houses(
            TimeInput::from_ut1_jd(JD_J2000),
            Observer::new(51.5, 0.0, 0.0).expect("Greenwich"),
            HouseOptions {
                system: HouseSystem::KrusinskiPisaGoelzer,
                include_speed: true,
                ..HouseOptions::default()
            },
        )
        .expect("Krusinski-Pisa-Goelzer native speeds");
    assert!(u_native.cusp_speeds_deg_day[1] > 100.0);

    // 5. chart_snapshot house_pos consistency in sidereal Placidus and non-3D house systems ('W', 'E', 'V')
    let london = Observer::new(51.5074, -0.1278, 0.0).expect("London");
    for sys in [
        HouseSystem::Placidus,
        HouseSystem::WholeSign,
        HouseSystem::EqualFromAsc,
        HouseSystem::VehlowEqual,
    ] {
        let h_opts = HouseOptions {
            system: sys,
            zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
            ..HouseOptions::default()
        };
        let c_opts = CalcOptions {
            zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
            ..CalcOptions::default()
        };
        let snap = eng
            .chart_snapshot(
                TimeInput::from_ut1_jd(JD_J2000),
                london,
                &[BodyId::Planet(PlanetId::Sun), BodyId::Planet(PlanetId::Moon)],
                h_opts,
                c_opts,
                AspectPreset::Major,
                1.0,
            )
            .expect("chart_snapshot sidereal");
        if !matches!(sys, HouseSystem::Placidus) {
            for b in &snap.bodies {
                assert_eq!(b.house_pos_3d.floor() as u8, b.house_number);
            }
        }
    }

    // 6. Heliacal setting: 5h before last visible setting (2451747.0 -> ~2451747.22) vs after last visible setting (2451748.0 -> next year > 2452100)
    let (reg_hel_before, _) = swe
        .heliacal_ut(
            2_451_747.0,
            [31.2, 30.0, 0.0],
            [1013.25, 15.0, 50.0, 0.2],
            [36.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            "Regulus",
            2,
            0,
        )
        .expect("Regulus heliacal setting from 2451747.0");
    assert!(
        (2_451_747.0..2_451_748.0).contains(&reg_hel_before[0]),
        "Expected heliacal setting on 2451747.22, got {}",
        reg_hel_before[0]
    );
    let (reg_hel_after, _) = swe
        .heliacal_ut(
            2_451_748.0,
            [31.2, 30.0, 0.0],
            [1013.25, 15.0, 50.0, 0.2],
            [36.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            "Regulus",
            2,
            0,
        )
        .expect("Regulus heliacal setting from 2451748.0");
    assert!(
        reg_hel_after[0] > 2_452_100.0,
        "Expected next year's heliacal setting > 2452100, got {}",
        reg_hel_after[0]
    );

    // 7. 1927-01-03 Hybrid Solar Eclipse & twilight rise_trans flags
    let (hyb_rf, hyb_tret) = swe
        .sol_eclipse_when_glob(2_424_880.0, SEFLG_SWIEPH, 0, false)
        .expect("1927 hybrid solar eclipse");
    assert_ne!(hyb_rf & SE_ECL_ANNULAR_TOTAL, 0);
    assert!(hyb_tret[1] > 0.0 && hyb_tret[6] > 0.0 && hyb_tret[7] > 0.0);

    let r_std = swe
        .rise_trans(
            JD_J2000,
            SE_SUN,
            None,
            SEFLG_SWIEPH,
            SE_CALC_RISE,
            [0.0, 51.5, 0.0],
            1013.25,
            15.0,
        )
        .expect("std rise")
        .expect("some std rise");
    let r_civ = swe
        .rise_trans(
            JD_J2000,
            SE_SUN,
            None,
            SEFLG_SWIEPH,
            SE_CALC_RISE | SE_BIT_CIVIL_TWILIGHT,
            [0.0, 51.5, 0.0],
            1013.25,
            15.0,
        )
        .expect("civil dawn")
        .expect("some civil dawn");
    let r_nau = swe
        .rise_trans(
            JD_J2000,
            SE_SUN,
            None,
            SEFLG_SWIEPH,
            SE_CALC_RISE | SE_BIT_NAUTIC_TWILIGHT,
            [0.0, 51.5, 0.0],
            1013.25,
            15.0,
        )
        .expect("nautic dawn")
        .expect("some nautic dawn");
    let r_ast = swe
        .rise_trans(
            JD_J2000,
            SE_SUN,
            None,
            SEFLG_SWIEPH,
            SE_CALC_RISE | SE_BIT_ASTRO_TWILIGHT,
            [0.0, 51.5, 0.0],
            1013.25,
            15.0,
        )
        .expect("astro dawn")
        .expect("some astro dawn");
    assert!(r_ast < r_nau && r_nau < r_civ && r_civ < r_std);
    let _ = SE_EVENING_FIRST;
}
