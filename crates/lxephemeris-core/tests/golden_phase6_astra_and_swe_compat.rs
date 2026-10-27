//! Phase 6 Verification Suite: Astra Requirements Closure (`R1.1..R11.4`) and
//! Generalized Swiss Ephemeris Migration Layer (`SweAdapter`).

#![cfg(feature = "swe-compat")]
#![allow(clippy::indexing_slicing)]

use lxephemeris_core::{
    house_of_lon, is_cusp_ring_monotonic, AtmosphereParams, BodyId, CalcOptions, DiscMode,
    EclipseOptions, Engine, EngineConfig, HeliacalOptions, HipId, HouseOptions, HouseSystem,
    HouseWarning, NutationMode, Observer, PlanetId, PolarPolicy, RiseTransEvent, RiseTransOptions,
    SolarEclipseKind, SweAdapter, TimeInput, Zodiac, AU_KM, JD_J2000, SEFLG_EQUATORIAL,
    SEFLG_HELCTR, SEFLG_J2000, SEFLG_MOSEPH, SEFLG_NONUT, SEFLG_RADIANS, SEFLG_SIDEREAL,
    SEFLG_SPEED, SEFLG_SWIEPH, SEFLG_TOPOCTR, SEFLG_TRUEPOS, SE_AST_OFFSET, SE_BIT_DISC_BOTTOM,
    SE_BIT_DISC_CENTER, SE_BIT_NO_REFRACTION, SE_CALC_RISE, SE_CERES, SE_CHIRON, SE_CUPIDO,
    SE_ECL2HOR, SE_ECL_NUT, SE_ECL_TOTAL, SE_EQU2HOR, SE_HELFLAG_HIGH_PRECISION,
    SE_HELIACAL_RISING, SE_MEAN_APOG, SE_MEAN_NODE, SE_MOON, SE_OSCU_APOG, SE_SIDM_KRISHNAMURTI,
    SE_SIDM_LAHIRI, SE_SUN, SE_TRUE_NODE,
};
use std::path::Path;
use std::sync::Arc;

fn load_engine() -> Arc<Engine> {
    let data_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .unwrap_or_else(|| Path::new("."))
        .join("data");
    Arc::new(
        Engine::from_data_dir(
            data_dir,
            EngineConfig {
                max_threads: Some(4),
                pin_time_table_version: Some("embedded".to_owned()),
                verify_manifest_sha256: false,
            },
        )
        .expect("load Engine"),
    )
}

#[test]
fn test_phase6_r1_delta_t_and_sidtime() {
    let engine = load_engine();
    let t_j2000_ut = TimeInput::from_ut1_jd(JD_J2000);

    // R1.2: delta_t at J2000.0
    let dt = engine.delta_t(t_j2000_ut);
    assert!(
        (dt.delta_t_sec - 63.8285).abs() < 0.05,
        "Unexpected delta_t at J2000: {}",
        dt.delta_t_sec
    );
    assert!(
        dt.tt_minus_tdb_sec.abs() < 0.002,
        "|TT - TDB| must be bounded by ~1.7 ms, got {}",
        dt.tt_minus_tdb_sec
    );
    assert!(!dt.time_extrapolated);

    // R1.4: sidtime at J2000.0 (Greenwich and local at +15 deg East = +1h)
    let st_greenwich = engine.sidtime(t_j2000_ut, NutationMode::True, 0.0);
    let st_15e = engine.sidtime(t_j2000_ut, NutationMode::True, 15.0);

    assert!(
        (st_greenwich.gmst_hours - 18.697_374_558).abs() < 1e-4,
        "GMST at J2000 should be ~18.69737h, got {}",
        st_greenwich.gmst_hours
    );
    let diff_lst = (st_15e.local_sidereal_hours - st_greenwich.gast_hours).rem_euclid(24.0);
    assert!((diff_lst - 1.0).abs() < 1e-12);
    assert!((st_15e.armc_deg - st_15e.local_sidereal_hours * 15.0).abs() < 1e-12);
}

#[test]
fn test_phase6_r2_ecl_nut_and_dist_km() {
    let engine = load_engine();
    let t_j2000_tt = TimeInput::from_tt_jd(JD_J2000);

    // R2.2: ecl_nut via Engine::ecl_nut and BodyId::EclNut
    let en = engine.ecl_nut(t_j2000_tt).expect("ecl_nut");
    assert!((en.mean_obliquity_deg - 23.439_291_1).abs() < 1e-4);
    assert!((en.true_obliquity_deg - (en.mean_obliquity_deg + en.nutation_obl_deg)).abs() < 1e-12);

    let en_body = engine
        .calc(BodyId::EclNut, t_j2000_tt, CalcOptions::default())
        .expect("calc BodyId::EclNut");
    assert_eq!(en_body.target_id, -1);
    assert!((en_body.lon_deg - en.true_obliquity_deg).abs() < 1e-12);
    assert!((en_body.lat_deg - en.mean_obliquity_deg).abs() < 1e-12);
    assert!((en_body.dist_au - en.nutation_lon_deg).abs() < 1e-12);
    assert!((en_body.lon_speed_deg_day - en.nutation_obl_deg).abs() < 1e-12);

    // R2.5: dist_km and dist_speed_km_day on CalcResult
    let moon = engine
        .calc(
            BodyId::Planet(PlanetId::Moon),
            t_j2000_tt,
            CalcOptions::default(),
        )
        .expect("calc Moon");
    assert!((moon.dist_km - moon.dist_au * AU_KM).abs() < 1e-6);
    assert!((moon.dist_speed_km_day - moon.dist_speed_au_day * AU_KM).abs() < 1e-6);
    assert!(
        moon.dist_km > 350_000.0 && moon.dist_km < 410_000.0,
        "Moon dist_km out of range: {}",
        moon.dist_km
    );
}

#[test]
fn test_phase6_r2_6_and_r7_3_azalt_ecliptic_and_near_horizon() {
    let engine = load_engine();
    let t_j2000_tt = TimeInput::from_tt_jd(JD_J2000);
    let obs = Observer::new(51.5074, -0.1278, 0.0).expect("London");
    let atmo = AtmosphereParams::default();

    // Compare calc_azalt(Sun) against azalt_ecliptic(sun.lon_deg, sun.lat_deg)
    let topocentric_ecl = CalcOptions {
        origin: lxephemeris_core::Origin::Topocentric(obs),
        plane: lxephemeris_core::Plane::Ecliptic,
        ..CalcOptions::default()
    };
    let sun_ecl = engine
        .calc(BodyId::Planet(PlanetId::Sun), t_j2000_tt, topocentric_ecl)
        .expect("sun topocentric ecl");

    let hz_direct = engine
        .calc_azalt(BodyId::Planet(PlanetId::Sun), t_j2000_tt, obs, atmo)
        .expect("calc_azalt Sun");
    let hz_from_ecl =
        engine.azalt_ecliptic(obs, t_j2000_tt, sun_ecl.lon_deg, sun_ecl.lat_deg, atmo);

    assert!((hz_direct.azimuth_deg - hz_from_ecl.azimuth_deg).abs() < 1e-9);
    assert!((hz_direct.altitude_true_deg - hz_from_ecl.altitude_true_deg).abs() < 1e-9);
    assert!((hz_direct.altitude_apparent_deg - hz_from_ecl.altitude_apparent_deg).abs() < 1e-9);
    assert_eq!(
        hz_direct.near_horizon,
        hz_direct.altitude_true_deg.abs() < 5.0
    );
}

#[test]
fn test_phase6_r4_4_and_r4_5_polar_fallback_whole_sign_and_house_of_lon() {
    let engine = load_engine();
    let t_j2000_tt = TimeInput::from_tt_jd(JD_J2000);
    // Tromso, Norway (inside Arctic Circle: 69.6492 N)
    let tromso = Observer::new(69.6492, 18.9553, 0.0).expect("Tromso");

    // 1. Placidus ('P', HardFail) with PolarPolicy::FallbackToWholeSign
    let h_placidus_ws = engine
        .houses(
            t_j2000_tt,
            tromso,
            HouseOptions {
                system: HouseSystem::Placidus,
                polar_policy: PolarPolicy::FallbackToWholeSign,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::True,
                include_speed: true,
            },
        )
        .expect("Placidus with FallbackToWholeSign at Tromso");

    assert_eq!(h_placidus_ws.requested_house_system, HouseSystem::Placidus);
    assert_eq!(h_placidus_ws.effective_house_system, HouseSystem::WholeSign);
    assert!(h_placidus_ws.polar_house_fallback);
    assert!(h_placidus_ws
        .warnings
        .contains(&HouseWarning::PolarFallbackToWholeSign));
    assert!(is_cusp_ring_monotonic(&h_placidus_ws.cusps_deg));

    // Verify Whole Sign cusps match true Ascendant sign (no latitude clamping!)
    let expected_c1 = (h_placidus_ws.angles.asc_deg / 30.0).floor() * 30.0;
    assert!((h_placidus_ws.cusps_deg[0] - expected_c1).abs() < 1e-12);

    // 2. Campanus ('C', SoftFail) at high latitude: check monotonic cusp ring or fallback
    let svalbard = Observer::new(78.2232, 15.6267, 0.0).expect("Svalbard");
    for step in 0..24 {
        let t = TimeInput::from_tt_jd(JD_J2000 + f64::from(step) / 24.0);
        let h_camp_ws = engine
            .houses(
                t,
                svalbard,
                HouseOptions {
                    system: HouseSystem::Campanus,
                    polar_policy: PolarPolicy::FallbackToWholeSign,
                    zodiac: Zodiac::Tropical,
                    nutation: NutationMode::True,
                    include_speed: false,
                },
            )
            .expect("Campanus at Svalbard with FallbackToWholeSign");
        assert!(
            is_cusp_ring_monotonic(&h_camp_ws.cusps_deg),
            "Cusp ring must always be strictly monotonic under FallbackToWholeSign"
        );
    }

    // 3. R4.5: house_of_lon with eps = 1e-9 deg and wrap-around across 360 -> 0 deg
    let cusps = [
        350.0, 20.0, 50.0, 80.0, 110.0, 140.0, 170.0, 200.0, 230.0, 260.0, 290.0, 320.0,
    ];
    // House 1 spans 350..20 across 0 deg
    assert_eq!(house_of_lon(355.0, &cusps, 1e-9), 1);
    assert_eq!(house_of_lon(0.0, &cusps, 1e-9), 1);
    assert_eq!(house_of_lon(19.999, &cusps, 1e-9), 1);
    // Boundary epsilon: 19.999_999_999_5 is within 1e-9 of cusp 2 (20.0), so snaps to House 2
    assert_eq!(house_of_lon(20.0 - 5e-10, &cusps, 1e-9), 2);
    assert_eq!(house_of_lon(20.0, &cusps, 1e-9), 2);
    assert_eq!(house_of_lon(349.999_999_999_5, &cusps, 1e-9), 1);
}

#[test]
fn test_phase6_r6_1_r6_2_r7_1_stars_vmag_and_rise_trans() {
    let engine = load_engine();

    // R6.1: star_info and list_stars
    let all_stars = engine.list_stars();
    assert_eq!(all_stars.len(), 52);

    let sirius = engine.star_info("Sirius").expect("Sirius metadata");
    assert_eq!(sirius.hip, 32349);
    assert!((sirius.vmag - (-1.46)).abs() < 0.05);

    let regulus = engine.star_info("Regulus").expect("Regulus metadata");
    assert_eq!(regulus.hip, 49669);
    assert!((regulus.vmag - 1.36).abs() < 0.05);

    // CalcResult for FixedStar also populates magnitude_v
    let regulus_pos = engine
        .calc_star_by_name(
            "Regulus",
            TimeInput::from_tt_jd(JD_J2000),
            CalcOptions::default(),
        )
        .expect("calc Regulus");
    assert_eq!(regulus_pos.magnitude_v, Some(regulus.vmag));

    // R6.2 & R7.1: rise_trans for FixedStar (Sirius) with proper motion and DiscMode options
    let greenwich = Observer::new(51.4769, 0.0, 0.0).expect("Greenwich");
    let sirius_rise = engine
        .rise_trans(
            BodyId::FixedStar(HipId(32349)),
            greenwich,
            JD_J2000,
            RiseTransEvent::Rise,
            RiseTransOptions::default(),
        )
        .expect("Sirius rise_trans")
        .expect("Sirius rises at Greenwich");
    assert!(sirius_rise > JD_J2000 && sirius_rise < JD_J2000 + 1.0);

    // DiscMode::Center vs DiscMode::UpperLimb for the Sun
    let sun_rise_upper = engine
        .rise_trans(
            BodyId::Planet(PlanetId::Sun),
            greenwich,
            JD_J2000,
            RiseTransEvent::Rise,
            RiseTransOptions {
                disc_mode: DiscMode::UpperLimb,
                ..RiseTransOptions::default()
            },
        )
        .expect("Sun rise upper")
        .expect("Sun rises");
    let sun_rise_center = engine
        .rise_trans(
            BodyId::Planet(PlanetId::Sun),
            greenwich,
            JD_J2000,
            RiseTransEvent::Rise,
            RiseTransOptions {
                disc_mode: DiscMode::Center,
                ..RiseTransOptions::default()
            },
        )
        .expect("Sun rise center")
        .expect("Sun rises");
    // Upper limb rises ~2-3 minutes earlier than disc center
    let diff_min = (sun_rise_center - sun_rise_upper) * 1440.0;
    assert!(
        diff_min > 1.0 && diff_min < 5.0,
        "Expected disc center to rise 1..5 min after upper limb, got {diff_min:.2} min"
    );
}

#[test]
fn test_phase6_r6_3_heliacal_phenomena_schaefer() {
    let engine = load_engine();
    // Cairo / Memphis (~30.0 N, 31.2 E) where Sirius heliacal rising occurs in mid-to-late July / early August
    let cairo = Observer::new(30.0444, 31.2357, 20.0).expect("Cairo");
    let jd_2024_jan1 = 2_460_310.5; // 2024-01-01 00:00 UT

    let sirius_hel = engine
        .heliacal_phenomena(
            BodyId::FixedStar(HipId(32349)),
            jd_2024_jan1,
            cairo,
            HeliacalOptions::default(),
        )
        .expect("Sirius heliacal phenomena at Cairo");

    assert_eq!(sirius_hel.target_name, "Sirius");
    assert!((sirius_hel.visual_magnitude - (-1.46)).abs() < 0.05);
    assert!(sirius_hel.arcus_visionis_deg > 5.0 && sirius_hel.arcus_visionis_deg < 12.0);
    let hel_rise = sirius_hel
        .heliacal_rise_jd_ut1
        .expect("Sirius heliacal rise found in 2024");
    let hel_set = sirius_hel
        .heliacal_set_jd_ut1
        .expect("Sirius heliacal set found in 2024");
    assert!(hel_rise > jd_2024_jan1 && hel_rise < jd_2024_jan1 + 366.0);
    assert!(hel_set > jd_2024_jan1 && hel_set < jd_2024_jan1 + 366.0);

    // Also test a planet (Venus or Mars)
    let mars_hel = engine
        .heliacal_phenomena(
            BodyId::Planet(PlanetId::Mars),
            jd_2024_jan1,
            cairo,
            HeliacalOptions {
                max_search_days: 800.0,
                ..HeliacalOptions::default()
            },
        )
        .expect("Mars heliacal phenomena");
    assert_eq!(mars_hel.target_name, "Mars");
    assert!(mars_hel.arcus_visionis_deg > 8.0 && mars_hel.arcus_visionis_deg < 18.0);
}

#[test]
fn test_phase6_r8_saros_path_width_local_and_range_eclipses() {
    let engine = load_engine();

    // 1. 1999-08-11 Total Solar Eclipse: Saros 145, Member 21
    let sol_1999 = engine
        .next_solar_eclipse(
            TimeInput::from_tt_jd(2_451_400.0),
            false,
            None,
            EclipseOptions::default(),
        )
        .expect("1999-08-11 solar eclipse");
    assert_eq!(sol_1999.kind, SolarEclipseKind::Total);
    assert_eq!(sol_1999.saros_series, 145);
    assert_eq!(sol_1999.saros_member, 21);
    assert!(sol_1999.central_duration_sec.unwrap_or(0.0) > 100.0);
    assert!(sol_1999.path_width_km.unwrap_or(0.0) > 80.0);

    // 2. 2017-08-21 Total Solar Eclipse (1 Saros = 223 synodic months after 1999-08-11): Saros 145, Member 22
    let sol_2017 = engine
        .next_solar_eclipse(
            TimeInput::from_tt_jd(2_457_980.0),
            false,
            None,
            EclipseOptions::default(),
        )
        .expect("2017-08-21 solar eclipse");
    assert_eq!(sol_2017.kind, SolarEclipseKind::Total);
    assert_eq!(sol_2017.saros_series, 145);
    assert_eq!(sol_2017.saros_member, 22);

    // 3. 2024-04-08 Total Solar Eclipse: Saros 139, Member 30
    let sol_2024 = engine
        .next_solar_eclipse(
            TimeInput::from_tt_jd(2_460_400.0),
            false,
            None,
            EclipseOptions::default(),
        )
        .expect("2024-04-08 solar eclipse");
    assert_eq!(sol_2024.kind, SolarEclipseKind::Total);
    assert_eq!(sol_2024.saros_series, 139);
    assert_eq!(sol_2024.saros_member, 30);
    let dur_2024 = sol_2024
        .central_duration_sec
        .expect("2024 central duration");
    let width_2024 = sol_2024.path_width_km.expect("2024 path width");
    assert!(
        dur_2024 > 240.0 && dur_2024 < 290.0,
        "2024-04-08 max central duration should be ~268s, got {dur_2024:.1}s"
    );
    assert!(
        width_2024 > 170.0 && width_2024 < 220.0,
        "2024-04-08 path width should be ~185-199 km, got {width_2024:.1} km"
    );

    // 4. R8.1: next_local_solar_eclipse for Dallas, TX (32.7767 N, -96.7970 E) starting in 2023
    let dallas = Observer::new(32.7767, -96.7970, 130.0).expect("Dallas");
    let local_dallas = engine
        .next_local_solar_eclipse(
            TimeInput::from_tt_jd(2_460_000.0), // early 2023
            dallas,
            false,
            EclipseOptions::default(),
        )
        .expect("next local solar eclipse at Dallas");
    let loc = local_dallas.local.expect("local circumstances populated");
    assert!(loc.visible_at_max);
    assert!(loc.magnitude > 0.70);

    // 5. R8.3: solar_eclipses_in_range and lunar_eclipses_in_range over 2024..2026
    let sols = engine
        .solar_eclipses_in_range(
            TimeInput::from_tt_jd(2_460_310.5),
            TimeInput::from_tt_jd(2_461_041.5),
            None,
            EclipseOptions::default(),
        )
        .expect("solar_eclipses_in_range");
    assert!(sols.len() >= 4);

    let luns = engine
        .lunar_eclipses_in_range(
            TimeInput::from_tt_jd(2_460_310.5),
            TimeInput::from_tt_jd(2_461_041.5),
            EclipseOptions::default(),
        )
        .expect("lunar_eclipses_in_range");
    assert!(luns.len() >= 4);
    for l in &luns {
        assert!(l.saros_series > 100 && l.saros_series < 180);
        assert!(l.saros_member > 0 && l.saros_member < 90);
    }
}

#[test]
fn test_phase6_r9_3_houses_batch() {
    let engine = load_engine();
    let london = Observer::new(51.5074, -0.1278, 0.0).expect("London");
    let epochs: Vec<TimeInput> = (0..50)
        .map(|i| TimeInput::from_tt_jd(JD_J2000 + f64::from(i) * 0.25))
        .collect();
    let opts = HouseOptions::default();

    let batch = engine
        .houses_batch(&epochs, london, opts)
        .expect("houses_batch");
    assert_eq!(batch.len(), 50);

    for (ep, b_res) in epochs.iter().zip(batch.iter()) {
        let single = engine.houses(*ep, london, opts).expect("single houses");
        assert!((b_res.angles.asc_deg - single.angles.asc_deg).abs() < 1e-12);
        assert!((b_res.cusps_deg[0] - single.cusps_deg[0]).abs() < 1e-12);
    }
}

#[test]
fn test_phase6_swe_adapter_migration_layer() {
    let engine = load_engine();
    let mut swe = SweAdapter::new(engine);

    // 1. Major planets, lunar nodes/apogees, Chiron, Ceres, Hamburg Cupido, and SE_ECL_NUT
    let (sun, flag) = swe
        .calc_ut(JD_J2000, SE_SUN, SEFLG_SPEED)
        .expect("SE_SUN");
    assert_eq!(flag.cast_unsigned() & SEFLG_SPEED, SEFLG_SPEED);
    assert!((sun[0] - 280.3689).abs() < 0.02);

    let (moon, _) = swe.calc_ut(JD_J2000, SE_MOON, SEFLG_SPEED).expect("SE_MOON");
    assert!(moon[3] > 11.0 && moon[3] < 16.0);

    let (mean_node, _) = swe
        .calc_ut(JD_J2000, SE_MEAN_NODE, SEFLG_SPEED)
        .expect("SE_MEAN_NODE");
    let (true_node, _) = swe
        .calc_ut(JD_J2000, SE_TRUE_NODE, SEFLG_SPEED)
        .expect("SE_TRUE_NODE");
    assert!((mean_node[0] - true_node[0]).abs() < 2.5);

    let (mean_apog, _) = swe
        .calc_ut(JD_J2000, SE_MEAN_APOG, SEFLG_SPEED)
        .expect("SE_MEAN_APOG");
    let (oscu_apog, _) = swe
        .calc_ut(JD_J2000, SE_OSCU_APOG, SEFLG_SPEED)
        .expect("SE_OSCU_APOG");
    assert!(mean_apog[0] >= 0.0 && oscu_apog[0] >= 0.0);

    let (chiron, _) = swe
        .calc_ut(JD_J2000, SE_CHIRON, SEFLG_SPEED)
        .expect("SE_CHIRON");
    let (ceres_legacy, _) = swe
        .calc_ut(JD_J2000, SE_CERES, SEFLG_SPEED)
        .expect("SE_CERES");
    let (ceres_offset, _) = swe
        .calc_ut(JD_J2000, SE_AST_OFFSET + 1, SEFLG_SPEED)
        .expect("SE_AST_OFFSET + 1");
    assert!(chiron[0] >= 0.0);
    assert!((ceres_legacy[0] - ceres_offset[0]).abs() < 1e-12);

    let (cupido, _) = swe
        .calc_ut(JD_J2000, SE_CUPIDO, SEFLG_SPEED)
        .expect("SE_CUPIDO");
    assert!(cupido[0] >= 0.0 && cupido[0] < 360.0);

    let (ecl_nut, _) = swe.calc_ut(JD_J2000, SE_ECL_NUT, 0).expect("SE_ECL_NUT");
    assert!((ecl_nut[0] - 23.439).abs() < 0.01);

    // 2. Flags: SEFLG_EQUATORIAL, SEFLG_HELCTR, SEFLG_TRUEPOS, SEFLG_J2000, SEFLG_NONUT, SEFLG_RADIANS, SEFLG_SIDEREAL, SEFLG_TOPOCTR
    let (sun_eq, _) = swe
        .calc(JD_J2000, SE_SUN, SEFLG_SPEED | SEFLG_EQUATORIAL)
        .expect("SEFLG_EQUATORIAL");
    assert!((sun_eq[1] - (-23.03)).abs() < 0.2);

    let (sun_rad, _) = swe
        .calc(JD_J2000, SE_SUN, SEFLG_SPEED | SEFLG_RADIANS)
        .expect("SEFLG_RADIANS");
    let (sun_deg, _) = swe
        .calc(JD_J2000, SE_SUN, SEFLG_SPEED)
        .expect("SEFLG_SPEED");
    assert!((sun_rad[0] - sun_deg[0].to_radians()).abs() < 1e-12);

    let _mars_hel = swe
        .calc(
            JD_J2000,
            lxephemeris_core::SE_MARS,
            SEFLG_SPEED | SEFLG_HELCTR | SEFLG_J2000 | SEFLG_TRUEPOS | SEFLG_NONUT,
        )
        .expect("SEFLG_HELCTR | SEFLG_J2000 | SEFLG_TRUEPOS | SEFLG_NONUT");

    swe.set_sid_mode(SE_SIDM_LAHIRI).expect("set_sid_mode Lahiri");
    let ayan_lahiri = swe.get_ayanamsa(JD_J2000).expect("get_ayanamsa");
    let (sun_mean, _) = swe
        .calc(JD_J2000, SE_SUN, SEFLG_SPEED | SEFLG_NONUT)
        .expect("SEFLG_NONUT");
    let (sun_sid, _) = swe
        .calc(JD_J2000, SE_SUN, SEFLG_SPEED | SEFLG_SIDEREAL)
        .expect("SEFLG_SIDEREAL");
    let expected_sid = (sun_mean[0] - ayan_lahiri).rem_euclid(360.0);
    assert!((sun_sid[0] - expected_sid).abs() < 1e-12);

    swe.set_sid_mode(SE_SIDM_KRISHNAMURTI)
        .expect("set_sid_mode KP");
    let ayan_kp = swe.get_ayanamsa_ut(JD_J2000).expect("get_ayanamsa_ut KP");
    assert!((ayan_kp - ayan_lahiri).abs() < 0.2);

    swe.set_topo(-0.1278, 51.5074, 0.0).expect("set_topo London");
    let (moon_topo, _) = swe
        .calc_ut(JD_J2000, SE_MOON, SEFLG_SPEED | SEFLG_TOPOCTR)
        .expect("Moon topocentric");
    assert!((moon_topo[0] - moon[0]).abs() > 0.05);

    // 3. Fixed stars, houses_ex, azalt, and rise_trans via SweAdapter
    let (star_out, star_name, _) = swe
        .fixstar2_ut("Regulus", JD_J2000, SEFLG_SPEED)
        .expect("fixstar2_ut Regulus");
    assert_eq!(star_name, "Regulus");
    assert!((star_out[0] - 149.829).abs() < 0.01);
    let mag = swe.fixstar2_mag("Regulus").expect("fixstar2_mag Regulus");
    assert!((mag - 1.36).abs() < 0.05);

    let (_, ascmc_trop) = swe
        .houses_ex(JD_J2000, 51.5074, -0.1278, 'P', 0)
        .expect("houses_ex tropical");
    let (_, ascmc_sid) = swe
        .houses_ex(JD_J2000, 51.5074, -0.1278, 'P', SEFLG_SIDEREAL)
        .expect("houses_ex sidereal");
    let (retflag_kp, ayan_kp_true) = swe
        .get_ayanamsa_ex_ut(JD_J2000, 0)
        .expect("get_ayanamsa_ex_ut KP");
    assert!(retflag_kp >= 0);
    let diff_asc = (ascmc_trop[0] - ascmc_sid[0]).rem_euclid(360.0);
    assert!((diff_asc - ayan_kp_true).abs() < 1e-12);
    assert!((diff_asc - ayan_kp).abs() < 0.005);

    // Moshier fallback on OutOfCoverage even when SEFLG_SWIEPH is passed
    let (sun_deep_past, retflag_moshier) = swe
        .calc_ut(1721425.5, SE_SUN, SEFLG_SWIEPH | SEFLG_SPEED)
        .expect("swe.calc_ut Moshier fallback for 1 CE");
    assert_ne!(retflag_moshier & (SEFLG_MOSEPH as i32), 0);
    assert!((0.0..360.0).contains(&sun_deep_past[0]));

    // Utility and eclipse/heliacal methods on SweAdapter
    swe.set_ephe_path("/tmp/ephe");
    swe.close();
    assert_eq!(SweAdapter::get_planet_name(SE_SUN), "Sun");
    assert_eq!(SweAdapter::get_planet_name(SE_TRUE_NODE), "True_Node");
    assert_eq!(SweAdapter::house_name('P'), "Placidus");

    let (sol_flags, sol_tret) = swe
        .sol_eclipse_when_glob(2451395.0, SEFLG_SWIEPH, 0, false)
        .expect("swe.sol_eclipse_when_glob");
    assert_ne!(sol_flags & SE_ECL_TOTAL, 0);
    assert!(sol_tret[0] > 2451395.0 && sol_tret[0] < 2451410.0);

    let (_, _, sol_attr) = swe
        .sol_eclipse_where(sol_tret[0], SEFLG_SWIEPH)
        .expect("swe.sol_eclipse_where");
    assert!(sol_attr[0] > 0.0);

    let (_, loc_tret, loc_attr) = swe
        .sol_eclipse_when_loc(2451395.0, [11.58, 48.14, 500.0], SEFLG_SWIEPH, false)
        .expect("swe.sol_eclipse_when_loc");
    assert!((loc_tret[0] - sol_tret[0]).abs() < 0.5);
    assert!(loc_attr[0] > 0.9);

    let (lun_flags, lun_tret) = swe
        .lun_eclipse_when(2451555.0, SEFLG_SWIEPH, 0, false)
        .expect("swe.lun_eclipse_when");
    assert_ne!(lun_flags & SE_ECL_TOTAL, 0);
    let (_, lun_attr) = swe
        .lun_eclipse_how(lun_tret[0], [-0.1278, 51.5074, 0.0], SEFLG_SWIEPH)
        .expect("swe.lun_eclipse_how");
    assert!(lun_attr[0] > 1.0);

    let (hel_tret, _) = swe
        .heliacal_ut(
            JD_J2000,
            [31.2357, 30.0444, 0.0],
            [1013.25, 20.0, 40.0, 0.25],
            [25.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            "Venus",
            SE_HELIACAL_RISING,
            SE_HELFLAG_HIGH_PRECISION,
        )
        .expect("swe.heliacal_ut");
    assert!(hel_tret[0] > JD_J2000);

    let hz_ecl = swe
        .azalt(
            JD_J2000,
            SE_ECL2HOR,
            [-0.1278, 51.5074, 0.0],
            1013.25,
            15.0,
            [sun[0], sun[1]],
        )
        .expect("swe.azalt SE_ECL2HOR");
    let hz_equ = swe
        .azalt(
            JD_J2000,
            SE_EQU2HOR,
            [-0.1278, 51.5074, 0.0],
            1013.25,
            15.0,
            [sun_eq[0], sun_eq[1]],
        )
        .expect("swe.azalt SE_EQU2HOR");
    assert!((hz_ecl.azimuth_deg - hz_equ.azimuth_deg).abs() < 0.05);

    let rt_upper = swe
        .rise_trans(
            JD_J2000,
            SE_SUN,
            None,
            0,
            SE_CALC_RISE,
            [-0.1278, 51.5074, 0.0],
            1013.25,
            15.0,
        )
        .expect("swe.rise_trans upper")
        .expect("sun rises");
    let rt_center = swe
        .rise_trans(
            JD_J2000,
            SE_SUN,
            None,
            0,
            SE_CALC_RISE | SE_BIT_DISC_CENTER,
            [-0.1278, 51.5074, 0.0],
            1013.25,
            15.0,
        )
        .expect("swe.rise_trans center")
        .expect("sun rises");
    let rt_bottom = swe
        .rise_trans(
            JD_J2000,
            SE_SUN,
            None,
            0,
            SE_CALC_RISE | SE_BIT_DISC_BOTTOM | SE_BIT_NO_REFRACTION,
            [-0.1278, 51.5074, 0.0],
            1013.25,
            15.0,
        )
        .expect("swe.rise_trans bottom no refraction")
        .expect("sun rises");
    assert!(rt_upper < rt_center && rt_center < rt_bottom);
}
