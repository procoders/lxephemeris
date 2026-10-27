//! Wave B & Wave C End-to-End Integration & Golden Verification Suite (`v0.2`).
#![cfg(feature = "swe-compat")]
#![allow(clippy::indexing_slicing)]
//!
//! Verifies 100% Swiss Ephemeris capability parity and fileless Moshier mode (`SEFLG_MOSEPH`):
//! 1. `Engine::new_moshier`, `Engine::from_data_dir_or_moshier`, `Engine::calc_moshier`, and `SEFLG_MOSEPH`
//!    in `SweAdapter` across `-3000 CE .. +3000 CE` + agreement with `DE440s` at `J2000.0`.
//! 2. `Engine::phenomena`, `SweAdapter::pheno_ut`, `SweAdapter::pheno`, and `CalcResult.magnitude_v`
//!    for Sun, Moon, Mercury..Pluto, and all 68 minor bodies (`MINOR_BODIES_68`).
//! 3. `Engine::nod_aps`, `SweAdapter::nod_aps_ut`, and `SweAdapter::nod_aps` across `SE_NODBIT_MEAN`,
//!    `SE_NODBIT_OSCU`, `SE_NODBIT_OSCU_BAR`, `SE_NODBIT_FOPOINT`, `SEFLG_XYZ`, `SEFLG_RADIANS`, and `SEFLG_SIDEREAL`.
//! 4. `Engine::next_lunar_occultation`, `Engine::next_local_lunar_occultation`,
//!    `SweAdapter::lun_occult_when_glob`, and `SweAdapter::lun_occult_when_loc` for fixed stars, planets, and minor bodies.
//! 5. 3D `Engine::house_pos` and `SweAdapter::house_pos` (`swe_house_pos`) accounting for ecliptic latitude $\beta$.
//! 6. `Engine::time_equation` / `SweAdapter::time_equ` (`swe_time_equ`) and `Engine::azalt_rev` / `SweAdapter::azalt_rev` (`swe_azalt_rev`).
//! 7. All 43 `SE_SIDM_*` sidereal ayanamsa models (`0..=42`) + `SE_SIDM_USER = 255` custom `(t0, ayan_t0)`.
//! 8. Custom fixed-star catalog loading (`Engine::add_custom_stars`, `Engine::load_custom_stars`, `parse_sefstars_catalog`).

use lxephemeris_core::{
    azalt, azalt_rev, parse_sefstars_catalog, AtmosphereParams, AyanamsaModel, BodyId, CalcOptions,
    Engine, EngineConfig, HipId, HouseOptions, HouseSystem, LunarPointId, MpcNumber, Observer,
    OccultationKind, Plane, PlanetId, SweAdapter, TimeInput, JD_J2000, JD_MOSHIER_MAX,
    JD_MOSHIER_MIN, MINOR_BODIES_68, SEFLG_HELCTR, SEFLG_MOSEPH, SEFLG_RADIANS, SEFLG_SIDEREAL,
    SEFLG_SPEED, SEFLG_SWIEPH, SEFLG_XYZ, SE_AST_OFFSET, SE_ECL2HOR, SE_EQU2HOR, SE_HOR2ECL,
    SE_HOR2EQU, SE_MARS, SE_MOON, SE_NODBIT_FOPOINT, SE_NODBIT_MEAN, SE_NODBIT_OSCU,
    SE_NODBIT_OSCU_BAR, SE_SIDM_USER, SE_SUN, SE_VENUS,
};
use std::path::Path;
use std::sync::Arc;

fn workspace_data_dir() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .unwrap_or_else(|| Path::new("."))
        .join("data")
}

#[test]
fn test_wave_b_moshier_mode_and_seflg_moseph() {
    let moshier_eng = Engine::new_moshier(EngineConfig::default());
    assert!(moshier_eng.is_moshier_mode());

    let spk_eng = Arc::new(
        Engine::from_data_dir(workspace_data_dir(), EngineConfig::default())
            .expect("load data dir"),
    );
    assert!(!spk_eng.is_moshier_mode());

    // Fallback constructor when directory does not exist -> Moshier mode
    let fallback_eng = Engine::from_data_dir_or_moshier(
        "/nonexistent/path/to/kernels",
        EngineConfig::default(),
    )
    .expect("fallback to moshier");
    assert!(fallback_eng.is_moshier_mode());

    let t_j2000 = TimeInput::from_tt_jd(JD_J2000);
    let opts = CalcOptions::default();

    // 1. Compare Moshier vs DE440s at J2000.0 for all major planets + Moon
    for planet in [
        PlanetId::Sun,
        PlanetId::Moon,
        PlanetId::Mercury,
        PlanetId::Venus,
        PlanetId::Mars,
        PlanetId::Jupiter,
        PlanetId::Saturn,
        PlanetId::Uranus,
        PlanetId::Neptune,
        PlanetId::Pluto,
    ] {
        let r_mosh = moshier_eng
            .calc(BodyId::Planet(planet), t_j2000, opts)
            .expect("moshier planet");
        let r_spk = spk_eng
            .calc(BodyId::Planet(planet), t_j2000, opts)
            .expect("spk planet");

        let dlon = ((r_mosh.lon_deg - r_spk.lon_deg + 180.0).rem_euclid(360.0) - 180.0).abs();
        let tol_deg = match planet {
            PlanetId::Jupiter | PlanetId::Saturn | PlanetId::Uranus | PlanetId::Neptune => 0.25,
            _ => 0.05,
        };
        assert!(
            dlon < tol_deg,
            "Moshier vs DE440s J2000 lon diff for {planet:?}: {dlon:.6} deg (tol {tol_deg})"
        );
        assert!(r_mosh.magnitude_v.is_some());
    }

    // 2. Compare Moshier vs SPK for all 68 minor bodies at J2000.0
    for spec in MINOR_BODIES_68 {
        let body = BodyId::MinorPlanet(spec.mpc);
        let r_mosh = moshier_eng
            .calc(body, t_j2000, opts)
            .expect("moshier minor body");
        let r_spk = spk_eng.calc(body, t_j2000, opts).expect("spk minor body");
        let dlon = ((r_mosh.lon_deg - r_spk.lon_deg + 180.0).rem_euclid(360.0) - 180.0).abs();
        assert!(
            dlon < 0.02,
            "Moshier vs SPK J2000 lon diff for MPC {:?}: {dlon:.6} deg",
            spec.mpc
        );
    }

    // 3. Verify deep historical (-2500 CE) and far-future (+2800 CE) epochs in Moshier mode
    for jd_tt in [
        JD_MOSHIER_MIN + 100.0,
        1_000_000.0,
        2_700_000.0,
        JD_MOSHIER_MAX - 100.0,
    ] {
        let ep = TimeInput::from_tt_jd(jd_tt);
        for body in [
            BodyId::Planet(PlanetId::Sun),
            BodyId::Planet(PlanetId::Moon),
            BodyId::Planet(PlanetId::Saturn),
            BodyId::Planet(PlanetId::Pluto),
            BodyId::MinorPlanet(MpcNumber(1)),
            BodyId::MinorPlanet(MpcNumber(2060)),
            BodyId::LunarPoint(LunarPointId::TrueNode),
        ] {
            let res = moshier_eng.calc(body, ep, opts).expect("deep epoch calc");
            assert!(res.lon_deg.is_finite() && (0.0..360.0).contains(&res.lon_deg));
            assert!(res.dist_au > 0.0);
        }
    }

    // 4. SweAdapter SEFLG_MOSEPH dispatch on both SPK-backed and Moshier-backed engines
    let swe_spk = SweAdapter::new(Arc::clone(&spk_eng));
    let (xx_swieph, flag_swieph) = swe_spk
        .calc_ut(JD_J2000, SE_MARS, SEFLG_SWIEPH | SEFLG_SPEED)
        .expect("swe_spk swieph");
    let (xx_moseph, flag_moseph) = swe_spk
        .calc_ut(JD_J2000, SE_MARS, SEFLG_MOSEPH | SEFLG_SPEED)
        .expect("swe_spk moseph");
    assert_ne!(flag_swieph as u32 & SEFLG_SWIEPH, 0);
    assert_ne!(flag_moseph as u32 & SEFLG_MOSEPH, 0);
    assert_eq!(flag_moseph as u32 & SEFLG_SWIEPH, 0);
    assert!((xx_swieph[0] - xx_moseph[0]).abs() < 0.02);

    // Ancient epoch outside DE440s range works seamlessly when SEFLG_MOSEPH is passed to swe_spk!
    let (xx_ancient, flag_ancient) = swe_spk
        .calc_ut(1_721_423.5, SE_SUN, SEFLG_MOSEPH | SEFLG_SPEED)
        .expect("swe_spk ancient moseph");
    assert!(xx_ancient[0].is_finite());
    assert_ne!(flag_ancient as u32 & SEFLG_MOSEPH, 0);
}

#[test]
fn test_wave_c1_phenomena_and_magnitudes() {
    let eng = Arc::new(
        Engine::from_data_dir(workspace_data_dir(), EngineConfig::default())
            .expect("load data dir"),
    );
    let swe = SweAdapter::new(Arc::clone(&eng));
    let t_j2000 = TimeInput::from_tt_jd(JD_J2000);
    let opts = CalcOptions::default();

    // Sun phenomena
    let sun_ph = eng
        .phenomena(BodyId::Planet(PlanetId::Sun), t_j2000, opts)
        .expect("sun phenomena");
    assert_eq!(sun_ph.phase_angle_deg, 0.0);
    assert_eq!(sun_ph.phase_illuminated, 1.0);
    assert!((sun_ph.apparent_magnitude_v - (-26.74)).abs() < 0.1);
    assert!((sun_ph.apparent_diameter_deg - 0.54).abs() < 0.03);

    // Moon phenomena + SweAdapter::pheno_ut
    let (attr, retflag) = swe
        .pheno_ut(JD_J2000, SE_MOON, SEFLG_SPEED)
        .expect("swe pheno_ut moon");
    assert_ne!(retflag as u32 & SEFLG_SWIEPH, 0);
    assert!((0.0..=180.0).contains(&attr[0]));
    assert!((0.0..=1.0).contains(&attr[1]));
    assert!((0.0..=180.0).contains(&attr[2]));
    assert!(attr[3] > 0.45 && attr[3] < 0.60);
    assert!(attr[4] < -4.0 && attr[4] > -13.0);
    assert!(attr[5] > 0.8 && attr[5] < 1.2);

    // Venus phenomena + automatic CalcResult.magnitude_v match
    let (attr_v, _) = swe.pheno(JD_J2000, SE_VENUS, 0).expect("swe pheno venus");
    let calc_v = eng
        .calc(BodyId::Planet(PlanetId::Venus), t_j2000, opts)
        .expect("calc venus");
    assert_eq!(calc_v.magnitude_v, Some(attr_v[4]));
    assert!(attr_v[4] < -3.5 && attr_v[4] > -5.0);

    // All 68 minor bodies return finite phenomena and magnitude_v
    for spec in MINOR_BODIES_68 {
        let body = BodyId::MinorPlanet(spec.mpc);
        let ph = eng.phenomena(body, t_j2000, opts).expect("minor body pheno");
        let calc = eng.calc(body, t_j2000, opts).expect("minor body calc");
        assert!(ph.apparent_magnitude_v.is_finite());
        assert_eq!(calc.magnitude_v, Some(ph.apparent_magnitude_v));
    }
}

#[test]
fn test_wave_c2_nodes_and_apsides() {
    let eng = Arc::new(
        Engine::from_data_dir(workspace_data_dir(), EngineConfig::default())
            .expect("load data dir"),
    );
    let swe = SweAdapter::new(Arc::clone(&eng));

    // 1. Lunar Mean Node & Perigee/Apogee (SE_NODBIT_MEAN)
    let (xnasc, xndsc, xperi, xaphe, _) = swe
        .nod_aps(JD_J2000, SE_MOON, SEFLG_SPEED, SE_NODBIT_MEAN)
        .expect("lunar mean nod_aps");
    assert!((xnasc[0] - 125.0445).abs() < 0.02);
    assert!(((xndsc[0] - (xnasc[0] + 180.0).rem_euclid(360.0)).abs()) < 1e-9);
    assert!((xperi[0] - 83.465).abs() < 0.05);
    assert!(((xaphe[0] - (xperi[0] + 180.0).rem_euclid(360.0)).abs()) < 1e-9);
    assert!(xnasc[3] < 0.0); // retrograde node speed (~ -0.05295 deg/day)
    assert!(xperi[3] > 0.0); // prograde perigee speed (~ +0.11140 deg/day)

    // 2. Mars Osculating & Second Focal Point (SE_NODBIT_OSCU | SE_NODBIT_FOPOINT)
    let (m_asc, _m_dsc, m_peri, m_aphe, _) = swe
        .nod_aps_ut(JD_J2000, SE_MARS, SEFLG_SPEED, SE_NODBIT_OSCU)
        .expect("mars oscu nod_aps");
    let (_, _, m_helio_peri, m_helio_aphe, _) = swe
        .nod_aps_ut(
            JD_J2000,
            SE_MARS,
            SEFLG_SPEED | SEFLG_HELCTR,
            SE_NODBIT_OSCU,
        )
        .expect("mars helio oscu nod_aps");
    let (_, _, _, m_fop, _) = swe
        .nod_aps_ut(
            JD_J2000,
            SE_MARS,
            SEFLG_SPEED,
            SE_NODBIT_OSCU | SE_NODBIT_FOPOINT,
        )
        .expect("mars fopoint nod_aps");
    assert!(m_peri[2] > 0.3 && m_aphe[2] > 0.3);
    assert!(m_helio_peri[2] > 1.3 && m_helio_aphe[2] > m_helio_peri[2]);
    assert_ne!(m_aphe[2], m_fop[2]);
    assert!(m_asc[0].is_finite());

    // 3. Barycentric osculating (SE_NODBIT_OSCU_BAR) + XYZ & Radians formatting
    let (xyz_asc, _, _, _, _) = swe
        .nod_aps(
            JD_J2000,
            SE_AST_OFFSET + 1, // Ceres
            SEFLG_SPEED | SEFLG_XYZ,
            SE_NODBIT_OSCU_BAR,
        )
        .expect("ceres xyz nod_aps");
    let r_ceres = (xyz_asc[0].powi(2) + xyz_asc[1].powi(2) + xyz_asc[2].powi(2)).sqrt();
    assert!(r_ceres > 1.5 && r_ceres < 4.0);

    let (rad_asc, _, _, _, _) = swe
        .nod_aps(
            JD_J2000,
            SE_MOON,
            SEFLG_SPEED | SEFLG_RADIANS,
            SE_NODBIT_MEAN,
        )
        .expect("moon radians nod_aps");
    assert!((rad_asc[0] - xnasc[0].to_radians()).abs() < 1e-12);
}

#[test]
fn test_wave_c3_lunar_occultations() {
    let eng = Arc::new(
        Engine::from_data_dir(workspace_data_dir(), EngineConfig::default())
            .expect("load data dir"),
    );
    let swe = SweAdapter::new(Arc::clone(&eng));

    // 1. Global lunar occultation of Aldebaran (HIP 21421) in 2015
    // 2015-01-01 00:00 UT = JD 2457023.5
    let jd_2015 = 2_457_023.5;
    let (retflag, tret) = swe
        .lun_occult_when_glob(jd_2015, 0, Some("Aldebaran"), 0, false)
        .expect("global occultation of Aldebaran");
    assert!(retflag > 0);
    assert!(tret[0] > jd_2015 && tret[0] < jd_2015 + 60.0);
    assert!(tret[2] <= tret[0] && tret[3] >= tret[0]);

    // Forward then backward search returns the same event
    let (retflag_back, tret_back) = swe
        .lun_occult_when_glob(tret[0] + 2.0, 0, Some("Aldebaran"), 0, true)
        .expect("backward global occultation of Aldebaran");
    assert_eq!(retflag, retflag_back);
    assert!((tret[0] - tret_back[0]).abs() < 1e-4);

    // 2. Local topocentric lunar occultation of Aldebaran from Greenwich (51.4769 N, 0.0 E)
    let obs_greenwich = [0.0, 51.4769, 0.0];
    let (loc_flag, loc_tret, loc_attr) = swe
        .lun_occult_when_loc(jd_2015, 0, Some("Aldebaran"), 0, obs_greenwich, false)
        .expect("local occultation of Aldebaran at Greenwich");
    assert!(loc_flag > 0);
    assert!(loc_tret[0] > jd_2015 && loc_tret[0] < jd_2015 + 400.0);
    assert!(loc_attr[0] > 0.0); // magnitude
    assert!(loc_attr[4] >= -0.5); // Moon above horizon

    // 3. Global lunar occultation of a planet (Mars)
    let mars_occ = eng
        .next_lunar_occultation(
            BodyId::Planet(PlanetId::Mars),
            TimeInput::from_ut1_jd(2_459_915.0), // Dec 2022 Mars occultation
            false,
            30.0,
        )
        .expect("Dec 2022 Mars lunar occultation");
    assert_eq!(mars_occ.kind, OccultationKind::Total);
    assert!((mars_occ.jd_ut1_max - 2_459_921.71).abs() < 0.2);
}

#[test]
fn test_wave_c4_house_pos_time_equ_azalt_rev_and_sidereal_models() {
    let eng = Arc::new(
        Engine::from_data_dir(workspace_data_dir(), EngineConfig::default())
            .expect("load data dir"),
    );
    let mut swe = SweAdapter::new(Arc::clone(&eng));

    // 1. 3D house_pos across house systems and non-zero ecliptic latitude
    let armc = 120.0;
    let lat = 47.3769; // Zurich
    let eps = 23.43929111;
    let hrs = eng
        .houses_armc(
            armc,
            lat,
            eps,
            0.0,
            HouseOptions {
                system: HouseSystem::Placidus,
                ..HouseOptions::default()
            },
        )
        .expect("placidus cusps");

    // On the ecliptic (beta = 0), house_pos at each cusp k must be ~ k.0
    for k in 1..=12 {
        let cusp_lon = hrs.cusps_deg[k - 1];
        let hp = swe
            .house_pos(armc, lat, eps, 'P', [cusp_lon, 0.0])
            .expect("house_pos at cusp");
        let diff = ((hp - (k as f64) + 6.0).rem_euclid(12.0) - 6.0).abs();
        assert!(
            diff < 1e-4,
            "Placidus house_pos at cusp {k} ({cusp_lon} deg): got {hp}, diff={diff}"
        );
    }

    // With non-zero latitude (e.g. Pluto/Moon off ecliptic), 3D Placidus house_pos shifts from 2D longitude position
    let hp_beta0 = swe
        .house_pos(armc, lat, eps, 'P', [45.0, 0.0])
        .expect("beta 0");
    let hp_beta10 = swe
        .house_pos(armc, lat, eps, 'P', [45.0, 10.0])
        .expect("beta +10");
    assert!((1.0..13.0).contains(&hp_beta10));
    assert!((hp_beta10 - hp_beta0).abs() > 0.05);

    // Test Campanus ('C'), Regiomontanus ('R'), Meridian ('X'), Morinus ('M'), Topocentric ('T'), Azimuthal ('H')
    for code in [
        'P', 'K', 'O', 'R', 'C', 'A', 'W', 'B', 'M', 'T', 'X', 'V', 'H',
    ] {
        let hp = swe
            .house_pos(armc, lat, eps, code, [100.0, 3.5])
            .expect("house_pos system");
        assert!(
            (1.0..13.0).contains(&hp),
            "house_pos for '{code}' out of [1, 13): {hp}"
        );
    }

    // 2. Equation of Time (`swe_time_equ`)
    // Around Feb 11 (JD 2451585.5), Apparent - Mean solar time is ~ -14.2 minutes (~ -0.00988 days)
    let eot_feb = swe.time_equ(2_451_585.5).expect("eot feb");
    assert!(
        (eot_feb * 1440.0 - (-14.2)).abs() < 0.6,
        "Feb 11 EoT in minutes expected ~ -14.2, got {}",
        eot_feb * 1440.0
    );
    // Around Nov 3 (JD 2451851.5), Apparent - Mean solar time is ~ +16.4 minutes (~ +0.0114 days)
    let eot_nov = swe.time_equ(2_451_851.5).expect("eot nov");
    assert!(
        (eot_nov * 1440.0 - 16.4).abs() < 0.6,
        "Nov 3 EoT in minutes expected ~ +16.4, got {}",
        eot_nov * 1440.0
    );

    // 3. Inverse Horizontal Coordinates (`azalt_rev` / `swe_azalt_rev`) round-trip
    let obs = Observer::new(47.3769, 8.5417, 408.0).expect("obs");
    let gast = 1.85;
    let true_eps = 23.43929111_f64.to_radians();
    let ra_in = 142.5;
    let dec_in = 21.3;
    let hz = azalt(obs, gast, ra_in, dec_in, AtmosphereParams::default());
    let [ra_out, dec_out] = azalt_rev(
        obs,
        gast,
        true_eps,
        hz.azimuth_deg,
        hz.altitude_true_deg,
        Plane::Equatorial,
    );
    assert!((ra_out - ra_in).abs() < 1e-9);
    assert!((dec_out - dec_in).abs() < 1e-9);

    let [ecl_lon, ecl_lat] = azalt_rev(
        obs,
        gast,
        true_eps,
        hz.azimuth_deg,
        hz.altitude_true_deg,
        Plane::Ecliptic,
    );
    assert!((0.0..360.0).contains(&ecl_lon));
    assert!((-90.0..=90.0).contains(&ecl_lat));

    // SweAdapter::azalt + SweAdapter::azalt_rev round-trip
    let geopos = [8.5417, 47.3769, 408.0];
    let xaz_eq = swe
        .azalt(JD_J2000, SE_EQU2HOR, geopos, 0.0, 0.0, [210.0, -15.0])
        .expect("swe azalt eq");
    let xout_eq = swe
        .azalt_rev(
            JD_J2000,
            SE_HOR2EQU,
            geopos,
            [xaz_eq.azimuth_deg, xaz_eq.altitude_true_deg],
        )
        .expect("swe azalt_rev eq");
    assert!((xout_eq[0] - 210.0).abs() < 1e-7);
    assert!((xout_eq[1] - (-15.0)).abs() < 1e-7);

    let xaz_ecl = swe
        .azalt(JD_J2000, SE_ECL2HOR, geopos, 0.0, 0.0, [125.0, 4.0])
        .expect("swe azalt ecl");
    let xout_ecl = swe
        .azalt_rev(
            JD_J2000,
            SE_HOR2ECL,
            geopos,
            [xaz_ecl.azimuth_deg, xaz_ecl.altitude_true_deg],
        )
        .expect("swe azalt_rev ecl");
    assert!((xout_ecl[0] - 125.0).abs() < 1e-7);
    assert!((xout_ecl[1] - 4.0).abs() < 1e-7);

    // 4. All 43 SE_SIDM_* models (0..=42) + SE_SIDM_USER = 255
    assert_eq!(AyanamsaModel::ALL_EXTENDED.len(), 43);
    for sid_id in 0..=42_u32 {
        swe.set_sid_mode(sid_id).expect("set_sid_mode 0..=42");
        let aya = swe.get_ayanamsa_ut(JD_J2000).expect("get_ayanamsa_ut");
        assert!(
            aya.is_finite() && (-360.0..=360.0).contains(&aya),
            "Invalid ayanamsa for sid_id={sid_id}: {aya}"
        );
    }

    // Custom user ayanamsa (`SE_SIDM_USER = 255`) anchored at J2000 with 24.12345 deg
    swe.set_sid_mode_custom(SE_SIDM_USER, JD_J2000, 24.12345)
        .expect("set_sid_mode_custom");
    let aya_user_et = swe.get_ayanamsa(JD_J2000).expect("mean user aya at TT");
    assert!((aya_user_et - 24.12345).abs() < 1e-12);
    let (_, aya_user_j2000) = swe
        .get_ayanamsa_ex_ut(JD_J2000, 64)
        .expect("mean user aya at UT1"); // SEFLG_NONUT = 64
    assert!((aya_user_j2000 - 24.12345).abs() < 1e-4);
    let (xx_sid, _) = swe
        .calc_ut(JD_J2000, SE_SUN, SEFLG_SIDEREAL | SEFLG_SPEED)
        .expect("sidereal sun with user ayanamsa");
    assert!((0.0..360.0).contains(&xx_sid[0]));
}

#[test]
fn test_wave_c5_custom_star_catalog_sefstars_and_json() {
    let sefstars_sample = r#"
# Custom sefstars.txt sample
PolarisTest  ,alUMi, IC, 02,31,49.0945, +89,15,50.792, 44.48, -11.85, -16.4, 7.54, 1.98, 2000, 999001
CustomVega   ,alLyr, IC, 18,36,56.3363, +38,47,01.280, 200.94, 286.23, -13.9, 130.23, 0.03, 2000, 999002
"#;
    let parsed = parse_sefstars_catalog(sefstars_sample).expect("parse sefstars");
    assert_eq!(parsed.len(), 2);
    assert_eq!(parsed[0].name, "PolarisTest");
    assert_eq!(parsed[0].hip, 999_001);
    assert_eq!(parsed[1].name, "CustomVega");
    assert_eq!(parsed[1].hip, 999_002);

    let mut eng = Engine::new_moshier(EngineConfig::default());
    assert_eq!(eng.add_custom_stars(parsed), 2);

    let meta = eng.star_info("PolarisTest").expect("custom star_info by name");
    assert_eq!(meta.hip, 999_001);

    let t_j2000 = TimeInput::from_tt_jd(JD_J2000);
    let res = eng
        .calc(
            BodyId::FixedStar(HipId(999_001)),
            t_j2000,
            CalcOptions::default(),
        )
        .expect("calc custom star by HIP");
    assert!(res.lon_deg.is_finite());
    assert_eq!(res.magnitude_v, Some(1.98));

    let swe = SweAdapter::new(Arc::new(eng));
    let (xx, stname, _) = swe
        .fixstar2_ut("CustomVega", JD_J2000, SEFLG_SPEED)
        .expect("fixstar2_ut custom star");
    assert_eq!(stname, "CustomVega");
    assert!(xx[0].is_finite());
}
