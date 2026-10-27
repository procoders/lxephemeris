//! Phase 3.2 (23 House Systems + PolarPolicy + Gauquelin) & Phase 4 (Solar and Lunar Eclipses)
//! + Codex Gate 1 Surgical Fixes Verification Suite.

#![allow(clippy::indexing_slicing, clippy::panic)]

use lxephemeris_core::{
    julday, AyanamsaModel, BodyId, CalcOptions, CalendarKind, EclipseOptions, Engine, EngineConfig,
    EphemerisError, HouseOptions, HouseSystem, HouseWarning, LunarEclipseKind, NutationMode,
    Observer, Origin, Plane, PlanetId, PolarCategory, PolarPolicy, ShadowEnlargementRule,
    SolarEclipseKind, TimeInput, TimeTable, Zodiac, JD_J2000,
};
use serde::Deserialize;
use std::fs;
use std::path::{Path, PathBuf};

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf()
}

fn load_de440s_engine() -> Engine {
    let root = workspace_root();
    Engine::from_bsp_paths(&[root.join("data/de440s.bsp")], EngineConfig::default())
        .unwrap_or_else(|e| panic!("Failed to load de440s.bsp: {e}"))
}

fn wrap_diff_deg(a: f64, b: f64) -> f64 {
    ((a - b + 180.0).rem_euclid(360.0) - 180.0).abs()
}

// ---------------------------------------------------------------------------
// Part A: Codex Gate 1 Surgical Fixes Verification
// ---------------------------------------------------------------------------

#[test]
fn test_codex_gate1_leapseconds_dynamic_parsing_and_horizontal_frame() {
    let root = workspace_root();
    let leap_text = fs::read_to_string(root.join("data/iers_leapseconds.txt"))
        .unwrap_or_else(|e| panic!("read iers_leapseconds.txt: {e}"));
    let parsed = TimeTable::parse_leapseconds_table(&leap_text);
    // From 1972 JAN 1 (10.0 s) through 2017 JAN 1 (37.0 s) there are 28 entries
    assert_eq!(parsed.len(), 28);
    assert_eq!(parsed.first().copied(), Some((41_317.0, 10.0)));
    assert_eq!(parsed.last().copied(), Some((57_754.0, 37.0)));

    // Synthetic future leap second in 2028
    let mut extended_text = leap_text.clone();
    extended_text.push_str("\n 2028 JAN  1 =JD 2461771.5  TAI-UTC=  38.0       S + (MJD - 41317.) X 0.0      S\n");
    let parsed_ext = TimeTable::parse_leapseconds_table(&extended_text);
    assert_eq!(parsed_ext.len(), 29);
    assert_eq!(parsed_ext.last().copied(), Some((61_771.0, 38.0)));

    // Verify Plane::Horizontal returns horizontal Cartesian vector and non-zero horizontal speeds
    let engine = load_de440s_engine();
    let obs = Observer::new(51.5074, -0.1276, 25.0).expect("London observer");
    let hz_opts = CalcOptions {
        origin: Origin::Topocentric(obs),
        plane: Plane::Horizontal,
        include_speed: true,
        ..Default::default()
    };
    let sun_hz = engine
        .calc(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_tt_jd(JD_J2000),
            hz_opts,
        )
        .unwrap_or_else(|e| panic!("calc horizontal sun: {e}"));

    // At J2000 (12:00 TT in London), Sun is near Upper Transit (Azimuth ~ 180 deg, Altitude ~ 15.5 deg)
    assert!(
        (sun_hz.lon_deg - 180.0).abs() < 2.0,
        "Expected Sun azimuth near 180 deg at noon in London, got {}",
        sun_hz.lon_deg
    );
    assert!(
        (sun_hz.lat_deg - 15.5).abs() < 1.5,
        "Expected Sun altitude near 15.5 deg at winter noon in London, got {}",
        sun_hz.lat_deg
    );
    // Horizontal azimuth speed at upper transit in northern hemisphere is positive (moving East -> South -> West, ~360+ deg/day)
    assert!(
        sun_hz.lon_speed_deg_day > 200.0,
        "Expected positive horizontal azimuth rate > 200 deg/day, got {}",
        sun_hz.lon_speed_deg_day
    );
    // Cartesian vector in Horizontal frame must match (lon_deg=Az, lat_deg=Alt, dist_au)
    let sph_from_cart = lxephemeris_core::cart_to_spherical(&sun_hz.cartesian_au);
    assert!(wrap_diff_deg(sph_from_cart.lon_deg, sun_hz.lon_deg) < 1e-10);
    assert!((sph_from_cart.lat_deg - sun_hz.lat_deg).abs() < 1e-10);
}

// ---------------------------------------------------------------------------
// Part B: Phase 3.2 — All 23 House Systems, PolarPolicy & Gauquelin Sectors
// ---------------------------------------------------------------------------

#[test]
fn test_all_23_house_systems_temperate_southern_and_dateline() {
    let engine = load_de440s_engine();
    // Summer solstice 2024: Sun declination ~ +23.44 deg
    let jd_solstice =
        julday(2024, 6, 20, 20.85, CalendarKind::Gregorian).expect("solstice jd");
    let t_in = TimeInput::from_tt_jd(jd_solstice);

    let london = Observer::new(51.5074, -0.1276, 0.0).expect("London");
    let sydney = Observer::new(-33.8688, 151.2093, 0.0).expect("Sydney");
    let dateline_pos = Observer::new(35.0, 180.0, 0.0).expect("+180 lon");
    let dateline_neg = Observer::new(35.0, -180.0, 0.0).expect("-180 lon");

    for sys in HouseSystem::ALL {
        let opts = HouseOptions {
            system: sys,
            polar_policy: PolarPolicy::Strict,
            zodiac: Zodiac::Tropical,
            nutation: NutationMode::True,
            include_speed: true,
        };

        // 1. Northern mid-latitude (London)
        let r_lon = engine
            .houses(t_in, london, opts)
            .unwrap_or_else(|e| panic!("London house system '{}' ({sys:?}) failed: {e}", sys.code()));
        assert_eq!(r_lon.requested_house_system, sys);
        assert_eq!(r_lon.effective_house_system, sys);
        assert!(r_lon.fallback_reason.is_none());
        assert!(r_lon.warnings.is_empty());
        for c in r_lon.cusps_deg {
            assert!((0.0..360.0).contains(&c), "invalid cusp {c} for {sys:?}");
        }

        // Verify distinctness of specific house system invariants:
        match sys {
            HouseSystem::WholeSignAries => {
                assert_eq!(r_lon.cusps_deg[0], 0.0);
                assert_eq!(r_lon.cusps_deg[9], 270.0);
            }
            HouseSystem::WholeSign => {
                assert_eq!(
                    r_lon.cusps_deg[0],
                    (r_lon.angles.asc_deg / 30.0).floor() * 30.0
                );
            }
            HouseSystem::EqualFromAsc => {
                assert!(wrap_diff_deg(r_lon.cusps_deg[0], r_lon.angles.asc_deg) < 1e-9);
                assert!(wrap_diff_deg(r_lon.cusps_deg[1], (r_lon.angles.asc_deg + 30.0) % 360.0) < 1e-9);
            }
            HouseSystem::EqualFromMc => {
                assert!(wrap_diff_deg(r_lon.cusps_deg[9], r_lon.angles.mc_deg) < 1e-9);
                assert!(wrap_diff_deg(r_lon.cusps_deg[0], (r_lon.angles.mc_deg + 90.0) % 360.0) < 1e-9);
            }
            HouseSystem::VehlowEqual => {
                assert!(
                    wrap_diff_deg(r_lon.cusps_deg[0], (r_lon.angles.asc_deg - 15.0).rem_euclid(360.0))
                        < 1e-9
                );
            }
            _ => {}
        }

        // 2. Southern hemisphere (Sydney)
        let r_syd = engine
            .houses(t_in, sydney, opts)
            .unwrap_or_else(|e| panic!("Sydney house system '{}' ({sys:?}) failed: {e}", sys.code()));
        for c in r_syd.cusps_deg {
            assert!((0.0..360.0).contains(&c));
        }

        // 3. International Date Line (+180° vs -180° must be bit-identical within 1e-10 deg)
        let r_p180 = engine
            .houses(t_in, dateline_pos, opts)
            .unwrap_or_else(|e| panic!("+180 house system {sys:?}: {e}"));
        let r_m180 = engine
            .houses(t_in, dateline_neg, opts)
            .unwrap_or_else(|e| panic!("-180 house system {sys:?}: {e}"));
        for i in 0..12 {
            assert!(
                wrap_diff_deg(r_p180.cusps_deg[i], r_m180.cusps_deg[i]) < 1e-10,
                "Date line mismatch on cusp {} for {sys:?}",
                i + 1
            );
        }
    }

    // Verify 'I' (Sunshine Treindl) and 'i' (Sunshine Makransky) match at temperate latitude (MODEL_PASSPORTS.md §5.3.16-17)
    let r_big_i = engine
        .houses(
            t_in,
            london,
            HouseOptions {
                system: HouseSystem::SunshineTreindl,
                ..Default::default()
            },
        )
        .expect("Sunshine Treindl");
    let r_small_i = engine
        .houses(
            t_in,
            london,
            HouseOptions {
                system: HouseSystem::SunshineMakransky,
                ..Default::default()
            },
        )
        .expect("Sunshine Makransky");
    assert!(
        wrap_diff_deg(r_big_i.cusps_deg[10], r_small_i.cusps_deg[10]) < 1e-12,
        "Sunshine 'I' and 'i' must agree at mid-latitude where |tan(dec_sun)*tan(lat)| <= 1"
    );
}

#[derive(Debug, Deserialize)]
struct L1bHousesFile {
    case_count: usize,
    system_count: usize,
    cases: Vec<L1bHouseCase>,
}

#[derive(Debug, Deserialize)]
struct L1bHouseCase {
    case_id: String,
    lat_deg: f64,
    ramc_deg: f64,
    eps_true_deg: f64,
    sun_decl_deg: f64,
    angles: Option<L1bAngles>,
    gauquelin_36_sectors: Option<Vec<f64>>,
    systems: std::collections::BTreeMap<String, L1bSystemEntry>,
}

#[derive(Debug, Deserialize)]
#[allow(non_snake_case)]
struct L1bAngles {
    ASC: f64,
    MC: f64,
    ARMC: f64,
    Vertex: f64,
    EqASC: f64,
    CoASC_Koch: f64,
    CoASC_Munkasey: f64,
    PolarASC: f64,
}

#[derive(Debug, Deserialize)]
struct L1bSystemEntry {
    status: String,
    cusps: Option<Vec<f64>>,
}

#[test]
fn test_l1b_houses_23_golden_json() {
    use lxephemeris_core::{gauquelin_sectors, houses_armc};

    let root = workspace_root();
    let text = std::fs::read_to_string(root.join("testdata/golden/l1b_houses_23.json"))
        .expect("read l1b_houses_23.json");
    let golden: L1bHousesFile = serde_json::from_str(&text).expect("parse l1b_houses_23.json");
    assert_eq!(golden.case_count, 11);
    assert_eq!(golden.system_count, 23);

    let mut max_cusp_err_arcsec = 0.0_f64;
    let mut max_angle_err_arcsec = 0.0_f64;
    let mut max_gauq_err_arcsec = 0.0_f64;

    for c in &golden.cases {
        for (code_str, sys_exp) in &c.systems {
            let ch = code_str.chars().next().expect("char");
            let sys = HouseSystem::from_code(ch).expect("valid house code");
            let res = houses_armc(
                c.ramc_deg,
                c.lat_deg,
                c.eps_true_deg,
                c.sun_decl_deg,
                HouseOptions {
                    system: sys,
                    polar_policy: PolarPolicy::Strict,
                    include_speed: false,
                    ..Default::default()
                },
            );

            if sys_exp.status.starts_with("Err") {
                assert!(
                    res.is_err(),
                    "Case {} system '{}' expected {}, got Ok",
                    c.case_id,
                    ch,
                    sys_exp.status
                );
            } else {
                let r = res.unwrap_or_else(|e| {
                    panic!("Case {} system '{}' failed: {e}", c.case_id, ch)
                });
                let exp_cusps = sys_exp.cusps.as_ref().expect("cusps present");
                for (i, (&actual, &expected)) in
                    r.cusps_deg.iter().zip(exp_cusps.iter()).enumerate()
                {
                    let err_arcsec = wrap_diff_deg(actual, expected) * 3_600.0;
                    max_cusp_err_arcsec = max_cusp_err_arcsec.max(err_arcsec);
                    assert!(
                        err_arcsec < 1e-4,
                        "Case {} system '{}' cusp {} mismatch: {err_arcsec:.8}\"",
                        c.case_id,
                        ch,
                        i + 1
                    );
                }
                if let Some(ang) = &c.angles {
                    for (actual, expected) in [
                        (r.angles.asc_deg, ang.ASC),
                        (r.angles.mc_deg, ang.MC),
                        (r.angles.armc_deg, ang.ARMC),
                        (r.angles.vertex_deg, ang.Vertex),
                        (r.angles.equatorial_asc_deg, ang.EqASC),
                        (r.angles.co_asc_koch_deg, ang.CoASC_Koch),
                        (r.angles.co_asc_munkasey_deg, ang.CoASC_Munkasey),
                        (r.angles.polar_asc_deg, ang.PolarASC),
                    ] {
                        let a_err = wrap_diff_deg(actual, expected) * 3_600.0;
                        max_angle_err_arcsec = max_angle_err_arcsec.max(a_err);
                        assert!(
                            a_err < 1e-4,
                            "Case {} angle mismatch: {a_err:.8}\"",
                            c.case_id
                        );
                    }
                }
            }
        }

        if let Some(exp_gauq) = &c.gauquelin_36_sectors {
            let g_res = gauquelin_sectors(
                c.ramc_deg,
                c.lat_deg,
                c.eps_true_deg,
                PolarPolicy::Strict,
            )
            .expect("gauquelin");
            for (i, (&actual, &expected)) in
                g_res.sectors_deg.iter().zip(exp_gauq.iter()).enumerate()
            {
                let g_err = wrap_diff_deg(actual, expected) * 3_600.0;
                max_gauq_err_arcsec = max_gauq_err_arcsec.max(g_err);
                assert!(
                    g_err < 1e-4,
                    "Case {} gauquelin sector {} mismatch: {g_err:.8}\"",
                    c.case_id,
                    i + 1
                );
            }
        }
    }

    println!(
        "L1B 23 House Systems (11 cases x 23 = 253 evals): max_cusp_err = {max_cusp_err_arcsec:.6}\", max_angle_err = {max_angle_err_arcsec:.6}\", max_gauq_err = {max_gauq_err_arcsec:.6}\""
    );
}

#[test]
fn test_polar_policy_strict_vs_fallback_near_pole_and_exact_poles() {
    let engine = load_de440s_engine();
    // Summer solstice so Sun declination is +23.44° (circumpolar at 68°N!)
    let jd_solstice =
        julday(2024, 6, 20, 20.85, CalendarKind::Gregorian).expect("solstice jd");
    let t_in = TimeInput::from_tt_jd(jd_solstice);

    // 1. Inside Arctic Circle: Tromsø, Norway (69.6492°N, 18.9553°E)
    let tromso = Observer::new(69.6492, 18.9553, 0.0).expect("Tromso");
    let ref_porphyry = engine
        .houses(
            t_in,
            tromso,
            HouseOptions {
                system: HouseSystem::Porphyry,
                ..Default::default()
            },
        )
        .expect("Porphyry at Tromso");

    for sys in HouseSystem::ALL {
        let strict_res = engine.houses(
            t_in,
            tromso,
            HouseOptions {
                system: sys,
                polar_policy: PolarPolicy::Strict,
                ..Default::default()
            },
        );
        let fallback_res = engine.houses(
            t_in,
            tromso,
            HouseOptions {
                system: sys,
                polar_policy: PolarPolicy::FallbackToPorphyry,
                ..Default::default()
            },
        );

        match sys.polar_category() {
            PolarCategory::HardFail => {
                // P, K, i MUST fail under Strict inside polar circle
                assert!(
                    matches!(strict_res, Err(EphemerisError::UndefinedHouseSystem { .. })),
                    "HardFail system '{}' ({sys:?}) must return UndefinedHouseSystem under Strict at 69.65°N",
                    sys.code()
                );
                // Under FallbackToPorphyry, must succeed and match Porphyry
                let fb = fallback_res.expect("FallbackToPorphyry must succeed");
                assert_eq!(fb.requested_house_system, sys);
                assert_eq!(fb.effective_house_system, HouseSystem::Porphyry);
                assert!(fb.fallback_reason.is_some());
                assert!(fb.warnings.contains(&HouseWarning::PolarFallbackToPorphyry));
                for i in 0..12 {
                    assert!(wrap_diff_deg(fb.cusps_deg[i], ref_porphyry.cusps_deg[i]) < 1e-10);
                }
            }
            PolarCategory::SoftFail => {
                let st = strict_res.expect("SoftFail system must succeed at 69.65°N");
                assert!(
                    st.warnings.contains(&HouseWarning::PolarUnreliable),
                    "SoftFail system '{}' ({sys:?}) must emit PolarUnreliable at 69.65°N",
                    sys.code()
                );
            }
            PolarCategory::PolarSafe => {
                let st = strict_res.expect("PolarSafe system must succeed at 69.65°N");
                assert!(!st.warnings.contains(&HouseWarning::PolarUnreliable));
                assert!(!st.warnings.contains(&HouseWarning::PolarFallbackToPorphyry));
            }
            _ => {}
        }
    }

    // 2. Near-pole (+89.9°N): emits NearPoleBoundary warning
    let near_pole = Observer::new(89.9, 0.0, 0.0).expect("89.9N");
    let np_res = engine
        .houses(
            t_in,
            near_pole,
            HouseOptions {
                system: HouseSystem::Porphyry,
                ..Default::default()
            },
        )
        .expect("Porphyry at 89.9N");
    assert!(np_res.warnings.contains(&HouseWarning::NearPoleBoundary));

    // 3. Exact North Pole (+90.0°N): N (WholeSignAries) succeeds, meridian-based systems return UndefinedHouseSystem
    let north_pole = Observer::new(90.0, 0.0, 0.0).expect("90N");
    let aries_pole = engine
        .houses(
            t_in,
            north_pole,
            HouseOptions {
                system: HouseSystem::WholeSignAries,
                ..Default::default()
            },
        )
        .expect("WholeSignAries at 90N");
    assert_eq!(aries_pole.cusps_deg[0], 0.0);

    let porph_pole = engine.houses(
        t_in,
        north_pole,
        HouseOptions {
            system: HouseSystem::Porphyry,
            ..Default::default()
        },
    );
    assert!(matches!(
        porph_pole,
        Err(EphemerisError::UndefinedHouseSystem { .. })
    ));

    // 4. Sidereal Houses & 36 Gauquelin Sectors
    let sid_houses = engine
        .houses(
            t_in,
            london_obs(),
            HouseOptions {
                system: HouseSystem::Placidus,
                zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
                ..Default::default()
            },
        )
        .expect("sidereal Placidus");
    let trop_houses = engine
        .houses(t_in, london_obs(), HouseOptions::default())
        .expect("tropical Placidus");
    let lahiri_val = engine
        .ayanamsa(t_in, AyanamsaModel::Lahiri, NutationMode::True)
        .expect("lahiri");
    assert!(
        wrap_diff_deg(
            (trop_houses.cusps_deg[0] - lahiri_val.effective_deg).rem_euclid(360.0),
            sid_houses.cusps_deg[0]
        ) < 1e-10
    );

    let gauq = engine
        .gauquelin_sectors(t_in, london_obs(), PolarPolicy::Strict)
        .expect("36 Gauquelin sectors");
    assert_eq!(gauq.sectors_deg.len(), 36);
    assert!(wrap_diff_deg(gauq.sectors_deg[0], trop_houses.angles.asc_deg) < 1e-9);
    assert!(wrap_diff_deg(gauq.sectors_deg[9], trop_houses.angles.mc_deg) < 1e-9);
}

fn london_obs() -> Observer {
    Observer::new(51.5074, -0.1276, 0.0).expect("London")
}

// ---------------------------------------------------------------------------
// Part C: Phase 4 — Solar & Lunar Eclipses (NASA Five Millennium Canon)
// ---------------------------------------------------------------------------

#[test]
fn test_solar_eclipses_nasa_canon_and_local_circumstances() {
    let engine = load_de440s_engine();
    let opts = EclipseOptions::default();

    // 1. 1999-08-11 Total Solar Eclipse
    // NASA Five Millennium Canon: Greatest Eclipse = 1999 Aug 11 11:04:09 TD (JD_TT = 2451401.961215),
    // gamma = +0.5062, mag = 1.0286, lat = 45.1°N, lon = 24.3°E (Romania)
    let jd_1999_expected = julday(1999, 8, 11, 11.0 + 4.0 / 60.0 + 9.0 / 3600.0, CalendarKind::Gregorian)
        .expect("1999 jd");
    let bucharest = Observer::new(44.4268, 26.1025, 80.0).expect("Bucharest");
    let sol_1999 = engine
        .solar_eclipse_at(
            TimeInput::from_tt_jd(jd_1999_expected - 0.1),
            Some(bucharest),
            opts,
        )
        .expect("solar_eclipse_at 1999")
        .expect("1999 must be a solar eclipse");

    assert_eq!(sol_1999.kind, SolarEclipseKind::Total);
    let dt_1999_sec = (sol_1999.max_eclipse_jd_tt - jd_1999_expected).abs() * 86_400.0;
    println!(
        "1999-08-11 Total Solar: dt_max = {dt_1999_sec:.2} s, gamma = {:.4}, mag = {:.4}, lat = {:?}, lon = {:?}",
        sol_1999.gamma, sol_1999.magnitude, sol_1999.central_latitude_deg, sol_1999.central_longitude_deg
    );
    assert!(dt_1999_sec < 5.0, "1999 max eclipse time error {dt_1999_sec:.2}s");
    assert!((sol_1999.gamma - 0.5062).abs() < 0.002);
    assert!((sol_1999.magnitude - 1.0286).abs() < 0.002);
    assert!((sol_1999.central_latitude_deg.unwrap_or(0.0) - 45.1).abs() < 0.2);
    assert!((sol_1999.central_longitude_deg.unwrap_or(0.0) - 24.3).abs() < 0.2);
    assert!(sol_1999.p1_jd_tt.is_some() && sol_1999.u1_jd_tt.is_some() && sol_1999.u2_jd_tt.is_some());
    assert!(sol_1999.u3_jd_tt.is_some() && sol_1999.u4_jd_tt.is_some() && sol_1999.p4_jd_tt.is_some());

    // Verify local totality in Bucharest (1999-08-11 passed directly over Bucharest!)
    let loc_1999 = sol_1999.local.expect("local Bucharest eclipse");
    assert!(loc_1999.visible_at_max);
    assert!(loc_1999.is_central_locally);
    assert_eq!(loc_1999.local_kind, SolarEclipseKind::Total);
    assert!(loc_1999.central_duration_sec > 120.0 && loc_1999.central_duration_sec < 150.0);
    assert!(loc_1999.c1_jd_tt.is_some() && loc_1999.c2_jd_tt.is_some());
    assert!(loc_1999.c3_jd_tt.is_some() && loc_1999.c4_jd_tt.is_some());

    // 2. 2017-08-21 Total Solar Eclipse ("Great American Eclipse")
    // NASA Canon: Greatest Eclipse = 2017 Aug 21 18:26:40 TD, gamma = +0.4367, mag = 1.0306
    let jd_2017_expected = julday(2017, 8, 21, 18.0 + 26.0 / 60.0 + 40.0 / 3600.0, CalendarKind::Gregorian)
        .expect("2017 jd");
    let hopkinsville = Observer::new(36.8656, -87.4886, 160.0).expect("Hopkinsville KY");
    let sol_2017 = engine
        .solar_eclipse_at(
            TimeInput::from_tt_jd(jd_2017_expected),
            Some(hopkinsville),
            opts,
        )
        .expect("solar_eclipse_at 2017")
        .expect("2017 eclipse");
    assert_eq!(sol_2017.kind, SolarEclipseKind::Total);
    assert!((sol_2017.max_eclipse_jd_tt - jd_2017_expected).abs() * 86_400.0 < 5.0);
    assert!((sol_2017.gamma - 0.4367).abs() < 0.002);
    let loc_2017 = sol_2017.local.expect("Hopkinsville local");
    assert!(loc_2017.is_central_locally);
    assert!((loc_2017.central_duration_sec - 160.0).abs() < 5.0);

    // 3. 2024-04-08 Total Solar Eclipse
    // NASA Canon: Greatest Eclipse = 2024 Apr 08 18:18:29 TD, gamma = +0.3431, mag = 1.0566
    let jd_2024_expected = julday(2024, 4, 8, 18.0 + 18.0 / 60.0 + 29.0 / 3600.0, CalendarKind::Gregorian)
        .expect("2024 jd");
    let sol_2024 = engine
        .next_solar_eclipse(
            TimeInput::from_tt_jd(jd_2024_expected - 20.0),
            false,
            None,
            opts,
        )
        .expect("next_solar_eclipse 2024");
    assert_eq!(sol_2024.kind, SolarEclipseKind::Total);
    assert!((sol_2024.max_eclipse_jd_tt - jd_2024_expected).abs() * 86_400.0 < 5.0);
    assert!((sol_2024.gamma - 0.3431).abs() < 0.002);

    // 4. 2023-10-14 Annular Solar Eclipse (found via backward search from 2024-04-01)
    let sol_2023 = engine
        .next_solar_eclipse(
            TimeInput::from_tt_jd(jd_2024_expected - 5.0),
            true,
            None,
            opts,
        )
        .expect("backward search to 2023-10-14");
    assert_eq!(sol_2023.kind, SolarEclipseKind::Annular);
    assert!(sol_2023.magnitude < 1.0 && sol_2023.magnitude > 0.94);

    // 5. 2013-11-03 Hybrid Solar Eclipse
    let jd_2013 = julday(2013, 11, 3, 12.0, CalendarKind::Gregorian).expect("2013 jd");
    let sol_2013 = engine
        .solar_eclipse_at(TimeInput::from_tt_jd(jd_2013), None, opts)
        .expect("2013 solar")
        .expect("2013 hybrid");
    assert_eq!(sol_2013.kind, SolarEclipseKind::Hybrid);

    // 6. 2022-10-25 Partial Solar Eclipse
    let jd_2022 = julday(2022, 10, 25, 11.0, CalendarKind::Gregorian).expect("2022 jd");
    let sol_2022 = engine
        .solar_eclipse_at(TimeInput::from_tt_jd(jd_2022), None, opts)
        .expect("2022 solar")
        .expect("2022 partial");
    assert_eq!(sol_2022.kind, SolarEclipseKind::Partial);
    assert!(sol_2022.central_latitude_deg.is_none());
}

#[test]
fn test_lunar_eclipses_nasa_canon_danjon_and_traditional() {
    let engine = load_de440s_engine();
    let opts = EclipseOptions::default();

    // 1. 2003-11-09 Total Lunar Eclipse
    // NASA Five Millennium Canon: Greatest Eclipse = 2003 Nov 09 01:19:38 TD, gamma = -0.4319, umbral mag = 1.022
    let jd_2003_expected = julday(2003, 11, 9, 1.0 + 19.0 / 60.0 + 38.0 / 3600.0, CalendarKind::Gregorian)
        .expect("2003 jd");
    let lun_2003 = engine
        .lunar_eclipse_at(TimeInput::from_tt_jd(jd_2003_expected), opts)
        .expect("lunar 2003")
        .expect("2003 total lunar");
    assert_eq!(lun_2003.kind, LunarEclipseKind::Total);
    let dt_2003_sec = (lun_2003.max_eclipse_jd_tt - jd_2003_expected).abs() * 86_400.0;
    println!(
        "2003-11-09 Total Lunar: dt_max = {dt_2003_sec:.2} s, gamma = {:.4}, umb_mag = {:.4}, pen_mag = {:.4}",
        lun_2003.gamma, lun_2003.umbral_magnitude, lun_2003.penumbral_magnitude
    );
    assert!(dt_2003_sec < 15.0, "2003 lunar max time error {dt_2003_sec:.2}s");
    assert!((lun_2003.gamma - (-0.4319)).abs() < 0.005);
    assert!((lun_2003.umbral_magnitude - 1.022).abs() < 0.01);
    assert!(lun_2003.p1_jd_tt.is_some() && lun_2003.u1_jd_tt.is_some() && lun_2003.u2_jd_tt.is_some());
    assert!(lun_2003.u3_jd_tt.is_some() && lun_2003.u4_jd_tt.is_some() && lun_2003.p4_jd_tt.is_some());

    // Compare Danjon (1/85) vs Traditional (2%) shadow enlargement rules
    let lun_2003_trad = engine
        .lunar_eclipse_at(
            TimeInput::from_tt_jd(jd_2003_expected),
            EclipseOptions {
                shadow_rule: ShadowEnlargementRule::Traditional2Percent,
                ..Default::default()
            },
        )
        .expect("lunar 2003 trad")
        .expect("2003 trad");
    assert!(lun_2003_trad.umbral_magnitude > lun_2003.umbral_magnitude);

    // 2. 2025-03-14 Total Lunar Eclipse (forward search from 2025-03-01)
    let jd_2025_start = julday(2025, 3, 1, 0.0, CalendarKind::Gregorian).expect("2025 start");
    let jd_2025_expected = julday(2025, 3, 14, 6.0 + 59.0 / 60.0 + 56.0 / 3600.0, CalendarKind::Gregorian)
        .expect("2025 jd");
    let lun_2025 = engine
        .next_lunar_eclipse(TimeInput::from_tt_jd(jd_2025_start), false, opts)
        .expect("next_lunar_eclipse 2025");
    assert_eq!(lun_2025.kind, LunarEclipseKind::Total);
    assert!((lun_2025.max_eclipse_jd_tt - jd_2025_expected).abs() * 86_400.0 < 15.0);
    assert!((lun_2025.gamma - 0.348).abs() < 0.005);

    // 3. 2024-03-25 Penumbral Lunar Eclipse
    let jd_2024_pen = julday(2024, 3, 25, 7.0, CalendarKind::Gregorian).expect("2024 pen");
    let lun_2024_pen = engine
        .lunar_eclipse_at(TimeInput::from_tt_jd(jd_2024_pen), opts)
        .expect("2024 pen")
        .expect("2024 penumbral");
    assert_eq!(lun_2024_pen.kind, LunarEclipseKind::Penumbral);
    assert!(lun_2024_pen.penumbral_magnitude > 0.90);
    assert!(lun_2024_pen.umbral_magnitude < 0.0);
    assert!(lun_2024_pen.u1_jd_tt.is_none() && lun_2024_pen.u2_jd_tt.is_none());
}
