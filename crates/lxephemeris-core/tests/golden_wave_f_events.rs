//! Golden verification suite for Wave F Astronomical & Astrological Event Solvers (`events.rs`):
//! 1. 2026 Chinese 24 Solar Terms & Winter Solstice (`R12` regression test against `ephem.Sun().hlong` $180^\circ$ bug)
//! 2. Human Design `-88°` Solar Arc prenatal epoch (`R9.2`)
//! 3. Planetary retrograde triple-pass longitude returns & Tropical/Sidereal (`Lahiri`) stations (`R9.1`, `R9.2`)
//! 4. 2025/2026 Major Lunar Standstill declination extrema / lunistices (`R9.2`)
//! 5. Principal lunar phases (`NewMoon`, `FirstQuarter`, `FullMoon`, `LastQuarter`) and `lunar_phase_summary` (`R0`, `R2.2`, `R9.2`)
//! 6. Bernadette Brady diurnal fixed-star parans (`R6.2`, Subject D London test + circumpolar stars + orb monotonicity)

#![forbid(unsafe_code)]

use lxephemeris_core::{
    julday_split, revjul, wrap_pm180, AyanamsaModel, BodyId, CalcOptions, CalendarKind,
    DeclinationExtremumKind, DiscMode, Engine, EngineConfig, Equinox, LunarPhaseKind, NutationMode,
    Observer, Origin, ParanAngle, ParanOptions, Plane, PlanetId, RiseTransOptions, StationKind,
    TimeInput, Zodiac,
};
use std::path::PathBuf;

fn load_spk_engine() -> Engine {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .map_or_else(|| PathBuf::from("data"), |p| p.join("data"));
    Engine::from_data_dir(data_dir, EngineConfig::default()).expect("load SPK engine from data/")
}

#[test]
fn test_r7_1_and_r6_2_rise_trans_options_presets() {
    let almanac = RiseTransOptions::almanac_upper_limb();
    assert_eq!(almanac.max_search_days, 2.0);
    assert_eq!(almanac.custom_horizon_alt_deg, None);
    assert_eq!(almanac.disc_mode, DiscMode::UpperLimb);
    assert!((almanac.pressure_hpa - 1013.25).abs() < 1e-12);
    assert!((almanac.temperature_c - 15.0).abs() < 1e-12);

    let geom = RiseTransOptions::geometric_center();
    assert_eq!(geom.max_search_days, 2.0);
    assert_eq!(geom.custom_horizon_alt_deg, None);
    assert_eq!(geom.disc_mode, DiscMode::Center);
    assert!((geom.pressure_hpa - 0.0).abs() < 1e-12);
    assert!((geom.temperature_c - 15.0).abs() < 1e-12);
}

#[test]
fn test_r12_chinese_24_solar_terms_and_winter_solstice_2026() {
    let spk_engine = load_spk_engine();
    let moshier_engine = Engine::new_moshier(EngineConfig::default());

    let start_2026 =
        TimeInput::Utc(julday_split(2026, 1, 1, 0.0, CalendarKind::Gregorian).expect("start 2026"));
    let end_2026 = TimeInput::Utc(
        julday_split(2026, 12, 31, 23.99, CalendarKind::Gregorian).expect("end 2026"),
    );

    for (label, engine) in [("SPK-DE440s", &spk_engine), ("Moshier", &moshier_engine)] {
        let terms = engine
            .solar_terms_in_range(start_2026, end_2026)
            .expect("2026 solar terms");
        assert_eq!(
            terms.len(),
            24,
            "{label}: calendar year 2026 must contain all 24 Chinese solar terms"
        );

        // Verify strict chronological ordering and exact geocentric apparent tropical solar longitude
        for pair in terms.windows(2) {
            if let [a, b] = pair {
                assert!(
                    b.jd_tt > a.jd_tt + 14.0 && b.jd_tt < a.jd_tt + 16.5,
                    "{label}: consecutive solar terms must be ~15.2 days apart: {a:?} vs {b:?}"
                );
            }
        }

        for t in &terms {
            let sun_res = engine
                .calc(
                    BodyId::Planet(PlanetId::Sun),
                    TimeInput::from_tt_jd(t.jd_tt),
                    CalcOptions::default(),
                )
                .expect("eval Sun at solar term");
            let diff = wrap_pm180(sun_res.lon_deg - t.target_lon_deg).abs();
            assert!(
                diff < 1.0e-7,
                "{label}: Solar term {} ({}) target={}, actual={}, diff={diff:e}",
                t.pinyin_name,
                t.term_index,
                t.target_lon_deg,
                sun_res.lon_deg
            );
            assert_eq!(t.is_major_term, t.term_index % 2 == 0);
        }

        // First solar term of Gregorian 2026 is Xiaohan (285°, around Jan 5, 2026)
        let first_term = terms.first().expect("at least one solar term");
        assert_eq!(first_term.pinyin_name, "Xiaohan");
        assert_eq!(first_term.term_index, 19);
        assert!((first_term.target_lon_deg - 285.0).abs() < 1e-12);

        // Find Dongzhi (Winter Solstice, 270°, term_index == 18): must occur on 2026-12-21 ~20:50 UTC
        // (never in June! Permanently prevents the ephem.Sun().hlong 180° heliocentric inversion bug R12)
        let dongzhi = terms
            .iter()
            .find(|t| t.pinyin_name == "Dongzhi")
            .expect("Dongzhi in 2026");
        assert_eq!(dongzhi.term_index, 18);
        assert!((dongzhi.target_lon_deg - 270.0).abs() < 1e-12);
        assert!(dongzhi.is_major_term);

        let cal = revjul(dongzhi.jd_ut1, CalendarKind::Gregorian).expect("revjul dongzhi");
        assert_eq!(cal.year, 2026);
        assert_eq!(cal.month, 12, "{label}: Dongzhi must be in December, not June");
        assert_eq!(cal.day, 21);
        assert!(
            (cal.hour - (20.0 + 50.0 / 60.0)).abs() < 0.05,
            "{label}: 2026 Dongzhi expected around 20:50 UT, got hour={}",
            cal.hour
        );

        // Cross-check against seasons_in_year(2026)
        let seasons = engine.seasons_in_year(2026).expect("seasons 2026");
        let chunfen = terms
            .iter()
            .find(|t| t.term_index == 0)
            .expect("Chunfen 2026");
        let xiazhi = terms
            .iter()
            .find(|t| t.term_index == 6)
            .expect("Xiazhi 2026");
        let qiufen = terms
            .iter()
            .find(|t| t.term_index == 12)
            .expect("Qiufen 2026");

        assert!((seasons.vernal_equinox_jd_tt - chunfen.jd_tt).abs() < 1.0e-8);
        assert!((seasons.summer_solstice_jd_tt - xiazhi.jd_tt).abs() < 1.0e-8);
        assert!((seasons.autumnal_equinox_jd_tt - qiufen.jd_tt).abs() < 1.0e-8);
        assert!((seasons.winter_solstice_jd_tt - dongzhi.jd_tt).abs() < 1.0e-8);
    }
}

#[test]
fn test_r9_2_human_design_88_deg_solar_arc_epoch() {
    let spk_engine = load_spk_engine();
    let moshier_engine = Engine::new_moshier(EngineConfig::default());

    let birth = TimeInput::Utc(
        julday_split(1990, 6, 15, 14.5, CalendarKind::Gregorian).expect("birth epoch"),
    );

    for (label, engine) in [("SPK-DE440s", &spk_engine), ("Moshier", &moshier_engine)] {
        let sun_birth = engine
            .calc(
                BodyId::Planet(PlanetId::Sun),
                birth,
                CalcOptions::default(),
            )
            .expect("sun at birth");
        let hd_event = engine
            .human_design_epoch(birth)
            .expect("human design epoch");

        let dt_days = sun_birth.jd_tt - hd_event.jd_tt;
        assert!(
            (87.0..=91.0).contains(&dt_days),
            "{label}: Human Design epoch must be ~88-90 days before birth, got {dt_days:.4} d"
        );

        let arc_deg = (sun_birth.lon_deg - hd_event.actual_lon_deg).rem_euclid(360.0);
        assert!(
            (arc_deg - 88.0).abs() < 1.0e-7,
            "{label}: Solar arc must be 88.0° +/- 1e-7°, got {arc_deg:.10}°"
        );
        assert!(!hd_event.is_retrograde);
    }
}

#[test]
fn test_r9_2_retrograde_triple_return_and_stations_tropical_and_sidereal() {
    let engine = load_spk_engine();

    // 1. Find Saturn's retrograde stations in 2026
    let start_2026 =
        TimeInput::Utc(julday_split(2026, 1, 1, 0.0, CalendarKind::Gregorian).expect("start 2026"));
    let end_2026 = TimeInput::Utc(
        julday_split(2026, 12, 31, 0.0, CalendarKind::Gregorian).expect("end 2026"),
    );
    let saturn = BodyId::Planet(PlanetId::Saturn);

    let sat_stations = engine
        .search_stations_in_range(saturn, start_2026, end_2026, &CalcOptions::default())
        .expect("Saturn stations 2026");
    assert_eq!(
        sat_stations.len(),
        2,
        "Saturn must have 1 retrograde station and 1 direct station in 2026"
    );
    let st0 = sat_stations.first().expect("station 0");
    let st1 = sat_stations.get(1).expect("station 1");
    assert_eq!(st0.kind, StationKind::StationRetrograde);
    assert_eq!(st1.kind, StationKind::StationDirect);
    assert!(st0.jd_tt < st1.jd_tt);

    // Midpoint longitude between the two stations is crossed 3 times in 2026 (direct -> retro -> direct)!
    let mid_lon_deg = 0.5 * (st0.lon_deg + st1.lon_deg);
    let start_2027 = TimeInput::Utc(
        julday_split(2027, 4, 1, 0.0, CalendarKind::Gregorian).expect("end triple window"),
    );
    let crossings = engine
        .body_longitude_crossings_in_range(
            saturn,
            start_2026,
            start_2027,
            mid_lon_deg,
            &CalcOptions::default(),
        )
        .expect("Saturn triple return crossings");
    assert_eq!(
        crossings.len(),
        3,
        "Target longitude inside Saturn's retrograde loop must be crossed 3 times: {crossings:?}"
    );
    let c0 = crossings.first().expect("pass 1");
    let c1 = crossings.get(1).expect("pass 2");
    let c2 = crossings.get(2).expect("pass 3");
    assert!(!c0.is_retrograde, "Pass 1 must be direct");
    assert!(c1.is_retrograde, "Pass 2 must be retrograde");
    assert!(!c2.is_retrograde, "Pass 3 must be direct");
    for c in &crossings {
        assert!(wrap_pm180(c.actual_lon_deg - mid_lon_deg).abs() < 1.0e-7);
    }

    // 2. Verify Mercury, Mars, and Saturn stations in both Tropical and Sidereal (Lahiri) modes
    let sid_opts = CalcOptions {
        zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
        ..CalcOptions::default()
    };
    let end_2027 = TimeInput::Utc(
        julday_split(2027, 6, 1, 0.0, CalendarKind::Gregorian).expect("2027-06-01"),
    );

    for body in [
        BodyId::Planet(PlanetId::Mercury),
        BodyId::Planet(PlanetId::Mars),
        BodyId::Planet(PlanetId::Saturn),
    ] {
        for (mode_label, opts) in [
            ("Tropical", CalcOptions::default()),
            ("Sidereal-Lahiri", sid_opts),
        ] {
            let stations = engine
                .search_stations_in_range(body, start_2026, end_2027, &opts)
                .expect("stations");
            assert!(
                stations.len() >= 2,
                "{mode_label} {body:?} must have at least 2 stations in 2026..2027-06"
            );
            for st in &stations {
                let check = engine
                    .calc(body, TimeInput::from_tt_jd(st.jd_tt), opts)
                    .expect("eval at station");
                assert!(
                    check.lon_speed_deg_day.abs() < 1.0e-6,
                    "{mode_label} {body:?} station residual speed must be < 1e-6 deg/day, got {:e}",
                    check.lon_speed_deg_day
                );
            }
        }
    }

    // Sun, Moon, and Mean_Node have no retrograde stations
    for monotonic_body in [
        BodyId::Planet(PlanetId::Sun),
        BodyId::Planet(PlanetId::Moon),
    ] {
        let empty = engine
            .search_stations_in_range(
                monotonic_body,
                start_2026,
                end_2026,
                &CalcOptions::default(),
            )
            .expect("monotonic body stations");
        assert!(empty.is_empty());
    }
}

#[test]
fn test_r9_2_lunar_declination_standstills_2025() {
    let engine = load_spk_engine();
    let start = TimeInput::Utc(
        julday_split(2025, 3, 1, 0.0, CalendarKind::Gregorian).expect("2025-03-01"),
    );
    let end = TimeInput::Utc(
        julday_split(2025, 3, 31, 0.0, CalendarKind::Gregorian).expect("2025-03-31"),
    );

    let extrema = engine
        .search_declination_extrema_in_range(
            BodyId::Planet(PlanetId::Moon),
            start,
            end,
            &CalcOptions::default(),
        )
        .expect("Moon declination extrema March 2025");

    assert!(
        extrema.len() >= 2,
        "30 days covers >2 draconic half-months (~13.6d each), expected >= 2 extrema, got {}",
        extrema.len()
    );

    let eq_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Equatorial,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: lxephemeris_core::AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: true,
    };

    for pair in extrema.windows(2) {
        if let [a, b] = pair {
            assert_ne!(
                a.kind, b.kind,
                "NorthernMaximum and SouthernMinimum must alternate"
            );
            let dt = b.jd_tt - a.jd_tt;
            assert!(
                (12.0..=15.5).contains(&dt),
                "Lunar standstill extrema must be ~13.6 days apart, got {dt:.3} d"
            );
        }
    }

    for ex in &extrema {
        // 2025 is the Major Lunar Standstill peak: |dec_deg| exceeds 28°!
        assert!(
            ex.dec_deg.abs() > 28.0,
            "2025 major lunar standstill extremum must exceed 28° declination, got {}",
            ex.dec_deg
        );
        match ex.kind {
            DeclinationExtremumKind::NorthernMaximum => assert!(ex.dec_deg > 28.0),
            DeclinationExtremumKind::SouthernMinimum => assert!(ex.dec_deg < -28.0),
        }

        let check = engine
            .calc(
                BodyId::Planet(PlanetId::Moon),
                TimeInput::from_tt_jd(ex.jd_tt),
                eq_opts,
            )
            .expect("eval Moon equatorial at extremum");
        assert!(
            check.lat_speed_deg_day.abs() < 1.0e-5,
            "Declination rate at extremum must be < 1e-5 deg/day, got {:e}",
            check.lat_speed_deg_day
        );
    }
}

#[test]
fn test_r0_r2_2_r9_2_lunar_phases_and_kerykeion_summary() {
    let spk_engine = load_spk_engine();
    let moshier_engine = Engine::new_moshier(EngineConfig::default());

    let start = TimeInput::Utc(
        julday_split(2026, 1, 1, 0.0, CalendarKind::Gregorian).expect("2026-01-01"),
    );
    let end = TimeInput::Utc(
        julday_split(2026, 2, 15, 0.0, CalendarKind::Gregorian).expect("2026-02-15"),
    );

    for (label, engine) in [("SPK-DE440s", &spk_engine), ("Moshier", &moshier_engine)] {
        let phases = engine
            .lunar_phases_in_range(start, end)
            .expect("lunar phases in range");
        assert!(
            phases.len() >= 6,
            "{label}: 45 days must contain at least 6 principal lunar phases, got {}",
            phases.len()
        );

        for p in &phases {
            assert!(p.kind.is_some());
            let actual_elong = (p.moon_lon_deg - p.sun_lon_deg).rem_euclid(360.0);
            let diff = wrap_pm180(actual_elong - p.target_elongation_deg).abs();
            assert!(
                diff < 1.0e-7,
                "{label}: Lunar phase {:?} elongation error {diff:e}",
                p.kind
            );
            assert!(
                (356_000.0..=406_800.0).contains(&p.moon_dist_km),
                "{label}: Moon distance km out of physical bounds: {}",
                p.moon_dist_km
            );
        }

        // Verify lunar_phase_summary at an intermediate epoch (2026-01-25 12:00 UTC)
        let probe = TimeInput::Utc(
            julday_split(2026, 1, 25, 12.0, CalendarKind::Gregorian).expect("2026-01-25"),
        );
        let summary = engine
            .lunar_phase_summary(probe)
            .expect("lunar_phase_summary");

        assert!((0.0..360.0).contains(&summary.elongation_deg));
        assert!((0.0..=1.0).contains(&summary.illumination_fraction));
        assert!((0.0..=180.0).contains(&summary.phase_angle_deg));
        assert!((356_000.0..=406_800.0).contains(&summary.moon_dist_km));
        assert!(summary.age_days > 0.0 && summary.age_days < 30.0);
        assert!(summary.synodic_fraction > 0.0 && summary.synodic_fraction < 1.0);
        assert!(summary.phase_8_index < 8);
        assert!(!summary.phase_8_name.is_empty());
        assert!((1..=28).contains(&summary.phase_28_index));
        assert!(summary.prev_new_moon_jd_ut1 < summary.jd_ut1);
        assert!(summary.next_new_moon_jd_ut1 > summary.jd_ut1);
        assert!(summary.next_first_quarter_jd_ut1 > summary.jd_ut1);
        assert!(summary.next_full_moon_jd_ut1 > summary.jd_ut1);
        assert!(summary.next_last_quarter_jd_ut1 > summary.jd_ut1);

        // Verify LunarPhaseKind mapping
        assert_eq!(
            LunarPhaseKind::from_elongation_deg(0.0),
            Some(LunarPhaseKind::NewMoon)
        );
        assert_eq!(
            LunarPhaseKind::from_elongation_deg(90.0),
            Some(LunarPhaseKind::FirstQuarter)
        );
        assert_eq!(
            LunarPhaseKind::from_elongation_deg(180.0),
            Some(LunarPhaseKind::FullMoon)
        );
        assert_eq!(
            LunarPhaseKind::from_elongation_deg(270.0),
            Some(LunarPhaseKind::LastQuarter)
        );
    }
}

#[test]
fn test_r6_2_brady_diurnal_fixed_star_parans_subject_d_london() {
    let spk_engine = load_spk_engine();
    let moshier_engine = Engine::new_moshier(EngineConfig::default());

    // Subject D: London (51.5074°N, -0.1278°E), 1969-05-17 14:00 UT
    let london = Observer::new(51.5074, -0.1278, 0.0).expect("London observer");
    let epoch = TimeInput::Ut1(
        julday_split(1969, 5, 17, 14.0, CalendarKind::Gregorian).expect("1969-05-17 14:00 UT"),
    );

    assert_eq!(ParanAngle::Rising.as_str(), "rising");
    assert_eq!(ParanAngle::Culminating.as_str(), "culminating");
    assert_eq!(ParanAngle::Setting.as_str(), "setting");
    assert_eq!(ParanAngle::AntiCulminating.as_str(), "anti-culminating");

    for (label, engine) in [("SPK-DE440s", &spk_engine), ("Moshier", &moshier_engine)] {
        let opts_4m = ParanOptions {
            orb_minutes: 4.0,
            ..ParanOptions::default()
        };
        let parans_4m = engine
            .find_parans(epoch, london, &opts_4m)
            .expect("find_parans 4.0m");

        assert!(
            !parans_4m.is_empty(),
            "{label}: Subject D (London 1969-05-17) must have at least 1 paran within 4.0m"
        );

        for pair in parans_4m.windows(2) {
            if let [a, b] = pair {
                assert!(
                    a.time_orb_minutes <= b.time_orb_minutes,
                    "{label}: Parans must be sorted ascending by time_orb_minutes: {} vs {}",
                    a.time_orb_minutes,
                    b.time_orb_minutes
                );
            }
        }

        for p in &parans_4m {
            assert!(p.time_orb_minutes >= 0.0 && p.time_orb_minutes <= 4.0);
            assert_eq!(
                p.paran_type, p.star_angle,
                "{label}: paran_type must always equal star_angle"
            );
            assert!(!p.star_name.is_empty());
            assert!(!p.star_bayer.is_empty());
            assert!(p.star_hip > 0);
            assert!(!p.planet_name.is_empty());
        }

        // Verify monotonicity of orb counts: 1.0 min <= 4.0 min <= 10.0 min
        let parans_1m = engine
            .find_parans(
                epoch,
                london,
                &ParanOptions {
                    orb_minutes: 1.0,
                    ..ParanOptions::default()
                },
            )
            .expect("find_parans 1.0m");
        let parans_10m = engine
            .find_parans(
                epoch,
                london,
                &ParanOptions {
                    orb_minutes: 10.0,
                    ..ParanOptions::default()
                },
            )
            .expect("find_parans 10.0m");

        assert!(
            parans_1m.len() <= parans_4m.len() && parans_4m.len() <= parans_10m.len(),
            "{label}: Paran count must be monotonic with orb: 1m={}, 4m={}, 10m={}",
            parans_1m.len(),
            parans_4m.len(),
            parans_10m.len()
        );

        // Verify circumpolar stars at 51.5°N (Dubhe, Polaris, Capella) never rise/set
        // but still culminate and anti-culminate cleanly without error!
        let circumpolar_parans = engine
            .find_parans(
                epoch,
                london,
                &ParanOptions {
                    orb_minutes: 720.0,
                    bodies: vec![BodyId::Planet(PlanetId::Sun)],
                    star_names: Some(vec![
                        "Dubhe".to_owned(),
                        "Polaris".to_owned(),
                        "Capella".to_owned(),
                    ]),
                },
            )
            .expect("circumpolar stars paran scan");

        // Each circumpolar star has 2 meridian transits (Culminating, AntiCulminating),
        // and the Sun has 4 diurnal events -> 3 * 2 * 4 = 24 pairs across the full 720m half-day!
        assert_eq!(
            circumpolar_parans.len(),
            24,
            "{label}: 3 circumpolar stars (2 transits each) x 4 Sun events must yield 24 pairs"
        );
        for cp in &circumpolar_parans {
            assert!(
                matches!(
                    cp.star_angle,
                    ParanAngle::Culminating | ParanAngle::AntiCulminating
                ),
                "{label}: Circumpolar star {} at 51.5°N cannot have Rising/Setting: {:?}",
                cp.star_name,
                cp.star_angle
            );
        }
    }
}
