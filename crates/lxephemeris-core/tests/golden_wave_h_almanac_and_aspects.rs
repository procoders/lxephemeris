//! Wave H Golden Integration Test Suite:
//! Verifies Zero-Threshold Initialization (`Engine::default()`, `Engine::from_env_or_moshier`),
//! Planetary Returns (`planetary_return`, `planetary_returns_in_range` including tropical/sidereal
//! Solar Returns and 3-pass retrograde Saturn Returns), Aspects (`calc_chart_aspects`,
//! `calc_synastry_aspects`, `calc_declination_aspects`, `search_exact_aspect`,
//! `mutual_aspects_in_range`, `natal_transits_in_range`, `void_of_course_moon_in_range`),
//! Vedic Jyotish (`panchanga`, `vimshottari_dasha`, `calc_varga`, `calc_shodashavarga`),
//! Traditional Western & Human Design (`chart_sect`, `planetary_hours`, `calc_lots`,
//! `map_zodiac_placement`), and one-call `chart_snapshot`.

#![allow(
    clippy::indexing_slicing,
    clippy::expect_used,
    clippy::unwrap_used,
    clippy::panic,
    dead_code
)]

use lxephemeris_core::{
    julday, map_zodiac_placement, vimshottari_dasha_from_lon, wrap_pm180, AspectPreset,
    AyanamsaModel, BodyId, CalcOptions, CalendarKind, DeclinationAspectKind, Engine, EngineConfig,
    HouseOptions, HouseSystem, Observer, PlanetId, TimeInput, Zodiac, HD_MANDALA_GATES_64,
    JD_J2000,
};
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Arc;

fn workspace_data_dir() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
        .join("data")
}

#[test]
fn test_1_zero_threshold_engine_default_and_from_env_or_moshier() {
    // 1. Engine::default() works out-of-the-box with zero external files required
    let default_engine = Engine::default();
    let sun_default = default_engine
        .calc(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_tt_jd(JD_J2000),
            CalcOptions::default(),
        )
        .expect("Engine::default() calc Sun");
    assert!(
        (sun_default.lon_deg - 280.368_165_28).abs() < 1.0e-3,
        "Default engine Sun lon_deg={}",
        sun_default.lon_deg
    );

    // 2. Engine::from_env_or_moshier initializes cleanly
    let env_engine = Engine::from_env_or_moshier(EngineConfig::default());
    let moon_env = env_engine
        .calc(
            BodyId::Planet(PlanetId::Moon),
            TimeInput::from_tt_jd(JD_J2000),
            CalcOptions::default(),
        )
        .expect("Engine::from_env_or_moshier calc Moon");
    assert!((0.0..360.0).contains(&moon_env.lon_deg));

    // 3. Engine::from_data_dir_or_moshier falls back gracefully on non-existent path
    let fallback_engine =
        Engine::from_data_dir_or_moshier("/nonexistent/lxephemeris/path", EngineConfig::default())
            .expect("from_data_dir_or_moshier fallback");
    assert!(fallback_engine.is_moshier_mode());
}

#[test]
fn test_2_planetary_returns_tropical_sidereal_and_saturn_triple_pass() {
    let engine = Engine::from_data_dir_or_moshier(workspace_data_dir(), EngineConfig::default())
        .expect("load engine");

    let natal_jd = julday(1990, 6, 21, 12.0, CalendarKind::Gregorian).expect("1990-06-21");
    let search_after_2024 = julday(2024, 1, 1, 0.0, CalendarKind::Gregorian).expect("2024-01-01");

    // 2a. Tropical Solar Return
    let trop_opts = CalcOptions::default();
    let natal_sun_trop = engine
        .calc(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_ut1_jd(natal_jd),
            trop_opts,
        )
        .expect("natal Sun tropical");
    let sr_trop = engine
        .planetary_return(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_ut1_jd(natal_jd),
            TimeInput::from_ut1_jd(search_after_2024),
            trop_opts,
        )
        .expect("tropical solar return");
    assert_eq!(sr_trop.body, BodyId::Planet(PlanetId::Sun));
    assert_eq!(sr_trop.pass_index, 1);
    assert!(!sr_trop.is_retrograde);
    assert!(
        wrap_pm180(sr_trop.actual_lon_deg - natal_sun_trop.lon_deg).abs() < 1.0e-7,
        "Tropical SR residual too large: {} vs {}",
        sr_trop.actual_lon_deg,
        natal_sun_trop.lon_deg
    );

    // 2b. Sidereal (Lahiri) Solar Return — occurs ~20 minutes later than Tropical SR due to precession over 34 years
    let sid_opts = CalcOptions {
        zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
        ..CalcOptions::default()
    };
    let natal_sun_sid = engine
        .calc(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_ut1_jd(natal_jd),
            sid_opts,
        )
        .expect("natal Sun sidereal");
    let sr_sid = engine
        .planetary_return(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_ut1_jd(natal_jd),
            TimeInput::from_ut1_jd(search_after_2024),
            sid_opts,
        )
        .expect("sidereal solar return");
    assert!(
        wrap_pm180(sr_sid.actual_lon_deg - natal_sun_sid.lon_deg).abs() < 1.0e-7,
        "Sidereal SR residual too large"
    );
    // Over 34 years, precession is ~0.47°, which takes the Sun ~0.48 days to traverse
    let sid_minus_trop_days = sr_sid.jd_ut1 - sr_trop.jd_ut1;
    assert!(
        sid_minus_trop_days > 0.40 && sid_minus_trop_days < 0.60,
        "Expected Sidereal SR ~0.48 days after Tropical SR over 34 years, got {sid_minus_trop_days}"
    );

    // 2c. Saturn triple-pass retrograde return (pass_index = 1, 2, 3)
    // In 2024-2025, Saturn retrogrades from ~19.4° Pisces (349.4°) back to ~12.7° Pisces (342.7°).
    // Find when Saturn was at 345.0° (15° Pisces) in 1995:
    let jd_1995 = julday(1995, 1, 1, 0.0, CalendarKind::Gregorian).expect("1995-01-01");
    let saturn_1995 = engine
        .search_body_longitude(
            BodyId::Planet(PlanetId::Saturn),
            TimeInput::from_ut1_jd(jd_1995),
            345.0,
            true,
            400.0,
            &trop_opts,
        )
        .expect("Saturn at 345° in 1995");
    let jd_2024_start = julday(2024, 1, 1, 0.0, CalendarKind::Gregorian).expect("2024-01-01");
    let jd_2025_mid = julday(2025, 5, 1, 0.0, CalendarKind::Gregorian).expect("2025-05-01");

    let saturn_returns = engine
        .planetary_returns_in_range(
            BodyId::Planet(PlanetId::Saturn),
            TimeInput::from_ut1_jd(saturn_1995.jd_ut1),
            TimeInput::from_ut1_jd(jd_2024_start),
            TimeInput::from_ut1_jd(jd_2025_mid),
            trop_opts,
        )
        .expect("Saturn returns in 2024-2025");
    assert_eq!(
        saturn_returns.len(),
        3,
        "Expected 3-pass Saturn return in 2024-2025, got {}",
        saturn_returns.len()
    );
    assert_eq!(saturn_returns[0].pass_index, 1);
    assert_eq!(saturn_returns[1].pass_index, 2);
    assert_eq!(saturn_returns[2].pass_index, 3);
    assert!(!saturn_returns[0].is_retrograde, "Pass 1 must be direct");
    assert!(saturn_returns[1].is_retrograde, "Pass 2 must be retrograde");
    assert!(!saturn_returns[2].is_retrograde, "Pass 3 must be direct");
    for sr in &saturn_returns {
        assert!(wrap_pm180(sr.actual_lon_deg - 345.0).abs() < 1.0e-6);
    }
}

#[test]
fn test_3_aspects_synastry_declination_transits_and_voc_moon() {
    let engine = Engine::from_data_dir_or_moshier(workspace_data_dir(), EngineConfig::default())
        .expect("load engine");

    let bodies = [
        BodyId::Planet(PlanetId::Sun),
        BodyId::Planet(PlanetId::Moon),
        BodyId::Planet(PlanetId::Mercury),
        BodyId::Planet(PlanetId::Venus),
        BodyId::Planet(PlanetId::Mars),
        BodyId::Planet(PlanetId::Jupiter),
        BodyId::Planet(PlanetId::Saturn),
        BodyId::Planet(PlanetId::Uranus),
        BodyId::Planet(PlanetId::Neptune),
        BodyId::Planet(PlanetId::Pluto),
    ];

    // 3a. calc_chart_aspects
    let chart_hits = engine
        .calc_chart_aspects(
            TimeInput::from_ut1_jd(JD_J2000),
            &bodies,
            AspectPreset::All,
            1.0,
            CalcOptions::default(),
        )
        .expect("calc_chart_aspects");
    assert!(!chart_hits.is_empty());
    for hit in &chart_hits {
        assert!(hit.orb_deg <= hit.max_orb_deg + 1.0e-12);
        assert!((0.0..=100.0).contains(&hit.tightness_pct));
    }

    // 3b. calc_synastry_aspects
    let syn_hits = engine
        .calc_synastry_aspects(
            TimeInput::from_ut1_jd(JD_J2000),
            &bodies[..5],
            TimeInput::from_ut1_jd(JD_J2000 + 90.0),
            &bodies[..5],
            AspectPreset::Standard,
            1.0,
            CalcOptions::default(),
        )
        .expect("calc_synastry_aspects");
    assert!(!syn_hits.is_empty());

    // 3c. calc_declination_aspects
    let dec_hits = engine
        .calc_declination_aspects(
            TimeInput::from_ut1_jd(JD_J2000),
            &bodies,
            2.0,
            CalcOptions::default(),
        )
        .expect("calc_declination_aspects");
    assert!(!dec_hits.is_empty());
    for dh in &dec_hits {
        assert!(dh.orb_deg <= 2.0 + 1.0e-12);
        match dh.kind {
            DeclinationAspectKind::Parallel => {
                assert!(((dh.dec1_deg - dh.dec2_deg).abs() - dh.orb_deg).abs() < 1.0e-12);
            }
            DeclinationAspectKind::Contraparallel => {
                assert!(((dh.dec1_deg + dh.dec2_deg).abs() - dh.orb_deg).abs() < 1.0e-12);
            }
        }
    }

    // 3d. search_exact_aspect (Moon-Sun 120° Trine)
    let jd_2024 = julday(2024, 1, 1, 0.0, CalendarKind::Gregorian).expect("2024-01-01");
    let exact_trine = engine
        .search_exact_aspect(
            BodyId::Planet(PlanetId::Moon),
            BodyId::Planet(PlanetId::Sun),
            120.0,
            TimeInput::from_ut1_jd(jd_2024),
            true,
            CalcOptions::default(),
        )
        .expect("search_exact_aspect Moon-Sun 120");
    assert_eq!(exact_trine.aspect_name, "Trine");
    let sep = wrap_pm180(exact_trine.body1_lon_deg - exact_trine.body2_lon_deg).abs();
    assert!(
        (sep - 120.0).abs() < 1.0e-6,
        "Expected exact 120° separation, got {sep}"
    );

    // 3e. mutual_aspects_in_range
    let mutual = engine
        .mutual_aspects_in_range(
            &[
                BodyId::Planet(PlanetId::Sun),
                BodyId::Planet(PlanetId::Moon),
                BodyId::Planet(PlanetId::Mercury),
            ],
            AspectPreset::Major,
            TimeInput::from_ut1_jd(jd_2024),
            TimeInput::from_ut1_jd(jd_2024 + 10.0),
            CalcOptions::default(),
        )
        .expect("mutual_aspects_in_range");
    assert!(!mutual.is_empty());
    for w in mutual.windows(2) {
        assert!(w[0].jd_tt <= w[1].jd_tt);
    }

    // 3f. natal_transits_in_range
    let transits = engine
        .natal_transits_in_range(
            &[BodyId::Planet(PlanetId::Moon), BodyId::Planet(PlanetId::Sun)],
            TimeInput::from_ut1_jd(JD_J2000),
            &[BodyId::Planet(PlanetId::Sun), BodyId::Planet(PlanetId::Moon)],
            AspectPreset::Major,
            TimeInput::from_ut1_jd(jd_2024),
            TimeInput::from_ut1_jd(jd_2024 + 7.0),
            CalcOptions::default(),
        )
        .expect("natal_transits_in_range");
    assert!(!transits.is_empty());

    // 3g. void_of_course_moon_in_range
    let voc = engine
        .void_of_course_moon_in_range(
            TimeInput::from_ut1_jd(jd_2024),
            TimeInput::from_ut1_jd(jd_2024 + 10.0),
            None,
            CalcOptions::default(),
        )
        .expect("void_of_course_moon_in_range");
    assert!(voc.len() >= 3, "Expected >= 3 VoC intervals in 10 days, got {}", voc.len());
    for v in &voc {
        assert!(v.end_jd_ut1 >= v.start_jd_ut1);
        assert!((v.duration_hours - (v.end_jd_ut1 - v.start_jd_ut1) * 24.0).abs() < 1.0e-9);
        assert_eq!((v.moon_sign_before + 1) % 12, v.moon_sign_after);
    }
}

#[test]
fn test_4_vedic_panchanga_vimshottari_dasha_and_shodashavarga() {
    let engine = Engine::from_data_dir_or_moshier(workspace_data_dir(), EngineConfig::default())
        .expect("load engine");
    let ujjain = Observer::new(23.1765, 75.7885, 491.0).expect("Ujjain observer");

    // 4a. Panchanga
    let panch = engine
        .panchanga(
            TimeInput::from_ut1_jd(JD_J2000),
            ujjain,
            AyanamsaModel::Lahiri,
        )
        .expect("panchanga");
    assert!(panch.sunrise_jd_ut1 <= panch.jd_ut1);
    assert!(panch.sunset_jd_ut1 > panch.sunrise_jd_ut1);
    assert!(panch.next_sunrise_jd_ut1 > panch.sunset_jd_ut1);
    assert!((1..=30).contains(&panch.tithi.tithi_index));
    assert!(panch.tithi.start_jd_ut1 <= panch.jd_ut1 + 1.0e-4);
    assert!(panch.tithi.end_jd_ut1 >= panch.jd_ut1 - 1.0e-4);
    assert!((1..=27).contains(&panch.nakshatra.nakshatra_index));
    assert!((1..=4).contains(&panch.nakshatra.pada));
    assert!((1..=27).contains(&panch.yoga.yoga_index));
    assert!((1..=60).contains(&panch.karana.karana_index));
    assert!(panch.rahukalam.end_jd_ut1 > panch.rahukalam.start_jd_ut1);
    assert!(panch.abhijit_muhurta.end_jd_ut1 > panch.abhijit_muhurta.start_jd_ut1);

    // 4b. Vimshottari Dasha
    let dasha = engine
        .vimshottari_dasha(
            TimeInput::from_ut1_jd(JD_J2000),
            AyanamsaModel::Lahiri,
            Some(365.2425),
        )
        .expect("vimshottari_dasha");
    assert_eq!(dasha.periods.len(), 9);
    assert!((0.0..=1.0).contains(&dasha.balance_fraction));
    for (idx, maha) in dasha.periods.iter().enumerate() {
        if idx == 0 {
            assert!((1..=9).contains(&maha.sub_periods.len()));
        } else {
            assert_eq!(maha.sub_periods.len(), 9);
        }
        assert!((maha.sub_periods.first().expect("first").start_jd_ut1 - maha.start_jd_ut1).abs() < 1.0e-9);
        assert!((maha.sub_periods.last().expect("last").end_jd_ut1 - maha.end_jd_ut1).abs() < 1.0e-9);
    }

    // 4c. calc_varga & calc_shodashavarga
    let v_d9 = engine.calc_varga(45.0, 9).expect("D9 Navamsha");
    assert_eq!(v_d9.division, 9);
    assert_eq!(v_d9.varga_code, "D9");
    assert_eq!(v_d9.varga_name, "Navamsha");

    let shodasha = engine
        .calc_shodashavarga(123.456)
        .expect("calc_shodashavarga");
    assert_eq!(shodasha.len(), 16);
    let expected_divs: [u16; 16] = [1, 2, 3, 4, 7, 9, 10, 12, 16, 20, 24, 27, 30, 40, 45, 60];
    for (v, &expected_div) in shodasha.iter().zip(expected_divs.iter()) {
        assert_eq!(v.division, expected_div);
        assert!((0..12).contains(&v.sign_index));
    }
}

#[test]
fn test_5_sect_planetary_hours_hermetic_lots_and_human_design_64_gates() {
    let engine = Engine::from_data_dir_or_moshier(workspace_data_dir(), EngineConfig::default())
        .expect("load engine");
    let greenwich = Observer::new(51.4769, 0.0, 0.0).expect("Greenwich observer");

    // 5a. chart_sect at Greenwich Noon (JD_J2000 = 2000-01-01 12:00 TT) vs Midnight (JD_J2000 + 0.5)
    let sect_day = engine
        .chart_sect(TimeInput::from_ut1_jd(JD_J2000), greenwich)
        .expect("sect_day");
    assert!(sect_day.is_day_chart);
    assert!(sect_day.sun_altitude_true_deg > 10.0);

    let sect_night = engine
        .chart_sect(TimeInput::from_ut1_jd(JD_J2000 + 0.5), greenwich)
        .expect("sect_night");
    assert!(!sect_night.is_day_chart);
    assert!(sect_night.sun_altitude_true_deg < -40.0);

    // 5b. planetary_hours
    let ph = engine
        .planetary_hours(TimeInput::from_ut1_jd(JD_J2000), greenwich)
        .expect("planetary_hours");
    assert_eq!(ph.hours.len(), 24);
    assert_eq!(ph.hours[0].lord, ph.day_lord);
    assert!(ph.hours.iter().take(12).all(|h| h.is_daytime));
    assert!(ph.hours.iter().skip(12).all(|h| !h.is_daytime));

    // 5c. calc_lots (7 Hermetic Lots with diurnal/nocturnal reversal)
    let lots_day = engine
        .calc_lots(
            TimeInput::from_ut1_jd(JD_J2000),
            greenwich,
            HouseOptions::default(),
        )
        .expect("calc_lots day");
    assert_eq!(lots_day.len(), 7);
    let expected_lot_names = [
        "Fortune",
        "Spirit",
        "Eros",
        "Necessity",
        "Courage",
        "Victory",
        "Nemesis",
    ];
    for (lot, &expected_name) in lots_day.iter().zip(expected_lot_names.iter()) {
        assert_eq!(lot.name, expected_name);
        assert!((0.0..360.0).contains(&lot.lon_deg));
        assert!((1..=12).contains(&lot.house_number));
    }

    // 5d. map_zodiac_placement: 64 Human Design Gates & 6 Lines starting at 302°00' Aquarius = Gate 41 Line 1
    let origin = map_zodiac_placement(302.0);
    assert_eq!(origin.sign_name, "Aquarius");
    assert_eq!(origin.dms_deg, 2);
    assert_eq!(origin.dms_min, 0);
    assert_eq!(origin.hd_gate, 41);
    assert_eq!(origin.hd_line, 1);
    assert_eq!(origin.hd_color, 1);
    assert_eq!(origin.hd_tone, 1);
    assert_eq!(origin.hd_base, 1);

    let mut seen_gates = HashSet::new();
    for (idx, &gate) in HD_MANDALA_GATES_64.iter().enumerate() {
        seen_gates.insert(gate);
        for line in 1..=6u8 {
            let lon = 302.0 + (idx as f64) * 5.625 + f64::from(line - 1) * 0.9375 + 1.0e-5;
            let zp = map_zodiac_placement(lon);
            assert_eq!(zp.hd_gate, gate);
            assert_eq!(zp.hd_line, line);
        }
    }
    assert_eq!(seen_gates.len(), 64);
}

#[test]
fn test_6_chart_snapshot_single_call() {
    let engine = Engine::from_data_dir_or_moshier(workspace_data_dir(), EngineConfig::default())
        .expect("load engine");
    let observer = Observer::new(40.7128, -74.0060, 10.0).expect("NYC observer");
    let bodies = [
        BodyId::Planet(PlanetId::Sun),
        BodyId::Planet(PlanetId::Moon),
        BodyId::Planet(PlanetId::Mercury),
        BodyId::Planet(PlanetId::Venus),
        BodyId::Planet(PlanetId::Mars),
        BodyId::Planet(PlanetId::Jupiter),
        BodyId::Planet(PlanetId::Saturn),
        BodyId::Planet(PlanetId::Uranus),
        BodyId::Planet(PlanetId::Neptune),
        BodyId::Planet(PlanetId::Pluto),
    ];
    let house_opts = HouseOptions {
        system: HouseSystem::Placidus,
        ..HouseOptions::default()
    };

    let snap = engine
        .chart_snapshot(
            TimeInput::from_ut1_jd(JD_J2000),
            observer,
            &bodies,
            house_opts,
            CalcOptions::default(),
            AspectPreset::Major,
            1.0,
        )
        .expect("chart_snapshot");

    assert_eq!(snap.bodies.len(), 10);
    for b in &snap.bodies {
        assert!((1..=12).contains(&b.house_number));
        assert!((1.0..13.0).contains(&b.house_pos_3d));
        assert!((1..=64).contains(&b.placement.hd_gate));
        assert!((1..=6).contains(&b.placement.hd_line));
        assert!((1..=27).contains(&b.placement.nakshatra_index));
    }
    assert_eq!(snap.houses.cusps_deg.len(), 12);
    assert!(!snap.aspects.is_empty());
    assert_eq!(snap.lots.len(), 7);
    assert!((1..=28).contains(&snap.lunar_phase.phase_28_index));
}

// ============================================================================
// Wave H5: 5,000-Vector Golden Oracle Verification (`l3_wave_h_almanac_aspects_5k.json`)
// ============================================================================

#[derive(Debug, Deserialize)]
struct GoldenWaveH5k {
    total_vector_count: usize,
    section_counts: SectionCountsH,
    section_h1_panchanga: SectionH1,
    section_h2_vargas_and_dasha: SectionH2,
    section_h3_traditional_almanac: SectionH3,
    section_h4_aspects_returns_voc: SectionH4,
    section_h5_hd_nakshatra_snapshot: SectionH5,
}

#[derive(Debug, Deserialize)]
struct SectionCountsH {
    section_h1_panchanga: usize,
    section_h2_vargas_and_dasha: usize,
    section_h3_traditional_almanac: usize,
    section_h4_aspects_returns_voc: usize,
    section_h5_hd_nakshatra_snapshot: usize,
}

#[derive(Debug, Deserialize)]
struct SectionH1 {
    vector_count: usize,
    panchanga_1000: Vec<GoldenPanchangaVec>,
}

#[derive(Debug, Deserialize)]
struct GoldenTimeWindow {
    start_jd_ut1: f64,
    end_jd_ut1: f64,
}

#[derive(Debug, Deserialize)]
struct GoldenTithi {
    tithi_index: u8,
    tithi_name: String,
    paksha: String,
    elapsed_fraction: f64,
    start_jd_ut1: f64,
    end_jd_ut1: f64,
}

#[derive(Debug, Deserialize)]
struct GoldenNakshatra {
    nakshatra_index: u8,
    nakshatra_name: String,
    vimshottari_lord: String,
    pada: u8,
    elapsed_fraction: f64,
    start_jd_ut1: f64,
    end_jd_ut1: f64,
}

#[derive(Debug, Deserialize)]
struct GoldenYoga {
    yoga_index: u8,
    yoga_name: String,
    sum_lon_deg: f64,
    elapsed_fraction: f64,
    start_jd_ut1: f64,
    end_jd_ut1: f64,
}

#[derive(Debug, Deserialize)]
struct GoldenKarana {
    karana_index: u8,
    karana_name: String,
    is_vishti_bhadra: bool,
    elapsed_fraction: f64,
    start_jd_ut1: f64,
    end_jd_ut1: f64,
}

#[derive(Debug, Deserialize)]
struct GoldenPanchangaVec {
    jd_ut1: f64,
    city: String,
    lat_deg: f64,
    lon_deg: f64,
    alt_m: f64,
    ayanamsa: String,
    ayanamsa_deg: f64,
    sunrise_jd_ut1: f64,
    sunset_jd_ut1: f64,
    next_sunrise_jd_ut1: f64,
    vara_index: u8,
    vara_name: String,
    vara_lord: String,
    tithi: GoldenTithi,
    nakshatra: GoldenNakshatra,
    yoga: GoldenYoga,
    karana: GoldenKarana,
    rahukalam: GoldenTimeWindow,
    yamagandam: GoldenTimeWindow,
    gulikakalam: GoldenTimeWindow,
    abhijit_muhurta: GoldenTimeWindow,
}

#[derive(Debug, Deserialize)]
struct SectionH2 {
    vector_count: usize,
    varga_1000: Vec<GoldenVargaVec>,
    vimshottari_dasha_200: Vec<GoldenDashaVec>,
}

#[derive(Debug, Deserialize)]
struct GoldenVargaVec {
    division: u16,
    varga_code: String,
    varga_name: String,
    source_sidereal_lon_deg: f64,
    sign_index: u8,
    sign_name: String,
    sign_lord: String,
    amsha_index: u16,
    varga_lon_deg: f64,
    deg_in_varga_sign: f64,
}

#[derive(Debug, Deserialize)]
struct GoldenSubPeriod {
    maha_lord: String,
    lord: String,
    start_jd_ut1: f64,
    end_jd_ut1: f64,
    duration_years: f64,
}

#[derive(Debug, Deserialize)]
struct GoldenDashaPeriod {
    lord: String,
    start_jd_ut1: f64,
    end_jd_ut1: f64,
    duration_years: f64,
    first_antardasha: GoldenSubPeriod,
    mid_antardasha: GoldenSubPeriod,
    last_antardasha: GoldenSubPeriod,
}

#[derive(Debug, Deserialize)]
struct GoldenDashaVec {
    birth_jd_ut1: f64,
    ayanamsa: String,
    year_length_days: f64,
    moon_sidereal_lon_deg: f64,
    birth_nakshatra_index: u8,
    birth_nakshatra_name: String,
    birth_pada: u8,
    balance_fraction: f64,
    balance_years: f64,
    cycle_end_jd_ut1: f64,
    periods: Vec<GoldenDashaPeriod>,
}

#[derive(Debug, Deserialize)]
struct SectionH3 {
    vector_count: usize,
    traditional_almanac_1000: Vec<GoldenTradAlmanacVec>,
}

#[derive(Debug, Deserialize)]
struct GoldenSect {
    is_day_chart: bool,
    sun_altitude_true_deg: f64,
    sun_altitude_apparent_deg: f64,
    near_horizon: bool,
}

#[derive(Debug, Deserialize)]
struct GoldenPlanetaryHours {
    sunrise_jd_ut1: f64,
    sunset_jd_ut1: f64,
    next_sunrise_jd_ut1: f64,
    day_lord: String,
    current_hour_index: u8,
    current_hour_lord: String,
    is_daytime_hour: bool,
    hour_1_lord: String,
    hour_12_end_jd_ut1: f64,
    hour_13_lord: String,
    hour_24_lord: String,
}

#[derive(Debug, Deserialize)]
struct GoldenLot {
    name: String,
    lon_deg: f64,
    sign_index: u8,
    sign_name: String,
    house_number: u8,
}

#[derive(Debug, Deserialize)]
struct GoldenTradAlmanacVec {
    jd_ut1: f64,
    city: String,
    lat_deg: f64,
    lon_deg: f64,
    alt_m: f64,
    house_system: String,
    sect: GoldenSect,
    planetary_hours: GoldenPlanetaryHours,
    lots: Vec<GoldenLot>,
}

#[derive(Debug, Deserialize)]
struct SectionH4 {
    vector_count: usize,
    static_aspects_300: Vec<serde_json::Value>,
    dynamic_exact_aspects_300: Vec<serde_json::Value>,
    planetary_returns_200: Vec<serde_json::Value>,
    void_of_course_moon_200: Vec<GoldenVocVec>,
}

#[derive(Debug, Deserialize)]
struct GoldenVocVec {
    zodiac: String,
    window_start_jd_ut1: f64,
    window_end_jd_ut1: f64,
    start_jd_ut1: f64,
    end_jd_ut1: f64,
    duration_hours: f64,
    last_aspect_body: String,
    last_aspect_name: String,
    last_aspect_angle_deg: f64,
    current_sign_index: u8,
    next_sign_index: u8,
    ingress_lon_deg: f64,
}

#[derive(Debug, Deserialize)]
struct SectionH5 {
    vector_count: usize,
    human_design_384: Vec<GoldenHdVec>,
    nakshatra_padas_108: Vec<GoldenNakshatraPadaVec>,
    chart_snapshots_308: Vec<serde_json::Value>,
}

#[derive(Debug, Deserialize)]
struct GoldenHdVec {
    mandala_slot_index: usize,
    lon_deg: f64,
    sign_index: u8,
    sign_name: String,
    hd_gate: u8,
    hd_line: u8,
    hd_color: u8,
    hd_tone: u8,
    hd_base: u8,
}

#[derive(Debug, Deserialize)]
struct GoldenNakshatraPadaVec {
    global_pada_index: usize,
    lon_deg: f64,
    sign_index: u8,
    sign_name: String,
    decan: u8,
    nakshatra_index: u8,
    nakshatra_name: String,
    nakshatra_lord: String,
    pada: u8,
    d9_navamsa_sign_index: u8,
    d9_navamsa_sign_name: String,
}

fn parse_ayanamsa(name: &str) -> AyanamsaModel {
    match name {
        "lahiri" => AyanamsaModel::Lahiri,
        "raman" => AyanamsaModel::Raman,
        "krishnamurti" => AyanamsaModel::Krishnamurti,
        "true_citra" => AyanamsaModel::TrueCitra,
        other => panic!("Unknown ayanamsa model: {other}"),
    }
}

fn parse_house_system(code: &str) -> HouseSystem {
    let ch = code.chars().next().expect("house code char");
    HouseSystem::from_code(ch).unwrap_or_else(|_| panic!("Unknown house code: {code}"))
}

fn parse_planet_body(name: &str) -> BodyId {
    match name {
        "Sun" => BodyId::Planet(PlanetId::Sun),
        "Moon" => BodyId::Planet(PlanetId::Moon),
        "Mercury" => BodyId::Planet(PlanetId::Mercury),
        "Venus" => BodyId::Planet(PlanetId::Venus),
        "Mars" => BodyId::Planet(PlanetId::Mars),
        "Jupiter" => BodyId::Planet(PlanetId::Jupiter),
        "Saturn" => BodyId::Planet(PlanetId::Saturn),
        "Uranus" => BodyId::Planet(PlanetId::Uranus),
        "Neptune" => BodyId::Planet(PlanetId::Neptune),
        "Pluto" => BodyId::Planet(PlanetId::Pluto),
        other => panic!("Unknown planet name: {other}"),
    }
}

fn parse_zodiac(name: &str) -> Zodiac {
    match name {
        "tropical" => Zodiac::Tropical,
        other => Zodiac::Sidereal(parse_ayanamsa(other)),
    }
}

#[test]
fn test_7_golden_wave_h_5000_vectors_all_sections() {
    let golden_path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
        .join("testdata")
        .join("golden")
        .join("l3_wave_h_almanac_aspects_5k.json");
    let raw = std::fs::read_to_string(&golden_path).expect("read l3_wave_h_almanac_aspects_5k.json");
    let data: GoldenWaveH5k =
        serde_json::from_str(&raw).expect("parse l3_wave_h_almanac_aspects_5k.json");

    assert_eq!(data.total_vector_count, 5_000);
    assert_eq!(data.section_counts.section_h1_panchanga, 1_000);
    assert_eq!(data.section_counts.section_h2_vargas_and_dasha, 1_200);
    assert_eq!(data.section_counts.section_h3_traditional_almanac, 1_000);
    assert_eq!(data.section_counts.section_h4_aspects_returns_voc, 1_000);
    assert_eq!(data.section_counts.section_h5_hd_nakshatra_snapshot, 800);

    assert_eq!(data.section_h1_panchanga.panchanga_1000.len(), 1_000);
    assert_eq!(data.section_h2_vargas_and_dasha.varga_1000.len(), 1_000);
    assert_eq!(
        data.section_h2_vargas_and_dasha.vimshottari_dasha_200.len(),
        200
    );
    assert_eq!(
        data.section_h3_traditional_almanac
            .traditional_almanac_1000
            .len(),
        1_000
    );
    assert_eq!(
        data.section_h4_aspects_returns_voc.static_aspects_300.len(),
        300
    );
    assert_eq!(
        data.section_h4_aspects_returns_voc
            .dynamic_exact_aspects_300
            .len(),
        300
    );
    assert_eq!(
        data.section_h4_aspects_returns_voc
            .planetary_returns_200
            .len(),
        200
    );
    assert_eq!(
        data.section_h4_aspects_returns_voc
            .void_of_course_moon_200
            .len(),
        200
    );
    assert_eq!(
        data.section_h5_hd_nakshatra_snapshot.human_design_384.len(),
        384
    );
    assert_eq!(
        data.section_h5_hd_nakshatra_snapshot
            .nakshatra_padas_108
            .len(),
        108
    );
    assert_eq!(
        data.section_h5_hd_nakshatra_snapshot
            .chart_snapshots_308
            .len(),
        308
    );

    let engine = Arc::new(
        Engine::from_data_dir_or_moshier(workspace_data_dir(), EngineConfig::default())
            .expect("load full DE440s engine"),
    );

    let (h1_metrics, h2_metrics, h3_metrics, h4_metrics, h5_metrics) = std::thread::scope(|s| {
        // --------------------------------------------------------------------
        // Section H1: 1,000 Vedic Panchanga & Muhurta vectors
        // --------------------------------------------------------------------
        let eng_h1 = Arc::clone(&engine);
        let sec_h1 = &data.section_h1_panchanga.panchanga_1000;
        let t_h1 = s.spawn(move || {
            let mut max_tk_days = 0.0_f64;
            let mut sum_sq_tk_days = 0.0_f64;
            let mut tk_count = 0usize;
            let mut max_sr_days = 0.0_f64;

            for (idx, v) in sec_h1.iter().enumerate() {
                let obs = Observer::new(v.lat_deg, v.lon_deg, v.alt_m).expect("observer");
                let aya = parse_ayanamsa(&v.ayanamsa);
                let p = eng_h1
                    .panchanga(TimeInput::from_ut1_jd(v.jd_ut1), obs, aya)
                    .unwrap_or_else(|e| panic!("H1[{idx}] {} {} failed: {e}", v.city, v.jd_ut1));

                assert_eq!(p.vara.vara_index, v.vara_index, "H1[{idx}] vara_index");
                assert_eq!(p.vara.vara_name, v.vara_name, "H1[{idx}] vara_name");
                assert_eq!(p.vara.vara_lord, v.vara_lord, "H1[{idx}] vara_lord");

                assert_eq!(
                    p.tithi.tithi_index, v.tithi.tithi_index,
                    "H1[{idx}] tithi_index"
                );
                assert_eq!(
                    p.tithi.tithi_name, v.tithi.tithi_name,
                    "H1[{idx}] tithi_name"
                );
                assert_eq!(p.tithi.paksha, v.tithi.paksha, "H1[{idx}] paksha");

                let dt_ts = (p.tithi.start_jd_ut1 - v.tithi.start_jd_ut1).abs();
                let dt_te = (p.tithi.end_jd_ut1 - v.tithi.end_jd_ut1).abs();
                let dt_ks = (p.karana.start_jd_ut1 - v.karana.start_jd_ut1).abs();
                let dt_ke = (p.karana.end_jd_ut1 - v.karana.end_jd_ut1).abs();
                for d in [dt_ts, dt_te, dt_ks, dt_ke] {
                    assert!(
                        d < 2.0e-5,
                        "H1[{idx}] tithi/karana boundary diff {d:.2e} days >= 2e-5 days"
                    );
                    max_tk_days = max_tk_days.max(d);
                    sum_sq_tk_days += d * d;
                    tk_count += 1;
                }

                assert_eq!(
                    p.nakshatra.nakshatra_index, v.nakshatra.nakshatra_index,
                    "H1[{idx}] nakshatra_index"
                );
                assert_eq!(
                    p.nakshatra.nakshatra_name, v.nakshatra.nakshatra_name,
                    "H1[{idx}] nakshatra_name"
                );
                assert_eq!(
                    p.nakshatra.vimshottari_lord, v.nakshatra.vimshottari_lord,
                    "H1[{idx}] vimshottari_lord"
                );
                assert_eq!(p.nakshatra.pada, v.nakshatra.pada, "H1[{idx}] pada");
                assert!((p.nakshatra.start_jd_ut1 - v.nakshatra.start_jd_ut1).abs() < 4.0e-4);
                assert!((p.nakshatra.end_jd_ut1 - v.nakshatra.end_jd_ut1).abs() < 4.0e-4);

                assert_eq!(p.yoga.yoga_index, v.yoga.yoga_index, "H1[{idx}] yoga_index");
                assert_eq!(p.yoga.yoga_name, v.yoga.yoga_name, "H1[{idx}] yoga_name");
                assert!((p.yoga.start_jd_ut1 - v.yoga.start_jd_ut1).abs() < 4.0e-4);
                assert!((p.yoga.end_jd_ut1 - v.yoga.end_jd_ut1).abs() < 4.0e-4);

                assert_eq!(
                    p.karana.karana_index, v.karana.karana_index,
                    "H1[{idx}] karana_index"
                );
                assert_eq!(
                    p.karana.karana_name, v.karana.karana_name,
                    "H1[{idx}] karana_name"
                );
                assert_eq!(
                    p.karana.is_vishti_bhadra, v.karana.is_vishti_bhadra,
                    "H1[{idx}] is_vishti_bhadra"
                );

                let d_sr = (p.sunrise_jd_ut1 - v.sunrise_jd_ut1).abs();
                let d_ss = (p.sunset_jd_ut1 - v.sunset_jd_ut1).abs();
                let d_nsr = (p.next_sunrise_jd_ut1 - v.next_sunrise_jd_ut1).abs();
                max_sr_days = max_sr_days.max(d_sr).max(d_ss).max(d_nsr);
                assert!(
                    d_sr < 3.5e-4 && d_ss < 3.5e-4 && d_nsr < 3.5e-4,
                    "H1[{idx}] sunrise/sunset diff: sr={d_sr:.2e}, ss={d_ss:.2e}, nsr={d_nsr:.2e}"
                );
                assert!(
                    (p.rahukalam.start_jd_ut1 - v.rahukalam.start_jd_ut1).abs() < 3.5e-4
                        && (p.rahukalam.end_jd_ut1 - v.rahukalam.end_jd_ut1).abs() < 3.5e-4
                );
                assert!(
                    (p.yamagandam.start_jd_ut1 - v.yamagandam.start_jd_ut1).abs() < 3.5e-4
                        && (p.yamagandam.end_jd_ut1 - v.yamagandam.end_jd_ut1).abs() < 3.5e-4
                );
                assert!(
                    (p.gulikakalam.start_jd_ut1 - v.gulikakalam.start_jd_ut1).abs() < 3.5e-4
                        && (p.gulikakalam.end_jd_ut1 - v.gulikakalam.end_jd_ut1).abs() < 3.5e-4
                );
                assert!(
                    (p.abhijit_muhurta.start_jd_ut1 - v.abhijit_muhurta.start_jd_ut1).abs()
                        < 3.5e-4
                        && (p.abhijit_muhurta.end_jd_ut1 - v.abhijit_muhurta.end_jd_ut1).abs()
                            < 3.5e-4
                );
            }

            let rms_tk_days = (sum_sq_tk_days / (tk_count as f64)).sqrt();
            (max_tk_days, rms_tk_days, max_sr_days)
        });

        // --------------------------------------------------------------------
        // Section H2: 1,200 Shodashavarga (1,000) & Vimshottari Dasha (200)
        // --------------------------------------------------------------------
        let eng_h2 = Arc::clone(&engine);
        let sec_h2 = &data.section_h2_vargas_and_dasha;
        let t_h2 = s.spawn(move || {
            let mut max_varga_lon_err = 0.0_f64;
            let mut sum_sq_varga_lon = 0.0_f64;

            for (idx, v) in sec_h2.varga_1000.iter().enumerate() {
                let res = eng_h2
                    .calc_varga(v.source_sidereal_lon_deg, v.division)
                    .unwrap_or_else(|e| panic!("H2A[{idx}] calc_varga failed: {e}"));
                assert_eq!(res.division, v.division);
                assert_eq!(res.varga_code, v.varga_code);
                assert_eq!(res.varga_name, v.varga_name);
                assert_eq!(res.sign_index, v.sign_index, "H2A[{idx}] sign_index");
                assert_eq!(res.sign_name, v.sign_name, "H2A[{idx}] sign_name");
                assert_eq!(res.sign_lord, v.sign_lord, "H2A[{idx}] sign_lord");
                assert_eq!(res.amsha_index, v.amsha_index, "H2A[{idx}] amsha_index");

                let err = wrap_pm180(res.varga_lon_deg - v.varga_lon_deg).abs();
                assert!(err < 1.0e-9, "H2A[{idx}] varga_lon_deg err={err:.3e}");
                max_varga_lon_err = max_varga_lon_err.max(err);
                sum_sq_varga_lon += err * err;

                if v.division == 1 {
                    let all16 = eng_h2
                        .calc_shodashavarga(v.source_sidereal_lon_deg)
                        .expect("calc_shodashavarga");
                    assert_eq!(all16.len(), 16);
                    assert_eq!(all16[0].sign_index, v.sign_index);
                }
            }

            let mut max_dasha_exact_jd_err = 0.0_f64;
            let mut max_dasha_engine_jd_err = 0.0_f64;
            for (idx, v) in sec_h2.vimshottari_dasha_200.iter().enumerate() {
                // 1. Exact formula verification from oracle Moon sidereal longitude
                let rep_exact = vimshottari_dasha_from_lon(
                    v.moon_sidereal_lon_deg,
                    v.birth_jd_ut1,
                    Some(v.year_length_days),
                )
                .expect("vimshottari_dasha_from_lon");
                assert_eq!(rep_exact.birth_nakshatra_index, v.birth_nakshatra_index);
                assert_eq!(rep_exact.birth_nakshatra_name, v.birth_nakshatra_name);
                assert_eq!(rep_exact.birth_pada, v.birth_pada);
                assert!((rep_exact.balance_fraction - v.balance_fraction).abs() < 1.0e-12);
                assert!((rep_exact.balance_years - v.balance_years).abs() < 1.0e-11);
                assert_eq!(rep_exact.periods.len(), 9);
                for (p_got, p_ref) in rep_exact.periods.iter().zip(v.periods.iter()) {
                    assert_eq!(p_got.lord, p_ref.lord);
                    let e_start = (p_got.start_jd_ut1 - p_ref.start_jd_ut1).abs();
                    let e_end = (p_got.end_jd_ut1 - p_ref.end_jd_ut1).abs();
                    max_dasha_exact_jd_err = max_dasha_exact_jd_err.max(e_start).max(e_end);
                    assert!(e_start < 1.0e-8 && e_end < 1.0e-8);
                    let mid_idx = p_got.sub_periods.len() / 2;
                    let last_idx = p_got.sub_periods.len() - 1;
                    assert_eq!(p_got.sub_periods[0].lord, p_ref.first_antardasha.lord);
                    assert_eq!(p_got.sub_periods[mid_idx].lord, p_ref.mid_antardasha.lord);
                    assert_eq!(p_got.sub_periods[last_idx].lord, p_ref.last_antardasha.lord);
                    assert!(
                        (p_got.sub_periods[mid_idx].end_jd_ut1 - p_ref.mid_antardasha.end_jd_ut1)
                            .abs()
                            < 1.0e-8
                    );
                }

                // 2. Full Engine DE440s Moon evaluation at birth_jd_ut1
                let aya = parse_ayanamsa(&v.ayanamsa);
                let rep_eng = eng_h2
                    .vimshottari_dasha(
                        TimeInput::from_ut1_jd(v.birth_jd_ut1),
                        aya,
                        Some(v.year_length_days),
                    )
                    .expect("engine.vimshottari_dasha");
                assert_eq!(rep_eng.birth_nakshatra_index, v.birth_nakshatra_index);
                assert_eq!(rep_eng.birth_nakshatra_name, v.birth_nakshatra_name);
                assert_eq!(rep_eng.birth_pada, v.birth_pada);
                let end_err = (rep_eng.periods[8].end_jd_ut1 - v.cycle_end_jd_ut1).abs();
                max_dasha_engine_jd_err = max_dasha_engine_jd_err.max(end_err);
                assert!(
                    end_err < 0.25,
                    "H2B[{idx}] 120-year cycle end diff={end_err:.4e} days"
                );
            }

            let rms_varga_lon = (sum_sq_varga_lon / 1000.0).sqrt();
            (
                max_varga_lon_err,
                rms_varga_lon,
                max_dasha_exact_jd_err,
                max_dasha_engine_jd_err,
            )
        });

        // --------------------------------------------------------------------
        // Section H3: 1,000 Traditional Western Almanac vectors
        // --------------------------------------------------------------------
        let eng_h3 = Arc::clone(&engine);
        let sec_h3 = &data.section_h3_traditional_almanac.traditional_almanac_1000;
        let t_h3 = s.spawn(move || {
            let mut max_lot_lon_err = 0.0_f64;
            let mut sum_sq_lot_lon = 0.0_f64;
            let mut lot_count = 0usize;
            let mut max_alt_err = 0.0_f64;

            for (idx, v) in sec_h3.iter().enumerate() {
                let obs = Observer::new(v.lat_deg, v.lon_deg, v.alt_m).expect("observer");
                let input = TimeInput::from_ut1_jd(v.jd_ut1);

                // 1. chart_sect
                let sect = eng_h3
                    .chart_sect(input, obs)
                    .unwrap_or_else(|e| panic!("H3[{idx}] chart_sect failed: {e}"));
                assert_eq!(sect.is_day_chart, v.sect.is_day_chart, "H3[{idx}] is_day_chart");
                assert_eq!(sect.near_horizon, v.sect.near_horizon, "H3[{idx}] near_horizon");
                let alt_err = (sect.sun_altitude_true_deg - v.sect.sun_altitude_true_deg).abs();
                max_alt_err = max_alt_err.max(alt_err);
                assert!(alt_err < 0.01, "H3[{idx}] sun_altitude_true_deg diff={alt_err}");

                // 2. planetary_hours
                let ph = eng_h3
                    .planetary_hours(input, obs)
                    .unwrap_or_else(|e| panic!("H3[{idx}] planetary_hours failed: {e}"));
                assert_eq!(ph.day_lord, v.planetary_hours.day_lord);
                assert_eq!(ph.current_hour_index, v.planetary_hours.current_hour_index);
                assert_eq!(ph.current_hour_lord, v.planetary_hours.current_hour_lord);
                assert_eq!(ph.is_daytime_hour, v.planetary_hours.is_daytime_hour);
                assert_eq!(ph.hours[0].lord, v.planetary_hours.hour_1_lord);
                assert_eq!(ph.hours[12].lord, v.planetary_hours.hour_13_lord);
                assert_eq!(ph.hours[23].lord, v.planetary_hours.hour_24_lord);
                assert!((ph.hours[11].end_jd_ut1 - v.planetary_hours.hour_12_end_jd_ut1).abs() < 3.5e-4);

                // 3. calc_lots (all 7 Hermetic Lots)
                let h_opts = HouseOptions {
                    system: parse_house_system(&v.house_system),
                    ..HouseOptions::default()
                };
                let lots = eng_h3
                    .calc_lots(input, obs, h_opts)
                    .unwrap_or_else(|e| panic!("H3[{idx}] calc_lots failed: {e}"));
                assert_eq!(lots.len(), 7);
                for (l_got, l_ref) in lots.iter().zip(v.lots.iter()) {
                    assert_eq!(l_got.name, l_ref.name);
                    assert_eq!(l_got.placement.sign_index, l_ref.sign_index);
                    assert_eq!(l_got.placement.sign_name, l_ref.sign_name);
                    assert_eq!(l_got.house_number, l_ref.house_number);
                    let l_err = wrap_pm180(l_got.lon_deg - l_ref.lon_deg).abs();
                    assert!(
                        l_err < 0.005,
                        "H3[{idx}] lot {} lon diff={l_err:.5}",
                        l_ref.name
                    );
                    max_lot_lon_err = max_lot_lon_err.max(l_err);
                    sum_sq_lot_lon += l_err * l_err;
                    lot_count += 1;
                }
            }

            let rms_lot_lon = (sum_sq_lot_lon / (lot_count as f64)).sqrt();
            (max_lot_lon_err, rms_lot_lon, max_alt_err)
        });

        // --------------------------------------------------------------------
        // Section H4: 1,000 Aspects, Transits, Returns & Void-of-Course Moon
        // --------------------------------------------------------------------
        let eng_h4 = Arc::clone(&engine);
        let sec_h4 = &data.section_h4_aspects_returns_voc;
        let t_h4 = s.spawn(move || {
            let bodies_10 = [
                BodyId::Planet(PlanetId::Sun),
                BodyId::Planet(PlanetId::Moon),
                BodyId::Planet(PlanetId::Mercury),
                BodyId::Planet(PlanetId::Venus),
                BodyId::Planet(PlanetId::Mars),
                BodyId::Planet(PlanetId::Jupiter),
                BodyId::Planet(PlanetId::Saturn),
                BodyId::Planet(PlanetId::Uranus),
                BodyId::Planet(PlanetId::Neptune),
                BodyId::Planet(PlanetId::Pluto),
            ];
            let bodies_7 = &bodies_10[..7];

            // H4A: 300 Static Aspects
            let mut max_static_orb_err = 0.0_f64;
            for (idx, val) in sec_h4.static_aspects_300.iter().enumerate() {
                let kind = val["kind"].as_str().expect("kind");
                match kind {
                    "chart_aspect" => {
                        let jd_ut1 = val["jd_ut1"].as_f64().expect("jd_ut1");
                        let b1 = parse_planet_body(val["body1"].as_str().expect("b1"));
                        let b2 = parse_planet_body(val["body2"].as_str().expect("b2"));
                        let aname = val["aspect_name"].as_str().expect("aspect_name");
                        let ref_orb = val["orb_deg"].as_f64().expect("orb_deg");
                        let ref_apply = val["is_applying"].as_bool().expect("is_applying");

                        let hits = eng_h4
                            .calc_chart_aspects(
                                TimeInput::from_ut1_jd(jd_ut1),
                                &bodies_10,
                                AspectPreset::All,
                                1.0,
                                CalcOptions::default(),
                            )
                            .expect("calc_chart_aspects");
                        let hit = hits
                            .iter()
                            .find(|h| {
                                ((h.body1 == b1 && h.body2 == b2)
                                    || (h.body1 == b2 && h.body2 == b1))
                                    && h.aspect_name == aname
                            })
                            .unwrap_or_else(|| {
                                panic!("H4A[{idx}] chart_aspect {b1:?}-{b2:?} {aname} not found")
                            });
                        assert_eq!(hit.is_applying, ref_apply, "H4A[{idx}] is_applying");
                        let orb_err = (hit.orb_deg - ref_orb).abs();
                        max_static_orb_err = max_static_orb_err.max(orb_err);
                        assert!(orb_err < 0.005, "H4A[{idx}] orb_err={orb_err}");
                    }
                    "synastry_aspect" => {
                        let jd1 = val["jd_ut1_1"].as_f64().expect("jd1");
                        let jd2 = val["jd_ut1_2"].as_f64().expect("jd2");
                        let b1 = parse_planet_body(val["body1"].as_str().expect("b1"));
                        let b2 = parse_planet_body(val["body2"].as_str().expect("b2"));
                        let aname = val["aspect_name"].as_str().expect("aspect_name");
                        let ref_orb = val["orb_deg"].as_f64().expect("orb_deg");
                        let ref_apply = val["is_applying"].as_bool().expect("is_applying");

                        let hits = eng_h4
                            .calc_synastry_aspects(
                                TimeInput::from_ut1_jd(jd1),
                                bodies_7,
                                TimeInput::from_ut1_jd(jd2),
                                bodies_7,
                                AspectPreset::Major,
                                1.0,
                                CalcOptions::default(),
                            )
                            .expect("calc_synastry_aspects");
                        let hit = hits
                            .iter()
                            .find(|h| h.body1 == b1 && h.body2 == b2 && h.aspect_name == aname)
                            .unwrap_or_else(|| {
                                panic!("H4A[{idx}] synastry_aspect {b1:?}-{b2:?} {aname} not found")
                            });
                        assert_eq!(hit.is_applying, ref_apply, "H4A[{idx}] synastry is_applying");
                        let orb_err = (hit.orb_deg - ref_orb).abs();
                        max_static_orb_err = max_static_orb_err.max(orb_err);
                        assert!(orb_err < 0.005, "H4A[{idx}] synastry orb_err={orb_err}");
                    }
                    "declination_aspect" => {
                        let jd_ut1 = val["jd_ut1"].as_f64().expect("jd_ut1");
                        let b1 = parse_planet_body(val["body1"].as_str().expect("b1"));
                        let b2 = parse_planet_body(val["body2"].as_str().expect("b2"));
                        let is_parallel = val["is_parallel"].as_bool().expect("is_parallel");
                        let ref_orb = val["orb_deg"].as_f64().expect("orb_deg");
                        let ref_apply = val["is_applying"].as_bool().expect("is_applying");
                        let exp_kind = if is_parallel {
                            DeclinationAspectKind::Parallel
                        } else {
                            DeclinationAspectKind::Contraparallel
                        };

                        let hits = eng_h4
                            .calc_declination_aspects(
                                TimeInput::from_ut1_jd(jd_ut1),
                                &bodies_10,
                                1.2,
                                CalcOptions::default(),
                            )
                            .expect("calc_declination_aspects");
                        let hit = hits
                            .iter()
                            .find(|h| {
                                ((h.body1 == b1 && h.body2 == b2)
                                    || (h.body1 == b2 && h.body2 == b1))
                                    && h.kind == exp_kind
                            })
                            .unwrap_or_else(|| {
                                panic!("H4A[{idx}] declination {b1:?}-{b2:?} {exp_kind:?} not found")
                            });
                        assert_eq!(hit.is_applying, ref_apply, "H4A[{idx}] decl is_applying");
                        let orb_err = (hit.orb_deg - ref_orb).abs();
                        max_static_orb_err = max_static_orb_err.max(orb_err);
                        assert!(orb_err < 0.005, "H4A[{idx}] decl orb_err={orb_err}");
                    }
                    other => panic!("Unknown H4A kind: {other}"),
                }
            }

            // H4B: 300 Dynamic Exact Aspects (incl. retrograde triple passes) — strictly < 1e-5 days!
            let mut max_dyn_aspect_days = 0.0_f64;
            let mut sum_sq_dyn_days = 0.0_f64;
            for (idx, val) in sec_h4.dynamic_exact_aspects_300.iter().enumerate() {
                let kind = val["kind"].as_str().expect("kind");
                let b1 = parse_planet_body(val["body1"].as_str().expect("b1"));
                let b2 = parse_planet_body(val["body2"].as_str().expect("b2"));
                let angle = val["exact_angle_deg"].as_f64().expect("angle");
                let ref_jd = val["jd_ut1"].as_f64().expect("jd_ut1");
                let ref_pass = val["pass_index"].as_u64().expect("pass_index") as u32;

                let got_jd = match kind {
                    "search_exact_aspect" => {
                        let s_start = val["search_start_jd_ut1"].as_f64().expect("search_start");
                        let ev = eng_h4
                            .search_exact_aspect(
                                b1,
                                b2,
                                angle,
                                TimeInput::from_ut1_jd(s_start),
                                true,
                                CalcOptions::default(),
                            )
                            .unwrap_or_else(|e| panic!("H4B[{idx}] search_exact_aspect failed: {e}"));
                        ev.jd_ut1
                    }
                    "mutual_aspect" => {
                        let w_start = val["window_start_jd_ut1"].as_f64().expect("w_start");
                        let w_end = val["window_end_jd_ut1"].as_f64().expect("w_end");
                        let evs = eng_h4
                            .mutual_aspects_in_range(
                                &[b1, b2],
                                AspectPreset::Major,
                                TimeInput::from_ut1_jd(w_start),
                                TimeInput::from_ut1_jd(w_end),
                                CalcOptions::default(),
                            )
                            .unwrap_or_else(|e| panic!("H4B[{idx}] mutual_aspects_in_range failed: {e}"));
                        let matching: Vec<_> = evs
                            .iter()
                            .filter(|e| (e.target_angle_deg - angle).abs() < 1.0e-6)
                            .collect();
                        let ev = matching
                            .get((ref_pass - 1) as usize)
                            .unwrap_or_else(|| panic!("H4B[{idx}] missing pass {ref_pass}"));
                        ev.jd_ut1
                    }
                    "natal_transit" => {
                        let nat_jd = val["natal_jd_ut1"].as_f64().expect("natal_jd_ut1");
                        let w_start = val["window_start_jd_ut1"].as_f64().expect("w_start");
                        let w_end = val["window_end_jd_ut1"].as_f64().expect("w_end");
                        let evs = eng_h4
                            .natal_transits_in_range(
                                &[b1],
                                TimeInput::from_ut1_jd(nat_jd),
                                &[b2],
                                AspectPreset::Major,
                                TimeInput::from_ut1_jd(w_start),
                                TimeInput::from_ut1_jd(w_end),
                                CalcOptions::default(),
                            )
                            .unwrap_or_else(|e| panic!("H4B[{idx}] natal_transits_in_range failed: {e}"));
                        let matching: Vec<_> = evs
                            .iter()
                            .filter(|e| (e.target_angle_deg - angle).abs() < 1.0e-6)
                            .collect();
                        let ev = matching
                            .get((ref_pass - 1) as usize)
                            .unwrap_or_else(|| panic!("H4B[{idx}] natal_transit missing pass {ref_pass}"));
                        ev.jd_ut1
                    }
                    other => panic!("Unknown H4B kind: {other}"),
                };

                let dt = (got_jd - ref_jd).abs();
                assert!(
                    dt < 1.0e-5,
                    "H4B[{idx}] {kind} {b1:?}-{b2:?} {angle}° jd_ut1 diff={dt:.2e} days >= 1e-5 days"
                );
                max_dyn_aspect_days = max_dyn_aspect_days.max(dt);
                sum_sq_dyn_days += dt * dt;
            }

            // H4C: 200 Planetary Returns (100 single + 100 range incl. 3-pass retrograde returns)
            let mut max_return_days = 0.0_f64;
            for (idx, val) in sec_h4.planetary_returns_200.iter().enumerate() {
                let kind = val["kind"].as_str().expect("kind");
                let body = parse_planet_body(val["body"].as_str().expect("body"));
                let zodiac = parse_zodiac(val["zodiac"].as_str().expect("zodiac"));
                let nat_jd = val["natal_jd_ut1"].as_f64().expect("natal_jd_ut1");
                let ref_jd = val["jd_ut1"].as_f64().expect("jd_ut1");
                let ref_retro = val["is_retrograde"].as_bool().expect("is_retrograde");
                let ref_pass = val["pass_index"].as_u64().expect("pass_index") as u32;

                let opts = CalcOptions {
                    zodiac,
                    ..CalcOptions::default()
                };

                let got_jd = match kind {
                    "single_return" => {
                        let s_after = val["search_after_jd_ut1"].as_f64().expect("search_after");
                        let ev = eng_h4
                            .planetary_return(
                                body,
                                TimeInput::from_ut1_jd(nat_jd),
                                TimeInput::from_ut1_jd(s_after),
                                opts,
                            )
                            .unwrap_or_else(|e| panic!("H4C[{idx}] planetary_return failed: {e}"));
                        assert_eq!(ev.is_retrograde, ref_retro);
                        assert_eq!(ev.pass_index, 1);
                        ev.jd_ut1
                    }
                    "range_return" => {
                        let w_start = val["window_start_jd_ut1"].as_f64().expect("w_start");
                        let w_end = val["window_end_jd_ut1"].as_f64().expect("w_end");
                        let evs = eng_h4
                            .planetary_returns_in_range(
                                body,
                                TimeInput::from_ut1_jd(nat_jd),
                                TimeInput::from_ut1_jd(w_start),
                                TimeInput::from_ut1_jd(w_end),
                                opts,
                            )
                            .unwrap_or_else(|e| {
                                panic!("H4C[{idx}] planetary_returns_in_range failed: {e}")
                            });
                        let ev = evs
                            .get((ref_pass - 1) as usize)
                            .unwrap_or_else(|| panic!("H4C[{idx}] missing return pass {ref_pass}"));
                        assert_eq!(ev.pass_index, ref_pass);
                        assert_eq!(ev.is_retrograde, ref_retro);
                        ev.jd_ut1
                    }
                    other => panic!("Unknown H4C kind: {other}"),
                };

                let dt = (got_jd - ref_jd).abs();
                assert!(
                    dt < 1.0e-5,
                    "H4C[{idx}] {kind} {body:?} return diff={dt:.2e} days >= 1e-5 days"
                );
                max_return_days = max_return_days.max(dt);
            }

            // H4D: 200 Void-of-Course Moon intervals
            let voc_planets_6 = [
                BodyId::Planet(PlanetId::Sun),
                BodyId::Planet(PlanetId::Mercury),
                BodyId::Planet(PlanetId::Venus),
                BodyId::Planet(PlanetId::Mars),
                BodyId::Planet(PlanetId::Jupiter),
                BodyId::Planet(PlanetId::Saturn),
            ];
            let mut max_voc_days = 0.0_f64;
            for (idx, v) in sec_h4.void_of_course_moon_200.iter().enumerate() {
                let opts = CalcOptions {
                    zodiac: parse_zodiac(&v.zodiac),
                    ..CalcOptions::default()
                };
                let vocs = eng_h4
                    .void_of_course_moon_in_range(
                        TimeInput::from_ut1_jd(v.window_start_jd_ut1),
                        TimeInput::from_ut1_jd(v.window_end_jd_ut1),
                        Some(&voc_planets_6),
                        opts,
                    )
                    .unwrap_or_else(|e| panic!("H4D[{idx}] void_of_course_moon failed: {e}"));
                let voc = vocs
                    .iter()
                    .find(|item| {
                        item.moon_sign_before == v.current_sign_index
                            && item.moon_sign_after == v.next_sign_index
                    })
                    .unwrap_or_else(|| panic!("H4D[{idx}] VoC interval not found"));
                assert_eq!(voc.last_aspect_body, parse_planet_body(&v.last_aspect_body));
                assert_eq!(voc.last_aspect_name, v.last_aspect_name);
                let dt_s = (voc.start_jd_ut1 - v.start_jd_ut1).abs();
                let dt_e = (voc.end_jd_ut1 - v.end_jd_ut1).abs();
                max_voc_days = max_voc_days.max(dt_s).max(dt_e);
                assert!(
                    dt_s < 3.0e-4 && dt_e < 3.0e-4,
                    "H4D[{idx}] VoC diff: start={dt_s:.2e}, end={dt_e:.2e}"
                );
            }

            let rms_dyn_days = (sum_sq_dyn_days / 300.0).sqrt();
            (
                max_static_orb_err,
                max_dyn_aspect_days,
                rms_dyn_days,
                max_return_days,
                max_voc_days,
            )
        });

        // --------------------------------------------------------------------
        // Section H5: 800 HD Gates (384), Nakshatras (108) & ChartSnapshots (308)
        // --------------------------------------------------------------------
        let eng_h5 = Arc::clone(&engine);
        let sec_h5 = &data.section_h5_hd_nakshatra_snapshot;
        let t_h5 = s.spawn(move || {
            // H5A: 384 Human Design Rave Mandala vectors (all 64 Gates x 6 Lines)
            for (idx, v) in sec_h5.human_design_384.iter().enumerate() {
                assert_eq!(v.mandala_slot_index, idx);
                let pl = map_zodiac_placement(v.lon_deg);
                assert_eq!(pl.sign_index, v.sign_index);
                assert_eq!(pl.sign_name, v.sign_name);
                assert_eq!(pl.hd_gate, v.hd_gate, "H5A[{idx}] hd_gate");
                assert_eq!(pl.hd_line, v.hd_line, "H5A[{idx}] hd_line");
                assert_eq!(pl.hd_color, v.hd_color, "H5A[{idx}] hd_color");
                assert_eq!(pl.hd_tone, v.hd_tone, "H5A[{idx}] hd_tone");
                assert_eq!(pl.hd_base, v.hd_base, "H5A[{idx}] hd_base");
            }

            // H5B: 108 Nakshatra x Pada vectors (all 27 Nakshatras x 4 Padas + D9 bijection)
            for (idx, v) in sec_h5.nakshatra_padas_108.iter().enumerate() {
                assert_eq!(v.global_pada_index, idx);
                let pl = map_zodiac_placement(v.lon_deg);
                let d9 = eng_h5.calc_varga(v.lon_deg, 9).expect("D9 Navamsha");
                assert_eq!(pl.sign_index, v.sign_index);
                assert_eq!(pl.sign_name, v.sign_name);
                assert_eq!(pl.decan, v.decan);
                assert_eq!(pl.nakshatra_index, v.nakshatra_index);
                assert_eq!(pl.nakshatra_name, v.nakshatra_name);
                assert_eq!(pl.nakshatra_lord, v.nakshatra_lord);
                assert_eq!(pl.pada, v.pada);
                assert_eq!(d9.sign_index, v.d9_navamsa_sign_index);
                assert_eq!(d9.sign_name, v.d9_navamsa_sign_name);
                // Canonical Jyotish identity: global pada index (0..108) mod 12 == D9 Navamsha sign_index!
                assert_eq!((idx % 12) as u8, d9.sign_index);
            }

            // H5C: 308 ChartSnapshot records (28 chart summaries + 280 body snapshots)
            let bodies_10 = [
                BodyId::Planet(PlanetId::Sun),
                BodyId::Planet(PlanetId::Moon),
                BodyId::Planet(PlanetId::Mercury),
                BodyId::Planet(PlanetId::Venus),
                BodyId::Planet(PlanetId::Mars),
                BodyId::Planet(PlanetId::Jupiter),
                BodyId::Planet(PlanetId::Saturn),
                BodyId::Planet(PlanetId::Uranus),
                BodyId::Planet(PlanetId::Neptune),
                BodyId::Planet(PlanetId::Pluto),
            ];
            let mut snapshots_cache = HashMap::new();
            let mut max_snap_lon_err = 0.0_f64;
            let mut max_snap_hpos_err = 0.0_f64;

            for (idx, rec) in sec_h5.chart_snapshots_308.iter().enumerate() {
                let c_idx = rec["chart_index"].as_u64().expect("chart_index") as usize;
                let jd_ut1 = rec["jd_ut1"].as_f64().expect("jd_ut1");
                let lat_deg = rec["lat_deg"].as_f64().expect("lat_deg");
                let lon_deg = rec["lon_deg"].as_f64().expect("lon_deg");
                let alt_m = rec["alt_m"].as_f64().expect("alt_m");
                let hsys = rec["house_system"].as_str().expect("house_system");

                let snap = snapshots_cache.entry(c_idx).or_insert_with(|| {
                    let obs = Observer::new(lat_deg, lon_deg, alt_m).expect("obs");
                    let h_opts = HouseOptions {
                        system: parse_house_system(hsys),
                        ..HouseOptions::default()
                    };
                    eng_h5
                        .chart_snapshot(
                            TimeInput::from_ut1_jd(jd_ut1),
                            obs,
                            &bodies_10,
                            h_opts,
                            CalcOptions::default(),
                            AspectPreset::Major,
                            1.0,
                        )
                        .expect("chart_snapshot")
                });

                let rtype = rec["record_type"].as_str().expect("record_type");
                if rtype == "chart_summary" {
                    let ref_asc = rec["asc_deg"].as_f64().expect("asc_deg");
                    let ref_mc = rec["mc_deg"].as_f64().expect("mc_deg");
                    let ref_day = rec["is_day_chart"].as_bool().expect("is_day_chart");
                    let ref_fortune = rec["fortune_lon_deg"].as_f64().expect("fortune_lon_deg");
                    let ref_spirit = rec["spirit_lon_deg"].as_f64().expect("spirit_lon_deg");

                    assert!(wrap_pm180(snap.houses.angles.asc_deg - ref_asc).abs() < 0.005);
                    assert!(wrap_pm180(snap.houses.angles.mc_deg - ref_mc).abs() < 0.005);
                    assert_eq!(snap.sect.is_day_chart, ref_day);
                    assert!(wrap_pm180(snap.lots[0].lon_deg - ref_fortune).abs() < 0.005);
                    assert!(wrap_pm180(snap.lots[1].lon_deg - ref_spirit).abs() < 0.005);
                } else {
                    let bname = rec["body"].as_str().expect("body");
                    let target_body = parse_planet_body(bname);
                    let b_snap = snap
                        .bodies
                        .iter()
                        .find(|b| b.body == target_body)
                        .expect("body in snapshot");

                    let ref_lon = rec["body_lon_deg"].as_f64().expect("body_lon_deg");
                    let ref_retro = rec["is_retrograde"].as_bool().expect("is_retrograde");
                    let ref_sign = rec["sign_index"].as_u64().expect("sign_index") as u8;
                    let ref_decan = rec["decan"].as_u64().expect("decan") as u8;
                    let ref_nak = rec["nakshatra_index"].as_u64().expect("nakshatra_index") as u8;
                    let ref_pada = rec["pada"].as_u64().expect("pada") as u8;
                    let ref_gate = rec["hd_gate"].as_u64().expect("hd_gate") as u8;
                    let ref_line = rec["hd_line"].as_u64().expect("hd_line") as u8;
                    let ref_hnum = rec["house_number"].as_u64().expect("house_number") as u8;
                    let ref_hpos = rec["house_pos_3d"].as_f64().expect("house_pos_3d");

                    let lon_err = wrap_pm180(b_snap.calc.lon_deg - ref_lon).abs();
                    max_snap_lon_err = max_snap_lon_err.max(lon_err);
                    assert!(lon_err < 0.005, "H5C[{idx}] {bname} lon_err={lon_err}");
                    assert_eq!(b_snap.is_retrograde, ref_retro, "H5C[{idx}] {bname} retro");
                    assert_eq!(b_snap.placement.sign_index, ref_sign);
                    assert_eq!(b_snap.placement.decan, ref_decan);
                    assert_eq!(b_snap.placement.nakshatra_index, ref_nak);
                    assert_eq!(b_snap.placement.pada, ref_pada);
                    assert_eq!(b_snap.placement.hd_gate, ref_gate);
                    assert_eq!(b_snap.placement.hd_line, ref_line);
                    assert_eq!(b_snap.house_number, ref_hnum);

                    let hpos_err = (b_snap.house_pos_3d - ref_hpos).abs();
                    let hpos_circ_err = hpos_err.min((12.0 - hpos_err).abs());
                    max_snap_hpos_err = max_snap_hpos_err.max(hpos_circ_err);
                    assert!(
                        hpos_circ_err < 0.005,
                        "H5C[{idx}] {bname} house_pos_3d diff={hpos_circ_err}"
                    );
                }
            }
            assert_eq!(snapshots_cache.len(), 28);
            (max_snap_lon_err, max_snap_hpos_err)
        });

        (
            t_h1.join().expect("H1 thread"),
            t_h2.join().expect("H2 thread"),
            t_h3.join().expect("H3 thread"),
            t_h4.join().expect("H4 thread"),
            t_h5.join().expect("H5 thread"),
        )
    });

    println!(
        "\n================================================================================"
    );
    println!("Wave H5 5,000-Vector Golden Oracle Verification Summary (100% Pass)");
    println!(
        "================================================================================"
    );
    println!(
        "  H1 (1,000 Panchanga & Muhurta):   max_tithi_karana_err = {:.2e} d ({:.2} s), rms = {:.2e} d, max_sunrise_err = {:.2e} d",
        h1_metrics.0,
        h1_metrics.0 * 86400.0,
        h1_metrics.1,
        h1_metrics.2
    );
    println!(
        "  H2 (1,200 Vargas & Vimshottari):  max_varga_lon_err = {:.2e} deg (rms = {:.2e} deg), max_dasha_formula_err = {:.2e} d, max_dasha_120yr_de440s_err = {:.4} d",
        h2_metrics.0, h2_metrics.1, h2_metrics.2, h2_metrics.3
    );
    println!(
        "  H3 (1,000 Sect, Hours & 7 Lots):  max_lot_lon_err = {:.5} deg (rms = {:.5} deg), max_sun_alt_err = {:.5} deg",
        h3_metrics.0, h3_metrics.1, h3_metrics.2
    );
    println!(
        "  H4 (1,000 Aspects, Returns, VoC): max_dyn_aspect_err = {:.2e} d ({:.2} s, rms = {:.2e} d), max_return_err = {:.2e} d, max_voc_err = {:.2e} d, max_static_orb_err = {:.5} deg",
        h4_metrics.1,
        h4_metrics.1 * 86400.0,
        h4_metrics.2,
        h4_metrics.3,
        h4_metrics.4,
        h4_metrics.0
    );
    println!(
        "  H5 (  800 HD, Nakshatra, Snap):   64x6 HD Gates/Lines 100% exact, 27x4 Nakshatras/Padas 100% exact, max_snap_lon_err = {:.5} deg, max_hpos_3d_err = {:.5}",
        h5_metrics.0, h5_metrics.1
    );
    println!(
        "================================================================================\n"
    );
}

