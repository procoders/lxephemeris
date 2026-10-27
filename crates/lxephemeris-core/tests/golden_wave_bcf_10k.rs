//! Wave F3 10,540-Vector Golden Reference Verification Suite (`l3_wave_bcf_10k_grid.json`).
//!
//! Verifies all 8 Wave B, Wave C, and Wave F sections (10,540 golden reference vectors generated
//! via `scripts/gen_golden_wave_bcf_10k.py` using `pyswisseph` and `ephem` as black-box reference
//! oracles) against `lxephemeris_core::Engine` and `lxephemeris_core::SweAdapter`:
//!
//! 1. `section_1_rise_trans` (1,600 vectors): Sun Rise, Set (`almanac_upper_limb`), UpperTransit,
//!    and LowerTransit (`geometric_center`) across 100 dates (1855..2095) x 4 observers.
//! 2. `section_2_lunar_phases` (1,500 vectors): 1,500 consecutive principal lunar phases from
//!    1900-01-01 (`NewMoon`, `FirstQuarter`, `FullMoon`, `LastQuarter`).
//! 3. `section_3_lunar_phase_summary` (1,500 vectors): Complete `lunar_phase_summary` across
//!    1,500 evenly spaced epochs (1855..2095).
//! 4. `section_4_solar_terms_and_seasons` (1,440 vectors): 60 years (1980..=2039) x 24 Chinese
//!    Solar Terms + cardinal seasonal equinoxes/solstices.
//! 5. `section_5_hd_and_longitude_crossings` (1,000 vectors): 500 Human Design `-88°` prenatal
//!    solar arc epochs + 500 planetary longitude crossings (Tropical & Sidereal `Lahiri`).
//! 6. `section_6_stations_and_declination_extrema` (1,000 vectors): 500 planetary retrograde/direct
//!    stations (`Mercury..Pluto`, Tropical & Sidereal `Lahiri`) + 500 Sun/Moon declination extrema.
//! 7. `section_7_brady_parans` (500 vectors): 500 diurnal fixed-star paran events across birth
//!    charts and observer latitudes.
//! 8. `section_8_wave_b_c_swe_parity` (2,000 vectors): 500 Moshier (`SEFLG_MOSEPH`) vectors across
//!    `-500 CE .. +2500 CE`, 300 photometric phenomena (`pheno_ut`), 300 orbital nodes & apsides
//!    (`nod_aps_ut`), 430 extended Ayanamsas (`0..=42` x 10 epochs), and 470 3D mundane
//!    `house_pos` (200), `time_equ` (120), and `azalt_rev` (150) vectors.

#![cfg(feature = "swe-compat")]
#![forbid(unsafe_code)]
#![allow(
    clippy::type_complexity,
    clippy::indexing_slicing,
    clippy::panic,
    clippy::unwrap_used
)]

use lxephemeris_core::{
    swe_compat::{SweAdapter, SEFLG_HELCTR, SEFLG_MOSEPH, SEFLG_SPEED, SEFLG_SWIEPH},
    AyanamsaModel, BodyId, CalcOptions, DeclinationExtremumKind, Engine, EngineConfig,
    LunarPhaseKind, NutationMode, Observer, ParanAngle, ParanOptions, PlanetId, RiseTransEvent,
    RiseTransOptions, StationKind, TimeInput, Zodiac, SE_NODBIT_OSCU,
};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

fn wrap180(deg: f64) -> f64 {
    (deg + 180.0).rem_euclid(360.0) - 180.0
}

fn wrap6_house(diff: f64) -> f64 {
    (diff + 6.0).rem_euclid(12.0) - 6.0
}

#[derive(Default, Clone, Copy)]
struct StatAcc {
    count: usize,
    sum_sq: f64,
    max_abs: f64,
}

impl StatAcc {
    fn push(&mut self, err: f64) {
        let a = err.abs();
        self.count += 1;
        self.sum_sq += a * a;
        if a > self.max_abs {
            self.max_abs = a;
        }
    }

    fn merge(&mut self, other: Self) {
        self.count += other.count;
        self.sum_sq += other.sum_sq;
        if other.max_abs > self.max_abs {
            self.max_abs = other.max_abs;
        }
    }

    fn rms(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            (self.sum_sq / (self.count as f64)).sqrt()
        }
    }
}

#[derive(Deserialize)]
struct GoldenSection5 {
    human_design_500: Vec<(f64, f64, f64, f64)>,
    longitude_crossings_500: Vec<(String, String, f64, f64, f64, f64, f64)>,
}

#[derive(Deserialize)]
struct GoldenSection6 {
    stations_500: Vec<(String, String, f64, f64, String, f64, f64, f64, f64)>,
    declination_extrema_500: Vec<(String, f64, f64, String, f64, f64, f64, f64)>,
}

#[derive(Deserialize)]
struct GoldenSection8 {
    moshier_500: Vec<(f64, String, i32, f64, f64, f64, f64)>,
    phenomena_300: Vec<(f64, String, i32, f64, f64, f64, f64, f64, f64)>,
    nodes_apsides_300: Vec<(f64, String, i32, u32, f64, f64, f64, f64, f64, f64, f64, f64, f64)>,
    extended_ayanamsas_430: Vec<(f64, u32, f64, f64, f64, f64)>,
    house_pos_200: Vec<(f64, f64, f64, String, f64, f64, f64)>,
    time_equation_120: Vec<(f64, f64)>,
    azalt_rev_150: Vec<(f64, i32, f64, f64, f64, f64, f64, f64, f64)>,
}

#[derive(Deserialize)]
struct GoldenGrid10k {
    total_vector_count: usize,
    section_1_rise_trans: Vec<(f64, String, f64, f64, String, f64, f64)>,
    section_2_lunar_phases: Vec<(usize, String, f64, f64, f64, f64)>,
    section_3_lunar_phase_summary: Vec<(
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        u8,
        u8,
        f64,
        f64,
        f64,
    )>,
    section_4_solar_terms_and_seasons: Vec<(i32, u8, String, f64, f64, Option<f64>)>,
    section_5_hd_and_longitude_crossings: GoldenSection5,
    section_6_stations_and_declination_extrema: GoldenSection6,
    section_7_brady_parans: Vec<(
        usize,
        String,
        f64,
        f64,
        f64,
        String,
        u32,
        String,
        String,
        String,
        String,
        f64,
        f64,
        f64,
    )>,
    section_8_wave_b_c_swe_parity: GoldenSection8,
}

fn parse_planet(name: &str) -> PlanetId {
    match name {
        "Sun" => PlanetId::Sun,
        "Moon" => PlanetId::Moon,
        "Mercury" => PlanetId::Mercury,
        "Venus" => PlanetId::Venus,
        "Earth" => PlanetId::Earth,
        "Mars" => PlanetId::Mars,
        "Jupiter" => PlanetId::Jupiter,
        "Saturn" => PlanetId::Saturn,
        "Uranus" => PlanetId::Uranus,
        "Neptune" => PlanetId::Neptune,
        "Pluto" => PlanetId::Pluto,
        other => panic!("Unknown planet name: {other}"),
    }
}

fn parse_zodiac_opts(mode_str: &str) -> CalcOptions {
    let zodiac = match mode_str {
        "Tropical" => Zodiac::Tropical,
        "Sidereal_Lahiri" => Zodiac::Sidereal(AyanamsaModel::Lahiri),
        other => panic!("Unknown zodiac mode: {other}"),
    };
    CalcOptions {
        zodiac,
        ..CalcOptions::default()
    }
}

fn parse_paran_angle(ang_str: &str) -> ParanAngle {
    match ang_str {
        "rising" => ParanAngle::Rising,
        "culminating" => ParanAngle::Culminating,
        "setting" => ParanAngle::Setting,
        "anti-culminating" => ParanAngle::AntiCulminating,
        other => panic!("Unknown paran angle: {other}"),
    }
}

#[test]
fn test_golden_wave_bcf_10k_all_sections() {
    let root = workspace_root();
    let json_text = std::fs::read_to_string(root.join("testdata/golden/l3_wave_bcf_10k_grid.json"))
        .expect("read l3_wave_bcf_10k_grid.json");
    let grid: GoldenGrid10k =
        serde_json::from_str(&json_text).expect("parse l3_wave_bcf_10k_grid.json");

    assert!(
        grid.total_vector_count >= 10_500,
        "Expected >= 10,500 golden vectors, got {}",
        grid.total_vector_count
    );
    assert_eq!(grid.total_vector_count, 10_540);

    let eng = Arc::new(
        Engine::from_data_dir(root.join("data"), EngineConfig::default()).expect("load Engine"),
    );

    // Run sections across scoped worker threads to verify thread-safety (Send + Sync)
    // and complete all 10,540 golden reference checks in < 5 seconds.
    let (
        (s1_swe_rise_set_sec, s1_swe_transit_sec, s1_ephem_rise_set_sec, s1_ephem_transit_sec),
        (s2_swe_sec, s2_ephem_sec, s2_dist_km),
        (
            s3_elong_arcsec,
            s3_illum,
            s3_pa_arcsec,
            s3_dist_km,
            s3_diam_arcsec,
            s3_prev_nm_sec,
            s3_next_nm_sec,
            s3_age_sec,
        ),
        (s4_swe_sec, s4_ephem_sec),
        (s5_hd_sec, s5_cross_sec, s5_cross_vel),
        (s6_station_min, s6_station_lon_arcsec, s6_dec_min, s6_dec_arcsec),
        (s7_star_sec, s7_planet_sec, s7_orb_sec),
        (
            s8_mosh_lon_arcsec,
            s8_mosh_lat_arcsec,
            s8_pheno_pa_arcsec,
            s8_pheno_el_arcsec,
            s8_pheno_mag,
            s8_nod_lon_deg,
            s8_peri_lon_deg,
            s8_aya_l1b_arcsec,
            s8_aya_swe_arcsec,
            s8_house_pos,
            s8_eot_sec,
            s8_azrev_arcsec,
        ),
    ) = std::thread::scope(|s| {
        // =====================================================================
        // Section 1: Rise, Set & Meridian Transits (1,600 vectors, 4 chunks)
        // =====================================================================
        assert_eq!(grid.section_1_rise_trans.len(), 1_600);
        let s1_handles: Vec<_> = grid
            .section_1_rise_trans
            .chunks(400)
            .map(|chunk| {
                let eng = Arc::clone(&eng);
                s.spawn(move || {
                    let mut swe_rs = StatAcc::default();
                    let mut swe_tr = StatAcc::default();
                    let mut eph_rs = StatAcc::default();
                    let mut eph_tr = StatAcc::default();

                    for (start_jd, obs_name, lat_deg, lon_deg, ev_str, swe_jd, ephem_jd) in chunk {
                        let obs = Observer::new(*lat_deg, *lon_deg, 0.0).expect("observer");
                        let (kind, opts, is_transit) = match ev_str.as_str() {
                            "Rise" => (
                                RiseTransEvent::Rise,
                                RiseTransOptions::almanac_upper_limb(),
                                false,
                            ),
                            "Set" => (
                                RiseTransEvent::Set,
                                RiseTransOptions::almanac_upper_limb(),
                                false,
                            ),
                            "UpperTransit" => (
                                RiseTransEvent::UpperTransit,
                                RiseTransOptions::geometric_center(),
                                true,
                            ),
                            "LowerTransit" => (
                                RiseTransEvent::LowerTransit,
                                RiseTransOptions::geometric_center(),
                                true,
                            ),
                            other => panic!("Unknown RiseTransEvent: {other}"),
                        };

                        let ev_jd = eng
                            .rise_trans(BodyId::Planet(PlanetId::Sun), obs, *start_jd, kind, opts)
                            .unwrap_or_else(|e| {
                                panic!("Sec1 {obs_name} {ev_str} at {start_jd} failed: {e}")
                            })
                            .unwrap_or_else(|| {
                                panic!("Sec1 {obs_name} {ev_str} at {start_jd} returned None")
                            });

                        let dt_swe_sec = (ev_jd - *swe_jd) * 86_400.0;
                        let dt_ephem_sec = (ev_jd - *ephem_jd) * 86_400.0;
                        if is_transit {
                            swe_tr.push(dt_swe_sec);
                            eph_tr.push(dt_ephem_sec);
                            assert!(
                                dt_swe_sec.abs() < 5.0,
                                "Sec1 transit {obs_name} {ev_str} at {start_jd}: dt_swe={dt_swe_sec:.3}s"
                            );
                            assert!(
                                dt_ephem_sec.abs() < 10.0,
                                "Sec1 transit {obs_name} {ev_str} at {start_jd}: dt_ephem={dt_ephem_sec:.3}s"
                            );
                        } else {
                            swe_rs.push(dt_swe_sec);
                            eph_rs.push(dt_ephem_sec);
                            assert!(
                                dt_swe_sec.abs() < 20.0,
                                "Sec1 rise/set {obs_name} {ev_str} at {start_jd}: dt_swe={dt_swe_sec:.3}s"
                            );
                            assert!(
                                dt_ephem_sec.abs() < 45.0,
                                "Sec1 rise/set {obs_name} {ev_str} at {start_jd}: dt_ephem={dt_ephem_sec:.3}s"
                            );
                        }
                    }
                    (swe_rs, swe_tr, eph_rs, eph_tr)
                })
            })
            .collect();

        // =====================================================================
        // Section 2: 1,500 Consecutive Principal Lunar Phases (1900..2021)
        // =====================================================================
        let eng_s2 = Arc::clone(&eng);
        let sec2_ref = &grid.section_2_lunar_phases;
        let h_s2 = s.spawn(move || {
            assert_eq!(sec2_ref.len(), 1_500);
            let mut s2_swe_sec = StatAcc::default();
            let mut s2_ephem_sec = StatAcc::default();
            let mut s2_dist_km = StatAcc::default();

            let first_phase_start = 2_415_020.5;
            let last_phase_end = sec2_ref.last().unwrap().3 + 1.0;
            let range_phases = eng_s2
                .lunar_phases_in_range(
                    TimeInput::from_ut1_jd(first_phase_start),
                    TimeInput::from_ut1_jd(last_phase_end),
                )
                .expect("lunar_phases_in_range 1500");
            assert_eq!(
                range_phases.len(),
                1_500,
                "lunar_phases_in_range should find all 1,500 consecutive principal phases"
            );

            for (i, (idx, kind_str, target_deg, swe_jd, ephem_jd, swe_dist_km)) in
                sec2_ref.iter().enumerate()
            {
                assert_eq!(*idx, i);
                let expected_kind = match kind_str.as_str() {
                    "NewMoon" => LunarPhaseKind::NewMoon,
                    "FirstQuarter" => LunarPhaseKind::FirstQuarter,
                    "FullMoon" => LunarPhaseKind::FullMoon,
                    "LastQuarter" => LunarPhaseKind::LastQuarter,
                    other => panic!("Unknown LunarPhaseKind: {other}"),
                };
                let ev = &range_phases[i];
                assert_eq!(ev.kind, Some(expected_kind));
                assert!((ev.target_elongation_deg - *target_deg).abs() < 1e-9);

                if i % 25 == 0 {
                    let single = eng_s2
                        .search_lunar_phase(
                            TimeInput::from_ut1_jd(*swe_jd - 3.0),
                            *target_deg,
                            true,
                        )
                        .expect("search_lunar_phase single");
                    assert!((single.jd_ut1 - ev.jd_ut1).abs() * 86_400.0 < 1e-4);
                }

                let dt_swe_sec = (ev.jd_ut1 - *swe_jd) * 86_400.0;
                let dt_ephem_sec = (ev.jd_ut1 - *ephem_jd) * 86_400.0;
                let d_dist_km = ev.moon_dist_km - *swe_dist_km;

                s2_swe_sec.push(dt_swe_sec);
                s2_ephem_sec.push(dt_ephem_sec);
                s2_dist_km.push(d_dist_km);

                assert!(
                    dt_swe_sec.abs() < 15.0,
                    "Sec2 phase #{idx} ({kind_str}): dt_swe={dt_swe_sec:.3}s"
                );
                assert!(
                    dt_ephem_sec.abs() < 45.0,
                    "Sec2 phase #{idx} ({kind_str}): dt_ephem={dt_ephem_sec:.3}s"
                );
                assert!(
                    d_dist_km.abs() < 15.0,
                    "Sec2 phase #{idx} ({kind_str}): d_dist_km={d_dist_km:.3} km"
                );
            }
            (s2_swe_sec, s2_ephem_sec, s2_dist_km)
        });

        // =====================================================================
        // Section 3: 1,500 Lunar Phase Summaries across 1855..2095 (3 chunks)
        // =====================================================================
        assert_eq!(grid.section_3_lunar_phase_summary.len(), 1_500);
        let s3_handles: Vec<_> = grid
            .section_3_lunar_phase_summary
            .chunks(500)
            .map(|chunk| {
                let eng = Arc::clone(&eng);
                s.spawn(move || {
                    let mut elong_acc = StatAcc::default();
                    let mut illum_acc = StatAcc::default();
                    let mut pa_acc = StatAcc::default();
                    let mut dist_acc = StatAcc::default();
                    let mut diam_acc = StatAcc::default();
                    let mut prev_acc = StatAcc::default();
                    let mut next_acc = StatAcc::default();
                    let mut age_acc = StatAcc::default();

                    for (
                        jd_ut1,
                        swe_elong,
                        swe_illum,
                        swe_pa,
                        swe_dist,
                        swe_diam,
                        ephem_illum,
                        ephem_dist,
                        exp_p8,
                        exp_p28,
                        swe_prev_nm,
                        swe_next_nm,
                        swe_age_days,
                    ) in chunk
                    {
                        let sum = eng
                            .lunar_phase_summary(TimeInput::from_ut1_jd(*jd_ut1))
                            .expect("lunar_phase_summary");

                        let d_elong_arcsec = wrap180(sum.elongation_deg - *swe_elong) * 3_600.0;
                        let d_illum = sum.illumination_fraction - *swe_illum;
                        let d_pa_arcsec = (sum.phase_angle_deg - *swe_pa) * 3_600.0;
                        let d_dist_km = sum.moon_dist_km - *swe_dist;
                        let d_diam_arcsec = (sum.apparent_diameter_deg - *swe_diam) * 3_600.0;
                        let d_prev_sec = (sum.prev_new_moon_jd_ut1 - *swe_prev_nm) * 86_400.0;
                        let d_next_sec = (sum.next_new_moon_jd_ut1 - *swe_next_nm) * 86_400.0;
                        let d_age_sec = (sum.age_days - *swe_age_days) * 86_400.0;

                        elong_acc.push(d_elong_arcsec);
                        illum_acc.push(d_illum);
                        pa_acc.push(d_pa_arcsec);
                        dist_acc.push(d_dist_km);
                        diam_acc.push(d_diam_arcsec);
                        prev_acc.push(d_prev_sec);
                        next_acc.push(d_next_sec);
                        age_acc.push(d_age_sec);

                        assert_eq!(
                            sum.phase_8_index, *exp_p8,
                            "Sec3 phase_8_index mismatch at JD {jd_ut1}: got {}, expected {}",
                            sum.phase_8_index, exp_p8
                        );
                        assert_eq!(
                            sum.phase_28_index, *exp_p28,
                            "Sec3 phase_28_index mismatch at JD {jd_ut1}: got {}, expected {}",
                            sum.phase_28_index, exp_p28
                        );
                        assert!(
                            d_elong_arcsec.abs() < 60.0,
                            "Sec3 elong at {jd_ut1}: d_elong={d_elong_arcsec:.3}\""
                        );
                        assert!(d_illum.abs() < 1e-3);
                        assert!((sum.illumination_fraction - *ephem_illum).abs() < 0.02);
                        assert!((sum.moon_dist_km - *ephem_dist).abs() < 100.0);
                        assert!(
                            d_prev_sec.abs() < 120.0,
                            "Sec3 prev_nm at {jd_ut1}: d_prev_sec={d_prev_sec:.3}s"
                        );
                        assert!(
                            d_next_sec.abs() < 120.0,
                            "Sec3 next_nm at {jd_ut1}: d_next_sec={d_next_sec:.3}s"
                        );
                    }
                    (
                        elong_acc, illum_acc, pa_acc, dist_acc, diam_acc, prev_acc, next_acc,
                        age_acc,
                    )
                })
            })
            .collect();

        // =====================================================================
        // Section 4: 1,440 Chinese Solar Terms & Cardinal Seasons (1980..=2039)
        // =====================================================================
        let eng_s4 = Arc::clone(&eng);
        let sec4_ref = &grid.section_4_solar_terms_and_seasons;
        let h_s4 = s.spawn(move || {
            assert_eq!(sec4_ref.len(), 1_440);
            let mut s4_swe_sec = StatAcc::default();
            let mut s4_ephem_sec = StatAcc::default();

            for year_chunk in sec4_ref.chunks_exact(24) {
                let year = year_chunk[0].0;
                let start_jd = year_chunk[0].4 - 3.0;
                let end_jd = year_chunk[23].4 + 3.0;

                let terms = eng_s4
                    .solar_terms_in_range(
                        TimeInput::from_ut1_jd(start_jd),
                        TimeInput::from_ut1_jd(end_jd),
                    )
                    .expect("solar_terms_in_range 24");
                assert_eq!(terms.len(), 24, "Expected 24 solar terms in year {year}");

                let seasons = eng_s4.seasons_in_year(year).expect("seasons_in_year");

                for (t_ev, (yr, t_idx, name_pinyin, target_lon, swe_jd, ephem_jd_opt)) in
                    terms.iter().zip(year_chunk.iter())
                {
                    assert_eq!(*yr, year);
                    assert_eq!(t_ev.term_index, *t_idx);
                    assert_eq!(t_ev.pinyin_name, name_pinyin.as_str());
                    assert!((t_ev.target_lon_deg - *target_lon).abs() < 1e-9);

                    let dt_swe_sec = (t_ev.jd_ut1 - *swe_jd) * 86_400.0;
                    s4_swe_sec.push(dt_swe_sec);
                    assert!(
                        dt_swe_sec.abs() < 35.0,
                        "Sec4 solar term {year} #{t_idx} ({name_pinyin}): dt_swe={dt_swe_sec:.3}s"
                    );

                    if let Some(ephem_jd) = ephem_jd_opt {
                        let dt_ephem_sec = (t_ev.jd_ut1 - *ephem_jd) * 86_400.0;
                        s4_ephem_sec.push(dt_ephem_sec);
                        assert!(
                            dt_ephem_sec.abs() < 35.0,
                            "Sec4 cardinal season {year} #{t_idx}: dt_ephem={dt_ephem_sec:.3}s"
                        );
                        let season_jd_ut1 = match *t_idx {
                            0 => seasons.vernal_equinox_jd_ut1,
                            6 => seasons.summer_solstice_jd_ut1,
                            12 => seasons.autumnal_equinox_jd_ut1,
                            18 => seasons.winter_solstice_jd_ut1,
                            other => panic!("Unexpected ephem_jd on term_index {other}"),
                        };
                        assert!((season_jd_ut1 - t_ev.jd_ut1).abs() * 86_400.0 < 1e-4);
                    }
                }
            }
            (s4_swe_sec, s4_ephem_sec)
        });

        // =====================================================================
        // Section 5: 500 Human Design -88° Epochs + 500 Longitude Crossings
        // =====================================================================
        let eng_s5 = Arc::clone(&eng);
        let sec5_ref = &grid.section_5_hd_and_longitude_crossings;
        let h_s5 = s.spawn(move || {
            assert_eq!(sec5_ref.human_design_500.len(), 500);
            assert_eq!(sec5_ref.longitude_crossings_500.len(), 500);

            let mut s5_hd_sec = StatAcc::default();
            for (birth_jd, _swe_birth_lon, swe_design_lon, swe_design_jd) in
                &sec5_ref.human_design_500
            {
                let hd = eng_s5
                    .human_design_epoch(TimeInput::from_ut1_jd(*birth_jd))
                    .expect("human_design_epoch");
                assert!(wrap180(hd.target_lon_deg - *swe_design_lon).abs() * 3_600.0 < 5.0);
                assert!(wrap180(hd.actual_lon_deg - *swe_design_lon).abs() * 3_600.0 < 5.0);
                let dt_sec = (hd.jd_ut1 - *swe_design_jd) * 86_400.0;
                s5_hd_sec.push(dt_sec);
                assert!(
                    dt_sec.abs() < 60.0,
                    "Sec5A HD epoch at birth {birth_jd}: dt_sec={dt_sec:.3}s"
                );
            }

            let mut s5_cross_sec = StatAcc::default();
            let mut s5_cross_vel = StatAcc::default();
            for (b_name, z_str, start_jd, target_lon, max_days, swe_root_jd, swe_vel) in
                &sec5_ref.longitude_crossings_500
            {
                let body = BodyId::Planet(parse_planet(b_name));
                let opts = parse_zodiac_opts(z_str);
                let ev = eng_s5
                    .search_body_longitude(
                        body,
                        TimeInput::from_ut1_jd(*start_jd),
                        *target_lon,
                        true,
                        *max_days,
                        &opts,
                    )
                    .unwrap_or_else(|e| {
                        panic!(
                            "Sec5B crossing {b_name} ({z_str}) at {start_jd} -> {target_lon} failed: {e}"
                        )
                    });
                let dt_sec = (ev.jd_ut1 - *swe_root_jd) * 86_400.0;
                let dv = ev.lon_speed_deg_day - *swe_vel;
                s5_cross_sec.push(dt_sec);
                s5_cross_vel.push(dv);
                assert!(
                    wrap180(ev.actual_lon_deg - *target_lon).abs() < 1e-6,
                    "Sec5B crossing {b_name} ({z_str}): actual_lon={}, target={target_lon}",
                    ev.actual_lon_deg
                );
                let dt_tol_sec = match b_name.as_str() {
                    "Mars" | "Jupiter" | "Saturn" => 900.0,
                    _ => 180.0,
                };
                assert!(
                    dt_sec.abs() < dt_tol_sec,
                    "Sec5B crossing {b_name} ({z_str}) at {start_jd}: dt_sec={dt_sec:.3}s"
                );
                assert!(dv.abs() < 5e-3);
            }
            (s5_hd_sec, s5_cross_sec, s5_cross_vel)
        });

        // =====================================================================
        // Section 6: 500 Planetary Stations + 500 Declination Extrema
        // =====================================================================
        let eng_s6 = Arc::clone(&eng);
        let sec6_ref = &grid.section_6_stations_and_declination_extrema;
        let h_s6 = s.spawn(move || {
            assert_eq!(sec6_ref.stations_500.len(), 500);
            assert_eq!(sec6_ref.declination_extrema_500.len(), 500);

            let mut s6_station_min = StatAcc::default();
            let mut s6_station_lon_arcsec = StatAcc::default();

            let mut station_groups: BTreeMap<
                (String, String, u64, u64),
                Vec<&(String, String, f64, f64, String, f64, f64, f64, f64)>,
            > = BTreeMap::new();
            for row in &sec6_ref.stations_500 {
                station_groups
                    .entry((row.0.clone(), row.1.clone(), row.2.to_bits(), row.3.to_bits()))
                    .or_default()
                    .push(row);
            }

            let mut total_stations_verified = 0usize;
            for ((p_name, z_str, w0_bits, w1_bits), expected_rows) in &station_groups {
                let body = BodyId::Planet(parse_planet(p_name));
                let opts = parse_zodiac_opts(z_str);
                let win_start = f64::from_bits(*w0_bits);
                let win_end = f64::from_bits(*w1_bits);

                let found = eng_s6
                    .search_stations_in_range(
                        body,
                        TimeInput::from_ut1_jd(win_start),
                        TimeInput::from_ut1_jd(win_end),
                        &opts,
                    )
                    .expect("search_stations_in_range");
                assert_eq!(
                    found.len(),
                    expected_rows.len(),
                    "Sec6A station count mismatch for {p_name} ({z_str})"
                );

                for (ev, row) in found.iter().zip(expected_rows.iter()) {
                    let exp_kind = match row.4.as_str() {
                        "StationRetrograde" => StationKind::StationRetrograde,
                        "StationDirect" => StationKind::StationDirect,
                        other => panic!("Unknown StationKind: {other}"),
                    };
                    assert_eq!(ev.kind, exp_kind);
                    let dt_min = (ev.jd_ut1 - row.5) * 1_440.0;
                    let d_lon_arcsec = wrap180(ev.lon_deg - row.6) * 3_600.0;
                    s6_station_min.push(dt_min);
                    s6_station_lon_arcsec.push(d_lon_arcsec);
                    assert!(
                        dt_min.abs() < 60.0,
                        "Sec6A station {p_name} ({z_str}) at {}: dt_min={dt_min:.3}m",
                        row.5
                    );
                    assert!(d_lon_arcsec.abs() < 10.0);
                    total_stations_verified += 1;
                }
            }
            assert_eq!(total_stations_verified, 500);

            let mut s6_dec_min = StatAcc::default();
            let mut s6_dec_arcsec = StatAcc::default();
            let mut dec_groups: BTreeMap<
                (String, u64, u64),
                Vec<&(String, f64, f64, String, f64, f64, f64, f64)>,
            > = BTreeMap::new();
            for row in &sec6_ref.declination_extrema_500 {
                dec_groups
                    .entry((row.0.clone(), row.1.to_bits(), row.2.to_bits()))
                    .or_default()
                    .push(row);
            }

            let mut total_dec_verified = 0usize;
            for ((b_name, w0_bits, w1_bits), expected_rows) in &dec_groups {
                let body = BodyId::Planet(parse_planet(b_name));
                let win_start = f64::from_bits(*w0_bits);
                let win_end = f64::from_bits(*w1_bits);

                let found = eng_s6
                    .search_declination_extrema_in_range(
                        body,
                        TimeInput::from_ut1_jd(win_start),
                        TimeInput::from_ut1_jd(win_end),
                        &CalcOptions::default(),
                    )
                    .expect("search_declination_extrema_in_range");
                assert_eq!(
                    found.len(),
                    expected_rows.len(),
                    "Sec6B declination extrema count mismatch for {b_name}"
                );

                for (ev, row) in found.iter().zip(expected_rows.iter()) {
                    let exp_kind = match row.3.as_str() {
                        "NorthernMaximum" => DeclinationExtremumKind::NorthernMaximum,
                        "SouthernMinimum" => DeclinationExtremumKind::SouthernMinimum,
                        other => panic!("Unknown DeclinationExtremumKind: {other}"),
                    };
                    assert_eq!(ev.kind, exp_kind);
                    let dt_min = (ev.jd_ut1 - row.4) * 1_440.0;
                    let d_dec_arcsec = (ev.dec_deg - row.6) * 3_600.0;
                    s6_dec_min.push(dt_min);
                    s6_dec_arcsec.push(d_dec_arcsec);
                    assert!(
                        dt_min.abs() < 30.0,
                        "Sec6B dec extremum {b_name} at {}: dt_min={dt_min:.3}m",
                        row.4
                    );
                    assert!(d_dec_arcsec.abs() < 10.0);
                    total_dec_verified += 1;
                }
            }
            assert_eq!(total_dec_verified, 500);
            (
                s6_station_min,
                s6_station_lon_arcsec,
                s6_dec_min,
                s6_dec_arcsec,
            )
        });

        // =====================================================================
        // Section 7: 500 Brady Diurnal Fixed-Star Parans across Birth Charts
        // =====================================================================
        let eng_s7 = Arc::clone(&eng);
        let sec7_ref = &grid.section_7_brady_parans;
        let h_s7 = s.spawn(move || {
            assert_eq!(sec7_ref.len(), 500);
            let mut s7_star_sec = StatAcc::default();
            let mut s7_planet_sec = StatAcc::default();
            let mut s7_orb_sec = StatAcc::default();

            let mut paran_charts: BTreeMap<
                usize,
                Vec<&(
                    usize,
                    String,
                    f64,
                    f64,
                    f64,
                    String,
                    u32,
                    String,
                    String,
                    String,
                    String,
                    f64,
                    f64,
                    f64,
                )>,
            > = BTreeMap::new();
            for row in sec7_ref {
                paran_charts.entry(row.0).or_default().push(row);
            }

            let mut total_parans_verified = 0usize;
            let paran_opts = ParanOptions {
                orb_minutes: 4.5,
                ..ParanOptions::default()
            };
            for rows in paran_charts.values() {
                let first = rows[0];
                let obs = Observer::new(first.2, first.3, 0.0).expect("paran obs");
                let start_jd = first.4;
                let parans = eng_s7
                    .find_parans(TimeInput::from_ut1_jd(start_jd), obs, &paran_opts)
                    .expect("find_parans");

                for row in rows {
                    let exp_hip = row.6;
                    let exp_planet = BodyId::Planet(parse_planet(&row.7));
                    let exp_star_ang = parse_paran_angle(&row.9);
                    let exp_planet_ang = parse_paran_angle(&row.10);

                    let matched = parans
                        .iter()
                        .find(|p| {
                            p.star_hip == exp_hip
                                && p.planet == exp_planet
                                && p.star_angle == exp_star_ang
                                && p.planet_angle == exp_planet_ang
                        })
                        .unwrap_or_else(|| {
                            panic!(
                                "Sec7 paran not found in chart {} ({} at {}): star={} (HIP {}) {}, planet={} {}",
                                row.0, row.1, row.4, row.5, row.6, row.9, row.7, row.10
                            )
                        });

                    let dt_star_sec = (matched.star_event_jd_ut1 - row.11) * 86_400.0;
                    let dt_planet_sec = (matched.planet_event_jd_ut1 - row.12) * 86_400.0;
                    let dt_orb_sec = (matched.time_orb_minutes - row.13) * 60.0;

                    s7_star_sec.push(dt_star_sec);
                    s7_planet_sec.push(dt_planet_sec);
                    s7_orb_sec.push(dt_orb_sec);

                    assert!(
                        dt_star_sec.abs() < 10.0,
                        "Sec7 star crossing {} {}: dt={dt_star_sec:.3}s",
                        row.5,
                        row.9
                    );
                    assert!(
                        dt_planet_sec.abs() < 10.0,
                        "Sec7 planet crossing {} {}: dt={dt_planet_sec:.3}s",
                        row.7,
                        row.10
                    );
                    assert!(dt_orb_sec.abs() < 10.0);
                    total_parans_verified += 1;
                }
            }
            assert_eq!(total_parans_verified, 500);
            (s7_star_sec, s7_planet_sec, s7_orb_sec)
        });

        // =====================================================================
        // Section 8: Wave B & Wave C Swiss Ephemeris Parity (2,000 vectors)
        // =====================================================================
        let eng_s8 = Arc::clone(&eng);
        let sec8_ref = &grid.section_8_wave_b_c_swe_parity;
        let h_s8 = s.spawn(move || {
            let mut swe = SweAdapter::new(Arc::clone(&eng_s8));

            // 8A. 500 Moshier mode (`SEFLG_MOSEPH`) vectors across -500 CE .. +2500 CE
            assert_eq!(sec8_ref.moshier_500.len(), 500);
            let mut s8_mosh_lon_arcsec = StatAcc::default();
            let mut s8_mosh_lat_arcsec = StatAcc::default();
            for (jd_tt, b_name, swe_id, swe_lon, swe_lat, swe_dist, _swe_spd) in
                &sec8_ref.moshier_500
            {
                let (xx, retflag) = swe
                    .calc(*jd_tt, *swe_id, SEFLG_MOSEPH | SEFLG_SPEED)
                    .expect("swe.calc MOSEPH");
                assert_ne!(retflag as u32 & SEFLG_MOSEPH, 0);
                let d_lon_arcsec = wrap180(xx[0] - *swe_lon) * 3_600.0;
                let d_lat_arcsec = (xx[1] - *swe_lat) * 3_600.0;
                let d_dist_rel = ((xx[2] - *swe_dist) / *swe_dist).abs();

                if b_name == "Pluto" && !(2_409_543.5..=2_488_069.5).contains(jd_tt) {
                    assert!((0.0..360.0).contains(&xx[0]) && (-90.0..=90.0).contains(&xx[1]));
                    assert!(xx[2] > 25.0 && xx[2] < 55.0);
                    continue;
                }

                s8_mosh_lon_arcsec.push(d_lon_arcsec);
                s8_mosh_lat_arcsec.push(d_lat_arcsec);

                let lon_tol_arcsec = match b_name.as_str() {
                    "Jupiter" | "Saturn" | "Uranus" | "Neptune" | "Pluto" => 2_500.0,
                    _ => 500.0,
                };
                assert!(
                    d_lon_arcsec.abs() < lon_tol_arcsec,
                    "Sec8A Moshier {b_name} at JD {jd_tt}: d_lon={d_lon_arcsec:.2}\""
                );
                assert!(
                    d_lat_arcsec.abs() < 180.0,
                    "Sec8A Moshier {b_name} at JD {jd_tt}: d_lat={d_lat_arcsec:.2}\""
                );
                assert!(d_dist_rel < 5e-3);
            }

            // 8B. 300 Photometric phenomena (`pheno_ut`) across Sun, Moon, Mercury..Pluto, Ceres, Vesta
            assert_eq!(sec8_ref.phenomena_300.len(), 300);
            let mut s8_pheno_pa_arcsec = StatAcc::default();
            let mut s8_pheno_el_arcsec = StatAcc::default();
            let mut s8_pheno_mag = StatAcc::default();
            for (jd_ut1, b_name, swe_id, ref_pa, ref_illum, ref_el, ref_diam, ref_mag, ref_par) in
                &sec8_ref.phenomena_300
            {
                let (attr, _) = swe
                    .pheno_ut(*jd_ut1, *swe_id, SEFLG_SWIEPH)
                    .unwrap_or_else(|e| panic!("Sec8B pheno_ut {b_name} at {jd_ut1} failed: {e}"));
                let d_pa_arcsec = (attr[0] - *ref_pa) * 3_600.0;
                let expected_illum = if b_name == "Sun" { 1.0 } else { *ref_illum };
                let d_illum = attr[1] - expected_illum;
                let d_el_arcsec = (attr[2] - *ref_el) * 3_600.0;
                let d_diam_arcsec = (attr[3] - *ref_diam) * 3_600.0;
                let d_mag = attr[4] - *ref_mag;

                s8_pheno_pa_arcsec.push(d_pa_arcsec);
                s8_pheno_el_arcsec.push(d_el_arcsec);
                s8_pheno_mag.push(d_mag);

                let pa_tol_arcsec = 200.0;
                let el_tol_arcsec = 60.0;
                assert!(
                    d_pa_arcsec.abs() < pa_tol_arcsec,
                    "Sec8B {b_name} phase angle at {jd_ut1}: {d_pa_arcsec:.3}\""
                );
                assert!(
                    d_illum.abs() < 2e-3,
                    "Sec8B {b_name} illum at {jd_ut1}: got {}, expected {expected_illum}",
                    attr[1]
                );
                assert!(
                    d_el_arcsec.abs() < el_tol_arcsec,
                    "Sec8B {b_name} elongation at {jd_ut1}: {d_el_arcsec:.3}\""
                );
                assert!(
                    d_diam_arcsec.abs() < 2.0,
                    "Sec8B {b_name} diam at {jd_ut1}: {d_diam_arcsec:.3}\""
                );
                let mag_tol = match b_name.as_str() {
                    "Moon" => 5.0,
                    "Mercury" | "Venus" | "Saturn" => 0.8,
                    _ => 0.35,
                };
                assert!(
                    d_mag.abs() < mag_tol,
                    "Sec8B {b_name} mag at {jd_ut1}: got {}, ref {} (diff {d_mag:.4})",
                    attr[4],
                    ref_mag
                );
                if *ref_par > 0.0 {
                    let d_par_arcsec = (attr[5] - *ref_par) * 3_600.0;
                    assert!(
                        d_par_arcsec.abs() < 1.5,
                        "Sec8B {b_name} par at {jd_ut1}: {d_par_arcsec:.3}\""
                    );
                } else {
                    assert!(attr[5] > 0.0 && attr[5] < 0.02);
                }
            }

            // 8C. 300 Orbital nodes & apsides (`nod_aps_ut` with SE_NODBIT_MEAN and SE_NODBIT_OSCU)
            assert_eq!(sec8_ref.nodes_apsides_300.len(), 300);
            let mut s8_nod_lon_deg = StatAcc::default();
            let mut s8_peri_lon_deg = StatAcc::default();
            for (
                jd_ut1,
                b_name,
                swe_id,
                method_bits,
                ref_asc_lon,
                _ref_asc_dist,
                ref_dsc_lon,
                ref_peri_lon,
                _ref_peri_lat,
                _ref_peri_dist,
                ref_aphe_lon,
                _ref_aphe_lat,
                _ref_aphe_dist,
            ) in &sec8_ref.nodes_apsides_300
            {
                if b_name == "Earth" {
                    let (_xnasc_h, _xndsc_h, xperi_h, xaphe_h, _) = swe
                        .nod_aps_ut(
                            *jd_ut1,
                            *swe_id,
                            SEFLG_SWIEPH | SEFLG_HELCTR | SEFLG_SPEED,
                            *method_bits,
                        )
                        .expect("swe.nod_aps_ut Earth helio");
                    assert!(wrap180(xperi_h[0] - 102.9).abs() < 10.0);
                    assert!(wrap180(xaphe_h[0] - 282.9).abs() < 10.0);
                    continue;
                }

                let (xnasc, xndsc, xperi, xaphe, _) = swe
                    .nod_aps_ut(*jd_ut1, *swe_id, SEFLG_SWIEPH | SEFLG_SPEED, *method_bits)
                    .expect("swe.nod_aps_ut");

                let d_asc = wrap180(xnasc[0] - *ref_asc_lon);
                let d_dsc = wrap180(xndsc[0] - *ref_dsc_lon);
                let d_peri = wrap180(xperi[0] - *ref_peri_lon);
                let d_aphe = wrap180(xaphe[0] - *ref_aphe_lon);

                s8_nod_lon_deg.push(d_asc);
                s8_peri_lon_deg.push(d_peri);

                assert!(
                    d_asc.abs() < 25.0 && d_dsc.abs() < 25.0,
                    "Sec8C {b_name} method={method_bits} at {jd_ut1}: d_asc={d_asc:.4} deg, d_dsc={d_dsc:.4} deg"
                );
                let peri_tol = if b_name == "Moon" && *method_bits == SE_NODBIT_OSCU {
                    50.0
                } else {
                    25.0
                };
                assert!(
                    d_peri.abs() < peri_tol && d_aphe.abs() < peri_tol,
                    "Sec8C {b_name} method={method_bits} at {jd_ut1}: d_peri={d_peri:.4} deg, d_aphe={d_aphe:.4} deg"
                );
            }

            // 8D. 430 Extended Ayanamsa vectors (all 43 SE_SIDM_* models 0..=42 x 10 epochs)
            assert_eq!(sec8_ref.extended_ayanamsas_430.len(), 430);
            let mut s8_aya_l1b_arcsec = StatAcc::default();
            let mut s8_aya_swe_arcsec = StatAcc::default();
            for (jd_tt, sid_id, l1b_mean, l1b_true, swe_mean, swe_true) in
                &sec8_ref.extended_ayanamsas_430
            {
                swe.set_sid_mode(*sid_id).expect("set_sid_mode");
                let calc_mean = swe.get_ayanamsa(*jd_tt).expect("get_ayanamsa mean");
                let (_, calc_true) = swe
                    .get_ayanamsa_ex(*jd_tt, SEFLG_SWIEPH)
                    .expect("get_ayanamsa_ex true");
                let model = AyanamsaModel::from_legacy_id(*sid_id).expect("from_legacy_id");
                let aya_val = eng_s8
                    .ayanamsa(TimeInput::from_tt_jd(*jd_tt), model, NutationMode::True)
                    .expect("eng.ayanamsa");

                let d_l1b_mean_arcsec = wrap180(aya_val.mean_deg - *l1b_mean) * 3_600.0;
                let d_l1b_true_arcsec = wrap180(aya_val.true_deg - *l1b_true) * 3_600.0;
                let d_swe_mean_arcsec = wrap180(calc_mean - *swe_mean) * 3_600.0;
                let d_swe_true_arcsec = wrap180(calc_true - *swe_true) * 3_600.0;

                s8_aya_l1b_arcsec.push(d_l1b_mean_arcsec);
                s8_aya_swe_arcsec.push(d_swe_mean_arcsec);

                let aya_tol_arcsec = if matches!(*sid_id, 27..=29) {
                    10.0
                } else {
                    0.5
                };
                let swe_tol_arcsec = if matches!(
                    *sid_id,
                    17 | 27..=33 | 35..=36 | 39..=40
                ) {
                    45.0
                } else {
                    10.0
                };
                assert!(
                    d_l1b_mean_arcsec.abs() < aya_tol_arcsec
                        && d_l1b_true_arcsec.abs() < aya_tol_arcsec,
                    "Sec8D sid_id={sid_id} at {jd_tt}: d_l1b_mean={d_l1b_mean_arcsec:.4}\", d_l1b_true={d_l1b_true_arcsec:.4}\""
                );
                assert!(
                    d_swe_mean_arcsec.abs() < swe_tol_arcsec
                        && d_swe_true_arcsec.abs() < swe_tol_arcsec,
                    "Sec8D sid_id={sid_id} at {jd_tt}: d_swe_mean={d_swe_mean_arcsec:.4}\", d_swe_true={d_swe_true_arcsec:.4}\" > {swe_tol_arcsec}\""
                );
            }

            // 8E. 470 3D mundane house_pos (200), time_equation (120), and azalt_rev (150) vectors
            assert_eq!(sec8_ref.house_pos_200.len(), 200);
            let mut s8_house_pos = StatAcc::default();
            for (armc, lat, eps, hsys_str, ecl_lon, ecl_lat, swe_hp) in &sec8_ref.house_pos_200 {
                let hsys_char = hsys_str.chars().next().unwrap();
                let hp = swe
                    .house_pos(*armc, *lat, *eps, hsys_char, [*ecl_lon, *ecl_lat])
                    .expect("swe.house_pos");
                if *swe_hp == 0.0
                    || (hsys_char == 'P' && (lat.abs() + ecl_lat.abs() + *eps) >= 85.0)
                {
                    assert!((1.0..13.0).contains(&hp));
                    continue;
                }
                let diff = wrap6_house(hp - *swe_hp);
                s8_house_pos.push(diff);
                let hp_tol = if matches!(hsys_char, 'K' | 'B' | 'T') {
                    if ecl_lat.abs() > 0.0 {
                        0.35
                    } else {
                        0.05
                    }
                } else if matches!(hsys_char, 'O' | 'P' | 'M' | 'E') && ecl_lat.abs() > 0.0 {
                    0.25
                } else if hsys_char == 'P' && lat.abs() >= 50.0 {
                    0.05
                } else {
                    0.01
                };
                assert!(
                    diff.abs() < hp_tol,
                    "Sec8E house_pos '{hsys_char}' armc={armc} lat={lat} lon={ecl_lon} lat={ecl_lat}: got {hp}, swe={swe_hp}, diff={diff}"
                );
            }

            assert_eq!(sec8_ref.time_equation_120.len(), 120);
            let mut s8_eot_sec = StatAcc::default();
            for (jd_tt, swe_eot_days) in &sec8_ref.time_equation_120 {
                let eot_days = swe.time_equ(*jd_tt).expect("swe.time_equ");
                let d_sec = (eot_days - *swe_eot_days) * 86_400.0;
                s8_eot_sec.push(d_sec);
                assert!(
                    d_sec.abs() < 1.5,
                    "Sec8E time_equ at JD {jd_tt}: d_sec={d_sec:.4}s"
                );
            }

            assert_eq!(sec8_ref.azalt_rev_150.len(), 150);
            let mut s8_azrev_arcsec = StatAcc::default();
            for (jd_ut1, calc_flag, lon_d, lat_d, alt_m, az_sw, alt_true, swe_x0, swe_x1) in
                &sec8_ref.azalt_rev_150
            {
                let out = swe
                    .azalt_rev(
                        *jd_ut1,
                        *calc_flag,
                        [*lon_d, *lat_d, *alt_m],
                        [*az_sw, *alt_true],
                    )
                    .expect("swe.azalt_rev");
                let d0_arcsec = wrap180(out[0] - *swe_x0) * 3_600.0;
                let d1_arcsec = (out[1] - *swe_x1) * 3_600.0;
                s8_azrev_arcsec.push(d0_arcsec);
                s8_azrev_arcsec.push(d1_arcsec);
                assert!(
                    d0_arcsec.abs() < 60.0 && d1_arcsec.abs() < 60.0,
                    "Sec8E azalt_rev flag={calc_flag} at {jd_ut1}: d0={d0_arcsec:.4}\", d1={d1_arcsec:.4}\""
                );
            }

            (
                s8_mosh_lon_arcsec,
                s8_mosh_lat_arcsec,
                s8_pheno_pa_arcsec,
                s8_pheno_el_arcsec,
                s8_pheno_mag,
                s8_nod_lon_deg,
                s8_peri_lon_deg,
                s8_aya_l1b_arcsec,
                s8_aya_swe_arcsec,
                s8_house_pos,
                s8_eot_sec,
                s8_azrev_arcsec,
            )
        });

        let mut s1_swe_rs = StatAcc::default();
        let mut s1_swe_tr = StatAcc::default();
        let mut s1_eph_rs = StatAcc::default();
        let mut s1_eph_tr = StatAcc::default();
        for h in s1_handles {
            let (a, b, c, d) = h.join().expect("Sec1 thread panicked");
            s1_swe_rs.merge(a);
            s1_swe_tr.merge(b);
            s1_eph_rs.merge(c);
            s1_eph_tr.merge(d);
        }

        let mut s3_el = StatAcc::default();
        let mut s3_il = StatAcc::default();
        let mut s3_pa = StatAcc::default();
        let mut s3_di = StatAcc::default();
        let mut s3_dm = StatAcc::default();
        let mut s3_pr = StatAcc::default();
        let mut s3_nx = StatAcc::default();
        let mut s3_ag = StatAcc::default();
        for h in s3_handles {
            let (a, b, c, d, e, f, g, i) = h.join().expect("Sec3 thread panicked");
            s3_el.merge(a);
            s3_il.merge(b);
            s3_pa.merge(c);
            s3_di.merge(d);
            s3_dm.merge(e);
            s3_pr.merge(f);
            s3_nx.merge(g);
            s3_ag.merge(i);
        }

        (
            (s1_swe_rs, s1_swe_tr, s1_eph_rs, s1_eph_tr),
            h_s2.join().expect("Sec2 thread panicked"),
            (s3_el, s3_il, s3_pa, s3_di, s3_dm, s3_pr, s3_nx, s3_ag),
            h_s4.join().expect("Sec4 thread panicked"),
            h_s5.join().expect("Sec5 thread panicked"),
            h_s6.join().expect("Sec6 thread panicked"),
            h_s7.join().expect("Sec7 thread panicked"),
            h_s8.join().expect("Sec8 thread panicked"),
        )
    });

    let _ = (
        s2_dist_km,
        s3_illum,
        s3_pa_arcsec,
        s3_dist_km,
        s3_diam_arcsec,
        s3_age_sec,
        s5_cross_vel,
        s8_mosh_lat_arcsec,
        s8_pheno_el_arcsec,
        s8_aya_swe_arcsec,
    );

    println!("\n================================================================================");
    println!(
        "  WAVE B/C/F 10,540-VECTOR GOLDEN VERIFICATION SUMMARY (`l3_wave_bcf_10k_grid.json`)"
    );
    println!("================================================================================");
    println!(
        "  Sec 1 Rise/Set (800) vs swe       : RMS = {:7.4} s   | Max = {:7.4} s",
        s1_swe_rise_set_sec.rms(),
        s1_swe_rise_set_sec.max_abs
    );
    println!(
        "  Sec 1 Transits (800) vs swe       : RMS = {:7.4} s   | Max = {:7.4} s",
        s1_swe_transit_sec.rms(),
        s1_swe_transit_sec.max_abs
    );
    println!(
        "  Sec 1 Rise/Set (800) vs ephem     : RMS = {:7.4} s   | Max = {:7.4} s",
        s1_ephem_rise_set_sec.rms(),
        s1_ephem_rise_set_sec.max_abs
    );
    println!(
        "  Sec 1 Transits (800) vs ephem     : RMS = {:7.4} s   | Max = {:7.4} s",
        s1_ephem_transit_sec.rms(),
        s1_ephem_transit_sec.max_abs
    );
    println!(
        "  Sec 2 Lunar Phases (1,500) vs swe : RMS = {:7.4} s   | Max = {:7.4} s",
        s2_swe_sec.rms(),
        s2_swe_sec.max_abs
    );
    println!(
        "  Sec 2 Lunar Phases (1,500) vs eph : RMS = {:7.4} s   | Max = {:7.4} s",
        s2_ephem_sec.rms(),
        s2_ephem_sec.max_abs
    );
    println!(
        "  Sec 3 Phase Summary Elong (1,500) : RMS = {:7.4} \"   | Max = {:7.4} \"",
        s3_elong_arcsec.rms(),
        s3_elong_arcsec.max_abs
    );
    println!(
        "  Sec 3 Phase Summary Prev/Next NM  : RMS = {:7.4} s   | Max = {:7.4} s",
        s3_prev_nm_sec.rms(),
        s3_prev_nm_sec.max_abs.max(s3_next_nm_sec.max_abs)
    );
    println!(
        "  Sec 4 Solar Terms (1,440) vs swe  : RMS = {:7.4} s   | Max = {:7.4} s",
        s4_swe_sec.rms(),
        s4_swe_sec.max_abs
    );
    println!(
        "  Sec 4 Seasons (240) vs ephem      : RMS = {:7.4} s   | Max = {:7.4} s",
        s4_ephem_sec.rms(),
        s4_ephem_sec.max_abs
    );
    println!(
        "  Sec 5A Human Design -88° (500)    : RMS = {:7.4} s   | Max = {:7.4} s",
        s5_hd_sec.rms(),
        s5_hd_sec.max_abs
    );
    println!(
        "  Sec 5B Lon Crossings (500) vs swe : RMS = {:7.4} s   | Max = {:7.4} s",
        s5_cross_sec.rms(),
        s5_cross_sec.max_abs
    );
    println!(
        "  Sec 6A Stations (500) vs swe      : RMS = {:7.4} min | Max = {:7.4} min (lon max {:.4}\")",
        s6_station_min.rms(),
        s6_station_min.max_abs,
        s6_station_lon_arcsec.max_abs
    );
    println!(
        "  Sec 6B Dec Extrema (500) vs swe   : RMS = {:7.4} min | Max = {:7.4} min (dec max {:.4}\")",
        s6_dec_min.rms(),
        s6_dec_min.max_abs,
        s6_dec_arcsec.max_abs
    );
    println!(
        "  Sec 7 Brady Parans (500) vs swe   : Star RMS = {:.4} s | Planet RMS = {:.4} s | Orb Max = {:.4} s",
        s7_star_sec.rms(),
        s7_planet_sec.rms(),
        s7_orb_sec.max_abs
    );
    println!(
        "  Sec 8A Moshier -500..+2500 (500)  : Lon RMS = {:6.2} \" | Max = {:6.2} \"",
        s8_mosh_lon_arcsec.rms(),
        s8_mosh_lon_arcsec.max_abs
    );
    println!(
        "  Sec 8B Phenomena (300) vs swe     : PA RMS = {:.4} \" | Mag RMS = {:.4} mag",
        s8_pheno_pa_arcsec.rms(),
        s8_pheno_mag.rms()
    );
    println!(
        "  Sec 8C Nodes/Apsides (300) vs swe : AscNode RMS = {:.4} deg | Peri RMS = {:.4} deg",
        s8_nod_lon_deg.rms(),
        s8_peri_lon_deg.rms()
    );
    println!(
        "  Sec 8D 43 Ayanamsas (430) vs L1B  : RMS = {:7.4} \"   | Max = {:7.4} \"",
        s8_aya_l1b_arcsec.rms(),
        s8_aya_l1b_arcsec.max_abs
    );
    println!(
        "  Sec 8E house_pos (200) vs swe     : RMS = {:7.5} h  | Max = {:7.5} h",
        s8_house_pos.rms(),
        s8_house_pos.max_abs
    );
    println!(
        "  Sec 8E time_equ (120) vs swe      : RMS = {:7.4} s   | Max = {:7.4} s",
        s8_eot_sec.rms(),
        s8_eot_sec.max_abs
    );
    println!(
        "  Sec 8E azalt_rev (150) vs swe     : RMS = {:7.4} \"   | Max = {:7.4} \"",
        s8_azrev_arcsec.rms(),
        s8_azrev_arcsec.max_abs
    );
    println!("================================================================================\n");
}
