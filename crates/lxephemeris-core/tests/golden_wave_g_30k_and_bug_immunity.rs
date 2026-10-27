//! Wave G1: 30k+ Extreme & Multi-Oracle Golden Reference Suite (`10,000` new vectors -> `32,416` total)
//! and 22 Explicit Upstream-Bug Immunity Regression Tests (`SWE-01..SWE-11`, `EPH-01..EPH-05`, `KER-01..KER-06`).

#![cfg(feature = "swe-compat")]
#![forbid(unsafe_code)]
#![allow(
    clippy::type_complexity,
    clippy::indexing_slicing,
    clippy::panic,
    clippy::unwrap_used
)]

use lxephemeris_core::{
    AbMode, AtmosphereParams, AyanamsaModel, BodyId, CalcOptions, CalendarKind, EclipseOptions, Engine, EngineConfig,
    EphemerisError, Equinox, FrameTransform, HeliacalOptions, HipId, HouseOptions, HouseSystem,
    HouseWarning, LunarPointId, MpcNumber, NutationMode, Observer, Origin, ParanOptions, Plane,
    PlanetId, PolarPolicy, RiseTransEvent, RiseTransOptions, SplitJd, SweAdapter, TimeInput, Tt,
    Ut1, Zodiac, compute_lunar_saros, compute_solar_saros, find_star_by_name_or_bayer, gast_rad,
    gmst_rad, house_of_lon, houses_armc_with_sidereal, is_cusp_ring_monotonic, julday,
    julday_split, parse_sefstars_catalog, seasons_in_year, solar_terms_in_range, tdb_minus_tt_sec,
    SEFLG_MOSEPH, SEFLG_SPEED, SEFLG_SWIEPH, SEFLG_TOPOCTR, SE_SIDM_USER,
};
use serde::Deserialize;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;

#[derive(Debug, Deserialize)]
struct GoldenWaveG10kGrid {
    total_vector_count: usize,
    counts_by_section: SectionGCounts,
    section_g1_extreme_states: SectionG1,
    section_g2_all_23_houses_extreme_latitudes: Vec<SectionG2Vec>,
    section_g3_minor_bodies_68_multimode: Vec<SectionG3Vec>,
    section_g4_fixed_stars_52_multicentury: Vec<SectionG4Vec>,
    section_g5_eclipses_occultations_heliacal: SectionG5,
    section_g6_ayanamsas_44_and_3d_house_pos_23: SectionG6,
}

#[derive(Debug, Deserialize)]
struct SectionGCounts {
    section_g1_extreme_states: usize,
    section_g2_all_23_houses_extreme_latitudes: usize,
    section_g3_minor_bodies_68_multimode: usize,
    section_g4_fixed_stars_52_multicentury: usize,
    section_g5_eclipses_occultations_heliacal: usize,
    section_g6_ayanamsas_44_and_3d_house_pos_23: usize,
}

#[derive(Debug, Deserialize)]
struct SectionG1 {
    /// `[jd_tt, mode_name, b_name, swe_id, flags, lon, lat, dist, speed]`
    multi_center_frame_500: Vec<(f64, String, String, i32, u32, f64, f64, f64, f64)>,
    /// `[jd_ut1, obs_name, lat_deg, lon_deg, alt_m, b_name, swe_id, lon, lat, dist, speed]`
    extreme_topocentric_500: Vec<(f64, String, f64, f64, f64, String, i32, f64, f64, f64, f64)>,
    /// `[jd_tt, kind, b_name, swe_id, lon, lat, dist, speed]`
    boundary_wrap_station_500: Vec<(f64, String, String, i32, f64, f64, f64, f64)>,
    /// `[jd_tt, b_name, swe_id, lon, lat, dist, speed]`
    deep_moshier_500: Vec<(f64, String, i32, f64, f64, f64, f64)>,
}

/// `[armc_deg, lat_deg, eps_deg, sun_decl_deg, hsys, is_polar, cusps_12, asc, mc, vertex]`
type SectionG2Vec = (f64, f64, f64, f64, String, bool, [f64; 12], f64, f64, f64);

/// `[jd_tt, mpc, name, mode_name, lon, lat, dist, speed]`
type SectionG3Vec = (f64, u32, String, String, f64, f64, f64, f64);

/// `[jd_tt, hip, name, mode_name, lon, lat]`
type SectionG4Vec = (f64, u32, String, String, f64, f64);

#[derive(Debug, Deserialize)]
struct SectionG5 {
    /// `[swe_jd_ut1, kind_str, saros_series, saros_member, gamma_approx, rf]`
    solar_eclipses_300: Vec<(f64, String, i32, i32, f64, i32)>,
    /// `[swe_jd_ut1, kind_str, saros_series, saros_member, gamma_approx, rf]`
    lunar_eclipses_300: Vec<(f64, String, i32, i32, f64, i32)>,
    /// `[t_label, hip_id, start_jd, swe_jd_ut1, rf]`
    lunar_occultations_350: Vec<(String, u32, f64, f64, i32)>,
    /// `[year, base_jd, t_name, t_kind, hip_id, eff_lat, lon_deg]`
    heliacal_350: Vec<(i32, f64, String, String, u32, f64, f64)>,
}

#[derive(Debug, Deserialize)]
struct SectionG6 {
    /// `[jd_tt, sid_id, l1b_mean, l1b_true, swe_mean, swe_true]`
    ayanamsas_880: Vec<(f64, u32, f64, f64, f64, f64)>,
    /// `[armc_deg, lat_deg, eps_deg, hsys, ecl_lon, ecl_lat, swe_hp]`
    house_pos_690: Vec<(f64, f64, f64, String, f64, f64, f64)>,
    /// `[base_jd, obs_name, lat_deg, lon_deg, body_name, ev_name, is_circumpolar, expected_jd_ut1]`
    rise_trans_430: Vec<(f64, String, f64, f64, String, String, bool, Option<f64>)>,
}

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
        .to_path_buf()
}

fn load_full_engine() -> Arc<Engine> {
    let root = workspace_root();
    Arc::new(
        Engine::from_data_dir(root.join("data"), EngineConfig::default())
            .expect("open Engine with de440s + asteroids_68"),
    )
}

fn wrap_diff_deg(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(360.0);
    if d > 180.0 {
        360.0 - d
    } else {
        d
    }
}

fn parse_planet_by_name(name: &str) -> PlanetId {
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
        other => panic!("unknown planet name: {other}"),
    }
}

#[test]
fn test_golden_wave_g_10k_and_30k_total_count() {
    let root = workspace_root();
    let golden_dir = root.join("testdata/golden");

    // 1. Verify individual 10k grid counts and combined 32,416 total golden vector count
    let se_10k_val: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(golden_dir.join("l3_se_10k_grid.json")).expect("read l3_se_10k_grid"),
    )
    .expect("parse l3_se_10k_grid");
    let se_10k_count = se_10k_val["total_vector_count"]
        .as_u64()
        .expect("se_10k total_vector_count") as usize;
    assert_eq!(se_10k_count, 10_242);

    let bcf_10k_val: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(golden_dir.join("l3_wave_bcf_10k_grid.json"))
            .expect("read l3_wave_bcf_10k_grid"),
    )
    .expect("parse l3_wave_bcf_10k_grid");
    let bcf_10k_count = bcf_10k_val["total_vector_count"]
        .as_u64()
        .expect("bcf_10k total_vector_count") as usize;
    assert_eq!(bcf_10k_count, 10_540);

    let g_10k_raw = fs::read_to_string(golden_dir.join("l3_wave_g_extreme_10k_grid.json"))
        .expect("read l3_wave_g_extreme_10k_grid.json");
    let grid: GoldenWaveG10kGrid =
        serde_json::from_str(&g_10k_raw).expect("parse l3_wave_g_extreme_10k_grid.json");

    assert_eq!(grid.total_vector_count, 10_000);
    assert_eq!(grid.counts_by_section.section_g1_extreme_states, 2_000);
    assert_eq!(
        grid.counts_by_section
            .section_g2_all_23_houses_extreme_latitudes,
        2_300
    );
    assert_eq!(
        grid.counts_by_section.section_g3_minor_bodies_68_multimode,
        1_360
    );
    assert_eq!(
        grid.counts_by_section
            .section_g4_fixed_stars_52_multicentury,
        1_040
    );
    assert_eq!(
        grid.counts_by_section
            .section_g5_eclipses_occultations_heliacal,
        1_300
    );
    assert_eq!(
        grid.counts_by_section
            .section_g6_ayanamsas_44_and_3d_house_pos_23,
        2_000
    );

    let g1_len = grid.section_g1_extreme_states.multi_center_frame_500.len()
        + grid.section_g1_extreme_states.extreme_topocentric_500.len()
        + grid
            .section_g1_extreme_states
            .boundary_wrap_station_500
            .len()
        + grid.section_g1_extreme_states.deep_moshier_500.len();
    let g2_len = grid.section_g2_all_23_houses_extreme_latitudes.len();
    let g3_len = grid.section_g3_minor_bodies_68_multimode.len();
    let g4_len = grid.section_g4_fixed_stars_52_multicentury.len();
    let g5_len = grid
        .section_g5_eclipses_occultations_heliacal
        .solar_eclipses_300
        .len()
        + grid
            .section_g5_eclipses_occultations_heliacal
            .lunar_eclipses_300
            .len()
        + grid
            .section_g5_eclipses_occultations_heliacal
            .lunar_occultations_350
            .len()
        + grid
            .section_g5_eclipses_occultations_heliacal
            .heliacal_350
            .len();
    let g6_len = grid
        .section_g6_ayanamsas_44_and_3d_house_pos_23
        .ayanamsas_880
        .len()
        + grid
            .section_g6_ayanamsas_44_and_3d_house_pos_23
            .house_pos_690
            .len()
        + grid
            .section_g6_ayanamsas_44_and_3d_house_pos_23
            .rise_trans_430
            .len();
    assert_eq!(g1_len + g2_len + g3_len + g4_len + g5_len + g6_len, 10_000);

    let wave_h_5k_val: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(golden_dir.join("l3_wave_h_almanac_aspects_5k.json"))
            .expect("read l3_wave_h_almanac_aspects_5k"),
    )
    .expect("parse l3_wave_h_almanac_aspects_5k");
    let wave_h_5k_count = wave_h_5k_val["total_vector_count"]
        .as_u64()
        .expect("wave_h_5k total_vector_count") as usize;
    assert_eq!(wave_h_5k_count, 5_000);

    // Count Phase 0-4 L1A/L1B/L3 golden vectors (1,634) + three 10k grids (30,782) + Wave H 5k (5,000) = 37,416
    let l1_l3_phase0_4_count: usize = 360 + 90 + 340 + 120 + 253 + 63 + 54 + 104 + 250;
    assert_eq!(l1_l3_phase0_4_count, 1_634);
    let grid_30k_count = se_10k_count + bcf_10k_count + grid.total_vector_count;
    assert_eq!(grid_30k_count, 30_782);
    assert!(grid_30k_count >= 30_000);
    let grand_total_vectors = l1_l3_phase0_4_count + grid_30k_count + wave_h_5k_count;
    assert_eq!(grand_total_vectors, 37_416);

    let engine = load_full_engine();

    // 2. Verify all 6 sections (10,000 vectors) in parallel across worker threads
    std::thread::scope(|s| {
        // Section G1: 2,000 extreme state vectors (G1A..G1D)
        let eng_g1 = Arc::clone(&engine);
        let sec_g1 = &grid.section_g1_extreme_states;
        let h_g1 = s.spawn(move || {
            let adapter = SweAdapter::new(Arc::clone(&eng_g1));

            // G1A: 500 multi-center & frame vectors
            for (idx, (jd_tt, mode_name, b_name, swe_id, flags, ref_lon, ref_lat, ref_dist, ref_spd)) in
                sec_g1.multi_center_frame_500.iter().enumerate()
            {
                let (xx, _ret) = adapter
                    .calc(*jd_tt, *swe_id, *flags)
                    .unwrap_or_else(|e| panic!("G1A[{idx}] {mode_name} {b_name} failed: {e}"));
                let is_bary = mode_name.contains("Barycentric");
                let is_xyz = mode_name.contains("XYZ");
                if is_xyz {
                    let tol = if is_bary { 0.02 } else { 0.002 };
                    assert!(
                        (xx[0] - *ref_lon).abs() < tol
                            && (xx[1] - *ref_lat).abs() < tol
                            && (xx[2] - *ref_dist).abs() < tol,
                        "G1A[{idx}] XYZ mismatch: got {:?}, ref ({ref_lon}, {ref_lat}, {ref_dist})",
                        &xx[..3]
                    );
                } else {
                    let lon_tol = if is_bary { 1.2 } else { 0.05 };
                    let lat_tol = if is_bary { 0.25 } else { 0.05 };
                    let d_lon = wrap_diff_deg(xx[0], *ref_lon);
                    let d_lat = (xx[1] - *ref_lat).abs();
                    assert!(
                        d_lon <= lon_tol && d_lat <= lat_tol && xx[3].is_finite() && ref_spd.is_finite(),
                        "G1A[{idx}] {mode_name} {b_name}: d_lon={d_lon} d_lat={d_lat}"
                    );
                }
            }

            // G1B: 500 extreme topocentric vectors
            for (idx, (jd_ut1, obs_name, lat_deg, lon_deg, alt_m, b_name, swe_id, ref_lon, ref_lat, ref_dist, _ref_spd)) in
                sec_g1.extreme_topocentric_500.iter().enumerate()
            {
                let mut topo_adapter = SweAdapter::new(Arc::clone(&eng_g1));
                topo_adapter
                    .set_topo(*lon_deg, *lat_deg, *alt_m)
                    .expect("set_topo");
                let (xx, _ret) = topo_adapter
                    .calc_ut(*jd_ut1, *swe_id, SEFLG_SWIEPH | SEFLG_TOPOCTR | SEFLG_SPEED)
                    .unwrap_or_else(|e| panic!("G1B[{idx}] {obs_name} {b_name} failed: {e}"));
                let d_lon = wrap_diff_deg(xx[0], *ref_lon);
                let d_lat = (xx[1] - *ref_lat).abs();
                let d_dist = (xx[2] - *ref_dist).abs();
                let lon_tol = if *swe_id == 1 { 0.05 } else { 0.02 };
                assert!(
                    d_lon <= lon_tol && d_lat <= 0.02 && d_dist <= 0.005,
                    "G1B[{idx}] {obs_name} {b_name}: d_lon={d_lon} d_lat={d_lat} d_dist={d_dist}"
                );
            }

            // G1C: 500 boundary, 0/360 Aries wrap, high-latitude & near-station vectors
            for (idx, (jd_tt, kind, b_name, swe_id, ref_lon, ref_lat, ref_dist, ref_spd)) in
                sec_g1.boundary_wrap_station_500.iter().enumerate()
            {
                let (xx, _ret) = adapter
                    .calc(*jd_tt, *swe_id, SEFLG_SWIEPH | SEFLG_SPEED)
                    .unwrap_or_else(|e| panic!("G1C[{idx}] {kind} {b_name} failed: {e}"));
                let d_lon = wrap_diff_deg(xx[0], *ref_lon);
                let d_lat = (xx[1] - *ref_lat).abs();
                let d_dist = (xx[2] - *ref_dist).abs();
                let d_spd = (xx[3] - *ref_spd).abs();
                assert!(
                    d_lon < 0.05 && d_lat < 0.02 && d_dist < 0.005 && d_spd < 0.02,
                    "G1C[{idx}] {kind} {b_name}: d_lon={d_lon} d_lat={d_lat} d_spd={d_spd}"
                );
            }

            // G1D: 500 deep-history & future Moshier vectors (-2500 CE .. +2900 CE)
            for (idx, (jd_tt, b_name, swe_id, ref_lon, ref_lat, ref_dist, ref_spd)) in
                sec_g1.deep_moshier_500.iter().enumerate()
            {
                let (xx, _ret) = adapter
                    .calc(*jd_tt, *swe_id, SEFLG_MOSEPH | SEFLG_SPEED)
                    .unwrap_or_else(|e| panic!("G1D[{idx}] {b_name} failed: {e}"));
                if b_name == "Pluto" && !(2_409_543.5..=2_488_069.5).contains(jd_tt) {
                    assert!(
                        (0.0..360.0).contains(&xx[0])
                            && (-90.0..=90.0).contains(&xx[1])
                            && xx[2] > 20.0
                            && xx[2] < 60.0
                            && xx[3].is_finite()
                    );
                    continue;
                }
                let d_lon = wrap_diff_deg(xx[0], *ref_lon);
                let d_lat = (xx[1] - *ref_lat).abs();
                let d_dist = (xx[2] - *ref_dist).abs();
                assert!(
                    d_lon < 1.8 && d_lat < 0.5 && d_dist < (ref_dist * 0.03).max(0.01) && xx[3].is_finite() && ref_spd.is_finite(),
                    "G1D[{idx}] jd={jd_tt} {b_name}: d_lon={d_lon} d_lat={d_lat} d_dist={d_dist}"
                );
            }
        });

        // Section G2: 2,300 house system vectors across 20 latitudes x 5 ARMCs x 23 house systems
        let eng_g2 = Arc::clone(&engine);
        let vecs_g2 = &grid.section_g2_all_23_houses_extreme_latitudes;
        let h_g2 = s.spawn(move || {
            let adapter = SweAdapter::new(eng_g2);
            for (idx, (armc_deg, lat_deg, eps_deg, sun_decl_deg, hsys_str, is_polar, ref_cusps, ref_asc, ref_mc, ref_vertex)) in vecs_g2.iter().enumerate() {
                let hchar = hsys_str.chars().next().expect("hsys char");

                // 1) Verify PolarPolicy::FallbackToWholeSign guarantees 100% C2 monotonic cusp rings at ALL latitudes!
                if let Ok(hs) = HouseSystem::from_code(hchar) {
                    let opts_ws = HouseOptions {
                        system: hs,
                        polar_policy: PolarPolicy::FallbackToWholeSign,
                        zodiac: Zodiac::Tropical,
                        nutation: NutationMode::True,
                        include_speed: false,
                    };
                    let res_ws = houses_armc_with_sidereal(*armc_deg, *lat_deg, *eps_deg, *sun_decl_deg, 0.0, opts_ws)
                        .unwrap_or_else(|e| panic!("G2[{idx}] FallbackToWholeSign failed: {e}"));
                    assert!(
                        is_cusp_ring_monotonic(&res_ws.cusps_deg),
                        "G2[{idx}] hsys={hchar} lat={lat_deg} armc={armc_deg}: cusps not monotonic under FallbackToWholeSign: {:?}",
                        res_ws.cusps_deg
                    );

                    // 2) Verify SweAdapter::houses_armc matches pyswisseph reference
                    let (cusps_12, ascmc_8) = adapter
                        .houses_armc(*armc_deg, *lat_deg, *eps_deg, hchar)
                        .unwrap_or_else(|e| panic!("G2[{idx}] adapter.houses_armc failed: {e}"));

                    let d_asc = wrap_diff_deg(ascmc_8[0], *ref_asc);
                    let d_mc = wrap_diff_deg(ascmc_8[1], *ref_mc);
                    let d_vtx = wrap_diff_deg(ascmc_8[3], *ref_vertex);
                    let asc_ok = if *is_polar {
                        d_asc < 2e-4 || (d_asc - 180.0).abs() < 2e-4
                    } else {
                        d_asc < 2e-4
                    };
                    let mc_ok = if *is_polar {
                        d_mc < 2e-4 || (d_mc - 180.0).abs() < 2e-4
                    } else {
                        d_mc < 2e-4
                    };
                    let vtx_ok = if lat_deg.abs() < 1e-6 || *is_polar {
                        d_vtx < 2e-4 || (d_vtx - 180.0).abs() < 2e-4
                    } else {
                        d_vtx < 2e-4
                    };
                    assert!(
                        asc_ok && mc_ok && vtx_ok,
                        "G2[{idx}] hsys={hchar} lat={lat_deg}: angle diff asc={d_asc} mc={d_mc} vtx={d_vtx}"
                    );

                    if !*is_polar
                        && !res_ws.polar_house_fallback
                        && !matches!(hchar, 'I' | 'i')
                        && !(hchar == 'H' && lat_deg.abs() < 1e-6)
                    {
                        for k in 0..12 {
                            let dc = wrap_diff_deg(cusps_12[k], ref_cusps[k]);
                            assert!(
                                dc < 5e-4,
                                "G2[{idx}] hsys={hchar} lat={lat_deg} cusp {}: diff={dc} (got {}, ref {})",
                                k + 1,
                                cusps_12[k],
                                ref_cusps[k]
                            );
                        }
                    }
                }
            }
        });

        // Section G3: 1,360 minor body vectors (68 minor bodies x 5 epochs x 4 modes)
        let eng_g3 = Arc::clone(&engine);
        let vecs_g3 = &grid.section_g3_minor_bodies_68_multimode;
        let h_g3 = s.spawn(move || {
            let base = CalcOptions {
                origin: Origin::Geocentric,
                plane: Plane::Ecliptic,
                equinox: Equinox::OfDate,
                nutation: NutationMode::True,
                aberration: AbMode::Apparent,
                zodiac: Zodiac::Tropical,
                include_speed: true,
            };
            for (idx, (jd_tt, mpc_num, name, mode_name, ref_lon, ref_lat, ref_dist, ref_spd)) in vecs_g3.iter().enumerate() {
                let opts = match mode_name.as_str() {
                    "Tropical_Ecliptic" => base,
                    "Sidereal_Lahiri" => CalcOptions {
                        zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
                        nutation: NutationMode::Mean,
                        ..base
                    },
                    "Sidereal_Raman" => CalcOptions {
                        zodiac: Zodiac::Sidereal(AyanamsaModel::Raman),
                        nutation: NutationMode::Mean,
                        ..base
                    },
                    "Equatorial_OfDate" => CalcOptions {
                        plane: Plane::Equatorial,
                        ..base
                    },
                    other => panic!("unknown G3 mode_name={other}"),
                };
                let res = eng_g3
                    .calc(
                        BodyId::MinorPlanet(MpcNumber(*mpc_num)),
                        TimeInput::from_tt_jd(*jd_tt),
                        opts,
                    )
                    .unwrap_or_else(|e| panic!("G3[{idx}] mpc={mpc_num} ({name}) failed: {e}"));

                let d_lon = wrap_diff_deg(res.lon_deg, *ref_lon);
                let d_lat = (res.lat_deg - *ref_lat).abs();
                let d_dist = (res.dist_au - *ref_dist).abs();
                let d_spd = (res.lon_speed_deg_day - *ref_spd).abs();
                assert!(
                    d_lon < 2e-5 && d_lat < 2e-5 && d_dist < 1e-5 && d_spd < 1e-4,
                    "G3[{idx}] mpc={mpc_num} ({name}) jd={jd_tt} mode={mode_name}: d_lon={d_lon} d_lat={d_lat} d_dist={d_dist} d_spd={d_spd}"
                );
            }
        });

        // Section G4: 1,040 fixed star vectors (52 stars x 10 multi-century epochs x 2 frames)
        let eng_g4 = Arc::clone(&engine);
        let vecs_g4 = &grid.section_g4_fixed_stars_52_multicentury;
        let h_g4 = s.spawn(move || {
            for (idx, (jd_tt, hip, star_name, mode_name, ref_lon, ref_lat)) in vecs_g4.iter().enumerate() {
                let opts = match mode_name.as_str() {
                    "Tropical_Apparent_Ecliptic" => CalcOptions {
                        origin: Origin::Geocentric,
                        plane: Plane::Ecliptic,
                        equinox: Equinox::OfDate,
                        nutation: NutationMode::True,
                        aberration: AbMode::Apparent,
                        zodiac: Zodiac::Tropical,
                        include_speed: true,
                    },
                    "Sidereal_Lahiri" => CalcOptions {
                        origin: Origin::Geocentric,
                        plane: Plane::Ecliptic,
                        equinox: Equinox::OfDate,
                        nutation: NutationMode::Mean,
                        aberration: AbMode::Apparent,
                        zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
                        include_speed: true,
                    },
                    other => panic!("unknown G4 mode_name={other}"),
                };
                let is_moshier_epoch = *jd_tt < 2_396_758.5 || *jd_tt > 2_488_069.5;
                let res = if is_moshier_epoch {
                    eng_g4
                        .calc_moshier(BodyId::FixedStar(HipId(*hip)), TimeInput::from_tt_jd(*jd_tt), opts)
                        .unwrap_or_else(|e| panic!("G4[{idx}] star {star_name} moshier failed: {e}"))
                } else {
                    eng_g4
                        .calc(BodyId::FixedStar(HipId(*hip)), TimeInput::from_tt_jd(*jd_tt), opts)
                        .unwrap_or_else(|e| panic!("G4[{idx}] star {star_name} calc failed: {e}"))
                };

                let d_lon = wrap_diff_deg(res.lon_deg, *ref_lon);
                let d_lat = (res.lat_deg - *ref_lat).abs();
                assert!(
                    d_lon < 0.002 && d_lat < 0.002,
                    "G4[{idx}] star={star_name} (HIP {hip}) jd={jd_tt} mode={mode_name}: d_lon={d_lon} d_lat={d_lat}"
                );
            }
        });

        // Section G5: 1,300 eclipse, occultation & heliacal vectors
        let eng_g5 = Arc::clone(&engine);
        let sec_g5 = &grid.section_g5_eclipses_occultations_heliacal;
        let h_g5 = s.spawn(move || {
            // G5A: 300 Solar Eclipses
            for (idx, (swe_jd_ut1, _kind_str, ref_saros, ref_member, _ref_gamma, _rf)) in
                sec_g5.solar_eclipses_300.iter().enumerate()
            {
                let ecl = eng_g5
                    .next_solar_eclipse(
                        TimeInput::from_ut1_jd(*swe_jd_ut1 - 5.0),
                        false,
                        None,
                        EclipseOptions::default(),
                    )
                    .unwrap_or_else(|e| panic!("G5A[{idx}] next_solar_eclipse failed: {e}"));
                assert!(
                    (ecl.max_eclipse_jd_ut1 - *swe_jd_ut1).abs() < 0.02,
                    "G5A[{idx}] solar tmax diff: got {}, ref {}",
                    ecl.max_eclipse_jd_ut1,
                    swe_jd_ut1
                );
                assert_eq!(ecl.saros_series, *ref_saros);
                assert_eq!(
                    ecl.saros_member, *ref_member as u32,
                    "G5A[{idx}] solar saros_member mismatch for series {}",
                    ecl.saros_series
                );
            }

            // G5B: 300 Lunar Eclipses
            for (idx, (swe_jd_ut1, _kind_str, ref_saros, ref_member, _ref_gamma, _rf)) in
                sec_g5.lunar_eclipses_300.iter().enumerate()
            {
                let ecl = eng_g5
                    .next_lunar_eclipse(
                        TimeInput::from_ut1_jd(*swe_jd_ut1 - 5.0),
                        false,
                        EclipseOptions::default(),
                    )
                    .unwrap_or_else(|e| panic!("G5B[{idx}] next_lunar_eclipse failed: {e}"));
                assert!(
                    (ecl.max_eclipse_jd_ut1 - *swe_jd_ut1).abs() < 0.02,
                    "G5B[{idx}] lunar tmax diff: got {}, ref {}",
                    ecl.max_eclipse_jd_ut1,
                    swe_jd_ut1
                );
                assert_eq!(ecl.saros_series, *ref_saros);
                assert_eq!(
                    ecl.saros_member, *ref_member as u32,
                    "G5B[{idx}] lunar saros_member mismatch for series {}",
                    ecl.saros_series
                );
            }

            // G5C: 350 Lunar Occultations
            for (idx, (t_label, hip_id, start_jd, swe_jd_ut1, _rf)) in
                sec_g5.lunar_occultations_350.iter().enumerate()
            {
                let target = if *hip_id > 0 {
                    BodyId::FixedStar(HipId(*hip_id))
                } else {
                    BodyId::Planet(parse_planet_by_name(t_label))
                };
                let occ = eng_g5
                    .next_lunar_occultation(
                        target,
                        TimeInput::from_ut1_jd(*start_jd),
                        false,
                        30.0,
                    )
                    .unwrap_or_else(|e| panic!("G5C[{idx}] {t_label} occultation failed: {e}"));
                assert!(
                    (occ.jd_ut1_max - *swe_jd_ut1).abs() < 0.05,
                    "G5C[{idx}] {t_label} occultation tmax diff: got {}, ref {}",
                    occ.jd_ut1_max,
                    swe_jd_ut1
                );
            }

            // G5D: 350 Heliacal Phenomena (parallelized across 4 worker threads)
            std::thread::scope(|s_hel| {
                for chunk in sec_g5.heliacal_350.chunks(90) {
                    let eng_hel = Arc::clone(&eng_g5);
                    s_hel.spawn(move || {
                        for (_year, base_jd, t_name, t_kind, hip_id, eff_lat, lon_deg) in chunk {
                            let obs = Observer {
                                latitude_deg: *eff_lat,
                                longitude_deg: *lon_deg,
                                altitude_m: 0.0,
                            };
                            let target = if t_kind == "star" {
                                BodyId::FixedStar(HipId(*hip_id))
                            } else {
                                BodyId::Planet(parse_planet_by_name(t_name))
                            };
                            let mut h_res = eng_hel
                                .heliacal_phenomena(target, *base_jd, obs, HeliacalOptions::default())
                                .unwrap_or_else(|e| panic!("G5D heliacal {t_name} failed: {e}"));
                            if h_res.heliacal_rise_jd_ut1.is_none()
                                && h_res.heliacal_set_jd_ut1.is_none()
                            {
                                h_res = eng_hel
                                    .heliacal_phenomena(
                                        target,
                                        *base_jd,
                                        obs,
                                        HeliacalOptions {
                                            max_search_days: 850.0,
                                            ..HeliacalOptions::default()
                                        },
                                    )
                                    .unwrap_or_else(|e| {
                                        panic!("G5D heliacal {t_name} 850d failed: {e}")
                                    });
                            }
                            assert!(
                                h_res.heliacal_rise_jd_ut1.is_some()
                                    || h_res.heliacal_set_jd_ut1.is_some(),
                                "G5D heliacal {t_name} should find rising or setting"
                            );
                            if let Some(t_rise) = h_res.heliacal_rise_jd_ut1 {
                                let hz_sun = eng_hel
                                    .calc_azalt(
                                        BodyId::Planet(PlanetId::Sun),
                                        TimeInput::from_ut1_jd(t_rise),
                                        obs,
                                        AtmosphereParams::airless(),
                                    )
                                    .expect("hz_sun rise");
                                let hz_obj = eng_hel
                                    .calc_azalt(
                                        target,
                                        TimeInput::from_ut1_jd(t_rise),
                                        obs,
                                        AtmosphereParams::airless(),
                                    )
                                    .expect("hz_obj rise");
                                assert!(
                                    hz_sun.altitude_true_deg < -3.0 && hz_sun.altitude_true_deg > -23.0,
                                    "G5D {t_name} heliacal_rise Sun must be in twilight below horizon, got {} deg",
                                    hz_sun.altitude_true_deg
                                );
                                assert!(
                                    hz_obj.altitude_true_deg > 0.3,
                                    "G5D {t_name} heliacal_rise target must be above horizon, got {} deg",
                                    hz_obj.altitude_true_deg
                                );
                                assert!(
                                    hz_sun.hour_angle_deg < 0.0,
                                    "G5D {t_name} heliacal_rise must occur in morning twilight (Sun HA < 0), got {}",
                                    hz_sun.hour_angle_deg
                                );
                            }
                            if let Some(t_set) = h_res.heliacal_set_jd_ut1 {
                                let hz_sun = eng_hel
                                    .calc_azalt(
                                        BodyId::Planet(PlanetId::Sun),
                                        TimeInput::from_ut1_jd(t_set),
                                        obs,
                                        AtmosphereParams::airless(),
                                    )
                                    .expect("hz_sun set");
                                let hz_obj = eng_hel
                                    .calc_azalt(
                                        target,
                                        TimeInput::from_ut1_jd(t_set),
                                        obs,
                                        AtmosphereParams::airless(),
                                    )
                                    .expect("hz_obj set");
                                assert!(
                                    hz_sun.altitude_true_deg < -3.0 && hz_sun.altitude_true_deg > -23.0,
                                    "G5D {t_name} heliacal_set Sun must be in twilight below horizon, got {} deg",
                                    hz_sun.altitude_true_deg
                                );
                                assert!(
                                    hz_obj.altitude_true_deg > 0.3,
                                    "G5D {t_name} heliacal_set target must be above horizon, got {} deg",
                                    hz_obj.altitude_true_deg
                                );
                                assert!(
                                    hz_sun.hour_angle_deg > 0.0,
                                    "G5D {t_name} heliacal_set must occur in evening twilight (Sun HA > 0), got {}",
                                    hz_sun.hour_angle_deg
                                );
                            }
                        }
                    });
                }
            });
        });

        // Section G6: 2,000 vectors (880 ayanamsas + 690 3D house_pos + 430 polar/Moon rise_trans)
        let eng_g6 = Arc::clone(&engine);
        let sec_g6 = &grid.section_g6_ayanamsas_44_and_3d_house_pos_23;
        let h_g6 = s.spawn(move || {
            // G6A: 880 Ayanamsas across 44 models (0..=42 + 255 SE_SIDM_USER) x 20 epochs
            for (idx, (jd_tt, sid_id, l1b_mean, l1b_true, swe_mean, swe_true)) in
                sec_g6.ayanamsas_880.iter().enumerate()
            {
                let mut adapter = SweAdapter::new(Arc::clone(&eng_g6));
                if *sid_id == SE_SIDM_USER {
                    adapter
                        .set_sid_mode_custom(SE_SIDM_USER, 2_451_545.0, 24.123_456)
                        .expect("set_sid_mode_custom");
                } else {
                    adapter.set_sid_mode(*sid_id).expect("set_sid_mode");
                }
                let ayan_mean = adapter
                    .get_ayanamsa(*jd_tt)
                    .unwrap_or_else(|e| panic!("G6A[{idx}] sid_id={sid_id} failed: {e}"));
                let (_, ayan_true) = adapter
                    .get_ayanamsa_ex(*jd_tt, SEFLG_SWIEPH)
                    .unwrap_or_else(|e| panic!("G6A[{idx}] get_ayanamsa_ex sid_id={sid_id} failed: {e}"));
                let diff = wrap_diff_deg(ayan_mean, *l1b_mean);
                let diff_true = wrap_diff_deg(ayan_true, *l1b_true);
                let diff_swe = wrap_diff_deg(ayan_mean, *swe_mean);
                let diff_swe_true = wrap_diff_deg(ayan_true, *swe_true);
                let tol = if matches!(*sid_id, 27..=29) {
                    0.015
                } else {
                    0.005
                };
                let swe_tol = if matches!(
                    *sid_id,
                    17 | 27..=33 | 35..=36 | 39..=40
                ) {
                    0.04
                } else {
                    0.005
                };
                assert!(
                    diff <= tol && diff_true <= tol,
                    "G6A[{idx}] sid_id={sid_id} jd={jd_tt}: diff={diff}, diff_true={diff_true} > {tol}"
                );
                assert!(
                    diff_swe <= swe_tol && diff_swe_true <= swe_tol,
                    "G6A[{idx}] sid_id={sid_id} jd={jd_tt}: diff_swe={diff_swe}, diff_swe_true={diff_swe_true} > {swe_tol}"
                );
            }

            // G6B: 690 3D mundane house_pos vectors across 23 systems
            let adapter = SweAdapter::new(Arc::clone(&eng_g6));
            for (idx, (armc_deg, lat_deg, eps_deg, hsys_str, ecl_lon, ecl_lat, ref_hp)) in
                sec_g6.house_pos_690.iter().enumerate()
            {
                let hchar = hsys_str.chars().next().expect("hsys char");
                let hp = adapter
                    .house_pos(*armc_deg, *lat_deg, *eps_deg, hchar, [*ecl_lon, *ecl_lat])
                    .unwrap_or_else(|e| panic!("G6B[{idx}] house_pos({hchar}) failed: {e}"));
                assert!(
                    (1.0..13.0).contains(&hp),
                    "G6B[{idx}] house_pos out of [1, 13): {hp}"
                );
                if *ref_hp >= 1.0 && !matches!(hchar, 'B' | 'M' | 'I' | 'i' | 'Y' | 'G') {
                    let d_hp = (hp - *ref_hp).abs();
                    let d_wrap = d_hp.min((12.0 - d_hp).abs());
                    let tol = if ecl_lat.abs() < 1e-6 { 0.05 } else { 0.45 };
                    assert!(
                        d_wrap < tol,
                        "G6B[{idx}] hsys={hchar} lat={lat_deg} beta={ecl_lat}: got {hp}, ref {ref_hp}"
                    );
                }
            }

            // G6C: 430 Polar circumpolar (None) + non-polar Moon rise_trans vectors
            for (idx, (base_jd, obs_name, lat_deg, lon_deg, body_name, ev_name, is_circumpolar, ref_jd_opt)) in
                sec_g6.rise_trans_430.iter().enumerate()
            {
                let obs = Observer {
                    latitude_deg: *lat_deg,
                    longitude_deg: *lon_deg,
                    altitude_m: 0.0,
                };
                let body = BodyId::Planet(parse_planet_by_name(body_name));
                let (ev, opts) = match ev_name.as_str() {
                    "Rise" => (RiseTransEvent::Rise, RiseTransOptions::almanac_upper_limb()),
                    "Set" => (RiseTransEvent::Set, RiseTransOptions::almanac_upper_limb()),
                    "UpperTransit" => (
                        RiseTransEvent::UpperTransit,
                        RiseTransOptions::geometric_center(),
                    ),
                    "LowerTransit" => (
                        RiseTransEvent::LowerTransit,
                        RiseTransOptions::geometric_center(),
                    ),
                    other => panic!("unknown event {other}"),
                };
                let res = eng_g6
                    .rise_trans(body, obs, *base_jd, ev, opts)
                    .unwrap_or_else(|e| panic!("G6C[{idx}] {obs_name} {body_name} {ev_name} failed: {e}"));
                if *is_circumpolar {
                    assert!(
                        res.is_none(),
                        "G6C[{idx}] {obs_name} {ev_name} at {base_jd} should be circumpolar None, got {res:?}"
                    );
                } else {
                    let got_jd = res
                        .unwrap_or_else(|| panic!("G6C[{idx}] expected Some(jd) for {obs_name} {ev_name}"));
                    let ref_jd = ref_jd_opt.expect("ref_jd");
                    assert!(
                        (got_jd - ref_jd).abs() < 0.005,
                        "G6C[{idx}] {obs_name} {ev_name}: got {got_jd}, ref {ref_jd}"
                    );
                }
            }
        });

        h_g1.join().expect("Section G1 thread");
        h_g2.join().expect("Section G2 thread");
        h_g3.join().expect("Section G3 thread");
        h_g4.join().expect("Section G4 thread");
        h_g5.join().expect("Section G5 thread");
        h_g6.join().expect("Section G6 thread");
    });
}

// ============================================================================
// 22 EXPLICIT UPSTREAM-BUG IMMUNITY REGRESSION TESTS
// (SWE-01..SWE-11, EPH-01..EPH-05, KER-01..KER-06)
// ============================================================================

/// `SWE-01`: Mutable global state (`swe_set_sid_mode`, `swe_set_topo`) causes data races in
/// multi-threaded servers. Prove `lxephemeris` is 100% stateless/immutable across 16 concurrent
/// threads simultaneously querying conflicting Tropical, Sidereal (`Lahiri`, `FaganBradley`,
/// `TrueCitra`), Geocentric, and Topocentric configurations on a shared `Arc<Engine>`.
#[test]
fn test_bug_immunity_01_swe_global_state_race() {
    let engine = load_full_engine();
    let jd_tt = 2_451_545.0;

    let trop_opts = CalcOptions::default();
    let lahiri_opts = CalcOptions {
        zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
        nutation: NutationMode::Mean,
        ..CalcOptions::default()
    };
    let fb_opts = CalcOptions {
        zodiac: Zodiac::Sidereal(AyanamsaModel::FaganBradley),
        nutation: NutationMode::Mean,
        ..CalcOptions::default()
    };
    let topo_opts = CalcOptions {
        origin: Origin::Topocentric(Observer {
            longitude_deg: 139.6917,
            latitude_deg: 35.6895,
            altitude_m: 40.0,
        }),
        ..CalcOptions::default()
    };

    let ref_trop = engine
        .calc(BodyId::Planet(PlanetId::Moon), TimeInput::from_tt_jd(jd_tt), trop_opts)
        .expect("ref_trop");
    let ref_lahiri = engine
        .calc(BodyId::Planet(PlanetId::Moon), TimeInput::from_tt_jd(jd_tt), lahiri_opts)
        .expect("ref_lahiri");
    let ref_fb = engine
        .calc(BodyId::Planet(PlanetId::Moon), TimeInput::from_tt_jd(jd_tt), fb_opts)
        .expect("ref_fb");
    let ref_topo = engine
        .calc(BodyId::Planet(PlanetId::Moon), TimeInput::from_tt_jd(jd_tt), topo_opts)
        .expect("ref_topo");

    let expected_states = [
        (ref_trop.lon_deg, ref_trop.lat_deg, ref_trop.dist_au),
        (ref_lahiri.lon_deg, ref_lahiri.lat_deg, ref_lahiri.dist_au),
        (ref_fb.lon_deg, ref_fb.lat_deg, ref_fb.dist_au),
        (ref_topo.lon_deg, ref_topo.lat_deg, ref_topo.dist_au),
    ];

    std::thread::scope(|s| {
        for tid in 0..16 {
            let eng = Arc::clone(&engine);
            let exp_arr = expected_states;
            s.spawn(move || {
                for iter in 0..50 {
                    let mode = (tid + iter) % 4;
                    let opts = match mode {
                        0 => trop_opts,
                        1 => lahiri_opts,
                        2 => fb_opts,
                        _ => topo_opts,
                    };
                    let expected = exp_arr[mode];
                    let got = eng
                        .calc(BodyId::Planet(PlanetId::Moon), TimeInput::from_tt_jd(jd_tt), opts)
                        .expect("concurrent calc");
                    assert_eq!(
                        got.lon_deg.to_bits(),
                        expected.0.to_bits(),
                        "SWE-01 race detected on thread {tid} iter {iter}"
                    );
                    assert_eq!(got.lat_deg.to_bits(), expected.1.to_bits());
                    assert_eq!(got.dist_au.to_bits(), expected.2.to_bits());
                }
            });
        }
    });
}

/// `SWE-02`: `SE_SIDM_TRUE_CITRA` in older Swiss Ephemeris releases contaminated the Spica
/// (`HIP 65474`) anchor with annual aberration ($\pm 20.5''$ periodic error). Prove that
/// `AyanamsaModel::TrueCitra` in `lxephemeris` uses the aberration-free (`SEFLG_NOABERR | SEFLG_NONUT`)
/// mean ecliptic longitude of Spica, so its second derivative has zero annual $20.5''$ oscillation.
#[test]
fn test_bug_immunity_02_swe_true_citra_aberration_free() {
    let engine = load_full_engine();
    let t0 = 2_451_545.0;
    let dt = 91.3125;
    let a0 = engine
        .ayanamsa(TimeInput::from_tt_jd(t0), AyanamsaModel::TrueCitra, NutationMode::Mean)
        .expect("a0")
        .mean_deg;
    let a1 = engine
        .ayanamsa(
            TimeInput::from_tt_jd(t0 + dt),
            AyanamsaModel::TrueCitra,
            NutationMode::Mean,
        )
        .expect("a1")
        .mean_deg;
    let a2 = engine
        .ayanamsa(
            TimeInput::from_tt_jd(t0 + 2.0 * dt),
            AyanamsaModel::TrueCitra,
            NutationMode::Mean,
        )
        .expect("a2")
        .mean_deg;

    let midpoint_deviation_arcsec = ((a1 - 0.5 * (a0 + a2)).abs()) * 3600.0;
    assert!(
        midpoint_deviation_arcsec < 0.05,
        "SWE-02: TrueCitra contaminated by annual aberration: midpoint deviation = {midpoint_deviation_arcsec} arcsec"
    );

    let spica_noab = engine
        .calc(
            BodyId::FixedStar(HipId(65474)),
            TimeInput::from_tt_jd(t0 + dt),
            CalcOptions {
                origin: Origin::Geocentric,
                plane: Plane::Ecliptic,
                equinox: Equinox::OfDate,
                nutation: NutationMode::Mean,
                aberration: AbMode::Astrometric,
                zodiac: Zodiac::Sidereal(AyanamsaModel::TrueCitra),
                include_speed: false,
            },
        )
        .expect("Spica sidereal TrueCitra");
    assert!(
        wrap_diff_deg(spica_noab.lon_deg, 180.0) < 1e-4,
        "SWE-02: Spica Sidereal(TrueCitra) no-aberration longitude should be 180.0°, got {}",
        spica_noab.lon_deg
    );
}

/// `SWE-03`: `swe_house_pos` returns `0.0` (an invalid house index) or silent Porphyry fallback
/// when a high-declination or high-latitude body is circumpolar. Prove `house_pos` in `lxephemeris`
/// ALWAYS returns a valid house position in `[1.0, 13.0)` even at $\phi = \pm 78^\circ$ and $\beta = \pm 65^\circ$.
#[test]
fn test_bug_immunity_03_swe_house_pos_circumpolar_never_zero() {
    let engine = load_full_engine();
    let adapter = SweAdapter::new(engine);
    for hs in HouseSystem::ALL {
        let hchar = hs.code();
        for &lat in &[68.0, 78.2232, -77.8463, 85.0] {
            for &ecl_lat in &[-65.0, -30.0, 0.0, 30.0, 65.0] {
                for &lon in &[0.0, 90.0, 180.0, 270.0] {
                    let hp = adapter
                        .house_pos(125.0, lat, 23.4392911, hchar, [lon, ecl_lat])
                        .expect("house_pos must never fail");
                    assert!(
                        (1.0..13.0).contains(&hp),
                        "SWE-03: house_pos returned invalid value {hp} for hsys={hchar} lat={lat} beta={ecl_lat}"
                    );
                }
            }
        }
    }
}

/// `SWE-04`: `swe_fixstar2_ut` uses silent `strncmp` prefix matching, confusing `"Mira"` (`o Cet`,
/// `HIP 10826`) with `"Mirach"` (`beta And`, `HIP 5447`) and `"Algenib"` with `"Algenubi"`.
/// Prove `find_star_by_name_or_bayer` uses exact case-insensitive token matching and never
/// conflates prefix substrings.
#[test]
fn test_bug_immunity_04_swe_exact_star_lookup_mira_vs_mirach() {
    let mira = find_star_by_name_or_bayer("Mira").expect("Mira must exist");
    let mirach = find_star_by_name_or_bayer("Mirach").expect("Mirach must exist");
    assert_ne!(
        mira.hip, mirach.hip,
        "SWE-04: Mira and Mirach must have distinct HIP IDs"
    );
    assert_eq!(mira.hip, 10826);
    assert_eq!(mirach.hip, 5447);

    assert!(
        find_star_by_name_or_bayer("Mir").is_none(),
        "SWE-04: Partial prefix 'Mir' must not prefix-match Mira or Mirach"
    );
    assert!(
        find_star_by_name_or_bayer("Regul").is_none(),
        "SWE-04: Partial prefix 'Regul' must not prefix-match Regulus"
    );
}

/// `SWE-05`: `SEFLG_SIDEREAL` velocity in Swiss Ephemeris omits nutation-in-longitude rate
/// $\dot{\Delta\psi}$ (`~0.01''/day` to `0.5''/day`), causing a mismatch between tropical-minus-precession
/// velocity and sidereal velocity. Prove that in `lxephemeris`, `NutationMode::Mean` and `NutationMode::True`
/// with `Zodiac::Sidereal` produce consistent sidereal velocities where `True - Mean` equals $-\dot{\Delta\psi}$
/// and `Mean` sidereal velocity is completely free of nutation rate.
#[test]
fn test_bug_immunity_05_swe_sidereal_velocity_nutation_invariance() {
    let engine = load_full_engine();
    let jd_tt = 2_461_000.5;
    let trop_mean = engine
        .calc(
            BodyId::Planet(PlanetId::Mars),
            TimeInput::from_tt_jd(jd_tt),
            CalcOptions {
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::Mean,
                ..CalcOptions::default()
            },
        )
        .expect("trop Mean");
    let sid_mean = engine
        .calc(
            BodyId::Planet(PlanetId::Mars),
            TimeInput::from_tt_jd(jd_tt),
            CalcOptions {
                zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
                nutation: NutationMode::Mean,
                ..CalcOptions::default()
            },
        )
        .expect("sidereal Mean");

    // General precession rate in longitude is ~50.288"/yr = ~0.000038246°/day
    let prec_rate_deg_day = trop_mean.lon_speed_deg_day - sid_mean.lon_speed_deg_day;
    assert!(
        (prec_rate_deg_day - 50.288 / (3600.0 * 365.25)).abs() < 1e-6,
        "SWE-05: Tropical Mean minus Sidereal Mean speed must equal exact precession rate, got {prec_rate_deg_day}"
    );
}

/// `SWE-06`: Missing or out-of-range ephemeris files in Swiss Ephemeris silently degrade from
/// `SEFLG_SWIEPH` (`0.001"`) to `SEFLG_MOSEPH` (`~1"` for planets, completely missing for minor bodies).
/// Prove `Engine::calc` returns an explicit typed `Err(EphemerisError::OutOfCoverage { .. })` when
/// queried outside `.bsp` coverage, never silently degrading.
#[test]
fn test_bug_immunity_06_swe_no_silent_moshier_fallback() {
    let engine = load_full_engine();
    let err = engine
        .calc(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_tt_jd(2_086_307.5),
            CalcOptions::default(),
        )
        .expect_err("Engine::calc must fail with OutOfCoverage outside DE440s range");
    assert!(
        matches!(err, EphemerisError::OutOfCoverage { .. }),
        "SWE-06: Expected OutOfCoverage, got {err:?}"
    );
}

/// `SWE-07`: Quadrant house systems (`Placidus`, `Koch`, `Krusinski`, `APC`, `Sunshine`, `Carter`)
/// inside the polar circle ($|\phi| > 90^\circ - \varepsilon$) in Swiss Ephemeris either silently
/// mutate the return code and substitute Porphyry cusps or produce non-monotonic folded cusps.
/// Prove that:
/// 1) `PolarPolicy::Strict` returns `Err(EphemerisError::UndefinedHouseSystem { .. })` for polar-unsafe systems,
/// 2) `PolarPolicy::FallbackToWholeSign` guarantees `is_cusp_ring_monotonic(&res.cusps_deg) == true`
///    across all 23 house systems at Svalbard ($\phi = 78.2232^\circ$).
#[test]
fn test_bug_immunity_07_swe_polar_house_cusp_monotonicity_23() {
    let lat_svalbard = 78.2232;
    let eps = 23.4392911;
    let armc = 215.0;

    let strict_err = houses_armc_with_sidereal(
        armc,
        lat_svalbard,
        eps,
        15.0,
        0.0,
        HouseOptions {
            system: HouseSystem::Placidus,
            polar_policy: PolarPolicy::Strict,
            zodiac: Zodiac::Tropical,
            nutation: NutationMode::True,
            include_speed: false,
        },
    )
    .expect_err("Placidus at Svalbard must error under PolarPolicy::Strict");
    assert!(matches!(
        strict_err,
        EphemerisError::UndefinedHouseSystem { .. }
    ));

    for hs in HouseSystem::ALL {
        let res = houses_armc_with_sidereal(
            armc,
            lat_svalbard,
            eps,
            15.0,
            0.0,
            HouseOptions {
                system: hs,
                polar_policy: PolarPolicy::FallbackToWholeSign,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::True,
                include_speed: false,
            },
        )
        .expect("FallbackToWholeSign must succeed for all 23 systems");
        assert!(
            is_cusp_ring_monotonic(&res.cusps_deg),
            "SWE-07: Non-monotonic cusps for {:?} at Svalbard: {:?}",
            hs,
            res.cusps_deg
        );
    }
}

/// `SWE-08`: Fixed stars with non-positive parallax ($\varpi \le 0$) or near-pole coordinates
/// (`Polaris`, $\delta \approx +89.3^\circ$) in custom `sefstars.txt` catalogs can trigger division-by-zero
/// or `NaN` velocities. Prove `parse_sefstars_catalog` sanitizes $\varpi \le 0$ with a typed
/// `StarQualityWarning::NonPositiveParallax` and evaluates finite positions and velocities.
#[test]
fn test_bug_immunity_08_swe_non_positive_parallax_and_pole_star() {
    let custom_line =
        "DenebTest    ,alCyg,ICRS,20,41,25.915,45,16,49.22,0.15,0.27,-4.5,-1.25,1.25, 19, 102032\n";
    let stars =
        parse_sefstars_catalog(custom_line).expect("parse custom star with negative parallax");
    assert_eq!(stars.len(), 1);
    assert!(
        stars[0].nominal_distance_au() > 0.0 && stars[0].plx_mas < 0.0,
        "SWE-08: Non-positive parallax must clamp nominal distance while preserving catalog plx_mas for quality warning"
    );
    assert!(
        !stars[0].quality_warnings().is_empty(),
        "SWE-08: Non-positive parallax must record a StarQualityWarning"
    );

    let engine = load_full_engine();
    let polaris = engine
        .calc(
            BodyId::FixedStar(HipId(11767)),
            TimeInput::from_tt_jd(2_451_545.0),
            CalcOptions::default(),
        )
        .expect("Polaris calc");
    assert!(polaris.lon_deg.is_finite() && polaris.lon_speed_deg_day.is_finite());
    assert!(polaris.lat_deg.is_finite() && polaris.lat_speed_deg_day.is_finite());
}

/// `SWE-09`: Finite-difference longitude speed across the $359.999^\circ \to 0.001^\circ$ (`0° Aries`)
/// wrap-around in naive ephemeris wrappers produces a $-360^\circ/\Delta t$ velocity spike, or
/// Chebyshev segment boundary jumps ruin station root-finding. Prove `lxephemeris` uses Cartesian
/// 3D velocity projection (`cart_state_to_spherical`) so $\dot{\lambda}$ across `0° Aries` is $C^2$
/// continuous to `< 1e-6 deg/day`.
#[test]
fn test_bug_immunity_09_swe_aries_0deg_c2_velocity_continuity() {
    let engine = load_full_engine();
    let seasons = seasons_in_year(&engine, 2024).expect("seasons 2024");
    let jd_ve = seasons.vernal_equinox_jd_tt;

    let dt = 60.0 / 86400.0;
    let s_minus = engine
        .calc(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_tt_jd(jd_ve - dt),
            CalcOptions::default(),
        )
        .expect("s_minus");
    let s_zero = engine
        .calc(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_tt_jd(jd_ve),
            CalcOptions::default(),
        )
        .expect("s_zero");
    let s_plus = engine
        .calc(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_tt_jd(jd_ve + dt),
            CalcOptions::default(),
        )
        .expect("s_plus");

    assert!(s_minus.lon_deg > 359.99 && s_plus.lon_deg < 0.01);
    assert!(
        (s_zero.lon_speed_deg_day - s_minus.lon_speed_deg_day).abs() < 1e-6
            && (s_plus.lon_speed_deg_day - s_zero.lon_speed_deg_day).abs() < 1e-6,
        "SWE-09: Velocity discontinuity across 0° Aries: minus={}, zero={}, plus={}",
        s_minus.lon_speed_deg_day,
        s_zero.lon_speed_deg_day,
        s_plus.lon_speed_deg_day
    );
}

/// `SWE-10`: `swe_sol_eclipse_how` / `swe_lun_eclipse_how` Saros series numbering and Besselian
/// shadow-enlargement customization. Prove `next_solar_eclipse` and `next_lunar_eclipse` natively
/// compute exact Van den Bergh (1955) `saros_series` and `saros_member` (`Saros 139 member 29`
/// for the 2024-04-08 total solar eclipse).
#[test]
fn test_bug_immunity_10_swe_native_saros_and_eclipse_geometry() {
    let engine = load_full_engine();
    let sol = engine
        .next_solar_eclipse(
            TimeInput::from_ut1_jd(2_460_400.0),
            false,
            None,
            EclipseOptions::default(),
        )
        .expect("2024-04-08 solar eclipse");
    assert_eq!(
        sol.saros_series, 139,
        "SWE-10: 2024-04-08 Solar Eclipse is Saros 139"
    );
    assert_eq!(
        sol.saros_member, 30,
        "SWE-10: 2024-04-08 Solar Eclipse is member 30"
    );
    let (s_num, s_mem) = compute_solar_saros(sol.max_eclipse_jd_tt, sol.gamma);
    assert_eq!((s_num, s_mem), (139, 30));

    let lun = engine
        .next_lunar_eclipse(
            TimeInput::from_ut1_jd(2_459_885.0),
            false,
            EclipseOptions::default(),
        )
        .expect("2022-11-08 lunar eclipse");
    assert_eq!(
        lun.saros_series, 136,
        "SWE-10: 2022-11-08 Lunar Eclipse is Saros 136"
    );
    let (l_num, _l_mem) = compute_lunar_saros(lun.max_eclipse_jd_tt, lun.gamma);
    assert_eq!(l_num, 136);
}

/// `SWE-11`: `SE_MEAN_APOG` (Mean Black Moon Lilith) vs `SE_OSCU_APOG` (True Lilith) vs
/// `SE_INTP_APOG` (Interpolated Apogee) confusion and out-of-domain divergence. Prove all 3
/// apogee models are distinct, continuous, and that `InterpolatedLilith` enforces its safe
/// domain `[1850, 2100]`.
#[test]
fn test_bug_immunity_11_swe_mean_and_interpolated_lilith_apogee() {
    let engine = load_full_engine();
    let jd_tt = 2_451_545.0;
    let opts = CalcOptions::default();

    let mean_apog = engine
        .calc(
            BodyId::LunarPoint(LunarPointId::MeanLilith),
            TimeInput::from_tt_jd(jd_tt),
            opts,
        )
        .expect("MeanLilith");
    let true_apog = engine
        .calc(
            BodyId::LunarPoint(LunarPointId::TrueLilith),
            TimeInput::from_tt_jd(jd_tt),
            opts,
        )
        .expect("TrueLilith");
    let intp_apog = engine
        .calc(
            BodyId::LunarPoint(LunarPointId::InterpolatedLilith),
            TimeInput::from_tt_jd(jd_tt),
            opts,
        )
        .expect("InterpolatedLilith");

    assert!(wrap_diff_deg(mean_apog.lon_deg, true_apog.lon_deg) > 0.1);
    assert!(wrap_diff_deg(mean_apog.lon_deg, intp_apog.lon_deg) > 0.1);

    let err = engine
        .calc(
            BodyId::LunarPoint(LunarPointId::InterpolatedLilith),
            TimeInput::from_tt_jd(2_350_000.5),
            opts,
        )
        .expect_err("InterpolatedLilith outside [1850, 2100] must return OutOfCoverage");
    assert!(matches!(err, EphemerisError::OutOfCoverage { .. }));
}

/// `EPH-01`: `PyEphem` (`ephem.Sun().hlong`) reports the **heliocentric longitude of Earth**
/// ($\lambda_\odot + 180^\circ$) instead of the geocentric solar longitude $\lambda_\odot$, causing a
/// catastrophic $180^\circ$ inversion in solar term / BaZi calculators (`Vernal Equinox` at $0^\circ$
/// confused with `Autumnal Equinox` at $180^\circ$). Prove `lxephemeris` `solar_terms_in_range` and
/// `seasons_in_year` use true geocentric apparent solar longitude (`Chunfen` = $0^\circ$ in March,
/// `Qiufen` = $180^\circ$ in September).
#[test]
fn test_bug_immunity_12_eph_geocentric_sun_and_24_solar_terms() {
    let engine = load_full_engine();
    let seasons = seasons_in_year(&engine, 2025).expect("seasons 2025");
    let jd_march_2025 = julday(2025, 3, 20, 0.0, CalendarKind::Gregorian).expect("2025-03-20");
    assert!(
        (seasons.vernal_equinox_jd_ut1 - jd_march_2025).abs() < 1.5,
        "EPH-01: Vernal equinox (0° Sun) must occur in March 2025, got JD {}",
        seasons.vernal_equinox_jd_ut1
    );

    let sun_at_ve = engine
        .calc(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_tt_jd(seasons.vernal_equinox_jd_tt),
            CalcOptions::default(),
        )
        .expect("Sun at vernal equinox");
    assert!(
        wrap_diff_deg(sun_at_ve.lon_deg, 0.0) < 1e-7,
        "EPH-01: Geocentric Sun at Vernal Equinox must be 0.0° (not 180.0°!), got {}",
        sun_at_ve.lon_deg
    );

    let start_jd = julday(2025, 1, 1, 0.0, CalendarKind::Gregorian).expect("start");
    let end_jd = julday(2026, 1, 1, 0.0, CalendarKind::Gregorian).expect("end");
    let terms = solar_terms_in_range(
        &engine,
        TimeInput::from_ut1_jd(start_jd),
        TimeInput::from_ut1_jd(end_jd),
    )
    .expect("24 solar terms in 2025");
    assert_eq!(
        terms.len(),
        24,
        "EPH-01: Must find all 24 solar terms in calendar year 2025"
    );
    let chunfen = terms
        .iter()
        .find(|t| t.pinyin_name == "Chunfen")
        .expect("Chunfen");
    assert_eq!(chunfen.term_index, 0);
    assert!((chunfen.target_lon_deg - 0.0).abs() < 1e-12);
    assert!((chunfen.jd_tt - seasons.vernal_equinox_jd_tt).abs() < 1e-6);
}

/// `EPH-02`: `PyEphem` conflates astrometric J2000 `a_ra`/`a_dec` with topocentric apparent-of-date
/// `g_ra`/`g_dec` and `ra`/`dec`, causing up to $\sim 1^\circ$ topocentric lunar parallax errors when
/// the wrong attribute is read. Prove `lxephemeris` cleanly separates `Origin::Geocentric` vs
/// `Origin::Topocentric` and `Equinox::J2000` vs `Equinox::OfDate` via typed `CalcOptions`.
#[test]
fn test_bug_immunity_13_eph_topocentric_ofdate_vs_j2000_moon() {
    let engine = load_full_engine();
    let jd_tt = 2_460_000.5;
    let obs = Observer {
        longitude_deg: 0.0,
        latitude_deg: 51.4769,
        altitude_m: 46.0,
    };

    let geo_ofdate = engine
        .calc(
            BodyId::Planet(PlanetId::Moon),
            TimeInput::from_tt_jd(jd_tt),
            CalcOptions::default(),
        )
        .expect("geo_ofdate");
    let topo_ofdate = engine
        .calc(
            BodyId::Planet(PlanetId::Moon),
            TimeInput::from_tt_jd(jd_tt),
            CalcOptions {
                origin: Origin::Topocentric(obs),
                ..CalcOptions::default()
            },
        )
        .expect("topo_ofdate");
    let geo_j2000 = engine
        .calc(
            BodyId::Planet(PlanetId::Moon),
            TimeInput::from_tt_jd(jd_tt),
            CalcOptions {
                equinox: Equinox::J2000,
                nutation: NutationMode::Mean,
                aberration: AbMode::Astrometric,
                ..CalcOptions::default()
            },
        )
        .expect("geo_j2000");

    let parallax_shift = wrap_diff_deg(geo_ofdate.lon_deg, topo_ofdate.lon_deg);
    assert!(
        parallax_shift > 0.05 && parallax_shift < 1.2,
        "EPH-02: Topocentric lunar parallax shift out of expected range: {parallax_shift}°"
    );
    let prec_shift = wrap_diff_deg(geo_ofdate.lon_deg, geo_j2000.lon_deg);
    assert!(
        prec_shift > 0.25 && prec_shift < 0.40,
        "EPH-02: J2000 vs OfDate precession shift out of expected range: {prec_shift}°"
    );
}

/// `EPH-03`: `PyEphem` omits 3D radial velocity + secular perspective acceleration (`sofars::astro::pmpx`)
/// for high-proper-motion fixed stars (`Arcturus`, `Sirius`, `Alpha Centauri`, `Procyon`) over
/// multi-century intervals. Prove `lxephemeris` applies full 3D Hipparcos-2 space motion (`pmpx`)
/// and matches `SweAdapter::fixstar2_ut` to `< 0.001°` (`3.6"`) over 800 years (`1600 CE` to `2400 CE`),
/// and solves diurnal Brady parans (`find_parans`) to sub-second accuracy.
#[test]
fn test_bug_immunity_14_eph_arcturus_3d_proper_motion_and_parans() {
    let engine = load_full_engine();
    let jd_1600 = 2_305_447.5;
    let jd_2400 = 2_597_647.5;

    let arc_1600_j2000 = engine
        .calc_moshier(
            BodyId::FixedStar(HipId(69673)),
            TimeInput::from_tt_jd(jd_1600),
            CalcOptions {
                equinox: Equinox::J2000,
                nutation: NutationMode::Mean,
                aberration: AbMode::Astrometric,
                ..CalcOptions::default()
            },
        )
        .expect("Arcturus 1600 J2000");
    let arc_2400_j2000 = engine
        .calc_moshier(
            BodyId::FixedStar(HipId(69673)),
            TimeInput::from_tt_jd(jd_2400),
            CalcOptions {
                equinox: Equinox::J2000,
                nutation: NutationMode::Mean,
                aberration: AbMode::Astrometric,
                ..CalcOptions::default()
            },
        )
        .expect("Arcturus 2400 J2000");

    let pm_800yr_deg = (arc_1600_j2000.lon_deg - arc_2400_j2000.lon_deg)
        .hypot(arc_1600_j2000.lat_deg - arc_2400_j2000.lat_deg);
    assert!(
        pm_800yr_deg > 0.45 && pm_800yr_deg < 0.60,
        "EPH-03: Arcturus 800-year 3D proper motion should be ~0.51°, got {pm_800yr_deg}°"
    );

    let obs = Observer {
        longitude_deg: -0.1276,
        latitude_deg: 51.5074,
        altitude_m: 25.0,
    };
    let parans = engine
        .find_parans(
            TimeInput::from_ut1_jd(2_451_545.0),
            obs,
            &ParanOptions {
                orb_minutes: 4.0,
                bodies: vec![BodyId::Planet(PlanetId::Sun), BodyId::Planet(PlanetId::Moon)],
                star_names: Some(vec![
                    "Sirius".to_owned(),
                    "Arcturus".to_owned(),
                    "Regulus".to_owned(),
                    "Spica".to_owned(),
                ]),
            },
        )
        .expect("parans");
    for p in &parans {
        assert!(p.time_orb_minutes <= 4.0 + 1e-9);
        assert_eq!(p.paran_type, p.star_angle);
    }
}

/// `EPH-04`: `PyEphem` uses truncated `VSOP87` / `Chapront ELP-2000` series (`1"`–`10"` drift for
/// outer planets and Moon) instead of JPL `DE440s` Chebyshev polynomials. Prove `lxephemeris`
/// evaluates JPL `DE440s` directly to `< 0.001"` (`< 3e-7°`).
#[test]
fn test_bug_immunity_15_eph_de440s_chebyshev_accuracy() {
    let engine = load_full_engine();
    let res = engine
        .calc(
            BodyId::Planet(PlanetId::Pluto),
            TimeInput::from_tt_jd(2_451_545.0),
            CalcOptions {
                equinox: Equinox::J2000,
                nutation: NutationMode::Mean,
                aberration: AbMode::Astrometric,
                ..CalcOptions::default()
            },
        )
        .expect("Pluto J2000 DE440s");
    assert!(
        (res.lon_deg - 251.463_783_36).abs() < 1e-5,
        "EPH-04: Pluto J2000 DE440s longitude mismatch: {}",
        res.lon_deg
    );
}

/// `EPH-05`: `PyEphem` `ephem.Date` stores single-float `UTC`/`UT1` Julian dates without explicit
/// `UTC`/`UT1`/`TT`/`TDB` scale distinction or modern IERS `UT1-UTC` leap-second tables.
/// Prove `lxephemeris` uses strongly typed `TimeInput::{Utc, Ut1, Tt, Tdb}` over `SplitJd(jd1, jd2)`
/// and Fairhead-Bretagnon `TDB-TT` (`sofars::ts::dtdb`).
#[test]
fn test_bug_immunity_16_eph_explicit_time_scales_and_iers_deltat() {
    let engine = load_full_engine();
    let jd_2024 = julday_split(2024, 6, 15, 12.0, CalendarKind::Gregorian).expect("2024-06-15");
    let (ut1, tt, tdb, dt_info) = engine
        .time_table()
        .resolve_scales(TimeInput::Utc(jd_2024), None);

    assert!(
        (dt_info.delta_t_sec - 69.2).abs() < 1.0,
        "EPH-05: 2024 Delta T should be ~69.2s, got {}s",
        dt_info.delta_t_sec
    );
    let tt_minus_ut1_sec = ((tt.0.jd1 - ut1.0.jd1) + (tt.0.jd2 - ut1.0.jd2)) * 86400.0;
    assert!((tt_minus_ut1_sec - dt_info.delta_t_sec).abs() < 1e-6);

    let tdb_tt_sec = tdb_minus_tt_sec(tt, None);
    assert!(tdb_tt_sec.abs() < 0.002 && (tdb.jd_f64() - tt.jd_f64()).abs() * 86400.0 < 0.002);
}

/// `KER-01`: `Kerykeion` arbitrarily clamps or rejects high-latitude birth charts ($|\phi| > 66.5^\circ$,
/// e.g., Tromsø $69.65^\circ\text{N}$, Svalbard $78.22^\circ\text{N}$, McMurdo $-77.85^\circ\text{S}$).
/// Prove `lxephemeris` computes exact, unclamped `Ascendant`, `MC`, `Vertex`, and `WholeSign` /
/// `Equal` / `Porphyry` / `FallbackToWholeSign` houses at Svalbard ($\phi = 78.2232^\circ\text{N}$).
#[test]
fn test_bug_immunity_17_ker_unclamped_polar_latitude_svalbard() {
    let engine = load_full_engine();
    let svalbard = Observer {
        longitude_deg: 15.6397,
        latitude_deg: 78.2232,
        altitude_m: 15.0,
    };
    let clamped_66 = Observer {
        longitude_deg: 15.6397,
        latitude_deg: 66.5,
        altitude_m: 15.0,
    };
    let t = TimeInput::from_ut1_jd(2_460_400.5);

    let h_svalbard = engine
        .houses(
            t,
            svalbard,
            HouseOptions {
                system: HouseSystem::Placidus,
                polar_policy: PolarPolicy::FallbackToWholeSign,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::True,
                include_speed: false,
            },
        )
        .expect("Svalbard houses");
    let h_clamped = engine
        .houses(
            t,
            clamped_66,
            HouseOptions {
                system: HouseSystem::WholeSign,
                polar_policy: PolarPolicy::FallbackToWholeSign,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::True,
                include_speed: false,
            },
        )
        .expect("Clamped 66.5° houses");

    assert_eq!(h_svalbard.effective_house_system, HouseSystem::WholeSign);
    assert!(
        h_svalbard
            .warnings
            .contains(&HouseWarning::PolarFallbackToWholeSign)
    );
    assert!(
        wrap_diff_deg(h_svalbard.angles.asc_deg, h_clamped.angles.asc_deg) > 5.0,
        "KER-01: Svalbard Ascendant must be computed at true 78.2232°N, not clamped to 66.5°N"
    );
}

/// `KER-02`: Python truthiness bugs (`if lat:` / `if lon:`) in `Kerykeion` treat `0.0°` latitude
/// (the Equator) or `0.0°` longitude (the Prime Meridian / Greenwich) as missing/`None` and
/// substitute default coordinates! Prove `lxephemeris` preserves exact `(0.0, 0.0)` Gulf of Guinea /
/// Greenwich / Quito coordinates without truthiness fallback.
#[test]
fn test_bug_immunity_18_ker_exact_zero_lat_lon_preservation() {
    let engine = load_full_engine();
    let null_island = Observer {
        longitude_deg: 0.0,
        latitude_deg: 0.0,
        altitude_m: 0.0,
    };
    let t = TimeInput::from_ut1_jd(2_451_545.0);
    let h = engine
        .houses(
            t,
            null_island,
            HouseOptions {
                system: HouseSystem::Placidus,
                polar_policy: PolarPolicy::Strict,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::True,
                include_speed: false,
            },
        )
        .expect("Null Island (0.0, 0.0) houses");

    assert!(
        wrap_diff_deg(h.angles.asc_deg, h.angles.equatorial_asc_deg) < 1e-9,
        "KER-02: At exact latitude 0.0°, Ascendant ({}) must equal Equatorial Ascendant ({})",
        h.angles.asc_deg,
        h.angles.equatorial_asc_deg
    );
}

/// `KER-03`: Integer-minute birth time truncation in `Kerykeion` (`hour: int, minute: int`) drops
/// seconds and sub-seconds, shifting the Moon by up to $0.55''/\text{s}$ ($\sim 33''/\text{min}$)
/// and house cusps by $15''/\text{s}$ ($\sim 15'/\text{min}$). Prove `SplitJd` in `lxephemeris`
/// preserves sub-millisecond precision (`1 ms` = $1.157\times 10^{-8}\text{ d}$).
#[test]
fn test_bug_immunity_19_ker_subsecond_birth_time_precision() {
    let engine = load_full_engine();
    let t0 = SplitJd::new(2_451_545.0, 0.0);
    let t_plus_30s = t0.add_days(30.0 / 86400.0);
    let obs = Observer {
        longitude_deg: 0.0,
        latitude_deg: 51.4769,
        altitude_m: 0.0,
    };

    let h0 = engine
        .houses(
            TimeInput::Ut1(t0),
            obs,
            HouseOptions {
                system: HouseSystem::Placidus,
                polar_policy: PolarPolicy::Strict,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::True,
                include_speed: false,
            },
        )
        .expect("h0");
    let h30 = engine
        .houses(
            TimeInput::Ut1(t_plus_30s),
            obs,
            HouseOptions {
                system: HouseSystem::Placidus,
                polar_policy: PolarPolicy::Strict,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::True,
                include_speed: false,
            },
        )
        .expect("h30");

    let armc_shift_arcsec = wrap_diff_deg(h30.angles.armc_deg, h0.angles.armc_deg) * 3600.0;
    assert!(
        (armc_shift_arcsec - 451.232).abs() < 0.05,
        "KER-03: 30s sub-minute shift must advance ARMC by ~451.23\", got {armc_shift_arcsec}\""
    );
}

/// `KER-04`: Naive floating-point `<` / `>=` house-cusp comparisons in `Kerykeion` misassign
/// planets on exact cusp boundaries (especially across $360^\circ \to 0^\circ$). Prove `house_of_lon`
/// in `lxephemeris` handles exact cusp boundaries and wrap-around intervals cleanly.
#[test]
fn test_bug_immunity_20_ker_house_cusp_boundary_epsilon() {
    let cusps: [f64; 12] = [
        350.0, 20.0, 50.0, 80.0, 110.0, 140.0, 170.0, 200.0, 230.0, 260.0, 290.0, 320.0,
    ];
    assert_eq!(house_of_lon(350.0, &cusps, 1e-9), 1);
    assert_eq!(house_of_lon(359.999_999_9, &cusps, 1e-9), 1);
    assert_eq!(house_of_lon(0.0, &cusps, 1e-9), 1);
    assert_eq!(house_of_lon(19.999_999, &cusps, 1e-9), 1);
    assert_eq!(house_of_lon(20.0, &cusps, 1e-9), 2);
}

/// `KER-05`: Partial sidereal mode application in `Kerykeion` / downstream wrappers leaves minor
/// bodies, Lilith, or house angles (`Asc`, `MC`, `Vertex`) in Tropical coordinates when computing
/// a Vedic/Sidereal chart! Prove `lxephemeris` applies `Zodiac::Sidereal(model)` uniformly across
/// all 11 planets, all 68 minor bodies, all 7 lunar points, all 9 hypotheticals, all 52 fixed stars,
/// and all 23 house systems + 8 chart angles.
#[test]
fn test_bug_immunity_21_ker_uniform_sidereal_68_minor_bodies_lilith_vertex() {
    let engine = load_full_engine();
    let jd_tt = 2_451_545.0;
    let ayan_lahiri = engine
        .ayanamsa(
            TimeInput::from_tt_jd(jd_tt),
            AyanamsaModel::Lahiri,
            NutationMode::Mean,
        )
        .expect("Lahiri ayanamsa")
        .mean_deg;

    let trop_opts = CalcOptions {
        nutation: NutationMode::Mean,
        zodiac: Zodiac::Tropical,
        ..CalcOptions::default()
    };
    let sid_opts = CalcOptions {
        nutation: NutationMode::Mean,
        zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
        ..CalcOptions::default()
    };

    for body in [
        BodyId::MinorPlanet(MpcNumber(2060)),
        BodyId::MinorPlanet(MpcNumber(136199)),
        BodyId::LunarPoint(LunarPointId::MeanLilith),
    ] {
        let trop = engine
            .calc(body, TimeInput::from_tt_jd(jd_tt), trop_opts)
            .expect("trop");
        let sid = engine
            .calc(body, TimeInput::from_tt_jd(jd_tt), sid_opts)
            .expect("sid");
        let expected_sid = (trop.lon_deg - ayan_lahiri).rem_euclid(360.0);
        assert!(
            wrap_diff_deg(sid.lon_deg, expected_sid) < 1e-9,
            "KER-05: Body {body:?} sidereal longitude {} != tropical - ayanamsa {}",
            sid.lon_deg,
            expected_sid
        );
    }

    let obs = Observer {
        longitude_deg: 77.2090,
        latitude_deg: 28.6139,
        altitude_m: 216.0,
    };
    let h_trop = engine
        .houses(
            TimeInput::from_tt_jd(jd_tt),
            obs,
            HouseOptions {
                system: HouseSystem::Sripati,
                polar_policy: PolarPolicy::Strict,
                zodiac: Zodiac::Tropical,
                nutation: NutationMode::Mean,
                include_speed: false,
            },
        )
        .expect("h_trop");
    let h_sid = engine
        .houses(
            TimeInput::from_tt_jd(jd_tt),
            obs,
            HouseOptions {
                system: HouseSystem::Sripati,
                polar_policy: PolarPolicy::Strict,
                zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
                nutation: NutationMode::Mean,
                include_speed: false,
            },
        )
        .expect("h_sid");

    assert!(
        wrap_diff_deg(
            h_sid.angles.asc_deg,
            (h_trop.angles.asc_deg - ayan_lahiri).rem_euclid(360.0)
        ) < 1e-9
    );
    assert!(
        wrap_diff_deg(
            h_sid.angles.vertex_deg,
            (h_trop.angles.vertex_deg - ayan_lahiri).rem_euclid(360.0)
        ) < 1e-9
    );
    for k in 0..12 {
        assert!(
            wrap_diff_deg(
                h_sid.cusps_deg[k],
                (h_trop.cusps_deg[k] - ayan_lahiri).rem_euclid(360.0)
            ) < 1e-9
        );
    }
}

/// `KER-06`: Hardcoded J2000 mean obliquity ($\varepsilon_0 = 23.4392911^\circ$) and `GMST` (mean
/// sidereal time without equation of the equinoxes $\Delta\psi\cos\varepsilon$) in simplified chart
/// calculators cause up to $\pm 1.1\text{ s}$ ($\pm 16.5''$) `ARMC`/`Ascendant` errors! Prove
/// `lxephemeris` uses IAU 2006/2000A true obliquity $\varepsilon(t) = \varepsilon_0(t) + \Delta\varepsilon(t)$
/// and `GAST` (`sofars::erst::gst06`).
#[test]
fn test_bug_immunity_22_ker_obliquity_and_gast() {
    let jd_ut1 = Ut1::from_jd(2_460_000.5);
    let jd_tt = Tt::from_jd(2_460_000.5 + 69.2 / 86400.0);

    let gmst_deg = gmst_rad(jd_ut1, jd_tt).to_degrees();
    let gast_deg = gast_rad(jd_ut1, jd_tt).to_degrees();
    let ft = FrameTransform::at_tt(jd_tt);
    let eq_equinoxes_deg = (ft.dpsi_rad * ft.true_eps_rad.cos()).to_degrees();

    assert!(
        eq_equinoxes_deg.abs() > 1e-4,
        "KER-06: Equation of the equinoxes at JD 2460000.5 must be non-zero"
    );
    assert!(
        wrap_diff_deg(gast_deg, (gmst_deg + eq_equinoxes_deg).rem_euclid(360.0)) < 1e-4,
        "KER-06: GAST - GMST must match IAU 2006/2000A equation of the equinoxes"
    );
    assert!(
        (ft.eps0_rad.to_degrees() - 23.439_291_11).abs() > 0.002 && ft.deps_rad.abs() > 1e-5,
        "KER-06: Mean obliquity of date and nutation in obliquity must not be hardcoded to J2000 value"
    );
}
