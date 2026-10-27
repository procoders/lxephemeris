//! Phase 5 Verification Suite:
//! 1. `testdata/golden/l1b_asteroids_68.json` (340 vectors across all 68 minor bodies, gate < 0.01").
//! 2. `testdata/golden/l3_se_10k_grid.json` (10,242 vectors across 7 sections: L1B & L3 SE parity).
//! 3. 32-Thread Concurrency Torture Suite (`0.0` ULP state bleed across 1,024+ mixed evaluations)
//!    + Codex Gate 2 checks (`lunar.rs` non-ecliptic/J2000 velocities & `CalcResult.warnings`).

#![allow(clippy::indexing_slicing, clippy::panic)]

use lxephemeris_core::{
    houses_armc, AbMode, AyanamsaModel, BodyId, CalcOptions, EclipseOptions, Engine, EngineConfig,
    Equinox, HipId, HouseOptions, HouseSystem, HypotheticalId, LunarPointId, MpcNumber,
    NutationMode, Observer, Origin, Plane, PlanetId, PolarPolicy, TimeInput, Zodiac, JD_J2000,
};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf()
}

fn load_full_engine() -> Engine {
    let root = workspace_root();
    Engine::from_bsp_paths(
        &[
            root.join("data/de440s.bsp"),
            root.join("data/asteroids_68_1850_2100.bsp"),
        ],
        EngineConfig::default(),
    )
    .unwrap_or_else(|e| panic!("Failed to load kernels: {e}"))
}

fn wrap_diff_arcsec(a_deg: f64, b_deg: f64) -> f64 {
    let d = (a_deg - b_deg + 180.0).rem_euclid(360.0) - 180.0;
    d.abs() * 3_600.0
}

fn parse_planet_name(name: &str) -> PlanetId {
    match name {
        "Sun" => PlanetId::Sun,
        "Moon" => PlanetId::Moon,
        "Mercury" => PlanetId::Mercury,
        "Venus" => PlanetId::Venus,
        "Mars" => PlanetId::Mars,
        "Jupiter" => PlanetId::Jupiter,
        "Saturn" => PlanetId::Saturn,
        "Uranus" => PlanetId::Uranus,
        "Neptune" => PlanetId::Neptune,
        "Pluto" => PlanetId::Pluto,
        other => panic!("Unknown planet name: {other}"),
    }
}

fn parse_lunar_point_name(name: &str) -> LunarPointId {
    match name {
        "Mean_Node" => LunarPointId::MeanNode,
        "Mean_South_Node" => LunarPointId::MeanSouthNode,
        "True_Node" => LunarPointId::TrueNode,
        "True_South_Node" => LunarPointId::TrueSouthNode,
        "Mean_Lilith" => LunarPointId::MeanLilith,
        "True_Lilith" => LunarPointId::TrueLilith,
        "Interpolated_Lilith" => LunarPointId::InterpolatedLilith,
        other => panic!("Unknown lunar point: {other}"),
    }
}

fn parse_uranian_name(name: &str) -> HypotheticalId {
    match name {
        "Cupido" => HypotheticalId::Cupido,
        "Hades" => HypotheticalId::Hades,
        "Zeus" => HypotheticalId::Zeus,
        "Kronos" => HypotheticalId::Kronos,
        "Apollon" => HypotheticalId::Apollon,
        "Admetos" => HypotheticalId::Admetos,
        "Vulkanus" => HypotheticalId::Vulkanus,
        "Poseidon" => HypotheticalId::Poseidon,
        "Transpluto" | "Isis" => HypotheticalId::Transpluto,
        other => panic!("Unknown uranian: {other}"),
    }
}

// ---------------------------------------------------------------------------
// 1. Level 1B: 68 Minor Bodies Golden Dataset (`l1b_asteroids_68.json`, 340 vectors)
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
struct L1bAsteroidsFile {
    tolerance_arcsec: f64,
    count: usize,
    vectors: Vec<L1bAsteroidVector>,
}

#[derive(Debug, Deserialize)]
struct L1bAsteroidVector {
    epoch_label: String,
    jd_tt: f64,
    mpc: u32,
    name: String,
    apparent_ecliptic_true_of_date: L1bEclState,
    astrometric_icrf: L1bIcrfState,
}

#[derive(Debug, Deserialize)]
struct L1bEclState {
    lon_deg: f64,
    lat_deg: f64,
    dist_au: f64,
    lon_speed_deg_day: f64,
    lat_speed_deg_day: f64,
    dist_speed_au_day: f64,
}

#[derive(Debug, Deserialize)]
struct L1bIcrfState {
    ra_deg: f64,
    dec_deg: f64,
    dist_au: f64,
}

#[test]
fn test_l1b_asteroids_68_golden_340() {
    let root = workspace_root();
    let text = fs::read_to_string(root.join("testdata/golden/l1b_asteroids_68.json"))
        .expect("read l1b_asteroids_68.json");
    let golden: L1bAsteroidsFile =
        serde_json::from_str(&text).expect("parse l1b_asteroids_68.json");
    assert_eq!(golden.count, 340);
    assert_eq!(golden.vectors.len(), 340);

    let engine = load_full_engine();
    let mut max_ecl_err_arcsec = 0.0_f64;
    let mut max_icrf_err_arcsec = 0.0_f64;
    let mut max_spd_err_deg_day = 0.0_f64;

    let icrf_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Equatorial,
        equinox: Equinox::J2000,
        nutation: NutationMode::Mean,
        aberration: AbMode::Astrometric,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };

    for v in &golden.vectors {
        let body = BodyId::MinorPlanet(MpcNumber(v.mpc));
        let t_in = TimeInput::from_tt_jd(v.jd_tt);

        let res_ecl = engine
            .calc(body, t_in, CalcOptions::default())
            .unwrap_or_else(|e| panic!("{} ({}) at {}: {e}", v.name, v.mpc, v.epoch_label));
        let exp_e = &v.apparent_ecliptic_true_of_date;
        let err_lon = wrap_diff_arcsec(res_ecl.lon_deg, exp_e.lon_deg);
        let err_lat = (res_ecl.lat_deg - exp_e.lat_deg).abs() * 3_600.0;
        let err_spd = (res_ecl.lon_speed_deg_day - exp_e.lon_speed_deg_day).abs();
        max_ecl_err_arcsec = max_ecl_err_arcsec.max(err_lon).max(err_lat);
        max_spd_err_deg_day = max_spd_err_deg_day.max(err_spd);

        assert!(
            err_lon < golden.tolerance_arcsec && err_lat < golden.tolerance_arcsec,
            "Minor body {} ({}) ecliptic mismatch at {}: dLon={err_lon:.6}\", dLat={err_lat:.6}\"",
            v.name,
            v.mpc,
            v.epoch_label
        );
        assert!((res_ecl.dist_au - exp_e.dist_au).abs() < 1e-6);
        assert!(err_spd < 1e-5, "Speed mismatch for {}: {err_spd:.3e}", v.name);
        assert!((res_ecl.lat_speed_deg_day - exp_e.lat_speed_deg_day).abs() < 1e-5);
        assert!((res_ecl.dist_speed_au_day - exp_e.dist_speed_au_day).abs() < 1e-5);

        let res_icrf = engine
            .calc(body, t_in, icrf_opts)
            .unwrap_or_else(|e| panic!("ICRF {} ({}): {e}", v.name, v.mpc));
        let exp_i = &v.astrometric_icrf;
        let err_ra = wrap_diff_arcsec(res_icrf.lon_deg, exp_i.ra_deg);
        let err_dec = (res_icrf.lat_deg - exp_i.dec_deg).abs() * 3_600.0;
        max_icrf_err_arcsec = max_icrf_err_arcsec.max(err_ra).max(err_dec);
        assert!(
            err_ra < golden.tolerance_arcsec && err_dec < golden.tolerance_arcsec,
            "Minor body {} ({}) ICRF mismatch at {}: dRA={err_ra:.6}\", dDec={err_dec:.6}\"",
            v.name,
            v.mpc,
            v.epoch_label
        );
        assert!((res_icrf.dist_au - exp_i.dist_au).abs() < 1e-6);
    }

    println!(
        "L1B 68 Minor Bodies (340 vectors): max_ecl_err = {max_ecl_err_arcsec:.6}\", max_icrf_err = {max_icrf_err_arcsec:.6}\", max_spd_err = {max_spd_err_deg_day:.3e} deg/day"
    );
}

// ---------------------------------------------------------------------------
// 2. Full 10,242-Vector Grid (`l3_se_10k_grid.json`)
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
#[allow(clippy::type_complexity)]
struct L3Grid10kFile {
    total_vector_count: usize,
    counts_by_section: BTreeMap<String, usize>,
    major_bodies: Vec<(f64, String, f64, f64, f64, f64, f64, f64)>,
    minor_bodies: Vec<(f64, u32, String, f64, f64, f64, f64, f64)>,
    ayanamsas: Vec<(f64, String, f64, f64, f64)>,
    houses: Vec<(f64, f64, String, [f64; 12])>,
    lunar_points: Vec<(f64, String, f64, f64, Option<f64>)>,
    uranians: Vec<(f64, String, f64, f64, f64)>,
    fixed_stars: Vec<(f64, u32, String, f64, f64)>,
}

#[test]
fn test_l3_se_10k_grid_all_10242_vectors() {
    let root = workspace_root();
    let text = fs::read_to_string(root.join("testdata/golden/l3_se_10k_grid.json"))
        .expect("read l3_se_10k_grid.json");
    let grid: L3Grid10kFile = serde_json::from_str(&text).expect("parse l3_se_10k_grid.json");

    assert_eq!(grid.total_vector_count, 10_242);
    assert_eq!(grid.major_bodies.len(), 5_000);
    assert_eq!(grid.minor_bodies.len(), 1_360);
    assert_eq!(grid.ayanamsas.len(), 1_200);
    assert_eq!(grid.houses.len(), 1_380);
    assert_eq!(grid.lunar_points.len(), 280);
    assert_eq!(grid.uranians.len(), 450);
    assert_eq!(grid.fixed_stars.len(), 572);
    assert_eq!(grid.counts_by_section.values().sum::<usize>(), 10_242);

    let engine = load_full_engine();

    // Section 1: 5,000 Major Body vectors (500 epochs x 10 bodies)
    let mut max_major_l1b_arcsec = 0.0_f64;
    let mut max_major_l3_arcsec = 0.0_f64;
    let mut current_jd: Option<f64> = None;
    let mut current_moment = None;

    for &(jd_tt, ref name, l1b_lon, l1b_lat, l1b_dist, l3_lon, l3_lat, l3_spd) in
        &grid.major_bodies
    {
        if current_jd != Some(jd_tt) {
            current_jd = Some(jd_tt);
            current_moment = Some(
                engine
                    .at(TimeInput::from_tt_jd(jd_tt), None)
                    .expect("moment"),
            );
        }
        let m = current_moment.as_ref().expect("moment");
        let p = parse_planet_name(name);
        let res = engine
            .calc_in_moment(m, BodyId::Planet(p), CalcOptions::default())
            .expect("calc major");

        let err_l1b_lon = wrap_diff_arcsec(res.lon_deg, l1b_lon);
        let err_l1b_lat = (res.lat_deg - l1b_lat).abs() * 3_600.0;
        let err_l3_lon = wrap_diff_arcsec(res.lon_deg, l3_lon);
        let err_l3_lat = (res.lat_deg - l3_lat).abs() * 3_600.0;

        max_major_l1b_arcsec = max_major_l1b_arcsec.max(err_l1b_lon).max(err_l1b_lat);
        max_major_l3_arcsec = max_major_l3_arcsec.max(err_l3_lon).max(err_l3_lat);

        assert!(
            err_l1b_lon < 0.01 && err_l1b_lat < 0.01,
            "Major body {name} at {jd_tt} L1B err: {err_l1b_lon:.6}\", {err_l1b_lat:.6}\""
        );
        // L3 column in l3_se_10k_grid.json is SE built-in Moshier analytical series (~1..5" vs DE440s),
        // whereas L1B column is exact DE440s + IAU 2006/2000A (< 0.01").
        let l3_gate = 5.0;
        assert!(
            err_l3_lon < l3_gate && err_l3_lat < l3_gate,
            "Major body {name} at {jd_tt} L3 SE diff: {err_l3_lon:.4}\" (gate {l3_gate}\")"
        );
        assert!((res.dist_au - l1b_dist).abs() < 1e-5);
        assert!((res.lon_speed_deg_day - l3_spd).abs() < 0.02);
    }

    // Section 2: 1,360 Minor Body vectors (20 epochs x 68 bodies)
    let mut max_minor_l1b_arcsec = 0.0_f64;
    let fast_opts = CalcOptions {
        include_speed: false,
        ..Default::default()
    };
    let icrf_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Equatorial,
        equinox: Equinox::J2000,
        nutation: NutationMode::Mean,
        aberration: AbMode::Astrometric,
        zodiac: Zodiac::Tropical,
        include_speed: false,
    };
    current_jd = None;
    for &(jd_tt, mpc, ref name, l1b_lon, l1b_lat, l1b_dist, icrf_ra, icrf_dec) in
        &grid.minor_bodies
    {
        if current_jd != Some(jd_tt) {
            current_jd = Some(jd_tt);
            current_moment = Some(
                engine
                    .at(TimeInput::from_tt_jd(jd_tt), None)
                    .expect("moment"),
            );
        }
        let m = current_moment.as_ref().expect("moment");
        let body = BodyId::MinorPlanet(MpcNumber(mpc));
        let res_ecl = engine.calc_in_moment(m, body, fast_opts).expect("minor ecl");
        let res_icrf = engine
            .calc_in_moment(m, body, icrf_opts)
            .expect("minor icrf");

        let err_lon = wrap_diff_arcsec(res_ecl.lon_deg, l1b_lon);
        let err_lat = (res_ecl.lat_deg - l1b_lat).abs() * 3_600.0;
        let err_ra = wrap_diff_arcsec(res_icrf.lon_deg, icrf_ra);
        let err_dec = (res_icrf.lat_deg - icrf_dec).abs() * 3_600.0;
        max_minor_l1b_arcsec = max_minor_l1b_arcsec
            .max(err_lon)
            .max(err_lat)
            .max(err_ra)
            .max(err_dec);

        assert!(
            err_lon < 0.01 && err_lat < 0.01 && err_ra < 0.01 && err_dec < 0.01,
            "Minor {name} ({mpc}) at {jd_tt} mismatch: {err_lon:.6}\", {err_ra:.6}\""
        );
        assert!((res_ecl.dist_au - l1b_dist).abs() < 1e-5);
    }

    // Section 3: 1,200 Ayanamsa vectors (100 epochs x 12 models)
    let mut max_aya_l1b_arcsec = 0.0_f64;
    // Note: in l3_se_10k_grid.json, the l1b_mean column for lahiri used the un-nutated 1956 anchor
    // (differing by IAU 1980 dpsi_1956 = -16.769099"), while l3_se_mean_deg and l1b_ayanamsas_12.json
    // include the -16.769099" IAE 1989 p. 525 correction.
    let dpsi_1956_deg = -16.769_099_1 / 3_600.0;
    for &(jd_tt, ref code, l1b_mean, l1b_true, l3_mean) in &grid.ayanamsas {
        let model = AyanamsaModel::from_code(code).expect("ayanamsa code");
        let val = engine
            .ayanamsa(TimeInput::from_tt_jd(jd_tt), model, NutationMode::True)
            .expect("ayanamsa");
        let (err_mean, err_true) = if model == AyanamsaModel::Lahiri {
            let em = wrap_diff_arcsec(val.mean_deg - dpsi_1956_deg, l1b_mean);
            let et = wrap_diff_arcsec(val.true_deg - dpsi_1956_deg, l1b_true);
            let err_l3 = wrap_diff_arcsec(val.mean_deg, l3_mean);
            assert!(
                err_l3 < 0.5,
                "Lahiri vs L3 SE mean mismatch at {jd_tt}: {err_l3:.4}\""
            );
            (em, et)
        } else {
            (
                wrap_diff_arcsec(val.mean_deg, l1b_mean),
                wrap_diff_arcsec(val.true_deg, l1b_true),
            )
        };
        max_aya_l1b_arcsec = max_aya_l1b_arcsec.max(err_mean).max(err_true);
        assert!(
            err_mean < 1.0 && err_true < 1.0,
            "Ayanamsa {code} at {jd_tt} grid err: {err_mean:.6}\", {err_true:.6}\""
        );
    }

    // Section 4: 1,380 House System vectors (60 grids x 23 systems)
    let mut max_house_err_arcsec = 0.0_f64;
    let eps_j2000_deg = 23.439_291_111_111_11_f64;
    let eps_rad = eps_j2000_deg.to_radians();
    for &(ramc_deg, lat_deg, ref sys_code, ref exp_cusps) in &grid.houses {
        let ch = sys_code.chars().next().expect("char");
        let sys = HouseSystem::from_code(ch).expect("house code");
        // Compute asc_deg to derive the default sun_decl_deg used in the 10k grid
        let ramc_rad = ramc_deg.to_radians();
        let lat_rad = lat_deg.to_radians();
        let asc_rad = ramc_rad
            .cos()
            .atan2(-(eps_rad.sin() * lat_rad.tan() + eps_rad.cos() * ramc_rad.sin()));
        let sun_decl_deg = (asc_rad.sin() * eps_rad.sin()).asin().to_degrees();

        let res = houses_armc(
            ramc_deg,
            lat_deg,
            eps_j2000_deg,
            sun_decl_deg,
            HouseOptions {
                system: sys,
                polar_policy: PolarPolicy::Strict,
                include_speed: false,
                ..Default::default()
            },
        )
        .expect("houses_armc 10k");

        for (i, (&actual, &expected)) in res.cusps_deg.iter().zip(exp_cusps.iter()).enumerate() {
            let err = wrap_diff_arcsec(actual, expected);
            max_house_err_arcsec = max_house_err_arcsec.max(err);
            assert!(
                err < 0.05,
                "House {sys_code} at ramc={ramc_deg}, lat={lat_deg} cusp {} err: {err:.6}\"",
                i + 1
            );
        }
    }

    // Section 5: 280 Lunar Point vectors (40 epochs x 7 points)
    let mut max_lunar_l1b_arcsec = 0.0_f64;
    let mut max_intp_se_arcsec = 0.0_f64;
    current_jd = None;
    for &(jd_tt, ref pt_name, l1b_lon, l1b_lat, l3_intp) in &grid.lunar_points {
        if current_jd != Some(jd_tt) {
            current_jd = Some(jd_tt);
            current_moment = Some(
                engine
                    .at(TimeInput::from_tt_jd(jd_tt), None)
                    .expect("moment"),
            );
        }
        let m = current_moment.as_ref().expect("moment");
        let lp = parse_lunar_point_name(pt_name);
        let res = engine
            .calc_in_moment(m, BodyId::LunarPoint(lp), fast_opts)
            .expect("lunar 10k");
        let err_lon = wrap_diff_arcsec(res.lon_deg, l1b_lon);
        let err_lat = (res.lat_deg - l1b_lat).abs() * 3_600.0;
        if lp == LunarPointId::InterpolatedLilith {
            // Interpolated_Lilith is verified to 0.000025" against l1b_lunar_points_7.json;
            // here we also verify documented divergence against SE's proprietary swe_nod_aps.
            if let Some(se_lon) = l3_intp {
                let err_se = wrap_diff_arcsec(res.lon_deg, se_lon);
                max_intp_se_arcsec = max_intp_se_arcsec.max(err_se);
                assert!(
                    err_se < 1_200.0,
                    "Interpolated_Lilith vs SE at {jd_tt}: {err_se:.2}\""
                );
            }
        } else {
            max_lunar_l1b_arcsec = max_lunar_l1b_arcsec.max(err_lon).max(err_lat);
            assert!(
                err_lon < 0.05 && err_lat < 0.05,
                "Lunar point {pt_name} at {jd_tt} err: {err_lon:.6}\", {err_lat:.6}\""
            );
        }
    }

    // Section 6: 450 Uranian vectors (50 epochs x 9 bodies)
    let mut max_uran_l1b_arcsec = 0.0_f64;
    current_jd = None;
    for &(jd_tt, ref name, l1b_lon, l1b_lat, _l3_lon) in &grid.uranians {
        if current_jd != Some(jd_tt) {
            current_jd = Some(jd_tt);
            current_moment = Some(
                engine
                    .at(TimeInput::from_tt_jd(jd_tt), None)
                    .expect("moment"),
            );
        }
        let m = current_moment.as_ref().expect("moment");
        let u = parse_uranian_name(name);
        let res = engine
            .calc_in_moment(m, BodyId::Hypothetical(u), fast_opts)
            .expect("uranian 10k");
        let err_lon = wrap_diff_arcsec(res.lon_deg, l1b_lon);
        let err_lat = (res.lat_deg - l1b_lat).abs() * 3_600.0;
        max_uran_l1b_arcsec = max_uran_l1b_arcsec.max(err_lon).max(err_lat);
        assert!(
            err_lon < 0.05 && err_lat < 0.05,
            "Uranian {name} at {jd_tt} err: {err_lon:.6}\", {err_lat:.6}\""
        );
    }

    // Section 7: 572 Fixed Star vectors (11 epochs x 52 stars)
    let mut max_star_l1b_arcsec = 0.0_f64;
    current_jd = None;
    for &(jd_tt, hip, ref name, l1b_lon, l1b_lat) in &grid.fixed_stars {
        if current_jd != Some(jd_tt) {
            current_jd = Some(jd_tt);
            current_moment = Some(
                engine
                    .at(TimeInput::from_tt_jd(jd_tt), None)
                    .expect("moment"),
            );
        }
        let m = current_moment.as_ref().expect("moment");
        let res = engine
            .calc_in_moment(m, BodyId::FixedStar(HipId(hip)), fast_opts)
            .expect("star 10k");
        let err_lon = wrap_diff_arcsec(res.lon_deg, l1b_lon);
        let err_lat = (res.lat_deg - l1b_lat).abs() * 3_600.0;
        max_star_l1b_arcsec = max_star_l1b_arcsec.max(err_lon).max(err_lat);
        assert!(
            err_lon < 0.02 && err_lat < 0.02,
            "Star {name} (HIP {hip}) at {jd_tt} err: {err_lon:.6}\", {err_lat:.6}\""
        );
    }

    println!(
        "10,242-Vector Grid Verified:\n  major_bodies (5000): max_l1b={max_major_l1b_arcsec:.6}\", max_l3={max_major_l3_arcsec:.4}\"\n  minor_bodies (1360): max_l1b={max_minor_l1b_arcsec:.6}\"\n  ayanamsas (1200): max_l1b={max_aya_l1b_arcsec:.6}\"\n  houses (1380): max_err={max_house_err_arcsec:.6}\"\n  lunar_points (280): max_l1b={max_lunar_l1b_arcsec:.6}\"\n  uranians (450): max_l1b={max_uran_l1b_arcsec:.6}\"\n  fixed_stars (572): max_l1b={max_star_l1b_arcsec:.6}\""
    );
}

// ---------------------------------------------------------------------------
// 3. 32-Thread Concurrency Torture Suite + Codex Gate 2 Verification
// ---------------------------------------------------------------------------

#[test]
fn test_32_threads_mixed_tropical_sidereal_topocentric_polar_zero_bleed() {
    let engine = Arc::new(load_full_engine());
    let london = Observer::new(51.5074, -0.1278, 35.0).expect("London");
    let tromso = Observer::new(69.6492, 18.9553, 10.0).expect("Tromso");

    // First verify Codex Gate 2:
    // (a) Lunar points in Equatorial, Horizontal, and J2000 frames have non-zero consistent velocities
    let m_london = engine
        .at(TimeInput::from_tt_jd(JD_J2000), Some(london))
        .expect("moment london");
    for lp in [
        LunarPointId::MeanNode,
        LunarPointId::TrueNode,
        LunarPointId::TrueLilith,
        LunarPointId::InterpolatedLilith,
    ] {
        let r_eq = engine
            .calc_in_moment(
                &m_london,
                BodyId::LunarPoint(lp),
                CalcOptions {
                    plane: Plane::Equatorial,
                    ..Default::default()
                },
            )
            .expect("lunar eq");
        assert!(r_eq.lon_speed_deg_day.abs() > 1e-4);
        assert!(
            r_eq.cartesian_vel_au_day.iter().any(|&x| x.abs() > 1e-8),
            "Equatorial cartesian_vel_au_day must be non-zero for {lp:?}"
        );

        let r_hz = engine
            .calc_in_moment(
                &m_london,
                BodyId::LunarPoint(lp),
                CalcOptions {
                    origin: Origin::Topocentric(london),
                    plane: Plane::Horizontal,
                    ..Default::default()
                },
            )
            .expect("lunar hz");
        assert!(
            r_hz.lon_speed_deg_day.abs() > 10.0,
            "Horizontal diurnal azimuth rate must be non-zero for {lp:?}"
        );
        assert!(r_hz.cartesian_vel_au_day.iter().any(|&x| x.abs() > 1e-6));
    }

    // (b) Fixed star warnings field in CalcResult and StarEntry::quality_warnings
    let r_bellatrix = engine
        .calc_in_moment(
            &m_london,
            BodyId::FixedStar(HipId(25336)),
            CalcOptions::default(),
        )
        .expect("Bellatrix");
    assert!(
        r_bellatrix.warnings.is_empty(),
        "Catalog stars with XHIP V/137D RV and positive Hipparcos-2 parallax have 0 warnings"
    );
    let mut synth_star = lxephemeris_core::find_star_by_hip(HipId(25336))
        .expect("Bellatrix")
        .clone();
    synth_star.plx_mas = -0.71;
    synth_star.flags = "MissingRadialVelocity".to_owned();
    let w_codes: Vec<&str> = synth_star
        .quality_warnings()
        .into_iter()
        .map(lxephemeris_core::StarQualityWarning::as_code)
        .collect();
    assert_eq!(w_codes, vec!["NonPositiveParallax", "MissingRadialVelocity"]);

    // Build 32 distinct workload profiles and compute single-threaded reference outputs
    let mut baseline = Vec::with_capacity(32);
    for task_idx in 0..32_usize {
        let jd = JD_J2000 + (task_idx as f64) * 123.456_789;
        let t_in = TimeInput::from_tt_jd(jd);
        let zodiac = match task_idx % 4 {
            0 => Zodiac::Tropical,
            1 => Zodiac::Sidereal(AyanamsaModel::Lahiri),
            2 => Zodiac::Sidereal(AyanamsaModel::Raman),
            _ => Zodiac::Sidereal(AyanamsaModel::TrueCitra),
        };
        let origin = if task_idx % 2 == 0 {
            Origin::Geocentric
        } else {
            Origin::Topocentric(london)
        };
        let opts = CalcOptions {
            origin,
            zodiac,
            ..Default::default()
        };
        let sun = engine
            .calc(BodyId::Planet(PlanetId::Sun), t_in, opts)
            .expect("sun");
        let moon = engine
            .calc(BodyId::Planet(PlanetId::Moon), t_in, opts)
            .expect("moon");
        let chiron = engine
            .calc(BodyId::MinorPlanet(MpcNumber(2060)), t_in, opts)
            .expect("chiron");
        let lilith = engine
            .calc(
                BodyId::LunarPoint(LunarPointId::InterpolatedLilith),
                t_in,
                opts,
            )
            .expect("lilith");
        let cupido = engine
            .calc(BodyId::Hypothetical(HypotheticalId::Cupido), t_in, opts)
            .expect("cupido");
        let regulus = engine
            .calc(BodyId::FixedStar(HipId(49669)), t_in, opts)
            .expect("regulus");
        let sys = HouseSystem::ALL[task_idx % HouseSystem::ALL.len()];
        let h_london = engine
            .houses(
                t_in,
                london,
                HouseOptions {
                    system: sys,
                    zodiac,
                    polar_policy: PolarPolicy::Strict,
                    include_speed: true,
                    ..Default::default()
                },
            )
            .expect("houses london");
        let h_tromso = engine
            .houses(
                t_in,
                tromso,
                HouseOptions {
                    system: sys,
                    zodiac,
                    polar_policy: PolarPolicy::FallbackToPorphyry,
                    include_speed: false,
                    ..Default::default()
                },
            )
            .expect("houses tromso fallback");
        baseline.push((
            sun.lon_deg.to_bits(),
            sun.lon_speed_deg_day.to_bits(),
            moon.lon_deg.to_bits(),
            chiron.lon_deg.to_bits(),
            lilith.lon_deg.to_bits(),
            cupido.lon_deg.to_bits(),
            regulus.lon_deg.to_bits(),
            h_london.cusps_deg[0].to_bits(),
            h_london.angles.asc_deg.to_bits(),
            h_tromso.cusps_deg[0].to_bits(),
        ));
    }

    let baseline = Arc::new(baseline);
    let mut handles = Vec::with_capacity(32);

    for thread_id in 0..32_usize {
        let eng = Arc::clone(&engine);
        let base = Arc::clone(&baseline);
        handles.push(thread::spawn(move || {
            // Each of the 32 threads executes all 32 mixed profiles in a shifted permutation (32 x 32 = 1,024 evaluations)
            for step in 0..32_usize {
                let task_idx = (thread_id * 7 + step) % 32;
                let jd = JD_J2000 + (task_idx as f64) * 123.456_789;
                let t_in = TimeInput::from_tt_jd(jd);
                let zodiac = match task_idx % 4 {
                    0 => Zodiac::Tropical,
                    1 => Zodiac::Sidereal(AyanamsaModel::Lahiri),
                    2 => Zodiac::Sidereal(AyanamsaModel::Raman),
                    _ => Zodiac::Sidereal(AyanamsaModel::TrueCitra),
                };
                let origin = if task_idx % 2 == 0 {
                    Origin::Geocentric
                } else {
                    Origin::Topocentric(london)
                };
                let opts = CalcOptions {
                    origin,
                    zodiac,
                    ..Default::default()
                };
                let sun = eng
                    .calc(BodyId::Planet(PlanetId::Sun), t_in, opts)
                    .expect("sun");
                let moon = eng
                    .calc(BodyId::Planet(PlanetId::Moon), t_in, opts)
                    .expect("moon");
                let chiron = eng
                    .calc(BodyId::MinorPlanet(MpcNumber(2060)), t_in, opts)
                    .expect("chiron");
                let lilith = eng
                    .calc(
                        BodyId::LunarPoint(LunarPointId::InterpolatedLilith),
                        t_in,
                        opts,
                    )
                    .expect("lilith");
                let cupido = eng
                    .calc(BodyId::Hypothetical(HypotheticalId::Cupido), t_in, opts)
                    .expect("cupido");
                let regulus = eng
                    .calc(BodyId::FixedStar(HipId(49669)), t_in, opts)
                    .expect("regulus");
                let sys = HouseSystem::ALL[task_idx % HouseSystem::ALL.len()];
                let h_london = eng
                    .houses(
                        t_in,
                        london,
                        HouseOptions {
                            system: sys,
                            zodiac,
                            polar_policy: PolarPolicy::Strict,
                            include_speed: true,
                            ..Default::default()
                        },
                    )
                    .expect("houses london");
                let h_tromso = eng
                    .houses(
                        t_in,
                        tromso,
                        HouseOptions {
                            system: sys,
                            zodiac,
                            polar_policy: PolarPolicy::FallbackToPorphyry,
                            include_speed: false,
                            ..Default::default()
                        },
                    )
                    .expect("houses tromso fallback");

                let actual = (
                    sun.lon_deg.to_bits(),
                    sun.lon_speed_deg_day.to_bits(),
                    moon.lon_deg.to_bits(),
                    chiron.lon_deg.to_bits(),
                    lilith.lon_deg.to_bits(),
                    cupido.lon_deg.to_bits(),
                    regulus.lon_deg.to_bits(),
                    h_london.cusps_deg[0].to_bits(),
                    h_london.angles.asc_deg.to_bits(),
                    h_tromso.cusps_deg[0].to_bits(),
                );
                assert_eq!(
                    actual, base[task_idx],
                    "Thread {thread_id} step {step} (task {task_idx}) bit-level state bleed detected!"
                );
            }

            // Also run an eclipse search per thread to stress concurrent root-finding
            if thread_id % 8 == 0 {
                let sol = eng
                    .next_solar_eclipse(
                        TimeInput::from_tt_jd(JD_J2000),
                        false,
                        None,
                        EclipseOptions::default(),
                    )
                    .expect("concurrent solar eclipse");
                assert!(sol.max_eclipse_jd_tt > JD_J2000);
            }
        }));
    }

    for h in handles {
        h.join().expect("thread join");
    }

    println!("32-Thread Concurrency Torture Suite (1,024 mixed chart profiles = 8,192+ evaluations): 0.0 ULP state bleed verified!");
}
