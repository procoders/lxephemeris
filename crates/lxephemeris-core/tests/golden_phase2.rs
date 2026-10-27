//! Phase 2 + Task 3.1 Golden Verification Suite:
//! - Level 1B 12 Vedic Ayanamsas (`testdata/golden/l1b_ayanamsas_12.json`, gate `< 0.005"`)
//! - Level 1B 7 Lunar Points (`testdata/golden/l1b_lunar_points_7.json`, gate `< 0.05"`)
//! - Level 1B 9 Uranian / Hypothetical Bodies (`testdata/golden/l1b_uranians_9.json`, gate `< 0.05"`)
//! - Level 1B 52 Fixed Stars (`testdata/golden/l1b_stars_50.json`, gate `< 0.01"`)
//! - All 68 Minor Bodies (`MINOR_BODIES_68`) via `Moment` + `safe_interval` enforcement
//! - Sidereal Zodiac (`Zodiac::Sidereal`) nutation invariance and `manifest.json` size validation

#![allow(clippy::indexing_slicing, clippy::panic)]

use lxephemeris_core::{
    find_star_by_hip, star_catalog_52, AbMode, AyanamsaModel, BodyId, CalcOptions, Engine,
    EngineConfig, EphemerisError, Equinox, HipId, HypotheticalId, LunarPointId, MpcNumber,
    NutationMode, Origin, Plane, PlanetId, StarQualityWarning, TimeInput, Zodiac,
    JD_1850_01_01, JD_2100_01_01, JD_INTERP_LILITH_SAFE_MAX, JD_INTERP_LILITH_SAFE_MIN,
    JD_J2000, MINOR_BODIES_68,
};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf()
}

fn wrap_diff_arcsec(a_deg: f64, b_deg: f64) -> f64 {
    let diff = (a_deg - b_deg + 180.0).rem_euclid(360.0) - 180.0;
    diff.abs() * 3_600.0
}

fn load_full_engine() -> Engine {
    let root = workspace_root();
    Engine::from_data_dir(root.join("data"), EngineConfig::default())
        .unwrap_or_else(|e| panic!("Failed to load Engine from data/: {e}"))
}

// ---------------------------------------------------------------------------
// 1. Level 1B: 12 Vedic Ayanamsas (120 vectors, gate < 0.005")
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
struct L1bAyanamsaFile {
    tolerance_arcsec: f64,
    count: usize,
    vectors: Vec<L1bAyanamsaVector>,
}

#[derive(Debug, Deserialize)]
struct L1bAyanamsaVector {
    epoch_label: String,
    jd_tt: f64,
    ayanamsa_code: String,
    ayanamsa_mean_deg: f64,
    ayanamsa_true_deg: f64,
    rate_deg_day: f64,
    dpsi_arcsec: f64,
}

#[test]
fn test_l1b_ayanamsas_12_golden_120() {
    let root = workspace_root();
    let json_text = fs::read_to_string(root.join("testdata/golden/l1b_ayanamsas_12.json"))
        .unwrap_or_else(|e| panic!("read l1b_ayanamsas_12.json: {e}"));
    let golden: L1bAyanamsaFile = serde_json::from_str(&json_text)
        .unwrap_or_else(|e| panic!("parse l1b_ayanamsas_12.json: {e}"));
    assert_eq!(golden.count, 120);
    assert_eq!(golden.vectors.len(), 120);

    let engine = load_full_engine();
    let mut max_mean_err_arcsec = 0.0_f64;
    let mut max_true_err_arcsec = 0.0_f64;
    let mut max_rate_err_deg_day = 0.0_f64;

    for v in &golden.vectors {
        let model = AyanamsaModel::from_code(&v.ayanamsa_code)
            .unwrap_or_else(|e| panic!("unknown ayanamsa {}: {e}", v.ayanamsa_code));
        let val = engine
            .ayanamsa(TimeInput::from_tt_jd(v.jd_tt), model, NutationMode::True)
            .unwrap_or_else(|e| panic!("ayanamsa failed at {}: {e}", v.epoch_label));

        let err_mean = wrap_diff_arcsec(val.mean_deg, v.ayanamsa_mean_deg);
        let err_true = wrap_diff_arcsec(val.true_deg, v.ayanamsa_true_deg);
        let err_dpsi = (val.dpsi_arcsec - v.dpsi_arcsec).abs();
        let err_rate = (val.rate_deg_day - v.rate_deg_day).abs();

        max_mean_err_arcsec = max_mean_err_arcsec.max(err_mean);
        max_true_err_arcsec = max_true_err_arcsec.max(err_true);
        max_rate_err_deg_day = max_rate_err_deg_day.max(err_rate);

        assert!(
            err_mean < golden.tolerance_arcsec,
            "L1B Ayanamsa mean mismatch for {} at {}: {:.6}\" (gate {:.4}\")",
            v.ayanamsa_code,
            v.epoch_label,
            err_mean,
            golden.tolerance_arcsec
        );
        assert!(
            err_true < golden.tolerance_arcsec,
            "L1B Ayanamsa true mismatch for {} at {}: {:.6}\" (gate {:.4}\")",
            v.ayanamsa_code,
            v.epoch_label,
            err_true,
            golden.tolerance_arcsec
        );
        assert!(err_dpsi < 1e-6, "dpsi mismatch: {err_dpsi}");
        assert!(err_rate < 1e-9, "rate mismatch: {err_rate}");
    }

    println!(
        "L1B 12 Ayanamsas (120 vectors): max_mean_err = {max_mean_err_arcsec:.6}\", max_true_err = {max_true_err_arcsec:.6}\", max_rate_err = {max_rate_err_deg_day:.3e} deg/day"
    );
}

#[test]
fn test_sidereal_zodiac_nutation_invariance() {
    let engine = load_full_engine();
    let moment = engine
        .at(TimeInput::from_tt_jd(JD_J2000), None)
        .unwrap_or_else(|e| panic!("moment: {e}"));

    for body in [
        BodyId::Planet(PlanetId::Sun),
        BodyId::Planet(PlanetId::Moon),
        BodyId::MinorPlanet(MpcNumber(1)),
        BodyId::LunarPoint(LunarPointId::TrueNode),
        BodyId::Hypothetical(HypotheticalId::Cupido),
    ] {
        let opts_true = CalcOptions {
            zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
            nutation: NutationMode::True,
            ..Default::default()
        };
        let opts_mean = CalcOptions {
            nutation: NutationMode::Mean,
            ..opts_true
        };

        let r_true = engine
            .calc_in_moment(&moment, body, opts_true)
            .unwrap_or_else(|e| panic!("calc true: {e}"));
        let r_mean = engine
            .calc_in_moment(&moment, body, opts_mean)
            .unwrap_or_else(|e| panic!("calc mean: {e}"));

        // First-order nutation in longitude (~17") cancels out in sidereal ecliptic longitude;
        // residual is only second-order nutation in obliquity (~deps * tan(lat) < 0.8" even at 5 deg lat, < 0.001" on ecliptic)
        let diff_arcsec = wrap_diff_arcsec(r_true.lon_deg, r_mean.lon_deg);
        let limit_arcsec = if r_true.lat_deg.abs() < 0.01 {
            0.01
        } else {
            1.0
        };
        assert!(
            diff_arcsec < limit_arcsec,
            "Sidereal nutation invariance failed for {body:?}: {diff_arcsec:.6}\" (limit {limit_arcsec}\")"
        );
        assert!(
            (r_true.lon_speed_deg_day - r_mean.lon_speed_deg_day).abs() < 1e-6,
            "Sidereal speed nutation invariance failed for {body:?}: {} vs {}",
            r_true.lon_speed_deg_day,
            r_mean.lon_speed_deg_day
        );
    }
}

// ---------------------------------------------------------------------------
// 2. Level 1B: 7 Lunar Points (63 vectors, gate < 0.05")
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
struct L1bLunarFile {
    tolerance_arcsec: f64,
    count: usize,
    vectors: Vec<L1bLunarEpoch>,
}

#[derive(Debug, Deserialize)]
struct L1bLunarEpoch {
    epoch_label: String,
    jd_tt: f64,
    points: BTreeMap<String, L1bLunarPointExpected>,
}

#[derive(Debug, Deserialize)]
struct L1bLunarPointExpected {
    lon_deg: f64,
    lon_mean_deg: Option<f64>,
    lat_deg: f64,
    lon_speed_deg_day: f64,
    lat_speed_deg_day: Option<f64>,
}

fn parse_lunar_point_key(key: &str) -> LunarPointId {
    match key {
        "Mean_Node" => LunarPointId::MeanNode,
        "Mean_South_Node" => LunarPointId::MeanSouthNode,
        "True_Node" => LunarPointId::TrueNode,
        "True_South_Node" => LunarPointId::TrueSouthNode,
        "Mean_Lilith" => LunarPointId::MeanLilith,
        "True_Lilith" => LunarPointId::TrueLilith,
        "Interpolated_Lilith" => LunarPointId::InterpolatedLilith,
        other => panic!("Unknown lunar point key: {other}"),
    }
}

#[test]
fn test_l1b_lunar_points_7_golden() {
    let root = workspace_root();
    let json_text = fs::read_to_string(root.join("testdata/golden/l1b_lunar_points_7.json"))
        .unwrap_or_else(|e| panic!("read l1b_lunar_points_7.json: {e}"));
    let golden: L1bLunarFile = serde_json::from_str(&json_text)
        .unwrap_or_else(|e| panic!("parse l1b_lunar_points_7.json: {e}"));
    assert_eq!(golden.count, 9);

    let engine = load_full_engine();
    let mut max_lon_err_arcsec = 0.0_f64;
    let mut max_lat_err_arcsec = 0.0_f64;
    let mut max_speed_err_deg_day = 0.0_f64;

    for ep in &golden.vectors {
        let moment = engine
            .at(TimeInput::from_tt_jd(ep.jd_tt), None)
            .unwrap_or_else(|e| panic!("moment at {}: {e}", ep.epoch_label));

        for (key, exp) in &ep.points {
            let lp = parse_lunar_point_key(key);
            let res = engine
                .calc_in_moment(&moment, BodyId::LunarPoint(lp), CalcOptions::default())
                .unwrap_or_else(|e| panic!("lunar point {key} at {}: {e}", ep.epoch_label));

            let err_lon = wrap_diff_arcsec(res.lon_deg, exp.lon_deg);
            let err_lat = (res.lat_deg - exp.lat_deg).abs() * 3_600.0;
            let err_spd = (res.lon_speed_deg_day - exp.lon_speed_deg_day).abs();

            max_lon_err_arcsec = max_lon_err_arcsec.max(err_lon);
            max_lat_err_arcsec = max_lat_err_arcsec.max(err_lat);
            max_speed_err_deg_day = max_speed_err_deg_day.max(err_spd);

            assert!(
                err_lon < golden.tolerance_arcsec,
                "Lunar point {key} lon mismatch at {}: {err_lon:.6}\" (gate {:.4}\")",
                ep.epoch_label,
                golden.tolerance_arcsec
            );
            assert!(
                err_lat < golden.tolerance_arcsec,
                "Lunar point {key} lat mismatch at {}: {err_lat:.6}\" (gate {:.4}\")",
                ep.epoch_label,
                golden.tolerance_arcsec
            );
            assert!(
                err_spd < 1e-6,
                "Lunar point {key} speed mismatch at {}: {err_spd:.3e}",
                ep.epoch_label
            );

            if let Some(exp_lat_spd) = exp.lat_speed_deg_day {
                assert!(
                    (res.lat_speed_deg_day - exp_lat_spd).abs() < 1e-6,
                    "Lunar point {key} lat speed mismatch at {}",
                    ep.epoch_label
                );
            }

            if let Some(exp_mean_lon) = exp.lon_mean_deg {
                let mean_opts = CalcOptions {
                    nutation: NutationMode::Mean,
                    ..Default::default()
                };
                let res_mean = engine
                    .calc_in_moment(&moment, BodyId::LunarPoint(lp), mean_opts)
                    .unwrap_or_else(|e| panic!("mean lunar point: {e}"));
                let err_mean_lon = wrap_diff_arcsec(res_mean.lon_deg, exp_mean_lon);
                assert!(
                    err_mean_lon < golden.tolerance_arcsec,
                    "Lunar point {key} mean lon mismatch at {}: {err_mean_lon:.6}\"",
                    ep.epoch_label
                );
            }
        }
    }

    // Verify InterpolatedLilith safe_interval [1850-03-01, 2149-11-01)
    let ok_start = engine.calc(
        BodyId::LunarPoint(LunarPointId::InterpolatedLilith),
        TimeInput::from_tt_jd(JD_INTERP_LILITH_SAFE_MIN),
        CalcOptions::default(),
    );
    assert!(ok_start.is_ok(), "InterpolatedLilith at safe min must succeed");

    let err_before = engine.calc(
        BodyId::LunarPoint(LunarPointId::InterpolatedLilith),
        TimeInput::from_tt_jd(JD_INTERP_LILITH_SAFE_MIN - 0.1),
        CalcOptions::default(),
    );
    assert!(matches!(err_before, Err(EphemerisError::OutOfCoverage { .. })));

    let err_at_end = engine.calc(
        BodyId::LunarPoint(LunarPointId::InterpolatedLilith),
        TimeInput::from_tt_jd(JD_INTERP_LILITH_SAFE_MAX),
        CalcOptions::default(),
    );
    assert!(matches!(err_at_end, Err(EphemerisError::OutOfCoverage { .. })));

    println!(
        "L1B 7 Lunar Points (63 vectors): max_lon_err = {max_lon_err_arcsec:.6}\", max_lat_err = {max_lat_err_arcsec:.6}\", max_speed_err = {max_speed_err_deg_day:.3e} deg/day"
    );
}

// ---------------------------------------------------------------------------
// 3. Level 1B: 9 Uranians / Hypotheticals (54 vectors, gate < 0.05")
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
struct L1bUranianFile {
    tolerance_arcsec: f64,
    count: usize,
    vectors: Vec<L1bUranianVector>,
}

#[derive(Debug, Deserialize)]
struct L1bUranianVector {
    epoch_label: String,
    jd_tt: f64,
    target_name: String,
    apparent_ecliptic_true_of_date: L1bUranianState,
}

#[derive(Debug, Deserialize)]
struct L1bUranianState {
    lon_deg: f64,
    lat_deg: f64,
    dist_au: f64,
    lon_speed_deg_day: f64,
    lat_speed_deg_day: f64,
    dist_speed_au_day: f64,
}

fn parse_hypothetical_name(name: &str) -> HypotheticalId {
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
        other => panic!("Unknown hypothetical name: {other}"),
    }
}

#[test]
fn test_l1b_uranians_9_golden_54() {
    let root = workspace_root();
    let json_text = fs::read_to_string(root.join("testdata/golden/l1b_uranians_9.json"))
        .unwrap_or_else(|e| panic!("read l1b_uranians_9.json: {e}"));
    let golden: L1bUranianFile = serde_json::from_str(&json_text)
        .unwrap_or_else(|e| panic!("parse l1b_uranians_9.json: {e}"));
    assert_eq!(golden.count, 54);

    let engine = load_full_engine();
    let mut max_lon_err_arcsec = 0.0_f64;
    let mut max_lat_err_arcsec = 0.0_f64;
    let mut max_spd_err_deg_day = 0.0_f64;

    for v in &golden.vectors {
        let id = parse_hypothetical_name(&v.target_name);
        let res = engine
            .calc(
                BodyId::Hypothetical(id),
                TimeInput::from_tt_jd(v.jd_tt),
                CalcOptions::default(),
            )
            .unwrap_or_else(|e| panic!("uranian {} at {}: {e}", v.target_name, v.epoch_label));

        let exp = &v.apparent_ecliptic_true_of_date;
        let err_lon = wrap_diff_arcsec(res.lon_deg, exp.lon_deg);
        let err_lat = (res.lat_deg - exp.lat_deg).abs() * 3_600.0;
        let err_dist = (res.dist_au - exp.dist_au).abs();
        let err_spd = (res.lon_speed_deg_day - exp.lon_speed_deg_day).abs();

        max_lon_err_arcsec = max_lon_err_arcsec.max(err_lon);
        max_lat_err_arcsec = max_lat_err_arcsec.max(err_lat);
        max_spd_err_deg_day = max_spd_err_deg_day.max(err_spd);

        assert!(
            err_lon < golden.tolerance_arcsec,
            "Uranian {} lon mismatch at {}: {err_lon:.6}\" (gate {:.4}\")",
            v.target_name,
            v.epoch_label,
            golden.tolerance_arcsec
        );
        assert!(
            err_lat < golden.tolerance_arcsec,
            "Uranian {} lat mismatch at {}: {err_lat:.6}\" (gate {:.4}\")",
            v.target_name,
            v.epoch_label,
            golden.tolerance_arcsec
        );
        assert!(err_dist < 1e-7, "Uranian {} dist mismatch: {err_dist}", v.target_name);
        assert!(
            err_spd < 1e-5,
            "Uranian {} speed mismatch at {}: {err_spd:.3e}",
            v.target_name,
            v.epoch_label
        );
        assert!((res.lat_speed_deg_day - exp.lat_speed_deg_day).abs() < 1e-5);
        assert!((res.dist_speed_au_day - exp.dist_speed_au_day).abs() < 1e-5);
    }

    println!(
        "L1B 9 Uranians (54 vectors): max_lon_err = {max_lon_err_arcsec:.6}\", max_lat_err = {max_lat_err_arcsec:.6}\", max_spd_err = {max_spd_err_deg_day:.3e} deg/day"
    );
}

// ---------------------------------------------------------------------------
// 4. Level 1B: 52 Fixed Stars (104 vectors, gate < 0.01")
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
struct L1bStarsFile {
    tolerance_arcsec: f64,
    count: usize,
    vectors: Vec<L1bStarVector>,
}

#[derive(Debug, Deserialize)]
struct L1bStarVector {
    epoch_label: String,
    jd_tt: f64,
    name: String,
    bayer: String,
    hip: u32,
    icrs_j2000_0: L1bStarEqCoords,
    apparent_of_date: L1bStarApparentCoords,
}

#[derive(Debug, Deserialize)]
struct L1bStarEqCoords {
    ra_deg: f64,
    dec_deg: f64,
}

#[derive(Debug, Deserialize)]
struct L1bStarApparentCoords {
    lon_deg: f64,
    lat_deg: f64,
    ra_deg: f64,
    dec_deg: f64,
}

#[test]
fn test_l1b_stars_50_golden_104() {
    let root = workspace_root();
    let json_text = fs::read_to_string(root.join("testdata/golden/l1b_stars_50.json"))
        .unwrap_or_else(|e| panic!("read l1b_stars_50.json: {e}"));
    let golden: L1bStarsFile = serde_json::from_str(&json_text)
        .unwrap_or_else(|e| panic!("parse l1b_stars_50.json: {e}"));
    assert_eq!(golden.count, 104);
    assert_eq!(star_catalog_52().len(), 52);

    let engine = load_full_engine();
    let mut max_j2000_err_arcsec = 0.0_f64;
    let mut max_ecl_err_arcsec = 0.0_f64;
    let mut max_eq_err_arcsec = 0.0_f64;

    for v in &golden.vectors {
        let hip = HipId(v.hip);
        let moment = engine
            .at(TimeInput::from_tt_jd(v.jd_tt), None)
            .unwrap_or_else(|e| panic!("moment at {}: {e}", v.epoch_label));

        // 1. ICRS J2000.0 (Geometric, Equatorial, J2000)
        let j2000_opts = CalcOptions {
            origin: Origin::Barycentric,
            plane: Plane::Equatorial,
            equinox: Equinox::J2000,
            nutation: NutationMode::Mean,
            aberration: AbMode::Geometric,
            zodiac: Zodiac::Tropical,
            include_speed: false,
        };
        let res_j2000 = engine
            .calc_in_moment(&moment, BodyId::FixedStar(hip), j2000_opts)
            .unwrap_or_else(|e| panic!("star j2000 {}: {e}", v.name));
        let err_ra_2000 = wrap_diff_arcsec(res_j2000.lon_deg, v.icrs_j2000_0.ra_deg);
        let err_dec_2000 = (res_j2000.lat_deg - v.icrs_j2000_0.dec_deg).abs() * 3_600.0;
        max_j2000_err_arcsec = max_j2000_err_arcsec.max(err_ra_2000).max(err_dec_2000);
        assert!(
            err_ra_2000 < golden.tolerance_arcsec && err_dec_2000 < golden.tolerance_arcsec,
            "Star {} J2000 mismatch: dRA={err_ra_2000:.6}\", dDec={err_dec_2000:.6}\"",
            v.name
        );

        // 2. Apparent Ecliptic of Date
        let res_ecl = engine
            .calc_in_moment(&moment, BodyId::FixedStar(hip), CalcOptions::default())
            .unwrap_or_else(|e| panic!("star ecl {}: {e}", v.name));
        let err_lon = wrap_diff_arcsec(res_ecl.lon_deg, v.apparent_of_date.lon_deg);
        let err_lat = (res_ecl.lat_deg - v.apparent_of_date.lat_deg).abs() * 3_600.0;
        max_ecl_err_arcsec = max_ecl_err_arcsec.max(err_lon).max(err_lat);
        assert!(
            err_lon < golden.tolerance_arcsec && err_lat < golden.tolerance_arcsec,
            "Star {} apparent ecliptic mismatch at {}: dLon={err_lon:.6}\", dLat={err_lat:.6}\"",
            v.name,
            v.epoch_label
        );

        // 3. Apparent Equatorial of Date
        let eq_opts = CalcOptions {
            plane: Plane::Equatorial,
            ..Default::default()
        };
        let res_eq = engine
            .calc_in_moment(&moment, BodyId::FixedStar(hip), eq_opts)
            .unwrap_or_else(|e| panic!("star eq {}: {e}", v.name));
        let err_ra = wrap_diff_arcsec(res_eq.lon_deg, v.apparent_of_date.ra_deg);
        let err_dec = (res_eq.lat_deg - v.apparent_of_date.dec_deg).abs() * 3_600.0;
        max_eq_err_arcsec = max_eq_err_arcsec.max(err_ra).max(err_dec);
        assert!(
            err_ra < golden.tolerance_arcsec && err_dec < golden.tolerance_arcsec,
            "Star {} apparent equatorial mismatch at {}: dRA={err_ra:.6}\", dDec={err_dec:.6}\"",
            v.name,
            v.epoch_label
        );

        // 4. Lookup by canonical name, Bayer designation, and "hip:<id>"
        let fast_opts = CalcOptions {
            include_speed: false,
            ..Default::default()
        };
        let by_name = engine
            .calc_star_by_name(&v.name, TimeInput::from_tt_jd(v.jd_tt), fast_opts)
            .unwrap_or_else(|e| panic!("by_name {}: {e}", v.name));
        let by_bayer = engine
            .calc_star_by_name(&v.bayer, TimeInput::from_tt_jd(v.jd_tt), fast_opts)
            .unwrap_or_else(|e| panic!("by_bayer {}: {e}", v.bayer));
        let by_hip_str = engine
            .calc_star_by_name(
                &format!("hip:{}", v.hip),
                TimeInput::from_tt_jd(v.jd_tt),
                fast_opts,
            )
            .unwrap_or_else(|e| panic!("by_hip_str {}: {e}", v.hip));
        assert_eq!(by_name.lon_deg, res_ecl.lon_deg);
        assert_eq!(by_bayer.lon_deg, res_ecl.lon_deg);
        assert_eq!(by_hip_str.lon_deg, res_ecl.lon_deg);
    }

    // Verify quality warning and clamping for non-positive parallax (e.g., Hipparcos-1 Mintaka plx_mas = -0.71)
    let mut synthetic_nonpos = find_star_by_hip(HipId(21421))
        .expect("Aldebaran in catalog")
        .clone();
    synthetic_nonpos.plx_mas = -0.71;
    assert_eq!(
        synthetic_nonpos.quality_warnings(),
        vec![StarQualityWarning::NonPositiveParallax]
    );
    assert_eq!(synthetic_nonpos.nominal_distance_au(), 1.0e9);

    println!(
        "L1B 52 Fixed Stars (104 vectors): max_j2000_err = {max_j2000_err_arcsec:.6}\", max_ecl_err = {max_ecl_err_arcsec:.6}\", max_eq_err = {max_eq_err_arcsec:.6}\""
    );
}

// ---------------------------------------------------------------------------
// 5. Task 2.1: All 68 Minor Bodies + Safe Interval [1850-01-01, 2100-01-01)
// ---------------------------------------------------------------------------

#[test]
fn test_all_68_minor_bodies_in_moment_and_safe_interval() {
    let engine = load_full_engine();
    assert_eq!(engine.capabilities().minor_bodies_loaded_count, 68);

    // Evaluate all 68 minor bodies at J2000.0 and at the exact inclusive start 1850-01-01T00:00:00 TT
    for jd_tt in [JD_1850_01_01, JD_J2000, JD_2100_01_01 - 1.0] {
        let moment = engine
            .at(TimeInput::from_tt_jd(jd_tt), None)
            .unwrap_or_else(|e| panic!("moment at {jd_tt}: {e}"));
        for spec in &MINOR_BODIES_68 {
            let res = engine
                .calc_in_moment(&moment, BodyId::MinorPlanet(spec.mpc), CalcOptions::default())
                .unwrap_or_else(|e| panic!("minor body ({}) {} at {jd_tt}: {e}", spec.mpc.0, spec.name));
            assert!(res.lon_deg >= 0.0 && res.lon_deg < 360.0);
            assert!(res.lat_deg.abs() <= 90.0);
            assert!(res.dist_au > 0.0);
            assert!(res.lon_speed_deg_day.is_finite());
        }
    }

    // Half-open bound enforcement: JD_1850_01_01 - 0.01 and JD_2100_01_01 must return OutOfCoverage
    let ceres = BodyId::MinorPlanet(MpcNumber(1));
    let before_start = engine.calc(
        ceres,
        TimeInput::from_tt_jd(JD_1850_01_01 - 0.01),
        CalcOptions::default(),
    );
    assert!(matches!(before_start, Err(EphemerisError::OutOfCoverage { .. })));

    let at_end = engine.calc(
        ceres,
        TimeInput::from_tt_jd(JD_2100_01_01),
        CalcOptions::default(),
    );
    assert!(matches!(at_end, Err(EphemerisError::OutOfCoverage { .. })));
}
