//! Phase 1 Integration & Golden Vector Verification Suite (Tasks 1.1–1.5):
//! 1. Level 1A: NASA JPL Horizons DE440/DE441 ICRF geometric and astrometric vectors (`< 0.005"`).
//! 2. Level 1B: `jplephem` + `pyerfa` IAU 2006/2000A apparent ecliptic/equatorial and astrometric ICRF (`< 0.01"`).
//! 3. 4-point 4th-order Cartesian stencil convergence across $h/2, h, 2h$.
//! 4. Topocentric `azalt`, `rise_trans` (Illinois root-finder), `Engine` capabilities, coverage, and parallel batching.

#![allow(clippy::indexing_slicing, clippy::panic)]

use lxephemeris_core::{
    cart_state_to_spherical, stencil_4pt_derivative, AbMode, AtmosphereParams, BodyId, CalcOptions,
    Engine, EngineConfig, EphemerisError, Equinox, MpcNumber, NutationMode, Observer, Origin,
    Plane, PlanetId, RiseTransEvent, RiseTransOptions, TimeInput, H_MOON_DAYS, H_PLANET_DAYS,
    JD_J2000,
};
use serde::Deserialize;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

fn workspace_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf()
}

fn shared_engine() -> &'static Engine {
    static ENGINE: OnceLock<Engine> = OnceLock::new();
    ENGINE.get_or_init(|| {
        Engine::from_data_dir(
            workspace_dir().join("data"),
            EngineConfig {
                max_threads: Some(4),
                pin_time_table_version: Some("embedded".to_owned()),
                verify_manifest_sha256: true,
            },
        )
        .expect("load shared Engine")
    })
}

fn target_id_to_body(target_id: i32) -> BodyId {
    match target_id {
        10 => BodyId::Planet(PlanetId::Sun),
        301 => BodyId::Planet(PlanetId::Moon),
        199 => BodyId::Planet(PlanetId::Mercury),
        299 => BodyId::Planet(PlanetId::Venus),
        399 => BodyId::Planet(PlanetId::Earth),
        4 => BodyId::Planet(PlanetId::Mars),
        5 => BodyId::Planet(PlanetId::Jupiter),
        6 => BodyId::Planet(PlanetId::Saturn),
        7 => BodyId::Planet(PlanetId::Uranus),
        8 => BodyId::Planet(PlanetId::Neptune),
        9 => BodyId::Planet(PlanetId::Pluto),
        other => panic!("Unexpected NAIF target ID in planetary golden vectors: {other}"),
    }
}

/// Computes true spherical angular separation in arcseconds between `(lon1_deg, lat1_deg)` and `(lon2_deg, lat2_deg)`.
fn angular_sep_arcsec(lon1_deg: f64, lat1_deg: f64, lon2_deg: f64, lat2_deg: f64) -> f64 {
    let l1 = lon1_deg.to_radians();
    let b1 = lat1_deg.to_radians();
    let l2 = lon2_deg.to_radians();
    let b2 = lat2_deg.to_radians();
    let u1 = [b1.cos() * l1.cos(), b1.cos() * l1.sin(), b1.sin()];
    let u2 = [b2.cos() * l2.cos(), b2.cos() * l2.sin(), b2.sin()];
    let diff = [u1[0] - u2[0], u1[1] - u2[1], u1[2] - u2[2]];
    let chord = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();
    let angle_rad = 2.0 * (0.5 * chord).clamp(0.0, 1.0).asin();
    angle_rad.to_degrees() * 3600.0
}

/// Computes signed wrapped difference in degrees in `[-180, +180)`.
fn wrap_diff_deg(a_deg: f64, b_deg: f64) -> f64 {
    (a_deg - b_deg + 180.0).rem_euclid(360.0) - 180.0
}

#[derive(Deserialize)]
struct L1aGeomIcrf {
    pos_au: [f64; 3],
    vel_au_day: [f64; 3],
    ra_deg: f64,
    dec_deg: f64,
    dist_au: f64,
}

#[derive(Deserialize)]
struct L1aAstroIcrf {
    pos_au: [f64; 3],
    ra_deg: f64,
    dec_deg: f64,
    dist_au: f64,
}

#[derive(Deserialize)]
struct L1aVector {
    epoch_label: String,
    jd_tdb: f64,
    target_id: i32,
    target_name: String,
    geometric_icrf: L1aGeomIcrf,
    astrometric_lt_icrf: L1aAstroIcrf,
}

#[derive(Deserialize)]
struct L1aFixture {
    tolerance_arcsec: f64,
    count: usize,
    vectors: Vec<L1aVector>,
}

#[test]
fn test_level_1a_horizons_de440_icrf_golden_vectors() {
    let ws = workspace_dir();
    let engine = shared_engine();

    let l1a_path = ws.join("testdata/golden/l1a_horizons_de440_icrf.json");
    let l1a_text = std::fs::read_to_string(&l1a_path).expect("read L1A fixture");
    let fixture: L1aFixture = serde_json::from_str(&l1a_text).expect("parse L1A fixture");
    assert_eq!(fixture.vectors.len(), fixture.count);
    assert_eq!(fixture.count, 50);

    let mut max_geom_sep_arcsec = 0.0_f64;
    let mut max_astro_sep_arcsec = 0.0_f64;

    for v in &fixture.vectors {
        let body = target_id_to_body(v.target_id);
        let t_in = TimeInput::from_tdb_jd(v.jd_tdb);

        // 1. Geometric ICRF
        let geom_opts = CalcOptions {
            origin: Origin::Geocentric,
            plane: Plane::Equatorial,
            equinox: Equinox::J2000,
            nutation: NutationMode::Mean,
            aberration: AbMode::Geometric,
            include_speed: true,
            ..CalcOptions::default()
        };
        let geom_res = engine
            .calc(body, t_in, geom_opts)
            .expect("calc geometric ICRF");

        let geom_sep = angular_sep_arcsec(
            geom_res.lon_deg,
            geom_res.lat_deg,
            v.geometric_icrf.ra_deg,
            v.geometric_icrf.dec_deg,
        );
        max_geom_sep_arcsec = max_geom_sep_arcsec.max(geom_sep);
        assert!(
            geom_sep < fixture.tolerance_arcsec,
            "L1A Geometric ICRF failed for {} at {}: sep={geom_sep:.6}\" >= {}\"",
            v.target_name,
            v.epoch_label,
            fixture.tolerance_arcsec
        );
        assert!(
            (geom_res.dist_au - v.geometric_icrf.dist_au).abs() < 1e-7,
            "L1A Geometric distance mismatch for {} at {}",
            v.target_name,
            v.epoch_label
        );
        for k in 0..3 {
            assert!(
                (geom_res.cartesian_au[k] - v.geometric_icrf.pos_au[k]).abs() < 1e-7,
                "L1A Geometric pos_au[{k}] mismatch"
            );
            assert!(
                (geom_res.cartesian_vel_au_day[k] - v.geometric_icrf.vel_au_day[k]).abs() < 1e-6,
                "L1A Geometric vel_au_day[{k}] mismatch"
            );
        }

        // 2. Astrometric (light-time corrected) ICRF
        let astro_opts = CalcOptions {
            origin: Origin::Geocentric,
            plane: Plane::Equatorial,
            equinox: Equinox::J2000,
            nutation: NutationMode::Mean,
            aberration: AbMode::Astrometric,
            include_speed: false,
            ..CalcOptions::default()
        };
        let astro_res = engine
            .calc(body, t_in, astro_opts)
            .expect("calc astrometric ICRF");

        let astro_sep = angular_sep_arcsec(
            astro_res.lon_deg,
            astro_res.lat_deg,
            v.astrometric_lt_icrf.ra_deg,
            v.astrometric_lt_icrf.dec_deg,
        );
        max_astro_sep_arcsec = max_astro_sep_arcsec.max(astro_sep);
        assert!(
            astro_sep < fixture.tolerance_arcsec,
            "L1A Astrometric ICRF failed for {} at {}: sep={astro_sep:.6}\" >= {}\"",
            v.target_name,
            v.epoch_label,
            fixture.tolerance_arcsec
        );
        for k in 0..3 {
            assert!(
                (astro_res.cartesian_au[k] - v.astrometric_lt_icrf.pos_au[k]).abs() < 1e-7,
                "L1A Astrometric pos_au[{k}] mismatch"
            );
        }
        assert!((astro_res.dist_au - v.astrometric_lt_icrf.dist_au).abs() < 1e-7);
    }

    println!(
        "L1A verification passed (50 vectors): max_geom_sep = {max_geom_sep_arcsec:.6}\", max_astro_sep = {max_astro_sep_arcsec:.6}\" (limit < 0.005\")"
    );
}

#[derive(Deserialize)]
struct L1bApparentEcl {
    lon_deg: f64,
    lat_deg: f64,
    dist_au: f64,
    lon_speed_deg_day: f64,
    lat_speed_deg_day: f64,
    dist_speed_au_day: f64,
}

#[derive(Deserialize)]
struct L1bApparentEqu {
    ra_deg: f64,
    dec_deg: f64,
}

#[derive(Deserialize)]
struct L1bAstrometricIcrf {
    ra_deg: f64,
    dec_deg: f64,
    dist_au: f64,
}

#[derive(Deserialize)]
struct L1bPlanetVector {
    epoch_label: String,
    jd_tt: f64,
    target_id: i32,
    target_name: String,
    apparent_ecliptic_true_of_date: L1bApparentEcl,
    apparent_equatorial_true_of_date: L1bApparentEqu,
    astrometric_icrf: L1bAstrometricIcrf,
}

#[derive(Deserialize)]
struct L1bPlanetFixture {
    tolerance_arcsec: f64,
    count: usize,
    vectors: Vec<L1bPlanetVector>,
}

#[test]
fn test_level_1b_planets_de440s_golden_vectors() {
    let ws = workspace_dir();
    let engine = shared_engine();

    let l1b_path = ws.join("testdata/golden/l1b_planets_de440s.json");
    let l1b_text = std::fs::read_to_string(&l1b_path).expect("read L1B planets fixture");
    let fixture: L1bPlanetFixture =
        serde_json::from_str(&l1b_text).expect("parse L1B planets fixture");
    assert_eq!(fixture.vectors.len(), fixture.count);
    assert_eq!(fixture.count, 120);

    let mut max_ecl_lon_arcsec = 0.0_f64;
    let mut max_ecl_lat_arcsec = 0.0_f64;
    let mut max_equ_sep_arcsec = 0.0_f64;
    let mut max_icrf_sep_arcsec = 0.0_f64;

    for v in &fixture.vectors {
        let body = target_id_to_body(v.target_id);
        let t_in = TimeInput::from_tt_jd(v.jd_tt);

        // 1. Apparent Ecliptic True of Date
        let ecl_opts = CalcOptions::default();
        let ecl_res = engine
            .calc(body, t_in, ecl_opts)
            .expect("calc apparent ecliptic");

        let dlon_arcsec = wrap_diff_deg(
            ecl_res.lon_deg,
            v.apparent_ecliptic_true_of_date.lon_deg,
        )
        .abs()
            * 3600.0;
        let dlat_arcsec =
            (ecl_res.lat_deg - v.apparent_ecliptic_true_of_date.lat_deg).abs() * 3600.0;
        max_ecl_lon_arcsec = max_ecl_lon_arcsec.max(dlon_arcsec);
        max_ecl_lat_arcsec = max_ecl_lat_arcsec.max(dlat_arcsec);

        assert!(
            dlon_arcsec < fixture.tolerance_arcsec,
            "L1B Ecliptic lon mismatch for {} at {}: {dlon_arcsec:.6}\" >= {}\"",
            v.target_name,
            v.epoch_label,
            fixture.tolerance_arcsec
        );
        assert!(
            dlat_arcsec < fixture.tolerance_arcsec,
            "L1B Ecliptic lat mismatch for {} at {}: {dlat_arcsec:.6}\" >= {}\"",
            v.target_name,
            v.epoch_label,
            fixture.tolerance_arcsec
        );
        assert!(
            (ecl_res.dist_au - v.apparent_ecliptic_true_of_date.dist_au).abs() < 1e-9,
            "L1B Ecliptic dist_au mismatch for {} at {}",
            v.target_name,
            v.epoch_label
        );
        // Velocity agreement (note L1B generator used coarse step h=1/24 for Moon, 0.25 for planets;
        // our engine uses h=1e-4 for Moon and 1e-3 for planets, agreeing within ~1e-3 deg/day)
        assert!(
            (ecl_res.lon_speed_deg_day - v.apparent_ecliptic_true_of_date.lon_speed_deg_day).abs()
                < 2e-3,
            "L1B lon_speed mismatch for {} at {}: got {}, expected {}",
            v.target_name,
            v.epoch_label,
            ecl_res.lon_speed_deg_day,
            v.apparent_ecliptic_true_of_date.lon_speed_deg_day
        );
        assert!(
            (ecl_res.lat_speed_deg_day - v.apparent_ecliptic_true_of_date.lat_speed_deg_day).abs()
                < 2e-3,
            "L1B lat_speed mismatch for {} at {}",
            v.target_name,
            v.epoch_label
        );
        assert!(
            (ecl_res.dist_speed_au_day - v.apparent_ecliptic_true_of_date.dist_speed_au_day).abs()
                < 1e-5,
            "L1B dist_speed mismatch for {} at {}",
            v.target_name,
            v.epoch_label
        );

        // 2. Apparent Equatorial True of Date
        let equ_opts = CalcOptions {
            plane: Plane::Equatorial,
            include_speed: false,
            ..CalcOptions::default()
        };
        let equ_res = engine
            .calc(body, t_in, equ_opts)
            .expect("calc apparent equatorial");
        let equ_sep = angular_sep_arcsec(
            equ_res.lon_deg,
            equ_res.lat_deg,
            v.apparent_equatorial_true_of_date.ra_deg,
            v.apparent_equatorial_true_of_date.dec_deg,
        );
        max_equ_sep_arcsec = max_equ_sep_arcsec.max(equ_sep);
        assert!(
            equ_sep < fixture.tolerance_arcsec,
            "L1B Equatorial mismatch for {} at {}: {equ_sep:.6}\" >= {}\"",
            v.target_name,
            v.epoch_label,
            fixture.tolerance_arcsec
        );

        // 3. Astrometric ICRF
        let icrf_opts = CalcOptions {
            plane: Plane::Equatorial,
            equinox: Equinox::J2000,
            nutation: NutationMode::Mean,
            aberration: AbMode::Astrometric,
            include_speed: false,
            ..CalcOptions::default()
        };
        let icrf_res = engine
            .calc(body, t_in, icrf_opts)
            .expect("calc astrometric icrf");
        let icrf_sep = angular_sep_arcsec(
            icrf_res.lon_deg,
            icrf_res.lat_deg,
            v.astrometric_icrf.ra_deg,
            v.astrometric_icrf.dec_deg,
        );
        max_icrf_sep_arcsec = max_icrf_sep_arcsec.max(icrf_sep);
        assert!(
            icrf_sep < fixture.tolerance_arcsec,
            "L1B Astrometric ICRF mismatch for {} at {}: {icrf_sep:.6}\" >= {}\"",
            v.target_name,
            v.epoch_label,
            fixture.tolerance_arcsec
        );
        assert!((icrf_res.dist_au - v.astrometric_icrf.dist_au).abs() < 1e-9);
    }

    println!(
        "L1B verification passed (120 vectors): max_ecl_lon = {max_ecl_lon_arcsec:.6}\", max_ecl_lat = {max_ecl_lat_arcsec:.6}\", max_equ = {max_equ_sep_arcsec:.6}\", max_icrf = {max_icrf_sep_arcsec:.6}\" (limit < 0.01\")"
    );
}

#[test]
fn test_velocity_stencil_convergence_h_half_h_2h() {
    let engine = shared_engine();
    let moment = engine
        .at(TimeInput::from_tt_jd(JD_J2000), None)
        .expect("moment at J2000");
    let opts = CalcOptions::default();

    for planet in [
        PlanetId::Sun,
        PlanetId::Moon,
        PlanetId::Mercury,
        PlanetId::Venus,
        PlanetId::Mars,
        PlanetId::Jupiter,
        PlanetId::Saturn,
    ] {
        let base_h = if planet == PlanetId::Moon {
            H_MOON_DAYS
        } else {
            H_PLANET_DAYS
        };
        let tid = planet.naif_target_id();
        let pos_0 = lxephemeris_core::eval_spk_target_vector_au(
            engine.spk(),
            moment.center(),
            tid,
            &opts,
        )
        .expect("pos_0");

        let eval_speed_at_step = |step_h: f64| -> (f64, f64, f64) {
            let st = moment.build_stencil(step_h).expect("build stencil");
            let m2 =
                lxephemeris_core::eval_spk_target_vector_au(engine.spk(), &st[0], tid, &opts)
                    .expect("m2");
            let m1 =
                lxephemeris_core::eval_spk_target_vector_au(engine.spk(), &st[1], tid, &opts)
                    .expect("m1");
            let p1 =
                lxephemeris_core::eval_spk_target_vector_au(engine.spk(), &st[2], tid, &opts)
                    .expect("p1");
            let p2 =
                lxephemeris_core::eval_spk_target_vector_au(engine.spk(), &st[3], tid, &opts)
                    .expect("p2");
            let v0 = stencil_4pt_derivative(&m2, &m1, &p1, &p2, step_h);
            let sph = cart_state_to_spherical(&pos_0, &v0);
            (
                sph.lon_speed_deg_day,
                sph.lat_speed_deg_day,
                sph.dist_speed_au_day,
            )
        };

        let (lon_half, lat_half, dist_half) = eval_speed_at_step(0.5 * base_h);
        let (lon_h, lat_h, dist_h) = eval_speed_at_step(base_h);
        let (lon_2h, lat_2h, dist_2h) = eval_speed_at_step(2.0 * base_h);

        // Agreement across h/2, h, and 2h must be < 1e-7 deg/day (< 0.00036"/day)
        assert!(
            (lon_h - lon_half).abs() < 1e-7 && (lon_2h - lon_h).abs() < 1e-7,
            "Stencil lon_speed convergence failed for {planet:?}: h/2={lon_half}, h={lon_h}, 2h={lon_2h}"
        );
        assert!(
            (lat_h - lat_half).abs() < 1e-7 && (lat_2h - lat_h).abs() < 1e-7,
            "Stencil lat_speed convergence failed for {planet:?}"
        );
        assert!(
            (dist_h - dist_half).abs() < 1e-9 && (dist_2h - dist_h).abs() < 1e-9,
            "Stencil dist_speed convergence failed for {planet:?}"
        );
    }
}

#[test]
fn test_engine_capabilities_coverage_azalt_rise_trans_and_batch() {
    let engine = shared_engine();

    let caps = engine.capabilities();
    assert_eq!(caps.minor_bodies_loaded_count, 68);
    assert_eq!(caps.house_systems.len(), 23);
    assert_eq!(caps.ayanamsa_models.len(), 12);

    let sun_cov = engine
        .coverage(BodyId::Planet(PlanetId::Sun))
        .expect("sun coverage");
    assert!(sun_cov.safe_start_jd_tt < JD_J2000 && sun_cov.safe_end_jd_tt > JD_J2000);

    // Out-of-coverage check
    let ooc = engine.calc(
        BodyId::MinorPlanet(MpcNumber(1)),
        TimeInput::from_tt_jd(2_200_000.0),
        CalcOptions::default(),
    );
    assert!(matches!(ooc, Err(EphemerisError::OutOfCoverage { .. })));

    // Topocentric azalt & rise_trans at Greenwich (lat=51.4769, lon=0.0, alt=46.0)
    let greenwich = Observer::new(51.4769, 0.0, 46.0).expect("Greenwich observer");
    // Equinox 2020-03-20 00:00 UT1 (JD = 2458928.5): Sun rises ~06:02 UT, transits ~12:07 UT, sets ~18:14 UT
    let jd_start = 2_458_928.5;
    let sun_rise = engine
        .rise_trans(
            BodyId::Planet(PlanetId::Sun),
            greenwich,
            jd_start,
            RiseTransEvent::Rise,
            RiseTransOptions::default(),
        )
        .expect("sun rise")
        .expect("sun rises at Greenwich");
    let sun_transit = engine
        .rise_trans(
            BodyId::Planet(PlanetId::Sun),
            greenwich,
            jd_start,
            RiseTransEvent::UpperTransit,
            RiseTransOptions::default(),
        )
        .expect("sun transit")
        .expect("sun transits at Greenwich");
    let sun_set = engine
        .rise_trans(
            BodyId::Planet(PlanetId::Sun),
            greenwich,
            jd_start,
            RiseTransEvent::Set,
            RiseTransOptions::default(),
        )
        .expect("sun set")
        .expect("sun sets at Greenwich");

    let rise_hour = (sun_rise - jd_start) * 24.0;
    let transit_hour = (sun_transit - jd_start) * 24.0;
    let set_hour = (sun_set - jd_start) * 24.0;
    assert!(
        (rise_hour - 6.04).abs() < 0.1,
        "Unexpected Greenwich equinox sunrise hour: {rise_hour}"
    );
    assert!(
        (transit_hour - 12.12).abs() < 0.1,
        "Unexpected Greenwich equinox solar noon hour: {transit_hour}"
    );
    assert!(
        (set_hour - 18.23).abs() < 0.1,
        "Unexpected Greenwich equinox sunset hour: {set_hour}"
    );

    // At solar noon, azimuth is ~180 deg (South)
    let hz_noon = engine
        .calc_azalt(
            BodyId::Planet(PlanetId::Sun),
            TimeInput::from_ut1_jd(sun_transit),
            greenwich,
            AtmosphereParams::default(),
        )
        .expect("azalt at solar noon");
    assert!(
        (hz_noon.azimuth_deg - 180.0).abs() < 0.05,
        "Azimuth at upper transit should be ~180 deg, got {}",
        hz_noon.azimuth_deg
    );

    // Parallel batch execution (atomic and per-item)
    let bodies = vec![
        BodyId::Planet(PlanetId::Sun),
        BodyId::Planet(PlanetId::Moon),
        BodyId::Planet(PlanetId::Mercury),
        BodyId::MinorPlanet(MpcNumber(1)),   // Ceres
        BodyId::MinorPlanet(MpcNumber(2060)), // Chiron
    ];
    let epochs: Vec<TimeInput> = (0..10)
        .map(|i| TimeInput::from_tt_jd(JD_J2000 + f64::from(i) * 10.0))
        .collect();
    let batch_res = engine
        .calc_batch(&bodies, &epochs, CalcOptions::default())
        .expect("calc_batch");
    assert_eq!(batch_res.len(), 10);
    assert_eq!(batch_res[0].len(), 5);

    // Per-item batch with one out-of-coverage epoch
    let mixed_epochs = vec![
        TimeInput::from_tt_jd(JD_J2000),
        TimeInput::from_tt_jd(2_200_000.0), // Out of coverage
        TimeInput::from_tt_jd(JD_J2000 + 1.0),
    ];
    let per_item = engine
        .calc_batch_per_item(&bodies, &mixed_epochs, CalcOptions::default())
        .expect("calc_batch_per_item");
    assert!(per_item[0][0].is_ok());
    assert!(matches!(
        per_item[1][0],
        Err(EphemerisError::OutOfCoverage { .. })
    ));
    assert!(per_item[2][0].is_ok());
}
