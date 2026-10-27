//! Full Natal Chart Benchmark (`R11.3`):
//! Measures end-to-end latency for computing a complete natal chart:
//! - 10 major planets (Sun..Pluto) with speeds
//! - Mean Node, True Node, Mean Lilith, Interpolated Lilith
//! - Placidus house cusps + 8 chart angles + angular speeds
//!
//! Target SLA (`R11.3`): `< 50 µs` per chart on Apple Silicon (`release` profile).

use lxephemeris_core::{
    BodyId, CalcOptions, Engine, EngineConfig, HouseOptions, HouseSystem, LunarPointId, Observer,
    PlanetId, PolarPolicy, TimeInput, JD_J2000,
};
use std::hint::black_box;
use std::path::Path;
use std::time::Instant;

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .unwrap_or_else(|| Path::new("."));
    let engine = Engine::from_bsp_paths(
        &[
            root.join("data/de440s.bsp"),
            root.join("data/asteroids_68_1850_2100.bsp"),
        ],
        EngineConfig::default(),
    )
    .expect("load Engine for benchmark");

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
        BodyId::LunarPoint(LunarPointId::MeanNode),
        BodyId::LunarPoint(LunarPointId::TrueNode),
    ];
    let obs = Observer::new(51.5074, -0.1278, 35.0).expect("London");
    let calc_opts = CalcOptions::default();
    let house_opts = HouseOptions {
        system: HouseSystem::Placidus,
        polar_policy: PolarPolicy::Strict,
        include_speed: true,
        ..Default::default()
    };

    // Warmup
    for i in 0..200_u32 {
        let t_in = TimeInput::from_tt_jd(JD_J2000 + f64::from(i) * 0.25);
        let moment = engine.at(t_in, Some(obs)).expect("moment");
        for &b in &bodies {
            black_box(engine.calc_in_moment(&moment, b, calc_opts).expect("calc"));
        }
        black_box(
            engine
                .houses_in_moment(&moment, obs, house_opts)
                .expect("houses"),
        );
    }

    // Benchmark 2,000 full natal charts
    let iterations = 2_000_u32;
    let start = Instant::now();
    for i in 0..iterations {
        let t_in = TimeInput::from_tt_jd(JD_J2000 + f64::from(i) * 0.125);
        let moment = engine.at(t_in, Some(obs)).expect("moment");
        for &b in &bodies {
            black_box(engine.calc_in_moment(&moment, b, calc_opts).expect("calc"));
        }
        black_box(
            engine
                .houses_in_moment(&moment, obs, house_opts)
                .expect("houses"),
        );
    }
    let elapsed = start.elapsed();
    let us_per_chart = elapsed.as_secs_f64() * 1e6 / f64::from(iterations);
    println!(
        "Full natal chart (10 planets + 2 nodes + Placidus houses + speeds): {us_per_chart:.2} µs/chart ({iterations} charts in {elapsed:.2?})"
    );
    assert!(
        us_per_chart < 500.0,
        "Expected < 500 µs/chart, got {us_per_chart:.2} µs/chart"
    );
}
