# lxephemeris

**The modern astronomical & astrological computation engine for Python, Node.js, and Rust.**  
*100% memory-safe Rust core (`#![forbid(unsafe_code)]`). Lock-free multi-threading. Zero external files required to start — with native NASA JPL `DE440s` support for `< 0.003"` precision.*

[![CI](https://github.com/procoders/lxephemeris/actions/workflows/ci.yml/badge.svg)](https://github.com/procoders/lxephemeris/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/lxephemeris.svg)](https://pypi.org/project/lxephemeris/)
[![npm](https://img.shields.io/npm/v/lxephemeris.svg)](https://www.npmjs.com/package/lxephemeris)
[![Crates.io](https://img.shields.io/crates/v/lxephemeris-core.svg)](https://crates.io/crates/lxephemeris-core)
[![License](https://img.shields.io/badge/License-PolyForm_Noncommercial_%7C_Commercial-purple.svg)](#licensing--commercial-use)

```python
import lxephemeris

# Works out of the box with zero external files (or loads NASA JPL DE440s if present)
engine = lxephemeris.Engine()

jd_utc = lxephemeris.utc_from_local(2026, 10, 27, 14, 30, 0.0, tz_name="America/New_York")
obs = lxephemeris.Observer(latitude_deg=40.7128, longitude_deg=-74.0060, altitude_m=10.0)

# Full natal snapshot: planets, 23 house systems, aspects, Hermetic lots, lunar phase & sect in ~150 µs
chart = engine.chart_snapshot(jd_utc, obs)
```

---

## Why `lxephemeris`?

Legacy backends stitch together 30-year-old C bindings (`pyswisseph`), `PyEphem`, and `Kerykeion` — inheriting process-global mutable C state, thread race conditions, silent polar house errors, and slow FFI loops. `lxephemeris` replaces the entire stack with a single clean-room Rust engine:

- **One Engine Replaces 6+ Libraries**: Ships 10 planets, **68 JPL minor bodies** (`Chiron`, `Ceres`, `Eris`, `Sedna`, `Lilith (1181)`...), 7 Lunar nodes/apogees, 9 Uranians, **23 house systems** (with explicit polar policies), **43+1 sidereal ayanamsas**, **52 Hipparcos-2/XHIP fixed stars**, solar/lunar **eclipses & occultations**, **Vedic Panchanga / 16 Vargas / Vimshottari Dasha**, **Traditional Sect / Planetary Hours / 7 Hermetic Lots**, **Human Design 64 Gates & `-88°` arc**, **24 Chinese Solar Terms (*Jieqi*)**, and native **Transit / Return / Station / VoC Moon** solvers.
- **Lock-Free Thread Safety & Memory Safety**: `#![forbid(unsafe_code)]`, `Send + Sync`, zero global C state, zero locks required across concurrent Tropical, Sidereal, and Topocentric requests.
- **5–15× Faster**: Single-pass analytical velocity Jacobians (`12.8 µs`/body with speed; `~150 µs`/full natal chart) + GIL-free Rayon multi-core batching.
- **22 Upstream Bugs Eliminated**: Fixes global C races, `±66°` latitude clamps (**7.85° Svalbard ASC error**), `Mira`/`Mirach` star collisions, star-anchor aberration wobble, and `Sun().hlong` $180^\circ$ inversions.
- **Drop-In `SweAdapter`**: Migrate existing `pyswisseph` codebases in minutes without rewriting call sites.

👉 **[Read the full architectural deep-dive & 22 eliminated legacy bugs →](docs/WHY_LXEPHEMERIS.md)**

### ☁️ Want a Turnkey Cloud API Instead of Self-Hosting?

`lxephemeris` is the low-level astronomical & astrological math engine that powers **[AstroAPI (`astrology-api.io`)](https://astrology-api.io/)** and **[Astrolium (`astrolium.com`)](https://astrolium.com/)**.

If you are building a commercial astrology app, SaaS, or mobile product and don't want to maintain your own ephemeris infrastructure — **or purchase a self-hosted commercial engine license** — use **[AstroAPI](https://astrology-api.io/)** directly:
- **Turnkey REST API** powered by `lxephemeris` in production
- **Built-in City Geocoding & Historical Timezone / DST Resolution** (`birth city + date/time` → full chart)
- **SVG Chart Wheel Rendering**, **Synastry Compatibility Scoring**, and **Structured Interpretations**
- **Commercial Rights Included**: Your **[AstroAPI subscription](https://astrology-api.io/)** includes commercial production usage — no separate `lxephemeris` commercial license required.

---

## How `lxephemeris` Compares

| Feature | **`lxephemeris`** | **Swiss Ephemeris (`pyswisseph`)** | **Kerykeion** | **PyEphem (`ephem`)** | **Skyfield** |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Core Architecture** | **100% Safe Rust (`#![forbid(unsafe_code)]`)** | Legacy C (`libswe`, 1990s) | Python wrapper over `pyswisseph` | Legacy C (`libastro`) | Pure Python / NumPy |
| **Thread Safety** | **Lock-free `Send + Sync` (0 state bleed)** | ❌ Global mutable C `static` state | ❌ Inherits `pyswisseph` globals | ⚠️ Mutable `Body` objects | ✅ Thread-safe |
| **Works with 0 External Files** | **✅ Built-in analytical engine (`-3000..+3000`)** | ⚠️ Silent fallback if `.se1` missing | ✅ Bundles `.se1` files | ✅ Analytical only | ❌ Requires downloading `.bsp` |
| **JPL `DE440s` / `DE441` Precision** | **✅ Native DAF/SPK (`< 0.003"`)** | Compressed `.se1` re-fits | Via `pyswisseph` | ❌ Truncated VSOP87 (multi-arcsec drift) | ✅ Native SPK |
| **Single Body + Speed Latency** | **`12.8 µs` (analytical Jacobian)** | `~45–90 µs` (3-pt finite diff) | `~25–80 ms` per subject | `~80–150 µs` | `~600–2,000 µs` |
| **23 Houses & Arctic/Polar Safety** | **✅ All 23 + explicit `PolarPolicy`** | ⚠️ Silent polar fallback & `house_pos=0.0` bug | ❌ Clamps `\|lat\| > 66°` (**7.85° ASC error**) | ❌ No astrological houses | ❌ No astrological houses |
| **68 Astrological Minor Bodies** | **✅ Single verified `.bsp` + Keplerian fallback** | Manage 68 separate `.se1` files | Partial | ❌ Manual orbital elements | Manual SPK segments |
| **43+1 Sidereal Ayanamsas** | **✅ Aberration-free star anchors** | ⚠️ `±20.5"` annual aberration wobble on star anchors | ⚠️ Leaves secondary points tropical | ❌ None | ❌ None |
| **Transit, Return & VoC Moon Solvers** | **✅ Native Rust Newton-Illinois solvers** | ❌ Write your own FFI loops | ❌ Static charts only | ❌ Basic rising/phase only | Discrete `almanac` only |
| **Vedic (Panchanga, Dasha, D1–D60)** | **✅ Built-in** | ❌ None | ❌ None | ❌ None | ❌ None |
| **Traditional (Sect, Hours, 7 Lots)** | **✅ Built-in** | ❌ None | ❌ None | ❌ None | ❌ None |
| **Human Design (64 Gates & `-88°` Arc)** | **✅ Built-in** | ❌ None | ❌ None | ❌ None | ❌ None |
| **SDKs with Identical Math** | **Python, Node.js/TS, Rust, CLI** | C + separate 3rd-party bindings | Python only | Python only | Python only |

---

## Quickstart

### Installation

```bash
# Python (PyPI)
pip install lxephemeris

# Node.js / TypeScript (npm)
npm install lxephemeris

# Rust Library & CLI (crates.io)
cargo add lxephemeris-core
cargo install lxephemeris-cli
```

### Python (`lxephemeris`)

```python
import lxephemeris

# 1. Start immediately (zero external files needed) or load NASA JPL DE440s from ./data
engine = lxephemeris.Engine()
# engine = lxephemeris.Engine.from_data_dir("./data", verify_manifest_sha256=True)

# 2. Local wall-clock time -> UTC Julian Date (preserves fractional seconds & PEP 495 DST fold)
jd_utc = lxephemeris.utc_from_local(2026, 10, 24, 14, 0, 32.5, tz_name="Europe/Berlin", fold=0)
obs = lxephemeris.Observer(51.5074, -0.1278, 35.0)

# 3. Single body position + analytical speed + visual magnitude
sun = engine.calc("Sun", jd_utc, scale="UTC")
chiron = engine.calc("Chiron", jd_utc, scale="UTC")
print(f"Sun: {sun.lon_deg:.6f}° (speed={sun.lon_speed_deg_day:.6f}°/d, Vmag={sun.magnitude_v:.2f})")

# 4. Polar-safe houses & cusp boundary guard (eps = 1e-9°)
houses = engine.houses(
    jd_utc,
    latitude_deg=51.5074,
    longitude_deg=-0.1278,
    system="P",
    scale="UTC",
    polar_policy="fallback_to_whole_sign",
)
sun_house = lxephemeris.house_of_lon(sun.lon_deg, houses.cusps_deg, eps_deg=1e-9)

# 5. All-in-one Natal Snapshot, Solar Return, Vedic Panchanga, & Human Design Gate
snapshot = engine.chart_snapshot(jd_utc, obs)
solar_return = engine.planetary_return("Sun", jd_utc, jd_utc + 350.0)
panchanga = engine.panchanga(jd_utc, obs, ayanamsa="lahiri")
hd = lxephemeris.map_zodiac_placement(sun.lon_deg)
print(f"ASC: {houses.asc_deg:.4f}° | Sun in House {sun_house} | HD Gate {hd.hd_gate}.{hd.hd_line}")
```

### Node.js / TypeScript (`lxephemeris`)

```typescript
import { Engine, utcFromLocal, houseOfLon, mapZodiacPlacement } from "lxephemeris";

const engine = new Engine();
// const engine = Engine.fromDataDir("./data", 4, null, true);

const jdUtc = utcFromLocal(2026, 10, 24, 14, 0, 32.5, "Europe/Berlin");
const jdSplit = { jd1: jdUtc, jd2: 0.0 };
const obs = { latitudeDeg: 51.5074, longitudeDeg: -0.1278, altitudeM: 35 };

const sun = engine.calc("Sun", jdSplit, "UTC");
const houses = engine.houses(jdSplit, obs, "P", "UTC", "fallback_to_whole_sign");
const sunHouse = houseOfLon(sun.lonDeg, houses.cuspsDeg, 1e-9);
const snapshot = engine.chartSnapshot(jdUtc, obs);
const hd = mapZodiacPlacement(sun.lonDeg);

console.log(
  `Sun: ${sun.lonDeg.toFixed(6)}° in House ${sunHouse} | HD Gate ${hd.hdGate}.${hd.hdLine} | Lots: ${snapshot.lots.length}`
);
```

### Rust (`lxephemeris-core`)

```rust
use lxephemeris_core::{
    AspectPreset, BodyId, CalcOptions, Engine, HouseOptions, HouseSystem, Observer,
    PlanetId, PolarPolicy, TimeInput, JD_J2000,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let engine = Engine::default();
    let t = TimeInput::from_tt_jd(JD_J2000);
    let sun = engine.calc(BodyId::Planet(PlanetId::Sun), t, CalcOptions::default())?;

    let obs = Observer::new(51.5074, -0.1278, 35.0)?;
    let houses = engine.houses(
        t,
        obs,
        HouseOptions {
            system: HouseSystem::Placidus,
            polar_policy: PolarPolicy::FallbackToWholeSign,
            ..Default::default()
        },
    )?;
    let snapshot = engine.chart_snapshot(
        t,
        obs,
        &[BodyId::Planet(PlanetId::Sun), BodyId::Planet(PlanetId::Moon)],
        HouseOptions::default(),
        CalcOptions::default(),
        AspectPreset::Major,
        1.0,
    )?;
    println!("Sun={:.6}°, ASC={:.4}°, Lots={}", sun.lon_deg, houses.angles.asc_deg, snapshot.lots.len());
    Ok(())
}
```

### CLI (`lxephemeris-cli`)

```bash
lxephemeris-cli calc --jd 2460676.5 --bodies Sun,Moon,Chiron
lxephemeris-cli houses --jd 2460676.5 --lat 51.5074 --lon -0.1278 --system P
lxephemeris-cli panchanga --jd 2460676.5 --lat 28.6139 --lon 77.2090 --ayanamsa lahiri
lxephemeris-cli fetch-kernel --kernel de440s --data-dir ./data
```

---

## Verification & Benchmarks

Verified against **37,416 golden reference vectors** from **NASA JPL Horizons (`DE440` / `SB441-N16`)**, **IAU SOFA**, and **USNO** across **145+ Rust tests**, **18 Python test suites**, and **12 Node.js test suites** (`100% PASS`). Single-body position + analytical velocity runs in **`12.8 µs`**; a full natal chart runs in **`~150–240 µs`**.

👉 **[Full Verification & Benchmark Report (`ACCEPTANCE_REPORT.md`) →](ACCEPTANCE_REPORT.md)**

---

## Documentation

- **[Why `lxephemeris` & 22 Eliminated Legacy Bugs (`docs/WHY_LXEPHEMERIS.md`)](docs/WHY_LXEPHEMERIS.md)** — Architecture overview and legacy bug breakdown.
- **[Quickstart & Complete API Guide (`docs/QUICKSTART_AND_API.md`)](docs/QUICKSTART_AND_API.md)** — Full examples for Python, Node.js/TypeScript, Rust, and CLI.
- **[Migration Guide from `pyswisseph`, `PyEphem` & `Kerykeion` (`docs/SWE_MIGRATION_GUIDE.md`)](docs/SWE_MIGRATION_GUIDE.md)** — Side-by-side before/after migration patterns.
- **[CLI Reference (`crates/lxephemeris-cli/README.md`)](crates/lxephemeris-cli/README.md)** — All subcommands and flags for `lxephemeris-cli`.
- **[Intentional Divergences & Bug Proofs (`docs/SE_DIVERGENCES.md`)](docs/SE_DIVERGENCES.md)** — Mathematical proofs of all 22 upstream bug fixes.
- **[Mathematical Model Passports (`docs/MODEL_PASSPORTS.md`)](docs/MODEL_PASSPORTS.md)** — Equations, constants, and peer-reviewed citations.
- **[Verification & Benchmark Report (`ACCEPTANCE_REPORT.md`)](ACCEPTANCE_REPORT.md)** — 37,416-vector verification metrics and latency benchmarks.

---

## Licensing & Commercial Use

`lxephemeris` is dual-licensed ([`LICENSE.md`](LICENSE.md)). Choose the option that fits your deployment:

| Option | Best For | Commercial Rights | How to Start |
| :--- | :--- | :--- | :--- |
| **☁️ [AstroAPI (Managed Cloud API)](https://astrology-api.io/)** | Web/mobile apps, SaaS startups & teams wanting turnkey REST endpoints, city geocoding, SVG chart wheels & interpretations | ✅ **Included** in API subscription | **[Start at `astrology-api.io` →](https://astrology-api.io/)** |
| **🏢 `lxephemeris` Commercial License (`LicenseRef-Commercial`)** | Self-hosted backends, high-throughput batch pipelines & on-premise services running `lxephemeris` on your own servers | ✅ **Included** via offline `Ed25519` key (`LXEPHEMERIS_LICENSE_KEY`) | **[Request License (`procoders.tech/contacts`) →](https://procoders.tech/contacts/)** or **[`lxephemeris@procoders.tech`](mailto:lxephemeris@procoders.tech)** |
| **🎓 `lxephemeris` Non-Commercial ([`PolyForm-Noncommercial-1.0.0`](LICENSE-NONCOMMERCIAL.md))** | Personal projects, academic research, education, non-profits & 30-day pre-purchase commercial evaluation | ❌ Non-commercial & 30-day staging eval only | `pip install lxephemeris` / `npm i lxephemeris` |

👉 **[Full Dual-Licensing Policy & Offline Key Activation (`LICENSE.md`) →](LICENSE.md)**
