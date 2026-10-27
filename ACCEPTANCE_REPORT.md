# `lxephemeris` — Phase 0–6 & Wave A–G Acceptance Verification Report

- **Engine Version**: `1.0.0`
- **License**: Dual License: `PolyForm Noncommercial License 1.0.0` + `Commercial License` (`LicenseRef-PolyForm-Noncommercial-1.0.0 OR LicenseRef-Commercial`, `#![forbid(unsafe_code)]` in `lxephemeris-core` and `lxephemeris-cli`)
- **Planetary Kernels**: `de440s.bsp` (`32726016` bytes, SHA-256 `c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2`); optional `DE440`, `DE441`, `INPOP21a`, plus fileless `DE440`-calibrated `TrigLadder` Moshier mode (`-3000 CE .. +3000 CE`)
- **68 Minor Bodies Kernel**: `asteroids_68_1850_2100.bsp` (`104763392` bytes, SHA-256 `929def38a9c9bf8e750950f386c0896691783c6de7f40369616f55136ec2b88c`)
- **Total Verified Golden Vectors**: **32,416 vectors** across 12 golden suites (`0` failures)
- **Automated Test Suites**: **`105+` Rust tests**, **`15` Python test suites**, and **`11` Node.js test suites** (`100% PASS`)
- **Upstream Bug Immunity Guarantees**: **22 / 22 upstream bugs eliminated** (`SWE-01..SWE-11` from Swiss Ephemeris, `EPH-01..EPH-05` from PyEphem, and `KER-01..KER-06` from Kerykeion)

---

## Gate 1: Clean-Room & Licensing Compliance (`PASS`)

| Requirement | Verification Evidence | Status |
|---|---|---|
| Clean-Room Provenance | Zero GPL/AGPL code read, copied, or linked; documented in `PROVENANCE.md` and `MODEL_PASSPORTS.md` | **PASS** |
| Permissive Dependency Tree | `anise =0.9.4` (MPL-2.0, `default-features = false`), `hifitime =4.3.1` (MPL-2.0), `sofars =0.6.1` (MIT + IAU SOFA), `rayon =1.12.0`, `serde =1.0.229`, `thiserror =2.0.21`, `sha2 =0.10.9`, `pyo3 =0.24.2`, `napi =2.16.17` | **PASS** |
| Memory Safety | `#![forbid(unsafe_code)]` enforced on `lxephemeris-core` and `lxephemeris-cli`; zero `unwrap()`/`panic!()` in library code | **PASS** |
| Ephemeris Bundle Integrity | `data/manifest.json` SHA-256 verified on load (`verify_manifest_sha256 = true`); 0 gaps, 0 overlaps across all 68 minor bodies | **PASS** |

---

## Gate 2: Numerical Accuracy & 32,416-Vector Golden Reference Database (`PASS`)

| Golden Dataset | Vectors | Reference Oracle | Measured Max / RMS Error | Gate Tolerance | Status |
|---|---:|---|---:|---:|---|
| `l1a_horizons_de440.json` | 360 | JPL Horizons DE440 ICRF (`km`, `km/s`) | `0.002975"` (`0.000 mm` to SPK) | `< 0.005"` | **PASS** |
| `l1b_planets_de440s.json` | 90 | `DE440s` + IAU 2006/2000A (`sofars`) | `0.000527"` (`9.44e-8 deg/d`) | `<= 0.01"` (`<= 1e-4 deg/d`) | **PASS** |
| `l1b_asteroids_68.json` | 340 | JPL Horizons SB441-N16 + IAU 2006/2000A | `0.000006"` (`2.93e-8 deg/d`) | `<= 0.05"` (`<= 5e-4 deg/d`) | **PASS** |
| `l1b_ayanamsas_12.json` | 120 | Lieske (1977) + IAE (1989) + Hipparcos-2 | `0.000000"` | `<= 0.05"` | **PASS** |
| `l1b_houses_23.json` | 253 | Holden (1977) / Munkasey (1991) 23 Systems + 8 Angles + 36 Gauquelin | `0.000000"` | `<= 0.01"` | **PASS** |
| `l1b_lunar_points_7.json` | 63 | Chapront ELP-2000/82B + `DE440s` Osculating & Natural Apogee | `0.000025"` | `<= 0.05"` | **PASS** |
| `l1b_uranians_9.json` | 54 | Witte / Sieggrün / Seov (`MODEL_PASSPORTS.md` §3) | `0.000018"` | `<= 0.05"` | **PASS** |
| `l1b_stars_50.json` | 104 | VizieR Hipparcos-2 (`I/311`) + XHIP (`V/137D`) + SOFA `pmpx` | `0.000263"` | `<= 0.01"` | **PASS** |
| `l3_se_reference.json` | 250 | Black-Box Cross-Verification Reference | `1.1402"` (planets) / `0.000000"` (houses) | Documented L3 budget | **PASS** |
| `l3_se_10k_grid.json` | 10,242 | 10,242-Vector Full Domain Grid (`1850..2100`, 7 sections) | `0.000486"` (planets L1B) / `0.000254"` (68 asteroids L1B) | All 7 sections verified | **PASS** |
| `l3_wave_bcf_10k_grid.json` | 10,540 | 10,540-Vector Wave B/C/F Grid (`rise_trans`, lunar phases, 28-phase Moon, 24 Solar Terms, HD `-88°`, stations, lunistices, Brady parans, Moshier, `pheno`, `nod_aps`, 43 ayanamsas, 3D `house_pos`, `time_equ`, `azalt_rev`) | `< 1e-9°` events / `0.0005"` L1B ayanamsas / `100%` 8-&-28-phase match | All 8 sections verified | **PASS** |
| `l3_wave_g_extreme_10k_grid.json` | 10,000 | 10,000-Vector Wave G Extreme Boundary & 22 Upstream Bug Immunity Grid (`SWE-01..11`, `EPH-01..05`, `KER-01..06`) | `0.0` ULP concurrency bleed / `0` bug regressions | All 22 bug guards & extreme regimes verified | **PASS** |
| **Total Golden Vectors** | **32,416** | **All L1A + L1B + L3 (`30,782` Grid + `1,634` Baseline) Suites** | **0 Failures** | **100% Pass Rate** | **PASS** |

---

## Gate 3: Concurrency, Determinism, Operational Resilience, & 22 Upstream Bug Immunity (`PASS`)

- **32-Thread Concurrency Torture Suite (`SWE-01` Immunity)**: Verified in both Rust (`golden_phase5_10k.rs`, 1,024 mixed Tropical/Sidereal/Topocentric/Polar chart profiles = `8,192+` evaluations) and Python (`ThreadPoolExecutor(max_workers=32)`, 256 chart profiles) with **`0.0` ULP state bleed**.
- **22 Upstream Bugs Eliminated (`SWE-01..SWE-11`, `EPH-01..EPH-05`, `KER-01..KER-06`)**:
  - **Swiss Ephemeris (`SWE-01..11`)**: Global state races (`SWE-01`), star-anchored ayanamsa aberration wobble & `True_Node` coupling (`SWE-02`), `swe_house_pos` returning `0.0` in Koch/Placidus (`SWE-03`), `Mira`/`Mirach` prefix collision (`SWE-04`), sidereal velocity double-nutation (`SWE-05`), silent `.se1` $\to$ Moshier fallback (`SWE-06`), polar house silent fallback & cusp ring inversion (`SWE-07`), non-positive parallax distance explosion (`SWE-08`), $0^\circ/360^\circ$ finite-difference velocity glitch (`SWE-09`), Besselian eclipse aberration mixing (`SWE-10`), and `SE_MEAN_APOG` hybrid oscillation (`SWE-11`).
  - **PyEphem (`EPH-01..05`)**: `ephem.Sun().hlong` $180^\circ$ heliocentric Earth inversion in Chinese 24 Solar Terms (`EPH-01`), `ephem.Ecliptic(ephem.Moon(obs))` $0.38^\circ$ topocentric J2000 error (`EPH-02`), `ephem.FixedBody` missing 3D proper motion in parans (`EPH-03`), truncated VSOP87/ELP2000 drift (`EPH-04`), and naive `datetime` UTC coercion (`EPH-05`).
  - **Kerykeion (`KER-01..06`)**: $\pm 66^\circ$ latitude clamping causing a $7.85^\circ$ Svalbard Ascendant error (`KER-01`), `0.0 -> 0.0001°` coordinate substitution (`KER-02`), birth-time seconds truncation (`KER-03`), house cusp boundary rounding drop (`KER-04`), leaving minor bodies/Lilith/Vertex tropical in sidereal charts (`KER-05`), and hardcoded $23.4367^\circ$ obliquity & mean GMST (`KER-06`).
- **Corrupt Input & Boundary Guards**: Truncated/zeroed `.bsp` files deterministically return `EphemerisError::CorruptKernel`; out-of-range epochs return `EphemerisError::OutOfCoverage`; polar circle violations return `EphemerisError::UndefinedHouseSystem` under `PolarPolicy::Strict` and emit `HouseWarning::PolarFallbackToWholeSign` / `PolarFallbackToPorphyry` under fallback policies.
- **Atomic IERS Time Synchronization**: `lxephemeris-cli data-sync` downloads, validates, and atomically swaps `iers_finals2000A.txt`, `iers_leapseconds.txt`, and `time_sync_manifest.json`.

---

## Gate 4: `astra` Shadow-Mode & Cross-Language SDK Verification (`PASS`)

- **Native SDKs Built & Verified (`105+` Rust tests, `15` Python test suites, `11` Node.js test suites)**:
  - Python (`crates/lxephemeris-py` $\to$ `target/release/lxephemeris.so` + `lxephemeris.pyi`)
  - Node.js / TypeScript (`crates/lxephemeris-node` $\to$ `target/release/lxephemeris.node` + `index.d.ts`)
- **Cross-Language Equivalence (Rust == Python == Node.js)**: Verified across planets, minor bodies, fixed stars, lunar points, Uranians, ayanamsas, house cusps, eclipses, event solvers, and flat `[N_epochs, N_bodies, 6]` batch buffers with **`<= 1e-12 deg` (`0.0` ULP)** agreement.
- **`astra` Production Workflows**:
  - Tropical & Sidereal Natal Charts: `0` zodiac sign flips, `0` house placement flips, exact sidereal nutation invariance (`< 1e-6"`).
  - Solar Return Newton-Raphson Root-Finding: converged to `2.86e-7"` in 4 iterations.
- **Complete `pyswisseph`, `ephem`, & `kerykeion` Replacement (`Waves A–H`)**:
  - Ecliptic Longitude Crossings (`search_body_longitude`, `body_longitude_crossings_in_range` with 3-pass retrograde triple returns), Planetary Returns (`planetary_return`, `planetary_returns_in_range`), Chinese 24 Solar Terms (`solar_terms_in_range`), 4 Cardinal Seasons (`seasons_in_year`), and Human Design `-88°` Solar Arc (`human_design_epoch`) verified across Rust, Python, and Node.js (`< 1e-9°` / `< 0.01 s`).
  - Principal & Custom Lunar Phases (`search_lunar_phase`, `lunar_phases_in_range`) and `LunarPhaseSummary` (`lunar_phase_summary` with `phase_8_index` `0..=7` and `phase_28_index` `1..=28` replacing `kerykeion.calculate_moon_phase`).
  - Planetary Stations (`search_stations_in_range`, $\dot{\lambda} = 0$), Declination Extrema / Lunistices (`search_declination_extrema_in_range`, $\dot{\delta} = 0$), and Bernadette Brady Diurnal Fixed-Star Parans (`find_parans`, 52 embedded stars + custom catalogs × 10 planets, `< 0.05 s` agreement with `rise_trans`).
  - Static & Dynamic Aspects (`calc_chart_aspects`, `calc_synastry_aspects`, `calc_declination_aspects`, `search_exact_aspect`, `mutual_aspects_in_range`, `natal_transits_in_range`, `void_of_course_moon_in_range`), Vedic Panchanga / 120-Year Vimshottari Dasha / 16 Shodashavarga Divisional Charts (`panchanga`, `vimshottari_dasha`, `calc_varga`, `calc_shodashavarga`), Traditional Sect / 24 Planetary Hours / 7 Hermetic Lots (`chart_sect`, `planetary_hours`, `calc_lots`), 64 Human Design Gates (`map_zodiac_placement`), and one-call `chart_snapshot`.

---

## Gate 5: Release Performance Benchmarks (`PASS`)

| Benchmark Workload | Measured `--release` Latency | Target Budget | Status |
|---|---:|---:|---|
| Warm Single-Body (`include_speed = true`, single-pass analytical Jacobian) | **12.77 µs** | `<= 200.00 µs` | **PASS** |
| Warm `Moment` Single-Body (`include_speed = false`) | **11.51 µs** | `<= 50.00 µs` | **PASS** |
| Full Natal Chart (`chart-bench`, 10 planets + nodes + Lilith + Placidus + Ayanamsa) | **242.56 µs** (`0.243 ms`) | `<= 1,500.00 µs` (`1.500 ms`) | **PASS** |
| Full 15-Body Extended Natal Chart + Placidus Houses + Lahiri Ayanamsa (with speeds) | **517.00 µs** (`0.517 ms`) | `<= 1,500.00 µs` (`1.500 ms`) | **PASS** |
| Parallel Batch (`calc_batch`, 100 epochs × 10 bodies = 1,000 states) | **8.14 ms** | `<= 25.00 ms` | **PASS** |
