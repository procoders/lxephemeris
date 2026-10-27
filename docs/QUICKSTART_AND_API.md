# `lxephemeris` — Quickstart & Complete API Reference (Python, Node.js, Rust)

This guide provides self-contained, copy-pasteable examples for every capability of `lxephemeris` across Python (`lxephemeris`), Node.js / TypeScript (`lxephemeris`), and Rust (`lxephemeris-core`).

`lxephemeris` completely replaces **Swiss Ephemeris (`pyswisseph`)**, **PyEphem (`ephem`)**, and **Kerykeion (`kerykeion`'s astronomical calculation layer)** with:
- **`#![forbid(unsafe_code)]` & `Send + Sync` thread safety** (`0.0` ULP state bleed across 32 concurrent threads).
- **Dual Backend**: NASA JPL `DE440s`/`DE440`/`DE441`/`INPOP21a` SPK + fileless `DE440`-calibrated `TrigLadder` Moshier mode (`-3000 CE .. +3000 CE`).
- **Single-Pass Analytical Velocity Jacobian (`12.77 µs` per body with speeds; `242.56 µs` per full natal chart)**.
- **32,416-Vector Golden Reference Database** (`l3_se_10k_grid.json` [10,242], `l3_wave_bcf_10k_grid.json` [10,540], `l3_wave_g_extreme_10k_grid.json` [10,000], plus `1,634` L1A/L1B/L3 vectors) verified across **`105+` Rust tests**, **`15` Python test suites**, and **`11` Node.js test suites**.
- **22 Upstream Bug Immunity Guarantees**: Eliminates `SWE-01..SWE-11` (`libswe`), `EPH-01..EPH-05` (`PyEphem`), and `KER-01..KER-06` (`Kerykeion`).

---

## 1. Engine Initialization & Data Bundle Integrity

Initialize `Engine` with zero arguments for instant fileless operation (or automatic `LXEPHEMERIS_DATA_DIR` / `SE_EPHE_PATH` discovery), or pass a data directory with `verify_manifest_sha256=True` for `< 0.003"` NASA JPL `DE440s` + 68 Minor Bodies SPK precision.

### Python: Zero-File Instant Start vs. High-Precision SPK Data Directory

```python
import json
import lxephemeris

# 1a. Zero-File Instant Start (works immediately with 0 external files, or auto-loads LXEPHEMERIS_DATA_DIR / SE_EPHE_PATH)
engine = lxephemeris.Engine()

# 1b. Or load NASA JPL DE440s + 68 Minor Bodies (< 0.003" precision) with SHA-256 verification:
engine = lxephemeris.Engine.from_data_dir(
    "./data",
    max_threads=8,
    verify_manifest_sha256=True,
)

# Inspect loaded capabilities and per-body safe coverage intervals
caps = json.loads(engine.capabilities_json())
chiron_cov = json.loads(engine.coverage_json("Chiron"))
print("Minor bodies loaded:", caps["minor_bodies_loaded_count"], "Source:", caps["ephemeris_source"])
print("Chiron safe interval (JD TT):", chiron_cov["safe_start_jd_tt"], "to", chiron_cov["safe_end_jd_tt"])
```

### TypeScript / Node.js: Load `Engine` (`JsEngine`)

```typescript
import { Engine } from "lxephemeris";

// 1a. Zero-File Instant Start (0 external files required, or auto-loads LXEPHEMERIS_DATA_DIR / SE_EPHE_PATH)
const defaultEngine = new Engine();

// 1b. Or load NASA JPL DE440s + 68 Minor Bodies (< 0.003" precision):
const engine = Engine.fromDataDir("./data", 8, null, true);
const chironCoverage = JSON.parse(engine.coverageJson("Chiron"));
console.log("Chiron safe coverage:", chironCoverage);
```

### Rust: Load `lxephemeris_core::Engine`

```rust
use lxephemeris_core::{
    BodyId, CalcOptions, Engine, EngineConfig, HouseOptions, HouseSystem, Observer, PlanetId,
    PolarPolicy, TimeInput, JD_J2000,
};

// 1a. Zero-File Instant Start (auto-detects LXEPHEMERIS_DATA_DIR / SE_EPHE_PATH or uses embedded Moshier)
let default_engine = Engine::default();

// 1b. Or load NASA JPL DE440s + 68 Minor Bodies (< 0.003" precision):
let engine = Engine::from_data_dir(
    "./data",
    EngineConfig {
        max_threads: Some(8),
        pin_time_table_version: None,
        verify_manifest_sha256: true,
    },
)?;

let sun = engine.calc(
    BodyId::Planet(PlanetId::Sun),
    TimeInput::from_tt_jd(JD_J2000),
    CalcOptions::default(),
)?;
let obs = Observer::new(51.5074, -0.1278, 35.0)?;
let houses = engine.houses(
    TimeInput::from_tt_jd(JD_J2000),
    obs,
    HouseOptions {
        system: HouseSystem::Placidus,
        polar_policy: PolarPolicy::FallbackToWholeSign,
        ..Default::default()
    },
)?;
```

---

## 2. Time Scales, Local Timezones (`utc_from_local`), Calendar (`julday` / `revjul`), $\Delta T$, and Sidereal Time (`sidtime`)

`lxephemeris` supports `"UT1"`, `"UTC"`, `"TT"`, and `"TDB"` time scales, IANA local time conversion with DST fold disambiguation (`utc_from_local`), high-precision two-part Julian Dates (`SplitJd`), standalone $\Delta T = \text{TT} - \text{UT1}$ (`delta_t`), and Greenwich/Local Apparent Sidereal Time (`sidtime`).

### Python: Local Timezone Conversion, Calendar Conversion, $\Delta T$, and Sidereal Time / ARMC

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")

# Convert local wall-clock time (with IANA timezone and PEP 495 DST fold) to UTC Julian Date
jd_utc = lxephemeris.utc_from_local(2026, 10, 25, 2, 30, 0.0, tz_name="Europe/Berlin", fold=1)

# Convert Gregorian or Julian date to Julian Day
jd_ut = lxephemeris.julday(2026, 10, 24, 14.5, calendar="gregorian")
year, month, day, hour_ut = lxephemeris.revjul(jd_ut, calendar="gregorian")

# Query Delta-T (TT - UT1) and TT - TDB directly (replaces swe.deltat)
dt_info = engine.delta_t(jd_ut, scale="UT1")
print(f"Delta-T: {dt_info.delta_t_sec:.4f} s, TT-TDB: {dt_info.tt_minus_tdb_sec:.6f} s")

# Query GMST, GAST, Local Sidereal Time, and ARMC (replaces swe.sidtime)
st = engine.sidtime(jd_ut, scale="UT1", nutation="true", longitude_deg=37.6173)
print(f"GMST: {st.gmst_hours:.6f} h, GAST: {st.gast_hours:.6f} h, ARMC: {st.armc_deg:.4f} deg")
```

---

## 3. Planetary, Minor Body, Lunar Point, and Uranian Positions (`calc`)

`engine.calc(body, jd1, jd2=0.0, scale="TT", options=None)` computes coordinates and 4-point analytical speeds for:
- **10 Major Bodies**: `"Sun"`, `"Moon"`, `"Mercury"`, `"Venus"`, `"Mars"`, `"Jupiter"`, `"Saturn"`, `"Uranus"`, `"Neptune"`, `"Pluto"` (plus `"Earth"` in heliocentric/barycentric frames).
- **68 Minor Bodies**: By name (`"Ceres"`, `"Pallas"`, `"Juno"`, `"Vesta"`, `"Chiron"`, `"Pholus"`, `"Nessus"`, `"Eris"`, `"Sedna"`, `"Haumea"`, `"Makemake"`, etc.) or MPC number (`"mpc:2060"`, `"mpc:433"`).
- **7 Lunar Points**: `"Mean_Node"`, `"True_Node"`, `"Mean_South_Node"`, `"True_South_Node"`, `"Mean_Lilith"`, `"True_Lilith"` (osculating apogee), `"Interpolated_Lilith"` ($C^2$ cubic spline over natural lunar apogees).
- **9 Hypothetical Bodies**: `"Cupido"`, `"Hades"`, `"Zeus"`, `"Kronos"`, `"Apollon"`, `"Admetos"`, `"Vulkanus"`, `"Poseidon"`, `"Transpluto"` (`"Isis"`).
- **Ecliptic Obliquity & Nutation**: `"ECL_NUT"` or `engine.ecl_nut(jd)`.

### Python: Tropical, Sidereal, Equatorial (RA/Dec), Heliocentric, and `ECL_NUT`

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")
jd_ut = lxephemeris.julday(2026, 6, 21, 12.0)

# 1. Default Geocentric Apparent Tropical Ecliptic of Date (with speeds & km distance)
moon = engine.calc("Moon", jd_ut, scale="UT1")
print(f"Moon lon={moon.lon_deg:.6f}°, lat={moon.lat_deg:.6f}°, dist={moon.dist_km:.1f} km")
print(f"Moon speeds: lon={moon.lon_speed_deg_day:.6f} °/d, dist={moon.dist_speed_km_day:.2f} km/d")

# 2. Equatorial Right Ascension (RA), Declination (Dec), and Declination speed (for OOB & parallel aspects)
eq_opts = lxephemeris.CalcOptions(plane="equatorial", equinox="of_date", nutation="true", aberration="apparent")
mars_eq = engine.calc("Mars", jd_ut, scale="UT1", options=eq_opts)
print(f"Mars RA={mars_eq.lon_deg:.6f}°, Dec={mars_eq.lat_deg:.6f}°, dDec/dt={mars_eq.lat_speed_deg_day:.6f} °/d")

# 3. Sidereal (Vedic) calculation without nutation for station finding
sid_mean_opts = lxephemeris.CalcOptions(zodiac="sidereal", ayanamsa="lahiri", nutation="mean")
saturn_sid = engine.calc("Saturn", jd_ut, scale="UT1", options=sid_mean_opts)

# 4. Heliocentric Earth & Planets
helio_opts = lxephemeris.CalcOptions(origin="heliocentric", aberration="geometric")
earth_helio = engine.calc("Earth", jd_ut, scale="UT1", options=helio_opts)

# 5. True & Mean Obliquity and Nutation (ECL_NUT)
ecl = engine.ecl_nut(jd_ut, scale="UT1")
print(f"True eps={ecl.true_obliquity_deg:.6f}°, Mean eps0={ecl.mean_obliquity_deg:.6f}°, dpsi={ecl.nutation_lon_deg * 3600:.3f}\"")
```

---

## 4. Vectorized Parallel Batch Calculations (`calc_batch`, `calc_batch_flat_f64`, `houses_batch`)

For transits, station finding, void-of-course Moon scans, and 1-minute angular triggers (`1,440` steps/day), use batch methods to compute thousands of epochs in parallel across all CPU cores without crossing the Python/JS FFI boundary in a loop.

### Python: Parallel Body Batch & 1-Minute Angular Trigger Houses Batch

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")
jd0 = lxephemeris.julday(2026, 10, 24, 0.0)

# 1. Evaluate 10 bodies across 365 daily epochs in parallel (returns flat [365 * 10 * 6] float list)
bodies = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
daily_epochs = [jd0 + i for i in range(365)]
flat_states = engine.calc_batch_flat_f64(bodies, daily_epochs, scale="UT1")

# 2. Evaluate 1,440 1-minute house & angle steps in a single parallel batch call
minute_epochs = [jd0 + m / 1440.0 for m in range(1440)]
minute_houses = engine.houses_batch(
    minute_epochs,
    latitude_deg=51.5074,
    longitude_deg=-0.1278,
    system="P",
    scale="UT1",
    polar_policy="fallback_to_whole_sign",
)
print(f"Minute 0 ASC={minute_houses[0].asc_deg:.4f}°, Minute 1439 ASC={minute_houses[-1].asc_deg:.4f}°")
```

---

## 5. All 23 House Systems, `PolarPolicy`, `houses_armc`, and `house_of_lon`

Supported case-sensitive house system codes:
- `'P'` (Placidus), `'K'` (Koch), `'O'` (Porphyry), `'R'` (Regiomontanus), `'C'` (Campanus), `'A'` (Equal from ASC), `'W'` (Whole Sign), `'B'` (Alcabitius), `'M'` (Morinus), `'T'` (Topocentric Polich-Page), `'U'` (Krusinski-Pisa-Goelzer), `'X'` (Meridian / Axial), `'V'` (Vehlow Equal), `'L'` (Pullen Sinusoidal Delta), `'Q'` (Pullen Sinusoidal Ratio), `'I'` (Sunshine Treindl), `'i'` (Sunshine Makransky), `'F'` (Carter Poli-Equatorial), `'N'` (Whole Sign $0^\circ$ Aries), `'S'` (Sripati), `'D'` (Equal from MC), `'Y'` (APC), `'H'` (Azimuthal / Horizontal).

### Python: Polar Fallback (`fallback_to_whole_sign`), Composite `houses_armc`, and Boundary-Safe `house_of_lon`

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")
jd_ut = lxephemeris.julday(2026, 6, 21, 12.0)

# High-latitude calculation (Tromsø, 69.65°N) with automatic Whole Sign fallback on true latitude
h_polar = engine.houses(
    jd_ut,
    latitude_deg=69.6492,
    longitude_deg=18.9553,
    system="P",
    scale="UT1",
    polar_policy="fallback_to_whole_sign",
)
print("Requested:", h_polar.system_code, "Effective:", h_polar.effective_system_code)
print("Polar fallback triggered:", h_polar.polar_house_fallback)

# Composite chart houses directly from ARMC without a Julian Day
ecl = engine.ecl_nut(jd_ut, scale="UT1")
h_comp = lxephemeris.houses_armc(
    armc_deg=142.5,
    eps_true_deg=ecl.true_obliquity_deg,
    latitude_deg=40.7128,
    system="P",
    polar_policy="fallback_to_whole_sign",
)

# Assign a planet to a house (1..12) with 1e-9 deg half-open boundary tolerance
house_num = lxephemeris.house_of_lon(h_comp.cusps_deg[9] - 1e-11, h_comp.cusps_deg, eps_deg=1e-9)
assert house_num == 10  # Snaps cleanly to House 10 instead of falling into House 9
```

---

## 6. Fixed Stars, Parans (`rise_trans`), Horizontal `azalt`, and Heliacal Phenomena (`heliacal_phenomena`)

### Python: Fixed Star Catalog, Parans with Proper Motion, Sunrise/Sunset, and Heliacal Rise/Set

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")
jd_ut = lxephemeris.julday(2026, 5, 1, 0.0)

# 1. Query star metadata (visual magnitude V, Bayer code, HIP number, parallax, radial velocity)
sirius_meta = engine.star_info("Sirius")  # or by Bayer code "alCMa"
sirius_pos = engine.calc_star_by_name("Sirius", jd_ut, scale="UT1")
print(f"{sirius_meta.canonical_name} ({sirius_meta.bayer}): V={sirius_meta.vmag}, lon={sirius_pos.lon_deg:.6f}°")

# 2. Star paran (rise/set/upper/lower transit) with full Hipparcos-2 proper motion, disc center, and 0 pressure
sirius_rise_jd = engine.rise_trans(
    "Sirius",
    "rise",
    jd_ut,
    latitude_deg=30.0444,
    longitude_deg=31.2357,
    scale="UT1",
    disc_mode="center",
    pressure_hpa=0.0,
)

# 3. Standard Almanac Sunrise (upper limb, standard refraction -0°34' at sea level)
sun_rise_jd = engine.rise_trans(
    "Sun",
    "rise",
    jd_ut,
    latitude_deg=30.0444,
    longitude_deg=31.2357,
    scale="UT1",
    disc_mode="upper_limb",
    pressure_hpa=1013.25,
)

# 4. Solar altitude for Traditional Sect (day/night chart & |alt| < 5° near_horizon flag)
sun_hz = engine.calc_azalt("Sun", jd_ut, latitude_deg=30.0444, longitude_deg=31.2357, scale="UT1")
print(f"Sun geometric alt={sun_hz.altitude_true_deg:.4f}°, apparent={sun_hz.altitude_apparent_deg:.4f}°, near_horizon={sun_hz.near_horizon}")

# 5. Heliacal rising & setting by Bradley Schaefer's arcus visionis model (age 36, Snellen 1.0)
hel = engine.heliacal_phenomena(
    "Sirius",
    jd_ut,
    latitude_deg=30.0444,
    longitude_deg=31.2357,
    scale="UT1",
    observer_age_years=36.0,
    snellen_ratio=1.0,
)
print("Sirius Heliacal Rise JD(UT1):", hel.heliacal_rise_jd_ut1, "Heliacal Set JD(UT1):", hel.heliacal_set_jd_ut1)
```

---

## 7. Solar & Lunar Eclipses (Global, Local, Saros, Path Width, and Range Search)

### Python: Next Global Eclipse, Next Local Eclipse, and 6-Year Eclipse Range

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")
jd_start = lxephemeris.julday(2024, 1, 1, 0.0)

# 1. Next solar eclipse with Saros series/member, path width (km), central duration (s), and local circumstances
obs_dallas = lxephemeris.Observer(32.7767, -96.7970, 130.0)
sol = engine.next_solar_eclipse(jd_start, scale="TT", observer=obs_dallas)
print(f"{sol.kind} Solar Eclipse: Saros {sol.saros_series} #{sol.saros_member}, width={sol.path_width_km:.1f} km")
print(f"Dallas local: kind={sol.local_kind}, mag={sol.local_magnitude:.4f}, totality={sol.local_central_duration_sec:.1f} s")

# 2. Next solar eclipse strictly visible at a specific observer location
next_visible = engine.next_local_solar_eclipse(jd_start, 32.7767, -96.7970, scale="TT", altitude_m=130.0)

# 3. Search all solar and lunar eclipses over a 6-year window (/eclipses/range)
jd_end = jd_start + 6 * 365.25
solar_list = engine.solar_eclipses_in_range(jd_start, jd_end, observer=obs_dallas)
lunar_list = engine.lunar_eclipses_in_range(jd_start, jd_end, shadow_rule="danjon")
print(f"Found {len(solar_list)} solar and {len(lunar_list)} lunar eclipses in 6-year window")
```

---

## 8. Fileless Moshier Mode (`new_moshier`), Phenomena (`phenomena`), Nodes & Apsides (`nod_aps`), Lunar Occultations, 3D `house_pos`, Equation of Time (`time_equation`), and Custom Star Catalogs

```python
import lxephemeris

# 1. Fileless Analytical Moshier Mode (-3000 CE .. +3000 CE, 0 bytes of external .bsp files required)
moshier_engine = lxephemeris.Engine.new_moshier()
jd_ut = lxephemeris.julday(2026, 10, 24, 12.0)
sun_mosh = moshier_engine.calc("Sun", jd_ut, scale="UT1")

# 2. Planetary Phenomena & Apparent Visual Magnitudes (Mallama & Hilton 2018 / Bowell H-G)
engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")
venus_pheno = engine.phenomena("Venus", jd_ut, scale="UT1")
print(
    f"Venus phase={venus_pheno.phase_illuminated:.4f}, "
    f"elongation={venus_pheno.elongation_deg:.2f}°, Vmag={venus_pheno.apparent_magnitude_v:.2f}"
)

# 3. Planetary, Lunar & Minor-Body Nodes and Apsides (Mean, Osculating, Barycentric, 2nd Focal Point)
mars_nodes = engine.nod_aps("Mars", jd_ut, scale="UT1", method=lxephemeris.SE_NODBIT_OSCU)
print("Mars Ascending Node lon:", mars_nodes.ascending_node[0], "Perihelion lon:", mars_nodes.periapsis[0])

# 4. Global & Local Lunar Occultations of Fixed Stars, Planets & Minor Bodies
occ = engine.next_lunar_occultation("Regulus", jd_ut, scale="UT1")
print(f"Next Regulus occultation max JD(TT): {occ.max_occultation_jd_tt:.5f}, kind={occ.kind}")

# 5. 3D Mundane House Position (accounting for ecliptic latitude beta) & Equation of Time
ecl = engine.ecl_nut(jd_ut, scale="UT1")
st = engine.sidtime(jd_ut, scale="UT1", longitude_deg=-0.1278)
pluto = engine.calc("Pluto", jd_ut, scale="UT1")
pluto_hpos = engine.house_pos(st.armc_deg, 51.5074, ecl.true_obliquity_deg, pluto.lon_deg, pluto.lat_deg, hsys="P")
eot_minutes = engine.time_equation(jd_ut, scale="UT1") * 1440.0
print(f"Pluto 3D Placidus house pos={pluto_hpos:.4f}, Equation of Time={eot_minutes:+.2f} min")
```

---

## 9. Astronomical Event Solvers, Chinese 24 Solar Terms, Human Design `-88°` Epoch, Lunar Phase Summary, Stations, Declination Extrema, & Brady Parans (`Wave F`)

`lxephemeris` includes built-in safeguarded Newton-Illinois root solvers and diurnal paran engines that completely replace `ephem` (`next_vernal_equinox`, `next_new_moon`, `Sun().hlong` solar terms) and `kerykeion` (`calculate_moon_phase`):

### Python: Event Solvers & Brady Fixed-Star Parans

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")
jd_2026 = lxephemeris.julday(2026, 1, 1, 0.0)

# 1. Exact Ecliptic Longitude Crossing & 3-Pass Retrograde Returns
aries_ingress = engine.search_body_longitude("Sun", jd_2026, target_lon_deg=0.0, forward=True, max_search_days=120.0)
merc_returns = engine.body_longitude_crossings_in_range("Mercury", jd_2026, jd_2026 + 365.25, target_lon_deg=15.0)

# 2. Chinese 24 Solar Terms (Jieqi) & 4 Cardinal Seasons (eliminates ephem.Sun().hlong 180° bug)
solar_terms = engine.solar_terms_in_range(jd_2026, jd_2026 + 365.25)
seasons = engine.seasons_in_year(2026)

# 3. Human Design -88° Prenatal Solar Arc Epoch
hd_design = engine.human_design_epoch(jd_2026 + 100.0)

# 4. Lunar Phases & Complete LunarPhaseSummary (8-phase & Kerykeion 28-phase indices)
next_full = engine.search_lunar_phase(jd_2026, target_elongation_deg=180.0, forward=True)
year_phases = engine.lunar_phases_in_range(jd_2026, jd_2026 + 365.25)
moon_summary = engine.lunar_phase_summary(jd_2026 + 10.0)
print(
    f"Moon Phase: {moon_summary.phase_8_name} (8-idx={moon_summary.phase_8_index}, "
    f"28-idx={moon_summary.phase_28_index}, age={moon_summary.age_days:.2f} d)"
)

# 5. Planetary Stations (d_lon/dt = 0) & Declination Extrema / Standstills (d_dec/dt = 0)
merc_stations = engine.search_stations_in_range("Mercury", jd_2026, jd_2026 + 365.25)
moon_standstills = engine.search_declination_extrema_in_range("Moon", jd_2026, jd_2026 + 30.0)

# 6. Bernadette Brady Diurnal Fixed-Star Parans (airless geometric center, 52 embedded stars + custom stars × 10 planets)
obs_london = lxephemeris.Observer(51.5074, -0.1278, 0.0)
parans = engine.find_parans(jd_2026, obs_london, orb_minutes=4.0)
for p in parans[:5]:
    print(f"{p.star_name} ({p.star_angle}) + {p.planet_name} ({p.planet_angle}): orb={p.time_orb_minutes:.2f} min")
```

### TypeScript / Node.js: Event Solvers & Brady Fixed-Star Parans

```typescript
import { Engine, julday } from "lxephemeris";

const engine = new Engine();
const jd2026 = julday(2026, 1, 1, 0.0, "gregorian");

const seasons = engine.seasonsInYear(2026);
const solarTerms = engine.solarTermsInRange(jd2026, jd2026 + 365.25);
const hdEpoch = engine.humanDesignEpoch(jd2026 + 100.0);
const moonSummary = engine.lunarPhaseSummary(jd2026 + 10.0);
const mercStations = engine.searchStationsInRange("Mercury", jd2026, jd2026 + 365.25);
const parans = engine.findParans(
  jd2026,
  { latitudeDeg: 51.5074, longitudeDeg: -0.1278, altitudeM: 0 },
  4.0
);
console.log(`Found ${solarTerms.length} solar terms, ${mercStations.length} Mercury stations, ${parans.length} parans`);
```

---

## 10. Planetary Returns, Aspects & Transits, Vedic Panchanga / Dasha / Vargas, Hermetic Lots, 64 HD Gates, & `chart_snapshot` (`Wave H`)

### Python: Returns, Aspects, VoC Moon, Panchanga, Vimshottari Dasha, Shodashavarga, Lots, & `chart_snapshot`

```python
import lxephemeris

engine = lxephemeris.Engine()
jd_birth = lxephemeris.utc_from_local(1990, 6, 21, 12, 0, 0.0, tz_name="Europe/Berlin")
jd_2026 = lxephemeris.utc_from_local(2026, 1, 1, 0, 0, 0.0, tz_name="UTC")
obs = lxephemeris.Observer(51.5074, -0.1278, 0.0)

# 1. Planetary Returns (Tropical & Sidereal Solar Return, plus 3-pass retrograde Saturn Return)
solar_return = engine.planetary_return("Sun", jd_birth, jd_2026)
saturn_returns = engine.planetary_returns_in_range("Saturn", jd_birth, jd_birth + 28 * 365.25, jd_birth + 31 * 365.25)

# 2. Static & Dynamic Aspects, Synastry, Declination Parallels, and Void-of-Course Moon
chart_aspects = engine.calc_chart_aspects(jd_birth, preset="major", orb_scale=1.0)
synastry = engine.calc_synastry_aspects(jd_birth, ["Sun", "Moon", "Venus"], jd_2026, ["Sun", "Moon", "Mars"])
decl_parallels = engine.calc_declination_aspects(jd_birth, max_orb_deg=1.2)
exact_conj = engine.search_exact_aspect("Moon", "Sun", 0.0, jd_2026, forward=True)
mutual_hits = engine.mutual_aspects_in_range(jd_2026, jd_2026 + 7.0, ["Sun", "Moon", "Mercury"], preset="major")
natal_transits = engine.natal_transits_in_range(jd_birth, jd_2026, jd_2026 + 7.0, ["Moon"], ["Sun", "Moon"])
voc_periods = engine.void_of_course_moon_in_range(jd_2026, jd_2026 + 7.0)

# 3. Vedic Jyotish: Panchanga, 120-Year Vimshottari Dasha, & 16 Shodashavarga Divisional Charts (D1..D60)
panch = engine.panchanga(jd_2026, obs, ayanamsa="lahiri")
dasha = engine.vimshottari_dasha(jd_birth, ayanamsa="lahiri")
d9_navamsha = engine.calc_varga(123.45, 9)
all_16_vargas = engine.calc_shodashavarga(123.45)

# 4. Traditional Sect, 24 Unequal Planetary Hours, 7 Hermetic Lots, 64 HD Gates, & Single-Call ChartSnapshot
sect = engine.chart_sect(jd_birth, obs)
hours_24 = engine.planetary_hours(jd_birth, obs)
lots_7 = engine.calc_lots(jd_birth, obs)
hd_placement = lxephemeris.map_zodiac_placement(302.0)  # 2°00' Aquarius -> Gate 41, Line 1
snapshot = engine.chart_snapshot(jd_birth, obs)
```

### TypeScript / Node.js: Returns, Aspects, Panchanga, Lots, & `chartSnapshot`

```typescript
import { Engine, utcFromLocal, mapZodiacPlacement } from "lxephemeris";

const engine = new Engine();
const jdBirth = utcFromLocal(1990, 6, 21, 12, 0, 0.0, "Europe/Berlin");
const jd2026 = utcFromLocal(2026, 1, 1, 0, 0, 0.0, "UTC");
const obs = { latitudeDeg: 51.5074, longitudeDeg: -0.1278, altitudeM: 0 };

const solarReturn = engine.planetaryReturn("Sun", jdBirth, jd2026);
const chartAspects = engine.calcChartAspects(jdBirth, null, "major", 1.0);
const vocMoon = engine.voidOfCourseMoonInRange(jd2026, jd2026 + 7.0);
const panch = engine.panchanga(jd2026, obs, "lahiri");
const dasha = engine.vimshottariDasha(jdBirth, "lahiri");
const shodasha = engine.calcShodashavarga(123.45);
const ph = engine.planetaryHours(jdBirth, obs);
const lots = engine.calcLots(jdBirth, obs, { system: "P" });
const hdGate41 = mapZodiacPlacement(302.0);
const snapshot = engine.chartSnapshot(jdBirth, obs, null, { system: "P" });
```

---

## 11. Verification Suite (`32,416` Golden Vectors) & 22 Upstream Bug Immunity Guarantees

Every API shown in Sections 1–10 is continuously verified against **32,416 independent golden reference vectors** across **110+ Rust tests**, **16 Python test suites**, and **12 Node.js test suites**:
- **`l3_se_10k_grid.json` (`10,242` vectors)**: 10 planets, 68 minor bodies, 7 lunar points, 9 Uranians, 52 embedded fixed stars (+ custom catalogs), 12 core ayanamsas, and 23 house systems across `1850..2100`.
- **`l3_wave_bcf_10k_grid.json` (`10,540` vectors)**: `rise_trans` (`almanac_upper_limb` & `geometric_center`), principal lunar phases, 8-phase + Kerykeion 28-phase `lunar_phase_summary`, 24 Chinese Solar Terms (*Jieqi*) & 4 cardinal seasons, Human Design `-88°` prenatal epochs, planetary stations ($\dot{\lambda}=0$), declination extrema / lunistices ($\dot{\delta}=0$), Bernadette Brady diurnal fixed-star parans, Moshier (`-500..+2500 CE`), `pheno_ut`, `nod_aps_ut`, 43 `SE_SIDM_*` ayanamsas, 3D `house_pos`, `time_equ`, and `azalt_rev`.
- **`l3_wave_g_extreme_10k_grid.json` (`10,000` vectors)**: Extreme boundary conditions and regression guards for all **22 eliminated upstream bugs** (`SWE-01..SWE-11` from Swiss Ephemeris, `EPH-01..EPH-05` from PyEphem, and `KER-01..KER-06` from Kerykeion).
- **L1A / L1B / L3 Baseline Suites (`1,634` vectors)**: Direct NASA JPL Horizons `DE440` ICRF (`360`), `DE440s` IAU 2006/2000A (`90`), `SB441-N16` 68 minor bodies (`340`), 12 ayanamsas (`120`), 23 house systems (`253`), 7 lunar points (`63`), 9 Uranians (`54`), Hipparcos-2 + XHIP stars (`104`), and L3 reference (`250`).

See **[`docs/SE_DIVERGENCES.md`](SE_DIVERGENCES.md)** for the mathematical breakdown of all 22 upstream bugs and **[`docs/SWE_MIGRATION_GUIDE.md`](SWE_MIGRATION_GUIDE.md)** for side-by-side migration patterns.


