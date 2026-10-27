# Complete Migration Guide: Replacing `pyswisseph`, `PyEphem` (`ephem`), and `Kerykeion` with `lxephemeris`

`lxephemeris` is engineered to **completely replace all three** legacy astronomical and astrological dependencies in modern Python, Node.js / TypeScript, and Rust backends:
1. **Swiss Ephemeris (`pyswisseph` / `swisseph` / `libswe`)** — eliminating global state race conditions (`SWE-01`), silent `.se1` $\to$ Moshier fallback (`SWE-06`), `Mira`/`Mirach` collisions (`SWE-04`), and `swe_house_pos` `0.0` failures (`SWE-03`).
2. **PyEphem (`ephem`)** — eliminating the $180^\circ$ `ephem.Sun().hlong` Chinese solar-term inversion (`EPH-01`), the $0.38^\circ$ topocentric `ephem.Ecliptic(ephem.Moon(obs))` J2000 bug (`EPH-02`), static `FixedBody` star parans (`EPH-03`), truncated analytical series drift (`EPH-04`), and naive `datetime` UTC coercion (`EPH-05`).
3. **Kerykeion (`kerykeion` astronomical calculation layer)** — eliminating $\pm 66^\circ$ silent latitude clamping (`KER-01`, which causes a $7.85^\circ$ Ascendant error in Svalbard), `0.0 -> 0.0001°` coordinate substitution (`KER-02`), birth-time seconds truncation (`KER-03`), house cusp rounding drops (`KER-04`), partial sidereal conversion (`KER-05`), and hardcoded $23.4367^\circ$ obliquity (`KER-06`).

---

## Two Migration Paths

1. **Path A — Drop-In `SweAdapter` (`lxephemeris.SweAdapter` / `JsSweAdapter`)**:
   Preserves standard Swiss Ephemeris integer body constants (`SE_SUN..SE_INTP_APOG`, `SE_CUPIDO..SE_ISIS`, `SE_AST_OFFSET + mpc`, `SE_ECL_NUT`), `SEFLG_*` bitflags (including `SEFLG_MOSEPH`), `SE_SIDM_*` ayanamsa IDs (`0..=42` + `SE_SIDM_USER = 255`), and `(values6, retflag)` return tuples — **scoped strictly to the adapter instance with zero process-global state**.
2. **Path B — Native Stateless `Engine` API (`lxephemeris.Engine` / `JsEngine`)**:
   Strongly typed, self-documenting methods (`calc`, `calc_batch`, `houses`, `houses_batch`, `ayanamsa`, `rise_trans`, `heliacal_phenomena`, `solar_eclipses_in_range`, `solar_terms_in_range`, `lunar_phase_summary`, `search_stations_in_range`, `find_parans`) that eliminate bitmask flags and return frozen structured objects.

---

## 1. Migrating from `pyswisseph` (`swisseph`)

### 1.1 Before vs. After (Python — Drop-In `SweAdapter` & Native `Engine`)

#### ❌ Before (`pyswisseph` — Global Mutable State & Process-Wide Lock Required)

```python
import threading
import swisseph as swe

# Global state in C `libswe` mutates across threads unless guarded by a lock!
_SWE_LOCK = threading.RLock()

with _SWE_LOCK:
    swe.set_ephe_path("/path/to/ephe")
    jd_ut = swe.julday(2026, 10, 24, 12.0, swe.GREG_CAL)
    flags = swe.FLG_SWIEPH | swe.FLG_SPEED

    # Tropical Sun & Asteroid 433 Eros
    xx_sun, ret = swe.calc_ut(jd_ut, swe.SUN, flags)
    xx_eros, _ = swe.calc_ut(jd_ut, swe.AST_OFFSET + 433, flags)

    # Sidereal Moon & Houses (mutates global sidereal mode for the entire process!)
    swe.set_sid_mode(swe.SIDM_LAHIRI)
    ay_deg = swe.get_ayanamsa_ut(jd_ut)
    xx_moon_sid, _ = swe.calc_ut(jd_ut, swe.MOON, flags | swe.FLG_SIDEREAL)
    cusps_sid, ascmc_sid = swe.houses_ex(jd_ut, 28.6139, 77.2090, b"P", swe.FLG_SIDEREAL)
```

#### ✅ After — Option A: Drop-In `lxephemeris.SweAdapter` (Zero Global State, No Lock Needed)

```python
import lxephemeris

# Load immutable Engine once at startup; create per-request or per-worker SweAdapter
engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data", verify_manifest_sha256=True)
swe = lxephemeris.SweAdapter(engine)

jd_ut = swe.julday(2026, 10, 24, 12.0, gregflag=1)
flags = lxephemeris.SEFLG_SWIEPH | lxephemeris.SEFLG_SPEED

# 1. Tropical Sun, Asteroid 433 Eros, & Ecliptic Obliquity/Nutation (SE_ECL_NUT = -1)
xx_sun, retflag = swe.calc_ut(jd_ut, lxephemeris.SE_SUN, flags)
xx_eros, _ = swe.calc_ut(jd_ut, lxephemeris.SE_AST_OFFSET + 433, flags)
xx_nut, _ = swe.calc_ut(jd_ut, lxephemeris.SE_ECL_NUT, 0)

# 2. Fixed Stars (exact lookup — immune to SWE-04 Mira/Mirach prefix collision!)
xx_mira, star_name, _ = swe.fixstar2_ut("Mira", jd_ut, flags)

# 3. Sidereal Mode (scoped ONLY to this `swe` adapter instance — 0.0 ULP bleed across threads!)
swe.set_sid_mode(lxephemeris.SE_SIDM_LAHIRI)
ay_deg, _ = swe.get_ayanamsa_ex_ut(jd_ut, lxephemeris.SEFLG_SWIEPH)
xx_moon_sid, _ = swe.calc_ut(jd_ut, lxephemeris.SE_MOON, flags | lxephemeris.SEFLG_SIDEREAL)
cusps_sid, ascmc_sid = swe.houses_ex(jd_ut, 28.6139, 77.2090, "P", lxephemeris.SEFLG_SIDEREAL)

# 4. Phenomena, Nodes/Apsides, 3D House Position (immune to SWE-03 0.0 bug), & Occultations
attr_20, _ = swe.pheno_ut(jd_ut, lxephemeris.SE_VENUS, flags)
xnasc, xndsc, xperi, xaphe, _ = swe.nod_aps_ut(jd_ut, lxephemeris.SE_MARS, flags, lxephemeris.SE_NODBIT_OSCU)
hpos_3d = swe.house_pos(ascmc_sid[2], 28.6139, xx_nut[0], [xx_sun[0], xx_sun[1]], "P")
occ_flag, occ_tret = swe.lun_occult_when_glob(jd_ut, 0, starname="Regulus", iflag=flags)

# 5. Fileless Analytical Moshier Mode (-3000 CE .. +3000 CE, 0 bytes of external .bsp files)
swe_mosh = lxephemeris.SweAdapter.new_moshier()
xx_sun_mosh, _ = swe_mosh.calc_ut(jd_ut, lxephemeris.SE_SUN, lxephemeris.SEFLG_MOSEPH | lxephemeris.SEFLG_SPEED)
```

#### ✅ After — Option B: Native Stateless `lxephemeris.Engine` (Recommended for New Code)

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data", verify_manifest_sha256=True)
jd_ut = lxephemeris.julday(2026, 10, 24, 12.0)

# Tropical & Sidereal queries can run concurrently across 32+ threads on the same `engine`
sun = engine.calc("Sun", jd_ut, scale="UT1")
eros = engine.calc("Eros", jd_ut, scale="UT1")  # or "mpc:433"

sid_opts = lxephemeris.CalcOptions(zodiac="sidereal", ayanamsa="lahiri", include_speed=True)
moon_sid = engine.calc("Moon", jd_ut, scale="UT1", options=sid_opts)
houses_sid = engine.houses(
    jd_ut, 28.6139, 77.2090, system="P", scale="UT1", zodiac="sidereal", ayanamsa="lahiri"
)
```

### 1.2 Before vs. After (Node.js / TypeScript)

```typescript
import {
  JsEngine,
  JsSweAdapter,
  SE_SUN,
  SE_MOON,
  SEFLG_SWIEPH,
  SEFLG_SPEED,
  SEFLG_SIDEREAL,
  SE_SIDM_LAHIRI,
} from "lxephemeris";

const engine = JsEngine.fromDataDir("/path/to/lxephemeris/data", 4, null, true);

// Path A: Instance-scoped JsSweAdapter
const swe = new JsSweAdapter(engine);
const jdUt = JsSweAdapter.julday(2026, 10, 24, 12.0, 1);
const sunRes = swe.calcUt(jdUt, SE_SUN, SEFLG_SWIEPH | SEFLG_SPEED);

swe.setSidMode(SE_SIDM_LAHIRI);
const moonSid = swe.calcUt(jdUt, SE_MOON, SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_SIDEREAL);
const housesSid = swe.housesEx(jdUt, 28.6139, 77.2090, "P", SEFLG_SIDEREAL);

// Path B: Native Stateless JsEngine
const sunNative = engine.calc("Sun", { jd1: jdUt, jd2: 0 }, "UT1");
const moonNativeSid = engine.calc("Moon", { jd1: jdUt, jd2: 0 }, "UT1", {
  zodiac: "sidereal",
  ayanamsa: "lahiri",
  includeSpeed: true,
});
```

---

## 2. Migrating from `PyEphem` (`ephem`)

`PyEphem` is frequently used alongside `pyswisseph` for almanac sunrise/sunset, solstice/equinox dates, lunar phases, Chinese 24 Solar Terms (*Jieqi*), and fixed-star diurnal parans. However, `PyEphem` suffers from 5 severe upstream bugs (`EPH-01..EPH-05`) that `lxephemeris` eliminates out of the box.

### 2.1 Chinese 24 Solar Terms (*Jieqi*) & Seasons — Eliminating `EPH-01` (`Sun().hlong` $180^\circ$ Bug)

#### ❌ Before (`ephem` — `Sun().hlong` Returns Heliocentric Earth Longitude $\lambda_\odot + 180^\circ$!)

```python
import ephem

# BUG EPH-01: ephem.Sun().hlong is the heliocentric longitude of EARTH (shifted by 180°
# and missing solar aberration & nutation), inverting Spring/Autumn and shifting terms by ~8 min!
s = ephem.Sun(ephem.Date("2026/03/20"))
wrong_lon_rad = s.hlong  # ~180° (Pi radians) at Vernal Equinox instead of 0°!
```

#### ✅ After (`lxephemeris` — Exact Geocentric Apparent Tropical Solar Longitude to $< 0.01\text{ s}$)

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")
jd_start = lxephemeris.julday(2026, 1, 1, 0.0)

# 1. All 24 Chinese Solar Terms (Jieqi, 15° multiples of apparent geocentric solar longitude)
solar_terms = engine.solar_terms_in_range(jd_start, jd_start + 365.25)
for term in solar_terms[:4]:
    print(f"{term.pinyin_name} ({term.english_name}): target={term.target_lon_deg}°, JD(TT)={term.jd_tt:.6f}")

# 2. Exact 4 Cardinal Seasons (Vernal Equinox 0°, Summer Solstice 90°, Autumnal 180°, Winter 270°)
seasons = engine.seasons_in_year(2026)
print("2026 Vernal Equinox JD(TT):", seasons.spring_equinox_jd_tt)
```

### 2.2 Topocentric Moon Ecliptic Coordinates — Eliminating `EPH-02` ($0.38^\circ$ J2000 Error)

#### ❌ Before (`ephem` — `ephem.Ecliptic(ephem.Moon(obs))` Mixes Topocentric Of-Date with J2000 Obliquity)

```python
import ephem

obs = ephem.Observer()
obs.lat, obs.lon, obs.elevation = "51.5074", "-0.1278", 35
obs.date = "2026/10/24 12:00:00"
moon = ephem.Moon(obs)
# BUG EPH-02: Rotates topocentric RA/Dec by J2000 obliquity -> up to 0.38° (23 arcmin) error!
ecl_moon = ephem.Ecliptic(moon)
```

#### ✅ After (`lxephemeris` — Exact WGS84 Topocentric Of-Date Frame Rotation)

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")
jd_ut = lxephemeris.julday(2026, 10, 24, 12.0)

topo_opts = lxephemeris.CalcOptions(
    origin="topocentric",
    observer=lxephemeris.Observer(51.5074, -0.1278, 35.0),
    plane="ecliptic",
    equinox="of_date",
    nutation="true",
    aberration="apparent",
)
moon_topo = engine.calc("Moon", jd_ut, scale="UT1", options=topo_opts)
print(f"Exact Topocentric Moon: lon={moon_topo.lon_deg:.6f}°, lat={moon_topo.lat_deg:.6f}°")
```

### 2.3 Almanac Sunrise/Sunset & Bernadette Brady Fixed-Star Parans — Eliminating `EPH-03`

#### ❌ Before (`ephem` — Static `FixedBody` Without 3D Proper Motion & Manual Paran Loops)

```python
import ephem

# BUG EPH-03: ephem.FixedBody has zero proper motion, parallax, or radial velocity!
# Fast stars like Sirius and Arcturus drift by arcminutes to degrees over centuries.
star = ephem.FixedBody()
star._ra, star._dec, star._epoch = "06:45:08.9", "-16:42:58", ephem.J2000
```

#### ✅ After (`lxephemeris` — Almanac `rise_trans` & Full 52-Star (+ Custom Catalogs) × 10-Planet `find_parans`)

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data")
jd_ut = lxephemeris.julday(2026, 10, 24, 0.0)

# 1. Standard Almanac Sunrise (upper limb + standard atmospheric refraction 1013.25 hPa, 15°C)
sunrise_jd = engine.rise_trans(
    "Sun", "rise", jd_ut, 51.5074, -0.1278, scale="UT1", disc_mode="upper_limb", pressure_hpa=1013.25
)

# 2. Bernadette Brady Diurnal Fixed-Star Parans (airless geometric center, 3D Hipparcos-2 space motion)
obs_london = lxephemeris.Observer(51.5074, -0.1278, 0.0)
parans = engine.find_parans(jd_ut, obs_london, orb_minutes=4.0)
for p in parans[:5]:
    print(f"{p.star_name} ({p.star_angle}) + {p.planet_name} ({p.planet_angle}): orb={p.time_orb_minutes:.2f} min")
```

---

## 3. Migrating from `Kerykeion` (`kerykeion`)

`lxephemeris` replaces Kerykeion's entire astronomical calculation layer while eliminating all 6 Kerykeion calculation defects (`KER-01..KER-06`).

### 3.1 Before vs. After (Python & Node.js — High-Latitude Charts, Seconds Precision, Boundary Guard, & 28-Phase Moon)

#### ❌ Before (`kerykeion` — Clamps Latitude at $\pm 66^\circ$, Drops Seconds, & Drops Cusp-Boundary Planets)

```python
from kerykeion import AstrologicalSubject

# BUG KER-01: Longyearbyen, Svalbard (78.22°N) is silently clamped to 66.0°N -> 7.85° ASC error!
# BUG KER-02: 0.0° lat/lon is silently mutated to 0.0001° (0.36" error).
# BUG KER-03: Birth seconds (e.g. 45.8s) cannot be passed -> up to 15 arcmin ASC error!
# BUG KER-04: Planets within 1e-4° of a house cusp drop into the previous house!
# BUG KER-05: Secondary points / Lilith / Vertex can remain Tropical in a Sidereal chart!
subject = AstrologicalSubject("Svalbard Native", 2026, 6, 21, 14, 30, lat=78.2232, lng=15.6267, tz_str="Arctic/Longyearbyen")
```

#### ✅ After (`lxephemeris` — Python)

```python
import lxephemeris

engine = lxephemeris.Engine.from_data_dir("/path/to/lxephemeris/data", verify_manifest_sha256=True)

# 1. Exact IANA timezone conversion preserving full fractional seconds (45.8s) and PEP 495 DST fold (Fixes KER-03 & EPH-05)
jd_utc = lxephemeris.utc_from_local(
    2026, 6, 21, 14, 30, 45.8, tz_name="Arctic/Longyearbyen", fold=0
)

# 2. Polar-safe houses at the TRUE geographic latitude 78.2232°N (Fixes KER-01 & KER-02)
houses = engine.houses(
    jd_utc,
    latitude_deg=78.2232,
    longitude_deg=15.6267,
    system="P",
    scale="UTC",
    polar_policy="fallback_to_whole_sign",
)
print(f"True Svalbard ASC={houses.asc_deg:.4f}° (effective system={houses.effective_system_code}, fallback={houses.polar_house_fallback})")

# 3. Boundary-safe house assignment with 1e-9° guard (Fixes KER-04)
sun = engine.calc("Sun", jd_utc, scale="UTC")
sun_house = lxephemeris.house_of_lon(sun.lon_deg, houses.cusps_deg, eps_deg=1e-9)

# 4. Complete 8-phase and Kerykeion 28-phase Lunar Phase Summary (Replaces kerykeion.calculate_moon_phase)
moon_phase = engine.lunar_phase_summary(jd_utc, scale="UTC")
print(
    f"Moon Phase: {moon_phase.phase_8_name} (8-idx={moon_phase.phase_8_index}, "
    f"Kerykeion 28-idx={moon_phase.phase_28_index}, illumination={moon_phase.illuminated_fraction:.4f})"
)
```

#### ✅ After (`lxephemeris` — Node.js / TypeScript)

```typescript
import { JsEngine, utcFromLocal, houseOfLon } from "lxephemeris";

const engine = JsEngine.fromDataDir("/path/to/lxephemeris/data", 4, null, true);

// 1. Exact IANA timezone + fractional seconds + DST fold
const jdUtc = utcFromLocal(2026, 6, 21, 14, 30, 45.8, "Arctic/Longyearbyen", 0);

// 2. Polar-safe houses at true latitude 78.2232°N
const houses = engine.houses(
  { jd1: jdUtc, jd2: 0 },
  { latitudeDeg: 78.2232, longitudeDeg: 15.6267, altitudeM: 0 },
  "P",
  "UTC",
  "fallback_to_whole_sign"
);

// 3. Boundary-safe house placement & 28-phase Moon summary
const sun = engine.calc("Sun", { jd1: jdUtc, jd2: 0 }, "UTC");
const sunHouse = houseOfLon(sun.lonDeg, houses.cuspsDeg, 1e-9);
const moonSummary = engine.lunarPhaseSummary(jdUtc);
console.log(`True Svalbard ASC=${houses.ascDeg.toFixed(4)}°, Sun in House ${sunHouse}, Moon 28-phase=${moonSummary.phase28Index}`);
```

---

## 4. Complete Side-by-Side API Crosswalk (`pyswisseph` / `ephem` / `kerykeion` $\to$ `lxephemeris`)

| Legacy Call (`pyswisseph` / `ephem` / `kerykeion`) | Python `lxephemeris.Engine` Call | Node.js `JsEngine` Call |
|---|---|---|
| `swe.set_ephe_path(path)` / Zero-arg init | `engine = lxephemeris.Engine()` or `Engine.from_data_dir(path)` | `const engine = new Engine()` or `Engine.fromDataDir(path)` |
| `swe.julday(y, m, d, h)` / `swe.revjul(jd)` | `lxephemeris.julday(y, m, d, h)` / `lxephemeris.revjul(jd)` | `julday(y, m, d, h)` / `revjul(jd)` |
| `ephem.Date(local_dt)` / `kerykeion` tz | `lxephemeris.utc_from_local(y, m, d, h, min, s, tz_name, fold=0)` | `utcFromLocal(y, m, d, h, min, s, tzName)` |
| `swe.deltat(jd_ut)` | `engine.delta_t(jd_ut, scale="UT1").delta_t_sec / 86400.0` | `engine.deltaT(jdUt, "UT1").deltaTSec / 86400.0` |
| `swe.sidtime(jd_ut)` | `engine.sidtime(jd_ut, scale="UT1").gast_hours` | `engine.sidtime(jdUt, "UT1").gastHours` |
| `swe.calc_ut(jd_ut, swe.ECL_NUT)` | `engine.ecl_nut(jd_ut, scale="UT1")` | `engine.eclNut(jdUt, "UT1")` |
| `swe.calc_ut(jd_ut, body_id, SEFLG_SWIEPH \| SEFLG_SPEED)` | `engine.calc(body_name, jd_ut, scale="UT1")` | `engine.calc(bodyName, jdUt, "UT1")` |
| `swe.calc_ut(jd_ut, body_id, SEFLG_MOSEPH \| SEFLG_SPEED)` | `Engine.new_moshier()` or `engine.calc_moshier(body, jd_ut, scale="UT1")` | `JsEngine.newMoshier()` or `engine.calcMoshier(body, jdUt, "UT1")` |
| `swe.calc_ut(jd_ut, body_id, SEFLG_EQUATORIAL \| SEFLG_SPEED)` | `engine.calc(body, jd_ut, scale="UT1", options=CalcOptions(plane="equatorial"))` | `engine.calc(body, jdUt, "UT1", { plane: "equatorial" })` |
| `swe.calc_ut(jd_ut, body_id, SEFLG_SIDEREAL)` | `engine.calc(body, jd_ut, scale="UT1", options=CalcOptions(zodiac="sidereal", ayanamsa="lahiri"))` | `engine.calc(body, jdUt, "UT1", { zodiac: "sidereal", ayanamsa: "lahiri" })` |
| `ephem.Ecliptic(ephem.Moon(obs))` | `engine.calc("Moon", jd_ut, scale="UT1", options=CalcOptions(origin="topocentric", observer=obs))` | `engine.calc("Moon", jdUt, "UT1", { origin: "topocentric", observer: obs })` |
| `swe.fixstar2_ut(star, jd_ut, flags)` | `engine.calc_star_by_name(star, jd_ut, scale="UT1")` | `engine.calcStarByName(star, jdUt, "UT1")` |
| `swe.houses_ex(jd_ut, lat, lon, hsys)` | `engine.houses(jd_ut, lat, lon, system=hsys, scale="UT1", polar_policy="fallback_to_whole_sign")` | `engine.houses(jdUt, obs, hsys, "UT1", "fallback_to_whole_sign")` |
| `swe.houses_armc(armc, lat, eps, hsys)` | `lxephemeris.houses_armc(armc, eps, lat, system=hsys)` | `housesArmc(armc, eps, lat, hsys)` |
| `kerykeion` house assignment | `lxephemeris.house_of_lon(lon_deg, houses.cusps_deg, eps_deg=1e-9)` | `houseOfLon(lonDeg, houses.cuspsDeg, 1e-9)` |
| `swe.house_pos(armc, lat, eps, hsys, [lon, lat])` | `engine.house_pos(armc, lat, eps, lon, lat_deg=lat, hsys=hsys)` | `engine.housePos(armc, lat, eps, lon, lat, hsys)` |
| `ephem.Observer.next_rising(Sun)` | `engine.rise_trans("Sun", "rise", jd_ut, lat, lon, scale="UT1", disc_mode="upper_limb", pressure_hpa=1013.25)` | `engine.riseTrans("Sun", "rise", jdUt, obs, "UT1", "upper_limb", 1013.25)` |
| `ephem` star parans (`FixedBody`) | `engine.find_parans(jd_ut, obs, orb_minutes=4.0)` | `engine.findParans(jdUt, obs, 4.0)` |
| `swe.heliacal_ut(...)` | `engine.heliacal_phenomena(star_or_planet, jd_ut, lat, lon, scale="UT1")` | `engine.heliacalPhenomena(starOrPlanet, jdUt, obs, "UT1")` |
| `swe.pheno_ut(jd_ut, body_id, flags)` | `engine.phenomena(body_name, jd_ut, scale="UT1")` | `engine.phenomena(bodyName, jdUt, "UT1")` |
| `swe.nod_aps_ut(jd_ut, body_id, flags, method)` | `engine.nod_aps(body_name, jd_ut, scale="UT1", method=method)` | `engine.nodAps(bodyName, jdUt, "UT1", method)` |
| `swe.sol_eclipse_when_glob` / `loc` | `engine.next_solar_eclipse` / `engine.next_local_solar_eclipse` / `engine.solar_eclipses_in_range` | `engine.nextSolarEclipse` / `engine.nextLocalSolarEclipse` / `engine.solarEclipsesInRange` |
| `swe.lun_eclipse_when` | `engine.next_lunar_eclipse` / `engine.lunar_eclipses_in_range` | `engine.nextLunarEclipse` / `engine.lunarEclipsesInRange` |
| `swe.lun_occult_when_glob` / `loc` | `engine.next_lunar_occultation` / `engine.next_local_lunar_occultation` | `engine.nextLunarOccultation` / `engine.nextLocalLunarOccultation` |
| `swe.time_equ(jd_ut)` | `engine.time_equation(jd_ut, scale="UT1")` | `engine.timeEquation(jdUt, "UT1")` |
| `swe.azalt_rev(jd_ut, flag, geopos, [az, alt])` | `engine.azalt_rev(jd_ut, flag, lat, lon, az, alt, scale="UT1")` | `engine.azaltRev(jdUt, flag, lat, lon, az, alt, "UT1")` |
| `ephem.Sun().hlong` Chinese Solar Terms | `engine.solar_terms_in_range(start_jd, end_jd)` | `engine.solarTermsInRange(startJd, endJd)` |
| `ephem.next_vernal_equinox` / `solstice` | `engine.seasons_in_year(year)` | `engine.seasonsInYear(year)` |
| Solar Return / Ingress / Retrograde 3-Pass | `engine.search_body_longitude` / `engine.body_longitude_crossings_in_range` | `engine.searchBodyLongitude` / `engine.bodyLongitudeCrossingsInRange` |
| Human Design `-88°` Solar Arc | `engine.human_design_epoch(birth_jd_tt)` | `engine.humanDesignEpoch(birthJdTt)` |
| `ephem.next_new_moon` / `next_full_moon` | `engine.search_lunar_phase` / `engine.lunar_phases_in_range` | `engine.searchLunarPhase` / `engine.lunarPhasesInRange` |
| `kerykeion.calculate_moon_phase` | `engine.lunar_phase_summary(jd_tt)` (`phase_8_index` + `phase_28_index`) | `engine.lunarPhaseSummary(jdTt)` (`phase8Index` + `phase28Index`) |
| Retrograde / Direct Station Search | `engine.search_stations_in_range(body, start_jd, end_jd)` | `engine.searchStationsInRange(body, startJd, endJd)` |
| Lunar Standstills / Declination Extrema | `engine.search_declination_extrema_in_range(body, start_jd, end_jd)` | `engine.searchDeclinationExtremaInRange(body, startJd, endJd)` |
| `swe.set_sid_mode(SE_SIDM_USER, t0, ayan_t0)` | `swe.set_sid_mode(SE_SIDM_USER, t0, ayan_t0)` | `swe.setSidMode(SE_SIDM_USER, t0, ayanT0)` |
| Custom `sefstars.txt` | `engine.load_custom_stars("/path/to/sefstars.txt")` | `engine.loadCustomStars("/path/to/sefstars.txt")` |
| Solar / Lunar / Saturn Returns (3-Pass) | `engine.planetary_return` / `engine.planetary_returns_in_range` | `engine.planetaryReturn` / `engine.planetaryReturnsInRange` |
| `kerykeion` Natal / Synastry / Transit Aspects | `engine.calc_chart_aspects` / `engine.calc_synastry_aspects` / `engine.natal_transits_in_range` | `engine.calcChartAspects` / `engine.calcSynastryAspects` / `engine.natalTransitsInRange` |
| Void-of-Course Moon | `engine.void_of_course_moon_in_range(start_jd, end_jd)` | `engine.voidOfCourseMoonInRange(startJd, endJd)` |
| Vedic Panchanga, Vimshottari Dasha & Vargas | `engine.panchanga` / `engine.vimshottari_dasha` / `engine.calc_shodashavarga` | `engine.panchanga` / `engine.vimshottariDasha` / `engine.calcShodashavarga` |
| Sect, Planetary Hours, Hermetic Lots & Snapshot | `engine.chart_sect` / `engine.planetary_hours` / `engine.calc_lots` / `engine.chart_snapshot` | `engine.chartSect` / `engine.planetaryHours` / `engine.calcLots` / `engine.chartSnapshot` |
