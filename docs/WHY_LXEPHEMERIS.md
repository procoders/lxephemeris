# Why We Built `lxephemeris` — Architecture, Performance & Eliminated Legacy Bugs

> **Production Proven**: `lxephemeris` was built as the high-precision astronomical and astrological computation engine powering **[AstroAPI (`https://astrology-api.io/`)](https://astrology-api.io/)** and **[Astrolium (`https://astrolium.com/`)](https://astrolium.com/)**. If you want a turnkey cloud REST API with built-in city geocoding, historical timezones, SVG chart rendering, and interpretations (without needing a self-hosted engine license), use **[AstroAPI](https://astrology-api.io/)** directly.

For nearly three decades, almost every astrology, horary, Vedic (*Jyotish*), Human Design, and astronomical backend in Python or Node.js has been forced to glue together **Swiss Ephemeris (`libswe` / `pyswisseph`)** — a C library designed in the 1990s for single-threaded desktop apps — alongside **`PyEphem`** and **`Kerykeion`**.

When running that legacy stack under high concurrency in modern production services (FastAPI, Node.js, Rust), five systemic problems emerge:

---

## 1. The 5 Pain Points of Legacy Ephemeris Stacks

1. **Global Mutable C State & Thread-Safety Bugs**: `libswe` stores ephemeris paths (`swe_set_ephe_path`), sidereal modes (`swe_set_sid_mode`), and topocentric observers (`swe_set_topo`) in process-global `static` C variables. If Thread A computes a Vedic (Lahiri) chart while Thread B computes a Western (Tropical) chart, **they silently corrupt each other's math**. Fixing this in Python required hundreds of lines of global `RLock` wrappers that serialized multi-core servers.
2. **A Mountain of Monkey-Patches in High-Level Wrappers**:
   - `Kerykeion` silently clamps latitudes above `±66°`, causing a **7.85° Ascendant error** in Northern Europe, Alaska, and Svalbard, and truncates birth-time seconds (`second=0`), shifting the Ascendant by up to **15 arcminutes**.
   - Rounding planet longitudes (`round(lon, 4)`) before house assignment routinely drops cusp planets into the wrong house.
   - `PyEphem`'s `Sun().hlong` returns Earth's heliocentric longitude shifted by `180°` (breaking Chinese *Jieqi* solar terms), and `ephem.Ecliptic(ephem.Moon(obs))` rotates topocentric coordinates using static J2000 obliquity (**23 arcminutes of error**).
3. **22 Upstream Calculation & C Buffer Bugs**:
   - Missing `.se1` ephemeris files cause `libswe` to **silently downgrade** to a lower-precision fallback without raising an error.
   - `swe_house_pos` returns an invalid `0.0` at high latitudes or circumpolar declinations.
   - `swe_fixstar2_ut("Mira")` matches `Mirach` (**>17° error**) due to naive prefix matching in C.
   - Star-anchored sidereal modes (`True Citra`, `Revati`, `Pushya`) wobble by **±20.5 arcseconds** every year because annual stellar aberration is mixed into the anchor star.
4. **Slow FFI Loops for Transits, Returns & Timelines**: Finding a year of exact transit windows, retrograde stations, Void-of-Course Moon intervals, or a 3-pass Saturn Return meant running tens of thousands of Python/JS bisection calls across the C FFI boundary — taking hundreds of milliseconds per request.
5. **Dependency Sprawl**: Teams had to maintain 5–6 different libraries — one for planetary positions, one for lunar phases, one for natal wheels, another for Vedic Panchanga/Dasha, another for Hellenistic Lots, and custom scripts for Human Design gates or Chinese calendars.

---

## 2. How `lxephemeris` Solves This

### One Engine Replaces 6+ Libraries (Batteries Included)
- **Solar System & 68 Minor Bodies**: `Sun..Pluto`, 7 Lunar points (`Mean_Node`, `True_Node`, `Mean_South_Node`, `True_South_Node`, `Mean_Lilith`, `True_Lilith`, and smooth $C^2$ `Interpolated_Lilith`), 9 Uranian hypotheticals (`Cupido..Poseidon`, `Isis`), and **68 JPL `SB441-N16` asteroids, centaurs, and dwarf planets** (`Chiron`, `Ceres`, `Pallas`, `Juno`, `Vesta`, `Pholus`, `Nessus`, `Eris`, `Sedna`, `Haumea`, `Makemake`, `Lilith (1181)`, `Proserpina (26)`, etc.).
- **23 House Systems + Polar Safety**: Placidus, Koch, Porphyry, Regiomontanus, Campanus, Equal, Whole Sign, Alcabitius, Morinus, Topocentric (Polich-Page), Krusinski, APC, Carter, Sunshine, Gauquelin (36 sectors), and more — with explicit `PolarPolicy` (`Strict`, `FallbackToWholeSign`, `FallbackToPorphyry`).
- **Native Event & Timing Solvers**:
  - **Aspects & Transits**: Static natal & synastry aspects (`calc_chart_aspects`, `calc_synastry_aspects`), declination parallels/contraparallels, exact sky-to-sky aspect windows (`mutual_aspects_in_range`), natal transit timelines (`natal_transits_in_range`), and **Void-of-Course Moon** intervals (`void_of_course_moon_in_range`).
  - **Planetary Returns & Stations**: Exact Solar, Lunar, Saturn, and Chiron returns (`planetary_return`, `planetary_returns_in_range`) that automatically detect all **3 passes** of a retrograde return loop, plus exact retrograde/direct stations (`search_stations_in_range`).
- **Vedic / Jyotish Suite**:
  - **43+1 Sidereal Ayanamsas**: All 43 standard models (`Lahiri`, `True Citra`, `Krishnamurti (KP)`, `Raman`, `Fagan-Bradley`, `Suryasiddhanta`, `Yukteshwar`, etc.) plus custom user-defined ayanamsas — free of annual aberration wobble.
  - **Panchanga, Muhurta & Vargas**: Full 5-limb `panchanga` (`Tithi`, `Nakshatra + Pada`, `Yoga`, `Karana`, `Vara`, plus `Rahukalam`, `Yamagandam`, `Gulika Kalam`, `Abhijit Muhurta`), all **16 Parashara Shodashavarga (`D1..D60`)**, and **120-year Vimshottari Dasha** (`Maha`, `Antar`, `Pratyantar`).
- **Traditional & Hellenistic Suite**: Day/Night chart sect (`chart_sect`), **24 Chaldean Planetary Hours** (`planetary_hours`), and **7 Hermetic Lots** (`Fortune`, `Spirit`, `Eros`, `Necessity`, `Courage`, `Victory`, `Nemesis` with diurnal/nocturnal reversal).
- **Human Design & Chinese Solar Calendar**:
  - **Human Design**: Exact `-88°` prenatal Design solar arc solver (`human_design_epoch`) and 64 I-Ching Mandala Gates × 6 Lines × Color × Tone × Base (`map_zodiac_placement`).
  - **Chinese Almanac**: All 24 Solar Terms (*Jieqi* / *Zhongqi*, `solar_terms_in_range`), 4 cardinal seasons (`seasons_in_year`), and lunar phase solvers (`lunar_phases_in_range`, 8-phase & 28-phase `lunar_phase_summary`).
- **Fixed Stars, Brady Parans, Eclipses & Occultations**:
  - **52 Embedded Fixed Stars** (VizieR Hipparcos-2 `I/311` + XHIP `V/137D` with full 3D space motion) + custom `sefstars.txt` / JSON catalog loader.
  - **Bernadette Brady Parans** (`find_parans`) & **Schaefer Heliacal Risings/Settings**.
  - **Solar & Lunar Eclipses** (global & local Besselian solvers with Saros series/member and path width) & **Lunar Occultations** of stars, planets, and asteroids.

### 5–15× Faster via Single-Pass Analytical Velocities
Legacy ephemeris libraries compute planet speed (`°/day`) by evaluating the ephemeris 3 separate times ($t - \Delta t, t, t + \Delta t$) and taking a finite difference — which is both 3× slower and prone to $0^\circ / 360^\circ$ wrap-around glitches.
- `lxephemeris` differentiates Chebyshev polynomials, relativistic light-time, aberration, and the IAU 2006/2000A frame rotation matrix **analytically in a single pass**:
  - **`12.8 µs`** per planetary position + instantaneous velocity.
  - **`~150–240 µs`** for a complete natal chart (`chart_snapshot`).
  - Releases the **Python GIL** and runs multi-core **Rayon** batch queries (`calc_batch`, `calc_batch_flat_f64`, `houses_batch`) across all CPU cores.

---

## 3. Summary of the 22 Eliminated Upstream Bugs

For the full mathematical equations and regression tests, see **[`docs/SE_DIVERGENCES.md`](SE_DIVERGENCES.md)**.

### Swiss Ephemeris (`SWE-01..SWE-11`)
- **`SWE-01` (Global State Race Conditions)**: `set_ephe_path`, `set_sid_mode`, and `set_topo` mutate process-global C state. In `lxephemeris`, all configuration is immutable or passed per call (`Send + Sync`).
- **`SWE-02` (Star-Anchored Ayanamsa Wobble)**: `True Citra`, `Revati`, and `Pushya` in `libswe` include annual stellar aberration ($\pm 20.5''$ yearly oscillation). `lxephemeris` evaluates star anchors strictly in the geometric frame.
- **`SWE-03` (`swe_house_pos` Returning `0.0`)**: `libswe` returns an invalid `0.0` house position at circumpolar declinations or high latitudes. `lxephemeris` uses a closed-form 3D horizon/prime-vertical projection valid on $[1.0, 13.0)$.
- **`SWE-04` (`Mira` / `Mirach` Star Collision)**: `swe_fixstar2_ut("Mira")` prefix-matches `Mirach` ($>17^\circ$ error). `lxephemeris` uses exact normalized name, Bayer, and Hipparcos ID matching.
- **`SWE-05` (Sidereal Velocity Double-Nutation)**: Fixes frame-rate and nutation-rate inconsistencies in sidereal speeds.
- **`SWE-06` (Silent Ephemeris Fallback)**: Missing `.se1` files silently downgrade accuracy in `libswe`. `lxephemeris` validates SHA-256 manifests explicitly and exposes `engine.is_moshier_mode()`.
- **`SWE-07` (Silent Polar House Flips)**: Replaces silent Porphyry substitution and $180^\circ$ Ascendant flips with explicit `PolarPolicy` (`Strict`, `FallbackToWholeSign`, `FallbackToPorphyry`).
- **`SWE-08` (Zero/Negative Stellar Parallax Explosion)**: Guards Hipparcos-2 parallax $\varpi \le 0$ against division-by-zero and $180^\circ$ sky inversions.
- **`SWE-09` ($0^\circ / 360^\circ$ Velocity Glitch)**: Single-pass analytical Jacobians eliminate finite-difference branch-cut spikes at $0^\circ \ Aries$.
- **`SWE-10` (Besselian Eclipse Aberration Mixing)**: Separates solar annual aberration from lunar shadow-plane coordinates (`0.2–0.6 s` contact improvement).
- **`SWE-11` (`SE_MEAN_APOG` Oscillation)**: Provides true secular IERS 2003 `Mean_Lilith` alongside osculating `True_Lilith` and smooth $C^2$ `Interpolated_Lilith`.

### PyEphem (`EPH-01..EPH-05`)
- **`EPH-01` (`Sun().hlong` $180^\circ$ Inversion)**: Fixes the heliocentric Earth / geocentric Sun confusion that breaks Chinese *Jieqi* solar terms.
- **`EPH-02` (Topocentric Moon $0.38^\circ$ Obliquity Error)**: Rotates topocentric coordinates using true obliquity of date $\varepsilon(t)$ instead of static J2000 obliquity.
- **`EPH-03` (Static `FixedBody` Stars)**: Applies full 3D Hipparcos-2 proper motion, parallax, and radial velocity in star parans.
- **`EPH-04` (Truncated VSOP87 / ELP2000 Drift)**: Replaces `libastro` series with JPL `DE440s` and `DE440`-calibrated series.
- **`EPH-05` (Naive Datetime & DST Fold Loss)**: `utc_from_local` enforces explicit IANA timezone + PEP 495 `fold` disambiguation.

### Kerykeion (`KER-01..KER-06`)
- **`KER-01` ($\pm 66^\circ$ Latitude Clamping)**: Eliminates the silent latitude clamp that caused a **$7.85^\circ$ Ascendant error** at Svalbard ($78.22^\circ\text{N}$).
- **`KER-02` ($0.0 \to 0.0001^\circ$ Equator/Greenwich Hack)**: Handles $(0^\circ, 0^\circ)$ exactly with zero coordinate falsification.
- **`KER-03` (Birth-Time Seconds Truncation)**: Preserves full `f64` sub-second precision (preventing up to $15'$ of Ascendant drift).
- **`KER-04` (Cusp Rounding House Drop)**: `house_of_lon` uses a formal `1e-9°` boundary epsilon instead of `round(lon, 4)`.
- **`KER-05` (Mixed Tropical Points in Sidereal Charts)**: Applies sidereal ayanamsa uniformly across all planets, minor bodies, nodes, Liliths, lots, and fixed stars.
- **`KER-06` (Hardcoded $23.4367^\circ$ Obliquity)**: Uses exact IAU 2006/2000A true obliquity and GAST everywhere.
