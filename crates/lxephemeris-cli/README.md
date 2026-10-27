# `lxephemeris-cli` — Command-Line Interface for `lxephemeris`

[![Crates.io](https://img.shields.io/crates/v/lxephemeris-cli.svg)](https://crates.io/crates/lxephemeris-cli)
[![Rust](https://img.shields.io/badge/Rust-1.80%2B%20%7C%20forbid(unsafe__code)-orange.svg)](https://github.com/procoders/lxephemeris)
[![License](https://img.shields.io/badge/License-PolyForm_Noncommercial_%7C_Commercial-purple.svg)](https://github.com/procoders/lxephemeris/blob/main/LICENSE.md)

`lxephemeris-cli` is the standalone command-line tool for the **[`lxephemeris`](https://github.com/procoders/lxephemeris)** astronomical and astrological engine.

It works **out of the box with zero external files** (using the built-in `DE440s`-calibrated analytical engine covering `-3000 CE .. +3000 CE`), or with high-precision NASA JPL `DE440s` / `DE440` / `DE441` `.bsp` SPK kernels (`< 0.003"` precision).

---

## Installation

```bash
cargo install lxephemeris-cli
```

---

## Commands & Examples

### 1. `calc` — Planetary, Asteroid, Lunar & Fixed Star Positions
Calculate apparent geocentric, topocentric, or heliocentric coordinates, zodiac placements, speeds, and Human Design gates:

```bash
# Instant calculation (zero external files required)
lxephemeris-cli calc --jd 2460676.5 --bodies Sun,Moon,Mercury,Venus,Mars,Jupiter,Saturn,Chiron,True_Node,Mean_Lilith

# Sidereal (Vedic Lahiri) calculation with JSON output
lxephemeris-cli calc --jd 2460676.5 --bodies Sun,Moon,Mars --zodiac sidereal --ayanamsa lahiri --json
```

### 2. `houses` — 23 Astrological House Systems & Chart Angles
Compute 12 house cusps and 8 angles (`ASC`, `MC`, `DSC`, `IC`, `Vertex`, `Equatorial ASC`, `Co-ASC`, `Polar ASC`) with explicit polar safety:

```bash
# Placidus houses for London
lxephemeris-cli houses --jd 2460676.5 --lat 51.5074 --lon -0.1278 --system P

# Polar-safe calculation for Tromsø / Svalbard with Whole Sign fallback
lxephemeris-cli houses --jd 2460676.5 --lat 78.2232 --lon 15.6267 --system P --polar-policy fallback_to_whole_sign
```

### 3. `panchanga` — Vedic Panchanga & Muhurta Timings
Compute `Tithi`, `Nakshatra + Pada`, `Yoga`, `Karana`, `Vara`, sunrise/sunset, `Rahukalam`, `Yamagandam`, `Gulika Kalam`, and `Abhijit Muhurta`:

```bash
lxephemeris-cli panchanga --jd 2460676.5 --lat 28.6139 --lon 77.2090 --ayanamsa lahiri
```

### 4. `swetest` — Legacy `swetest`-Compatible Output
Drop-in CLI compatibility for scripts expecting Swiss Ephemeris `swetest` flags:

```bash
lxephemeris-cli swetest -b1.1.2000 -p0123456789 -fPlbR
```

### 5. `fetch-kernel`, `data-sync` & `validate-bundle` — Ephemeris & IERS Data Management
Download official NASA JPL / IMCCE SPK kernels, synchronize IERS `finals2000A` ($\Delta T$ / `UT1-UTC`) and leap-second tables, and verify SHA-256 manifest checksums:

```bash
# Fetch NASA JPL DE440s (1849–2150), DE440 (1550–2650), DE441 (-13200..+17191), or INPOP21a
lxephemeris-cli fetch-kernel --kernel de440s --data-dir ./data

# Synchronize IERS finals2000A Earth orientation & leap-second tables atomically
lxephemeris-cli data-sync --data-dir ./data

# Verify SHA-256 integrity of ./data/manifest.json and inspect SPK segments
lxephemeris-cli validate-bundle --data-dir ./data
```

### 6. `chart-bench` — Latency Benchmark
Benchmark single-thread and multi-thread full natal chart calculation speed on your hardware:

```bash
lxephemeris-cli chart-bench --data-dir ./data --iterations 1000
```

### 7. `license-verify` — Offline Commercial License Verification
Verify an offline `Ed25519`-signed commercial license key (`LXEPHEMERIS_LICENSE_KEY` or `LXEPHEMERIS_LICENSE_FILE`) in CI/CD or production containers:

```bash
lxephemeris-cli license-verify
lxephemeris-cli license-verify --key "lxeph_live_v1...." --require-commercial
```

---

## Licensing & Cloud API

Distributed under a Dual License (**`PolyForm-Noncommercial-1.0.0`** for non-commercial use / **`LicenseRef-Commercial`** for self-hosted commercial production).
- **Turnkey Cloud API**: Need ready-to-use REST endpoints with city geocoding, historical timezones, SVG chart wheels & interpretations? Use **[AstroAPI (`https://astrology-api.io/`)](https://astrology-api.io/)** (powered by `lxephemeris`, commercial rights included).
- **Self-Hosted Commercial License**: Request at **[https://procoders.tech/contacts/](https://procoders.tech/contacts/)** or **[lxephemeris@procoders.tech](mailto:lxephemeris@procoders.tech)**.
- **Full Documentation**: **[https://github.com/procoders/lxephemeris](https://github.com/procoders/lxephemeris)**
