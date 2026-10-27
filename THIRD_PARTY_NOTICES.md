# Third-Party Notices & Clean-Room Provenance Record

`lxephemeris` is engineered under strict clean-room rules. No source code, coefficient tables, or files from Swiss Ephemeris (`swisseph`/`pyswisseph`), `Kerykeion`, `libephemeris`, `Anonyfox/celestine`, or any GPL/AGPL-licensed software have been read, translated, copied, or linked into this repository.

This file documents the provenance audit (Phase `0.0`) and reproduces all required copyright, attribution, and license notices for third-party dependencies, algorithmic references, and astronomical datasets used by `lxephemeris`.

---

## 1. Provenance Gate Audit Summary (Phase `0.0`)

| Component / Reference | Version / Commit | License | Clean-Room Audit Status |
|---|---|---|---|
| `nyx-space/anise` | `0.9.x` | `MPL-2.0` | **Approved** (Crate dependency; unmodified MPL-2.0 library isolated behind internal `SpkSource` trait) |
| `nyx-space/hifitime` | `4.x` | `MPL-2.0` | **Approved** (Crate dependency; nanosecond TAI/TT/UTC/UT1/TDB epoch arithmetic) |
| `astro-xao/sofars` | `0.6.x` | `MIT` + IAU SOFA License | **Approved** (Pure-Rust port of IAU SOFA C library without `unsafe`; mandatory 6-clause IAU SOFA notice included below) |
| `liberfa/erfa` | `2.0.x` | `BSD-3-Clause` | **Approved** (Algorithmic reference & Python `pyerfa` Level 1B golden vector oracle) |
| `jayeshmepani/jpl-ephemeris` (JME) | MIT | `MIT` | **Audited & Approved as Architectural Reference Only** — All mathematical constants in `lxephemeris` are sourced directly from primary peer-reviewed literature (*Lieske et al. 1977*, *Simon et al. 1994*, *Holden 1977*, *Munkasey 1991*, *Espenak & Meeus 2006*, *IAU 2000A/2006*). No unverified constant tables are copied from JME. |
| `skyfielders/python-skyfield` | `1.x` | `MIT` | **Approved** (Level 1B golden vector verification oracle & astrometric pipeline reference) |
| `cosinekitty/astronomy` | `2.x` | `MIT` | **Approved** (Algorithmic reference for rise/set/culmination & eclipse root-finding) |
| `TheDaniel166/moira` | MIT | `MIT` | **Approved** (Reference for house geometry cross-checks alongside *Munkasey 1991* & *Holden 1977*) |
| NASA JPL `DE440s` & Horizons SPK | DE440 / SB441-N16 | US Gov Public Work / JPL Terms | **Approved** (Binary DAF/SPK ephemeris kernels) |
| CDS VizieR `I/311`, `I/239`, `V/137D` | van Leeuwen 2007 / ESA 1997 / Anderson & Francis 2012 | CDS / Academic Open Access | **Approved** (50 fixed stars astrometric cross-match with attribution below) |
| IERS EOP `finals2000A` & `leapseconds` | IERS / USNO | Public Scientific Data | **Approved** (`UT1-UTC` and TAI-UTC leap second tables) |

---

## 2. `astro-xao/sofars` and IAU SOFA Software License Notice

`lxephemeris-core` depends on [`sofars`](https://github.com/astro-xao/sofars) (`MIT` license), a pure-Rust translation of the International Astronomical Union (IAU) Standards of Fundamental Astronomy (SOFA) library. In compliance with the IAU SOFA Software License (clauses 3(a)–3(e)), the following notice applies to all routines in `sofars` and any algorithms in `lxephemeris` derived from IAU SOFA:

```text
Copyright (C) 2023-2026 astro-xao contributors (sofars Rust port, MIT License)
Copyright (C) International Astronomical Union Standards of Fundamental Astronomy (IAU SOFA)

NOTICE TO USER:

1. This software incorporates and/or is derived from software provided by the
   International Astronomical Union's SOFA (Standards of Fundamental Astronomy)
   initiative (http://www.iausofa.org).
2. The algorithms and source code in `sofars` and `lxephemeris` differ from the
   original ANSI C / Fortran 77 code distributed by the IAU SOFA Center (specifically,
   they are implemented in safe Rust).
3. Neither the IAU SOFA Center nor the IAU SOFA Board has endorsed, reviewed, or
   certified `sofars` or `lxephemeris`.
4. These routines may not be redistributed under the pristine `iau*` or `sofa*`
   C/Fortran routine names in a manner that could be confused with official SOFA
   distributions.

SOFA Software License (Standard 6-Clause Terms):

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this
   list of conditions and the following disclaimer in the documentation and/or
   other materials provided with the distribution.
3. Neither the name of the International Astronomical Union nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.
4. Any modified version of SOFA software must carry a prominent notice stating
   (a) that it is derived from software provided by IAU SOFA, (b) the nature of
   the modifications, and (c) that it is not endorsed by the IAU SOFA Board.
5. Users are requested to cite the IAU SOFA service in publications benefiting
   from the use of SOFA software.
6. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

---

## 3. `nyx-space/anise` and `nyx-space/hifitime` (MPL-2.0)

`lxephemeris` uses [`anise`](https://github.com/nyx-space/anise) and [`hifitime`](https://github.com/nyx-space/hifitime) as unmodified Cargo crate dependencies under the **Mozilla Public License 2.0 (MPL-2.0)**.

- **Copyright**: Copyright (c) Christopher Rabotin and Nyx Space contributors.
- **License**: Mozilla Public License Version 2.0 (`https://www.mozilla.org/en-US/MPL/2.0/`).
- **Source Code Availability**: Both `anise` and `hifitime` are linked as unmodified upstream crates from `crates.io` (`https://crates.io/crates/anise`, `https://crates.io/crates/hifitime`), and their complete source code is publicly available at `https://github.com/nyx-space/anise` and `https://github.com/nyx-space/hifitime`.

---

## 4. `liberfa/erfa` (BSD-3-Clause)

```text
Copyright (C) 2013-2023, NumFOCUS Foundation.
Derived, with permission, from the SOFA library.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the Standards Of Fundamental Astronomy (SOFA),
   the International Astronomical Union (IAU), nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.
```

---

## 5. MIT-Licensed Algorithmic & Verification References

The following open-source libraries are acknowledged as clean-room architectural and verification references under the MIT License:

- **`jayeshmepani/jpl-ephemeris` (JME)** — Copyright (c) Jayesh Mepani (`https://github.com/jayeshmepani/jpl-ephemeris`, MIT License).
- **`skyfielders/python-skyfield`** — Copyright (c) Brandon Rhodes (`https://github.com/skyfielders/python-skyfield`, MIT License).
- **`cosinekitty/astronomy` (Astronomy Engine)** — Copyright (c) Don Cross (`https://github.com/cosinekitty/astronomy`, MIT License).
- **`TheDaniel166/moira`** — Copyright (c) Daniel (`https://github.com/TheDaniel166/moira`, MIT License).

```text
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
```

---

## 6. Astronomical Ephemerides, Star Catalogs, and Earth Orientation Data

### 6.1. NASA JPL `DE440s` and JPL Horizons Small-Body SPK Kernels
- **Source**: Jet Propulsion Laboratory (JPL), California Institute of Technology, under contract with the National Aeronautics and Space Administration (NASA).
- **References**:
  - Park, R. S., Folkner, W. M., Williams, J. G., & Boggs, D. H. (2021). *The JPL Planetary and Lunar Ephemerides DE440 and DE441*. The Astronomical Journal, 161(3), 105.
  - JPL Solar System Dynamics Group, Horizons On-Line Ephemeris System (`https://ssd.jpl.nasa.gov/horizons/`), Navigation and Ancillary Information Facility (`https://naif.jpl.nasa.gov/`).

### 6.2. CDS VizieR Star Catalogs (`I/311`, `I/239`, `V/137D`)
This work has made use of the VizieR catalogue access tool, CDS, Strasbourg, France (DOI: `10.26093/cds/vizier`). The original description of the VizieR service was published in *A&AS 143, 23 (2000)*.

1. **VizieR `I/311` (Hipparcos, the New Reduction)**:
   - van Leeuwen, F. (2007). *Validation of the new Hipparcos reduction*. Astronomy & Astrophysics, 474(2), 653–664.
2. **VizieR `I/239` (The Hipparcos and Tycho Catalogues)**:
   - European Space Agency (ESA) (1997). *The Hipparcos and Tycho Catalogues*, ESA SP-1200.
3. **VizieR `V/137D` (Extended Hipparcos Compilation — XHIP)**:
   - Anderson, E., & Francis, Ch. (2012). *XHIP: An extended hipparcos compilation*. Astronomy Letters, 38(5), 331–346.

### 6.3. IERS Earth Orientation Parameters & Historical $\Delta T$
- **IERS / USNO `finals2000A` and Bulletin C**: International Earth Rotation and Reference Systems Service (IERS) Earth Orientation Center / U.S. Naval Observatory (USNO).
- **Historical $\Delta T$ Series**: Stephenson, F. R., Morrison, L. V., & Hohenkerk, C. Y. (2016). *Measurement of the Earth's rotation: 720 BC to AD 2015*. Proceedings of the Royal Society A, 472(2196), 20160404.
