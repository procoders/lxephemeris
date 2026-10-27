# `lxephemeris` — Mathematical Model Passports (`MODEL_PASSPORTS.md`)

**Status:** Phase 0.3 Clean-Room Mathematical Specification (Verified)  
**Target Crate:** `lxephemeris-core`  
**Permissive Reference Sources:**
- IAU SOFA / ERFA (`liberfa/erfa`, BSD-3; `astro-xao/sofars`, MIT + IAU SOFA License)
- Lieske, Lederle, Fricke & Morando (1977), *A&A* 58, 1–16 (Table 5 two-epoch precession)
- *The Indian Astronomical Ephemeris 1989* (Controller of Publications, Delhi, pp. 525, 556)
- Simon, Bretagnon, Chapront, Chapront-Touzé, Francou & Laskar (1994), *A&A* 282, 663–683
- IERS Conventions (2003), Technical Note 32, §5.7.2 (Fundamental Arguments `fal03`, `faf03`, `faom03`)
- Park, Folkner, Williams & Boggs (2021), *AJ* 161, 105 (JPL DE440/DE441 constants)
- Neely, James (1980), "New Ephemerides for the Uranian Planets," *Matrix Journal* VII
- Sevin, Émile (1946), *C. R. Acad. Sci. Paris* / Strubell, W. (1952), *Die Sterne* 3/1952 (`Transpluto` / `Isis`)
- CDS VizieR `I/311` (Hipparcos-2, van Leeuwen 2007, *A&A* 474, 653), `V/137D` (XHIP, Anderson & Francis 2012, *Astron. Lett.* 38, 331), and `V/50` (Bright Star Catalogue 5th Ed., Hoffleit & Warren 1991)
- Permissive reference implementations: `jayeshmepani/jpl-ephemeris` (MIT) and `TheDaniel166/moira` (MIT)

---

## 1. The 12 Ayanamsas (`AyanamsaModel`)

### 1.1 Fundamental Frame & Nutation Rules

For any evaluation epoch $t$ given as Terrestrial Time Julian Date $\text{JD}_{\text{TT}}(t)$:

1. **Mean Ayanamsa $A_{\text{mean}}(t)$**:
   The longitude of the sidereal origin (`0° Aries sidereal`) measured eastward along the **mean ecliptic of date** from the **mean equinox of date** $\Upsilon_{\text{mean}}(t)$, normalized to $[0^\circ, 360^\circ)$.
2. **True Ayanamsa $A_{\text{true}}(t)$**:
   The longitude of the sidereal origin measured from the **true equinox of date** $\Upsilon_{\text{true}}(t)$:
   $$A_{\text{true}}(t) = \bigl(A_{\text{mean}}(t) + \Delta\psi(t)\bigr) \bmod 360^\circ$$
   where $\Delta\psi(t)$ is the IAU 2000A nutation in longitude (`sofars::pnm::nut00a`) converted from radians to degrees.
3. **Tropical-to-Sidereal Longitude Invariance**:
   - If the tropical longitude $\lambda_{\text{trop, true}}(t)$ is referred to the **true equinox of date** (`NutationMode::True`, the default for `Frame::Ecliptic`), subtract $A_{\text{true}}(t)$:
     $$\lambda_{\text{sid}}(t) = \bigl(\lambda_{\text{trop, true}}(t) - A_{\text{true}}(t)\bigr) \bmod 360^\circ$$
   - If the tropical longitude $\lambda_{\text{trop, mean}}(t)$ is referred to the **mean equinox of date** (`NutationMode::Mean`), subtract $A_{\text{mean}}(t)$:
     $$\lambda_{\text{sid}}(t) = \bigl(\lambda_{\text{trop, mean}}(t) - A_{\text{mean}}(t)\bigr) \bmod 360^\circ$$
   Because $\lambda_{\text{trop, true}}(t) = \lambda_{\text{trop, mean}}(t) + \Delta\psi(t)$ to first order, $\Delta\psi(t)$ cancels in the sidereal longitude, keeping the sidereal zodiac inertially anchored rather than oscillating by $\pm 17''$ with lunar node nutation.

---

### 1.2 Category A: 9 Epoch-Defined Models (Lieske 1977 Two-Epoch Precession)

Historical sidereal traditions define a fixed anchor value $A_{\text{mean}}(t_0)$ at a reference epoch $t_0$ ($\text{JD}_{\text{TT}}(t_0)$) and propagate it via the IAU 1976 two-epoch general precession in longitude $p_A(t_0 \to t)$ (Lieske et al. 1977, *A&A* 58, 1–16, Table 5; *Indian Astronomical Ephemeris 1989*, p. 525):

$$\begin{aligned}
T &= \frac{\text{JD}_{\text{TT}}(t_0) - 2451545.0}{36525.0} \quad (\text{Julian centuries from } J2000.0 \text{ to } t_0) \\
t &= \frac{\text{JD}_{\text{TT}} - \text{JD}_{\text{TT}}(t_0)}{36525.0} \quad (\text{Julian centuries from } t_0 \text{ to target } \text{JD}_{\text{TT}})
\end{aligned}$$

$$p_A(t_0 \to t) = \bigl(5029.0966'' + 2.22226''\,T - 0.000042''\,T^2\bigr)\,t + \bigl(1.11161'' - 0.000127''\,T\bigr)\,t^2 - 0.000113''\,t^3 \quad [\text{arcsec}]$$

The mean ayanamsa at $\text{JD}_{\text{TT}}$ is:

$$A_{\text{mean}}(\text{JD}_{\text{TT}}) = \left(A_{\text{mean}}(t_0) + \frac{p_A(t_0 \to t)}{3600.0}\right) \bmod 360^\circ$$

The analytic secular rate $\dot{A}_{\text{mean}}(\text{JD}_{\text{TT}})$ in **degrees per TT day** (`deg/day_TT`) is:

$$\dot{A}_{\text{mean}}(\text{JD}_{\text{TT}}) = \frac{\bigl(5029.0966 + 2.22226\,T - 0.000042\,T^2\bigr) + 2\bigl(1.11161 - 0.000127\,T\bigr)\,t - 3(0.000113)\,t^2}{3600.0 \times 36525.0}$$

*(Note: Even when `NutationMode::True` is requested for `AyanamsaResult`, `rate_deg_day` reports the secular rate $\dot{A}_{\text{mean}}$ unless the caller explicitly differentiates $A_{\text{true}}(t)$).*

#### Complete Passport Table for the 9 Epoch-Defined Models

| # | `astra` Key | Legacy ID | Enum Variant | Reference Epoch $t_0$ ($\text{JD}_{\text{TT}}$) | Calendar Instant ($t_0$) | Mean Anchor $A_{\text{mean}}(t_0)$ (DMS) | Mean Anchor $A_{\text{mean}}(t_0)$ (deg) | Mean Value at $J2000.0$ ($\text{JD}_{\text{TT}} = 2451545.0$) | Primary Literature & Provenance Notes |
|---:|---|---:|---|---:|---|---|---:|---:|---|
| 1 | `fagan_bradley` | `0` | `FaganBradley` | `2433282.42345905` | `B1950.0` (`1949-12-31 22:09:46.86 TT`) | `24°02'31.3600"` (`24°02'30.9474"` in legacy SE) | `24.042044444444444` (`24.041929833333333` SE) | `24.7404145675` (`24.7402999564` SE) | Fagan & Firebrace, *A Primer of Sidereal Astrology* (1971): "synetic vernal point" with Spica at $29^\circ 06' 05''$ Virgo (no aberration) at `1950.0` (`B1950.0`). |
| 2 | `lahiri` | `1` | `Lahiri` | `2435553.50000000` | `1956-03-21 00:00:00 TT` | `23°14'43.8891"` (from $A_{\text{true}}(t_0) = 23^\circ 15' 00.658'' - \Delta\psi_{1980}(t_0)$) | `23.245524741262587` (`23.2455609515` SE) | `23.8570561072` (`23.8570923175` SE) | *Indian Astronomical Ephemeris 1989*, p. 556 footnote: $A_{\text{true}}(t_0) = 23^\circ 15' 00.658''$ at `1956-03-21 0h ET` with $\Delta\psi_{1980}(t_0) = +16.76893''$ (`eraNut80`). |
| 3 | `raman` | `3` | `Raman` | `2415020.00000000` | `J1900.0` (`1899-12-31 12:00:00 TT`) | `21°00'51.1560"` | `21.014210000000000` | `22.4107910119` | B. V. Raman, *Hindu Predictive Astrology* (1938), pp. 378–379 (`21.01421°` at `J1900.0`; zero year 397 CE). |
| 4 | `ushashashi` | `4` | `Ushashashi` | `2415020.00000000` | `J1900.0` (`1899-12-31 12:00:00 TT`) | `18°39'39.4560"` | `18.660960000000000` | `20.0575410119` | Ushā-Shashī Revati-based zodiac at `1900.0` (Mercier 1976/1977). |
| 5 | `krishnamurti` | `5` | `Krishnamurti` | `2415020.00000000` | `J1900.0` (`1899-12-31 12:00:00 TT`) | `22°21'49.1724"` | `22.363659000000000` | `23.7602400119` | K. S. Krishnamurti, *KP Reader I* (1963), pp. 55–59 (`22.363659°` at `J1900.0`). |
| 6 | `yukteshwar` | `7` | `Yukteshwar` | `2415020.00000000` | `J1900.0` (`1899-12-31 12:00:00 TT`) | `21°04'56.0000"` | `21.082222222222222` | `22.4788032341` | Sri Yukteswar, *The Holy Science* (1894): zero in 499 CE (`21°04'56"` at `J1900.0`). |
| 7 | `jn_bhasin` | `8` | `JnBhasin` | `2415020.00000000` | `J1900.0` (`1899-12-31 12:00:00 TT`) | `21°21'56.0000"` | `21.365555555555556` | `22.7621365675` | J. N. Bhasin sidereal tables (`21°21'56"` at `J1900.0`; zero at $\text{JD}_{\text{TT}} \approx 1854239.3$). |
| 8 | `suryasiddhanta` | `21` | `Suryasiddhanta` | `1903396.87706980` | `0499-03-21 07:30:31.57 UT` (Julian; Kali 3600 midnight) | `0°00'00.0000"` | `0.000000000000000` | `20.8950588452` | E. Burgess, *Translation of the Sūrya-Siddhānta* (1860): Kaliyuga 3600 mean vernal equinox ($\text{JD}_{\text{UT}} = 1903396.8128654$, $\Delta T \approx +5547.26\text{ s}$). |
| 9 | `ss_citra` | `26` | `SsCitra` | `1903396.87706980` | `0499-03-21 07:30:31.57 UT` (Julian; Kali 3600 midnight) | `+2°06'38.5360"` | `+2.110704444444444` | `23.0057632896` | Surya Siddhanta polar longitude of Citra (Spica) fixed at $180^\circ$ at the 499 CE Kali 3600 epoch ($A_{\text{mean}}(t_0) = +2.11070444^\circ$). |

> **Implementation Note on Canonical `J2000.0` Form vs Two-Epoch Form:**
> Because `Lieske (1977)` Table 5 is a cubic polynomial on the tangent bundle of the 1-parameter precession group, $p_A(t_0 \to t) - [p_A(t_0 \to J2000) + p_A(J2000 \to t)]$ is at most $0.0002''$ over `1850–2150` and $0.08''$ back to `499 CE`. Storing the exact $(t_0, A_{\text{mean}}(t_0))$ pairs from the table above and evaluating `iae_general_precession_arcsec(t0, jd_tt)` reproduces the historical epoch definitions to `< 0.0001"`.

---

### 1.3 Category B: 3 Star-Anchored Models (`true_citra`, `true_revati`, `true_pushya`)

A star-anchored sidereal model locks a specific fiducial star $S$ to a constant sidereal ecliptic longitude $\lambda_{\text{target}}$ at every epoch $t$:

#### Exact Mathematical Algorithm (`StarAnchored`)
Given target epoch $\text{JD}_{\text{TT}}$:
1. **Space-Motion Propagation (`J1991.25` $\to \text{JD}_{\text{TT}}$)**:
   Starting from the star's Hipparcos-2 (`I/311`, epoch $\text{JD}_1 = 2448349.0625$) ICRS coordinates $(\alpha_1, \delta_1)$, proper motions:
   $$\mu_\alpha = \frac{\mu_{\alpha*} \cdot \text{MAS2RAD}}{\cos\delta_1}, \qquad \mu_\delta = \mu_\delta \cdot \text{MAS2RAD}$$
   parallax $\varpi$ (arcseconds), and radial velocity $v_r$ (km/s from `XHIP` `V/137D`), propagate the star's direction to $\text{JD}_2 = \text{JD}_{\text{TT}}$ including **geocentric annual parallax** (using Earth's barycentric position $\mathbf{p}_{\text{EB}}(t)$ in AU from DE440s via `sofars::astro::pmpx`), with **no light deflection** and **no annual aberration**:
   $$\hat{\mathbf{u}}_{\text{ICRS}}(t) = \text{pmpx}(\alpha_1, \delta_1, \mu_\alpha, \mu_\delta, \varpi, v_r, \Delta t_{\text{J1991.25}\to t}, \mathbf{p}_{\text{EB}}(t))$$
   *(Why no annual aberration? As documented in `jayeshmepani/jpl-ephemeris` `AYANAMSA_MODELS.md` §1.3, the sidereal origin is a direction in space, not an apparent place; applying annual aberration to the anchor star causes the entire sidereal zodiac origin to oscillate seasonally by $\pm 20.49552''$.)*
2. **ICRS to Mean Equator and Equinox of Date (`bp06`)**:
   Compute the IAU 2006 frame-bias + precession matrix $\mathbf{R}_{\text{bp}}(\text{JD}_{\text{TT}})$ via `sofars::pnm::bp06(jd_tt, 0.0).2`:
   $$\hat{\mathbf{u}}_{\text{mean\_eq}}(t) = \mathbf{R}_{\text{bp}}(\text{JD}_{\text{TT}}) \, \hat{\mathbf{u}}_{\text{ICRS}}(t)$$
3. **Mean Equator of Date to Mean Ecliptic of Date (`obl06`)**:
   Compute the IAU 2006 mean obliquity $\varepsilon_0(\text{JD}_{\text{TT}})$ via `sofars::pnm::obl06(jd_tt, 0.0)` and rotate around the $+X$ axis by $+\varepsilon_0$:
   $$\begin{pmatrix} x_{\text{ecl}} \\ y_{\text{ecl}} \\ z_{\text{ecl}} \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \cos\varepsilon_0 & \sin\varepsilon_0 \\ 0 & -\sin\varepsilon_0 & \cos\varepsilon_0 \end{pmatrix} \hat{\mathbf{u}}_{\text{mean\_eq}}(t)$$
   $$\lambda_{S,\text{mean}}(t) = \operatorname{atan2}(y_{\text{ecl}}, x_{\text{ecl}}) \bmod 2\pi \quad \longrightarrow \quad \text{degrees in } [0^\circ, 360^\circ)$$
4. **Mean & True Ayanamsa**:
   $$A_{\text{mean}}(t) = \bigl(\lambda_{S,\text{mean}}(t) - \lambda_{\text{target}}\bigr) \bmod 360^\circ$$
   $$A_{\text{true}}(t) = \bigl(A_{\text{mean}}(t) + \Delta\psi(t)\bigr) \bmod 360^\circ$$
5. **Secular Rate $\dot{A}_{\text{mean}}(t)$**:
   Computed via central difference of $A_{\text{mean}}$ over $\pm 0.5\text{ d}_{\text{TT}}$:
   $$\dot{A}_{\text{mean}}(t) = \operatorname{wrap}_{[-180^\circ, 180^\circ)}\bigl(A_{\text{mean}}(t + 0.5) - A_{\text{mean}}(t - 0.5)\bigr) \quad [\text{deg/day}_{\text{TT}}]$$

#### Anchor Star Passport Table (Catalog Epoch `J1991.25` = `JD 2448349.0625`, `I/311` + `XHIP V/137D`)

| # | `astra` Key | Legacy ID | Enum Variant | Fiducial Star $S$ | HIP ID | HR / HD | Comp | $\alpha_{\text{ICRS}}(J1991.25)$ (deg) | $\delta_{\text{ICRS}}(J1991.25)$ (deg) | $\mu_{\alpha*}$ (mas/yr) | $\mu_\delta$ (mas/yr) | $\varpi$ (mas) | $v_r$ (km/s) | Target Sidereal Longitude $\lambda_{\text{target}}$ | $A_{\text{mean}}(J2000.0)$ (deg) |
|---:|---|---:|---|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | `true_citra` | `27` | `TrueCitra` | Spica ($\alpha$ Vir) | `65474` | `5056` / `116658` | `A` | `201.29835228` | `-11.16124494` | `-42.35` | `-30.67` | `13.06` | `0.00` (`1.0` BSC5) | `180°00'00"` (`180.0°` = `0° Libra`) | `23.84134212` |
| 11 | `true_revati` | `28` | `TrueRevati` | Revati ($\zeta$ Psc A) | `5737` | `361` / `7344` | `A` | `18.43250842` | `+7.57548938` | `+145.00` | `-55.69` | `18.76` | `+15.00` (`9.0` BSC5) | `359°50'00"` (`359.8333333333333°` = `29°50' Pisces`) | `20.04419835` |
| 12 | `true_pushya` | `29` | `TruePushya` | Pushya ($\delta$ Cnc, Asellus Australis) | `42911` | `3461` / `74442` | `A` | `131.17129191` | `+18.15486373` | `-17.67` | `-229.26` | `24.98` | `+17.14` (`17.0` BSC5) | `106°00'00"` (`106.0°` = `16°00' Cancer`) | `22.72197510` |

---

## 2. The 4 Lunar Nodes & 3 Black Moon Lilith Variants

### 2.1 Shared Physical & Mathematical Constants
- **Time Argument**: $T = (\text{JD}_{\text{TT}} - 2451545.0) / 36525.0$ (Julian centuries of TT from `J2000.0`; per IERS Conventions 2003 §5.7.2, TT is used directly in the lunar fundamental arguments).
- **Earth–Moon Gravitational Parameter (`DE440` Header Constants, Park et al. 2021 Table 2)**:
  $$\mu_{\text{EM}} = GM_\oplus + GM_{\text{Moon}} = 398600.435507 + 4902.800118 = 403503.235625\text{ km}^3/\text{s}^2$$
- **Synthetic Point Output Semantics**:
  For analytical mean points (`Mean_Node`, `Mean_South_Node`, `Mean_Lilith`), ecliptic latitude is identically $\beta = 0^\circ$, latitude rate is $\dot{\beta} = 0^\circ/\text{d}$, and nominal distance is $r = 1.0\text{ AU}$ ($\dot{r} = 0$) unless actual lunar apogee/perigee distance in AU is requested.

---

### 2.2 The 4 Lunar Nodes (`Mean_Node`, `Mean_South_Node`, `True_Node`, `True_South_Node`)

#### Model 1: `Mean_Node` & `Mean_South_Node` (Simon et al. 1994 / IERS Conventions 2003 `faom03`)
The Moon's mean longitude of the ascending node $\Omega(T)$ referred to the **mean ecliptic and mean equinox of date** is given by the IERS 2003 / Simon et al. (1994) polynomial (`sofars::fundargs::faom03`):

$$\Omega_{\text{mean}}(T) = 450160.398036'' - 6962890.5431''\,T + 7.4722''\,T^2 + 0.007702''\,T^3 - 0.00005939''\,T^4 \quad [\text{arcsec}]$$

*(Note: $450160.398036'' = 125^\circ 02' 40.398036'' = 125.04455501^\circ$. In Meeus Ch. 47 / Simon 1994 §3.5, the linear coefficient is $-6962890.2665''\,T = -1934.1361851^\circ\,T$, differing by only $0.2766''/\text{cy}$ from `faom03`.)*

The analytic secular rate in **degrees per TT day** (`deg/day_TT`) is:
$$\dot{\Omega}_{\text{mean}}(T) = \frac{-6962890.5431 + 2(7.4722)\,T + 3(0.007702)\,T^2 - 4(0.00005939)\,T^3}{3600.0 \times 36525.0} \approx -0.05295376\text{ deg/day}$$

- **`Mean_Node` (`LunarPoint::MeanNode`)**:
  - In `NutationMode::Mean`: $\lambda = \Omega_{\text{mean}}(T) \bmod 360^\circ$.
  - In `NutationMode::True` (default tropical): $\lambda = \bigl(\Omega_{\text{mean}}(T) + \Delta\psi(T)\bigr) \bmod 360^\circ$.
  - In `ZodiacMode::Sidereal(model)`: $\lambda_{\text{sid}} = \bigl(\Omega_{\text{mean}}(T) - A_{\text{mean}}(T)\bigr) \bmod 360^\circ$, $\dot{\lambda}_{\text{sid}} = \dot{\Omega}_{\text{mean}}(T) - \dot{A}_{\text{mean}}(T)$.
- **`Mean_South_Node` (`LunarPoint::MeanSouthNode`)**:
  - Exact antipodal point:
    $$\lambda_{\text{south}}(T) = \bigl(\lambda_{\text{north}}(T) + 180.0^\circ\bigr) \bmod 360^\circ, \qquad \dot{\lambda}_{\text{south}}(T) = \dot{\lambda}_{\text{north}}(T)$$

#### Model 2: `True_Node` & `True_South_Node` (Osculating Orbital Plane Intersection from DE440s)
Given target epoch $t$ ($\text{JD}_{\text{TDB}} \approx \text{JD}_{\text{TT}} + \Delta_{\text{TDB-TT}}$):
1. Query DE440s for the instantaneous **geometric** geocentric state vector of the Moon (`NAIF 301`) relative to the Earth (`NAIF 399`) in the ICRF (`J2000`) frame:
   $$\mathbf{r}_{\text{ICRF}} = \mathbf{r}_{301,\text{SSB}}(t) - \mathbf{r}_{399,\text{SSB}}(t) \quad [\text{km}], \qquad \mathbf{v}_{\text{ICRF}} = \mathbf{v}_{301,\text{SSB}}(t) - \mathbf{v}_{399,\text{SSB}}(t) \quad [\text{km/s}]$$
2. Rotate $\mathbf{r}_{\text{ICRF}}$ and $\mathbf{v}_{\text{ICRF}}$ into the target ecliptic frame of date:
   - For `NutationMode::Mean`: $\mathbf{R}_{\text{ecl}}(t) = \mathbf{R}_x(\varepsilon_0(t)) \, \mathbf{R}_{\text{bp}}(t)$
   - For `NutationMode::True`: $\mathbf{R}_{\text{ecl}}(t) = \mathbf{R}_x(\varepsilon_0(t) + \Delta\varepsilon(t)) \, \mathbf{R}_{\text{npb}}(t)$
   $$\mathbf{r} = \mathbf{R}_{\text{ecl}}(t)\,\mathbf{r}_{\text{ICRF}}, \qquad \mathbf{v} = \mathbf{R}_{\text{ecl}}(t)\,\mathbf{v}_{\text{ICRF}}$$
3. Compute the specific angular momentum vector (orbital plane normal):
   $$\mathbf{h} = \mathbf{r} \times \mathbf{v} = (h_x, h_y, h_z)^T$$
4. The line of nodes is the intersection of the lunar orbital plane with the $XY$ ecliptic plane:
   $$\mathbf{n} = \hat{\mathbf{z}}_{\text{ecl}} \times \mathbf{h} = (-h_y, h_x, 0)^T$$
5. Therefore, the osculating ascending node longitude is:
   $$\lambda_{\text{TrueNode}}(t) = \operatorname{atan2}(h_x, -h_y) \bmod 2\pi \quad \longrightarrow \quad \text{degrees in } [0^\circ, 360^\circ)$$
   $$\lambda_{\text{TrueSouthNode}}(t) = \bigl(\lambda_{\text{TrueNode}}(t) + 180.0^\circ\bigr) \bmod 360^\circ$$
6. **Speed $\dot{\lambda}_{\text{TrueNode}}(t)$**:
   Evaluated via the 4-point 4th-order central finite-difference stencil (`h = 1/24 day = 1 hour` or `h = 0.002 day`) on the unwrapped/unit-circle vector $\hat{\mathbf{n}}(t) = (\cos\lambda(t), \sin\lambda(t))$:
   $$\dot{\lambda}(t) = \cos\lambda(t)\,\dot{n}_y(t) - \sin\lambda(t)\,\dot{n}_x(t)$$

---

### 2.3 The 3 Black Moon Lilith Variants (`Mean_Lilith`, `True_Lilith`, `Interpolated_Lilith`)

#### Variant 1: `Mean_Lilith` (`LunarPoint::MeanLilith` — Simon et al. 1994 / IERS 2003 Mean Lunar Apogee)
The Moon's mean longitude of perigee $\varpi_{\text{mean}}(T)$ referred to the **mean ecliptic and mean equinox of date** is given by two equivalent forms that agree to $< 0.14''$ over `1900–2100`:

1. **Direct Simon et al. (1994, *A&A* 282, 663, §3.5) / Chapront ELP2000-82B Secular Perigee Polynomial**:
   $$\varpi_{\text{mean}}(T) = 83^\circ 21' 11.67475'' + 14643420.2669''\,T - 38.2776''\,T^2 - 0.045047''\,T^3 + 0.00021301''\,T^4$$
   *(In decimal degrees: $\varpi_{\text{mean}}(T) = 83.3532430^\circ + 4067.6167408^\circ\,T - 0.0106327^\circ\,T^2 - 1.2513 \times 10^{-5\circ}\,T^3 + 5.917 \times 10^{-8\circ}\,T^4$.)*
2. **IERS Conventions (2003) Fundamental Argument Identity** (`faf03 + faom03 - fal03`):
   $$\varpi_{\text{mean, IERS}}(T) = F(T) + \Omega(T) - l(T)$$
   where:
   $$\begin{aligned}
   l(T) &= 485868.249036'' + 1717915923.2178''\,T + 31.8792''\,T^2 + 0.051635''\,T^3 - 0.00024470''\,T^4 \\
   F(T) &= 335779.526232'' + 1739527262.8478''\,T - 12.7512''\,T^2 - 0.001037''\,T^3 + 0.00000417''\,T^4 \\
   \Omega(T) &= 450160.398036'' - 6962890.5431''\,T + 7.4722''\,T^2 + 0.007702''\,T^3 - 0.00005939''\,T^4
   \end{aligned}$$
   Combining the IERS 2003 coefficients term-by-term yields:
   $$\varpi_{\text{mean, IERS}}(T) = 300071.675232'' + 14648449.0869''\,T - 37.1582''\,T^2 - 0.044970''\,T^3 + 0.00018948''\,T^4 \quad [\text{arcsec}]$$
   *(Note: $300071.675232'' = 83^\circ 21' 11.675232'' = 83.35324312^\circ$, and $14648449.0869'' \bmod 1296000'' = 4069.01363525^\circ/\text{cy}$, matching Simon 1994 / Meeus Ch. 47 to $0.0005''$ at $J2000.0$ and $0.14''$ across a century!)*

**Mean Lunar Apogee (`Mean_Lilith`)**:
$$\lambda_{\text{MeanLilith, mean}}(T) = \bigl(\varpi_{\text{mean, IERS}}(T) + 180.0^\circ\bigr) \bmod 360^\circ$$
- At $J2000.0$ ($T = 0$): $\lambda_{\text{MeanLilith, mean}}(0) = 263^\circ 21' 11.675232'' = 263.35324312^\circ$.
- Analytic secular speed (`deg/day_TT`):
  $$\dot{\lambda}_{\text{MeanLilith}}(T) = \frac{14648449.0869 - 2(37.1582)\,T - 3(0.044970)\,T^2 + 4(0.00018948)\,T^3}{3600.0 \times 36525.0} \approx +0.11140408\text{ deg/day}$$
- In `NutationMode::True`: add $\Delta\psi(T)$.
- *(Note on legacy Swiss Ephemeris `SE_MEAN_APOG`: Legacy SE adds a subset of medium-period solar-perturbation harmonic terms from ELP2000 on top of the secular polynomial while still calling it "Mean Apogee", causing a $\pm 4'$ to $\pm 7'$ oscillation relative to pure secular `Mean_Lilith`. In `lxephemeris`, `Mean_Lilith` is the pure IERS/Simon secular mean apogee, while the solar-perturbed natural apogee is provided cleanly by `Interpolated_Lilith`.)*

---

#### Variant 2: `True_Lilith` (`LunarPoint::TrueLilith`, alias `osculating_apogee`)
The instantaneous **osculating lunar apogee** derived from the Laplace–Runge–Lenz eccentricity vector of the Earth–Moon two-body system using DE440s state vectors:
1. Obtain the Moon's geometric geocentric position $\mathbf{r}$ ($\text{km}$) and velocity $\mathbf{v}$ ($\text{km/s}$) rotated into the target ecliptic frame of date (identical to Step 2 of `True_Node`).
2. Compute specific angular momentum $\mathbf{h} = \mathbf{r} \times \mathbf{v}$ and the Laplace–Runge–Lenz eccentricity vector (which points toward osculating **perigee**):
   $$\mathbf{e} = \frac{\mathbf{v} \times \mathbf{h}}{\mu_{\text{EM}}} - \frac{\mathbf{r}}{\|\mathbf{r}\|}, \qquad \mu_{\text{EM}} = 403503.235625\text{ km}^3/\text{s}^2$$
3. Osculating **apogee** (`True_Lilith`) points in the opposite direction:
   $$\mathbf{a} = -\mathbf{e} = (a_x, a_y, a_z)^T$$
4. Extract ecliptic longitude and latitude:
   $$\lambda_{\text{TrueLilith}}(t) = \operatorname{atan2}(a_y, a_x) \bmod 2\pi \quad \longrightarrow \quad [0^\circ, 360^\circ)$$
   $$\beta_{\text{TrueLilith}}(t) = \operatorname{atan2}\!\left(a_z, \sqrt{a_x^2 + a_y^2}\right)$$
   and osculating apogee distance:
   $$r_{\text{apogee}}(t) = a_{\text{sma}}(t)\,\bigl(1 + \|\mathbf{e}\|\bigr), \qquad a_{\text{sma}}(t) = \left(\frac{2}{\|\mathbf{r}\|} - \frac{\|\mathbf{v}\|^2}{\mu_{\text{EM}}}\right)^{-1}$$
5. **Speed**: 4-point 4th-order central finite difference on the unit vector $\hat{\mathbf{a}}_{xy}(t) = (\cos\lambda(t), \sin\lambda(t))$ with step $h = 1\text{ hour}$ (`1/24 day`). Because the Moon's eccentricity $e \approx 0.0549$ is small and strongly perturbed by the Sun's evection, $\lambda_{\text{TrueLilith}}(t)$ oscillates by $\pm 30^\circ$ over the 31.8-day evection cycle and its true instantaneous speed ranges from $-20^\circ/\text{d}$ to $+15^\circ/\text{d}$.

---

#### Variant 3: `Interpolated_Lilith` (`LunarPoint::InterpolatedLilith` — Natural Lunar Apogee via $C^2$ Cubic Spline)
Because osculating `True_Lilith` treats the strongly 3-body Sun–Earth–Moon system as an unperturbed 2-body Keplerian ellipse at every instant, it exhibits large artificial oscillations ($\pm 30^\circ$) even at the moments when the Moon actually passes through apogee or perigee. Conversely, at the actual physical apsides where radial velocity vanishes ($\dot{r}(t_k) = \mathbf{r}(t_k)\cdot\mathbf{v}(t_k)/\|\mathbf{r}(t_k)\| = 0$), the Moon's true position coincides with the physical line of apsides (Dieter Koch 2000; Chapront-Touzé & Chapront 1983).

`lxephemeris` defines **`Interpolated_Lilith` (`natural_apogee_c2_spline_v1`)** as a completely self-contained, deterministic mathematical algorithm on the DE440s geocentric lunar state:

##### Complete Deterministic Algorithm (`natural_apogee_c2_spline_v1`)
Given target epoch $t \in \text{safe\_interval} = [\text{1850-03-01T00:00:00 TT}, \text{2149-11-01T00:00:00 TT})$:
1. **Apsidal Node Search Window ($[t - 60\text{ d}, t + 60\text{ d}]$)**:
   Define the geocentric radial rate function in $\text{km/s}$:
   $$g(\tau) = \frac{\mathbf{r}_{\text{ICRF}}(\tau) \cdot \mathbf{v}_{\text{ICRF}}(\tau)}{\|\mathbf{r}_{\text{ICRF}}(\tau)\|}$$
   Since the lunar anomalistic half-period is $\frac{1}{2} P_{\text{anom}} \approx 13.777\text{ d}$, the $120\text{ d}$ window $[t - 60\text{ d}, t + 60\text{ d}]$ contains $8$ to $9$ alternating apogee and perigee extrema.
2. **Bracket & Refine Roots of $g(\tau) = 0$**:
   - Step across $[t - 60\text{ d}, t + 60\text{ d}]$ with uniform step $\Delta\tau = 2.0\text{ d}$ ($61$ grid evaluations). Since successive apsides are separated by at least $12.5\text{ d}$, every root of $g(\tau) = 0$ is strictly bracketed by a sign change $g(\tau_j)\cdot g(\tau_{j+1}) \le 0$.
   - Refine each bracket $[\tau_j, \tau_{j+1}]$ to machine precision ($|\Delta\tau| < 10^{-11}\text{ d} \approx 1\,\mu\text{s}$) using Brent's method (`lxephemeris_core::math::roots::brent`).
   - Classify each root $\tau_k$:
     - If $g(\tau)$ changes from $+$ to $-$ ($g'(\tau_k) < 0$, distance maximum $\ddot{r}(\tau_k) < 0$): $\tau_k$ is a **True Apogee**.
     - If $g(\tau)$ changes from $-$ to $+$ ($g'(\tau_k) > 0$, distance minimum $\ddot{r}(\tau_k) > 0$): $\tau_k$ is a **True Perigee**.
3. **Node Longitudes & Distances in Target Frame of Epoch $\tau_k$**:
   At each extremum $\tau_k$ ($k = 0, \dots, M-1$, where $M \in \{8, 9\}$), evaluate the Moon's geometric geocentric position $\mathbf{r}_{\text{ecl}}(\tau_k)$ in the mean ecliptic of date $\tau_k$, with spherical coordinates $(\lambda_{\text{Moon}}(\tau_k), \beta_{\text{Moon}}(\tau_k), r_{\text{Moon}}(\tau_k))$:
   - If $\tau_k$ is a **True Apogee**:
     $$\Lambda_k = \lambda_{\text{Moon}}(\tau_k), \qquad B_k = \beta_{\text{Moon}}(\tau_k)$$
   - If $\tau_k$ is a **True Perigee** (projected to the apogee side by $180^\circ$ opposition):
     $$\Lambda_k = \bigl(\lambda_{\text{Moon}}(\tau_k) + 180.0^\circ\bigr) \bmod 360^\circ, \qquad B_k = -\beta_{\text{Moon}}(\tau_k)$$
4. **Phase Unwrapping**:
   Unwrap the sequence $\{\Lambda_0, \Lambda_1, \dots, \Lambda_{M-1}\}$ into a continuous real-valued sequence $\{\tilde{\Lambda}_0, \dots, \tilde{\Lambda}_{M-1}\}$ so that $|\tilde{\Lambda}_{k+1} - \tilde{\Lambda}_k| < 180^\circ$ for all $k$:
   $$\tilde{\Lambda}_0 = \Lambda_0, \qquad \tilde{\Lambda}_{k+1} = \tilde{\Lambda}_k + \operatorname{wrap}_{[-180^\circ, 180^\circ)}\bigl(\Lambda_{k+1} - \tilde{\Lambda}_k\bigr)$$
5. **Natural $C^2$ Cubic Spline Interpolation**:
   Construct the unique natural $C^2$ cubic spline $S_\Lambda(\tau)$ (`S''(\tau_0) = S''(\tau_{M-1}) = 0`, solved via the $O(M)$ Thomas tridiagonal algorithm) through the knots $\{(\tau_k, \tilde{\Lambda}_k)\}_{k=0}^{M-1}$ (and similarly $S_B(\tau)$ through $\{(\tau_k, B_k)\}_{k=0}^{M-1}$).
6. **Evaluate at Target Epoch $t$**:
   $$\lambda_{\text{InterpLilith, mean}}(t) = S_\Lambda(t) \bmod 360^\circ, \qquad \dot{\lambda}_{\text{InterpLilith}}(t) = S'_\Lambda(t) \quad [\text{deg/day}_{\text{TT}}]$$
   $$\beta_{\text{InterpLilith}}(t) = S_B(t), \qquad \dot{\beta}_{\text{InterpLilith}}(t) = S'_B(t)$$
   If `NutationMode::True` is requested, add $\Delta\psi(t)$ to $\lambda_{\text{InterpLilith, mean}}(t)$.
7. **Coverage & Boundary Safety**:
   Because the $\pm 60\text{ d}$ window requires valid DE440s Moon states on $[t - 60\text{ d}, t + 60\text{ d}]$, `Interpolated_Lilith` declares `safe_interval = [1850-03-01T00:00:00 TT, 2149-11-01T00:00:00 TT)` on `de440s.bsp` (`1849-12-25` to `2150-01-22`). Queries outside `safe_interval` return `Err(LxError::OutOfRange { .. })`.

---

## 3. The 8 Hamburg Uranian Planets + `Transpluto` (`Isis`)

### 3.1 Provenance & Classification
- **8 Classical Hamburg Uranians** (`Cupido`, `Hades`, `Zeus`, `Kronos`, `Apollon`, `Admetos`, `Vulkanus`, `Poseidon`):
  Introduced by Alfred Witte (1919–1928) and Friedrich Sieggrün (1930s) of the Hamburg School; canonical Keplerian orbital elements reconstructed and published by **James Neely** in *Matrix Journal* VII (1980) and Michael Erlewine (1980).
- **`Transpluto` (alias `Isis`)**:
  Hypothetical trans-Plutonian planet proposed by **Émile Sevin** (*Comptes Rendus de l'Académie des Sciences*, Paris, 1946) from gravitational residuals of Uranus/Neptune, with orbital elements published by **W. Strubell** in *Die Sterne* 3/1952 (and Ebertin's *Kosmobiшеs Jahrbuch* 1953). *(Distinct from asteroid `(42) Isis`.)*
- **Permissive Reference Implementation**: `TheDaniel166/moira` (`moira/uranian.py` & `URANIAN_BACKEND_STANDARD.md`, MIT License), verified against `uranian_reference_matrix.json`.

### 3.2 Complete Orbital Element Table

Gaussian gravitational constant in degrees/day:
$$k^\circ = 0.01720209895 \times \frac{180}{\pi} = 0.9856076686^\circ/\text{day}$$
Mean daily motion:
$$n = \frac{k^\circ}{a^{3/2}} \quad [\text{deg/day}_{\text{TT}}]$$

| # | `astra` ID | Alias | Source Reference Epoch $t_0$ ($\text{JD}_{\text{TT}}$) | Source Equinox $\text{JD}_{\text{eq}}$ ($\text{JD}_{\text{TT}}$) | Mean Anomaly $M_0$ (deg) | Semi-Major Axis $a$ (AU) | Eccentricity $e$ | Arg. of Perihelion $\omega$ (deg) | Ascending Node $\Omega$ (deg) | Inclination $i$ (deg) | Mean Motion $n = k^\circ / a^{3/2}$ ($\text{deg/day}$) | Orbital Period $P$ (Julian yr) |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | `Cupido` | — | `2415020.0` (`J1900.0`) | `2415020.0` (`J1900.0`) | `163.7409` | `40.99837` | `0.00460` | `171.4333` | `129.8325` | `1.0833` | `0.0037543866` | `262.536` |
| 2 | `Hades` | — | `2415020.0` (`J1900.0`) | `2415020.0` (`J1900.0`) | `27.6496` | `50.66744` | `0.00245` | `148.1796` | `161.3339` | `1.0500` | `0.0027329447` | `360.659` |
| 3 | `Zeus` | — | `2415020.0` (`J1900.0`) | `2415020.0` (`J1900.0`) | `165.1232` | `59.21436` | `0.00120` | `299.0440` | `0.0000` | `0.0000` | `0.0021630106` | `455.689` |
| 4 | `Kronos` | — | `2415020.0` (`J1900.0`) | `2415020.0` (`J1900.0`) | `169.0193` | `64.81690` | `0.00305` | `208.8801` | `0.0000` | `0.0000` | `0.0018886875` | `521.876` |
| 5 | `Apollon` | — | `2415020.0` (`J1900.0`) | `2415020.0` (`J1900.0`) | `138.0533` | `70.29949` | `0.00000` | `0.0000` | `0.0000` | `0.0000` | `0.0016726243` | `589.289` |
| 6 | `Admetos` | — | `2415020.0` (`J1900.0`) | `2415020.0` (`J1900.0`) | `351.3350` | `73.62765` | `0.00000` | `0.0000` | `0.0000` | `0.0000` | `0.0015604618` | `631.646` |
| 7 | `Vulkanus` | — | `2415020.0` (`J1900.0`) | `2415020.0` (`J1900.0`) | `55.8983` | `77.25568` | `0.00000` | `0.0000` | `0.0000` | `0.0000` | `0.0014518442` | `678.895` |
| 8 | `Poseidon` | — | `2415020.0` (`J1900.0`) | `2415020.0` (`J1900.0`) | `165.5163` | `83.66907` | `0.00000` | `0.0000` | `0.0000` | `0.0000` | `0.0012881224` | `765.184` |
| 9 | `Transpluto` | `Isis` | `2368547.66` (`1772-08-29 03:50:24 TT`) | `2431456.5` (`1945-01-01 00:00:00 TT`) | `0.0000` | `77.77500` | `0.30000` | `0.7000` | `0.0000` | `0.0000` | `0.0014373284` | `685.751` |

*(Note on `Kronos` semi-major axis: Some secondary prose summaries contain a digit transposition `64.81960 AU`, whereas the active Neely/seorbel element table uses `64.81690 AU`, verified against `uranian_reference_matrix.json` at both `J2000.0` and `1980-06-15`.)*

### 3.3 Exact Mathematical Reduction Pipeline for Uranians & `Transpluto`

Given evaluation epoch $t$ ($\text{JD}_{\text{TT}}$):

1. **Heliocentric Keplerian State at Light-Time Retarded Epoch $\tau = t - \Delta t_{\text{LT}}$**:
   - Propagate mean anomaly (in radians):
     $$M(\tau) = \bigl(M_0 + n\,(\tau - t_0)\bigr) \cdot \frac{\pi}{180}$$
   - Solve Kepler's equation $E - e\sin E = M(\tau)$ via Newton–Raphson to $|E_{k+1} - E_k| < 10^{-14}$:
     $$E_0 = M, \qquad E_{k+1} = E_k - \frac{E_k - e\sin E_k - M}{1 - e\cos E_k}$$
   - Compute orbital-plane position $(x_{\text{orb}}, y_{\text{orb}})$ in AU:
     $$x_{\text{orb}} = a(\cos E - e), \qquad y_{\text{orb}} = a\sqrt{1 - e^2}\,\sin E$$
   - Rotate from the orbital plane to the **source ecliptic frame** of $\text{JD}_{\text{eq}}$ via $(\omega, \Omega, i)$:
     $$\begin{aligned}
     x_{\text{ecl,eq}} &= (\cos\Omega\cos\omega - \sin\Omega\sin\omega\cos i)\,x_{\text{orb}} + (-\cos\Omega\sin\omega - \sin\Omega\cos\omega\cos i)\,y_{\text{orb}} \\
     y_{\text{ecl,eq}} &= (\sin\Omega\cos\omega + \cos\Omega\sin\omega\cos i)\,x_{\text{orb}} + (-\sin\Omega\sin\omega + \cos\Omega\cos\omega\cos i)\,y_{\text{orb}} \\
     z_{\text{ecl,eq}} &= (\sin\omega\sin i)\,x_{\text{orb}} + (\cos\omega\sin i)\,y_{\text{orb}}
     \end{aligned}$$
2. **Transform Source Ecliptic ($\text{JD}_{\text{eq}}$) $\to$ ICRF (`J2000`)**:
   - Rotate from source ecliptic of $\text{JD}_{\text{eq}}$ to mean equator of $\text{JD}_{\text{eq}}$ using mean obliquity $\varepsilon_0(\text{JD}_{\text{eq}})$ (`obl06`):
     $$\mathbf{r}_{\text{eq}}(\text{JD}_{\text{eq}}) = \mathbf{R}_x\bigl(-\varepsilon_0(\text{JD}_{\text{eq}})\bigr)\,\mathbf{r}_{\text{ecl,eq}}$$
   - Rotate from mean equator of $\text{JD}_{\text{eq}}$ to ICRF via the transpose of the IAU 2006 bias-precession matrix $\mathbf{R}_{\text{bp}}(\text{JD}_{\text{eq}})$ (`bp06`):
     $$\mathbf{r}_{\text{helio, ICRF}}(\tau) = \mathbf{R}_{\text{bp}}(\text{JD}_{\text{eq}})^T\,\mathbf{r}_{\text{eq}}(\text{JD}_{\text{eq}})$$
3. **Geocentric Light-Time Iteration & Apparent Reduction**:
   - Let $\mathbf{r}_{\text{Earth, helio, ICRF}}(t) = \mathbf{r}_{399,\text{SSB}}(t) - \mathbf{r}_{10,\text{SSB}}(t)$ be Earth's heliocentric ICRF position at epoch $t$ from DE440s (in AU), and $\mathbf{v}_{\text{Earth, SSB}}(t)$ be Earth's barycentric velocity (in AU/day).
   - Iterate light time $\Delta t_{\text{LT}}$ (2 iterations suffice to $< 10^{-12}\text{ d}$):
     $$\Delta t_{\text{LT}}^{(0)} = 0, \qquad \boldsymbol{\rho}^{(k)} = \mathbf{r}_{\text{helio, ICRF}}\!\bigl(t - \Delta t_{\text{LT}}^{(k)}\bigr) - \mathbf{r}_{\text{Earth, helio, ICRF}}(t), \qquad \Delta t_{\text{LT}}^{(k+1)} = \frac{\|\boldsymbol{\rho}^{(k)}\|}{c_{\text{AU/day}}}$$
     where $c_{\text{AU/day}} = 173.14463267424034\text{ AU/day}$.
   - Apply stellar annual aberration (`sofars::astro::ab` or first-order $\hat{\boldsymbol{\rho}}' = \text{normalize}(\hat{\boldsymbol{\rho}} + \mathbf{v}_{\text{Earth, SSB}}(t)/c)$).
   - Rotate from ICRF to true (or mean) equator and ecliptic of date $t$ via $\mathbf{R}_{\text{npb}}(t)$ and $\varepsilon(t) = \varepsilon_0(t) + \Delta\varepsilon(t)$.
4. **Verification Benchmarks (`J2000.0` = `JD 2451545.0 TT`, Apparent Geocentric Tropical Ecliptic Longitude)**:
   - `Cupido`: $258.5906^\circ$ ($18^\circ 35'\text{ Sagittarius}$)
   - `Hades`: $82.7712^\circ$ ($22^\circ 46'\text{ Gemini}$)
   - `Zeus`: $186.8091^\circ$ ($6^\circ 49'\text{ Libra}$)
   - `Kronos`: $99.1369^\circ$ ($9^\circ 08'\text{ Cancer}$)
   - `Apollon`: $200.6290^\circ$ ($20^\circ 38'\text{ Libra}$)
   - `Admetos`: $50.0085^\circ$ ($20^\circ 01'\text{ Taurus}$)
   - `Vulkanus`: $110.2286^\circ$ ($20^\circ 14'\text{ Cancer}$)
   - `Poseidon`: $213.9721^\circ$ ($3^\circ 58'\text{ Scorpio}$)
   - `Transpluto`: $144.8965^\circ$ ($24^\circ 54'\text{ Leo}$)

---

## 4. The 50 Fixed Stars of `astra` (`STAR_CATALOG_50`)

### 4.1 Catalog Provenance, Licensing, and Astrometric Reduction

1. **Primary Astrometric Source**:
   - **Hipparcos, the New Reduction (`VizieR I/311`, van Leeuwen 2007, *A&A* 474, 653)**: Provides ICRS right ascension $\alpha_{\text{ICRS}}$ and declination $\delta_{\text{ICRS}}$ at catalog epoch **`J1991.25` ($\text{JD}_{\text{TT}} = 2448349.0625$)**, trigonometric parallax $\varpi$ (mas), and proper motions $\mu_{\alpha*} = \mu_\alpha\cos\delta$ (mas/yr) and $\mu_\delta$ (mas/yr).
   - **Extended Hipparcos Compilation (`VizieR V/137D`, XHIP, Anderson & Francis 2012, *Astron. Lett.* 38, 331)** & **Yale Bright Star Catalogue 5th Ed. (`VizieR V/50`, Hoffleit & Warren 1991)**: Provides heliocentric radial velocity $v_r$ (km/s), component multiplicity flags (`Comp`), `HR` / `HD` cross-identifications, and visual Johnson $V$ magnitude (`Vmag`).
   - **Licensing & Clean-Room Status**: All 50 entries are physical measurements of natural stars extracted from public scientific astronomical catalogs (`I/311`, `V/137D`, `V/50`) and cross-verified against the MIT-licensed `jayeshmepani/jpl-ephemeris` (`fixstar_hipparcos.inc`). Zero data or strings were read from `sefstars.txt`.
2. **Astrometric Reduction Pipeline (`J1991.25` $\to$ Target Epoch $t$)**:
   - **Never treat `J1991.25` coordinates as `J2000.0` without proper-motion propagation!** Over the $8.75\text{ yr}$ interval from `J1991.25` to `J2000.0`, high-proper-motion stars move by up to $20''$ (e.g., Arcturus moves by $19.9''$, Sirius by $11.7''$, Procyon by $11.0''$).
   - At target epoch $t$ ($\text{JD}_{\text{TDB}}$):
     1. Convert `I/311` `J1991.25` angles and rates to radians and radians/Julian year:
        $$\alpha_1 = \alpha_{\text{J1991.25}}\,\frac{\pi}{180}, \quad \delta_1 = \delta_{\text{J1991.25}}\,\frac{\pi}{180}, \quad \mu_\alpha = \frac{\mu_{\alpha*} \cdot \text{MAS2RAD}}{\cos\delta_1}, \quad \mu_\delta = \mu_\delta \cdot \text{MAS2RAD}$$
     2. Propagate from $\text{JD}_1 = 2448349.0625$ (`J1991.25`) to $\text{JD}_2 = \text{JD}_{\text{TDB}}$ via `sofars::astro::pmpx` (for geocentric observer $\mathbf{p}_{\text{EB}}(t)$, applying 3D space motion + annual parallax).
     3. For `AbMode::Apparent` (default): apply relativistic solar gravitational light deflection (`sofars::astro::ldsun`) and stellar annual aberration (`sofars::astro::ab`).
     4. Rotate from ICRS to the requested output frame (`Frame::Ecliptic` or `Frame::Equatorial`, `NutationMode::True` or `NutationMode::Mean`, `ZodiacMode::Tropical` or `ZodiacMode::Sidereal`).
3. **Quality Flags & Edge-Case Policies**:
   - **Parallax Policy**: All 50 `astra` stars have strictly positive trigonometric parallaxes in `I/311` ($\varpi_{\min} = 1.65\text{ mas}$ for Alnilam, $\varpi_{\max} = 379.21\text{ mas}$ for Sirius). For any general catalog entry with $\varpi \le 0$, clamp $\varpi_{\text{eff}} = 0.0$ in `pmpx` (angular-only space motion) and emit quality flag `PLX_LE_ZERO`.
   - **Radial Velocity Policy**: All 50 `astra` stars have measured radial velocities in both `XHIP V/137D` and `BSC5 V/50`. Both values are recorded below; `XHIP V/137D` is used as primary.
   - **Multi-Component Systems (`Comp` Flag)**: Multi-component systems are explicitly disambiguated by their Hipparcos component identifier (`COMP_A`, `MULTICOMP(AB)`, `MULTICOMP(AP)`, `MULTICOMP(ABC)`):
     - `Castor` ($\alpha$ Gem): `HIP 36850` (`HR 2891` / `HD 60179`, component `AB` joint Hipparcos solution).
     - `Zubenelgenubi` ($\alpha^2$ Lib): `HIP 72622` (`HR 5531` / `HD 130841`, primary component `A`).
     - `Acrab` ($\beta^1$ Sco): `HIP 78820` (`HR 5984` / `HD 144217`, primary component `AB`).
     - `Rasalgethi` ($\alpha^1$ Her): `HIP 84345` (`HR 6406` / `HD 156014`, primary component `AB`).
     - `Alcyone` ($\eta$ Tau): `HIP 17702` (`HR 1165` / `HD 23630`, primary star $\eta$ Tauri `A`, not a cluster centroid).

### 4.2 Complete 50-Star Astrometric Cross-Match Table (+ 2 Additional Ayanamsa Anchor Stars)

| # | `astra` Name | Bayer Key | HIP ID | HR | HD | Comp | $\alpha_{\text{ICRS}}(J1991.25)$ (deg) | $\delta_{\text{ICRS}}(J1991.25)$ (deg) | $\alpha_{\text{ICRS}}(J2000.0)$ (deg) | $\delta_{\text{ICRS}}(J2000.0)$ (deg) | $\mu_{\alpha*}$ (mas/yr) | $\mu_\delta$ (mas/yr) | $\varpi$ (mas) | $v_r$ XHIP (km/s) | $v_r$ BSC5 (km/s) | $V$ (mag) | Quality / Component Flags |
|---:|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | `Aldebaran` | `alTau` | `21421` | `1457` | `29139` | `A` | `68.98000194` | `+16.50976158` | `68.98016279` | `+16.50930236` | `+63.45` | `-188.94` | `48.94` | `+54.26` | `+54.0` | `0.85` | `COMP_A` |
| 2 | `Regulus` | `alLeo` | `49669` | `3982` | `87901` | `A` | `152.09358042` | `+11.96719519` | `152.09296244` | `+11.96720878` | `-248.73` | `+5.59` | `41.13` | `+5.90` | `+6.0` | `1.35` | `COMP_A` |
| 3 | `Antares` | `alSco` | `80763` | `6134` | `148478` | `AB` | `247.35194829` | `-26.43194598` | `247.35191542` | `-26.43200261` | `-12.11` | `-23.30` | `5.89` | `-3.50` | `-3.0` | `0.96` | `MULTICOMP(AB)` |
| 4 | `Fomalhaut` | `alPsA` | `113368` | `8728` | `216956` | `A` | `344.41177299` | `-29.62183680` | `344.41269272` | `-29.62223703` | `+328.95` | `-164.67` | `129.81` | `+6.50` | `+7.0` | `1.16` | `OK` |
| 5 | `Sirius` | `alCMa` | `32349` | `2491` | `48915` | `ABC` | `101.28854105` | `-16.71314306` | `101.28715536` | `-16.71611585` | `-546.01` | `-1223.07` | `379.21` | `-5.50` | `-8.0` | `-1.46` | `MULTICOMP(ABC)` |
| 6 | `Spica` | `alVir` | `65474` | `5056` | `116658` | `A` | `201.29835228` | `-11.16124494` | `201.29824736` | `-11.16131949` | `-42.35` | `-30.67` | `13.06` | `0.00` | `+1.0` | `0.98` | `COMP_A` (Anchor: `true_citra`) |
| 7 | `Algol` | `bePer` | `14576` | `936` | `19356` | `A` | `47.04220893` | `+40.95565071` | `47.04221855` | `+40.95564668` | `+2.99` | `-1.66` | `36.27` | `+3.70` | `+4.0` | `2.12` | `COMP_A` |
| 8 | `Arcturus` | `alBoo` | `69673` | `5340` | `124897` | `AB` | `213.91811408` | `+19.18727046` | `213.91530029` | `+19.18240916` | `-1093.39` | `-2000.06` | `88.83` | `-5.24` | `-5.0` | `-0.04` | `MULTICOMP(AB)` |
| 9 | `Vega` | `alLyr` | `91262` | `7001` | `172167` | `A` | `279.23410825` | `+38.78299326` | `279.23473480` | `+38.78368897` | `+200.94` | `+286.23` | `130.23` | `-20.60` | `-14.0` | `0.03` | `COMP_A` |
| 10 | `Capella` | `alAur` | `24608` | `1708` | `34029` | `AP` | `79.17206466` | `+45.99902905` | `79.17232794` | `+45.99799149` | `+75.25` | `-426.89` | `76.20` | `+29.19` | `+30.0` | `0.08` | `MULTICOMP(AP)` |
| 11 | `Procyon` | `alCMi` | `37279` | `2943` | `61421` | `AB` | `114.82724202` | `+5.22750758` | `114.82549791` | `+5.22498756` | `-714.59` | `-1036.80` | `284.56` | `-2.70` | `-3.0` | `0.38` | `MULTICOMP(AB)` |
| 12 | `Rigel` | `beOri` | `24436` | `1713` | `34085` | `ABC` | `78.63446385` | `-8.20163958` | `78.63446707` | `-8.20163836` | `+1.31` | `+0.50` | `3.78` | `+17.80` | `+21.0` | `0.12` | `MULTICOMP(ABC)` |
| 13 | `Betelgeuse` | `alOri` | `27989` | `2061` | `39801` | `AP` | `88.79287149` | `+7.40703653` | `88.79293899` | `+7.40706400` | `+27.54` | `+11.30` | `6.55` | `+21.91` | `+21.0` | `0.50` | `MULTICOMP(AP)` |
| 14 | `Altair` | `alAql` | `97649` | `7557` | `187642` | `A` | `297.69450819` | `+8.86738473` | `297.69582736` | `+8.86832124` | `+536.23` | `+385.29` | `194.95` | `-26.60` | `-26.0` | `0.77` | `COMP_A` |
| 15 | `Pollux` | `beGem` | `37826` | `2990` | `62509` | `A` | `116.33068294` | `+28.02631022` | `116.32895777` | `+28.02619889` | `-626.55` | `-45.80` | `96.54` | `+3.23` | `+3.0` | `1.14` | `COMP_A` |
| 16 | `Castor` | `alGem` | `36850` | `2891` | `60179` | `AB` | `113.65001968` | `+31.88863511` | `113.64947164` | `+31.88828222` | `-191.45` | `-145.19` | `64.12` | `+5.40` | `+6.0` | `1.98` | `MULTICOMP(AB)` |
| 17 | `Bellatrix` | `gaOri` | `25336` | `1790` | `35468` | `A` | `81.28278339` | `+6.34973457` | `81.28276356` | `+6.34970326` | `-8.11` | `-12.88` | `12.92` | `+18.20` | `+18.0` | `1.64` | `COMP_A` |
| 18 | `Alnilam` | `epOri` | `26311` | `1903` | `37128` | `A` | `84.05338544` | `-1.20191724` | `84.05338894` | `-1.20191914` | `+1.44` | `-0.78` | `1.65` | `+27.30` | `+26.0` | `1.70` | `COMP_A` |
| 19 | `Deneb` | `alCyg` | `102098` | `7924` | `197345` | `A` | `310.35797281` | `+45.28033431` | `310.35797975` | `+45.28033881` | `+2.01` | `+1.85` | `2.31` | `-4.90` | `-5.0` | `1.25` | `COMP_A` |
| 20 | `Achernar` | `alEri` | `7588` | `472` | `10144` | `A` | `24.42813208` | `-57.23665985` | `24.42852282` | `-57.23675279` | `+87.00` | `-38.24` | `23.39` | `+18.60` | `+16.0` | `0.46` | `OK` |
| 21 | `Hamal` | `alAri` | `9884` | `617` | `12929` | `A` | `31.79285751` | `+23.46277747` | `31.79335710` | `+23.46241755` | `+188.55` | `-148.08` | `49.56` | `-14.64` | `-14.0` | `2.00` | `OK` |
| 22 | `Mirach` | `beAnd` | `5447` | `337` | `6860` | `A` | `17.43249023` | `+35.62083036` | `17.43301617` | `+35.62055765` | `+175.90` | `-112.20` | `16.52` | `+0.06` | `+3.0` | `2.06` | `COMP_A` |
| 23 | `Alioth` | `epUMa` | `62956` | `4905` | `112185` | `AB` | `193.50680405` | `+55.95984299` | `193.50728997` | `+55.95982296` | `+111.91` | `-8.24` | `39.51` | `-12.70` | `-9.0` | `1.77` | `MULTICOMP(AB)` |
| 24 | `Alcyone` | `etTau` | `17702` | `1165` | `23630` | `A` | `56.87110081` | `+24.10524179` | `56.87115231` | `+24.10513565` | `+19.34` | `-43.67` | `8.09` | `+5.40` | `+10.0` | `2.87` | `COMP_A` |
| 25 | `Alphecca` | `alCrB` | `76267` | `5793` | `139006` | `A` | `233.67162276` | `+26.71491051` | `233.67195002` | `+26.71469278` | `+120.27` | `-89.58` | `43.46` | `+1.40` | `+2.0` | `2.23` | `COMP_A` |
| 26 | `Markab` | `alPeg` | `113963` | `8781` | `218045` | `A` | `346.19007056` | `+15.20536753` | `346.19022269` | `+15.20526715` | `+60.40` | `-41.30` | `24.46` | `-2.70` | `-4.0` | `2.49` | `OK` |
| 27 | `Zubenelgenubi` | `al-2Lib` | `72622` | `5531` | `130841` | `A` | `222.71990516` | `-16.04161027` | `222.71963789` | `-16.04177652` | `-105.68` | `-68.40` | `43.03` | `-22.00` | `-10.0` | `2.75` | `COMP_A` |
| 28 | `Zubeneschamali` | `beLib` | `74785` | `5685` | `135742` | `*` | `229.25196592` | `-9.38286665` | `229.25172425` | `-9.38291441` | `-98.10` | `-19.65` | `17.62` | `-35.60` | `-35.0` | `2.61` | `MULTICOMP(*)` |
| 29 | `Acrab` | `be-1Sco` | `78820` | `5984` | `144217` | `AB` | `241.35931336` | `-19.80539435` | `241.35929993` | `-19.80545278` | `-5.20` | `-24.04` | `8.07` | `-1.00` | `-1.0` | `2.62` | `MULTICOMP(AB)` |
| 30 | `Dschubba` | `deSco` | `78401` | `5953` | `143275` | `AB` | `240.08338223` | `-22.62162036` | `240.08335535` | `-22.62170643` | `-10.21` | `-35.41` | `6.64` | `-6.00` | `-7.0` | `2.32` | `MULTICOMP(AB)` |
| 31 | `Ras Alhague` | `alOph` | `86032` | `6556` | `159561` | `A` | `263.73335361` | `+12.56057593` | `263.73362272` | `+12.56003740` | `+108.07` | `-221.57` | `67.13` | `+11.70` | `+13.0` | `2.08` | `COMP_A` |
| 32 | `Rasalgethi` | `alHer` | `84345` | `6406` | `156014` | `AB` | `258.66192746` | `+14.39025294` | `258.66190909` | `+14.39034061` | `-7.32` | `+36.07` | `9.07` | `-32.09` | `-33.0` | `3.48` | `MULTICOMP(AB)` |
| 33 | `Canopus` | `alCar` | `30438` | `2326` | `45348` | `A` | `95.98787790` | `-52.69571787` | `95.98795783` | `-52.69566138` | `+19.93` | `+23.24` | `10.55` | `+20.30` | `+21.0` | `-0.72` | `OK` |
| 34 | `Alphard` | `alHya` | `46390` | `3748` | `81797` | `A` | `141.89688204` | `-8.65868307` | `141.89684460` | `-8.65859953` | `-15.23` | `+34.37` | `18.09` | `-4.70` | `-4.0` | `1.98` | `COMP_A` |
| 35 | `Dubhe` | `alUMa` | `54061` | `4301` | `95689` | `AB` | `165.93265337` | `+61.75111903` | `165.93196467` | `+61.75103469` | `-134.11` | `-34.70` | `26.54` | `-9.40` | `-9.0` | `1.79` | `MULTICOMP(AB)` |
| 36 | `Merak` | `beUMa` | `53910` | `4295` | `95418` | `A` | `165.45996142` | `+56.38234469` | `165.46031891` | `+56.38242609` | `+81.43` | `+33.49` | `40.90` | `-13.10` | `-12.0` | `2.37` | `OK` |
| 37 | `Phecda` | `gaUMa` | `58001` | `4554` | `103287` | `A` | `178.45725512` | `+53.69473297` | `178.45769716` | `+53.69475973` | `+107.68` | `+11.01` | `39.21` | `-11.90` | `-13.0` | `2.44` | `OK` |
| 38 | `Polaris` | `alUMi` | `11767` | `424` | `8890` | `A` | `37.94614300` | `+89.26413778` | `37.95456067` | `+89.26410897` | `+44.48` | `-11.85` | `7.54` | `-16.42` | `-17.0` | `2.02` | `COMP_A` |
| 39 | `Scheat` | `bePeg` | `113881` | `8775` | `217906` | `A` | `345.94305579` | `+28.08245431` | `345.94357275` | `+28.08278712` | `+187.65` | `+136.93` | `16.64` | `+7.99` | `+9.0` | `2.42` | `COMP_A` |
| 40 | `Algenib` | `gaPeg` | `1067` | `39` | `886` | `A` | `3.30895847` | `+15.18361610` | `3.30896346` | `+15.18359354` | `+1.98` | `-9.28` | `8.33` | `+3.20` | `+4.0` | `2.83` | `COMP_A` |
| 41 | `Alpheratz` | `alAnd` | `677` | `15` | `358` | `A` | `2.09653385` | `+29.09082837` | `2.09691619` | `+29.09043112` | `+137.46` | `-163.44` | `33.62` | `-10.10` | `-12.0` | `2.06` | `COMP_A` |
| 42 | `Menkar` | `alCet` | `14135` | `911` | `18884` | `A` | `45.56991317` | `+4.08992556` | `45.56988780` | `+4.08973877` | `-10.41` | `-76.85` | `13.09` | `-26.08` | `-26.0` | `2.53` | `OK` |
| 43 | `Mira` | `omiCet` | `10826` | `681` | `14386` | `AP` | `34.83661105` | `-2.97706076` | `34.83663376` | `-2.97763767` | `+9.33` | `-237.36` | `10.91` | `+63.50` | `+64.0` | `3.04` | `MULTICOMP(AP)` |
| 44 | `Baten Kaitos` | `zeCet` | `8645` | `539` | `11353` | `A` | `27.86504430` | `-10.33494527` | `27.86514510` | `-10.33503581` | `+40.80` | `-37.25` | `13.88` | `+7.80` | `+9.0` | `3.73` | `COMP_A` |
| 45 | `Difda` | `beCet` | `3419` | `188` | `4128` | `A` | `10.89678447` | `-17.98668407` | `10.89737874` | `-17.98660632` | `+232.55` | `+31.99` | `33.86` | `+13.10` | `+13.0` | `2.04` | `OK` |
| 46 | `Enif` | `epPeg` | `107315` | `8308` | `206778` | `A` | `326.04641750` | `+9.87500758` | `326.04648391` | `+9.87500865` | `+26.92` | `+0.44` | `4.73` | `+3.39` | `+5.0` | `2.39` | `COMP_A` |
| 47 | `Sadalmelik` | `alAqr` | `109074` | `8414` | `209750` | `A` | `331.44593869` | `-0.31982628` | `331.44598305` | `-0.31984910` | `+18.25` | `-9.39` | `6.23` | `+6.85` | `+8.0` | `2.96` | `COMP_A` |
| 48 | `Sadalsuud` | `beAqr` | `106278` | `8232` | `204867` | `A` | `322.88966964` | `-5.57115560` | `322.88971548` | `-5.57117555` | `+18.77` | `-8.21` | `6.07` | `+6.28` | `+7.0` | `2.91` | `COMP_A` |
| 49 | `Deneb Algedi` | `deCap` | `107556` | `8322` | `207098` | `A` | `326.75952220` | `-16.12656594` | `326.76018434` | `-16.12728709` | `+261.70` | `-296.70` | `84.27` | `-3.40` | `-6.0` | `2.87` | `COMP_A` |
| 50 | `Nashira` | `gaCap` | `106985` | `8278` | `206088` | `A` | `325.02225947` | `-16.66225299` | `325.02273533` | `-16.66230756` | `+187.56` | `-22.45` | `20.77` | `-31.20` | `-31.0` | `3.68` | `COMP_A` |
| A1 | `Revati (zeta Psc A)` | `zePsc` | `5737` | `361` | `7344` | `A` | `18.43250842` | `+7.57548938` | `18.43286395` | `+7.57535402` | `+145.00` | `-55.69` | `18.76` | `+15.00` | `+9.0` | `5.24` | `COMP_A` (Anchor: `true_revati`) |
| A2 | `Pushya (delta Cnc)` | `deCnc` | `42911` | `3461` | `74442` | `A` | `131.17129191` | `+18.15486373` | `131.17124671` | `+18.15430650` | `-17.67` | `-229.26` | `24.98` | `+17.14` | `+17.0` | `3.94` | `COMP_A` (Anchor: `true_pushya`) |

---

## 5. The 23 House Systems, 8 Chart Angles, and 36 Gauquelin Sectors (`HouseSystem`)

### 5.1 Fundamental Spherical-Trigonometric Conventions & Inputs

Given an evaluation epoch $t$ ($\text{JD}_{\text{UT1}}$ and $\text{JD}_{\text{TT}}$) and geodetic observer coordinates $(\lambda_{\text{obs}}, \phi)$ with east-positive longitude $\lambda_{\text{obs}} \in (-180^\circ, +180^\circ]$ and north-positive geodetic latitude $\phi \in [-90^\circ, +90^\circ]$:

1. **Right Ascension of the Midheaven (`RAMC` / `ARMC`)**:
   - In **Tropical** mode (`NutationMode::True`, standard for apparent tropical houses):
     $$\text{RAMC} = \bigl(\text{GAST}(\text{JD}_{\text{UT1}}, \text{JD}_{\text{TT}}) + \lambda_{\text{obs}}\bigr) \bmod 360^\circ$$
     where $\text{GAST}$ is Greenwich Apparent Sidereal Time (`sofars::st::gst06a`) and the obliquity of the ecliptic is the true obliquity $\varepsilon = \varepsilon_0(\text{JD}_{\text{TT}}) + \Delta\varepsilon(\text{JD}_{\text{TT}})$.
   - In **Mean** mode (`NutationMode::Mean`):
     $$\text{RAMC} = \bigl(\text{GMST}(\text{JD}_{\text{UT1}}, \text{JD}_{\text{TT}}) + \lambda_{\text{obs}}\bigr) \bmod 360^\circ, \qquad \varepsilon = \varepsilon_0(\text{JD}_{\text{TT}})$$
   - In **Sidereal** mode (`ZodiacMode::Sidereal(model)`):
     Compute tropical cusps and angles first using $(\text{RAMC}, \phi, \varepsilon_{\text{true}})$, then subtract $A_{\text{true}}(t)$ from all ecliptic longitudes (`cusps[1..=12]` and all ecliptic chart angles, except `ARMC` which is an equatorial right ascension).

2. **Fundamental Projection Primitives**:
   - **Meridian-to-Ecliptic Projection** $\text{MC}(\alpha, \varepsilon)$: Given right ascension $\alpha$, the ecliptic longitude of the intersection of the hour circle $\text{RA} = \alpha$ with the ecliptic is:
     $$\text{MC}(\alpha, \varepsilon) = \operatorname{atan2}\bigl(\sin\alpha, \cos\alpha\cos\varepsilon\bigr) \bmod 360^\circ$$
   - **Horizon-to-Ecliptic Projection** $\text{ASC}(r, \psi, \varepsilon)$: Given meridian right ascension $r$ and pole height (latitude) $\psi$, the rising intersection of the horizon of pole height $\psi$ with the ecliptic is:
     $$\text{ASC}(r, \psi, \varepsilon) = \operatorname{atan2}\bigl(\cos r, -\sin r\cos\varepsilon - \tan\psi\sin\varepsilon\bigr) \bmod 360^\circ$$
   - **Equatorial Frame Unit Vectors (of Date)**:
     Let $r = \text{RAMC}$ and $\phi$ be the observer's geodetic latitude. In the true equatorial triad of date $(\hat{\mathbf{x}} = \Upsilon, \hat{\mathbf{y}}, \hat{\mathbf{z}} = \text{NCP})$:
     $$\begin{aligned}
     \mathbf{e} &= (-\sin r,\, \cos r,\, 0) && \text{(East horizon unit vector)} \\
     \mathbf{z} &= (\cos\phi\cos r,\, \cos\phi\sin r,\, \sin\phi) && \text{(Zenith unit vector)} \\
     \mathbf{n}_{\text{north}} &= \mathbf{z} \times \mathbf{e} = (-\sin\phi\cos r,\, -\sin\phi\sin r,\, \cos\phi) && \text{(North horizon unit vector)} \\
     \mathbf{v}(\lambda) &= (\cos\lambda,\, \sin\lambda\cos\varepsilon,\, \sin\lambda\sin\varepsilon) && \text{(Ecliptic unit vector at longitude } \lambda\text{)}
     \end{aligned}$$

---

### 5.2 The 8 Standard Chart Angles (`ChartAngles`)

All 8 standard chart angles are well-defined across both Northern ($\phi > 0$) and Southern ($\phi < 0$) hemispheres (verified to `0.000000"` against L3 reference across Zurich, London, Sydney, São Paulo, and $\pm 180^\circ$ antimeridian):

| Index | Angle Name | Formula (in terms of $r = \text{RAMC}$, $\phi$, $\varepsilon$) |
| :--- | :--- | :--- |
| `0` | `ASC` (Ascendant) | $\operatorname{atan2}\bigl(\cos r, -\sin r\cos\varepsilon - \tan\phi\sin\varepsilon\bigr) \bmod 360^\circ$ |
| `1` | `MC` (Midheaven / Medium Coeli) | $\operatorname{atan2}\bigl(\sin r, \cos r\cos\varepsilon\bigr) \bmod 360^\circ$ |
| `2` | `ARMC` (Sidereal Time in Degrees) | $r = \text{RAMC} \bmod 360^\circ$ |
| `3` | `Vertex` | $\operatorname{atan2}\bigl(-\cos r, \sin r\cos\varepsilon - \cot\phi\sin\varepsilon\bigr) \bmod 360^\circ$ |
| `4` | `EqASC` (Equatorial Ascendant / East Point) | $\operatorname{atan2}\bigl(\cos r, -\sin r\cos\varepsilon\bigr) \bmod 360^\circ = \text{ASC}(r, 0^\circ, \varepsilon)$ |
| `5` | `CoASC_Koch` (Koch Co-Ascendant) | $\bigl(\text{PolarASC} + 180^\circ\bigr) \bmod 360^\circ$ |
| `6` | `CoASC_Munkasey` (Munkasey Co-Ascendant) | $\operatorname{atan2}\bigl(\cos r, -\sin r\cos\varepsilon - \cot\phi\sin\varepsilon\bigr) \bmod 360^\circ$ |
| `7` | `PolarASC` (Munkasey Polar Ascendant) | $\text{ASC}(r + 180^\circ, \phi, \varepsilon) = \operatorname{atan2}\bigl(-\cos r, \sin r\cos\varepsilon - \tan\phi\sin\varepsilon\bigr) \bmod 360^\circ$ |

- **Opposing Angles**: $\text{DSC} = (\text{ASC} + 180^\circ) \bmod 360^\circ$, $\text{IC} = (\text{MC} + 180^\circ) \bmod 360^\circ$, $\text{AntiVertex} = (\text{Vertex} + 180^\circ) \bmod 360^\circ$.
- **Critical Southern Hemisphere Invariant for `CoASC_Munkasey`**:
  Evaluating $\text{ASC}(r + 180^\circ, 90^\circ - \phi, \varepsilon)$ directly with $90^\circ + \phi$ for $\phi < 0$ negates both arguments of `atan2` and introduces a spurious $180^\circ$ flip in the Southern Hemisphere. Always evaluate $\operatorname{atan2}(\cos r, -\sin r\cos\varepsilon - \cot\phi\sin\varepsilon) \bmod 360^\circ$ directly; it is continuous across both hemispheres (at $\phi = 0^\circ$, $\cot\phi \to \infty$ so `CoASC_Munkasey` and `Vertex` approach $180^\circ$).

---

### 5.3 Complete Catalog of the 23 House Systems

For quadrant systems, Cusp $1 = \text{ASC}$, Cusp $10 = \text{MC}$, Cusp $7 = \text{DSC} = (\text{ASC} + 180^\circ) \bmod 360^\circ$, Cusp $4 = \text{IC} = (\text{MC} + 180^\circ) \bmod 360^\circ$, and opposing intermediate cusps satisfy $C_{k+6} = (C_k + 180^\circ) \bmod 360^\circ$ for $k \in \{11, 12, 2, 3\} \mapsto \{5, 6, 8, 9\}$.

#### 1. Placidus (`'P'`, `HouseSystem::Placidus`)
- **Classification**: Diurnal/nocturnal semi-arc trisection (iterative).
- **Equations**:
  For each diurnal semi-arc fraction $f_k \in \{+1/3, +2/3, -2/3, -1/3\}$ corresponding to Cusps $(11, 12, 8, 9)$:
  - Cusp 11: $f_{11} = +1/3$ of rising diurnal semi-arc ($\Delta\alpha \in (0, \pi)$ east of MC)
  - Cusp 12: $f_{12} = +2/3$ of rising diurnal semi-arc ($\Delta\alpha \in (0, \pi)$ east of MC)
  - Cusp 8: $f_8 = -2/3$ of setting diurnal semi-arc ($\Delta\alpha \in (-\pi, 0)$ west of MC), with $C_2 = (C_8 + 180^\circ) \bmod 360^\circ$
  - Cusp 9: $f_9 = -1/3$ of setting diurnal semi-arc ($\Delta\alpha \in (-\pi, 0)$ west of MC), with $C_3 = (C_9 + 180^\circ) \bmod 360^\circ$
  Solve the transcendental equation for meridian offset $\Delta\alpha = \alpha - \text{RAMC}$:
  $$g(\Delta\alpha) = \Delta\alpha - f_k \arccos\bigl(-\sin(\text{RAMC} + \Delta\alpha)\tan\phi\tan\varepsilon\bigr) = 0$$
  (via bisection on $[0, \pi f_k]$ or $[\pi f_k, 0]$ to $|\Delta\alpha| < 10^{-14}\text{ rad}$), then project the solution onto the ecliptic:
  $$C_k = \text{MC}(\text{RAMC} + \Delta\alpha_k, \varepsilon)$$
- **Polar Behavior (`hard_fail`)**:
  Undefined whenever $|\tan\phi\tan\varepsilon| \ge 1$, i.e. $|\phi| \ge 90^\circ - \varepsilon \approx 66.56^\circ$.
  - `HighLatPolicy::StrictError` (default): Returns `Err(EphemerisError::HighLatitudeCircumpolar { system: Placidus, latitude: phi, max_valid_latitude: 90.0 - eps })`.
  - `HighLatPolicy::FallbackPorphyry`: Returns Porphyry cusps with `PolarDegradationWarning::FallbackToPorphyry`.

#### 2. Koch (`'K'`, `HouseSystem::Koch`)
- **Classification**: Birthplace system (trisection of the MC's diurnal semi-arc projected via oblique ascension).
- **Equations**:
  $$\delta_{\text{MC}} = \arcsin(\sin\text{MC}\sin\varepsilon), \qquad \text{AD}_{\text{MC}} = \arcsin(\tan\delta_{\text{MC}}\tan\phi), \qquad \text{SDA}_{\text{MC}} = 90^\circ + \text{AD}_{\text{MC}}$$
  For $k \in \{11, 12, 2, 3\}$ with $f_k \in \{-2/3, -1/3, +1/3, +2/3\}$:
  $$C_k = \text{ASC}\bigl(\text{RAMC} + f_k \cdot \text{SDA}_{\text{MC}},\, \phi,\, \varepsilon\bigr)$$
- **Polar Behavior (`hard_fail`)**:
  Undefined when $|\tan\delta_{\text{MC}}\tan\phi| > 1$ (and globally restricted when $|\phi| \ge 90^\circ - \varepsilon$). Under `StrictError`, returns `Err(HighLatitudeCircumpolar)`.

#### 3. Porphyry (`'O'`, `HouseSystem::Porphyry`)
- **Classification**: Ecliptic quadrant trisection.
- **Equations**:
  Let $Q_4 = (\text{ASC} - \text{MC}) \bmod 360^\circ$ (4th quadrant arc) and $Q_1 = (\text{IC} - \text{ASC}) \bmod 360^\circ = 180^\circ - Q_4$ (1st quadrant arc):
  $$C_{11} = \left(\text{MC} + \frac{Q_4}{3}\right) \bmod 360^\circ, \quad C_{12} = \left(\text{MC} + \frac{2Q_4}{3}\right) \bmod 360^\circ$$
  $$C_2 = \left(\text{ASC} + \frac{Q_1}{3}\right) \bmod 360^\circ, \quad C_3 = \left(\text{ASC} + \frac{2Q_1}{3}\right) \bmod 360^\circ$$
- **Polar Behavior (`polar_safe`)**: Well-defined for all $|\phi| < 90^\circ$. At $|\phi| = 90^\circ$, returns `Err(SingularPole)` because the horizon coincides with the equator and `RAMC` is undefined.

#### 4. Regiomontanus (`'R'`, `HouseSystem::Regiomontanus`)
- **Classification**: Equatorial mundi trisection projected via North/South horizon poles.
- **Equations**:
  For $k \in \{11, 12, 2, 3\}$ with equatorial hour angle $h_k \in \{30^\circ, 60^\circ, 120^\circ, 150^\circ\}$:
  $$C_k = \operatorname{atan2}\bigl(\sin(\text{RAMC} + h_k),\, \cos(\text{RAMC} + h_k)\cos\varepsilon - \sin h_k \tan\phi \sin\varepsilon\bigr) \bmod 360^\circ$$
- **Polar Behavior (`soft_fail`)**: Mathematically continuous up to $|\phi| < 90^\circ$, though quadrants can become highly skewed above the Arctic/Antarctic circle.

#### 5. Campanus (`'C'`, `HouseSystem::Campanus`)
- **Classification**: Prime vertical trisection projected via North/South horizon poles.
- **Equations**:
  For $k \in \{11, 12, 2, 3\}$ with prime-vertical angle from the meridian $h_k \in \{30^\circ, 60^\circ, 120^\circ, 150^\circ\}$ (so prime-vertical altitude is $\theta_k = 90^\circ - h_k \in \{+60^\circ, +30^\circ, -30^\circ, -60^\circ\}$):
  $$\mathbf{p}_k = \cos\theta_k\,\mathbf{e} + \sin\theta_k\,\mathbf{z}, \qquad \mathbf{n}_k = \mathbf{n}_{\text{north}} \times \mathbf{p}_k$$
  The great circle with normal $\mathbf{n}_k = (n_{k,x}, n_{k,y}, n_{k,z})$ intersects the ecliptic at:
  $$\lambda_k = \operatorname{atan2}\bigl(-n_{k,x},\, n_{k,y}\cos\varepsilon + n_{k,z}\sin\varepsilon\bigr) \bmod 360^\circ$$
  choosing the branch $\lambda_k$ or $(\lambda_k + 180^\circ) \bmod 360^\circ$ such that $\mathbf{v}(\lambda_k) \cdot \mathbf{p}_k \ge 0$.
  Equivalently, in closed spherical-trigonometric form (let $H_k = 90^\circ - \theta_k = h_k$):
  $$\text{KO}_k = \arctan(\sin\phi\tan h_k), \quad \text{DN}_k = \arctan\bigl(\tan(\arccos(\sin\phi\sin h_k))\cos(\text{RAMC} + \text{KO}_k)\bigr)$$
  $$C_k = \operatorname{atan2}\bigl(\sin\text{DN}_k,\, \cos\text{DN}_k\cos\varepsilon\sin(\text{RAMC} + \text{KO}_k) - \sin\varepsilon\cos(\text{RAMC} + \text{KO}_k)\bigr) \bmod 360^\circ$$
- **Polar Behavior (`soft_fail`)**: Defined for all $|\phi| < 90^\circ$; above $|\phi| \ge 90^\circ - \varepsilon$, prime vertical circles can cross the ecliptic in retrograde order when the ecliptic pole passes through the zenith.

#### 6. Equal from ASC (`'A'`, `HouseSystem::Equal`)
- **Equations**: For $k \in 1..=12$:
  $$C_k = \bigl(\text{ASC} + 30^\circ(k - 1)\bigr) \bmod 360^\circ$$
- **Polar Behavior (`polar_safe`)**: Defined for all $|\phi| < 90^\circ$.

#### 7. Whole Sign (`'W'`, `HouseSystem::WholeSign`)
- **Equations**: For $k \in 1..=12$:
  $$C_k = \left(\left\lfloor \frac{\text{ASC}}{30^\circ} \right\rfloor \times 30^\circ + 30^\circ(k - 1)\right) \bmod 360^\circ$$
- **Polar Behavior (`polar_safe`)**: Defined for all $|\phi| < 90^\circ$.

#### 8. Alcabitius (`'B'`, `HouseSystem::Alcabitius`)
- **Classification**: Diurnal/nocturnal semi-arc of the Ascendant trisected along the equator and projected via hour circles (`MC` projection).
- **Equations**:
  $$\delta_{\text{ASC}} = \arcsin(\sin\varepsilon\sin\text{ASC}), \qquad s = \operatorname{clamp}(\tan\phi\tan\delta_{\text{ASC}}, -1, +1)$$
  $$\text{SDA}_{\text{ASC}} = 90^\circ + \arcsin(s), \qquad \text{SNA}_{\text{ASC}} = 180^\circ - \text{SDA}_{\text{ASC}}$$
  $$C_{11} = \text{MC}\!\left(\text{RAMC} + \frac{1}{3}\text{SDA}_{\text{ASC}}, \varepsilon\right), \quad C_{12} = \text{MC}\!\left(\text{RAMC} + \frac{2}{3}\text{SDA}_{\text{ASC}}, \varepsilon\right)$$
  $$C_2 = \text{MC}\!\left(\text{RAMC} + \text{SDA}_{\text{ASC}} + \frac{1}{3}\text{SNA}_{\text{ASC}}, \varepsilon\right), \quad C_3 = \text{MC}\!\left(\text{RAMC} + \text{SDA}_{\text{ASC}} + \frac{2}{3}\text{SNA}_{\text{ASC}}, \varepsilon\right)$$
- **Polar Behavior (`soft_fail`)**: When $|\tan\phi\tan\delta_{\text{ASC}}| > 1$, clamping $s \in [-1, 1]$ collapses $\text{SDA}_{\text{ASC}}$ to $0^\circ$ or $180^\circ$ (flagged under `StrictError` when $|\phi| \ge 90^\circ - \varepsilon$).

#### 9. Morinus (`'M'`, `HouseSystem::Morinus`)
- **Classification**: Equatorial sectors orthogonal to the ecliptic.
- **Equations**: For $k \in 1..=12$, let $\alpha_k = \text{RAMC} + 60^\circ + 30^\circ k$:
  $$C_k = \operatorname{atan2}\bigl(\sin\alpha_k\cos\varepsilon,\, \cos\alpha_k\bigr) \bmod 360^\circ$$
- **Polar Behavior (`polar_safe`)**: Independent of latitude $\phi$; valid even at $|\phi| = 90^\circ$ as long as local meridian longitude is specified.

#### 10. Polich-Page Topocentric (`'T'`, `HouseSystem::Topocentric`)
- **Classification**: Tangent-latitude cone projection of Placidus semi-arcs.
- **Equations**:
  For $k \in \{11, 12, 2, 3\}$ with $(h_k, f_k) \in \{(30^\circ, 1/3), (60^\circ, 2/3), (120^\circ, 2/3), (150^\circ, 1/3)\}$:
  $$\psi_k = \arctan(f_k \tan\phi), \qquad C_k = \text{ASC}(\text{RAMC} + h_k - 90^\circ,\, \psi_k,\, \varepsilon)$$
- **Polar Behavior (`soft_fail`)**: Defined for all $|\phi| < 90^\circ$ because $\arctan(f_k\tan\phi) \in (-90^\circ, +90^\circ)$ never diverges for $|\phi| < 90^\circ$.

#### 11. Krusinski-Pisa-Goelzer (`'U'`, `HouseSystem::Krusinski`)
- **Classification**: Great circle through Ascendant and Zenith trisected and projected via hour circles (`MC` projection).
- **Equations**:
  Let $\mathbf{v}_{\text{ASC}} = (\cos\text{ASC},\, \sin\text{ASC}\cos\varepsilon,\, \sin\text{ASC}\sin\varepsilon)$ and $\mathbf{z} = (\cos\phi\cos\text{RAMC},\, \cos\phi\sin\text{RAMC},\, \sin\phi)$.
  Using spherical linear interpolation $\operatorname{slerp}(\mathbf{a}, \mathbf{b}, t) = \frac{\sin((1-t)\Omega)}{\sin\Omega}\mathbf{a} + \frac{\sin(t\Omega)}{\sin\Omega}\mathbf{b}$ where $\Omega = \arccos(\mathbf{a}\cdot\mathbf{b})$:
  - Above horizon (Ascendant $\to$ Zenith): $\mathbf{p}_{12} = \operatorname{slerp}(\mathbf{v}_{\text{ASC}}, \mathbf{z}, 1/3)$, $\mathbf{p}_{11} = \operatorname{slerp}(\mathbf{v}_{\text{ASC}}, \mathbf{z}, 2/3)$
  - Below horizon (Ascendant $\to$ Nadir $-\mathbf{z}$): $\mathbf{p}_2 = \operatorname{slerp}(\mathbf{v}_{\text{ASC}}, -\mathbf{z}, 1/3)$, $\mathbf{p}_3 = \operatorname{slerp}(\mathbf{v}_{\text{ASC}}, -\mathbf{z}, 2/3)$
  For each $k \in \{11, 12, 2, 3\}$, compute right ascension $\alpha_k = \operatorname{atan2}(p_{k,y}, p_{k,x})$ and project onto the ecliptic:
  $$C_k = \text{MC}(\alpha_k, \varepsilon)$$
- **Polar Behavior (`polar_safe`)**: Defined for all $|\phi| < 90^\circ$.

#### 12. Meridian / Axial Rotation (`'X'`, `HouseSystem::Meridian`)
- **Classification**: Equal $30^\circ$ right ascension intervals from `RAMC` projected via hour circles.
- **Equations**: For $k \in 1..=12$, let $\alpha_k = \text{RAMC} + 60^\circ + 30^\circ k$:
  $$C_k = \text{MC}(\alpha_k, \varepsilon) = \operatorname{atan2}\bigl(\sin\alpha_k,\, \cos\alpha_k\cos\varepsilon\bigr) \bmod 360^\circ$$
  *(Note: $C_{10} = \text{MC}$, $C_4 = \text{IC}$, while $C_1 = \text{EqASC} = \text{MC}(\text{RAMC} + 90^\circ, \varepsilon)$ rather than the horizon Ascendant).*
- **Polar Behavior (`polar_safe`)**: Independent of latitude $\phi$.

#### 13. Vehlow Equal (`'V'`, `HouseSystem::Vehlow`)
- **Equations**: For $k \in 1..=12$:
  $$C_k = \bigl(\text{ASC} - 15^\circ + 30^\circ(k - 1)\bigr) \bmod 360^\circ$$
- **Polar Behavior (`polar_safe`)**: Defined for all $|\phi| < 90^\circ$.

#### 14. Pullen SD / Sinusoidal Delta (`'L'`, `HouseSystem::PullenSd` / `NeoPorphyry`)
- **Classification**: Sinusoidal quadrant partitioning by additive delta.
- **Equations**:
  For each quadrant width $Q \in \{Q_4, Q_1\}$ where $Q_4 = (\text{ASC} - \text{MC}) \bmod 360^\circ$ and $Q_1 = 180^\circ - Q_4$:
  - If $Q \ge 30^\circ$: the middle house has width $w_{\text{mid}} = \frac{Q - 30^\circ}{2}$ and the two flanking houses each have width $w_{\text{flank}} = \frac{Q + 30^\circ}{4}$ (note $2w_{\text{flank}} + w_{\text{mid}} = Q$).
  - If $Q < 30^\circ$ (extreme high-latitude quadrant compression): $w_{\text{mid}} = 0^\circ$ and $w_{\text{flank}} = \frac{Q}{2}$.
  Thus:
  $$C_{11} = (\text{MC} + w_{\text{flank}}(Q_4)) \bmod 360^\circ, \qquad C_{12} = (\text{ASC} - w_{\text{flank}}(Q_4)) \bmod 360^\circ$$
  $$C_2 = (\text{ASC} + w_{\text{flank}}(Q_1)) \bmod 360^\circ, \qquad C_3 = (\text{IC} - w_{\text{flank}}(Q_1)) \bmod 360^\circ$$
- **Polar Behavior (`soft_fail`)**: Continuous for all $|\phi| < 90^\circ$, with house 11 or 2 pinching to $0^\circ$ width when $Q \le 30^\circ$.

#### 15. Pullen SR / Sinusoidal Ratio (`'Q'`, `HouseSystem::PullenSr`)
- **Classification**: Sinusoidal quadrant partitioning by multiplicative ratio $r \ge 1$.
- **Equations**:
  Let $Q_4 = (\text{ASC} - \text{MC}) \bmod 360^\circ$, $Q_1 = 180^\circ - Q_4$, and let $q = \min(Q_4, Q_1) \in (0^\circ, 90^\circ]$ be the narrower quadrant and $Q_{\text{large}} = 180^\circ - q$ be the wider quadrant, with ratio $R_q = Q_{\text{large}} / q$.
  Adjacent house sizes oscillate multiplicatively by factor $r \ge 1$ such that the narrow quadrant has house widths $(r x,\, x,\, r x)$ and the wide quadrant has house widths $(r^3 x,\, r^4 x,\, r^3 x)$.
  Taking the ratio of total quadrant widths gives the quartic equation in $r \ge 1$:
  $$f(r) = r^3(r + 2) - R_q(2r + 1) = 0, \qquad f'(r) = 4r^3 + 6r^2 - 2R_q$$
  Solve $f(r) = 0$ via Newton-Raphson (starting from $r_0 = R_q^{1/4}$, converging in $\le 6$ iterations to $10^{-15}$), set $x = \frac{q}{2r + 1}$, and assign widths $(r x, x, r x)$ to the narrower quadrant and $(r^3 x, r^4 x, r^3 x)$ to the wider quadrant.
- **Polar Behavior (`polar_safe`)**: strictly positive house widths for all $q \in (0^\circ, 90^\circ]$.

#### 16 & 17. Sunshine — Treindl (`'I'`, `HouseSystem::SunshineTreindl`) & Makransky (`'i'`, `HouseSystem::SunshineMakransky`)
- **Classification**: Solar parallel of declination trisection (Bob Makransky 1988, *Primary Directions*; Alois Treindl general spherical solution).
- **Additional Input**: True apparent equatorial declination of the Sun $\delta_\odot(\text{JD}_{\text{TT}})$.
- **Equations**:
  Compute the solar ascensional difference parameter:
  $$s_\odot = \tan\delta_\odot \tan\phi$$
  - **Treindl formulation (`'I'`, `soft_fail`)**: Clamp $s_\odot \in [-1, +1]$ and set $k_v = \arcsin(\operatorname{clamp}(s_\odot, -1, 1))$.
    Cusp $1 = \text{ASC}$, Cusp $4 = \text{IC}$, Cusp $7 = \text{DSC}$, Cusp $10 = \text{MC}$.
    For each intermediate house $k \in \{2, 3, 5, 6, 8, 9, 11, 12\}$:
    $$a_k = \begin{cases}
    k_v + \text{RAMC} + \dfrac{\pi}{2} + \dfrac{k - 1}{3}\left(\dfrac{\pi}{2} - k_v\right) & \text{if } k < 8 \\[8pt]
    k_v + \text{RAMC} + \dfrac{\pi}{2} + \dfrac{k - 13}{3}\left(\dfrac{\pi}{2} + k_v\right) & \text{if } k \ge 8
    \end{cases}$$
    $$\begin{aligned}
    y_k &= \tan\delta_\odot\tan\phi\sin\text{RAMC} + \sin a_k \\
    x_k &= \cos\varepsilon\bigl(\tan\delta_\odot\tan\phi\cos\text{RAMC} + \cos a_k\bigr) + \sin\varepsilon\tan\phi\sin(\text{RAMC} - a_k) \\
    C_k &= \operatorname{atan2}(y_k, x_k) \bmod 360^\circ
    \end{aligned}$$
  - **Makransky formulation (`'i'`, `hard_fail`)**:
    Requires $|s_\odot| \le 1$ (returns `Err(HighLatitudeCircumpolar)` if $|\tan\delta_\odot\tan\phi| > 1$, i.e. when the Sun is circumpolar or never rises).
    Let $\text{AD}_\odot = \arcsin(\tan\delta_\odot\tan\phi)$, $\text{SDA}_\odot = 90^\circ + \text{AD}_\odot$, $\text{SNA}_\odot = 90^\circ - \text{AD}_\odot$.
    For all latitudes and solar declinations where the trisected meridian distances satisfy $m_k \le 90^\circ$ ($\frac{2}{3}\max(\text{SDA}_\odot, \text{SNA}_\odot) \le 90^\circ$, which holds whenever $|\tan\delta_\odot\tan\phi| \le 1/\sqrt{2}$ — covering all $|\phi| \le 58.48^\circ$ at solstice and higher latitudes nearer the equinoxes), `'I'` and `'i'` are **mathematically identical to `0.000000"`**.
    When $m_k > 90^\circ$ ($1/\sqrt{2} < |\tan\delta_\odot\tan\phi| \le 1$), Makransky's 1988 right-triangle `atan` expressions wrap across $m_k = 90^\circ$, producing a discontinuous branch jump in cusps 12/8 (summer) or 2/6 (winter); `'I'` (Treindl) uses the continuous 3D great-circle formula above and is the recommended default Sunshine variant.

#### 18. Carter Poli-Equatorial (`'F'`, `HouseSystem::Carter`)
- **Classification**: Equal $30^\circ$ right ascension steps from the Ascendant's right ascension $\alpha_{\text{ASC}}$, projected via hour circles.
- **Equations**:
  $$\alpha_{\text{ASC}} = \operatorname{atan2}\bigl(\sin\text{ASC}\cos\varepsilon,\, \cos\text{ASC}\bigr)$$
  For $k \in 1..=12$:
  $$C_k = \text{MC}\bigl(\alpha_{\text{ASC}} + 30^\circ(k - 1),\, \varepsilon\bigr)$$
- **Polar Behavior (`polar_safe`)**: Defined for all $|\phi| < 90^\circ$.

#### 19. Whole Sign Aries / Equal $0^\circ$ Aries (`'N'`, `HouseSystem::EqualAries`)
- **Equations**: For $k \in 1..=12$:
  $$C_k = 30^\circ(k - 1)$$
  (In Sidereal mode, $C_k = (30^\circ(k - 1) - A_{\text{true}}(t)) \bmod 360^\circ$ or $30^\circ(k-1)$ depending on whether sidereal $0^\circ$ Aries is the origin.)
- **Polar Behavior (`polar_safe`)**: Independent of latitude $\phi$.

#### 20. Sripati (`'S'`, `HouseSystem::Sripati`)
- **Classification**: Hindu Bhava Chalit system (Porphyrycusps are the *centers* of the 12 houses; Sripati cusps are the *midpoints* between consecutive Porphyry cusps).
- **Equations**:
  Let $P_1, \dots, P_{12}$ be the 12 Porphyry cusps (`'O'`). For each $k \in 1..=12$, let $P_{\text{prev}} = P_{((k - 2) \bmod 12) + 1}$:
  $$C_k = \left(P_{\text{prev}} + \frac{(P_k - P_{\text{prev}}) \bmod 360^\circ}{2}\right) \bmod 360^\circ$$
- **Polar Behavior (`polar_safe`)**: Defined for all $|\phi| < 90^\circ$.

#### 21. Equal from MC (`'D'`, `HouseSystem::EqualMc`)
- **Equations**: For $k \in 1..=12$:
  $$C_k = \bigl(\text{MC} + 90^\circ + 30^\circ(k - 1)\bigr) \bmod 360^\circ$$
- **Polar Behavior (`polar_safe`)**: Defined for all $|\phi| < 90^\circ$.

#### 22. APC / Apparent Pole System (`'Y'`, `HouseSystem::Apc`)
- **Classification**: Parallel of declination of the Ascendant trisected and projected via North/South horizon points.
- **Equations**:
  Compute the ascensional difference $k_v$ of the Ascendant directly from $(\text{RAMC}, \phi, \varepsilon)$:
  $$k_v = \arctan\!\left(\frac{\tan\phi\tan\varepsilon\cos\text{RAMC}}{1 + \tan\phi\tan\varepsilon\sin\text{RAMC}}\right), \qquad \tan\delta_{\text{ASC}} = \frac{\sin k_v}{\tan\phi}$$
  Cusp $1 = \text{ASC}$, Cusp $4 = \text{IC}$, Cusp $7 = \text{DSC}$, Cusp $10 = \text{MC}$.
  For each intermediate house $k \in \{2, 3, 5, 6, 8, 9, 11, 12\}$:
  $$a_k = \begin{cases}
  k_v + \text{RAMC} + \dfrac{\pi}{2} + \dfrac{k - 1}{3}\left(\dfrac{\pi}{2} - k_v\right) & \text{if } k < 8 \\[8pt]
  k_v + \text{RAMC} + \dfrac{\pi}{2} + \dfrac{k - 13}{3}\left(\dfrac{\pi}{2} + k_v\right) & \text{if } k \ge 8
  \end{cases}$$
  $$\begin{aligned}
  y_k &= \sin k_v \sin\text{RAMC} + \sin a_k \\
  x_k &= \cos\varepsilon\bigl(\sin k_v \cos\text{RAMC} + \cos a_k\bigr) + \sin\varepsilon\tan\phi\sin(\text{RAMC} - a_k) \\
  C_k &= \operatorname{atan2}(y_k, x_k) \bmod 360^\circ
  \end{aligned}$$
  *(Note: Cusps 5, 6, 8, 9 in APC are **not** $180^\circ$ opposite Cusps 11, 12, 2, 3 because the parallel of declination of the Ascendant is not a great circle when $\delta_{\text{ASC}} \ne 0$!)*
- **Polar Behavior (`soft_fail`)**: Defined for all $|\phi| < 90^\circ$.

#### 23. Azimuthal / Horizontal (`'H'`, `HouseSystem::Horizontal`)
- **Classification**: Horizon divided into twelve $30^\circ$ azimuth sectors from the meridian, projected onto the ecliptic via vertical circles through the Zenith and Nadir.
- **Equations**:
  Let $\mathbf{e} = (-\sin r, \cos r, 0)$, $\mathbf{z} = (\cos\phi\cos r, \cos\phi\sin r, \sin\phi)$, and define the **IC-ward horizon unit vector**:
  $$\mathbf{h}_{\text{IC}} = \operatorname{sgn}(\phi)\,\mathbf{n}_{\text{north}} = \operatorname{sgn}(\phi)\,(\mathbf{z} \times \mathbf{e})$$
  (where $\operatorname{sgn}(\phi) = +1$ for $\phi \ge 0$ and $-1$ for $\phi < 0$).
  Cusp $10 = \text{MC}$, Cusp $4 = \text{IC}$.
  For $k \in \{11, 12, 1, 2, 3\}$ with horizon angle from the MC $h_k \in \{30^\circ, 60^\circ, 90^\circ, 120^\circ, 150^\circ\}$ and $\theta_k = h_k - 90^\circ \in \{-60^\circ, -30^\circ, 0^\circ, +30^\circ, +60^\circ\}$:
  $$\mathbf{p}_k = \cos\theta_k\,\mathbf{e} + \sin\theta_k\,\mathbf{h}_{\text{IC}}, \qquad \mathbf{n}_k = \mathbf{z} \times \mathbf{p}_k$$
  $$\lambda_k = \operatorname{atan2}\bigl(-n_{k,x},\, n_{k,y}\cos\varepsilon + n_{k,z}\sin\varepsilon\bigr) \bmod 360^\circ$$
  choosing the branch $\lambda_k$ or $(\lambda_k + 180^\circ) \bmod 360^\circ$ such that $\mathbf{v}(\lambda_k) \cdot \mathbf{p}_k \ge 0$, and setting $C_{k \pm 6} = (C_k + 180^\circ) \bmod 360^\circ$.
  *(Note: For $k = 1$, $\theta_1 = 0^\circ \implies \mathbf{p}_1 = \mathbf{e}$, so $C_1$ is the intersection of the Prime Vertical with the ecliptic, i.e. $C_1 = (\text{Vertex} + 180^\circ) \bmod 360^\circ = \text{AntiVertex}$ in the Northern Hemisphere and $\text{Vertex}$ in the Southern Hemisphere).*
  - **Bug Fix vs `jayeshmepani/jpl-ephemeris`**: JME's `horizontal_cusp` omits $\operatorname{sgn}(\phi)$ on $\mathbf{n}_{\text{north}}$, causing Southern Hemisphere cusps $11, 12$ and $2, 3$ to swap across the Prime Vertical. Including $\operatorname{sgn}(\phi)$ makes `'H'` match L3 to `0.000000"` in both hemispheres.
- **Polar Behavior (`polar_safe`)**: Defined for all $|\phi| < 90^\circ$.

---

### 5.4 Summary of High-Latitude Behavior across All 23 Systems

| Code | Enum Variant | Name | Polar Class | Behavior for $\|\phi\| \ge 90^\circ - \varepsilon$ (`StrictError` Default) |
| :---: | :--- | :--- | :---: | :--- |
| `'P'` | `Placidus` | Placidus | `hard_fail` | Returns `Err(HighLatitudeCircumpolar)` when $\|\tan\phi\tan\varepsilon\| \ge 1$ |
| `'K'` | `Koch` | Koch | `hard_fail` | Returns `Err(HighLatitudeCircumpolar)` when $\|\tan\delta_{\text{MC}}\tan\phi\| > 1$ or $\|\phi\| \ge 90^\circ - \varepsilon$ |
| `'i'` | `SunshineMakransky` | Sunshine (Makransky) | `hard_fail` | Returns `Err(HighLatitudeCircumpolar)` when $\|\tan\delta_\odot\tan\phi\| > 1$ |
| `'R'` | `Regiomontanus` | Regiomontanus | `soft_fail` | Continuous up to $\|\phi\| < 90^\circ$; emits `HighLatitudeQuadrantSkew` warning above $66.56^\circ$ |
| `'C'` | `Campanus` | Campanus | `soft_fail` | Continuous up to $\|\phi\| < 90^\circ$; emits warning when ecliptic pole crosses zenith |
| `'B'` | `Alcabitius` | Alcabitius | `soft_fail` | Clamps $\tan\phi\tan\delta_{\text{ASC}} \in [-1, 1]$ (`SoftClamp` warning) or `Err` in strict mode |
| `'T'` | `Topocentric` | Polich-Page Topocentric | `soft_fail` | Tangent cone $\arctan(f\tan\phi)$ remains finite for all $\|\phi\| < 90^\circ$ |
| `'L'` | `PullenSd` | Pullen Sinusoidal Delta | `soft_fail` | Clamps middle house width to $0^\circ$ when quadrant $Q < 30^\circ$ |
| `'I'` | `SunshineTreindl` | Sunshine (Treindl) | `soft_fail` | Clamps $\tan\delta_\odot\tan\phi \in [-1, 1]$; continuous for all $\|\phi\| < 90^\circ$ |
| `'Y'` | `Apc` | APC (Apparent Pole) | `soft_fail` | Continuous for all $\|\phi\| < 90^\circ$ |
| `'O'` | `Porphyry` | Porphyry | `polar_safe` | Valid for all $\|\phi\| < 90^\circ$ |
| `'A'` | `Equal` | Equal from ASC | `polar_safe` | Valid for all $\|\phi\| < 90^\circ$ |
| `'W'` | `WholeSign` | Whole Sign | `polar_safe` | Valid for all $\|\phi\| < 90^\circ$ |
| `'M'` | `Morinus` | Morinus | `polar_safe` | Latitude-independent |
| `'U'` | `Krusinski` | Krusinski-Pisa-Goelzer | `polar_safe` | Valid for all $\|\phi\| < 90^\circ$ |
| `'X'` | `Meridian` | Meridian / Axial | `polar_safe` | Latitude-independent |
| `'V'` | `Vehlow` | Vehlow Equal | `polar_safe` | Valid for all $\|\phi\| < 90^\circ$ |
| `'Q'` | `PullenSr` | Pullen Sinusoidal Ratio | `polar_safe` | Strictly positive widths for all $\|\phi\| < 90^\circ$ |
| `'F'` | `Carter` | Carter Poli-Equatorial | `polar_safe` | Valid for all $\|\phi\| < 90^\circ$ |
| `'N'` | `EqualAries` | Whole Sign Aries | `polar_safe` | Latitude- and time-independent |
| `'S'` | `Sripati` | Sripati | `polar_safe` | Valid for all $\|\phi\| < 90^\circ$ |
| `'D'` | `EqualMc` | Equal from MC | `polar_safe` | Valid for all $\|\phi\| < 90^\circ$ |
| `'H'` | `Horizontal` | Azimuthal / Horizontal | `polar_safe` | Valid for all $\|\phi\| < 90^\circ$ |

*(At the exact geographic poles $|\phi| = 90^\circ$, all horizon-dependent systems return `Err(EphemerisError::SingularPole { latitude: phi })` because the horizon coincides with the celestial equator and local meridian `RAMC` is undefined).*

---

### 5.5 The 36 Gauquelin Sectors (`GauquelinSectors`)

Michel Gauquelin's 36 diurnal/nocturnal sectors divide the diurnal and nocturnal semi-arcs into 18 sectors each, numbered **clockwise** (in the direction of diurnal motion: Rise $\to$ MC $\to$ Set $\to$ IC $\to$ Rise):
- **Sector 1**: Starts at the eastern horizon (Ascendant) and rises toward the MC.
- **Sector 9 $\to$ 10 boundary**: Upper culmination (MC = Cusp of Sector 10).
- **Sector 18 $\to$ 19 boundary**: Western horizon (Descendant = Cusp of Sector 19).
- **Sector 27 $\to$ 28 boundary**: Lower culmination (IC = Cusp of Sector 28).

#### Ecliptic Cusps of the 36 Gauquelin Sectors
Using Placidus diurnal/nocturnal semi-arc subdivision into 18 intervals per semi-arc ($f_j = \frac{9 - j}{9}$ of the diurnal semi-arc for $j = 0..18$, and nocturnal semi-arc for $j = 18..36$):
- For each sector boundary $s \in 1..=36$:
  - $S_1 = \text{ASC}$, $S_{10} = \text{MC}$, $S_{19} = \text{DSC}$, $S_{28} = \text{IC}$.
  - For $s \in 2..=9$ (rising diurnal sectors between ASC and MC): Placidus root with diurnal fraction $f_s = \frac{10 - s}{9} \in \{8/9, 7/9, \dots, 1/9\}$ east of the meridian ($\Delta\alpha \in (0, \pi)$).
  - For $s \in 11..=18$ (setting diurnal sectors between MC and DSC): $(S_{s + 18} + 180^\circ) \bmod 360^\circ$, corresponding to nocturnal fraction $f = \frac{s - 10}{9}$ east of the IC.
  - For $s \in 20..=27$ (nocturnal sectors between DSC and IC): $(S_{s - 18} + 180^\circ) \bmod 360^\circ$.
  - For $s \in 29..=36$ (rising nocturnal sectors between IC and ASC): Placidus root with nocturnal fraction $f_s = -\frac{s - 28}{9} \in \{-1/9, -2/9, \dots, -8/9\}$ east of the meridian ($\Delta\alpha \in (-\pi, 0)$).
- **Body Sector Position (`gauquelin_sector_for_body`)**:
  Given a body's apparent equatorial coordinates $(\alpha, \delta)$, its diurnal semi-arc is $\text{SDA} = \arccos(\operatorname{clamp}(-\tan\phi\tan\delta, -1, 1))$ and nocturnal semi-arc is $\text{SNA} = \pi - \text{SDA}$. With meridian hour angle $H = (\text{RAMC} - \alpha) \bmod 2\pi \in [-\pi, \pi)$:
  - If the body is above the horizon ($|H| \le \text{SDA}$):
    $$\text{sector\_float} = 10.0 + 9.0 \frac{H}{\text{SDA}} \in [1.0, 19.0]$$
  - If the body is below the horizon ($|H| > \text{SDA}$):
    $$\text{sector\_float} = \begin{cases}
    19.0 + 9.0 \dfrac{H - \text{SDA}}{\text{SNA}} \in (19.0, 28.0] & \text{if } H > 0 \\[8pt]
    28.0 + 9.0 \dfrac{\pi + H}{\text{SNA}} \in (28.0, 37.0) & \text{if } H < 0
    \end{cases}$$

---

## 6. Solar & Lunar Eclipses (`eclipses`)

### 6.1 Fundamental Constants & Physical Radii

All eclipse calculations evaluate apparent geocentric (or topocentric) positions of the Sun and Moon from the active SPK kernel (`DE440s` / `DE440`) using `AberrationMode::Apparent` (light-time + stellar aberration for the Sun; light-time for the Moon):

| Symbol | Constant Name | Value | Reference |
| :--- | :--- | :--- | :--- |
| $R_\odot$ | IAU Nominal Solar Photospheric Radius | $696\,000.0\text{ km}$ | IAU 2015 Resolution B3 (`au = 149597870.700 km`, $\sigma_\odot = 959.63''$ at $1\text{ AU}$) |
| $R_E$ | WGS84 / IERS Earth Equatorial Radius | $6\,378.1366\text{ km}$ | IERS Conventions (2010) |
| $f_E$ | WGS84 Earth Flattening | $1 / 298.257223563$ | WGS84 |
| $R_{\text{Moon}, k}$ | Lunar Mean Radius Ratio (Solar Eclipses) | $k_{\text{mean}} = 0.2725076\,R_E = 1\,738.1\text{ km}$ | IAU 1982 / Espenak & Meeus (NASA TP-2006-214141) |
| $k_{\text{umbra}}$ | Lunar Valley-Limb Ratio (Total/Annular Centrality) | $k_{\text{valley}} = 0.2722810\,R_E = 1\,736.65\text{ km}$ | Danjon / Espenak (used for eclipse magnitude & duration) |
| $\alpha_{\text{atm}}$ | Earth Atmospheric Shadow Enlargement (Chauvenet / Danjon) | $1.0100$ ($1\%$ enlargement) or Danjon $1/85$ rule | Espenak & Meeus (2009), *Five Millennium Canon of Lunar Eclipses* |

---

### 6.2 Lunar Eclipses (Penumbral, Partial, Total)

1. **Coarse Search & Refine (Full Moon / Opposition)**:
   Step across the search interval in $\Delta t = 1.0\text{ day}$ increments monitoring the ecliptic elongation $\Delta\lambda(t) = (\lambda_{\text{Moon}}(t) - \lambda_\odot(t)) \bmod 360^\circ$. Whenever $\Delta\lambda(t)$ crosses $180^\circ$ and $|\beta_{\text{Moon}}| < 1.65^\circ$, bracket and minimize the geocentric angular separation $\psi(t)$ between the Moon's apparent center $\hat{\mathbf{u}}_{\text{Moon}}(t)$ and the anti-solar (shadow-axis) unit vector $\hat{\mathbf{u}}_{\text{shadow}}(t) = -\hat{\mathbf{u}}_\odot(t)$ via Brent's method (`brent_minimize`) to $|\Delta t| < 10^{-8}\text{ day}$ ($< 1\text{ ms}$).

2. **Equatorial Parallax & Shadow Radii at Lunar Distance**:
   At any instant $t$ with geocentric distances $d_\odot = \|\mathbf{r}_\odot\|$ and $d_{\text{Moon}} = \|\mathbf{r}_{\text{Moon}}\|$:
   $$\begin{aligned}
   \pi_\odot &= \arcsin\!\left(\frac{R_E}{d_\odot}\right), & s_\odot &= \arcsin\!\left(\frac{R_\odot}{d_\odot}\right) \\
   \pi_{\text{Moon}} &= \arcsin\!\left(\frac{R_E}{d_{\text{Moon}}}\right), & s_{\text{Moon}} &= \arcsin\!\left(\frac{k_{\text{valley}}\,R_E}{d_{\text{Moon}}}\right)
   \end{aligned}$$
   Using the standard Chauvenet/Espenak $1\%$ atmospheric shadow enlargement factor $f_{\text{sh}} = 1.01$ (with mean oblateness factor $\pi_1 = 0.998340\,\pi_{\text{Moon}}$ per Espenak & Meeus 2009):
   $$\begin{aligned}
   \rho_{\text{umbra}}(t) &= f_{\text{sh}}\bigl(\pi_1 + \pi_\odot - s_\odot\bigr) && \text{(Angular radius of Earth's umbra at lunar distance)} \\
   \rho_{\text{penumbra}}(t) &= f_{\text{sh}}\bigl(\pi_1 + \pi_\odot + s_\odot\bigr) && \text{(Angular radius of Earth's penumbra at lunar distance)}
   \end{aligned}$$

3. **Classification & Magnitudes at Greatest Eclipse ($t_{\max}$)**:
   Let $\gamma_{\text{ang}} = \psi(t_{\max})$ be the minimum angular distance between the Moon's center and the shadow axis:
   $$\begin{aligned}
   M_{\text{umbra}} &= \frac{\rho_{\text{umbra}}(t_{\max}) + s_{\text{Moon}}(t_{\max}) - \gamma_{\text{ang}}}{2\,s_{\text{Moon}}(t_{\max})} && \text{(Umbral magnitude)} \\
   M_{\text{penumbra}} &= \frac{\rho_{\text{penumbra}}(t_{\max}) + s_{\text{Moon}}(t_{\max}) - \gamma_{\text{ang}}}{2\,s_{\text{Moon}}(t_{\max})} && \text{(Penumbral magnitude)}
   \end{aligned}$$
   - If $M_{\text{umbra}} \ge 1.0$: **`LunarEclipseKind::Total`** ($\gamma_{\text{ang}} \le \rho_{\text{umbra}} - s_{\text{Moon}}$).
   - Else if $M_{\text{umbra}} > 0.0$: **`LunarEclipseKind::Partial`** ($\rho_{\text{umbra}} - s_{\text{Moon}} < \gamma_{\text{ang}} < \rho_{\text{umbra}} + s_{\text{Moon}}$).
   - Else if $M_{\text{penumbra}} > 0.0$: **`LunarEclipseKind::Penumbral`** ($\rho_{\text{umbra}} + s_{\text{Moon}} \le \gamma_{\text{ang}} < \rho_{\text{penumbra}} + s_{\text{Moon}}$).
   - Else: No lunar eclipse at this syzygy.

4. **Contact Times ($P_1, U_1, U_2, t_{\max}, U_3, U_4, P_4$)**:
   Using `brent_zero` on $f_C(t) = \psi(t) - R_C(t)$ over $[t_{\max} - 0.30\text{ d},\, t_{\max}]$ (ingress) and $[t_{\max},\, t_{\max} + 0.30\text{ d}]$ (egress):
   - **Penumbral contacts $P_1, P_4$**: $R_{P}(t) = \rho_{\text{penumbra}}(t) + s_{\text{Moon}}(t)$.
   - **Partial (umbral) contacts $U_1, U_4$** (when $M_{\text{umbra}} > 0$): $R_{U,\text{ext}}(t) = \rho_{\text{umbra}}(t) + s_{\text{Moon}}(t)$.
   - **Totality contacts $U_2, U_3$** (when $M_{\text{umbra}} \ge 1$): $R_{U,\text{int}}(t) = \rho_{\text{umbra}}(t) - s_{\text{Moon}}(t)$.

---

### 6.3 Solar Eclipses (Global & Local via Fundamental Plane / Besselian Elements)

1. **Fundamental Plane Reference Frame (Besselian Coordinates)**:
   At any instant $t$, let $\mathbf{r}_\odot(t)$ and $\mathbf{r}_{\text{Moon}}(t)$ be the apparent geocentric equatorial vectors of the Sun and Moon, and define the shadow-axis vector from Moon to Sun:
   $$\mathbf{g}(t) = \mathbf{r}_\odot(t) - \mathbf{r}_{\text{Moon}}(t), \qquad \hat{\mathbf{k}}(t) = \frac{\mathbf{g}(t)}{\|\mathbf{g}(t)\|}$$
   With right ascension $a(t)$ and declination $d(t)$ of $\hat{\mathbf{k}}(t)$, construct the orthonormal Besselian triad $(\hat{\mathbf{i}}, \hat{\mathbf{j}}, \hat{\mathbf{k}})$ where $\hat{\mathbf{k}}$ points toward the shadow axis (Sunward), $\hat{\mathbf{i}} = (-\sin a, \cos a, 0)$ points eastward in the Fundamental Plane, and $\hat{\mathbf{j}} = \hat{\mathbf{k}} \times \hat{\mathbf{i}}$ points northward.
   To account for the WGS84 oblateness of the Earth ($f_E = 1/298.257223563$, $e^2 = 2f_E - f_E^2$), scale $z$-coordinates into the equivalent spherical Earth of radius $R_E$ via $\rho_z = 1 / \sqrt{1 - e^2}$.
   In units of $R_E$, the coordinates of the shadow axis in the Fundamental Plane are:
   $$x(t) = \frac{\mathbf{r}_{\text{Moon}}(t) \cdot \hat{\mathbf{i}}}{R_E}, \qquad y(t) = \frac{\mathbf{r}_{\text{Moon}}(t) \cdot \hat{\mathbf{j}}}{R_E}, \qquad z(t) = \frac{\mathbf{r}_{\text{Moon}}(t) \cdot \hat{\mathbf{k}}}{R_E}$$

2. **Umbral and Penumbral Cone Angles ($f_1, f_2$) and Fundamental Plane Radii ($l_1, l_2$)**:
   Let $D_{MS} = \|\mathbf{g}(t)\| / R_E$ be the Moon–Sun distance in Earth radii, $r_S = R_\odot / R_E = 109.12278$, and $r_M = k_{\text{valley}} = 0.2722810$:
   $$\sin f_1 = \frac{r_S + r_M}{D_{MS}} \quad \text{(Penumbral cone)}, \qquad \sin f_2 = \frac{r_S - r_M}{D_{MS}} \quad \text{(Umbral cone)}$$
   $$c_1(t) = z(t) + \frac{r_M}{\sin f_1}, \qquad c_2(t) = z(t) - \frac{r_M}{\sin f_2}$$
   $$l_1(t) = c_1(t)\tan f_1 \quad \text{(Penumbral radius on Fundamental Plane)}, \qquad l_2(t) = c_2(t)\tan f_2 \quad \text{(Umbral radius on Fundamental Plane)}$$
   *(Convention: $l_2(t) < 0$ when the umbral vertex reaches the surface of the Earth — **Total**; $l_2(t) > 0$ when the umbral vertex falls short of the fundamental plane — **Annular**).*

3. **Global Solar Eclipse Classification & Gamma ($\gamma$)**:
   At greatest eclipse $t_{\max}$ where $\gamma_{\text{abs}} = \sqrt{x(t_{\max})^2 + y(t_{\max})^2}$ is minimized ($\gamma = \operatorname{sgn}(y(t_{\max}))\sqrt{x^2 + y^2}$):
   - Let $\zeta_0 = \sqrt{\max(0,\, 1 - \gamma^2)}$ be the height of the surface point closest to the shadow axis above the Fundamental Plane, with surface shadow radii $L_1 = l_1 - \zeta_0\tan f_1$ and $L_2 = l_2 - \zeta_0\tan f_2$.
   - **Central vs Non-Central**: Central if $|\gamma| \le 1.0$ (or $|\gamma| \le 0.9972 + |L_2|$ grazing limit).
   - **Type (`SolarEclipseKind`)**:
     - If $|\gamma| > 1.0 + L_1$: No solar eclipse.
     - If $|\gamma| > 1.0 + |L_2|$: **`SolarEclipseKind::Partial`** with magnitude $M = \frac{1.0 + L_1 - |\gamma|}{L_1 + L_2}$.
     - If central/umbral ($|\gamma| \le 1.0 + |L_2|$):
       - Inspect $L_2(\zeta)$ along the central path ($x^2(t) + y^2(t) \le 1$):
         - If $L_2 < 0$ everywhere along the central path: **`SolarEclipseKind::Total`**.
         - If $L_2 > 0$ everywhere along the central path: **`SolarEclipseKind::Annular`**.
         - If $L_2$ changes sign along the central path ($L_2 < 0$ near $t_{\max}$ where $\zeta \approx 1$, and $L_2 > 0$ near the sunrise/sunset limbs where $\zeta \approx 0$): **`SolarEclipseKind::Hybrid`** (Annular-Total).

4. **Local (Topocentric) Solar Eclipse Circumstances**:
   For an observer at geodetic $(\lambda_{\text{obs}}, \phi, h)$ with Fundamental Plane coordinates $(\xi(t), \eta(t), \zeta(t))$:
   - Observer shadow-plane radii: $L_1'(t) = l_1(t) - \zeta(t)\tan f_1$, $L_2'(t) = l_2(t) - \zeta(t)\tan f_2$.
   - Shadow-axis distance in the observer's plane: $\Delta(t) = \sqrt{(x(t) - \xi(t))^2 + (y(t) - \eta(t))^2}$.
   - Maximum local eclipse $t_{\max,\text{loc}}$ minimizes $\Delta(t)$ (subject to topocentric solar altitude $> 0^\circ$).
   - Local magnitude: $M_{\text{loc}} = \frac{L_1'(t_{\max,\text{loc}}) - \Delta(t_{\max,\text{loc}})}{L_1'(t_{\max,\text{loc}}) + L_2'(t_{\max,\text{loc}})}$.
   - Local contacts $C_1, C_4$ solve $\Delta(t) = L_1'(t)$; inner contacts $C_2, C_3$ solve $\Delta(t) = |L_2'(t)|$.
