#!/usr/bin/env python3
"""
Level 1B Golden Vector Generator for `lxephemeris`.

Generates independent scientific reference fixtures in `testdata/golden/` using:
  - `jplephem` + `skyfield` (reading `data/de440s.bsp`)
  - `pyerfa` (IAU SOFA / ERFA v2.0.1+)
  - `numpy` / `scipy`-free pure numerical routines (Brent root finder + Thomas C2 cubic spline)

Fixtures produced:
  1. `testdata/golden/l1b_planets_de440s.json`
  2. `testdata/golden/l1b_ayanamsas_12.json`
  3. `testdata/golden/l1b_lunar_points_7.json`
  4. `testdata/golden/l1b_uranians_9.json`
  5. `testdata/golden/l1b_stars_50.json`
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import erfa
import numpy as np
from jplephem.spk import SPK

ROOT = Path(__file__).resolve().parent.parent
DE440S_PATH = ROOT / "data" / "de440s.bsp"
GOLDEN_DIR = ROOT / "testdata" / "golden"

AU_KM = 149597870.700
C_KM_S = 299792.458
C_AU_DAY = (C_KM_S * 86400.0) / AU_KM
DAS2R = 4.848136811095359935899141e-6
DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi

# DE440 header GM values in km^3/s^2 (Park et al. 2021, Table 2)
GM_SUN_KM3_S2 = 132712440041.279419
GM_EARTH_KM3_S2 = 398600.435507
GM_MOON_KM3_S2 = 4902.800118
MU_EM_KM3_S2 = GM_EARTH_KM3_S2 + GM_MOON_KM3_S2  # 403503.235625 km^3/s^2


def tt_to_tdb_jd(jd_tt: float) -> float:
    """Convert JD(TT) to JD(TDB) using SOFA/ERFA dtdb at geocenter."""
    g = math.radians((357.53 + 0.98560028 * (jd_tt - 2451545.0)) % 360.0)
    dt_sec = erfa.dtdb(jd_tt, 0.0, 0.0, 0.0, 0.0, g)
    return jd_tt + dt_sec / 86400.0


def spk_barycentric_state_km_s(spk: SPK, target_id: int, jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
    """Return barycentric (SSB=0) position [km] and velocity [km/s] in ICRF from DE440s."""
    if target_id == 0:
        return np.zeros(3), np.zeros(3)
    if target_id in (199, 299):
        center = 1 if target_id == 199 else 2
        p1, v1 = spk[0, center].compute_and_differentiate(jd_tdb)
        p2, v2 = spk[center, target_id].compute_and_differentiate(jd_tdb)
        return p1 + p2, (v1 + v2) / 86400.0
    if target_id in (301, 399):
        p1, v1 = spk[0, 3].compute_and_differentiate(jd_tdb)
        p2, v2 = spk[3, target_id].compute_and_differentiate(jd_tdb)
        return p1 + p2, (v1 + v2) / 86400.0
    p, v = spk[0, target_id].compute_and_differentiate(jd_tdb)
    return p, v / 86400.0


def icrf_to_ecliptic_matrix(jd_tt: float, nutation: bool = True) -> np.ndarray:
    """Return 3x3 rotation matrix from ICRF to mean or true ecliptic of date."""
    eps0 = float(erfa.obl06(jd_tt, 0.0))
    if nutation:
        dpsi, deps = erfa.nut00a(jd_tt, 0.0)
        eps = eps0 + float(deps)
        r_eq = erfa.num06a(jd_tt, 0.0) @ erfa.bp06(jd_tt, 0.0)[2]
    else:
        eps = eps0
        r_eq = erfa.bp06(jd_tt, 0.0)[2]
    ce, se = math.cos(eps), math.sin(eps)
    r_x = np.array([[1.0, 0.0, 0.0], [0.0, ce, se], [0.0, -se, ce]])
    return r_x @ r_eq


def icrf_to_equatorial_matrix(jd_tt: float, nutation: bool = True) -> np.ndarray:
    if nutation:
        return erfa.num06a(jd_tt, 0.0) @ erfa.bp06(jd_tt, 0.0)[2]
    return erfa.bp06(jd_tt, 0.0)[2]


def cart_to_spherical(u: np.ndarray) -> tuple[float, float, float]:
    """Return (lon_deg in [0, 360), lat_deg in [-90, 90], r)."""
    x, y, z = float(u[0]), float(u[1]), float(u[2])
    r = math.hypot(x, y, z)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(max(-1.0, min(1.0, z / r)))) if r > 0.0 else 0.0
    return lon, lat, r


def cart_state_to_spherical_rates(u: np.ndarray, du: np.ndarray) -> tuple[float, float, float]:
    """Convert 3D position u and velocity du (per day) to (dlon_deg_day, dlat_deg_day, dr_day)."""
    x, y, z = float(u[0]), float(u[1]), float(u[2])
    dx, dy, dz = float(du[0]), float(du[1]), float(du[2])
    r2 = x * x + y * y + z * z
    r = math.sqrt(r2)
    rho2 = x * x + y * y
    rho = math.sqrt(rho2)
    dr = (x * dx + y * dy + z * dz) / r if r > 0.0 else 0.0
    dlon = math.degrees((x * dy - y * dx) / rho2) if rho2 > 0.0 else 0.0
    dlat = math.degrees((rho * dz - z * (x * dx + y * dy) / rho) / r2) if (r2 > 0.0 and rho > 0.0) else 0.0
    return dlon, dlat, dr


def planet_geocentric_vector_au(
    spk: SPK, target_id: int, jd_tt: float, apparent: bool = True, nutation: bool = True
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute geocentric vector (in target ecliptic and equatorial frames, in AU)
    with light-time iteration, solar gravitational light deflection (for planets),
    and relativistic stellar aberration (when apparent=True).
    """
    jd_tdb = tt_to_tdb_jd(jd_tt)
    p_earth_km, v_earth_kms = spk_barycentric_state_km_s(spk, 399, jd_tdb)
    p_earth_au = p_earth_km / AU_KM

    # Light-time iteration
    dt_days = 0.0
    for _ in range(6):
        p_targ_km, _ = spk_barycentric_state_km_s(spk, target_id, jd_tdb - dt_days)
        rho_au = (p_targ_km - p_earth_km) / AU_KM
        dist_au = float(np.linalg.norm(rho_au))
        dt_days = dist_au / C_AU_DAY

    u_unit = rho_au / dist_au
    if apparent:
        # Relativistic solar light deflection (except for Sun=10 and Moon=301)
        if target_id not in (10, 301):
            p_sun_km, _ = spk_barycentric_state_km_s(spk, 10, jd_tdb - dt_days)
            q_vec = (p_targ_km - p_sun_km) / AU_KM
            e_vec = (p_earth_km - p_sun_km) / AU_KM
            em = float(np.linalg.norm(e_vec))
            q_unit = q_vec / float(np.linalg.norm(q_vec))
            e_unit = e_vec / em
            u_unit = np.asarray(erfa.ld(1.0, u_unit, q_unit, e_unit, em, 0.0), dtype=float)
            u_unit = u_unit / float(np.linalg.norm(u_unit))

        # Relativistic aberration via erfa.ab
        p_sun_now_km, _ = spk_barycentric_state_km_s(spk, 10, jd_tdb)
        bm1 = math.sqrt(1.0 - min(float(np.dot(v_earth_kms, v_earth_kms)) / (C_KM_S * C_KM_S), 0.999999))
        v_ab = v_earth_kms / C_KM_S
        s_dist_au = float(np.linalg.norm((p_earth_km - p_sun_now_km) / AU_KM))
        u_unit = np.asarray(erfa.ab(u_unit, v_ab, s_dist_au, bm1), dtype=float)

    vec_icrf_au = u_unit * dist_au
    r_ecl = icrf_to_ecliptic_matrix(jd_tt, nutation=nutation)
    r_equ = icrf_to_equatorial_matrix(jd_tt, nutation=nutation)
    return r_ecl @ vec_icrf_au, r_equ @ vec_icrf_au


def generate_l1b_planets(spk: SPK) -> dict:
    bodies = [
        (10, "Sun", "planet"),
        (301, "Moon", "planet"),
        (199, "Mercury", "planet"),
        (299, "Venus", "planet"),
        (4, "Mars_BC", "barycenter"),
        (5, "Jupiter_BC", "barycenter"),
        (6, "Saturn_BC", "barycenter"),
        (7, "Uranus_BC", "barycenter"),
        (8, "Neptune_BC", "barycenter"),
        (9, "Pluto_BC", "barycenter"),
    ]
    epochs = [
        ("1855-01-01T00:00:00", 2398585.5),
        ("1900-01-01T12:00:00", 2415021.0),
        ("1950-01-01T00:00:00", 2433282.5),
        ("1956-03-21T00:00:00", 2435553.5),
        ("1980-06-15T12:00:00", 2444406.0),
        ("2000-01-01T12:00:00", 2451545.0),
        ("2020-03-20T03:50:00", 2458928.659722222),  # Sun near 0 deg Aries wrap
        ("2024-04-08T18:18:00", 2460409.2625),       # Solar eclipse syzygy
        ("2026-06-15T00:00:00", 2461206.5),
        ("2050-01-01T00:00:00", 2469807.5),
        ("2100-01-01T00:00:00", 2488069.5),
        ("2145-01-01T00:00:00", 2504505.5),
    ]
    cases = []
    for label, jd_tt in epochs:
        for tid, tname, tkind in bodies:
            h = (1.0 / 24.0) if tid == 301 else 0.25
            # 4-point 4th-order Cartesian stencil [-2h, -h, +h, +2h]
            ecl_m2, _ = planet_geocentric_vector_au(spk, tid, jd_tt - 2 * h, apparent=True, nutation=True)
            ecl_m1, _ = planet_geocentric_vector_au(spk, tid, jd_tt - h, apparent=True, nutation=True)
            ecl_0, equ_0 = planet_geocentric_vector_au(spk, tid, jd_tt, apparent=True, nutation=True)
            ecl_p1, _ = planet_geocentric_vector_au(spk, tid, jd_tt + h, apparent=True, nutation=True)
            ecl_p2, _ = planet_geocentric_vector_au(spk, tid, jd_tt + 2 * h, apparent=True, nutation=True)
            decl = (ecl_m2 - 8.0 * ecl_m1 + 8.0 * ecl_p1 - ecl_p2) / (12.0 * h)

            lon_deg, lat_deg, dist_au = cart_to_spherical(ecl_0)
            dlon, dlat, ddist = cart_state_to_spherical_rates(ecl_0, decl)
            ra_deg, dec_deg, _ = cart_to_spherical(equ_0)

            # Astrometric ICRF (no deflection, no aberration, no precession)
            jd_tdb = tt_to_tdb_jd(jd_tt)
            p_earth_km, _ = spk_barycentric_state_km_s(spk, 399, jd_tdb)
            dt = 0.0
            for _ in range(6):
                p_targ_km, _ = spk_barycentric_state_km_s(spk, tid, jd_tdb - dt)
                rho_km = p_targ_km - p_earth_km
                dt = float(np.linalg.norm(rho_km)) / (C_KM_S * 86400.0)
            ra_icrf, dec_icrf, dist_icrf_km = cart_to_spherical(rho_km)

            cases.append({
                "epoch_label": label,
                "jd_tt": jd_tt,
                "jd_tdb": jd_tdb,
                "target_id": tid,
                "target_name": tname,
                "target_kind": tkind,
                "apparent_ecliptic_true_of_date": {
                    "lon_deg": lon_deg,
                    "lat_deg": lat_deg,
                    "dist_au": dist_au,
                    "lon_speed_deg_day": dlon,
                    "lat_speed_deg_day": dlat,
                    "dist_speed_au_day": ddist,
                },
                "apparent_equatorial_true_of_date": {
                    "ra_deg": ra_deg,
                    "dec_deg": dec_deg,
                },
                "astrometric_icrf": {
                    "ra_deg": ra_icrf,
                    "dec_deg": dec_icrf,
                    "dist_au": dist_icrf_km / AU_KM,
                },
            })

    return {
        "schema_version": "1.0.0",
        "level": "L1B",
        "source": "jplephem + pyerfa (IAU 2006/2000A) + 4th-order Cartesian stencil",
        "versions": {
            "erfa": erfa.__version__,
            "numpy": np.__version__,
        },
        "ephemeris": "DE440s (de440s.bsp)",
        "time_scale": "TT",
        "tolerance_arcsec": 0.01,
        "count": len(cases),
        "vectors": cases,
    }


def iae_general_precession_arcsec(jd_t0: float, jd_t: float) -> tuple[float, float]:
    """Lieske (1977) Table 5 / IAE (1989) p. 525 two-epoch precession in arcsec and rate in arcsec/day."""
    big_t = (jd_t0 - 2451545.0) / 36525.0
    t = (jd_t - jd_t0) / 36525.0
    p_arcsec = (
        (5029.0966 + 2.22226 * big_t - 0.000042 * big_t * big_t) * t
        + (1.11161 - 0.000127 * big_t) * t * t
        - 0.000113 * t * t * t
    )
    rate_arcsec_day = (
        (5029.0966 + 2.22226 * big_t - 0.000042 * big_t * big_t)
        + 2.0 * (1.11161 - 0.000127 * big_t) * t
        - 3.0 * 0.000113 * t * t
    ) / 36525.0
    return p_arcsec, rate_arcsec_day


def star_mean_ecliptic_lon_no_ab(spk: SPK, star: dict, jd_tt: float) -> float:
    """Propagate Hipparcos-2 J1991.25 star via pmpx (with geocentric parallax, no aberration) to mean ecliptic of date."""
    jd_tdb = tt_to_tdb_jd(jd_tt)
    p_earth_km, _ = spk_barycentric_state_km_s(spk, 399, jd_tdb)
    p_eb_au = p_earth_km / AU_KM

    ra1 = math.radians(star["ra_1991_deg"])
    de1 = math.radians(star["de_1991_deg"])
    pmr = (star["pm_ra_cosdec_mas_yr"] * 1e-3 * DAS2R) / math.cos(de1)
    pmd = star["pm_dec_mas_yr"] * 1e-3 * DAS2R
    px = max(star["plx_mas"] * 1e-3, 0.0)
    pmt = (jd_tdb - 2448349.0625) / 365.25
    u_icrs = np.asarray(erfa.pmpx(ra1, de1, pmr, pmd, px, star["rv_kms"], pmt, p_eb_au), dtype=float)
    r_ecl_mean = icrf_to_ecliptic_matrix(jd_tt, nutation=False)
    u_ecl = r_ecl_mean @ u_icrs
    return math.degrees(math.atan2(float(u_ecl[1]), float(u_ecl[0]))) % 360.0


def generate_l1b_ayanamsas(spk: SPK, stars_by_hip: dict[int, dict]) -> dict:
    # 9 Epoch-defined models + 3 Star-anchored models
    dpsi_1956, _ = erfa.nut80(2435553.5, 0.0)
    lahiri_mean_1956 = (23.0 * 3600.0 + 15.0 * 60.0 + 0.658 - math.degrees(float(dpsi_1956)) * 3600.0) / 3600.0

    epoch_models = [
        ("fagan_bradley", 0, 2433282.42345905, 24.0 + 2.0 / 60.0 + 31.36 / 3600.0),
        ("lahiri", 1, 2435553.5, lahiri_mean_1956),
        ("raman", 3, 2415020.0, 21.014210000000000),
        ("ushashashi", 4, 2415020.0, 18.660960000000000),
        ("krishnamurti", 5, 2415020.0, 22.363659000000000),
        ("yukteshwar", 7, 2415020.0, 21.082222222222222),
        ("jn_bhasin", 8, 2415020.0, 21.365555555555556),
        ("suryasiddhanta", 21, 1903396.8770698, 0.0),
        ("ss_citra", 26, 1903396.8770698, 2.110704444444444),
    ]
    star_models = [
        ("true_citra", 27, 65474, 180.0),
        ("true_revati", 28, 5737, 359.0 + 50.0 / 60.0),
        ("true_pushya", 29, 42911, 106.0),
    ]
    epochs = [
        ("1855-01-01T00:00:00", 2398585.5),
        ("1900-01-01T00:00:00", 2415020.0),
        ("1949-12-31T22:09:47", 2433282.42345905),
        ("1956-03-21T00:00:00", 2435553.5),
        ("1980-06-15T12:00:00", 2444406.0),
        ("2000-01-01T00:00:00", 2451544.500738756),
        ("2000-01-01T12:00:00", 2451545.0),
        ("2026-06-15T00:00:00", 2461206.5),
        ("2050-01-01T00:00:00", 2469807.5),
        ("2100-01-01T00:00:00", 2488069.5),
    ]
    vectors = []
    for label, jd_tt in epochs:
        dpsi_rad, _ = erfa.nut00a(jd_tt, 0.0)
        dpsi_deg = math.degrees(float(dpsi_rad))
        for code, legacy_id, t0, a0_deg in epoch_models:
            p_arcsec, rate_arcsec_day = iae_general_precession_arcsec(t0, jd_tt)
            mean_deg = (a0_deg + p_arcsec / 3600.0) % 360.0
            true_deg = (mean_deg + dpsi_deg) % 360.0
            vectors.append({
                "epoch_label": label,
                "jd_tt": jd_tt,
                "ayanamsa_code": code,
                "legacy_id": legacy_id,
                "category": "epoch_defined",
                "ayanamsa_mean_deg": mean_deg,
                "ayanamsa_true_deg": true_deg,
                "rate_deg_day": rate_arcsec_day / 3600.0,
                "dpsi_arcsec": dpsi_deg * 3600.0,
            })
        for code, legacy_id, hip_id, target_lon in star_models:
            star = stars_by_hip[hip_id]
            lon_0 = star_mean_ecliptic_lon_no_ab(spk, star, jd_tt)
            lon_m = star_mean_ecliptic_lon_no_ab(spk, star, jd_tt - 0.5)
            lon_p = star_mean_ecliptic_lon_no_ab(spk, star, jd_tt + 0.5)
            mean_deg = (lon_0 - target_lon) % 360.0
            true_deg = (mean_deg + dpsi_deg) % 360.0
            rate_deg_day = ((lon_p - lon_m + 180.0) % 360.0) - 180.0
            vectors.append({
                "epoch_label": label,
                "jd_tt": jd_tt,
                "ayanamsa_code": code,
                "legacy_id": legacy_id,
                "category": "star_anchored",
                "anchor_hip": hip_id,
                "target_sidereal_lon_deg": target_lon,
                "ayanamsa_mean_deg": mean_deg,
                "ayanamsa_true_deg": true_deg,
                "rate_deg_day": rate_deg_day,
                "dpsi_arcsec": dpsi_deg * 3600.0,
            })

    return {
        "schema_version": "1.0.0",
        "level": "L1B",
        "source": "Lieske (1977) Table 5 + IAE (1989) + Hipparcos-2 I/311 + XHIP V/137D + pyerfa",
        "ephemeris": "DE440s (de440s.bsp)",
        "time_scale": "TT",
        "tolerance_arcsec": 0.005,
        "count": len(vectors),
        "vectors": vectors,
    }


def brent_root(f, a: float, b: float, tol: float = 1e-11, max_iter: int = 100) -> float:
    """Deterministic Brent's root-finding algorithm on [a, b]."""
    fa = f(a)
    fb = f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    c, fc = a, fa
    d = e = b - a
    for _ in range(max_iter):
        if (fb > 0.0 and fc > 0.0) or (fb < 0.0 and fc < 0.0):
            c, fc = a, fa
            d = e = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol1 = 2.0 * 2.220446049250313e-16 * abs(b) + 0.5 * tol
        xm = 0.5 * (c - b)
        if abs(xm) <= tol1 or fb == 0.0:
            return b
        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p = 2.0 * xm * s
                q = 1.0 - s
            else:
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0.0:
                q = -q
            p = abs(p)
            if 2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q)):
                e = d
                d = p / q
            else:
                d = xm
                e = d
        else:
            d = xm
            e = d
        a, fa = b, fb
        if abs(d) > tol1:
            b += d
        else:
            b += math.copysign(tol1, xm)
        fb = f(b)
    return b


def natural_cubic_spline_eval(xs: list[float], ys: list[float], x_eval: float) -> tuple[float, float]:
    """Natural C2 cubic spline (S''(x_0) = S''(x_{n-1}) = 0) solved via Thomas algorithm."""
    n = len(xs)
    h = [xs[i + 1] - xs[i] for i in range(n - 1)]
    alpha = [0.0] * n
    for i in range(1, n - 1):
        alpha[i] = (3.0 / h[i]) * (ys[i + 1] - ys[i]) - (3.0 / h[i - 1]) * (ys[i] - ys[i - 1])
    l = [1.0] + [0.0] * (n - 1)
    mu = [0.0] * n
    z = [0.0] * n
    for i in range(1, n - 1):
        l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]
    l[n - 1] = 1.0
    c = [0.0] * n
    b = [0.0] * (n - 1)
    d = [0.0] * (n - 1)
    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (ys[j + 1] - ys[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j])
    # Find interval containing x_eval
    idx = 0
    for i in range(n - 1):
        if xs[i] <= x_eval <= xs[i + 1]:
            idx = i
            break
    dx = x_eval - xs[idx]
    val = ys[idx] + b[idx] * dx + c[idx] * dx * dx + d[idx] * dx * dx * dx
    deriv = b[idx] + 2.0 * c[idx] * dx + 3.0 * d[idx] * dx * dx
    return val, deriv


def lunar_osculating_node_and_lilith(spk: SPK, jd_tt: float, nutation: bool = True) -> tuple[float, float, float]:
    """Return (true_node_lon_deg, true_lilith_lon_deg, true_lilith_lat_deg) at jd_tt."""
    jd_tdb = tt_to_tdb_jd(jd_tt)
    p_m, v_m = spk_barycentric_state_km_s(spk, 301, jd_tdb)
    p_e, v_e = spk_barycentric_state_km_s(spk, 399, jd_tdb)
    r_ecl = icrf_to_ecliptic_matrix(jd_tt, nutation=nutation)
    r = r_ecl @ (p_m - p_e)
    v = r_ecl @ (v_m - v_e)
    h_vec = np.cross(r, v)
    node_lon = math.degrees(math.atan2(float(h_vec[0]), -float(h_vec[1]))) % 360.0

    r_norm = float(np.linalg.norm(r))
    e_vec = np.cross(v, h_vec) / MU_EM_KM3_S2 - (r / r_norm)
    a_vec = -e_vec
    lilith_lon, lilith_lat, _ = cart_to_spherical(a_vec)
    return node_lon, lilith_lon, lilith_lat


def interpolated_lilith_c2(spk: SPK, jd_tt: float, nutation: bool = True) -> tuple[float, float, float, float]:
    """Compute Interpolated_Lilith (`natural_apogee_c2_spline_v1`) over [jd_tt - 60, jd_tt + 60]."""
    def radial_rate(t_tt: float) -> float:
        t_tdb = tt_to_tdb_jd(t_tt)
        p_m, v_m = spk_barycentric_state_km_s(spk, 301, t_tdb)
        p_e, v_e = spk_barycentric_state_km_s(spk, 399, t_tdb)
        r = p_m - p_e
        v = v_m - v_e
        return float(np.dot(r, v)) / float(np.linalg.norm(r))

    t_start = jd_tt - 60.0
    step = 2.0
    grid_t = [t_start + i * step for i in range(61)]
    grid_g = [radial_rate(t) for t in grid_t]

    knots_t = []
    knots_lon = []
    knots_lat = []
    for i in range(60):
        if grid_g[i] * grid_g[i + 1] <= 0.0:
            is_apogee = grid_g[i] > grid_g[i + 1]  # g changes + to - at distance maximum
            if is_apogee:
                t_root = brent_root(radial_rate, grid_t[i], grid_t[i + 1])
                t_tdb = tt_to_tdb_jd(t_root)
                p_m, _ = spk_barycentric_state_km_s(spk, 301, t_tdb)
                p_e, _ = spk_barycentric_state_km_s(spk, 399, t_tdb)
                r_ecl = icrf_to_ecliptic_matrix(t_root, nutation=False)
                lon_m, lat_m, _ = cart_to_spherical(r_ecl @ (p_m - p_e))
                knots_lon.append(lon_m)
                knots_lat.append(lat_m)
                knots_t.append(t_root)

    # Phase unwrap longitudes
    unwrapped_lon = [knots_lon[0]]
    for k in range(1, len(knots_lon)):
        diff = ((knots_lon[k] - unwrapped_lon[-1] + 180.0) % 360.0) - 180.0
        unwrapped_lon.append(unwrapped_lon[-1] + diff)

    lon_mean, dlon = natural_cubic_spline_eval(knots_t, unwrapped_lon, jd_tt)
    lat_mean, dlat = natural_cubic_spline_eval(knots_t, knots_lat, jd_tt)
    lon_out = lon_mean % 360.0
    if nutation:
        dpsi_rad, _ = erfa.nut00a(jd_tt, 0.0)
        lon_out = (lon_out + math.degrees(float(dpsi_rad))) % 360.0
    return lon_out, lat_mean, dlon, dlat


def generate_l1b_lunar_points(spk: SPK) -> dict:
    epochs = [
        ("1855-06-01T00:00:00", 2398736.5),
        ("1900-01-01T12:00:00", 2415021.0),
        ("1956-03-21T00:00:00", 2435553.5),
        ("1980-06-15T12:00:00", 2444406.0),
        ("2000-01-01T12:00:00", 2451545.0),
        ("2026-06-15T00:00:00", 2461206.5),
        ("2050-01-01T00:00:00", 2469807.5),
        ("2100-01-01T00:00:00", 2488069.5),
        ("2145-06-01T00:00:00", 2504656.5),
    ]
    vectors = []
    for label, jd_tt in epochs:
        T = (jd_tt - 2451545.0) / 36525.0
        dpsi_rad, _ = erfa.nut00a(jd_tt, 0.0)
        dpsi_deg = math.degrees(float(dpsi_rad))

        # 1. Mean_Node & Mean_South_Node (IERS 2003 faom03)
        om_mean_deg = math.degrees(float(erfa.faom03(T))) % 360.0
        om_true_deg = (om_mean_deg + dpsi_deg) % 360.0
        om_rate_deg_day = (
            -6962890.5431 + 2.0 * 7.4722 * T + 3.0 * 0.007702 * T**2 - 4.0 * 0.00005939 * T**3
        ) / (3600.0 * 36525.0)

        # 2. Mean_Lilith (IERS 2003 F + Om - l + 180)
        f_rad = float(erfa.faf03(T))
        om_rad = float(erfa.faom03(T))
        l_rad = float(erfa.fal03(T))
        mlil_mean_deg = (math.degrees(f_rad + om_rad - l_rad) + 180.0) % 360.0
        mlil_true_deg = (mlil_mean_deg + dpsi_deg) % 360.0
        mlil_rate_deg_day = (
            14648449.0869 - 2.0 * 37.1582 * T - 3.0 * 0.044970 * T**2 + 4.0 * 0.00018948 * T**3
        ) / (3600.0 * 36525.0)

        # 3. True_Node, True_South_Node, True_Lilith with 4-point stencil (h = 1/24 day)
        h = 1.0 / 24.0
        tn_0, tl_0, tlat_0 = lunar_osculating_node_and_lilith(spk, jd_tt, nutation=True)
        tn_m2, tl_m2, _ = lunar_osculating_node_and_lilith(spk, jd_tt - 2 * h, nutation=True)
        tn_m1, tl_m1, _ = lunar_osculating_node_and_lilith(spk, jd_tt - h, nutation=True)
        tn_p1, tl_p1, _ = lunar_osculating_node_and_lilith(spk, jd_tt + h, nutation=True)
        tn_p2, tl_p2, _ = lunar_osculating_node_and_lilith(spk, jd_tt + 2 * h, nutation=True)

        def stencil_angle_rate(a_m2, a_m1, a_p1, a_p2, step):
            d1 = ((a_p1 - a_m1 + 180.0) % 360.0) - 180.0
            d2 = ((a_p2 - a_m2 + 180.0) % 360.0) - 180.0
            return (8.0 * d1 - d2) / (12.0 * step)

        tn_speed = stencil_angle_rate(tn_m2, tn_m1, tn_p1, tn_p2, h)
        tl_speed = stencil_angle_rate(tl_m2, tl_m1, tl_p1, tl_p2, h)

        # 4. Interpolated_Lilith (natural_apogee_c2_spline_v1)
        il_lon, il_lat, il_dlon, il_dlat = interpolated_lilith_c2(spk, jd_tt, nutation=True)

        vectors.append({
            "epoch_label": label,
            "jd_tt": jd_tt,
            "points": {
                "Mean_Node": {"lon_deg": om_true_deg, "lon_mean_deg": om_mean_deg, "lat_deg": 0.0, "lon_speed_deg_day": om_rate_deg_day},
                "Mean_South_Node": {"lon_deg": (om_true_deg + 180.0) % 360.0, "lon_mean_deg": (om_mean_deg + 180.0) % 360.0, "lat_deg": 0.0, "lon_speed_deg_day": om_rate_deg_day},
                "True_Node": {"lon_deg": tn_0, "lat_deg": 0.0, "lon_speed_deg_day": tn_speed},
                "True_South_Node": {"lon_deg": (tn_0 + 180.0) % 360.0, "lat_deg": 0.0, "lon_speed_deg_day": tn_speed},
                "Mean_Lilith": {"lon_deg": mlil_true_deg, "lon_mean_deg": mlil_mean_deg, "lat_deg": 0.0, "lon_speed_deg_day": mlil_rate_deg_day},
                "True_Lilith": {"lon_deg": tl_0, "lat_deg": tlat_0, "lon_speed_deg_day": tl_speed},
                "Interpolated_Lilith": {"lon_deg": il_lon, "lat_deg": il_lat, "lon_speed_deg_day": il_dlon, "lat_speed_deg_day": il_dlat},
            },
        })

    return {
        "schema_version": "1.0.0",
        "level": "L1B",
        "source": "IERS 2003 (faom03/faf03/fal03) + DE440s osculating state + natural_apogee_c2_spline_v1",
        "ephemeris": "DE440s (de440s.bsp)",
        "mu_earth_moon_km3_s2": MU_EM_KM3_S2,
        "time_scale": "TT",
        "tolerance_arcsec": 0.05,
        "count": len(vectors),
        "vectors": vectors,
    }


def uranian_heliocentric_icrf_au(el: dict, jd_tt: float) -> np.ndarray:
    k_deg = 0.9856076686
    a = el["a"]
    e = el["e"]
    n = k_deg / (a ** 1.5)
    m_rad = math.radians((el["M0"] + n * (jd_tt - el["epoch_jd"])) % 360.0)
    E = m_rad
    for _ in range(15):
        dE = (E - e * math.sin(E) - m_rad) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < 1e-14:
            break
    x_orb = a * (math.cos(E) - e)
    y_orb = a * math.sqrt(1.0 - e * e) * math.sin(E)

    w = math.radians(el["w"])
    om = math.radians(el["om"])
    inc = math.radians(el["i"])
    cw, sw = math.cos(w), math.sin(w)
    com, som = math.cos(om), math.sin(om)
    ci, si = math.cos(inc), math.sin(inc)

    x_ecl = (com * cw - som * sw * ci) * x_orb + (-com * sw - som * cw * ci) * y_orb
    y_ecl = (som * cw + com * sw * ci) * x_orb + (-som * sw + com * cw * ci) * y_orb
    z_ecl = (sw * si) * x_orb + (cw * si) * y_orb
    v_ecl_eq = np.array([x_ecl, y_ecl, z_ecl])

    # Source ecliptic of equinox_jd -> ICRF
    r_ecl_eq = icrf_to_ecliptic_matrix(el["equinox_jd"], nutation=False)
    return r_ecl_eq.T @ v_ecl_eq


def generate_l1b_uranians(spk: SPK) -> dict:
    elements = [
        {"name": "Cupido", "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 163.7409, "a": 40.99837, "e": 0.00460, "w": 171.4333, "om": 129.8325, "i": 1.0833},
        {"name": "Hades", "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 27.6496, "a": 50.66744, "e": 0.00245, "w": 148.1796, "om": 161.3339, "i": 1.0500},
        {"name": "Zeus", "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 165.1232, "a": 59.21436, "e": 0.00120, "w": 299.0440, "om": 0.0, "i": 0.0},
        {"name": "Kronos", "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 169.0193, "a": 64.81690, "e": 0.00305, "w": 208.8801, "om": 0.0, "i": 0.0},
        {"name": "Apollon", "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 138.0533, "a": 70.29949, "e": 0.0, "w": 0.0, "om": 0.0, "i": 0.0},
        {"name": "Admetos", "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 351.3350, "a": 73.62765, "e": 0.0, "w": 0.0, "om": 0.0, "i": 0.0},
        {"name": "Vulkanus", "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 55.8983, "a": 77.25568, "e": 0.0, "w": 0.0, "om": 0.0, "i": 0.0},
        {"name": "Poseidon", "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 165.5163, "a": 83.66907, "e": 0.0, "w": 0.0, "om": 0.0, "i": 0.0},
        {"name": "Transpluto", "alias": "Isis", "epoch_jd": 2368547.66, "equinox_jd": 2431456.5, "M0": 0.0, "a": 77.77500, "e": 0.30000, "w": 0.7000, "om": 0.0, "i": 0.0},
    ]
    epochs = [
        ("1900-01-01T00:00:00", 2415020.0),
        ("1945-01-01T00:00:00", 2431456.5),
        ("1980-06-15T12:00:00", 2444406.0),
        ("2000-01-01T12:00:00", 2451545.0),
        ("2026-06-15T00:00:00", 2461206.5),
        ("2050-01-01T00:00:00", 2469807.5),
    ]

    def uranian_apparent_ecl(el: dict, t_tt: float) -> np.ndarray:
        t_tdb = tt_to_tdb_jd(t_tt)
        p_earth_km, v_earth_kms = spk_barycentric_state_km_s(spk, 399, t_tdb)
        p_sun_km, _ = spk_barycentric_state_km_s(spk, 10, t_tdb)
        r_earth_helio_au = (p_earth_km - p_sun_km) / AU_KM
        dt = 0.0
        for _ in range(3):
            r_obj_helio_au = uranian_heliocentric_icrf_au(el, t_tt - dt)
            rho_au = r_obj_helio_au - r_earth_helio_au
            dist_au = float(np.linalg.norm(rho_au))
            dt = dist_au / C_AU_DAY
        u_unit = rho_au / dist_au
        bm1 = math.sqrt(1.0 - min(float(np.dot(v_earth_kms, v_earth_kms)) / (C_KM_S * C_KM_S), 0.999999))
        u_unit = np.asarray(erfa.ab(u_unit, v_earth_kms / C_KM_S, float(np.linalg.norm(r_earth_helio_au)), bm1), dtype=float)
        r_ecl = icrf_to_ecliptic_matrix(t_tt, nutation=True)
        return r_ecl @ (u_unit * dist_au)

    vectors = []
    h = 0.25
    for label, jd_tt in epochs:
        for el in elements:
            v_m2 = uranian_apparent_ecl(el, jd_tt - 2 * h)
            v_m1 = uranian_apparent_ecl(el, jd_tt - h)
            v_0 = uranian_apparent_ecl(el, jd_tt)
            v_p1 = uranian_apparent_ecl(el, jd_tt + h)
            v_p2 = uranian_apparent_ecl(el, jd_tt + 2 * h)
            dv = (v_m2 - 8.0 * v_m1 + 8.0 * v_p1 - v_p2) / (12.0 * h)
            lon_deg, lat_deg, dist_au = cart_to_spherical(v_0)
            dlon, dlat, ddist = cart_state_to_spherical_rates(v_0, dv)
            vectors.append({
                "epoch_label": label,
                "jd_tt": jd_tt,
                "target_name": el["name"],
                "target_kind": "hypothetical",
                "apparent_ecliptic_true_of_date": {
                    "lon_deg": lon_deg,
                    "lat_deg": lat_deg,
                    "dist_au": dist_au,
                    "lon_speed_deg_day": dlon,
                    "lat_speed_deg_day": dlat,
                    "dist_speed_au_day": ddist,
                },
            })

    return {
        "schema_version": "1.0.0",
        "level": "L1B",
        "source": "James Neely (1980) + Sevin/Strubell (1952) Keplerian elements + DE440s Earth state",
        "ephemeris": "DE440s (de440s.bsp)",
        "time_scale": "TT",
        "tolerance_arcsec": 0.05,
        "count": len(vectors),
        "vectors": vectors,
    }


def generate_l1b_stars(spk: SPK, stars_list: list[dict]) -> dict:
    epochs = [
        ("2000-01-01T12:00:00", 2451545.0),
        ("2026-06-15T00:00:00", 2461206.5),
    ]
    vectors = []
    for label, jd_tt in epochs:
        jd_tdb = tt_to_tdb_jd(jd_tt)
        p_earth_km, v_earth_kms = spk_barycentric_state_km_s(spk, 399, jd_tdb)
        p_sun_km, _ = spk_barycentric_state_km_s(spk, 10, jd_tdb)
        p_eb_au = p_earth_km / AU_KM
        e_sun_vec = (p_earth_km - p_sun_km) / AU_KM
        em = float(np.linalg.norm(e_sun_vec))
        e_sun_unit = e_sun_vec / em
        bm1 = math.sqrt(1.0 - min(float(np.dot(v_earth_kms, v_earth_kms)) / (C_KM_S * C_KM_S), 0.999999))
        v_ab = v_earth_kms / C_KM_S
        r_ecl_true = icrf_to_ecliptic_matrix(jd_tt, nutation=True)
        r_equ_true = icrf_to_equatorial_matrix(jd_tt, nutation=True)

        pmt = (jd_tdb - 2448349.0625) / 365.25
        for s in stars_list:
            ra1 = math.radians(s["ra_1991_deg"])
            de1 = math.radians(s["de_1991_deg"])
            pmr = (s["pm_ra_cosdec_mas_yr"] * 1e-3 * DAS2R) / math.cos(de1)
            pmd = s["pm_dec_mas_yr"] * 1e-3 * DAS2R
            px = max(s["plx_mas"] * 1e-3, 0.0)
            u_ppx = np.asarray(erfa.pmpx(ra1, de1, pmr, pmd, px, s["rv_kms"], pmt, p_eb_au), dtype=float)
            u_ld = np.asarray(erfa.ldsun(u_ppx, e_sun_unit, em), dtype=float)
            u_ab = np.asarray(erfa.ab(u_ld, v_ab, em, bm1), dtype=float)

            lon_deg, lat_deg, _ = cart_to_spherical(r_ecl_true @ u_ab)
            ra_deg, dec_deg, _ = cart_to_spherical(r_equ_true @ u_ab)
            vectors.append({
                "epoch_label": label,
                "jd_tt": jd_tt,
                "name": s["name"],
                "bayer": s["bayer"],
                "hip": s["hip"],
                "hr": s["hr"],
                "hd": s["hd"],
                "comp": s["comp"],
                "is_astra_50": s["is_astra_50"],
                "vmag": s["vmag"],
                "flags": s["flags"],
                "icrs_j1991_25": {
                    "ra_deg": s["ra_1991_deg"],
                    "dec_deg": s["de_1991_deg"],
                    "pm_ra_cosdec_mas_yr": s["pm_ra_cosdec_mas_yr"],
                    "pm_dec_mas_yr": s["pm_dec_mas_yr"],
                    "plx_mas": s["plx_mas"],
                    "rv_kms": s["rv_kms"],
                },
                "icrs_j2000_0": {
                    "ra_deg": s["ra_2000_deg"],
                    "dec_deg": s["de_2000_deg"],
                },
                "apparent_of_date": {
                    "lon_deg": lon_deg,
                    "lat_deg": lat_deg,
                    "ra_deg": ra_deg,
                    "dec_deg": dec_deg,
                },
            })

    return {
        "schema_version": "1.0.0",
        "level": "L1B",
        "source": "VizieR I/311 (Hipparcos-2, J1991.25) + XHIP V/137D RV + pyerfa (pmpx -> ldsun -> ab -> npb06)",
        "ephemeris": "DE440s (de440s.bsp)",
        "time_scale": "TT",
        "tolerance_arcsec": 0.01,
        "count": len(vectors),
        "vectors": vectors,
    }


def main() -> None:
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    spk = SPK.open(str(DE440S_PATH))

    # Load verified 52-star catalog from testdata/golden/stars_catalog_52.json
    stars_cat_path = GOLDEN_DIR / "stars_catalog_52.json"
    with open(stars_cat_path) as f:
        stars_list = json.load(f)
    stars_by_hip = {s["hip"]: s for s in stars_list}

    p_data = generate_l1b_planets(spk)
    with open(GOLDEN_DIR / "l1b_planets_de440s.json", "w") as f:
        json.dump(p_data, f, indent=2)

    a_data = generate_l1b_ayanamsas(spk, stars_by_hip)
    with open(GOLDEN_DIR / "l1b_ayanamsas_12.json", "w") as f:
        json.dump(a_data, f, indent=2)

    l_data = generate_l1b_lunar_points(spk)
    with open(GOLDEN_DIR / "l1b_lunar_points_7.json", "w") as f:
        json.dump(l_data, f, indent=2)

    u_data = generate_l1b_uranians(spk)
    with open(GOLDEN_DIR / "l1b_uranians_9.json", "w") as f:
        json.dump(u_data, f, indent=2)

    s_data = generate_l1b_stars(spk, stars_list)
    with open(GOLDEN_DIR / "l1b_stars_50.json", "w") as f:
        json.dump(s_data, f, indent=2)

    print(
        f"Generated L1B golden fixtures:\n"
        f"  - l1b_planets_de440s.json : {p_data['count']} vectors\n"
        f"  - l1b_ayanamsas_12.json   : {a_data['count']} vectors\n"
        f"  - l1b_lunar_points_7.json : {l_data['count']} epochs x 7 points = {l_data['count']*7} vectors\n"
        f"  - l1b_uranians_9.json     : {u_data['count']} vectors\n"
        f"  - l1b_stars_50.json       : {s_data['count']} vectors"
    )


if __name__ == "__main__":
    main()
