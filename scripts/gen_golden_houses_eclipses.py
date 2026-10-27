#!/usr/bin/env python3
"""
Clean-Room Golden Vector Generator for:
  1. testdata/golden/l1b_houses_23.json
     - All 23 house systems, 8 chart angles, and 36 Gauquelin sectors
     - Mid-latitude N & S hemispheres, +/-180 deg antimeridian, equator, and polar latitudes
  2. testdata/golden/l1b_eclipses.json
     - Canonical Solar (Total, Annular, Hybrid, Partial) and Lunar (Total, Partial, Penumbral)
       eclipses computed from JPL DE440s using first-principles Besselian / Fundamental Plane
       geometry (MODEL_PASSPORTS.md §6) and cross-referenced against the NASA Espenak & Meeus
       Five Millennium Canon of Solar and Lunar Eclipses (TD/TT).
"""

import json
import math
from pathlib import Path
import erfa
import numpy as np
from jplephem.spk import SPK

ROOT = Path(__file__).resolve().parent.parent
GOLDEN_DIR = ROOT / "testdata" / "golden"
BSP_PATH = ROOT / "data" / "de440s.bsp"

DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi
AU_KM = 149597870.700
C_KM_DAY = 299792.458 * 86400.0
R_EARTH_KM = 6378.1366
F_EARTH = 1.0 / 298.257223563
R_SUN_KM = 696000.0
K_MOON_VALLEY = 0.2722810
K_MOON_MEAN = 0.2725076


def norm360(x: float) -> float:
    return x % 360.0


def sind(x: float) -> float:
    return math.sin(x * DEG2RAD)


def cosd(x: float) -> float:
    return math.cos(x * DEG2RAD)


def tand(x: float) -> float:
    return math.tan(x * DEG2RAD)


def asind(x: float) -> float:
    return math.asin(max(-1.0, min(1.0, x))) * RAD2DEG


def acosd(x: float) -> float:
    return math.acos(max(-1.0, min(1.0, x))) * RAD2DEG


def atand(x: float) -> float:
    return math.atan(x) * RAD2DEG


def atan2d(y: float, x: float) -> float:
    return norm360(math.atan2(y, x) * RAD2DEG)


# ---------------------------------------------------------------------------
# Section 5: 8 Chart Angles, 23 House Systems, and 36 Gauquelin Sectors
# ---------------------------------------------------------------------------

def calc_mc(ramc: float, eps: float) -> float:
    return atan2d(sind(ramc), cosd(ramc) * cosd(eps))


def calc_asc(ramc: float, lat: float, eps: float) -> float:
    y = cosd(ramc)
    x = -sind(ramc) * cosd(eps) - tand(lat) * sind(eps)
    return atan2d(y, x)


def compute_8_angles(ramc: float, lat: float, eps: float) -> dict:
    asc = calc_asc(ramc, lat, eps)
    mc = calc_mc(ramc, eps)
    armc = norm360(ramc)
    if abs(lat) > 1e-12:
        cot_lat = 1.0 / tand(lat)
        vertex = atan2d(-cosd(ramc), sind(ramc) * cosd(eps) - cot_lat * sind(eps))
        coasc_munk = atan2d(cosd(ramc), -sind(ramc) * cosd(eps) - cot_lat * sind(eps))
    else:
        vertex = 180.0
        coasc_munk = 180.0
    if norm360(asc - mc) >= 180.0:
        asc = norm360(asc + 180.0)
    if 1e-12 < norm360(vertex - mc) < 180.0:
        vertex = norm360(vertex + 180.0)
    eqasc = calc_asc(ramc, 0.0, eps)
    polar_asc = calc_asc(norm360(ramc + 180.0), lat, eps)
    coasc_koch = norm360(polar_asc + 180.0)
    return {
        "ASC": asc,
        "MC": mc,
        "ARMC": armc,
        "Vertex": vertex,
        "EqASC": eqasc,
        "CoASC_Koch": coasc_koch,
        "CoASC_Munkasey": coasc_munk,
        "PolarASC": polar_asc,
        "DSC": norm360(asc + 180.0),
        "IC": norm360(mc + 180.0),
    }


def placidus_offset(ramc: float, lat: float, eps: float, f: float) -> float:
    """Solve g(d) = d - f * acos(-sin(ramc + d)*tan(lat)*tan(eps)) = 0 via bisection."""
    tt = tand(lat) * tand(eps)
    r = ramc * DEG2RAD
    lo, hi = (0.0, math.pi) if f > 0 else (-math.pi, 0.0)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        arg = -math.sin(r + mid) * tt
        if abs(arg) >= 1.0:
            raise ValueError("Circumpolar ecliptic point in Placidus")
        g_mid = mid - f * math.acos(arg)
        arg_lo = -math.sin(r + lo) * tt
        g_lo = lo - f * math.acos(max(-1.0, min(1.0, arg_lo)))
        if g_lo * g_mid <= 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi) * RAD2DEG


def compute_gauquelin_36(ramc: float, lat: float, eps: float) -> list | None:
    """Compute the 36 Gauquelin sector ecliptic cusps (numbered 1..36 clockwise)."""
    if abs(lat) >= 90.0 - eps:
        return None
    asc = calc_asc(ramc, lat, eps)
    mc = calc_mc(ramc, eps)
    dsc = norm360(asc + 180.0)
    ic = norm360(mc + 180.0)
    sectors = [0.0] * 36
    sectors[0] = asc   # Sector 1 cusp = ASC
    sectors[9] = mc    # Sector 10 cusp = MC
    sectors[18] = dsc  # Sector 19 cusp = DSC
    sectors[27] = ic   # Sector 28 cusp = IC
    # Sectors 2..9: rising diurnal between ASC and MC (fractions 8/9 .. 1/9 east of MC)
    for idx in range(1, 9):
        f = (9 - idx) / 9.0
        d = placidus_offset(ramc, lat, eps, f)
        c = calc_mc(ramc + d, eps)
        sectors[idx] = c
        sectors[idx + 18] = norm360(c + 180.0)
    # Sectors 29..36: rising nocturnal between IC and ASC (opposed to setting diurnal Sectors 11..18)
    for idx in range(1, 9):
        f = -idx / 9.0
        d = placidus_offset(ramc, lat, eps, f)
        c_setting = calc_mc(ramc + d, eps)
        sectors[9 + idx] = c_setting
        sectors[27 + idx] = norm360(c_setting + 180.0)
    return sectors


def _slerp(u: np.ndarray, v: np.ndarray, t: float) -> np.ndarray:
    dot = max(-1.0, min(1.0, float(np.dot(u, v))))
    omega = math.acos(dot)
    if abs(omega) < 1e-12:
        return u
    return (math.sin((1.0 - t) * omega) * u + math.sin(t * omega) * v) / math.sin(omega)


def _great_circle_ecl_lon(normal: np.ndarray, target_dir: np.ndarray, eps_rad: float) -> float:
    nx, ny, nz = normal
    lam = math.atan2(-nx, ny * math.cos(eps_rad) + nz * math.sin(eps_rad))
    v = np.array([
        math.cos(lam),
        math.sin(lam) * math.cos(eps_rad),
        math.sin(lam) * math.sin(eps_rad),
    ])
    if float(np.dot(v, target_dir)) < 0.0:
        lam += math.pi
    return norm360(lam * RAD2DEG)


def compute_all_23_houses(ramc: float, lat: float, eps: float, sun_decl: float) -> dict:
    """
    Compute all 23 house systems according to MODEL_PASSPORTS.md §5.
    Returns a dict mapping single-char code -> {name, polar_class, status, cusps}.
    """
    if abs(lat) >= 90.0:
        # At exact geographic pole, only latitude-independent systems are defined
        pass

    angles = compute_8_angles(ramc, lat, eps) if abs(lat) < 90.0 else None
    asc = angles["ASC"] if angles else 0.0
    mc = angles["MC"] if angles else calc_mc(ramc, eps)
    dsc = norm360(asc + 180.0)
    ic = norm360(mc + 180.0)
    q4 = norm360(asc - mc)
    q1 = norm360(ic - asc)

    results = {}

    def quad_from_4(c11, c12, c2, c3):
        return [
            asc,
            norm360(c2),
            norm360(c3),
            ic,
            norm360(c11 + 180.0),
            norm360(c12 + 180.0),
            dsc,
            norm360(c2 + 180.0),
            norm360(c3 + 180.0),
            mc,
            norm360(c11),
            norm360(c12),
        ]

    is_polar = abs(lat) >= (90.0 - eps)
    is_pole = abs(lat) >= 90.0

    # 1. P - Placidus (hard_fail)
    if is_pole or is_polar:
        results["P"] = {
            "name": "Placidus",
            "polar_class": "hard_fail",
            "status": "Err(SingularPole)" if is_pole else "Err(HighLatitudeCircumpolar)",
            "cusps": None,
        }
    else:
        c11 = calc_mc(ramc + placidus_offset(ramc, lat, eps, 1.0 / 3.0), eps)
        c12 = calc_mc(ramc + placidus_offset(ramc, lat, eps, 2.0 / 3.0), eps)
        c8 = calc_mc(ramc + placidus_offset(ramc, lat, eps, -2.0 / 3.0), eps)
        c9 = calc_mc(ramc + placidus_offset(ramc, lat, eps, -1.0 / 3.0), eps)
        c2 = norm360(c8 + 180.0)
        c3 = norm360(c9 + 180.0)
        results["P"] = {
            "name": "Placidus",
            "polar_class": "hard_fail",
            "status": "Ok",
            "cusps": quad_from_4(c11, c12, c2, c3),
        }

    # 2. K - Koch (hard_fail)
    decl_mc = asind(sind(mc) * sind(eps))
    koch_arg = tand(decl_mc) * tand(lat) if not is_pole else 999.0
    if is_pole or is_polar or abs(koch_arg) >= 1.0:
        results["K"] = {
            "name": "Koch",
            "polar_class": "hard_fail",
            "status": "Err(SingularPole)" if is_pole else "Err(HighLatitudeCircumpolar)",
            "cusps": None,
        }
    else:
        ad_mc = asind(koch_arg)
        sda_mc = 90.0 + ad_mc
        c11 = calc_asc(ramc - 2.0 * sda_mc / 3.0, lat, eps)
        c12 = calc_asc(ramc - 1.0 * sda_mc / 3.0, lat, eps)
        c2 = calc_asc(ramc + 1.0 * sda_mc / 3.0, lat, eps)
        c3 = calc_asc(ramc + 2.0 * sda_mc / 3.0, lat, eps)
        results["K"] = {
            "name": "Koch",
            "polar_class": "hard_fail",
            "status": "Ok",
            "cusps": quad_from_4(c11, c12, c2, c3),
        }

    # 3. O - Porphyry (polar_safe)
    if is_pole:
        results["O"] = {"name": "Porphyry", "polar_class": "polar_safe", "status": "Err(SingularPole)", "cusps": None}
    else:
        porph = quad_from_4(mc + q4 / 3.0, mc + 2.0 * q4 / 3.0, asc + q1 / 3.0, asc + 2.0 * q1 / 3.0)
        results["O"] = {"name": "Porphyry", "polar_class": "polar_safe", "status": "Ok", "cusps": porph}

    # 4. R - Regiomontanus (soft_fail)
    if is_pole:
        results["R"] = {"name": "Regiomontanus", "polar_class": "soft_fail", "status": "Err(SingularPole)", "cusps": None}
    else:
        def regio_cusp(h):
            rh = ramc + h
            return atan2d(sind(rh), cosd(rh) * cosd(eps) - sind(h) * tand(lat) * sind(eps))
        results["R"] = {
            "name": "Regiomontanus",
            "polar_class": "soft_fail",
            "status": "OkWithWarning(HighLatitudeQuadrantSkew)" if is_polar else "Ok",
            "cusps": quad_from_4(regio_cusp(30.0), regio_cusp(60.0), regio_cusp(120.0), regio_cusp(150.0)),
        }

    # Unit vectors for 3D geometric systems (C, U, H)
    if not is_pole:
        r_rad = ramc * DEG2RAD
        p_rad = lat * DEG2RAD
        e_rad = eps * DEG2RAD
        east = np.array([-math.sin(r_rad), math.cos(r_rad), 0.0])
        zenith = np.array([math.cos(p_rad) * math.cos(r_rad), math.cos(p_rad) * math.sin(r_rad), math.sin(p_rad)])
        north = np.cross(zenith, east)

    # 5. C - Campanus (soft_fail)
    if is_pole:
        results["C"] = {"name": "Campanus", "polar_class": "soft_fail", "status": "Err(SingularPole)", "cusps": None}
    else:
        def camp_cusp(h):
            elev = (90.0 - h) * DEG2RAD
            pt = math.cos(elev) * east + math.sin(elev) * zenith
            normal = np.cross(north, pt)
            return _great_circle_ecl_lon(normal, pt, e_rad)
        results["C"] = {
            "name": "Campanus",
            "polar_class": "soft_fail",
            "status": "OkWithWarning(HighLatitudeQuadrantSkew)" if is_polar else "Ok",
            "cusps": quad_from_4(camp_cusp(30.0), camp_cusp(60.0), camp_cusp(120.0), camp_cusp(150.0)),
        }

    # 6. A - Equal from ASC (polar_safe)
    results["A"] = {
        "name": "Equal (ASC)",
        "polar_class": "polar_safe",
        "status": "Err(SingularPole)" if is_pole else "Ok",
        "cusps": None if is_pole else [norm360(asc + 30.0 * k) for k in range(12)],
    }

    # 7. W - Whole Sign (polar_safe)
    ws_base = math.floor(asc / 30.0) * 30.0
    results["W"] = {
        "name": "Whole Sign",
        "polar_class": "polar_safe",
        "status": "Err(SingularPole)" if is_pole else "Ok",
        "cusps": None if is_pole else [norm360(ws_base + 30.0 * k) for k in range(12)],
    }

    # 8. B - Alcabitius (soft_fail)
    if is_pole:
        results["B"] = {"name": "Alcabitius", "polar_class": "soft_fail", "status": "Err(SingularPole)", "cusps": None}
    else:
        decl_asc = asind(sind(eps) * sind(asc))
        s_alc = tand(lat) * tand(decl_asc)
        clamped = max(-1.0, min(1.0, s_alc))
        sday = 90.0 + asind(clamped)
        snight = 180.0 - sday
        results["B"] = {
            "name": "Alcabitius",
            "polar_class": "soft_fail",
            "status": "OkWithWarning(SoftClamp)" if abs(s_alc) > 1.0 or is_polar else "Ok",
            "cusps": quad_from_4(
                calc_mc(ramc + sday / 3.0, eps),
                calc_mc(ramc + 2.0 * sday / 3.0, eps),
                calc_mc(ramc + sday + snight / 3.0, eps),
                calc_mc(ramc + sday + 2.0 * snight / 3.0, eps),
            ),
        }

    # 9. M - Morinus (polar_safe, latitude-independent)
    results["M"] = {
        "name": "Morinus",
        "polar_class": "polar_safe",
        "status": "Ok",
        "cusps": [atan2d(sind(ramc + 60.0 + 30.0 * (k + 1)) * cosd(eps), cosd(ramc + 60.0 + 30.0 * (k + 1))) for k in range(12)],
    }

    # 10. T - Polich-Page Topocentric (soft_fail)
    if is_pole:
        results["T"] = {"name": "Topocentric (Polich-Page)", "polar_class": "soft_fail", "status": "Err(SingularPole)", "cusps": None}
    else:
        t1 = atand(tand(lat) / 3.0)
        t2 = atand(2.0 * tand(lat) / 3.0)
        results["T"] = {
            "name": "Topocentric (Polich-Page)",
            "polar_class": "soft_fail",
            "status": "OkWithWarning(HighLatitudeQuadrantSkew)" if is_polar else "Ok",
            "cusps": quad_from_4(
                calc_asc(ramc - 60.0, t1, eps),
                calc_asc(ramc - 30.0, t2, eps),
                calc_asc(ramc + 30.0, t2, eps),
                calc_asc(ramc + 60.0, t1, eps),
            ),
        }

    # 11. U - Krusinski-Pisa-Goelzer (polar_safe)
    if is_pole:
        results["U"] = {"name": "Krusinski-Pisa-Goelzer", "polar_class": "polar_safe", "status": "Err(SingularPole)", "cusps": None}
    else:
        asc_rad = asc * DEG2RAD
        v_asc = np.array([math.cos(asc_rad), math.sin(asc_rad) * math.cos(e_rad), math.sin(asc_rad) * math.sin(e_rad)])
        def krus_pt(target, frac):
            p = _slerp(v_asc, target, frac)
            ra = math.atan2(p[1], p[0]) * RAD2DEG
            return calc_mc(ra, eps)
        results["U"] = {
            "name": "Krusinski-Pisa-Goelzer",
            "polar_class": "polar_safe",
            "status": "Ok",
            "cusps": quad_from_4(
                krus_pt(zenith, 2.0 / 3.0),
                krus_pt(zenith, 1.0 / 3.0),
                krus_pt(-zenith, 1.0 / 3.0),
                krus_pt(-zenith, 2.0 / 3.0),
            ),
        }

    # 12. X - Meridian / Axial Rotation (polar_safe, latitude-independent)
    results["X"] = {
        "name": "Meridian (Axial Rotation)",
        "polar_class": "polar_safe",
        "status": "Ok",
        "cusps": [calc_mc(ramc + 60.0 + 30.0 * (k + 1), eps) for k in range(12)],
    }

    # 13. V - Vehlow Equal (polar_safe)
    results["V"] = {
        "name": "Vehlow Equal",
        "polar_class": "polar_safe",
        "status": "Err(SingularPole)" if is_pole else "Ok",
        "cusps": None if is_pole else [norm360(asc - 15.0 + 30.0 * k) for k in range(12)],
    }

    # 14. L - Pullen SD / Sinusoidal Delta (soft_fail)
    if is_pole:
        results["L"] = {"name": "Pullen SD (Sinusoidal Delta)", "polar_class": "soft_fail", "status": "Err(SingularPole)", "cusps": None}
    else:
        w4 = (q4 + 30.0) / 4.0 if q4 >= 30.0 else q4 / 2.0
        w1 = (q1 + 30.0) / 4.0 if q1 >= 30.0 else q1 / 2.0
        results["L"] = {
            "name": "Pullen SD (Sinusoidal Delta)",
            "polar_class": "soft_fail",
            "status": "OkWithWarning(PinchedHouse)" if (q4 < 30.0 or q1 < 30.0) else "Ok",
            "cusps": quad_from_4(mc + w4, asc - w4, asc + w1, ic - w1),
        }

    # 15. Q - Pullen SR / Sinusoidal Ratio (polar_safe)
    if is_pole:
        results["Q"] = {"name": "Pullen SR (Sinusoidal Ratio)", "polar_class": "polar_safe", "status": "Err(SingularPole)", "cusps": None}
    else:
        q_small = min(q4, q1)
        ratio_q = (180.0 - q_small) / max(q_small, 1e-12)
        r = max(1.0, (2.0 * ratio_q) ** (1.0 / 3.0) + 1.0)
        for _ in range(30):
            f_val = (r ** 3) * (r + 2.0) - ratio_q * (2.0 * r + 1.0)
            df_val = 4.0 * (r ** 3) + 6.0 * (r ** 2) - 2.0 * ratio_q
            r -= f_val / df_val
        x_base = q_small / (2.0 * r + 1.0)
        w_small = r * x_base
        w_large = (r ** 3) * x_base
        if q4 <= q1:
            w4_sr, w1_sr = w_small, w_large
        else:
            w4_sr, w1_sr = w_large, w_small
        results["Q"] = {
            "name": "Pullen SR (Sinusoidal Ratio)",
            "polar_class": "polar_safe",
            "status": "Ok",
            "cusps": quad_from_4(mc + w4_sr, asc - w4_sr, asc + w1_sr, ic - w1_sr),
        }

    # 16 & 17. I & i - Sunshine (Treindl & Makransky)
    if is_pole:
        results["I"] = {"name": "Sunshine (Treindl)", "polar_class": "soft_fail", "status": "Err(SingularPole)", "cusps": None}
        results["i"] = {"name": "Sunshine (Makransky)", "polar_class": "hard_fail", "status": "Err(SingularPole)", "cusps": None}
    else:
        sin_kv = tand(sun_decl) * tand(lat)
        kv_clamped = math.asin(max(-1.0, min(1.0, sin_kv)))
        cusps_I = [0.0] * 12
        cusps_I[0], cusps_I[3], cusps_I[6], cusps_I[9] = asc, ic, dsc, mc
        for h in [2, 3, 5, 6, 8, 9, 11, 12]:
            if h < 8:
                a = kv_clamped + r_rad + 0.5 * math.pi + ((h - 1.0) / 3.0) * (0.5 * math.pi - kv_clamped)
            else:
                a = kv_clamped + r_rad + 0.5 * math.pi + ((h - 13.0) / 3.0) * (0.5 * math.pi + kv_clamped)
            y_val = sin_kv * math.sin(r_rad) + math.sin(a)
            x_val = math.cos(e_rad) * (sin_kv * math.cos(r_rad) + math.cos(a)) + math.sin(e_rad) * math.tan(p_rad) * math.sin(r_rad - a)
            cusps_I[h - 1] = norm360(math.atan2(y_val, x_val) * RAD2DEG)
        results["I"] = {
            "name": "Sunshine (Treindl)",
            "polar_class": "soft_fail",
            "status": "OkWithWarning(SoftClamp)" if abs(sin_kv) > 1.0 else "Ok",
            "cusps": cusps_I,
        }
        if abs(sin_kv) > 1.0:
            results["i"] = {
                "name": "Sunshine (Makransky)",
                "polar_class": "hard_fail",
                "status": "Err(HighLatitudeCircumpolar)",
                "cusps": None,
            }
        else:
            results["i"] = {
                "name": "Sunshine (Makransky)",
                "polar_class": "hard_fail",
                "status": "Ok",
                "cusps": list(cusps_I),
            }

    # 18. F - Carter Poli-Equatorial (polar_safe)
    if is_pole:
        results["F"] = {"name": "Carter Poli-Equatorial", "polar_class": "polar_safe", "status": "Err(SingularPole)", "cusps": None}
    else:
        ra_asc = atan2d(sind(asc) * cosd(eps), cosd(asc))
        results["F"] = {
            "name": "Carter Poli-Equatorial",
            "polar_class": "polar_safe",
            "status": "Ok",
            "cusps": [calc_mc(ra_asc + 30.0 * k, eps) for k in range(12)],
        }

    # 19. N - Whole Sign Aries (polar_safe, latitude-independent)
    results["N"] = {
        "name": "Whole Sign (0 deg Aries)",
        "polar_class": "polar_safe",
        "status": "Ok",
        "cusps": [30.0 * k for k in range(12)],
    }

    # 20. S - Sripati (polar_safe)
    if is_pole:
        results["S"] = {"name": "Sripati", "polar_class": "polar_safe", "status": "Err(SingularPole)", "cusps": None}
    else:
        porph = results["O"]["cusps"]
        sri = []
        for k in range(12):
            prev_c = porph[(k - 1) % 12]
            curr_c = porph[k]
            sri.append(norm360(prev_c + norm360(curr_c - prev_c) / 2.0))
        results["S"] = {"name": "Sripati", "polar_class": "polar_safe", "status": "Ok", "cusps": sri}

    # 21. D - Equal from MC (polar_safe)
    results["D"] = {
        "name": "Equal (MC)",
        "polar_class": "polar_safe",
        "status": "Ok",
        "cusps": [norm360(mc + 90.0 + 30.0 * k) for k in range(12)],
    }

    # 22. Y - APC / Apparent Pole (soft_fail)
    if is_pole:
        results["Y"] = {"name": "APC (Apparent Pole)", "polar_class": "soft_fail", "status": "Err(SingularPole)", "cusps": None}
    else:
        kv_apc = math.atan(
            math.tan(p_rad) * math.tan(e_rad) * math.cos(r_rad)
            / (1.0 + math.tan(p_rad) * math.tan(e_rad) * math.sin(r_rad))
        )
        sin_kv_apc = math.sin(kv_apc)
        cusps_Y = [0.0] * 12
        cusps_Y[0], cusps_Y[3], cusps_Y[6], cusps_Y[9] = asc, ic, dsc, mc
        for h in [2, 3, 5, 6, 8, 9, 11, 12]:
            if h < 8:
                a = kv_apc + r_rad + 0.5 * math.pi + ((h - 1.0) / 3.0) * (0.5 * math.pi - kv_apc)
            else:
                a = kv_apc + r_rad + 0.5 * math.pi + ((h - 13.0) / 3.0) * (0.5 * math.pi + kv_apc)
            y_val = sin_kv_apc * math.sin(r_rad) + math.sin(a)
            x_val = math.cos(e_rad) * (sin_kv_apc * math.cos(r_rad) + math.cos(a)) + math.sin(e_rad) * math.tan(p_rad) * math.sin(r_rad - a)
            cusps_Y[h - 1] = norm360(math.atan2(y_val, x_val) * RAD2DEG)
        results["Y"] = {
            "name": "APC (Apparent Pole)",
            "polar_class": "soft_fail",
            "status": "OkWithWarning(HighLatitudeQuadrantSkew)" if is_polar else "Ok",
            "cusps": cusps_Y,
        }

    # 23. H - Azimuthal / Horizontal (polar_safe)
    if is_pole:
        results["H"] = {"name": "Azimuthal (Horizontal)", "polar_class": "polar_safe", "status": "Err(SingularPole)", "cusps": None}
    else:
        ic_dir = north if lat >= 0.0 else -north
        def horiz_cusp(h):
            theta = (h - 90.0) * DEG2RAD
            pt = math.cos(theta) * east + math.sin(theta) * ic_dir
            normal = np.cross(zenith, pt)
            return _great_circle_ecl_lon(normal, pt, e_rad)
        c11_h = horiz_cusp(30.0)
        c12_h = horiz_cusp(60.0)
        c1_h = horiz_cusp(90.0)
        c2_h = horiz_cusp(120.0)
        c3_h = horiz_cusp(150.0)
        results["H"] = {
            "name": "Azimuthal (Horizontal)",
            "polar_class": "polar_safe",
            "status": "Ok",
            "cusps": [
                c1_h,
                c2_h,
                c3_h,
                ic,
                norm360(c11_h + 180.0),
                norm360(c12_h + 180.0),
                norm360(c1_h + 180.0),
                norm360(c2_h + 180.0),
                norm360(c3_h + 180.0),
                mc,
                c11_h,
                c12_h,
            ],
        }

    return results


# ---------------------------------------------------------------------------
# Section 6: Solar & Lunar Eclipses via JPL DE440s + Besselian Geometry
# ---------------------------------------------------------------------------

class De440sEclipseEvaluator:
    def __init__(self, bsp_path: Path):
        self.kernel = SPK.open(str(bsp_path))

    def _earth_icrf(self, jd_tdb: float):
        p_emb, v_emb = self.kernel[0, 3].compute_and_differentiate(jd_tdb)
        p_e, v_e = self.kernel[3, 399].compute_and_differentiate(jd_tdb)
        return p_emb + p_e, v_emb + v_e

    def _moon_icrf(self, jd_tdb: float):
        p_emb, v_emb = self.kernel[0, 3].compute_and_differentiate(jd_tdb)
        p_m, v_m = self.kernel[3, 301].compute_and_differentiate(jd_tdb)
        return p_emb + p_m, v_emb + v_m

    def _sun_icrf(self, jd_tdb: float):
        return self.kernel[0, 10].compute_and_differentiate(jd_tdb)

    def apparent_sun_moon_eq_of_date(self, jd_tt: float):
        """
        Return apparent geocentric vectors (r_sun_km, r_moon_km) rotated to the
        mean/true equatorial frame of date using Lieske/IAU precession.
        """
        p_e, v_e = self._earth_icrf(jd_tt)
        # Light-time to Sun
        p_s0, _ = self._sun_icrf(jd_tt)
        tau_s = float(np.linalg.norm(p_s0 - p_e)) / C_KM_DAY
        p_s, _ = self._sun_icrf(jd_tt - tau_s)
        r_sun = p_s - p_e
        # Annual stellar aberration for Sun
        d_s = float(np.linalg.norm(r_sun))
        r_sun_app = r_sun + d_s * (v_e / C_KM_DAY)
        r_sun_app = r_sun_app * (d_s / float(np.linalg.norm(r_sun_app)))

        # Light-time to Moon (geocentric)
        p_m0, _ = self._moon_icrf(jd_tt)
        tau_m = float(np.linalg.norm(p_m0 - p_e)) / C_KM_DAY
        p_m, _ = self._moon_icrf(jd_tt - tau_m)
        r_moon_app = p_m - p_e

        # Rotate ICRF -> Equatorial of Date via IAU 1976/2000 precession angles (zeta_A, z_A, theta_A)
        t = (jd_tt - 2451545.0) / 36525.0
        zeta = ((2306.2181 + (0.30188 - 0.017998 * t) * t) * t) / 3600.0 * DEG2RAD
        z = ((2306.2181 + (1.09468 + 0.018203 * t) * t) * t) / 3600.0 * DEG2RAD
        theta = ((2004.3109 - (0.42665 + 0.041833 * t) * t) * t) / 3600.0 * DEG2RAD

        cz, sz = math.cos(zeta), math.sin(zeta)
        cz2, sz2 = math.cos(z), math.sin(z)
        ct, st = math.cos(theta), math.sin(theta)
        p_mat = np.array([
            [cz2 * ct * cz - sz2 * sz, -cz2 * ct * sz - sz2 * cz, -cz2 * st],
            [sz2 * ct * cz + cz2 * sz, -sz2 * ct * sz + cz2 * cz, -sz2 * st],
            [st * cz, -st * sz, ct],
        ])
        return p_mat @ r_sun_app, p_mat @ r_moon_app

    def brent_minimize(self, f, a: float, b: float, tol: float = 1e-8):
        gr = 0.3819660112501051
        x = w = v = a + gr * (b - a)
        fx = fw = fv = f(x)
        d = e = 0.0
        for _ in range(80):
            xm = 0.5 * (a + b)
            tol1 = tol * abs(x) + 1e-12
            tol2 = 2.0 * tol1
            if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
                return x, fx
            if abs(e) > tol1:
                r = (x - w) * (fx - fv)
                q = (x - v) * (fx - fw)
                p = (x - v) * q - (x - w) * r
                q = 2.0 * (q - r)
                if q > 0.0:
                    p = -p
                q = abs(q)
                etemp = e
                e = d
                if abs(p) < abs(0.5 * q * etemp) and p > q * (a - x) and p < q * (b - x):
                    d = p / q
                    u = x + d
                    if (u - a) < tol2 or (b - u) < tol2:
                        d = tol1 if xm >= x else -tol1
                else:
                    e = (a - x) if x >= xm else (b - x)
                    d = gr * e
            else:
                e = (a - x) if x >= xm else (b - x)
                d = gr * e
            u = x + (d if abs(d) >= tol1 else (tol1 if d > 0 else -tol1))
            fu = f(u)
            if fu <= fx:
                if u >= x:
                    a = x
                else:
                    b = x
                v, fv = w, fw
                w, fw = x, fx
                x, fx = u, fu
            else:
                if u < x:
                    a = u
                else:
                    b = u
                if fu <= fw or w == x:
                    v, fv = w, fw
                    w, fw = u, fu
                elif fu <= fv or v == x or v == w:
                    v, fv = u, fu
        return x, fx

    def besselian_solar_state(self, jd_tt: float):
        r_sun, r_moon = self.apparent_sun_moon_eq_of_date(jd_tt)
        # Oblateness z-scale
        e2 = 2.0 * F_EARTH - F_EARTH ** 2
        rho_z = 1.0 / math.sqrt(1.0 - e2)
        r_s = np.array([r_sun[0], r_sun[1], r_sun[2] * rho_z]) / R_EARTH_KM
        r_m = np.array([r_moon[0], r_moon[1], r_moon[2] * rho_z]) / R_EARTH_KM

        g = r_s - r_m
        d_ms = float(np.linalg.norm(g))
        k_hat = g / d_ms
        a = math.atan2(k_hat[1], k_hat[0])
        i_hat = np.array([-math.sin(a), math.cos(a), 0.0])
        j_hat = np.cross(k_hat, i_hat)

        x = float(np.dot(r_m, i_hat))
        y = float(np.dot(r_m, j_hat))
        z = float(np.dot(r_m, k_hat))

        rs_er = R_SUN_KM / R_EARTH_KM
        rm_er = K_MOON_VALLEY
        sin_f1 = (rs_er + rm_er) / d_ms
        sin_f2 = (rs_er - rm_er) / d_ms
        tan_f1 = sin_f1 / math.sqrt(1.0 - sin_f1 ** 2)
        tan_f2 = sin_f2 / math.sqrt(1.0 - sin_f2 ** 2)

        c1 = z + rm_er / sin_f1
        c2 = z - rm_er / sin_f2
        l1 = c1 * tan_f1
        l2 = c2 * tan_f2

        gamma_abs = math.hypot(x, y)
        gamma = gamma_abs if y >= 0.0 else -gamma_abs
        zeta0 = math.sqrt(max(0.0, 1.0 - gamma_abs ** 2))
        L1 = l1 - zeta0 * tan_f1
        L2 = l2 - zeta0 * tan_f2

        return {
            "x": x,
            "y": y,
            "z": z,
            "gamma": gamma,
            "gamma_abs": gamma_abs,
            "l1": l1,
            "l2": l2,
            "L1": L1,
            "L2": L2,
            "tan_f1": tan_f1,
            "tan_f2": tan_f2,
        }

    def analyze_solar_eclipse(self, approx_jd_tt: float, label: str, expected_kind: str, nasa_td_str: str):
        jd_max, _ = self.brent_minimize(
            lambda jd: self.besselian_solar_state(jd)["gamma_abs"],
            approx_jd_tt - 0.5,
            approx_jd_tt + 0.5,
            tol=1e-9,
        )
        st = self.besselian_solar_state(jd_max)
        gamma = st["gamma"]
        g_abs = st["gamma_abs"]
        L1, L2 = st["L1"], st["L2"]

        if g_abs > 1.0 + abs(L2):
            kind = "Partial"
            mag = (1.0 + L1 - g_abs) / (L1 + L2)
        else:
            # Check L2 at limbs (zeta = 0 => l2) vs at greatest eclipse (L2)
            if st["l2"] > 0.0 and L2 < 0.0:
                kind = "Hybrid"
            elif L2 < 0.0:
                kind = "Total"
            else:
                kind = "Annular"
            mag = (L1 - L2) / (L1 + L2)

        assert kind == expected_kind, f"Expected {expected_kind} for {label}, got {kind}"
        return {
            "label": label,
            "nasa_canon_date_td": nasa_td_str,
            "kind": kind,
            "jd_tt_greatest": round(jd_max, 7),
            "gamma_earth_radii": round(gamma, 6),
            "magnitude": round(mag, 6),
            "fundamental_plane_l1_er": round(st["l1"], 6),
            "fundamental_plane_l2_er": round(st["l2"], 6),
            "surface_L1_er": round(L1, 6),
            "surface_L2_er": round(L2, 6),
        }

    def lunar_eclipse_state(self, jd_tt: float):
        r_sun, r_moon = self.apparent_sun_moon_eq_of_date(jd_tt)
        d_s = float(np.linalg.norm(r_sun))
        d_m = float(np.linalg.norm(r_moon))
        u_shadow = -r_sun / d_s
        u_moon = r_moon / d_m
        dot = max(-1.0, min(1.0, float(np.dot(u_shadow, u_moon))))
        sep_rad = math.acos(dot)

        pi_s = math.asin(R_EARTH_KM / d_s)
        s_s = math.asin(R_SUN_KM / d_s)
        pi_m = math.asin(R_EARTH_KM / d_m)
        s_m = math.asin(K_MOON_VALLEY * R_EARTH_KM / d_m)
        pi_1 = 0.998340 * pi_m

        # Chauvenet / Espenak 1% atmospheric enlargement
        f_sh = 1.01
        r_umbra = f_sh * (pi_1 + pi_s - s_s)
        r_penumbra = f_sh * (pi_1 + pi_s + s_s)

        mag_u = (r_umbra + s_m - sep_rad) / (2.0 * s_m)
        mag_p = (r_penumbra + s_m - sep_rad) / (2.0 * s_m)
        return {
            "sep_rad": sep_rad,
            "sep_deg": sep_rad * RAD2DEG,
            "r_umbra_deg": r_umbra * RAD2DEG,
            "r_penumbra_deg": r_penumbra * RAD2DEG,
            "s_moon_deg": s_m * RAD2DEG,
            "mag_umbra": mag_u,
            "mag_penumbra": mag_p,
        }

    def brent_root(self, f, a: float, b: float, tol: float = 1e-8) -> float:
        fa, fb = f(a), f(b)
        for _ in range(80):
            mid = 0.5 * (a + b)
            if abs(b - a) < tol:
                return mid
            fm = f(mid)
            if fa * fm <= 0.0:
                b, fb = mid, fm
            else:
                a, fa = mid, fm
        return 0.5 * (a + b)

    def analyze_lunar_eclipse(self, approx_jd_tt: float, label: str, expected_kind: str, nasa_td_str: str):
        jd_max, _ = self.brent_minimize(
            lambda jd: self.lunar_eclipse_state(jd)["sep_rad"],
            approx_jd_tt - 0.5,
            approx_jd_tt + 0.5,
            tol=1e-9,
        )
        st = self.lunar_eclipse_state(jd_max)
        mu = st["mag_umbra"]
        mp = st["mag_penumbra"]
        if mu >= 1.0:
            kind = "Total"
        elif mu > 0.0:
            kind = "Partial"
        else:
            kind = "Penumbral"
        assert kind == expected_kind, f"Expected {expected_kind} for {label}, got {kind}"

        # Contact times (P1, U1, U2, U3, U4, P4)
        def f_pen(jd):
            s = self.lunar_eclipse_state(jd)
            return s["sep_deg"] - (s["r_penumbra_deg"] + s["s_moon_deg"])

        def f_umb_ext(jd):
            s = self.lunar_eclipse_state(jd)
            return s["sep_deg"] - (s["r_umbra_deg"] + s["s_moon_deg"])

        def f_umb_int(jd):
            s = self.lunar_eclipse_state(jd)
            return s["sep_deg"] - (s["r_umbra_deg"] - s["s_moon_deg"])

        p1 = self.brent_root(f_pen, jd_max - 0.30, jd_max)
        p4 = self.brent_root(f_pen, jd_max, jd_max + 0.30)
        u1 = self.brent_root(f_umb_ext, jd_max - 0.25, jd_max) if mu > 0.0 else None
        u4 = self.brent_root(f_umb_ext, jd_max, jd_max + 0.25) if mu > 0.0 else None
        u2 = self.brent_root(f_umb_int, jd_max - 0.12, jd_max) if mu >= 1.0 else None
        u3 = self.brent_root(f_umb_int, jd_max, jd_max + 0.12) if mu >= 1.0 else None

        return {
            "label": label,
            "nasa_canon_date_td": nasa_td_str,
            "kind": kind,
            "jd_tt_greatest": round(jd_max, 7),
            "umbral_magnitude": round(mu, 6),
            "penumbral_magnitude": round(mp, 6),
            "min_separation_deg": round(st["sep_deg"], 6),
            "umbra_radius_deg": round(st["r_umbra_deg"], 6),
            "penumbra_radius_deg": round(st["r_penumbra_deg"], 6),
            "moon_semidiameter_deg": round(st["s_moon_deg"], 6),
            "contacts_jd_tt": {
                "P1": round(p1, 7),
                "U1": round(u1, 7) if u1 is not None else None,
                "U2": round(u2, 7) if u2 is not None else None,
                "greatest": round(jd_max, 7),
                "U3": round(u3, 7) if u3 is not None else None,
                "U4": round(u4, 7) if u4 is not None else None,
                "P4": round(p4, 7),
            },
        }


def main():
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)

    ev = De440sEclipseEvaluator(BSP_PATH)

    # 1. Generate l1b_houses_23.json with exact GAST (erfa.gst06a), True Obliquity, and DE440s Solar Declination
    raw_house_cases = [
        ("zurich_midlat_north", "Zurich, Switzerland (+47.37N, +8.55E) at J2000.0", 2451545.0, 47.37, 8.55),
        ("sydney_midlat_south", "Sydney, Australia (-33.87S, +151.21E) at 2024-06-21 06:00 UT1", 2460482.75, -33.87, 151.21),
        ("london_midlat_north", "London, UK (+51.51N, -0.13W) at 2010-03-20 12:00 UT1", 2455276.0, 51.51, -0.13),
        ("sao_paulo_south", "Sao Paulo, Brazil (-23.55S, -46.63W) at 1985-11-05 18:00 UT1", 2446375.25, -23.55, -46.63),
        ("antimeridian_pos180_north", "Antimeridian East (+45.00N, +180.00E) at J2000.0", 2451545.0, 45.0, 180.0),
        ("antimeridian_neg180_south", "Antimeridian West (-45.00S, -180.00W) at J2000.0", 2451545.0, -45.0, -180.0),
        ("quito_equator", "Quito, Ecuador (0.00N, -78.50W) at J2000.0", 2451545.0, 0.0, -78.50),
        ("tromso_arctic_north", "Tromso, Norway (+69.65N, +18.96E) above Arctic Circle at J2000.0", 2451545.0, 69.65, 18.96),
        ("mcmurdo_antarctic_south", "McMurdo Station, Antarctica (-77.85S, +166.67E) at 2024-06-21 06:00 UT1", 2460482.75, -77.85, 166.67),
        ("near_north_pole_89_9", "Near Geographic North Pole (+89.90N, 0.00E) at J2000.0", 2451545.0, 89.90, 0.0),
        ("exact_north_pole_90_0", "Exact Geographic North Pole (+90.00N, 0.00E) Singularity Test", 2451545.0, 90.0, 0.0),
    ]

    house_cases = []
    for cid, desc, jd_ut1, lat_deg, lon_deg in raw_house_cases:
        # Approximate delta-T (~64-69 s) for TT
        jd_tt = jd_ut1 + 65.0 / 86400.0
        gast_rad = float(erfa.gst06a(jd_ut1, 0.0, jd_tt, 0.0))
        ramc_deg = norm360(gast_rad * RAD2DEG + lon_deg)
        eps0_rad = float(erfa.obl06(jd_tt, 0.0))
        _, deps_rad = erfa.nut00a(jd_tt, 0.0)
        eps_true_deg = (eps0_rad + float(deps_rad)) * RAD2DEG
        r_sun, _ = ev.apparent_sun_moon_eq_of_date(jd_tt)
        sun_decl_deg = math.asin(float(r_sun[2]) / float(np.linalg.norm(r_sun))) * RAD2DEG
        house_cases.append({
            "case_id": cid,
            "description": desc,
            "jd_ut1": jd_ut1,
            "jd_tt": round(jd_tt, 8),
            "lat_deg": lat_deg,
            "lon_deg": lon_deg,
            "ramc_deg": round(ramc_deg, 10),
            "eps_true_deg": round(eps_true_deg, 10),
            "sun_decl_deg": round(sun_decl_deg, 10),
        })

    house_records = []
    for c in house_cases:
        is_pole = abs(c["lat_deg"]) >= 90.0
        angles = None if is_pole else compute_8_angles(c["ramc_deg"], c["lat_deg"], c["eps_true_deg"])
        systems = compute_all_23_houses(c["ramc_deg"], c["lat_deg"], c["eps_true_deg"], c["sun_decl_deg"])
        gauquelin = compute_gauquelin_36(c["ramc_deg"], c["lat_deg"], c["eps_true_deg"])
        house_records.append({
            **c,
            "angles": angles,
            "gauquelin_36_sectors": gauquelin,
            "systems": systems,
        })

    houses_out = {
        "schema_version": "1.0.0",
        "tier": "L1B_spherical_trigonometry",
        "source": "MODEL_PASSPORTS.md §5 (23 House Systems, 8 Chart Angles, 36 Gauquelin Sectors)",
        "case_count": len(house_records),
        "system_count": 23,
        "cases": house_records,
    }
    houses_path = GOLDEN_DIR / "l1b_houses_23.json"
    houses_path.write_text(json.dumps(houses_out, indent=2) + "\n")
    print(f"Wrote {houses_path} ({len(house_records)} cases x 23 house systems = {len(house_records)*23} system evaluations)")

    # 2. Generate l1b_eclipses.json
    ev = De440sEclipseEvaluator(BSP_PATH)
    solar_targets = [
        (2451401.964, "1999-08-11 Total Solar Eclipse (Europe)", "Total", "1999 Aug 11 11:04:09 TD"),
        (2453469.359, "2005-04-08 Hybrid Solar Eclipse (Pacific / Americas)", "Hybrid", "2005 Apr 08 20:36:51 TD"),
        (2456599.969, "2013-11-03 Hybrid Solar Eclipse (Atlantic / Africa)", "Hybrid", "2013 Nov 03 12:47:36 TD"),
        (2457987.269, "2017-08-21 Total Solar Eclipse (North America)", "Total", "2017 Aug 21 18:26:40 TD"),
        (2460054.679, "2023-04-20 Hybrid Solar Eclipse (Exmouth / Indonesia)", "Hybrid", "2023 Apr 20 04:17:56 TD"),
        (2460232.250, "2023-10-14 Annular Solar Eclipse (Americas)", "Annular", "2023 Oct 14 18:00:41 TD"),
        (2460409.263, "2024-04-08 Total Solar Eclipse (North America)", "Total", "2024 Apr 08 18:18:29 TD"),
        (2460763.950, "2025-03-29 Partial Solar Eclipse (North Atlantic / Europe)", "Partial", "2025 Mar 29 10:48:36 TD"),
    ]
    lunar_targets = [
        (2451742.081, "2000-07-16 Total Lunar Eclipse (Century-Max Duration)", "Total", "2000 Jul 16 13:56:39 TD"),
        (2458327.349, "2018-07-27 Total Lunar Eclipse (21st-Century Longest)", "Total", "2018 Jul 27 20:22:54 TD"),
        (2459537.878, "2021-11-19 Near-Total Partial Lunar Eclipse", "Partial", "2021 Nov 19 09:04:06 TD"),
        (2459891.958, "2022-11-08 Total Lunar Eclipse", "Total", "2022 Nov 08 11:00:22 TD"),
        (2460070.225, "2023-05-05 Deep Penumbral Lunar Eclipse", "Penumbral", "2023 May 05 17:24:05 TD"),
        (2460748.791, "2025-03-14 Total Lunar Eclipse", "Total", "2025 Mar 14 06:59:56 TD"),
    ]

    solar_records = [ev.analyze_solar_eclipse(*t) for t in solar_targets]
    lunar_records = [ev.analyze_lunar_eclipse(*t) for t in lunar_targets]

    eclipses_out = {
        "schema_version": "1.0.0",
        "tier": "L1B_besselian_eclipses",
        "source": "MODEL_PASSPORTS.md §6 + JPL DE440s + NASA Espenak & Meeus Five Millennium Canon",
        "solar_eclipse_count": len(solar_records),
        "lunar_eclipse_count": len(lunar_records),
        "solar_eclipses": solar_records,
        "lunar_eclipses": lunar_records,
    }
    eclipses_path = GOLDEN_DIR / "l1b_eclipses.json"
    eclipses_path.write_text(json.dumps(eclipses_out, indent=2) + "\n")
    print(f"Wrote {eclipses_path} ({len(solar_records)} solar + {len(lunar_records)} lunar = {len(solar_records)+len(lunar_records)} eclipses)")


if __name__ == "__main__":
    main()
