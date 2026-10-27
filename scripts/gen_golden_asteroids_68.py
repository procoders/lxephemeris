#!/usr/bin/env python3
"""
Level 1B Golden Vector Generator for all 68 Minor Bodies (`l1b_asteroids_68.json`).

Evaluates all 68 minor bodies (MPC 1 Ceres .. MPC 136472 Makemake) from
`data/asteroids_68_1850_2100.bsp` + `data/de440s.bsp` across 5 epochs in `[1850, 2100)`:
  - 1850-01-01T00:00:00 TT (JD 2396758.5)
  - 1900-01-01T00:00:00 TT (JD 2415020.5)
  - 2000-01-01T12:00:00 TT (JD 2451545.0)
  - 2025-01-01T00:00:00 TT (JD 2460676.5)
  - 2099-12-31T00:00:00 TT (JD 2488068.5)
Total: 68 x 5 = 340 golden vectors.
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
AST_BSP_PATH = ROOT / "data" / "asteroids_68_1850_2100.bsp"
MANIFEST_PATH = ROOT / "data" / "manifest.json"
GOLDEN_DIR = ROOT / "testdata" / "golden"

AU_KM = 149597870.700
C_KM_S = 299792.458
C_AU_DAY = (C_KM_S * 86400.0) / AU_KM


def tt_to_tdb_jd(jd_tt: float) -> float:
    g = math.radians((357.53 + 0.98560028 * (jd_tt - 2451545.0)) % 360.0)
    dt_sec = erfa.dtdb(jd_tt, 0.0, 0.0, 0.0, 0.0, g)
    return jd_tt + dt_sec / 86400.0


def spk_earth_and_sun_bary_km_s(spk_pl: SPK, jd_tdb: float):
    p_emb, v_emb = spk_pl[0, 3].compute_and_differentiate(jd_tdb)
    p_e, v_e = spk_pl[3, 399].compute_and_differentiate(jd_tdb)
    p_s, v_s = spk_pl[0, 10].compute_and_differentiate(jd_tdb)
    return (p_emb + p_e), (v_emb + v_e) / 86400.0, p_s, v_s / 86400.0


def spk_minor_body_bary_km(spk_pl: SPK, spk_ast: SPK, naif_id: int, center_id: int, jd_tdb: float) -> np.ndarray:
    p_rel = spk_ast[center_id, naif_id].compute(jd_tdb)
    if center_id == 0:
        return p_rel
    p_sun = spk_pl[0, center_id].compute(jd_tdb)
    return p_sun + p_rel


def icrf_to_ecliptic_matrix(jd_tt: float, nutation: bool = True) -> np.ndarray:
    eps0 = float(erfa.obl06(jd_tt, 0.0))
    if nutation:
        _, deps = erfa.nut00a(jd_tt, 0.0)
        eps = eps0 + float(deps)
        r_eq = erfa.num06a(jd_tt, 0.0) @ erfa.bp06(jd_tt, 0.0)[2]
    else:
        eps = eps0
        r_eq = erfa.bp06(jd_tt, 0.0)[2]
    ce, se = math.cos(eps), math.sin(eps)
    r_x = np.array([[1.0, 0.0, 0.0], [0.0, ce, se], [0.0, -se, ce]])
    return r_x @ r_eq


def cart_to_spherical(u: np.ndarray) -> tuple[float, float, float]:
    x, y, z = float(u[0]), float(u[1]), float(u[2])
    r = math.hypot(x, y, z)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(max(-1.0, min(1.0, z / r)))) if r > 0.0 else 0.0
    return lon, lat, r


def cart_state_to_spherical_rates(u: np.ndarray, du: np.ndarray) -> tuple[float, float, float]:
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


def minor_body_apparent_ecl_au(
    spk_pl: SPK, spk_ast: SPK, naif_id: int, center_id: int, jd_tt: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return (apparent_ecl_true_au, astrometric_icrf_au) for a minor body."""
    jd_tdb = tt_to_tdb_jd(jd_tt)
    p_earth_km, v_earth_kms, p_sun_now_km, _ = spk_earth_and_sun_bary_km_s(spk_pl, jd_tdb)

    dt_days = 0.0
    for _ in range(5):
        p_targ_km = spk_minor_body_bary_km(spk_pl, spk_ast, naif_id, center_id, jd_tdb - dt_days)
        rho_au = (p_targ_km - p_earth_km) / AU_KM
        dist_au = float(np.linalg.norm(rho_au))
        dt_days = dist_au / C_AU_DAY

    astro_icrf_au = rho_au.copy()
    u_unit = rho_au / dist_au

    # Relativistic solar light deflection + stellar aberration
    p_sun_ret_km = spk_pl[0, 10].compute(jd_tdb - dt_days)
    q_vec = (p_targ_km - p_sun_ret_km) / AU_KM
    e_vec = (p_earth_km - p_sun_ret_km) / AU_KM
    em = float(np.linalg.norm(e_vec))
    q_unit = q_vec / float(np.linalg.norm(q_vec))
    e_unit = e_vec / em
    u_unit = np.asarray(erfa.ld(1.0, u_unit, q_unit, e_unit, em, 0.0), dtype=float)
    u_unit = u_unit / float(np.linalg.norm(u_unit))

    bm1 = math.sqrt(1.0 - min(float(np.dot(v_earth_kms, v_earth_kms)) / (C_KM_S * C_KM_S), 0.999999))
    v_ab = v_earth_kms / C_KM_S
    s_dist_au = float(np.linalg.norm((p_earth_km - p_sun_now_km) / AU_KM))
    u_unit = np.asarray(erfa.ab(u_unit, v_ab, s_dist_au, bm1), dtype=float)

    r_ecl = icrf_to_ecliptic_matrix(jd_tt, nutation=True)
    return r_ecl @ (u_unit * dist_au), astro_icrf_au


def main() -> None:
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST_PATH.read_text())
    bodies_meta = manifest["minor_bodies_kernel"]["bodies"]
    assert len(bodies_meta) == 68, f"Expected 68 minor bodies, found {len(bodies_meta)}"

    spk_pl = SPK.open(str(DE440S_PATH))
    spk_ast = SPK.open(str(AST_BSP_PATH))

    epochs = [
        ("1850-01-01T00:00:00", 2396758.5),
        ("1900-01-01T00:00:00", 2415020.5),
        ("2000-01-01T12:00:00", 2451545.0),
        ("2025-01-01T00:00:00", 2460676.5),
        ("2099-12-31T00:00:00", 2488068.5),
    ]

    h = 0.25  # 6-hour stencil step (stays well within raw kernel [2396752.5, 2488073.5])
    vectors = []

    for ep_label, jd_tt in epochs:
        for b in bodies_meta:
            naif_id = b["naif_spk_id"]
            center_id = b["center_id"]
            ecl_m2, _ = minor_body_apparent_ecl_au(spk_pl, spk_ast, naif_id, center_id, jd_tt - 2 * h)
            ecl_m1, _ = minor_body_apparent_ecl_au(spk_pl, spk_ast, naif_id, center_id, jd_tt - h)
            ecl_0, astro_0 = minor_body_apparent_ecl_au(spk_pl, spk_ast, naif_id, center_id, jd_tt)
            ecl_p1, _ = minor_body_apparent_ecl_au(spk_pl, spk_ast, naif_id, center_id, jd_tt + h)
            ecl_p2, _ = minor_body_apparent_ecl_au(spk_pl, spk_ast, naif_id, center_id, jd_tt + 2 * h)
            decl = (ecl_m2 - 8.0 * ecl_m1 + 8.0 * ecl_p1 - ecl_p2) / (12.0 * h)

            lon_deg, lat_deg, dist_au = cart_to_spherical(ecl_0)
            dlon, dlat, ddist = cart_state_to_spherical_rates(ecl_0, decl)
            ra_icrf, dec_icrf, dist_icrf_au = cart_to_spherical(astro_0)

            vectors.append({
                "epoch_label": ep_label,
                "jd_tt": jd_tt,
                "mpc": b["mpc"],
                "name": b["name"],
                "category": b["category"],
                "naif_spk_id": naif_id,
                "apparent_ecliptic_true_of_date": {
                    "lon_deg": round(lon_deg, 9),
                    "lat_deg": round(lat_deg, 9),
                    "dist_au": round(dist_au, 9),
                    "lon_speed_deg_day": round(dlon, 9),
                    "lat_speed_deg_day": round(dlat, 9),
                    "dist_speed_au_day": round(ddist, 9),
                },
                "astrometric_icrf": {
                    "ra_deg": round(ra_icrf, 9),
                    "dec_deg": round(dec_icrf, 9),
                    "dist_au": round(dist_icrf_au, 9),
                },
            })

    out = {
        "schema_version": "1.0.0",
        "level": "L1B",
        "source": "jplephem + pyerfa (IAU 2006/2000A) + 4th-order Cartesian stencil on asteroids_68_1850_2100.bsp + de440s.bsp",
        "ephemeris": "asteroids_68_1850_2100.bsp + de440s.bsp",
        "time_scale": "TT",
        "tolerance_arcsec": 0.01,
        "body_count": len(bodies_meta),
        "epoch_count": len(epochs),
        "count": len(vectors),
        "vectors": vectors,
    }
    out_path = GOLDEN_DIR / "l1b_asteroids_68.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Wrote {out_path} ({len(vectors)} vectors across {len(bodies_meta)} minor bodies x {len(epochs)} epochs)")


if __name__ == "__main__":
    main()
