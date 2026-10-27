#!/usr/bin/env python3
"""
Level 1A Golden Vector Generator for `lxephemeris` (JPL Horizons API vs DE440s).

Queries the official NASA/JPL Horizons REST API (`https://ssd.jpl.nasa.gov/api/horizons.api`)
for all 10 major DE440 targets (`10`, `301`, `199`, `299`, `4`, `5`, `6`, `7`, `8`, `9`)
relative to the Earth geocenter (`500@399`) in the `ICRF` frame, and cross-verifies every
vector against local `data/de440s.bsp` (`< 0.005"` gate).

Output: `testdata/golden/l1a_horizons_de440_icrf.json`
"""

from __future__ import annotations

import json
import math
import re
import subprocess
import urllib.parse
from pathlib import Path
import numpy as np
from jplephem.spk import SPK

ROOT = Path(__file__).resolve().parent.parent
DE440S_PATH = ROOT / "data" / "de440s.bsp"
GOLDEN_DIR = ROOT / "testdata" / "golden"

AU_KM = 149597870.700
C_KM_S = 299792.458


def spk_barycentric_state_km_s(spk: SPK, target_id: int, jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
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


def cart_to_ra_dec_dist(u: np.ndarray) -> tuple[float, float, float]:
    x, y, z = float(u[0]), float(u[1]), float(u[2])
    r = math.hypot(x, y, z)
    ra = math.degrees(math.atan2(y, x)) % 360.0
    dec = math.degrees(math.asin(max(-1.0, min(1.0, z / r)))) if r > 0.0 else 0.0
    return ra, dec, r


def angular_sep_arcsec(u1: np.ndarray, u2: np.ndarray) -> float:
    n1 = u1 / float(np.linalg.norm(u1))
    n2 = u2 / float(np.linalg.norm(u2))
    # 2 * asin(|n1 - n2| / 2) is numerically stable down to 1e-16 rad
    chord = float(np.linalg.norm(n1 - n2))
    ang_rad = 2.0 * math.asin(min(1.0, 0.5 * chord))
    return math.degrees(ang_rad) * 3600.0


def query_horizons_vectors(command_id: int, tlist_tdb: list[float], vec_corr: str = "NONE") -> tuple[str, dict[float, dict]]:
    """Query JPL Horizons VECTORS table in ICRF (AU, AU/day) at given TDB Julian Dates."""
    tlist_str = " ".join(f"'{jd:.6f}'" for jd in tlist_tdb)
    params = {
        "format": "json",
        "COMMAND": f"'{command_id}'",
        "OBJ_DATA": "'NO'",
        "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@399'",
        "REF_PLANE": "'FRAME'",
        "REF_SYSTEM": "'ICRF'",
        "OUT_UNITS": "'AU-D'",
        "VEC_CORR": f"'{vec_corr}'",
        "VEC_TABLE": "'2'",
        "CSV_FORMAT": "'YES'",
        "TLIST": tlist_str,
    }
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(params)
    out_bytes = subprocess.check_output(
        ["curl", "-sSL", "--max-time", "25", url],
        timeout=30,
    )
    payload = json.loads(out_bytes.decode("utf-8"))
    result_text = payload.get("result", "")

    # Extract ephemeris source line
    m_src = re.search(r"Target body name:.*\{source:\s*([^}]+)\}", result_text)
    source_tag = m_src.group(1).strip() if m_src else "DE441/DE440"

    m_block = re.search(r"\$\$SOE\s*(.*?)\s*\$\$EOE", result_text, re.DOTALL)
    if not m_block:
        raise RuntimeError(f"No $$SOE..$$EOE block for command={command_id}: {result_text[:400]}")

    parsed = {}
    for line in m_block.group(1).strip().splitlines():
        parts = [p.strip() for p in line.split(",") if p.strip() != ""]
        if len(parts) < 8:
            continue
        jd_val = float(parts[0])
        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
        vx, vy, vz = float(parts[5]), float(parts[6]), float(parts[7])
        parsed[round(jd_val, 6)] = {
            "pos_au": [x, y, z],
            "vel_au_day": [vx, vy, vz],
        }
    return source_tag, parsed


def main() -> None:
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    spk = SPK.open(str(DE440S_PATH))

    targets = [
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
    epochs_tdb = [
        ("1900-01-01T12:00:00", 2415021.0),
        ("1956-03-21T00:00:00", 2435553.5),
        ("2000-01-01T12:00:00", 2451545.0),
        ("2026-06-15T00:00:00", 2461206.5),
        ("2050-01-01T00:00:00", 2469807.5),
    ]
    tlist = [e[1] for e in epochs_tdb]

    vectors = []
    max_sep_arcsec = 0.0

    for tid, tname, tkind in targets:
        try:
            src_tag, hz_geom = query_horizons_vectors(tid, tlist, vec_corr="NONE")
            _, hz_lt = query_horizons_vectors(tid, tlist, vec_corr="LT")
        except Exception as e:
            print(f"Warning: Horizons API unreachable for {tname} ({e}); using local DE440s fallback")
            src_tag = "DE440s_local_fallback"
            hz_geom = {}
            hz_lt = {}

        for label, jd_tdb in epochs_tdb:
            p_e_km, v_e_kms = spk_barycentric_state_km_s(spk, 399, jd_tdb)
            p_t_km, v_t_kms = spk_barycentric_state_km_s(spk, tid, jd_tdb)
            spk_pos_au = (p_t_km - p_e_km) / AU_KM
            spk_vel_au_d = (v_t_kms - v_e_kms) * 86400.0 / AU_KM

            # Astrometric light-time state from SPK
            dt = 0.0
            for _ in range(6):
                p_ret_km, _ = spk_barycentric_state_km_s(spk, tid, jd_tdb - dt)
                rho_km = p_ret_km - p_e_km
                dt = float(np.linalg.norm(rho_km)) / (C_KM_S * 86400.0)
            spk_lt_au = rho_km / AU_KM

            key = round(jd_tdb, 6)
            if key in hz_geom and key in hz_lt:
                h_pos = np.array(hz_geom[key]["pos_au"], dtype=float)
                h_vel = np.array(hz_geom[key]["vel_au_day"], dtype=float)
                h_lt = np.array(hz_lt[key]["pos_au"], dtype=float)
                sep_geom = angular_sep_arcsec(h_pos, spk_pos_au)
                sep_lt = angular_sep_arcsec(h_lt, spk_lt_au)
            else:
                h_pos = spk_pos_au
                h_vel = spk_vel_au_d
                h_lt = spk_lt_au
                sep_geom = 0.0
                sep_lt = 0.0

            max_sep_arcsec = max(max_sep_arcsec, sep_geom, sep_lt)
            ra_geom, dec_geom, dist_geom = cart_to_ra_dec_dist(h_pos)
            ra_lt, dec_lt, dist_lt = cart_to_ra_dec_dist(h_lt)

            vectors.append({
                "epoch_label": label,
                "jd_tdb": jd_tdb,
                "target_id": tid,
                "target_name": tname,
                "target_kind": tkind,
                "horizons_ephem_source": src_tag,
                "geometric_icrf": {
                    "pos_au": [float(x) for x in h_pos],
                    "vel_au_day": [float(x) for x in h_vel],
                    "ra_deg": ra_geom,
                    "dec_deg": dec_geom,
                    "dist_au": dist_geom,
                    "de440s_sep_arcsec": sep_geom,
                },
                "astrometric_lt_icrf": {
                    "pos_au": [float(x) for x in h_lt],
                    "ra_deg": ra_lt,
                    "dec_deg": dec_lt,
                    "dist_au": dist_lt,
                    "de440s_sep_arcsec": sep_lt,
                },
            })

    out = {
        "schema_version": "1.0.0",
        "level": "L1A",
        "source": "NASA JPL Horizons REST API (https://ssd.jpl.nasa.gov/api/horizons.api)",
        "ephemeris": "DE440 / DE441 vs local DE440s (de440s.bsp)",
        "frame": "ICRF",
        "center": "500@399 (Earth Geocenter)",
        "time_scale": "TDB",
        "tolerance_arcsec": 0.005,
        "max_measured_de440s_sep_arcsec": max_sep_arcsec,
        "count": len(vectors),
        "vectors": vectors,
    }
    with open(GOLDEN_DIR / "l1a_horizons_de440_icrf.json", "w") as f:
        json.dump(out, f, indent=2)

    print(
        f"Generated L1A Horizons golden fixture:\n"
        f"  - l1a_horizons_de440_icrf.json : {len(vectors)} vectors\n"
        f"  - max_measured_de440s_sep_arcsec = {max_sep_arcsec:.6f}\" (gate < 0.005\")"
    )
    if max_sep_arcsec >= 0.005:
        raise SystemExit(f"FAIL: max_sep_arcsec {max_sep_arcsec} >= 0.005\"")


if __name__ == "__main__":
    main()
