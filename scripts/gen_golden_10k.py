#!/usr/bin/env python3
"""
10,000+ Vector Suite Generator (`testdata/golden/l3_se_10k_grid.json`).

Generates 10,242 golden verification vectors across all domains:
  1. 10 Major Bodies (Sun..Pluto) x 500 epochs (1850..2100)       = 5,000 vectors
  2. 68 Minor Bodies (MPC 1..136472) x 20 epochs (1850..2100)     = 1,360 vectors
  3. 12 Ayanamsas x 100 epochs (1850..2100)                       = 1,200 vectors
  4. 23 House Systems x 60 (RAMC, lat, eps, sun_decl) states      = 1,380 vectors
  5. 7 Lunar Points x 40 epochs (1855..2145)                      =   280 vectors
  6. 9 Uranians (8 Hamburg + Transpluto) x 50 epochs (1850..2100) =   450 vectors
  7. 52 Fixed Stars (50 Astra + 2 anchors) x 11 epochs            =   572 vectors
Total in `l3_se_10k_grid.json`: 10,242 vectors.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import erfa
import numpy as np
from jplephem.spk import SPK
import swisseph as swe

from gen_golden_l1b import (
    AU_KM,
    C_AU_DAY,
    C_KM_S,
    DAS2R,
    MU_EM_KM3_S2,
    cart_to_spherical,
    iae_general_precession_arcsec,
    icrf_to_ecliptic_matrix,
    icrf_to_equatorial_matrix,
    interpolated_lilith_c2,
    lunar_osculating_node_and_lilith,
    planet_geocentric_vector_au,
    spk_barycentric_state_km_s,
    star_mean_ecliptic_lon_no_ab,
    tt_to_tdb_jd,
    uranian_heliocentric_icrf_au,
)
from gen_golden_houses_eclipses import compute_8_angles, compute_all_23_houses
from gen_golden_asteroids_68 import minor_body_apparent_ecl_au

ROOT = Path(__file__).resolve().parent.parent
DE440S_PATH = ROOT / "data" / "de440s.bsp"
AST_BSP_PATH = ROOT / "data" / "asteroids_68_1850_2100.bsp"
MANIFEST_PATH = ROOT / "data" / "manifest.json"
GOLDEN_DIR = ROOT / "testdata" / "golden"


def linspace(start: float, stop: float, n: int) -> list[float]:
    if n == 1:
        return [start]
    step = (stop - start) / (n - 1)
    return [round(start + i * step, 6) for i in range(n)]


def main() -> None:
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    spk_pl = SPK.open(str(DE440S_PATH))
    spk_ast = SPK.open(str(AST_BSP_PATH))
    manifest = json.loads(MANIFEST_PATH.read_text())
    bodies_meta = manifest["minor_bodies_kernel"]["bodies"]
    stars_list = json.loads((GOLDEN_DIR / "stars_catalog_52.json").read_text())
    stars_by_hip = {s["hip"]: s for s in stars_list}

    # -----------------------------------------------------------------------
    # 1. 10 Major Bodies x 500 epochs (1850..2100) = 5,000 vectors
    # -----------------------------------------------------------------------
    major_bodies = [
        (10, swe.SUN, "Sun"),
        (301, swe.MOON, "Moon"),
        (199, swe.MERCURY, "Mercury"),
        (299, swe.VENUS, "Venus"),
        (4, swe.MARS, "Mars"),
        (5, swe.JUPITER, "Jupiter"),
        (6, swe.SATURN, "Saturn"),
        (7, swe.URANUS, "Uranus"),
        (8, swe.NEPTUNE, "Neptune"),
        (9, swe.PLUTO, "Pluto"),
    ]
    epochs_500 = linspace(2396760.5, 2488065.5, 500)
    planetary_vectors = []
    se_flags = swe.FLG_SWIEPH | swe.FLG_SPEED

    for jd_tt in epochs_500:
        for naif_id, swe_id, name in major_bodies:
            ecl_0, _ = planet_geocentric_vector_au(spk_pl, naif_id, jd_tt, apparent=True, nutation=True)
            lon_l1b, lat_l1b, dist_l1b = cart_to_spherical(ecl_0)
            se_vals, _ = swe.calc(jd_tt, swe_id, se_flags)
            planetary_vectors.append([
                jd_tt,
                name,
                round(lon_l1b, 7),
                round(lat_l1b, 7),
                round(dist_l1b, 8),
                round(se_vals[0], 7),
                round(se_vals[1], 7),
                round(se_vals[3], 7),
            ])

    # -----------------------------------------------------------------------
    # 2. 68 Minor Bodies x 20 epochs (1850..2100) = 1,360 vectors
    # -----------------------------------------------------------------------
    epochs_20 = linspace(2396760.5, 2488065.5, 20)
    asteroid_vectors = []
    for jd_tt in epochs_20:
        for b in bodies_meta:
            ecl_0, astro_0 = minor_body_apparent_ecl_au(
                spk_pl, spk_ast, b["naif_spk_id"], b["center_id"], jd_tt
            )
            lon_deg, lat_deg, dist_au = cart_to_spherical(ecl_0)
            ra_icrf, dec_icrf, _ = cart_to_spherical(astro_0)
            asteroid_vectors.append([
                jd_tt,
                b["mpc"],
                b["name"],
                round(lon_deg, 7),
                round(lat_deg, 7),
                round(dist_au, 8),
                round(ra_icrf, 7),
                round(dec_icrf, 7),
            ])

    # -----------------------------------------------------------------------
    # 3. 12 Ayanamsas x 100 epochs (1850..2100) = 1,200 vectors
    # -----------------------------------------------------------------------
    epoch_models = [
        ("fagan_bradley", 0, 2433282.42345905, 24.042044444444444),
        ("lahiri", 1, 2435553.5, 23.250182777777778),
        ("raman", 3, 2415020.0, 21.013972222222223),
        ("ushashashi", 4, 2415020.0, 18.660722222222222),
        ("krishnamurti", 5, 2415020.0, 22.363421222222222),
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
    epochs_100 = linspace(2396760.5, 2488065.5, 100)
    ayanamsa_vectors = []
    for jd_tt in epochs_100:
        dpsi_rad, _ = erfa.nut00a(jd_tt, 0.0)
        dpsi_deg = math.degrees(float(dpsi_rad))
        for code, legacy_id, t0, a0_deg in epoch_models:
            p_arcsec, _ = iae_general_precession_arcsec(t0, jd_tt)
            mean_deg = (a0_deg + p_arcsec / 3600.0) % 360.0
            true_deg = (mean_deg + dpsi_deg) % 360.0
            swe.set_sid_mode(legacy_id, 0, 0)
            se_mean = swe.get_ayanamsa_ex(jd_tt, swe.FLG_SWIEPH | swe.FLG_NONUT)[1]
            ayanamsa_vectors.append([
                jd_tt,
                code,
                round(mean_deg, 7),
                round(true_deg, 7),
                round(se_mean, 7),
            ])
        for code, legacy_id, hip_id, target_lon in star_models:
            star = stars_by_hip[hip_id]
            lon_0 = star_mean_ecliptic_lon_no_ab(spk_pl, star, jd_tt)
            mean_deg = (lon_0 - target_lon) % 360.0
            true_deg = (mean_deg + dpsi_deg) % 360.0
            swe.set_sid_mode(legacy_id, 0, 0)
            se_mean = swe.get_ayanamsa_ex(jd_tt, swe.FLG_SWIEPH | swe.FLG_NONUT)[1]
            ayanamsa_vectors.append([
                jd_tt,
                code,
                round(mean_deg, 7),
                round(true_deg, 7),
                round(se_mean, 7),
            ])

    # -----------------------------------------------------------------------
    # 4. 23 House Systems x 60 (RAMC, lat, eps, sun_decl) states = 1,380 vectors
    # -----------------------------------------------------------------------
    # 12 RAMC steps (0..330 deg) x 5 latitudes (-50, -25, +15, +35, +55 deg) = 60 states
    ramc_grid = [30.0 * k for k in range(12)]
    lat_grid = [-50.0, -25.0, 15.0, 35.0, 55.0]
    eps_val = 23.4392911
    sun_decl_val = 12.5
    house_vectors = []
    for lat in lat_grid:
        for ramc in ramc_grid:
            all_23 = compute_all_23_houses(ramc, lat, eps_val, sun_decl_val)
            for code, sys_data in all_23.items():
                house_vectors.append([
                    ramc,
                    lat,
                    code,
                    [round(c, 6) for c in sys_data["cusps"]],
                ])

    # -----------------------------------------------------------------------
    # 5. 7 Lunar Points x 40 epochs (1855..2145) = 280 vectors
    # -----------------------------------------------------------------------
    epochs_40 = linspace(2398736.5, 2504656.5, 40)
    lunar_vectors = []
    for jd_tt in epochs_40:
        T = (jd_tt - 2451545.0) / 36525.0
        dpsi_rad, _ = erfa.nut00a(jd_tt, 0.0)
        dpsi_deg = math.degrees(float(dpsi_rad))
        om_true = (math.degrees(float(erfa.faom03(T))) + dpsi_deg) % 360.0
        mlil_true = (
            math.degrees(float(erfa.faf03(T) + erfa.faom03(T) - erfa.fal03(T))) + 180.0 + dpsi_deg
        ) % 360.0
        tn_0, tl_0, tlat_0 = lunar_osculating_node_and_lilith(spk_pl, jd_tt, nutation=True)
        il_lon, il_lat, il_dlon, _ = interpolated_lilith_c2(spk_pl, jd_tt, nutation=True)
        se_intp, _ = swe.calc(jd_tt, swe.INTP_APOG, se_flags)

        pts = [
            ("Mean_Node", om_true, 0.0),
            ("Mean_South_Node", (om_true + 180.0) % 360.0, 0.0),
            ("True_Node", tn_0, 0.0),
            ("True_South_Node", (tn_0 + 180.0) % 360.0, 0.0),
            ("Mean_Lilith", mlil_true, 0.0),
            ("True_Lilith", tl_0, tlat_0),
            ("Interpolated_Lilith", il_lon, il_lat),
        ]
        for pname, plon, plat in pts:
            lunar_vectors.append([
                jd_tt,
                pname,
                round(plon, 6),
                round(plat, 6),
                round(se_intp[0], 6) if pname == "Interpolated_Lilith" else None,
            ])

    # -----------------------------------------------------------------------
    # 6. 9 Uranians x 50 epochs (1850..2100) = 450 vectors
    # -----------------------------------------------------------------------
    uranian_elements = [
        {"name": "Cupido", "swe_id": swe.CUPIDO, "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 163.7409, "a": 40.99837, "e": 0.00460, "w": 171.4333, "om": 129.8325, "i": 1.0833},
        {"name": "Hades", "swe_id": swe.HADES, "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 27.6496, "a": 50.66744, "e": 0.00245, "w": 148.1796, "om": 161.3339, "i": 1.0500},
        {"name": "Zeus", "swe_id": swe.ZEUS, "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 165.1232, "a": 59.21436, "e": 0.00120, "w": 299.0440, "om": 0.0, "i": 0.0},
        {"name": "Kronos", "swe_id": swe.KRONOS, "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 169.0193, "a": 64.81690, "e": 0.00305, "w": 208.8801, "om": 0.0, "i": 0.0},
        {"name": "Apollon", "swe_id": swe.APOLLON, "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 138.0533, "a": 70.29949, "e": 0.0, "w": 0.0, "om": 0.0, "i": 0.0},
        {"name": "Admetos", "swe_id": swe.ADMETOS, "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 351.3350, "a": 73.62765, "e": 0.0, "w": 0.0, "om": 0.0, "i": 0.0},
        {"name": "Vulkanus", "swe_id": swe.VULKANUS, "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 55.8983, "a": 77.25568, "e": 0.0, "w": 0.0, "om": 0.0, "i": 0.0},
        {"name": "Poseidon", "swe_id": swe.POSEIDON, "epoch_jd": 2415020.0, "equinox_jd": 2415020.0, "M0": 165.5163, "a": 83.66907, "e": 0.0, "w": 0.0, "om": 0.0, "i": 0.0},
        {"name": "Transpluto", "swe_id": swe.ISIS, "epoch_jd": 2368547.66, "equinox_jd": 2431456.5, "M0": 0.0, "a": 77.77500, "e": 0.30000, "w": 0.7000, "om": 0.0, "i": 0.0},
    ]
    epochs_50 = linspace(2396760.5, 2488065.5, 50)
    uranian_vectors = []
    for jd_tt in epochs_50:
        t_tdb = tt_to_tdb_jd(jd_tt)
        p_earth_km, v_earth_kms = spk_barycentric_state_km_s(spk_pl, 399, t_tdb)
        p_sun_km, _ = spk_barycentric_state_km_s(spk_pl, 10, t_tdb)
        r_earth_helio_au = (p_earth_km - p_sun_km) / AU_KM
        bm1 = math.sqrt(1.0 - min(float(np.dot(v_earth_kms, v_earth_kms)) / (C_KM_S * C_KM_S), 0.999999))
        r_ecl = icrf_to_ecliptic_matrix(jd_tt, nutation=True)

        for el in uranian_elements:
            dt = 0.0
            for _ in range(3):
                r_obj_helio_au = uranian_heliocentric_icrf_au(el, jd_tt - dt)
                rho_au = r_obj_helio_au - r_earth_helio_au
                dist_au = float(np.linalg.norm(rho_au))
                dt = dist_au / C_AU_DAY
            u_unit = rho_au / dist_au
            u_unit = np.asarray(erfa.ab(u_unit, v_earth_kms / C_KM_S, float(np.linalg.norm(r_earth_helio_au)), bm1), dtype=float)
            lon_l1b, lat_l1b, _ = cart_to_spherical(r_ecl @ (u_unit * dist_au))
            se_u, _ = swe.calc(jd_tt, el["swe_id"], se_flags)
            uranian_vectors.append([
                jd_tt,
                el["name"],
                round(lon_l1b, 6),
                round(lat_l1b, 6),
                round(se_u[0], 6),
            ])

    # -----------------------------------------------------------------------
    # 7. 52 Fixed Stars x 11 epochs (1850..2100) = 572 vectors
    # -----------------------------------------------------------------------
    epochs_11 = linspace(2396760.5, 2488065.5, 11)
    star_vectors = []
    for jd_tt in epochs_11:
        jd_tdb = tt_to_tdb_jd(jd_tt)
        p_earth_km, v_earth_kms = spk_barycentric_state_km_s(spk_pl, 399, jd_tdb)
        p_sun_km, _ = spk_barycentric_state_km_s(spk_pl, 10, jd_tdb)
        p_eb_au = p_earth_km / AU_KM
        e_sun_vec = (p_earth_km - p_sun_km) / AU_KM
        em = float(np.linalg.norm(e_sun_vec))
        e_sun_unit = e_sun_vec / em
        bm1 = math.sqrt(1.0 - min(float(np.dot(v_earth_kms, v_earth_kms)) / (C_KM_S * C_KM_S), 0.999999))
        v_ab = v_earth_kms / C_KM_S
        r_ecl_true = icrf_to_ecliptic_matrix(jd_tt, nutation=True)
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
            star_vectors.append([
                jd_tt,
                s["hip"],
                s["name"],
                round(lon_deg, 6),
                round(lat_deg, 6),
            ])

    total_count = (
        len(planetary_vectors)
        + len(asteroid_vectors)
        + len(ayanamsa_vectors)
        + len(house_vectors)
        + len(lunar_vectors)
        + len(uranian_vectors)
        + len(star_vectors)
    )

    out = {
        "schema_version": "1.0.0",
        "tier": "L1B_and_L3_10k_grid",
        "total_vector_count": total_count,
        "counts_by_section": {
            "major_bodies_10x500": len(planetary_vectors),
            "minor_bodies_68x20": len(asteroid_vectors),
            "ayanamsas_12x100": len(ayanamsa_vectors),
            "houses_23x60": len(house_vectors),
            "lunar_points_7x40": len(lunar_vectors),
            "uranians_9x50": len(uranian_vectors),
            "fixed_stars_52x11": len(star_vectors),
        },
        "columns": {
            "major_bodies": ["jd_tt", "body", "l1b_lon_deg", "l1b_lat_deg", "l1b_dist_au", "l3_se_lon_deg", "l3_se_lat_deg", "l3_se_speed_deg_day"],
            "minor_bodies": ["jd_tt", "mpc", "name", "l1b_lon_deg", "l1b_lat_deg", "l1b_dist_au", "icrf_ra_deg", "icrf_dec_deg"],
            "ayanamsas": ["jd_tt", "code", "l1b_mean_deg", "l1b_true_deg", "l3_se_mean_deg"],
            "houses": ["ramc_deg", "lat_deg", "system_code", "cusps_1_to_12_deg"],
            "lunar_points": ["jd_tt", "point", "l1b_lon_deg", "l1b_lat_deg", "l3_se_intp_lon_deg_if_applicable"],
            "uranians": ["jd_tt", "name", "l1b_lon_deg", "l1b_lat_deg", "l3_se_lon_deg"],
            "fixed_stars": ["jd_tt", "hip", "name", "l1b_lon_deg", "l1b_lat_deg"],
        },
        "major_bodies": planetary_vectors,
        "minor_bodies": asteroid_vectors,
        "ayanamsas": ayanamsa_vectors,
        "houses": house_vectors,
        "lunar_points": lunar_vectors,
        "uranians": uranian_vectors,
        "fixed_stars": star_vectors,
    }

    out_path = GOLDEN_DIR / "l3_se_10k_grid.json"
    out_path.write_text(json.dumps(out, separators=(",", ":")) + "\n")
    size_kb = out_path.stat().st_size / 1024.0
    print(f"Wrote {out_path} ({total_count} vectors, {size_kb:.1f} KB)")
    for k, v in out["counts_by_section"].items():
        print(f"  - {k}: {v}")


if __name__ == "__main__":
    main()
