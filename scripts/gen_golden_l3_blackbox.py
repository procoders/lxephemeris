#!/usr/bin/env python3
"""
Black-Box Reference Oracle Generator (L3 Tier — ARCHITECTURE.md §1.2.2b)
Outputs: testdata/golden/l3_se_reference.json

Strictly invokes `pyswisseph` as an external black-box binary via its public Python API
to record numerical outputs for cross-validation against our Clean-Room L1A/L1B models
and to verify all documented tolerances in `docs/SE_DIVERGENCES.md`.
"""

import json
import math
from pathlib import Path
import swisseph as swe

ROOT = Path(__file__).resolve().parent.parent
GOLDEN_DIR = ROOT / "testdata" / "golden"


def norm360(x: float) -> float:
    return x % 360.0


def wrap180(x: float) -> float:
    return (x + 180.0) % 360.0 - 180.0


def main():
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)

    l1b_ayanamsas = json.loads((GOLDEN_DIR / "l1b_ayanamsas_12.json").read_text())
    l1b_nodes = json.loads((GOLDEN_DIR / "l1b_lunar_points_7.json").read_text())
    l1b_houses = json.loads((GOLDEN_DIR / "l1b_houses_23.json").read_text())
    l1b_eclipses = json.loads((GOLDEN_DIR / "l1b_eclipses.json").read_text())

    # 1. Major Bodies (Apparent Tropical Ecliptic of Date + Speed)
    epochs = [
        ("1900-01-01T12:00:00_TT", 2415021.0),
        ("1950-01-01T00:00:00_TT", 2433282.5),
        ("2000-01-01T12:00:00_TT", 2451545.0),
        ("2025-06-21T00:00:00_TT", 2460847.5),
        ("2050-01-01T12:00:00_TT", 2469807.0),
    ]
    bodies = [
        ("Sun", swe.SUN),
        ("Moon", swe.MOON),
        ("Mercury", swe.MERCURY),
        ("Venus", swe.VENUS),
        ("Mars", swe.MARS),
        ("Jupiter", swe.JUPITER),
        ("Saturn", swe.SATURN),
        ("Uranus", swe.URANUS),
        ("Neptune", swe.NEPTUNE),
        ("Pluto", swe.PLUTO),
    ]
    major_body_vectors = []
    flags = swe.FLG_SWIEPH | swe.FLG_SPEED
    for ep_label, jd_tt in epochs:
        for b_name, b_id in bodies:
            vals, ret_flag = swe.calc(jd_tt, b_id, flags)
            major_body_vectors.append({
                "epoch_label": ep_label,
                "jd_tt": jd_tt,
                "body": b_name,
                "lon_deg": round(vals[0], 10),
                "lat_deg": round(vals[1], 10),
                "dist_au": round(vals[2], 10),
                "lon_speed_deg_day": round(vals[3], 10),
                "lat_speed_deg_day": round(vals[4], 10),
                "dist_speed_au_day": round(vals[5], 10),
                "ret_flag": int(ret_flag),
            })

    # 2. Lunar Nodes & 3 Liliths (Mean_Node, True_Node, Mean_Lilith, True_Lilith, Interpolated_Lilith)
    node_points = [
        ("Mean_Node", swe.MEAN_NODE),
        ("True_Node", swe.TRUE_NODE),
        ("Mean_Lilith", swe.MEAN_APOG),
        ("True_Lilith", swe.OSCU_APOG),
        ("Interpolated_Lilith", swe.INTP_APOG),
    ]
    l3_nodes_liliths = []
    intp_diffs_arcmin = []
    for rec in l1b_nodes["vectors"]:
        jd_tt = rec["jd_tt"]
        entry = {"epoch_label": rec["epoch_label"], "jd_tt": jd_tt, "points": {}}
        for p_name, p_id in node_points:
            vals, _ = swe.calc(jd_tt, p_id, flags)
            l1b_pt = rec["points"][p_name]
            diff_arcsec = wrap180(l1b_pt["lon_deg"] - vals[0]) * 3600.0
            if p_name == "Interpolated_Lilith":
                intp_diffs_arcmin.append(diff_arcsec / 60.0)
            entry["points"][p_name] = {
                "l3_se_lon_deg": round(vals[0], 9),
                "l3_se_lat_deg": round(vals[1], 9),
                "l3_se_dist_au": round(vals[2], 9),
                "l3_se_lon_speed_deg_day": round(vals[3], 9),
                "l1b_lon_deg": round(l1b_pt["lon_deg"], 9),
                "diff_lon_arcsec": round(diff_arcsec, 4),
            }
        l3_nodes_liliths.append(entry)

    rms_intp_arcmin = math.sqrt(sum(d * d for d in intp_diffs_arcmin) / len(intp_diffs_arcmin))
    max_intp_arcmin = max(abs(d) for d in intp_diffs_arcmin)

    # 3. 12 Ayanamsas (Mean & True) vs L1B
    l3_ayanamsas = []
    for rec in l1b_ayanamsas["vectors"]:
        code = rec["ayanamsa_code"]
        legacy_id = rec["legacy_id"]
        jd_tt = rec["jd_tt"]
        swe.set_sid_mode(legacy_id, 0, 0)
        se_mean = swe.get_ayanamsa_ex(jd_tt, swe.FLG_SWIEPH | swe.FLG_NONUT)[1]
        se_true = swe.get_ayanamsa_ex(jd_tt, swe.FLG_SWIEPH)[1]
        diff_mean_arcsec = wrap180(rec["ayanamsa_mean_deg"] - se_mean) * 3600.0
        diff_true_arcsec = wrap180(rec["ayanamsa_true_deg"] - se_true) * 3600.0
        l3_ayanamsas.append({
            "ayanamsa_code": code,
            "legacy_id": legacy_id,
            "epoch_label": rec["epoch_label"],
            "jd_tt": jd_tt,
            "l1b_mean_deg": round(rec["ayanamsa_mean_deg"], 9),
            "l3_se_mean_deg": round(se_mean, 9),
            "diff_mean_arcsec": round(diff_mean_arcsec, 4),
            "l1b_true_deg": round(rec["ayanamsa_true_deg"], 9),
            "l3_se_true_deg": round(se_true, 9),
            "diff_true_arcsec": round(diff_true_arcsec, 4),
        })

    # 4. 8 Hamburg Uranians + Transpluto
    uranian_map = [
        ("Cupido", swe.CUPIDO),
        ("Hades", swe.HADES),
        ("Zeus", swe.ZEUS),
        ("Kronos", swe.KRONOS),
        ("Apollon", swe.APOLLON),
        ("Admetos", swe.ADMETOS),
        ("Vulkanus", swe.VULKANUS),
        ("Poseidon", swe.POSEIDON),
        ("Transpluto", swe.ISIS),
    ]
    l3_uranians = []
    for ep_label, jd_tt in epochs:
        for u_name, u_id in uranian_map:
            geo_vals, _ = swe.calc(jd_tt, u_id, swe.FLG_SWIEPH | swe.FLG_SPEED)
            helio_vals, _ = swe.calc(jd_tt, u_id, swe.FLG_SWIEPH | swe.FLG_HELCTR | swe.FLG_SPEED)
            l3_uranians.append({
                "epoch_label": ep_label,
                "jd_tt": jd_tt,
                "point": u_name,
                "geo_lon_deg": round(geo_vals[0], 8),
                "helio_lon_deg": round(helio_vals[0], 8),
                "geo_speed_deg_day": round(geo_vals[3], 8),
            })

    # 5. All 23 House Systems + 8 Chart Angles Verification vs L1B
    l3_houses_comparison = []
    max_midlat_house_diff_arcsec = 0.0
    max_midlat_angle_diff_arcsec = 0.0

    angle_key_order = [
        "ASC",
        "MC",
        "ARMC",
        "Vertex",
        "EqASC",
        "CoASC_Koch",
        "CoASC_Munkasey",
        "PolarASC",
    ]

    for case in l1b_houses["cases"]:
        if abs(case["lat_deg"]) >= 66.5 or abs(case["lat_deg"]) < 1e-6:
            continue
        ramc = case["ramc_deg"]
        lat = case["lat_deg"]
        eps = case["eps_true_deg"]

        # Verify 8 angles using P
        _, ascmc = swe.houses_armc(ramc, lat, eps, b"P")
        angle_diffs = {}
        for idx, k in enumerate(angle_key_order):
            d_arcsec = abs(wrap180(case["angles"][k] - ascmc[idx])) * 3600.0
            angle_diffs[k] = round(d_arcsec, 6)
            max_midlat_angle_diff_arcsec = max(max_midlat_angle_diff_arcsec, d_arcsec)

        sys_comparison = {}
        for code_str, sys_data in case["systems"].items():
            bcode = code_str.encode("ascii")
            if code_str in ("I", "i"):
                cusps_se, _ = swe.houses_ex(case["jd_ut1"], lat, case["lon_deg"], bcode, 0)
            else:
                cusps_se, _ = swe.houses_armc(ramc, lat, eps, bcode)
            l1b_cusps = sys_data["cusps"]
            diffs_arcsec = [
                abs(wrap180(l1b_cusps[i] - cusps_se[i])) * 3600.0
                for i in range(12)
            ]
            max_d = max(diffs_arcsec)
            max_midlat_house_diff_arcsec = max(max_midlat_house_diff_arcsec, max_d)
            sys_comparison[code_str] = {
                "l3_se_cusps": [round(c, 9) for c in cusps_se[:12]],
                "max_diff_vs_l1b_arcsec": round(max_d, 6),
            }

        l3_houses_comparison.append({
            "case_id": case["case_id"],
            "lat_deg": lat,
            "lon_deg": case["lon_deg"],
            "l3_se_angles_8": {k: round(ascmc[idx], 9) for idx, k in enumerate(angle_key_order)},
            "angle_diffs_arcsec": angle_diffs,
            "systems": sys_comparison,
        })

    # 6. Solar & Lunar Eclipses Verification
    l3_solar_eclipses = []
    for sol in l1b_eclipses["solar_eclipses"]:
        jd_start = sol["jd_tt_greatest"] - 5.0
        ret_flag, tret = swe.sol_eclipse_when_glob(jd_start, swe.FLG_SWIEPH, 0, False)
        dt_ut = swe.deltat(tret[0])
        jd_tt_se = tret[0] + dt_ut
        diff_sec = (sol["jd_tt_greatest"] - jd_tt_se) * 86400.0
        l3_solar_eclipses.append({
            "label": sol["label"],
            "kind": sol["kind"],
            "l1b_jd_tt_greatest": sol["jd_tt_greatest"],
            "l3_se_jd_tt_greatest": round(jd_tt_se, 7),
            "diff_seconds": round(diff_sec, 2),
            "se_ret_flag": int(ret_flag),
        })

    l3_lunar_eclipses = []
    for lun in l1b_eclipses["lunar_eclipses"]:
        jd_start = lun["jd_tt_greatest"] - 5.0
        ret_flag, tret = swe.lun_eclipse_when(jd_start, swe.FLG_SWIEPH, 0, False)
        dt_ut = swe.deltat(tret[0])
        jd_tt_se = tret[0] + dt_ut
        diff_sec = (lun["jd_tt_greatest"] - jd_tt_se) * 86400.0
        l3_lunar_eclipses.append({
            "label": lun["label"],
            "kind": lun["kind"],
            "l1b_jd_tt_greatest": lun["jd_tt_greatest"],
            "l3_se_jd_tt_greatest": round(jd_tt_se, 7),
            "diff_seconds": round(diff_sec, 2),
            "se_ret_flag": int(ret_flag),
        })

    out = {
        "schema_version": "1.0.0",
        "tier": "L3_blackbox_se_reference",
        "oracle": "pyswisseph (strictly black-box CLI/API invocation per ARCHITECTURE.md §1.2.2b)",
        "summary_metrics": {
            "max_midlat_23_houses_diff_arcsec": round(max_midlat_house_diff_arcsec, 6),
            "max_midlat_8_angles_diff_arcsec": round(max_midlat_angle_diff_arcsec, 6),
            "interpolated_lilith_vs_se_intp_apog_rms_arcmin": round(rms_intp_arcmin, 2),
            "interpolated_lilith_vs_se_intp_apog_max_arcmin": round(max_intp_arcmin, 2),
        },
        "major_bodies": major_body_vectors,
        "nodes_and_liliths": l3_nodes_liliths,
        "ayanamsas_12": l3_ayanamsas,
        "uranians_and_transpluto": l3_uranians,
        "houses_23_comparison": l3_houses_comparison,
        "solar_eclipses_comparison": l3_solar_eclipses,
        "lunar_eclipses_comparison": l3_lunar_eclipses,
    }

    out_path = GOLDEN_DIR / "l3_se_reference.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Wrote {out_path}")
    print(f"  Max mid-latitude 23-houses diff vs L1B: {max_midlat_house_diff_arcsec:.6f} arcsec")
    print(f"  Max mid-latitude 8-angles diff vs L1B:  {max_midlat_angle_diff_arcsec:.6f} arcsec")
    print(f"  Interpolated_Lilith C2 spline vs SE_INTP_APOG: RMS = {rms_intp_arcmin:.2f}', max = {max_intp_arcmin:.2f}'")


if __name__ == "__main__":
    main()
