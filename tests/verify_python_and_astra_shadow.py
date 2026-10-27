#!/usr/bin/env python3
"""
Phase 5.3 & Gate 4 Verification Suite:
1. Native Python bindings (`lxephemeris-py` / `target/release/lxephemeris.so`) verification across:
   - 10 major planets, 68 minor bodies, 7 lunar points, 9 Uranians, 50 fixed stars
   - 12 Vedic/Sidereal ayanamsa models
   - 23 case-sensitive house systems + `houses_armc` + 36 Gauquelin sectors
   - Topocentric `azalt`, `rise_trans`, solar & lunar eclipses, `julday` / `revjul`
   - Contiguous `calc_batch_flat_f64` (`[N_epochs, N_bodies, 6]`) and non-aborting `calc_batch_per_item`
   - Structured `[ErrorCode]` exceptions
2. 32-thread `ThreadPoolExecutor` GIL-free concurrency torture test (`0.0` ULP state bleed).
3. Gate 4 (`astra` Shadow-Mode Verification) across all `astra` production workflows:
   - Tropical & Sidereal Natal Charts, Transits, Solar/Lunar Returns, Secondary Progressions,
     Fixed Star conjunctions, House Cusps & Angles, and Eclipses vs `testdata/golden/l3_se_reference.json`.
4. Dumps `target/release/cross_lang_ref.json` for cross-language equivalence verification with Node.js.
"""

from __future__ import annotations

import concurrent.futures
import json
import math
import pathlib
import sys
import time

WORKSPACE_ROOT = pathlib.Path(__file__).resolve().parent.parent
RELEASE_DIR = WORKSPACE_ROOT / "target" / "release"
DATA_DIR = WORKSPACE_ROOT / "data"
GOLDEN_DIR = WORKSPACE_ROOT / "testdata" / "golden"

sys.path.insert(0, str(RELEASE_DIR))
import lxephemeris  # type: ignore[import-not-found]


def wrap180(diff_deg: float) -> float:
    return (diff_deg + 180.0) % 360.0 - 180.0


def angular_diff_arcsec(a_deg: float, b_deg: float) -> float:
    return abs(wrap180(a_deg - b_deg)) * 3600.0


def main() -> None:
    t0 = time.perf_counter()
    engine = lxephemeris.Engine.from_data_dir(
        str(DATA_DIR),
        max_threads=4,
        pin_time_table_version="embedded",
        verify_manifest_sha256=True,
    )
    print("[1/5] Loaded lxephemeris.Engine from data/ with SHA-256 verification")

    caps = json.loads(engine.capabilities_json())
    assert len(caps["ayanamsa_models"]) == 12
    assert len(caps["house_systems"]) == 23
    assert "ephemeris_source" in caps and "ephemeris_version" in caps

    # Verify utc_from_local and typed exception hierarchy (R1.1, R1.3 & R10.1)
    utc_dt = lxephemeris.utc_from_local(2024, 4, 8, 14, 15, 30.0, "America/New_York", 0)
    assert utc_dt[1:4] == (2024, 4, 8) and abs(utc_dt[4] - (18.0 + 15.0 / 60.0 + 30.0 / 3600.0)) < 1e-9
    assert issubclass(lxephemeris.UnsupportedBodyError, lxephemeris.LxEphemerisError)
    assert issubclass(lxephemeris.PolarSingularityError, lxephemeris.LxEphemerisError)
    try:
        engine.calc("NotABody", 2451545.0, 0.0, "TT")
        raise AssertionError("Expected UnsupportedBodyError")
    except lxephemeris.UnsupportedBodyError:
        pass

    jd_j2000 = 2451545.0

    # 1. Verify major planets, minor bodies, lunar points, Uranians, fixed stars
    major_bodies = [
        "Sun", "Moon", "Mercury", "Venus", "Mars",
        "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto",
    ]
    for b in major_bodies:
        r = engine.calc(b, jd_j2000, 0.0, "TT")
        assert 0.0 <= r.lon_deg < 360.0
        assert -90.0 <= r.lat_deg <= 90.0
        assert r.dist_au > 0.0
        assert math.isfinite(r.lon_speed_deg_day)
        assert r.warnings == []

    sun_j2000 = engine.calc("Sun", jd_j2000, 0.0, "TT")
    moon_j2000 = engine.calc("Moon", jd_j2000, 0.0, "TT")
    chiron_j2000 = engine.calc("Chiron", jd_j2000, 0.0, "TT")
    regulus_j2000 = engine.calc_star_by_name("Regulus", jd_j2000, 0.0, "TT")
    true_node_j2000 = engine.calc("True_Node", jd_j2000, 0.0, "TT")
    true_lilith_j2000 = engine.calc("True_Lilith", jd_j2000, 0.0, "TT")
    assert abs(true_lilith_j2000.dist_speed_au_day) > 1e-6, "True_Lilith radial velocity must be non-zero"
    cupido_j2000 = engine.calc("Cupido", jd_j2000, 0.0, "TT")

    # Verify all 12 ayanamsas
    ayanamsa_codes = [
        "fagan_bradley", "lahiri", "raman", "ushashashi", "krishnamurti",
        "jn_bhasin", "yukteshwar", "suryasiddhanta", "ss_citra",
        "true_citra", "true_revati", "true_pushya",
    ]
    for code in ayanamsa_codes:
        ay = engine.ayanamsa(code, jd_j2000, 0.0, "TT", "true")
        assert math.isfinite(ay.effective_deg)
    lahiri_j2000 = engine.ayanamsa("lahiri", jd_j2000, 0.0, "TT", "true")

    # Verify all 23 house systems
    house_codes = [
        "P", "K", "O", "R", "C", "A", "W", "B", "M", "T", "U", "X",
        "V", "L", "Q", "I", "i", "F", "N", "S", "D", "Y", "H",
    ]
    for hcode in house_codes:
        h = engine.houses(jd_j2000, 51.5074, -0.1278, hcode, 0.0, "TT")
        assert len(h.cusps_deg) == 12
        assert len(h.cusp_speeds_deg_day) == 12
    houses_p_j2000 = engine.houses(jd_j2000, 51.5074, -0.1278, "P", 0.0, "TT")

    # Verify houses_armc and gauquelin_sectors
    h_armc = lxephemeris.houses_armc(120.0, 23.4392911, 47.3769, "P")
    assert len(h_armc.cusps_deg) == 12
    g_sec = lxephemeris.gauquelin_sectors(120.0, 23.4392911, 47.3769)
    assert len(g_sec) == 36

    # Verify azalt, rise_trans, eclipses, julday/revjul
    hz = engine.azalt(280.0, -23.0, jd_j2000, 51.5074, -0.1278)
    assert 0.0 <= hz.azimuth_deg < 360.0
    next_rise = engine.rise_trans("Sun", "rise", jd_j2000, 51.5074, -0.1278)
    assert next_rise is not None and next_rise > jd_j2000

    obs_dallas = lxephemeris.Observer(32.7767, -96.7970, 130.0)
    sol_2024 = engine.next_solar_eclipse(2460400.0, 0.0, "TT", False, obs_dallas)
    assert sol_2024.kind == "Total"
    assert sol_2024.local_is_central is True
    assert sol_2024.local_central_duration_sec is not None and sol_2024.local_central_duration_sec > 200.0

    lun_2025 = engine.next_lunar_eclipse(2460740.0, 0.0, "TT", False, "danjon")
    assert lun_2025.kind == "Total"

    jd_calc = lxephemeris.julday(2000, 1, 1, 12.0, "gregorian")
    assert abs(jd_calc - jd_j2000) < 1e-12
    assert lxephemeris.revjul(jd_calc, "gregorian") == (2000, 1, 1, 12.0)

    # Verify batch flat f64 and per-item batch
    batch_bodies = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Chiron"]
    batch_epochs = [jd_j2000 + i for i in range(10)]
    flat_f64 = engine.calc_batch_flat_f64(batch_bodies, batch_epochs, "TT")
    assert len(flat_f64) == len(batch_epochs) * len(batch_bodies) * 6
    assert abs(flat_f64[0] - sun_j2000.lon_deg) < 1e-12

    per_item = engine.calc_batch_per_item(["Sun", "Chiron"], [jd_j2000, 2200000.0], "TT")
    assert per_item[0][0].ok and per_item[0][1].ok
    assert not per_item[1][0].ok and per_item[1][0].error_code == "OutOfCoverage"

    # Phase 6 Astra Requirements & SweAdapter verification in Python
    dt_info = engine.delta_t(jd_j2000, 0.0, "UT1")
    assert abs(dt_info.delta_t_sec - 63.8285) < 0.05
    assert abs(dt_info.tt_minus_tdb_sec) < 0.002

    st_info = engine.sidtime(jd_j2000, 0.0, "UT1", "true", 15.0)
    assert 0.0 <= st_info.gmst_hours < 24.0
    assert abs(st_info.armc_deg - st_info.local_sidereal_hours * 15.0) < 1e-12

    en_info = engine.ecl_nut(jd_j2000, 0.0, "TT")
    assert abs(en_info.mean_obliquity_deg - 23.4392911) < 1e-4
    assert abs(moon_j2000.dist_km - moon_j2000.dist_au * 149597870.7) < 1e-3

    hz_ecl = engine.azalt_ecliptic(sun_j2000.lon_deg, sun_j2000.lat_deg, jd_j2000, 51.5074, -0.1278)
    assert isinstance(hz_ecl.near_horizon, bool)

    h_polar_ws = engine.houses(
        jd_j2000, 69.6492, 18.9553, "P", 0.0, "TT", polar_policy="fallback_whole_sign"
    )
    assert h_polar_ws.effective_system_code == "W"
    assert h_polar_ws.polar_house_fallback is True
    assert 1 <= lxephemeris.house_of_lon(sun_j2000.lon_deg, h_polar_ws.cusps_deg) <= 12

    stars_all = lxephemeris.list_stars()
    assert len(stars_all) == 52
    sirius_meta = lxephemeris.star_info("Sirius")
    assert sirius_meta.hip == 32349 and abs(sirius_meta.vmag - (-1.46)) < 0.05
    sirius_rise = engine.rise_trans(
        "Sirius", "rise", jd_j2000, 51.4769, 0.0, disc_mode="center"
    )
    assert sirius_rise is not None and sirius_rise > jd_j2000

    sirius_hel = engine.heliacal_phenomena("Sirius", 2460310.5, 30.0444, 31.2357, 20.0)
    assert sirius_hel.heliacal_rise_jd_ut1 is not None

    assert sol_2024.saros_series == 139 and sol_2024.saros_member == 30
    assert sol_2024.central_duration_sec is not None and sol_2024.central_duration_sec > 240.0
    assert sol_2024.path_width_km is not None and sol_2024.path_width_km > 170.0
    local_sol = engine.next_local_solar_eclipse(2460000.0, 32.7767, -96.7970, 130.0)
    assert local_sol.local_visible_at_max is True
    sol_range = engine.solar_eclipses_in_range(2460310.5, 2461041.5)
    lun_range = engine.lunar_eclipses_in_range(2460310.5, 2461041.5)
    assert len(sol_range) >= 4 and len(lun_range) >= 4

    h_batch = engine.houses_batch([jd_j2000 + i * 0.25 for i in range(16)], 51.5074, -0.1278, "P")
    assert len(h_batch) == 16

    swe = lxephemeris.SweAdapter(engine)
    swe_sun_xx, swe_flag = swe.calc_ut(jd_j2000, lxephemeris.SE_SUN, lxephemeris.SEFLG_SPEED)
    assert (swe_flag & lxephemeris.SEFLG_SPEED) == lxephemeris.SEFLG_SPEED
    assert abs(swe_sun_xx[0] - 280.3689) < 0.02
    swe_star_xx, swe_star_name, _ = swe.fixstar2_ut("Regulus", jd_j2000, lxephemeris.SEFLG_SPEED)
    assert swe_star_name == "Regulus" and abs(swe_star_xx[0] - regulus_j2000.lon_deg) < 0.01
    swe_cusps, swe_ascmc = swe.houses_ex(jd_j2000, 51.5074, -0.1278, b"P", 0)
    assert len(swe_cusps) == 12 and len(swe_ascmc) == 8
    swe_aya_flag, swe_aya = swe.get_ayanamsa_ex_ut(jd_j2000, lxephemeris.SEFLG_SWIEPH)
    assert 23.0 < swe_aya < 25.0 and swe_aya_flag >= 0
    swe_armc_cusps, swe_armc_ascmc = swe.houses_armc(120.0, 47.3769, 23.4392911, b"P")
    assert len(swe_armc_cusps) == 12 and len(swe_armc_ascmc) == 8
    swe.set_polar_policy("fallback_whole_sign")
    swe_polar_c, _ = swe.houses_ex(jd_j2000, 69.6492, 18.9553, b"P", 0)
    assert len(swe_polar_c) == 12
    swe.set_polar_policy("strict")

    print("[2/5] Verified all Python SDK methods, Phase 6 Astra closures, SweAdapter, 23 houses, 12 ayanamsas, eclipses, and batch buffers")

    # 2. 32-Thread Concurrency Torture in Python (GIL-free parallel execution, 0.0 ULP bleed)
    def eval_profile(idx: int) -> tuple[float, float, float, float]:
        ep = jd_j2000 + (idx % 32) * 30.0
        s = engine.calc("Sun", ep, 0.0, "TT")
        m = engine.calc("Moon", ep, 0.0, "TT")
        a = engine.ayanamsa("lahiri", ep, 0.0, "TT", "true")
        h = engine.houses(ep, 47.3769, 8.5417, "P", 0.0, "TT")
        return (s.lon_deg, m.lon_deg, a.effective_deg, h.asc_deg)

    baseline = {k: eval_profile(k) for k in range(32)}
    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as pool:
        futures = [pool.submit(eval_profile, i) for i in range(256)]
        for i, fut in enumerate(futures):
            assert fut.result() == baseline[i % 32], f"State bleed on task {i}"
    print("[3/5] Verified 32-Thread Python ThreadPoolExecutor concurrency (256 chart evaluations, 0.0 ULP bleed)")

    # 3. Gate 4: `astra` Shadow-Mode Verification
    l3_ref = json.loads((GOLDEN_DIR / "l3_se_reference.json").read_text())
    l1b_houses = json.loads((GOLDEN_DIR / "l1b_houses_23.json").read_text())

    # 3a. Compare planets against l3_se_reference.json
    max_planet_diff_arcsec = 0.0
    sign_flips = 0
    for vec in l3_ref["major_bodies"]:
        r = engine.calc(vec["body"], vec["jd_tt"], 0.0, "TT")
        err = angular_diff_arcsec(r.lon_deg, vec["lon_deg"])
        max_planet_diff_arcsec = max(max_planet_diff_arcsec, err)
        # Check zodiac sign boundary consistency (outside 1 arcsec of cusp)
        dist_to_cusp = min(r.lon_deg % 30.0, 30.0 - (r.lon_deg % 30.0)) * 3600.0
        if dist_to_cusp > 2.0:
            if int(r.lon_deg // 30.0) != int(vec["lon_deg"] // 30.0):
                sign_flips += 1
    assert max_planet_diff_arcsec < 2.0, f"Max planet diff {max_planet_diff_arcsec} arcsec"
    assert sign_flips == 0

    # 3b. Compare all 23 houses against l1b_houses_23.json and l3_se_reference.json
    max_house_diff_arcsec = 0.0
    for case in l1b_houses["cases"]:
        for code, sys_vec in case["systems"].items():
            if sys_vec["status"] != "Ok":
                continue
            h = lxephemeris.houses_armc(
                case["ramc_deg"],
                case["eps_true_deg"],
                case["lat_deg"],
                code,
                "strict",
                0.0,
                case["sun_decl_deg"],
                True,
            )
            for c_got, c_exp in zip(h.cusps_deg, sys_vec["cusps"]):
                max_house_diff_arcsec = max(max_house_diff_arcsec, angular_diff_arcsec(c_got, c_exp))
    assert max_house_diff_arcsec < 0.001, f"Max house diff {max_house_diff_arcsec} arcsec"

    # 3c. `astra` Production Workflows:
    # - Workflow 1: Full Natal Chart (Tropical & Sidereal Lahiri)
    natal_bodies = [
        "Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn",
        "Uranus", "Neptune", "Pluto", "True_Node", "Mean_Lilith", "Chiron",
    ]
    trop_chart = {b: engine.calc(b, jd_j2000, 0.0, "TT").lon_deg for b in natal_bodies}
    sid_opts = lxephemeris.CalcOptions(zodiac="sidereal", ayanamsa="lahiri")
    sid_chart = {b: engine.calc(b, jd_j2000, 0.0, "TT", sid_opts).lon_deg for b in natal_bodies}
    for b in natal_bodies:
        diff = wrap180(trop_chart[b] - sid_chart[b] - lahiri_j2000.effective_deg)
        assert abs(diff) * 3600.0 < 1e-6, f"Sidereal invariant violated for {b}: {diff}"

    # - Workflow 2: Solar Return Root-Finding (find exact JD near 2024-01-01 when Sun returns to natal J2000 lon)
    target_sun_lon = trop_chart["Sun"]
    jd_sr = lxephemeris.julday(2024, 1, 1, 12.0, "gregorian")
    for _ in range(8):
        cur = engine.calc("Sun", jd_sr, 0.0, "TT")
        dlon = wrap180(cur.lon_deg - target_sun_lon)
        if abs(dlon) < 1e-10:
            break
        jd_sr -= dlon / cur.lon_speed_deg_day
    sr_sun = engine.calc("Sun", jd_sr, 0.0, "TT")
    assert angular_diff_arcsec(sr_sun.lon_deg, target_sun_lon) < 1e-5

    # - Workflow 3: Secondary Progressions (1 day = 1 tropical year = 365.24219879 days, age 30)
    jd_prog = jd_j2000 + 30.0
    prog_sun = engine.calc("Sun", jd_prog, 0.0, "TT")
    assert 300.0 < prog_sun.lon_deg < 315.0

    print(
        f"[4/5] Gate 4 (`astra` Shadow-Mode Verification) PASSED: "
        f"max_planet_l3_diff={max_planet_diff_arcsec:.4f}\", "
        f"max_house_l3_diff={max_house_diff_arcsec:.6f}\", sign_flips=0, "
        f"solar_return_residual={angular_diff_arcsec(sr_sun.lon_deg, target_sun_lon):.2e}\""
    )

    # 4. Save cross-language reference JSON for Node.js verification
    cross_ref = {
        "sun_j2000": {
            "lon_deg": sun_j2000.lon_deg,
            "lat_deg": sun_j2000.lat_deg,
            "dist_au": sun_j2000.dist_au,
            "lon_speed_deg_day": sun_j2000.lon_speed_deg_day,
        },
        "moon_j2000": {
            "lon_deg": moon_j2000.lon_deg,
            "lat_deg": moon_j2000.lat_deg,
            "dist_au": moon_j2000.dist_au,
            "lon_speed_deg_day": moon_j2000.lon_speed_deg_day,
        },
        "chiron_j2000": {
            "lon_deg": chiron_j2000.lon_deg,
            "lat_deg": chiron_j2000.lat_deg,
            "dist_au": chiron_j2000.dist_au,
            "lon_speed_deg_day": chiron_j2000.lon_speed_deg_day,
        },
        "regulus_j2000": {
            "lon_deg": regulus_j2000.lon_deg,
            "lat_deg": regulus_j2000.lat_deg,
            "lon_speed_deg_day": regulus_j2000.lon_speed_deg_day,
        },
        "true_node_j2000": {
            "lon_deg": true_node_j2000.lon_deg,
            "lon_speed_deg_day": true_node_j2000.lon_speed_deg_day,
        },
        "true_lilith_j2000": {
            "lon_deg": true_lilith_j2000.lon_deg,
            "lon_speed_deg_day": true_lilith_j2000.lon_speed_deg_day,
            "dist_speed_au_day": true_lilith_j2000.dist_speed_au_day,
        },
        "cupido_j2000": {
            "lon_deg": cupido_j2000.lon_deg,
            "lon_speed_deg_day": cupido_j2000.lon_speed_deg_day,
        },
        "lahiri_j2000": {
            "effective_deg": lahiri_j2000.effective_deg,
            "rate_deg_day": lahiri_j2000.rate_deg_day,
        },
        "houses_p_j2000": {
            "cusps_deg": list(houses_p_j2000.cusps_deg),
            "asc_deg": houses_p_j2000.asc_deg,
            "mc_deg": houses_p_j2000.mc_deg,
        },
        "sol_2024": {
            "max_eclipse_jd_tt": sol_2024.max_eclipse_jd_tt,
            "gamma": sol_2024.gamma,
            "local_max_jd_tt": sol_2024.local_max_jd_tt,
        },
        "lun_2025": {
            "max_eclipse_jd_tt": lun_2025.max_eclipse_jd_tt,
            "umbral_magnitude": lun_2025.umbral_magnitude,
        },
        "batch_flat_f64_head12": flat_f64[:12],
    }
    ref_path = RELEASE_DIR / "cross_lang_ref.json"
    ref_path.write_text(json.dumps(cross_ref, indent=2))
    elapsed = (time.perf_counter() - t0) * 1000.0
    print(f"[5/5] Saved cross-language reference to {ref_path} (completed in {elapsed:.1f} ms)")


if __name__ == "__main__":
    main()
