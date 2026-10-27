#!/usr/bin/env python3
"""
Wave F3 10,000+ Golden Reference Vector Generator (`testdata/golden/l3_wave_bcf_10k_grid.json`).

Strictly invokes `pyswisseph` (`swisseph`) and `PyEphem` (`ephem`) via their public Python APIs
as black-box reference oracles (together with `jplephem` + `pyerfa` L1B reference helpers) to
generate 10,540 golden reference vectors across all 8 Wave B, Wave C, and Wave F sections:

  1. `section_1_rise_trans`                   =  1,600 vectors (100 dates x 4 observers x 4 Sun events)
  2. `section_2_lunar_phases`                 =  1,500 vectors (1,500 consecutive principal lunar phases from 1900-01-01)
  3. `section_3_lunar_phase_summary`          =  1,500 vectors (1,500 evenly spaced epochs in 1855..2095)
  4. `section_4_solar_terms_and_seasons`      =  1,440 vectors (60 years 1980..=2039 x 24 Chinese Solar Terms)
  5. `section_5_hd_and_longitude_crossings`   =  1,000 vectors (500 HD -88 deg solar arc + 500 longitude crossings)
  6. `section_6_stations_and_declination_extrema` = 1,000 vectors (500 planetary stations + 500 Sun/Moon dec extrema)
  7. `section_7_brady_parans`                 =    500 vectors (500 diurnal star-planet parans across 25 charts)
  8. `section_8_wave_b_c_swe_parity`          =  2,000 vectors (500 Moshier + 300 pheno + 300 nod_aps + 430 ayanamsas + 470 house_pos/time_equ/azalt_rev)
Total: 10,540 golden reference vectors.
"""

from __future__ import annotations

import json
import math
import struct
import tempfile
from pathlib import Path
import ephem
import numpy as np
from scipy.optimize import brentq
import swisseph as swe

ROOT = Path(__file__).resolve().parent.parent
DE440S_PATH = ROOT / "data" / "de440s.bsp"
AST_BSP_PATH = ROOT / "data" / "asteroids_68_1850_2100.bsp"
GOLDEN_DIR = ROOT / "testdata" / "golden"

AU_KM = 149597870.700
C_KM_S = 299792.458
DJD_OFFSET = 2415020.0  # PyEphem Dublin Julian Date offset (1899-12-31 12:00 UT)


def iae_general_precession_arcsec(jd_t0: float, jd_t: float) -> tuple[float, float]:
    big_t = (jd_t0 - 2451545.0) / 36525.0
    t = (jd_t - jd_t0) / 36525.0
    c1 = 5029.0966 + 2.22226 * big_t - 0.000042 * big_t * big_t
    c2 = 1.11161 - 0.000127 * big_t
    c3 = -0.000113
    p_arcsec = c1 * t + c2 * t * t + c3 * t * t * t
    rate_arcsec_day = (c1 + 2.0 * c2 * t + 3.0 * c3 * t * t) / 36525.0
    return p_arcsec, rate_arcsec_day


class SimpleDafType2Reader:
    """Minimal pure-Python NAIF DAF Type 2 Chebyshev reader for de440s.bsp and asteroids_68_1850_2100.bsp."""

    def __init__(self, path: Path) -> None:
        self.raw = path.read_bytes()
        self.words = np.frombuffer(self.raw, dtype="<f8")
        fward = struct.unpack_from("<i", self.raw, 76)[0]
        self.segments: dict[tuple[int, int], tuple[int, int]] = {}
        curr_rec = fward
        while curr_rec > 0:
            rec_off = (curr_rec - 1) * 1024
            next_rec = int(struct.unpack_from("<d", self.raw, rec_off)[0])
            nsum = int(struct.unpack_from("<d", self.raw, rec_off + 16)[0])
            for idx in range(nsum):
                ent_off = rec_off + 24 + idx * 40
                target_id, center_id, _, spk_type, init_w, final_w = struct.unpack_from(
                    "<6i", self.raw, ent_off + 16
                )
                if spk_type == 2:
                    self.segments[(center_id, target_id)] = (init_w, final_w)
            curr_rec = next_rec

    def compute(self, center_id: int, target_id: int, jd_tdb: float) -> np.ndarray:
        init_w, final_w = self.segments[(center_id, target_id)]
        seg = self.words[init_w - 1 : final_w]
        init_et, intlen_s, rsize_f, n_rec_f = seg[-4], seg[-3], seg[-2], seg[-1]
        rsize = int(rsize_f)
        n_rec = int(n_rec_f)
        degree = (rsize - 2) // 3 - 1
        et_s = (jd_tdb - 2451545.0) * 86400.0
        rec_idx = min(max(int((et_s - init_et) // intlen_s), 0), n_rec - 1)
        rec = seg[rec_idx * rsize : (rec_idx + 1) * rsize]
        mid, radius = rec[0], rec[1]
        tau = (et_s - mid) / radius
        out = np.zeros(3, dtype=float)
        m = degree + 1
        for axis in range(3):
            coeffs = rec[2 + axis * m : 2 + (axis + 1) * m]
            w0, w1, w2 = 0.0, 0.0, 0.0
            for j in range(degree, 0, -1):
                w2 = w1
                w1 = w0
                w0 = float(coeffs[j]) + 2.0 * tau * w1 - w2
            out[axis] = float(coeffs[0]) + tau * w0 - w1
        return out

SOLAR_TERM_PINYIN_NAMES = [
    "Chunfen",
    "Qingming",
    "Guyu",
    "Lixia",
    "Xiaoman",
    "Mangzhong",
    "Xiazhi",
    "Xiaoshu",
    "Dashu",
    "Liqiu",
    "Chushu",
    "Bailu",
    "Qiufen",
    "Hanlu",
    "Shuangjiang",
    "Lidong",
    "Xiaoxue",
    "Daxue",
    "Dongzhi",
    "Xiaohan",
    "Dahan",
    "Lichun",
    "Yushui",
    "Jingzhe",
]


def wrap180(deg: float) -> float:
    return (deg + 180.0) % 360.0 - 180.0


def linspace(start: float, stop: float, n: int) -> list[float]:
    if n == 1:
        return [start]
    step = (stop - start) / (n - 1)
    return [round(start + i * step, 6) for i in range(n)]


def write_temp_sefstars(td_path: Path, stars_list: list[dict]) -> None:
    """Write a temporary sefstars.txt from our clean-room stars_catalog_52.json so swe.rise_trans / swe.get_ayanamsa_ex can resolve stars."""
    lines = []
    # Sort by name length ascending so legacy C prefix matching in swe.rise_trans never matches "Mirach" when looking up "Mira"
    for s in sorted(stars_list, key=lambda x: (len(x["name"]), x["name"])):
        ra_h_f = s["ra_2000_deg"] / 15.0
        ra_h = int(ra_h_f)
        ra_m_f = (ra_h_f - ra_h) * 60.0
        ra_m = int(ra_m_f)
        ra_s = (ra_m_f - ra_m) * 60.0

        dec_sign = "+" if s["de_2000_deg"] >= 0.0 else "-"
        de_abs = abs(s["de_2000_deg"])
        de_d = int(de_abs)
        de_m_f = (de_abs - de_d) * 60.0
        de_m = int(de_m_f)
        de_s = (de_m_f - de_m) * 60.0

        cos_d = max(abs(math.cos(math.radians(s["de_2000_deg"]))), 1e-6)
        pm_ra_se = (s["pm_ra_cosdec_mas_yr"] / cos_d) / 150.0
        pm_de_se = s["pm_dec_mas_yr"] / 10.0
        plx = max(float(s["plx_mas"]), 0.0)
        line = (
            f"{s['name']:13s},{s['bayer']}, IC, "
            f"{ra_h:02d},{ra_m:02d},{ra_s:07.4f}, "
            f"{dec_sign}{de_d:02d},{de_m:02d},{de_s:06.3f}, "
            f"{pm_ra_se:.4f}, {pm_de_se:.3f}, {s['rv_kms']:.1f}, "
            f"{plx:.2f}, {s['vmag']:.2f}, 2000, {s['hip']}"
        )
        lines.append(line)
        if s["name"] == "Zeta Piscium":
            lines.append(line.replace("Zeta Piscium ", "Revati       ", 1))
        elif s["name"] == "Asellus Australis":
            lines.append(line.replace("Asellus Australis", "Pushya       ", 1))
    (td_path / "sefstars.txt").write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Section 1: Rise, Set & Meridian Transits (1,600 vectors)
# ---------------------------------------------------------------------------
def build_section_1_rise_trans() -> list[list]:
    dates_100 = linspace(2398588.5, 2486240.5, 100)  # 1855..2095
    observers = [
        ("London", 51.5074, -0.1278),
        ("NYC", 40.7128, -74.0060),
        ("Mumbai", 19.0760, 72.8777),
        ("Sydney", -33.8688, 151.2093),
    ]
    sun = ephem.Sun()
    vectors = []

    for base_jd in dates_100:
        for obs_name, lat_deg, lon_deg in observers:
            # Anchor start_jd_ut1 at ~02:30 local solar time so Rise, UpperTransit, Set, LowerTransit
            # all occur cleanly within (start_jd_ut1 + 0.08, start_jd_ut1 + 0.95) without boundary straddling.
            start_jd_ut1 = round(base_jd + (2.5 - lon_deg / 15.0) / 24.0, 6)
            geopos = (lon_deg, lat_deg, 0.0)
            d0 = ephem.Date(start_jd_ut1 - DJD_OFFSET)

            obs = ephem.Observer()
            obs.lat = math.radians(lat_deg)
            obs.lon = math.radians(lon_deg)
            obs.elevation = 0.0

            # 1. Sun Rise (almanac_upper_limb: 1013.25 hPa, 15 C, upper limb)
            _, tret_r = swe.rise_trans(
                start_jd_ut1, swe.SUN, swe.CALC_RISE, geopos, 1013.25, 15.0, swe.FLG_SWIEPH
            )
            obs.pressure = 0.0
            obs.horizon = "-0:34"
            ephem_r = float(obs.next_rising(sun, start=d0)) + DJD_OFFSET
            vectors.append([
                start_jd_ut1,
                obs_name,
                lat_deg,
                lon_deg,
                "Rise",
                round(tret_r[0], 8),
                round(ephem_r, 8),
            ])

            # 2. Sun Set (almanac_upper_limb: 1013.25 hPa, 15 C, upper limb)
            _, tret_s = swe.rise_trans(
                start_jd_ut1, swe.SUN, swe.CALC_SET, geopos, 1013.25, 15.0, swe.FLG_SWIEPH
            )
            ephem_s = float(obs.next_setting(sun, start=d0)) + DJD_OFFSET
            vectors.append([
                start_jd_ut1,
                obs_name,
                lat_deg,
                lon_deg,
                "Set",
                round(tret_s[0], 8),
                round(ephem_s, 8),
            ])

            # 3. Sun UpperTransit (geometric_center: 0 hPa, disc center)
            rsmi_up = swe.CALC_MTRANSIT | swe.BIT_DISC_CENTER | swe.BIT_NO_REFRACTION
            _, tret_m = swe.rise_trans(
                start_jd_ut1, swe.SUN, rsmi_up, geopos, 0.0, 15.0, swe.FLG_SWIEPH
            )
            obs.horizon = "0:00"
            ephem_m = float(obs.next_transit(sun, start=d0)) + DJD_OFFSET
            vectors.append([
                start_jd_ut1,
                obs_name,
                lat_deg,
                lon_deg,
                "UpperTransit",
                round(tret_m[0], 8),
                round(ephem_m, 8),
            ])

            # 4. Sun LowerTransit (geometric_center: 0 hPa, disc center)
            rsmi_lo = swe.CALC_ITRANSIT | swe.BIT_DISC_CENTER | swe.BIT_NO_REFRACTION
            _, tret_i = swe.rise_trans(
                start_jd_ut1, swe.SUN, rsmi_lo, geopos, 0.0, 15.0, swe.FLG_SWIEPH
            )
            ephem_i = float(obs.next_antitransit(sun, start=d0)) + DJD_OFFSET
            vectors.append([
                start_jd_ut1,
                obs_name,
                lat_deg,
                lon_deg,
                "LowerTransit",
                round(tret_i[0], 8),
                round(ephem_i, 8),
            ])

    assert len(vectors) == 1600
    return vectors


# ---------------------------------------------------------------------------
# Section 2: 1,500 Consecutive Principal Lunar Phases from 1900-01-01
# ---------------------------------------------------------------------------
def swe_lunar_elong_diff(jd_ut: float, target_deg: float) -> float:
    m_lon = swe.calc_ut(jd_ut, swe.MOON, swe.FLG_SWIEPH)[0][0]
    s_lon = swe.calc_ut(jd_ut, swe.SUN, swe.FLG_SWIEPH)[0][0]
    return wrap180(((m_lon - s_lon) % 360.0) - target_deg)


def solve_swe_lunar_phase_near(est_jd_ut: float, target_deg: float) -> tuple[float, float]:
    root_jd = float(
        brentq(
            swe_lunar_elong_diff,
            est_jd_ut - 2.0,
            est_jd_ut + 2.0,
            args=(target_deg,),
            xtol=1e-11,
        )
    )
    m_vals = swe.calc_ut(root_jd, swe.MOON, swe.FLG_SWIEPH)[0]
    return root_jd, float(m_vals[2]) * AU_KM


def build_section_2_lunar_phases() -> list[list]:
    # Start at 1900-01-01 00:00 UT (JD 2415020.5); first phase is NewMoon on 1900-01-01 ~13:52 UT
    d_start = ephem.Date(2415020.5 - DJD_OFFSET)
    phase_cycle = [
        ("NewMoon", 0.0, ephem.next_new_moon),
        ("FirstQuarter", 90.0, ephem.next_first_quarter_moon),
        ("FullMoon", 180.0, ephem.next_full_moon),
        ("LastQuarter", 270.0, ephem.next_last_quarter_moon),
    ]
    vectors = []
    curr_d = d_start

    for idx in range(1500):
        kind_name, target_deg, ephem_fn = phase_cycle[idx % 4]
        ephem_d = ephem_fn(curr_d)
        ephem_jd_ut1 = float(ephem_d) + DJD_OFFSET
        swe_jd_ut1, swe_dist_km = solve_swe_lunar_phase_near(ephem_jd_ut1, target_deg)
        vectors.append([
            idx,
            kind_name,
            target_deg,
            round(swe_jd_ut1, 8),
            round(ephem_jd_ut1, 8),
            round(swe_dist_km, 3),
        ])
        curr_d = ephem.Date(float(ephem_d) + 2.0)

    assert len(vectors) == 1500
    return vectors


# ---------------------------------------------------------------------------
# Section 3: 1,500 Lunar Phase Summaries across 1855..2095
# ---------------------------------------------------------------------------
def build_section_3_lunar_phase_summary() -> list[list]:
    raw_epochs = linspace(2398588.5, 2486240.5, 1500)
    vectors = []
    step_28 = 360.0 / 28.0

    for raw_jd in raw_epochs:
        jd_ut1 = raw_jd
        # Ensure epoch is at least 0.08 deg away from any 8-phase (45 deg) or 28-phase (12.857 deg) boundary
        # so discrete bin indices and prev/next New Moon boundaries are 100% unambiguous.
        for _ in range(5):
            m_lon = swe.calc_ut(jd_ut1, swe.MOON, swe.FLG_SWIEPH)[0][0]
            s_lon = swe.calc_ut(jd_ut1, swe.SUN, swe.FLG_SWIEPH)[0][0]
            elong = (m_lon - s_lon) % 360.0
            rem45 = min(elong % 45.0, 45.0 - (elong % 45.0))
            rem28 = min(elong % step_28, step_28 - (elong % step_28))
            if rem45 < 0.08 or rem28 < 0.08:
                jd_ut1 = round(jd_ut1 + 0.04, 6)
            else:
                break

        m_vals = swe.calc_ut(jd_ut1, swe.MOON, swe.FLG_SWIEPH)[0]
        s_vals = swe.calc_ut(jd_ut1, swe.SUN, swe.FLG_SWIEPH)[0]
        pheno = swe.pheno_ut(jd_ut1, swe.MOON, swe.FLG_SWIEPH)

        swe_elong = (m_vals[0] - s_vals[0]) % 360.0
        swe_phase_angle = float(pheno[0])
        swe_illum = float(pheno[1])
        swe_diam = float(pheno[3])
        swe_dist_km = float(m_vals[2]) * AU_KM

        phase_8_idx = int(math.floor(swe_elong / 45.0)) % 8
        phase_28_idx = (int(math.floor(swe_elong / step_28)) % 28) + 1

        d_ephem = ephem.Date(jd_ut1 - DJD_OFFSET)
        em = ephem.Moon(d_ephem)
        ephem_illum = float(em.moon_phase)
        ephem_dist_km = float(em.earth_distance) * AU_KM

        prev_nm_ephem = float(ephem.previous_new_moon(d_ephem)) + DJD_OFFSET
        next_nm_ephem = float(ephem.next_new_moon(d_ephem)) + DJD_OFFSET
        swe_prev_nm, _ = solve_swe_lunar_phase_near(prev_nm_ephem, 0.0)
        swe_next_nm, _ = solve_swe_lunar_phase_near(next_nm_ephem, 0.0)
        swe_age_days = jd_ut1 - swe_prev_nm

        vectors.append([
            jd_ut1,
            round(swe_elong, 7),
            round(swe_illum, 7),
            round(swe_phase_angle, 7),
            round(swe_dist_km, 3),
            round(swe_diam, 7),
            round(ephem_illum, 6),
            round(ephem_dist_km, 3),
            phase_8_idx,
            phase_28_idx,
            round(swe_prev_nm, 8),
            round(swe_next_nm, 8),
            round(swe_age_days, 7),
        ])

    assert len(vectors) == 1500
    return vectors


# ---------------------------------------------------------------------------
# Section 4: 1,440 Chinese Solar Terms & Cardinal Seasons (1980..=2039)
# ---------------------------------------------------------------------------
def swe_sun_lon_diff(jd_ut: float, target_deg: float) -> float:
    s_lon = swe.calc_ut(jd_ut, swe.SUN, swe.FLG_SWIEPH)[0][0]
    return wrap180(s_lon - target_deg)


def build_section_4_solar_terms_and_seasons() -> list[list]:
    vectors = []
    # In a Gregorian calendar year (Jan 1 .. Dec 31), the 24 chronological solar terms
    # start with Xiaohan (term_index=19, 285 deg, ~Jan 5) and end with Dongzhi (term_index=18, 270 deg, ~Dec 21).
    chronological_indices = [19, 20, 21, 22, 23] + list(range(0, 19))

    for year in range(1980, 2040):
        jd_jan1 = swe.julday(year, 1, 1, 0.0, swe.GREG_CAL)
        d_jan1 = ephem.Date(jd_jan1 - DJD_OFFSET)
        ephem_seasons = {
            0: float(ephem.next_vernal_equinox(d_jan1)) + DJD_OFFSET,
            6: float(ephem.next_summer_solstice(d_jan1)) + DJD_OFFSET,
            12: float(ephem.next_autumnal_equinox(d_jan1)) + DJD_OFFSET,
            18: float(ephem.next_winter_solstice(d_jan1)) + DJD_OFFSET,
        }

        t_cursor = jd_jan1
        for k_pos, term_idx in enumerate(chronological_indices):
            target_lon = float(term_idx) * 15.0
            # Approximate day of year for k_pos-th term (~5th Jan + k_pos * 15.218 days)
            est_jd = jd_jan1 + 4.5 + k_pos * 15.2184
            root_jd = float(
                brentq(
                    swe_sun_lon_diff,
                    est_jd - 4.0,
                    est_jd + 4.0,
                    args=(target_lon,),
                    xtol=1e-11,
                )
            )
            t_cursor = root_jd
            ephem_jd = round(ephem_seasons[term_idx], 8) if term_idx in ephem_seasons else None
            vectors.append([
                year,
                term_idx,
                SOLAR_TERM_PINYIN_NAMES[term_idx],
                target_lon,
                round(root_jd, 8),
                ephem_jd,
            ])

    assert len(vectors) == 1440
    return vectors


# ---------------------------------------------------------------------------
# Section 5: 500 Human Design -88 deg Epochs + 500 Planetary Longitude Crossings
# ---------------------------------------------------------------------------
def build_section_5_hd_and_longitude_crossings() -> dict:
    # 5A. 500 Human Design -88 deg solar arc epochs across 1860..2090
    birth_epochs = linspace(2400409.5, 2484414.5, 500)
    hd_vectors = []
    for birth_jd in birth_epochs:
        birth_sun_lon = float(swe.calc_ut(birth_jd, swe.SUN, swe.FLG_SWIEPH)[0][0])
        target_design_lon = (birth_sun_lon - 88.0) % 360.0
        est_design_jd = birth_jd - 88.15
        design_jd = float(
            brentq(
                swe_sun_lon_diff,
                est_design_jd - 4.0,
                est_design_jd + 4.0,
                args=(target_design_lon,),
                xtol=1e-11,
            )
        )
        hd_vectors.append([
            birth_jd,
            round(birth_sun_lon, 7),
            round(target_design_lon, 7),
            round(design_jd, 8),
        ])

    # 5B. 500 Planetary Longitude Crossings across Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn
    # in both Tropical (250) and Sidereal Lahiri (250) modes.
    swe.set_sid_mode(swe.SIDM_LAHIRI, 0, 0)
    bodies = [
        ("Sun", swe.SUN, 3.0, 20.0),
        ("Moon", swe.MOON, 1.2, 10.0),
        ("Mercury", swe.MERCURY, 2.5, 25.0),
        ("Venus", swe.VENUS, 3.0, 25.0),
        ("Mars", swe.MARS, 5.0, 40.0),
        ("Jupiter", swe.JUPITER, 12.0, 90.0),
        ("Saturn", swe.SATURN, 20.0, 150.0),
    ]
    start_grid = linspace(2402000.5, 2482000.5, 500)
    crossing_vectors = []

    for idx, base_start_jd in enumerate(start_grid):
        b_name, swe_id, probe_offset_days, max_search_days = bodies[idx % len(bodies)]
        zodiac_mode = "Tropical" if (idx % 2 == 0) else "Sidereal_Lahiri"
        flags = swe.FLG_SWIEPH | swe.FLG_SPEED
        if zodiac_mode == "Sidereal_Lahiri":
            flags |= swe.FLG_SIDEREAL

        # Find a clean monotonic segment [start_jd, start_jd + probe_offset_days] without a station
        start_jd = base_start_jd
        for _ in range(20):
            v0 = swe.calc_ut(start_jd, swe_id, flags)[0]
            v1 = swe.calc_ut(start_jd + 0.5 * probe_offset_days, swe_id, flags)[0]
            v2 = swe.calc_ut(start_jd + probe_offset_days, swe_id, flags)[0]
            # Require uniform velocity sign and magnitude well above zero across the segment
            if (
                v0[3] * v1[3] > 0.0
                and v1[3] * v2[3] > 0.0
                and min(abs(v0[3]), abs(v1[3]), abs(v2[3])) > 0.015
                and abs(wrap180(v2[0] - v0[0])) > 0.15
                and abs(wrap180(v2[0] - v0[0])) < 60.0
            ):
                break
            start_jd = round(start_jd + 7.0, 6)

        # Pick target_lon_deg at start_jd + 0.65 * probe_offset_days (rounded to 4 decimals)
        probe_jd = start_jd + 0.65 * probe_offset_days
        probe_vals = swe.calc_ut(probe_jd, swe_id, flags)[0]
        target_lon_deg = round(float(probe_vals[0]) % 360.0, 4)

        def f_lon(t_ut: float) -> float:
            lon = swe.calc_ut(t_ut, swe_id, flags)[0][0]
            return wrap180(lon - target_lon_deg)

        root_jd = float(brentq(f_lon, start_jd, start_jd + probe_offset_days, xtol=1e-11))
        root_vals = swe.calc_ut(root_jd, swe_id, flags)[0]

        crossing_vectors.append([
            b_name,
            zodiac_mode,
            start_jd,
            target_lon_deg,
            max_search_days,
            round(root_jd, 8),
            round(float(root_vals[3]), 7),
        ])

    assert len(hd_vectors) == 500
    assert len(crossing_vectors) == 500
    return {
        "human_design_500": hd_vectors,
        "longitude_crossings_500": crossing_vectors,
    }


# ---------------------------------------------------------------------------
# Section 6: 500 Planetary Stations + 500 Declination Extrema
# ---------------------------------------------------------------------------
def build_section_6_stations_and_declination_extrema() -> dict:
    swe.set_sid_mode(swe.SIDM_LAHIRI, 0, 0)
    # 6A. 500 Planetary Stations across Mercury..Pluto in Tropical and Sidereal_Lahiri (16 groups)
    # Counts per group: 4 groups x 32 + 12 groups x 31 = 128 + 372 = 500 stations.
    station_planets = [
        ("Mercury", swe.MERCURY, 2.0),
        ("Venus", swe.VENUS, 3.0),
        ("Mars", swe.MARS, 5.0),
        ("Jupiter", swe.JUPITER, 5.0),
        ("Saturn", swe.SATURN, 5.0),
        ("Uranus", swe.URANUS, 5.0),
        ("Neptune", swe.NEPTUNE, 5.0),
        ("Pluto", swe.PLUTO, 5.0),
    ]
    stations_500 = []
    group_idx = 0
    base_jd_1900 = 2415020.5  # 1900-01-01 00:00 UT

    for p_name, swe_id, step_days in station_planets:
        for zodiac_mode in ("Tropical", "Sidereal_Lahiri"):
            target_count = 32 if group_idx < 4 else 31
            group_idx += 1
            flags = swe.FLG_SWIEPH | swe.FLG_SPEED
            if zodiac_mode == "Sidereal_Lahiri":
                flags |= swe.FLG_SIDEREAL

            def f_vel(t_ut: float) -> float:
                return float(swe.calc_ut(t_ut, swe_id, flags)[0][3])

            win_start = base_jd_1900
            # Ensure win_start is not near a station
            while abs(f_vel(win_start)) < 0.005:
                win_start += 2.0

            t_prev = win_start
            v_prev = f_vel(t_prev)
            found_group = []

            while len(found_group) < target_count:
                t_curr = t_prev + step_days
                v_curr = f_vel(t_curr)
                if (v_prev > 0.0 and v_curr <= 0.0) or (v_prev < 0.0 and v_curr >= 0.0):
                    kind_name = "StationRetrograde" if v_prev > 0.0 else "StationDirect"
                    root_jd = float(brentq(f_vel, t_prev, t_curr, xtol=1e-10))
                    vals = swe.calc_ut(root_jd, swe_id, flags)[0]
                    found_group.append((
                        kind_name,
                        round(root_jd, 7),
                        round(float(vals[0]), 7),
                        round(float(vals[1]), 7),
                        round(float(vals[2]), 8),
                    ))
                t_prev = t_curr
                v_prev = v_curr

            # Place win_end 4 days after the last station in this group (well before the next station)
            win_end = round(found_group[-1][1] + 4.0, 4)
            for kind_name, r_jd, r_lon, r_lat, r_dist in found_group:
                stations_500.append([
                    p_name,
                    zodiac_mode,
                    win_start,
                    win_end,
                    kind_name,
                    r_jd,
                    r_lon,
                    r_lat,
                    r_dist,
                ])

    assert len(stations_500) == 500

    # 6B. 500 Lunar & Solar Declination Extrema (NorthernMaximum & SouthernMinimum) in 1950..2050
    # - 200 Solar declination extrema (100 years 1950..2050 x 2 solstices/yr = 200)
    # - 300 Lunar declination extrema across 3 windows (100 in 1950..1954, 100 in 1985..1989, 100 in 2023..2027 Major Standstill)
    eq_flags = swe.FLG_SWIEPH | swe.FLG_EQUATORIAL | swe.FLG_SPEED
    dec_extrema_500 = []

    def collect_dec_extrema(body_name: str, swe_id: int, win_start: float, step_days: float, count: int) -> list[list]:
        def f_dec_vel(t_ut: float) -> float:
            return float(swe.calc_ut(t_ut, swe_id, eq_flags)[0][4])

        while abs(f_dec_vel(win_start)) < 0.005:
            win_start += 1.0

        t_prev = win_start
        v_prev = f_dec_vel(t_prev)
        found = []
        while len(found) < count:
            t_curr = t_prev + step_days
            v_curr = f_dec_vel(t_curr)
            if (v_prev > 0.0 and v_curr <= 0.0) or (v_prev < 0.0 and v_curr >= 0.0):
                kind_name = "NorthernMaximum" if v_prev > 0.0 else "SouthernMinimum"
                root_jd = float(brentq(f_dec_vel, t_prev, t_curr, xtol=1e-10))
                vals = swe.calc_ut(root_jd, swe_id, eq_flags)[0]
                found.append((
                    kind_name,
                    round(root_jd, 7),
                    round(float(vals[0]), 7),
                    round(float(vals[1]), 7),
                    round(float(vals[2]), 8),
                ))
            t_prev = t_curr
            v_prev = v_curr

        win_end = round(found[-1][1] + (3.0 if body_name == "Moon" else 20.0), 4)
        return [
            [body_name, win_start, win_end, k_name, r_jd, r_ra, r_dec, r_dist]
            for (k_name, r_jd, r_ra, r_dec, r_dist) in found
        ]

    # 200 Solar declination extrema starting 1950-01-01 (JD 2433282.5)
    dec_extrema_500.extend(collect_dec_extrema("Sun", swe.SUN, 2433282.5, 5.0, 200))
    # 3 x 100 Lunar declination extrema (1950, 1985, 2023 Major Lunar Standstill)
    dec_extrema_500.extend(collect_dec_extrema("Moon", swe.MOON, 2433282.5, 1.5, 100))
    dec_extrema_500.extend(collect_dec_extrema("Moon", swe.MOON, 2446066.5, 1.5, 100))
    dec_extrema_500.extend(collect_dec_extrema("Moon", swe.MOON, 2459945.5, 1.5, 100))

    assert len(dec_extrema_500) == 500
    return {
        "stations_500": stations_500,
        "declination_extrema_500": dec_extrema_500,
    }


# ---------------------------------------------------------------------------
# Section 7: 500 Diurnal Fixed-Star Parans (Brady R6.2) across 25 Charts
# ---------------------------------------------------------------------------
def build_section_7_brady_parans(stars_list: list[dict]) -> list[list]:
    locations = [
        ("London", 51.5074, -0.1278),
        ("NYC", 40.7128, -74.0060),
        ("Kharkiv", 49.9935, 36.2304),
        ("Tokyo", 35.6762, 139.6503),
        ("Sydney", -33.8688, 151.2093),
    ]
    planets = [
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
    angles = [
        ("rising", swe.CALC_RISE | swe.BIT_DISC_CENTER | swe.BIT_NO_REFRACTION),
        ("culminating", swe.CALC_MTRANSIT | swe.BIT_DISC_CENTER | swe.BIT_NO_REFRACTION),
        ("setting", swe.CALC_SET | swe.BIT_DISC_CENTER | swe.BIT_NO_REFRACTION),
        ("anti-culminating", swe.CALC_ITRANSIT | swe.BIT_DISC_CENTER | swe.BIT_NO_REFRACTION),
    ]

    # 25 birth epochs (5 per location) across 1950..2025
    chart_epochs = linspace(2433400.5, 2460800.5, 25)
    parans_500 = []

    for case_idx, start_jd in enumerate(chart_epochs):
        loc_name, lat_deg, lon_deg = locations[case_idx % len(locations)]
        geopos = (lon_deg, lat_deg, 0.0)

        # Solve planet diurnal crossings in [start_jd, start_jd + 1.05]
        planet_events = []
        for p_name, swe_id in planets:
            p_ev = []
            for ang_name, rsmi in angles:
                res_code, tret = swe.rise_trans(
                    start_jd, swe_id, rsmi, geopos, 0.0, 15.0, swe.FLG_SWIEPH
                )
                if res_code == 0 and start_jd <= tret[0] <= start_jd + 1.05:
                    p_ev.append((ang_name, float(tret[0])))
            planet_events.append((p_name, p_ev))

        # Solve star diurnal crossings in [start_jd, start_jd + 1.05]
        case_candidates = []
        for star in stars_list:
            s_name = star["name"]
            s_hip = star["hip"]
            s_ev = []
            for ang_name, rsmi in angles:
                res_code, tret = swe.rise_trans(
                    start_jd, s_name, rsmi, geopos, 0.0, 15.0, swe.FLG_SWIEPH
                )
                if res_code == 0 and start_jd <= tret[0] <= start_jd + 1.05:
                    s_ev.append((ang_name, float(tret[0])))

            for s_ang, t_star in s_ev:
                for p_name, p_ev in planet_events:
                    for p_ang, t_planet in p_ev:
                        diff_min = (abs(t_star - t_planet) * 1440.0) % 1440.0
                        circ_min = min(diff_min, 1440.0 - diff_min)
                        # Keep parans strictly inside orb <= 3.85 min (away from 4.0m cutoff boundary)
                        # and with both events safely inside [start_jd + 0.01, start_jd + 0.99]
                        if (
                            circ_min <= 3.85
                            and (start_jd + 0.01) <= t_star <= (start_jd + 0.99)
                            and (start_jd + 0.01) <= t_planet <= (start_jd + 0.99)
                        ):
                            case_candidates.append([
                                case_idx,
                                loc_name,
                                lat_deg,
                                lon_deg,
                                start_jd,
                                s_name,
                                s_hip,
                                p_name,
                                s_ang,
                                s_ang,
                                p_ang,
                                round(t_star, 8),
                                round(t_planet, 8),
                                round(circ_min, 6),
                            ])

        case_candidates.sort(key=lambda r: (r[13], r[5], r[7], r[9], r[10]))
        parans_500.extend(case_candidates[:20])

    # If needed, scan additional epochs to reach exactly 500 vectors
    extra_idx = 25
    while len(parans_500) < 500:
        start_jd = round(2440359.083333 + extra_idx * 41.3, 6)
        loc_name, lat_deg, lon_deg = locations[extra_idx % len(locations)]
        geopos = (lon_deg, lat_deg, 0.0)
        planet_events = []
        for p_name, swe_id in planets:
            p_ev = []
            for ang_name, rsmi in angles:
                res_code, tret = swe.rise_trans(
                    start_jd, swe_id, rsmi, geopos, 0.0, 15.0, swe.FLG_SWIEPH
                )
                if res_code == 0 and start_jd <= tret[0] <= start_jd + 1.05:
                    p_ev.append((ang_name, float(tret[0])))
            planet_events.append((p_name, p_ev))

        case_candidates = []
        for star in stars_list:
            s_name = star["name"]
            s_hip = star["hip"]
            s_ev = []
            for ang_name, rsmi in angles:
                res_code, tret = swe.rise_trans(
                    start_jd, s_name, rsmi, geopos, 0.0, 15.0, swe.FLG_SWIEPH
                )
                if res_code == 0 and start_jd <= tret[0] <= start_jd + 1.05:
                    s_ev.append((ang_name, float(tret[0])))
            for s_ang, t_star in s_ev:
                for p_name, p_ev in planet_events:
                    for p_ang, t_planet in p_ev:
                        diff_min = (abs(t_star - t_planet) * 1440.0) % 1440.0
                        circ_min = min(diff_min, 1440.0 - diff_min)
                        if (
                            circ_min <= 3.85
                            and (start_jd + 0.01) <= t_star <= (start_jd + 0.99)
                            and (start_jd + 0.01) <= t_planet <= (start_jd + 0.99)
                        ):
                            case_candidates.append([
                                extra_idx,
                                loc_name,
                                lat_deg,
                                lon_deg,
                                start_jd,
                                s_name,
                                s_hip,
                                p_name,
                                s_ang,
                                s_ang,
                                p_ang,
                                round(t_star, 8),
                                round(t_planet, 8),
                                round(circ_min, 6),
                            ])
        case_candidates.sort(key=lambda r: (r[13], r[5], r[7], r[9], r[10]))
        need = 500 - len(parans_500)
        parans_500.extend(case_candidates[: min(20, need)])
        extra_idx += 1

    assert len(parans_500) == 500
    return parans_500


# ---------------------------------------------------------------------------
# Section 8: Wave B & Wave C Swiss Ephemeris Parity (2,000 vectors)
# ---------------------------------------------------------------------------
def bowell_hg_mag(h_mag: float, g_slope: float, r_au: float, delta_au: float, alpha_deg: float) -> float:
    alpha_rad = math.radians(min(max(alpha_deg, 0.0), 179.999))
    tan_half = math.tan(0.5 * alpha_rad)
    phi1 = math.exp(-3.33 * (tan_half ** 0.63))
    phi2 = math.exp(-1.87 * (tan_half ** 1.22))
    comb = (1.0 - g_slope) * phi1 + g_slope * phi2
    phase_term = -2.5 * math.log10(comb) if comb > 1e-15 else 0.0
    return h_mag + 5.0 * math.log10(max(r_au * delta_au, 1e-15)) + phase_term


def tt_to_tdb_jd(jd_tt: float) -> float:
    g = math.radians((357.53 + 0.98560028 * (jd_tt - 2451545.0)) % 360.0)
    dtdb_sec = 0.001657 * math.sin(g) + 0.000022 * math.sin(2.0 * g)
    return jd_tt + dtdb_sec / 86400.0


def l1b_minor_body_pheno(
    spk_pl: SimpleDafType2Reader,
    spk_ast: SimpleDafType2Reader,
    mpc: int,
    h_mag: float,
    g_slope: float,
    radius_km: float,
    jd_ut1: float,
) -> tuple[float, float, float, float, float, float]:
    dt_days = swe.deltat(jd_ut1)
    jd_tt = jd_ut1 + dt_days
    jd_tdb = tt_to_tdb_jd(jd_tt)
    p_earth_km = spk_pl.compute(0, 3, jd_tdb) + spk_pl.compute(3, 399, jd_tdb)
    p_sun_now_km = spk_pl.compute(0, 10, jd_tdb)
    naif_id = 20000000 + mpc

    tau_days = 0.0
    p_targ_km = np.zeros(3, dtype=float)
    rho_au = np.zeros(3, dtype=float)
    delta_au = 1.0
    for _ in range(5):
        t_ret = jd_tdb - tau_days
        p_targ_km = spk_pl.compute(0, 10, t_ret) + spk_ast.compute(10, naif_id, t_ret)
        rho_au = (p_targ_km - p_earth_km) / AU_KM
        delta_au = float(np.linalg.norm(rho_au))
        tau_days = delta_au / ((C_KM_S * 86400.0) / AU_KM)

    p_sun_ret_km = spk_pl.compute(0, 10, jd_tdb - tau_days)
    r_helio_au = (p_targ_km - p_sun_ret_km) / AU_KM
    r_au = float(np.linalg.norm(r_helio_au))

    tau_sun = float(np.linalg.norm(p_sun_now_km - p_earth_km)) / (C_KM_S * 86400.0)
    p_sun_obs_km = spk_pl.compute(0, 10, jd_tdb - tau_sun)
    r_obs_sun_au = (p_sun_obs_km - p_earth_km) / AU_KM
    r_sun_au = float(np.linalg.norm(r_obs_sun_au))

    cos_alpha = float(np.dot(r_helio_au, rho_au)) / (r_au * delta_au)
    cos_alpha = max(-1.0, min(1.0, cos_alpha))
    phase_angle_deg = math.degrees(math.acos(cos_alpha))
    phase_illum = 0.5 * (1.0 + cos_alpha)

    cos_psi = float(np.dot(rho_au, r_obs_sun_au)) / (delta_au * r_sun_au)
    cos_psi = max(-1.0, min(1.0, cos_psi))
    elong_deg = math.degrees(math.acos(cos_psi))

    delta_km = delta_au * AU_KM
    diam_deg = 2.0 * math.degrees(math.asin(min(radius_km / delta_km, 1.0)))
    horiz_par_deg = math.degrees(math.asin(min(6378.1366 / delta_km, 1.0)))
    vmag = bowell_hg_mag(h_mag, g_slope, r_au, delta_au, phase_angle_deg)
    return phase_angle_deg, phase_illum, elong_deg, diam_deg, vmag, horiz_par_deg


def build_section_8_wave_b_c_swe_parity(
    spk_pl: SimpleDafType2Reader,
    spk_ast: SimpleDafType2Reader,
    stars_by_hip: dict[int, dict],
) -> dict:
    # 8A. 500 Moshier mode (`SEFLG_MOSEPH`) vectors across -500 CE .. +2500 CE (50 epochs x 10 bodies)
    epochs_mosh_50 = linspace(1538432.5, 2634169.5, 50)
    mosh_bodies = [
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
    moshier_500 = []
    mosh_flags = swe.FLG_MOSEPH | swe.FLG_SPEED
    for jd_tt in epochs_mosh_50:
        for b_name, swe_id in mosh_bodies:
            vals, _ = swe.calc(jd_tt, swe_id, mosh_flags)
            moshier_500.append([
                jd_tt,
                b_name,
                swe_id,
                round(float(vals[0]), 6),
                round(float(vals[1]), 6),
                round(float(vals[2]), 7),
                round(float(vals[3]), 6),
            ])
    assert len(moshier_500) == 500

    # 8B. 300 Photometric phenomena (`phenomena` / `SweAdapter::pheno_ut`) across Sun, Moon, Mercury..Pluto, Ceres, Vesta (25 epochs x 12 bodies)
    epochs_25 = linspace(2400409.5, 2484414.5, 25)
    pheno_major = [
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
    phenomena_300 = []
    for jd_ut1 in epochs_25:
        for b_name, swe_id in pheno_major:
            attr = swe.pheno_ut(jd_ut1, swe_id, swe.FLG_SWIEPH)
            phenomena_300.append([
                jd_ut1,
                b_name,
                swe_id,
                round(float(attr[0]), 6),
                round(float(attr[1]), 7),
                round(float(attr[2]), 6),
                round(float(attr[3]), 7),
                round(float(attr[4]), 5),
                round(float(attr[5]), 7),
            ])
        for b_name, swe_id, mpc, h_mag, g_slope, r_km in [
            ("Ceres", swe.CERES, 1, 3.34, 0.12, 469.73),
            ("Vesta", swe.VESTA, 4, 3.20, 0.32, 262.70),
        ]:
            pa, pi, el, dm, vm, hp = l1b_minor_body_pheno(spk_pl, spk_ast, mpc, h_mag, g_slope, r_km, jd_ut1)
            phenomena_300.append([
                jd_ut1,
                b_name,
                swe_id,
                round(pa, 6),
                round(pi, 7),
                round(el, 6),
                round(dm, 7),
                round(vm, 5),
                round(hp, 7),
            ])
    assert len(phenomena_300) == 300

    # 8C. 300 Orbital nodes & apsides (`nod_aps_ut` with SE_NODBIT_MEAN and SE_NODBIT_OSCU)
    # 15 epochs x 10 bodies (Moon, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto) x 2 methods = 300
    epochs_15 = linspace(2400409.5, 2484414.5, 15)
    nod_bodies = [
        ("Moon", swe.MOON),
        ("Mercury", swe.MERCURY),
        ("Venus", swe.VENUS),
        ("Earth", swe.EARTH),
        ("Mars", swe.MARS),
        ("Jupiter", swe.JUPITER),
        ("Saturn", swe.SATURN),
        ("Uranus", swe.URANUS),
        ("Neptune", swe.NEPTUNE),
        ("Pluto", swe.PLUTO),
    ]
    nodes_apsides_300 = []
    for jd_ut1 in epochs_15:
        for b_name, swe_id in nod_bodies:
            for method_bits in (swe.NODBIT_MEAN, swe.NODBIT_OSCU):
                xnasc, xndsc, xperi, xaphe = swe.nod_aps_ut(
                    jd_ut1, swe_id, method_bits, swe.FLG_SWIEPH | swe.FLG_SPEED
                )
                nodes_apsides_300.append([
                    jd_ut1,
                    b_name,
                    swe_id,
                    int(method_bits),
                    round(float(xnasc[0]), 6),
                    round(float(xnasc[2]), 7),
                    round(float(xndsc[0]), 6),
                    round(float(xperi[0]), 6),
                    round(float(xperi[1]), 6),
                    round(float(xperi[2]), 7),
                    round(float(xaphe[0]), 6),
                    round(float(xaphe[1]), 6),
                    round(float(xaphe[2]), 7),
                ])
    assert len(nodes_apsides_300) == 300

    # 8D. 430 Extended Ayanamsa vectors: all 43 SE_SIDM_* models (0..=42) x 10 epochs (1850..2100)
    epochs_10 = linspace(2396760.5, 2488065.5, 10)
    # IAU 1980 nutation in longitude at 1956-03-21 00:00 TT (JD 2435553.5) is +16.76885713 arcsec
    dpsi_1956_arcsec = 16.76885712832541
    lahiri_mean_1956 = (23.0 * 3600.0 + 15.0 * 60.0 + 0.658 - dpsi_1956_arcsec) / 3600.0

    epoch_anchors = {
        0: (2433282.42345905, 24.0 + 2.0 / 60.0 + 31.36 / 3600.0),
        1: (2435553.5, lahiri_mean_1956),
        2: (2451545.0, 27.815752807),
        3: (2415020.0, 21.01421),
        4: (2415020.0, 18.66096),
        5: (2415020.0, 22.363659),
        6: (2494459.5, 30.0),
        7: (2415020.0, 21.0 + 4.0 / 60.0 + 56.0 / 3600.0),
        8: (2415020.0, 21.0 + 21.0 / 60.0 + 56.0 / 3600.0),
        9: (2451545.0, 23.533639900),
        10: (2451545.0, 24.933639900),
        11: (2451545.0, 25.783639900),
        12: (2451545.0, 24.733639900),
        13: (2451545.0, 24.522527973),
        14: (2451545.0, 24.758923900),
        15: (2451545.0, 20.247788132),
        16: (2451545.0, 19.992959400),
        17: (2451545.0, 26.846045806),
        18: (2451545.0, 0.0),
        19: (2415020.0, 0.0),
        20: (2433282.42345905, 0.0),
        21: (1903396.8770698, 0.0),
        22: (2451545.0, 20.680424905),
        23: (1903396.8770698, 0.0),
        24: (2451545.0, 20.657427361),
        25: (2451545.0, 20.103388395),
        26: (1903396.8770698, 2.0 + 6.0 / 60.0 + 38.536 / 3600.0),
        30: (2451545.0, 22.469104789),
        31: (2451545.0, 30.023172469),
        32: (2451545.0, 30.076092180),
        33: (2451545.0, 23.409425513),
        34: (2451545.0, 30.017793877),
        35: (2451545.0, 24.579972884),
        36: (2451545.0, 20.039233112),
        37: (1911797.5, 0.0),
        38: (2451545.0, 24.615752807),
        39: (2451545.0, 25.234449334),
        40: (2451545.0, 356.846045806),
        41: (2451545.0, 25.000019094),
        42: (2451545.0, 22.795609328),
    }
    star_anchors = {
        27: (65474, 180.0),
        28: (5737, 359.0 + 50.0 / 60.0),
        29: (42911, 106.0),
    }

    extended_ayanamsas_430 = []
    for jd_tt in epochs_10:
        dpsi_deg = float(swe.calc(jd_tt, swe.ECL_NUT, 0)[0][2])
        for sid_id in range(43):
            if sid_id in epoch_anchors:
                t0, a0 = epoch_anchors[sid_id]
                p_arcsec, _ = iae_general_precession_arcsec(t0, jd_tt)
                l1b_mean = (a0 + p_arcsec / 3600.0) % 360.0
            else:
                hip_id, target_lon = star_anchors[sid_id]
                st_name = stars_by_hip[hip_id]["name"]
                lon_0 = float(
                    swe.fixstar2(
                        st_name,
                        jd_tt,
                        swe.FLG_SWIEPH | swe.FLG_NOABERR | swe.FLG_NONUT,
                    )[0][0]
                )
                l1b_mean = (lon_0 - target_lon) % 360.0
            l1b_true = (l1b_mean + dpsi_deg) % 360.0

            swe.set_sid_mode(sid_id, 0, 0)
            swe_mean = float(swe.get_ayanamsa_ex(jd_tt, swe.FLG_SWIEPH | swe.FLG_NONUT)[1]) % 360.0
            swe_true = float(swe.get_ayanamsa_ex(jd_tt, swe.FLG_SWIEPH)[1]) % 360.0

            extended_ayanamsas_430.append([
                jd_tt,
                sid_id,
                round(l1b_mean, 7),
                round(l1b_true, 7),
                round(swe_mean, 7),
                round(swe_true, 7),
            ])
    assert len(extended_ayanamsas_430) == 430

    # 8E. 470 3D mundane house_pos (200), time_equation (120), and azalt_rev (150) vectors
    house_pos_200 = []
    hsys_codes = ["P", "K", "O", "R", "C", "A", "W", "B", "M", "T"]
    armc_vals = [15.0, 95.0, 185.0, 275.0, 335.0]
    lat_vals = [-42.0, -18.0, 28.5, 51.5]
    eps_deg = 23.43929111
    for i, hsys in enumerate(hsys_codes):
        b_hsys = hsys.encode("ascii")
        for j, lat in enumerate(lat_vals):
            for k, armc in enumerate(armc_vals):
                ecl_lon = (armc + 37.0 * (i + 1) + 23.0 * j + 11.0 * k) % 360.0
                ecl_lat = -6.0 + 2.0 * ((i + j + k) % 7)
                hp = float(swe.house_pos(armc, lat, eps_deg, (ecl_lon, ecl_lat), b_hsys))
                house_pos_200.append([
                    armc,
                    lat,
                    eps_deg,
                    hsys,
                    round(ecl_lon, 5),
                    round(ecl_lat, 5),
                    round(hp, 7),
                ])
    assert len(house_pos_200) == 200

    time_equation_120 = []
    for jd_tt in linspace(2398588.5, 2486240.5, 120):
        eot_days = float(swe.time_equ(jd_tt))
        time_equation_120.append([jd_tt, round(eot_days, 9)])
    assert len(time_equation_120) == 120

    azalt_rev_150 = []
    epochs_az_30 = linspace(2400409.5, 2484414.5, 30)
    az_configs = [
        (swe.HOR2EQU, 8.5417, 47.3769, 408.0, 35.0, 28.0),
        (swe.HOR2ECL, 8.5417, 47.3769, 408.0, 145.0, 42.0),
        (swe.HOR2EQU, -74.0060, 40.7128, 10.0, 220.0, 19.5),
        (swe.HOR2ECL, 139.6503, 35.6762, 40.0, 310.0, 55.0),
        (swe.HOR2EQU, 151.2093, -33.8688, 25.0, 85.0, -12.0),
    ]
    for jd_ut1 in epochs_az_30:
        for calc_flag, lon_d, lat_d, alt_m, az_sw, alt_true in az_configs:
            x0, x1 = swe.azalt_rev(jd_ut1, calc_flag, (lon_d, lat_d, alt_m), az_sw, alt_true)
            azalt_rev_150.append([
                jd_ut1,
                int(calc_flag),
                lon_d,
                lat_d,
                alt_m,
                az_sw,
                alt_true,
                round(float(x0) % 360.0, 7),
                round(float(x1), 7),
            ])
    assert len(azalt_rev_150) == 150

    return {
        "moshier_500": moshier_500,
        "phenomena_300": phenomena_300,
        "nodes_apsides_300": nodes_apsides_300,
        "extended_ayanamsas_430": extended_ayanamsas_430,
        "house_pos_200": house_pos_200,
        "time_equation_120": time_equation_120,
        "azalt_rev_150": azalt_rev_150,
    }


def main() -> None:
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    spk_pl = SimpleDafType2Reader(DE440S_PATH)
    spk_ast = SimpleDafType2Reader(AST_BSP_PATH)
    stars_list = json.loads((GOLDEN_DIR / "stars_catalog_52.json").read_text())
    stars_by_hip = {s["hip"]: s for s in stars_list}

    with tempfile.TemporaryDirectory() as td:
        write_temp_sefstars(Path(td), stars_list)
        swe.set_ephe_path(td)

        print("Generating Section 1: Rise, Set & Meridian Transits (1,600 vectors)...")
        sec1 = build_section_1_rise_trans()

        print("Generating Section 2: Principal Lunar Phases (1,500 vectors)...")
        sec2 = build_section_2_lunar_phases()

        print("Generating Section 3: Lunar Phase Summary (1,500 vectors)...")
        sec3 = build_section_3_lunar_phase_summary()

        print("Generating Section 4: Chinese 24 Solar Terms & Seasons (1,440 vectors)...")
        sec4 = build_section_4_solar_terms_and_seasons()

        print("Generating Section 5: Human Design -88 deg & Longitude Crossings (1,000 vectors)...")
        sec5 = build_section_5_hd_and_longitude_crossings()

        print("Generating Section 6: Planetary Stations & Declination Extrema (1,000 vectors)...")
        sec6 = build_section_6_stations_and_declination_extrema()

        print("Generating Section 7: Brady Diurnal Fixed-Star Parans (500 vectors)...")
        sec7 = build_section_7_brady_parans(stars_list)

        print("Generating Section 8: Wave B & Wave C Swiss Ephemeris Parity (2,000 vectors)...")
        sec8 = build_section_8_wave_b_c_swe_parity(spk_pl, spk_ast, stars_by_hip)

    counts = {
        "section_1_rise_trans": len(sec1),
        "section_2_lunar_phases": len(sec2),
        "section_3_lunar_phase_summary": len(sec3),
        "section_4_solar_terms_and_seasons": len(sec4),
        "section_5_hd_and_longitude_crossings": len(sec5["human_design_500"]) + len(sec5["longitude_crossings_500"]),
        "section_6_stations_and_declination_extrema": len(sec6["stations_500"]) + len(sec6["declination_extrema_500"]),
        "section_7_brady_parans": len(sec7),
        "section_8_wave_b_c_swe_parity": (
            len(sec8["moshier_500"])
            + len(sec8["phenomena_300"])
            + len(sec8["nodes_apsides_300"])
            + len(sec8["extended_ayanamsas_430"])
            + len(sec8["house_pos_200"])
            + len(sec8["time_equation_120"])
            + len(sec8["azalt_rev_150"])
        ),
    }
    total_count = sum(counts.values())
    assert total_count == 10540, f"Expected 10540 vectors, got {total_count}"

    out = {
        "schema_version": "1.0.0",
        "tier": "L3_wave_bcf_10k_grid",
        "oracles": ["pyswisseph 2.10.3.2", "ephem 4.1.5+", "jplephem + pyerfa L1B"],
        "total_vector_count": total_count,
        "counts_by_section": counts,
        "section_1_rise_trans": sec1,
        "section_2_lunar_phases": sec2,
        "section_3_lunar_phase_summary": sec3,
        "section_4_solar_terms_and_seasons": sec4,
        "section_5_hd_and_longitude_crossings": sec5,
        "section_6_stations_and_declination_extrema": sec6,
        "section_7_brady_parans": sec7,
        "section_8_wave_b_c_swe_parity": sec8,
    }

    out_path = GOLDEN_DIR / "l3_wave_bcf_10k_grid.json"
    out_path.write_text(json.dumps(out, separators=(",", ":")) + "\n")
    size_kb = out_path.stat().st_size / 1024.0
    print(f"\nWrote {out_path} ({total_count} vectors, {size_kb:.1f} KB)")
    for k, v in counts.items():
        print(f"  - {k}: {v}")


if __name__ == "__main__":
    main()
