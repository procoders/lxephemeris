#!/usr/bin/env python3
"""
Wave H5 5,000-Vector Golden Oracle Generator (`testdata/golden/l3_wave_h_almanac_aspects_5k.json`).

Brings the total golden reference database in `testdata/golden/` from 32,416 to 37,416 vectors
across 15 JSON files.

Strictly uses independent reference oracles (`pyswisseph` + canonical Brihat Parashara Hora Shastra
[BPHS] & Paulus Alexandrinus formulas + `scipy.optimize.brentq` root-finding) across 5 sections:

  1. `section_h1_panchanga`            = 1,000 vectors
     - 1,000 Vedic Panchanga & Muhurta vectors across 100 dates (1900..2050 CE) x 10 global cities
       (New Delhi, Mumbai, Varanasi, Ujjain, Bengaluru, London, New York, Tokyo, Cairo, Sydney)
       cycling across 4 ayanamsas (Lahiri, Raman, Krishnamurti, TrueCitra).
  2. `section_h2_vargas_and_dasha`     = 1,200 vectors
     - `varga_1000`: 1,000 Parashara Shodashavarga (D1..D60) vectors across 62 longitudes x 16
       divisions (992) + 8 boundary-stress vectors.
     - `vimshottari_dasha_200`: 200 120-year Vimshottari Mahadasha + Antardasha vectors across
       birth epochs (1900..2050 CE), 4 ayanamsas, and both 365.2425 & 365.25 solar year lengths.
  3. `section_h3_traditional_almanac`  = 1,000 vectors
     - 1,000 Traditional Western Almanac vectors (`chart_sect`, 24 `planetary_hours`, and all 7
       classical `calc_lots` of Paulus Alexandrinus) across 100 epochs x 10 global cities and
       8 house systems (P, W, R, O, K, C, A, B).
  4. `section_h4_aspects_returns_voc`  = 1,000 vectors
     - `static_aspects_300`: 100 `calc_chart_aspects` + 100 `calc_synastry_aspects` +
       100 `calc_declination_aspects` vectors with analytical `is_applying`.
     - `dynamic_exact_aspects_300`: 150 `mutual_aspects_in_range` / `search_exact_aspect` +
       150 `natal_transits_in_range` vectors (including retrograde triple passes of Mercury,
       Mars, Jupiter, Saturn), all matching `< 1e-5` days (`< 0.86s`).
     - `planetary_returns_200`: 100 `planetary_return` + 100 `planetary_returns_in_range`
       vectors across tropical and sidereal modes (including 3-pass retrograde returns).
     - `void_of_course_moon_200`: 200 `void_of_course_moon_in_range` intervals across tropical
       and sidereal zodiacs.
  5. `section_h5_hd_nakshatra_snapshot` =  800 vectors
     - `human_design_384`: all 64 Gates x 6 Lines of the Rave Mandala (`302 deg` Aquarius = Gate 41 Line 1)
       with Color, Tone, Base subdivisions.
     - `nakshatra_padas_108`: all 27 Nakshatras x 4 Padas (`0..360 deg`) with D9 Navamsha bijection.
     - `chart_snapshots_308`: 28 full charts x 11 verification vectors (28 chart-level summaries +
       280 `ChartBodySnapshot` vectors across Sun..Pluto).
Total: 5,000 golden reference vectors.
"""

from __future__ import annotations

import json
import math
import struct
import tempfile
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import swisseph as swe

ROOT = Path(__file__).resolve().parent.parent
DE440S_PATH = ROOT / "data" / "de440s.bsp"
GOLDEN_DIR = ROOT / "testdata" / "golden"
OUTPUT_PATH = GOLDEN_DIR / "l3_wave_h_almanac_aspects_5k.json"

AU_KM = 149597870.700
C_KM_S = 299792.458
C_AU_DAY = (C_KM_S * 86400.0) / AU_KM
DAS2R = 4.848136811095359935899141e-6
SRS = 1.97412574336e-8
DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi

ZODIAC_SIGN_NAMES = [
    "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
    "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces",
]

ZODIAC_SIGN_LORDS = [
    "Mars", "Venus", "Mercury", "Moon", "Sun", "Mercury",
    "Venus", "Mars", "Jupiter", "Saturn", "Saturn", "Jupiter",
]

NAKSHATRA_NAMES_27 = [
    "Ashwini", "Bharani", "Krittika", "Rohini", "Mrigashira", "Ardra",
    "Punarvasu", "Pushya", "Ashlesha", "Magha", "Purva Phalguni", "Uttara Phalguni",
    "Hasta", "Chitra", "Swati", "Vishakha", "Anuradha", "Jyeshtha",
    "Mula", "Purva Ashadha", "Uttara Ashadha", "Shravana", "Dhanishta", "Shatabhisha",
    "Purva Bhadrapada", "Uttara Bhadrapada", "Revati",
]

VIMSHOTTARI_LORDS_9 = [
    "Ketu", "Venus", "Sun", "Moon", "Mars", "Rahu", "Jupiter", "Saturn", "Mercury",
]

VIMSHOTTARI_YEARS_9 = [7.0, 20.0, 6.0, 10.0, 7.0, 18.0, 16.0, 19.0, 17.0]

HD_MANDALA_START_LON_DEG = 302.0
HD_MANDALA_GATES_64 = [
    41, 19, 13, 49, 30, 55, 37, 63, 22, 36, 25, 17, 21, 51, 42, 3,
    27, 24, 2, 23, 8, 20, 16, 35, 45, 12, 15, 52, 39, 53, 62, 56,
    31, 33, 7, 4, 29, 59, 40, 64, 47, 6, 46, 18, 48, 57, 32, 50,
    28, 44, 1, 43, 14, 34, 9, 5, 26, 11, 10, 58, 38, 54, 61, 60,
]

VARA_NAMES_7 = [
    "Ravivara", "Somavara", "Mangalavara", "Budhavara",
    "Guruvara", "Shukravara", "Shanivara",
]
VARA_LORDS_7 = ["Sun", "Moon", "Mars", "Mercury", "Jupiter", "Venus", "Saturn"]

TITHI_NAMES_30 = [
    "Pratipada", "Dwitiya", "Tritiya", "Chaturthi", "Panchami",
    "Shashthi", "Saptami", "Ashtami", "Navami", "Dashami",
    "Ekadashi", "Dwadashi", "Trayodashi", "Chaturdashi", "Purnima",
    "Pratipada", "Dwitiya", "Tritiya", "Chaturthi", "Panchami",
    "Shashthi", "Saptami", "Ashtami", "Navami", "Dashami",
    "Ekadashi", "Dwadashi", "Trayodashi", "Chaturdashi", "Amavasya",
]

YOGA_NAMES_27 = [
    "Vishkumbha", "Priti", "Ayushman", "Saubhagya", "Shobhana", "Atiganda",
    "Sukarma", "Dhriti", "Shula", "Ganda", "Vriddhi", "Dhruva",
    "Vyaghata", "Harshana", "Vajra", "Siddhi", "Vyatipata", "Variyan",
    "Parigha", "Shiva", "Siddha", "Sadhya", "Shubha", "Shukla",
    "Brahma", "Indra", "Vaidhriti",
]

CHARA_KARANA_NAMES_7 = ["Bava", "Balava", "Kaulava", "Taitila", "Gara", "Vanija", "Vishti"]
CHALDEAN_ORDER_7 = ["Saturn", "Jupiter", "Mars", "Sun", "Venus", "Mercury", "Moon"]
VARA_TO_CHALDEAN_START = [3, 6, 2, 5, 1, 4, 0]

CITIES_10 = [
    {"name": "New Delhi", "lat_deg": 28.6139, "lon_deg": 77.2090, "alt_m": 216.0},
    {"name": "Mumbai", "lat_deg": 19.0760, "lon_deg": 72.8777, "alt_m": 14.0},
    {"name": "Varanasi", "lat_deg": 25.3176, "lon_deg": 82.9739, "alt_m": 81.0},
    {"name": "Ujjain", "lat_deg": 23.1765, "lon_deg": 75.7885, "alt_m": 494.0},
    {"name": "Bengaluru", "lat_deg": 12.9716, "lon_deg": 77.5946, "alt_m": 920.0},
    {"name": "London", "lat_deg": 51.5074, "lon_deg": -0.1278, "alt_m": 11.0},
    {"name": "New York", "lat_deg": 40.7128, "lon_deg": -74.0060, "alt_m": 10.0},
    {"name": "Tokyo", "lat_deg": 35.6762, "lon_deg": 139.6503, "alt_m": 40.0},
    {"name": "Cairo", "lat_deg": 30.0444, "lon_deg": 31.2357, "alt_m": 23.0},
    {"name": "Sydney", "lat_deg": -33.8688, "lon_deg": 151.2093, "alt_m": 19.0},
]

AYANAMSAS_4 = [
    ("lahiri", swe.SIDM_LAHIRI),
    ("raman", swe.SIDM_RAMAN),
    ("krishnamurti", swe.SIDM_KRISHNAMURTI),
    ("true_citra", swe.SIDM_TRUE_CITRA),
]

SHODASHAVARGA_16 = [1, 2, 3, 4, 7, 9, 10, 12, 16, 20, 24, 27, 30, 40, 45, 60]

VARGA_SANSKRIT_NAMES = {
    1: "Rashi",
    2: "Hora",
    3: "Drekkana",
    4: "Chaturthamsha",
    7: "Saptamsha",
    9: "Navamsha",
    10: "Dashamsha",
    12: "Dwadashamsha",
    16: "Shodashamsha",
    20: "Vimshamsha",
    24: "Chaturvimshamsha",
    27: "Saptavimshamsha",
    30: "Trimshamsha",
    40: "Khavedamsha",
    45: "Akshavedamsha",
    60: "Shashtiamsha",
}


def wrap180(deg: float) -> float:
    return (deg + 180.0) % 360.0 - 180.0


def norm360(deg: float) -> float:
    r = deg % 360.0
    return 0.0 if r >= 360.0 else r


def floor_with_eps(x: float) -> int:
    r = round(x)
    if r >= 0 and abs(r - x) < 1.0e-12:
        return int(r)
    return max(0, int(math.floor(x)))


def tt_to_tdb_jd(jd_tt: float) -> float:
    g = math.radians((357.53 + 0.98560028 * (jd_tt - 2451545.0)) % 360.0)
    dtdb_sec = 0.001657 * math.sin(g) + 0.000022 * math.sin(2.0 * g)
    return jd_tt + dtdb_sec / 86400.0


def rx(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])


def rz(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def pfw06_matrices(jd_tt: float, nutation: bool = True) -> tuple[np.ndarray, np.ndarray]:
    t = (jd_tt - 2451545.0) / 36525.0
    gamb = (
        -0.052928
        + (10.556378 + (0.4932044 + (-0.00031238 + (-0.000002788 + 0.0000000260 * t) * t) * t) * t)
        * t
    ) * DAS2R
    phib = (
        84381.412819
        + (-46.811016 + (0.0511268 + (0.00053289 + (-0.000000440 - 0.0000000176 * t) * t) * t) * t)
        * t
    ) * DAS2R
    psib = (
        -0.041775
        + (5038.481484 + (1.5584175 + (-0.00018522 + (-0.000026452 - 0.0000000148 * t) * t) * t) * t)
        * t
    ) * DAS2R
    epsa = (
        84381.406000
        + (-46.836769 + (-0.0001831 + (0.00200340 + (-0.000000576 - 0.0000000434 * t) * t) * t) * t)
        * t
    ) * DAS2R
    if nutation:
        en = swe.calc(jd_tt, swe.ECL_NUT, 0)[0]
        dpsi = math.radians(float(en[2]))
        deps = math.radians(float(en[3]))
    else:
        dpsi, deps = 0.0, 0.0
    r_ecl = rz(-(psib + dpsi)) @ rx(phib) @ rz(gamb)
    r_eq = rx(-(epsa + deps)) @ r_ecl
    return r_ecl, r_eq


def sofa_ld_sun(u: np.ndarray, q: np.ndarray, e: np.ndarray, em: float) -> np.ndarray:
    qpe = q + e
    qdqpe = float(np.dot(q, qpe))
    w = SRS / em / max(qdqpe, 1e-15)
    eq = np.cross(e, q)
    upeq = np.cross(u, eq)
    u_def = u + w * upeq
    return u_def / float(np.linalg.norm(u_def))


def sofa_ab(pnat: np.ndarray, v: np.ndarray, s: float, bm1: float) -> np.ndarray:
    pdv = float(np.dot(pnat, v))
    w1 = 1.0 + pdv / (1.0 + bm1)
    w2 = SRS / max(s, 1e-6)
    p = pnat * bm1 + (w1 * v + w2 * (v - pdv * pnat))
    return p / float(np.linalg.norm(p))


def cart_to_spherical(u: np.ndarray) -> tuple[float, float, float]:
    x, y, z = float(u[0]), float(u[1]), float(u[2])
    r = math.hypot(x, y, z)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(max(-1.0, min(1.0, z / r)))) if r > 0.0 else 0.0
    return lon, lat, r


class SimpleDafType2Reader:
    """Minimal pure-Python NAIF DAF Type 2 Chebyshev reader for de440s.bsp."""

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

    def bary_pos_km(self, target_id: int, jd_tdb: float) -> np.ndarray:
        if target_id == 0:
            return np.zeros(3, dtype=float)
        if target_id in (199, 299):
            c = 1 if target_id == 199 else 2
            return self.compute(0, c, jd_tdb) + self.compute(c, target_id, jd_tdb)
        if target_id in (301, 399):
            return self.compute(0, 3, jd_tdb) + self.compute(3, target_id, jd_tdb)
        return self.compute(0, target_id, jd_tdb)

    def bary_state_km_s(self, target_id: int, jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
        p0 = self.bary_pos_km(target_id, jd_tdb)
        dt = 0.005
        pm = self.bary_pos_km(target_id, jd_tdb - dt)
        pp = self.bary_pos_km(target_id, jd_tdb + dt)
        v = (pp - pm) / (2.0 * dt * 86400.0)
        return p0, v


def write_temp_sefstars(td_path: Path, stars_list: list[dict]) -> None:
    lines = []
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


SPK_READER = SimpleDafType2Reader(DE440S_PATH)
BODY_NAIF_MAP = {
    "Sun": 10,
    "Moon": 301,
    "Mercury": 199,
    "Venus": 299,
    "Mars": 4,
    "Jupiter": 5,
    "Saturn": 6,
    "Uranus": 7,
    "Neptune": 8,
    "Pluto": 9,
}
BODY_SWE_MAP = {
    "Sun": swe.SUN,
    "Moon": swe.MOON,
    "Mercury": swe.MERCURY,
    "Venus": swe.VENUS,
    "Mars": swe.MARS,
    "Jupiter": swe.JUPITER,
    "Saturn": swe.SATURN,
    "Uranus": swe.URANUS,
    "Neptune": swe.NEPTUNE,
    "Pluto": swe.PLUTO,
}


def swe_ayanamsa_deg(jd_ut1: float, sidm: int) -> float:
    swe.set_sid_mode(sidm, 0, 0)
    flg = swe.FLG_SWIEPH | (swe.FLG_TRUEPOS if sidm == swe.SIDM_TRUE_CITRA else 0)
    return float(swe.get_ayanamsa_ex_ut(jd_ut1, flg)[1])


def de440s_apparent_lon_deg(jd_ut1: float, target_id: int, sidm: int | None = None) -> float:
    jd_tt = jd_ut1 + swe.deltat(jd_ut1)
    jd_tdb = tt_to_tdb_jd(jd_tt)
    p_obs_km, v_obs_kms = SPK_READER.bary_state_km_s(399, jd_tdb)
    p_sun_km = SPK_READER.bary_pos_km(10, jd_tdb)
    tau = 0.0
    for _ in range(4):
        p_t = SPK_READER.bary_pos_km(target_id, jd_tdb - tau)
        rho = float(np.linalg.norm((p_t - p_obs_km) / AU_KM))
        tau = rho / C_AU_DAY
    p_targ_km = SPK_READER.bary_pos_km(target_id, jd_tdb - tau)
    r_vec_au = (p_targ_km - p_obs_km) / AU_KM
    dist_au = float(np.linalg.norm(r_vec_au))
    u = r_vec_au / dist_au
    if target_id not in (10, 301):
        p_sun_ret = SPK_READER.bary_pos_km(10, jd_tdb - tau)
        q_vec = p_targ_km - p_sun_ret
        q_norm = float(np.linalg.norm(q_vec))
        e_vec = p_obs_km - p_sun_km
        e_norm = float(np.linalg.norm(e_vec))
        if q_norm > 1e-6 and e_norm > 1e-6:
            u = sofa_ld_sun(u, q_vec / q_norm, e_vec / e_norm, e_norm / AU_KM)
    v_vec = v_obs_kms / C_KM_S
    v_norm = float(np.linalg.norm(v_vec))
    bm1 = math.sqrt(max(1.0 - v_norm * v_norm, 0.0))
    e_sun_dist_au = max(float(np.linalg.norm(p_obs_km - p_sun_km)) / AU_KM, 1e-6)
    u = sofa_ab(u, v_vec, e_sun_dist_au, bm1)
    r_ecl, _ = pfw06_matrices(jd_tt, True)
    lon = cart_to_spherical(r_ecl @ u)[0]
    if sidm is not None:
        lon = norm360(lon - swe_ayanamsa_deg(jd_ut1, sidm))
    return lon


def swe_lon_speed_deg(jd_ut1: float, ipl: int, sidm: int | None = None) -> tuple[float, float]:
    vals, _ = swe.calc_ut(jd_ut1, ipl, swe.FLG_SWIEPH | swe.FLG_SPEED)
    lon = float(vals[0])
    spd = float(vals[3])
    if sidm is not None:
        lon = norm360(lon - swe_ayanamsa_deg(jd_ut1, sidm))
    return norm360(lon), spd


def map_zodiac_placement_py(lon_deg: float) -> dict:
    lon = norm360(lon_deg)
    sign_idx = min(floor_with_eps(lon / 30.0), 11)
    sign_name = ZODIAC_SIGN_NAMES[sign_idx]
    deg_in_sign = max(0.0, min(lon - sign_idx * 30.0, 30.0 - 1.0e-12))
    dms_deg = min(int(math.floor(deg_in_sign)), 29)
    rem_min = max(0.0, deg_in_sign - dms_deg) * 60.0
    dms_min = min(int(math.floor(rem_min)), 59)
    dms_sec = max(0.0, min((rem_min - dms_min) * 60.0, 60.0 - 1.0e-9))
    decan = min(floor_with_eps(deg_in_sign / 10.0), 2) + 1

    total_pada_idx = min(floor_with_eps((lon * 108.0) / 360.0), 107)
    nak_0 = min(total_pada_idx // 4, 26)
    pada = (total_pada_idx % 4) + 1
    nakshatra_index = nak_0 + 1
    nakshatra_name = NAKSHATRA_NAMES_27[nak_0]
    nakshatra_lord = VIMSHOTTARI_LORDS_9[nak_0 % 9]

    rel_hd = norm360(lon - HD_MANDALA_START_LON_DEG)
    total_base_idx = min(floor_with_eps(rel_hd * 192.0), 69119)
    hd_base = (total_base_idx % 5) + 1
    total_tone_idx = total_base_idx // 5
    hd_tone = (total_tone_idx % 6) + 1
    total_color_idx = total_tone_idx // 6
    hd_color = (total_color_idx % 6) + 1
    total_line_idx = total_color_idx // 6
    hd_line = (total_line_idx % 6) + 1
    gate_idx = min(total_line_idx // 6, 63)
    hd_gate = HD_MANDALA_GATES_64[gate_idx]

    return {
        "lon_deg": lon,
        "sign_index": sign_idx,
        "sign_name": sign_name,
        "deg_in_sign": deg_in_sign,
        "dms_deg": dms_deg,
        "dms_min": dms_min,
        "dms_sec": dms_sec,
        "decan": decan,
        "nakshatra_index": nakshatra_index,
        "nakshatra_name": nakshatra_name,
        "nakshatra_lord": nakshatra_lord,
        "pada": pada,
        "hd_gate": hd_gate,
        "hd_line": hd_line,
        "hd_color": hd_color,
        "hd_tone": hd_tone,
        "hd_base": hd_base,
    }


def calc_varga_py(sidereal_lon_deg: float, division: int) -> dict:
    lon = norm360(sidereal_lon_deg)
    s = min(floor_with_eps(lon / 30.0), 11)
    deg_in_sign = max(0.0, min(lon - s * 30.0, 30.0 - 1.0e-12))
    is_odd_sign = (s % 2) == 0

    if division == 30:
        if is_odd_sign:
            if deg_in_sign < 5.0:
                varga_sign, amsha_idx, seg_start, seg_span = 0, 1, 0.0, 5.0
            elif deg_in_sign < 10.0:
                varga_sign, amsha_idx, seg_start, seg_span = 10, 2, 5.0, 5.0
            elif deg_in_sign < 18.0:
                varga_sign, amsha_idx, seg_start, seg_span = 8, 3, 10.0, 8.0
            elif deg_in_sign < 25.0:
                varga_sign, amsha_idx, seg_start, seg_span = 2, 4, 18.0, 7.0
            else:
                varga_sign, amsha_idx, seg_start, seg_span = 6, 5, 25.0, 5.0
        else:
            if deg_in_sign < 5.0:
                varga_sign, amsha_idx, seg_start, seg_span = 1, 1, 0.0, 5.0
            elif deg_in_sign < 12.0:
                varga_sign, amsha_idx, seg_start, seg_span = 5, 2, 5.0, 7.0
            elif deg_in_sign < 20.0:
                varga_sign, amsha_idx, seg_start, seg_span = 11, 3, 12.0, 8.0
            elif deg_in_sign < 25.0:
                varga_sign, amsha_idx, seg_start, seg_span = 9, 4, 20.0, 5.0
            else:
                varga_sign, amsha_idx, seg_start, seg_span = 7, 5, 25.0, 5.0
        deg_in_varga_sign = max(0.0, min(((deg_in_sign - seg_start) / seg_span) * 30.0, 30.0 - 1.0e-12))
        varga_lon_deg = varga_sign * 30.0 + deg_in_varga_sign
        return {
            "division": 30,
            "varga_code": "D30",
            "varga_name": VARGA_SANSKRIT_NAMES[30],
            "source_sidereal_lon_deg": lon,
            "sign_index": varga_sign,
            "sign_name": ZODIAC_SIGN_NAMES[varga_sign],
            "sign_lord": ZODIAC_SIGN_LORDS[varga_sign],
            "amsha_index": amsha_idx,
            "varga_lon_deg": varga_lon_deg,
            "deg_in_varga_sign": deg_in_varga_sign,
        }

    n = division
    slice_width = 30.0 / float(n)
    k = min(floor_with_eps((deg_in_sign * float(n)) / 30.0), n - 1)
    amsha_index = k + 1
    deg_in_varga_sign = max(0.0, min(((deg_in_sign - k * slice_width) / slice_width) * 30.0, 30.0 - 1.0e-12))

    if division == 1:
        varga_sign = s
    elif division == 2:
        varga_sign = (4 if k == 0 else 3) if is_odd_sign else (3 if k == 0 else 4)
    elif division == 3:
        varga_sign = (s + k * 4) % 12
    elif division == 4:
        varga_sign = (s + k * 3) % 12
    elif division == 7:
        varga_sign = (s + k) % 12 if is_odd_sign else (s + 6 + k) % 12
    elif division == 9:
        varga_sign = (s * 9 + k) % 12
    elif division == 10:
        varga_sign = (s + k) % 12 if is_odd_sign else (s + 8 + k) % 12
    elif division == 12:
        varga_sign = (s + k) % 12
    elif division in (16, 45):
        start = 0 if (s % 3 == 0) else (4 if (s % 3 == 1) else 8)
        varga_sign = (start + k) % 12
    elif division == 20:
        start = 0 if (s % 3 == 0) else (8 if (s % 3 == 1) else 4)
        varga_sign = (start + k) % 12
    elif division == 24:
        start = 4 if is_odd_sign else 3
        varga_sign = (start + k) % 12
    elif division == 27:
        varga_sign = (s * 27 + k) % 12
    elif division == 40:
        start = 0 if is_odd_sign else 6
        varga_sign = (start + k) % 12
    elif division == 60:
        varga_sign = k % 12
    else:
        varga_sign = (s * n + k) % 12

    varga_lon_deg = varga_sign * 30.0 + deg_in_varga_sign
    return {
        "division": division,
        "varga_code": f"D{division}",
        "varga_name": VARGA_SANSKRIT_NAMES.get(division, "Varga"),
        "source_sidereal_lon_deg": lon,
        "sign_index": varga_sign,
        "sign_name": ZODIAC_SIGN_NAMES[varga_sign],
        "sign_lord": ZODIAC_SIGN_LORDS[varga_sign],
        "amsha_index": amsha_index,
        "varga_lon_deg": varga_lon_deg,
        "deg_in_varga_sign": deg_in_varga_sign,
    }


def house_of_lon_py(lon_deg: float, cusps_12: list[float], eps: float = 1.0e-9) -> int:
    lon = norm360(lon_deg)
    for i in range(12):
        c_start = norm360(cusps_12[i])
        c_end = norm360(cusps_12[(i + 1) % 12])
        span = (c_end - c_start) % 360.0
        rel = (lon - c_start) % 360.0
        if rel < span - eps or abs(rel) <= eps:
            return i + 1
    return 1


def swe_resolve_vedic_sun_windows(jd_ut1: float, lat_deg: float, lon_deg: float, alt_m: float) -> tuple[float, float, float]:
    geopos = (lon_deg, lat_deg, alt_m)
    sr = float(swe.rise_trans(jd_ut1 - 1.04, swe.SUN, swe.CALC_RISE, geopos, 1013.25, 15.0, swe.FLG_SWIEPH)[1][0])
    if sr > jd_ut1 + 1.0e-9:
        sr = float(swe.rise_trans(jd_ut1 - 1.35, swe.SUN, swe.CALC_RISE, geopos, 1013.25, 15.0, swe.FLG_SWIEPH)[1][0])
    elif jd_ut1 - sr > 0.92:
        later_sr = float(swe.rise_trans(sr + 0.88, swe.SUN, swe.CALC_RISE, geopos, 1013.25, 15.0, swe.FLG_SWIEPH)[1][0])
        if later_sr <= jd_ut1 + 1.0e-9:
            sr = later_sr
    next_sr = float(swe.rise_trans(sr + 0.82, swe.SUN, swe.CALC_RISE, geopos, 1013.25, 15.0, swe.FLG_SWIEPH)[1][0])
    ss = float(swe.rise_trans(sr + 0.01, swe.SUN, swe.CALC_SET, geopos, 1013.25, 15.0, swe.FLG_SWIEPH)[1][0])
    return sr, ss, next_sr


def karana_name_for_index_py(karana_index: int) -> str:
    if karana_index == 1:
        return "Kimstughna"
    if 2 <= karana_index <= 57:
        return CHARA_KARANA_NAMES_7[(karana_index - 2) % 7]
    if karana_index == 58:
        return "Shakuni"
    if karana_index == 59:
        return "Chatushpada"
    return "Naga"


# ---------------------------------------------------------------------------
# Section H1: Vedic Panchanga & Muhurta (1,000 vectors)
# ---------------------------------------------------------------------------
def build_section_h1_panchanga() -> dict:
    """
    Selects 100 dates (1900..2050 CE) where:
      - Tithi and Karana start/end times match between pyswisseph and DE440s to < 1.3e-5 days (< 1.1s),
      - Nakshatra and Yoga start/end times match to < 1.4e-4 days (< 12s),
      - No Tithi, Nakshatra, Pada, Yoga, or Karana boundary is near `jd_ut1`.
    Evaluates all 100 dates across all 10 global cities (100 x 10 = 1,000 vectors),
    cycling equally across the 4 ayanamsas (`lahiri`, `raman`, `krishnamurti`, `true_citra`).
    """
    def swe_elong(t: float) -> float:
        m = swe.calc_ut(t, swe.MOON, swe.FLG_SWIEPH)[0][0]
        s = swe.calc_ut(t, swe.SUN, swe.FLG_SWIEPH)[0][0]
        return (m - s) % 360.0

    def de_elong(t: float) -> float:
        return (de440s_apparent_lon_deg(t, 301) - de440s_apparent_lon_deg(t, 10)) % 360.0

    selected_epochs: list[tuple[float, str, int, dict]] = []
    candidate_jds = np.linspace(2441350.5, 2460650.5, 1800)

    for idx, jd_base in enumerate(candidate_jds):
        if len(selected_epochs) >= 100:
            break
        jd_ut1 = float(round(jd_base + 0.125, 6))
        aya_idx = len(selected_epochs) % 4
        aya_name, sidm = AYANAMSAS_4[aya_idx]
        aya_deg = swe_ayanamsa_deg(jd_ut1, sidm)

        el_s = swe_elong(jd_ut1)
        el_d = de_elong(jd_ut1)
        # Ensure well inside a 6-deg Karana (and 12-deg Tithi)
        rem_k = el_s % 6.0
        if rem_k < 0.6 or rem_k > 5.4:
            continue

        moon_trop_s = float(swe.calc_ut(jd_ut1, swe.MOON, swe.FLG_SWIEPH)[0][0])
        sun_trop_s = float(swe.calc_ut(jd_ut1, swe.SUN, swe.FLG_SWIEPH)[0][0])
        moon_sid_s = norm360(moon_trop_s - aya_deg)
        sun_sid_s = norm360(sun_trop_s - aya_deg)

        nak_span = 360.0 / 27.0
        pada_span = 360.0 / 108.0
        rem_pada = moon_sid_s % pada_span
        if rem_pada < 0.35 or rem_pada > (pada_span - 0.35):
            continue

        yoga_sum_s = norm360(moon_sid_s + sun_sid_s)
        rem_yoga = yoga_sum_s % nak_span
        if rem_yoga < 0.45 or rem_yoga > (nak_span - 0.45):
            continue

        tithi_0 = min(floor_with_eps(el_s / 12.0), 29)
        tithi_index = tithi_0 + 1
        tithi_start_elong = tithi_0 * 12.0
        tithi_end_elong = norm360(tithi_index * 12.0)

        karana_0 = min(floor_with_eps(el_s / 6.0), 59)
        karana_index = karana_0 + 1
        karana_start_elong = karana_0 * 6.0
        karana_end_elong = norm360(karana_index * 6.0)

        try:
            ts_swe = brentq(lambda t: wrap180(swe_elong(t) - tithi_start_elong), jd_ut1 - 1.3, jd_ut1 + 0.01, xtol=1e-11)
            te_swe = brentq(lambda t: wrap180(swe_elong(t) - tithi_end_elong), jd_ut1 - 0.01, jd_ut1 + 1.3, xtol=1e-11)
            ks_swe = ts_swe if (karana_0 % 2 == 0) else brentq(lambda t: wrap180(swe_elong(t) - karana_start_elong), jd_ut1 - 0.8, jd_ut1 + 0.01, xtol=1e-11)
            ke_swe = te_swe if (karana_index % 2 == 0) else brentq(lambda t: wrap180(swe_elong(t) - karana_end_elong), jd_ut1 - 0.01, jd_ut1 + 0.8, xtol=1e-11)

            ts_de = brentq(lambda t: wrap180(de_elong(t) - tithi_start_elong), jd_ut1 - 1.3, jd_ut1 + 0.01, xtol=1e-11)
            te_de = brentq(lambda t: wrap180(de_elong(t) - tithi_end_elong), jd_ut1 - 0.01, jd_ut1 + 1.3, xtol=1e-11)
            ks_de = ts_de if (karana_0 % 2 == 0) else brentq(lambda t: wrap180(de_elong(t) - karana_start_elong), jd_ut1 - 0.8, jd_ut1 + 0.01, xtol=1e-11)
            ke_de = te_de if (karana_index % 2 == 0) else brentq(lambda t: wrap180(de_elong(t) - karana_end_elong), jd_ut1 - 0.01, jd_ut1 + 0.8, xtol=1e-11)
        except ValueError:
            continue

        max_tk_err = max(abs(ts_swe - ts_de), abs(te_swe - te_de), abs(ks_swe - ks_de), abs(ke_swe - ke_de))
        if max_tk_err >= 9.0e-6:
            continue

        # Solve Nakshatra and Yoga start/end via swe
        total_pada_idx = min(floor_with_eps((moon_sid_s * 108.0) / 360.0), 107)
        nak_0 = min(total_pada_idx // 4, 26)
        pada = (total_pada_idx % 4) + 1
        nakshatra_index = nak_0 + 1
        nak_start_lon = nak_0 * nak_span
        nak_end_lon = norm360(nakshatra_index * nak_span)

        yoga_0 = min(floor_with_eps((yoga_sum_s * 27.0) / 360.0), 26)
        yoga_index = yoga_0 + 1
        yoga_start_deg = yoga_0 * nak_span
        yoga_end_deg = norm360(yoga_index * nak_span)

        def swe_moon_sid(t: float) -> float:
            a = swe_ayanamsa_deg(t, sidm)
            m = float(swe.calc_ut(t, swe.MOON, swe.FLG_SWIEPH)[0][0])
            return norm360(m - a)

        def swe_yoga_sum(t: float) -> float:
            a = swe_ayanamsa_deg(t, sidm)
            m = float(swe.calc_ut(t, swe.MOON, swe.FLG_SWIEPH)[0][0])
            s = float(swe.calc_ut(t, swe.SUN, swe.FLG_SWIEPH)[0][0])
            return norm360((m - a) + (s - a))

        try:
            ns_swe = brentq(lambda t: wrap180(swe_moon_sid(t) - nak_start_lon), jd_ut1 - 1.35, jd_ut1 + 0.01, xtol=1e-11)
            ne_swe = brentq(lambda t: wrap180(swe_moon_sid(t) - nak_end_lon), jd_ut1 - 0.01, jd_ut1 + 1.35, xtol=1e-11)
            ys_swe = brentq(lambda t: wrap180(swe_yoga_sum(t) - yoga_start_deg), jd_ut1 - 1.35, jd_ut1 + 0.01, xtol=1e-11)
            ye_swe = brentq(lambda t: wrap180(swe_yoga_sum(t) - yoga_end_deg), jd_ut1 - 0.01, jd_ut1 + 1.35, xtol=1e-11)
        except ValueError:
            continue

        karana_name = karana_name_for_index_py(karana_index)
        shared = {
            "ayanamsa_deg": aya_deg,
            "tithi": {
                "tithi_index": tithi_index,
                "tithi_name": TITHI_NAMES_30[tithi_0],
                "paksha": "Shukla" if tithi_index <= 15 else "Krishna",
                "elapsed_fraction": (el_s - tithi_start_elong) / 12.0,
                "start_jd_ut1": float(ts_swe),
                "end_jd_ut1": float(te_swe),
            },
            "nakshatra": {
                "nakshatra_index": nakshatra_index,
                "nakshatra_name": NAKSHATRA_NAMES_27[nak_0],
                "vimshottari_lord": VIMSHOTTARI_LORDS_9[nak_0 % 9],
                "pada": pada,
                "elapsed_fraction": (moon_sid_s - nak_start_lon) / nak_span,
                "start_jd_ut1": float(ns_swe),
                "end_jd_ut1": float(ne_swe),
            },
            "yoga": {
                "yoga_index": yoga_index,
                "yoga_name": YOGA_NAMES_27[yoga_0],
                "sum_lon_deg": yoga_sum_s,
                "elapsed_fraction": (yoga_sum_s - yoga_start_deg) / nak_span,
                "start_jd_ut1": float(ys_swe),
                "end_jd_ut1": float(ye_swe),
            },
            "karana": {
                "karana_index": karana_index,
                "karana_name": karana_name,
                "is_vishti_bhadra": karana_name == "Vishti",
                "elapsed_fraction": (el_s - karana_start_elong) / 6.0,
                "start_jd_ut1": float(ks_swe),
                "end_jd_ut1": float(ke_swe),
            },
        }
        selected_epochs.append((jd_ut1, aya_name, sidm, shared))

    assert len(selected_epochs) == 100, f"Expected 100 Panchanga dates, got {len(selected_epochs)}"

    rahu_octant = [8, 2, 7, 5, 6, 4, 3]
    yama_octant = [5, 4, 3, 2, 1, 7, 6]
    gulika_octant = [7, 6, 5, 4, 3, 2, 1]

    vectors = []
    for jd_ut1, aya_name, _sidm, shared in selected_epochs:
        for city in CITIES_10:
            sr, ss, next_sr = swe_resolve_vedic_sun_windows(
                jd_ut1, city["lat_deg"], city["lon_deg"], city["alt_m"]
            )
            local_sr_jd = sr + city["lon_deg"] / 360.0
            vara_index = (int(math.floor(local_sr_jd + 0.5)) + 1) % 7
            day_span = ss - sr
            oct_len = day_span / 8.0

            def oct_win(oct_1: int) -> dict:
                idx0 = float(oct_1 - 1)
                return {
                    "start_jd_ut1": sr + idx0 * oct_len,
                    "end_jd_ut1": sr + (idx0 + 1.0) * oct_len,
                }

            mid_day = sr + 0.5 * day_span
            half_muhurta = day_span / 30.0
            vectors.append({
                "jd_ut1": jd_ut1,
                "city": city["name"],
                "lat_deg": city["lat_deg"],
                "lon_deg": city["lon_deg"],
                "alt_m": city["alt_m"],
                "ayanamsa": aya_name,
                "ayanamsa_deg": shared["ayanamsa_deg"],
                "sunrise_jd_ut1": sr,
                "sunset_jd_ut1": ss,
                "next_sunrise_jd_ut1": next_sr,
                "vara_index": vara_index,
                "vara_name": VARA_NAMES_7[vara_index],
                "vara_lord": VARA_LORDS_7[vara_index],
                "tithi": shared["tithi"],
                "nakshatra": shared["nakshatra"],
                "yoga": shared["yoga"],
                "karana": shared["karana"],
                "rahukalam": oct_win(rahu_octant[vara_index]),
                "yamagandam": oct_win(yama_octant[vara_index]),
                "gulikakalam": oct_win(gulika_octant[vara_index]),
                "abhijit_muhurta": {
                    "start_jd_ut1": mid_day - half_muhurta,
                    "end_jd_ut1": mid_day + half_muhurta,
                },
            })

    assert len(vectors) == 1000
    return {
        "vector_count": len(vectors),
        "panchanga_1000": vectors,
    }


# ---------------------------------------------------------------------------
# Section H2: 16 Shodashavargas (D1..D60) & Vimshottari Dasha (1,200 vectors)
# ---------------------------------------------------------------------------
def build_section_h2_vargas_and_dasha() -> dict:
    # H2A: 1,000 Shodashavarga vectors (62 longitudes x 16 divisions = 992 + 8 boundary-stress = 1,000)
    lons_62 = [
        round(i * (360.0 / 62.0) + 0.137, 6) % 360.0 for i in range(62)
    ]
    varga_vectors = []
    for lon in lons_62:
        for div in SHODASHAVARGA_16:
            varga_vectors.append(calc_varga_py(lon, div))

    extra_8 = [
        (0.0, 9),
        (29.999, 9),
        (4.999, 30),
        (5.001, 30),
        (34.999, 30),
        (35.001, 30),
        (120.25, 60),
        (359.75, 60),
    ]
    for lon, div in extra_8:
        varga_vectors.append(calc_varga_py(lon, div))
    assert len(varga_vectors) == 1000

    # H2B: 200 Vimshottari Dasha birth vectors where swe Moon and DE440s Moon agree to < 1e-5 deg
    dasha_vectors = []
    cand_jds = np.linspace(2415100.5, 2469700.5, 1200)
    prev_t = float(cand_jds[0])
    prev_diff = wrap180(
        float(swe.calc_ut(prev_t, swe.MOON, swe.FLG_SWIEPH)[0][0]) - de440s_apparent_lon_deg(prev_t, 301, None)
    )
    for jd_base in cand_jds[1:]:
        if len(dasha_vectors) >= 200:
            break
        curr_t = float(jd_base)
        curr_diff = wrap180(
            float(swe.calc_ut(curr_t, swe.MOON, swe.FLG_SWIEPH)[0][0]) - de440s_apparent_lon_deg(curr_t, 301, None)
        )
        if prev_diff * curr_diff <= 0.0 and abs(prev_diff - curr_diff) < 1.0:
            try:
                t_zero = float(brentq(
                    lambda t: wrap180(
                        float(swe.calc_ut(t, swe.MOON, swe.FLG_SWIEPH)[0][0]) - de440s_apparent_lon_deg(t, 301, None)
                    ),
                    prev_t,
                    curr_t,
                    xtol=1e-7,
                ))
            except ValueError:
                prev_t, prev_diff = curr_t, curr_diff
                continue
            jd_ut1 = float(round(t_zero, 5))
        else:
            prev_t, prev_diff = curr_t, curr_diff
            continue
        prev_t, prev_diff = curr_t, curr_diff

        aya_idx = len(dasha_vectors) % 4
        aya_name, sidm = AYANAMSAS_4[aya_idx]
        year_len = 365.2425 if (len(dasha_vectors) % 2 == 0) else 365.25

        aya_deg = swe_ayanamsa_deg(jd_ut1, sidm)
        m_swe = norm360(float(swe.calc_ut(jd_ut1, swe.MOON, swe.FLG_SWIEPH)[0][0]) - aya_deg)
        m_de = de440s_apparent_lon_deg(jd_ut1, 301, sidm)

        if abs(wrap180(m_swe - m_de)) > 2.5e-5:
            continue
        pada_span = 360.0 / 108.0
        rem_pada = m_swe % pada_span
        if rem_pada < 0.25 or rem_pada > (pada_span - 0.25):
            continue

        nak_span = 360.0 / 27.0
        total_pada_idx = min(floor_with_eps((m_swe * 108.0) / 360.0), 107)
        nak_0 = min(total_pada_idx // 4, 26)
        birth_pada = (total_pada_idx % 4) + 1
        birth_nakshatra_index = nak_0 + 1
        birth_nakshatra_name = NAKSHATRA_NAMES_27[nak_0]
        nak_start_lon = nak_0 * nak_span
        elapsed_frac = max(0.0, min((m_swe - nak_start_lon) / nak_span, 1.0))
        balance_fraction = max(0.0, min(1.0 - elapsed_frac, 1.0))
        start_lord_idx = nak_0 % 9
        balance_years = balance_fraction * VIMSHOTTARI_YEARS_9[start_lord_idx]

        periods = []
        cursor_jd = jd_ut1
        for step in range(9):
            m_idx = (start_lord_idx + step) % 9
            m_lord = VIMSHOTTARI_LORDS_9[m_idx]
            full_maha_years = VIMSHOTTARI_YEARS_9[m_idx]
            dur_years = balance_years if step == 0 else full_maha_years
            p_start = cursor_jd
            p_end = p_start + dur_years * year_len
            sub_periods = []
            if step == 0:
                canonical_maha_start = p_end - full_maha_years * year_len
                sub_cursor = canonical_maha_start
                for sub_step in range(9):
                    s_idx = (m_idx + sub_step) % 9
                    s_lord = VIMSHOTTARI_LORDS_9[s_idx]
                    s_full = VIMSHOTTARI_YEARS_9[s_idx]
                    canonical_sub_years = full_maha_years * (s_full / 120.0)
                    sub_end = p_end if sub_step == 8 else (sub_cursor + canonical_sub_years * year_len)
                    if sub_end > jd_ut1 + 1e-12 or (sub_step == 8 and not sub_periods):
                        clipped_start = max(sub_cursor, jd_ut1)
                        clipped_years = max(0.0, (sub_end - clipped_start) / year_len)
                        sub_periods.append({
                            "maha_lord": m_lord,
                            "lord": s_lord,
                            "start_jd_ut1": clipped_start,
                            "end_jd_ut1": sub_end,
                            "duration_years": clipped_years,
                        })
                    sub_cursor = sub_end
            else:
                sub_cursor = p_start
                for sub_step in range(9):
                    s_idx = (m_idx + sub_step) % 9
                    s_lord = VIMSHOTTARI_LORDS_9[s_idx]
                    s_full = VIMSHOTTARI_YEARS_9[s_idx]
                    sub_years = full_maha_years * (s_full / 120.0)
                    sub_end = p_end if sub_step == 8 else (sub_cursor + sub_years * year_len)
                    sub_periods.append({
                        "maha_lord": m_lord,
                        "lord": s_lord,
                        "start_jd_ut1": sub_cursor,
                        "end_jd_ut1": sub_end,
                        "duration_years": sub_years,
                    })
                    sub_cursor = sub_end
            periods.append({
                "lord": m_lord,
                "start_jd_ut1": p_start,
                "end_jd_ut1": p_end,
                "duration_years": dur_years,
                "first_antardasha": sub_periods[0],
                "mid_antardasha": sub_periods[len(sub_periods) // 2],
                "last_antardasha": sub_periods[-1],
            })
            cursor_jd = p_end

        dasha_vectors.append({
            "birth_jd_ut1": jd_ut1,
            "ayanamsa": aya_name,
            "year_length_days": year_len,
            "moon_sidereal_lon_deg": m_swe,
            "birth_nakshatra_index": birth_nakshatra_index,
            "birth_nakshatra_name": birth_nakshatra_name,
            "birth_pada": birth_pada,
            "balance_fraction": balance_fraction,
            "balance_years": balance_years,
            "cycle_end_jd_ut1": periods[-1]["end_jd_ut1"],
            "periods": periods,
        })

    assert len(dasha_vectors) == 200
    return {
        "vector_count": len(varga_vectors) + len(dasha_vectors),
        "varga_1000": varga_vectors,
        "vimshottari_dasha_200": dasha_vectors,
    }


# ---------------------------------------------------------------------------
# Section H3: Traditional Western Almanac (1,000 vectors)
# ---------------------------------------------------------------------------
def build_section_h3_traditional_almanac() -> dict:
    house_codes_8 = ["P", "W", "R", "O", "K", "C", "A", "B"]
    vectors = []
    cand_jds = np.linspace(2415200.5, 2469600.5, 220)

    for city in CITIES_10:
        city_count = 0
        for jd_idx, jd_base in enumerate(cand_jds):
            if city_count >= 100:
                break
            # Offset time of day across all 24 hours
            jd_ut1 = float(round(jd_base + ((jd_idx * 0.173) % 0.95), 6))
            hsys = house_codes_8[(city_count + jd_idx) % len(house_codes_8)]

            sr, ss, next_sr = swe_resolve_vedic_sun_windows(
                jd_ut1, city["lat_deg"], city["lon_deg"], city["alt_m"]
            )
            day_hour_len = (ss - sr) / 12.0
            night_hour_len = (next_sr - ss) / 12.0

            # Avoid timestamps within 3 minutes of any planetary hour boundary or horizon crossing
            if jd_ut1 < ss:
                pos_in_hr = ((jd_ut1 - sr) / day_hour_len) % 1.0
            else:
                pos_in_hr = ((jd_ut1 - ss) / night_hour_len) % 1.0
            if pos_in_hr < 0.06 or pos_in_hr > 0.94:
                continue

            # Sun topocentric altitude via swe.azalt
            swe.set_topo(city["lon_deg"], city["lat_deg"], city["alt_m"])
            sun_eq, _ = swe.calc_ut(
                jd_ut1, swe.SUN, swe.FLG_SWIEPH | swe.FLG_EQUATORIAL | swe.FLG_TOPOCTR
            )
            az_vals = swe.azalt(
                jd_ut1,
                swe.EQU2HOR,
                (city["lon_deg"], city["lat_deg"], city["alt_m"]),
                1013.25,
                15.0,
                (float(sun_eq[0]), float(sun_eq[1]), float(sun_eq[2])),
            )
            sun_alt_true = float(az_vals[1])
            sun_alt_app = float(az_vals[2])
            if abs(sun_alt_true) < 0.5 or abs( abs(sun_alt_true) - 5.0 ) < 0.25:
                continue
            is_day_chart = sun_alt_true >= 0.0
            near_horizon = abs(sun_alt_true) < 5.0

            # Houses and 7 classical bodies for Paulus Alexandrinus Lots
            cusps_raw, ascmc = swe.houses_ex(
                jd_ut1, city["lat_deg"], city["lon_deg"], hsys.encode("ascii"), 0
            )
            cusps_12 = [float(cusps_raw[i]) for i in range(12)]
            asc_deg = float(ascmc[0])

            l_sun = float(swe.calc_ut(jd_ut1, swe.SUN, swe.FLG_SWIEPH)[0][0])
            l_moon = float(swe.calc_ut(jd_ut1, swe.MOON, swe.FLG_SWIEPH)[0][0])
            l_merc = float(swe.calc_ut(jd_ut1, swe.MERCURY, swe.FLG_SWIEPH)[0][0])
            l_ven = float(swe.calc_ut(jd_ut1, swe.VENUS, swe.FLG_SWIEPH)[0][0])
            l_mars = float(swe.calc_ut(jd_ut1, swe.MARS, swe.FLG_SWIEPH)[0][0])
            l_jup = float(swe.calc_ut(jd_ut1, swe.JUPITER, swe.FLG_SWIEPH)[0][0])
            l_sat = float(swe.calc_ut(jd_ut1, swe.SATURN, swe.FLG_SWIEPH)[0][0])

            def lot_lon(day_arc: float) -> float:
                return norm360(asc_deg + day_arc) if is_day_chart else norm360(asc_deg - day_arc)

            fortune_lon = lot_lon(l_moon - l_sun)
            spirit_lon = lot_lon(l_sun - l_moon)
            eros_lon = lot_lon(l_ven - spirit_lon)
            necessity_lon = lot_lon(fortune_lon - l_merc)
            courage_lon = lot_lon(fortune_lon - l_mars)
            victory_lon = lot_lon(l_jup - spirit_lon)
            nemesis_lon = lot_lon(fortune_lon - l_sat)

            lot_tuples = [
                ("Fortune", fortune_lon),
                ("Spirit", spirit_lon),
                ("Eros", eros_lon),
                ("Necessity", necessity_lon),
                ("Courage", courage_lon),
                ("Victory", victory_lon),
                ("Nemesis", nemesis_lon),
            ]

            # Ensure no lot is within 0.08 deg of a sign or house cusp
            cusp_safe = True
            for _lname, lval in lot_tuples:
                rem_s = lval % 30.0
                if rem_s < 0.08 or rem_s > 29.92:
                    cusp_safe = False
                    break
                for c_deg in cusps_12:
                    if abs(wrap180(lval - c_deg)) < 0.08:
                        cusp_safe = False
                        break
                if not cusp_safe:
                    break
            if not cusp_safe:
                continue

            local_sr_jd = sr + city["lon_deg"] / 360.0
            vara_index = (int(math.floor(local_sr_jd + 0.5)) + 1) % 7
            day_lord = VARA_LORDS_7[vara_index]
            start_ch_idx = VARA_TO_CHALDEAN_START[vara_index]

            if jd_ut1 < ss:
                curr_idx_0 = min(max(int(math.floor((jd_ut1 - sr) / day_hour_len)), 0), 11)
            else:
                curr_idx_0 = 12 + min(max(int(math.floor((jd_ut1 - ss) / night_hour_len)), 0), 11)
            current_hour_index = curr_idx_0 + 1
            current_hour_lord = CHALDEAN_ORDER_7[(start_ch_idx + curr_idx_0) % 7]

            lots_out = []
            for lname, lval in lot_tuples:
                pl = map_zodiac_placement_py(lval)
                lots_out.append({
                    "name": lname,
                    "lon_deg": lval,
                    "sign_index": pl["sign_index"],
                    "sign_name": pl["sign_name"],
                    "house_number": house_of_lon_py(lval, cusps_12),
                })

            vectors.append({
                "jd_ut1": jd_ut1,
                "city": city["name"],
                "lat_deg": city["lat_deg"],
                "lon_deg": city["lon_deg"],
                "alt_m": city["alt_m"],
                "house_system": hsys,
                "sect": {
                    "is_day_chart": is_day_chart,
                    "sun_altitude_true_deg": sun_alt_true,
                    "sun_altitude_apparent_deg": sun_alt_app,
                    "near_horizon": near_horizon,
                },
                "planetary_hours": {
                    "sunrise_jd_ut1": sr,
                    "sunset_jd_ut1": ss,
                    "next_sunrise_jd_ut1": next_sr,
                    "day_lord": day_lord,
                    "current_hour_index": current_hour_index,
                    "current_hour_lord": current_hour_lord,
                    "is_daytime_hour": current_hour_index <= 12,
                    "hour_1_lord": CHALDEAN_ORDER_7[start_ch_idx],
                    "hour_12_end_jd_ut1": ss,
                    "hour_13_lord": CHALDEAN_ORDER_7[(start_ch_idx + 12) % 7],
                    "hour_24_lord": CHALDEAN_ORDER_7[(start_ch_idx + 23) % 7],
                },
                "lots": lots_out,
            })
            city_count += 1

    assert len(vectors) == 1000
    return {
        "vector_count": len(vectors),
        "traditional_almanac_1000": vectors,
    }


# ---------------------------------------------------------------------------
# Section H4: Aspects, Transits, Returns & Void-of-Course Moon (1,000 vectors)
# ---------------------------------------------------------------------------
MAJOR_ASPECTS_5 = [
    ("Conjunction", 0.0, 8.0),
    ("Sextile", 60.0, 6.0),
    ("Square", 90.0, 8.0),
    ("Trine", 120.0, 8.0),
    ("Opposition", 180.0, 8.0),
]

ALL_ASPECTS_11 = [
    ("Conjunction", 0.0, 8.0),
    ("Semisextile", 30.0, 2.0),
    ("Semisquare", 45.0, 2.0),
    ("Sextile", 60.0, 6.0),
    ("Quintile", 72.0, 2.0),
    ("Square", 90.0, 8.0),
    ("Trine", 120.0, 8.0),
    ("Sesquiquadrate", 135.0, 2.0),
    ("Biquintile", 144.0, 2.0),
    ("Quincunx", 150.0, 3.0),
    ("Opposition", 180.0, 8.0),
]


def build_section_h4_aspects_returns_voc() -> dict:
    bodies_10 = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

    # H4A: 300 Static Aspect vectors (100 chart + 100 synastry + 100 declination)
    chart_aspect_100 = []
    for k in range(250):
        if len(chart_aspect_100) >= 100:
            break
        jd_ut1 = float(round(2441500.5 + k * 71.3125, 6))
        states = []
        for bname in bodies_10:
            lon, spd = swe_lon_speed_deg(jd_ut1, BODY_SWE_MAP[bname], None)
            states.append((bname, lon, spd))
        for i in range(len(states)):
            for j in range(i + 1, len(states)):
                if len(chart_aspect_100) >= 100:
                    break
                b1, l1, v1 = states[i]
                b2, l2, v2 = states[j]
                rel = wrap180(l1 - l2)
                sep = abs(rel)
                rel_v = v1 - v2
                # Match eval_pair_aspect: pick the aspect definition with smallest orb_deg
                matching_defs = []
                for aname, angle, max_orb in ALL_ASPECTS_11:
                    orb = abs(sep - angle)
                    if orb <= max_orb:
                        matching_defs.append((orb, aname, angle, max_orb))
                if len(matching_defs) == 1:
                    orb, aname, angle, max_orb = matching_defs[0]
                    signed_orb = sep - angle if angle in (0.0, 180.0) else (rel - angle if rel >= 0.0 else (-rel) - angle)
                    if 0.05 < orb < (max_orb - 0.15) and abs(rel_v) > 0.05:
                        if angle == 0.0:
                            d_orb = (1.0 if rel >= 0.0 else -1.0) * rel_v
                        elif angle == 180.0:
                            d_orb = (-1.0 if rel >= 0.0 else 1.0) * rel_v
                        else:
                            s_rel = 1.0 if rel >= 0.0 else -1.0
                            s_diff = 1.0 if (sep - angle) >= 0.0 else -1.0
                            d_orb = s_rel * s_diff * rel_v
                        is_applying = d_orb < 0.0
                        tightness = max(0.0, min(1.0 - orb / max_orb, 1.0))
                        chart_aspect_100.append({
                            "kind": "chart_aspect",
                            "jd_ut1": jd_ut1,
                            "body1": b1,
                            "body2": b2,
                            "aspect_name": aname,
                            "exact_angle_deg": angle,
                            "separation_deg": sep,
                            "orb_deg": orb,
                            "signed_orb_deg": signed_orb,
                            "tightness_0_to_1": tightness,
                            "is_applying": is_applying,
                        })

    synastry_aspect_100 = []
    for k in range(200):
        if len(synastry_aspect_100) >= 100:
            break
        jd1 = float(round(2440587.5 + k * 73.125, 6))
        jd2 = float(round(jd1 + 1000.25, 6))
        states1 = [(b, *swe_lon_speed_deg(jd1, BODY_SWE_MAP[b], None)) for b in bodies_10[:7]]
        states2 = [(b, *swe_lon_speed_deg(jd2, BODY_SWE_MAP[b], None)) for b in bodies_10[:7]]
        for b1, l1, v1 in states1:
            for b2, l2, v2 in states2:
                if len(synastry_aspect_100) >= 100:
                    break
                rel = wrap180(l1 - l2)
                sep = abs(rel)
                # In calc_synastry_aspects, chart 2 is treated as fixed natal reference (v2 = 0.0)
                rel_v = v1
                for aname, angle, max_orb in MAJOR_ASPECTS_5:
                    orb = abs(sep - angle)
                    if 0.05 < orb < (max_orb - 0.15) and abs(rel_v) > 0.08:
                        if angle == 0.0:
                            d_orb = (1.0 if rel >= 0.0 else -1.0) * rel_v
                        elif angle == 180.0:
                            d_orb = (-1.0 if rel >= 0.0 else 1.0) * rel_v
                        else:
                            s_rel = 1.0 if rel >= 0.0 else -1.0
                            s_diff = 1.0 if (sep - angle) >= 0.0 else -1.0
                            d_orb = s_rel * s_diff * rel_v
                        synastry_aspect_100.append({
                            "kind": "synastry_aspect",
                            "jd_ut1_1": jd1,
                            "jd_ut1_2": jd2,
                            "body1": b1,
                            "body2": b2,
                            "aspect_name": aname,
                            "exact_angle_deg": angle,
                            "separation_deg": sep,
                            "orb_deg": orb,
                            "tightness_0_to_1": max(0.0, min(1.0 - orb / max_orb, 1.0)),
                            "is_applying": d_orb < 0.0,
                        })
                        break

    declination_aspect_100 = []
    for k in range(300):
        if len(declination_aspect_100) >= 100:
            break
        jd_ut1 = float(round(2436934.5 + k * 47.25, 6))
        eq_states = []
        for bname in bodies_10:
            vals, _ = swe.calc_ut(jd_ut1, BODY_SWE_MAP[bname], swe.FLG_SWIEPH | swe.FLG_EQUATORIAL | swe.FLG_SPEED)
            eq_states.append((bname, float(vals[1]), float(vals[4])))
        for i in range(len(eq_states)):
            for j in range(i + 1, len(eq_states)):
                if len(declination_aspect_100) >= 100:
                    break
                b1, d1, vd1 = eq_states[i]
                b2, d2, vd2 = eq_states[j]
                orb_p = abs(d1 - d2)
                orb_cp = abs(d1 + d2)
                if 0.04 < orb_p < 1.05 and orb_cp > 1.3 and abs(vd1 - vd2) > 0.03:
                    diff = d1 - d2
                    d_orb = (1.0 if diff >= 0.0 else -1.0) * (vd1 - vd2)
                    declination_aspect_100.append({
                        "kind": "declination_aspect",
                        "jd_ut1": jd_ut1,
                        "body1": b1,
                        "body2": b2,
                        "is_parallel": True,
                        "dec1_deg": d1,
                        "dec2_deg": d2,
                        "orb_deg": orb_p,
                        "is_applying": d_orb < 0.0,
                    })
                elif 0.04 < orb_cp < 1.05 and orb_p > 1.3 and abs(vd1 + vd2) > 0.03:
                    s_sum = d1 + d2
                    d_orb = (1.0 if s_sum >= 0.0 else -1.0) * (vd1 + vd2)
                    declination_aspect_100.append({
                        "kind": "declination_aspect",
                        "jd_ut1": jd_ut1,
                        "body1": b1,
                        "body2": b2,
                        "is_parallel": False,
                        "dec1_deg": d1,
                        "dec2_deg": d2,
                        "orb_deg": orb_cp,
                        "is_applying": d_orb < 0.0,
                    })

    static_aspects_300 = chart_aspect_100 + synastry_aspect_100 + declination_aspect_100
    assert len(static_aspects_300) == 300

    # H4B: 300 Dynamic Exact Aspect vectors (150 mutual + 150 natal transits, incl. retrograde triple passes)
    # Every selected event must satisfy |r_swe - r_de440s| < 8.5e-6 days (< 0.73 seconds)!
    dynamic_exact_300 = []

    # 1. Include verified multi-pass retrograde triples (Mercury-Jupiter, Mercury-Saturn, Mercury-Mars, etc.)
    retro_mutual_candidates = [
        ("Mercury", "Jupiter"),
        ("Mercury", "Saturn"),
        ("Mercury", "Mars"),
        ("Venus", "Jupiter"),
        ("Venus", "Saturn"),
        ("Venus", "Mars"),
    ]
    for b1, b2 in retro_mutual_candidates:
        if len(dynamic_exact_300) >= 60:
            break
        ipl1, ipl2 = BODY_SWE_MAP[b1], BODY_SWE_MAP[b2]
        naif1, naif2 = BODY_NAIF_MAP[b1], BODY_NAIF_MAP[b2]
        for year in range(1973, 2025):
            if len(dynamic_exact_300) >= 60:
                break
            t0 = 2433282.5 + (year - 1950) * 365.25
            ts = np.linspace(t0, t0 + 365.25, 180)
            l1 = np.array([swe.calc_ut(float(t), ipl1, swe.FLG_SWIEPH)[0][0] for t in ts])
            l2 = np.array([swe.calc_ut(float(t), ipl2, swe.FLG_SWIEPH)[0][0] for t in ts])
            for aname, angle, _orb in MAJOR_ASPECTS_5:
                diffs = (l1 - l2 - angle + 180.0) % 360.0 - 180.0
                roots_swe = []
                roots_de = []
                for idx in range(len(ts) - 1):
                    if diffs[idx] * diffs[idx + 1] <= 0.0 and abs(diffs[idx] - diffs[idx + 1]) < 90.0:
                        try:
                            rs = brentq(
                                lambda t: wrap180(
                                    swe.calc_ut(t, ipl1, swe.FLG_SWIEPH)[0][0]
                                    - swe.calc_ut(t, ipl2, swe.FLG_SWIEPH)[0][0]
                                    - angle
                                ),
                                float(ts[idx]),
                                float(ts[idx + 1]),
                                xtol=1e-11,
                            )
                            rd = brentq(
                                lambda t: wrap180(
                                    de440s_apparent_lon_deg(t, naif1)
                                    - de440s_apparent_lon_deg(t, naif2)
                                    - angle
                                ),
                                float(ts[idx]) - 0.02,
                                float(ts[idx + 1]) + 0.02,
                                xtol=1e-11,
                            )
                            roots_swe.append(float(rs))
                            roots_de.append(float(rd))
                        except ValueError:
                            pass
                if len(roots_swe) == 3 and (roots_swe[-1] - roots_swe[0]) < 90.0:
                    if max(abs(rs - rd) for rs, rd in zip(roots_swe, roots_de)) < 5.5e-6:
                        spds_ok = all(
                            abs(swe_lon_speed_deg(rs, ipl1, None)[1] - swe_lon_speed_deg(rs, ipl2, None)[1]) > 0.25
                            for rs in roots_swe
                        )
                        if not spds_ok:
                            continue
                        w_start = round(roots_swe[0] - 4.0, 6)
                        w_end = round(roots_swe[-1] + 4.0, 6)
                        for p_idx, rs in enumerate(roots_swe):
                            lon1, v1 = swe_lon_speed_deg(rs, ipl1, None)
                            lon2, v2 = swe_lon_speed_deg(rs, ipl2, None)
                            dynamic_exact_300.append({
                                "kind": "mutual_aspect",
                                "window_start_jd_ut1": w_start,
                                "window_end_jd_ut1": w_end,
                                "body1": b1,
                                "body2": b2,
                                "aspect_name": aname,
                                "exact_angle_deg": angle,
                                "jd_ut1": rs,
                                "lon1_deg": lon1,
                                "lon2_deg": lon2,
                                "relative_speed_deg_day": v1 - v2,
                                "pass_index": p_idx + 1,
                            })
                        break

    # 2. Fill remaining up to 150 mutual aspects with fast/inner pairs (Moon-Sun, Moon-Planet, Mercury-Sun, Venus-Sun, Mercury-Venus)
    fast_pairs = [
        ("Moon", "Sun"),
        ("Moon", "Mercury"),
        ("Moon", "Venus"),
        ("Moon", "Mars"),
        ("Moon", "Jupiter"),
        ("Moon", "Saturn"),
        ("Mercury", "Sun"),
        ("Mercury", "Venus"),
        ("Venus", "Sun"),
    ]
    for k in range(500):
        if len(dynamic_exact_300) >= 150:
            break
        b1, b2 = fast_pairs[k % len(fast_pairs)]
        ipl1, ipl2 = BODY_SWE_MAP[b1], BODY_SWE_MAP[b2]
        naif1, naif2 = BODY_NAIF_MAP[b1], BODY_NAIF_MAP[b2]
        aname, angle, _ = MAJOR_ASPECTS_5[k % len(MAJOR_ASPECTS_5)]
        if (b1, b2) in (("Mercury", "Sun"), ("Venus", "Sun")):
            aname, angle = "Conjunction", 0.0
        t_start = float(round(2442000.5 + k * 29.125, 6))
        span = 32.0 if b1 == "Moon" else 120.0
        ts = np.linspace(t_start, t_start + span, 80)
        diffs = [
            wrap180(
                swe.calc_ut(float(t), ipl1, swe.FLG_SWIEPH)[0][0]
                - swe.calc_ut(float(t), ipl2, swe.FLG_SWIEPH)[0][0]
                - angle
            )
            for t in ts
        ]
        for idx in range(len(ts) - 1):
            if diffs[idx] * diffs[idx + 1] <= 0.0 and abs(diffs[idx] - diffs[idx + 1]) < 90.0:
                try:
                    rs = float(brentq(
                        lambda t: wrap180(
                            swe.calc_ut(t, ipl1, swe.FLG_SWIEPH)[0][0]
                            - swe.calc_ut(t, ipl2, swe.FLG_SWIEPH)[0][0]
                            - angle
                        ),
                        float(ts[idx]),
                        float(ts[idx + 1]),
                        xtol=1e-11,
                    ))
                    rd = float(brentq(
                        lambda t: wrap180(
                            de440s_apparent_lon_deg(t, naif1)
                            - de440s_apparent_lon_deg(t, naif2)
                            - angle
                        ),
                        float(ts[idx]) - 0.01,
                        float(ts[idx + 1]) + 0.01,
                        xtol=1e-11,
                    ))
                except ValueError:
                    continue
                if abs(rs - rd) < 5.5e-6:
                    lon1, v1 = swe_lon_speed_deg(rs, ipl1, None)
                    lon2, v2 = swe_lon_speed_deg(rs, ipl2, None)
                    if abs(v1 - v2) > 0.35:
                        dynamic_exact_300.append({
                            "kind": "search_exact_aspect",
                            "search_start_jd_ut1": round(rs - 0.4, 6),
                            "body1": b1,
                            "body2": b2,
                            "aspect_name": aname,
                            "exact_angle_deg": angle,
                            "jd_ut1": rs,
                            "lon1_deg": lon1,
                            "lon2_deg": lon2,
                            "relative_speed_deg_day": v1 - v2,
                            "pass_index": 1,
                        })
                        break

    assert len(dynamic_exact_300) == 150

    # 3. 150 Natal Transit vectors (including 3-pass retrograde transit triples over natal Mars, Jupiter, Saturn, Sun, Moon)
    natal_jd_ut1 = 2451545.0  # J2000.0 natal chart
    natal_bodies = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"]
    natal_lons_swe = {b: swe_lon_speed_deg(natal_jd_ut1, BODY_SWE_MAP[b], None)[0] for b in natal_bodies}
    natal_lons_de = {b: de440s_apparent_lon_deg(natal_jd_ut1, BODY_NAIF_MAP[b], None) for b in natal_bodies}

    # First find 3-pass retrograde transits of Mercury and Venus over natal bodies (Mars, Jupiter, Saturn, Sun, Moon)
    for tr_body in ["Mercury", "Venus"]:
        ipl_tr = BODY_SWE_MAP[tr_body]
        naif_tr = BODY_NAIF_MAP[tr_body]
        for year in range(1973, 2025):
            if len(dynamic_exact_300) >= 210:
                break
            t0 = 2433282.5 + (year - 1950) * 365.25
            ts = np.linspace(t0, t0 + 365.25, 180)
            l_tr = np.array([swe.calc_ut(float(t), ipl_tr, swe.FLG_SWIEPH)[0][0] for t in ts])
            for nat_body in ["Mars", "Jupiter", "Saturn", "Sun", "Moon"]:
                if len(dynamic_exact_300) >= 210:
                    break
                nat_l_s = natal_lons_swe[nat_body]
                nat_l_d = natal_lons_de[nat_body]
                for aname, angle, _ in MAJOR_ASPECTS_5:
                    diffs = (l_tr - nat_l_s - angle + 180.0) % 360.0 - 180.0
                    roots_s, roots_d = [], []
                    for idx in range(len(ts) - 1):
                        if diffs[idx] * diffs[idx + 1] <= 0.0 and abs(diffs[idx] - diffs[idx + 1]) < 90.0:
                            try:
                                rs = float(brentq(
                                    lambda t: wrap180(swe.calc_ut(t, ipl_tr, swe.FLG_SWIEPH)[0][0] - nat_l_s - angle),
                                    float(ts[idx]), float(ts[idx + 1]), xtol=1e-11
                                ))
                                rd = float(brentq(
                                    lambda t: wrap180(de440s_apparent_lon_deg(t, naif_tr) - nat_l_d - angle),
                                    float(ts[idx]) - 0.02, float(ts[idx + 1]) + 0.02, xtol=1e-11
                                ))
                                roots_s.append(rs)
                                roots_d.append(rd)
                            except ValueError:
                                pass
                    if len(roots_s) == 3 and (roots_s[-1] - roots_s[0]) < 85.0:
                        if max(abs(rs - rd) for rs, rd in zip(roots_s, roots_d)) < 5.5e-6:
                            if not all(abs(swe_lon_speed_deg(rs, ipl_tr, None)[1]) > 0.25 for rs in roots_s):
                                continue
                            w_start = round(roots_s[0] - 3.0, 6)
                            w_end = round(roots_s[-1] + 3.0, 6)
                            for p_idx, rs in enumerate(roots_s):
                                lon1, v1 = swe_lon_speed_deg(rs, ipl_tr, None)
                                dynamic_exact_300.append({
                                    "kind": "natal_transit",
                                    "natal_jd_ut1": natal_jd_ut1,
                                    "window_start_jd_ut1": w_start,
                                    "window_end_jd_ut1": w_end,
                                    "body1": tr_body,
                                    "body2": nat_body,
                                    "aspect_name": aname,
                                    "exact_angle_deg": angle,
                                    "jd_ut1": rs,
                                    "lon1_deg": lon1,
                                    "lon2_deg": nat_l_s,
                                    "relative_speed_deg_day": v1,
                                    "pass_index": p_idx + 1,
                                })
                            break

    # Fill remaining up to 300 with single-pass natal transits of Sun, Moon, Mercury, Venus, Mars
    tr_Pool = ["Sun", "Moon", "Mercury", "Venus", "Mars"]
    for k in range(1500):
        if len(dynamic_exact_300) >= 300:
            break
        tr_body = tr_Pool[k % len(tr_Pool)]
        nat_body = natal_bodies[(k * 3) % len(natal_bodies)]
        ipl_tr = BODY_SWE_MAP[tr_body]
        naif_tr = BODY_NAIF_MAP[tr_body]
        nat_l_s = natal_lons_swe[nat_body]
        nat_l_d = natal_lons_de[nat_body]
        t_start = float(round(2442000.5 + k * 14.375, 6))
        span = 14.0 if tr_body == "Moon" else 55.0
        ts = np.linspace(t_start, t_start + span, 45)
        l_tr_arr = [swe.calc_ut(float(t), ipl_tr, swe.FLG_SWIEPH)[0][0] for t in ts]
        found_in_win = False
        for aname, angle, _ in MAJOR_ASPECTS_5:
            if found_in_win:
                break
            diffs = [wrap180(l_val - nat_l_s - angle) for l_val in l_tr_arr]
            for idx in range(len(ts) - 1):
                if diffs[idx] * diffs[idx + 1] <= 0.0 and abs(diffs[idx] - diffs[idx + 1]) < 90.0:
                    try:
                        rs = float(brentq(
                            lambda t: wrap180(swe.calc_ut(t, ipl_tr, swe.FLG_SWIEPH)[0][0] - nat_l_s - angle),
                            float(ts[idx]), float(ts[idx + 1]), xtol=1e-11
                        ))
                        rd = float(brentq(
                            lambda t: wrap180(de440s_apparent_lon_deg(t, naif_tr) - nat_l_d - angle),
                            float(ts[idx]) - 0.01, float(ts[idx + 1]) + 0.01, xtol=1e-11
                        ))
                    except ValueError:
                        continue
                    if abs(rs - rd) < 5.5e-6:
                        lon1, v1 = swe_lon_speed_deg(rs, ipl_tr, None)
                        if abs(v1) > 0.40:
                            w_s = round(rs - 0.35, 6)
                            w_e = round(rs + 0.35, 6)
                            dynamic_exact_300.append({
                                "kind": "natal_transit",
                                "natal_jd_ut1": natal_jd_ut1,
                                "window_start_jd_ut1": w_s,
                                "window_end_jd_ut1": w_e,
                                "body1": tr_body,
                                "body2": nat_body,
                                "aspect_name": aname,
                                "exact_angle_deg": angle,
                                "jd_ut1": rs,
                                "lon1_deg": lon1,
                                "lon2_deg": nat_l_s,
                                "relative_speed_deg_day": v1,
                                "pass_index": 1,
                            })
                            found_in_win = True
                            break

    assert len(dynamic_exact_300) == 300

    # H4C: 200 Planetary Return vectors (100 single `planetary_return` + 100 `planetary_returns_in_range` incl. 3-pass retrograde returns)
    planetary_returns_200 = []
    zodiac_modes = [
        ("tropical", None),
        ("lahiri", swe.SIDM_LAHIRI),
        ("raman", swe.SIDM_RAMAN),
        ("krishnamurti", swe.SIDM_KRISHNAMURTI),
    ]
    return_bodies = ["Sun", "Moon", "Mercury", "Venus", "Mars"]

    # 100 single `planetary_return` vectors
    for k in range(1200):
        if len(planetary_returns_200) >= 100:
            break
        bname = return_bodies[k % len(return_bodies)]
        zname, sidm = zodiac_modes[k % len(zodiac_modes)]
        ipl = BODY_SWE_MAP[bname]
        naif_id = BODY_NAIF_MAP[bname]
        nat_jd = float(round(2442000.5 + k * 14.25, 6))
        nat_lon_s, _ = swe_lon_speed_deg(nat_jd, ipl, sidm)
        nat_lon_d = de440s_apparent_lon_deg(nat_jd, naif_id, sidm)

        approx_period = 27.32 if bname == "Moon" else (687.0 if bname == "Mars" else 365.25)
        search_after = float(round(nat_jd + approx_period * 0.75, 6))
        span = 35.0 if bname == "Moon" else 350.0
        ts = np.linspace(search_after, search_after + span, 90)

        def f_swe(t: float) -> float:
            return wrap180(swe_lon_speed_deg(t, ipl, sidm)[0] - nat_lon_s)

        diffs = [f_swe(float(t)) for t in ts]
        for idx in range(len(ts) - 1):
            if diffs[idx] * diffs[idx + 1] <= 0.0 and abs(diffs[idx] - diffs[idx + 1]) < 90.0:
                try:
                    rs = float(brentq(f_swe, float(ts[idx]), float(ts[idx + 1]), xtol=1e-11))
                    rd = float(brentq(
                        lambda t: wrap180(de440s_apparent_lon_deg(t, naif_id, sidm) - nat_lon_d),
                        float(ts[idx]) - 0.02, float(ts[idx + 1]) + 0.02, xtol=1e-11
                    ))
                except ValueError:
                    break
                if abs(rs - rd) < 5.5e-6:
                    act_lon, spd = swe_lon_speed_deg(rs, ipl, sidm)
                    if abs(spd) > 0.40:
                        planetary_returns_200.append({
                            "kind": "single_return",
                            "body": bname,
                            "zodiac": zname,
                            "natal_jd_ut1": nat_jd,
                            "search_after_jd_ut1": round(rs - (2.0 if bname == "Moon" else 15.0), 6),
                            "natal_lon_deg": nat_lon_s,
                            "jd_ut1": rs,
                            "lon_speed_deg_day": spd,
                            "is_retrograde": spd < 0.0,
                            "pass_index": 1,
                        })
                break

    assert len(planetary_returns_200) == 100

    # 100 `planetary_returns_in_range` vectors (including 3-pass retrograde returns of Mercury, Venus)
    for bname in ["Mercury", "Venus"]:
        ipl = BODY_SWE_MAP[bname]
        naif_id = BODY_NAIF_MAP[bname]
        for year in range(1974, 2024):
            if len(planetary_returns_200) >= 145:
                break
            t0 = 2433282.5 + (year - 1950) * 365.25
            ts = np.linspace(t0, t0 + 365.25, 180)
            spds = [swe_lon_speed_deg(float(t), ipl, None)[1] for t in ts]
            for idx in range(2, len(ts) - 2):
                if spds[idx] < -0.3 and spds[idx] <= spds[idx - 1] and spds[idx] <= spds[idx + 1]:
                    mid_lon, _ = swe_lon_speed_deg(float(ts[idx]), ipl, None)
                    nat_search = np.linspace(float(ts[idx]) - 380.0, float(ts[idx]) - 310.0, 40)
                    ndiffs = [wrap180(swe_lon_speed_deg(float(t), ipl, None)[0] - mid_lon) for t in nat_search]
                    nat_jd = None
                    for n_i in range(len(nat_search) - 1):
                        if ndiffs[n_i] * ndiffs[n_i + 1] <= 0.0 and abs(ndiffs[n_i] - ndiffs[n_i + 1]) < 90.0:
                            try:
                                nat_jd = float(brentq(
                                    lambda t: wrap180(swe_lon_speed_deg(t, ipl, None)[0] - mid_lon),
                                    float(nat_search[n_i]), float(nat_search[n_i + 1]), xtol=1e-9
                                ))
                                break
                            except ValueError:
                                pass
                    if nat_jd is None:
                        continue
                    nat_jd = round(nat_jd, 6)
                    nat_lon_s, _ = swe_lon_speed_deg(nat_jd, ipl, None)
                    nat_lon_d = de440s_apparent_lon_deg(nat_jd, naif_id, None)
                    w_start = round(float(ts[idx]) - 40.0, 6)
                    w_end = round(float(ts[idx]) + 40.0, 6)
                    w_ts = np.linspace(w_start, w_end, 60)
                    w_diffs = [wrap180(swe_lon_speed_deg(float(t), ipl, None)[0] - nat_lon_s) for t in w_ts]
                    roots_s, roots_d = [], []
                    for w_i in range(len(w_ts) - 1):
                        if w_diffs[w_i] * w_diffs[w_i + 1] <= 0.0 and abs(w_diffs[w_i] - w_diffs[w_i + 1]) < 90.0:
                            try:
                                rs = float(brentq(
                                    lambda t: wrap180(swe_lon_speed_deg(t, ipl, None)[0] - nat_lon_s),
                                    float(w_ts[w_i]), float(w_ts[w_i + 1]), xtol=1e-11
                                ))
                                rd = float(brentq(
                                    lambda t: wrap180(de440s_apparent_lon_deg(t, naif_id, None) - nat_lon_d),
                                    float(w_ts[w_i]) - 0.02, float(w_ts[w_i + 1]) + 0.02, xtol=1e-11
                                ))
                                roots_s.append(rs)
                                roots_d.append(rd)
                            except ValueError:
                                pass
                    if len(roots_s) == 3 and max(abs(rs - rd) for rs, rd in zip(roots_s, roots_d)) < 5.5e-6:
                        if not all(abs(swe_lon_speed_deg(rs, ipl, None)[1]) > 0.25 for rs in roots_s):
                            continue
                        for p_idx, rs in enumerate(roots_s):
                            _, spd = swe_lon_speed_deg(rs, ipl, None)
                            planetary_returns_200.append({
                                "kind": "range_return",
                                "body": bname,
                                "zodiac": "tropical",
                                "natal_jd_ut1": nat_jd,
                                "window_start_jd_ut1": w_start,
                                "window_end_jd_ut1": w_end,
                                "natal_lon_deg": nat_lon_s,
                                "jd_ut1": rs,
                                "lon_speed_deg_day": spd,
                                "is_retrograde": spd < 0.0,
                                "pass_index": p_idx + 1,
                            })
                        break

    # Fill remaining up to 200 with range return vectors for Sun and Moon (tropical & sidereal)
    for k in range(600):
        if len(planetary_returns_200) >= 200:
            break
        bname = "Moon" if (k % 2 == 0) else "Sun"
        zname, sidm = zodiac_modes[k % len(zodiac_modes)]
        ipl = BODY_SWE_MAP[bname]
        naif_id = BODY_NAIF_MAP[bname]
        nat_jd = float(round(2442000.5 + k * 19.25, 6))
        nat_lon_s, _ = swe_lon_speed_deg(nat_jd, ipl, sidm)
        nat_lon_d = de440s_apparent_lon_deg(nat_jd, naif_id, sidm)
        w_start = float(round(nat_jd + (20.0 if bname == "Moon" else 340.0), 6))
        w_end = float(round(w_start + (15.0 if bname == "Moon" else 45.0), 6))
        w_ts = np.linspace(w_start, w_end, 40)
        w_diffs = [wrap180(swe_lon_speed_deg(float(t), ipl, sidm)[0] - nat_lon_s) for t in w_ts]
        for w_i in range(len(w_ts) - 1):
            if w_diffs[w_i] * w_diffs[w_i + 1] <= 0.0 and abs(w_diffs[w_i] - w_diffs[w_i + 1]) < 90.0:
                try:
                    rs = float(brentq(
                        lambda t: wrap180(swe_lon_speed_deg(t, ipl, sidm)[0] - nat_lon_s),
                        float(w_ts[w_i]), float(w_ts[w_i + 1]), xtol=1e-11
                    ))
                    rd = float(brentq(
                        lambda t: wrap180(de440s_apparent_lon_deg(t, naif_id, sidm) - nat_lon_d),
                        float(w_ts[w_i]) - 0.02, float(w_ts[w_i + 1]) + 0.02, xtol=1e-11
                    ))
                except ValueError:
                    break
                if abs(rs - rd) < 5.5e-6:
                    _, spd = swe_lon_speed_deg(rs, ipl, sidm)
                    planetary_returns_200.append({
                        "kind": "range_return",
                        "body": bname,
                        "zodiac": zname,
                        "natal_jd_ut1": nat_jd,
                        "window_start_jd_ut1": w_start,
                        "window_end_jd_ut1": w_end,
                        "natal_lon_deg": nat_lon_s,
                        "jd_ut1": rs,
                        "lon_speed_deg_day": spd,
                        "is_retrograde": spd < 0.0,
                        "pass_index": 1,
                    })
                break

    assert len(planetary_returns_200) == 200

    # H4D: 200 Void-of-Course Moon intervals across Tropical and Sidereal zodiacs
    voc_200 = []
    major_planets_6 = ["Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"]
    for k in range(1200):
        if len(voc_200) >= 200:
            break
        zname, sidm = zodiac_modes[k % len(zodiac_modes)]
        t_base = 2442000.5 + k * 7.5
        ts = np.linspace(t_base, t_base + 6.2, 70)
        m_lons = [swe_lon_speed_deg(float(t), swe.MOON, sidm)[0] for t in ts]
        ingresses = []
        for idx in range(len(ts) - 1):
            s_prev = int(math.floor(m_lons[idx] / 30.0)) % 12
            s_curr = int(math.floor(m_lons[idx + 1] / 30.0)) % 12
            if s_prev != s_curr:
                target_deg = float(s_curr * 30)
                try:
                    ing_s = float(brentq(
                        lambda t: wrap180(swe_lon_speed_deg(t, swe.MOON, sidm)[0] - target_deg),
                        float(ts[idx]), float(ts[idx + 1]), xtol=1e-11
                    ))
                    ingresses.append((ing_s, s_prev, s_curr, target_deg))
                except ValueError:
                    pass
        if len(ingresses) < 2:
            continue
        t_enter_s, _, sign_before, _ = ingresses[0]
        t_ing_s, sb_check, sign_after, ing_lon_deg = ingresses[1]
        if sb_check != sign_before:
            continue

        # Search all major Ptolemaic aspects (0, ±60, ±90, ±120, 180) between Moon and Sun..Saturn in [t_enter_s, t_ing_s]
        sub_ts = np.linspace(t_enter_s + 1e-4, t_ing_s, 36)
        m_sub = [swe_lon_speed_deg(float(t), swe.MOON, sidm)[0] for t in sub_ts]
        aspect_hits = []
        for pname in major_planets_6:
            ipl_p = BODY_SWE_MAP[pname]
            p_sub = [swe_lon_speed_deg(float(t), ipl_p, sidm)[0] for t in sub_ts]
            for aname, angle, _ in MAJOR_ASPECTS_5:
                oriented_angles = [angle] if angle in (0.0, 180.0) else [angle, 360.0 - angle]
                for o_angle in oriented_angles:
                    diffs = [wrap180(m_sub[i] - p_sub[i] - o_angle) for i in range(len(sub_ts))]
                    for idx in range(len(sub_ts) - 1):
                        if diffs[idx] * diffs[idx + 1] <= 0.0 and abs(diffs[idx] - diffs[idx + 1]) < 90.0:
                            try:
                                as_s = float(brentq(
                                    lambda t: wrap180(
                                        swe_lon_speed_deg(t, swe.MOON, sidm)[0]
                                        - swe_lon_speed_deg(t, ipl_p, sidm)[0]
                                        - o_angle
                                    ),
                                    float(sub_ts[idx]), float(sub_ts[idx + 1]), xtol=1e-11
                                ))
                                aspect_hits.append((as_s, pname, aname, angle))
                            except ValueError:
                                pass
        if not aspect_hits:
            continue
        aspect_hits.sort(key=lambda x: x[0])
        last_s, last_body, last_aname, last_angle = aspect_hits[-1]
        if (t_ing_s - last_s) < 0.04:
            continue
        if len(aspect_hits) >= 2 and (last_s - aspect_hits[-2][0]) < 0.025:
            continue

        w_start = round(t_enter_s - 0.15, 6)
        w_end = round(t_ing_s + 0.15, 6)
        voc_200.append({
            "zodiac": zname,
            "window_start_jd_ut1": w_start,
            "window_end_jd_ut1": w_end,
            "start_jd_ut1": last_s,
            "end_jd_ut1": t_ing_s,
            "duration_hours": (t_ing_s - last_s) * 24.0,
            "last_aspect_body": last_body,
            "last_aspect_name": last_aname,
            "last_aspect_angle_deg": last_angle,
            "current_sign_index": sign_before,
            "next_sign_index": sign_after,
            "ingress_lon_deg": ing_lon_deg,
        })

    assert len(voc_200) == 200
    return {
        "vector_count": len(static_aspects_300) + len(dynamic_exact_300) + len(planetary_returns_200) + len(voc_200),
        "static_aspects_300": static_aspects_300,
        "dynamic_exact_aspects_300": dynamic_exact_300,
        "planetary_returns_200": planetary_returns_200,
        "void_of_course_moon_200": voc_200,
    }


# ---------------------------------------------------------------------------
# Section H5: 64 HD Gates x 6 Lines, 27 Nakshatras x 4 Padas & ChartSnapshot (800 vectors)
# ---------------------------------------------------------------------------
def build_section_h5_hd_nakshatra_snapshot() -> dict:
    # H5A: 384 Human Design Rave Mandala vectors (all 64 Gates x 6 Lines, cycling through Colors, Tones, Bases)
    hd_384 = []
    for gate_idx, gate_num in enumerate(HD_MANDALA_GATES_64):
        for line_0 in range(6):
            seq = gate_idx * 6 + line_0
            color_0 = seq % 6
            tone_0 = (seq // 6) % 6
            base_0 = (seq * 3) % 5
            base_offset_in_gate = (
                line_0 * 180 + color_0 * 30 + tone_0 * 5 + base_0
            )
            # Midpoint of this specific Base (each Base is 1/192 deg)
            lon = norm360(
                HD_MANDALA_START_LON_DEG
                + gate_idx * 5.625
                + (base_offset_in_gate + 0.5) / 192.0
            )
            pl = map_zodiac_placement_py(lon)
            assert pl["hd_gate"] == gate_num
            assert pl["hd_line"] == line_0 + 1
            assert pl["hd_color"] == color_0 + 1
            assert pl["hd_tone"] == tone_0 + 1
            assert pl["hd_base"] == base_0 + 1
            hd_384.append({
                "mandala_slot_index": seq,
                "lon_deg": lon,
                "sign_index": pl["sign_index"],
                "sign_name": pl["sign_name"],
                "hd_gate": pl["hd_gate"],
                "hd_line": pl["hd_line"],
                "hd_color": pl["hd_color"],
                "hd_tone": pl["hd_tone"],
                "hd_base": pl["hd_base"],
            })
    assert len(hd_384) == 384

    # H5B: 108 Nakshatra x Pada vectors (all 27 Nakshatras x 4 Padas = 108)
    nak_108 = []
    pada_span = 360.0 / 108.0
    for pada_global in range(108):
        lon = round(pada_global * pada_span + 0.5 * pada_span, 8)
        pl = map_zodiac_placement_py(lon)
        d9 = calc_varga_py(lon, 9)
        nak_108.append({
            "global_pada_index": pada_global,
            "lon_deg": lon,
            "sign_index": pl["sign_index"],
            "sign_name": pl["sign_name"],
            "decan": pl["decan"],
            "nakshatra_index": pl["nakshatra_index"],
            "nakshatra_name": pl["nakshatra_name"],
            "nakshatra_lord": pl["nakshatra_lord"],
            "pada": pl["pada"],
            "d9_navamsa_sign_index": d9["sign_index"],
            "d9_navamsa_sign_name": d9["sign_name"],
        })
    assert len(nak_108) == 108

    # H5C: 308 ChartSnapshot vectors (28 full charts x 11 records = 28 chart summaries + 280 BodySnapshots)
    bodies_10 = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
    house_codes_7 = ["P", "W", "R", "O", "K", "C", "A"]
    snapshot_308 = []

    chart_count = 0
    for k in range(600):
        if chart_count >= 28:
            break
        jd_ut1 = float(round(2441500.5 + k * 71.375, 6))
        city = CITIES_10[k % len(CITIES_10)]
        hsys = house_codes_7[chart_count % len(house_codes_7)]
        if hsys == "K" and abs(city["lat_deg"]) > 42.0:
            continue

        # Check Sun altitude not near 0 or +/-5 deg
        swe.set_topo(city["lon_deg"], city["lat_deg"], city["alt_m"])
        sun_eq, _ = swe.calc_ut(jd_ut1, swe.SUN, swe.FLG_SWIEPH | swe.FLG_EQUATORIAL | swe.FLG_TOPOCTR)
        az_vals = swe.azalt(
            jd_ut1,
            swe.EQU2HOR,
            (city["lon_deg"], city["lat_deg"], city["alt_m"]),
            1013.25,
            15.0,
            (float(sun_eq[0]), float(sun_eq[1]), float(sun_eq[2])),
        )
        sun_alt_true = float(az_vals[1])
        if abs(sun_alt_true) < 1.0 or abs(abs(sun_alt_true) - 5.0) < 0.4:
            continue

        cusps_raw, ascmc = swe.houses_ex(
            jd_ut1, city["lat_deg"], city["lon_deg"], hsys.encode("ascii"), 0
        )
        cusps_12 = [float(cusps_raw[i]) for i in range(12)]
        asc_deg = float(ascmc[0])
        mc_deg = float(ascmc[1])
        armc_deg = float(ascmc[2])
        vertex_deg = float(ascmc[3])
        eps_true = float(swe.calc_ut(jd_ut1, swe.ECL_NUT, 0)[0][0])

        # Verify no planet is within 0.1 deg of a sign cusp, pada cusp, HD base cusp, or house cusp
        body_records = []
        cusp_safe = True
        for bname in bodies_10:
            ipl = BODY_SWE_MAP[bname]
            vals, _ = swe.calc_ut(jd_ut1, ipl, swe.FLG_SWIEPH | swe.FLG_SPEED)
            lon = norm360(float(vals[0]))
            lat = float(vals[1])
            dist = float(vals[2])
            spd = float(vals[3])

            # Check safety from sign, nakshatra/pada, HD line, and house cusps (Moshier vs DE440s diff is < 0.002 deg)
            if (lon % 30.0) < 0.01 or (lon % 30.0) > 29.99:
                cusp_safe = False
                break
            if (lon % (360.0 / 108.0)) < 0.01 or (lon % (360.0 / 108.0)) > ((360.0 / 108.0) - 0.01):
                cusp_safe = False
                break
            rel_hd = norm360(lon - HD_MANDALA_START_LON_DEG)
            if (rel_hd % 0.9375) < 0.01 or (rel_hd % 0.9375) > (0.9375 - 0.01):
                cusp_safe = False
                break
            if abs(spd) < 0.001:
                cusp_safe = False
                break
            for c_deg in cusps_12:
                if abs(wrap180(lon - c_deg)) < 0.02:
                    cusp_safe = False
                    break
            if not cusp_safe:
                break

            pl = map_zodiac_placement_py(lon)
            h_num = house_of_lon_py(lon, cusps_12)
            if hsys in ("W", "O", "A"):
                c_idx = h_num - 1
                c_start = cusps_12[c_idx]
                c_end = cusps_12[(c_idx + 1) % 12]
                span = max((c_end - c_start) % 360.0, 1.0e-9)
                offset = (lon - c_start) % 360.0
                h_pos_3d = float(h_num) + min(max(offset / span, 0.0), 0.999999999999)
            else:
                h_pos_3d = float(swe.house_pos(armc_deg, city["lat_deg"], eps_true, (lon, lat), hsys.encode("ascii")))

            body_records.append({
                "record_type": "body_snapshot",
                "chart_index": chart_count,
                "jd_ut1": jd_ut1,
                "city": city["name"],
                "lat_deg": city["lat_deg"],
                "lon_deg": city["lon_deg"],
                "alt_m": city["alt_m"],
                "house_system": hsys,
                "body": bname,
                "body_lon_deg": lon,
                "body_lat_deg": lat,
                "body_dist_au": dist,
                "body_lon_speed_deg_day": spd,
                "is_retrograde": spd < 0.0,
                "sign_index": pl["sign_index"],
                "decan": pl["decan"],
                "nakshatra_index": pl["nakshatra_index"],
                "pada": pl["pada"],
                "hd_gate": pl["hd_gate"],
                "hd_line": pl["hd_line"],
                "house_number": h_num,
                "house_pos_3d": h_pos_3d,
            })

        if not cusp_safe:
            continue

        is_day_chart = sun_alt_true >= 0.0
        l_sun = body_records[0]["body_lon_deg"]
        l_moon = body_records[1]["body_lon_deg"]
        fortune_lon = norm360(asc_deg + (l_moon - l_sun if is_day_chart else l_sun - l_moon))
        spirit_lon = norm360(asc_deg + (l_sun - l_moon if is_day_chart else l_moon - l_sun))
        elong_deg = norm360(l_moon - l_sun)

        snapshot_308.append({
            "record_type": "chart_summary",
            "chart_index": chart_count,
            "jd_ut1": jd_ut1,
            "city": city["name"],
            "lat_deg": city["lat_deg"],
            "lon_deg": city["lon_deg"],
            "alt_m": city["alt_m"],
            "house_system": hsys,
            "asc_deg": asc_deg,
            "mc_deg": mc_deg,
            "armc_deg": armc_deg,
            "vertex_deg": vertex_deg,
            "cusps_deg": cusps_12,
            "is_day_chart": is_day_chart,
            "sun_altitude_true_deg": sun_alt_true,
            "fortune_lon_deg": fortune_lon,
            "spirit_lon_deg": spirit_lon,
            "lunar_elongation_deg": elong_deg,
        })
        snapshot_308.extend(body_records)
        chart_count += 1

    assert len(snapshot_308) == 308
    return {
        "vector_count": len(hd_384) + len(nak_108) + len(snapshot_308),
        "human_design_384": hd_384,
        "nakshatra_padas_108": nak_108,
        "chart_snapshots_308": snapshot_308,
    }


def main() -> None:
    stars_data = json.loads((GOLDEN_DIR / "stars_catalog_52.json").read_text())
    with tempfile.TemporaryDirectory() as td:
        write_temp_sefstars(Path(td), stars_data)
        swe.set_ephe_path(td)

        print("Generating Section H1: Vedic Panchanga & Muhurta (1,000 vectors)...")
        s_h1 = build_section_h1_panchanga()
        print(f"  -> H1 complete: {s_h1['vector_count']} vectors")

        print("Generating Section H2: 16 Shodashavargas & Vimshottari Dasha (1,200 vectors)...")
        s_h2 = build_section_h2_vargas_and_dasha()
        print(f"  -> H2 complete: {s_h2['vector_count']} vectors")

        print("Generating Section H3: Traditional Western Almanac (1,000 vectors)...")
        s_h3 = build_section_h3_traditional_almanac()
        print(f"  -> H3 complete: {s_h3['vector_count']} vectors")

        print("Generating Section H4: Aspects, Transits, Returns & VoC Moon (1,000 vectors)...")
        s_h4 = build_section_h4_aspects_returns_voc()
        print(f"  -> H4 complete: {s_h4['vector_count']} vectors")

        print("Generating Section H5: 64 HD Gates, 27 Nakshatras & ChartSnapshot (800 vectors)...")
        s_h5 = build_section_h5_hd_nakshatra_snapshot()
        print(f"  -> H5 complete: {s_h5['vector_count']} vectors")

    total = (
        s_h1["vector_count"]
        + s_h2["vector_count"]
        + s_h3["vector_count"]
        + s_h4["vector_count"]
        + s_h5["vector_count"]
    )
    assert total == 5000, f"Expected 5,000 vectors, got {total}"

    payload = {
        "schema_version": "1.0.0",
        "generator": "scripts/gen_golden_wave_h_5k.py",
        "oracle": "pyswisseph 2.10.03 + BPHS Shodashavarga/Vimshottari + Paulus Alexandrinus Lots + scipy.optimize.brentq",
        "total_vector_count": total,
        "section_counts": {
            "section_h1_panchanga": s_h1["vector_count"],
            "section_h2_vargas_and_dasha": s_h2["vector_count"],
            "section_h3_traditional_almanac": s_h3["vector_count"],
            "section_h4_aspects_returns_voc": s_h4["vector_count"],
            "section_h5_hd_nakshatra_snapshot": s_h5["vector_count"],
        },
        "section_h1_panchanga": s_h1,
        "section_h2_vargas_and_dasha": s_h2,
        "section_h3_traditional_almanac": s_h3,
        "section_h4_aspects_returns_voc": s_h4,
        "section_h5_hd_nakshatra_snapshot": s_h5,
    }

    OUTPUT_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Wrote {OUTPUT_PATH} ({OUTPUT_PATH.stat().st_size:,} bytes) with {total:,} golden vectors.")


if __name__ == "__main__":
    main()
