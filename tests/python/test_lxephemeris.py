#!/usr/bin/env python3
"""Comprehensive test suite for `lxephemeris` Python bindings (`PyO3`).

Tests:
1. Moshier engine creation (`new_moshier`), mode check (`is_moshier_mode`), and analytical evaluation (`calc_moshier`).
2. Fallback loading (`from_data_dir_or_moshier`).
3. Custom stars runtime loading (`load_custom_stars`).
4. Planetary/lunar phenomena (`phenomena`) and Swiss Ephemeris 20-attribute formatting (`to_swe_attr20`).
5. Orbital nodes and apsides (`nod_aps`) for Moon and major planets across mean/osculating bitmasks (`SE_NODBIT_*`).
6. Lunar occultations (`next_lunar_occultation`, `next_local_lunar_occultation`, `lunar_occultations_in_range`).
7. 3D House position (`house_pos`), Equation of Time (`time_equation`), and Inverse Horizontal transformation (`azalt_rev`).
8. `SweAdapter` extensions: `new_moshier`, `set_sid_mode` (`SE_SIDM_USER` with `(t0, ayan_t0)`), `pheno_ut`, `pheno`,
   `nod_aps_ut`, `nod_aps`, `lun_occult_when_glob`, `lun_occult_when_loc`, `house_pos`, `time_equ`, `azalt_rev`.
9. Constants: `SE_NODBIT_*`, `SE_HOR2*`, all 43 `SE_SIDM_*` + `SE_SIDM_USER`.
"""

from __future__ import annotations

import math
import os
import pathlib
import sys
import tempfile
import unittest

WORKSPACE_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
RELEASE_DIR = WORKSPACE_ROOT / "target" / "release"
DATA_DIR = WORKSPACE_ROOT / "data"

sys.path.insert(0, str(WORKSPACE_ROOT))
sys.path.insert(0, str(RELEASE_DIR))
import lxephemeris  # type: ignore[import-not-found]

JD_J2000 = 2_451_545.0
JD_2015_ALDEBARAN = 2_457_023.5


class TestLxEphemerisPy(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.kernel_engine = lxephemeris.Engine.from_data_dir(
            str(DATA_DIR),
            max_threads=2,
            pin_time_table_version="embedded",
            verify_manifest_sha256=False,
        )

    def test_01_moshier_engine_and_calc_moshier(self) -> None:
        """Verify fileless Moshier engine creation and calc_moshier."""
        mosh_eng = lxephemeris.Engine.new_moshier(max_threads=2)
        self.assertTrue(mosh_eng.is_moshier_mode())

        # Test calc_moshier for Sun at J2000
        sun_mosh = mosh_eng.calc_moshier("Sun", JD_J2000, 0.0, "TT")
        self.assertAlmostEqual(sun_mosh.lon_deg, 280.368, delta=0.1)

        # Standard calc() delegates to Moshier in pure analytical mode
        sun_calc = mosh_eng.calc("Sun", JD_J2000, 0.0, "TT")
        self.assertAlmostEqual(sun_calc.lon_deg, sun_mosh.lon_deg, delta=1e-7)

    def test_02_from_data_dir_or_moshier_fallback(self) -> None:
        """Verify fallback to Moshier when data dir does not exist or has no kernels."""
        fallback = lxephemeris.Engine.from_data_dir_or_moshier("/nonexistent_test_directory_xyz")
        self.assertTrue(fallback.is_moshier_mode())

        loaded = lxephemeris.Engine.from_data_dir_or_moshier(str(DATA_DIR))
        self.assertFalse(loaded.is_moshier_mode())

    def test_03_load_custom_stars(self) -> None:
        """Verify dynamic loading of custom stars from sefstars format."""
        mosh_eng = lxephemeris.Engine.new_moshier()
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write("PyTestStar,alp Test,1991,12,0,0,10,0,0,0,0,0,10,2.5,test\n")
            tmp_path = f.name

        try:
            count = mosh_eng.load_custom_stars(tmp_path)
            self.assertGreaterEqual(count, 1)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_04_phenomena(self) -> None:
        """Verify phenomena calculation for major planets and attr20 formatting."""
        ph = self.kernel_engine.phenomena("Mars", JD_J2000, 0.0, "TT")
        self.assertGreaterEqual(ph.phase_angle_deg, 0.0)
        self.assertLessEqual(ph.phase_angle_deg, 180.0)
        self.assertGreaterEqual(ph.phase_illuminated, 0.0)
        self.assertLessEqual(ph.phase_illuminated, 1.0)
        self.assertGreater(ph.geocentric_distance_au, 0.0)
        self.assertGreater(ph.heliocentric_distance_au, 0.0)

        attr20 = ph.to_swe_attr20()
        self.assertEqual(len(attr20), 20)
        self.assertAlmostEqual(attr20[0], ph.phase_angle_deg, delta=1e-9)
        self.assertAlmostEqual(attr20[1], ph.phase_illuminated, delta=1e-9)
        self.assertAlmostEqual(attr20[2], ph.elongation_deg, delta=1e-9)
        self.assertAlmostEqual(attr20[3], ph.apparent_diameter_deg, delta=1e-9)
        self.assertAlmostEqual(attr20[4], ph.apparent_magnitude_v, delta=1e-9)
        self.assertAlmostEqual(attr20[5], ph.horizontal_parallax_deg, delta=1e-9)

    def test_05_nod_aps(self) -> None:
        """Verify nodes and apsides calculation for Moon and major planets."""
        # Moon mean node
        nod_mean = self.kernel_engine.nod_aps("Moon", JD_J2000, 0.0, "TT", method=lxephemeris.SE_NODBIT_MEAN)
        self.assertEqual(len(nod_mean.ascending_node), 6)
        self.assertEqual(len(nod_mean.descending_node), 6)
        self.assertEqual(len(nod_mean.periapsis), 6)
        self.assertEqual(len(nod_mean.apoapsis), 6)
        self.assertEqual(len(nod_mean.asc_node), 6)
        self.assertEqual(len(nod_mean.desc_node), 6)
        self.assertGreaterEqual(nod_mean.ascending_node[0], 0.0)
        self.assertLessEqual(nod_mean.ascending_node[0], 360.0)

        # Jupiter osculating node
        nod_oscu = self.kernel_engine.nod_aps("Jupiter", JD_J2000, 0.0, "TT", method=lxephemeris.SE_NODBIT_OSCU)
        self.assertGreaterEqual(nod_oscu.ascending_node[0], 0.0)
        self.assertLessEqual(nod_oscu.ascending_node[0], 360.0)

    def test_06_lunar_occultations(self) -> None:
        """Verify lunar occultation search forward and locally."""
        occ = self.kernel_engine.next_lunar_occultation("Aldebaran", JD_2015_ALDEBARAN, max_synodic_months=12)
        self.assertIn(occ.kind, ["Total", "Grazing", "Partial", "Annular"])
        self.assertGreater(occ.jd_ut1_max, JD_2015_ALDEBARAN)
        self.assertEqual(len(occ.to_swe_glob_tret()), 10)
        self.assertEqual(len(occ.to_swe_tret10()), 10)
        self.assertGreaterEqual(occ.swe_retflag(), 0)

        # Local occultation from Greenwich
        local_occ = self.kernel_engine.next_local_lunar_occultation(
            "Aldebaran",
            JD_2015_ALDEBARAN,
            latitude_deg=51.4769,
            longitude_deg=0.0,
            max_synodic_months=12,
        )
        self.assertGreater(local_occ.jd_ut1_max, JD_2015_ALDEBARAN)

        # Occultations in range
        occ_list = self.kernel_engine.lunar_occultations_in_range(
            "Aldebaran",
            JD_2015_ALDEBARAN,
            JD_2015_ALDEBARAN + 60.0,
        )
        self.assertGreaterEqual(len(occ_list), 1)

    def test_07_house_pos_time_equation_azalt_rev(self) -> None:
        """Verify 3D house_pos, time_equation, and azalt_rev on Engine."""
        hp = self.kernel_engine.house_pos(180.0, 51.5074, 23.439, 280.368, 0.0, "P")
        self.assertGreaterEqual(hp, 1.0)
        self.assertLess(hp, 13.0)

        teq = self.kernel_engine.time_equation(JD_J2000)
        self.assertLess(abs(teq), 60.0)

        coords = self.kernel_engine.azalt_rev(
            JD_J2000,
            lxephemeris.SE_HOR2EQU,
            51.5074,
            -0.1278,
            180.0,
            45.0,
        )
        self.assertEqual(len(coords), 2)
        self.assertGreaterEqual(coords[0], 0.0)
        self.assertLessEqual(coords[0], 360.0)
        self.assertGreaterEqual(coords[1], -90.0)
        self.assertLessEqual(coords[1], 90.0)

    def test_08_swe_adapter_extensions(self) -> None:
        """Verify SweAdapter extensions including Moshier, USER sidereal mode, pheno, nod_aps, occultations."""
        swe_mosh = lxephemeris.SweAdapter.new_moshier()
        sun_arr, flag = swe_mosh.calc(JD_J2000, lxephemeris.SE_SUN, lxephemeris.SEFLG_SPEED)
        self.assertAlmostEqual(sun_arr[0], 280.368, delta=0.1)
        self.assertTrue(bool(flag & lxephemeris.SEFLG_SPEED))

        # USER sidereal mode with custom (t0, ayan_t0)
        swe_mosh.set_sid_mode(lxephemeris.SE_SIDM_USER, 2451545.0, 24.0)
        swe_mosh.set_sid_mode(lxephemeris.SE_SIDM_LAHIRI)

        # pheno_ut and pheno
        attr20_ut, retflag_ut = swe_mosh.pheno_ut(JD_J2000, lxephemeris.SE_MARS, 0)
        self.assertEqual(len(attr20_ut), 20)
        self.assertEqual(retflag_ut, lxephemeris.SEFLG_MOSEPH)

        attr20_tt, retflag_tt = swe_mosh.pheno(JD_J2000, lxephemeris.SE_MARS, 0)
        self.assertEqual(len(attr20_tt), 20)

        # nod_aps_ut and nod_aps
        asc, desc, peri, apo, ret_nod = swe_mosh.nod_aps_ut(
            JD_J2000, lxephemeris.SE_MOON, 258, lxephemeris.SE_NODBIT_MEAN
        )
        self.assertEqual(len(asc), 6)
        self.assertEqual(len(desc), 6)
        self.assertEqual(len(peri), 6)
        self.assertEqual(len(apo), 6)

        asc2, desc2, peri2, apo2, _ = swe_mosh.nod_aps(
            JD_J2000, lxephemeris.SE_MOON, 258, lxephemeris.SE_NODBIT_MEAN
        )
        self.assertEqual(len(asc2), 6)

        # house_pos, time_equ, azalt_rev
        shp = swe_mosh.house_pos(180.0, 51.5074, 23.439, [280.368, 0.0], "P")
        self.assertGreaterEqual(shp, 1.0)
        self.assertLess(shp, 13.0)

        ste = swe_mosh.time_equ(JD_J2000)
        self.assertLess(abs(ste), 0.1)  # days

        s_coords = swe_mosh.azalt_rev(
            JD_J2000,
            lxephemeris.SE_HOR2EQU,
            [-0.1278, 51.5074, 0.0],
            [180.0, 45.0],
        )
        self.assertEqual(len(s_coords), 2)
        self.assertGreaterEqual(s_coords[0], 0.0)
        self.assertLessEqual(s_coords[0], 360.0)

        # SweAdapter occultations with kernel engine
        swe_kernel = lxephemeris.SweAdapter(self.kernel_engine)
        ret_glob, tret_glob = swe_kernel.lun_occult_when_glob(
            JD_2015_ALDEBARAN, 0, starname="Aldebaran", iflag=258, backward=False
        )
        self.assertEqual(len(tret_glob), 10)
        self.assertGreaterEqual(ret_glob, 0)

        ret_loc, tret_loc, attr_loc = swe_kernel.lun_occult_when_loc(
            JD_2015_ALDEBARAN,
            0,
            geopos=[0.0, 51.5, 0.0],
            starname="Aldebaran",
            iflag=258,
            backward=False,
        )
        self.assertEqual(len(tret_loc), 10)
        self.assertEqual(len(attr_loc), 20)
        self.assertGreaterEqual(ret_loc, 0)

    def test_09_constants(self) -> None:
        """Verify presence and exact values of all required constants."""
        self.assertEqual(lxephemeris.SE_NODBIT_MEAN, 1)
        self.assertEqual(lxephemeris.SE_NODBIT_OSCU, 2)
        self.assertEqual(lxephemeris.SE_NODBIT_OSCU_BAR, 4)
        self.assertEqual(lxephemeris.SE_NODBIT_FOPOINT, 256)

        self.assertEqual(lxephemeris.SE_HOR2ECL, 0)
        self.assertEqual(lxephemeris.SE_HOR2EQU, 1)

        self.assertEqual(lxephemeris.SE_SIDM_FAGAN_BRADLEY, 0)
        self.assertEqual(lxephemeris.SE_SIDM_LAHIRI, 1)
        self.assertEqual(lxephemeris.SE_SIDM_DELUCE, 2)
        self.assertEqual(lxephemeris.SE_SIDM_RAMAN, 3)
        self.assertEqual(lxephemeris.SE_SIDM_USHASHASHI, 4)
        self.assertEqual(lxephemeris.SE_SIDM_KRISHNAMURTI, 5)
        self.assertEqual(lxephemeris.SE_SIDM_DJWHAL_KHUL, 6)
        self.assertEqual(lxephemeris.SE_SIDM_YUKTESHWAR, 7)
        self.assertEqual(lxephemeris.SE_SIDM_JN_BHASIN, 8)
        self.assertEqual(lxephemeris.SE_SIDM_BABYL_KUGLER1, 9)
        self.assertEqual(lxephemeris.SE_SIDM_BABYL_KUGLER2, 10)
        self.assertEqual(lxephemeris.SE_SIDM_BABYL_KUGLER3, 11)
        self.assertEqual(lxephemeris.SE_SIDM_BABYL_HUBER, 12)
        self.assertEqual(lxephemeris.SE_SIDM_BABYL_ETPSC, 13)
        self.assertEqual(lxephemeris.SE_SIDM_ALDEBARAN_15TAU, 14)
        self.assertEqual(lxephemeris.SE_SIDM_HIPPARCHOS, 15)
        self.assertEqual(lxephemeris.SE_SIDM_SASSANIAN, 16)
        self.assertEqual(lxephemeris.SE_SIDM_GALCENT_0SAG, 17)
        self.assertEqual(lxephemeris.SE_SIDM_J2000, 18)
        self.assertEqual(lxephemeris.SE_SIDM_J1900, 19)
        self.assertEqual(lxephemeris.SE_SIDM_B1950, 20)
        self.assertEqual(lxephemeris.SE_SIDM_SURYASIDDHANTA, 21)
        self.assertEqual(lxephemeris.SE_SIDM_SURYASIDDHANTA_MSUN, 22)
        self.assertEqual(lxephemeris.SE_SIDM_ARYABHATA, 23)
        self.assertEqual(lxephemeris.SE_SIDM_ARYABHATA_MSUN, 24)
        self.assertEqual(lxephemeris.SE_SIDM_SS_REVATI, 25)
        self.assertEqual(lxephemeris.SE_SIDM_SS_CITRA, 26)
        self.assertEqual(lxephemeris.SE_SIDM_TRUE_CITRA, 27)
        self.assertEqual(lxephemeris.SE_SIDM_TRUE_REVATI, 28)
        self.assertEqual(lxephemeris.SE_SIDM_TRUE_PUSHYA, 29)
        self.assertEqual(lxephemeris.SE_SIDM_GALCENT_RGILBRAND, 30)
        self.assertEqual(lxephemeris.SE_SIDM_GALEQ_IAU1958, 31)
        self.assertEqual(lxephemeris.SE_SIDM_GALEQ_TRUE, 32)
        self.assertEqual(lxephemeris.SE_SIDM_GALEQ_MULA, 33)
        self.assertEqual(lxephemeris.SE_SIDM_GALALIGN_MARDYKS, 34)
        self.assertEqual(lxephemeris.SE_SIDM_TRUE_MULA, 35)
        self.assertEqual(lxephemeris.SE_SIDM_GALCENT_MULA_WILHELM, 36)
        self.assertEqual(lxephemeris.SE_SIDM_ARYABHATA_522, 37)
        self.assertEqual(lxephemeris.SE_SIDM_BABYL_BRITTON, 38)
        self.assertEqual(lxephemeris.SE_SIDM_TRUE_SHERATAN, 39)
        self.assertEqual(lxephemeris.SE_SIDM_TRUE_ALDEBARAN, 40)
        self.assertEqual(lxephemeris.SE_SIDM_GALEQU_FIORENZA, 41)
        self.assertEqual(lxephemeris.SE_SIDM_VALENS_MOON, 42)
        self.assertEqual(lxephemeris.SE_SIDM_USER, 255)

    def test_10_solar_terms_2026_and_seasons(self) -> None:
        """Verify 2026 Chinese 24 Solar Terms (R12) and seasons_in_year."""
        jd_start_2026 = lxephemeris.julday(2026, 1, 1, 0.0, "gregorian")
        jd_end_2026 = lxephemeris.julday(2026, 12, 31, 23.99, "gregorian")
        mosh_eng = lxephemeris.Engine.new_moshier(max_threads=2)

        for label, eng in [("SPK", self.kernel_engine), ("Moshier", mosh_eng)]:
            terms = eng.solar_terms_in_range(jd_start_2026, jd_end_2026)
            self.assertEqual(len(terms), 24, f"{label}: expected 24 solar terms in 2026")

            # Check first term of 2026 is Xiaohan (285 deg, index 19)
            self.assertEqual(terms[0]["pinyin_name"], "Xiaohan")
            self.assertEqual(terms[0]["term_index"], 19)
            self.assertAlmostEqual(terms[0]["target_lon_deg"], 285.0, delta=1e-12)
            self.assertFalse(terms[0]["is_major_term"])

            # Verify chronological spacing and keys
            for i, term in enumerate(terms):
                for key in (
                    "term_index",
                    "target_lon_deg",
                    "pinyin_name",
                    "is_major_term",
                    "jd_tt",
                    "jd_ut1",
                    "delta_t_sec",
                ):
                    self.assertIn(key, term)
                self.assertEqual(term["is_major_term"], term["term_index"] % 2 == 0)
                if i > 0:
                    dt = term["jd_tt"] - terms[i - 1]["jd_tt"]
                    self.assertGreater(dt, 14.0)
                    self.assertLess(dt, 16.5)

            # Find Dongzhi (Winter Solstice, 270 deg, index 18) -> 2026-12-21 ~20:50 UT (never June!)
            dongzhi = next(t for t in terms if t["pinyin_name"] == "Dongzhi")
            self.assertEqual(dongzhi["term_index"], 18)
            self.assertAlmostEqual(dongzhi["target_lon_deg"], 270.0, delta=1e-12)
            self.assertTrue(dongzhi["is_major_term"])
            y, m, d, hr = lxephemeris.revjul(dongzhi["jd_ut1"], "gregorian")
            self.assertEqual((y, m, d), (2026, 12, 21))
            self.assertAlmostEqual(hr, 20.0 + 50.0 / 60.0, delta=0.05)

            # Cross-check seasons_in_year(2026)
            seasons = eng.seasons_in_year(2026)
            self.assertEqual(seasons["year"], 2026)
            chunfen = next(t for t in terms if t["term_index"] == 0)
            xiazhi = next(t for t in terms if t["term_index"] == 6)
            qiufen = next(t for t in terms if t["term_index"] == 12)
            self.assertAlmostEqual(seasons["vernal_equinox_jd_tt"], chunfen["jd_tt"], delta=1e-8)
            self.assertAlmostEqual(seasons["vernal_equinox_jd_ut1"], chunfen["jd_ut1"], delta=1e-8)
            self.assertAlmostEqual(seasons["summer_solstice_jd_tt"], xiazhi["jd_tt"], delta=1e-8)
            self.assertAlmostEqual(seasons["summer_solstice_jd_ut1"], xiazhi["jd_ut1"], delta=1e-8)
            self.assertAlmostEqual(seasons["autumnal_equinox_jd_tt"], qiufen["jd_tt"], delta=1e-8)
            self.assertAlmostEqual(seasons["autumnal_equinox_jd_ut1"], qiufen["jd_ut1"], delta=1e-8)
            self.assertAlmostEqual(seasons["winter_solstice_jd_tt"], dongzhi["jd_tt"], delta=1e-8)
            self.assertAlmostEqual(seasons["winter_solstice_jd_ut1"], dongzhi["jd_ut1"], delta=1e-8)

    def test_11_search_body_longitude_and_human_design_epoch(self) -> None:
        """Verify search_body_longitude (forward/backward/sidereal) and human_design_epoch (-88 deg)."""
        jd_start = lxephemeris.julday(2026, 3, 1, 0.0, "gregorian")
        vernal = self.kernel_engine.search_body_longitude("Sun", jd_start, 0.0, forward=True, max_search_days=40.0)
        for key in (
            "body",
            "target_lon_deg",
            "jd_tt",
            "jd_ut1",
            "delta_t_sec",
            "actual_lon_deg",
            "lat_deg",
            "dist_au",
            "dist_km",
            "lon_speed_deg_day",
            "is_retrograde",
        ):
            self.assertIn(key, vernal)
        self.assertEqual(vernal["body"], "Sun")
        self.assertAlmostEqual(vernal["target_lon_deg"], 0.0, delta=1e-12)
        lon_err = ((vernal["actual_lon_deg"] + 180.0) % 360.0) - 180.0
        self.assertLess(abs(lon_err), 1e-7)
        self.assertFalse(vernal["is_retrograde"])

        # Sidereal longitude crossing (Lahiri)
        sid_opts = lxephemeris.CalcOptions(zodiac="sidereal", ayanamsa="lahiri")
        sid_ing = self.kernel_engine.search_body_longitude(
            "Sun", jd_start, 0.0, forward=True, max_search_days=60.0, options=sid_opts
        )
        self.assertGreater(sid_ing["jd_ut1"], vernal["jd_ut1"] + 20.0)

        # Human Design -88 deg solar arc prenatal epoch
        birth_jd = lxephemeris.julday(1990, 6, 15, 14.5, "gregorian")
        sun_birth = self.kernel_engine.calc("Sun", birth_jd, 0.0, "UT1")
        hd = self.kernel_engine.human_design_epoch(birth_jd)
        dt_days = sun_birth.jd_tt - hd["jd_tt"]
        self.assertGreaterEqual(dt_days, 87.0)
        self.assertLessEqual(dt_days, 91.0)
        arc_deg = (sun_birth.lon_deg - hd["actual_lon_deg"]) % 360.0
        self.assertAlmostEqual(arc_deg, 88.0, delta=1e-7)
        self.assertFalse(hd["is_retrograde"])

    def test_12_lunar_phases_and_28_phase_summary(self) -> None:
        """Verify search_lunar_phase, lunar_phases_in_range, and lunar_phase_summary (8/28-phase Kerykeion replacement)."""
        jd_start = lxephemeris.julday(2026, 1, 1, 0.0, "gregorian")
        jd_end = lxephemeris.julday(2026, 2, 15, 0.0, "gregorian")

        # Principal phase search
        for elong, expected_kind in [
            (0.0, "new_moon"),
            (90.0, "first_quarter"),
            (180.0, "full_moon"),
            (270.0, "last_quarter"),
            (45.0, None),
        ]:
            ev = self.kernel_engine.search_lunar_phase(jd_start, elong, forward=True)
            self.assertEqual(ev["kind"], expected_kind)
            self.assertAlmostEqual(ev["target_elongation_deg"], elong, delta=1e-9)
            actual_elong = (ev["moon_lon_deg"] - ev["sun_lon_deg"]) % 360.0
            diff = ((actual_elong - elong + 180.0) % 360.0) - 180.0
            self.assertLess(abs(diff), 1e-7)
            self.assertGreaterEqual(ev["moon_dist_km"], 356_000.0)
            self.assertLessEqual(ev["moon_dist_km"], 406_800.0)

        # Lunar phases in range
        phases = self.kernel_engine.lunar_phases_in_range(jd_start, jd_end)
        self.assertGreaterEqual(len(phases), 6)
        for p in phases:
            self.assertIn(p["kind"], ("new_moon", "first_quarter", "full_moon", "last_quarter"))
            self.assertGreaterEqual(p["moon_dist_km"], 356_000.0)
            self.assertLessEqual(p["moon_dist_km"], 406_800.0)

        # Lunar phase summary at 2026-01-25 12:00 UT
        probe_jd = lxephemeris.julday(2026, 1, 25, 12.0, "gregorian")
        summary = self.kernel_engine.lunar_phase_summary(probe_jd)
        expected_keys = (
            "jd_tt",
            "jd_ut1",
            "elongation_deg",
            "illumination_fraction",
            "phase_angle_deg",
            "moon_lon_deg",
            "moon_lat_deg",
            "moon_dist_au",
            "moon_dist_km",
            "apparent_diameter_deg",
            "magnitude_v",
            "age_days",
            "synodic_fraction",
            "phase_8_index",
            "phase_8_name",
            "phase_28_index",
            "prev_new_moon_jd_ut1",
            "next_new_moon_jd_ut1",
            "next_first_quarter_jd_ut1",
            "next_full_moon_jd_ut1",
            "next_last_quarter_jd_ut1",
        )
        for key in expected_keys:
            self.assertIn(key, summary)
        self.assertGreaterEqual(summary["elongation_deg"], 0.0)
        self.assertLess(summary["elongation_deg"], 360.0)
        self.assertGreaterEqual(summary["illumination_fraction"], 0.0)
        self.assertLessEqual(summary["illumination_fraction"], 1.0)
        self.assertGreaterEqual(summary["moon_dist_km"], 356_000.0)
        self.assertLessEqual(summary["moon_dist_km"], 406_800.0)
        self.assertGreater(summary["age_days"], 0.0)
        self.assertLess(summary["age_days"], 30.0)
        self.assertGreater(summary["synodic_fraction"], 0.0)
        self.assertLess(summary["synodic_fraction"], 1.0)
        self.assertGreaterEqual(summary["phase_8_index"], 0)
        self.assertLess(summary["phase_8_index"], 8)
        self.assertTrue(bool(summary["phase_8_name"]))
        self.assertGreaterEqual(summary["phase_28_index"], 1)
        self.assertLessEqual(summary["phase_28_index"], 28)
        self.assertLess(summary["prev_new_moon_jd_ut1"], summary["jd_ut1"])
        self.assertGreater(summary["next_new_moon_jd_ut1"], summary["jd_ut1"])

    def test_13_stations_and_retrograde_triple_return_crossings(self) -> None:
        """Verify tropical/sidereal stations and 3-pass retrograde longitude crossings."""
        jd_start_2026 = lxephemeris.julday(2026, 1, 1, 0.0, "gregorian")
        jd_end_2026 = lxephemeris.julday(2026, 12, 31, 0.0, "gregorian")
        jd_end_triple = lxephemeris.julday(2027, 4, 1, 0.0, "gregorian")

        sat_stations = self.kernel_engine.search_stations_in_range("Saturn", jd_start_2026, jd_end_2026)
        self.assertEqual(len(sat_stations), 2)
        self.assertEqual(sat_stations[0]["body"], "Saturn")
        self.assertEqual(sat_stations[0]["kind"], "station_retrograde")
        self.assertEqual(sat_stations[1]["kind"], "station_direct")
        self.assertLess(sat_stations[0]["jd_tt"], sat_stations[1]["jd_tt"])

        # Midpoint longitude between the two stations is crossed 3 times (direct -> retro -> direct)
        mid_lon = 0.5 * (sat_stations[0]["lon_deg"] + sat_stations[1]["lon_deg"])
        crossings = self.kernel_engine.body_longitude_crossings_in_range(
            "Saturn", jd_start_2026, jd_end_triple, mid_lon
        )
        self.assertEqual(len(crossings), 3)
        self.assertFalse(crossings[0]["is_retrograde"])
        self.assertTrue(crossings[1]["is_retrograde"])
        self.assertFalse(crossings[2]["is_retrograde"])
        for c in crossings:
            err = ((c["actual_lon_deg"] - mid_lon + 180.0) % 360.0) - 180.0
            self.assertLess(abs(err), 1e-7)

        # Sidereal stations (Lahiri) for Mercury, Mars, Saturn
        sid_opts = lxephemeris.CalcOptions(zodiac="sidereal", ayanamsa="lahiri")
        jd_2027_mid = lxephemeris.julday(2027, 6, 1, 0.0, "gregorian")
        for body in ("Mercury", "Mars", "Saturn"):
            st_sid = self.kernel_engine.search_stations_in_range(
                body, jd_start_2026, jd_2027_mid, options=sid_opts
            )
            self.assertGreaterEqual(len(st_sid), 2)
            for st in st_sid:
                self.assertIn(st["kind"], ("station_retrograde", "station_direct"))
                check = self.kernel_engine.calc(body, st["jd_tt"], 0.0, "TT", options=sid_opts)
                self.assertLess(abs(check.lon_speed_deg_day), 1e-6)

    def test_14_lunar_declination_standstills_2025(self) -> None:
        """Verify 2025 Major Lunar Standstill declination extrema (|dec| > 28 deg)."""
        jd_start = lxephemeris.julday(2025, 3, 1, 0.0, "gregorian")
        jd_end = lxephemeris.julday(2025, 3, 31, 0.0, "gregorian")

        extrema = self.kernel_engine.search_declination_extrema_in_range("Moon", jd_start, jd_end)
        self.assertGreaterEqual(len(extrema), 2)
        for i, ex in enumerate(extrema):
            for key in (
                "body",
                "kind",
                "jd_tt",
                "jd_ut1",
                "delta_t_sec",
                "ra_deg",
                "dec_deg",
                "dist_au",
                "dist_km",
            ):
                self.assertIn(key, ex)
            self.assertEqual(ex["body"], "Moon")
            self.assertGreater(abs(ex["dec_deg"]), 28.0)
            if ex["kind"] == "northern_maximum":
                self.assertGreater(ex["dec_deg"], 28.0)
            else:
                self.assertEqual(ex["kind"], "southern_minimum")
                self.assertLess(ex["dec_deg"], -28.0)
            if i > 0:
                self.assertNotEqual(ex["kind"], extrema[i - 1]["kind"])

    def test_15_brady_diurnal_fixed_star_parans_subject_d_london(self) -> None:
        """Verify Bernadette Brady diurnal fixed-star parans for Subject D (London 1969-05-17 14:00 UT)."""
        london = lxephemeris.Observer(51.5074, -0.1278, 0.0)
        jd_subj_d = lxephemeris.julday(1969, 5, 17, 14.0, "gregorian")

        parans_4m = self.kernel_engine.find_parans(jd_subj_d, london, orb_minutes=4.0)
        self.assertGreater(len(parans_4m), 0)

        valid_angles = ("rising", "culminating", "setting", "anti-culminating")
        for i, p in enumerate(parans_4m):
            for key in (
                "star_name",
                "star_bayer",
                "star_hip",
                "planet",
                "paran_type",
                "star_angle",
                "planet_angle",
                "star_event_jd_ut1",
                "planet_event_jd_ut1",
                "time_orb_minutes",
            ):
                self.assertIn(key, p)
            self.assertEqual(p["paran_type"], p["star_angle"])
            self.assertIn(p["star_angle"], valid_angles)
            self.assertIn(p["planet_angle"], valid_angles)
            self.assertGreaterEqual(p["time_orb_minutes"], 0.0)
            self.assertLessEqual(p["time_orb_minutes"], 4.0)
            if i > 0:
                self.assertGreaterEqual(p["time_orb_minutes"], parans_4m[i - 1]["time_orb_minutes"])

        # Monotonicity of orb counts: 1m <= 4m <= 10m
        parans_1m = self.kernel_engine.find_parans(jd_subj_d, london, orb_minutes=1.0)
        parans_10m = self.kernel_engine.find_parans(jd_subj_d, london, orb_minutes=10.0)
        self.assertLessEqual(len(parans_1m), len(parans_4m))
        self.assertLessEqual(len(parans_4m), len(parans_10m))

        # Circumpolar stars at 51.5N (Dubhe, Polaris, Capella) with Sun over 720m -> 24 meridian-transit parans
        circumpolar = self.kernel_engine.find_parans(
            jd_subj_d,
            london,
            orb_minutes=720.0,
            bodies=["Sun"],
            star_names=["Dubhe", "Polaris", "Capella"],
        )
        self.assertEqual(len(circumpolar), 24)
        for cp in circumpolar:
            self.assertIn(cp["star_angle"], ("culminating", "anti-culminating"))

    def test_16_wave_h_constructors_time_and_serialization(self) -> None:
        """Verify Engine(), SweAdapter(), jd_utc_from_local, AttrDict, .to_dict(), .to_json(), and repr()."""
        import json

        # 1. Default constructors
        eng_default = lxephemeris.Engine()
        sun_def = eng_default.calc("Sun", JD_J2000)
        self.assertAlmostEqual(sun_def.lon_deg, 280.368, delta=0.1)

        eng_dir = lxephemeris.Engine(str(DATA_DIR))
        self.assertFalse(eng_dir.is_moshier_mode())

        swe_default = lxephemeris.SweAdapter()
        xx_def, _ = swe_default.calc_ut(JD_J2000, lxephemeris.SE_SUN, 258)
        self.assertAlmostEqual(xx_def[0], 280.368, delta=0.1)

        swe_dir = lxephemeris.SweAdapter(str(DATA_DIR))
        xx_dir, _ = swe_dir.calc_ut(JD_J2000, lxephemeris.SE_SUN, 258)
        self.assertAlmostEqual(xx_dir[0], xx_def[0], delta=0.05)

        swe_from_eng = lxephemeris.SweAdapter.from_engine(self.kernel_engine)
        xx_fe, _ = swe_from_eng.calc_ut(JD_J2000, lxephemeris.SE_SUN, 258)
        self.assertAlmostEqual(xx_fe[0], xx_dir[0], delta=1e-6)

        # 2. jd_utc_from_local vs utc_from_local
        tup = lxephemeris.utc_from_local(2026, 3, 20, 12, 0, 0.0, "Europe/London", 0)
        self.assertEqual(len(tup), 5)
        jd_only = lxephemeris.jd_utc_from_local(2026, 3, 20, 12, 0, 0.0, "Europe/London", 0)
        self.assertIsInstance(jd_only, float)
        self.assertAlmostEqual(jd_only, tup[0], delta=1e-12)

        # 3. .to_dict(), .to_json(), and __repr__() on pyclass result types
        calc_res = self.kernel_engine.calc("Sun", JD_J2000)
        calc_dict = calc_res.to_dict()
        self.assertIsInstance(calc_dict, dict)
        self.assertIsInstance(calc_dict, lxephemeris.AttrDict)
        self.assertAlmostEqual(calc_dict.lon_deg, calc_res.lon_deg, delta=1e-12)
        self.assertAlmostEqual(calc_dict["lon_deg"], calc_res.lon_deg, delta=1e-12)
        calc_json = json.loads(calc_res.to_json())
        self.assertAlmostEqual(calc_json["lon_deg"], calc_res.lon_deg, delta=1e-12)
        self.assertIn("CalcResult", repr(calc_res))

        h_res = self.kernel_engine.houses(JD_J2000, 51.5074, -0.1278, "P")
        h_dict = h_res.to_dict()
        self.assertEqual(len(h_dict.cusps_deg), 12)
        self.assertIn("HouseResult", repr(h_res))
        self.assertEqual(len(json.loads(h_res.to_json())["cusps_deg"]), 12)

        ayan_res = self.kernel_engine.ayanamsa("lahiri", JD_J2000)
        self.assertAlmostEqual(ayan_res.to_dict().effective_deg, ayan_res.effective_deg, delta=1e-12)
        self.assertIn("AyanamsaValue", repr(ayan_res))

        dt_res = self.kernel_engine.delta_t(JD_J2000)
        self.assertAlmostEqual(dt_res.to_dict().delta_t_sec, dt_res.delta_t_sec, delta=1e-12)
        self.assertIn("DeltaTInfo", repr(dt_res))

        st_res = self.kernel_engine.sidtime(JD_J2000)
        self.assertAlmostEqual(st_res.to_dict().gast_hours, st_res.gast_hours, delta=1e-12)
        self.assertIn("SiderealTimeInfo", repr(st_res))

        en_res = self.kernel_engine.ecl_nut(JD_J2000)
        self.assertAlmostEqual(en_res.to_dict().true_obliquity_deg, en_res.true_obliquity_deg, delta=1e-12)
        self.assertIn("EclNutInfo", repr(en_res))

        star_meta = self.kernel_engine.star_info("Regulus")
        self.assertEqual(star_meta.to_dict().canonical_name, "Regulus")
        self.assertIn("StarMetadata", repr(star_meta))

        obs = lxephemeris.Observer(51.5074, -0.1278, 0.0)
        self.assertAlmostEqual(obs.to_dict().latitude_deg, 51.5074, delta=1e-12)
        self.assertIn("Observer", repr(obs))

        sjd = lxephemeris.SplitJd(JD_J2000, 0.25)
        self.assertAlmostEqual(sjd.to_dict().jd1, sjd.jd1, delta=1e-12)
        self.assertAlmostEqual(sjd.to_dict().jd2, sjd.jd2, delta=1e-12)
        self.assertIn("SplitJd", repr(sjd))

        # 4. AttrDict attribute access and AttributeError on missing key
        summary = self.kernel_engine.lunar_phase_summary(JD_J2000, scale="TT")
        self.assertIsInstance(summary, dict)
        self.assertIsInstance(summary, lxephemeris.AttrDict)
        self.assertEqual(summary.phase_8_name, summary["phase_8_name"])
        with self.assertRaises(AttributeError):
            _ = summary.non_existent_attribute_xyz

    def test_17_wave_h_returns_aspects_transits_and_voc_moon(self) -> None:
        """Verify planetary_return, aspects, search_exact_aspect, mutual_aspects_in_range, transits, and VoC Moon."""
        natal_jd = lxephemeris.julday(1990, 6, 15, 12.0, "gregorian")
        search_jd = lxephemeris.julday(2026, 6, 1, 0.0, "gregorian")

        # 1. Solar Return
        sr = self.kernel_engine.planetary_return("Sun", natal_jd, search_jd)
        self.assertIsInstance(sr, lxephemeris.AttrDict)
        self.assertEqual(sr.body, "Sun")
        self.assertAlmostEqual(sr.actual_lon_deg, sr.natal_lon_deg, delta=1e-6)
        self.assertGreater(sr.jd_ut1, search_jd)

        # Lunar returns in a 60-day window
        lrs = self.kernel_engine.planetary_returns_in_range("Moon", natal_jd, search_jd, search_jd + 60.0)
        self.assertGreaterEqual(len(lrs), 2)
        for lr in lrs:
            self.assertAlmostEqual(lr.actual_lon_deg, lr.natal_lon_deg, delta=1e-6)

        # 2. Chart aspects, synastry aspects, declination aspects
        chart_asps = self.kernel_engine.calc_chart_aspects(
            JD_J2000,
            bodies=["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"],
            preset="all",
        )
        self.assertGreater(len(chart_asps), 0)
        self.assertTrue(hasattr(chart_asps[0], "aspect_name"))
        self.assertTrue(hasattr(chart_asps[0], "orb_deg"))

        syn_asps = self.kernel_engine.calc_synastry_aspects(
            natal_jd,
            ["Sun", "Moon", "Venus"],
            JD_J2000,
            ["Sun", "Moon", "Mars"],
            preset="all",
        )
        self.assertGreater(len(syn_asps), 0)
        self.assertTrue(hasattr(syn_asps[0], "aspect_name"))

        dec_asps = self.kernel_engine.calc_declination_aspects(
            JD_J2000,
            bodies=["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"],
            max_orb_deg=3.0,
        )
        self.assertGreater(len(dec_asps), 0)
        self.assertIn(dec_asps[0].kind, ("Parallel", "Contraparallel"))

        # 3. Exact aspect search, mutual aspects in range, natal transits, and VoC Moon
        exact_conj = self.kernel_engine.search_exact_aspect("Moon", "Sun", 0.0, search_jd, forward=True)
        self.assertAlmostEqual(exact_conj.target_angle_deg, 0.0, delta=1e-9)
        self.assertGreater(exact_conj.jd_ut1, search_jd)

        mutual = self.kernel_engine.mutual_aspects_in_range(
            search_jd,
            search_jd + 5.0,
            bodies=["Moon", "Sun", "Mars"],
            preset="major",
        )
        self.assertGreater(len(mutual), 0)

        transits = self.kernel_engine.natal_transits_in_range(
            natal_jd,
            search_jd,
            search_jd + 5.0,
            transiting_bodies=["Moon"],
            natal_bodies=["Sun", "Moon"],
            preset="major",
        )
        self.assertGreater(len(transits), 0)

        voc = self.kernel_engine.void_of_course_moon_in_range(search_jd, search_jd + 7.0)
        self.assertGreaterEqual(len(voc), 2)
        for v in voc:
            self.assertLessEqual(v.start_jd_ut1, v.end_jd_ut1)
            self.assertGreater(v.duration_hours, 0.0)

    def test_18_wave_h_panchanga_dasha_vargas_lots_and_snapshot(self) -> None:
        """Verify Vedic Panchanga, Vimshottari Dasha, Vargas, Sect, Planetary Hours, Lots, and Chart Snapshot."""
        jd_2026 = lxephemeris.julday(2026, 3, 20, 12.0, "gregorian")

        # 1. Zodiac placement & Vargas
        zp_mod = lxephemeris.map_zodiac_placement(125.5)
        zp_eng = self.kernel_engine.map_zodiac_placement(125.5)
        self.assertEqual(zp_mod.sign_name, "Leo")
        self.assertEqual(zp_eng.sign_name, "Leo")
        self.assertGreaterEqual(zp_mod.hd_gate, 1)
        self.assertLessEqual(zp_mod.hd_gate, 64)

        v9 = lxephemeris.calc_varga(125.5, 9)
        self.assertEqual(v9.division, 9)
        self.assertEqual(v9.varga_name, "Navamsha")
        v16 = self.kernel_engine.calc_shodashavarga(125.5)
        self.assertEqual(len(v16), 16)

        defs_major = lxephemeris.aspect_definitions("major")
        self.assertEqual(len(defs_major), 5)
        self.assertEqual(defs_major[0].name, "Conjunction")
        defs_all = lxephemeris.aspect_definitions("all")
        self.assertEqual(len(defs_all), 11)

        # 2. Panchanga & Vimshottari Dasha
        panch = self.kernel_engine.panchanga(jd_2026, 28.6139, 77.2090)
        self.assertGreaterEqual(panch.tithi.tithi_index, 1)
        self.assertLessEqual(panch.tithi.tithi_index, 30)
        self.assertGreaterEqual(panch.nakshatra.nakshatra_index, 1)
        self.assertLessEqual(panch.nakshatra.nakshatra_index, 27)
        self.assertIsNotNone(panch.rahukalam)

        dasha = self.kernel_engine.vimshottari_dasha(jd_2026)
        self.assertEqual(len(dasha.periods), 9)
        self.assertGreaterEqual(len(dasha.periods[0].sub_periods), 1)
        self.assertLessEqual(len(dasha.periods[0].sub_periods), 9)
        self.assertEqual(len(dasha.periods[1].sub_periods), 9)

        # 3. Chart Sect, Planetary Hours, Hermetic Lots, and Chart Snapshot
        sect = self.kernel_engine.chart_sect(jd_2026, 51.5074, -0.1278)
        self.assertIsInstance(sect.is_day_chart, bool)

        ph = self.kernel_engine.planetary_hours(jd_2026, 51.5074, -0.1278)
        self.assertEqual(len(ph.hours), 24)
        self.assertGreaterEqual(ph.current_hour_index, 1)
        self.assertLessEqual(ph.current_hour_index, 24)

        lots = self.kernel_engine.calc_lots(jd_2026, 51.5074, -0.1278)
        self.assertEqual(len(lots), 7)
        self.assertEqual(lots[0].name, "Fortune")

        snap = self.kernel_engine.chart_snapshot(
            jd_2026,
            51.5074,
            -0.1278,
            bodies=["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"],
        )
        self.assertEqual(len(snap.bodies), 7)
        self.assertEqual(len(snap.lots), 7)
        self.assertEqual(snap.bodies[0].body_name, "Sun")
        self.assertGreaterEqual(snap.bodies[0].house_number, 1)
        self.assertLessEqual(snap.bodies[0].house_number, 12)

    def test_19_swe_adapter_drop_in_compatibility_enhancements(self) -> None:
        """Verify all 8 SweAdapter / pyswisseph drop-in compatibility enhancements."""
        swe = lxephemeris.SweAdapter(self.kernel_engine)

        # 1. get_ayanamsa_ex_ut returns (retflag: int, ayanamsa_deg: float)
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        retflag, ayan_deg = swe.get_ayanamsa_ex_ut(JD_J2000, swe.FLG_SWIEPH)
        self.assertIsInstance(retflag, int)
        self.assertIsInstance(ayan_deg, float)
        self.assertGreaterEqual(retflag, 0)
        self.assertAlmostEqual(ayan_deg, 23.85, delta=0.05)

        # 2. Automatic Moshier fallback in calc_ut / calc when epoch is outside SPK coverage even with SEFLG_SWIEPH
        jd_1ce = swe.julday(1, 1, 1, 12.0, swe.GREG_CAL)
        xx_ut, flg_ut = swe.calc_ut(jd_1ce, swe.SUN, swe.FLG_SWIEPH | swe.FLG_SPEED)
        self.assertNotEqual(flg_ut & swe.FLG_MOSEPH, 0)
        self.assertGreaterEqual(xx_ut[0], 0.0)
        self.assertLess(xx_ut[0], 360.0)

        xx_tt, flg_tt = swe.calc(jd_1ce, swe.MOON, swe.SEFLG_SWIEPH | swe.SEFLG_SPEED)
        self.assertNotEqual(flg_tt & swe.SEFLG_MOSEPH, 0)
        self.assertGreaterEqual(xx_tt[0], 0.0)
        self.assertLess(xx_tt[0], 360.0)

        # 3. azalt and azalt_rev accept both 2-element and 3-element xin sequences
        hz2 = swe.azalt(JD_J2000, swe.ECL2HOR, [-0.1278, 51.5074, 0.0], 1013.25, 15.0, [280.368, 0.0])
        hz3 = swe.azalt(JD_J2000, swe.ECL2HOR, [-0.1278, 51.5074, 0.0], 1013.25, 15.0, [280.368, 0.0, 0.983])
        self.assertAlmostEqual(hz2[0], hz3[0], delta=1e-12)
        self.assertAlmostEqual(hz2[1], hz3[1], delta=1e-12)
        self.assertAlmostEqual(hz2[2], hz3[2], delta=1e-12)

        rev2 = swe.azalt_rev(JD_J2000, swe.HOR2EQU, [-0.1278, 51.5074, 0.0], [180.0, 45.0])
        rev3 = swe.azalt_rev(JD_J2000, swe.HOR2EQU, [-0.1278, 51.5074, 0.0], [180.0, 45.0, 1.0])
        self.assertAlmostEqual(rev2[0], rev3[0], delta=1e-12)
        self.assertAlmostEqual(rev2[1], rev3[1], delta=1e-12)

        # 4. Utility methods: set_ephe_path, close, get_planet_name, house_name
        swe.set_ephe_path("/tmp/ephe")
        swe.set_ephe_path(None)
        swe.close()
        self.assertEqual(swe.get_planet_name(swe.SUN), "Sun")
        self.assertEqual(swe.get_planet_name(swe.TRUE_NODE), "True_Node")
        self.assertEqual(swe.get_planet_name(swe.CHIRON), "Chiron")
        self.assertEqual(swe.house_name(b"P"), "Placidus")
        self.assertEqual(swe.house_name("W"), "equal/ whole sign")

        # 5. Eclipse and Heliacal methods on SweAdapter
        sol_flags, sol_tret = swe.sol_eclipse_when_glob(2451395.0, swe.FLG_SWIEPH, 0, False)
        self.assertNotEqual(sol_flags & swe.ECL_TOTAL, 0)
        self.assertEqual(len(sol_tret), 10)
        self.assertGreater(sol_tret[0], 2451395.0)

        where_flags, geopos_max, sol_attr = swe.sol_eclipse_where(sol_tret[0], swe.FLG_SWIEPH)
        self.assertNotEqual(where_flags & swe.ECL_TOTAL, 0)
        self.assertEqual(len(geopos_max), 10)
        self.assertEqual(len(sol_attr), 20)
        self.assertGreater(sol_attr[0], 0.0)

        loc_flags, loc_tret, loc_attr = swe.sol_eclipse_when_loc(
            2451395.0, [11.58, 48.14, 500.0], swe.FLG_SWIEPH, False
        )
        self.assertEqual(len(loc_tret), 10)
        self.assertEqual(len(loc_attr), 20)
        self.assertAlmostEqual(loc_tret[0], sol_tret[0], delta=0.5)
        self.assertGreater(loc_attr[0], 0.9)

        lun_flags, lun_tret = swe.lun_eclipse_when(2451555.0, swe.FLG_SWIEPH, 0, False)
        self.assertNotEqual(lun_flags & swe.ECL_TOTAL, 0)
        self.assertEqual(len(lun_tret), 10)

        how_flags, lun_attr = swe.lun_eclipse_how(lun_tret[0], [-0.1278, 51.5074, 0.0], swe.FLG_SWIEPH)
        self.assertNotEqual(how_flags & swe.ECL_TOTAL, 0)
        self.assertEqual(len(lun_attr), 20)
        self.assertGreater(lun_attr[0], 1.0)

        hel_tret = swe.heliacal_ut(
            JD_J2000,
            [31.2357, 30.0444, 0.0],
            [1013.25, 20.0, 40.0, 0.25],
            [25.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            "Venus",
            swe.HELIACAL_RISING,
            swe.HELFLAG_HIGH_PRECISION,
        )
        self.assertIsInstance(hel_tret, tuple)
        self.assertEqual(len(hel_tret), 3)
        self.assertGreater(hel_tret[0], JD_J2000)
        self.assertLess(hel_tret[0], hel_tret[1])
        self.assertLess(hel_tret[1], hel_tret[2])

        # 6. Unprefixed constants, Error, and version on both module and SweAdapter instance
        for target in (lxephemeris, swe):
            self.assertEqual(target.SUN, 0)
            self.assertEqual(target.MOON, 1)
            self.assertEqual(target.FLG_JPLEPH, 1)
            self.assertEqual(target.SEFLG_JPLEPH, 1)
            self.assertEqual(target.FLG_SWIEPH, 2)
            self.assertEqual(target.FLG_MOSEPH, 4)
            self.assertEqual(target.FLG_SPEED, 256)
            self.assertEqual(target.FLG_SIDEREAL, 65536)
            self.assertEqual(target.SIDM_LAHIRI, 1)
            self.assertEqual(target.SIDM_RAMAN, 3)
            self.assertEqual(target.SIDM_KRISHNAMURTI, 5)
            self.assertEqual(target.ECL_CENTRAL, 1)
            self.assertEqual(target.ECL_NONCENTRAL, 2)
            self.assertEqual(target.ECL_TOTAL, 4)
            self.assertEqual(target.ECL_ANNULAR, 8)
            self.assertEqual(target.ECL_PARTIAL, 16)
            self.assertEqual(target.ECL_ANNULAR_TOTAL, 32)
            self.assertEqual(target.ECL_HYBRID, 32)
            self.assertEqual(target.ECL_PENUMBRAL, 64)
            self.assertEqual(target.HELIACAL_RISING, 1)
            self.assertEqual(target.HELIACAL_SETTING, 2)
            self.assertEqual(target.EVENING_FIRST, 3)
            self.assertEqual(target.MORNING_LAST, 4)
            self.assertEqual(target.HELFLAG_HIGH_PRECISION, 256)
            self.assertIs(target.Error, lxephemeris.LxEphemerisError)
            self.assertTrue(isinstance(target.version, str) and len(target.version) > 0)
            self.assertTrue(isinstance(target.__version__, str) and len(target.__version__) > 0)

        # 7. julday and revjul accept both integer calendar flags (1, 0) and string names
        jd_greg_int = lxephemeris.julday(2000, 1, 1, 12.0, lxephemeris.GREG_CAL)
        jd_greg_str = lxephemeris.julday(2000, 1, 1, 12.0, "gregorian")
        self.assertAlmostEqual(jd_greg_int, JD_J2000, delta=1e-12)
        self.assertAlmostEqual(jd_greg_str, JD_J2000, delta=1e-12)

        jd_jul_int = lxephemeris.julday(2000, 1, 1, 12.0, lxephemeris.JUL_CAL)
        jd_jul_str = lxephemeris.julday(2000, 1, 1, 12.0, "julian")
        self.assertAlmostEqual(jd_jul_int, jd_jul_str, delta=1e-12)
        self.assertAlmostEqual(jd_jul_int - jd_greg_int, 13.0, delta=1e-12)

        self.assertEqual(lxephemeris.revjul(JD_J2000, lxephemeris.GREG_CAL)[:3], (2000, 1, 1))
        self.assertEqual(lxephemeris.revjul(jd_jul_int, lxephemeris.JUL_CAL)[:3], (2000, 1, 1))
        self.assertEqual(swe.revjul(JD_J2000, "gregorian")[:3], (2000, 1, 1))

    def test_20_codex_audit_19_findings_regressions(self) -> None:
        """Verify all 19 Codex audit fixes across Python and SweAdapter."""
        swe = lxephemeris.SweAdapter(self.kernel_engine)

        # #1 & #16: heliacal_ut 3-tuple (t_start, t_opt, t_end), meteorological visibility (km), inner planets & stars
        for ev_type in (swe.HELIACAL_RISING, swe.HELIACAL_SETTING, swe.EVENING_FIRST, swe.MORNING_LAST):
            t_start, t_opt, t_end = swe.heliacal_ut(
                JD_J2000,
                [31.2357, 30.0444, 0.0],
                [1013.25, 20.0, 40.0, 25.0],  # 25 km meteorological range -> Koschmieder k_V
                [30.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                "Venus",
                ev_type,
                swe.HELFLAG_HIGH_PRECISION,
            )
            self.assertGreater(t_start, JD_J2000)
            self.assertLess(t_start, t_opt)
            self.assertLess(t_opt, t_end)

        sirius_t = swe.heliacal_ut(
            JD_J2000,
            [31.2357, 30.0444, 0.0],
            [1013.25, 20.0, 40.0, 25.0],
            [30.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            "Sirius",
            swe.HELIACAL_RISING,
            0,
        )
        self.assertEqual(len(sirius_t), 3)
        self.assertLess(sirius_t[0], sirius_t[1])

        # #2: All custom exceptions inherit from LxEphemerisError / swe.Error
        for exc_cls in (
            lxephemeris.InvalidInputError,
            lxephemeris.UnsupportedBodyError,
            lxephemeris.UnsupportedModelError,
            lxephemeris.OutOfRangeError,
            lxephemeris.PolarSingularityError,
            lxephemeris.ConvergenceError,
            lxephemeris.MissingDataError,
            lxephemeris.CorruptKernelError,
        ):
            self.assertTrue(issubclass(exc_cls, lxephemeris.LxEphemerisError))
            self.assertTrue(issubclass(exc_cls, swe.Error))

        # #5 & #6: Vimshottari Dasha birth Mahadasha skips elapsed Antardashas and clips active birth Antardasha
        birth_jd = lxephemeris.julday(1990, 6, 15, 12.0, "gregorian")
        dasha = self.kernel_engine.vimshottari_dasha(birth_jd)
        p0 = dasha.periods[0]
        self.assertGreater(len(p0.sub_periods), 0)
        self.assertLessEqual(len(p0.sub_periods), 9)
        self.assertAlmostEqual(p0.sub_periods[0].start_jd_ut1, birth_jd, delta=1e-9)
        self.assertAlmostEqual(p0.sub_periods[-1].end_jd_ut1, p0.end_jd_ut1, delta=1e-9)
        for sp in p0.sub_periods:
            self.assertGreater(sp.duration_years, 0.0)
            self.assertGreaterEqual(sp.start_jd_ut1, birth_jd - 1e-9)

        # #7: SEFLG_MOSEPH clears both SEFLG_JPLEPH and SEFLG_SWIEPH in retflag
        _, ret_mosh = swe.calc_ut(JD_J2000, swe.SUN, swe.FLG_MOSEPH | swe.FLG_SWIEPH | swe.FLG_JPLEPH)
        self.assertNotEqual(ret_mosh & swe.FLG_MOSEPH, 0)
        self.assertEqual(ret_mosh & swe.FLG_SWIEPH, 0)
        self.assertEqual(ret_mosh & swe.FLG_JPLEPH, 0)

        # #8: azalt_rev accepts both 4-arg and 5-arg forms
        rev_4 = swe.azalt_rev(JD_J2000, swe.HOR2EQU, [-0.1278, 51.5074, 0.0], [180.0, 45.0])
        rev_5 = swe.azalt_rev(JD_J2000, swe.HOR2EQU, [-0.1278, 51.5074, 0.0], 180.0, 45.0)
        self.assertAlmostEqual(rev_4[0], rev_5[0], delta=1e-12)
        self.assertAlmostEqual(rev_4[1], rev_5[1], delta=1e-12)

        # #9: rise_trans returns (res_code, tret_10_tuple) and supports int/str body + pyswisseph & C-style signatures
        res_sun, tret_sun = swe.rise_trans(JD_J2000, swe.SUN, swe.CALC_RISE, [-0.1278, 51.5074, 0.0])
        self.assertEqual(res_sun, 0)
        self.assertEqual(len(tret_sun), 10)
        self.assertGreater(tret_sun[0], JD_J2000)

        res_star, tret_star = swe.rise_trans(JD_J2000, "Regulus", swe.CALC_RISE, [-0.1278, 51.5074, 0.0])
        self.assertEqual(res_star, 0)
        self.assertGreater(tret_star[0], JD_J2000)

        res_c, tret_c = swe.rise_trans(
            JD_J2000, swe.SUN, None, swe.FLG_SWIEPH, swe.CALC_RISE, [-0.1278, 51.5074, 0.0], 1013.25, 15.0
        )
        self.assertEqual(res_c, 0)
        self.assertAlmostEqual(tret_c[0], tret_sun[0], delta=1e-6)

        # Circumpolar Sun at Svalbard in June returns (-2, (0.0, ...))
        jd_june = swe.julday(2026, 6, 21, 12.0)
        res_circ, tret_circ = swe.rise_trans(jd_june, swe.SUN, swe.CALC_RISE, [15.63, 78.22, 0.0])
        self.assertEqual(res_circ, -2)
        self.assertEqual(tret_circ[0], 0.0)

        # #10 & #11: sol_eclipse_when_glob and lun_eclipse_when ifltype filtering + SE_ECL_CENTRAL + tret[0..8]
        sol_tot_flags, sol_tot_tret = swe.sol_eclipse_when_glob(2451395.0, swe.FLG_SWIEPH, swe.ECL_TOTAL, False)
        self.assertNotEqual(sol_tot_flags & swe.ECL_TOTAL, 0)
        self.assertNotEqual(sol_tot_flags & swe.ECL_CENTRAL, 0)
        self.assertLess(sol_tot_tret[2], sol_tot_tret[4])  # P1 < U1
        self.assertLess(sol_tot_tret[4], sol_tot_tret[6])  # U1 < U2
        self.assertLess(sol_tot_tret[6], sol_tot_tret[0])  # U2 < Max
        self.assertLess(sol_tot_tret[0], sol_tot_tret[7])  # Max < U3
        self.assertLess(sol_tot_tret[7], sol_tot_tret[5])  # U3 < U4
        self.assertLess(sol_tot_tret[5], sol_tot_tret[3])  # U4 < P4

        sol_ann_flags, sol_ann_tret = swe.sol_eclipse_when_glob(2451395.0, swe.FLG_SWIEPH, swe.ECL_ANNULAR, False)
        self.assertNotEqual(sol_ann_flags & swe.ECL_ANNULAR, 0)
        self.assertGreater(sol_ann_tret[0], sol_tot_tret[0] + 100.0)

        lun_part_flags, _ = swe.lun_eclipse_when(2451555.0, swe.FLG_SWIEPH, swe.ECL_PARTIAL, False)
        self.assertNotEqual(lun_part_flags & swe.ECL_PARTIAL, 0)

        # #12: sol_eclipse_where and lun_eclipse_how evaluate at jd_ut (return 0 outside an eclipse)
        no_sol_flags, _, _ = swe.sol_eclipse_where(JD_J2000, swe.FLG_SWIEPH)
        self.assertEqual(no_sol_flags, 0)
        no_lun_flags, _ = swe.lun_eclipse_how(JD_J2000, [0.0, 51.5, 0.0], swe.FLG_SWIEPH)
        self.assertEqual(no_lun_flags, 0)

        # #13: Out-of-coverage Moshier fallback on SweAdapter methods (1000 CE)
        jd_1000ce = swe.julday(1000, 3, 20, 12.0)
        rt_code, rt_tret = swe.rise_trans(jd_1000ce, swe.SUN, swe.CALC_RISE, [0.0, 51.5, 0.0])
        self.assertEqual(rt_code, 0)
        self.assertGreater(rt_tret[0], jd_1000ce)
        _, ph_flg = swe.pheno_ut(jd_1000ce, swe.MARS, swe.FLG_SWIEPH)
        self.assertNotEqual(ph_flg & swe.FLG_MOSEPH, 0)
        _, _, _, _, na_flg = swe.nod_aps_ut(jd_1000ce, swe.MOON, swe.FLG_SWIEPH)
        self.assertNotEqual(na_flg & swe.FLG_MOSEPH, 0)
        cusps_1000, _ = swe.houses_ex(jd_1000ce, 51.5, 0.0, b"P", swe.FLG_SIDEREAL)
        self.assertEqual(len(cusps_1000), 12)

        # #14: SweAdapter default polar policy is FallbackToPorphyry (Placidus at Svalbard 78.22N succeeds)
        cusps_polar, ascmc_polar = swe.houses_ex(JD_J2000, 78.22, 15.63, b"P", 0)
        self.assertEqual(len(cusps_polar), 12)
        self.assertGreaterEqual(ascmc_polar[0], 0.0)

        # #15: void_of_course_moon_in_range includes ongoing VoC window whose ingress occurs after end_jd
        voc_full = self.kernel_engine.void_of_course_moon_in_range(JD_J2000, JD_J2000 + 6.0)
        self.assertGreaterEqual(len(voc_full), 2)
        first_voc = voc_full[0]
        if first_voc.end_jd_ut1 - first_voc.start_jd_ut1 > 1e-4:
            mid_jd = 0.5 * (first_voc.start_jd_ut1 + first_voc.end_jd_ut1)
            voc_clipped = self.kernel_engine.void_of_course_moon_in_range(JD_J2000, mid_jd)
            self.assertTrue(any(abs(w.start_jd_ut1 - first_voc.start_jd_ut1) < 1e-6 for w in voc_clipped))

        # #17: Engine native method aliases
        jd_2026_start = lxephemeris.julday(2026, 1, 1, 0.0)
        jd_2026_end = lxephemeris.julday(2026, 12, 31, 0.0)
        st_alias = self.kernel_engine.stations_in_range("Saturn", jd_2026_start, jd_2026_end)
        self.assertEqual(len(st_alias), 2)
        cross_alias = self.kernel_engine.longitude_crossings_in_range("Sun", jd_2026_start, jd_2026_start + 100.0, 0.0)
        self.assertEqual(len(cross_alias), 1)
        dyn_alias = self.kernel_engine.dynamic_aspects_in_range(
            jd_2026_start, jd_2026_start + 3.0, bodies=["Sun", "Moon"]
        )
        self.assertGreater(len(dyn_alias), 0)
        lp_detail = self.kernel_engine.lunar_phase_detail(JD_J2000)
        self.assertIn("phase_28_index", lp_detail)

        # Primary source & astrology-api enhancements: Cor Caroli (HIP 63125), paran same_angle, Bonatti mild signs
        cc_info = self.kernel_engine.star_info("Cor Caroli")
        self.assertEqual(cc_info.hip, 63125)
        self.assertEqual(cc_info.bayer, "al2CVn")
        cc_xx, cc_name, cc_flg = swe.fixstar2_ut("Cor Caroli", JD_J2000, swe.FLG_SWIEPH | swe.FLG_SPEED)
        self.assertEqual(cc_name, "Cor Caroli")
        self.assertGreaterEqual(cc_xx[0], 0.0)
        self.assertAlmostEqual(swe.fixstar2_mag("al2CVn"), 2.89, delta=1e-6)

        obs_london = lxephemeris.Observer(51.5074, -0.1278, 0.0)
        parans = self.kernel_engine.find_parans(JD_J2000, obs_london, orb_minutes=10.0, star_names=["Regulus", "Cor Caroli"])
        for p in parans:
            self.assertIn("same_angle", p)
            self.assertEqual(p.same_angle, p.star_angle == p.planet_angle)

        voc_month = self.kernel_engine.void_of_course_moon_in_range(JD_J2000, JD_J2000 + 30.0)
        self.assertGreater(len(voc_month), 8)
        for v in voc_month:
            self.assertIn("is_bonatti_mild_sign", v)
            self.assertEqual(v.is_bonatti_mild_sign, v.moon_sign_before in (1, 3, 8, 11))

        # Reviewer #2 regressions:
        # 1. Regulus heliacal rising in Cairo (30N, 31.2E, start 2026-01-01) -> 2026-09-13 02:41 UTC
        jd_2026_jan1 = swe.julday(2026, 1, 1, 0.0)
        reg_hel = swe.heliacal_ut(
            jd_2026_jan1,
            [31.2, 30.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "Regulus",
            swe.HELIACAL_RISING,
            swe.FLG_SWIEPH | swe.HELFLAG_HIGH_PRECISION,
        )
        y_r, m_r, d_r, hr_r = swe.revjul(reg_hel[0])
        self.assertEqual((y_r, m_r, d_r), (2026, 9, 13))
        self.assertAlmostEqual(hr_r, 2.0 + 41.0 / 60.0, delta=0.05)

        # 2. NASA Five Millennium Canon Saros series & member numbers:
        #    Solar 2026-08-12 -> 126/48; Lunar 2026-03-03 -> 133/27
        _, sol_2026_tret = swe.sol_eclipse_when_glob(swe.julday(2026, 8, 1, 0.0), swe.FLG_SWIEPH, 0, False)
        _, _, sol_2026_attr = swe.sol_eclipse_where(sol_2026_tret[0], swe.FLG_SWIEPH)
        self.assertEqual(int(sol_2026_attr[9]), 126)
        self.assertEqual(int(sol_2026_attr[10]), 48)

        _, lun_2026_tret = swe.lun_eclipse_when(swe.julday(2026, 3, 1, 0.0), swe.FLG_SWIEPH, 0, False)
        _, lun_2026_attr = swe.lun_eclipse_how(lun_2026_tret[0], [0.0, 0.0, 0.0], swe.FLG_SWIEPH)
        self.assertEqual(int(lun_2026_attr[9]), 133)
        self.assertEqual(int(lun_2026_attr[10]), 27)

        # Reviewer #4 Canonical Doctrine regressions:
        # 1. Sirius heliacal rising at Lisbon (38.72N, -9.14E, start 2026-01-01):
        #    Must occur in morning twilight before sunrise (Sun altitude < -6°, Sirius altitude > 0°),
        #    never in daylight with Sirius below the horizon.
        sir_hel = swe.heliacal_ut(
            jd_2026_jan1,
            [-9.14, 38.72, 0.0],
            [1013.25, 15.0, 50.0, 0.25],
            [36.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            "Sirius",
            swe.HELIACAL_RISING,
            swe.FLG_SWIEPH | swe.HELFLAG_HIGH_PRECISION,
        )
        sun_eq_lis = swe.calc_ut(sir_hel[0], swe.SUN, swe.FLG_EQUATORIAL)[0]
        sir_eq_lis = swe.fixstar2_ut("Sirius", sir_hel[0], swe.FLG_EQUATORIAL)[0]
        sun_hz_lis = swe.azalt(sir_hel[0], swe.EQU2HOR, [-9.14, 38.72, 0.0], 1013.25, 15.0, sun_eq_lis[:3])
        sir_hz_lis = swe.azalt(sir_hel[0], swe.EQU2HOR, [-9.14, 38.72, 0.0], 1013.25, 15.0, sir_eq_lis[:3])
        self.assertLess(sun_hz_lis[1], -6.0)
        self.assertGreater(sir_hz_lis[1], 0.0)

        # 2. Canonical BPHS Ch. 6 vv. 33–35 D60 (Shashtiamsha): ignores natal Rashi index and counts from Aries.
        v60_tau0 = self.kernel_engine.calc_varga(30.1, 60)
        self.assertEqual(v60_tau0.sign_index, 0)
        self.assertEqual(v60_tau0.sign_name, "Aries")
        self.assertEqual(v60_tau0.sign_lord, "Mars")
        v60_tau1 = self.kernel_engine.calc_varga(30.6, 60)
        self.assertEqual(v60_tau1.sign_index, 1)
        self.assertEqual(v60_tau1.sign_name, "Taurus")
        self.assertEqual(v60_tau1.sign_lord, "Venus")

        # 3. Native polar Ascendant & Hermetic Lots at Longyearbyen (78.2232N, 15.6469E) on 2026-10-25 15:00 UT1:
        #    Must select the physically rising point (~355.321°, dh/dt > 0), not the setting point (~175.321°).
        jd_sval = swe.julday(2026, 10, 25, 15.0)
        h_sval = self.kernel_engine.houses(jd_sval, latitude_deg=78.2232, longitude_deg=15.6469, system="W", scale="UT1")
        self.assertAlmostEqual(h_sval.asc_deg, 355.321421, delta=0.01)
        self.assertLess((h_sval.asc_deg - h_sval.mc_deg) % 360.0, 180.0)
        # Verify physical altitude of ASC point is increasing (+dh/dt > 0) over ±30 seconds
        dt_30s = 30.0 / 86400.0
        alt_minus = swe.azalt(jd_sval - dt_30s, swe.ECL2HOR, [15.6469, 78.2232, 0.0], 0.0, 0.0, [h_sval.asc_deg, 0.0])[1]
        alt_plus = swe.azalt(jd_sval + dt_30s, swe.ECL2HOR, [15.6469, 78.2232, 0.0], 0.0, 0.0, [h_sval.asc_deg, 0.0])[1]
        self.assertLess(alt_minus, 0.0)
        self.assertGreater(alt_plus, 0.0)

        lots_sval = self.kernel_engine.calc_lots(jd_sval, 78.2232, 15.6469, 0.0, "W")
        self.assertAlmostEqual(lots_sval[0].lon_deg, 182.634221, delta=0.01)
        self.assertAlmostEqual(lots_sval[1].lon_deg, 168.008621, delta=0.01)

        # 4. Pullen Sinusoidal Ratio ('Q') at polar latitude 70N (jd=2460482.0, lat=70.0, lon=25.0)
        q_cusps, q_ascmc = swe.houses_ex(2460482.0, 70.0, 25.0, b"Q", 0)
        self.assertEqual(len(q_cusps), 12)
        self.assertAlmostEqual(q_cusps[0], q_ascmc[0], delta=1e-9)
        self.assertAlmostEqual(q_cusps[9], q_ascmc[1], delta=1e-9)
        # Cusps 10 -> 11 -> 12 -> 1 must progress monotonically across the small quadrant
        q4_span = (q_ascmc[0] - q_ascmc[1]) % 360.0
        self.assertGreater((q_cusps[10] - q_cusps[9]) % 360.0, 0.0)
        self.assertLess((q_cusps[10] - q_cusps[9]) % 360.0, q4_span)
        self.assertGreater((q_cusps[11] - q_cusps[10]) % 360.0, 0.0)
        self.assertLess((q_cusps[11] - q_cusps[9]) % 360.0, q4_span)

        # 5. Native SweAdapter utilities: get_ayanamsa_name, difdeg2n, degnorm, cotrans, house_name, keyword houses_ex, minor-body range guard
        self.assertEqual(swe.get_ayanamsa_name(0), "Fagan/Bradley")
        self.assertEqual(swe.get_ayanamsa_name(1), "Lahiri")
        self.assertEqual(swe.get_ayanamsa_name(46), "Lahiri ICRC")
        self.assertEqual(lxephemeris.get_ayanamsa_name(1), "Lahiri")
        self.assertEqual(swe.house_name(bytearray(b"E")), "equal")
        self.assertEqual(swe.house_name(ord("G")), "Gauquelin sectors")
        self.assertAlmostEqual(swe.difdeg2n(10.0, 350.0), 20.0, delta=1e-12)
        self.assertAlmostEqual(swe.degnorm(-30.0), 330.0, delta=1e-12)
        c2 = swe.cotrans((120.0, 5.0, 1.0), -23.4392911)
        c4 = swe.cotrans(120.0, 5.0, 1.0, -23.4392911)
        self.assertAlmostEqual(c2[0], c4[0], delta=1e-12)
        self.assertAlmostEqual(c2[1], c4[1], delta=1e-12)
        kw_cusps, kw_ascmc = swe.houses_ex(tjdut=JD_J2000, lat=51.5, lon=-0.12, hsys=b"E", flags=0)
        self.assertEqual(len(kw_cusps), 12)
        self.assertAlmostEqual(kw_cusps[0], kw_ascmc[0], delta=1e-9)
        with self.assertRaises(lxephemeris.OutOfRangeError):
            swe.calc_ut(swe.julday(1700, 1, 1, 0.0), swe.CHIRON, swe.FLG_SWIEPH)

        # 6. Codex Round 2 regressions: calc(MEAN_APOG) 3D latitude, FLG_SPEED3, FLG_TOPOCTR|FLG_HELCTR,
        #    houses_ex2 speeds, houses_ex(FLG_RADIANS), Polich-Page house_pos('T'), outer planet heliacal_ut(3/4)
        apog_tt, _ = swe.calc(JD_J2000, swe.MEAN_APOG, swe.FLG_SWIEPH | swe.FLG_SPEED)
        apog_ut, _ = swe.calc_ut(JD_J2000, swe.MEAN_APOG, swe.FLG_SWIEPH | swe.FLG_SPEED)
        self.assertGreater(abs(apog_tt[1]), 3.4)
        self.assertAlmostEqual(apog_tt[1], apog_ut[1], delta=0.001)
        sun_spd3, _ = swe.calc_ut(JD_J2000, swe.SUN, swe.FLG_SWIEPH | swe.FLG_SPEED3)
        self.assertGreater(sun_spd3[3], 0.95)
        hel_node, _ = swe.calc_ut(JD_J2000, swe.MEAN_NODE, swe.FLG_SWIEPH | swe.FLG_HELCTR)
        self.assertEqual(list(hel_node), [0.0] * 6)

        cusps_d, ascmc_d, cspd_d, aspd_d = swe.houses_ex2(JD_J2000, 51.5, 0.0, b"P", 0)
        self.assertGreater(cspd_d[0], 200.0)
        self.assertAlmostEqual(aspd_d[2], 360.9856, delta=0.05)
        cusps_r, ascmc_r = swe.houses_ex(JD_J2000, 51.5, 0.0, b"P", swe.FLG_RADIANS)
        self.assertAlmostEqual(cusps_r[0], math.radians(cusps_d[0]), delta=1e-12)
        self.assertAlmostEqual(ascmc_r[0], math.radians(ascmc_d[0]), delta=1e-12)

        t_cusps, _ = swe.houses_armc(120.0, 51.5, 23.4392911, b"T")
        for idx, c_lon in enumerate(t_cusps):
            hp_t = swe.house_pos(120.0, 51.5, 23.4392911, [c_lon, 0.0], b"T")
            self.assertAlmostEqual((hp_t - (idx + 1.0) + 6.0) % 12.0 - 6.0, 0.0, delta=1e-6)
        with self.assertRaises(ValueError):
            swe.heliacal_ut(JD_J2000, [31.2, 30.0, 0.0], [1013.25, 15.0, 50.0, 0.2], [36.0, 1.0], "Mars", 3, 0)

        # 7. Codex Round 3 regressions: Moon heliacal_ut(3/4), ancient time_equ, houses_ex2(RADIANS) speeds,
        #    T cusps at -65°, global occultation tret layout + ifltype, lun_eclipse_how horizon check, difdeg2n(180, 0)
        with self.assertRaises(ValueError):
            swe.heliacal_ut(JD_J2000, [31.2, 30.0, 0.0], [1013.25, 15.0, 50.0, 0.2], [36.0, 1.0], "Moon", 1, 0)
        moon_ev3 = swe.heliacal_ut(JD_J2000, [31.2, 30.0, 0.0], [1013.25, 15.0, 50.0, 0.2], [36.0, 1.0], "Moon", 3, 0)
        moon_ev4 = swe.heliacal_ut(JD_J2000, [31.2, 30.0, 0.0], [1013.25, 15.0, 50.0, 0.2], [36.0, 1.0], "Moon", 4, 0)
        self.assertAlmostEqual(moon_ev4[1], 2451548.65, delta=1.2)
        self.assertAlmostEqual(moon_ev3[1], 2451552.17, delta=1.2)

        jd_2000bce = swe.julday(-2000, 1, 1, 12.0, swe.GREG_CAL)
        self.assertAlmostEqual(swe.time_equ(jd_2000bce) * 86400.0, -722.35, delta=35.0)

        _, _, cspd_r, aspd_r = swe.houses_ex2(JD_J2000, 51.5, 0.0, b"P", swe.FLG_RADIANS)
        self.assertAlmostEqual(cspd_r[0], cspd_d[0], delta=1e-12)
        self.assertAlmostEqual(aspd_r[2], 360.9856, delta=0.05)

        t65_cusps, _, _, _ = swe.houses_ex2(2461041.5, -65.0, 13.0, b"T", 0)
        self.assertAlmostEqual(t65_cusps[5], 342.1706, delta=0.01)

        occ_rf, occ_tret = swe.lun_occult_when_glob(JD_J2000, swe.VENUS, None, swe.FLG_SWIEPH, swe.ECL_TOTAL, False)
        self.assertNotEqual(occ_rf & swe.ECL_TOTAL, 0)
        self.assertLess(occ_tret[2], occ_tret[4])
        self.assertLess(occ_tret[4], occ_tret[0])
        self.assertLess(occ_tret[0], occ_tret[5])
        self.assertLess(occ_tret[5], occ_tret[3])

        lhow_rf, _ = swe.lun_eclipse_how(2460394.800615769, [0.0, 51.5, 0.0], swe.FLG_SWIEPH)
        self.assertEqual(lhow_rf, 0)
        sloc_rf, _, sloc_attr = swe.sol_eclipse_when_loc(2460310.5, [-96.8, 32.8, 0.0], swe.FLG_SWIEPH, False)
        self.assertEqual(sloc_rf, 8068)
        self.assertAlmostEqual(sloc_attr[4], 8.07, delta=0.5)
        self.assertEqual(swe.difdeg2n(180.0, 0.0), -180.0)

        # 8. Codex Round 4 (gpt-6-astra) regressions: SE_INTP_PERG, extended SE_SIDM_*, SWE_* aliases,
        #    5-arg lun_occult_when_loc, 3-arg nod_aps_ut + SE_SUN/SE_EARTH, heliacal_ut integer body ID,
        #    get_ayanamsa_ex, houses_armc sun_decl_deg, and chart_snapshot bodies=[]
        self.assertEqual(swe.INTP_PERG, 22)
        self.assertEqual(lxephemeris.SE_INTP_PERG, 22)
        self.assertEqual(lxephemeris.SWE_SUN, 0)
        self.assertEqual(lxephemeris.SWE_SIDM_LAHIRI, 1)
        perg_xx, _ = swe.calc_ut(JD_J2000, swe.INTP_PERG, swe.FLG_SWIEPH | swe.FLG_SPEED)
        self.assertTrue(90.0 <= perg_xx[0] <= 94.0)

        swe.set_sid_mode(swe.SIDM_LAHIRI_ICRC)
        _, aya_icrc = swe.get_ayanamsa_ex(JD_J2000, swe.FLG_SWIEPH)
        self.assertTrue(23.0 <= aya_icrc <= 25.0)
        swe.set_sid_mode(swe.SIDM_TRUE_CITRA)
        spica_sid, _, _ = swe.fixstar2_ut("Spica", JD_J2000, swe.FLG_SWIEPH | swe.FLG_SIDEREAL)
        self.assertAlmostEqual(spica_sid[0], 180.0, delta=1e-4)
        star1_xx, star1_name, _ = swe.fixstar2_ut("1", JD_J2000, swe.FLG_SWIEPH)
        self.assertTrue(len(star1_name) > 0)
        self.assertTrue(0.0 <= star1_xx[0] < 360.0)
        swe.set_sid_mode(swe.SIDM_LAHIRI)

        # 3-arg nod_aps_ut default method=0 and SE_SUN / SE_EARTH support
        n_moon3 = swe.nod_aps_ut(JD_J2000, swe.MOON, swe.FLG_SWIEPH)
        self.assertEqual(len(n_moon3), 5)
        n_sun = swe.nod_aps_ut(JD_J2000, swe.SUN, swe.FLG_SWIEPH, swe.NODBIT_MEAN)
        n_earth = swe.nod_aps_ut(JD_J2000, swe.EARTH, swe.FLG_SWIEPH, swe.NODBIT_MEAN)
        self.assertAlmostEqual(n_sun[2][0], 282.9124, delta=0.05)
        self.assertAlmostEqual(n_earth[2][0], 191.6875, delta=0.05)

        # 5-arg pyswisseph signature for lun_occult_when_loc(tjd_start, body, geopos, flags=..., backwards=...)
        loc_occ_rf, loc_occ_tret, loc_occ_attr = swe.lun_occult_when_loc(
            JD_J2000, swe.VENUS, [0.0, 51.5, 0.0], swe.FLG_SWIEPH, False
        )
        self.assertGreater(loc_occ_rf, 0)
        self.assertEqual(len(loc_occ_tret), 10)
        self.assertEqual(len(loc_occ_attr), 20)
        loc_occ_rf2, _, _ = swe.lun_occult_when_loc(
            JD_J2000, swe.VENUS, "", [0.0, 51.5, 0.0], swe.FLG_SWIEPH
        )
        self.assertEqual(loc_occ_rf2, loc_occ_rf)

        # heliacal_ut accepts integer body ID (e.g. swe.VENUS = 3)
        hel_int = swe.heliacal_ut(JD_J2000, [31.2, 30.0, 0.0], [1013.25, 15.0, 50.0, 0.2], [36.0, 1.0], swe.VENUS, 1, 0)
        self.assertEqual(len(hel_int), 3)
        self.assertGreater(hel_int[0], JD_J2000)

        # houses_armc optional sun_decl_deg for Sunshine 'i'
        sun_i_cusps, _ = swe.houses_armc(120.0, 51.5, 23.4392911, b"i", 10.0)
        self.assertEqual(len(sun_i_cusps), 12)

        # chart_snapshot with bodies=[] returns 0 bodies and sidereal house_pos_3d consistency
        empty_snap = self.kernel_engine.chart_snapshot(JD_J2000, 51.5, 0.0, bodies=[])
        self.assertEqual(len(empty_snap.bodies), 0)
        sid_snap = self.kernel_engine.chart_snapshot(
            JD_J2000,
            51.5,
            0.0,
            bodies=["Sun", "Moon"],
            system="W",
            zodiac="sidereal",
            ayanamsa="lahiri",
        )
        self.assertEqual(len(sid_snap.bodies), 2)
        for b in sid_snap.bodies:
            self.assertEqual(int(b.house_pos_3d), b.house_number)


if __name__ == "__main__":
    unittest.main()

