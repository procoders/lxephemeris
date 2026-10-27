#!/usr/bin/env node
/**
 * Test Suite for lxephemeris-node:
 * Exhaustively exercises Wave B and Wave C functionality:
 * 1. Moshier mode: JsEngine.newMoshier, JsEngine.fromDataDirOrMoshier, isMoshierMode, calcMoshier
 * 2. Custom stars: loadCustomStars
 * 3. Photometric phenomena: JsEngine.phenomena, JsSweAdapter.phenoUt, JsSweAdapter.pheno
 * 4. Orbital nodes & apsides: JsEngine.nodAps, JsSweAdapter.nodApsUt, JsSweAdapter.nodAps
 * 5. Lunar occultations: JsEngine.nextLunarOccultation, nextLocalLunarOccultation,
 *    lunarOccultationsInRange, localLunarOccultationsInRange,
 *    JsSweAdapter.lunOccultWhenGlob, JsSweAdapter.lunOccultWhenLoc
 * 6. 3D House position: JsEngine.housePos, JsSweAdapter.housePos
 * 7. Equation of Time: JsEngine.timeEquation, JsSweAdapter.timeEqu
 * 8. Inverse horizontal coordinates: JsEngine.azaltRev, JsSweAdapter.azaltRev
 * 9. Extended Sidereal Ayanamsa: SE_SIDM_USER, all 43 SE_SIDM_* models, JsSweAdapter.setSidMode with t0/ayanT0
 * 10. Node bits & horizon constants: SE_NODBIT_*, SE_HOR2ECL, SE_HOR2EQU
 * 11. Error handling & edge cases
 */

import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const WORKSPACE_ROOT = path.resolve(__dirname, "..", "..");
const RELEASE_DIR = path.join(WORKSPACE_ROOT, "target", "release");
const DATA_DIR = path.join(WORKSPACE_ROOT, "data");

const require = createRequire(import.meta.url);
const addon = require(path.join(WORKSPACE_ROOT, "crates", "lxephemeris-node", "index.js"));

function assertClose(actual, expected, tol, label) {
  const diff = Math.abs(actual - expected);
  assert.ok(
    diff <= tol,
    `${label}: actual=${actual}, expected=${expected}, diff=${diff} > tol=${tol}`
  );
}

const tStart = performance.now();
console.log("=== Starting lxephemeris-node Wave B & C Test Suite ===");

// --------------------------------------------------------------------------
// 1. Exported Constants Verification
// --------------------------------------------------------------------------
console.log("[1/10] Verifying exported constants...");
assert.equal(addon.SE_NODBIT_MEAN, 1);
assert.equal(addon.SE_NODBIT_OSCU, 2);
assert.equal(addon.SE_NODBIT_OSCU_BAR, 4);
assert.equal(addon.SE_NODBIT_FOPOINT, 256);

assert.equal(addon.SE_HOR2ECL, 0);
assert.equal(addon.SE_HOR2EQU, 1);
assert.equal(addon.SE_ECL2HOR, 0);
assert.equal(addon.SE_EQU2HOR, 1);

assert.equal(addon.SE_SIDM_USER, 255);
assert.equal(addon.SE_SIDM_FAGAN_BRADLEY, 0);
assert.equal(addon.SE_SIDM_LAHIRI, 1);
assert.equal(addon.SE_SIDM_DELUCE, 2);
assert.equal(addon.SE_SIDM_RAMAN, 3);
assert.equal(addon.SE_SIDM_USHASHASHI, 4);
assert.equal(addon.SE_SIDM_KRISHNAMURTI, 5);
assert.equal(addon.SE_SIDM_DJWHAL_KHUL, 6);
assert.equal(addon.SE_SIDM_YUKTESHWAR, 7);
assert.equal(addon.SE_SIDM_JN_BHASIN, 8);
assert.equal(addon.SE_SIDM_BABYL_KUGLER1, 9);
assert.equal(addon.SE_SIDM_BABYL_KUGLER2, 10);
assert.equal(addon.SE_SIDM_BABYL_KUGLER3, 11);
assert.equal(addon.SE_SIDM_BABYL_HUBER, 12);
assert.equal(addon.SE_SIDM_BABYL_ETPSC, 13);
assert.equal(addon.SE_SIDM_ALDEBARAN_15TAU, 14);
assert.equal(addon.SE_SIDM_HIPPARCHOS, 15);
assert.equal(addon.SE_SIDM_SASSANIAN, 16);
assert.equal(addon.SE_SIDM_GALCENT_0SAG, 17);
assert.equal(addon.SE_SIDM_J2000, 18);
assert.equal(addon.SE_SIDM_J1900, 19);
assert.equal(addon.SE_SIDM_B1950, 20);
assert.equal(addon.SE_SIDM_SURYASIDDHANTA, 21);
assert.equal(addon.SE_SIDM_SURYASIDDHANTA_MSUN, 22);
assert.equal(addon.SE_SIDM_ARYABHATA, 23);
assert.equal(addon.SE_SIDM_ARYABHATA_MSUN, 24);
assert.equal(addon.SE_SIDM_SS_REVATI, 25);
assert.equal(addon.SE_SIDM_SS_CITRA, 26);
assert.equal(addon.SE_SIDM_TRUE_CITRA, 27);
assert.equal(addon.SE_SIDM_TRUE_REVATI, 28);
assert.equal(addon.SE_SIDM_TRUE_PUSHYA, 29);
assert.equal(addon.SE_SIDM_GALCENT_RGILBRAND, 30);
assert.equal(addon.SE_SIDM_GALEQ_IAU1958, 31);
assert.equal(addon.SE_SIDM_GALEQ_TRUE, 32);
assert.equal(addon.SE_SIDM_GALEQ_MULA, 33);
assert.equal(addon.SE_SIDM_GALALIGN_MARDYKS, 34);
assert.equal(addon.SE_SIDM_TRUE_MULA, 35);
assert.equal(addon.SE_SIDM_GALCENT_MULA_WILHELM, 36);
assert.equal(addon.SE_SIDM_ARYABHATA_522, 37);
assert.equal(addon.SE_SIDM_BABYL_BRITTON, 38);
assert.equal(addon.SE_SIDM_TRUE_SHERATAN, 39);
assert.equal(addon.SE_SIDM_TRUE_ALDEBARAN, 40);
assert.equal(addon.SE_SIDM_GALEQU_FIORENZA, 41);
assert.equal(addon.SE_SIDM_VALENS_MOON, 42);

// --------------------------------------------------------------------------
// 2. Moshier Analytical Mode (JsEngine & JsSweAdapter)
// --------------------------------------------------------------------------
console.log("[2/10] Verifying Moshier analytical mode...");
const engineMoshier = addon.JsEngine.newMoshier();
assert.equal(engineMoshier.isMoshierMode(), true);

const engineSpk = addon.JsEngine.fromDataDir(DATA_DIR, 2, "embedded", false);
assert.equal(engineSpk.isMoshierMode(), false);

const engineAutoDir = addon.JsEngine.fromDataDirOrMoshier(DATA_DIR);
assert.equal(engineAutoDir.isMoshierMode(), false);

const engineAutoFallback = addon.JsEngine.fromDataDirOrMoshier(path.join(os.tmpdir(), "nonexistent_lx_dir_12345"));
assert.equal(engineAutoFallback.isMoshierMode(), true);

const jdJ2000 = { jd1: 2451545.0, jd2: 0.0 };
const sunMoshier = engineMoshier.calc("Sun", jdJ2000, "TT", null);
const sunViaCalcMoshier = engineSpk.calcMoshier("Sun", jdJ2000, "TT", null);
assertClose(sunMoshier.lonDeg, 280.368, 0.05, "Moshier Sun longitude");
assertClose(sunMoshier.lonDeg, sunViaCalcMoshier.lonDeg, 1e-10, "Moshier direct vs calcMoshier");

const sweMoshier = addon.JsSweAdapter.newMoshier();
const sweAuto = addon.JsSweAdapter.fromDataDirOrMoshier(DATA_DIR);
const sweCalcOut = sweMoshier.calcUt(2451545.0, addon.SE_SUN, null);
assert.equal(sweCalcOut.xx.length, 6);
assertClose(sweCalcOut.xx[0], 280.368, 0.05, "SweAdapter Moshier Sun longitude");

// --------------------------------------------------------------------------
// 3. Custom Fixed Stars Catalog Loading
// --------------------------------------------------------------------------
console.log("[3/10] Verifying custom fixed stars catalog...");
const tempStarFile = path.join(os.tmpdir(), `test_sefstars_${Date.now()}.txt`);
const starContent = `
PolarisTest  ,alUMi, IC, 02,31,49.0945, +89,15,50.792, 44.48, -11.85, -16.4, 7.54, 1.98, 2000, 999001
CustomVega   ,alLyr, IC, 18,36,56.3363, +38,47,01.280, 200.94, 286.23, -13.9, 130.23, 0.03, 2000, 999002
`;
fs.writeFileSync(tempStarFile, starContent.trim());

try {
  const loadedCount = engineSpk.loadCustomStars(tempStarFile);
  assert.equal(loadedCount, 2);

  const starInfo1 = engineSpk.starInfo("PolarisTest");
  assert.equal(starInfo1.hip, 999001);
  assert.equal(starInfo1.canonicalName, "PolarisTest");

  const starCalc = engineSpk.calcStarByName("PolarisTest", jdJ2000, "TT", null);
  assert.ok(starCalc.lonDeg >= 0 && starCalc.lonDeg < 360);

  const starList = engineSpk.listStars();
  assert.ok(starList.some((s) => s.canonicalName === "PolarisTest"));
  assert.ok(starList.some((s) => s.canonicalName === "CustomVega"));
} finally {
  fs.unlinkSync(tempStarFile);
}

// --------------------------------------------------------------------------
// 4. Photometric Phenomena (JsEngine & JsSweAdapter)
// --------------------------------------------------------------------------
console.log("[4/10] Verifying photometric phenomena...");
const phenoVenus = engineSpk.phenomena("Venus", jdJ2000, "TT", null);
assert.ok(phenoVenus.phaseAngleDeg >= 0.0 && phenoVenus.phaseAngleDeg <= 180.0);
assert.ok(phenoVenus.phaseIlluminated >= 0.0 && phenoVenus.phaseIlluminated <= 1.0);
assert.ok(phenoVenus.elongationDeg >= 0.0 && phenoVenus.elongationDeg <= 180.0);
assert.ok(phenoVenus.apparentDiameterDeg > 0.0);
assert.ok(phenoVenus.apparentMagnitudeV < 0.0); // Venus is very bright
assert.ok(phenoVenus.horizontalParallaxDeg > 0.0);
assert.equal(phenoVenus.attr20.length, 20);

const swePhenoEt = sweAuto.pheno(2451545.0, addon.SE_VENUS, null);
assert.equal(swePhenoEt.attr.length, 20);
assert.ok(swePhenoEt.retflag > 0);
assertClose(swePhenoEt.attr[0], phenoVenus.phaseAngleDeg, 1e-9, "swe.pheno phase_angle");
assertClose(swePhenoEt.attr[1], phenoVenus.phaseIlluminated, 1e-9, "swe.pheno phase_illuminated");

const swePhenoUt = sweAuto.phenoUt(2451545.0, addon.SE_VENUS, null);
assert.equal(swePhenoUt.attr.length, 20);
assert.ok(swePhenoUt.retflag > 0);
assertClose(swePhenoUt.attr[0], phenoVenus.phaseAngleDeg, 0.01, "swe.phenoUt phase_angle near TT");

// --------------------------------------------------------------------------
// 5. Planetary & Lunar Orbital Nodes & Apsides
// --------------------------------------------------------------------------
console.log("[5/10] Verifying orbital nodes and apsides...");
const nodMars = engineSpk.nodAps("Mars", jdJ2000, "TT", null, addon.SE_NODBIT_MEAN);
assert.equal(nodMars.ascendingNode.length, 6);
assert.equal(nodMars.descendingNode.length, 6);
assert.equal(nodMars.periapsis.length, 6);
assert.equal(nodMars.apoapsis.length, 6);
assert.ok(nodMars.ascendingNode[0] >= 0.0 && nodMars.ascendingNode[0] < 360.0);
assert.ok(nodMars.descendingNode[0] >= 0.0 && nodMars.descendingNode[0] < 360.0);

const sweNodEt = sweAuto.nodAps(2451545.0, addon.SE_MARS, null, addon.SE_NODBIT_MEAN);
assert.equal(sweNodEt.ascendingNode.length, 6);
assertClose(sweNodEt.ascendingNode[0], nodMars.ascendingNode[0], 1e-9, "swe.nodAps ascending node lon at TT");

const sweNodUt = sweAuto.nodApsUt(2451545.0, addon.SE_MARS, null, addon.SE_NODBIT_MEAN);
assert.equal(sweNodUt.ascendingNode.length, 6);
assert.equal(sweNodUt.descendingNode.length, 6);
assert.equal(sweNodUt.periapsis.length, 6);
assert.equal(sweNodUt.apoapsis.length, 6);
assertClose(sweNodUt.ascendingNode[0], nodMars.ascendingNode[0], 0.01, "swe.nodApsUt ascending node lon near TT");

// --------------------------------------------------------------------------
// 6. Lunar Occultations (Global & Local Circumstances)
// --------------------------------------------------------------------------
console.log("[6/10] Verifying lunar occultations...");
const jd2015 = { jd1: 2457023.5, jd2: 0.0 }; // 2015-01-01 UT
const obsGreenwich = { latitudeDeg: 51.4769, longitudeDeg: 0.0, altitudeM: 0.0 };

const occGlob = engineSpk.nextLunarOccultation("Aldebaran", jd2015, "UT1", false, 60.0);
assert.ok(occGlob.jdUt1Max > 2457023.5);
assert.ok(occGlob.jdUt1Max < 2457023.5 + 60.0);
assert.equal(occGlob.tret.length, 10);
assert.ok(occGlob.tret[1] <= occGlob.tret[0]);
assert.ok(occGlob.tret[4] >= occGlob.tret[0]);

const occLoc = engineSpk.nextLocalLunarOccultation("Aldebaran", obsGreenwich, jd2015, "UT1", false, 400.0);
assert.ok(occLoc.jdUt1Max > 2457023.5);
assert.equal(occLoc.tret.length, 10);
assert.ok(occLoc.moonAltitudeDeg !== null);

const occListGlob = engineSpk.lunarOccultationsInRange(
  "Aldebaran",
  { jd1: 2457023.5, jd2: 0.0 },
  { jd1: 2457100.0, jd2: 0.0 },
  "UT1"
);
assert.ok(occListGlob.length >= 1);
assert.equal(occListGlob[0].target.includes("21421") || occListGlob[0].target.includes("Aldebaran"), true);

const occListLoc = engineSpk.localLunarOccultationsInRange(
  "Aldebaran",
  obsGreenwich,
  { jd1: 2457023.5, jd2: 0.0 },
  { jd1: 2457400.0, jd2: 0.0 },
  "UT1"
);
assert.ok(occListLoc.length >= 1);

const sweOccGlob = sweAuto.lunOccultWhenGlob(2457023.5, 0, "Aldebaran", 0, false);
assert.ok(sweOccGlob.retflag > 0);
assert.equal(sweOccGlob.tret.length, 10);
assertClose(sweOccGlob.tret[0], occGlob.tret[0], 1e-4, "swe.lunOccultWhenGlob vs nextLunarOccultation");

const sweOccLoc = sweAuto.lunOccultWhenLoc(2457023.5, 0, "Aldebaran", 0, [0.0, 51.4769, 0.0], false);
assert.ok(sweOccLoc.retflag > 0);
assert.equal(sweOccLoc.tret.length, 10);
assert.equal(sweOccLoc.attr.length, 20);

// --------------------------------------------------------------------------
// 7. 3D Fractional Mundane House Position
// --------------------------------------------------------------------------
console.log("[7/10] Verifying 3D mundane house position...");
const armc = 100.0;
const geolat = 51.5;
const eps = 23.44;
const lonPoint = 120.0;
const latPoint = 5.0;

const hposEngine = engineSpk.housePos(armc, geolat, eps, "P", lonPoint, latPoint);
assert.ok(hposEngine >= 1.0 && hposEngine < 13.0);

const hposSwe = sweAuto.housePos(armc, geolat, eps, [lonPoint, latPoint], "P");
assertClose(hposSwe, hposEngine, 1e-9, "JsEngine.housePos vs JsSweAdapter.housePos");

// --------------------------------------------------------------------------
// 8. Equation of Time (EoT)
// --------------------------------------------------------------------------
console.log("[8/10] Verifying Equation of Time...");
// Around Feb 11 (JD 2451586.0), EoT in minutes is ~ -14.2 minutes
const jdFeb11 = { jd1: 2451586.0, jd2: 0.0 };
const eotFeb = engineSpk.timeEquation(jdFeb11, "UT1");
const eotFebMin = eotFeb * 1440.0;
assertClose(eotFebMin, -14.2, 0.6, "Feb 11 Equation of Time (minutes)");

const sweEotFeb = sweAuto.timeEqu(2451586.0);
assertClose(sweEotFeb * 1440.0, eotFebMin, 1e-6, "JsEngine vs JsSweAdapter EoT");

// Around Nov 3 (JD 2451851.5), EoT in minutes is ~ +16.4 minutes
const jdNov3 = { jd1: 2451851.5, jd2: 0.0 };
const eotNov = engineSpk.timeEquation(jdNov3, "TT");
assertClose(eotNov * 1440.0, 16.4, 0.6, "Nov 3 Equation of Time (minutes)");

// --------------------------------------------------------------------------
// 9. Inverse Horizontal Coordinates (azalt_rev)
// --------------------------------------------------------------------------
console.log("[9/10] Verifying inverse horizontal coordinate transformation...");
const obsZurich = { latitudeDeg: 47.3769, longitudeDeg: 8.5417, altitudeM: 408.0 };
const geoposZurich = [8.5417, 47.3769, 408.0];

// Round-trip equatorial
const raIn = 142.5;
const decIn = 21.3;
const hz = engineSpk.azalt(raIn, decIn, jdJ2000, obsZurich, "UT1", 1013.25, 15.0);
const revEqu = engineSpk.azaltRev(obsZurich, jdJ2000, "UT1", hz.azimuthDeg, hz.altitudeTrueDeg, "equatorial");
assertClose(revEqu[0], raIn, 1e-7, "azaltRev round-trip RA");
assertClose(revEqu[1], decIn, 1e-7, "azaltRev round-trip Dec");

// Round-trip ecliptic via SweAdapter
const lonIn = 125.0;
const latIn = 4.0;
const xaz = sweAuto.azalt(2451545.0, addon.SE_ECL2HOR, geoposZurich, 0.0, 0.0, [lonIn, latIn]);
const xout = sweAuto.azaltRev(2451545.0, addon.SE_HOR2ECL, geoposZurich, [xaz[0], xaz[1]]);
assertClose(xout[0], lonIn, 1e-7, "swe.azaltRev round-trip Lon");
assertClose(xout[1], latIn, 1e-7, "swe.azaltRev round-trip Lat");

// --------------------------------------------------------------------------
// 10. Extended Sidereal Ayanamsa & Error Handling
// --------------------------------------------------------------------------
console.log("[10/10] Verifying sidereal ayanamsa models and error handling...");
// All 47 legacy models (0..=46)
for (let m = 0; m <= 46; m++) {
  sweAuto.setSidMode(m, null, null);
  const val = sweAuto.getAyanamsaUt(2451545.0);
  assert.ok(Number.isFinite(val) && val >= -360.0 && val <= 360.0);
}

// User-defined ayanamsa mode (SE_SIDM_USER = 255)
const t0User = 2451545.0;
const ayanT0User = 24.123456;
sweAuto.setSidMode(addon.SE_SIDM_USER, t0User, ayanT0User);
const userAyaEt = sweAuto.getAyanamsa(t0User);
assertClose(userAyaEt, ayanT0User, 1e-12, "SE_SIDM_USER exact t0 value");

const userAyaUt = sweAuto.getAyanamsaExUt(t0User, 64); // SEFLG_NONUT = 64
assertClose(userAyaUt.ayanamsa, ayanT0User, 1e-4, "SE_SIDM_USER UT1 value");

// Error handling tests
assert.throws(() => {
  engineSpk.phenomena("NonExistentBodyXYZ", jdJ2000, "TT", null);
}, /UnsupportedBody/);

assert.throws(() => {
  sweAuto.housePos(armc, geolat, eps, [120.0], "P"); // xpin length < 2
}, /xpin must have at least 2 elements/);

assert.throws(() => {
  sweAuto.azaltRev(2451545.0, addon.SE_HOR2ECL, [8.5417, 47.3769], [180.0, 45.0]); // geopos length < 3
}, /geopos must have at least 3 elements/);

// --------------------------------------------------------------------------
// 11. Wave F Astronomical Event Solvers & Brady Parans
// --------------------------------------------------------------------------
console.log("[11/11] Verifying Wave F astronomical event solvers and Brady parans...");

const jd2024Jan1 = addon.julday(2024, 1, 1, 0.0, "gregorian");
const jd2024End = jd2024Jan1 + 366.0;

// 11.1 searchBodyLongitude (Vernal Equinox 2024 ~ 0 deg Aries ingress)
const vernal2024 = engineSpk.searchBodyLongitude("Sun", jd2024Jan1, 0.0, true, 120.0, null);
assert.equal(vernal2024.body, "Sun");
assertClose(vernal2024.targetLonDeg, 0.0, 1e-12, "Vernal targetLonDeg");
assert.ok(
  Math.abs(vernal2024.actualLonDeg - 0.0) < 1e-6 ||
    Math.abs(vernal2024.actualLonDeg - 360.0) < 1e-6,
  `Vernal actualLonDeg=${vernal2024.actualLonDeg}`
);
assert.equal(vernal2024.isRetrograde, false);
assert.ok(vernal2024.distKm > 147_000_000.0 && vernal2024.distKm < 153_000_000.0);

// 11.2 bodyLongitudeCrossingsInRange (Mercury retrograde triple pass in 2024)
const mercCrossings = engineSpk.bodyLongitudeCrossingsInRange(
  "Mercury",
  jd2024Jan1,
  jd2024End,
  20.0,
  null
);
assert.ok(mercCrossings.length >= 3, `Expected >= 3 Mercury crossings, got ${mercCrossings.length}`);
assert.ok(mercCrossings.some((c) => c.isRetrograde === true));
assert.ok(mercCrossings.some((c) => c.isRetrograde === false));

// 11.3 solarTermsInRange (24 Chinese Solar Terms in 2024)
const solarTerms2024 = engineSpk.solarTermsInRange(jd2024Jan1, jd2024End);
assert.equal(solarTerms2024.length, 24);
const chunfen = solarTerms2024.find((t) => t.pinyinName === "Chunfen");
assert.ok(chunfen);
assert.equal(chunfen.termIndex, 0);
assert.equal(chunfen.isMajorTerm, true);
assertClose(chunfen.targetLonDeg, 0.0, 1e-12, "Chunfen targetLonDeg");
assertClose(chunfen.jdUt1, vernal2024.jdUt1, 1e-6, "Chunfen vs searchBodyLongitude(Sun, 0)");

// 11.4 seasonsInYear
const seasons2024 = engineSpk.seasonsInYear(2024);
assert.equal(seasons2024.year, 2024);
assertClose(seasons2024.vernalEquinoxJdUt1, vernal2024.jdUt1, 1e-6, "seasonsInYear vernalEquinoxJdUt1");
assert.ok(seasons2024.summerSolsticeJdUt1 > seasons2024.vernalEquinoxJdUt1);
assert.ok(seasons2024.autumnalEquinoxJdUt1 > seasons2024.summerSolsticeJdUt1);
assert.ok(seasons2024.winterSolsticeJdUt1 > seasons2024.autumnalEquinoxJdUt1);

// 11.5 humanDesignEpoch (-88 deg solar arc)
const hdEpoch = engineSpk.humanDesignEpoch(2451545.0);
assert.equal(hdEpoch.body, "Sun");
const sunAtBirth = engineSpk.calc("Sun", { jd1: 2451545.0, jd2: 0.0 }, "UT1", null);
const expectedDesignLon = ((sunAtBirth.lonDeg - 88.0) % 360.0 + 360.0) % 360.0;
assertClose(hdEpoch.actualLonDeg, expectedDesignLon, 1e-6, "Human Design -88 deg solar longitude");
assert.ok(2451545.0 - hdEpoch.jdUt1 > 85.0 && 2451545.0 - hdEpoch.jdUt1 < 92.0);

// 11.6 searchLunarPhase
const nextFullMoon = engineSpk.searchLunarPhase(jd2024Jan1, 180.0, true);
assert.equal(nextFullMoon.kind, "FullMoon");
assertClose(nextFullMoon.targetElongationDeg, 180.0, 1e-9, "FullMoon targetElongationDeg");
assert.ok(nextFullMoon.moonDistKm > 350_000.0 && nextFullMoon.moonDistKm < 410_000.0);

const customPhase = engineSpk.searchLunarPhase(jd2024Jan1, 45.0, true);
assert.equal(customPhase.kind, undefined);
assertClose(customPhase.targetElongationDeg, 45.0, 1e-9, "Custom 45 deg lunar elongation");

// 11.7 lunarPhasesInRange
const janPhases = engineSpk.lunarPhasesInRange(jd2024Jan1, jd2024Jan1 + 32.0);
assert.ok(janPhases.length >= 4);
const kinds = new Set(janPhases.map((p) => p.kind));
assert.ok(kinds.has("NewMoon") && kinds.has("FirstQuarter") && kinds.has("FullMoon") && kinds.has("LastQuarter"));

// 11.8 lunarPhaseSummary (8-phase & Kerykeion 28-phase replacement)
const phaseSummary = engineSpk.lunarPhaseSummary(2451545.0);
assert.ok(phaseSummary.elongationDeg >= 0.0 && phaseSummary.elongationDeg < 360.0);
assert.ok(phaseSummary.illuminationFraction >= 0.0 && phaseSummary.illuminationFraction <= 1.0);
assert.ok(phaseSummary.phase8Index >= 0 && phaseSummary.phase8Index <= 7);
assert.equal(typeof phaseSummary.phase8Name, "string");
assert.ok(phaseSummary.phase28Index >= 1 && phaseSummary.phase28Index <= 28);
assert.ok(phaseSummary.ageDays >= 0.0 && phaseSummary.ageDays <= 30.0);
assert.ok(phaseSummary.synodicFraction >= 0.0 && phaseSummary.synodicFraction <= 1.0);
assert.ok(phaseSummary.prevNewMoonJdUt1 <= 2451545.0);
assert.ok(phaseSummary.nextNewMoonJdUt1 >= 2451545.0);

// 11.9 searchStationsInRange
const mercStations = engineSpk.searchStationsInRange("Mercury", jd2024Jan1, jd2024End, null);
assert.ok(mercStations.length >= 6);
assert.ok(mercStations.some((s) => s.kind === "StationRetrograde"));
assert.ok(mercStations.some((s) => s.kind === "StationDirect"));

// 11.10 searchDeclinationExtremaInRange
const moonExtrema = engineSpk.searchDeclinationExtremaInRange(
  "Moon",
  jd2024Jan1,
  jd2024Jan1 + 30.0,
  null
);
assert.ok(moonExtrema.length >= 2);
assert.ok(moonExtrema.some((e) => e.kind === "NorthernMaximum" && e.decDeg > 18.0));
assert.ok(moonExtrema.some((e) => e.kind === "SouthernMinimum" && e.decDeg < -18.0));

// 11.11 findParans (Bernadette Brady diurnal fixed-star parans)
const parans = engineSpk.findParans(
  2451545.0,
  obsGreenwich,
  12.0,
  ["Sun", "Moon", "Mercury", "Venus", "Mars"],
  ["Regulus", "Sirius", "Spica", "Aldebaran", "Polaris"]
);
assert.ok(Array.isArray(parans));
for (let i = 0; i < parans.length; i++) {
  const p = parans[i];
  assert.ok(p.timeOrbMinutes >= 0.0 && p.timeOrbMinutes <= 12.0);
  assert.equal(p.paranType, p.starAngle);
  assert.ok(["rising", "culminating", "setting", "anti-culminating"].includes(p.starAngle));
  assert.ok(["rising", "culminating", "setting", "anti-culminating"].includes(p.planetAngle));
  if (i > 0) {
    assert.ok(p.timeOrbMinutes >= parans[i - 1].timeOrbMinutes - 1e-12);
  }
}

// --------------------------------------------------------------------------
// 12. Wave H Zero-Threshold Onboarding, Almanac, Aspects & ChartSnapshot
// --------------------------------------------------------------------------
console.log("[12/12] Verifying Wave H zero-arg Engine, utcFromLocal (IANA), aspects, Panchanga & chartSnapshot...");

// 12.1 Zero-threshold constructor new Engine() and new JsEngine()
assert.equal(addon.Engine, addon.JsEngine);
assert.equal(addon.SweAdapter, addon.JsSweAdapter);
const zeroEngine = new addon.Engine();
const sunZero = zeroEngine.calc("Sun", { jd1: 2451545.0, jd2: 0.0 }, "TT", null);
assertClose(sunZero.lonDeg, 280.36816528, 1e-3, "new Engine() Sun longitude at J2000");

// 12.2 utcFromLocal (numeric offset + IANA timezone strings via Intl.DateTimeFormat)
const jdNumeric = addon.utcFromLocal(2000, 1, 1, 13, 0, 0.0, 1.0, true);
assertClose(jdNumeric, 2451545.0, 1e-9, "utcFromLocal numeric +1.0h");
const jdBerlinWinter = addon.utcFromLocal(2000, 1, 1, 13, 0, 0.0, "Europe/Berlin", true);
assertClose(jdBerlinWinter, 2451545.0, 1e-9, "utcFromLocal Europe/Berlin winter CET (+1h)");
const jdBerlinSummer = addon.utcFromLocal(2024, 7, 15, 14, 0, 0.0, "Europe/Berlin", true);
const jdExpectedSummer = addon.julday(2024, 7, 15, 12.0, "gregorian");
assertClose(jdBerlinSummer, jdExpectedSummer, 1e-9, "utcFromLocal Europe/Berlin summer CEST (+2h)");
const jdKolkata = addon.utcFromLocal(2024, 1, 1, 17, 30, 0.0, "Asia/Kolkata", true);
const jdExpectedKolkata = addon.julday(2024, 1, 1, 12.0, "gregorian");
assertClose(jdKolkata, jdExpectedKolkata, 1e-9, "utcFromLocal Asia/Kolkata IST (+5.5h)");

// 12.3 planetaryReturn & planetaryReturnsInRange
const solarReturn = engineSpk.planetaryReturn("Sun", 2451545.0, 2451545.0 + 350.0, null);
assert.equal(solarReturn.body, "Sun");
assert.ok(solarReturn.jdUt1 > 2451545.0 + 364.0 && solarReturn.jdUt1 < 2451545.0 + 367.0);
const solarReturns2Y = engineSpk.planetaryReturnsInRange(
  "Sun",
  2451545.0,
  2451545.0 + 300.0,
  2451545.0 + 800.0,
  null
);
assert.equal(solarReturns2Y.length, 2);

// 12.4 Aspects: calcChartAspects, calcSynastryAspects, calcDeclinationAspects, searchExactAspect, mutualAspectsInRange, natalTransitsInRange, voidOfCourseMoonInRange
const chartAspects = engineSpk.calcChartAspects(2451545.0, null, "major", 1.0, null);
assert.ok(Array.isArray(chartAspects) && chartAspects.length > 0);
assert.equal(typeof chartAspects[0].isApplying, "boolean");

const synastryAspects = engineSpk.calcSynastryAspects(
  2451545.0,
  ["Sun", "Moon", "Venus"],
  2451545.0 + 180.0,
  ["Sun", "Moon", "Mars"],
  "major",
  1.0,
  null
);
assert.ok(Array.isArray(synastryAspects));

const decAspects = engineSpk.calcDeclinationAspects(2451545.0, null, 1.5, null);
assert.ok(Array.isArray(decAspects));

const exactConj = engineSpk.searchExactAspect("Moon", "Sun", 0.0, jd2024Jan1, true, null);
assert.equal(exactConj.body1, "Moon");
assert.equal(exactConj.body2, "Sun");
assertClose(exactConj.targetAngleDeg, 0.0, 1e-9, "searchExactAspect targetAngleDeg");

const mutualAspects = engineSpk.mutualAspectsInRange(
  jd2024Jan1,
  jd2024Jan1 + 5.0,
  ["Sun", "Moon", "Mercury"],
  "major",
  null
);
assert.ok(Array.isArray(mutualAspects) && mutualAspects.length > 0);

const natalTransits = engineSpk.natalTransitsInRange(
  2451545.0,
  jd2024Jan1,
  jd2024Jan1 + 5.0,
  ["Moon"],
  ["Sun", "Moon"],
  "major",
  null
);
assert.ok(Array.isArray(natalTransits) && natalTransits.length > 0);

const vocPeriods = engineSpk.voidOfCourseMoonInRange(jd2024Jan1, jd2024Jan1 + 7.0, null);
assert.ok(Array.isArray(vocPeriods) && vocPeriods.length >= 2);
assert.ok(vocPeriods[0].endJdUt1 > vocPeriods[0].startJdUt1);

// 12.5 Vedic Panchanga, Vimshottari Dasha & Vargas
const obsDelhi = { latitudeDeg: 28.6139, longitudeDeg: 77.2090, altitudeM: 216.0 };
const panch = engineSpk.panchanga(2451545.0, obsDelhi, "lahiri");
assert.ok(panch.tithi.tithiIndex >= 1 && panch.tithi.tithiIndex <= 30);
assert.ok(panch.nakshatra.nakshatraIndex >= 1 && panch.nakshatra.nakshatraIndex <= 27);
assert.ok(panch.yoga.yogaIndex >= 1 && panch.yoga.yogaIndex <= 27);
assert.ok(panch.karana.karanaIndex >= 1 && panch.karana.karanaIndex <= 60);

const dasha = engineSpk.vimshottariDasha(2451545.0, "lahiri", 365.2425);
assert.equal(dasha.periods.length, 9);
assert.ok(dasha.periods[0].subPeriods.length >= 1 && dasha.periods[0].subPeriods.length <= 9);
assert.equal(dasha.periods[1].subPeriods.length, 9);

const d9 = engineSpk.calcVarga(45.0, 9);
assert.equal(d9.division, 9);
assert.equal(d9.vargaCode, "D9");
const shodasha = engineSpk.calcShodashavarga(45.0);
assert.equal(shodasha.length, 16);

// 12.6 Traditional Sect, Planetary Hours, Hermetic Lots, Zodiac Placement & ChartSnapshot
const sect = engineSpk.chartSect(2451545.0, obsGreenwich);
assert.equal(typeof sect.isDayChart, "boolean");

const ph = engineSpk.planetaryHours(2451545.0, obsGreenwich);
assert.equal(ph.hours.length, 24);
assert.equal(ph.hours[0].lord, ph.dayLord);

const lots = engineSpk.calcLots(2451545.0, obsGreenwich, { system: "P" });
assert.equal(lots.length, 7);
assert.equal(lots[0].name, "Fortune");
assert.equal(lots[1].name, "Spirit");

const zp = addon.mapZodiacPlacement(302.0);
assert.equal(zp.signName, "Aquarius");
assert.equal(zp.hdGate, 41);
assert.equal(zp.hdLine, 1);

const snapshot = engineSpk.chartSnapshot(2451545.0, obsGreenwich, null, { system: "P" }, null, "major", 1.0);
assert.equal(snapshot.bodies.length, 13);
assert.equal(snapshot.lots.length, 7);
assert.equal(snapshot.houses.cuspsDeg.length, 12);
assert.ok(Array.isArray(snapshot.aspects) && snapshot.aspects.length > 0);

const emptySnapshot = engineSpk.chartSnapshot(2451545.0, obsGreenwich, [], { system: "P" }, null, "major", 1.0);
assert.equal(emptySnapshot.bodies.length, 0);

// 12.7 Round 4 (gpt-6-astra) Node.js regressions
assert.equal(addon.SE_INTP_PERG, 22);
assert.equal(addon.SE_SIDM_LAHIRI_ICRC, 46);
const pergCalc = sweAuto.calcUt(2451545.0, addon.SE_INTP_PERG, addon.SEFLG_SWIEPH | addon.SEFLG_SPEED);
assert.ok(pergCalc.xx[0] >= 90.0 && pergCalc.xx[0] <= 94.0);
sweAuto.setSidMode(addon.SE_SIDM_LAHIRI_ICRC, null, null);
const ayaEx = sweAuto.getAyanamsaEx(2451545.0, addon.SEFLG_SWIEPH);
assert.ok(ayaEx.ayanamsa >= 23.0 && ayaEx.ayanamsa <= 25.0);
const sunICusps = sweAuto.housesArmc(120.0, 51.5, 23.4392911, "i", 10.0);
assert.equal(sunICusps.cusps.length, 12);

const tTotal = (performance.now() - tStart).toFixed(1);
console.log(`\nAll Wave B, C, F & H Node.js SDK tests PASSED in ${tTotal} ms!`);
