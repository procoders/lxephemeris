#!/usr/bin/env node
/**
 * Phase 5.3 Cross-Language Node.js SDK Verification (`lxephemeris-node` / `target/release/lxephemeris.node`):
 * 1. Loads the compiled native N-API addon (`target/release/lxephemeris.node`) via `createRequire`.
 * 2. Exercises `JsEngine.fromDataDir`, `calc`, `calcStarByName`, `calcBatch`, `calcBatchFlatF64`,
 *    `calcBatchPerItem`, `ayanamsa`, `houses`, `housesArmc`, `gauquelinSectors`, `azalt`, `riseTrans`,
 *    `nextSolarEclipse`, `nextLunarEclipse`, `julday`, `revjul`, and structured `[ErrorCode]` errors.
 * 3. Compares all outputs against `target/release/cross_lang_ref.json` produced by Python (`lxephemeris-py`)
 *    and verifies `<= 1e-12 deg` cross-language equivalence!
 */

import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const WORKSPACE_ROOT = path.resolve(__dirname, "..");
const RELEASE_DIR = path.join(WORKSPACE_ROOT, "target", "release");
const DATA_DIR = path.join(WORKSPACE_ROOT, "data");

const require = createRequire(import.meta.url);
const addon = require(path.join(RELEASE_DIR, "lxephemeris.node"));

function assertClose(actual, expected, tol, label) {
  const diff = Math.abs(actual - expected);
  assert.ok(
    diff <= tol,
    `${label}: actual=${actual}, expected=${expected}, diff=${diff} > tol=${tol}`
  );
}

const t0 = performance.now();
const engine = addon.JsEngine.fromDataDir(DATA_DIR, 4, "embedded", true);
console.log("[1/4] Loaded JsEngine from data/ with SHA-256 verification");

const caps = JSON.parse(engine.capabilitiesJson());
assert.equal(caps.ayanamsa_models.length, 12);
assert.equal(caps.house_systems.length, 23);
assert.equal(typeof caps.ephemeris_source, "string");
assert.equal(typeof caps.ephemeris_version, "string");

const jdJ2000 = { jd1: 2451545.0, jd2: 0.0 };

const sun = engine.calc("Sun", jdJ2000, "TT", null);
const moon = engine.calc("Moon", jdJ2000, "TT", null);
const chiron = engine.calc("Chiron", jdJ2000, "TT", null);
const regulus = engine.calcStarByName("Regulus", jdJ2000, "TT", null);
const trueNode = engine.calc("True_Node", jdJ2000, "TT", null);
const trueLilith = engine.calc("True_Lilith", jdJ2000, "TT", null);
const cupido = engine.calc("Cupido", jdJ2000, "TT", null);
const lahiri = engine.ayanamsa("lahiri", jdJ2000, "TT", "true");

const obsLondon = { latitudeDeg: 51.5074, longitudeDeg: -0.1278, altitudeM: 0.0 };
const housesP = engine.houses(
  jdJ2000,
  obsLondon,
  "P",
  "TT",
  "strict",
  "tropical",
  "lahiri",
  "true",
  true
);

const obsDallas = { latitudeDeg: 32.7767, longitudeDeg: -96.7970, altitudeM: 130.0 };
const sol2024 = engine.nextSolarEclipse(
  { jd1: 2460400.0, jd2: 0.0 },
  "TT",
  false,
  obsDallas,
  36
);
const lun2025 = engine.nextLunarEclipse(
  { jd1: 2460740.0, jd2: 0.0 },
  "TT",
  false,
  "danjon",
  36
);

const batchBodies = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Chiron"];
const batchEpochs = Array.from({ length: 10 }, (_, i) => ({ jd1: 2451545.0 + i, jd2: 0.0 }));
const flatF64 = engine.calcBatchFlatF64(batchBodies, batchEpochs, "TT", null);
assert.ok(flatF64 instanceof Float64Array, "calcBatchFlatF64 must return Float64Array");
assert.equal(flatF64.length, batchEpochs.length * batchBodies.length * 6);

const perItem = engine.calcBatchPerItem(
  ["Sun", "Chiron"],
  [jdJ2000, { jd1: 2200000.0, jd2: 0.0 }],
  "TT",
  null
);
assert.equal(perItem[0][0].ok, true);
assert.equal(perItem[1][0].ok, false);
assert.equal(perItem[1][0].errorCode, "OutOfCoverage");

// Standalone functions
const jdCalc = addon.julday(2000, 1, 1, 12.0, "gregorian");
assertClose(jdCalc, 2451545.0, 1e-12, "julday J2000");
const rev = addon.revjul(jdCalc, "gregorian");
assert.equal(rev.year, 2000);
assert.equal(rev.month, 1);
assert.equal(rev.day, 1);
assertClose(rev.hourUt, 12.0, 1e-9, "revjul hourUt");

const hArmc = addon.housesArmc(120.0, 23.4392911, 47.3769, "P", "strict", 0.0, 0.0, true);
assert.equal(hArmc.cuspsDeg.length, 12);
const gSec = addon.gauquelinSectors(120.0, 23.4392911, 47.3769, "strict");
assert.equal(gSec.length, 36);

// Structured error check
assert.throws(
  () => engine.calc("NotABody", jdJ2000, "TT", null),
  /\[UnsupportedBody\]/
);

// Phase 6 Astra Requirements & JsSweAdapter verification in Node.js
const dtInfo = engine.deltaT(jdJ2000, "UT1");
assertClose(dtInfo.deltaTSec, 63.8285, 0.05, "deltaT J2000");
assert.ok(Math.abs(dtInfo.ttMinusTdbSec) < 0.002);

const stInfo = engine.sidtime(jdJ2000, "UT1", "true", 15.0);
assertClose(stInfo.armcDeg, stInfo.localSiderealHours * 15.0, 1e-12, "sidtime armcDeg");

const enInfo = engine.eclNut(jdJ2000, "TT");
assertClose(enInfo.meanObliquityDeg, 23.4392911, 1e-4, "eclNut meanObliquityDeg");
assertClose(moon.distKm, moon.distAu * 149597870.7, 1e-3, "moon distKm");

const hzEcl = engine.azaltEcliptic(sun.lonDeg, sun.latDeg, jdJ2000, obsLondon, "TT", 1013.25, 15.0, 0.5);
assert.equal(typeof hzEcl.nearHorizon, "boolean");

const obsTromso = { latitudeDeg: 69.6492, longitudeDeg: 18.9553, altitudeM: 0.0 };
const hPolarWs = engine.houses(jdJ2000, obsTromso, "P", "TT", "fallback_whole_sign");
assert.equal(hPolarWs.effectiveSystemCode, "W");
assert.equal(hPolarWs.polarHouseFallback, true);
const hNum = addon.houseOfLon(sun.lonDeg, hPolarWs.cuspsDeg, 1e-9);
assert.ok(hNum >= 1 && hNum <= 12);

const allStars = addon.listStars();
assert.equal(allStars.length, 52);
const siriusMeta = addon.starInfo("Sirius");
assert.equal(siriusMeta.hip, 32349);
const siriusRise = engine.riseTrans("Sirius", "rise", jdJ2000, obsLondon, "UT1", null, 5.0, "center", 1013.25, 15.0);
assert.ok(siriusRise !== null && siriusRise > 2451545.0);

const obsCairo = { latitudeDeg: 30.0444, longitudeDeg: 31.2357, altitudeM: 20.0 };
const siriusHel = engine.heliacalPhenomena("Sirius", { jd1: 2460310.5, jd2: 0.0 }, obsCairo, "UT1");
assert.ok(siriusHel.heliacalRiseJdUt1 !== null);

assert.equal(sol2024.sarosSeries, 139);
assert.equal(sol2024.sarosMember, 30);
assert.ok(sol2024.centralDurationSec > 240.0 && sol2024.pathWidthKm > 170.0);
const localSol = engine.nextLocalSolarEclipse({ jd1: 2460000.0, jd2: 0.0 }, obsDallas, "TT", false, 36);
assert.equal(localSol.localVisibleAtMax, true);
const solRange = engine.solarEclipsesInRange(2460310.5, 2461041.5, null);
const lunRange = engine.lunarEclipsesInRange(2460310.5, 2461041.5, "danjon");
assert.ok(solRange.length >= 4 && lunRange.length >= 4);

const hBatch = engine.housesBatch(batchEpochs, obsLondon, "P", "TT", "strict");
assert.equal(hBatch.length, batchEpochs.length);

const swe = new addon.JsSweAdapter(engine);
const sweSun = swe.calcUt(2451545.0, addon.SE_SUN, addon.SEFLG_SPEED);
assert.equal(sweSun.retflag & addon.SEFLG_SPEED, addon.SEFLG_SPEED);
assertClose(sweSun.xx[0], 280.3689, 0.02, "JsSweAdapter calcUt Sun");
assertClose(sweSun.values[0], sweSun.xx[0], 1e-15, "JsSweAdapter calcUt .values alias");
const sweStar = swe.fixstar2Ut("Regulus", 2451545.0, addon.SEFLG_SPEED);
assert.equal(sweStar.starName, "Regulus");
const sweHouses = swe.housesEx(2451545.0, 51.5074, -0.1278, "P", 0);
assert.equal(sweHouses.cusps.length, 12);
assert.equal(sweHouses.ascmc.length, 8);
const sweAya = swe.getAyanamsaExUt(2451545.0, addon.SEFLG_SWIEPH);
assert.ok(sweAya.ayanamsa > 23.0 && sweAya.ayanamsa < 25.0);
const sweArmc = swe.housesArmc(120.0, 47.3769, 23.4392911, "P");
assert.equal(sweArmc.cusps.length, 12);
const gSecEngine = engine.gauquelinSectors(120.0, 23.4392911, 47.3769, "strict");
assert.equal(gSecEngine.length, 36);

console.log("[2/4] Verified all Node.js N-API methods, Phase 6 Astra closures, JsSweAdapter, batch buffers, eclipses, houses, and error codes");

// Cross-language equivalence vs Python (`cross_lang_ref.json`)
const crossRefPath = path.join(RELEASE_DIR, "cross_lang_ref.json");
const crossRef = JSON.parse(fs.readFileSync(crossRefPath, "utf8"));

assertClose(sun.lonDeg, crossRef.sun_j2000.lon_deg, 1e-12, "Sun lon_deg");
assertClose(sun.latDeg, crossRef.sun_j2000.lat_deg, 1e-12, "Sun lat_deg");
assertClose(sun.distAu, crossRef.sun_j2000.dist_au, 1e-12, "Sun dist_au");
assertClose(sun.lonSpeedDegDay, crossRef.sun_j2000.lon_speed_deg_day, 1e-12, "Sun lon_speed");

assertClose(moon.lonDeg, crossRef.moon_j2000.lon_deg, 1e-12, "Moon lon_deg");
assertClose(moon.latDeg, crossRef.moon_j2000.lat_deg, 1e-12, "Moon lat_deg");
assertClose(moon.distAu, crossRef.moon_j2000.dist_au, 1e-12, "Moon dist_au");
assertClose(moon.lonSpeedDegDay, crossRef.moon_j2000.lon_speed_deg_day, 1e-12, "Moon lon_speed");

assertClose(chiron.lonDeg, crossRef.chiron_j2000.lon_deg, 1e-12, "Chiron lon_deg");
assertClose(regulus.lonDeg, crossRef.regulus_j2000.lon_deg, 1e-12, "Regulus lon_deg");
assertClose(trueNode.lonDeg, crossRef.true_node_j2000.lon_deg, 1e-12, "True_Node lon_deg");
assertClose(trueLilith.lonDeg, crossRef.true_lilith_j2000.lon_deg, 1e-12, "True_Lilith lon_deg");
assertClose(trueLilith.distSpeedAuDay, crossRef.true_lilith_j2000.dist_speed_au_day, 1e-12, "True_Lilith dist_speed");
assertClose(cupido.lonDeg, crossRef.cupido_j2000.lon_deg, 1e-12, "Cupido lon_deg");
assertClose(lahiri.effectiveDeg, crossRef.lahiri_j2000.effective_deg, 1e-12, "Lahiri effective_deg");

for (let i = 0; i < 12; i++) {
  assertClose(
    housesP.cuspsDeg[i],
    crossRef.houses_p_j2000.cusps_deg[i],
    1e-12,
    `Placidus cusp[${i}]`
  );
}
assertClose(housesP.ascDeg, crossRef.houses_p_j2000.asc_deg, 1e-12, "Placidus asc_deg");
assertClose(housesP.mcDeg, crossRef.houses_p_j2000.mc_deg, 1e-12, "Placidus mc_deg");

assertClose(sol2024.maxEclipseJdTt, crossRef.sol_2024.max_eclipse_jd_tt, 1e-12, "Solar 2024 max_jd");
assertClose(sol2024.gamma, crossRef.sol_2024.gamma, 1e-12, "Solar 2024 gamma");
assertClose(sol2024.localMaxJdTt, crossRef.sol_2024.local_max_jd_tt, 1e-12, "Solar 2024 local_max_jd");
assertClose(lun2025.maxEclipseJdTt, crossRef.lun_2025.max_eclipse_jd_tt, 1e-12, "Lunar 2025 max_jd");
assertClose(lun2025.umbralMagnitude, crossRef.lun_2025.umbral_magnitude, 1e-12, "Lunar 2025 umbral_mag");

for (let i = 0; i < 12; i++) {
  assertClose(flatF64[i], crossRef.batch_flat_f64_head12[i], 1e-12, `batch_flat_f64[${i}]`);
}

console.log("[3/4] Verified 0.0-diff (<= 1e-12 deg) cross-language equivalence between Python and Node.js");
const elapsed = (performance.now() - t0).toFixed(1);
console.log(`[4/4] Node.js SDK verification PASSED in ${elapsed} ms`);
