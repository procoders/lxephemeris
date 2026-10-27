const { existsSync } = require('node:fs');
const { join } = require('node:path');

const platform = process.platform;
const arch = process.arch;

const candidates = [
  join(__dirname, `lxephemeris.${platform}-${arch}.node`),
  join(__dirname, `lxephemeris.${platform}-${arch}-gnu.node`),
  join(__dirname, `lxephemeris.${platform}-${arch}-msvc.node`),
  join(__dirname, 'lxephemeris.node'),
  join(__dirname, '..', '..', 'target', 'release', 'liblxephemeris_node.dylib'),
  join(__dirname, '..', '..', 'target', 'release', 'liblxephemeris_node.so'),
  join(__dirname, '..', '..', 'target', 'release', 'lxephemeris_node.dll'),
];

function loadNativeAddon(filePath) {
  if (filePath.endsWith('.node')) {
    return require(filePath);
  }
  const mod = { exports: {} };
  process.dlopen(mod, filePath);
  return mod.exports;
}

let nativeBinding = null;
for (const candidate of candidates) {
  if (existsSync(candidate)) {
    nativeBinding = loadNativeAddon(candidate);
    break;
  }
}

if (!nativeBinding) {
  throw new Error(
    'Failed to load native binding for lxephemeris-node. Looked in: ' + candidates.join(', ')
  );
}

function makeUtcMs(year, month, day, hour, minute, second) {
  const d = new Date(0);
  d.setUTCFullYear(year, month - 1, day);
  d.setUTCHours(hour, minute, Math.floor(second), Math.round((second % 1) * 1000));
  return d.getTime();
}

function getIanaOffsetMs(utcMs, timeZone) {
  const fmt = new Intl.DateTimeFormat('en-US', {
    timeZone,
    year: 'numeric',
    month: '2-digit',
    day: '2-digit',
    hour: '2-digit',
    minute: '2-digit',
    second: '2-digit',
    hourCycle: 'h23',
  });
  const parts = fmt.formatToParts(new Date(utcMs));
  const map = {};
  for (const p of parts) {
    if (p.type !== 'literal') {
      map[p.type] = parseInt(p.value, 10);
    }
  }
  const zonedMs = makeUtcMs(
    map.year,
    map.month,
    map.day,
    map.hour === 24 ? 0 : map.hour,
    map.minute,
    map.second
  );
  const roundedUtcMs = Math.floor(utcMs / 1000) * 1000;
  return zonedMs - roundedUtcMs;
}

function resolveTimezoneOffsetHours(year, month, day, hour, minute, second, tzOffsetOrIana) {
  if (tzOffsetOrIana === undefined || tzOffsetOrIana === null) {
    return 0.0;
  }
  if (typeof tzOffsetOrIana === 'number') {
    return tzOffsetOrIana;
  }
  if (typeof tzOffsetOrIana === 'string') {
    const trimmed = tzOffsetOrIana.trim();
    if (/^[+-]?\d+(\.\d+)?$/.test(trimmed)) {
      return parseFloat(trimmed);
    }
    const localWallMs = makeUtcMs(year, month, day, hour || 0, minute || 0, second || 0);
    const offset1 = getIanaOffsetMs(localWallMs, trimmed);
    const offset2 = getIanaOffsetMs(localWallMs - offset1, trimmed);
    return offset2 / 3600000.0;
  }
  throw new TypeError(`Invalid timezone argument: ${String(tzOffsetOrIana)}`);
}

const nativeUtcFromLocal = nativeBinding.utcFromLocal;

function utcFromLocal(
  year,
  month,
  day,
  hour = 0,
  minute = 0,
  second = 0.0,
  tzOffsetOrIana = 0.0,
  gregorian = true
) {
  const offsetHours = resolveTimezoneOffsetHours(
    year,
    month,
    day,
    hour ?? 0,
    minute ?? 0,
    second ?? 0.0,
    tzOffsetOrIana
  );
  return nativeUtcFromLocal(
    year,
    month,
    day,
    hour ?? 0,
    minute ?? 0,
    second ?? 0.0,
    offsetHours,
    gregorian ?? true
  );
}

function toScalarJd(input) {
  if (typeof input === 'number') {
    return input;
  }
  if (input instanceof Date) {
    const ms = input.getTime();
    if (Number.isNaN(ms)) {
      throw new TypeError('Invalid Date passed to lxephemeris');
    }
    return ms / 86400000.0 + 2440587.5;
  }
  if (typeof input === 'string') {
    const trimmed = input.trim();
    const d = new Date(trimmed);
    const ms = d.getTime();
    if (Number.isNaN(ms)) {
      throw new TypeError(`Invalid ISO date string passed to lxephemeris: ${input}`);
    }
    return ms / 86400000.0 + 2440587.5;
  }
  if (input && typeof input === 'object' && typeof input.jd1 === 'number') {
    return input.jd1 + (typeof input.jd2 === 'number' ? input.jd2 : 0.0);
  }
  throw new TypeError(`Unsupported time input: ${String(input)}`);
}

function toSplitJd(input) {
  if (input && typeof input === 'object' && !(input instanceof Date) && typeof input.jd1 === 'number') {
    return input;
  }
  return { jd1: toScalarJd(input), jd2: 0.0 };
}

function jdFromDate(dateOrIso) {
  return toScalarJd(dateOrIso);
}

function dateFromJd(jdInput) {
  const jd = toScalarJd(jdInput);
  return new Date((jd - 2440587.5) * 86400000.0);
}

const proto = nativeBinding.JsEngine.prototype;

// Wrap methods that expect JsSplitJd at a specific argument index so they accept number | Date | string | JsSplitJd
function wrapSplitJdArg(methodName, argIndex) {
  const orig = proto[methodName];
  if (typeof orig !== 'function') return;
  proto[methodName] = function (...args) {
    if (args.length > argIndex && args[argIndex] !== undefined && args[argIndex] !== null) {
      args[argIndex] = toSplitJd(args[argIndex]);
    }
    return orig.apply(this, args);
  };
}

// Wrap methods that expect Array<JsSplitJd> at a specific argument index
function wrapSplitJdArrayArg(methodName, argIndex) {
  const orig = proto[methodName];
  if (typeof orig !== 'function') return;
  proto[methodName] = function (...args) {
    if (Array.isArray(args[argIndex])) {
      args[argIndex] = args[argIndex].map(toSplitJd);
    }
    return orig.apply(this, args);
  };
}

// Wrap methods that expect two JsSplitJd arguments (e.g., startJd, endJd)
function wrapTwoSplitJdArgs(methodName, idx1, idx2) {
  const orig = proto[methodName];
  if (typeof orig !== 'function') return;
  proto[methodName] = function (...args) {
    if (args.length > idx1 && args[idx1] !== undefined && args[idx1] !== null) {
      args[idx1] = toSplitJd(args[idx1]);
    }
    if (args.length > idx2 && args[idx2] !== undefined && args[idx2] !== null) {
      args[idx2] = toSplitJd(args[idx2]);
    }
    return orig.apply(this, args);
  };
}

// Wrap methods that expect scalar JD (number) at specific argument indices so they also accept Date | string | JsSplitJd
function wrapScalarJdArgs(methodName, indices) {
  const orig = proto[methodName];
  if (typeof orig !== 'function') return;
  proto[methodName] = function (...args) {
    for (const idx of indices) {
      if (args.length > idx && args[idx] !== undefined && args[idx] !== null && typeof args[idx] !== 'number') {
        args[idx] = toScalarJd(args[idx]);
      }
    }
    return orig.apply(this, args);
  };
}

// SplitJd methods at index 0
for (const m of ['deltaT', 'sidtime', 'eclNut', 'houses', 'nextSolarEclipse', 'prevSolarEclipse', 'nextLocalSolarEclipse', 'nextLunarEclipse', 'prevLunarEclipse', 'timeEquation']) {
  wrapSplitJdArg(m, 0);
}
// SplitJd methods at index 1
for (const m of ['calc', 'calcStarByName', 'ayanamsa', 'calcAzalt', 'heliacalPhenomena', 'calcMoshier', 'phenomena', 'nodAps', 'nextLunarOccultation', 'azaltRev']) {
  wrapSplitJdArg(m, 1);
}
// SplitJd methods at index 2
for (const m of ['azalt', 'azaltEcliptic', 'riseTrans', 'nextLocalLunarOccultation']) {
  wrapSplitJdArg(m, 2);
}
// Array<JsSplitJd> methods
wrapSplitJdArrayArg('calcBatch', 1);
wrapSplitJdArrayArg('calcBatchFlatF64', 1);
wrapSplitJdArrayArg('calcBatchPerItem', 1);
wrapSplitJdArrayArg('housesBatch', 0);
// Two SplitJd methods
wrapTwoSplitJdArgs('lunarOccultationsInRange', 1, 2);
wrapTwoSplitJdArgs('localLunarOccultationsInRange', 2, 3);

// Scalar JD methods (accept number | Date | string | JsSplitJd)
for (const m of ['humanDesignEpoch', 'lunarPhaseSummary', 'findParans', 'calcChartAspects', 'calcDeclinationAspects', 'panchanga', 'vimshottariDasha', 'chartSect', 'planetaryHours', 'calcLots', 'chartSnapshot']) {
  wrapScalarJdArgs(m, [0]);
}
for (const m of ['solarEclipsesInRange', 'lunarEclipsesInRange', 'solarTermsInRange', 'lunarPhasesInRange', 'mutualAspectsInRange', 'voidOfCourseMoonInRange']) {
  wrapScalarJdArgs(m, [0, 1]);
}
wrapScalarJdArgs('searchBodyLongitude', [1]);
wrapScalarJdArgs('bodyLongitudeCrossingsInRange', [1, 2]);
wrapScalarJdArgs('searchLunarPhase', [0]);
wrapScalarJdArgs('searchStationsInRange', [1, 2]);
wrapScalarJdArgs('searchDeclinationExtremaInRange', [1, 2]);
wrapScalarJdArgs('planetaryReturn', [1, 2]);
wrapScalarJdArgs('planetaryReturnsInRange', [1, 2, 3]);
wrapScalarJdArgs('calcSynastryAspects', [0, 2]);
wrapScalarJdArgs('searchExactAspect', [3]);
wrapScalarJdArgs('natalTransitsInRange', [0, 1, 2]);

nativeBinding.Engine = nativeBinding.JsEngine;
nativeBinding.SweAdapter = nativeBinding.JsSweAdapter;
nativeBinding.utcFromLocal = utcFromLocal;
nativeBinding.toScalarJd = toScalarJd;
nativeBinding.toSplitJd = toSplitJd;
nativeBinding.jdFromDate = jdFromDate;
nativeBinding.dateFromJd = dateFromJd;

module.exports = nativeBinding;
