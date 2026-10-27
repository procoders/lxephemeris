/* TypeScript definitions for lxephemeris-node (napi-rs) */

export declare const SE_ECL_NUT: number;
export declare const SE_SUN: number;
export declare const SE_MOON: number;
export declare const SE_MERCURY: number;
export declare const SE_VENUS: number;
export declare const SE_MARS: number;
export declare const SE_JUPITER: number;
export declare const SE_SATURN: number;
export declare const SE_URANUS: number;
export declare const SE_NEPTUNE: number;
export declare const SE_PLUTO: number;
export declare const SE_MEAN_NODE: number;
export declare const SE_TRUE_NODE: number;
export declare const SE_MEAN_APOG: number;
export declare const SE_OSCU_APOG: number;
export declare const SE_EARTH: number;
export declare const SE_CHIRON: number;
export declare const SE_PHOLUS: number;
export declare const SE_CERES: number;
export declare const SE_PALLAS: number;
export declare const SE_JUNO: number;
export declare const SE_VESTA: number;
export declare const SE_INTP_APOG: number;
export declare const SE_FICT_OFFSET: number;
export declare const SE_CUPIDO: number;
export declare const SE_HADES: number;
export declare const SE_ZEUS: number;
export declare const SE_KRONOS: number;
export declare const SE_APOLLON: number;
export declare const SE_ADMETOS: number;
export declare const SE_VULKANUS: number;
export declare const SE_POSEIDON: number;
export declare const SE_ISIS: number;
export declare const SE_AST_OFFSET: number;

export declare const SEFLG_SWIEPH: number;
export declare const SEFLG_MOSEPH: number;
export declare const SEFLG_HELCTR: number;
export declare const SEFLG_TRUEPOS: number;
export declare const SEFLG_J2000: number;
export declare const SEFLG_NONUT: number;
export declare const SEFLG_SPEED: number;
export declare const SEFLG_NOGDEFL: number;
export declare const SEFLG_NOABERR: number;
export declare const SEFLG_BARYCTR: number;
export declare const SEFLG_TOPOCTR: number;
export declare const SEFLG_EQUATORIAL: number;
export declare const SEFLG_XYZ: number;
export declare const SEFLG_RADIANS: number;
export declare const SEFLG_SIDEREAL: number;

export declare const SE_SIDM_FAGAN_BRADLEY: number;
export declare const SE_SIDM_LAHIRI: number;
export declare const SE_SIDM_RAMAN: number;
export declare const SE_SIDM_USHASHASHI: number;
export declare const SE_SIDM_KRISHNAMURTI: number;
export declare const SE_SIDM_YUKTESHWAR: number;
export declare const SE_SIDM_JN_BHASIN: number;
export declare const SE_SIDM_SURYASIDDHANTA: number;
export declare const SE_SIDM_SS_CITRA: number;
export declare const SE_SIDM_TRUE_CITRA: number;
export declare const SE_SIDM_TRUE_REVATI: number;
export declare const SE_SIDM_TRUE_PUSHYA: number;
export declare const SE_SIDM_USER: number;
export declare const SE_SIDM_DELUCE: number;
export declare const SE_SIDM_ARYABHATA: number;
export declare const SE_SIDM_ARYABHATA_MSUN: number;
export declare const SE_SIDM_SS_REVATI: number;
export declare const SE_SIDM_SURYASIDDHANTA_MSUN: number;
export declare const SE_SIDM_HIPPARCHOS: number;
export declare const SE_SIDM_SASSANIAN: number;
export declare const SE_SIDM_J1900: number;
export declare const SE_SIDM_B1950: number;
export declare const SE_SIDM_J2000: number;
export declare const SE_SIDM_BABYL_KUGLER1: number;
export declare const SE_SIDM_BABYL_KUGLER2: number;
export declare const SE_SIDM_BABYL_KUGLER3: number;
export declare const SE_SIDM_BABYL_HUBER: number;
export declare const SE_SIDM_BABYL_ETPSC: number;
export declare const SE_SIDM_ALDEBARAN_15TAU: number;
export declare const SE_SIDM_GALCENT_RGILBRAND: number;
export declare const SE_SIDM_GALCENT_0SAG: number;
export declare const SE_SIDM_GALALIGN_MARDYKS: number;
export declare const SE_SIDM_TRUE_MULA: number;
export declare const SE_SIDM_GALCENT_MULA_WILHELM: number;
export declare const SE_SIDM_ARYABHATA_522: number;
export declare const SE_SIDM_BABYL_BRITTON: number;
export declare const SE_SIDM_TRUE_SHERATAN: number;
export declare const SE_SIDM_TRUE_ALDEBARAN: number;
export declare const SE_SIDM_DJWHAL_KHUL: number;
export declare const SE_SIDM_VALENS_MOON: number;
export declare const SE_SIDM_GALEQ_IAU1958: number;
export declare const SE_SIDM_GALEQ_TRUE: number;
export declare const SE_SIDM_GALEQU_FIORENZA: number;
export declare const SE_SIDM_GALEQ_MULA: number;

export declare const SE_CALC_RISE: number;
export declare const SE_CALC_SET: number;
export declare const SE_CALC_MTRANSIT: number;
export declare const SE_CALC_ITRANSIT: number;
export declare const SE_BIT_DISC_CENTER: number;
export declare const SE_BIT_DISC_BOTTOM: number;
export declare const SE_BIT_NO_REFRACTION: number;
export declare const SE_ECL2HOR: number;
export declare const SE_EQU2HOR: number;
export declare const SE_HOR2ECL: number;
export declare const SE_HOR2EQU: number;
export declare const SE_NODBIT_MEAN: number;
export declare const SE_NODBIT_OSCU: number;
export declare const SE_NODBIT_OSCU_BAR: number;
export declare const SE_NODBIT_FOPOINT: number;
export declare const SE_JUL_CAL: number;
export declare const SE_GREG_CAL: number;

export interface JsSplitJd {
  jd1: number;
  jd2?: number;
}

export interface JsObserver {
  latitudeDeg: number;
  longitudeDeg: number;
  altitudeM?: number;
}

export interface JsKernelBytes {
  name: string;
  data: Buffer;
}

export interface JsCalcOptions {
  origin?: string;
  plane?: string;
  equinox?: string;
  nutation?: string;
  aberration?: string;
  zodiac?: string;
  ayanamsa?: string;
  includeSpeed?: boolean;
  observer?: JsObserver;
}

export interface JsCalcResult {
  targetId: number;
  centerId: number;
  jdTt: number;
  jdUt1: number;
  jdTdb: number;
  deltaTSec: number;
  timeExtrapolated: boolean;
  timeTableVersion: string;
  lonDeg: number;
  latDeg: number;
  distAu: number;
  distKm: number;
  lonSpeedDegDay: number;
  latSpeedDegDay: number;
  distSpeedAuDay: number;
  distSpeedKmDay: number;
  magnitudeV?: number;
  cartesianAu: Array<number>;
  cartesianVelAuDay: Array<number>;
  warnings: Array<string>;
}

export interface JsBatchItemResult {
  ok: boolean;
  result?: JsCalcResult;
  errorCode?: string;
  errorMessage?: string;
}

export interface JsDeltaTInfo {
  deltaTSec: number;
  ttMinusTdbSec: number;
  extrapolated: boolean;
  tableVersion: string;
  warning?: string;
}

export interface JsSiderealTimeInfo {
  gmstHours: number;
  gastHours: number;
  localSiderealHours: number;
  armcDeg: number;
}

export interface JsEclNutInfo {
  trueObliquityDeg: number;
  meanObliquityDeg: number;
  nutationLonDeg: number;
  nutationOblDeg: number;
}

export interface JsStarMetadata {
  hip: number;
  canonicalName: string;
  bayer: string;
  vmag: number;
  raRad: number;
  decRad: number;
  pmRaMasYr: number;
  pmDecMasYr: number;
  plxMas: number;
  rvKms: number;
  rvMissing: boolean;
}

export interface JsHeliacalResult {
  targetName: string;
  visualMagnitude: number;
  arcusVisionisDeg: number;
  heliacalRiseJdUt1?: number;
  morningSetJdUt1?: number;
  eveningRiseJdUt1?: number;
  heliacalSetJdUt1?: number;
}

export interface JsAyanamsaValue {
  modelCode: string;
  jdTt: number;
  meanDeg: number;
  trueDeg: number;
  effectiveDeg: number;
  rateDegDay: number;
  dpsiArcsec: number;
}

export interface JsHouseResult {
  systemCode: string;
  effectiveSystemCode: string;
  polarHouseFallback: boolean;
  cuspsDeg: Array<number>;
  cuspSpeedsDegDay: Array<number>;
  ascDeg: number;
  mcDeg: number;
  armcDeg: number;
  vertexDeg: number;
  equatorialAscDeg: number;
  coAscKochDeg: number;
  coAscMunkaseyDeg: number;
  polarAscDeg: number;
  ascSpeedDegDay: number;
  mcSpeedDegDay: number;
  warnings: Array<string>;
}

export interface JsHorizontalCoords {
  azimuthDeg: number;
  altitudeTrueDeg: number;
  altitudeApparentDeg: number;
  hourAngleDeg: number;
  nearHorizon: boolean;
}

export interface JsSolarEclipseInfo {
  kind: string;
  sarosSeries: number;
  sarosMember: number;
  maxEclipseJdTt: number;
  maxEclipseJdUt1: number;
  deltaTSec: number;
  gamma: number;
  uRadiusEarth: number;
  penumbraRadiusEarth: number;
  magnitude: number;
  centralDurationSec?: number;
  pathWidthKm?: number;
  centralLatitudeDeg?: number;
  centralLongitudeDeg?: number;
  p1JdTt?: number;
  u1JdTt?: number;
  u2JdTt?: number;
  u3JdTt?: number;
  u4JdTt?: number;
  p4JdTt?: number;
  localVisibleAtMax?: boolean;
  localIsCentral?: boolean;
  localKind?: string;
  localC1JdTt?: number;
  localC2JdTt?: number;
  localMaxJdTt?: number;
  localC3JdTt?: number;
  localC4JdTt?: number;
  localMagnitude?: number;
  localObscuration?: number;
  localCentralDurationSec?: number;
  localSunAltitudeDeg?: number;
  localSunAzimuthDeg?: number;
}

export interface JsLunarEclipseInfo {
  kind: string;
  sarosSeries: number;
  sarosMember: number;
  shadowRule: string;
  maxEclipseJdTt: number;
  maxEclipseJdUt1: number;
  deltaTSec: number;
  gamma: number;
  penumbralMagnitude: number;
  umbralMagnitude: number;
  penumbralSemiDurationMin: number;
  partialSemiDurationMin: number;
  totalSemiDurationMin: number;
  p1JdTt?: number;
  u1JdTt?: number;
  u2JdTt?: number;
  u3JdTt?: number;
  u4JdTt?: number;
  p4JdTt?: number;
}

export interface JsCalendarDate {
  year: number;
  month: number;
  day: number;
  hourUt: number;
}

export interface JsSweCalcOut {
  xx: Array<number>;
  values: Array<number>;
  retflag: number;
}

export interface JsSweFixstarOut {
  xx: Array<number>;
  values: Array<number>;
  starName: string;
  retflag: number;
}

export interface JsSweAyanamsaOut {
  retflag: number;
  ayanamsa: number;
}

export interface JsSweHousesOut {
  cusps: Array<number>;
  ascmc: Array<number>;
}

export interface JsPhenomenaResult {
  phaseAngleDeg: number;
  phaseIlluminated: number;
  elongationDeg: number;
  apparentDiameterDeg: number;
  apparentMagnitudeV: number;
  horizontalParallaxDeg: number;
  attr20: Array<number>;
}

export interface JsNodApsResult {
  ascendingNode: Array<number>;
  descendingNode: Array<number>;
  periapsis: Array<number>;
  apoapsis: Array<number>;
  retflag: number;
}

export interface JsLunarOccultationInfo {
  kind: string;
  target: string;
  jdUt1Max: number;
  jdTtMax: number;
  t1FirstContactJdUt1?: number;
  t2InnerImmersionJdUt1?: number;
  t3InnerEmersionJdUt1?: number;
  t4LastContactJdUt1?: number;
  gamma: number;
  magnitude: number;
  moonAltitudeDeg?: number;
  tret: Array<number>;
}

export interface JsSwePhenoOut {
  attr: Array<number>;
  retflag: number;
}

export interface JsSweOccultGlobOut {
  retflag: number;
  tret: Array<number>;
}

export interface JsSweOccultLocOut {
  retflag: number;
  tret: Array<number>;
  attr: Array<number>;
}

export interface JsLongitudeCrossingEvent {
  body: string;
  targetLonDeg: number;
  jdTt: number;
  jdUt1: number;
  deltaTSec: number;
  actualLonDeg: number;
  latDeg: number;
  distAu: number;
  distKm: number;
  lonSpeedDegDay: number;
  isRetrograde: boolean;
}

export interface JsSolarTermEvent {
  termIndex: number;
  targetLonDeg: number;
  pinyinName: string;
  isMajorTerm: boolean;
  jdTt: number;
  jdUt1: number;
  deltaTSec: number;
}

export interface JsSeasonsOfYear {
  year: number;
  vernalEquinoxJdUt1: number;
  vernalEquinoxJdTt: number;
  summerSolsticeJdUt1: number;
  summerSolsticeJdTt: number;
  autumnalEquinoxJdUt1: number;
  autumnalEquinoxJdTt: number;
  winterSolsticeJdUt1: number;
  winterSolsticeJdTt: number;
}

export interface JsLunarPhaseEvent {
  kind?: string;
  targetElongationDeg: number;
  jdTt: number;
  jdUt1: number;
  deltaTSec: number;
  moonLonDeg: number;
  sunLonDeg: number;
  moonDistAu: number;
  moonDistKm: number;
}

export interface JsLunarPhaseSummary {
  jdTt: number;
  jdUt1: number;
  elongationDeg: number;
  illuminationFraction: number;
  phaseAngleDeg: number;
  moonLonDeg: number;
  moonLatDeg: number;
  moonDistAu: number;
  moonDistKm: number;
  apparentDiameterDeg: number;
  magnitudeV: number;
  ageDays: number;
  synodicFraction: number;
  phase8Index: number;
  phase8Name: string;
  phase28Index: number;
  prevNewMoonJdUt1: number;
  nextNewMoonJdUt1: number;
  nextFirstQuarterJdUt1: number;
  nextFullMoonJdUt1: number;
  nextLastQuarterJdUt1: number;
}

export interface JsStationEvent {
  body: string;
  kind: string;
  jdTt: number;
  jdUt1: number;
  deltaTSec: number;
  lonDeg: number;
  latDeg: number;
  distAu: number;
  distKm: number;
}

export interface JsDeclinationExtremumEvent {
  body: string;
  kind: string;
  jdTt: number;
  jdUt1: number;
  deltaTSec: number;
  raDeg: number;
  decDeg: number;
  distAu: number;
  distKm: number;
}

export interface JsParanEvent {
  starName: string;
  starBayer: string;
  starHip: number;
  planetName: string;
  paranType: string;
  starAngle: string;
  planetAngle: string;
  sameAngle: boolean;
  starEventJdUt1: number;
  planetEventJdUt1: number;
  timeOrbMinutes: number;
}

export interface JsHouseOptions {
  system?: HouseSystemCode | undefined | null;
  polarPolicy?: 'strict_error' | 'fallback_porphyry' | 'fallback_whole_sign' | (string & {}) | undefined | null;
  zodiac?: 'tropical' | 'sidereal' | (string & {}) | undefined | null;
  ayanamsa?: AyanamsaName | undefined | null;
  nutation?: 'true' | 'mean' | (string & {}) | undefined | null;
  includeSpeed?: boolean | undefined | null;
}

/**
 * Flexible time input accepted by all `Engine` methods:
 * - scalar Julian Day number (`2461338.0833`)
 * - native JavaScript `Date` (`new Date("2026-10-24T14:00:00Z")`)
 * - ISO-8601 string (`"2026-10-24T14:00:00+02:00"`)
 * - high-precision split Julian Day (`{ jd1: 2461338.0, jd2: 0.0833 }`)
 */
export type TimeLike = number | Date | string | JsSplitJd;

export type BodyName =
  | 'Sun'
  | 'Moon'
  | 'Mercury'
  | 'Venus'
  | 'Mars'
  | 'Jupiter'
  | 'Saturn'
  | 'Uranus'
  | 'Neptune'
  | 'Pluto'
  | 'Earth'
  | 'Mean_Node'
  | 'True_Node'
  | 'South_Mean_Node'
  | 'South_True_Node'
  | 'Mean_Lilith'
  | 'True_Lilith'
  | 'Interpolated_Lilith'
  | 'Chiron'
  | 'Pholus'
  | 'Ceres'
  | 'Pallas'
  | 'Juno'
  | 'Vesta'
  | 'Eris'
  | 'Sedna'
  | 'Cupido'
  | 'Hades'
  | 'Zeus'
  | 'Kronos'
  | 'Apollon'
  | 'Admetos'
  | 'Vulkanus'
  | 'Poseidon'
  | 'Proserpina'
  | (string & {});

export type HouseSystemCode =
  | 'P' | 'K' | 'O' | 'R' | 'C' | 'A' | 'E' | 'V' | 'W' | 'X' | 'M' | 'B'
  | 'F' | 'G' | 'H' | 'T' | 'U' | 'Y' | 'D' | 'I' | 'L' | 'N' | 'S'
  | (string & {});

export type AyanamsaName =
  | 'lahiri'
  | 'fagan_bradley'
  | 'raman'
  | 'krishnamurti'
  | 'yukteshwar'
  | 'jn_bhasin'
  | 'true_citra'
  | 'true_revati'
  | 'true_pushya'
  | 'suryasiddhanta'
  | 'ushashashi'
  | 'j2000'
  | (string & {});

export type AspectPresetName = 'major' | 'standard' | 'kepler' | 'all' | (string & {});

export type TimeScale = 'UTC' | 'UT1' | 'TT' | 'TDB' | (string & {});

export interface JsZodiacPlacement {
  lonDeg: number;
  signIndex: number;
  signName: string;
  degInSign: number;
  dmsDeg: number;
  dmsMin: number;
  dmsSec: number;
  decan: number;
  nakshatraIndex: number;
  nakshatraName: string;
  nakshatraLord: string;
  pada: number;
  hdGate: number;
  hdLine: number;
  hdColor: number;
  hdTone: number;
  hdBase: number;
}

export interface JsAspectHit {
  body1: string;
  body2: string;
  aspectName: string;
  exactAngleDeg: number;
  actualSeparationDeg: number;
  orbDeg: number;
  maxOrbDeg: number;
  tightnessPct: number;
  isApplying: boolean;
  isMajor: boolean;
  relativeSpeedDegDay: number;
}

export interface JsDeclinationAspectHit {
  body1: string;
  body2: string;
  kind: 'Parallel' | 'Contraparallel' | string;
  dec1Deg: number;
  dec2Deg: number;
  orbDeg: number;
  maxOrbDeg: number;
  tightnessPct: number;
  isApplying: boolean;
  relativeDecSpeedDegDay: number;
}

export interface JsExactAspectEvent {
  body1: string;
  body2: string;
  aspectName: string;
  targetAngleDeg: number;
  jdUt1: number;
  jdTt: number;
  body1LonDeg: number;
  body2LonDeg: number;
  relativeSpeedDegDay: number;
  body1Retrograde: boolean;
  body2Retrograde: boolean;
}

export interface JsPlanetaryReturnEvent {
  body: string;
  natalLonDeg: number;
  jdUt1: number;
  jdTt: number;
  actualLonDeg: number;
  lonSpeedDegDay: number;
  isRetrograde: boolean;
  passIndex: number;
}

export interface JsVoidOfCourseInterval {
  startJdUt1: number;
  endJdUt1: number;
  durationHours: number;
  lastAspectBody: string;
  lastAspectName: string;
  lastAspectAngleDeg: number;
  moonSignBefore: number;
  moonSignAfter: number;
  isBonattiMildSign: boolean;
}

export interface JsTimeWindow {
  startJdUt1: number;
  endJdUt1: number;
}

export interface JsPanchangaReport {
  jdUt1: number;
  sunriseJdUt1: number;
  sunsetJdUt1: number;
  nextSunriseJdUt1: number;
  polarDayFallback: boolean;
  ayanamsaDeg: number;
  vara: {
    varaIndex: number;
    name: string;
    lord: string;
  };
  tithi: {
    tithiIndex: number;
    name: string;
    paksha: 'Shukla' | 'Krishna' | string;
    elapsedFraction: number;
    startJdUt1: number;
    endJdUt1: number;
  };
  nakshatra: {
    nakshatraIndex: number;
    name: string;
    vimshottariLord: string;
    pada: number;
    elapsedFraction: number;
    startJdUt1: number;
    endJdUt1: number;
  };
  yoga: {
    yogaIndex: number;
    name: string;
    elapsedFraction: number;
    startJdUt1: number;
    endJdUt1: number;
  };
  karana: {
    karanaIndex: number;
    name: string;
    isVishtiBhadra: boolean;
    elapsedFraction: number;
    startJdUt1: number;
    endJdUt1: number;
  };
  rahukalam: JsTimeWindow;
  yamagandam: JsTimeWindow;
  gulikakalam: JsTimeWindow;
  abhijitMuhurta: JsTimeWindow;
}

export interface JsDashaSubPeriod {
  mahaLord: string;
  antarLord: string;
  startJdUt1: number;
  endJdUt1: number;
  durationDays: number;
}

export interface JsDashaPeriod {
  lord: string;
  startJdUt1: number;
  endJdUt1: number;
  durationYears: number;
  subPeriods: Array<JsDashaSubPeriod>;
}

export interface JsVimshottariDashaReport {
  birthJdUt1: number;
  moonSiderealLonDeg: number;
  birthNakshatraIndex: number;
  birthNakshatraName: string;
  birthPada: number;
  startingLord: string;
  balanceOfFirstDashaYears: number;
  yearLengthDays: number;
  periods: Array<JsDashaPeriod>;
}

export interface JsVargaPlacement {
  division: number;
  sanskritName: string;
  vargaLonDeg: number;
  signIndex: number;
  signName: string;
  signLord: string;
  degInSign: number;
}

export interface JsChartSectInfo {
  isDayChart: boolean;
  sunAltitudeTrueDeg: number;
  sunAltitudeApparentDeg: number;
  nearHorizon: boolean;
}

export interface JsPlanetaryHourSlot {
  hourIndex: number;
  isDaytime: boolean;
  lord: string;
  startJdUt1: number;
  endJdUt1: number;
}

export interface JsPlanetaryHoursReport {
  jdUt1: number;
  sunriseJdUt1: number;
  sunsetJdUt1: number;
  nextSunriseJdUt1: number;
  dayLord: string;
  currentHourIndex: number;
  currentHourLord: string;
  isDaytimeHour: boolean;
  hours: Array<JsPlanetaryHourSlot>;
}

export interface JsHermeticLot {
  name: string;
  lonDeg: number;
  placement: JsZodiacPlacement;
  houseNumber: number;
}

export interface JsChartBodySnapshot {
  body: string;
  bodyName: string;
  calc: JsCalcResult;
  placement: JsZodiacPlacement;
  houseNumber: number;
  housePos3d: number;
  isRetrograde: boolean;
}

export interface JsChartSnapshot {
  jdUt1: number;
  jdTt: number;
  observer: JsObserver;
  bodies: Array<JsChartBodySnapshot>;
  houses: JsHouseResult;
  aspects: Array<JsAspectHit>;
  lunarPhase: JsLunarPhaseSummary;
  sect: JsChartSectInfo;
  lots: Array<JsHermeticLot>;
}

export interface JsLicenseInfo {
  status: 'missing' | 'valid' | 'expired' | 'invalid_signature' | 'malformed' | 'unsupported_version';
  tier: string;
  isCommercialValid: boolean;
  licenseId: string | null;
  holder: string | null;
  keyId: string | null;
  issuedAt: string | null;
  updatesUntil: string | null;
  engineReleaseDate: string;
  message: string;
}

export declare class JsEngine {
  constructor(
    dataDir?: string | undefined | null,
    maxThreads?: number | undefined | null,
    pinTimeTableVersion?: string | undefined | null,
    verifyManifestSha256?: boolean | undefined | null
  );
  static fromDataDir(
    dataDir: string,
    maxThreads?: number | undefined | null,
    pinTimeTableVersion?: string | undefined | null,
    verifyManifestSha256?: boolean | undefined | null
  ): JsEngine;
  static fromBspPaths(
    paths: Array<string>,
    maxThreads?: number | undefined | null,
    pinTimeTableVersion?: string | undefined | null
  ): JsEngine;
  static fromBytes(
    kernels: Array<JsKernelBytes>,
    maxThreads?: number | undefined | null
  ): JsEngine;
  static newMoshier(
    maxThreads?: number | undefined | null,
    pinTimeTableVersion?: string | undefined | null
  ): JsEngine;
  static fromDataDirOrMoshier(
    dataDir: string,
    maxThreads?: number | undefined | null,
    pinTimeTableVersion?: string | undefined | null,
    verifyManifestSha256?: boolean | undefined | null
  ): JsEngine;
  capabilitiesJson(): string;
  licenseInfo(): JsLicenseInfo;
  coverageJson(body: BodyName): string;
  deltaT(jd: TimeLike, scale?: TimeScale | undefined | null): JsDeltaTInfo;
  sidtime(
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    nutation?: string | undefined | null,
    longitudeDeg?: number | undefined | null
  ): JsSiderealTimeInfo;
  eclNut(jd: TimeLike, scale?: TimeScale | undefined | null): JsEclNutInfo;
  starInfo(nameOrBayer: string): JsStarMetadata;
  listStars(): Array<JsStarMetadata>;
  calc(
    body: BodyName,
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): JsCalcResult;
  calcStarByName(
    nameOrBayer: string,
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): JsCalcResult;
  calcBatch(
    bodies: Array<BodyName>,
    epochs: Array<TimeLike>,
    scale?: TimeScale | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): Array<Array<JsCalcResult>>;
  calcBatchFlatF64(
    bodies: Array<BodyName>,
    epochs: Array<TimeLike>,
    scale?: TimeScale | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): Float64Array;
  calcBatchPerItem(
    bodies: Array<BodyName>,
    epochs: Array<TimeLike>,
    scale?: TimeScale | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): Array<Array<JsBatchItemResult>>;
  ayanamsa(
    model: AyanamsaName,
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    nutation?: string | undefined | null
  ): JsAyanamsaValue;
  houses(
    jd: TimeLike,
    observer: JsObserver,
    system?: HouseSystemCode | undefined | null,
    scale?: TimeScale | undefined | null,
    polarPolicy?: string | undefined | null,
    zodiac?: string | undefined | null,
    ayanamsa?: AyanamsaName | undefined | null,
    nutation?: string | undefined | null,
    includeSpeed?: boolean | undefined | null
  ): JsHouseResult;
  housesBatch(
    epochs: Array<TimeLike>,
    observer: JsObserver,
    system?: HouseSystemCode | undefined | null,
    scale?: TimeScale | undefined | null,
    polarPolicy?: string | undefined | null,
    zodiac?: string | undefined | null,
    ayanamsa?: AyanamsaName | undefined | null,
    nutation?: string | undefined | null,
    includeSpeed?: boolean | undefined | null
  ): Array<JsHouseResult>;
  azalt(
    raDeg: number,
    decDeg: number,
    jd: TimeLike,
    observer: JsObserver,
    scale?: TimeScale | undefined | null,
    pressureHpa?: number | undefined | null,
    temperatureC?: number | undefined | null
  ): JsHorizontalCoords;
  azaltEcliptic(
    eclLonDeg: number,
    eclLatDeg: number,
    jd: TimeLike,
    observer: JsObserver,
    scale?: TimeScale | undefined | null,
    pressureHpa?: number | undefined | null,
    temperatureC?: number | undefined | null
  ): JsHorizontalCoords;
  calcAzalt(
    body: BodyName,
    jd: TimeLike,
    observer: JsObserver,
    scale?: TimeScale | undefined | null,
    pressureHpa?: number | undefined | null,
    temperatureC?: number | undefined | null
  ): JsHorizontalCoords;
  riseTrans(
    body: BodyName,
    event: 'rise' | 'set' | 'upper_transit' | 'lower_transit' | (string & {}),
    jd: TimeLike,
    observer: JsObserver,
    scale?: TimeScale | undefined | null,
    customHorizonAltDeg?: number | undefined | null,
    maxSearchDays?: number | undefined | null,
    discMode?: string | undefined | null,
    pressureHpa?: number | undefined | null,
    temperatureC?: number | undefined | null
  ): number | null;
  heliacalPhenomena(
    body: BodyName,
    jd: TimeLike,
    observer: JsObserver,
    scale?: TimeScale | undefined | null,
    extinctionKV?: number | undefined | null,
    snellenFactor?: number | undefined | null,
    pressureHpa?: number | undefined | null,
    temperatureC?: number | undefined | null,
    relativeHumidity?: number | undefined | null,
    maxSearchDays?: number | undefined | null,
    observerAgeYears?: number | undefined | null,
    snellenRatio?: number | undefined | null
  ): JsHeliacalResult;
  nextSolarEclipse(
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    backward?: boolean | undefined | null,
    observer?: JsObserver | undefined | null,
    maxSynodicMonths?: number | undefined | null
  ): JsSolarEclipseInfo;
  prevSolarEclipse(
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    observer?: JsObserver | undefined | null,
    maxSynodicMonths?: number | undefined | null
  ): JsSolarEclipseInfo;
  nextLocalSolarEclipse(
    jd: TimeLike,
    observer: JsObserver,
    scale?: TimeScale | undefined | null,
    backward?: boolean | undefined | null,
    maxSynodicMonths?: number | undefined | null
  ): JsSolarEclipseInfo;
  solarEclipsesInRange(
    startJdTt: TimeLike,
    endJdTt: TimeLike,
    observer?: JsObserver | undefined | null
  ): Array<JsSolarEclipseInfo>;
  nextLunarEclipse(
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    backward?: boolean | undefined | null,
    shadowRule?: string | undefined | null,
    maxSynodicMonths?: number | undefined | null
  ): JsLunarEclipseInfo;
  prevLunarEclipse(
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    shadowRule?: string | undefined | null,
    maxSynodicMonths?: number | undefined | null
  ): JsLunarEclipseInfo;
  lunarEclipsesInRange(
    startJdTt: TimeLike,
    endJdTt: TimeLike,
    shadowRule?: string | undefined | null
  ): Array<JsLunarEclipseInfo>;
  gauquelinSectors(
    armcDeg: number,
    epsTrueDeg: number,
    latitudeDeg: number,
    polarPolicy?: string | undefined | null
  ): Array<number>;
  isMoshierMode(): boolean;
  calcMoshier(
    body: BodyName,
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): JsCalcResult;
  loadCustomStars(path: string): number;
  phenomena(
    body: BodyName,
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): JsPhenomenaResult;
  nodAps(
    body: BodyName,
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    options?: JsCalcOptions | undefined | null,
    methodBits?: number | undefined | null
  ): JsNodApsResult;
  nextLunarOccultation(
    target: string,
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    backward?: boolean | undefined | null,
    maxSearchDays?: number | undefined | null
  ): JsLunarOccultationInfo;
  nextLocalLunarOccultation(
    target: string,
    observer: JsObserver,
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    backward?: boolean | undefined | null,
    maxSearchDays?: number | undefined | null
  ): JsLunarOccultationInfo;
  lunarOccultationsInRange(
    target: string,
    startJd: TimeLike,
    endJd: TimeLike,
    scale?: TimeScale | undefined | null
  ): Array<JsLunarOccultationInfo>;
  localLunarOccultationsInRange(
    target: string,
    observer: JsObserver,
    startJd: TimeLike,
    endJd: TimeLike,
    scale?: TimeScale | undefined | null
  ): Array<JsLunarOccultationInfo>;
  housePos(
    armcDeg: number,
    geolatDeg: number,
    epsDeg: number,
    system: HouseSystemCode,
    eclLonDeg: number,
    eclLatDeg: number
  ): number;
  timeEquation(jd: TimeLike, scale?: TimeScale | undefined | null): number;
  azaltRev(
    observer: JsObserver,
    jd: TimeLike,
    scale?: TimeScale | undefined | null,
    azimuthNorthEastDeg: number,
    altitudeTrueDeg: number,
    targetPlane?: string | undefined | null
  ): Array<number>;
  searchBodyLongitude(
    body: BodyName,
    jdUt1: TimeLike,
    targetLonDeg: number,
    forward?: boolean | undefined | null,
    maxSearchDays?: number | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): JsLongitudeCrossingEvent;
  bodyLongitudeCrossingsInRange(
    body: BodyName,
    startJdUt1: TimeLike,
    endJdUt1: TimeLike,
    targetLonDeg: number,
    options?: JsCalcOptions | undefined | null
  ): Array<JsLongitudeCrossingEvent>;
  solarTermsInRange(
    startJdUt1: TimeLike,
    endJdUt1: TimeLike
  ): Array<JsSolarTermEvent>;
  seasonsInYear(year: number): JsSeasonsOfYear;
  humanDesignEpoch(birthJdUt1: TimeLike): JsLongitudeCrossingEvent;
  searchLunarPhase(
    startJdUt1: TimeLike,
    targetElongationDeg: number,
    forward?: boolean | undefined | null
  ): JsLunarPhaseEvent;
  lunarPhasesInRange(
    startJdUt1: TimeLike,
    endJdUt1: TimeLike
  ): Array<JsLunarPhaseEvent>;
  lunarPhaseSummary(jdUt1: TimeLike): JsLunarPhaseSummary;
  searchStationsInRange(
    body: BodyName,
    startJdUt1: TimeLike,
    endJdUt1: TimeLike,
    options?: JsCalcOptions | undefined | null
  ): Array<JsStationEvent>;
  searchDeclinationExtremaInRange(
    body: BodyName,
    startJdUt1: TimeLike,
    endJdUt1: TimeLike,
    options?: JsCalcOptions | undefined | null
  ): Array<JsDeclinationExtremumEvent>;
  findParans(
    jdUt1: TimeLike,
    observer: JsObserver,
    orbMinutes?: number | undefined | null,
    bodies?: Array<BodyName> | undefined | null,
    starNames?: Array<string> | undefined | null
  ): Array<JsParanEvent>;
  planetaryReturn(
    body: BodyName,
    natalJdUt1: TimeLike,
    searchAfterJdUt1: TimeLike,
    options?: JsCalcOptions | undefined | null
  ): JsPlanetaryReturnEvent;
  planetaryReturnsInRange(
    body: BodyName,
    natalJdUt1: TimeLike,
    startJdUt1: TimeLike,
    endJdUt1: TimeLike,
    options?: JsCalcOptions | undefined | null
  ): Array<JsPlanetaryReturnEvent>;
  calcChartAspects(
    jdUt1: TimeLike,
    bodies?: Array<BodyName> | undefined | null,
    preset?: AspectPresetName | undefined | null,
    orbScale?: number | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): Array<JsAspectHit>;
  calcSynastryAspects(
    jdUt1A: TimeLike,
    bodiesA: Array<BodyName>,
    jdUt1B: TimeLike,
    bodiesB: Array<BodyName>,
    preset?: AspectPresetName | undefined | null,
    orbScale?: number | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): Array<JsAspectHit>;
  calcDeclinationAspects(
    jdUt1: TimeLike,
    bodies?: Array<BodyName> | undefined | null,
    maxOrbDeg?: number | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): Array<JsDeclinationAspectHit>;
  searchExactAspect(
    body1: BodyName,
    body2: BodyName,
    targetAngleDeg: number,
    startJdUt1: TimeLike,
    forward?: boolean | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): JsExactAspectEvent;
  mutualAspectsInRange(
    startJdUt1: TimeLike,
    endJdUt1: TimeLike,
    bodies?: Array<BodyName> | undefined | null,
    preset?: AspectPresetName | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): Array<JsExactAspectEvent>;
  natalTransitsInRange(
    natalJdUt1: TimeLike,
    startJdUt1: TimeLike,
    endJdUt1: TimeLike,
    transitingBodies?: Array<BodyName> | undefined | null,
    natalBodies?: Array<BodyName> | undefined | null,
    preset?: AspectPresetName | undefined | null,
    options?: JsCalcOptions | undefined | null
  ): Array<JsExactAspectEvent>;
  voidOfCourseMoonInRange(
    startJdUt1: TimeLike,
    endJdUt1: TimeLike,
    options?: JsCalcOptions | undefined | null
  ): Array<JsVoidOfCourseInterval>;
  panchanga(
    jdUt1: TimeLike,
    observer: JsObserver,
    ayanamsa?: AyanamsaName | undefined | null
  ): JsPanchangaReport;
  vimshottariDasha(
    birthJdUt1: TimeLike,
    ayanamsa?: AyanamsaName | undefined | null,
    yearLengthDays?: number | undefined | null
  ): JsVimshottariDashaReport;
  calcVarga(siderealLonDeg: number, division: number): JsVargaPlacement;
  calcShodashavarga(siderealLonDeg: number): Array<JsVargaPlacement>;
  chartSect(jdUt1: TimeLike, observer: JsObserver): JsChartSectInfo;
  planetaryHours(jdUt1: TimeLike, observer: JsObserver): JsPlanetaryHoursReport;
  calcLots(
    jdUt1: TimeLike,
    observer: JsObserver,
    houseOptions?: JsHouseOptions | undefined | null
  ): Array<JsHermeticLot>;
  chartSnapshot(
    jdUt1: TimeLike,
    observer: JsObserver,
    bodies?: Array<BodyName> | undefined | null,
    houseOptions?: JsHouseOptions | undefined | null,
    calcOptions?: JsCalcOptions | undefined | null,
    preset?: AspectPresetName | undefined | null,
    orbScale?: number | undefined | null
  ): JsChartSnapshot;
  mapZodiacPlacement(lonDeg: number): JsZodiacPlacement;
}

export declare class JsSweAdapter {
  constructor(engine: JsEngine);
  static fromEngine(engine: JsEngine): JsSweAdapter;
  static newMoshier(): JsSweAdapter;
  static fromDataDir(dataDir: string): JsSweAdapter;
  static fromDataDirOrMoshier(dataDir: string): JsSweAdapter;
  setSidMode(
    sidMode: number,
    t0?: number | undefined | null,
    ayanT0?: number | undefined | null
  ): void;
  setTopo(lonDeg: number, latDeg: number, altM?: number | undefined | null): void;
  setPolarPolicy(policy: string): void;
  deltat(jdUt: number): number;
  deltatSec(jdUt: number): number;
  sidtime(jdUt: number): number;
  sidtime0(jdUt: number, epsDeg: number, nutLonDeg: number): number;
  getAyanamsaUt(jdUt: number): number;
  getAyanamsa(jdEt: number): number;
  getAyanamsaExUt(jdUt: number, iflag?: number | undefined | null): JsSweAyanamsaOut;
  calcUt(jdUt: number, ipl: number, iflag?: number | undefined | null): JsSweCalcOut;
  calc(jdEt: number, ipl: number, iflag?: number | undefined | null): JsSweCalcOut;
  fixstar2Ut(
    star: string,
    jdUt: number,
    iflag?: number | undefined | null
  ): JsSweFixstarOut;
  fixstar2Mag(star: string): number;
  housesEx(
    jdUt: number,
    latDeg: number,
    lonDeg: number,
    hsys?: string | undefined | null,
    iflag?: number | undefined | null
  ): JsSweHousesOut;
  houses(
    jdUt: number,
    latDeg: number,
    lonDeg: number,
    hsys?: string | undefined | null
  ): JsSweHousesOut;
  housesArmc(
    armcDeg: number,
    latDeg: number,
    epsDeg: number,
    hsys?: string | undefined | null
  ): JsSweHousesOut;
  azalt(
    jdUt: number,
    calcFlag: number,
    geopos: Array<number>,
    atpress: number,
    attemp: number,
    xin: Array<number>
  ): Array<number>;
  riseTrans(
    jdUt: number,
    ipl: number,
    rsmi: number,
    geopos: Array<number>,
    atpress?: number | undefined | null,
    attemp?: number | undefined | null,
    starname?: string | undefined | null,
    epheflag?: number | undefined | null
  ): number | null;
  phenoUt(
    jdUt: number,
    ipl: number,
    iflag?: number | undefined | null
  ): JsSwePhenoOut;
  pheno(
    jdEt: number,
    ipl: number,
    iflag?: number | undefined | null
  ): JsSwePhenoOut;
  nodApsUt(
    jdUt: number,
    ipl: number,
    iflag?: number | undefined | null,
    method?: number | undefined | null
  ): JsNodApsResult;
  nodAps(
    jdEt: number,
    ipl: number,
    iflag?: number | undefined | null,
    method?: number | undefined | null
  ): JsNodApsResult;
  lunOccultWhenGlob(
    jdUtStart: number,
    ipl: number,
    starname?: string | undefined | null,
    iflag?: number | undefined | null,
    backward?: boolean | undefined | null
  ): JsSweOccultGlobOut;
  lunOccultWhenLoc(
    jdUtStart: number,
    ipl: number,
    starname?: string | undefined | null,
    iflag?: number | undefined | null,
    geopos: Array<number>,
    backward?: boolean | undefined | null
  ): JsSweOccultLocOut;
  housePos(
    armcDeg: number,
    geolatDeg: number,
    epsDeg: number,
    xpin: Array<number>,
    hsys?: string | undefined | null
  ): number;
  timeEqu(jdEt: number): number;
  azaltRev(
    jdUt: number,
    calcFlag: number,
    geopos: Array<number>,
    xin: Array<number>
  ): Array<number>;
  static julday(
    year: number,
    month: number,
    day: number,
    hour?: number | undefined | null,
    gregflag?: number | undefined | null
  ): number;
  static revjul(jd: number, gregflag?: number | undefined | null): JsCalendarDate;
}

export declare function julday(
  year: number,
  month: number,
  day: number,
  hourUt?: number | undefined | null,
  calendar?: string | undefined | null
): number;

export declare function revjul(
  jd: number,
  calendar?: string | undefined | null
): JsCalendarDate;

export declare function housesArmc(
  armcDeg: number,
  epsTrueDeg: number,
  latitudeDeg: number,
  system?: string | undefined | null,
  polarPolicy?: string | undefined | null,
  ayanamsaDeg?: number | undefined | null,
  sunDeclDeg?: number | undefined | null,
  includeSpeed?: boolean | undefined | null
): JsHouseResult;

export declare function houseOfLon(
  lonDeg: number,
  cuspsDeg: Array<number>,
  epsDeg?: number | undefined | null
): number;

export declare function starInfo(nameOrBayer: string): JsStarMetadata;

export declare function listStars(): Array<JsStarMetadata>;

export declare function gauquelinSectors(
  armcDeg: number,
  epsTrueDeg: number,
  latitudeDeg: number,
  polarPolicy?: string | undefined | null
): Array<number>;

export declare function utcFromLocal(
  year: number,
  month: number,
  day: number,
  hour?: number | undefined | null,
  minute?: number | undefined | null,
  second?: number | undefined | null,
  tzOffsetHoursOrIana?: number | string | undefined | null,
  gregorian?: boolean | undefined | null
): number;

export declare function mapZodiacPlacement(lonDeg: number): JsZodiacPlacement;

export declare function verifyLicenseKey(token: string): JsLicenseInfo;
export declare function registerCommercialLicense(token: string): JsLicenseInfo;
export declare function clearCommercialLicense(): void;
export declare function licenseInfo(): JsLicenseInfo;

export declare function jdFromDate(dateOrIso: Date | string | number | JsSplitJd): number;
export declare function dateFromJd(jdInput: TimeLike): Date;
export declare function toScalarJd(input: TimeLike): number;
export declare function toSplitJd(input: TimeLike): JsSplitJd;

export { JsEngine as Engine, JsSweAdapter as SweAdapter };
