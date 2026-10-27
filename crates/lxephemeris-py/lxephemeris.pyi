"""PEP 484 Type Stubs for `lxephemeris` (`PyO3` extension module)."""

from typing import Literal, Optional, TypedDict, Union

# Typed Exception Hierarchy (R10.2)
class LxEphemerisError(Exception): ...
Error = LxEphemerisError
class InvalidInputError(LxEphemerisError, ValueError): ...
class UnsupportedBodyError(LxEphemerisError, ValueError): ...
class UnsupportedModelError(LxEphemerisError, ValueError): ...
class OutOfRangeError(LxEphemerisError, ValueError): ...
class PolarSingularityError(LxEphemerisError, ValueError): ...
class ConvergenceError(LxEphemerisError, RuntimeError): ...
class MissingDataError(LxEphemerisError, RuntimeError): ...
class CorruptKernelError(LxEphemerisError, RuntimeError): ...

version: str
__version__: str

# Swiss Ephemeris Migration Constants (SE_*/SEFLG_* and unprefixed aliases)
SE_ECL_NUT: int
ECL_NUT: int
SE_SUN: int
SUN: int
SE_MOON: int
MOON: int
SE_MERCURY: int
MERCURY: int
SE_VENUS: int
VENUS: int
SE_MARS: int
MARS: int
SE_JUPITER: int
JUPITER: int
SE_SATURN: int
SATURN: int
SE_URANUS: int
URANUS: int
SE_NEPTUNE: int
NEPTUNE: int
SE_PLUTO: int
PLUTO: int
SE_MEAN_NODE: int
MEAN_NODE: int
SE_TRUE_NODE: int
TRUE_NODE: int
SE_MEAN_APOG: int
MEAN_APOG: int
SE_OSCU_APOG: int
OSCU_APOG: int
SE_EARTH: int
EARTH: int
SE_CHIRON: int
CHIRON: int
SE_PHOLUS: int
PHOLUS: int
SE_CERES: int
CERES: int
SE_PALLAS: int
PALLAS: int
SE_JUNO: int
JUNO: int
SE_VESTA: int
VESTA: int
SE_INTP_APOG: int
INTP_APOG: int
SE_INTP_PERG: int
INTP_PERG: int
SE_FICT_OFFSET: int
FICT_OFFSET: int
SE_CUPIDO: int
CUPIDO: int
SE_HADES: int
HADES: int
SE_ZEUS: int
ZEUS: int
SE_KRONOS: int
KRONOS: int
SE_APOLLON: int
APOLLON: int
SE_ADMETOS: int
ADMETOS: int
SE_VULKANUS: int
VULKANUS: int
SE_POSEIDON: int
POSEIDON: int
SE_ISIS: int
ISIS: int
SE_AST_OFFSET: int
AST_OFFSET: int

SEFLG_JPLEPH: int
FLG_JPLEPH: int
SEFLG_SWIEPH: int
FLG_SWIEPH: int
SEFLG_MOSEPH: int
FLG_MOSEPH: int
SEFLG_HELCTR: int
FLG_HELCTR: int
SEFLG_TRUEPOS: int
FLG_TRUEPOS: int
SEFLG_J2000: int
FLG_J2000: int
SEFLG_NONUT: int
FLG_NONUT: int
SEFLG_SPEED: int
FLG_SPEED: int
SEFLG_NOGDEFL: int
FLG_NOGDEFL: int
SEFLG_NOABERR: int
FLG_NOABERR: int
SEFLG_BARYCTR: int
FLG_BARYCTR: int
SEFLG_TOPOCTR: int
FLG_TOPOCTR: int
SEFLG_EQUATORIAL: int
FLG_EQUATORIAL: int
SEFLG_XYZ: int
FLG_XYZ: int
SEFLG_RADIANS: int
FLG_RADIANS: int
SEFLG_SIDEREAL: int
FLG_SIDEREAL: int

SE_NODBIT_MEAN: int
NODBIT_MEAN: int
SE_NODBIT_OSCU: int
NODBIT_OSCU: int
SE_NODBIT_OSCU_BAR: int
NODBIT_OSCU_BAR: int
SE_NODBIT_FOPOINT: int
NODBIT_FOPOINT: int

SE_HOR2ECL: int
HOR2ECL: int
SE_HOR2EQU: int
HOR2EQU: int

SE_SIDM_FAGAN_BRADLEY: int
SIDM_FAGAN_BRADLEY: int
SE_SIDM_LAHIRI: int
SIDM_LAHIRI: int
SE_SIDM_DELUCE: int
SIDM_DELUCE: int
SE_SIDM_RAMAN: int
SIDM_RAMAN: int
SE_SIDM_USHASHASHI: int
SIDM_USHASHASHI: int
SE_SIDM_KRISHNAMURTI: int
SIDM_KRISHNAMURTI: int
SE_SIDM_DJWHAL_KHUL: int
SIDM_DJWHAL_KHUL: int
SE_SIDM_YUKTESHWAR: int
SIDM_YUKTESHWAR: int
SE_SIDM_JN_BHASIN: int
SIDM_JN_BHASIN: int
SE_SIDM_BABYL_KUGLER1: int
SIDM_BABYL_KUGLER1: int
SE_SIDM_BABYL_KUGLER2: int
SIDM_BABYL_KUGLER2: int
SE_SIDM_BABYL_KUGLER3: int
SIDM_BABYL_KUGLER3: int
SE_SIDM_BABYL_HUBER: int
SIDM_BABYL_HUBER: int
SE_SIDM_BABYL_ETPSC: int
SIDM_BABYL_ETPSC: int
SE_SIDM_ALDEBARAN_15TAU: int
SIDM_ALDEBARAN_15TAU: int
SE_SIDM_HIPPARCHOS: int
SIDM_HIPPARCHOS: int
SE_SIDM_SASSANIAN: int
SIDM_SASSANIAN: int
SE_SIDM_GALCENT_0SAG: int
SIDM_GALCENT_0SAG: int
SE_SIDM_J2000: int
SIDM_J2000: int
SE_SIDM_J1900: int
SIDM_J1900: int
SE_SIDM_B1950: int
SIDM_B1950: int
SE_SIDM_SURYASIDDHANTA: int
SIDM_SURYASIDDHANTA: int
SE_SIDM_SURYASIDDHANTA_MSUN: int
SIDM_SURYASIDDHANTA_MSUN: int
SE_SIDM_ARYABHATA: int
SIDM_ARYABHATA: int
SE_SIDM_ARYABHATA_MSUN: int
SIDM_ARYABHATA_MSUN: int
SE_SIDM_SS_REVATI: int
SIDM_SS_REVATI: int
SE_SIDM_SS_CITRA: int
SIDM_SS_CITRA: int
SE_SIDM_TRUE_CITRA: int
SIDM_TRUE_CITRA: int
SE_SIDM_TRUE_REVATI: int
SIDM_TRUE_REVATI: int
SE_SIDM_TRUE_PUSHYA: int
SIDM_TRUE_PUSHYA: int
SE_SIDM_GALCENT_RGILBRAND: int
SIDM_GALCENT_RGILBRAND: int
SE_SIDM_GALEQ_IAU1958: int
SIDM_GALEQ_IAU1958: int
SE_SIDM_GALEQ_TRUE: int
SIDM_GALEQ_TRUE: int
SE_SIDM_GALEQ_MULA: int
SIDM_GALEQ_MULA: int
SE_SIDM_GALALIGN_MARDYKS: int
SIDM_GALALIGN_MARDYKS: int
SE_SIDM_TRUE_MULA: int
SIDM_TRUE_MULA: int
SE_SIDM_GALCENT_MULA_WILHELM: int
SIDM_GALCENT_MULA_WILHELM: int
SE_SIDM_ARYABHATA_522: int
SIDM_ARYABHATA_522: int
SE_SIDM_BABYL_BRITTON: int
SIDM_BABYL_BRITTON: int
SE_SIDM_TRUE_SHERATAN: int
SIDM_TRUE_SHERATAN: int
SE_SIDM_TRUE_ALDEBARAN: int
SIDM_TRUE_ALDEBARAN: int
SE_SIDM_GALEQU_FIORENZA: int
SIDM_GALEQU_FIORENZA: int
SE_SIDM_VALENS_MOON: int
SIDM_VALENS_MOON: int
SE_SIDM_LAHIRI_1940: int
SIDM_LAHIRI_1940: int
SE_SIDM_LAHIRI_VP285: int
SIDM_LAHIRI_VP285: int
SE_SIDM_KP_SENTHILATHIBAN: int
SIDM_KP_SENTHILATHIBAN: int
SE_SIDM_LAHIRI_ICRC: int
SIDM_LAHIRI_ICRC: int
SE_SIDM_USER: int
SIDM_USER: int
SE_SIDBIT_ECL_T0: int
SIDBIT_ECL_T0: int
SE_SIDBIT_SSY_PLANE: int
SIDBIT_SSY_PLANE: int
SE_SIDBIT_USER_UT: int
SIDBIT_USER_UT: int
SE_SIDBIT_ECL_DATE: int
SIDBIT_ECL_DATE: int
SE_SIDBIT_NO_PREC_OFFSET: int
SIDBIT_NO_PREC_OFFSET: int
SE_SIDBIT_PREC_ORIG: int
SIDBIT_PREC_ORIG: int

SE_CALC_RISE: int
CALC_RISE: int
SWE_CALC_RISE: int
SE_CALC_SET: int
CALC_SET: int
SWE_CALC_SET: int
SE_CALC_MTRANSIT: int
CALC_MTRANSIT: int
SWE_CALC_MTRANSIT: int
SE_CALC_ITRANSIT: int
CALC_ITRANSIT: int
SWE_CALC_ITRANSIT: int
SE_BIT_DISC_CENTER: int
BIT_DISC_CENTER: int
SWE_BIT_DISC_CENTER: int
SE_BIT_DISC_BOTTOM: int
BIT_DISC_BOTTOM: int
SWE_BIT_DISC_BOTTOM: int
SE_BIT_GEOCTR_NO_ECL_LAT: int
BIT_GEOCTR_NO_ECL_LAT: int
SWE_BIT_GEOCTR_NO_ECL_LAT: int
SE_BIT_NO_REFRACTION: int
BIT_NO_REFRACTION: int
SWE_BIT_NO_REFRACTION: int
SE_BIT_CIVIL_TWILIGHT: int
BIT_CIVIL_TWILIGHT: int
SWE_BIT_CIVIL_TWILIGHT: int
SE_BIT_NAUTIC_TWILIGHT: int
BIT_NAUTIC_TWILIGHT: int
SWE_BIT_NAUTIC_TWILIGHT: int
SE_BIT_ASTRO_TWILIGHT: int
BIT_ASTRO_TWILIGHT: int
SWE_BIT_ASTRO_TWILIGHT: int
SE_BIT_FIXED_DISC_SIZE: int
BIT_FIXED_DISC_SIZE: int
SWE_BIT_FIXED_DISC_SIZE: int
SE_BIT_HINDU_RISING: int
BIT_HINDU_RISING: int
SWE_BIT_HINDU_RISING: int
SE_ECL2HOR: int
ECL2HOR: int
SWE_ECL2HOR: int
SE_EQU2HOR: int
EQU2HOR: int
SWE_EQU2HOR: int
SE_JUL_CAL: int
JUL_CAL: int
SWE_JUL_CAL: int
SE_GREG_CAL: int
GREG_CAL: int
SWE_GREG_CAL: int

SE_ECL_CENTRAL: int
ECL_CENTRAL: int
SWE_ECL_CENTRAL: int
SE_ECL_NONCENTRAL: int
ECL_NONCENTRAL: int
SWE_ECL_NONCENTRAL: int
SE_ECL_TOTAL: int
ECL_TOTAL: int
SWE_ECL_TOTAL: int
SE_ECL_ANNULAR: int
ECL_ANNULAR: int
SWE_ECL_ANNULAR: int
SE_ECL_PARTIAL: int
ECL_PARTIAL: int
SWE_ECL_PARTIAL: int
SE_ECL_ANNULAR_TOTAL: int
ECL_ANNULAR_TOTAL: int
SWE_ECL_ANNULAR_TOTAL: int
SE_ECL_HYBRID: int
ECL_HYBRID: int
SWE_ECL_HYBRID: int
SE_ECL_PENUMBRAL: int
ECL_PENUMBRAL: int
SWE_ECL_PENUMBRAL: int
SE_ECL_ALLTYPES_SOLAR: int
ECL_ALLTYPES_SOLAR: int
SWE_ECL_ALLTYPES_SOLAR: int
SE_ECL_ALLTYPES_LUNAR: int
ECL_ALLTYPES_LUNAR: int
SWE_ECL_ALLTYPES_LUNAR: int
SE_ECL_VISIBLE: int
ECL_VISIBLE: int
SWE_ECL_VISIBLE: int
SE_ECL_MAX_VISIBLE: int
ECL_MAX_VISIBLE: int
SWE_ECL_MAX_VISIBLE: int
SE_ECL_1ST_VISIBLE: int
ECL_1ST_VISIBLE: int
SWE_ECL_1ST_VISIBLE: int
SE_ECL_2ND_VISIBLE: int
ECL_2ND_VISIBLE: int
SWE_ECL_2ND_VISIBLE: int
SE_ECL_3RD_VISIBLE: int
ECL_3RD_VISIBLE: int
SWE_ECL_3RD_VISIBLE: int
SE_ECL_4TH_VISIBLE: int
ECL_4TH_VISIBLE: int
SWE_ECL_4TH_VISIBLE: int
SE_ECL_ONE_TRY: int
ECL_ONE_TRY: int
SWE_ECL_ONE_TRY: int
SE_HELIACAL_RISING: int
HELIACAL_RISING: int
SE_HELIACAL_SETTING: int
HELIACAL_SETTING: int
SE_EVENING_FIRST: int
EVENING_FIRST: int
SE_MORNING_LAST: int
MORNING_LAST: int
SE_HELFLAG_HIGH_PRECISION: int
HELFLAG_HIGH_PRECISION: int

class AttrDict(dict[str, object]):
    def __getattr__(self, key: str) -> object: ...
    def __setattr__(self, key: str, val: object) -> None: ...
    def __delattr__(self, key: str) -> None: ...

class SplitJd:
    @property
    def jd1(self) -> float: ...
    @property
    def jd2(self) -> float: ...
    def __init__(self, jd1: float, jd2: float = 0.0) -> None: ...
    def to_f64(self) -> float: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class Observer:
    @property
    def latitude_deg(self) -> float: ...
    @property
    def longitude_deg(self) -> float: ...
    @property
    def altitude_m(self) -> float: ...
    def __init__(self, latitude_deg: float, longitude_deg: float, altitude_m: float = 0.0) -> None: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class CalcOptions:
    def __init__(
        self,
        origin: str = "geocentric",
        plane: str = "ecliptic",
        equinox: str = "of_date",
        nutation: str = "true",
        aberration: str = "apparent",
        zodiac: str = "tropical",
        ayanamsa: str = "lahiri",
        include_speed: bool = True,
        observer: Optional[Observer] = None,
    ) -> None: ...

class CalcResult:
    @property
    def target_id(self) -> int: ...
    @property
    def center_id(self) -> int: ...
    @property
    def jd_tt(self) -> float: ...
    @property
    def jd_ut1(self) -> float: ...
    @property
    def jd_tdb(self) -> float: ...
    @property
    def delta_t_sec(self) -> float: ...
    @property
    def time_extrapolated(self) -> bool: ...
    @property
    def time_table_version(self) -> str: ...
    @property
    def lon_deg(self) -> float: ...
    @property
    def lat_deg(self) -> float: ...
    @property
    def dist_au(self) -> float: ...
    @property
    def dist_km(self) -> float: ...
    @property
    def lon_speed_deg_day(self) -> float: ...
    @property
    def lat_speed_deg_day(self) -> float: ...
    @property
    def dist_speed_au_day(self) -> float: ...
    @property
    def dist_speed_km_day(self) -> float: ...
    @property
    def magnitude_v(self) -> Optional[float]: ...
    @property
    def cartesian_au(self) -> list[float]: ...
    @property
    def cartesian_vel_au_day(self) -> list[float]: ...
    @property
    def warnings(self) -> list[str]: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class BatchItemResult:
    @property
    def ok(self) -> bool: ...
    @property
    def result(self) -> Optional[CalcResult]: ...
    @property
    def error_code(self) -> Optional[str]: ...
    @property
    def error_message(self) -> Optional[str]: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class DeltaTInfo:
    @property
    def delta_t_sec(self) -> float: ...
    @property
    def tt_minus_tdb_sec(self) -> float: ...
    @property
    def extrapolated(self) -> bool: ...
    @property
    def table_version(self) -> str: ...
    @property
    def warning(self) -> Optional[str]: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class SiderealTimeInfo:
    @property
    def gmst_hours(self) -> float: ...
    @property
    def gast_hours(self) -> float: ...
    @property
    def local_sidereal_hours(self) -> float: ...
    @property
    def armc_deg(self) -> float: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class EclNutInfo:
    @property
    def true_obliquity_deg(self) -> float: ...
    @property
    def mean_obliquity_deg(self) -> float: ...
    @property
    def nutation_lon_deg(self) -> float: ...
    @property
    def nutation_obl_deg(self) -> float: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class StarMetadata:
    @property
    def hip(self) -> int: ...
    @property
    def canonical_name(self) -> str: ...
    @property
    def bayer(self) -> str: ...
    @property
    def vmag(self) -> float: ...
    @property
    def ra_rad(self) -> float: ...
    @property
    def dec_rad(self) -> float: ...
    @property
    def pm_ra_mas_yr(self) -> float: ...
    @property
    def pm_dec_mas_yr(self) -> float: ...
    @property
    def plx_mas(self) -> float: ...
    @property
    def rv_kms(self) -> float: ...
    @property
    def rv_missing(self) -> bool: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class HeliacalResult:
    @property
    def target_name(self) -> str: ...
    @property
    def visual_magnitude(self) -> float: ...
    @property
    def arcus_visionis_deg(self) -> float: ...
    @property
    def heliacal_rise_jd_ut1(self) -> Optional[float]: ...
    @property
    def morning_set_jd_ut1(self) -> Optional[float]: ...
    @property
    def evening_rise_jd_ut1(self) -> Optional[float]: ...
    @property
    def heliacal_set_jd_ut1(self) -> Optional[float]: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class AyanamsaValue:
    @property
    def model_code(self) -> str: ...
    @property
    def jd_tt(self) -> float: ...
    @property
    def mean_deg(self) -> float: ...
    @property
    def true_deg(self) -> float: ...
    @property
    def effective_deg(self) -> float: ...
    @property
    def rate_deg_day(self) -> float: ...
    @property
    def dpsi_arcsec(self) -> float: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class HouseResult:
    @property
    def system_code(self) -> str: ...
    @property
    def effective_system_code(self) -> str: ...
    @property
    def polar_house_fallback(self) -> bool: ...
    @property
    def cusps_deg(self) -> list[float]: ...
    @property
    def cusp_speeds_deg_day(self) -> list[float]: ...
    @property
    def asc_deg(self) -> float: ...
    @property
    def mc_deg(self) -> float: ...
    @property
    def armc_deg(self) -> float: ...
    @property
    def vertex_deg(self) -> float: ...
    @property
    def equatorial_asc_deg(self) -> float: ...
    @property
    def co_asc_koch_deg(self) -> float: ...
    @property
    def co_asc_munkasey_deg(self) -> float: ...
    @property
    def polar_asc_deg(self) -> float: ...
    @property
    def asc_speed_deg_day(self) -> float: ...
    @property
    def mc_speed_deg_day(self) -> float: ...
    @property
    def warnings(self) -> list[str]: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class HorizontalCoords:
    @property
    def azimuth_deg(self) -> float: ...
    @property
    def altitude_true_deg(self) -> float: ...
    @property
    def altitude_apparent_deg(self) -> float: ...
    @property
    def hour_angle_deg(self) -> float: ...
    @property
    def near_horizon(self) -> bool: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class SolarEclipseInfo:
    @property
    def kind(self) -> str: ...
    @property
    def saros_series(self) -> int: ...
    @property
    def saros_member(self) -> int: ...
    @property
    def max_eclipse_jd_tt(self) -> float: ...
    @property
    def max_eclipse_jd_ut1(self) -> float: ...
    @property
    def delta_t_sec(self) -> float: ...
    @property
    def gamma(self) -> float: ...
    @property
    def u_radius_earth(self) -> float: ...
    @property
    def penumbra_radius_earth(self) -> float: ...
    @property
    def magnitude(self) -> float: ...
    @property
    def central_duration_sec(self) -> Optional[float]: ...
    @property
    def path_width_km(self) -> Optional[float]: ...
    @property
    def central_latitude_deg(self) -> Optional[float]: ...
    @property
    def central_longitude_deg(self) -> Optional[float]: ...
    @property
    def p1_jd_tt(self) -> Optional[float]: ...
    @property
    def u1_jd_tt(self) -> Optional[float]: ...
    @property
    def u2_jd_tt(self) -> Optional[float]: ...
    @property
    def u3_jd_tt(self) -> Optional[float]: ...
    @property
    def u4_jd_tt(self) -> Optional[float]: ...
    @property
    def p4_jd_tt(self) -> Optional[float]: ...
    @property
    def local_visible_at_max(self) -> Optional[bool]: ...
    @property
    def local_is_central(self) -> Optional[bool]: ...
    @property
    def local_kind(self) -> Optional[str]: ...
    @property
    def local_c1_jd_tt(self) -> Optional[float]: ...
    @property
    def local_c2_jd_tt(self) -> Optional[float]: ...
    @property
    def local_max_jd_tt(self) -> Optional[float]: ...
    @property
    def local_c3_jd_tt(self) -> Optional[float]: ...
    @property
    def local_c4_jd_tt(self) -> Optional[float]: ...
    @property
    def local_magnitude(self) -> Optional[float]: ...
    @property
    def local_obscuration(self) -> Optional[float]: ...
    @property
    def local_central_duration_sec(self) -> Optional[float]: ...
    @property
    def local_sun_altitude_deg(self) -> Optional[float]: ...
    @property
    def local_sun_azimuth_deg(self) -> Optional[float]: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class LunarEclipseInfo:
    @property
    def kind(self) -> str: ...
    @property
    def saros_series(self) -> int: ...
    @property
    def saros_member(self) -> int: ...
    @property
    def shadow_rule(self) -> str: ...
    @property
    def max_eclipse_jd_tt(self) -> float: ...
    @property
    def max_eclipse_jd_ut1(self) -> float: ...
    @property
    def delta_t_sec(self) -> float: ...
    @property
    def gamma(self) -> float: ...
    @property
    def penumbral_magnitude(self) -> float: ...
    @property
    def umbral_magnitude(self) -> float: ...
    @property
    def penumbral_semi_duration_min(self) -> float: ...
    @property
    def partial_semi_duration_min(self) -> float: ...
    @property
    def total_semi_duration_min(self) -> float: ...
    @property
    def p1_jd_tt(self) -> Optional[float]: ...
    @property
    def u1_jd_tt(self) -> Optional[float]: ...
    @property
    def u2_jd_tt(self) -> Optional[float]: ...
    @property
    def u3_jd_tt(self) -> Optional[float]: ...
    @property
    def u4_jd_tt(self) -> Optional[float]: ...
    @property
    def p4_jd_tt(self) -> Optional[float]: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class PhenomenaResult:
    @property
    def phase_angle_deg(self) -> float: ...
    @property
    def phase_illuminated(self) -> float: ...
    @property
    def elongation_deg(self) -> float: ...
    @property
    def apparent_diameter_deg(self) -> float: ...
    @property
    def apparent_magnitude_v(self) -> float: ...
    @property
    def horizontal_parallax_deg(self) -> float: ...
    @property
    def geocentric_distance_au(self) -> float: ...
    @property
    def heliocentric_distance_au(self) -> float: ...
    def to_swe_attr20(self) -> list[float]: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class NodApsResult:
    @property
    def ascending_node(self) -> list[float]: ...
    @property
    def descending_node(self) -> list[float]: ...
    @property
    def periapsis(self) -> list[float]: ...
    @property
    def apoapsis(self) -> list[float]: ...
    @property
    def asc_node(self) -> list[float]: ...
    @property
    def desc_node(self) -> list[float]: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class LunarOccultationInfo:
    @property
    def kind(self) -> str: ...
    @property
    def target(self) -> str: ...
    @property
    def jd_ut1_max(self) -> float: ...
    @property
    def jd_tt_max(self) -> float: ...
    @property
    def t1_first_contact_jd_ut1(self) -> Optional[float]: ...
    @property
    def t2_second_contact_jd_ut1(self) -> Optional[float]: ...
    @property
    def t3_third_contact_jd_ut1(self) -> Optional[float]: ...
    @property
    def t4_fourth_contact_jd_ut1(self) -> Optional[float]: ...
    @property
    def gamma(self) -> float: ...
    @property
    def magnitude(self) -> float: ...
    @property
    def moon_altitude_deg(self) -> Optional[float]: ...
    @property
    def local_moon_altitude_deg(self) -> Optional[float]: ...
    def to_swe_glob_tret(self) -> list[float]: ...
    def to_swe_tret10(self) -> list[float]: ...
    def swe_retflag(self) -> int: ...
    def to_dict(self) -> AttrDict: ...
    def to_json(self) -> str: ...

class LongitudeCrossingDict(TypedDict):
    body: str
    target_lon_deg: float
    jd_tt: float
    jd_ut1: float
    delta_t_sec: float
    actual_lon_deg: float
    lat_deg: float
    dist_au: float
    dist_km: float
    lon_speed_deg_day: float
    is_retrograde: bool

class SolarTermDict(TypedDict):
    term_index: int
    target_lon_deg: float
    pinyin_name: str
    is_major_term: bool
    jd_tt: float
    jd_ut1: float
    delta_t_sec: float

class SeasonsOfYearDict(TypedDict):
    year: int
    vernal_equinox_jd_ut1: float
    vernal_equinox_jd_tt: float
    summer_solstice_jd_ut1: float
    summer_solstice_jd_tt: float
    autumnal_equinox_jd_ut1: float
    autumnal_equinox_jd_tt: float
    winter_solstice_jd_ut1: float
    winter_solstice_jd_tt: float

class LunarPhaseEventDict(TypedDict):
    kind: Optional[Literal["new_moon", "first_quarter", "full_moon", "last_quarter"]]
    target_elongation_deg: float
    jd_tt: float
    jd_ut1: float
    delta_t_sec: float
    moon_lon_deg: float
    sun_lon_deg: float
    moon_dist_au: float
    moon_dist_km: float

class LunarPhaseSummaryDict(TypedDict):
    jd_tt: float
    jd_ut1: float
    elongation_deg: float
    illumination_fraction: float
    phase_angle_deg: float
    moon_lon_deg: float
    moon_lat_deg: float
    moon_dist_au: float
    moon_dist_km: float
    apparent_diameter_deg: float
    magnitude_v: float
    age_days: float
    synodic_fraction: float
    phase_8_index: int
    phase_8_name: str
    phase_28_index: int
    prev_new_moon_jd_ut1: float
    next_new_moon_jd_ut1: float
    next_first_quarter_jd_ut1: float
    next_full_moon_jd_ut1: float
    next_last_quarter_jd_ut1: float

class StationEventDict(TypedDict):
    body: str
    kind: Literal["station_retrograde", "station_direct"]
    jd_tt: float
    jd_ut1: float
    delta_t_sec: float
    lon_deg: float
    lat_deg: float
    dist_au: float
    dist_km: float

class DeclinationExtremumDict(TypedDict):
    body: str
    kind: Literal["northern_maximum", "southern_minimum"]
    jd_tt: float
    jd_ut1: float
    delta_t_sec: float
    ra_deg: float
    dec_deg: float
    dist_au: float
    dist_km: float

ParanAngleLiteral = Literal["rising", "culminating", "setting", "anti-culminating"]

class ParanEventDict(TypedDict):
    star_name: str
    star_bayer: str
    star_hip: int
    planet: str
    paran_type: ParanAngleLiteral
    star_angle: ParanAngleLiteral
    planet_angle: ParanAngleLiteral
    same_angle: bool
    star_event_jd_ut1: float
    planet_event_jd_ut1: float
    time_orb_minutes: float

class Engine:
    def __init__(
        self,
        data_dir: Optional[str] = None,
        max_threads: Optional[int] = None,
        pin_time_table_version: Optional[str] = None,
        verify_manifest_sha256: bool = False,
    ) -> None: ...
    @staticmethod
    def from_data_dir(
        data_dir: str,
        max_threads: Optional[int] = None,
        pin_time_table_version: Optional[str] = None,
        verify_manifest_sha256: bool = False,
    ) -> Engine: ...
    @staticmethod
    def from_bsp_paths(
        paths: list[str],
        max_threads: Optional[int] = None,
        pin_time_table_version: Optional[str] = None,
    ) -> Engine: ...
    @staticmethod
    def from_bytes(
        kernels: list[tuple[str, bytes]],
        max_threads: Optional[int] = None,
    ) -> Engine: ...
    @staticmethod
    def new_moshier(
        max_threads: Optional[int] = None,
        pin_time_table_version: Optional[str] = None,
    ) -> Engine: ...
    @staticmethod
    def from_data_dir_or_moshier(
        data_dir: str,
        max_threads: Optional[int] = None,
        pin_time_table_version: Optional[str] = None,
        verify_manifest_sha256: bool = False,
    ) -> Engine: ...
    def is_moshier_mode(self) -> bool: ...
    def load_custom_stars(self, path: str) -> int: ...
    def calc_moshier(
        self,
        body: str,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "TT",
        options: Optional[CalcOptions] = None,
    ) -> CalcResult: ...
    def phenomena(
        self,
        body: str,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "TT",
        options: Optional[CalcOptions] = None,
    ) -> PhenomenaResult: ...
    def nod_aps(
        self,
        body: str,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "TT",
        method: int = 0,
        options: Optional[CalcOptions] = None,
    ) -> NodApsResult: ...
    def next_lunar_occultation(
        self,
        body: str,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "TT",
        backward: bool = False,
        observer: Optional[Observer] = None,
        max_synodic_months: int = 36,
    ) -> LunarOccultationInfo: ...
    def next_local_lunar_occultation(
        self,
        body: str,
        jd1: float,
        latitude_deg: float,
        longitude_deg: float,
        jd2: float = 0.0,
        scale: str = "TT",
        backward: bool = False,
        altitude_m: float = 0.0,
        max_synodic_months: int = 120,
    ) -> LunarOccultationInfo: ...
    def lunar_occultations_in_range(
        self,
        body: str,
        start_jd_tt: float,
        end_jd_tt: float,
        observer: Optional[Observer] = None,
    ) -> list[LunarOccultationInfo]: ...
    def local_lunar_occultations_in_range(
        self,
        body: str,
        start_jd_tt: float,
        end_jd_tt: float,
        latitude_deg: float,
        longitude_deg: float,
        altitude_m: float = 0.0,
    ) -> list[LunarOccultationInfo]: ...
    def house_pos(
        self,
        armc_deg: float,
        geolat_deg: float,
        eps_deg: float,
        lon_deg: float,
        lat_deg: float = 0.0,
        hsys: Optional[Union[str, bytes]] = None,
    ) -> float: ...
    def time_equation(
        self,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "UT1",
    ) -> float: ...
    def azalt_rev(
        self,
        jd1: float,
        calc_flag: int,
        latitude_deg: float,
        longitude_deg: float,
        az_deg: float,
        true_alt_deg: float,
        jd2: float = 0.0,
        scale: str = "UT1",
        altitude_m: float = 0.0,
    ) -> list[float]: ...
    def capabilities_json(self) -> str: ...
    def license_info(self) -> AttrDict: ...
    def coverage_json(self, body: str) -> str: ...
    def delta_t(self, jd1: float, jd2: float = 0.0, scale: str = "UT1") -> DeltaTInfo: ...
    def sidtime(
        self,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "UT1",
        nutation: str = "true",
        longitude_deg: float = 0.0,
    ) -> SiderealTimeInfo: ...
    def ecl_nut(self, jd1: float, jd2: float = 0.0, scale: str = "TT") -> EclNutInfo: ...
    def star_info(self, name_or_bayer: str) -> StarMetadata: ...
    def list_stars(self) -> list[StarMetadata]: ...
    def calc(
        self,
        body: str,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "TT",
        options: Optional[CalcOptions] = None,
    ) -> CalcResult: ...
    def calc_star_by_name(
        self,
        name_or_bayer: str,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "TT",
        options: Optional[CalcOptions] = None,
    ) -> CalcResult: ...
    def calc_batch(
        self,
        bodies: list[str],
        epochs_jd: list[float],
        scale: str = "TT",
        options: Optional[CalcOptions] = None,
    ) -> list[list[CalcResult]]: ...
    def calc_batch_flat_f64(
        self,
        bodies: list[str],
        epochs_jd: list[float],
        scale: str = "TT",
        options: Optional[CalcOptions] = None,
    ) -> list[float]: ...
    def calc_batch_per_item(
        self,
        bodies: list[str],
        epochs_jd: list[float],
        scale: str = "TT",
        options: Optional[CalcOptions] = None,
    ) -> list[list[BatchItemResult]]: ...
    def ayanamsa(
        self,
        model: str = "lahiri",
        jd1: float = 2451545.0,
        jd2: float = 0.0,
        scale: str = "TT",
        nutation: str = "true",
    ) -> AyanamsaValue: ...
    def houses(
        self,
        jd1: float,
        latitude_deg: float,
        longitude_deg: float,
        system: Union[str, bytes] = "P",
        jd2: float = 0.0,
        scale: str = "TT",
        altitude_m: float = 0.0,
        polar_policy: str = "strict",
        zodiac: str = "tropical",
        ayanamsa: str = "lahiri",
        nutation: str = "true",
        include_speed: bool = True,
    ) -> HouseResult: ...
    def houses_batch(
        self,
        epochs_jd: list[float],
        latitude_deg: float,
        longitude_deg: float,
        system: Union[str, bytes] = "P",
        scale: str = "TT",
        altitude_m: float = 0.0,
        polar_policy: str = "strict",
        zodiac: str = "tropical",
        ayanamsa: str = "lahiri",
        nutation: str = "true",
        include_speed: bool = True,
    ) -> list[HouseResult]: ...
    def azalt(
        self,
        ra_deg: float,
        dec_deg: float,
        jd1: float,
        latitude_deg: float,
        longitude_deg: float,
        jd2: float = 0.0,
        scale: str = "TT",
        altitude_m: float = 0.0,
        pressure_hpa: float = 1013.25,
        temperature_c: float = 15.0,
    ) -> HorizontalCoords: ...
    def azalt_ecliptic(
        self,
        ecl_lon_deg: float,
        ecl_lat_deg: float,
        jd1: float,
        latitude_deg: float,
        longitude_deg: float,
        jd2: float = 0.0,
        scale: str = "TT",
        altitude_m: float = 0.0,
        pressure_hpa: float = 1013.25,
        temperature_c: float = 15.0,
    ) -> HorizontalCoords: ...
    def calc_azalt(
        self,
        body: str,
        jd1: float,
        latitude_deg: float,
        longitude_deg: float,
        jd2: float = 0.0,
        scale: str = "TT",
        altitude_m: float = 0.0,
        pressure_hpa: float = 1013.25,
        temperature_c: float = 15.0,
    ) -> HorizontalCoords: ...
    def rise_trans(
        self,
        body: str,
        event: str,
        jd1: float,
        latitude_deg: float,
        longitude_deg: float,
        jd2: float = 0.0,
        scale: str = "UT1",
        altitude_m: float = 0.0,
        custom_horizon_alt_deg: Optional[float] = None,
        max_search_days: float = 2.0,
        disc_mode: str = "upper_limb",
        pressure_hpa: float = 1013.25,
        temperature_c: float = 15.0,
    ) -> Optional[float]: ...
    def heliacal_phenomena(
        self,
        body: str,
        jd1: float,
        latitude_deg: float,
        longitude_deg: float,
        jd2: float = 0.0,
        scale: str = "UT1",
        altitude_m: float = 0.0,
        extinction_k_v: float = 0.20,
        snellen_factor: float = 1.0,
        observer_age_years: float = 36.0,
        snellen_ratio: Optional[float] = None,
        pressure_hpa: float = 1013.25,
        temperature_c: float = 15.0,
        relative_humidity: float = 0.50,
        max_search_days: float = 400.0,
    ) -> HeliacalResult: ...
    def next_solar_eclipse(
        self,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "TT",
        backward: bool = False,
        observer: Optional[Observer] = None,
        max_synodic_months: int = 36,
    ) -> SolarEclipseInfo: ...
    def prev_solar_eclipse(
        self,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "TT",
        observer: Optional[Observer] = None,
        max_synodic_months: int = 36,
    ) -> SolarEclipseInfo: ...
    def next_local_solar_eclipse(
        self,
        jd1: float,
        latitude_deg: float,
        longitude_deg: float,
        jd2: float = 0.0,
        scale: str = "TT",
        backward: bool = False,
        altitude_m: float = 0.0,
        max_synodic_months: int = 120,
    ) -> SolarEclipseInfo: ...
    def solar_eclipses_in_range(
        self,
        start_jd_tt: float,
        end_jd_tt: float,
        observer: Optional[Observer] = None,
    ) -> list[SolarEclipseInfo]: ...
    def next_lunar_eclipse(
        self,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "TT",
        backward: bool = False,
        shadow_rule: str = "danjon",
        max_synodic_months: int = 36,
    ) -> LunarEclipseInfo: ...
    def prev_lunar_eclipse(
        self,
        jd1: float,
        jd2: float = 0.0,
        scale: str = "TT",
        shadow_rule: str = "danjon",
        max_synodic_months: int = 36,
    ) -> LunarEclipseInfo: ...
    def lunar_eclipses_in_range(
        self,
        start_jd_tt: float,
        end_jd_tt: float,
        shadow_rule: str = "danjon",
    ) -> list[LunarEclipseInfo]: ...
    def search_body_longitude(
        self,
        body: str,
        jd_ut1: float,
        target_lon_deg: float,
        forward: bool = True,
        max_search_days: float = 400.0,
        options: Optional[CalcOptions] = None,
        scale: str = "UT1",
        jd2: float = 0.0,
    ) -> AttrDict: ...
    def body_longitude_crossings_in_range(
        self,
        body: str,
        start_jd_ut1: float,
        end_jd_ut1: float,
        target_lon_deg: float,
        options: Optional[CalcOptions] = None,
        scale: str = "UT1",
    ) -> list[AttrDict]: ...
    def longitude_crossings_in_range(
        self,
        body: str,
        start_jd_ut1: float,
        end_jd_ut1: float,
        target_lon_deg: float,
        options: Optional[CalcOptions] = None,
        scale: str = "UT1",
    ) -> list[AttrDict]: ...
    def solar_terms_in_range(
        self,
        start_jd_ut1: float,
        end_jd_ut1: float,
        scale: str = "UT1",
    ) -> list[AttrDict]: ...
    def seasons_in_year(self, year: int) -> AttrDict: ...
    def human_design_epoch(
        self,
        birth_jd_ut1: float,
        scale: str = "UT1",
        jd2: float = 0.0,
    ) -> AttrDict: ...
    def search_lunar_phase(
        self,
        start_jd_ut1: float,
        target_elongation_deg: float,
        forward: bool = True,
        scale: str = "UT1",
        jd2: float = 0.0,
    ) -> AttrDict: ...
    def lunar_phases_in_range(
        self,
        start_jd_ut1: float,
        end_jd_ut1: float,
        scale: str = "UT1",
    ) -> list[AttrDict]: ...
    def lunar_phase_summary(
        self,
        jd_ut1: float,
        scale: str = "UT1",
        jd2: float = 0.0,
    ) -> AttrDict: ...
    def lunar_phase_detail(
        self,
        jd_ut1: float,
        scale: str = "UT1",
        jd2: float = 0.0,
    ) -> AttrDict: ...
    def search_stations_in_range(
        self,
        body: str,
        start_jd_ut1: float,
        end_jd_ut1: float,
        options: Optional[CalcOptions] = None,
        scale: str = "UT1",
    ) -> list[AttrDict]: ...
    def stations_in_range(
        self,
        body: str,
        start_jd_ut1: float,
        end_jd_ut1: float,
        options: Optional[CalcOptions] = None,
        scale: str = "UT1",
    ) -> list[AttrDict]: ...
    def search_declination_extrema_in_range(
        self,
        body: str,
        start_jd_ut1: float,
        end_jd_ut1: float,
        options: Optional[CalcOptions] = None,
        scale: str = "UT1",
    ) -> list[AttrDict]: ...
    def find_parans(
        self,
        jd_ut1: float,
        observer: Observer,
        orb_minutes: float = 4.0,
        bodies: Optional[list[str]] = None,
        star_names: Optional[list[str]] = None,
        scale: str = "UT1",
        jd2: float = 0.0,
    ) -> list[AttrDict]: ...
    def map_zodiac_placement(self, lon_deg: float) -> AttrDict: ...
    def planetary_return(
        self,
        body: str,
        natal_jd: float,
        search_after_jd: float,
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
        options: Optional[CalcOptions] = None,
    ) -> AttrDict: ...
    def planetary_returns_in_range(
        self,
        body: str,
        natal_jd: float,
        start_jd: float,
        end_jd: float,
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
        options: Optional[CalcOptions] = None,
    ) -> list[AttrDict]: ...
    def calc_chart_aspects(
        self,
        jd: float,
        bodies: Optional[list[str]] = None,
        preset: str = "major",
        orb_scale: float = 1.0,
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
    ) -> list[AttrDict]: ...
    def calc_synastry_aspects(
        self,
        jd1: float,
        bodies1: Optional[list[str]],
        jd2: float,
        bodies2: Optional[list[str]],
        preset: str = "major",
        orb_scale: float = 1.0,
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
    ) -> list[AttrDict]: ...
    def calc_declination_aspects(
        self,
        jd: float,
        bodies: Optional[list[str]] = None,
        max_orb_deg: float = 1.2,
        scale: str = "UT1",
    ) -> list[AttrDict]: ...
    def search_exact_aspect(
        self,
        body1: str,
        body2: str,
        target_angle_deg: float,
        start_jd: float,
        forward: bool = True,
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
    ) -> AttrDict: ...
    def mutual_aspects_in_range(
        self,
        start_jd: float,
        end_jd: float,
        bodies: Optional[list[str]] = None,
        preset: str = "major",
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
    ) -> list[AttrDict]: ...
    def dynamic_aspects_in_range(
        self,
        start_jd: float,
        end_jd: float,
        bodies: Optional[list[str]] = None,
        preset: str = "major",
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
    ) -> list[AttrDict]: ...
    def natal_transits_in_range(
        self,
        natal_jd: float,
        start_jd: float,
        end_jd: float,
        transiting_bodies: Optional[list[str]] = None,
        natal_bodies: Optional[list[str]] = None,
        preset: str = "major",
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
    ) -> list[AttrDict]: ...
    def void_of_course_moon_in_range(
        self,
        start_jd: float,
        end_jd: float,
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
    ) -> list[AttrDict]: ...
    def panchanga(
        self,
        jd: float,
        lat_deg: float,
        lon_deg: float,
        alt_m: float = 0.0,
        ayanamsa: str = "lahiri",
        scale: str = "UT1",
    ) -> AttrDict: ...
    def vimshottari_dasha(
        self,
        birth_jd: float,
        ayanamsa: str = "lahiri",
        year_length_days: Optional[float] = None,
        scale: str = "UT1",
    ) -> AttrDict: ...
    def calc_varga(self, sidereal_lon_deg: float, division: int) -> AttrDict: ...
    def calc_shodashavarga(self, sidereal_lon_deg: float) -> list[AttrDict]: ...
    def chart_sect(
        self,
        jd: float,
        lat_deg: float,
        lon_deg: float,
        alt_m: float = 0.0,
        scale: str = "UT1",
    ) -> AttrDict: ...
    def planetary_hours(
        self,
        jd: float,
        lat_deg: float,
        lon_deg: float,
        alt_m: float = 0.0,
        scale: str = "UT1",
    ) -> AttrDict: ...
    def calc_lots(
        self,
        jd: float,
        lat_deg: float,
        lon_deg: float,
        alt_m: float = 0.0,
        system: Optional[Union[str, bytes]] = None,
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
    ) -> list[AttrDict]: ...
    def chart_snapshot(
        self,
        jd: float,
        lat_deg: float,
        lon_deg: float,
        alt_m: float = 0.0,
        bodies: Optional[list[str]] = None,
        system: Optional[Union[str, bytes]] = None,
        preset: str = "major",
        orb_scale: float = 1.0,
        scale: str = "UT1",
        zodiac: str = "tropical",
        ayanamsa: Optional[str] = None,
    ) -> AttrDict: ...

class SweAdapter:
    version: str
    __version__: str
    Error = LxEphemerisError
    def __init__(self, data_dir: Optional[Union[str, Engine]] = None) -> None: ...
    @staticmethod
    def from_engine(engine: Engine) -> SweAdapter: ...
    @staticmethod
    def from_data_dir(data_dir: str) -> SweAdapter: ...
    @staticmethod
    def new_moshier() -> SweAdapter: ...
    def set_ephe_path(self, path: Optional[str] = None) -> None: ...
    def close(self) -> None: ...
    @staticmethod
    def get_planet_name(ipl: int) -> str: ...
    @staticmethod
    def house_name(hsys: Optional[Union[str, bytes]] = None) -> str: ...
    def set_sid_mode(self, sid_mode: int, t0: float = 0.0, ayan_t0: float = 0.0, sid_bits: Optional[int] = None) -> None: ...
    def set_topo(self, lon_deg: float, lat_deg: float, alt_m: float = 0.0, lon: Optional[float] = None, lat: Optional[float] = None, alt: Optional[float] = None) -> None: ...
    def set_polar_policy(self, policy: str) -> None: ...
    def deltat(self, jd_ut: float) -> float: ...
    def deltat_sec(self, jd_ut: float) -> float: ...
    def sidtime(self, jd_ut: float) -> float: ...
    def sidtime0(self, jd_ut: float, eps_deg: float, nut_lon_deg: float) -> float: ...
    def get_ayanamsa_ex_ut(self, jd_ut: float, iflag: int = 2, flags: Optional[int] = None) -> tuple[int, float]: ...
    def get_ayanamsa_ex(self, jd_et: float, iflag: int = 2, flags: Optional[int] = None) -> tuple[int, float]: ...
    def get_ayanamsa_ut(self, jd_ut: float) -> float: ...
    def get_ayanamsa(self, jd_et: float) -> float: ...
    def pheno_ut(self, jd_ut: float, ipl: int, iflag: int = 0, flags: Optional[int] = None) -> tuple[list[float], int]: ...
    def pheno(self, jd_et: float, ipl: int, iflag: int = 0, flags: Optional[int] = None) -> tuple[list[float], int]: ...
    def nod_aps_ut(
        self,
        jd_ut: float,
        ipl: int,
        iflag: Optional[int] = None,
        method: Optional[int] = None,
        flags: Optional[int] = None,
        method_bits: Optional[int] = None,
    ) -> tuple[list[float], list[float], list[float], list[float], int]: ...
    def nod_aps(
        self,
        jd_et: float,
        ipl: int,
        iflag: Optional[int] = None,
        method: Optional[int] = None,
        flags: Optional[int] = None,
        method_bits: Optional[int] = None,
    ) -> tuple[list[float], list[float], list[float], list[float], int]: ...
    def lun_occult_when_glob(
        self,
        jd_ut: float,
        body: Union[int, str],
        starname: Optional[Union[str, int]] = None,
        iflag: Optional[Union[int, bool]] = None,
        ifltype: Optional[Union[int, bool]] = None,
        backward: Optional[bool] = None,
        flags: Optional[int] = None,
        backwards: Optional[bool] = None,
    ) -> tuple[int, list[float]]: ...
    def lun_occult_when_loc(
        self,
        jd_ut: float,
        body: Union[int, str],
        geopos: Union[list[float], int, str],
        starname: Optional[Union[str, int, list[float]]] = None,
        iflag: Optional[Union[int, bool]] = None,
        backward: Optional[bool] = None,
        flags: Optional[int] = None,
        backwards: Optional[bool] = None,
    ) -> tuple[int, list[float], list[float]]: ...
    def sol_eclipse_when_glob(
        self,
        jd_start: float,
        iflag: int = 2,
        ifltype: int = 0,
        backward: bool = False,
        flags: Optional[int] = None,
        backwards: Optional[bool] = None,
    ) -> tuple[int, list[float]]: ...
    def sol_eclipse_where(
        self,
        jd_ut: float,
        iflag: int = 2,
        flags: Optional[int] = None,
    ) -> tuple[int, list[float], list[float]]: ...
    def sol_eclipse_when_loc(
        self,
        jd_start: float,
        geopos: list[float],
        iflag: int = 2,
        backward: bool = False,
        flags: Optional[int] = None,
        backwards: Optional[bool] = None,
    ) -> tuple[int, list[float], list[float]]: ...
    def lun_eclipse_when(
        self,
        jd_start: float,
        iflag: int = 2,
        ifltype: int = 0,
        backward: bool = False,
        flags: Optional[int] = None,
        backwards: Optional[bool] = None,
    ) -> tuple[int, list[float]]: ...
    def lun_eclipse_how(
        self,
        jd_ut: float,
        geopos: list[float],
        iflag: int = 2,
        flags: Optional[int] = None,
    ) -> tuple[int, list[float]]: ...
    def heliacal_ut(
        self,
        jd_start: float,
        geopos: list[float],
        datm: list[float],
        dobs: list[float],
        object_name: Union[str, int],
        event_type: int,
        helflag: int = 0,
    ) -> tuple[float, float, float]: ...
    def calc_ut(self, jd_ut: float, ipl: int, iflag: int = 258, flags: Optional[int] = None) -> tuple[list[float], int]: ...
    def calc(self, jd_et: float, ipl: int, iflag: int = 258, flags: Optional[int] = None) -> tuple[list[float], int]: ...
    def fixstar2_ut(self, star: str, jd_ut: float, iflag: int = 258, flags: Optional[int] = None) -> tuple[list[float], str, int]: ...
    def fixstar2_mag(self, star: str) -> float: ...
    def house_pos(
        self,
        armc_deg: float,
        geolat_deg: float,
        eps_deg: float,
        xpin: list[float],
        hsys: Optional[Union[str, bytes]] = None,
    ) -> float: ...
    def houses_ex(
        self,
        jd_ut: float,
        lat_deg: float,
        lon_deg: float,
        hsys: Union[str, bytes] = "P",
        iflag: int = 0,
    ) -> tuple[list[float], list[float]]: ...
    def houses(
        self,
        jd_ut: float,
        lat_deg: float,
        lon_deg: float,
        hsys: Union[str, bytes] = "P",
    ) -> tuple[list[float], list[float]]: ...
    def houses_armc(
        self,
        armc_deg: float,
        lat_deg: float,
        eps_deg: float,
        hsys: Union[str, bytes] = "P",
        sun_decl_deg: Optional[float] = None,
    ) -> tuple[list[float], list[float]]: ...
    def time_equ(self, jd_ut: float) -> float: ...
    def azalt(
        self,
        jd_ut: float,
        calc_flag: int,
        geopos: list[float],
        atpress: float,
        attemp: float,
        xin: list[float],
    ) -> list[float]: ...
    def azalt_rev(
        self,
        jd_ut: float,
        calc_flag: int,
        geopos: list[float],
        xin: Union[list[float], float],
        true_alt_deg: Optional[float] = None,
    ) -> list[float]: ...
    def rise_trans(
        self,
        jd_ut: float,
        ipl: Union[int, str],
        rsmi: Optional[Union[int, str]] = None,
        geopos: Optional[Union[list[float], int]] = None,
        atpress: Union[float, int] = 1013.25,
        attemp: Union[float, list[float]] = 15.0,
        starname: Optional[Union[str, float]] = None,
        epheflag: Union[int, float] = 258,
    ) -> tuple[int, tuple[float, float, float, float, float, float, float, float, float, float]]: ...
    @staticmethod
    def julday(
        year: int,
        month: int,
        day: int,
        hour: float = 0.0,
        gregflag: Optional[Union[int, str]] = None,
    ) -> float: ...
    @staticmethod
    def revjul(
        jd: float,
        gregflag: Optional[Union[int, str]] = None,
    ) -> tuple[int, int, int, float]: ...

def julday(
    year: int,
    month: int,
    day: int,
    hour_ut: float = 0.0,
    calendar: Optional[Union[int, str]] = None,
) -> float: ...
def revjul(
    jd: float,
    calendar: Optional[Union[int, str]] = None,
) -> tuple[int, int, int, float]: ...
def utc_from_local(
    year: int,
    month: int,
    day: int,
    hour: int = 0,
    minute: int = 0,
    second: float = 0.0,
    tz_name: str = "UTC",
    fold: int = 0,
) -> tuple[float, int, int, int, float]: ...
def jd_utc_from_local(
    year: int,
    month: int,
    day: int,
    hour: int = 0,
    minute: int = 0,
    second: float = 0.0,
    tz_name: str = "UTC",
    fold: int = 0,
    gregorian: bool = True,
) -> float: ...
def map_zodiac_placement(lon_deg: float) -> AttrDict: ...
def calc_varga(sidereal_lon_deg: float, division: int) -> AttrDict: ...
def calc_shodashavarga(sidereal_lon_deg: float) -> list[AttrDict]: ...
def aspect_definitions(preset: str = "major") -> list[AttrDict]: ...
def houses_armc(
    armc_deg: float,
    eps_true_deg: float,
    latitude_deg: float,
    system: Union[str, bytes] = "P",
    polar_policy: str = "strict",
    ayanamsa_deg: float = 0.0,
    sun_decl_deg: float = 0.0,
    include_speed: bool = True,
) -> HouseResult: ...
def house_of_lon(lon_deg: float, cusps_deg: list[float], eps_deg: float = 1e-9) -> int: ...
def star_info(name_or_bayer: str) -> StarMetadata: ...
def list_stars() -> list[StarMetadata]: ...
def gauquelin_sectors(
    armc_deg: float,
    eps_true_deg: float,
    latitude_deg: float,
    polar_policy: str = "strict",
) -> list[float]: ...
def verify_license_key(token: str) -> AttrDict: ...
def register_commercial_license(token: str) -> AttrDict: ...
def clear_commercial_license() -> None: ...
def license_info() -> AttrDict: ...
