//! Time scales (`UTC`, `UT1`, `TT`, `TDB`), `SplitJd`, calendar conversions (`julday` / `revjul`),
//! two-tier $\Delta T$ / `UT1-UTC` / leap seconds tables, `TDB-TT` (`sofars::ts::dtdb`), and
//! sidereal time (`GMST`, `GAST`, `ERA` strictly from `UT1`).

use crate::error::EphemerisError;
use crate::manifest::{sha256_hex, TimeSyncManifest};
use crate::types::{Observer, JD_J2000, SECONDS_PER_DAY};
use serde::{Deserialize, Serialize};
use sofars::coords::gd2gc;
use sofars::erst::{era00, gmst06, gst06};
use sofars::ts::dtdb;
use std::fs;
use std::path::Path;

/// Two-part Julian Date `(jd1, jd2)` preserving sub-microsecond precision across centuries.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct SplitJd {
    pub jd1: f64,
    pub jd2: f64,
}

impl SplitJd {
    /// Creates a normalized [`SplitJd`] where `jd1` is on a `.5` or `.0` day boundary and `|jd2| <= 0.5`.
    #[must_use]
    pub fn new(jd1: f64, jd2: f64) -> Self {
        if !jd1.is_finite() || !jd2.is_finite() {
            return Self {
                jd1: f64::NAN,
                jd2: f64::NAN,
            };
        }
        let sum = jd1 + jd2;
        if !sum.is_finite() {
            return Self {
                jd1: f64::NAN,
                jd2: f64::NAN,
            };
        }
        let base = (sum - 0.5).floor() + 0.5;
        let rem = (jd1 - base) + jd2;
        Self {
            jd1: base,
            jd2: rem,
        }
    }

    /// Returns `true` if both `jd1` and `jd2` are finite and `|jd1 + jd2| <= 1e9`.
    #[must_use]
    pub fn is_valid_jd(self) -> bool {
        let jd = self.to_f64();
        self.jd1.is_finite() && self.jd2.is_finite() && jd.is_finite() && jd.abs() <= 1.0e9
    }

    /// Creates a [`SplitJd`] from a single `f64` Julian Date.
    #[must_use]
    pub fn from_f64(jd: f64) -> Self {
        Self::new(jd, 0.0)
    }

    /// Returns the combined `f64` Julian Date (`jd1 + jd2`).
    #[must_use]
    pub const fn to_f64(self) -> f64 {
        self.jd1 + self.jd2
    }

    /// Returns a new [`SplitJd`] shifted by `delta_days`.
    #[must_use]
    pub fn add_days(self, delta_days: f64) -> Self {
        Self::new(self.jd1, self.jd2 + delta_days)
    }
}

/// Coordinated Universal Time (`UTC`) Julian Date.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Utc(pub SplitJd);

impl Utc {
    #[must_use]
    pub fn from_jd(jd: f64) -> Self {
        Self(SplitJd::from_f64(jd))
    }

    #[must_use]
    pub const fn jd_f64(self) -> f64 {
        self.0.to_f64()
    }
}

/// Universal Time (`UT1`, Earth rotation angle time scale) Julian Date.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Ut1(pub SplitJd);

impl Ut1 {
    #[must_use]
    pub fn from_jd(jd: f64) -> Self {
        Self(SplitJd::from_f64(jd))
    }

    #[must_use]
    pub const fn jd_f64(self) -> f64 {
        self.0.to_f64()
    }
}

/// Terrestrial Time (`TT = TAI + 32.184 s`) Julian Date.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Tt(pub SplitJd);

impl Tt {
    #[must_use]
    pub fn from_jd(jd: f64) -> Self {
        Self(SplitJd::from_f64(jd))
    }

    #[must_use]
    pub const fn jd_f64(self) -> f64 {
        self.0.to_f64()
    }
}

/// Barycentric Dynamical Time (`TDB`) Julian Date (ephemeris argument for JPL DE kernels).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Tdb(pub SplitJd);

impl Tdb {
    #[must_use]
    pub fn from_jd(jd: f64) -> Self {
        Self(SplitJd::from_f64(jd))
    }

    #[must_use]
    pub const fn jd_f64(self) -> f64 {
        self.0.to_f64()
    }
}

/// Explicit time scale input for `Engine` / `Moment` queries.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum TimeInput {
    Tt(SplitJd),
    Ut1(SplitJd),
    Utc(SplitJd),
    Tdb(SplitJd),
}

impl TimeInput {
    #[must_use]
    pub fn from_tt_jd(jd_tt: f64) -> Self {
        Self::Tt(SplitJd::from_f64(jd_tt))
    }

    #[must_use]
    pub fn from_ut1_jd(jd_ut1: f64) -> Self {
        Self::Ut1(SplitJd::from_f64(jd_ut1))
    }

    #[must_use]
    pub fn from_utc_jd(jd_utc: f64) -> Self {
        Self::Utc(SplitJd::from_f64(jd_utc))
    }

    #[must_use]
    pub fn from_tdb_jd(jd_tdb: f64) -> Self {
        Self::Tdb(SplitJd::from_f64(jd_tdb))
    }
}

/// Calendar system convention for `julday` and `revjul`.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum CalendarKind {
    /// Proleptic Gregorian calendar for all dates.
    #[default]
    Gregorian,
    /// Proleptic Julian calendar for all dates.
    Julian,
    /// Historical changeover: Julian through `1582-10-04` (`JD < 2299160.5`),
    /// Gregorian from `1582-10-15` (`JD >= 2299160.5`).
    AstroAuto,
}

/// Decomposed calendar date and fractional hour.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CalendarDate {
    pub year: i32,
    pub month: u8,
    pub day: u8,
    pub hour: f64,
}

/// Converts a calendar date `(year, month, day, hour)` into a [`SplitJd`] using astronomical year
/// numbering (`0 = 1 BCE`, `-1 = 2 BCE`; Meeus *Astronomical Algorithms* Ch. 7).
pub fn julday_split(
    year: i32,
    month: u8,
    day: u8,
    hour: f64,
    calendar: CalendarKind,
) -> Result<SplitJd, EphemerisError> {
    if !(-1_000_000..=1_000_000).contains(&year)
        || !(1..=12).contains(&month)
        || !(1..=31).contains(&day)
        || !hour.is_finite()
        || !(-10_000.0..=10_000.0).contains(&hour)
    {
        return Err(EphemerisError::InvalidInput {
            message: format!("Invalid calendar date {year}-{month:02}-{day:02} hour={hour}"),
        });
    }
    if calendar == CalendarKind::AstroAuto && year == 1582 && month == 10 && (5..=14).contains(&day)
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Date 1582-10-{day:02} does not exist in the historical civil calendar (Papal Gregorian reform gap)"
            ),
        });
    }

    let use_gregorian = match calendar {
        CalendarKind::Gregorian => true,
        CalendarKind::Julian => false,
        CalendarKind::AstroAuto => {
            (year, month, day) >= (1582, 10, 15)
        }
    };

    let (y, m) = if month <= 2 {
        (i64::from(year) - 1, i64::from(month) + 12)
    } else {
        (i64::from(year), i64::from(month))
    };

    let b = if use_gregorian {
        let a = y.div_euclid(100);
        2 - a + a.div_euclid(4)
    } else {
        0
    };

    let jd_int = (365.25 * ((y + 4716) as f64)).floor()
        + (30.6001 * ((m + 1) as f64)).floor()
        + f64::from(day)
        + (b as f64)
        - 1524.5;

    Ok(SplitJd::new(jd_int, hour / 24.0))
}

/// Converts a calendar date `(year, month, day, hour)` into a single `f64` Julian Date.
pub fn julday(
    year: i32,
    month: u8,
    day: u8,
    hour: f64,
    calendar: CalendarKind,
) -> Result<f64, EphemerisError> {
    julday_split(year, month, day, hour, calendar).map(SplitJd::to_f64)
}

/// Inverse calendar conversion (`revjul`) from a Julian Date to [`CalendarDate`] (Meeus Ch. 7).
pub fn revjul(jd: f64, calendar: CalendarKind) -> Result<CalendarDate, EphemerisError> {
    if !jd.is_finite() || !(-1.0e8..=1.0e8).contains(&jd) {
        return Err(EphemerisError::InvalidInput {
            message: format!("Invalid Julian Date {jd} in revjul"),
        });
    }

    let shifted = jd + 0.5;
    let z_orig = shifted.floor() as i64;
    let f = shifted - (z_orig as f64);

    let use_gregorian = match calendar {
        CalendarKind::Gregorian => true,
        CalendarKind::Julian => false,
        CalendarKind::AstroAuto => z_orig >= 2_299_161,
    };

    // Fold negative Julian Dates forward by an integer multiple of 400-year calendar cycles
    // (146,097 days in Gregorian; 146,100 days in Julian) so Meeus Ch. 7 integer arithmetic
    // always operates in the strictly positive B >= 1524 regime.
    let (z, year_offset) = if z_orig < 0 {
        let cycle_days = if use_gregorian { 146_097_i64 } else { 146_100_i64 };
        let cycles = ((-z_orig) / cycle_days) + 1;
        (z_orig + cycles * cycle_days, cycles * 400)
    } else {
        (z_orig, 0)
    };

    let a = if use_gregorian {
        let alpha = (((z as f64) - 1_867_216.25) / 36_524.25).floor() as i64;
        z + 1 + alpha - alpha.div_euclid(4)
    } else {
        z
    };

    let b = a + 1524;
    let c = (((b as f64) - 122.1) / 365.25).floor() as i64;
    let d = (365.25 * (c as f64)).floor() as i64;
    let e = (((b - d) as f64) / 30.6001).floor() as i64;

    let day_num = b - d - ((30.6001 * (e as f64)).floor() as i64);
    let month_num = if e < 14 { e - 1 } else { e - 13 };
    let year_num = (if month_num > 2 { c - 4716 } else { c - 4715 }) - year_offset;

    Ok(CalendarDate {
        year: year_num as i32,
        month: month_num as u8,
        day: day_num as u8,
        hour: f * 24.0,
    })
}

/// Result of a $\Delta T = \text{TT} - \text{UT1}$ lookup, including provenance metadata (`R1.2`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DeltaTInfo {
    /// $\Delta T = \text{TT} - \text{UT1}$ in seconds.
    pub delta_t_sec: f64,
    /// $\text{TT} - \text{TDB}$ in seconds at this epoch (`-sofars::ts::dtdb`).
    pub tt_minus_tdb_sec: f64,
    /// `true` if the epoch fell outside empirical IERS/historical spline coverage and used
    /// Morrison & Stephenson parabolic extrapolation.
    pub time_extrapolated: bool,
    /// Identifier of the active time table (`"embedded-smh2016-iers2026"` or synced IERS version).
    pub time_table_version: String,
}

/// Greenwich and Local Sidereal Time information (`R1.4`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SiderealTimeInfo {
    /// Greenwich Mean Sidereal Time (`GMST`, IAU 2006) in hours `[0, 24)`.
    pub gmst_hours: f64,
    /// Greenwich Apparent Sidereal Time (`GAST`, IAU 2006/2000A) in hours `[0, 24)`.
    pub gast_hours: f64,
    /// Local Sidereal Time (`GAST` if `nutation == True`, `GMST` if `nutation == Mean`, plus `longitude_deg / 15`) in hours `[0, 24)`.
    pub local_sidereal_hours: f64,
    /// Right Ascension of the Midheaven (`ARMC = local_sidereal_hours * 15`) in degrees `[0, 360)`.
    pub armc_deg: f64,
}

/// Embedded historical & modern $\Delta T$ table knots `(decimal_year, delta_t_sec)`
/// from Stephenson, Morrison & Hohenkerk (2016, Proc. R. Soc. A 472, 20160404, Table 1) and USNO/IERS series.
const EMBEDDED_DELTAT_KNOTS: &[(f64, f64)] = &[
    (-700.0, 20270.0),
    (-600.0, 18450.0),
    (-500.0, 16800.0),
    (-400.0, 15300.0),
    (-300.0, 13890.0),
    (-200.0, 12560.0),
    (-100.0, 11320.0),
    (0.0, 10160.0),
    (100.0, 9070.0),
    (200.0, 8040.0),
    (300.0, 7020.0),
    (400.0, 6030.0),
    (500.0, 5080.0),
    (600.0, 4170.0),
    (700.0, 3300.0),
    (800.0, 2530.0),
    (900.0, 1900.0),
    (1000.0, 1440.0),
    (1100.0, 1090.0),
    (1200.0, 790.0),
    (1300.0, 530.0),
    (1400.0, 320.0),
    (1500.0, 175.0),
    (1600.0, 110.0),
    (1620.0, 92.0),
    (1640.0, 62.0),
    (1660.0, 22.0),
    (1680.0, 16.0),
    (1700.0, 9.0),
    (1720.0, 11.0),
    (1740.0, 12.0),
    (1760.0, 15.0),
    (1780.0, 17.0),
    (1800.0, 13.7),
    (1810.0, 12.5),
    (1820.0, 12.0),
    (1830.0, 7.5),
    (1840.0, 5.7),
    (1850.0, 7.1),
    (1860.0, 7.9),
    (1870.0, 1.6),
    (1880.0, -5.4),
    (1890.0, -6.3),
    (1900.0, -2.72),
    (1905.0, 3.86),
    (1910.0, 10.46),
    (1915.0, 17.20),
    (1920.0, 21.16),
    (1925.0, 23.62),
    (1930.0, 24.02),
    (1935.0, 23.93),
    (1940.0, 24.33),
    (1945.0, 26.77),
    (1950.0, 29.15),
    (1955.0, 31.07),
    (1960.0, 33.15),
    (1965.0, 35.73),
    (1970.0, 40.18),
    (1972.0, 42.23),
    (1975.0, 45.48),
    (1980.0, 50.54),
    (1985.0, 54.34),
    (1990.0, 56.86),
    (1995.0, 60.78),
    (2000.0, 63.83),
    (2005.0, 64.69),
    (2010.0, 66.07),
    (2015.0, 67.64),
    (2017.0, 68.80),
    (2020.0, 69.36),
    (2022.0, 69.29),
    (2024.0, 69.18),
    (2026.0, 69.15),
];

/// Embedded IERS leap-second steps `(mjd_utc, tai_minus_utc_sec)` since 1972-01-01.
const EMBEDDED_LEAP_SECONDS: &[(f64, f64)] = &[
    (41317.0, 10.0), // 1972-01-01
    (41499.0, 11.0), // 1972-07-01
    (41683.0, 12.0), // 1973-01-01
    (42048.0, 13.0), // 1974-01-01
    (42413.0, 14.0), // 1975-01-01
    (42778.0, 15.0), // 1976-01-01
    (43144.0, 16.0), // 1977-01-01
    (43509.0, 17.0), // 1978-01-01
    (43874.0, 18.0), // 1979-01-01
    (44239.0, 19.0), // 1980-01-01
    (44786.0, 20.0), // 1981-07-01
    (45151.0, 21.0), // 1982-07-01
    (45516.0, 22.0), // 1983-07-01
    (46247.0, 23.0), // 1985-07-01
    (47161.0, 24.0), // 1988-01-01
    (47892.0, 25.0), // 1990-01-01
    (48257.0, 26.0), // 1991-01-01
    (48804.0, 27.0), // 1992-07-01
    (49169.0, 28.0), // 1993-07-01
    (49534.0, 29.0), // 1994-07-01
    (50083.0, 30.0), // 1996-01-01
    (50630.0, 31.0), // 1997-07-01
    (51179.0, 32.0), // 1999-01-01
    (53736.0, 33.0), // 2006-01-01
    (54832.0, 34.0), // 2009-01-01
    (56109.0, 35.0), // 2012-07-01
    (57204.0, 36.0), // 2015-07-01
    (57754.0, 37.0), // 2017-01-01
];

/// Immutable, thread-safe (`Send + Sync`) time conversion table supporting embedded baseline
/// and atomic IERS updates (`ARCHITECTURE.md` §2.1).
#[derive(Debug, Clone)]
pub struct TimeTable {
    version: String,
    /// Daily `(mjd, ut1_minus_utc_sec, is_predicted)` rows loaded from `finals2000A` if available.
    finals_mjd_dut1: Vec<(f64, f64, bool)>,
    /// Leap second steps `(mjd_utc, tai_minus_utc_sec)`.
    leap_seconds: Vec<(f64, f64)>,
}

impl Default for TimeTable {
    fn default() -> Self {
        Self::embedded()
    }
}

/// Identifier for the built-in baseline $\Delta T$ and leap second table.
pub const EMBEDDED_TIME_TABLE_VERSION: &str = "embedded-smh2016-iers2026";

impl TimeTable {
    /// Returns the built-in baseline time table (`"embedded-smh2016-iers2026"`).
    #[must_use]
    pub fn embedded() -> Self {
        Self {
            version: EMBEDDED_TIME_TABLE_VERSION.to_owned(),
            finals_mjd_dut1: Vec::new(),
            leap_seconds: EMBEDDED_LEAP_SECONDS.to_vec(),
        }
    }

    /// Loads an updated [`TimeTable`] from `data_dir` if verified `time_sync_manifest.json`,
    /// `iers_finals2000A.txt`, and `iers_leapseconds.txt` are present; otherwise returns `embedded()`.
    pub fn load_or_embedded(
        data_dir: &Path,
        pin_version: Option<&str>,
    ) -> Result<Self, EphemerisError> {
        if let Some(pinned) = pin_version {
            if pinned == "embedded" || pinned == "embedded-smh2016-iers2026" {
                return Ok(Self::embedded());
            }
        }

        let sync_manifest_path = data_dir.join("time_sync_manifest.json");
        let finals_path = data_dir.join("iers_finals2000A.txt");
        let leap_path = data_dir.join("iers_leapseconds.txt");

        if !sync_manifest_path.exists() || !finals_path.exists() || !leap_path.exists() {
            if let Some(pinned) = pin_version {
                return Err(EphemerisError::MissingData {
                    detail: format!(
                        "Pinned time_table_version '{pinned}' requested, but synced IERS files not found in {}",
                        data_dir.display()
                    ),
                });
            }
            return Ok(Self::embedded());
        }

        let manifest_text =
            fs::read_to_string(&sync_manifest_path).map_err(|e| EphemerisError::MissingData {
                detail: format!("Failed to read {}: {e}", sync_manifest_path.display()),
            })?;
        let sync_manifest: TimeSyncManifest =
            serde_json::from_str(&manifest_text).map_err(|e| EphemerisError::CorruptKernel {
                path: sync_manifest_path.display().to_string(),
                reason: format!("Invalid time_sync_manifest.json: {e}"),
            })?;

        if let Some(pinned) = pin_version {
            if pinned != sync_manifest.iers_release_date {
                return Err(EphemerisError::InvalidInput {
                    message: format!(
                        "Pinned time_table_version '{pinned}' does not match synced IERS version '{}'",
                        sync_manifest.iers_release_date
                    ),
                });
            }
        }

        let finals_bytes = fs::read(&finals_path).map_err(|e| EphemerisError::MissingData {
            detail: format!("Failed to read {}: {e}", finals_path.display()),
        })?;
        if sha256_hex(&finals_bytes) != sync_manifest.finals2000a_sha256 {
            return Err(EphemerisError::CorruptKernel {
                path: finals_path.display().to_string(),
                reason: "SHA-256 mismatch on iers_finals2000A.txt".to_owned(),
            });
        }

        let leap_bytes = fs::read(&leap_path).map_err(|e| EphemerisError::MissingData {
            detail: format!("Failed to read {}: {e}", leap_path.display()),
        })?;
        if sha256_hex(&leap_bytes) != sync_manifest.leapseconds_sha256 {
            return Err(EphemerisError::CorruptKernel {
                path: leap_path.display().to_string(),
                reason: "SHA-256 mismatch on iers_leapseconds.txt".to_owned(),
            });
        }

        let finals_text = String::from_utf8_lossy(&finals_bytes);
        let mut finals_mjd_dut1 = Vec::new();
        for line in finals_text.lines() {
            // IERS finals2000A fixed-width columns:
            // Cols 8..15 (0-indexed 7..15): MJD
            // Col 58 (0-indexed 57..58): IERS flag ('I' = observed, 'P' = predicted)
            // Cols 59..68 (0-indexed 58..68): UT1-UTC (seconds)
            if let (Some(mjd_s), Some(dut1_s)) = (line.get(7..15), line.get(58..68)) {
                if let (Ok(mjd), Ok(dut1)) =
                    (mjd_s.trim().parse::<f64>(), dut1_s.trim().parse::<f64>())
                {
                    let is_predicted = line.get(57..58) == Some("P");
                    finals_mjd_dut1.push((mjd, dut1, is_predicted));
                }
            }
        }

        let leap_text = String::from_utf8_lossy(&leap_bytes);
        let parsed_leaps = Self::parse_leapseconds_table(&leap_text);
        let leap_seconds = if parsed_leaps.is_empty() {
            EMBEDDED_LEAP_SECONDS.to_vec()
        } else {
            parsed_leaps
        };

        Ok(Self {
            version: sync_manifest.iers_release_date,
            finals_mjd_dut1,
            leap_seconds,
        })
    }

    /// Parses an IERS `Leap_Second.dat` / `iers_leapseconds.txt` table into `(mjd, tai_minus_utc_sec)`
    /// step pairs from `1972-01-01` (`MJD >= 41317.0`) onward.
    #[must_use]
    pub fn parse_leapseconds_table(text: &str) -> Vec<(f64, f64)> {
        let mut out = Vec::new();
        for raw_line in text.lines() {
            let line = raw_line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if let (Some(jd_pos), Some(tai_pos)) = (line.find("=JD"), line.find("TAI-UTC=")) {
                if tai_pos > jd_pos + 3 {
                    let jd_part = &line[(jd_pos + 3)..tai_pos];
                    let after_tai = &line[(tai_pos + 8)..];
                    let tai_token = after_tai.split_whitespace().next().unwrap_or("");
                    if let (Ok(jd), Ok(tai_utc)) =
                        (jd_part.trim().parse::<f64>(), tai_token.parse::<f64>())
                    {
                        if jd.is_finite() && tai_utc.is_finite() {
                            let mjd = jd - 2_400_000.5;
                            if mjd >= 41_317.0 {
                                out.push((mjd, tai_utc));
                            }
                        }
                    }
                }
            } else {
                // Also support simple 2-column "<mjd> <tai_utc>" lines
                let mut parts = line.split_whitespace();
                if let (Some(m_s), Some(t_s)) = (parts.next(), parts.next()) {
                    if let (Ok(mjd), Ok(tai_utc)) = (m_s.parse::<f64>(), t_s.parse::<f64>()) {
                        if mjd.is_finite() && tai_utc.is_finite() && mjd >= 41_317.0 {
                            out.push((mjd, tai_utc));
                        }
                    }
                }
            }
        }
        out
    }

    /// Returns the active time table version identifier.
    #[must_use]
    pub fn version(&self) -> &str {
        &self.version
    }

    /// Returns `TAI - UTC` in seconds at the given MJD (`JD - 2400000.5`).
    #[must_use]
    pub fn tai_minus_utc_sec(&self, mjd: f64) -> f64 {
        let idx = self.leap_seconds.partition_point(|&(step_mjd, _)| step_mjd <= mjd);
        if idx == 0 {
            10.0
        } else {
            self.leap_seconds
                .get(idx - 1)
                .map_or(37.0, |&(_, ls)| ls)
        }
    }

    /// Evaluates $\Delta T = \text{TT} - \text{UT1}$ in seconds and extrapolation status at a given Julian Date.
    #[must_use]
    pub fn delta_t(&self, jd: f64) -> DeltaTInfo {
        let mjd = jd - 2_400_000.5;
        let tt_minus_tdb_sec = -tdb_minus_tt_sec(Tt::from_jd(jd), None);

        // 1. Check high-resolution IERS finals2000A daily table if loaded and within coverage.
        if let (Some(&(first_mjd, _, _)), Some(&(last_mjd, _, _))) =
            (self.finals_mjd_dut1.first(), self.finals_mjd_dut1.last())
        {
            if mjd >= first_mjd && mjd <= last_mjd {
                let idx = self
                    .finals_mjd_dut1
                    .partition_point(|&(m, _, _)| m < mjd)
                    .min(self.finals_mjd_dut1.len() - 1);
                let (dut1, is_predicted) = if idx == 0 {
                    self.finals_mjd_dut1
                        .first()
                        .map_or((0.0, false), |&(_, d, p)| (d, p))
                } else {
                    let (m0, d0, p0) = self
                        .finals_mjd_dut1
                        .get(idx - 1)
                        .copied()
                        .unwrap_or((mjd, 0.0, false));
                    let (m1, d1, p1) = self
                        .finals_mjd_dut1
                        .get(idx)
                        .copied()
                        .unwrap_or((mjd, 0.0, false));
                    let frac = if (m1 - m0).abs() > 1e-9 {
                        ((mjd - m0) / (m1 - m0)).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    (d0 + frac * (d1 - d0), p0 || p1)
                };
                let tai_utc = self.tai_minus_utc_sec(mjd);
                let delta_t_sec = 32.184 + tai_utc - dut1;
                return DeltaTInfo {
                    delta_t_sec,
                    tt_minus_tdb_sec,
                    time_extrapolated: is_predicted,
                    time_table_version: self.version.clone(),
                };
            }
        }

        // 2. Check embedded Stephenson–Morrison–Hohenkerk (2016) + USNO/IERS knots [-700 .. 2026].
        let dec_year = 2000.0 + (jd - JD_J2000) / 365.2425;
        let first_knot = EMBEDDED_DELTAT_KNOTS.first().copied().unwrap_or((-700.0, 20270.0));
        let last_knot = EMBEDDED_DELTAT_KNOTS.last().copied().unwrap_or((2026.0, 69.15));

        if dec_year >= first_knot.0 && dec_year <= last_knot.0 {
            let idx = EMBEDDED_DELTAT_KNOTS
                .partition_point(|&(y, _)| y < dec_year)
                .min(EMBEDDED_DELTAT_KNOTS.len() - 1);
            let dt = if idx == 0 {
                first_knot.1
            } else {
                let (y0, dt0) = EMBEDDED_DELTAT_KNOTS
                    .get(idx - 1)
                    .copied()
                    .unwrap_or(first_knot);
                let (y1, dt1) = EMBEDDED_DELTAT_KNOTS
                    .get(idx)
                    .copied()
                    .unwrap_or(last_knot);
                let frac = if (y1 - y0).abs() > 1e-9 {
                    ((dec_year - y0) / (y1 - y0)).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                dt0 + frac * (dt1 - dt0)
            };
            return DeltaTInfo {
                delta_t_sec: dt,
                tt_minus_tdb_sec,
                time_extrapolated: false,
                time_table_version: self.version.clone(),
            };
        }

        // 3. Stephenson, Morrison & Hohenkerk (2016) parabolic extrapolation outside empirical coverage:
        //    Delta T(t) = -320.0 + 32.5 * ((t - 1825.0) / 100.0)^2
        let u = (dec_year - 1825.0) / 100.0;
        let parabola = -320.0 + 32.5 * u * u;
        let delta_t_sec = if dec_year > last_knot.0 {
            // Smooth C0-continuous offset at last_knot.0
            let u_last = (last_knot.0 - 1825.0) / 100.0;
            let parabola_last = -320.0 + 32.5 * u_last * u_last;
            last_knot.1 + (parabola - parabola_last)
        } else {
            // Smooth C0-continuous offset at first_knot.0
            let u_first = (first_knot.0 - 1825.0) / 100.0;
            let parabola_first = -320.0 + 32.5 * u_first * u_first;
            first_knot.1 + (parabola - parabola_first)
        };

        DeltaTInfo {
            delta_t_sec,
            tt_minus_tdb_sec,
            time_extrapolated: true,
            time_table_version: self.version.clone(),
        }
    }

    /// Evaluates [`SiderealTimeInfo`] (`GMST`, `GAST`, local sidereal time, and `ARMC`) at `epoch`
    /// for a given `nutation` mode and geographic `longitude_deg` (East positive).
    #[must_use]
    pub fn sidtime(
        &self,
        epoch: TimeInput,
        nutation: crate::types::NutationMode,
        longitude_deg: f64,
    ) -> SiderealTimeInfo {
        let (ut1, tt, _, _) = self.resolve_scales(epoch, None);
        let gmst_hours = (gmst_rad(ut1, tt).to_degrees() / 15.0).rem_euclid(24.0);
        let gast_hours = (gast_rad(ut1, tt).to_degrees() / 15.0).rem_euclid(24.0);
        let base_hours = match nutation {
            crate::types::NutationMode::True => gast_hours,
            crate::types::NutationMode::Mean => gmst_hours,
        };
        let local_sidereal_hours = (base_hours + longitude_deg / 15.0).rem_euclid(24.0);
        let armc_deg = (local_sidereal_hours * 15.0).rem_euclid(360.0);
        SiderealTimeInfo {
            gmst_hours,
            gast_hours,
            local_sidereal_hours,
            armc_deg,
        }
    }

    /// Converts any [`TimeInput`] into `(Ut1, Tt, Tdb, DeltaTInfo)` for a given optional topocentric observer.
    #[must_use]
    pub fn resolve_scales(
        &self,
        input: TimeInput,
        observer: Option<Observer>,
    ) -> (Ut1, Tt, Tdb, DeltaTInfo) {
        match input {
            TimeInput::Tt(tt_split) => {
                let tt = Tt(tt_split);
                let dt_info = self.delta_t(tt_split.to_f64());
                let ut1 = Ut1(tt_split.add_days(-dt_info.delta_t_sec / SECONDS_PER_DAY));
                let tdb = tt_to_tdb(tt, observer);
                (ut1, tt, tdb, dt_info)
            }
            TimeInput::Ut1(ut1_split) => {
                let ut1 = Ut1(ut1_split);
                let dt_info = self.delta_t(ut1_split.to_f64());
                let tt = Tt(ut1_split.add_days(dt_info.delta_t_sec / SECONDS_PER_DAY));
                let tdb = tt_to_tdb(tt, observer);
                (ut1, tt, tdb, dt_info)
            }
            TimeInput::Utc(utc_split) => {
                let mjd = utc_split.to_f64() - 2_400_000.5;
                if mjd < 41_317.0 {
                    // Before 1972-01-01 (MJD 41317.0), atomic stepped UTC did not exist;
                    // civil/historical "UTC" timestamps represent UT1 (GMT/UT), so UTC == UT1.
                    let ut1 = Ut1(utc_split);
                    let dt_info = self.delta_t(utc_split.to_f64());
                    let tt = Tt(utc_split.add_days(dt_info.delta_t_sec / SECONDS_PER_DAY));
                    let tdb = tt_to_tdb(tt, observer);
                    (ut1, tt, tdb, dt_info)
                } else {
                    let tai_utc = self.tai_minus_utc_sec(mjd);
                    let tt_split = utc_split.add_days((tai_utc + 32.184) / SECONDS_PER_DAY);
                    let tt = Tt(tt_split);
                    let dt_info = self.delta_t(tt_split.to_f64());
                    let ut1 = Ut1(tt_split.add_days(-dt_info.delta_t_sec / SECONDS_PER_DAY));
                    let tdb = tt_to_tdb(tt, observer);
                    (ut1, tt, tdb, dt_info)
                }
            }
            TimeInput::Tdb(tdb_split) => {
                let tdb = Tdb(tdb_split);
                let tt = tdb_to_tt(tdb, observer);
                let dt_info = self.delta_t(tt.0.to_f64());
                let ut1 = Ut1(tt.0.add_days(-dt_info.delta_t_sec / SECONDS_PER_DAY));
                (ut1, tt, tdb, dt_info)
            }
        }
    }
}

/// Computes $\text{TDB} - \text{TT}$ in seconds via IAU SOFA `sofars::ts::dtdb` (Fairhead & Bretagnon 1990).
#[must_use]
pub fn tdb_minus_tt_sec(tt: Tt, observer: Option<Observer>) -> f64 {
    let (elong, u_km, v_km) = if let Some(obs) = observer {
        let lon_rad = obs.longitude_deg.to_radians();
        let lat_rad = obs.latitude_deg.to_radians();
        if let Ok(xyz_m) = gd2gc(1, lon_rad, lat_rad, obs.altitude_m) {
            let u = (xyz_m[0] * xyz_m[0] + xyz_m[1] * xyz_m[1]).sqrt() / 1000.0;
            let v = xyz_m[2] / 1000.0;
            (lon_rad, u, v)
        } else {
            (0.0, 0.0, 0.0)
        }
    } else {
        (0.0, 0.0, 0.0)
    };
    let ut_frac = (tt.0.jd2 + 0.5).rem_euclid(1.0);
    dtdb(tt.0.jd1, tt.0.jd2, ut_frac, elong, u_km, v_km)
}

/// Converts [`Tt`] to [`Tdb`] using `sofars::ts::dtdb`.
#[must_use]
pub fn tt_to_tdb(tt: Tt, observer: Option<Observer>) -> Tdb {
    let dt_sec = tdb_minus_tt_sec(tt, observer);
    Tdb(tt.0.add_days(dt_sec / SECONDS_PER_DAY))
}

/// Converts [`Tdb`] to [`Tt`] by inverting `sofars::ts::dtdb`.
#[must_use]
pub fn tdb_to_tt(tdb: Tdb, observer: Option<Observer>) -> Tt {
    let approx_tt = Tt(tdb.0);
    let dt_sec = tdb_minus_tt_sec(approx_tt, observer);
    let refined_tt = Tt(tdb.0.add_days(-dt_sec / SECONDS_PER_DAY));
    let dt_sec_2 = tdb_minus_tt_sec(refined_tt, observer);
    Tt(tdb.0.add_days(-dt_sec_2 / SECONDS_PER_DAY))
}

/// Computes Earth Rotation Angle (ERA, IAU 2000) in radians `[0, 2π)` strictly from [`Ut1`].
#[must_use]
pub fn earth_rotation_angle_rad(ut1: Ut1) -> f64 {
    era00(ut1.0.jd1, ut1.0.jd2)
}

/// Computes Greenwich Mean Sidereal Time (GMST, IAU 2006) in radians `[0, 2π)` from [`Ut1`] and [`Tt`].
#[must_use]
pub fn gmst_rad(ut1: Ut1, tt: Tt) -> f64 {
    gmst06(ut1.0.jd1, ut1.0.jd2, tt.0.jd1, tt.0.jd2)
}

/// Computes Greenwich Apparent Sidereal Time (GAST, IAU 2006/2000A) in radians `[0, 2π)` strictly from [`Ut1`] and [`Tt`].
#[must_use]
pub fn gast_rad(ut1: Ut1, tt: Tt) -> f64 {
    let ft = crate::frame::FrameTransform::at_tt(tt);
    gst06(ut1.0.jd1, ut1.0.jd2, tt.0.jd1, tt.0.jd2, &ft.r_bpn)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_julday_revjul_roundtrip_and_papal_reform() {
        // J2000.0: 2000-01-01 12:00
        let j2000 = julday(2000, 1, 1, 12.0, CalendarKind::Gregorian).expect("julday J2000");
        assert!((j2000 - JD_J2000).abs() < 1e-12);
        let dt = revjul(j2000, CalendarKind::Gregorian).expect("revjul J2000");
        assert_eq!((dt.year, dt.month, dt.day), (2000, 1, 1));
        assert!((dt.hour - 12.0).abs() < 1e-10);

        // Papal reform boundary: 1582-10-04 Julian == JD 2299159.5, 1582-10-15 Gregorian == JD 2299160.5
        let jd_jul = julday(1582, 10, 4, 0.0, CalendarKind::AstroAuto).expect("1582-10-04");
        let jd_greg = julday(1582, 10, 15, 0.0, CalendarKind::AstroAuto).expect("1582-10-15");
        assert!((jd_jul - 2_299_159.5).abs() < 1e-10);
        assert!((jd_greg - 2_299_160.5).abs() < 1e-10);
        assert!((jd_greg - jd_jul - 1.0).abs() < 1e-10);

        // Negative year (astronomical year 0 = 1 BCE, -4712 = 4713 BCE)
        let jd_zero = julday(-4712, 1, 1, 12.0, CalendarKind::Julian).expect("-4712-01-01");
        assert!(jd_zero.abs() < 1e-10);
        let rev_zero = revjul(jd_zero, CalendarKind::Julian).expect("revjul 0");
        assert_eq!((rev_zero.year, rev_zero.month, rev_zero.day), (-4712, 1, 1));
    }

    #[test]
    fn test_delta_t_embedded_and_extrapolated() {
        let table = TimeTable::embedded();
        assert_eq!(table.version(), EMBEDDED_TIME_TABLE_VERSION);

        // J2000.0 (2000.0): Delta T should be ~63.83 s and not extrapolated
        let dt_j2000 = table.delta_t(JD_J2000);
        assert!(!dt_j2000.time_extrapolated);
        assert!((dt_j2000.delta_t_sec - 63.83).abs() < 0.1);

        // Far future (2100.0): extrapolated = true
        let dt_2100 = table.delta_t(JD_J2000 + 36525.0);
        assert!(dt_2100.time_extrapolated);
        assert!(dt_2100.delta_t_sec > 80.0);

        // Far past (-1000 CE): extrapolated = true
        let jd_ancient = julday(-1000, 1, 1, 12.0, CalendarKind::Julian).expect("-1000");
        let dt_ancient = table.delta_t(jd_ancient);
        assert!(dt_ancient.time_extrapolated);
        assert!(dt_ancient.delta_t_sec > 20_000.0);

        // Verify strict monotonicity of EMBEDDED_DELTAT_KNOTS from -700 to 1700 CE (SMH2016)
        for w in EMBEDDED_DELTAT_KNOTS.windows(2) {
            if let [(y0, dt0), (y1, dt1)] = *w {
                if y1 <= 1700.0 {
                    assert!(
                        dt0 > dt1,
                        "Non-monotonic Delta T between {y0} ({dt0}) and {y1} ({dt1})"
                    );
                }
            }
        }

        // Pre-1972 UTC must resolve identically to UT1 (UTC == UT1 before MJD 41317)
        let jd_1900 = julday(1900, 1, 1, 12.0, CalendarKind::Gregorian).expect("1900");
        let (ut1_from_utc, tt_from_utc, _, _) =
            table.resolve_scales(TimeInput::from_utc_jd(jd_1900), None);
        let (ut1_from_ut1, tt_from_ut1, _, _) =
            table.resolve_scales(TimeInput::from_ut1_jd(jd_1900), None);
        assert!((ut1_from_utc.jd_f64() - ut1_from_ut1.jd_f64()).abs() < 1e-14);
        assert!((tt_from_utc.jd_f64() - tt_from_ut1.jd_f64()).abs() < 1e-14);
    }

    #[test]
    fn test_tt_tdb_roundtrip_and_sidereal_time() {
        let tt = Tt::from_jd(JD_J2000);
        let tdb = tt_to_tdb(tt, None);
        // At J2000, |TDB - TT| < 2 ms
        assert!((tdb.jd_f64() - tt.jd_f64()).abs() * SECONDS_PER_DAY < 0.002);
        let tt_back = tdb_to_tt(tdb, None);
        assert!((tt_back.jd_f64() - tt.jd_f64()).abs() * SECONDS_PER_DAY < 1e-9);

        let ut1 = Ut1::from_jd(JD_J2000);
        let era = earth_rotation_angle_rad(ut1);
        let gmst = gmst_rad(ut1, tt);
        let gast = gast_rad(ut1, tt);
        assert!((0.0..std::f64::consts::TAU).contains(&era));
        assert!((0.0..std::f64::consts::TAU).contains(&gmst));
        assert!((0.0..std::f64::consts::TAU).contains(&gast));
        // Equation of the equinoxes |GAST - GMST| is less than ~1.2 arcseconds (~20 arcsec on sky)
        assert!((gast - gmst).abs() < 1e-3);
    }
}
