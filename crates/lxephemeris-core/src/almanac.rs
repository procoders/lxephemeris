//! Vedic Jyotish Panchanga, Divisional Charts (`Varga` `D1..D60`), Vimshottari Dasha,
//! Human Design 64-Gate Rave Mandala (`Gate`, `Line`, `Color`, `Tone`, `Base`),
//! and Traditional Western Almanac (`ChartSectInfo`, `PlanetaryHoursReport`, and the
//! 7 Classical Hermetic Lots of Paulus Alexandrinus).
//!
//! Clean-room mathematical implementation adhering strictly to `#![forbid(unsafe_code)]`.

use crate::astrometry::{AtmosphereParams, RiseTransEvent, RiseTransOptions};
use crate::engine::Engine;
use crate::error::EphemerisError;
use crate::events::{search_body_longitude, search_lunar_phase, wrap_pm180};
use crate::houses::{house_of_lon, HouseOptions};
use crate::time::TimeInput;
use crate::types::{
    AbMode, AyanamsaModel, BodyId, CalcOptions, Equinox, NutationMode, Observer, Origin, Plane,
    PlanetId, Zodiac,
};
use serde::{Deserialize, Serialize};

/// Type alias for [`Observer`] used in almanac and Panchanga APIs.
pub type ObserverLocation = Observer;

/// Canonical names of the 12 Tropical / Sidereal zodiac signs (`0 = Aries .. 11 = Pisces`).
pub const ZODIAC_SIGN_NAMES: [&str; 12] = [
    "Aries",
    "Taurus",
    "Gemini",
    "Cancer",
    "Leo",
    "Virgo",
    "Libra",
    "Scorpio",
    "Sagittarius",
    "Capricorn",
    "Aquarius",
    "Pisces",
];

/// Traditional Parashara / Classical planetary rulers of the 12 zodiac signs (`0 = Aries .. 11 = Pisces`).
pub const ZODIAC_SIGN_LORDS: [&str; 12] = [
    "Mars",    // 0: Aries
    "Venus",   // 1: Taurus
    "Mercury", // 2: Gemini
    "Moon",    // 3: Cancer
    "Sun",     // 4: Leo
    "Mercury", // 5: Virgo
    "Venus",   // 6: Libra
    "Mars",    // 7: Scorpio
    "Jupiter", // 8: Sagittarius
    "Saturn",  // 9: Capricorn
    "Saturn",  // 10: Aquarius
    "Jupiter", // 11: Pisces
];

/// Canonical names of the 27 Vedic Nakshatras (`1..=27`, indexed `0..27`, each spanning $13^\circ 20' = 40^\circ/3$).
pub const NAKSHATRA_NAMES_27: [&str; 27] = [
    "Ashwini",
    "Bharani",
    "Krittika",
    "Rohini",
    "Mrigashira",
    "Ardra",
    "Punarvasu",
    "Pushya",
    "Ashlesha",
    "Magha",
    "Purva Phalguni",
    "Uttara Phalguni",
    "Hasta",
    "Chitra",
    "Swati",
    "Vishakha",
    "Anuradha",
    "Jyeshtha",
    "Mula",
    "Purva Ashadha",
    "Uttara Ashadha",
    "Shravana",
    "Dhanishta",
    "Shatabhisha",
    "Purva Bhadrapada",
    "Uttara Bhadrapada",
    "Revati",
];

/// Canonical 9 Vimshottari Dasha lords in sequence (`(nakshatra_index - 1) % 9`).
pub const VIMSHOTTARI_LORDS_9: [&str; 9] = [
    "Ketu",
    "Venus",
    "Sun",
    "Moon",
    "Mars",
    "Rahu",
    "Jupiter",
    "Saturn",
    "Mercury",
];

/// Canonical full Mahadasha durations in solar years for [`VIMSHOTTARI_LORDS_9`] (sum = `120.0` years).
pub const VIMSHOTTARI_YEARS_9: [f64; 9] = [7.0, 20.0, 6.0, 10.0, 7.0, 18.0, 16.0, 19.0, 17.0];

/// Default Gregorian mean solar year length in days (`365.2425`) for Vimshottari Dasha calculations.
pub const DEFAULT_DASHAS_YEAR_DAYS: f64 = 365.2425;

/// Exact tropical longitude (in degrees) of the start of **Gate 41, Line 1** on the Human Design Rave Mandala
/// ($302^\circ 00' 00'' = 02^\circ 00'$ Aquarius).
pub const HD_MANDALA_START_LON_DEG: f64 = 302.0;

/// Canonical counter-clockwise sequence of the 64 Human Design / I-Ching Rave Mandala gates
/// starting at $302^\circ 00' 00''$ tropical longitude (`Gate 41, Line 1`), each spanning $5.625^\circ$ ($5^\circ 37' 30''$).
pub const HD_MANDALA_GATES_64: [u8; 64] = [
    41, 19, 13, 49, 30, 55, 37, 63, 22, 36, 25, 17, 21, 51, 42, 3, 27, 24, 2, 23, 8, 20, 16, 35,
    45, 12, 15, 52, 39, 53, 62, 56, 31, 33, 7, 4, 29, 59, 40, 64, 47, 6, 46, 18, 48, 57, 32, 50,
    28, 44, 1, 43, 14, 34, 9, 5, 26, 11, 10, 58, 38, 54, 61, 60,
];

/// Canonical names of the 7 Vedic weekdays (*Vara*, `0 = Sunday / Ravivara .. 6 = Saturday / Shanivara`).
pub const VARA_NAMES_7: [&str; 7] = [
    "Ravivara",
    "Somavara",
    "Mangalavara",
    "Budhavara",
    "Guruvara",
    "Shukravara",
    "Shanivara",
];

/// Planetary lords of the 7 weekdays (`0 = Sunday .. 6 = Saturday`).
pub const VARA_LORDS_7: [&str; 7] = [
    "Sun",
    "Moon",
    "Mars",
    "Mercury",
    "Jupiter",
    "Venus",
    "Saturn",
];

/// Canonical names of the 30 lunar Tithis (`1..=30`, indexed `0..30`, each $12^\circ$ of lunar-solar elongation).
pub const TITHI_NAMES_30: [&str; 30] = [
    "Pratipada",
    "Dwitiya",
    "Tritiya",
    "Chaturthi",
    "Panchami",
    "Shashthi",
    "Saptami",
    "Ashtami",
    "Navami",
    "Dashami",
    "Ekadashi",
    "Dwadashi",
    "Trayodashi",
    "Chaturdashi",
    "Purnima",
    "Pratipada",
    "Dwitiya",
    "Tritiya",
    "Chaturthi",
    "Panchami",
    "Shashthi",
    "Saptami",
    "Ashtami",
    "Navami",
    "Dashami",
    "Ekadashi",
    "Dwadashi",
    "Trayodashi",
    "Chaturdashi",
    "Amavasya",
];

/// Canonical names of the 27 Vedic Yogas (`1..=27`, indexed `0..27`, each $13^\circ 20'$ of $(\lambda_{\mathbb{C},\text{sid}} + \lambda_{\odot,\text{sid}}) \bmod 360^\circ$).
pub const YOGA_NAMES_27: [&str; 27] = [
    "Vishkumbha",
    "Priti",
    "Ayushman",
    "Saubhagya",
    "Shobhana",
    "Atiganda",
    "Sukarma",
    "Dhriti",
    "Shula",
    "Ganda",
    "Vriddhi",
    "Dhruva",
    "Vyaghata",
    "Harshana",
    "Vajra",
    "Siddhi",
    "Vyatipata",
    "Variyan",
    "Parigha",
    "Shiva",
    "Siddha",
    "Sadhya",
    "Shubha",
    "Shukla",
    "Brahma",
    "Indra",
    "Vaidhriti",
];

/// Canonical names of the 7 movable (*Chara*) Karanas repeating across half-tithis `2..=57`.
pub const CHARA_KARANA_NAMES_7: [&str; 7] = [
    "Bava",
    "Balava",
    "Kaulava",
    "Taitila",
    "Gara",
    "Vanija",
    "Vishti",
];

/// Classical Chaldean descending order of the 7 planetary rulers used for planetary hours.
pub const CHALDEAN_ORDER_7: [&str; 7] = [
    "Saturn",
    "Jupiter",
    "Mars",
    "Sun",
    "Venus",
    "Mercury",
    "Moon",
];

/// Normalizes an angle in degrees into `[0.0, 360.0)`.
#[inline]
#[must_use]
pub fn normalize_lon_360(lon_deg: f64) -> f64 {
    if !lon_deg.is_finite() {
        return 0.0;
    }
    let r = lon_deg.rem_euclid(360.0);
    if r >= 360.0 {
        0.0
    } else {
        r
    }
}

/// Floors a non-negative floating-point value with a `1e-12` boundary guard against IEEE-754 rounding just below an integer.
#[inline]
fn floor_with_eps(x: f64) -> usize {
    let r = x.round();
    if r >= 0.0 && (r - x).abs() < 1.0e-12 {
        r as usize
    } else {
        x.floor().max(0.0) as usize
    }
}

// ============================================================================
// 1. ZodiacPlacement, 27 Nakshatras, 64 Human Design Gates, and Vargas (D1..D60)
// ============================================================================

/// Complete zodiacal, Nakshatra (`1..=27`, `Pada 1..=4`), and Human Design Rave Mandala
/// (`Gate 1..=64`, `Line 1..=6`, `Color 1..=6`, `Tone 1..=6`, `Base 1..=5`) placement for an ecliptic longitude.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ZodiacPlacement {
    /// Normalized ecliptic longitude in degrees `[0.0, 360.0)`.
    pub lon_deg: f64,
    /// 0-based zodiac sign index (`0 = Aries .. 11 = Pisces`).
    pub sign_index: u8,
    /// Canonical zodiac sign name (`"Aries" .. "Pisces"`).
    pub sign_name: String,
    /// Degrees within the sign `[0.0, 30.0)`.
    pub deg_in_sign: f64,
    /// Whole degrees within the sign (`0..=29`).
    pub dms_deg: u8,
    /// Arcminutes within the degree (`0..=59`).
    pub dms_min: u8,
    /// Arcseconds within the arcminute `[0.0, 60.0)`.
    pub dms_sec: f64,
    /// Decan (face) within the sign (`1..=3`, each $10^\circ$).
    pub decan: u8,
    /// 1-based Vedic Nakshatra index (`1..=27`, each $13^\circ 20' = 40^\circ/3$).
    pub nakshatra_index: u8,
    /// Canonical Vedic Nakshatra name (`"Ashwini" .. "Revati"`).
    pub nakshatra_name: String,
    /// Vimshottari planetary lord of the Nakshatra (`"Ketu"`, `"Venus"`, `"Sun"`, `"Moon"`, `"Mars"`, `"Rahu"`, `"Jupiter"`, `"Saturn"`, `"Mercury"`).
    pub nakshatra_lord: String,
    /// Quarter (*Pada*) within the Nakshatra (`1..=4`, each $3^\circ 20' = 10^\circ/3$).
    pub pada: u8,
    /// Human Design / I-Ching Rave Mandala gate (`1..=64`, each $5.625^\circ$ from $302^\circ 00'$ tropical).
    pub hd_gate: u8,
    /// Human Design line within the gate (`1..=6`, each $0.9375^\circ = 56'15''$).
    pub hd_line: u8,
    /// Human Design color within the line (`1..=6`, each $0.15625^\circ = 9'22.5''$).
    pub hd_color: u8,
    /// Human Design tone within the color (`1..=6`, each $0.026041\overline{6}^\circ = 1'33.75''$).
    pub hd_tone: u8,
    /// Human Design base within the tone (`1..=5`, each $0.005208\overline{3}^\circ = 18.75''$).
    pub hd_base: u8,
}

impl ZodiacPlacement {
    /// Maps any ecliptic longitude `lon_deg` into its [`ZodiacPlacement`].
    #[must_use]
    pub fn from_longitude(lon_deg: f64) -> Self {
        map_zodiac_placement(lon_deg)
    }
}

/// Maps an ecliptic longitude `lon_deg` into a [`ZodiacPlacement`] (Sign, DMS, Decan, Nakshatra/Pada, and HD Gate/Line/Color/Tone/Base).
#[must_use]
pub fn map_zodiac_placement(lon_deg: f64) -> ZodiacPlacement {
    let lon = normalize_lon_360(lon_deg);

    // Sign & DMS
    let sign_idx = floor_with_eps(lon / 30.0).min(11);
    let sign_index = sign_idx as u8;
    let sign_name = ZODIAC_SIGN_NAMES
        .get(sign_idx)
        .copied()
        .unwrap_or("Aries")
        .to_owned();
    let deg_in_sign = (lon - (sign_idx as f64) * 30.0).clamp(0.0, 30.0 - 1.0e-12);

    let dms_deg = (deg_in_sign.floor() as u8).min(29);
    let rem_min = (deg_in_sign - f64::from(dms_deg)).max(0.0) * 60.0;
    let dms_min = (rem_min.floor() as u8).min(59);
    let dms_sec = ((rem_min - f64::from(dms_min)).max(0.0) * 60.0).clamp(0.0, 60.0 - 1.0e-9);
    let decan = (floor_with_eps(deg_in_sign / 10.0).min(2) as u8) + 1;

    // 27 Nakshatras & 4 Padas (108 padas of 10/3° across 360°)
    let total_pada_idx = floor_with_eps((lon * 108.0) / 360.0).min(107);
    let nak_0 = (total_pada_idx / 4).min(26);
    let pada = ((total_pada_idx % 4) as u8) + 1;
    let nakshatra_index = (nak_0 as u8) + 1;
    let nakshatra_name = NAKSHATRA_NAMES_27
        .get(nak_0)
        .copied()
        .unwrap_or("Ashwini")
        .to_owned();
    let nakshatra_lord = VIMSHOTTARI_LORDS_9
        .get(nak_0 % 9)
        .copied()
        .unwrap_or("Ketu")
        .to_owned();

    // 64 Human Design Gates, 6 Lines, 6 Colors, 6 Tones, 5 Bases
    // 64 * 6 * 6 * 6 * 5 = 69,120 bases in 360° => exactly 192 bases per degree.
    let rel_hd = normalize_lon_360(lon - HD_MANDALA_START_LON_DEG);
    let total_base_idx = floor_with_eps(rel_hd * 192.0).min(69_119);
    let hd_base = ((total_base_idx % 5) as u8) + 1;
    let total_tone_idx = total_base_idx / 5;
    let hd_tone = ((total_tone_idx % 6) as u8) + 1;
    let total_color_idx = total_tone_idx / 6;
    let hd_color = ((total_color_idx % 6) as u8) + 1;
    let total_line_idx = total_color_idx / 6;
    let hd_line = ((total_line_idx % 6) as u8) + 1;
    let gate_idx = (total_line_idx / 6).min(63);
    let hd_gate = HD_MANDALA_GATES_64.get(gate_idx).copied().unwrap_or(41);

    ZodiacPlacement {
        lon_deg: lon,
        sign_index,
        sign_name,
        deg_in_sign,
        dms_deg,
        dms_min,
        dms_sec,
        decan,
        nakshatra_index,
        nakshatra_name,
        nakshatra_lord,
        pada,
        hd_gate,
        hd_line,
        hd_color,
        hd_tone,
        hd_base,
    }
}

/// Placement of a body or angle inside a Vedic divisional chart (*Varga*, `D1..D60`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VargaPlacement {
    /// Harmonic division factor (`1..=60`, e.g. `9` for Navamsha `D9`).
    pub division: u16,
    /// Short Varga code (`"D1" .. "D60"`).
    pub varga_code: String,
    /// Traditional Sanskrit name of the Varga (`"Rashi"`, `"Hora"`, `"Drekkana"`, `"Navamsha"`, etc.).
    pub varga_name: String,
    /// Input sidereal ecliptic longitude in degrees `[0.0, 360.0)`.
    pub source_sidereal_lon_deg: f64,
    /// Resulting Varga sign index (`0..12`, `0 = Aries .. 11 = Pisces`).
    pub sign_index: u8,
    /// Resulting Varga sign name (`"Aries" .. "Pisces"`).
    pub sign_name: String,
    /// Planetary ruler of the resulting Varga sign (`"Mars" .. "Jupiter"`).
    pub sign_lord: String,
    /// 1-based slice (*amsha*) index within the source sign (`1..=division`, or `1..=5` for unequal `D30`).
    pub amsha_index: u16,
    /// Projected longitude in the Varga chart `[0.0, 360.0)`.
    pub varga_lon_deg: f64,
    /// Projected degrees within the Varga sign `[0.0, 30.0)`.
    pub deg_in_varga_sign: f64,
}

/// Returns the canonical Sanskrit name for a Parashara Varga division (`1..=60`).
#[must_use]
pub fn varga_sanskrit_name(division: u16) -> &'static str {
    match division {
        1 => "Rashi",
        2 => "Hora",
        3 => "Drekkana",
        4 => "Chaturthamsha",
        5 => "Panchamsha",
        6 => "Shashthamsha",
        7 => "Saptamsha",
        8 => "Ashtamsha",
        9 => "Navamsha",
        10 => "Dashamsha",
        11 => "Rudramsha",
        12 => "Dwadashamsha",
        16 => "Shodashamsha",
        20 => "Vimshamsha",
        24 => "Chaturvimshamsha",
        27 => "Saptavimshamsha",
        30 => "Trimshamsha",
        40 => "Khavedamsha",
        45 => "Akshavedamsha",
        60 => "Shashtiamsha",
        _ => "Varga",
    }
}

/// Evaluates the Parashara divisional chart (`Varga`, `D1..D60`) placement for `sidereal_lon_deg`.
///
/// Implements the classical *Brihat Parashara Hora Shastra* rules for all 16 *Shodashavarga*
/// (`D1, D2, D3, D4, D7, D9, D10, D12, D16, D20, D24, D27, D30, D40, D45, D60`), plus the
/// continuous harmonic formula for any other `division` in `1..=60`.
pub fn calc_varga(sidereal_lon_deg: f64, division: u16) -> Result<VargaPlacement, EphemerisError> {
    if !sidereal_lon_deg.is_finite() || !(1..=60).contains(&division) {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid calc_varga parameters: sidereal_lon_deg={sidereal_lon_deg}, division={division} (expected 1..=60)"
            ),
        });
    }

    let lon = normalize_lon_360(sidereal_lon_deg);
    let s = floor_with_eps(lon / 30.0).min(11);
    let deg_in_sign = (lon - (s as f64) * 30.0).clamp(0.0, 30.0 - 1.0e-12);
    let is_odd_sign = s % 2 == 0; // 0 = Aries (1st sign, odd), 1 = Taurus (2nd sign, even)

    // Special unequal Parashara rule for D30 (Trimshamsha)
    if division == 30 {
        let (varga_sign, amsha_idx, seg_start, seg_span) = if is_odd_sign {
            // Odd signs: Mars 5° (Aries=0), Saturn 5° (Aquarius=10), Jupiter 8° (Sagittarius=8), Mercury 7° (Gemini=2), Venus 5° (Libra=6)
            if deg_in_sign < 5.0 {
                (0usize, 1u16, 0.0, 5.0)
            } else if deg_in_sign < 10.0 {
                (10usize, 2u16, 5.0, 5.0)
            } else if deg_in_sign < 18.0 {
                (8usize, 3u16, 10.0, 8.0)
            } else if deg_in_sign < 25.0 {
                (2usize, 4u16, 18.0, 7.0)
            } else {
                (6usize, 5u16, 25.0, 5.0)
            }
        } else {
            // Even signs: Venus 5° (Taurus=1), Mercury 7° (Virgo=5), Jupiter 8° (Pisces=11), Saturn 5° (Capricorn=9), Mars 5° (Scorpio=7)
            if deg_in_sign < 5.0 {
                (1usize, 1u16, 0.0, 5.0)
            } else if deg_in_sign < 12.0 {
                (5usize, 2u16, 5.0, 7.0)
            } else if deg_in_sign < 20.0 {
                (11usize, 3u16, 12.0, 8.0)
            } else if deg_in_sign < 25.0 {
                (9usize, 4u16, 20.0, 5.0)
            } else {
                (7usize, 5u16, 25.0, 5.0)
            }
        };
        let deg_in_varga_sign =
            (((deg_in_sign - seg_start) / seg_span) * 30.0).clamp(0.0, 30.0 - 1.0e-12);
        let varga_lon_deg = (varga_sign as f64) * 30.0 + deg_in_varga_sign;

        return Ok(VargaPlacement {
            division: 30,
            varga_code: "D30".to_owned(),
            varga_name: varga_sanskrit_name(30).to_owned(),
            source_sidereal_lon_deg: lon,
            sign_index: varga_sign as u8,
            sign_name: ZODIAC_SIGN_NAMES
                .get(varga_sign)
                .copied()
                .unwrap_or("Aries")
                .to_owned(),
            sign_lord: ZODIAC_SIGN_LORDS
                .get(varga_sign)
                .copied()
                .unwrap_or("Mars")
                .to_owned(),
            amsha_index: amsha_idx,
            varga_lon_deg,
            deg_in_varga_sign,
        });
    }

    let n = division as usize;
    let slice_width = 30.0 / (n as f64);
    let k = floor_with_eps((deg_in_sign * (n as f64)) / 30.0).min(n - 1);
    let amsha_index = (k as u16) + 1;
    let deg_in_varga_sign = (((deg_in_sign - (k as f64) * slice_width) / slice_width) * 30.0)
        .clamp(0.0, 30.0 - 1.0e-12);

    let varga_sign = match division {
        1 => s,
        2 => {
            // Parashara Hora: Odd sign -> [Leo(4), Cancer(3)]; Even sign -> [Cancer(3), Leo(4)]
            if is_odd_sign {
                if k == 0 { 4 } else { 3 }
            } else if k == 0 {
                3
            } else {
                4
            }
        }
        3 => (s + k * 4) % 12,
        4 => (s + k * 3) % 12,
        7 => {
            if is_odd_sign {
                (s + k) % 12
            } else {
                (s + 6 + k) % 12
            }
        }
        9 => {
            // Fire(0)->Aries(0), Earth(1)->Capricorn(9), Air(2)->Libra(6), Water(3)->Cancer(3)
            (s * 9 + k) % 12
        }
        10 => {
            if is_odd_sign {
                (s + k) % 12
            } else {
                (s + 8 + k) % 12
            }
        }
        12 => (s + k) % 12,
        16 => {
            // Movable(0)->Aries(0), Fixed(1)->Leo(4), Dual(2)->Sagittarius(8)
            let start = match s % 3 {
                0 => 0,
                1 => 4,
                _ => 8,
            };
            (start + k) % 12
        }
        20 => {
            // Movable(0)->Aries(0), Fixed(1)->Sagittarius(8), Dual(2)->Leo(4)
            let start = match s % 3 {
                0 => 0,
                1 => 8,
                _ => 4,
            };
            (start + k) % 12
        }
        24 => {
            // Odd->Leo(4), Even->Cancer(3)
            let start = if is_odd_sign { 4 } else { 3 };
            (start + k) % 12
        }
        27 => {
            // Fire(0)->Aries(0), Earth(1)->Cancer(3), Air(2)->Libra(6), Water(3)->Capricorn(9)
            (s * 27 + k) % 12
        }
        40 => {
            // Odd->Aries(0), Even->Libra(6)
            let start = if is_odd_sign { 0 } else { 6 };
            (start + k) % 12
        }
        45 => {
            // Movable(0)->Aries(0), Fixed(1)->Leo(4), Dual(2)->Sagittarius(8)
            let start = match s % 3 {
                0 => 0,
                1 => 4,
                _ => 8,
            };
            (start + k) % 12
        }
        // BPHS Ch. 6 vv. 33–35: "ignore the Rāśi position of a Grah, take the degrees etc.
        // it traversed in that Rāśi, multiply that figure by 2, divide the degrees by 12,
        // add 1 to the remainder; from Mesha (Aries) count..."
        60 => k % 12,
        _ => (s * n + k) % 12,
    };

    let varga_lon_deg = (varga_sign as f64) * 30.0 + deg_in_varga_sign;

    Ok(VargaPlacement {
        division,
        varga_code: format!("D{division}"),
        varga_name: varga_sanskrit_name(division).to_owned(),
        source_sidereal_lon_deg: lon,
        sign_index: varga_sign as u8,
        sign_name: ZODIAC_SIGN_NAMES
            .get(varga_sign)
            .copied()
            .unwrap_or("Aries")
            .to_owned(),
        sign_lord: ZODIAC_SIGN_LORDS
            .get(varga_sign)
            .copied()
            .unwrap_or("Mars")
            .to_owned(),
        amsha_index,
        varga_lon_deg,
        deg_in_varga_sign,
    })
}

/// Evaluates all 16 classical Parashara *Shodashavarga* divisions (`D1, D2, D3, D4, D7, D9, D10, D12, D16, D20, D24, D27, D30, D40, D45, D60`).
pub fn calc_shodashavarga(sidereal_lon_deg: f64) -> Result<Vec<VargaPlacement>, EphemerisError> {
    const DIVISIONS: [u16; 16] = [1, 2, 3, 4, 7, 9, 10, 12, 16, 20, 24, 27, 30, 40, 45, 60];
    DIVISIONS
        .iter()
        .map(|&d| calc_varga(sidereal_lon_deg, d))
        .collect()
}

// ============================================================================
// 2. Vedic Panchanga and Muhurtas of the Day (`calc_panchanga`)
// ============================================================================

/// Start and end timestamps (`JD(UT1)`) of an auspicious or inauspicious daily window (*Muhurta* / *Kalam*).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TimeWindow {
    /// Window start in `JD(UT1)`.
    pub start_jd_ut1: f64,
    /// Window end in `JD(UT1)`.
    pub end_jd_ut1: f64,
}

/// Vedic civil weekday (*Vara*, `0 = Ravivara / Sunday .. 6 = Shanivara / Saturday`) governed by local sunrise.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct VaraInfo {
    /// Weekday index `0..=6` (`0 = Sunday .. 6 = Saturday`).
    pub vara_index: u8,
    /// Sanskrit weekday name (`"Ravivara" .. "Shanivara"`).
    pub vara_name: String,
    /// Planetary lord of the day (`"Sun" .. "Saturn"`).
    pub vara_lord: String,
}

/// Lunar day (*Tithi*, `1..=30`) with exact start and end timestamps (`JD(UT1)`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TithiInfo {
    /// 1-based Tithi index (`1..=30`, where `15` = Purnima, `30` = Amavasya).
    pub tithi_index: u8,
    /// Canonical Sanskrit name of the Tithi (`"Pratipada" .. "Amavasya"`).
    pub tithi_name: String,
    /// Lunar fortnight (`"Shukla"` for `1..=15`, `"Krishna"` for `16..=30`).
    pub paksha: String,
    /// Fractional progress through the current Tithi in `[0.0, 1.0]`.
    pub elapsed_fraction: f64,
    /// Exact start of the current Tithi in `JD(UT1)`.
    pub start_jd_ut1: f64,
    /// Exact end of the current Tithi in `JD(UT1)`.
    pub end_jd_ut1: f64,
}

/// Lunar mansion (*Nakshatra*, `1..=27`) and *Pada* (`1..=4`) with exact start and end timestamps (`JD(UT1)`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NakshatraInfo {
    /// 1-based Nakshatra index (`1..=27`).
    pub nakshatra_index: u8,
    /// Canonical Sanskrit name (`"Ashwini" .. "Revati"`).
    pub nakshatra_name: String,
    /// Vimshottari planetary lord (`"Ketu" .. "Mercury"`).
    pub vimshottari_lord: String,
    /// Quarter (*Pada*) within the Nakshatra (`1..=4`).
    pub pada: u8,
    /// Fractional progress through the current Nakshatra in `[0.0, 1.0]`.
    pub elapsed_fraction: f64,
    /// Exact start of the current Nakshatra in `JD(UT1)`.
    pub start_jd_ut1: f64,
    /// Exact end of the current Nakshatra in `JD(UT1)`.
    pub end_jd_ut1: f64,
}

/// Luni-solar *Yoga* (`1..=27`) defined by $(\lambda_{\mathbb{C},\text{sid}} + \lambda_{\odot,\text{sid}}) \bmod 360^\circ$.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct YogaInfo {
    /// 1-based Yoga index (`1..=27`).
    pub yoga_index: u8,
    /// Canonical Sanskrit name (`"Vishkumbha" .. "Vaidhriti"`).
    pub yoga_name: String,
    /// Sum of sidereal lunar and solar longitudes in degrees `[0.0, 360.0)`.
    pub sum_lon_deg: f64,
    /// Fractional progress through the current Yoga in `[0.0, 1.0]`.
    pub elapsed_fraction: f64,
    /// Exact start of the current Yoga in `JD(UT1)`.
    pub start_jd_ut1: f64,
    /// Exact end of the current Yoga in `JD(UT1)`.
    pub end_jd_ut1: f64,
}

/// Half-tithi *Karana* (`1..=60`) with canonical name out of the 11 Karanas (4 fixed + 7 movable).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KaranaInfo {
    /// 1-based half-tithi index (`1..=60`, each $6^\circ$ of lunar-solar elongation).
    pub karana_index: u8,
    /// Canonical Karana name (`"Kimstughna"`, `"Bava".."Vishti"`, `"Shakuni"`, `"Chatushpada"`, `"Naga"`).
    pub karana_name: String,
    /// `true` for `"Vishti"` (*Bhadra* Karana).
    pub is_vishti_bhadra: bool,
    /// Fractional progress through the current Karana in `[0.0, 1.0]`.
    pub elapsed_fraction: f64,
    /// Exact start of the current Karana in `JD(UT1)`.
    pub start_jd_ut1: f64,
    /// Exact end of the current Karana in `JD(UT1)`.
    pub end_jd_ut1: f64,
}

/// Complete Vedic *Panchanga* (5 limbs: *Vara*, *Tithi*, *Nakshatra*, *Yoga*, *Karana*) and daily *Muhurta* windows.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PanchangaReport {
    /// Evaluated `JD(TT)` of the input instant.
    pub jd_tt: f64,
    /// Evaluated `JD(UT1)` of the input instant.
    pub jd_ut1: f64,
    /// Ayanamsa model used for sidereal Nakshatra and Yoga calculations.
    pub ayanamsa_model: AyanamsaModel,
    /// Evaluated true ayanamsa offset in degrees at `jd_tt`.
    pub ayanamsa_deg: f64,
    /// Start of the current Vedic civil day (most recent local sunrise at or before `jd_ut1`) in `JD(UT1)`.
    pub sunrise_jd_ut1: f64,
    /// Local sunset of the current Vedic civil day in `JD(UT1)`.
    pub sunset_jd_ut1: f64,
    /// Next local sunrise (end of the current Vedic civil day) in `JD(UT1)`.
    pub next_sunrise_jd_ut1: f64,
    /// `true` if high-latitude polar fallback (local solar noon $\mp 6\text{ h}$) was used because the Sun was circumpolar.
    pub polar_day_fallback: bool,
    /// 1st limb: *Vara* (weekday governed by `sunrise_jd_ut1`).
    pub vara: VaraInfo,
    /// 2nd limb: *Tithi* (lunar day `1..=30`).
    pub tithi: TithiInfo,
    /// 3rd limb: *Nakshatra* (`1..=27` and *Pada* `1..=4`).
    pub nakshatra: NakshatraInfo,
    /// 4th limb: *Yoga* (`1..=27`).
    pub yoga: YogaInfo,
    /// 5th limb: *Karana* (`1..=60`).
    pub karana: KaranaInfo,
    /// Inauspicious *Rahu Kalam* window (1/8th of daytime arc `[sunrise, sunset]`).
    pub rahukalam: TimeWindow,
    /// Inauspicious *Yamagandam* window (1/8th of daytime arc `[sunrise, sunset]`).
    pub yamagandam: TimeWindow,
    /// *Gulika Kalam* window (1/8th of daytime arc `[sunrise, sunset]`).
    pub gulikakalam: TimeWindow,
    /// Auspicious midday *Abhijit Muhurta* window (8th of 15 daytime muhurtas, `mid ± D/30`).
    pub abhijit_muhurta: TimeWindow,
}

/// Returns the canonical name of a Karana for `karana_index` in `1..=60`.
#[must_use]
pub fn karana_name_for_index(karana_index: u8) -> &'static str {
    match karana_index {
        1 => "Kimstughna",
        2..=57 => CHARA_KARANA_NAMES_7
            .get(((karana_index - 2) % 7) as usize)
            .copied()
            .unwrap_or("Bava"),
        58 => "Shakuni",
        59 => "Chatushpada",
        _ => "Naga",
    }
}

/// Resolves the Vedic civil day bounds `(sunrise_jd_ut1, sunset_jd_ut1, next_sunrise_jd_ut1, polar_day_fallback)`
/// containing `jd_ut1` at `observer`.
///
/// If `jd_ut1` occurs before local sunrise on a civil calendar date, the active Vedic day began at the
/// *preceding* sunrise (`sunrise_jd_ut1 <= jd_ut1 < next_sunrise_jd_ut1`).
/// At polar latitudes where the Sun does not rise or set within the diurnal cycle, falls back
/// deterministically to local solar noon $\mp 6\text{ h}$ (`12\text{ h}$ daytime arc).
pub fn resolve_vedic_day_sun_windows(
    engine: &Engine,
    jd_ut1: f64,
    observer: ObserverLocation,
) -> Result<(f64, f64, f64, bool), EphemerisError> {
    let sun = BodyId::Planet(PlanetId::Sun);
    let rt_opts = RiseTransOptions {
        max_search_days: 1.15,
        ..RiseTransOptions::almanac_upper_limb()
    };

    // Search from jd_ut1 - 1.04 days so the first found sunrise is normally the most recent sunrise <= jd_ut1.
    if let Some(mut sr) =
        engine.rise_trans(sun, observer, jd_ut1 - 1.04, RiseTransEvent::Rise, rt_opts)?
    {
        if sr > jd_ut1 + 1.0e-9 {
            // Edge case near high-latitude seasonal change where sunrise drifted > 1.04d
            if let Some(earlier) = engine.rise_trans(
                sun,
                observer,
                jd_ut1 - 1.35,
                RiseTransEvent::Rise,
                RiseTransOptions {
                    max_search_days: 0.4,
                    ..rt_opts
                },
            )? {
                if earlier <= jd_ut1 + 1.0e-9 {
                    sr = earlier;
                }
            }
        } else if jd_ut1 - sr > 0.92 {
            // Check if a second sunrise occurred right before jd_ut1 (around sr + 1.0d)
            if let Some(later_sr) = engine.rise_trans(
                sun,
                observer,
                sr + 0.88,
                RiseTransEvent::Rise,
                RiseTransOptions {
                    max_search_days: 0.25,
                    ..rt_opts
                },
            )? {
                if later_sr <= jd_ut1 + 1.0e-9 {
                    sr = later_sr;
                }
            }
        }

        if sr <= jd_ut1 + 1.0e-9 {
            // Next sunrise is ~1.0 day after sr; start search at sr + 0.82 to avoid scanning the full day
            let next_sr_opt = engine
                .rise_trans(
                    sun,
                    observer,
                    sr + 0.82,
                    RiseTransEvent::Rise,
                    RiseTransOptions {
                        max_search_days: 0.45,
                        ..rt_opts
                    },
                )?
                .or_else(|| {
                    engine
                        .rise_trans(sun, observer, sr + 0.15, RiseTransEvent::Rise, rt_opts)
                        .ok()
                        .flatten()
                });

            if let Some(ss) =
                engine.rise_trans(sun, observer, sr + 1.0e-5, RiseTransEvent::Set, rt_opts)?
            {
                let next_sr = next_sr_opt.unwrap_or(sr + 1.0);
                if ss > sr && ss < next_sr {
                    let is_fallback = next_sr_opt.is_none();
                    return Ok((sr, ss, next_sr, is_fallback));
                }
            }
        }
    }

    // Deterministic polar fallback: anchor on the most recent 06:00 Local Mean Solar Time <= jd_ut1
    // (i.e., local solar noon - 6h), with sunset = sunrise + 12h (0.5 d) and next_sunrise = sunrise + 24h (1.0 d).
    let lon_days = observer.longitude_deg / 360.0;
    let shifted = jd_ut1 + lon_days + 0.25;
    let sunrise_fallback = shifted.floor() - 0.25 - lon_days;
    let sunset_fallback = sunrise_fallback + 0.5;
    let next_sunrise_fallback = sunrise_fallback + 1.0;

    Ok((
        sunrise_fallback,
        sunset_fallback,
        next_sunrise_fallback,
        true,
    ))
}

/// Evaluates the sidereal Yoga sum $S(t) = (\lambda_{\mathbb{C},\text{sid}} + \lambda_{\odot,\text{sid}}) \bmod 360^\circ$,
/// its analytical derivative $\dot{S}(t) = \dot{\lambda}_{\mathbb{C},\text{sid}} + \dot{\lambda}_{\odot,\text{sid}}$ (`deg/day`),
/// and `jd_ut1` at `jd_tt`.
fn eval_yoga_sum_at_tt(
    engine: &Engine,
    jd_tt: f64,
    ayanamsa: AyanamsaModel,
) -> Result<(f64, f64, f64), EphemerisError> {
    let sid_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Ecliptic,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Sidereal(ayanamsa),
        include_speed: true,
    };
    let moment = engine.at(TimeInput::from_tt_jd(jd_tt), None)?;
    let moon = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Moon), sid_opts)?;
    let sun = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Sun), sid_opts)?;
    let sum_deg = normalize_lon_360(moon.lon_deg + sun.lon_deg);
    let sum_speed = moon.lon_speed_deg_day + sun.lon_speed_deg_day;
    Ok((sum_deg, sum_speed, moon.jd_ut1))
}

/// Solves for the exact `JD(UT1)` when the sidereal Yoga sum $(\lambda_{\mathbb{C},\text{sid}} + \lambda_{\odot,\text{sid}}) \bmod 360^\circ$
/// reaches `target_sum_deg` starting from `jd_tt_start` (moving forward if `forward == true`, or backward if `forward == false`)
/// using a safeguarded hybrid Newton-Illinois solver with analytical derivative $\dot{\lambda}_{\mathbb{C}} + \dot{\lambda}_\odot$.
fn search_yoga_boundary_ut1(
    engine: &Engine,
    jd_tt_start: f64,
    target_sum_deg: f64,
    forward: bool,
    ayanamsa: AyanamsaModel,
) -> Result<f64, EphemerisError> {
    let target = normalize_lon_360(target_sum_deg);
    let dir = if forward { 1.0 } else { -1.0 };
    let step_days = 0.5_f64;
    let max_steps = 8usize;

    let mut t_prev = jd_tt_start;
    let (mut sum_prev, mut speed_prev, mut ut1_prev) =
        eval_yoga_sum_at_tt(engine, t_prev, ayanamsa)?;

    for step in 1..=max_steps {
        let t_curr = jd_tt_start + dir * (step as f64) * step_days;
        let (sum_curr, speed_curr, ut1_curr) = eval_yoga_sum_at_tt(engine, t_curr, ayanamsa)?;

        let f_prev = wrap_pm180(sum_prev - target);
        let f_curr = wrap_pm180(sum_curr - target);

        if f_prev.abs() < 60.0 && f_curr.abs() < 60.0 && f_prev * f_curr <= 0.0 {
            let mut a = t_prev;
            let mut b = t_curr;
            let mut fa = f_prev;
            let mut fb = f_curr;
            let mut va = speed_prev;
            let mut vb = speed_curr;
            let mut best_ut1 = if fa.abs() <= fb.abs() {
                ut1_prev
            } else {
                ut1_curr
            };
            let mut side = 0i8;

            if fa.abs() < 1.0e-9 {
                return Ok(ut1_prev);
            }
            if fb.abs() < 1.0e-9 {
                return Ok(ut1_curr);
            }

            for _ in 0..40 {
                let span = (b - a).abs();
                if span < 1.0e-10 {
                    return Ok(best_ut1);
                }

                let (x_curr, f_c, v_c) = if fa.abs() <= fb.abs() {
                    (a, fa, va)
                } else {
                    (b, fb, vb)
                };
                let t_min = a.min(b);
                let t_max = a.max(b);
                let pad = 0.02 * span;

                let c_newton = if v_c.abs() > 1.0 {
                    Some(x_curr - f_c / v_c)
                } else {
                    None
                };
                let c = match c_newton {
                    Some(cn) if cn > t_min + pad && cn < t_max - pad => cn,
                    _ => {
                        let denom = fb - fa;
                        if denom.abs() > 1.0e-15 {
                            ((a * fb - b * fa) / denom).clamp(t_min + pad, t_max - pad)
                        } else {
                            0.5 * (a + b)
                        }
                    }
                };

                let (sum_c, vc, ut1_c) = eval_yoga_sum_at_tt(engine, c, ayanamsa)?;
                let fc = wrap_pm180(sum_c - target);

                if fc.abs() < fa.abs().min(fb.abs()) {
                    best_ut1 = ut1_c;
                }
                if fc.abs() < 1.0e-9 {
                    return Ok(ut1_c);
                }

                if (fc > 0.0) == (fb > 0.0) {
                    b = c;
                    fb = fc;
                    vb = vc;
                    if side == -1 {
                        fa *= 0.5;
                    }
                    side = -1;
                } else {
                    a = c;
                    fa = fc;
                    va = vc;
                    if side == 1 {
                        fb *= 0.5;
                    }
                    side = 1;
                }
            }

            return Ok(best_ut1);
        }

        t_prev = t_curr;
        sum_prev = sum_curr;
        speed_prev = speed_curr;
        ut1_prev = ut1_curr;
    }

    Err(EphemerisError::NonConvergence {
        solver: "search_yoga_boundary_ut1".to_owned(),
        iterations: max_steps as u32,
        residual: wrap_pm180(sum_prev - target).abs(),
    })
}

/// Computes the complete Vedic *Panchanga* (`Vara`, `Tithi`, `Nakshatra`, `Yoga`, `Karana`) and daily
/// *Muhurta* / *Kalam* windows (`Rahukalam`, `Yamagandam`, `Gulikakalam`, `Abhijit Muhurta`) at `input` for `observer`.
pub fn calc_panchanga(
    engine: &Engine,
    input: TimeInput,
    observer: ObserverLocation,
    ayanamsa: AyanamsaModel,
) -> Result<PanchangaReport, EphemerisError> {
    let moment = engine.at(input, Some(observer))?;
    let jd_tt = moment.center().tt.jd_f64();
    let jd_ut1 = moment.center().ut1.jd_f64();

    let ayan_val = engine.ayanamsa(TimeInput::from_tt_jd(jd_tt), ayanamsa, NutationMode::True)?;
    let ayanamsa_deg = ayan_val.effective_deg;

    // 1. Local Vedic day sunrise, sunset, next_sunrise & Vara
    let (sunrise_jd_ut1, sunset_jd_ut1, next_sunrise_jd_ut1, polar_day_fallback) =
        resolve_vedic_day_sun_windows(engine, jd_ut1, observer)?;

    let local_sunrise_jd = sunrise_jd_ut1 + observer.longitude_deg / 360.0;
    let vara_index = (((local_sunrise_jd + 0.5).floor() as i64) + 1).rem_euclid(7) as u8;
    let v_idx = vara_index as usize;
    let vara = VaraInfo {
        vara_index,
        vara_name: VARA_NAMES_7
            .get(v_idx)
            .copied()
            .unwrap_or("Ravivara")
            .to_owned(),
        vara_lord: VARA_LORDS_7
            .get(v_idx)
            .copied()
            .unwrap_or("Sun")
            .to_owned(),
    };

    // Evaluate Sun and Moon (tropical and sidereal)
    let trop_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Ecliptic,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Tropical,
        include_speed: true,
    };
    let sid_opts = CalcOptions {
        zodiac: Zodiac::Sidereal(ayanamsa),
        ..trop_opts
    };

    let sun_trop = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Sun), trop_opts)?;
    let moon_trop = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Moon), trop_opts)?;
    let sun_sid = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Sun), sid_opts)?;
    let moon_sid = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Moon), sid_opts)?;

    // 2. Tithi (1..=30, each 12° of lunar-solar elongation)
    let elong_deg = normalize_lon_360(moon_trop.lon_deg - sun_trop.lon_deg);
    let tithi_0 = floor_with_eps(elong_deg / 12.0).min(29);
    let tithi_index = (tithi_0 as u8) + 1;
    let tithi_start_elong = (tithi_0 as f64) * 12.0;
    let tithi_end_elong = normalize_lon_360((tithi_index as f64) * 12.0);
    let tithi_elapsed = ((elong_deg - tithi_start_elong) / 12.0).clamp(0.0, 1.0);

    let tithi_start_ev = search_lunar_phase(
        engine,
        TimeInput::from_tt_jd(jd_tt + 1.0e-4),
        tithi_start_elong,
        false,
    )?;
    let tithi_end_ev = search_lunar_phase(
        engine,
        TimeInput::from_tt_jd(jd_tt - 1.0e-4),
        tithi_end_elong,
        true,
    )?;

    let tithi = TithiInfo {
        tithi_index,
        tithi_name: TITHI_NAMES_30
            .get(tithi_0)
            .copied()
            .unwrap_or("Pratipada")
            .to_owned(),
        paksha: if tithi_index <= 15 {
            "Shukla".to_owned()
        } else {
            "Krishna".to_owned()
        },
        elapsed_fraction: tithi_elapsed,
        start_jd_ut1: tithi_start_ev.jd_ut1,
        end_jd_ut1: tithi_end_ev.jd_ut1,
    };

    // 3. Nakshatra (1..=27) and Pada (1..=4) from sidereal lunar longitude
    let moon_sid_lon = normalize_lon_360(moon_sid.lon_deg);
    let nak_span = 360.0 / 27.0;
    let total_pada_idx = floor_with_eps((moon_sid_lon * 108.0) / 360.0).min(107);
    let nak_0 = (total_pada_idx / 4).min(26);
    let pada = ((total_pada_idx % 4) as u8) + 1;
    let nakshatra_index = (nak_0 as u8) + 1;
    let nak_start_lon = (nak_0 as f64) * nak_span;
    let nak_end_lon = normalize_lon_360((nakshatra_index as f64) * nak_span);
    let nak_elapsed = ((moon_sid_lon - nak_start_lon) / nak_span).clamp(0.0, 1.0);

    let nak_start_ev = search_body_longitude(
        engine,
        BodyId::Planet(PlanetId::Moon),
        TimeInput::from_tt_jd(jd_tt + 1.0e-4),
        nak_start_lon,
        false,
        5.0,
        &sid_opts,
    )?;
    let nak_end_ev = search_body_longitude(
        engine,
        BodyId::Planet(PlanetId::Moon),
        TimeInput::from_tt_jd(jd_tt - 1.0e-4),
        nak_end_lon,
        true,
        5.0,
        &sid_opts,
    )?;

    let nakshatra = NakshatraInfo {
        nakshatra_index,
        nakshatra_name: NAKSHATRA_NAMES_27
            .get(nak_0)
            .copied()
            .unwrap_or("Ashwini")
            .to_owned(),
        vimshottari_lord: VIMSHOTTARI_LORDS_9
            .get(nak_0 % 9)
            .copied()
            .unwrap_or("Ketu")
            .to_owned(),
        pada,
        elapsed_fraction: nak_elapsed,
        start_jd_ut1: nak_start_ev.jd_ut1,
        end_jd_ut1: nak_end_ev.jd_ut1,
    };

    // 4. Yoga (1..=27) from sum of sidereal Moon and Sun longitudes
    let yoga_sum_deg = normalize_lon_360(moon_sid.lon_deg + sun_sid.lon_deg);
    let yoga_span = 360.0 / 27.0;
    let yoga_0 = floor_with_eps((yoga_sum_deg * 27.0) / 360.0).min(26);
    let yoga_index = (yoga_0 as u8) + 1;
    let yoga_start_deg = (yoga_0 as f64) * yoga_span;
    let yoga_end_deg = normalize_lon_360((yoga_index as f64) * yoga_span);
    let yoga_elapsed = ((yoga_sum_deg - yoga_start_deg) / yoga_span).clamp(0.0, 1.0);

    let yoga_start_ut1 =
        search_yoga_boundary_ut1(engine, jd_tt + 1.0e-4, yoga_start_deg, false, ayanamsa)?;
    let yoga_end_ut1 =
        search_yoga_boundary_ut1(engine, jd_tt - 1.0e-4, yoga_end_deg, true, ayanamsa)?;

    let yoga = YogaInfo {
        yoga_index,
        yoga_name: YOGA_NAMES_27
            .get(yoga_0)
            .copied()
            .unwrap_or("Vishkumbha")
            .to_owned(),
        sum_lon_deg: yoga_sum_deg,
        elapsed_fraction: yoga_elapsed,
        start_jd_ut1: yoga_start_ut1,
        end_jd_ut1: yoga_end_ut1,
    };

    // 5. Karana (1..=60, each 6° of lunar-solar elongation)
    let karana_0 = floor_with_eps(elong_deg / 6.0).min(59);
    let karana_index = (karana_0 as u8) + 1;
    let karana_start_elong = (karana_0 as f64) * 6.0;
    let karana_end_elong = normalize_lon_360((karana_index as f64) * 6.0);
    let karana_elapsed = ((elong_deg - karana_start_elong) / 6.0).clamp(0.0, 1.0);
    let karana_name = karana_name_for_index(karana_index).to_owned();

    let karana_start_ev = if karana_0 % 2 == 0 {
        tithi_start_ev
    } else {
        search_lunar_phase(
            engine,
            TimeInput::from_tt_jd(jd_tt + 1.0e-4),
            karana_start_elong,
            false,
        )?
    };
    let karana_end_ev = if karana_index % 2 == 0 {
        tithi_end_ev
    } else {
        search_lunar_phase(
            engine,
            TimeInput::from_tt_jd(jd_tt - 1.0e-4),
            karana_end_elong,
            true,
        )?
    };

    let karana = KaranaInfo {
        karana_index,
        is_vishti_bhadra: karana_name == "Vishti",
        karana_name,
        elapsed_fraction: karana_elapsed,
        start_jd_ut1: karana_start_ev.jd_ut1,
        end_jd_ut1: karana_end_ev.jd_ut1,
    };

    // 6. Muhurtas and Kalam windows of the day
    // Daytime arc D = sunset - sunrise divided into 8 equal segments (1-based octants):
    // Rahukalam:   Sun->8, Mon->2, Tue->7, Wed->5, Thu->6, Fri->4, Sat->3
    // Yamagandam:  Sun->5, Mon->4, Tue->3, Wed->2, Thu->1, Fri->7, Sat->6
    // Gulikakalam: Sun->7, Mon->6, Tue->5, Wed->4, Thu->3, Fri->2, Sat->1
    const RAHU_OCTANT: [u8; 7] = [8, 2, 7, 5, 6, 4, 3];
    const YAMA_OCTANT: [u8; 7] = [5, 4, 3, 2, 1, 7, 6];
    const GULIKA_OCTANT: [u8; 7] = [7, 6, 5, 4, 3, 2, 1];

    let day_span = sunset_jd_ut1 - sunrise_jd_ut1;
    let octant_len = day_span / 8.0;
    let make_octant_window = |octant_1based: u8| -> TimeWindow {
        let idx = f64::from(octant_1based - 1);
        TimeWindow {
            start_jd_ut1: sunrise_jd_ut1 + idx * octant_len,
            end_jd_ut1: sunrise_jd_ut1 + (idx + 1.0) * octant_len,
        }
    };

    let rahukalam = make_octant_window(RAHU_OCTANT.get(v_idx).copied().unwrap_or(8));
    let yamagandam = make_octant_window(YAMA_OCTANT.get(v_idx).copied().unwrap_or(5));
    let gulikakalam = make_octant_window(GULIKA_OCTANT.get(v_idx).copied().unwrap_or(7));

    let mid_day = sunrise_jd_ut1 + 0.5 * day_span;
    let half_muhurta = day_span / 30.0;
    let abhijit_muhurta = TimeWindow {
        start_jd_ut1: mid_day - half_muhurta,
        end_jd_ut1: mid_day + half_muhurta,
    };

    Ok(PanchangaReport {
        jd_tt,
        jd_ut1,
        ayanamsa_model: ayanamsa,
        ayanamsa_deg,
        sunrise_jd_ut1,
        sunset_jd_ut1,
        next_sunrise_jd_ut1,
        polar_day_fallback,
        vara,
        tithi,
        nakshatra,
        yoga,
        karana,
        rahukalam,
        yamagandam,
        gulikakalam,
        abhijit_muhurta,
    })
}

// ============================================================================
// 3. Vimshottari Dasha (`vimshottari_dasha` & `calc_vimshottari_dasha`)
// ============================================================================

/// Sub-period (*Antardasha* / *Bhukti*) within a Vimshottari *Mahadasha*.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DashaSubPeriod {
    /// Governing Mahadasha lord (`"Ketu" .. "Mercury"`).
    pub maha_lord: String,
    /// Antardasha sub-lord (`"Ketu" .. "Mercury"`).
    pub lord: String,
    /// Start timestamp in `JD(UT1)`.
    pub start_jd_ut1: f64,
    /// End timestamp in `JD(UT1)`.
    pub end_jd_ut1: f64,
    /// Duration of this sub-period in solar years.
    pub duration_years: f64,
}

/// Major planetary period (*Mahadasha*) in the 120-year Vimshottari cycle, containing its 9 *Antardashas*.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DashaPeriod {
    /// Mahadasha planetary lord (`"Ketu"`, `"Venus"`, `"Sun"`, `"Moon"`, `"Mars"`, `"Rahu"`, `"Jupiter"`, `"Saturn"`, `"Mercury"`).
    pub lord: String,
    /// Start timestamp in `JD(UT1)` (`birth_jd_ut1` for the opening Mahadasha balance).
    pub start_jd_ut1: f64,
    /// End timestamp in `JD(UT1)`.
    pub end_jd_ut1: f64,
    /// Effective duration of this Mahadasha in solar years (remaining balance at birth for the 1st period, full duration for periods `2..=9`).
    pub duration_years: f64,
    /// The 9 Antardasha sub-periods within this Mahadasha.
    pub sub_periods: Vec<DashaSubPeriod>,
}

/// Complete Vimshottari Dasha report evaluated for a birth epoch.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VimshottariDashaReport {
    /// Birth timestamp in `JD(UT1)`.
    pub birth_jd_ut1: f64,
    /// Sidereal lunar longitude at birth in degrees `[0.0, 360.0)`.
    pub moon_sidereal_lon_deg: f64,
    /// Birth Nakshatra index (`1..=27`).
    pub birth_nakshatra_index: u8,
    /// Birth Nakshatra name (`"Ashwini" .. "Revati"`).
    pub birth_nakshatra_name: String,
    /// Birth Pada (`1..=4`).
    pub birth_pada: u8,
    /// Fractional remaining balance of the birth Nakshatra in `(0.0, 1.0]`.
    pub balance_fraction: f64,
    /// Remaining years of the opening Mahadasha at birth.
    pub balance_years: f64,
    /// Length of 1 solar year in days used for `JD(UT1)` conversion (`365.2425` default).
    pub year_length_days: f64,
    /// The 9 Mahadasha periods of the 120-year cycle starting from `birth_jd_ut1`.
    pub periods: Vec<DashaPeriod>,
}

/// Computes the 9 [`DashaPeriod`] entries (opening Mahadasha balance + 8 full subsequent Mahadashas,
/// each with 9 [`DashaSubPeriod`] Antardashas) directly from `moon_sidereal_lon_deg` and `birth_jd_ut1`.
///
/// `year_length_days` defaults to [`DEFAULT_DASHAS_YEAR_DAYS`] (`365.2425`) when `None` (supports `365.25`, `365.24219878`, etc.).
pub fn vimshottari_dasha(
    moon_sidereal_lon_deg: f64,
    birth_jd_ut1: f64,
    year_length_days: Option<f64>,
) -> Result<Vec<DashaPeriod>, EphemerisError> {
    let report = vimshottari_dasha_from_lon(moon_sidereal_lon_deg, birth_jd_ut1, year_length_days)?;
    Ok(report.periods)
}

/// Computes the full [`VimshottariDashaReport`] directly from `moon_sidereal_lon_deg` and `birth_jd_ut1`.
pub fn vimshottari_dasha_from_lon(
    moon_sidereal_lon_deg: f64,
    birth_jd_ut1: f64,
    year_length_days: Option<f64>,
) -> Result<VimshottariDashaReport, EphemerisError> {
    let year_days = year_length_days.unwrap_or(DEFAULT_DASHAS_YEAR_DAYS);
    if !moon_sidereal_lon_deg.is_finite()
        || !birth_jd_ut1.is_finite()
        || !year_days.is_finite()
        || !(350.0..=370.0).contains(&year_days)
    {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid vimshottari_dasha parameters: moon_sidereal_lon_deg={moon_sidereal_lon_deg}, birth_jd_ut1={birth_jd_ut1}, year_days={year_days}"
            ),
        });
    }

    let lon = normalize_lon_360(moon_sidereal_lon_deg);
    let nak_span = 360.0 / 27.0;
    let total_pada_idx = floor_with_eps((lon * 108.0) / 360.0).min(107);
    let nak_0 = (total_pada_idx / 4).min(26);
    let birth_pada = ((total_pada_idx % 4) as u8) + 1;
    let birth_nakshatra_index = (nak_0 as u8) + 1;
    let birth_nakshatra_name = NAKSHATRA_NAMES_27
        .get(nak_0)
        .copied()
        .unwrap_or("Ashwini")
        .to_owned();

    let nak_start_lon = (nak_0 as f64) * nak_span;
    let elapsed_frac = ((lon - nak_start_lon) / nak_span).clamp(0.0, 1.0);
    let balance_fraction = (1.0 - elapsed_frac).clamp(0.0, 1.0);

    let start_lord_idx = nak_0 % 9;
    let start_lord_full_years = VIMSHOTTARI_YEARS_9
        .get(start_lord_idx)
        .copied()
        .unwrap_or(7.0);
    let balance_years = balance_fraction * start_lord_full_years;

    let mut periods = Vec::with_capacity(9);
    let mut cursor_jd = birth_jd_ut1;

    for step in 0..9usize {
        let m_idx = (start_lord_idx + step) % 9;
        let m_lord = VIMSHOTTARI_LORDS_9.get(m_idx).copied().unwrap_or("Ketu");
        let full_maha_years = VIMSHOTTARI_YEARS_9.get(m_idx).copied().unwrap_or(7.0);
        let duration_years = if step == 0 {
            balance_years
        } else {
            full_maha_years
        };
        let period_start_jd = cursor_jd;
        let period_end_jd = period_start_jd + duration_years * year_days;

        let mut sub_periods = Vec::with_capacity(9);
        if step == 0 {
            // Reconstruct the canonical unclipped birth Mahadasha start before birth_jd_ut1,
            // generate the 9 canonical Antardashas of full_maha_years * (s_full_years / 120.0),
            // skip Antardashas that already elapsed prior to birth_jd_ut1, and clip the active
            // birth Antardasha's start_jd_ut1 to birth_jd_ut1.
            let canonical_maha_start_jd = period_end_jd - full_maha_years * year_days;
            let mut sub_cursor_jd = canonical_maha_start_jd;
            for sub_step in 0..9usize {
                let s_idx = (m_idx + sub_step) % 9;
                let s_lord = VIMSHOTTARI_LORDS_9.get(s_idx).copied().unwrap_or("Ketu");
                let s_full_years = VIMSHOTTARI_YEARS_9.get(s_idx).copied().unwrap_or(7.0);
                let canonical_sub_years = full_maha_years * (s_full_years / 120.0);
                let sub_end_jd = if sub_step == 8 {
                    period_end_jd
                } else {
                    sub_cursor_jd + canonical_sub_years * year_days
                };
                if sub_end_jd > birth_jd_ut1 + 1.0e-12 || (sub_step == 8 && sub_periods.is_empty())
                {
                    let clipped_start_jd = sub_cursor_jd.max(birth_jd_ut1);
                    let clipped_years = ((sub_end_jd - clipped_start_jd) / year_days).max(0.0);
                    sub_periods.push(DashaSubPeriod {
                        maha_lord: m_lord.to_owned(),
                        lord: s_lord.to_owned(),
                        start_jd_ut1: clipped_start_jd,
                        end_jd_ut1: sub_end_jd,
                        duration_years: clipped_years,
                    });
                }
                sub_cursor_jd = sub_end_jd;
            }
        } else {
            let mut sub_cursor_jd = period_start_jd;
            for sub_step in 0..9usize {
                let s_idx = (m_idx + sub_step) % 9;
                let s_lord = VIMSHOTTARI_LORDS_9.get(s_idx).copied().unwrap_or("Ketu");
                let s_full_years = VIMSHOTTARI_YEARS_9.get(s_idx).copied().unwrap_or(7.0);
                let sub_years = full_maha_years * (s_full_years / 120.0);
                let sub_end_jd = if sub_step == 8 {
                    period_end_jd
                } else {
                    sub_cursor_jd + sub_years * year_days
                };
                sub_periods.push(DashaSubPeriod {
                    maha_lord: m_lord.to_owned(),
                    lord: s_lord.to_owned(),
                    start_jd_ut1: sub_cursor_jd,
                    end_jd_ut1: sub_end_jd,
                    duration_years: sub_years,
                });
                sub_cursor_jd = sub_end_jd;
            }
        }

        periods.push(DashaPeriod {
            lord: m_lord.to_owned(),
            start_jd_ut1: period_start_jd,
            end_jd_ut1: period_end_jd,
            duration_years,
            sub_periods,
        });
        cursor_jd = period_end_jd;
    }

    Ok(VimshottariDashaReport {
        birth_jd_ut1,
        moon_sidereal_lon_deg: lon,
        birth_nakshatra_index,
        birth_nakshatra_name,
        birth_pada,
        balance_fraction,
        balance_years,
        year_length_days: year_days,
        periods,
    })
}

/// Evaluates the Moon's geocentric apparent sidereal longitude at `birth_input` under `ayanamsa`
/// and computes the complete 120-year [`VimshottariDashaReport`].
pub fn calc_vimshottari_dasha(
    engine: &Engine,
    birth_input: TimeInput,
    ayanamsa: AyanamsaModel,
    year_length_days: Option<f64>,
) -> Result<VimshottariDashaReport, EphemerisError> {
    let sid_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Ecliptic,
        equinox: Equinox::OfDate,
        nutation: NutationMode::True,
        aberration: AbMode::Apparent,
        zodiac: Zodiac::Sidereal(ayanamsa),
        include_speed: false,
    };
    let moon = engine.calc(BodyId::Planet(PlanetId::Moon), birth_input, sid_opts)?;
    vimshottari_dasha_from_lon(moon.lon_deg, moon.jd_ut1, year_length_days)
}

// ============================================================================
// 4. Traditional Astrology: Chart Sect, Planetary Hours & 7 Hermetic Lots
// ============================================================================

/// Diurnal / nocturnal chart sect evaluation (`chart_sect`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ChartSectInfo {
    /// `true` if the Sun's true geometric altitude is on or above the horizon ($h_{\odot,\text{true}} \ge 0^\circ$).
    pub is_day_chart: bool,
    /// True (airless) topocentric altitude of the Sun in degrees `[-90, +90]`.
    pub sun_altitude_true_deg: f64,
    /// Apparent (refracted) topocentric altitude of the Sun in degrees `[-90, +90]`.
    pub sun_altitude_apparent_deg: f64,
    /// `true` when $|h_{\odot,\text{true}}| < 5.0^\circ$ (`R7.3`).
    pub near_horizon: bool,
}

/// Evaluates the Traditional [`ChartSectInfo`] at `input` for `observer`.
pub fn chart_sect(
    engine: &Engine,
    input: TimeInput,
    observer: ObserverLocation,
) -> Result<ChartSectInfo, EphemerisError> {
    let hz = engine.calc_azalt(
        BodyId::Planet(PlanetId::Sun),
        input,
        observer,
        AtmosphereParams::default(),
    )?;
    Ok(ChartSectInfo {
        is_day_chart: hz.altitude_true_deg >= 0.0,
        sun_altitude_true_deg: hz.altitude_true_deg,
        sun_altitude_apparent_deg: hz.altitude_apparent_deg,
        near_horizon: hz.near_horizon,
    })
}

/// Single unequal planetary hour slot (`1..=24`, where `1..=12` are daytime hours and `13..=24` are nighttime hours).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlanetaryHourSlot {
    /// 1-based hour index from sunrise (`1..=12` daytime, `13..=24` nighttime).
    pub hour_index: u8,
    /// Governing Chaldean planetary ruler (`"Saturn"`, `"Jupiter"`, `"Mars"`, `"Sun"`, `"Venus"`, `"Mercury"`, `"Moon"`).
    pub lord: String,
    /// `true` for daytime hours (`1..=12`), `false` for nighttime hours (`13..=24`).
    pub is_daytime: bool,
    /// Start timestamp of the planetary hour in `JD(UT1)`.
    pub start_jd_ut1: f64,
    /// End timestamp of the planetary hour in `JD(UT1)`.
    pub end_jd_ut1: f64,
}

/// Complete 24-hour Traditional Planetary Hours report from local sunrise to next local sunrise.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlanetaryHoursReport {
    /// Evaluated `JD(UT1)` of the query instant.
    pub jd_ut1: f64,
    /// Start of the planetary day (local sunrise) in `JD(UT1)`.
    pub sunrise_jd_ut1: f64,
    /// Local sunset in `JD(UT1)`.
    pub sunset_jd_ut1: f64,
    /// End of the planetary day (next local sunrise) in `JD(UT1)`.
    pub next_sunrise_jd_ut1: f64,
    /// `true` if high-latitude polar fallback was used.
    pub polar_day_fallback: bool,
    /// Planetary ruler of the day (and ruler of the 1st hour after sunrise).
    pub day_lord: String,
    /// Active planetary hour index (`1..=24`) containing `jd_ut1`.
    pub current_hour_index: u8,
    /// Planetary ruler of the active hour containing `jd_ut1`.
    pub current_hour_lord: String,
    /// `true` if `current_hour_index <= 12`.
    pub is_daytime_hour: bool,
    /// All 24 unequal planetary hour slots (`1..=12` daytime, `13..=24` nighttime).
    pub hours: Vec<PlanetaryHourSlot>,
}

/// Computes the 24 unequal Traditional Planetary Hours (`planetary_hours`) at `input` for `observer`.
///
/// Divides `[sunrise, sunset]` into 12 equal daytime hours and `[sunset, next_sunrise]` into 12 equal
/// nighttime hours. The 1st hour after sunrise is ruled by `day_lord` (governed by local *Vara*),
/// and subsequent hours follow the descending Chaldean order `["Saturn", "Jupiter", "Mars", "Sun", "Venus", "Mercury", "Moon"]`.
pub fn planetary_hours(
    engine: &Engine,
    input: TimeInput,
    observer: ObserverLocation,
) -> Result<PlanetaryHoursReport, EphemerisError> {
    let (ut1, _, _, _) = engine.time_table().resolve_scales(input, Some(observer));
    let jd_ut1 = ut1.jd_f64();

    let (sunrise_jd_ut1, sunset_jd_ut1, next_sunrise_jd_ut1, polar_day_fallback) =
        resolve_vedic_day_sun_windows(engine, jd_ut1, observer)?;

    let local_sunrise_jd = sunrise_jd_ut1 + observer.longitude_deg / 360.0;
    let vara_index = (((local_sunrise_jd + 0.5).floor() as i64) + 1).rem_euclid(7) as usize;
    let day_lord = VARA_LORDS_7
        .get(vara_index)
        .copied()
        .unwrap_or("Sun")
        .to_owned();

    // Chaldean index of the 1st hour ruler:
    // Sun->3, Moon->6, Mars->2, Mercury->5, Jupiter->1, Venus->4, Saturn->0
    const VARA_TO_CHALDEAN_START: [usize; 7] = [3, 6, 2, 5, 1, 4, 0];
    let start_chaldean_idx = VARA_TO_CHALDEAN_START
        .get(vara_index)
        .copied()
        .unwrap_or(3);

    let day_hour_len = (sunset_jd_ut1 - sunrise_jd_ut1) / 12.0;
    let night_hour_len = (next_sunrise_jd_ut1 - sunset_jd_ut1) / 12.0;

    let mut hours = Vec::with_capacity(24);
    for i in 0..24usize {
        let hour_index = (i as u8) + 1;
        let is_daytime = i < 12;
        let (start_jd, end_jd) = if is_daytime {
            let s = sunrise_jd_ut1 + (i as f64) * day_hour_len;
            let e = if i == 11 {
                sunset_jd_ut1
            } else {
                sunrise_jd_ut1 + ((i + 1) as f64) * day_hour_len
            };
            (s, e)
        } else {
            let ni = i - 12;
            let s = sunset_jd_ut1 + (ni as f64) * night_hour_len;
            let e = if ni == 11 {
                next_sunrise_jd_ut1
            } else {
                sunset_jd_ut1 + ((ni + 1) as f64) * night_hour_len
            };
            (s, e)
        };
        let lord = CHALDEAN_ORDER_7
            .get((start_chaldean_idx + i) % 7)
            .copied()
            .unwrap_or("Sun")
            .to_owned();
        hours.push(PlanetaryHourSlot {
            hour_index,
            lord,
            is_daytime,
            start_jd_ut1: start_jd,
            end_jd_ut1: end_jd,
        });
    }

    let current_idx_0 = if jd_ut1 < sunset_jd_ut1 {
        let frac = ((jd_ut1 - sunrise_jd_ut1) / day_hour_len.max(1.0e-12)).floor() as isize;
        frac.clamp(0, 11) as usize
    } else {
        let frac = ((jd_ut1 - sunset_jd_ut1) / night_hour_len.max(1.0e-12)).floor() as isize;
        12 + (frac.clamp(0, 11) as usize)
    };

    let current_hour_index = (current_idx_0 as u8) + 1;
    let current_hour_lord = hours
        .get(current_idx_0)
        .map_or_else(|| day_lord.clone(), |h| h.lord.clone());
    let is_daytime_hour = current_hour_index <= 12;

    Ok(PlanetaryHoursReport {
        jd_ut1,
        sunrise_jd_ut1,
        sunset_jd_ut1,
        next_sunrise_jd_ut1,
        polar_day_fallback,
        day_lord,
        current_hour_index,
        current_hour_lord,
        is_daytime_hour,
        hours,
    })
}

/// Classical Hellenistic / Hermetic Lot (*Klaros*) of Paulus Alexandrinus (`Fortune`, `Spirit`, `Eros`, `Necessity`, `Courage`, `Victory`, `Nemesis`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HermeticLot {
    /// Canonical Lot name (`"Fortune"`, `"Spirit"`, `"Eros"`, `"Necessity"`, `"Courage"`, `"Victory"`, `"Nemesis"`).
    pub name: String,
    /// Ecliptic longitude in degrees `[0.0, 360.0)`.
    pub lon_deg: f64,
    /// Complete [`ZodiacPlacement`] (Sign, DMS, Decan, Nakshatra/Pada, HD Gate/Line/Color/Tone/Base).
    pub placement: ZodiacPlacement,
    /// House placement (`1..=12`) evaluated via [`house_of_lon`].
    pub house_number: u8,
}

/// Pure mathematical evaluator for all 7 Classical Hermetic Lots of Paulus Alexandrinus
/// from precomputed `asc_deg`, 12 `cusps_deg`, `is_day_chart`, and the 7 classical planetary longitudes.
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn calc_lots_from_positions(
    asc_deg: f64,
    cusps_deg: &[f64; 12],
    is_day_chart: bool,
    sun_lon_deg: f64,
    moon_lon_deg: f64,
    mercury_lon_deg: f64,
    venus_lon_deg: f64,
    mars_lon_deg: f64,
    jupiter_lon_deg: f64,
    saturn_lon_deg: f64,
) -> [HermeticLot; 7] {
    let formula = |day_arc: f64| -> f64 {
        if is_day_chart {
            normalize_lon_360(asc_deg + day_arc)
        } else {
            normalize_lon_360(asc_deg - day_arc)
        }
    };

    // 1. Fortune: Day = ASC + Moon - Sun, Night = ASC + Sun - Moon
    let fortune_lon = formula(moon_lon_deg - sun_lon_deg);
    // 2. Spirit: Day = ASC + Sun - Moon, Night = ASC + Moon - Sun
    let spirit_lon = formula(sun_lon_deg - moon_lon_deg);
    // 3. Eros: Day = ASC + Venus - Spirit, Night = ASC + Spirit - Venus
    let eros_lon = formula(venus_lon_deg - spirit_lon);
    // 4. Necessity: Day = ASC + Fortune - Mercury, Night = ASC + Mercury - Fortune
    let necessity_lon = formula(fortune_lon - mercury_lon_deg);
    // 5. Courage: Day = ASC + Fortune - Mars, Night = ASC + Mars - Fortune
    let courage_lon = formula(fortune_lon - mars_lon_deg);
    // 6. Victory: Day = ASC + Jupiter - Spirit, Night = ASC + Spirit - Jupiter
    let victory_lon = formula(jupiter_lon_deg - spirit_lon);
    // 7. Nemesis: Day = ASC + Fortune - Saturn, Night = ASC + Saturn - Fortune
    let nemesis_lon = formula(fortune_lon - saturn_lon_deg);

    let make_lot = |name: &str, lon_deg: f64| -> HermeticLot {
        HermeticLot {
            name: name.to_owned(),
            lon_deg,
            placement: map_zodiac_placement(lon_deg),
            house_number: house_of_lon(lon_deg, cusps_deg, 1.0e-9),
        }
    };

    [
        make_lot("Fortune", fortune_lon),
        make_lot("Spirit", spirit_lon),
        make_lot("Eros", eros_lon),
        make_lot("Necessity", necessity_lon),
        make_lot("Courage", courage_lon),
        make_lot("Victory", victory_lon),
        make_lot("Nemesis", nemesis_lon),
    ]
}

/// Evaluates all 7 Classical Hermetic Lots of Paulus Alexandrinus (`Fortune`, `Spirit`, `Eros`,
/// `Necessity`, `Courage`, `Victory`, `Nemesis`) at `input` for `observer` and `house_options`.
pub fn calc_lots(
    engine: &Engine,
    input: TimeInput,
    observer: ObserverLocation,
    house_options: HouseOptions,
) -> Result<Vec<HermeticLot>, EphemerisError> {
    let moment = engine.at(input, Some(observer))?;
    let houses = engine.houses_in_moment(&moment, observer, house_options)?;
    let sect = chart_sect(engine, input, observer)?;

    let body_opts = CalcOptions {
        origin: Origin::Geocentric,
        plane: Plane::Ecliptic,
        equinox: Equinox::OfDate,
        nutation: house_options.nutation,
        aberration: AbMode::Apparent,
        zodiac: house_options.zodiac,
        include_speed: false,
    };

    let sun = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Sun), body_opts)?;
    let moon = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Moon), body_opts)?;
    let mercury = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Mercury), body_opts)?;
    let venus = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Venus), body_opts)?;
    let mars = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Mars), body_opts)?;
    let jupiter = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Jupiter), body_opts)?;
    let saturn = engine.calc_in_moment(&moment, BodyId::Planet(PlanetId::Saturn), body_opts)?;

    Ok(calc_lots_from_positions(
        houses.angles.asc_deg,
        &houses.cusps_deg,
        sect.is_day_chart,
        sun.lon_deg,
        moon.lon_deg,
        mercury.lon_deg,
        venus.lon_deg,
        mars.lon_deg,
        jupiter.lon_deg,
        saturn.lon_deg,
    )
    .to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::EngineConfig;
    use crate::types::{HouseSystem, PolarPolicy, JD_J2000};
    use std::collections::HashSet;

    #[test]
    fn test_hd_mandala_gates_64_bijection_and_boundaries() {
        assert_eq!(HD_MANDALA_GATES_64.len(), 64);
        let unique: HashSet<u8> = HD_MANDALA_GATES_64.iter().copied().collect();
        assert_eq!(unique.len(), 64);
        for g in 1..=64u8 {
            assert!(unique.contains(&g), "Missing gate {g} in HD_MANDALA_GATES_64");
        }

        // Exact origin: 302.0° (2°00' Aquarius) = Gate 41, Line 1, Color 1, Tone 1, Base 1
        let p_origin = map_zodiac_placement(302.0);
        assert_eq!(p_origin.sign_index, 10);
        assert_eq!(p_origin.sign_name, "Aquarius");
        assert_eq!(p_origin.dms_deg, 2);
        assert_eq!(p_origin.dms_min, 0);
        assert!(p_origin.dms_sec.abs() < 1.0e-9);
        assert_eq!(p_origin.decan, 1);
        assert_eq!(p_origin.hd_gate, 41);
        assert_eq!(p_origin.hd_line, 1);
        assert_eq!(p_origin.hd_color, 1);
        assert_eq!(p_origin.hd_tone, 1);
        assert_eq!(p_origin.hd_base, 1);

        // Just before 302.0° -> last gate (Gate 60), Line 6, Color 6, Tone 6, Base 5
        let p_end = map_zodiac_placement(302.0 - 1.0e-4);
        assert_eq!(p_end.hd_gate, 60);
        assert_eq!(p_end.hd_line, 6);
        assert_eq!(p_end.hd_color, 6);
        assert_eq!(p_end.hd_tone, 6);
        assert_eq!(p_end.hd_base, 5);

        // Second gate (Gate 19) starts at 302.0 + 5.625 = 307.625° (7°37'30" Aquarius)
        let p_g19 = map_zodiac_placement(307.625);
        assert_eq!(p_g19.hd_gate, 19);
        assert_eq!(p_g19.hd_line, 1);
        assert_eq!(p_g19.hd_color, 1);
        assert_eq!(p_g19.hd_tone, 1);
        assert_eq!(p_g19.hd_base, 1);
        assert_eq!(p_g19.dms_deg, 7);
        assert_eq!(p_g19.dms_min, 37);
        assert!((p_g19.dms_sec - 30.0).abs() < 1.0e-6);

        // Verify all 64 gates and 6 lines step cleanly
        for (gate_idx, &expected_gate) in HD_MANDALA_GATES_64.iter().enumerate() {
            for line in 1..=6u8 {
                let lon = 302.0
                    + (gate_idx as f64) * 5.625
                    + f64::from(line - 1) * 0.9375
                    + 1.0e-4;
                let pl = map_zodiac_placement(lon);
                assert_eq!(pl.hd_gate, expected_gate);
                assert_eq!(pl.hd_line, line);
            }
        }
    }

    #[test]
    fn test_nakshatra_27_and_pada_4_mapping() {
        // 0° Aries -> Ashwini (1), Ketu, Pada 1
        let p0 = map_zodiac_placement(0.0);
        assert_eq!(p0.nakshatra_index, 1);
        assert_eq!(p0.nakshatra_name, "Ashwini");
        assert_eq!(p0.nakshatra_lord, "Ketu");
        assert_eq!(p0.pada, 1);

        // 3°20' Aries -> Ashwini (1), Pada 2
        let p_pada2 = map_zodiac_placement(10.0 / 3.0);
        assert_eq!(p_pada2.nakshatra_index, 1);
        assert_eq!(p_pada2.pada, 2);

        // 13°20' Aries -> Bharani (2), Venus, Pada 1
        let p_bharani = map_zodiac_placement(40.0 / 3.0);
        assert_eq!(p_bharani.nakshatra_index, 2);
        assert_eq!(p_bharani.nakshatra_name, "Bharani");
        assert_eq!(p_bharani.nakshatra_lord, "Venus");
        assert_eq!(p_bharani.pada, 1);

        // 120°00' (0° Leo) -> Magha (10), Ketu, Pada 1
        let p_magha = map_zodiac_placement(120.0);
        assert_eq!(p_magha.sign_index, 4);
        assert_eq!(p_magha.sign_name, "Leo");
        assert_eq!(p_magha.nakshatra_index, 10);
        assert_eq!(p_magha.nakshatra_name, "Magha");
        assert_eq!(p_magha.nakshatra_lord, "Ketu");
        assert_eq!(p_magha.pada, 1);

        // 359.9° -> Revati (27), Mercury, Pada 4
        let p_revati = map_zodiac_placement(359.9);
        assert_eq!(p_revati.nakshatra_index, 27);
        assert_eq!(p_revati.nakshatra_name, "Revati");
        assert_eq!(p_revati.nakshatra_lord, "Mercury");
        assert_eq!(p_revati.pada, 4);
    }

    #[test]
    fn test_parashara_vargas_d1_to_d60() {
        // Reject invalid divisions
        assert!(calc_varga(10.0, 0).is_err());
        assert!(calc_varga(10.0, 61).is_err());

        // D1 (Rashi): 15° Taurus (45.0°) -> Taurus (1), Venus
        let d1 = calc_varga(45.0, 1).expect("D1");
        assert_eq!(d1.sign_index, 1);
        assert_eq!(d1.sign_name, "Taurus");
        assert_eq!(d1.sign_lord, "Venus");

        // D2 (Hora):
        // Odd sign (Aries 0): [0..15) -> Leo (4), [15..30) -> Cancer (3)
        assert_eq!(calc_varga(5.0, 2).expect("D2 odd 1").sign_index, 4);
        assert_eq!(calc_varga(20.0, 2).expect("D2 odd 2").sign_index, 3);
        // Even sign (Taurus 1): [0..15) -> Cancer (3), [15..30) -> Leo (4)
        assert_eq!(calc_varga(35.0, 2).expect("D2 even 1").sign_index, 3);
        assert_eq!(calc_varga(50.0, 2).expect("D2 even 2").sign_index, 4);

        // D3 (Drekkana): Aries 0 -> [Aries(0), Leo(4), Sagittarius(8)]
        assert_eq!(calc_varga(5.0, 3).expect("D3 1").sign_index, 0);
        assert_eq!(calc_varga(15.0, 3).expect("D3 2").sign_index, 4);
        assert_eq!(calc_varga(25.0, 3).expect("D3 3").sign_index, 8);

        // D9 (Navamsha):
        // 0° Aries -> Aries (0); 3°20' Taurus (33.334°) -> 2nd Navamsha of Taurus (starts Capricorn 9 -> Aquarius 10)
        assert_eq!(calc_varga(1.0, 9).expect("D9 Aries 1").sign_index, 0);
        assert_eq!(calc_varga(30.5, 9).expect("D9 Taurus 1").sign_index, 9);
        assert_eq!(calc_varga(34.0, 9).expect("D9 Taurus 2").sign_index, 10);

        // D30 (Trimshamsha) unequal segments:
        // Odd sign (Aries): [0..5)->Aries(0), [5..10)->Aquarius(10), [10..18)->Sagittarius(8), [18..25)->Gemini(2), [25..30)->Libra(6)
        assert_eq!(calc_varga(2.0, 30).expect("D30 odd 1").sign_index, 0);
        assert_eq!(calc_varga(7.0, 30).expect("D30 odd 2").sign_index, 10);
        assert_eq!(calc_varga(14.0, 30).expect("D30 odd 3").sign_index, 8);
        assert_eq!(calc_varga(21.0, 30).expect("D30 odd 4").sign_index, 2);
        assert_eq!(calc_varga(27.0, 30).expect("D30 odd 5").sign_index, 6);

        // Even sign (Taurus): [0..5)->Taurus(1), [5..12)->Virgo(5), [12..20)->Pisces(11), [20..25)->Capricorn(9), [25..30)->Scorpio(7)
        assert_eq!(calc_varga(32.0, 30).expect("D30 even 1").sign_index, 1);
        assert_eq!(calc_varga(38.0, 30).expect("D30 even 2").sign_index, 5);
        assert_eq!(calc_varga(46.0, 30).expect("D30 even 3").sign_index, 11);
        assert_eq!(calc_varga(52.0, 30).expect("D30 even 4").sign_index, 9);
        assert_eq!(calc_varga(58.0, 30).expect("D30 even 5").sign_index, 7);

        // Shodashavarga 16-chart bundle
        let all16 = calc_shodashavarga(123.456).expect("shodashavarga");
        assert_eq!(all16.len(), 16);
    }

    #[test]
    fn test_vimshottari_dasha_120_year_cycle_and_antardashas() {
        // At 0.0° Aries (start of Ashwini), full 7-year Ketu Mahadasha remains, and total 9 Mahadashas = 120 years
        let birth_jd = JD_J2000;
        let report =
            vimshottari_dasha_from_lon(0.0, birth_jd, Some(365.2425)).expect("vimshottari");
        assert_eq!(report.birth_nakshatra_index, 1);
        assert_eq!(report.birth_nakshatra_name, "Ashwini");
        assert_eq!(report.birth_pada, 1);
        assert!((report.balance_fraction - 1.0).abs() < 1.0e-12);
        assert!((report.balance_years - 7.0).abs() < 1.0e-12);
        assert_eq!(report.periods.len(), 9);

        let total_years: f64 = report.periods.iter().map(|p| p.duration_years).sum();
        assert!((total_years - 120.0).abs() < 1.0e-10);
        let last_period = report.periods.last().expect("9th period");
        assert!((last_period.end_jd_ut1 - (birth_jd + 120.0 * 365.2425)).abs() < 1.0e-7);

        for p in &report.periods {
            assert_eq!(p.sub_periods.len(), 9);
            let first_sub = p.sub_periods.first().expect("1st sub");
            let last_sub = p.sub_periods.last().expect("9th sub");
            assert!((first_sub.start_jd_ut1 - p.start_jd_ut1).abs() < 1.0e-10);
            assert!((last_sub.end_jd_ut1 - p.end_jd_ut1).abs() < 1.0e-10);
            let sub_sum: f64 = p.sub_periods.iter().map(|s| s.duration_years).sum();
            assert!((sub_sum - p.duration_years).abs() < 1.0e-10);
        }

        // At midpoint of Ashwini (6°40' = 20/3 deg), 50% (3.5 years) of Ketu Mahadasha remains.
        // The first 4 Antardashas of Ketu (Ketu 7/120*7, Venus 20/120*7, Sun 6/120*7, Moon 10/120*7 = 43/120*7 yr)
        // and part of Mars (7/120*7 yr, cumulative 50/120*7 < 60/120*7) have already elapsed before birth,
        // and birth falls inside Ketu-Rahu (cumulative (43+7)=50/120 to (50+18)=68/120 of Ketu)!
        // So 4 surviving Antardashas remain: Rahu (clipped), Jupiter, Saturn, Mercury!
        let mid_report =
            vimshottari_dasha_from_lon(20.0 / 3.0, birth_jd, Some(365.25)).expect("mid ashwini");
        assert!((mid_report.balance_fraction - 0.5).abs() < 1.0e-12);
        assert!((mid_report.balance_years - 3.5).abs() < 1.0e-12);
        let first_maha = mid_report.periods.first().expect("1st maha");
        assert_eq!(first_maha.lord, "Ketu");
        assert!((first_maha.duration_years - 3.5).abs() < 1.0e-12);
        assert_eq!(first_maha.sub_periods.len(), 4);
        let sp0 = first_maha.sub_periods.first().expect("sp0");
        let sp1 = first_maha.sub_periods.get(1).expect("sp1");
        let sp2 = first_maha.sub_periods.get(2).expect("sp2");
        let sp3 = first_maha.sub_periods.get(3).expect("sp3");
        assert_eq!(sp0.lord, "Rahu");
        assert!((sp0.start_jd_ut1 - birth_jd).abs() < 1.0e-10);
        // Unclipped subsequent Antardashas in the birth Mahadasha have full canonical durations (7 * Y / 120)
        assert_eq!(sp1.lord, "Jupiter");
        assert!((sp1.duration_years - 7.0 * 16.0 / 120.0).abs() < 1.0e-12);
        assert_eq!(sp2.lord, "Saturn");
        assert!((sp2.duration_years - 7.0 * 19.0 / 120.0).abs() < 1.0e-12);
        assert_eq!(sp3.lord, "Mercury");
        assert!((sp3.duration_years - 7.0 * 17.0 / 120.0).abs() < 1.0e-12);
        let sub_sum_mid: f64 = first_maha.sub_periods.iter().map(|s| s.duration_years).sum();
        assert!((sub_sum_mid - 3.5).abs() < 1.0e-10);
    }

    #[test]
    fn test_panchanga_and_muhurtas_new_delhi_and_polar_fallback() {
        let engine = Engine::new_moshier(EngineConfig::default());
        let new_delhi = Observer::new(28.6139, 77.2090, 216.0).expect("New Delhi");
        let input = TimeInput::from_ut1_jd(JD_J2000); // 2000-01-01 12:00 UT1 (17:30 IST Saturday)

        let p = calc_panchanga(&engine, input, new_delhi, AyanamsaModel::Lahiri)
            .expect("calc_panchanga");

        assert!(!p.polar_day_fallback);
        assert!(p.sunrise_jd_ut1 < p.sunset_jd_ut1);
        assert!(p.sunset_jd_ut1 < p.next_sunrise_jd_ut1);
        assert!(p.sunrise_jd_ut1 <= p.jd_ut1 && p.jd_ut1 < p.next_sunrise_jd_ut1);

        // 2000-01-01 in New Delhi was Saturday (Shanivara, 6, Saturn)
        assert_eq!(p.vara.vara_index, 6);
        assert_eq!(p.vara.vara_name, "Shanivara");
        assert_eq!(p.vara.vara_lord, "Saturn");

        // Pre-sunrise check: 2000-01-02 03:00 IST (1999-12-31 21:30 UT1 + 1.0d = JD_J2000 + 0.395833)
        // Before Sunday's sunrise in New Delhi, the active Vedic day is still Saturday (Shanivara, 6)!
        let pre_sunrise_sunday = TimeInput::from_ut1_jd(JD_J2000 + 0.3958333333333333);
        let p_pre = calc_panchanga(&engine, pre_sunrise_sunday, new_delhi, AyanamsaModel::Lahiri)
            .expect("pre-sunrise panchanga");
        assert_eq!(p_pre.vara.vara_index, 6);
        assert_eq!(p_pre.vara.vara_name, "Shanivara");
        assert!((p_pre.sunrise_jd_ut1 - p.sunrise_jd_ut1).abs() < 1.0e-6);

        // Tithi, Nakshatra, Yoga, Karana bounds and timestamps
        assert!((1..=30).contains(&p.tithi.tithi_index));
        assert!(p.tithi.start_jd_ut1 <= p.jd_ut1 && p.jd_ut1 <= p.tithi.end_jd_ut1);
        assert!((1..=27).contains(&p.nakshatra.nakshatra_index));
        assert!((1..=4).contains(&p.nakshatra.pada));
        assert!(p.nakshatra.start_jd_ut1 <= p.jd_ut1 && p.jd_ut1 <= p.nakshatra.end_jd_ut1);
        assert!((1..=27).contains(&p.yoga.yoga_index));
        assert!(p.yoga.start_jd_ut1 <= p.jd_ut1 && p.jd_ut1 <= p.yoga.end_jd_ut1);
        assert!((1..=60).contains(&p.karana.karana_index));
        assert!(p.karana.start_jd_ut1 <= p.jd_ut1 && p.jd_ut1 <= p.karana.end_jd_ut1);

        // Muhurtas within [sunrise, sunset]
        assert!(
            p.rahukalam.start_jd_ut1 >= p.sunrise_jd_ut1 - 1.0e-9
                && p.rahukalam.end_jd_ut1 <= p.sunset_jd_ut1 + 1.0e-9
        );
        assert!(
            p.abhijit_muhurta.start_jd_ut1 > p.sunrise_jd_ut1
                && p.abhijit_muhurta.end_jd_ut1 < p.sunset_jd_ut1
        );

        // Polar fallback at Longyearbyen (78.22° N) during polar night / midnight sun
        let svalbard = Observer::new(78.2232, 15.6267, 0.0).expect("Svalbard");
        let p_polar = calc_panchanga(&engine, input, svalbard, AyanamsaModel::Lahiri)
            .expect("polar panchanga must never fail");
        assert!(p_polar.polar_day_fallback);
        assert!((p_polar.sunset_jd_ut1 - p_polar.sunrise_jd_ut1 - 0.5).abs() < 1.0e-12);
        assert!((p_polar.next_sunrise_jd_ut1 - p_polar.sunrise_jd_ut1 - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_traditional_sect_planetary_hours_and_7_hermetic_lots() {
        let engine = Engine::new_moshier(EngineConfig::default());
        let london = Observer::new(51.5074, -0.1278, 0.0).expect("London");
        let noon_ut1 = TimeInput::from_ut1_jd(JD_J2000); // 2000-01-01 12:00 UT1 (Saturday noon in London)

        let sect = chart_sect(&engine, noon_ut1, london).expect("chart_sect");
        assert!(sect.is_day_chart);
        assert!(sect.sun_altitude_true_deg > 10.0);

        let ph = planetary_hours(&engine, noon_ut1, london).expect("planetary_hours");
        assert_eq!(ph.day_lord, "Saturn"); // Saturday
        assert_eq!(ph.hours.len(), 24);
        assert_eq!(ph.hours.first().expect("h1").lord, "Saturn"); // 1st hour of Saturday = Saturn
        assert_eq!(ph.hours.get(1).expect("h2").lord, "Jupiter"); // Chaldean sequence
        assert!(ph.is_daytime_hour);
        assert!((1..=12).contains(&ph.current_hour_index));

        // 25th hour (next day 1st hour) would be index (0 + 24) % 7 = 3 => Sun (Sunday)!
        assert_eq!(ph.hours.last().expect("h24").lord, "Mars");

        let h_opts = HouseOptions {
            system: HouseSystem::WholeSign,
            polar_policy: PolarPolicy::FallbackToWholeSign,
            zodiac: Zodiac::Tropical,
            nutation: NutationMode::True,
            include_speed: false,
        };
        let lots = calc_lots(&engine, noon_ut1, london, h_opts).expect("calc_lots");
        assert_eq!(lots.len(), 7);
        let expected_names = [
            "Fortune",
            "Spirit",
            "Eros",
            "Necessity",
            "Courage",
            "Victory",
            "Nemesis",
        ];
        for (lot, &exp_name) in lots.iter().zip(expected_names.iter()) {
            assert_eq!(lot.name, exp_name);
            assert!((0.0..360.0).contains(&lot.lon_deg));
            assert!((1..=12).contains(&lot.house_number));
        }

        // Verify diurnal vs nocturnal reversal between Fortune and Spirit
        let cusps = [
            0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0,
        ];
        let day_lots =
            calc_lots_from_positions(10.0, &cusps, true, 100.0, 140.0, 90.0, 80.0, 200.0, 30.0, 50.0);
        let night_lots = calc_lots_from_positions(
            10.0, &cusps, false, 100.0, 140.0, 90.0, 80.0, 200.0, 30.0, 50.0,
        );
        // Day Fortune (10 + 140 - 100 = 50°) == Night Spirit
        assert!((day_lots[0].lon_deg - 50.0).abs() < 1.0e-12);
        assert!((day_lots[0].lon_deg - night_lots[1].lon_deg).abs() < 1.0e-12);
        // Day Spirit (10 + 100 - 140 = 330°) == Night Fortune
        assert!((day_lots[1].lon_deg - 330.0).abs() < 1.0e-12);
        assert!((day_lots[1].lon_deg - night_lots[0].lon_deg).abs() < 1.0e-12);

        // Canonical D60 (BPHS Ch. 6 vv. 33–35): ignores natal Rashi index and counts from Mesha (Aries).
        // 30.1° (0°06' Taurus, k=0) -> Aries (0, Mars); 30.6° (0°36' Taurus, k=1) -> Taurus (1, Venus).
        let d60_tau_0 = calc_varga(30.1, 60).expect("D60 30.1");
        assert_eq!(d60_tau_0.sign_index, 0);
        assert_eq!(d60_tau_0.sign_name, "Aries");
        assert_eq!(d60_tau_0.sign_lord, "Mars");
        let d60_tau_1 = calc_varga(30.6, 60).expect("D60 30.6");
        assert_eq!(d60_tau_1.sign_index, 1);
        assert_eq!(d60_tau_1.sign_name, "Taurus");
        assert_eq!(d60_tau_1.sign_lord, "Venus");

        // Polar Ascendant & Hermetic Lots at Longyearbyen (78.2232°N, 15.6469°E) on 2026-10-25 15:00 UT1:
        // Must select the eastern rising horizon intersection (~355.61°, (ASC - MC) mod 360 < 180°),
        // NOT the western setting point (~175.61°), so Fortune (~182.63°) and Spirit (~168.01°) are not flipped 180°.
        let longyearbyen = Observer::new(78.2232, 15.6469, 0.0).expect("Longyearbyen");
        let jd_polar = crate::time::julday(2026, 10, 25, 15.0, crate::CalendarKind::Gregorian)
            .expect("julday");
        let t_polar = TimeInput::from_ut1_jd(jd_polar);
        let h_polar = engine.houses(t_polar, longyearbyen, h_opts).expect("polar houses");
        assert!(
            (h_polar.angles.asc_deg - h_polar.angles.mc_deg).rem_euclid(360.0) < 180.0,
            "Polar ASC must be in the eastern rising hemisphere: asc={}, mc={}",
            h_polar.angles.asc_deg,
            h_polar.angles.mc_deg
        );
        assert!((h_polar.angles.asc_deg - 355.321421).abs() < 0.01);
        let lots_polar = calc_lots(&engine, t_polar, longyearbyen, h_opts).expect("polar lots");
        assert!((lots_polar.first().expect("fortune").lon_deg - 182.634).abs() < 0.05);
        assert!((lots_polar.get(1).expect("spirit").lon_deg - 168.009).abs() < 0.05);
    }
}
