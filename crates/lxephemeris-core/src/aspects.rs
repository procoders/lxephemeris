//! High-precision astrological aspect, mutual/natal transit, declination parallel, and
//! Void-of-Course Moon solvers (`Wave H1`):
//! - Aspect presets (`AspectPreset::Major`, `Standard`, `Kepler`, `All`) and `AspectDefinition`.
//! - Static ecliptic and equatorial declination aspect evaluators (`calc_aspects`, `synastry_aspects`,
//!   `declination_aspects`, `calc_chart_aspects`, `calc_synastry_aspects`, `calc_declination_aspects`)
//!   with analytical applying/separating detection ($e \cdot \dot{e} < 0$).
//! - Dynamic exact mutual and natal transit aspect root-finders (`search_exact_aspect`,
//!   `mutual_aspects_in_range`, `natal_transits_in_range`) using hybrid Newton–Illinois polishing
//!   and relative stationary point ($v_{\text{rel}}(t) = 0$) interval splitting so grazing
//!   retrograde double/triple passes are never missed.
//! - Void-of-Course Moon interval solver (`void_of_course_moon_in_range`, `VoidOfCourseInterval`).

use crate::astrometry::CalcResult;
use crate::engine::Engine;
use crate::error::EphemerisError;
use crate::events::wrap_pm180;
use crate::time::TimeInput;
use crate::types::{
    BodyId, CalcOptions, Equinox, LunarPointId, Origin, Plane, PlanetId, Zodiac,
};
use serde::{Deserialize, Serialize};

// ============================================================================
// 1. Aspect Definitions & Presets
// ============================================================================

/// Definition of an astrological aspect angle, default orb, and major/minor classification.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AspectDefinition {
    /// Canonical aspect name (e.g., `"Conjunction"`, `"Sextile"`, `"Square"`, `"Trine"`, `"Opposition"`).
    pub name: String,
    /// Exact aspect separation angle in degrees $\in [0^\circ, 180^\circ]$.
    pub angle_deg: f64,
    /// Default maximum orb in degrees (before `orb_scale`).
    pub default_orb_deg: f64,
    /// `true` for the 5 Ptolemaic major aspects (`0°`, `60°`, `90°`, `120°`, `180°`).
    pub is_major: bool,
}

/// Standard aspect set presets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AspectPreset {
    /// Ptolemaic major aspects: `0° Conjunction`, `60° Sextile`, `90° Square`, `120° Trine`, `180° Opposition`.
    Major,
    /// `Major` + 4 minor aspects: `30° Semisextile`, `45° Semisquare`, `135° Sesquiquadrate`, `150° Quincunx`.
    Standard,
    /// `Standard` + Kepler 5th-harmonic aspects: `72° Quintile`, `144° Biquintile`.
    Kepler,
    /// Alias for the complete 11-aspect set (`Standard` + `72° Quintile`, `144° Biquintile`).
    All,
}

/// Default 9 non-lunar major planets used for Void-of-Course Moon calculations (`Sun..=Pluto` excluding `Moon`).
pub const DEFAULT_VOC_MAJOR_PLANETS: [BodyId; 9] = [
    BodyId::Planet(PlanetId::Sun),
    BodyId::Planet(PlanetId::Mercury),
    BodyId::Planet(PlanetId::Venus),
    BodyId::Planet(PlanetId::Mars),
    BodyId::Planet(PlanetId::Jupiter),
    BodyId::Planet(PlanetId::Saturn),
    BodyId::Planet(PlanetId::Uranus),
    BodyId::Planet(PlanetId::Neptune),
    BodyId::Planet(PlanetId::Pluto),
];

/// Returns the canonical list of [`AspectDefinition`]s for `preset`, sorted ascending by `angle_deg`.
#[must_use]
pub fn aspect_definitions(preset: AspectPreset) -> Vec<AspectDefinition> {
    let mut defs = vec![
        AspectDefinition {
            name: "Conjunction".to_owned(),
            angle_deg: 0.0,
            default_orb_deg: 8.0,
            is_major: true,
        },
        AspectDefinition {
            name: "Sextile".to_owned(),
            angle_deg: 60.0,
            default_orb_deg: 6.0,
            is_major: true,
        },
        AspectDefinition {
            name: "Square".to_owned(),
            angle_deg: 90.0,
            default_orb_deg: 8.0,
            is_major: true,
        },
        AspectDefinition {
            name: "Trine".to_owned(),
            angle_deg: 120.0,
            default_orb_deg: 8.0,
            is_major: true,
        },
        AspectDefinition {
            name: "Opposition".to_owned(),
            angle_deg: 180.0,
            default_orb_deg: 8.0,
            is_major: true,
        },
    ];

    if matches!(
        preset,
        AspectPreset::Standard | AspectPreset::Kepler | AspectPreset::All
    ) {
        defs.extend([
            AspectDefinition {
                name: "Semisextile".to_owned(),
                angle_deg: 30.0,
                default_orb_deg: 2.0,
                is_major: false,
            },
            AspectDefinition {
                name: "Semisquare".to_owned(),
                angle_deg: 45.0,
                default_orb_deg: 2.0,
                is_major: false,
            },
            AspectDefinition {
                name: "Sesquiquadrate".to_owned(),
                angle_deg: 135.0,
                default_orb_deg: 2.0,
                is_major: false,
            },
            AspectDefinition {
                name: "Quincunx".to_owned(),
                angle_deg: 150.0,
                default_orb_deg: 3.0,
                is_major: false,
            },
        ]);
    }

    if matches!(preset, AspectPreset::Kepler | AspectPreset::All) {
        defs.extend([
            AspectDefinition {
                name: "Quintile".to_owned(),
                angle_deg: 72.0,
                default_orb_deg: 2.0,
                is_major: false,
            },
            AspectDefinition {
                name: "Biquintile".to_owned(),
                angle_deg: 144.0,
                default_orb_deg: 2.0,
                is_major: false,
            },
        ]);
    }

    defs.sort_by(|a, b| a.angle_deg.total_cmp(&b.angle_deg));
    defs
}

/// Maps a target angle in degrees to its canonical aspect name (or formats custom angles).
#[must_use]
pub fn canonical_aspect_name(angle_deg: f64) -> String {
    let norm = angle_deg.rem_euclid(360.0);
    let sym = if norm > 180.0 { 360.0 - norm } else { norm };
    for def in aspect_definitions(AspectPreset::All) {
        if (sym - def.angle_deg).abs() < 1.0e-6 {
            return def.name;
        }
    }
    format!("{sym:.4}°")
}

// ============================================================================
// 2. Static Aspect Calculation (Chart, Synastry, Declination)
// ============================================================================

/// Sampled ecliptic (and optional equatorial declination) state of a body for static aspect evaluation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BodyPositionSample {
    /// Body identifier.
    pub body: BodyId,
    /// Ecliptic longitude $\lambda$ in degrees `[0, 360)`.
    pub lon_deg: f64,
    /// Ecliptic latitude $\beta$ in degrees `[-90, +90]`.
    pub lat_deg: f64,
    /// Ecliptic longitudinal velocity $\dot{\lambda}$ in degrees/day.
    pub lon_speed_deg_day: f64,
    /// Optional equatorial declination $\delta$ in degrees `[-90, +90]`.
    pub dec_deg: Option<f64>,
    /// Optional equatorial declination velocity $\dot{\delta}$ in degrees/day.
    pub dec_speed_deg_day: Option<f64>,
}

/// Detected static ecliptic aspect between `body1` and `body2`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AspectHit {
    /// First body identifier.
    pub body1: BodyId,
    /// Second body identifier.
    pub body2: BodyId,
    /// Canonical aspect name (`"Conjunction"`, `"Sextile"`, `"Square"`, `"Trine"`, `"Opposition"`, etc.).
    pub aspect_name: String,
    /// Exact aspect angle $A \in [0^\circ, 180^\circ]$.
    pub exact_angle_deg: f64,
    /// Actual unsigned angular separation $|\text{wrap}_{[-180^\circ, 180^\circ)}(\lambda_1 - \lambda_2)| \in [0^\circ, 180^\circ]$.
    pub actual_separation_deg: f64,
    /// Absolute orb $|\text{actual\_separation\_deg} - A|$ in degrees.
    pub orb_deg: f64,
    /// Maximum allowed orb (`default_orb_deg * orb_scale`) in degrees.
    pub max_orb_deg: f64,
    /// Aspect tightness percentage `((1.0 - orb_deg.abs() / max_orb_deg).clamp(0.0, 1.0)) * 100.0`.
    pub tightness_pct: f64,
    /// `true` if the angular error $e$ and its time derivative $\dot{e}$ satisfy $e \cdot \dot{e} < 0$.
    pub is_applying: bool,
    /// `true` if this is one of the 5 Ptolemaic major aspects.
    pub is_major: bool,
    /// Relative longitudinal speed $\dot{\lambda}_1 - \dot{\lambda}_2$ in degrees/day.
    pub relative_speed_deg_day: f64,
}

/// Evaluates whether `(s1, s2)` forms an aspect from `defs` (scaled by `orb_scale`).
/// If multiple definitions match, returns the one with the smallest `orb_deg`.
fn eval_pair_aspect(
    s1: &BodyPositionSample,
    s2: &BodyPositionSample,
    defs: &[AspectDefinition],
    orb_scale: f64,
) -> Option<AspectHit> {
    let scale = if orb_scale.is_finite() && orb_scale > 0.0 {
        orb_scale
    } else {
        1.0
    };
    let delta_lon = wrap_pm180(s1.lon_deg - s2.lon_deg);
    let actual_separation_deg = delta_lon.abs();
    let v_rel = s1.lon_speed_deg_day - s2.lon_speed_deg_day;

    let mut best_hit: Option<AspectHit> = None;

    for def in defs {
        let a = def.angle_deg;
        let max_orb_deg = def.default_orb_deg * scale;

        let e = if a <= 1.0e-9 {
            delta_lon
        } else if (a - 180.0).abs() <= 1.0e-9 {
            wrap_pm180(delta_lon - 180.0)
        } else {
            let sign = if delta_lon >= 0.0 { 1.0 } else { -1.0 };
            delta_lon - sign * a
        };

        let orb_deg = e.abs();
        if orb_deg <= max_orb_deg {
            let is_applying = (e * v_rel) < 0.0;
            let tightness_pct = if max_orb_deg > 0.0 {
                ((1.0 - orb_deg.abs() / max_orb_deg).clamp(0.0, 1.0)) * 100.0
            } else {
                100.0
            };

            let candidate = AspectHit {
                body1: s1.body,
                body2: s2.body,
                aspect_name: def.name.clone(),
                exact_angle_deg: a,
                actual_separation_deg,
                orb_deg,
                max_orb_deg,
                tightness_pct,
                is_applying,
                is_major: def.is_major,
                relative_speed_deg_day: v_rel,
            };

            let replace = match &best_hit {
                None => true,
                Some(prev) => candidate.orb_deg < prev.orb_deg,
            };
            if replace {
                best_hit = Some(candidate);
            }
        }
    }

    best_hit
}

/// Computes all static mutual aspects among unique pairs `(i < j)` in `samples`.
#[must_use]
pub fn calc_aspects(
    samples: &[BodyPositionSample],
    defs: &[AspectDefinition],
    orb_scale: f64,
) -> Vec<AspectHit> {
    let mut hits = Vec::new();
    for (i, s1) in samples.iter().enumerate() {
        for s2 in samples.iter().skip(i + 1) {
            if s1.body == s2.body {
                continue;
            }
            if let Some(hit) = eval_pair_aspect(s1, s2, defs, orb_scale) {
                hits.push(hit);
            }
        }
    }
    hits
}

/// Computes all cross-chart (synastry / transit-to-natal) aspects between `samples1` and `samples2`.
#[must_use]
pub fn synastry_aspects(
    samples1: &[BodyPositionSample],
    samples2: &[BodyPositionSample],
    defs: &[AspectDefinition],
    orb_scale: f64,
) -> Vec<AspectHit> {
    let mut hits = Vec::new();
    for s1 in samples1 {
        for s2 in samples2 {
            if let Some(hit) = eval_pair_aspect(s1, s2, defs, orb_scale) {
                hits.push(hit);
            }
        }
    }
    hits
}

/// Evaluates all `bodies` at `input` and computes all chart aspects for `preset` and `orb_scale`.
pub fn calc_chart_aspects(
    engine: &Engine,
    input: TimeInput,
    bodies: &[BodyId],
    preset: AspectPreset,
    orb_scale: f64,
    options: CalcOptions,
) -> Result<Vec<AspectHit>, EphemerisError> {
    if !orb_scale.is_finite() || orb_scale <= 0.0 {
        return Err(EphemerisError::InvalidInput {
            message: format!("Invalid orb_scale: {orb_scale} (must be finite and > 0)"),
        });
    }
    let mut eval_opts = options;
    eval_opts.plane = Plane::Ecliptic;
    eval_opts.include_speed = true;

    let obs = match eval_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let moment = engine.at(input, obs)?;
    let mut samples = Vec::with_capacity(bodies.len());
    for &body in bodies {
        let res = engine.calc_in_moment(&moment, body, eval_opts)?;
        samples.push(BodyPositionSample {
            body,
            lon_deg: res.lon_deg,
            lat_deg: res.lat_deg,
            lon_speed_deg_day: res.lon_speed_deg_day,
            dec_deg: None,
            dec_speed_deg_day: None,
        });
    }

    let defs = aspect_definitions(preset);
    Ok(calc_aspects(&samples, &defs, orb_scale))
}

/// Evaluates `bodies1` at `input1` (moving/transiting chart) and `bodies2` at `input2`
/// (reference/natal chart, with speeds set to `0.0` so `relative_speed_deg_day` reflects
/// `bodies1`'s motion against fixed natal positions) and returns all cross-chart [`AspectHit`]s.
#[allow(clippy::too_many_arguments)]
pub fn calc_synastry_aspects(
    engine: &Engine,
    input1: TimeInput,
    bodies1: &[BodyId],
    input2: TimeInput,
    bodies2: &[BodyId],
    preset: AspectPreset,
    orb_scale: f64,
    options: CalcOptions,
) -> Result<Vec<AspectHit>, EphemerisError> {
    if !orb_scale.is_finite() || orb_scale <= 0.0 {
        return Err(EphemerisError::InvalidInput {
            message: format!("Invalid orb_scale: {orb_scale} (must be finite and > 0)"),
        });
    }
    let mut eval_opts = options;
    eval_opts.plane = Plane::Ecliptic;
    eval_opts.include_speed = true;

    let obs = match eval_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let m1 = engine.at(input1, obs)?;
    let m2 = engine.at(input2, obs)?;

    let mut samples1 = Vec::with_capacity(bodies1.len());
    for &body in bodies1 {
        let res = engine.calc_in_moment(&m1, body, eval_opts)?;
        samples1.push(BodyPositionSample {
            body,
            lon_deg: res.lon_deg,
            lat_deg: res.lat_deg,
            lon_speed_deg_day: res.lon_speed_deg_day,
            dec_deg: None,
            dec_speed_deg_day: None,
        });
    }

    let mut samples2 = Vec::with_capacity(bodies2.len());
    for &body in bodies2 {
        let res = engine.calc_in_moment(&m2, body, eval_opts)?;
        samples2.push(BodyPositionSample {
            body,
            lon_deg: res.lon_deg,
            lat_deg: res.lat_deg,
            lon_speed_deg_day: 0.0,
            dec_deg: None,
            dec_speed_deg_day: None,
        });
    }

    let defs = aspect_definitions(preset);
    Ok(synastry_aspects(&samples1, &samples2, &defs, orb_scale))
}

/// Equatorial declination aspect kind (`Parallel` or `Contraparallel`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DeclinationAspectKind {
    /// Parallel of declination ($|\delta_1 - \delta_2| \le \text{orb}$).
    Parallel,
    /// Contraparallel of declination ($|\delta_1 + \delta_2| \le \text{orb}$).
    Contraparallel,
}

/// Detected equatorial declination aspect (`Parallel` or `Contraparallel`) between `body1` and `body2`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DeclinationAspectHit {
    /// First body identifier.
    pub body1: BodyId,
    /// Second body identifier.
    pub body2: BodyId,
    /// Canonical aspect name (`"Parallel"` or `"Contraparallel"`).
    pub aspect_name: String,
    /// Declination aspect classification (`Parallel` or `Contraparallel`).
    pub kind: DeclinationAspectKind,
    /// Declination $\delta_1$ of `body1` in degrees `[-90, +90]`.
    pub dec1_deg: f64,
    /// Declination $\delta_2$ of `body2` in degrees `[-90, +90]`.
    pub dec2_deg: f64,
    /// Unsigned declination orb in degrees ($|\delta_1 - \delta_2|$ for `Parallel`, $|\delta_1 + \delta_2|$ for `Contraparallel`).
    pub orb_deg: f64,
    /// Maximum allowed declination orb in degrees.
    pub max_orb_deg: f64,
    /// Tightness percentage `((1.0 - orb_deg.abs() / max_orb_deg).clamp(0.0, 1.0)) * 100.0`.
    pub tightness_pct: f64,
    /// `true` if the declination separation is decreasing with time ($e \cdot \dot{e} < 0$).
    pub is_applying: bool,
    /// Relative declination speed in degrees/day ($\dot{\delta}_1 - \dot{\delta}_2$ for `Parallel`,
    /// $\dot{\delta}_1 + \dot{\delta}_2$ for `Contraparallel`).
    pub relative_speed_deg_day: f64,
}

/// Computes all `Parallel` ($|\delta_1 - \delta_2| \le \text{max\_orb\_deg}$) and
/// `Contraparallel` ($|\delta_1 + \delta_2| \le \text{max\_orb\_deg}$) declination aspects among
/// unique pairs `(i < j)` in `samples`.
#[must_use]
pub fn declination_aspects(
    samples: &[BodyPositionSample],
    max_orb_deg: f64,
) -> Vec<DeclinationAspectHit> {
    if !max_orb_deg.is_finite() || max_orb_deg <= 0.0 {
        return Vec::new();
    }
    let mut hits = Vec::new();

    for (i, s1) in samples.iter().enumerate() {
        for s2 in samples.iter().skip(i + 1) {
            if s1.body == s2.body {
                continue;
            }
            let (Some(d1), Some(d2)) = (s1.dec_deg, s2.dec_deg) else {
                continue;
            };
            let v1 = s1.dec_speed_deg_day.unwrap_or(0.0);
            let v2 = s2.dec_speed_deg_day.unwrap_or(0.0);

            // Parallel: e = d1 - d2, de/dt = v1 - v2
            let e_par = d1 - d2;
            let orb_par = e_par.abs();
            if orb_par <= max_orb_deg {
                let v_par = v1 - v2;
                let tightness_pct = ((1.0 - orb_par / max_orb_deg).clamp(0.0, 1.0)) * 100.0;
                hits.push(DeclinationAspectHit {
                    body1: s1.body,
                    body2: s2.body,
                    aspect_name: "Parallel".to_owned(),
                    kind: DeclinationAspectKind::Parallel,
                    dec1_deg: d1,
                    dec2_deg: d2,
                    orb_deg: orb_par,
                    max_orb_deg,
                    tightness_pct,
                    is_applying: (e_par * v_par) < 0.0,
                    relative_speed_deg_day: v_par,
                });
            }

            // Contraparallel: e = d1 + d2, de/dt = v1 + v2
            let e_contra = d1 + d2;
            let orb_contra = e_contra.abs();
            if orb_contra <= max_orb_deg {
                let v_contra = v1 + v2;
                let tightness_pct = ((1.0 - orb_contra / max_orb_deg).clamp(0.0, 1.0)) * 100.0;
                hits.push(DeclinationAspectHit {
                    body1: s1.body,
                    body2: s2.body,
                    aspect_name: "Contraparallel".to_owned(),
                    kind: DeclinationAspectKind::Contraparallel,
                    dec1_deg: d1,
                    dec2_deg: d2,
                    orb_deg: orb_contra,
                    max_orb_deg,
                    tightness_pct,
                    is_applying: (e_contra * v_contra) < 0.0,
                    relative_speed_deg_day: v_contra,
                });
            }
        }
    }

    hits
}

/// Evaluates equatorial declination $\delta$ and speed $\dot{\delta}$ for `bodies` at `input`
/// and returns all `Parallel` and `Contraparallel` [`DeclinationAspectHit`]s within `max_orb_deg`.
pub fn calc_declination_aspects(
    engine: &Engine,
    input: TimeInput,
    bodies: &[BodyId],
    max_orb_deg: f64,
    options: CalcOptions,
) -> Result<Vec<DeclinationAspectHit>, EphemerisError> {
    if !max_orb_deg.is_finite() || max_orb_deg <= 0.0 {
        return Err(EphemerisError::InvalidInput {
            message: format!("Invalid max_orb_deg: {max_orb_deg} (must be finite and > 0)"),
        });
    }

    let mut eq_opts = options;
    eq_opts.plane = Plane::Equatorial;
    eq_opts.equinox = Equinox::OfDate;
    eq_opts.zodiac = Zodiac::Tropical;
    eq_opts.include_speed = true;

    let obs = match eq_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let moment = engine.at(input, obs)?;
    let mut samples = Vec::with_capacity(bodies.len());

    for &body in bodies {
        let eq_res = engine.calc_in_moment(&moment, body, eq_opts)?;
        samples.push(BodyPositionSample {
            body,
            lon_deg: eq_res.lon_deg,
            lat_deg: eq_res.lat_deg,
            lon_speed_deg_day: eq_res.lon_speed_deg_day,
            dec_deg: Some(eq_res.lat_deg),
            dec_speed_deg_day: Some(eq_res.lat_speed_deg_day),
        });
    }

    Ok(declination_aspects(&samples, max_orb_deg))
}

// ============================================================================
// 3. Dynamic Solver for Exact Mutual & Transit Aspects
// ============================================================================

/// Exact aspect event between two moving bodies (or a transiting body and a natal body).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExactAspectEvent {
    /// First (or transiting) body identifier.
    pub body1: BodyId,
    /// Second (or natal) body identifier.
    pub body2: BodyId,
    /// Canonical aspect name (`"Conjunction"`, `"Sextile"`, `"Square"`, `"Trine"`, `"Opposition"`, etc.).
    pub aspect_name: String,
    /// Canonical aspect angle in degrees (e.g., `0.0`, `60.0`, `90.0`, `120.0`, `180.0`).
    pub target_angle_deg: f64,
    /// Universal Time (`UT1`) Julian Date of exact aspect.
    pub jd_ut1: f64,
    /// Terrestrial Time (`TT`) Julian Date of exact aspect.
    pub jd_tt: f64,
    /// Ecliptic longitude of `body1` at the exact aspect in degrees `[0, 360)`.
    pub body1_lon_deg: f64,
    /// Ecliptic longitude of `body2` at the exact aspect in degrees `[0, 360)`.
    pub body2_lon_deg: f64,
    /// Relative longitudinal velocity $\dot{\lambda}_1 - \dot{\lambda}_2$ in degrees/day at the exact aspect.
    pub relative_speed_deg_day: f64,
    /// `true` if `body1` is in apparent retrograde motion ($\dot{\lambda}_1 < 0$) at the exact aspect.
    pub body1_retrograde: bool,
    /// `true` if `body2` is in apparent retrograde motion ($\dot{\lambda}_2 < 0$) at the exact aspect.
    pub body2_retrograde: bool,
}

#[derive(Debug, Clone, Copy)]
struct PairEvalState {
    jd_tt: f64,
    jd_ut1: f64,
    lon1_deg: f64,
    lon2_deg: f64,
    speed1_deg_day: f64,
    speed2_deg_day: f64,
}

impl PairEvalState {
    #[inline]
    fn rel_lon_deg(self) -> f64 {
        (self.lon1_deg - self.lon2_deg).rem_euclid(360.0)
    }

    #[inline]
    fn rel_speed_deg_day(self) -> f64 {
        self.speed1_deg_day - self.speed2_deg_day
    }

    #[inline]
    fn error_for_oriented_target(self, oriented_target_deg: f64) -> f64 {
        wrap_pm180(self.lon1_deg - self.lon2_deg - oriented_target_deg)
    }

    fn to_event(
        self,
        body1: BodyId,
        body2: BodyId,
        aspect_name: &str,
        target_angle_deg: f64,
    ) -> ExactAspectEvent {
        ExactAspectEvent {
            body1,
            body2,
            aspect_name: aspect_name.to_owned(),
            target_angle_deg,
            jd_ut1: self.jd_ut1,
            jd_tt: self.jd_tt,
            body1_lon_deg: self.lon1_deg,
            body2_lon_deg: self.lon2_deg,
            relative_speed_deg_day: self.rel_speed_deg_day(),
            body1_retrograde: self.speed1_deg_day < 0.0,
            body2_retrograde: self.speed2_deg_day < 0.0,
        }
    }
}

/// Evaluates two moving bodies at `jd_tt` using a single shared [` crate::moment::Moment`].
fn eval_moving_pair_at_tt(
    engine: &Engine,
    body1: BodyId,
    body2: BodyId,
    jd_tt: f64,
    eval_opts: CalcOptions,
) -> Result<PairEvalState, EphemerisError> {
    let obs = match eval_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let moment = engine.at(TimeInput::from_tt_jd(jd_tt), obs)?;
    let r1 = engine.calc_in_moment(&moment, body1, eval_opts)?;
    let r2 = engine.calc_in_moment(&moment, body2, eval_opts)?;
    Ok(PairEvalState {
        jd_tt: r1.jd_tt,
        jd_ut1: r1.jd_ut1,
        lon1_deg: r1.lon_deg,
        lon2_deg: r2.lon_deg,
        speed1_deg_day: r1.lon_speed_deg_day,
        speed2_deg_day: r2.lon_speed_deg_day,
    })
}

/// Evaluates one transiting body against a fixed natal position `(natal_lon_deg, natal_speed_deg_day)` at `jd_tt`.
fn eval_transit_pair_at_tt(
    engine: &Engine,
    transiting_body: BodyId,
    natal_lon_deg: f64,
    natal_retrograde: bool,
    jd_tt: f64,
    eval_opts: CalcOptions,
) -> Result<PairEvalState, EphemerisError> {
    let r1 = engine.calc(transiting_body, TimeInput::from_tt_jd(jd_tt), eval_opts)?;
    Ok(PairEvalState {
        jd_tt: r1.jd_tt,
        jd_ut1: r1.jd_ut1,
        lon1_deg: r1.lon_deg,
        lon2_deg: natal_lon_deg,
        speed1_deg_day: r1.lon_speed_deg_day,
        speed2_deg_day: if natal_retrograde { -0.0 } else { 0.0 },
    })
}

/// Refines a bracketed root of $f(t) = \text{wrap}_{[-180^\circ, 180^\circ)}(\lambda_1(t) - \lambda_2(t) - \Delta\theta) = 0$
/// on `[state_a.jd_tt, state_b.jd_tt]` using a safeguarded hybrid **Newton–Illinois** solver.
#[allow(clippy::too_many_arguments)]
fn refine_exact_aspect_root<F>(
    eval_at_tt: &F,
    body1: BodyId,
    body2: BodyId,
    aspect_name: &str,
    target_angle_deg: f64,
    oriented_target_deg: f64,
    state_a: PairEvalState,
    state_b: PairEvalState,
) -> Result<ExactAspectEvent, EphemerisError>
where
    F: Fn(f64) -> Result<PairEvalState, EphemerisError>,
{
    let mut a = state_a.jd_tt;
    let mut b = state_b.jd_tt;
    let mut fa = state_a.error_for_oriented_target(oriented_target_deg);
    let mut fb = state_b.error_for_oriented_target(oriented_target_deg);
    let mut va = state_a.rel_speed_deg_day();
    let mut vb = state_b.rel_speed_deg_day();

    if fa.abs() < 1.0e-9 {
        return Ok(state_a.to_event(body1, body2, aspect_name, target_angle_deg));
    }
    if fb.abs() < 1.0e-9 {
        return Ok(state_b.to_event(body1, body2, aspect_name, target_angle_deg));
    }

    let mut best_state = if fa.abs() <= fb.abs() {
        state_a
    } else {
        state_b
    };
    let mut side = 0i8;

    for _ in 0..50 {
        let span = (b - a).abs();
        if span < 1.0e-10 {
            return Ok(best_state.to_event(body1, body2, aspect_name, target_angle_deg));
        }

        let (x_curr, f_curr, v_curr) = if fa.abs() <= fb.abs() {
            (a, fa, va)
        } else {
            (b, fb, vb)
        };
        let t_min = a.min(b);
        let t_max = a.max(b);
        let pad = 0.02 * span;

        let c_newton = if v_curr.abs() > 1.0e-5 {
            Some(x_curr - f_curr / v_curr)
        } else {
            None
        };

        let c = match c_newton {
            Some(cn) if cn > t_min + pad && cn < t_max - pad && (cn - x_curr).abs() < 0.75 * span => {
                cn
            }
            _ => {
                let denom = fb - fa;
                let cf = if denom.abs() > 1.0e-15 {
                    (a * fb - b * fa) / denom
                } else {
                    0.5 * (a + b)
                };
                cf.clamp(t_min + 0.01 * span, t_max - 0.01 * span)
            }
        };

        let state_c = eval_at_tt(c)?;
        let fc = state_c.error_for_oriented_target(oriented_target_deg);
        let vc = state_c.rel_speed_deg_day();

        if fc.abs() < fa.abs().min(fb.abs()) {
            best_state = state_c;
        }
        if fc.abs() < 1.0e-9 {
            return Ok(state_c.to_event(body1, body2, aspect_name, target_angle_deg));
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

    if best_state
        .error_for_oriented_target(oriented_target_deg)
        .abs()
        < 1.0e-6
    {
        return Ok(best_state.to_event(body1, body2, aspect_name, target_angle_deg));
    }

    Err(EphemerisError::NonConvergence {
        solver: "refine_exact_aspect_root".to_owned(),
        iterations: 50,
        residual: (b - a).abs(),
    })
}

/// Solves for the point of relative station $v_{\text{rel}}(t_{\text{stat}}) = \dot{\lambda}_1(t) - \dot{\lambda}_2(t) = 0$
/// inside `[state_a.jd_tt, state_b.jd_tt]` when `va * vb <= 0` using safeguarded Illinois / bisection.
fn solve_relative_station<F>(
    eval_at_tt: &F,
    state_a: PairEvalState,
    state_b: PairEvalState,
) -> Result<PairEvalState, EphemerisError>
where
    F: Fn(f64) -> Result<PairEvalState, EphemerisError>,
{
    let mut a = state_a.jd_tt;
    let mut b = state_b.jd_tt;
    let mut va = state_a.rel_speed_deg_day();
    let mut vb = state_b.rel_speed_deg_day();

    if va.abs() < 1.0e-9 {
        return Ok(state_a);
    }
    if vb.abs() < 1.0e-9 {
        return Ok(state_b);
    }

    let mut best_state = if va.abs() <= vb.abs() {
        state_a
    } else {
        state_b
    };
    let mut side = 0i8;

    for _ in 0..35 {
        let span = (b - a).abs();
        if span < 1.0e-7 || va.abs().min(vb.abs()) < 1.0e-8 {
            break;
        }
        let t_min = a.min(b);
        let t_max = a.max(b);
        let denom = vb - va;
        let c = if denom.abs() > 1.0e-15 {
            ((a * vb - b * va) / denom).clamp(t_min + 0.02 * span, t_max - 0.02 * span)
        } else {
            0.5 * (a + b)
        };

        let state_c = eval_at_tt(c)?;
        let vc = state_c.rel_speed_deg_day();
        if vc.abs() < va.abs().min(vb.abs()) {
            best_state = state_c;
        }
        if vc.abs() < 1.0e-8 {
            return Ok(state_c);
        }

        if (vc > 0.0) == (vb > 0.0) {
            b = c;
            vb = vc;
            if side == -1 {
                va *= 0.5;
            }
            side = -1;
        } else {
            a = c;
            va = vc;
            if side == 1 {
                vb *= 0.5;
            }
            side = 1;
        }
    }

    Ok(best_state)
}

/// Checks a monotonic segment `[state_a, state_b]` for crossings of `(aspect_name, target_angle_deg, oriented_target_deg)`
/// and appends any found root to `out`.
#[allow(clippy::too_many_arguments)]
fn check_monotonic_segment_for_target<F>(
    eval_at_tt: &F,
    body1: BodyId,
    body2: BodyId,
    aspect_name: &str,
    target_angle_deg: f64,
    oriented_target_deg: f64,
    state_a: PairEvalState,
    state_b: PairEvalState,
    out: &mut Vec<ExactAspectEvent>,
) -> Result<(), EphemerisError>
where
    F: Fn(f64) -> Result<PairEvalState, EphemerisError>,
{
    let fa = state_a.error_for_oriented_target(oriented_target_deg);
    let fb = state_b.error_for_oriented_target(oriented_target_deg);

    if fa.abs() < 90.0 && fb.abs() < 90.0 && fa * fb <= 0.0 {
        let ev = refine_exact_aspect_root(
            eval_at_tt,
            body1,
            body2,
            aspect_name,
            target_angle_deg,
            oriented_target_deg,
            state_a,
            state_b,
        )?;
        let duplicate = out.iter().any(|prev| {
            prev.body1 == ev.body1
                && prev.body2 == ev.body2
                && (prev.target_angle_deg - ev.target_angle_deg).abs() < 1.0e-6
                && (prev.jd_tt - ev.jd_tt).abs() < 1.0e-5
        });
        if !duplicate {
            out.push(ev);
        }
    }
    Ok(())
}

/// Expands an aspect definition or angle into its 1 or 2 oriented targets $\Delta\theta \in [0^\circ, 360^\circ)$.
/// - For $A = 0^\circ$ or $A = 180^\circ$: 1 oriented target ($0^\circ$ or $180^\circ$).
/// - For $0^\circ < A < 180^\circ$: 2 oriented targets ($A$ and $360^\circ - A$).
/// - For $180^\circ < A < 360^\circ$: 1 explicit oriented target ($A$).
fn oriented_targets_for_angle(angle_deg: f64) -> (f64, Vec<f64>) {
    let norm = angle_deg.rem_euclid(360.0);
    if norm <= 1.0e-9 || (norm - 360.0).abs() <= 1.0e-9 {
        (0.0, vec![0.0])
    } else if (norm - 180.0).abs() <= 1.0e-9 {
        (180.0, vec![180.0])
    } else if norm < 180.0 {
        (norm, vec![norm, 360.0 - norm])
    } else {
        (norm, vec![norm])
    }
}

/// Inspects a single coarse step `[state_a, state_b]` for a list of `(aspect_name, target_angle_deg, oriented_targets)`.
/// If the relative velocity changes sign ($v_a \cdot v_b \le 0$), solves for the relative stationary epoch
/// $t_{\text{stat}} \in (t_a, t_b)$ and checks both $[t_a, t_{\text{stat}}]$ and $[t_{\text{stat}}, t_b]$
/// so grazing retrograde double/triple passes are never missed.
fn collect_pair_aspects_in_step<F>(
    eval_at_tt: &F,
    body1: BodyId,
    body2: BodyId,
    targets: &[(String, f64, Vec<f64>)],
    state_a: PairEvalState,
    state_b: PairEvalState,
    out: &mut Vec<ExactAspectEvent>,
) -> Result<(), EphemerisError>
where
    F: Fn(f64) -> Result<PairEvalState, EphemerisError>,
{
    let va = state_a.rel_speed_deg_day();
    let vb = state_b.rel_speed_deg_day();

    // Check whether any target is within 30° of this step before doing station root-finding
    let near_any_target = targets.iter().any(|(_, _, oriented_list)| {
        oriented_list.iter().any(|&ot| {
            state_a.error_for_oriented_target(ot).abs() < 30.0
                || state_b.error_for_oriented_target(ot).abs() < 30.0
        })
    });

    if va * vb <= 0.0 && near_any_target && (state_b.jd_tt - state_a.jd_tt).abs() > 1.0e-6 {
        let state_stat = solve_relative_station(eval_at_tt, state_a, state_b)?;
        for (aspect_name, target_angle_deg, oriented_list) in targets {
            for &oriented_target_deg in oriented_list {
                check_monotonic_segment_for_target(
                    eval_at_tt,
                    body1,
                    body2,
                    aspect_name,
                    *target_angle_deg,
                    oriented_target_deg,
                    state_a,
                    state_stat,
                    out,
                )?;
                check_monotonic_segment_for_target(
                    eval_at_tt,
                    body1,
                    body2,
                    aspect_name,
                    *target_angle_deg,
                    oriented_target_deg,
                    state_stat,
                    state_b,
                    out,
                )?;
            }
        }
        return Ok(());
    }

    for (aspect_name, target_angle_deg, oriented_list) in targets {
        for &oriented_target_deg in oriented_list {
            check_monotonic_segment_for_target(
                eval_at_tt,
                body1,
                body2,
                aspect_name,
                *target_angle_deg,
                oriented_target_deg,
                state_a,
                state_b,
                out,
            )?;
        }
    }

    Ok(())
}

/// Characteristic coarse step in days for a pair of bodies.
fn characteristic_pair_step_days(body1: BodyId, body2: BodyId) -> f64 {
    let is_fast_lunar = |b: BodyId| {
        matches!(
            b,
            BodyId::Planet(PlanetId::Moon)
                | BodyId::LunarPoint(
                    LunarPointId::TrueLilith
                        | LunarPointId::InterpolatedLilith
                        | LunarPointId::TrueNode
                        | LunarPointId::TrueSouthNode
                )
        )
    };
    if is_fast_lunar(body1) || is_fast_lunar(body2) {
        0.5
    } else if matches!(
        body1,
        BodyId::Planet(PlanetId::Sun | PlanetId::Mercury | PlanetId::Venus)
    ) || matches!(
        body2,
        BodyId::Planet(PlanetId::Sun | PlanetId::Mercury | PlanetId::Venus)
    ) {
        1.0
    } else {
        2.0
    }
}

/// Characteristic maximum search window in days for [`search_exact_aspect`].
fn characteristic_pair_max_search_days(body1: BodyId, body2: BodyId) -> f64 {
    if matches!(body1, BodyId::Planet(PlanetId::Moon))
        || matches!(body2, BodyId::Planet(PlanetId::Moon))
    {
        40.0
    } else if matches!(
        body1,
        BodyId::Planet(PlanetId::Sun | PlanetId::Mercury | PlanetId::Venus)
    ) || matches!(
        body2,
        BodyId::Planet(PlanetId::Sun | PlanetId::Mercury | PlanetId::Venus)
    ) {
        600.0
    } else if matches!(body1, BodyId::Planet(PlanetId::Mars | PlanetId::Jupiter))
        || matches!(body2, BodyId::Planet(PlanetId::Mars | PlanetId::Jupiter))
    {
        5_000.0
    } else {
        20_000.0
    }
}

/// Searches forward (`forward == true`) or backward (`forward == false`) from `start_input`
/// for the first exact aspect of `target_angle_deg` between `body1` and `body2` under `options`.
///
/// For $0^\circ < \text{target\_angle\_deg} < 180^\circ$, both $+\Delta\theta = A$ and
/// $-\Delta\theta = 360^\circ - A$ are searched and the earliest occurrence is returned.
pub fn search_exact_aspect(
    engine: &Engine,
    body1: BodyId,
    body2: BodyId,
    target_angle_deg: f64,
    start_input: TimeInput,
    forward: bool,
    options: CalcOptions,
) -> Result<ExactAspectEvent, EphemerisError> {
    if body1 == BodyId::EclNut || body2 == BodyId::EclNut {
        return Err(EphemerisError::UnsupportedBody {
            body: "EclNut".to_owned(),
        });
    }
    if body1 == body2 {
        return Err(EphemerisError::InvalidInput {
            message: "Cannot search_exact_aspect between a body and itself".to_owned(),
        });
    }
    if !target_angle_deg.is_finite() {
        return Err(EphemerisError::InvalidInput {
            message: format!("Invalid target_angle_deg: {target_angle_deg}"),
        });
    }

    let mut eval_opts = options;
    eval_opts.plane = Plane::Ecliptic;
    eval_opts.include_speed = true;

    let obs = match eval_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let (_, tt_start, _, _) = engine.time_table().resolve_scales(start_input, obs);
    let jd_tt_start = tt_start.jd_f64();

    let (canonical_angle, oriented_list) = oriented_targets_for_angle(target_angle_deg);
    let aspect_name = canonical_aspect_name(canonical_angle);
    let targets = vec![(aspect_name, canonical_angle, oriented_list)];

    let step_days = characteristic_pair_step_days(body1, body2);
    let max_search_days = characteristic_pair_max_search_days(body1, body2);
    let total_steps = ((max_search_days / step_days).ceil() as usize).max(2);
    let dir = if forward { 1.0 } else { -1.0 };

    let eval_closure =
        |t: f64| -> Result<PairEvalState, EphemerisError> { eval_moving_pair_at_tt(engine, body1, body2, t, eval_opts) };

    let mut state_prev = eval_closure(jd_tt_start)?;
    let mut found = Vec::new();

    for step in 1..=total_steps {
        let dt = ((step as f64) * step_days).min(max_search_days);
        let t_curr = jd_tt_start + dir * dt;
        let state_curr = eval_closure(t_curr)?;

        collect_pair_aspects_in_step(
            &eval_closure,
            body1,
            body2,
            &targets,
            state_prev,
            state_curr,
            &mut found,
        )?;

        if !found.is_empty() {
            if forward {
                found.sort_by(|a, b| a.jd_tt.total_cmp(&b.jd_tt));
            } else {
                found.sort_by(|a, b| b.jd_tt.total_cmp(&a.jd_tt));
            }
            if let Some(first) = found.first().cloned() {
                return Ok(first);
            }
        }

        state_prev = state_curr;
    }

    Err(EphemerisError::NonConvergence {
        solver: "search_exact_aspect".to_owned(),
        iterations: total_steps as u32,
        residual: wrap_pm180(state_prev.rel_lon_deg() - canonical_angle).abs(),
    })
}

/// Finds all exact mutual aspects for `preset` among unique pairs `(i < j)` of `bodies`
/// within `[start_input, end_input]`, sorted in chronological order.
///
/// Shares a single precomputed coarse [`CalcResult`] grid across all bodies and splits steps
/// at relative velocity stations ($v_{\text{rel}}(t) = 0$) so retrograde triple passes and
/// grazing double passes are captured with `< 1e-9°` accuracy.
pub fn mutual_aspects_in_range(
    engine: &Engine,
    bodies: &[BodyId],
    preset: AspectPreset,
    start_input: TimeInput,
    end_input: TimeInput,
    options: CalcOptions,
) -> Result<Vec<ExactAspectEvent>, EphemerisError> {
    for &b in bodies {
        if b == BodyId::EclNut {
            return Err(EphemerisError::UnsupportedBody {
                body: "EclNut".to_owned(),
            });
        }
    }

    let mut eval_opts = options;
    eval_opts.plane = Plane::Ecliptic;
    eval_opts.include_speed = true;

    let obs = match eval_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let (_, tt_start, _, _) = engine.time_table().resolve_scales(start_input, obs);
    let (_, tt_end, _, _) = engine.time_table().resolve_scales(end_input, obs);
    let jd_tt_start = tt_start.jd_f64();
    let jd_tt_end = tt_end.jd_f64();

    if !jd_tt_start.is_finite() || !jd_tt_end.is_finite() || jd_tt_end < jd_tt_start {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid time range in mutual_aspects_in_range: start={jd_tt_start}, end={jd_tt_end}"
            ),
        });
    }
    if bodies.len() < 2 || (jd_tt_end - jd_tt_start) < 1.0e-9 {
        return Ok(Vec::new());
    }

    let defs = aspect_definitions(preset);
    let targets: Vec<(String, f64, Vec<f64>)> = defs
        .into_iter()
        .map(|d| {
            let (_, oriented) = oriented_targets_for_angle(d.angle_deg);
            (d.name, d.angle_deg, oriented)
        })
        .collect();

    let has_fast_body = bodies.iter().any(|&b| {
        matches!(
            b,
            BodyId::Planet(PlanetId::Moon) | BodyId::LunarPoint(_)
        )
    });
    let span = jd_tt_end - jd_tt_start;
    let base_step = if has_fast_body { 0.5_f64 } else { 1.0_f64 };
    let step_days = base_step.min(span * 0.5).max(1.0e-3);
    let total_steps = ((span / step_days).ceil() as usize).max(1);

    // Precompute coarse grid once across all requested bodies
    let mut grid: Vec<Vec<CalcResult>> = Vec::with_capacity(total_steps + 1);
    for step in 0..=total_steps {
        let t = if step == total_steps {
            jd_tt_end
        } else {
            (jd_tt_start + (step as f64) * step_days).min(jd_tt_end)
        };
        let moment = engine.at(TimeInput::from_tt_jd(t), obs)?;
        let mut row = Vec::with_capacity(bodies.len());
        for &body in bodies {
            row.push(engine.calc_in_moment(&moment, body, eval_opts)?);
        }
        grid.push(row);
    }

    let mut events = Vec::new();

    for (i, &body1) in bodies.iter().enumerate() {
        for (j, &body2) in bodies.iter().enumerate().skip(i + 1) {
            if body1 == body2 {
                continue;
            }

            let eval_closure = |t: f64| -> Result<PairEvalState, EphemerisError> {
                eval_moving_pair_at_tt(engine, body1, body2, t, eval_opts)
            };

            for pair_rows in grid.windows(2) {
                let [row_a, row_b] = pair_rows else {
                    continue;
                };
                let (Some(r1_a), Some(r2_a), Some(r1_b), Some(r2_b)) =
                    (row_a.get(i), row_a.get(j), row_b.get(i), row_b.get(j))
                else {
                    continue;
                };

                let state_a = PairEvalState {
                    jd_tt: r1_a.jd_tt,
                    jd_ut1: r1_a.jd_ut1,
                    lon1_deg: r1_a.lon_deg,
                    lon2_deg: r2_a.lon_deg,
                    speed1_deg_day: r1_a.lon_speed_deg_day,
                    speed2_deg_day: r2_a.lon_speed_deg_day,
                };
                let state_b = PairEvalState {
                    jd_tt: r1_b.jd_tt,
                    jd_ut1: r1_b.jd_ut1,
                    lon1_deg: r1_b.lon_deg,
                    lon2_deg: r2_b.lon_deg,
                    speed1_deg_day: r1_b.lon_speed_deg_day,
                    speed2_deg_day: r2_b.lon_speed_deg_day,
                };

                collect_pair_aspects_in_step(
                    &eval_closure,
                    body1,
                    body2,
                    &targets,
                    state_a,
                    state_b,
                    &mut events,
                )?;
            }
        }
    }

    events.sort_by(|a, b| a.jd_tt.total_cmp(&b.jd_tt));
    Ok(events)
}

/// Finds all exact transit-to-natal aspects where `transiting_bodies` aspect the fixed positions
/// of `natal_bodies` (evaluated at `natal_input`) over `[start_input, end_input]` under `options`.
#[allow(clippy::too_many_arguments)]
pub fn natal_transits_in_range(
    engine: &Engine,
    transiting_bodies: &[BodyId],
    natal_input: TimeInput,
    natal_bodies: &[BodyId],
    preset: AspectPreset,
    start_input: TimeInput,
    end_input: TimeInput,
    options: CalcOptions,
) -> Result<Vec<ExactAspectEvent>, EphemerisError> {
    for &b in transiting_bodies.iter().chain(natal_bodies.iter()) {
        if b == BodyId::EclNut {
            return Err(EphemerisError::UnsupportedBody {
                body: "EclNut".to_owned(),
            });
        }
    }

    let mut eval_opts = options;
    eval_opts.plane = Plane::Ecliptic;
    eval_opts.include_speed = true;

    let obs = match eval_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let (_, tt_start, _, _) = engine.time_table().resolve_scales(start_input, obs);
    let (_, tt_end, _, _) = engine.time_table().resolve_scales(end_input, obs);
    let jd_tt_start = tt_start.jd_f64();
    let jd_tt_end = tt_end.jd_f64();

    if !jd_tt_start.is_finite() || !jd_tt_end.is_finite() || jd_tt_end < jd_tt_start {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid time range in natal_transits_in_range: start={jd_tt_start}, end={jd_tt_end}"
            ),
        });
    }
    if transiting_bodies.is_empty() || natal_bodies.is_empty() || (jd_tt_end - jd_tt_start) < 1.0e-9
    {
        return Ok(Vec::new());
    }

    // 1. Evaluate all natal bodies once at natal_input
    let natal_moment = engine.at(natal_input, obs)?;
    let mut natal_states: Vec<(BodyId, f64, bool)> = Vec::with_capacity(natal_bodies.len());
    for &n_body in natal_bodies {
        let res = engine.calc_in_moment(&natal_moment, n_body, eval_opts)?;
        natal_states.push((n_body, res.lon_deg, res.lon_speed_deg_day < 0.0));
    }

    let defs = aspect_definitions(preset);
    let targets: Vec<(String, f64, Vec<f64>)> = defs
        .into_iter()
        .map(|d| {
            let (_, oriented) = oriented_targets_for_angle(d.angle_deg);
            (d.name, d.angle_deg, oriented)
        })
        .collect();

    let has_fast_body = transiting_bodies.iter().any(|&b| {
        matches!(
            b,
            BodyId::Planet(PlanetId::Moon) | BodyId::LunarPoint(_)
        )
    });
    let span = jd_tt_end - jd_tt_start;
    let base_step = if has_fast_body { 0.5_f64 } else { 1.0_f64 };
    let step_days = base_step.min(span * 0.5).max(1.0e-3);
    let total_steps = ((span / step_days).ceil() as usize).max(1);

    // 2. Precompute coarse grid for transiting bodies across [jd_tt_start, jd_tt_end]
    let mut grid: Vec<Vec<CalcResult>> = Vec::with_capacity(total_steps + 1);
    for step in 0..=total_steps {
        let t = if step == total_steps {
            jd_tt_end
        } else {
            (jd_tt_start + (step as f64) * step_days).min(jd_tt_end)
        };
        let moment = engine.at(TimeInput::from_tt_jd(t), obs)?;
        let mut row = Vec::with_capacity(transiting_bodies.len());
        for &t_body in transiting_bodies {
            row.push(engine.calc_in_moment(&moment, t_body, eval_opts)?);
        }
        grid.push(row);
    }

    let mut events = Vec::new();

    for (t_idx, &t_body) in transiting_bodies.iter().enumerate() {
        for &(n_body, natal_lon_deg, natal_retrograde) in &natal_states {
            let eval_closure = |t: f64| -> Result<PairEvalState, EphemerisError> {
                eval_transit_pair_at_tt(
                    engine,
                    t_body,
                    natal_lon_deg,
                    natal_retrograde,
                    t,
                    eval_opts,
                )
            };

            for pair_rows in grid.windows(2) {
                let [row_a, row_b] = pair_rows else {
                    continue;
                };
                let (Some(r_a), Some(r_b)) = (row_a.get(t_idx), row_b.get(t_idx)) else {
                    continue;
                };

                let state_a = PairEvalState {
                    jd_tt: r_a.jd_tt,
                    jd_ut1: r_a.jd_ut1,
                    lon1_deg: r_a.lon_deg,
                    lon2_deg: natal_lon_deg,
                    speed1_deg_day: r_a.lon_speed_deg_day,
                    speed2_deg_day: 0.0,
                };
                let state_b = PairEvalState {
                    jd_tt: r_b.jd_tt,
                    jd_ut1: r_b.jd_ut1,
                    lon1_deg: r_b.lon_deg,
                    lon2_deg: natal_lon_deg,
                    speed1_deg_day: r_b.lon_speed_deg_day,
                    speed2_deg_day: 0.0,
                };

                let mut step_events = Vec::new();
                collect_pair_aspects_in_step(
                    &eval_closure,
                    t_body,
                    n_body,
                    &targets,
                    state_a,
                    state_b,
                    &mut step_events,
                )?;
                for mut ev in step_events {
                    ev.body2_retrograde = natal_retrograde;
                    let duplicate = events.iter().any(|prev: &ExactAspectEvent| {
                        prev.body1 == ev.body1
                            && prev.body2 == ev.body2
                            && (prev.target_angle_deg - ev.target_angle_deg).abs() < 1.0e-6
                            && (prev.jd_tt - ev.jd_tt).abs() < 1.0e-5
                    });
                    if !duplicate {
                        events.push(ev);
                    }
                }
            }
        }
    }

    events.sort_by(|a, b| a.jd_tt.total_cmp(&b.jd_tt));
    Ok(events)
}

// ============================================================================
// 4. Void of Course Moon Solver
// ============================================================================

/// Void-of-Course (VoC) Moon interval between the Moon's last exact major aspect in a zodiac sign
/// and its ingress into the next zodiac sign.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VoidOfCourseInterval {
    /// `JD(UT1)` of the last exact major aspect made by the Moon before leaving `moon_sign_before`.
    pub start_jd_ut1: f64,
    /// `JD(UT1)` when the Moon ingresses into `moon_sign_after` ($\lambda_{\text{Moon}} = k \cdot 30^\circ$).
    pub end_jd_ut1: f64,
    /// Duration of the Void-of-Course interval in hours (`(end_jd_ut1 - start_jd_ut1) * 24.0`).
    pub duration_hours: f64,
    /// Planet with which the Moon formed its last exact major aspect before sign ingress.
    pub last_aspect_body: BodyId,
    /// Canonical name of that last major aspect (`"Conjunction"`, `"Sextile"`, `"Square"`, `"Trine"`, `"Opposition"`).
    pub last_aspect_name: String,
    /// Angle of that last major aspect in degrees (`0.0`, `60.0`, `90.0`, `120.0`, `180.0`).
    pub last_aspect_angle_deg: f64,
    /// Zodiac sign index `0..=11` (`0 = Aries`, ..., `11 = Pisces`) before ingress.
    pub moon_sign_before: u8,
    /// Zodiac sign index `0..=11` (`0 = Aries`, ..., `11 = Pisces`) entered at `end_jd_ut1`.
    pub moon_sign_after: u8,
    /// `true` when `moon_sign_before` is Taurus (`1`), Cancer (`3`), Sagittarius (`8`), or Pisces (`11`),
    /// where the Moon "perfmeth somewhat" despite being Void of Course per Guido Bonatti
    /// (*Anima Astrologiae*, Considerations 113 & 141) and William Lilly (*Christian Astrology*, Book I, pp. 113, 141).
    pub is_bonatti_mild_sign: bool,
}

/// Finds all Void-of-Course Moon intervals in `[start_input, end_input]` under `options`:
/// - Scans lunar zodiac sign ingresses ($\lambda_{\text{Moon}} = k \cdot 30^\circ, k = 0..11$)
///   starting `2.8` days before `start_input` so the Moon's sign entry preceding the first
///   ingress in `[start_input, end_input]` is always captured.
/// - For each sign ingress in `[start_input, end_input]`, finds all exact major aspects
///   (`0°, 60°, 90°, 120°, 180°`) between the Moon and `major_planets` (defaulting to
///   [`DEFAULT_VOC_MAJOR_PLANETS`]: `Sun, Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto`)
///   during the Moon's transit of that sign, and selects the final exact aspect before ingress.
pub fn void_of_course_moon_in_range(
    engine: &Engine,
    start_input: TimeInput,
    end_input: TimeInput,
    major_planets: Option<&[BodyId]>,
    options: CalcOptions,
) -> Result<Vec<VoidOfCourseInterval>, EphemerisError> {
    let mut pos_opts = options;
    pos_opts.plane = Plane::Ecliptic;
    pos_opts.include_speed = false;

    let obs = match pos_opts.origin {
        Origin::Topocentric(o) => Some(o),
        _ => None,
    };
    let (ut1_start, tt_start, _, _) = engine.time_table().resolve_scales(start_input, obs);
    let (ut1_end, tt_end, _, _) = engine.time_table().resolve_scales(end_input, obs);
    let jd_tt_start = tt_start.jd_f64();
    let jd_tt_end = tt_end.jd_f64();
    let jd_ut1_start = ut1_start.jd_f64();
    let jd_ut1_end = ut1_end.jd_f64();

    if !jd_tt_start.is_finite() || !jd_tt_end.is_finite() || jd_tt_end < jd_tt_start {
        return Err(EphemerisError::InvalidInput {
            message: format!(
                "Invalid time range in void_of_course_moon_in_range: start={jd_tt_start}, end={jd_tt_end}"
            ),
        });
    }
    if (jd_tt_end - jd_tt_start) < 1.0e-9 {
        return Ok(Vec::new());
    }

    let moon = BodyId::Planet(PlanetId::Moon);
    let planet_slice = match major_planets {
        Some(list) if !list.is_empty() => list,
        _ => &DEFAULT_VOC_MAJOR_PLANETS,
    };
    let planets: Vec<BodyId> = planet_slice
        .iter()
        .copied()
        .filter(|&b| b != moon)
        .collect();
    if planets.is_empty() {
        return Ok(Vec::new());
    }

    // 1. Precompute coarse position-only grid once for [Moon] + planets from (jd_tt_start - 2.8)
    //    to (jd_tt_end + 2.8) so any ongoing Void-of-Course interval that starts before jd_ut1_end
    //    and finishes at the next sign ingress after jd_ut1_end is also captured.
    let lookback_days = 2.8_f64;
    let lookahead_days = 2.8_f64;
    let moon_cov = engine.coverage(moon)?;
    let scan_start_tt = (jd_tt_start - lookback_days).max(moon_cov.safe_start_jd_tt + 1.0e-6);
    let scan_end_tt = (jd_tt_end + lookahead_days).min(moon_cov.safe_end_jd_tt - 1.0e-6);
    let span = (scan_end_tt - scan_start_tt).max(1.0e-6);
    let step_days = 0.5_f64;
    let total_steps = ((span / step_days).ceil() as usize).max(1);

    let mut all_bodies = Vec::with_capacity(1 + planets.len());
    all_bodies.push(moon);
    all_bodies.extend_from_slice(&planets);

    let mut grid: Vec<Vec<CalcResult>> = Vec::with_capacity(total_steps + 1);
    for step in 0..=total_steps {
        let t = if step == total_steps {
            scan_end_tt
        } else {
            (scan_start_tt + (step as f64) * step_days).min(scan_end_tt)
        };
        let moment = engine.at(TimeInput::from_tt_jd(t), obs)?;
        let mut row = Vec::with_capacity(all_bodies.len());
        for &b in &all_bodies {
            row.push(engine.calc_in_moment(&moment, b, pos_opts)?);
        }
        grid.push(row);
    }

    // Fast position-only secant/Illinois refiner for Moon sign ingress on [m_a, m_b]
    let refine_moon_ingress_fast =
        |target_lon: f64, m_a: &CalcResult, m_b: &CalcResult| -> Result<(f64, f64), EphemerisError> {
            let mut a = m_a.jd_tt;
            let mut b = m_b.jd_tt;
            let mut fa = wrap_pm180(m_a.lon_deg - target_lon);
            let mut fb = wrap_pm180(m_b.lon_deg - target_lon);
            let mut best = if fa.abs() <= fb.abs() {
                (m_a.jd_tt, m_a.jd_ut1)
            } else {
                (m_b.jd_tt, m_b.jd_ut1)
            };
            let mut side = 0i8;

            for _ in 0..25 {
                let s = (b - a).abs();
                if s < 1.0e-10 || fa.abs().min(fb.abs()) < 1.0e-9 {
                    return Ok(best);
                }
                let denom = fb - fa;
                let c = if denom.abs() > 1.0e-15 {
                    ((a * fb - b * fa) / denom).clamp(a + 0.01 * s, b - 0.01 * s)
                } else {
                    0.5 * (a + b)
                };
                let rc = engine.calc(moon, TimeInput::from_tt_jd(c), pos_opts)?;
                let fc = wrap_pm180(rc.lon_deg - target_lon);
                if fc.abs() < fa.abs().min(fb.abs()) {
                    best = (rc.jd_tt, rc.jd_ut1);
                }
                if fc.abs() < 1.0e-9 {
                    return Ok((rc.jd_tt, rc.jd_ut1));
                }
                if (fc > 0.0) == (fb > 0.0) {
                    b = c;
                    fb = fc;
                    if side == -1 {
                        fa *= 0.5;
                    }
                    side = -1;
                } else {
                    a = c;
                    fa = fc;
                    if side == 1 {
                        fb *= 0.5;
                    }
                    side = 1;
                }
            }
            Ok(best)
        };

    // Fast position-only secant/Illinois refiner for Moon-planet aspect on [step, step+1]
    let refine_moon_aspect_fast = |planet: BodyId,
                                   oriented_target_deg: f64,
                                   m_a: &CalcResult,
                                   p_a: &CalcResult,
                                   m_b: &CalcResult,
                                   p_b: &CalcResult|
     -> Result<(f64, f64), EphemerisError> {
        let mut a = m_a.jd_tt;
        let mut b = m_b.jd_tt;
        let mut fa = wrap_pm180(m_a.lon_deg - p_a.lon_deg - oriented_target_deg);
        let mut fb = wrap_pm180(m_b.lon_deg - p_b.lon_deg - oriented_target_deg);
        let mut best = if fa.abs() <= fb.abs() {
            (m_a.jd_tt, m_a.jd_ut1)
        } else {
            (m_b.jd_tt, m_b.jd_ut1)
        };
        let mut side = 0i8;

        for _ in 0..25 {
            let s = (b - a).abs();
            if s < 1.0e-10 || fa.abs().min(fb.abs()) < 1.0e-9 {
                return Ok(best);
            }
            let denom = fb - fa;
            let c = if denom.abs() > 1.0e-15 {
                ((a * fb - b * fa) / denom).clamp(a + 0.01 * s, b - 0.01 * s)
            } else {
                0.5 * (a + b)
            };
            let moment = engine.at(TimeInput::from_tt_jd(c), obs)?;
            let mc = engine.calc_in_moment(&moment, moon, pos_opts)?;
            let pc = engine.calc_in_moment(&moment, planet, pos_opts)?;
            let fc = wrap_pm180(mc.lon_deg - pc.lon_deg - oriented_target_deg);
            if fc.abs() < fa.abs().min(fb.abs()) {
                best = (mc.jd_tt, mc.jd_ut1);
            }
            if fc.abs() < 1.0e-9 {
                return Ok((mc.jd_tt, mc.jd_ut1));
            }
            if (fc > 0.0) == (fb > 0.0) {
                b = c;
                fb = fc;
                if side == -1 {
                    fa *= 0.5;
                }
                side = -1;
            } else {
                a = c;
                fa = fc;
                if side == 1 {
                    fb *= 0.5;
                }
                side = 1;
            }
        }
        Ok(best)
    };

    // 2. Detect and refine all Moon sign ingresses using the already-computed Moon column (col 0)
    // Each ingress stores (sign_after: u8, jd_tt: f64, jd_ut1: f64, step_idx: usize)
    let mut ingresses: Vec<(u8, f64, f64, usize)> = Vec::new();
    for (step, pair_rows) in grid.windows(2).enumerate() {
        let [row_a, row_b] = pair_rows else {
            continue;
        };
        let (Some(m_a), Some(m_b)) = (row_a.first(), row_b.first()) else {
            continue;
        };
        for sign_idx in 0u8..12u8 {
            let target_lon = f64::from(sign_idx) * 30.0;
            let f_prev = wrap_pm180(m_a.lon_deg - target_lon);
            let f_curr = wrap_pm180(m_b.lon_deg - target_lon);

            if f_prev <= 0.0 && f_curr >= 0.0 && f_prev.abs() < 25.0 && f_curr.abs() < 25.0 {
                let (ing_tt, ing_ut1) = refine_moon_ingress_fast(target_lon, m_a, m_b)?;
                if ing_tt >= scan_start_tt - 1.0e-7
                    && ing_tt <= scan_end_tt + 1.0e-7
                    && ingresses
                        .last()
                        .map_or(true, |last| (ing_tt - last.1).abs() > 0.5)
                {
                    ingresses.push((sign_idx, ing_tt, ing_ut1, step));
                }
            }
        }
    }
    ingresses.sort_by(|a, b| a.1.total_cmp(&b.1));

    // 3. For each ingress with ingress_ut1 >= jd_ut1_start, scan coarse steps backward from
    //    the ingress step to the preceding ingress step, refining only the latest bracket(s),
    //    and include the VoC interval if [last_ut1, ingress_ut1] overlaps [jd_ut1_start, jd_ut1_end].
    let major_defs = aspect_definitions(AspectPreset::Major);
    let targets: Vec<(String, f64, Vec<f64>)> = major_defs
        .into_iter()
        .map(|d| {
            let (_, oriented) = oriented_targets_for_angle(d.angle_deg);
            (d.name, d.angle_deg, oriented)
        })
        .collect();

    let mut voc_intervals = Vec::new();

    for (idx, &(sign_after, ingress_tt, ingress_ut1, end_step)) in ingresses.iter().enumerate() {
        if ingress_ut1 < jd_ut1_start - 1.0e-9 {
            continue;
        }
        let is_post_window_ingress = ingress_ut1 > jd_ut1_end + 1.0e-9;
        let (sign_enter_tt, sign_enter_ut1, start_step) = idx
            .checked_sub(1)
            .and_then(|prev_idx| ingresses.get(prev_idx))
            .map_or((scan_start_tt, scan_start_tt, 0), |prev| (prev.1, prev.2, prev.3));
        let sign_before = (sign_after + 11) % 12;

        // Track best (jd_tt, jd_ut1, body2, aspect_name, target_angle_deg)
        let mut best_last: Option<(f64, f64, BodyId, String, f64)> = None;

        for step in (start_step..=end_step).rev() {
            let (Some(row_a), Some(row_b)) = (grid.get(step), grid.get(step + 1)) else {
                continue;
            };
            let (Some(m_a), Some(m_b)) = (row_a.first(), row_b.first()) else {
                continue;
            };

            // If we already found a valid aspect in a later step whose jd_tt >= m_b.jd_tt,
            // no earlier step can possibly have a larger jd_tt!
            if let Some((best_tt, _, _, _, _)) = &best_last {
                if *best_tt >= m_b.jd_tt - 1.0e-9 {
                    break;
                }
            }

            for (p_idx, &planet) in planets.iter().enumerate() {
                let col_idx = p_idx + 1;
                let (Some(p_a), Some(p_b)) = (row_a.get(col_idx), row_b.get(col_idx)) else {
                    continue;
                };

                for (aspect_name, target_angle_deg, oriented_list) in &targets {
                    for &oriented_target_deg in oriented_list {
                        let fa = wrap_pm180(m_a.lon_deg - p_a.lon_deg - oriented_target_deg);
                        let fb = wrap_pm180(m_b.lon_deg - p_b.lon_deg - oriented_target_deg);

                        if fa.abs() < 45.0 && fb.abs() < 45.0 && fa * fb <= 0.0 {
                            let (asp_tt, asp_ut1) = refine_moon_aspect_fast(
                                planet,
                                oriented_target_deg,
                                m_a,
                                p_a,
                                m_b,
                                p_b,
                            )?;
                            if asp_tt > sign_enter_tt + 1.0e-7 && asp_tt < ingress_tt - 1.0e-7 {
                                let replace = match &best_last {
                                    None => true,
                                    Some((prev_tt, _, _, _, _)) => asp_tt > *prev_tt,
                                };
                                if replace {
                                    best_last = Some((
                                        asp_tt,
                                        asp_ut1,
                                        planet,
                                        aspect_name.clone(),
                                        *target_angle_deg,
                                    ));
                                }
                            }
                        }
                    }
                }
            }
        }

        if let Some((_, last_ut1, last_body, last_name, last_angle)) = best_last {
            if last_ut1 <= jd_ut1_end + 1.0e-9 {
                let duration_hours = ((ingress_ut1 - last_ut1) * 24.0).max(0.0);
                voc_intervals.push(VoidOfCourseInterval {
                    start_jd_ut1: last_ut1,
                    end_jd_ut1: ingress_ut1,
                    duration_hours,
                    last_aspect_body: last_body,
                    last_aspect_name: last_name,
                    last_aspect_angle_deg: last_angle,
                    moon_sign_before: sign_before,
                    moon_sign_after: sign_after,
                    is_bonatti_mild_sign: matches!(sign_before, 1 | 3 | 8 | 11),
                });
            }
        } else if idx > 0 && sign_enter_ut1 <= jd_ut1_end + 1.0e-9 {
            let duration_hours = ((ingress_ut1 - sign_enter_ut1) * 24.0).max(0.0);
            voc_intervals.push(VoidOfCourseInterval {
                start_jd_ut1: sign_enter_ut1,
                end_jd_ut1: ingress_ut1,
                duration_hours,
                last_aspect_body: BodyId::Planet(PlanetId::Moon),
                last_aspect_name: "None".to_owned(),
                last_aspect_angle_deg: 0.0,
                moon_sign_before: sign_before,
                moon_sign_after: sign_after,
                is_bonatti_mild_sign: matches!(sign_before, 1 | 3 | 8 | 11),
            });
        }
        if is_post_window_ingress {
            break;
        }
    }

    Ok(voc_intervals)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::EngineConfig;
    use crate::types::{AyanamsaModel, MpcNumber};

    #[test]
    fn test_aspect_presets_counts_and_ordering() {
        let major = aspect_definitions(AspectPreset::Major);
        assert_eq!(major.len(), 5);
        assert!(major.iter().all(|d| d.is_major));
        let angles: Vec<f64> = major.iter().map(|d| d.angle_deg).collect();
        assert_eq!(angles, vec![0.0, 60.0, 90.0, 120.0, 180.0]);

        let standard = aspect_definitions(AspectPreset::Standard);
        assert_eq!(standard.len(), 9);

        let kepler = aspect_definitions(AspectPreset::Kepler);
        assert_eq!(kepler.len(), 11);

        let all = aspect_definitions(AspectPreset::All);
        assert_eq!(all, kepler);
    }

    #[test]
    fn test_static_aspect_applying_vs_separating_analytical() {
        let defs = aspect_definitions(AspectPreset::Major);

        // Case 1: Moon at 88° (speed +13°/d) and Sun at 0° (speed +1°/d) -> 88° separation,
        // Δλ = +88°, target = +90°, e = -2°, v_rel = +12°/d -> e * v_rel = -24 < 0 -> APPLYING Square!
        let samples_applying = vec![
            BodyPositionSample {
                body: BodyId::Planet(PlanetId::Moon),
                lon_deg: 88.0,
                lat_deg: 0.0,
                lon_speed_deg_day: 13.0,
                dec_deg: None,
                dec_speed_deg_day: None,
            },
            BodyPositionSample {
                body: BodyId::Planet(PlanetId::Sun),
                lon_deg: 0.0,
                lat_deg: 0.0,
                lon_speed_deg_day: 1.0,
                dec_deg: None,
                dec_speed_deg_day: None,
            },
        ];
        let hits = calc_aspects(&samples_applying, &defs, 1.0);
        assert_eq!(hits.len(), 1);
        let h = hits.first().expect("1 hit");
        assert_eq!(h.aspect_name, "Square");
        assert!((h.exact_angle_deg - 90.0).abs() < 1e-12);
        assert!((h.actual_separation_deg - 88.0).abs() < 1e-12);
        assert!((h.orb_deg - 2.0).abs() < 1e-12);
        assert!((h.tightness_pct - 75.0).abs() < 1e-12);
        assert!(h.is_applying, "88° moving toward 90° must be applying");

        // Case 2: Moon at 92° (speed +13°/d) and Sun at 0° (speed +1°/d) -> 92° separation,
        // e = +2°, v_rel = +12°/d -> e * v_rel = +24 > 0 -> SEPARATING Square!
        let samples_separating = vec![
            BodyPositionSample {
                body: BodyId::Planet(PlanetId::Moon),
                lon_deg: 92.0,
                lat_deg: 0.0,
                lon_speed_deg_day: 13.0,
                dec_deg: None,
                dec_speed_deg_day: None,
            },
            BodyPositionSample {
                body: BodyId::Planet(PlanetId::Sun),
                lon_deg: 0.0,
                lat_deg: 0.0,
                lon_speed_deg_day: 1.0,
                dec_deg: None,
                dec_speed_deg_day: None,
            },
        ];
        let hits_sep = calc_aspects(&samples_separating, &defs, 1.0);
        assert_eq!(hits_sep.len(), 1);
        assert!(
            !hits_sep.first().expect("1 sep hit").is_applying,
            "92° moving away from 90° must be separating"
        );

        // Case 3: Waning Square (Δλ = -88°, i.e., Moon at 272°, Sun at 0°):
        // Moon approaching 360° (Δλ increasing from -88° toward 0°, so separating from -90°!)
        // Wait: at Moon = 268° (Δλ = -92°), Moon is approaching 270° (Δλ = -90°):
        // e = -92° - (-1)*90° = -2°, v_rel = +12°/d -> e * v_rel = -24 < 0 -> APPLYING waning Square!
        let samples_waning_app = vec![
            BodyPositionSample {
                body: BodyId::Planet(PlanetId::Moon),
                lon_deg: 268.0,
                lat_deg: 0.0,
                lon_speed_deg_day: 13.0,
                dec_deg: None,
                dec_speed_deg_day: None,
            },
            BodyPositionSample {
                body: BodyId::Planet(PlanetId::Sun),
                lon_deg: 0.0,
                lat_deg: 0.0,
                lon_speed_deg_day: 1.0,
                dec_deg: None,
                dec_speed_deg_day: None,
            },
        ];
        let hits_waning = calc_aspects(&samples_waning_app, &defs, 1.0);
        assert_eq!(hits_waning.len(), 1);
        let hw = hits_waning.first().expect("1 waning hit");
        assert_eq!(hw.aspect_name, "Square");
        assert!(hw.is_applying, "Moon at 268° approaching 270° must be applying");
    }

    #[test]
    fn test_declination_parallel_and_contraparallel() {
        let engine = Engine::new_moshier(EngineConfig::default());
        let jd = TimeInput::from_ut1_jd(2_451_545.0);
        let bodies = [
            BodyId::Planet(PlanetId::Sun),
            BodyId::Planet(PlanetId::Moon),
            BodyId::Planet(PlanetId::Mercury),
            BodyId::Planet(PlanetId::Venus),
            BodyId::Planet(PlanetId::Mars),
            BodyId::Planet(PlanetId::Jupiter),
            BodyId::Planet(PlanetId::Saturn),
        ];
        let hits = calc_declination_aspects(&engine, jd, &bodies, 2.0, CalcOptions::default())
            .expect("calc_declination_aspects");
        for h in &hits {
            assert!(h.orb_deg <= 2.0);
            assert!((0.0..=100.0).contains(&h.tightness_pct));
            match h.kind {
                DeclinationAspectKind::Parallel => {
                    assert!(((h.dec1_deg - h.dec2_deg).abs() - h.orb_deg).abs() < 1e-12);
                }
                DeclinationAspectKind::Contraparallel => {
                    assert!(((h.dec1_deg + h.dec2_deg).abs() - h.orb_deg).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn test_search_exact_aspect_and_mutual_aspects_in_range() {
        let engine = Engine::new_moshier(EngineConfig::default());
        let start = TimeInput::from_ut1_jd(2_461_041.5); // 2026-01-01 00:00 UT
        let end = TimeInput::from_ut1_jd(2_461_071.5);   // 2026-01-31 00:00 UT
        let opts = CalcOptions::default();

        // Next Moon-Sun opposition (Full Moon = 180°) after 2026-01-01
        let opp = search_exact_aspect(
            &engine,
            BodyId::Planet(PlanetId::Moon),
            BodyId::Planet(PlanetId::Sun),
            180.0,
            start,
            true,
            opts,
        )
        .expect("search_exact_aspect Moon-Sun 180");
        assert_eq!(opp.aspect_name, "Opposition");
        let sep = wrap_pm180(opp.body1_lon_deg - opp.body2_lon_deg).abs();
        assert!((sep - 180.0).abs() < 1.0e-7, "sep={sep}");

        // Next Moon-Sun Square (90° or 270°) after 2026-01-01
        let sq = search_exact_aspect(
            &engine,
            BodyId::Planet(PlanetId::Moon),
            BodyId::Planet(PlanetId::Sun),
            90.0,
            start,
            true,
            opts,
        )
        .expect("search_exact_aspect Moon-Sun 90");
        assert_eq!(sq.aspect_name, "Square");
        let sq_sep = wrap_pm180(sq.body1_lon_deg - sq.body2_lon_deg).abs();
        assert!((sq_sep - 90.0).abs() < 1.0e-7, "sq_sep={sq_sep}");

        // Mutual aspects in range for Sun, Moon, Mercury, Venus, Mars over 30 days
        let bodies = [
            BodyId::Planet(PlanetId::Sun),
            BodyId::Planet(PlanetId::Moon),
            BodyId::Planet(PlanetId::Mercury),
            BodyId::Planet(PlanetId::Venus),
            BodyId::Planet(PlanetId::Mars),
        ];
        let events = mutual_aspects_in_range(
            &engine,
            &bodies,
            AspectPreset::Major,
            start,
            end,
            opts,
        )
        .expect("mutual_aspects_in_range");
        assert!(events.len() > 25, "Expected > 25 major aspects in 30 days, got {}", events.len());

        for pair in events.windows(2) {
            if let [a, b] = pair {
                assert!(a.jd_tt <= b.jd_tt, "Events must be sorted chronologically");
            }
        }
        for ev in &events {
            let actual_sep = wrap_pm180(ev.body1_lon_deg - ev.body2_lon_deg).abs();
            let err = wrap_pm180(actual_sep - ev.target_angle_deg).abs();
            assert!(
                err < 1.0e-7,
                "Aspect {:?}-{:?} {} target={} actual_sep={} err={err:e}",
                ev.body1,
                ev.body2,
                ev.aspect_name,
                ev.target_angle_deg,
                actual_sep
            );
        }
    }

    #[test]
    fn test_natal_transits_retrograde_triple_pass_station_splitting() {
        let engine = Engine::new_moshier(EngineConfig::default());
        let opts = CalcOptions::default();

        // Saturn retrograde loop in 2026: find Saturn's station retrograde in 2026 and place a
        // synthetic natal Sun slightly below Saturn's stationary longitude so Saturn crosses it
        // 3 times (direct -> retrograde -> direct).
        let start_2026 = TimeInput::from_ut1_jd(2_461_041.5); // 2026-01-01
        let end_2027 = TimeInput::from_ut1_jd(2_461_467.5);   // 2027-03-03
        let stations = engine
            .search_stations_in_range(
                BodyId::Planet(PlanetId::Saturn),
                start_2026,
                end_2027,
                &opts,
            )
            .expect("Saturn stations 2026");
        assert!(stations.len() >= 2);
        let sr = stations.first().expect("sr");
        let sd = stations.get(1).expect("sd");
        let mid_lon = 0.5 * (sr.lon_deg + sd.lon_deg);

        // Find when the Sun is at `mid_lon` in spring 2026 to use as `natal_input`
        let sun_at_mid = engine
            .search_body_longitude(
                BodyId::Planet(PlanetId::Sun),
                start_2026,
                mid_lon,
                true,
                200.0,
                &opts,
            )
            .expect("Sun at mid_lon");

        let transits = natal_transits_in_range(
            &engine,
            &[BodyId::Planet(PlanetId::Saturn)],
            TimeInput::from_tt_jd(sun_at_mid.jd_tt),
            &[BodyId::Planet(PlanetId::Sun)],
            AspectPreset::Major,
            start_2026,
            end_2027,
            opts,
        )
        .expect("Saturn transits to natal Sun");

        let conjunctions: Vec<&ExactAspectEvent> = transits
            .iter()
            .filter(|e| e.aspect_name == "Conjunction")
            .collect();
        assert_eq!(
            conjunctions.len(),
            3,
            "Saturn must make a 3-pass triple conjunction over the natal point inside its retrograde loop: {conjunctions:?}"
        );
        let c0 = conjunctions.first().expect("c0");
        let c1 = conjunctions.get(1).expect("c1");
        let c2 = conjunctions.get(2).expect("c2");
        assert!(!c0.body1_retrograde, "Pass 1 must be direct");
        assert!(c1.body1_retrograde, "Pass 2 must be retrograde");
        assert!(!c2.body1_retrograde, "Pass 3 must be direct");

        // Also test a grazing double pass very close to the stationary point (within 0.02° of station)
        let grazing_lon = sr.lon_deg - 0.02;
        let sun_at_grazing = engine
            .search_body_longitude(
                BodyId::Planet(PlanetId::Sun),
                start_2026,
                grazing_lon,
                true,
                200.0,
                &opts,
            )
            .expect("Sun at grazing_lon");
        let grazing_transits = natal_transits_in_range(
            &engine,
            &[BodyId::Planet(PlanetId::Saturn)],
            TimeInput::from_tt_jd(sun_at_grazing.jd_tt),
            &[BodyId::Planet(PlanetId::Sun)],
            AspectPreset::Major,
            TimeInput::from_tt_jd(sr.jd_tt - 15.0),
            TimeInput::from_tt_jd(sr.jd_tt + 15.0),
            opts,
        )
        .expect("Saturn grazing transits");
        let grazing_conj: Vec<&ExactAspectEvent> = grazing_transits
            .iter()
            .filter(|e| e.aspect_name == "Conjunction")
            .collect();
        assert_eq!(
            grazing_conj.len(),
            2,
            "Relative station splitting must capture both passes of a grazing 0.02° station turn: {grazing_conj:?}"
        );
    }

    #[test]
    fn test_void_of_course_moon_in_range_and_sidereal_support() {
        let engine = Engine::new_moshier(EngineConfig::default());
        let start = TimeInput::from_ut1_jd(2_461_041.5); // 2026-01-01
        let end = TimeInput::from_ut1_jd(2_461_071.5);   // 2026-01-31

        let voc_trop = void_of_course_moon_in_range(&engine, start, end, None, CalcOptions::default())
            .expect("VoC Moon tropical");
        // Moon changes signs every ~2.2-2.5 days -> ~12-14 ingresses in 30 days
        assert!(
            (11..=16).contains(&voc_trop.len()),
            "Expected 11..=16 VoC Moon intervals in 30 days, got {}",
            voc_trop.len()
        );
        for v in &voc_trop {
            assert!(v.end_jd_ut1 > v.start_jd_ut1);
            assert!(v.duration_hours > 0.0 && v.duration_hours < 72.0);
            assert_eq!((v.moon_sign_before + 1) % 12, v.moon_sign_after);
            assert!(
                matches!(
                    v.last_aspect_angle_deg as u32,
                    0 | 60 | 90 | 120 | 180
                ),
                "Unexpected last aspect angle: {}",
                v.last_aspect_angle_deg
            );
            // Verify serde_json round-trip (DeserializeOwned)
            let json = serde_json::to_string(v).expect("serialize VoC");
            let back: VoidOfCourseInterval = serde_json::from_str(&json).expect("deserialize VoC");
            assert_eq!(back.last_aspect_body, v.last_aspect_body);
            assert_eq!(back.last_aspect_name, v.last_aspect_name);
            assert_eq!(back.moon_sign_before, v.moon_sign_before);
            assert_eq!(back.moon_sign_after, v.moon_sign_after);
            assert!((back.start_jd_ut1 - v.start_jd_ut1).abs() < 1e-12);
            assert!((back.end_jd_ut1 - v.end_jd_ut1).abs() < 1e-12);
            assert!((back.duration_hours - v.duration_hours).abs() < 1e-12);
        }

        // Also verify Sidereal (Lahiri) VoC Moon works seamlessly
        let sid_opts = CalcOptions {
            zodiac: Zodiac::Sidereal(AyanamsaModel::Lahiri),
            ..Default::default()
        };
        let voc_sid = void_of_course_moon_in_range(&engine, start, end, None, sid_opts)
            .expect("VoC Moon sidereal Lahiri");
        assert!((11..=16).contains(&voc_sid.len()));

        // Verify a short query window strictly inside an ongoing VoC interval returns that VoC interval
        let v0 = voc_trop.first().expect("first VoC interval");
        let sub_start = v0.start_jd_ut1 + 0.25 * (v0.end_jd_ut1 - v0.start_jd_ut1);
        let sub_end = v0.start_jd_ut1 + 0.75 * (v0.end_jd_ut1 - v0.start_jd_ut1);
        let sub_voc = void_of_course_moon_in_range(
            &engine,
            TimeInput::from_ut1_jd(sub_start),
            TimeInput::from_ut1_jd(sub_end),
            None,
            CalcOptions::default(),
        )
        .expect("sub-window inside ongoing VoC");
        assert_eq!(sub_voc.len(), 1, "Ongoing VoC interval must be returned for interior sub-window");
        let sub_v0 = sub_voc.first().expect("sub_v0");
        assert!((sub_v0.start_jd_ut1 - v0.start_jd_ut1).abs() < 1.0e-7);
        assert!((sub_v0.end_jd_ut1 - v0.end_jd_ut1).abs() < 1.0e-7);

        // Verify chart & synastry aspects
        let chart_hits = calc_chart_aspects(
            &engine,
            start,
            &DEFAULT_VOC_MAJOR_PLANETS,
            AspectPreset::All,
            1.0,
            CalcOptions::default(),
        )
        .expect("calc_chart_aspects");
        assert!(!chart_hits.is_empty());

        let syn_hits = calc_synastry_aspects(
            &engine,
            end,
            &[BodyId::Planet(PlanetId::Sun), BodyId::Planet(PlanetId::Moon)],
            start,
            &DEFAULT_VOC_MAJOR_PLANETS,
            AspectPreset::Major,
            1.0,
            CalcOptions::default(),
        )
        .expect("calc_synastry_aspects");
        assert!(!syn_hits.is_empty());
        let _ = MpcNumber(2060);
    }
}
