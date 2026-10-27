//! `SpkSource` trait and `AniseSpkBackend` (`ARCHITECTURE.md` §1.2, §2.2.1).
//!
//! Isolates all DAF/SPK kernel access behind the internal [`SpkSource`] trait and statically
//! verifies `Send + Sync` on [`anise::almanac::Almanac`].

use crate::daf::{ eval_spk_type21, extract_segment_words, parse_daf_spk, DafSegmentDescriptor };
use crate::error::EphemerisError;
use crate::types::{AU_KM, JD_J2000, SECONDS_PER_DAY};
use anise::almanac::Almanac;
use anise::astro::Aberration;
use anise::constants::orientations::J2000;
use anise::naif::SPK;
use anise::prelude::{Epoch, Frame};
use serde::{Deserialize, Serialize};
use static_assertions::assert_impl_all;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

// Compile-time proof that `anise::almanac::Almanac` is `Send + Sync` (`PLAN.md` §0.1).
assert_impl_all!(Almanac: Send, Sync);

/// Raw geometric state vector in ICRF (`J2000` equatorial frame) returned by an [`SpkSource`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SpkState {
    /// Position vector `(x, y, z)` in kilometers in ICRF (`J2000`).
    pub pos_km: [f64; 3],
    /// Velocity vector `(vx, vy, vz)` in km/s in ICRF (`J2000`).
    pub vel_km_s: [f64; 3],
    /// Requested NAIF target ID.
    pub target_id: i32,
    /// Requested NAIF observer/center ID.
    pub observer_id: i32,
}

impl SpkState {
    /// Converts position from kilometers to Astronomical Units (AU).
    #[must_use]
    pub fn pos_au(&self) -> [f64; 3] {
        [
            self.pos_km[0] / AU_KM,
            self.pos_km[1] / AU_KM,
            self.pos_km[2] / AU_KM,
        ]
    }

    /// Converts velocity from km/s to AU/day.
    #[must_use]
    pub fn vel_au_day(&self) -> [f64; 3] {
        let scale = SECONDS_PER_DAY / AU_KM;
        [
            self.vel_km_s[0] * scale,
            self.vel_km_s[1] * scale,
            self.vel_km_s[2] * scale,
        ]
    }
}

/// Merged continuous coverage interval `[start_jd_tdb, end_jd_tdb]` for a NAIF target ID.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TargetCoverage {
    pub target_id: i32,
    pub center_id: i32,
    pub frame_id: i32,
    pub spk_type: i32,
    pub start_jd_tdb: f64,
    pub end_jd_tdb: f64,
    pub segment_count: usize,
    pub has_gaps: bool,
    pub has_overlaps: bool,
}

/// Abstraction over binary SPK ephemeris providers (`ARCHITECTURE.md` §1.2).
pub trait SpkSource: Send + Sync {
    /// Queries the instantaneous geometric (`Aberration::NONE`) state vector of `target_naif_id`
    /// relative to `observer_naif_id` in ICRF (`J2000`) at `jd_tdb`.
    fn state_icrf_km(
        &self,
        target_naif_id: i32,
        observer_naif_id: i32,
        jd_tdb: f64,
    ) -> Result<SpkState, EphemerisError> {
        self.state_icrf_km_split(target_naif_id, observer_naif_id, JD_J2000, jd_tdb - JD_J2000)
    }

    /// Queries the instantaneous geometric (`Aberration::NONE`) state vector of `target_naif_id`
    /// relative to `observer_naif_id` in ICRF (`J2000`) at two-part `(jd_tdb_1, jd_tdb_2)`,
    /// preserving full sub-nanosecond precision.
    fn state_icrf_km_split(
        &self,
        target_naif_id: i32,
        observer_naif_id: i32,
        jd_tdb_1: f64,
        jd_tdb_2: f64,
    ) -> Result<SpkState, EphemerisError>;

    /// Returns the merged continuous coverage metadata for `target_naif_id`.
    fn target_coverage(&self, target_naif_id: i32) -> Result<&TargetCoverage, EphemerisError>;

    /// Returns all parsed DAF/SPK segment descriptors loaded in this source.
    fn segments(&self) -> &[DafSegmentDescriptor];
}

struct Type21SegmentData {
    start_et_s: f64,
    end_et_s: f64,
    center_id: i32,
    words: Vec<f64>,
}

/// Production [`SpkSource`] backed by [`anise::almanac::Almanac`] and pre-validated DAF segment headers.
/// Also includes a clean-room NAIF `SPK Type 21` (`Extended Modified Difference Array`) evaluator for
/// raw JPL Horizons small-body kernels that `anise` 0.9.4 cannot evaluate natively.
pub struct AniseSpkBackend {
    almanac: Almanac,
    segments: Vec<DafSegmentDescriptor>,
    coverage_by_target: BTreeMap<i32, TargetCoverage>,
    type21_by_target: BTreeMap<i32, Vec<Type21SegmentData>>,
}

assert_impl_all!(AniseSpkBackend: Send, Sync);

impl AniseSpkBackend {
    /// Creates an [`AniseSpkBackend`] from one or more `.bsp` file paths on disk.
    pub fn from_bsp_paths<P: AsRef<Path>>(paths: &[P]) -> Result<Self, EphemerisError> {
        if paths.is_empty() {
            return Err(EphemerisError::MissingData {
                detail: "No BSP file paths provided to AniseSpkBackend::from_bsp_paths".to_owned(),
            });
        }
        let mut named_buffers = Vec::with_capacity(paths.len());
        for p in paths {
            let path_ref = p.as_ref();
            let label = path_ref.display().to_string();
            let bytes = fs::read(path_ref).map_err(|e| EphemerisError::MissingData {
                detail: format!("Failed to read BSP file '{label}': {e}"),
            })?;
            named_buffers.push((label, bytes));
        }
        Self::from_named_bytes(&named_buffers)
    }

    /// Creates an [`AniseSpkBackend`] from in-memory `(label, bytes)` DAF/SPK buffers.
    pub fn from_named_bytes(buffers: &[(String, Vec<u8>)]) -> Result<Self, EphemerisError> {
        if buffers.is_empty() {
            return Err(EphemerisError::MissingData {
                detail: "No BSP buffers provided to AniseSpkBackend".to_owned(),
            });
        }

        let mut almanac = Almanac::default();
        let mut all_segments = Vec::new();
        let mut grouped: BTreeMap<i32, Vec<DafSegmentDescriptor>> = BTreeMap::new();
        let mut type21_by_target: BTreeMap<i32, Vec<Type21SegmentData>> = BTreeMap::new();

        for (label, bytes) in buffers {
            let daf = parse_daf_spk(bytes, label)?;
            for seg in &daf.segments {
                if seg.spk_type != 2
                    && seg.spk_type != 3
                    && seg.spk_type != 9
                    && seg.spk_type != 13
                    && seg.spk_type != 21
                {
                    return Err(EphemerisError::UnsupportedSpkType {
                        target_id: seg.target_id,
                        spk_type: seg.spk_type,
                    });
                }
                if seg.spk_type == 21 {
                    let words = extract_segment_words(bytes, seg, daf.little_endian, label)?;
                    type21_by_target
                        .entry(seg.target_id)
                        .or_default()
                        .push(Type21SegmentData {
                            start_et_s: seg.start_et_s,
                            end_et_s: seg.end_et_s,
                            center_id: seg.center_id,
                            words,
                        });
                }
                grouped
                    .entry(seg.target_id)
                    .or_default()
                    .push(seg.clone());
                all_segments.push(seg.clone());
            }

            let spk = SPK::parse(bytes.clone()).map_err(|e| EphemerisError::CorruptKernel {
                path: label.clone(),
                reason: format!("ANISE SPK::parse failed: {e}"),
            })?;
            almanac = almanac.with_spk(spk);
        }

        for segs in type21_by_target.values_mut() {
            segs.sort_by(|a, b| {
                a.start_et_s
                    .partial_cmp(&b.start_et_s)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }

        let mut coverage_by_target = BTreeMap::new();
        for (target_id, mut segs) in grouped {
            segs.sort_by(|a, b| {
                a.start_et_s
                    .partial_cmp(&b.start_et_s)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let first = segs.first().ok_or_else(|| EphemerisError::MissingData {
                detail: format!("Empty segment list for target {target_id}"),
            })?;
            let center_id = first.center_id;
            let frame_id = first.frame_id;
            let spk_type = first.spk_type;
            let start_jd_tdb = first.start_jd_tdb;
            let mut end_jd_tdb = first.end_jd_tdb;
            let mut prev_end_et = first.end_et_s;
            let mut has_gaps = false;
            let mut has_overlaps = false;

            for seg in segs.iter().skip(1) {
                let diff = seg.start_et_s - prev_end_et;
                if diff > 1e-3 {
                    has_gaps = true;
                } else if diff < -1e-3 {
                    has_overlaps = true;
                }
                if seg.end_et_s > prev_end_et {
                    prev_end_et = seg.end_et_s;
                    end_jd_tdb = seg.end_jd_tdb;
                }
            }

            coverage_by_target.insert(
                target_id,
                TargetCoverage {
                    target_id,
                    center_id,
                    frame_id,
                    spk_type,
                    start_jd_tdb,
                    end_jd_tdb,
                    segment_count: segs.len(),
                    has_gaps,
                    has_overlaps,
                },
            );
        }

        Ok(Self {
            almanac,
            segments: all_segments,
            coverage_by_target,
            type21_by_target,
        })
    }

    /// Returns a reference to the underlying [`anise::almanac::Almanac`].
    #[must_use]
    pub const fn almanac(&self) -> &Almanac {
        &self.almanac
    }
}

impl SpkSource for AniseSpkBackend {
    fn state_icrf_km_split(
        &self,
        target_naif_id: i32,
        observer_naif_id: i32,
        jd_tdb_1: f64,
        jd_tdb_2: f64,
    ) -> Result<SpkState, EphemerisError> {
        let jd_tdb = jd_tdb_1 + jd_tdb_2;
        let et_s = ((jd_tdb_1 - JD_J2000) + jd_tdb_2) * SECONDS_PER_DAY;
        if !jd_tdb_1.is_finite() || !jd_tdb_2.is_finite() || !jd_tdb.is_finite() || !et_s.is_finite()
        {
            return Err(EphemerisError::InvalidInput {
                message: format!("Non-finite TDB epoch: jd1={jd_tdb_1}, jd2={jd_tdb_2}"),
            });
        }

        if target_naif_id != 0 {
            let cov = self.target_coverage(target_naif_id)?;
            if jd_tdb < cov.start_jd_tdb || jd_tdb > cov.end_jd_tdb {
                return Err(EphemerisError::OutOfCoverage {
                    body: format!("NAIF:{target_naif_id}"),
                    epoch_jd_tdb: jd_tdb,
                    start_jd_tdb: cov.start_jd_tdb,
                    end_jd_tdb: cov.end_jd_tdb,
                });
            }
        }
        if observer_naif_id != 0 {
            let cov = self.target_coverage(observer_naif_id)?;
            if jd_tdb < cov.start_jd_tdb || jd_tdb > cov.end_jd_tdb {
                return Err(EphemerisError::OutOfCoverage {
                    body: format!("NAIF:{observer_naif_id}"),
                    epoch_jd_tdb: jd_tdb,
                    start_jd_tdb: cov.start_jd_tdb,
                    end_jd_tdb: cov.end_jd_tdb,
                });
            }
        }

        // If target_naif_id is a raw SPK Type 21 segment (not supported natively by ANISE 0.9.4),
        // evaluate it via our clean-room NAIF SPKE21 evaluator and add the parent center state.
        if let Some(t21_segs) = self.type21_by_target.get(&target_naif_id) {
            let idx = t21_segs
                .partition_point(|s| s.end_et_s < et_s)
                .min(t21_segs.len().saturating_sub(1));
            let seg = t21_segs
                .get(idx)
                .ok_or_else(|| EphemerisError::MissingData {
                    detail: format!("Missing Type 21 segment for target {target_naif_id}"),
                })?;
            let (rel_pos, rel_vel) = eval_spk_type21(&seg.words, et_s, target_naif_id)?;
            let center_id = seg.center_id;
            if center_id == observer_naif_id {
                return Ok(SpkState {
                    pos_km: rel_pos,
                    vel_km_s: rel_vel,
                    target_id: target_naif_id,
                    observer_id: observer_naif_id,
                });
            }
            let center_state =
                self.state_icrf_km_split(center_id, observer_naif_id, jd_tdb_1, jd_tdb_2)?;
            return Ok(SpkState {
                pos_km: [
                    rel_pos[0] + center_state.pos_km[0],
                    rel_pos[1] + center_state.pos_km[1],
                    rel_pos[2] + center_state.pos_km[2],
                ],
                vel_km_s: [
                    rel_vel[0] + center_state.vel_km_s[0],
                    rel_vel[1] + center_state.vel_km_s[1],
                    rel_vel[2] + center_state.vel_km_s[2],
                ],
                target_id: target_naif_id,
                observer_id: observer_naif_id,
            });
        }

        let epoch = Epoch::from_et_seconds(et_s);
        let target_frame = Frame::new(target_naif_id, J2000);
        let observer_frame = Frame::new(observer_naif_id, J2000);

        let state = self
            .almanac
            .translate(target_frame, observer_frame, epoch, Aberration::NONE)
            .map_err(|e| {
                let msg = e.to_string();
                if msg.contains("data type") || msg.contains("Type21") {
                    if let Ok(cov) = self.target_coverage(target_naif_id) {
                        return EphemerisError::UnsupportedSpkType {
                            target_id: target_naif_id,
                            spk_type: cov.spk_type,
                        };
                    }
                }
                EphemerisError::CorruptKernel {
                    path: format!("NAIF:{target_naif_id}->{observer_naif_id}"),
                    reason: msg,
                }
            })?;

        Ok(SpkState {
            pos_km: [state.radius_km.x, state.radius_km.y, state.radius_km.z],
            vel_km_s: [
                state.velocity_km_s.x,
                state.velocity_km_s.y,
                state.velocity_km_s.z,
            ],
            target_id: target_naif_id,
            observer_id: observer_naif_id,
        })
    }

    fn target_coverage(&self, target_naif_id: i32) -> Result<&TargetCoverage, EphemerisError> {
        self.coverage_by_target
            .get(&target_naif_id)
            .ok_or_else(|| EphemerisError::UnsupportedBody {
                body: format!("NAIF:{target_naif_id}"),
            })
    }

    fn segments(&self) -> &[DafSegmentDescriptor] {
        &self.segments
    }
}
