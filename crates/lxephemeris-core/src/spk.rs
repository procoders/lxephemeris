//! `SpkSource` trait and `AniseSpkBackend` (`ARCHITECTURE.md` §1.2, §2.2.1).
//!
//! Isolates all DAF/SPK kernel access behind the internal [`SpkSource`] trait and statically
//! verifies `Send + Sync` on [`anise::almanac::Almanac`].

use crate::daf::{
    eval_spk_type21, extract_segment_words, parse_daf_spk, parse_daf_spk_file_headers,
    DafSegmentDescriptor, DafSpkFile,
};
use crate::error::EphemerisError;
use crate::types::{AU_KM, JD_J2000, SECONDS_PER_DAY};
use anise::almanac::Almanac;
use anise::astro::Aberration;
use anise::constants::orientations::J2000;
use anise::naif::SPK;
use anise::prelude::{Epoch, Frame};
use bytes::{BufMut as _, BytesMut};
use serde::{Deserialize, Serialize};
use static_assertions::assert_impl_all;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::marker::PhantomData;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

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

struct LazyMinorKernel {
    path: PathBuf,
    label: String,
    targets: BTreeSet<i32>,
    segments: Vec<DafSegmentDescriptor>,
    almanac: OnceLock<Result<Almanac, EphemerisError>>,
}

impl LazyMinorKernel {
    fn get_or_init_almanac(&self) -> Result<&Almanac, EphemerisError> {
        let res = self.almanac.get_or_init(|| {
            let buf = read_file_to_bytes_mut(&self.path, &self.label)?;
            let _ = parse_daf_spk(&buf, &self.label)?;
            let spk = build_spk_from_bytes_mut(buf, &self.label)?;
            Ok(Almanac::default().with_spk_as(spk, Some(self.label.clone())))
        });
        match res {
            Ok(alm) => Ok(alm),
            Err(err) => Err(err.clone()),
        }
    }
}

fn validate_spk_segment_type(seg: &DafSegmentDescriptor) -> Result<(), EphemerisError> {
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
    Ok(())
}

fn read_file_to_bytes_mut(path: &Path, label: &str) -> Result<BytesMut, EphemerisError> {
    let mut file = fs::File::open(path).map_err(|e| EphemerisError::MissingData {
        detail: format!("Failed to read BSP file '{label}': {e}"),
    })?;
    let file_len_u64 = file
        .metadata()
        .map_err(|e| EphemerisError::MissingData {
            detail: format!("Failed to stat BSP file '{label}': {e}"),
        })?
        .len();
    let file_len = usize::try_from(file_len_u64).map_err(|_| EphemerisError::CorruptKernel {
        path: label.to_owned(),
        reason: format!("File size {file_len_u64} exceeds addressable usize"),
    })?;
    let mut writer = BytesMut::with_capacity(file_len).writer();
    std::io::copy(&mut file, &mut writer).map_err(|e| EphemerisError::MissingData {
        detail: format!("Failed to read BSP file '{label}': {e}"),
    })?;
    Ok(writer.into_inner())
}

fn build_spk_from_bytes_mut(buf: BytesMut, label: &str) -> Result<SPK, EphemerisError> {
    let spk = SPK {
        bytes: buf,
        crc32: None,
        _daf_type: PhantomData,
    };
    spk.file_record().map_err(|e| EphemerisError::CorruptKernel {
        path: label.to_owned(),
        reason: format!("ANISE SPK::parse failed: {e}"),
    })?;
    spk.name_record(None)
        .map_err(|e| EphemerisError::CorruptKernel {
            path: label.to_owned(),
            reason: format!("ANISE SPK::parse failed: {e}"),
        })?;
    Ok(spk)
}

fn build_coverage_by_target(
    grouped: BTreeMap<i32, Vec<DafSegmentDescriptor>>,
) -> Result<BTreeMap<i32, TargetCoverage>, EphemerisError> {
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
    Ok(coverage_by_target)
}

/// Production [`SpkSource`] backed by [`anise::almanac::Almanac`] and pre-validated DAF segment headers.
/// Also includes a clean-room NAIF `SPK Type 21` (`Extended Modified Difference Array`) evaluator for
/// raw JPL Horizons small-body kernels that `anise` 0.9.4 cannot evaluate natively, and lazy `OnceLock`
/// loading for minor-body-only `.bsp` bundles.
pub struct AniseSpkBackend {
    almanac: Almanac,
    segments: Vec<DafSegmentDescriptor>,
    coverage_by_target: BTreeMap<i32, TargetCoverage>,
    type21_by_target: BTreeMap<i32, Vec<Type21SegmentData>>,
    lazy_kernels: Vec<LazyMinorKernel>,
    lazy_targets: BTreeSet<i32>,
}

assert_impl_all!(AniseSpkBackend: Send, Sync);

struct LazyParentTranslation {
    pos_km: [f64; 3],
    vel_km_s: [f64; 3],
    center_id: i32,
}

impl AniseSpkBackend {
    /// Creates an [`AniseSpkBackend`] from one or more `.bsp` file paths on disk.
    ///
    /// Inspects DAF headers first (~12 KB I/O) so segment metadata and target coverage are populated
    /// immediately; streams primary kernels directly into a single `BytesMut` allocation without
    /// intermediate `Vec<u8>` duplication; and defers loading minor-body-only `.bsp` payloads until
    /// a minor body is actually queried.
    pub fn from_bsp_paths<P: AsRef<Path>>(paths: &[P]) -> Result<Self, EphemerisError> {
        if paths.is_empty() {
            return Err(EphemerisError::MissingData {
                detail: "No BSP file paths provided to AniseSpkBackend::from_bsp_paths".to_owned(),
            });
        }

        struct InspectedKernel {
            path: PathBuf,
            label: String,
            daf: DafSpkFile,
            is_minor_candidate: bool,
        }

        let mut inspected = Vec::with_capacity(paths.len());

        for p in paths {
            let raw_path = p.as_ref();
            let label = raw_path.display().to_string();
            let daf = parse_daf_spk_file_headers(raw_path)?;
            let path = fs::canonicalize(raw_path).unwrap_or_else(|_| {
                if raw_path.is_absolute() {
                    raw_path.to_path_buf()
                } else {
                    std::env::current_dir()
                        .map(|cwd| cwd.join(raw_path))
                        .unwrap_or_else(|_| raw_path.to_path_buf())
                }
            });
            for seg in &daf.segments {
                validate_spk_segment_type(seg)?;
            }
            let is_minor_candidate = !daf.segments.is_empty()
                && daf
                    .segments
                    .iter()
                    .all(|seg| seg.target_id > 2_000_000 && seg.spk_type != 21);
            inspected.push(InspectedKernel {
                path,
                label,
                daf,
                is_minor_candidate,
            });
        }

        // Fixed-point closure for eager vs. lazy partitioning:
        // 1. Any non-minor-candidate kernel starts eager (or all kernels if no primary kernel is present).
        // 2. Any minor-candidate kernel that shares a `target_id` with an eager kernel, OR defines a
        //    `target_id` used as a `center_id` by an eager kernel, is transitively promoted to eager.
        // This guarantees every NAIF target is either 100% eager (preserving exact SPK override order)
        // or 100% lazy, and no eager body ever depends on a lazy intermediate center.
        let has_non_minor_kernel = inspected.iter().any(|item| !item.is_minor_candidate);
        let mut is_eager: Vec<bool> = inspected
            .iter()
            .map(|item| !has_non_minor_kernel || !item.is_minor_candidate)
            .collect();

        if has_non_minor_kernel {
            loop {
                let mut eager_targets = BTreeSet::new();
                let mut eager_centers = BTreeSet::new();
                for (idx, item) in inspected.iter().enumerate() {
                    if is_eager.get(idx).copied().unwrap_or(true) {
                        for seg in &item.daf.segments {
                            eager_targets.insert(seg.target_id);
                            eager_centers.insert(seg.center_id);
                        }
                    }
                }

                let mut changed = false;
                for (idx, item) in inspected.iter().enumerate() {
                    if !is_eager.get(idx).copied().unwrap_or(true) {
                        let must_promote = item.daf.segments.iter().any(|seg| {
                            eager_targets.contains(&seg.target_id)
                                || eager_centers.contains(&seg.target_id)
                        });
                        if must_promote {
                            if let Some(slot) = is_eager.get_mut(idx) {
                                *slot = true;
                                changed = true;
                            }
                        }
                    }
                }
                if !changed {
                    break;
                }
            }
        }

        let mut almanac = Almanac::default();
        let mut all_segments = Vec::new();
        let mut grouped: BTreeMap<i32, Vec<DafSegmentDescriptor>> = BTreeMap::new();
        let mut type21_by_target: BTreeMap<i32, Vec<Type21SegmentData>> = BTreeMap::new();
        let mut lazy_kernels = Vec::new();
        let mut lazy_targets = BTreeSet::new();

        for (idx, item) in inspected.into_iter().enumerate() {
            let load_lazily = !is_eager.get(idx).copied().unwrap_or(true);

            if load_lazily {
                let mut targets = BTreeSet::new();
                let mut kernel_segments = Vec::with_capacity(item.daf.segments.len());
                for seg in item.daf.segments {
                    targets.insert(seg.target_id);
                    lazy_targets.insert(seg.target_id);
                    grouped
                        .entry(seg.target_id)
                        .or_default()
                        .push(seg.clone());
                    kernel_segments.push(seg.clone());
                    all_segments.push(seg);
                }
                lazy_kernels.push(LazyMinorKernel {
                    path: item.path,
                    label: item.label,
                    targets,
                    segments: kernel_segments,
                    almanac: OnceLock::new(),
                });
            } else {
                let buf = read_file_to_bytes_mut(&item.path, &item.label)?;
                let daf = parse_daf_spk(&buf, &item.label)?;
                for seg in &daf.segments {
                    validate_spk_segment_type(seg)?;
                    if seg.spk_type == 21 {
                        let words =
                            extract_segment_words(&buf, seg, daf.little_endian, &item.label)?;
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
                let spk = build_spk_from_bytes_mut(buf, &item.label)?;
                almanac = almanac.with_spk_as(spk, Some(format!("{idx}:{}", item.label)));
            }
        }

        for segs in type21_by_target.values_mut() {
            segs.sort_by(|a, b| {
                a.start_et_s
                    .partial_cmp(&b.start_et_s)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }

        let coverage_by_target = build_coverage_by_target(grouped)?;

        Ok(Self {
            almanac,
            segments: all_segments,
            coverage_by_target,
            type21_by_target,
            lazy_kernels,
            lazy_targets,
        })
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

        for (idx, (label, bytes)) in buffers.iter().enumerate() {
            let daf = parse_daf_spk(bytes, label)?;
            for seg in &daf.segments {
                validate_spk_segment_type(seg)?;
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

            let buf = BytesMut::from(bytes.as_slice());
            let spk = build_spk_from_bytes_mut(buf, label)?;
            almanac = almanac.with_spk_as(spk, Some(format!("{idx}:{label}")));
        }

        for segs in type21_by_target.values_mut() {
            segs.sort_by(|a, b| {
                a.start_et_s
                    .partial_cmp(&b.start_et_s)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }

        let coverage_by_target = build_coverage_by_target(grouped)?;

        Ok(Self {
            almanac,
            segments: all_segments,
            coverage_by_target,
            type21_by_target,
            lazy_kernels: Vec::new(),
            lazy_targets: BTreeSet::new(),
        })
    }

    /// Returns a reference to the underlying [`anise::almanac::Almanac`].
    #[must_use]
    pub const fn almanac(&self) -> &Almanac {
        &self.almanac
    }

    /// ANISE 0.9.4 SPK summary epoch selector tolerance (`100 ns` = `1e-7 s`).
    const ANISE_EPOCH_TOL_S: f64 = 1e-7;

    fn segment_covers_et(seg: &DafSegmentDescriptor, et_s: f64) -> bool {
        (seg.start_et_s - Self::ANISE_EPOCH_TOL_S) <= et_s
            && et_s <= (seg.end_et_s + Self::ANISE_EPOCH_TOL_S)
    }

    fn active_center_at(&self, id: i32, et_s: f64) -> Option<i32> {
        if let Some(seg) = self
            .segments
            .iter()
            .rev()
            .find(|s| s.target_id == id && s.start_et_s <= et_s && et_s <= s.end_et_s)
        {
            return Some(seg.center_id);
        }
        if let Some(seg) = self
            .segments
            .iter()
            .rev()
            .find(|s| s.target_id == id && Self::segment_covers_et(s, et_s))
        {
            return Some(seg.center_id);
        }
        self.coverage_by_target.get(&id).map(|cov| cov.center_id)
    }

    fn ancestor_chain(&self, id: i32, et_s: f64) -> Vec<i32> {
        let mut chain = Vec::with_capacity(8);
        chain.push(id);
        let mut cur = id;
        for _ in 0..16 {
            if cur == 0 {
                break;
            }
            if let Some(center_id) = self.active_center_at(cur, et_s) {
                if chain.contains(&center_id) {
                    break;
                }
                chain.push(center_id);
                cur = center_id;
            } else {
                break;
            }
        }
        chain
    }

    fn find_common_ancestor(
        &self,
        target_id: i32,
        observer_id: i32,
        et_s: f64,
    ) -> Result<i32, EphemerisError> {
        let target_chain = self.ancestor_chain(target_id, et_s);
        let observer_chain = self.ancestor_chain(observer_id, et_s);
        for node in &target_chain {
            if observer_chain.contains(node) {
                return Ok(*node);
            }
        }
        Err(EphemerisError::CorruptKernel {
            path: format!("NAIF:{target_id}->{observer_id}"),
            reason: format!(
                "No common SPK ephemeris ancestor between NAIF {target_id} and {observer_id}"
            ),
        })
    }

    fn eval_lazy_target_if_any(
        &self,
        id: i32,
        jd_tdb: f64,
        et_s: f64,
        epoch: Epoch,
    ) -> Result<Option<LazyParentTranslation>, EphemerisError> {
        if !self.lazy_targets.contains(&id) {
            return Ok(None);
        }

        // Resolve applicable coverage from headers BEFORE any disk I/O so gap queries
        // return OutOfCoverage deterministically without opening or initializing lazy files.
        let any_lazy_covers = self.lazy_kernels.iter().any(|kernel| {
            kernel.targets.contains(&id)
                && kernel
                    .segments
                    .iter()
                    .any(|s| s.target_id == id && Self::segment_covers_et(s, et_s))
        });
        if !any_lazy_covers {
            if let Some(cov) = self.coverage_by_target.get(&id) {
                return Err(EphemerisError::OutOfCoverage {
                    body: format!("NAIF:{id}"),
                    epoch_jd_tdb: jd_tdb,
                    start_jd_tdb: cov.start_jd_tdb,
                    end_jd_tdb: cov.end_jd_tdb,
                });
            }
            return Err(EphemerisError::UnsupportedBody {
                body: format!("NAIF:{id}"),
            });
        }

        let mut last_err: Option<EphemerisError> = None;
        for kernel in self.lazy_kernels.iter().rev() {
            if !kernel.targets.contains(&id) {
                continue;
            }
            let covers_epoch = kernel
                .segments
                .iter()
                .any(|s| s.target_id == id && Self::segment_covers_et(s, et_s));
            if !covers_epoch {
                continue;
            }
            let alm = kernel.get_or_init_almanac()?;
            match alm.translate_to_parent(Frame::new(id, J2000), epoch) {
                Ok(state) => {
                    return Ok(Some(LazyParentTranslation {
                        pos_km: [state.radius_km.x, state.radius_km.y, state.radius_km.z],
                        vel_km_s: [
                            state.velocity_km_s.x,
                            state.velocity_km_s.y,
                            state.velocity_km_s.z,
                        ],
                        center_id: state.frame.ephemeris_id,
                    }));
                }
                Err(e) => {
                    last_err = Some(EphemerisError::CorruptKernel {
                        path: format!("NAIF:{id}"),
                        reason: e.to_string(),
                    });
                }
            }
        }

        let strictly_covered = self
            .segments
            .iter()
            .any(|s| s.target_id == id && s.start_et_s <= et_s && et_s <= s.end_et_s);
        if !strictly_covered {
            if let Some(cov) = self.coverage_by_target.get(&id) {
                return Err(EphemerisError::OutOfCoverage {
                    body: format!("NAIF:{id}"),
                    epoch_jd_tdb: jd_tdb,
                    start_jd_tdb: cov.start_jd_tdb,
                    end_jd_tdb: cov.end_jd_tdb,
                });
            }
        }

        Err(last_err.unwrap_or_else(|| EphemerisError::UnsupportedBody {
            body: format!("NAIF:{id}"),
        }))
    }

    fn state_to_ancestor(
        &self,
        id: i32,
        ancestor: i32,
        jd_tdb_1: f64,
        jd_tdb_2: f64,
        et_s: f64,
        epoch: Epoch,
    ) -> Result<([f64; 3], [f64; 3]), EphemerisError> {
        if id == ancestor {
            return Ok(([0.0; 3], [0.0; 3]));
        }
        let jd_tdb = jd_tdb_1 + jd_tdb_2;
        if let Some(rel) = self.eval_lazy_target_if_any(id, jd_tdb, et_s, epoch)? {
            if rel.center_id == ancestor {
                return Ok((rel.pos_km, rel.vel_km_s));
            }
            let parent = self.state_icrf_km_split(rel.center_id, ancestor, jd_tdb_1, jd_tdb_2)?;
            return Ok((
                [
                    rel.pos_km[0] + parent.pos_km[0],
                    rel.pos_km[1] + parent.pos_km[1],
                    rel.pos_km[2] + parent.pos_km[2],
                ],
                [
                    rel.vel_km_s[0] + parent.vel_km_s[0],
                    rel.vel_km_s[1] + parent.vel_km_s[1],
                    rel.vel_km_s[2] + parent.vel_km_s[2],
                ],
            ));
        }
        let st = self.state_icrf_km_split(id, ancestor, jd_tdb_1, jd_tdb_2)?;
        Ok((st.pos_km, st.vel_km_s))
    }

    #[cfg(test)]
    fn any_lazy_kernel_loaded(&self) -> bool {
        self.lazy_kernels
            .iter()
            .any(|k| k.almanac.get().is_some())
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

        if self.lazy_targets.contains(&target_naif_id)
            || self.lazy_targets.contains(&observer_naif_id)
        {
            if target_naif_id == observer_naif_id {
                return Ok(SpkState {
                    pos_km: [0.0; 3],
                    vel_km_s: [0.0; 3],
                    target_id: target_naif_id,
                    observer_id: observer_naif_id,
                });
            }
            let epoch = Epoch::from_et_seconds(et_s);
            let common_node = self.find_common_ancestor(target_naif_id, observer_naif_id, et_s)?;
            let (pos_bwrd, vel_bwrd) = self.state_to_ancestor(
                target_naif_id,
                common_node,
                jd_tdb_1,
                jd_tdb_2,
                et_s,
                epoch,
            )?;
            let (pos_fwrd, vel_fwrd) = self.state_to_ancestor(
                observer_naif_id,
                common_node,
                jd_tdb_1,
                jd_tdb_2,
                et_s,
                epoch,
            )?;
            return Ok(SpkState {
                pos_km: [
                    pos_bwrd[0] - pos_fwrd[0],
                    pos_bwrd[1] - pos_fwrd[1],
                    pos_bwrd[2] - pos_fwrd[2],
                ],
                vel_km_s: [
                    vel_bwrd[0] - vel_fwrd[0],
                    vel_bwrd[1] - vel_fwrd[1],
                    vel_bwrd[2] - vel_fwrd[2],
                ],
                target_id: target_naif_id,
                observer_id: observer_naif_id,
            });
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{Engine, EngineConfig};
    use crate::types::{
        BodyId, MpcNumber, PlanetId, JD_1850_01_01, JD_2100_01_01, MINOR_BODIES_68,
    };

    #[test]
    fn test_lazy_minor_bsp_coverage_and_bit_identical_states() {
        let workspace_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(|p| p.parent())
            .unwrap_or_else(|| Path::new("."));
        let data_dir = workspace_dir.join("data");
        let de440s_path = data_dir.join("de440s.bsp");
        let ast_path = data_dir.join("asteroids_68_1850_2100.bsp");
        if !de440s_path.exists() || !ast_path.exists() {
            return;
        }

        let de_bytes = fs::read(&de440s_path).expect("read de440s.bsp");
        let ast_bytes = fs::read(&ast_path).expect("read asteroids_68_1850_2100.bsp");

        // 1. Verify header-only DAF parser matches full-buffer DAF parser identically.
        let de_hdr = parse_daf_spk_file_headers(&de440s_path).expect("de440s headers");
        let de_full =
            parse_daf_spk(&de_bytes, &de440s_path.display().to_string()).expect("de440s full");
        assert_eq!(de_hdr, de_full);

        let ast_hdr = parse_daf_spk_file_headers(&ast_path).expect("asteroids headers");
        let ast_full =
            parse_daf_spk(&ast_bytes, &ast_path.display().to_string()).expect("asteroids full");
        assert_eq!(ast_hdr, ast_full);

        // 2. Construct lazy backend (from_bsp_paths & Engine::from_data_dir) vs eager backend (from_named_bytes).
        let named_buffers = vec![
            ("de440s.bsp".to_owned(), de_bytes),
            ("asteroids_68_1850_2100.bsp".to_owned(), ast_bytes),
        ];
        let eager_backend =
            AniseSpkBackend::from_named_bytes(&named_buffers).expect("eager from_named_bytes");
        let lazy_backend = AniseSpkBackend::from_bsp_paths(&[&de440s_path, &ast_path])
            .expect("lazy from_bsp_paths");

        assert_eq!(lazy_backend.lazy_kernels.len(), 1);
        assert!(
            !lazy_backend.any_lazy_kernel_loaded(),
            "minor-body .bsp payload must not be loaded before a minor body is queried"
        );

        // Segments and target_coverage must be 100% populated before minor-body evaluation.
        assert_eq!(lazy_backend.segments(), eager_backend.segments());
        for seg in eager_backend.segments() {
            assert_eq!(
                lazy_backend.target_coverage(seg.target_id),
                eager_backend.target_coverage(seg.target_id)
            );
        }

        // Engine::from_data_dir coverage & capabilities must match Engine::from_bytes BEFORE minor-body evaluation.
        let cfg = EngineConfig {
            pin_time_table_version: Some("embedded".to_owned()),
            ..EngineConfig::default()
        };
        let eng_dir = Engine::from_data_dir(&data_dir, cfg.clone()).expect("Engine::from_data_dir");
        let eng_bytes = Engine::from_bytes(&named_buffers, cfg).expect("Engine::from_bytes");

        assert_eq!(
            eng_dir.capabilities().loaded_segment_count,
            eng_bytes.capabilities().loaded_segment_count
        );
        assert_eq!(
            eng_dir.capabilities().minor_bodies_loaded_count,
            eng_bytes.capabilities().minor_bodies_loaded_count
        );
        for spec in MINOR_BODIES_68 {
            let body = BodyId::MinorPlanet(spec.mpc);
            assert_eq!(eng_dir.coverage(body), eng_bytes.coverage(body));
        }

        // Evaluating major planets must NOT trigger lazy minor-body kernel loading.
        for planet in [
            PlanetId::Sun,
            PlanetId::Moon,
            PlanetId::Earth,
            PlanetId::Saturn,
            PlanetId::Pluto,
        ] {
            let st_lazy = lazy_backend
                .state_icrf_km_split(planet.naif_target_id(), 399, JD_J2000, 0.25)
                .expect("lazy planet state");
            let st_eager = eager_backend
                .state_icrf_km_split(planet.naif_target_id(), 399, JD_J2000, 0.25)
                .expect("eager planet state");
            assert_eq!(st_lazy, st_eager);
        }
        assert!(
            !lazy_backend.any_lazy_kernel_loaded(),
            "major-planet queries must not load minor-body .bsp payload"
        );

        // 3. Evaluate Ceres (20000001) and Chiron (20002060) relative to Earth (399), Sun (10), and SSB (0)
        // (plus symmetric observer cases) and assert bit-identical SpkState.
        let ceres_id = MpcNumber(1).standard_naif_spk_id();
        let chiron_id = MpcNumber(2060).standard_naif_spk_id();
        assert_eq!(ceres_id, 20_000_001);
        assert_eq!(chiron_id, 20_002_060);

        let split_epochs = [
            (JD_1850_01_01, 10.25),
            (JD_J2000, 0.0),
            (JD_J2000, 3_652.123_456_789),
            (JD_2100_01_01, -10.75),
        ];
        for &(jd1, jd2) in &split_epochs {
            for target_id in [ceres_id, chiron_id] {
                for observer_id in [399, 10, 0] {
                    let st_lazy = lazy_backend
                        .state_icrf_km_split(target_id, observer_id, jd1, jd2)
                        .expect("lazy minor->observer state");
                    let st_eager = eager_backend
                        .state_icrf_km_split(target_id, observer_id, jd1, jd2)
                        .expect("eager minor->observer state");
                    assert_eq!(
                        st_lazy, st_eager,
                        "bit-identical mismatch for target={target_id}, obs={observer_id} at ({jd1}, {jd2})"
                    );

                    // Symmetric case: observer is minor body
                    let rev_lazy = lazy_backend
                        .state_icrf_km_split(observer_id, target_id, jd1, jd2)
                        .expect("lazy observer->minor state");
                    let rev_eager = eager_backend
                        .state_icrf_km_split(observer_id, target_id, jd1, jd2)
                        .expect("eager observer->minor state");
                    assert_eq!(
                        rev_lazy, rev_eager,
                        "bit-identical mismatch for target={observer_id}, obs={target_id} at ({jd1}, {jd2})"
                    );

                    // And via Engine::from_data_dir vs Engine::from_bytes
                    let eng_st_dir = eng_dir
                        .spk()
                        .state_icrf_km_split(target_id, observer_id, jd1, jd2)
                        .expect("eng_dir minor state");
                    assert_eq!(eng_st_dir, st_eager);
                }
            }

            // Minor-to-minor translation (Ceres relative to Chiron)
            let m2m_lazy = lazy_backend
                .state_icrf_km_split(ceres_id, chiron_id, jd1, jd2)
                .expect("lazy ceres->chiron");
            let m2m_eager = eager_backend
                .state_icrf_km_split(ceres_id, chiron_id, jd1, jd2)
                .expect("eager ceres->chiron");
            assert_eq!(m2m_lazy, m2m_eager);
        }

        assert!(
            lazy_backend.any_lazy_kernel_loaded(),
            "minor-body .bsp payload must be initialized on first minor-body query"
        );
    }

    fn const_type2_seg(
        name: &str,
        target_id: i32,
        center_id: i32,
        start_et_s: f64,
        end_et_s: f64,
        pos_km: [f64; 3],
    ) -> crate::daf::DafSegmentInput {
        let mid = 0.5 * (start_et_s + end_et_s);
        let radius = 0.5 * (end_et_s - start_et_s);
        let intlen = end_et_s - start_et_s;
        crate::daf::DafSegmentInput {
            name: name.to_owned(),
            target_id,
            center_id,
            frame_id: 1,
            spk_type: 2,
            start_et_s,
            end_et_s,
            data_words: vec![
                mid, radius, pos_km[0], 0.0, pos_km[1], 0.0, pos_km[2], 0.0, start_et_s, intlen,
                8.0, 1.0,
            ],
        }
    }

    #[test]
    fn test_lazy_spk_override_closure_parent_centers_epoch_fallback_and_truncated_names() {
        use crate::daf::build_daf_spk_bytes;

        let tmp_dir = std::env::temp_dir().join(format!("lx_spk_audit_{}", std::process::id()));
        fs::create_dir_all(&tmp_dir).expect("create temp dir");

        let x_id = 20_000_001;
        let y_id = 20_000_002;
        let day = SECONDS_PER_DAY;

        // 1. [P, B, A] transitive override closure:
        // P defines (10 -> 0) and (X -> 0); B defines (Y -> 0) = [100, 0, 0]; A defines (X -> 0) and (Y -> 0) = [999, 0, 0].
        // Because A overlaps X (eager), A becomes eager; fixed-point closure must also promote B so A overrides B on Y!
        let p_bytes = build_daf_spk_bytes(
            "P",
            "",
            &[
                const_type2_seg("SUN", 10, 0, 0.0, 100.0 * day, [10.0, 0.0, 0.0]),
                const_type2_seg("X_P", x_id, 0, 0.0, 100.0 * day, [50.0, 0.0, 0.0]),
            ],
        )
        .expect("build P");
        let b_bytes = build_daf_spk_bytes(
            "B",
            "",
            &[const_type2_seg(
                "Y_B",
                y_id,
                0,
                0.0,
                100.0 * day,
                [100.0, 0.0, 0.0],
            )],
        )
        .expect("build B");
        let a_bytes = build_daf_spk_bytes(
            "A",
            "",
            &[
                const_type2_seg("X_A", x_id, 0, 0.0, 100.0 * day, [55.0, 0.0, 0.0]),
                const_type2_seg("Y_A", y_id, 0, 0.0, 100.0 * day, [999.0, 0.0, 0.0]),
            ],
        )
        .expect("build A");

        let p_path = tmp_dir.join("p.bsp");
        let b_path = tmp_dir.join("b.bsp");
        let a_path = tmp_dir.join("a.bsp");
        fs::write(&p_path, &p_bytes).expect("write p");
        fs::write(&b_path, &b_bytes).expect("write b");
        fs::write(&a_path, &a_bytes).expect("write a");

        let pba_paths =
            AniseSpkBackend::from_bsp_paths(&[&p_path, &b_path, &a_path]).expect("pba paths");
        let pba_eager = AniseSpkBackend::from_named_bytes(&[
            ("p".to_owned(), p_bytes),
            ("b".to_owned(), b_bytes),
            ("a".to_owned(), a_bytes),
        ])
        .expect("pba eager");
        assert!(pba_paths.lazy_kernels.is_empty());
        let st_y_paths = pba_paths
            .state_icrf_km_split(y_id, 10, JD_J2000, 10.0)
            .expect("Y->10 paths");
        let st_y_eager = pba_eager
            .state_icrf_km_split(y_id, 10, JD_J2000, 10.0)
            .expect("Y->10 eager");
        assert_eq!(st_y_paths, st_y_eager);
        assert!((st_y_paths.pos_km[0] - (999.0 - 10.0)).abs() < 1e-9);

        // 2. [P_child, L_center]: eager body 301 has center_id = X where X lives in minor-only L_center.
        // Fixed-point closure over eager_centers must promote L_center so 301 -> 10 succeeds.
        let p_child_bytes = build_daf_spk_bytes(
            "P_CHILD",
            "",
            &[
                const_type2_seg("SUN", 10, 0, 0.0, 100.0 * day, [10.0, 0.0, 0.0]),
                const_type2_seg("PROBE", 301, x_id, 0.0, 100.0 * day, [7.0, 0.0, 0.0]),
            ],
        )
        .expect("build P_child");
        let l_center_bytes = build_daf_spk_bytes(
            "L_CENTER",
            "",
            &[const_type2_seg(
                "X_L",
                x_id,
                0,
                0.0,
                100.0 * day,
                [200.0, 0.0, 0.0],
            )],
        )
        .expect("build L_center");
        let pc_path = tmp_dir.join("p_child.bsp");
        let lc_path = tmp_dir.join("l_center.bsp");
        fs::write(&pc_path, &p_child_bytes).expect("write pc");
        fs::write(&lc_path, &l_center_bytes).expect("write lc");

        let pc_paths =
            AniseSpkBackend::from_bsp_paths(&[&pc_path, &lc_path]).expect("pc+lc from_bsp_paths");
        let pc_eager = AniseSpkBackend::from_named_bytes(&[
            ("pc".to_owned(), p_child_bytes),
            ("lc".to_owned(), l_center_bytes),
        ])
        .expect("pc+lc eager");
        let st_301_paths = pc_paths
            .state_icrf_km_split(301, 10, JD_J2000, 10.0)
            .expect("301->10 paths");
        let st_301_eager = pc_eager
            .state_icrf_km_split(301, 10, JD_J2000, 10.0)
            .expect("301->10 eager");
        assert_eq!(st_301_paths, st_301_eager);
        assert!((st_301_paths.pos_km[0] - (200.0 + 7.0 - 10.0)).abs() < 1e-9);

        // 3. Multi-lazy-kernel epoch-dependent center_id override & coverage gap:
        // Base defines (10 -> 0) over [0, 300d].
        // L1 defines (Y -> 0) over [0, 100d] = [300, 0, 0] and [200d, 300d] = [400, 0, 0] (gap at 100d..200d).
        // L2 defines (Y -> 10) over [220d, 300d] = [500, 0, 0] (different center_id = 10!).
        let base_bytes = build_daf_spk_bytes(
            "BASE",
            "",
            &[const_type2_seg(
                "SUN",
                10,
                0,
                0.0,
                300.0 * day,
                [25.0, 0.0, 0.0],
            )],
        )
        .expect("build base");
        let l1_bytes = build_daf_spk_bytes(
            "L1",
            "",
            &[
                const_type2_seg("Y_EARLY", y_id, 0, 0.0, 100.0 * day, [300.0, 0.0, 0.0]),
                const_type2_seg(
                    "Y_LATE_SSB",
                    y_id,
                    0,
                    200.0 * day,
                    300.0 * day,
                    [400.0, 0.0, 0.0],
                ),
            ],
        )
        .expect("build L1");
        let l2_bytes = build_daf_spk_bytes(
            "L2",
            "",
            &[const_type2_seg(
                "Y_LATE_HELIO",
                y_id,
                10,
                220.0 * day,
                300.0 * day,
                [500.0, 0.0, 0.0],
            )],
        )
        .expect("build L2");
        let base_path = tmp_dir.join("base.bsp");
        let l1_path = tmp_dir.join("l1.bsp");
        let l2_path = tmp_dir.join("l2.bsp");
        fs::write(&base_path, &base_bytes).expect("write base");
        fs::write(&l1_path, &l1_bytes).expect("write l1");
        fs::write(&l2_path, &l2_bytes).expect("write l2");

        let multi_lazy =
            AniseSpkBackend::from_bsp_paths(&[&base_path, &l1_path, &l2_path]).expect("multi lazy");
        let multi_eager = AniseSpkBackend::from_named_bytes(&[
            ("base".to_owned(), base_bytes.clone()),
            ("l1".to_owned(), l1_bytes.clone()),
            ("l2".to_owned(), l2_bytes.clone()),
        ])
        .expect("multi eager");
        assert_eq!(multi_lazy.lazy_kernels.len(), 2);
        assert!(
            multi_lazy
                .lazy_kernels
                .iter()
                .all(|k| k.path.is_absolute()),
            "lazy kernel paths must be canonicalized to absolute paths at construction"
        );

        // Sub-millisecond boundary gap (+100d + 0.0005s = 0.5 ms past L1's first segment end,
        // before ANY lazy kernel is loaded): must return OutOfCoverage WITHOUT disk I/O!
        assert!(matches!(
            multi_lazy.state_icrf_km_split(y_id, 10, JD_J2000, 100.0 + 0.0005 / day),
            Err(EphemerisError::OutOfCoverage { .. })
        ));
        assert!(
            !multi_lazy.any_lazy_kernel_loaded(),
            "0.5 ms gap query must reject from headers without initializing any lazy kernel"
        );

        // Also verify relative-path construction canonicalizes to an absolute path at construction
        // (thread-safe: no process-wide set_current_dir mutation in parallel test runner):
        let rel_dir_name = format!("../../target/tmp_lx_rel_{}", std::process::id());
        let rel_sub = Path::new(&rel_dir_name);
        fs::create_dir_all(rel_sub).expect("create relative target temp dir");
        let rel_base = rel_sub.join("base.bsp");
        let rel_l1 = rel_sub.join("l1.bsp");
        fs::write(&rel_base, &base_bytes).expect("write rel base");
        fs::write(&rel_l1, &l1_bytes).expect("write rel l1");
        assert!(!rel_l1.is_absolute());
        let rel_backend =
            AniseSpkBackend::from_bsp_paths(&[&rel_base, &rel_l1]).expect("construct relative");
        assert!(
            rel_backend
                .lazy_kernels
                .iter()
                .all(|k| k.path.is_absolute()),
            "relative paths passed to from_bsp_paths must be resolved to absolute paths eagerly"
        );
        assert!(!rel_backend.any_lazy_kernel_loaded());
        let rel_st = rel_backend
            .state_icrf_km_split(y_id, 10, JD_J2000, 50.0)
            .expect("query lazy relative-constructed backend");
        assert!((rel_st.pos_km[0] - (300.0 - 25.0)).abs() < 1e-9);
        let _ = fs::remove_dir_all(rel_sub);

        // At +50d: covered only by L1 (center 0); L2 must NOT even be loaded from disk!
        let st_50_lazy = multi_lazy
            .state_icrf_km_split(y_id, 10, JD_J2000, 50.0)
            .expect("Y->10 at +50d lazy");
        let st_50_eager = multi_eager
            .state_icrf_km_split(y_id, 10, JD_J2000, 50.0)
            .expect("Y->10 at +50d eager");
        assert_eq!(st_50_lazy, st_50_eager);
        assert!(
            multi_lazy
                .lazy_kernels
                .get(1)
                .is_some_and(|k| k.almanac.get().is_none()),
            "L2 must remain unloaded when querying epoch outside L2's segments"
        );

        // At +150d: inside outer [0, 300d] range, but falls in coverage gap (100d..200d) -> OutOfCoverage!
        assert!(matches!(
            multi_lazy.state_icrf_km_split(y_id, 10, JD_J2000, 150.0),
            Err(EphemerisError::OutOfCoverage { .. })
        ));

        // At +250d: covered by both L1 (center 0) and L2 (center 10); L2 wins with active center 10!
        assert_eq!(multi_lazy.active_center_at(y_id, 250.0 * day), Some(10));
        let st_250_lazy = multi_lazy
            .state_icrf_km_split(y_id, 10, JD_J2000, 250.0)
            .expect("Y->10 at +250d lazy");
        let st_250_eager = multi_eager
            .state_icrf_km_split(y_id, 10, JD_J2000, 250.0)
            .expect("Y->10 at +250d eager");
        assert_eq!(st_250_lazy, st_250_eager);
        assert!((st_250_lazy.pos_km[0] - 500.0).abs() < 1e-9);

        // 4. Truncated DAF name record rejection:
        // Build a DAF buffer with FWARD=2, NSUM=0 (no segment word references beyond EOF), and truncate
        // right after record 2 (2048 bytes, missing record 3 Name Record) or mid-record-3 (2500 bytes).
        let mut raw_no_name = l1_bytes.get(0..2048).expect("first 2 records").to_vec();
        if let Some(nsum_slot) = raw_no_name.get_mut(1024 + 16..1024 + 24) {
            nsum_slot.copy_from_slice(&0.0_f64.to_le_bytes());
        }
        let trunc_path = tmp_dir.join("truncated_name.bsp");
        fs::write(&trunc_path, &raw_no_name).expect("write truncated_name.bsp");
        assert!(matches!(
            parse_daf_spk(&raw_no_name, "truncated_name"),
            Err(EphemerisError::CorruptKernel { .. })
        ));
        assert!(matches!(
            parse_daf_spk_file_headers(&trunc_path),
            Err(EphemerisError::CorruptKernel { .. })
        ));

        let mut partial_name = raw_no_name;
        partial_name.extend_from_slice(&[b' '; 452]);
        let partial_path = tmp_dir.join("partial_name.bsp");
        fs::write(&partial_path, &partial_name).expect("write partial_name.bsp");
        assert!(matches!(
            parse_daf_spk(&partial_name, "partial_name"),
            Err(EphemerisError::CorruptKernel { .. })
        ));
        assert!(matches!(
            parse_daf_spk_file_headers(&partial_path),
            Err(EphemerisError::CorruptKernel { .. })
        ));

        let _ = fs::remove_dir_all(&tmp_dir);
    }
}
