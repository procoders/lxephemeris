//! Data bundle manifest (`manifest.json`) schema and validation (`PLAN.md` §0.2, `REVIEW_AND_SPEC_ASTRA.md` §3.2).

use crate::error::EphemerisError;
use crate::types::{JD_1850_01_01, JD_2100_01_01, MINOR_BODIES_68};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

/// Computes the lowercase hex-encoded SHA-256 digest of a byte slice.
#[must_use]
pub fn sha256_hex(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    let digest = hasher.finalize();
    let mut out = String::with_capacity(64);
    for b in digest {
        use std::fmt::Write;
        let _ = write!(&mut out, "{b:02x}");
    }
    out
}

/// Segment descriptor entry recorded in `manifest.json`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SegmentManifestEntry {
    pub name: String,
    pub target_id: i32,
    pub center_id: i32,
    pub frame_id: i32,
    pub spk_type: i32,
    pub start_et_s: f64,
    pub end_et_s: f64,
    pub start_jd_tdb: f64,
    pub end_jd_tdb: f64,
}

/// Operation-dependent half-open safe intervals `[start, end)` in TT (`ARCHITECTURE.md` §2.1).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SafeIntervalsManifest {
    pub planetary_and_minor_bodies_jd_tt: [f64; 2],
    pub planetary_and_minor_bodies_iso_tt: [String; 2],
    pub de440s_only_planetary_jd_tt: [f64; 2],
    pub de440s_only_planetary_iso_tt: [String; 2],
    pub interpolated_lilith_jd_tt: [f64; 2],
    pub interpolated_lilith_iso_tt: [String; 2],
}

/// Metadata for a base planetary SPK kernel (`de440s.bsp`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlanetaryKernelManifest {
    pub filename: String,
    pub source_url: String,
    pub sha256: String,
    pub size_bytes: u64,
    pub raw_interval_jd_tdb: [f64; 2],
    pub raw_interval_iso_tt: [String; 2],
    pub segment_count: usize,
    pub segments: Vec<SegmentManifestEntry>,
}

/// Per-body provenance and coverage entry for each of the 68 minor bodies (`REVIEW_AND_SPEC_ASTRA.md` §6).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MinorBodyManifestEntry {
    pub mpc: u32,
    pub name: String,
    pub category: String,
    pub naif_spk_id: i32,
    pub orbit_solution_id: String,
    pub orbit_solution_date: String,
    pub center_id: i32,
    pub frame_id: i32,
    pub spk_type: i32,
    pub segment_count: usize,
    pub has_gaps: bool,
    pub has_overlaps: bool,
    pub raw_interval_jd_tdb: [f64; 2],
    pub raw_interval_iso_tt: [String; 2],
    pub safe_interval_jd_tt: [f64; 2],
    pub safe_interval_iso_tt: [String; 2],
    pub source_spk_sha256: String,
    pub conversion_max_pos_err_m: Option<f64>,
    pub conversion_max_vel_err_mm_s: Option<f64>,
    pub segments: Vec<SegmentManifestEntry>,
}

/// Metadata for the stitched 68 minor bodies SPK bundle (`asteroids_68_1850_2100.bsp`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MinorBodiesKernelManifest {
    pub filename: String,
    pub source: String,
    pub sha256: String,
    pub size_bytes: u64,
    pub requested_interval_iso: [String; 2],
    pub raw_interval_jd_tdb: [f64; 2],
    pub raw_interval_iso_tt: [String; 2],
    pub safe_interval_jd_tt: [f64; 2],
    pub safe_interval_iso_tt: [String; 2],
    pub body_count: usize,
    pub bodies: Vec<MinorBodyManifestEntry>,
}

/// Top-level `manifest.json` structure for an `lxephemeris` data directory.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DataManifest {
    pub schema_version: String,
    pub builder_version: String,
    pub generated_at_utc: String,
    pub license_and_provenance: String,
    pub safe_intervals: SafeIntervalsManifest,
    pub planetary_kernel: PlanetaryKernelManifest,
    pub minor_bodies_kernel: MinorBodiesKernelManifest,
}

impl DataManifest {
    /// Strictly validates that the manifest contains all 68 required minor bodies without gaps or
    /// overlaps and covers the half-open public safe interval `[1850-01-01T00:00:00 TT, 2100-01-01T00:00:00 TT)`.
    pub fn validate(&self) -> Result<(), EphemerisError> {
        if self.minor_bodies_kernel.body_count != 68 || self.minor_bodies_kernel.bodies.len() != 68
        {
            return Err(EphemerisError::MissingData {
                detail: format!(
                    "Expected exactly 68 minor bodies in manifest, found {}",
                    self.minor_bodies_kernel.bodies.len()
                ),
            });
        }

        // Margin required for light-time (~0.5 day for outer TNOs), 4-point velocity stencil (+/- 2e-3 day),
        // and TT-TDB offset: at least 2 days on each side of [1850-01-01, 2100-01-01).
        let required_raw_start = JD_1850_01_01 - 2.0;
        let required_raw_end = JD_2100_01_01 + 2.0;

        for (expected, actual) in MINOR_BODIES_68
            .iter()
            .zip(self.minor_bodies_kernel.bodies.iter())
        {
            if actual.mpc != expected.mpc.0 {
                return Err(EphemerisError::InvalidInput {
                    message: format!(
                        "Manifest minor body ordering mismatch: expected MPC {} ({}), found MPC {} ({})",
                        expected.mpc.0, expected.name, actual.mpc, actual.name
                    ),
                });
            }
            if actual.has_gaps {
                return Err(EphemerisError::CorruptKernel {
                    path: self.minor_bodies_kernel.filename.clone(),
                    reason: format!("Gap detected in SPK segments for ({}) {}", actual.mpc, actual.name),
                });
            }
            if actual.has_overlaps {
                return Err(EphemerisError::CorruptKernel {
                    path: self.minor_bodies_kernel.filename.clone(),
                    reason: format!(
                        "Conflicting segment overlap detected for ({}) {}",
                        actual.mpc, actual.name
                    ),
                });
            }
            let raw_start = actual.raw_interval_jd_tdb[0];
            let raw_end = actual.raw_interval_jd_tdb[1];
            if raw_start > required_raw_start || raw_end < required_raw_end {
                return Err(EphemerisError::OutOfCoverage {
                    body: format!("({}) {}", actual.mpc, actual.name),
                    epoch_jd_tdb: JD_1850_01_01,
                    start_jd_tdb: raw_start,
                    end_jd_tdb: raw_end,
                });
            }
        }

        Ok(())
    }

    /// Validates both the manifest metadata (`validate()`) and that the kernel files on disk in
    /// `data_dir` exist and match `size_bytes` recorded in `manifest.json`.
    pub fn validate_with_dir(&self, data_dir: & std::path::Path) -> Result<(), EphemerisError> {
        self.validate()?;

        let p_path = data_dir.join(&self.planetary_kernel.filename);
        let p_meta = std::fs::metadata(&p_path).map_err(|e| EphemerisError::MissingData {
            detail: format!("Missing planetary kernel {}: {e}", p_path.display()),
        })?;
        if p_meta.len() != self.planetary_kernel.size_bytes {
            return Err(EphemerisError::CorruptKernel {
                path: p_path.display().to_string(),
                reason: format!(
                    "File size mismatch against manifest.json: expected {} bytes, found {} bytes",
                    self.planetary_kernel.size_bytes,
                    p_meta.len()
                ),
            });
        }

        let m_path = data_dir.join(&self.minor_bodies_kernel.filename);
        let m_meta = std::fs::metadata(&m_path).map_err(|e| EphemerisError::MissingData {
            detail: format!("Missing minor bodies kernel {}: {e}", m_path.display()),
        })?;
        if m_meta.len() != self.minor_bodies_kernel.size_bytes {
            return Err(EphemerisError::CorruptKernel {
                path: m_path.display().to_string(),
                reason: format!(
                    "File size mismatch against manifest.json: expected {} bytes, found {} bytes",
                    self.minor_bodies_kernel.size_bytes,
                    m_meta.len()
                ),
            });
        }

        Ok(())
    }
}

/// Metadata for atomic IERS time table synchronization (`lxephemeris-cli data-sync`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TimeSyncManifest {
    pub schema_version: String,
    pub synced_at_utc: String,
    pub iers_release_date: String,
    pub finals2000a_sha256: String,
    pub finals2000a_bytes: usize,
    pub leapseconds_sha256: String,
    pub leapseconds_bytes: usize,
    pub coverage_start_mjd: f64,
    pub coverage_end_mjd: f64,
}
