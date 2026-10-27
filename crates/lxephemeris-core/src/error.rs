//! Structured, non-exhaustive error types and stable error codes for `lxephemeris`.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Stable machine-readable error code returned across Rust, Python, and Node.js boundaries.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ErrorCode {
    InvalidInput,
    UnsupportedBody,
    UnsupportedModel,
    UnsupportedSpkType,
    OutOfCoverage,
    OutOfModelValidity,
    MissingData,
    CorruptKernel,
    UndefinedHouseSystem,
    NonConvergence,
}

impl ErrorCode {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::InvalidInput => "InvalidInput",
            Self::UnsupportedBody => "UnsupportedBody",
            Self::UnsupportedModel => "UnsupportedModel",
            Self::UnsupportedSpkType => "UnsupportedSpkType",
            Self::OutOfCoverage => "OutOfCoverage",
            Self::OutOfModelValidity => "OutOfModelValidity",
            Self::MissingData => "MissingData",
            Self::CorruptKernel => "CorruptKernel",
            Self::UndefinedHouseSystem => "UndefinedHouseSystem",
            Self::NonConvergence => "NonConvergence",
        }
    }
}

/// Primary error type for all `lxephemeris` operations.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Error, Serialize, Deserialize)]
pub enum EphemerisError {
    #[error("Invalid input: {message}")]
    InvalidInput { message: String },

    #[error("Unsupported body: {body}")]
    UnsupportedBody { body: String },

    #[error("Unsupported model: {model}")]
    UnsupportedModel { model: String },

    #[error("Unsupported SPK segment type {spk_type} for target NAIF ID {target_id}")]
    UnsupportedSpkType { target_id: i32, spk_type: i32 },

    #[error(
        "Epoch JD(TDB) {epoch_jd_tdb:.6} is out of coverage [{start_jd_tdb:.6}, {end_jd_tdb:.6}) for body {body}"
    )]
    OutOfCoverage {
        body: String,
        epoch_jd_tdb: f64,
        start_jd_tdb: f64,
        end_jd_tdb: f64,
    },

    #[error("Epoch JD(TT) {epoch_jd_tt:.6} is outside validity window [{min_jd_tt:.6}, {max_jd_tt:.6}] for model {model}")]
    OutOfModelValidity {
        model: String,
        epoch_jd_tt: f64,
        min_jd_tt: f64,
        max_jd_tt: f64,
    },

    #[error("Missing required ephemeris data: {detail}")]
    MissingData { detail: String },

    #[error("Corrupt or invalid DAF/SPK kernel ({path}): {reason}")]
    CorruptKernel { path: String, reason: String },

    #[error(
        "House system '{house_system}' is undefined at latitude {latitude_deg:.4}°: {reason}"
    )]
    UndefinedHouseSystem {
        house_system: char,
        latitude_deg: f64,
        reason: String,
    },

    #[error("Iterative solver '{solver}' failed to converge after {iterations} iterations (residual: {residual:.3e})")]
    NonConvergence {
        solver: String,
        iterations: u32,
        residual: f64,
    },
}

impl EphemerisError {
    /// Returns the stable machine-readable [`ErrorCode`] for this error variant.
    #[must_use]
    pub const fn code(&self) -> ErrorCode {
        match self {
            Self::InvalidInput { .. } => ErrorCode::InvalidInput,
            Self::UnsupportedBody { .. } => ErrorCode::UnsupportedBody,
            Self::UnsupportedModel { .. } => ErrorCode::UnsupportedModel,
            Self::UnsupportedSpkType { .. } => ErrorCode::UnsupportedSpkType,
            Self::OutOfCoverage { .. } => ErrorCode::OutOfCoverage,
            Self::OutOfModelValidity { .. } => ErrorCode::OutOfModelValidity,
            Self::MissingData { .. } => ErrorCode::MissingData,
            Self::CorruptKernel { .. } => ErrorCode::CorruptKernel,
            Self::UndefinedHouseSystem { .. } => ErrorCode::UndefinedHouseSystem,
            Self::NonConvergence { .. } => ErrorCode::NonConvergence,
        }
    }
}
