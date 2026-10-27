//! Offline commercial license verification (`Ed25519` / RFC 8032) for `lxephemeris`.
//!
//! Public packages contain **only** signature verification and the embedded 32-byte public key.
//! License key generation and signing are performed offline by the maintainer.
//!
//! # Token Format
//! ```text
//! LXEPH-1.<base64url_nopad(payload_bytes)>.<base64url_nopad(ed25519_sig_64_bytes)>
//! ```
//! The Ed25519 signature covers `b"lxephemeris-license-v1:" || payload_bytes`.
//! Signature verification is performed on the exact raw payload bytes **before** JSON parsing.

#![allow(clippy::indexing_slicing)]

use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{OnceLock, RwLock};

use base64::engine::general_purpose::URL_SAFE_NO_PAD;
use base64::Engine as _;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha512};

/// Release date (`YYYY-MM-DD`) of this `lxephemeris` build, compared against `updates_until`
/// in commercial license keys so that a running deployment never expires on its own.
pub const ENGINE_RELEASE_DATE: &str = "2026-10-27";

/// Maximum accepted license token length in bytes.
const MAX_TOKEN_LEN: usize = 4096;

/// Domain separator prepended to raw payload bytes prior to Ed25519 verification.
const DOMAIN_SEPARATOR: &[u8] = b"lxephemeris-license-v1:";

/// Embedded maintainer Ed25519 public key (`kid = "k1-2026"`).
pub const PUBLIC_KEY_K1_2026: [u8; 32] = [
    0xa5, 0x37, 0xef, 0xda, 0x2c, 0x9d, 0xe0, 0xc0, 0xb5, 0xa5, 0x59, 0x1b, 0x21, 0xc6, 0x7f, 0xce,
    0xe1, 0x85, 0xb3, 0xc5, 0x22, 0xac, 0x8f, 0x98, 0xd9, 0x68, 0x40, 0xe3, 0x7c, 0x2e, 0x68, 0x06,
];

static NOTICE_EMITTED: AtomicBool = AtomicBool::new(false);
static REGISTERED_LICENSE_KEY: OnceLock<RwLock<Option<String>>> = OnceLock::new();

fn registered_key_lock() -> &'static RwLock<Option<String>> {
    REGISTERED_LICENSE_KEY.get_or_init(|| RwLock::new(None))
}

/// Status of the configured `lxephemeris` license credential.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LicenseStatus {
    /// No commercial license key configured; running under `PolyForm-Noncommercial-1.0.0`.
    Missing,
    /// Valid Ed25519-signed commercial license covering this engine release.
    Valid,
    /// Valid Ed25519 signature, but `updates_until` precedes [`ENGINE_RELEASE_DATE`].
    Expired,
    /// Ed25519 signature verification failed or key ID is unknown.
    InvalidSignature,
    /// Token structure, base64 encoding, or JSON payload is malformed.
    Malformed,
    /// Token format version is unsupported by this engine version.
    UnsupportedVersion,
}

impl LicenseStatus {
    /// Returns the snake_case string representation of the license status.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Missing => "missing",
            Self::Valid => "valid",
            Self::Expired => "expired",
            Self::InvalidSignature => "invalid_signature",
            Self::Malformed => "malformed",
            Self::UnsupportedVersion => "unsupported_version",
        }
    }
}

/// Verified license metadata returned by [`verify_license_key`] and [`crate::Engine::license_info`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct LicenseInfo {
    /// Verification status.
    pub status: LicenseStatus,
    /// Active license tier (`"noncommercial"`, `"commercial"`, or `"enterprise"`).
    pub tier: String,
    /// Whether a valid, non-expired commercial license covering this release is active.
    pub is_commercial_valid: bool,
    /// Unique license identifier (e.g. `"LIC-2026-001"`), if present.
    pub license_id: Option<String>,
    /// Licensed organization or holder name, if present.
    pub holder: Option<String>,
    /// Signing key identifier (e.g. `"k1-2026"`), if present.
    pub key_id: Option<String>,
    /// Issuance date (`YYYY-MM-DD`), if present.
    pub issued_at: Option<String>,
    /// Latest engine release date (`YYYY-MM-DD`) covered by annual updates (`None` = perpetual).
    pub updates_until: Option<String>,
    /// Release date of this `lxephemeris` build (`YYYY-MM-DD`).
    pub engine_release_date: String,
    /// Human-readable diagnostic summary.
    pub message: String,
}

impl LicenseInfo {
    #[must_use]
    fn missing() -> Self {
        Self {
            status: LicenseStatus::Missing,
            tier: "noncommercial".to_owned(),
            is_commercial_valid: false,
            license_id: None,
            holder: None,
            key_id: None,
            issued_at: None,
            updates_until: None,
            engine_release_date: ENGINE_RELEASE_DATE.to_owned(),
            message: "No commercial license credential configured (running under PolyForm-Noncommercial-1.0.0).".to_owned(),
        }
    }

    #[must_use]
    fn error(status: LicenseStatus, msg: impl Into<String>) -> Self {
        Self {
            status,
            tier: "noncommercial".to_owned(),
            is_commercial_valid: false,
            license_id: None,
            holder: None,
            key_id: None,
            issued_at: None,
            updates_until: None,
            engine_release_date: ENGINE_RELEASE_DATE.to_owned(),
            message: msg.into(),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawLicensePayload {
    v: u32,
    kid: String,
    lid: String,
    holder: String,
    tier: String,
    issued_at: String,
    #[serde(default)]
    updates_until: Option<String>,
}

/// Normalizes a pasted token by stripping outer quotes (`"..."` / `'...'`) and removing
/// any accidental whitespace or line breaks inserted by email/chat clients.
fn normalize_token_str(raw: &str) -> String {
    let mut s = raw.trim();
    if (s.starts_with('"') && s.ends_with('"') && s.len() >= 2)
        || (s.starts_with('\'') && s.ends_with('\'') && s.len() >= 2)
    {
        s = &s[1..s.len() - 1];
    }
    s.chars().filter(|c| !c.is_ascii_whitespace()).collect()
}

/// Validates a calendar date string in strict `YYYY-MM-DD` format.
#[must_use]
pub fn is_valid_iso_date(s: &str) -> bool {
    if s.len() != 10 {
        return false;
    }
    let b = s.as_bytes();
    if b[4] != b'-' || b[7] != b'-' {
        return false;
    }
    let Ok(year) = s[0..4].parse::<u32>() else {
        return false;
    };
    let Ok(month) = s[5..7].parse::<u32>() else {
        return false;
    };
    let Ok(day) = s[8..10].parse::<u32>() else {
        return false;
    };
    if !(2000..=2999).contains(&year) || !(1..=12).contains(&month) || day < 1 {
        return false;
    }
    let is_leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    let max_day = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 => {
            if is_leap {
                29
            } else {
                28
            }
        }
        _ => return false,
    };
    day <= max_day
}

/// Registers a commercial license key globally for the current process and returns its [`LicenseInfo`].
///
/// Passing an empty string clears any previously registered key. Passing an invalid key returns
/// the diagnostic [`LicenseInfo`] **without** overwriting an already valid registered key.
pub fn register_commercial_license(key: &str) -> LicenseInfo {
    let normalized = normalize_token_str(key);
    if normalized.is_empty() {
        if let Ok(mut guard) = registered_key_lock().write() {
            *guard = None;
        }
        return LicenseInfo::missing();
    }
    let info = verify_license_key(&normalized);
    if info.is_commercial_valid {
        if let Ok(mut guard) = registered_key_lock().write() {
            *guard = Some(normalized);
        }
    }
    info
}

/// Clears any process-global license key set via [`register_commercial_license`].
pub fn clear_commercial_license() {
    if let Ok(mut guard) = registered_key_lock().write() {
        *guard = None;
    }
}

/// Resolves the active license credential from (1) explicit config key, (2) [`register_commercial_license`],
/// (3) `LXEPHEMERIS_LICENSE_KEY`, or (4) `LXEPHEMERIS_LICENSE_FILE`.
#[must_use]
pub fn resolve_license_info(explicit_key: Option<&str>) -> LicenseInfo {
    if let Some(k) = explicit_key {
        let norm = normalize_token_str(k);
        if !norm.is_empty() {
            return verify_license_key(&norm);
        }
    }
    if let Ok(guard) = registered_key_lock().read() {
        if let Some(ref k) = *guard {
            return verify_license_key(k);
        }
    }
    if let Ok(env_key) = std::env::var("LXEPHEMERIS_LICENSE_KEY") {
        let norm = normalize_token_str(&env_key);
        if !norm.is_empty() {
            return verify_license_key(&norm);
        }
    }
    if let Ok(file_path) = std::env::var("LXEPHEMERIS_LICENSE_FILE") {
        let trimmed_path = file_path.trim();
        if !trimmed_path.is_empty() {
            match std::fs::read_to_string(trimmed_path) {
                Ok(contents) => {
                    let norm = normalize_token_str(&contents);
                    if !norm.is_empty() {
                        return verify_license_key(&norm);
                    }
                }
                Err(e) => {
                    return LicenseInfo::error(
                        LicenseStatus::Malformed,
                        format!("Failed to read LXEPHEMERIS_LICENSE_FILE ({trimmed_path}): {e}"),
                    );
                }
            }
        }
    }
    LicenseInfo::missing()
}

/// Emits a one-time, non-fatal `stderr` reminder when running in a production environment
/// (`LXEPHEMERIS_ENV=production`, or `NODE_ENV=production` when `LXEPHEMERIS_ENV` is unset)
/// without a valid commercial license.
pub fn maybe_emit_production_license_notice(explicit_key: Option<&str>) {
    if NOTICE_EMITTED.load(Ordering::Relaxed) {
        return;
    }
    if std::env::var("LXEPHEMERIS_SUPPRESS_LICENSE_NOTICE")
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        .unwrap_or(false)
    {
        return;
    }
    let is_prod = match std::env::var("LXEPHEMERIS_ENV") {
        Ok(val) => val.trim().eq_ignore_ascii_case("production"),
        Err(_) => std::env::var("NODE_ENV")
            .map(|v| v.trim().eq_ignore_ascii_case("production"))
            .unwrap_or(false),
    };
    if !is_prod {
        return;
    }
    let info = resolve_license_info(explicit_key);
    if info.is_commercial_valid {
        return;
    }
    if !NOTICE_EMITTED.swap(true, Ordering::SeqCst) {
        let _ = writeln!(
            std::io::stderr(),
            "[lxephemeris] {} Set LXEPHEMERIS_LICENSE_KEY or LXEPHEMERIS_LICENSE_FILE for commercial deployments.",
            info.message
        );
    }
}

/// Verifies an offline `LXEPH-1.<payload>.<signature>` commercial license key against the embedded
/// `Ed25519` public key (`PUBLIC_KEY_K1_2026`).
#[must_use]
pub fn verify_license_key(token: &str) -> LicenseInfo {
    verify_license_key_with_pubkey(token, &PUBLIC_KEY_K1_2026)
}

/// Internal helper that verifies a license token against a provided 32-byte Ed25519 public key.
#[must_use]
fn verify_license_key_with_pubkey(token: &str, pub_key: &[u8; 32]) -> LicenseInfo {
    if token.len() > MAX_TOKEN_LEN {
        return LicenseInfo::error(
            LicenseStatus::Malformed,
            "License token exceeds maximum allowed length (4096 bytes).",
        );
    }
    let normalized = normalize_token_str(token);
    if normalized.is_empty() {
        return LicenseInfo::missing();
    }
    let parts: Vec<&str> = normalized.split('.').collect();
    if parts.len() != 3 {
        return LicenseInfo::error(
            LicenseStatus::Malformed,
            "License token must have format LXEPH-1.<payload>.<signature>.",
        );
    }
    let prefix = parts[0];
    if prefix != "LXEPH-1" {
        if prefix.starts_with("LXEPH-") {
            return LicenseInfo::error(
                LicenseStatus::UnsupportedVersion,
                format!("Unsupported license token prefix '{prefix}' (expected 'LXEPH-1')."),
            );
        }
        return LicenseInfo::error(
            LicenseStatus::Malformed,
            "License token must begin with 'LXEPH-1.'.",
        );
    }
    let Ok(payload_bytes) = URL_SAFE_NO_PAD.decode(parts[1]) else {
        return LicenseInfo::error(
            LicenseStatus::Malformed,
            "License payload is not valid base64url (no-pad).",
        );
    };
    let Ok(sig_bytes) = URL_SAFE_NO_PAD.decode(parts[2]) else {
        return LicenseInfo::error(
            LicenseStatus::Malformed,
            "License signature is not valid base64url (no-pad).",
        );
    };
    let Ok(sig_arr): Result<[u8; 64], _> = sig_bytes.try_into() else {
        return LicenseInfo::error(
            LicenseStatus::InvalidSignature,
            "License Ed25519 signature must be exactly 64 bytes.",
        );
    };

    // Build domain-separated signed message BEFORE parsing JSON
    let mut signed_msg = Vec::with_capacity(DOMAIN_SEPARATOR.len() + payload_bytes.len());
    signed_msg.extend_from_slice(DOMAIN_SEPARATOR);
    signed_msg.extend_from_slice(&payload_bytes);

    if !ed25519_verify(pub_key, &signed_msg, &sig_arr) {
        return LicenseInfo::error(
            LicenseStatus::InvalidSignature,
            "Ed25519 signature verification failed for license token.",
        );
    }

    let Ok(payload) = serde_json::from_slice::<RawLicensePayload>(&payload_bytes) else {
        return LicenseInfo::error(
            LicenseStatus::Malformed,
            "License JSON payload failed schema validation.",
        );
    };
    if payload.v != 1 {
        return LicenseInfo::error(
            LicenseStatus::UnsupportedVersion,
            format!("Unsupported license schema version v={} (expected 1).", payload.v),
        );
    }
    if payload.kid != "k1-2026" {
        return LicenseInfo::error(
            LicenseStatus::InvalidSignature,
            format!("Unknown signing key ID '{}'.", payload.kid),
        );
    }
    if !matches!(payload.tier.as_str(), "commercial" | "enterprise") {
        return LicenseInfo::error(
            LicenseStatus::Malformed,
            format!(
                "Invalid commercial license tier '{}' (expected 'commercial' or 'enterprise').",
                payload.tier
            ),
        );
    }
    if payload.lid.trim().is_empty()
        || payload.holder.trim().is_empty()
        || !is_valid_iso_date(&payload.issued_at)
    {
        return LicenseInfo::error(
            LicenseStatus::Malformed,
            "License payload contains empty holder/lid or invalid issued_at calendar date.",
        );
    }
    if let Some(ref u) = payload.updates_until {
        if !is_valid_iso_date(u) || u.as_str() < payload.issued_at.as_str() {
            return LicenseInfo::error(
                LicenseStatus::Malformed,
                "License updates_until field is not a valid calendar date >= issued_at.",
            );
        }
    }

    // Check release coverage (updates_until vs ENGINE_RELEASE_DATE)
    if let Some(ref u) = payload.updates_until {
        if ENGINE_RELEASE_DATE > u.as_str() {
            let reason = format!(
                "Commercial license '{}' for '{}' covers engine releases up to {}, but this engine release is {}. Please renew your annual license.",
                payload.lid, payload.holder, u, ENGINE_RELEASE_DATE
            );
            return LicenseInfo {
                status: LicenseStatus::Expired,
                tier: payload.tier,
                is_commercial_valid: false,
                license_id: Some(payload.lid),
                holder: Some(payload.holder),
                key_id: Some(payload.kid),
                issued_at: Some(payload.issued_at),
                updates_until: payload.updates_until,
                engine_release_date: ENGINE_RELEASE_DATE.to_owned(),
                message: reason,
            };
        }
    }

    let coverage = payload
        .updates_until
        .as_deref()
        .map_or_else(|| "perpetual".to_owned(), |d| format!("updates until {d}"));
    let msg = format!(
        "Valid commercial license '{}' for '{}' ({coverage}).",
        payload.lid, payload.holder
    );
    LicenseInfo {
        status: LicenseStatus::Valid,
        tier: payload.tier,
        is_commercial_valid: true,
        license_id: Some(payload.lid),
        holder: Some(payload.holder),
        key_id: Some(payload.kid),
        issued_at: Some(payload.issued_at),
        updates_until: payload.updates_until,
        engine_release_date: ENGINE_RELEASE_DATE.to_owned(),
        message: msg,
    }
}

// ============================================================================
// Pure-Rust RFC 8032 Ed25519 Signature Verification (#![forbid(unsafe_code)])
// ============================================================================

type Fe = [i64; 16];

const FE_ZERO: Fe = [0; 16];
const FE_ONE: Fe = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

// d = -121665 / 121666 mod p
const FE_D: Fe = [
    0x78a3, 0x1359, 0x4dca, 0x75eb, 0xd8ab, 0x4141, 0x0a4d, 0x0070, 0xe898, 0x7779, 0x4079, 0x8cc7,
    0xfe73, 0x2b6f, 0x6cee, 0x5203,
];

// Base point B = (X, Y) where Y = 4/5 mod p
const FE_BX: Fe = [
    0xd51a, 0x8f25, 0x2d60, 0xc956, 0xa7b2, 0x9525, 0xc760, 0x692c, 0xdc5c, 0xfdd6, 0xe231, 0xc0a4,
    0x53fe, 0xcd6e, 0x36d3, 0x2169,
];
const FE_BY: Fe = [
    0x6658, 0x6666, 0x6666, 0x6666, 0x6666, 0x6666, 0x6666, 0x6666, 0x6666, 0x6666, 0x6666, 0x6666,
    0x6666, 0x6666, 0x6666, 0x6666,
];

// Group order L = 2^252 + 27742317777372353535851937790883648493 in little-endian bytes
const ED25519_L: [u8; 32] = [
    0xed, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58, 0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10,
];

fn fe_car(o: &mut Fe) {
    for _ in 0..2 {
        for i in 0..16 {
            o[i] += 1 << 16;
            let c = o[i] >> 16;
            if i < 15 {
                o[i + 1] += c - 1;
            } else {
                o[0] += 38 * (c - 1);
            }
            o[i] -= c << 16;
        }
    }
}

fn fe_add(a: &Fe, b: &Fe) -> Fe {
    let mut o = FE_ZERO;
    for i in 0..16 {
        o[i] = a[i] + b[i];
    }
    o
}

fn fe_sub(a: &Fe, b: &Fe) -> Fe {
    let mut o = FE_ZERO;
    for i in 0..16 {
        o[i] = a[i] - b[i];
    }
    o
}

fn fe_mul(a: &Fe, b: &Fe) -> Fe {
    let mut t = [0i64; 31];
    for i in 0..16 {
        for j in 0..16 {
            t[i + j] += a[i] * b[j];
        }
    }
    for i in 0..15 {
        t[i] += 38 * t[i + 16];
    }
    let mut o = FE_ZERO;
    o.copy_from_slice(&t[..16]);
    fe_car(&mut o);
    fe_car(&mut o);
    o
}

fn fe_sq(a: &Fe) -> Fe {
    fe_mul(a, a)
}

fn fe_inv(i: &Fe) -> Fe {
    let mut c = *i;
    for a in (0..=253).rev() {
        c = fe_sq(&c);
        if a != 2 && a != 4 {
            c = fe_mul(&c, i);
        }
    }
    c
}

fn fe_pow2523(i: &Fe) -> Fe {
    let mut c = *i;
    for a in (0..=250).rev() {
        c = fe_sq(&c);
        if a != 1 {
            c = fe_mul(&c, i);
        }
    }
    c
}

fn fe_from_bytes(s: &[u8; 32]) -> Fe {
    let mut o = FE_ZERO;
    for i in 0..16 {
        o[i] = i64::from(s[2 * i]) + (i64::from(s[2 * i + 1]) << 8);
    }
    o[15] &= 0x7fff;
    o
}

fn fe_to_bytes(n: &Fe) -> [u8; 32] {
    let mut t = *n;
    fe_car(&mut t);
    fe_car(&mut t);
    fe_car(&mut t);
    let mut m = FE_ZERO;
    for _ in 0..2 {
        m[0] = t[0] - 0xffed;
        for i in 1..15 {
            m[i] = t[i] - 0xffff - ((m[i - 1] >> 16) & 1);
            m[i - 1] &= 0xffff;
        }
        m[15] = t[15] - 0x7fff - ((m[14] >> 16) & 1);
        let b = (m[15] >> 16) & 1;
        m[14] &= 0xffff;
        for i in 0..16 {
            t[i] = if b == 0 { m[i] } else { t[i] };
        }
    }
    let mut out = [0u8; 32];
    for i in 0..16 {
        out[2 * i] = (t[i] & 0xff) as u8;
        out[2 * i + 1] = ((t[i] >> 8) & 0xff) as u8;
    }
    out
}

fn fe_par(a: &Fe) -> u8 {
    let d = fe_to_bytes(a);
    d[0] & 1
}

type GePoint = [Fe; 4];

fn ge_cswap(p: &mut GePoint, q: &mut GePoint, b: i64) {
    let mask = !(b - 1);
    for k in 0..4 {
        for i in 0..16 {
            let t = mask & (p[k][i] ^ q[k][i]);
            p[k][i] ^= t;
            q[k][i] ^= t;
        }
    }
}

fn ge_add(p: &mut GePoint, q: &GePoint) {
    let a = fe_add(&p[1], &p[0]);
    let t = fe_add(&q[1], &q[0]);
    let a = fe_mul(&a, &t);
    let b = fe_sub(&p[1], &p[0]);
    let t = fe_sub(&q[1], &q[0]);
    let b = fe_mul(&b, &t);
    let c = fe_mul(&p[3], &q[3]);
    let c = fe_mul(&c, &FE_D);
    let c = fe_add(&c, &c);
    let d = fe_mul(&p[2], &q[2]);
    let d = fe_add(&d, &d);
    let e = fe_sub(&a, &b);
    let f = fe_sub(&d, &c);
    let g = fe_add(&d, &c);
    let h = fe_add(&a, &b);
    p[0] = fe_mul(&e, &f);
    p[1] = fe_mul(&h, &g);
    p[2] = fe_mul(&g, &f);
    p[3] = fe_mul(&e, &h);
}

fn ge_scalarmult(q: &GePoint, s: &[u8; 32]) -> GePoint {
    let mut p = [FE_ZERO, FE_ONE, FE_ONE, FE_ZERO];
    let mut r = *q;
    for i in (0..=255).rev() {
        let b = i64::from((s[i >> 3] >> (i & 7)) & 1);
        ge_cswap(&mut p, &mut r, b);
        ge_add(&mut r, &p);
        let mut p2 = p;
        ge_add(&mut p, &p2);
        p2.fill(FE_ZERO);
        ge_cswap(&mut p, &mut r, b);
    }
    p
}

fn ge_pack(p: &GePoint) -> [u8; 32] {
    let zi = fe_inv(&p[2]);
    let tx = fe_mul(&p[0], &zi);
    let ty = fe_mul(&p[1], &zi);
    let mut r = fe_to_bytes(&ty);
    r[31] ^= fe_par(&tx) << 7;
    r
}

fn ge_unpackneg(r: &[u8; 32]) -> Option<GePoint> {
    let mut p = [FE_ZERO, FE_ZERO, FE_ONE, FE_ZERO];
    p[1] = fe_from_bytes(r);
    let c = fe_sq(&p[1]);
    let d = fe_mul(&c, &FE_D);
    let a = fe_sub(&c, &FE_ONE);
    let b = fe_add(&d, &FE_ONE);

    let t = fe_sq(&b);
    let t = fe_mul(&t, &b);
    let mut x = fe_sq(&t);
    x = fe_mul(&x, &b);
    x = fe_mul(&x, &a);
    x = fe_pow2523(&x);
    x = fe_mul(&x, &t);
    x = fe_mul(&x, &a);

    let chk = fe_mul(&b, &fe_sq(&x));
    if fe_to_bytes(&chk) != fe_to_bytes(&a) {
        const I_SQRT_M1: Fe = [
            0xa0b0, 0x4a0e, 0x1b27, 0xc4ee, 0xe478, 0xad2f, 0x1806, 0x2f43, 0xd7a7, 0x3dfb, 0x0099,
            0x2b4d, 0xdf0b, 0x4fc1, 0x2480, 0x2b83,
        ];
        x = fe_mul(&x, &I_SQRT_M1);
        let chk2 = fe_mul(&b, &fe_sq(&x));
        if fe_to_bytes(&chk2) != fe_to_bytes(&a) {
            return None;
        }
    }
    if fe_par(&x) != (r[31] >> 7) {
        x = fe_sub(&FE_ZERO, &x);
    }
    p[0] = x;
    p[3] = fe_mul(&p[0], &p[1]);
    p[0] = fe_sub(&FE_ZERO, &p[0]);
    p[3] = fe_sub(&FE_ZERO, &p[3]);
    Some(p)
}

fn scalar_is_canonical(s: &[u8; 32]) -> bool {
    for i in (0..32).rev() {
        if s[i] < ED25519_L[i] {
            return true;
        }
        if s[i] > ED25519_L[i] {
            return false;
        }
    }
    false
}

fn reduce_mod_l(x: &[u8; 64]) -> [u8; 32] {
    let mut v = [0i64; 64];
    for i in 0..64 {
        v[i] = i64::from(x[i]);
    }
    for i in (32..64).rev() {
        let mut carry = 0i64;
        for j in (i - 32)..(i - 12) {
            v[j] -= v[i] * i64::from(ED25519_L[j - (i - 32)]) * 16;
            carry = v[j] >> 8;
            v[j + 1] += carry;
            v[j] &= 255;
        }
        v[i - 12] += carry;
        v[i] = 0;
    }
    let mut carry = 0i64;
    for j in 0..32 {
        v[j] += carry - (v[31] >> 4) * i64::from(ED25519_L[j]);
        carry = v[j] >> 8;
        v[j] &= 255;
    }
    for j in 0..32 {
        v[j] -= carry * i64::from(ED25519_L[j]);
    }
    let mut out = [0u8; 32];
    for i in 0..32 {
        v[i + 1] += v[i] >> 8;
        out[i] = (v[i] & 255) as u8;
    }
    out
}

/// Verifies an RFC 8032 `Ed25519` signature `sig` (`64` bytes) over `msg` against `pub_key` (`32` bytes).
#[must_use]
pub fn ed25519_verify(pub_key: &[u8; 32], msg: &[u8], sig: &[u8; 64]) -> bool {
    let mut r_bytes = [0u8; 32];
    let mut s_bytes = [0u8; 32];
    r_bytes.copy_from_slice(&sig[..32]);
    s_bytes.copy_from_slice(&sig[32..]);

    if !scalar_is_canonical(&s_bytes) {
        return false;
    }
    let Some(neg_a) = ge_unpackneg(pub_key) else {
        return false;
    };

    let mut hasher = Sha512::new();
    hasher.update(r_bytes);
    hasher.update(pub_key);
    hasher.update(msg);
    let h_digest = hasher.finalize();
    let mut h_raw = [0u8; 64];
    h_raw.copy_from_slice(&h_digest);
    let k = reduce_mod_l(&h_raw);

    let b_point = [FE_BX, FE_BY, FE_ONE, fe_mul(&FE_BX, &FE_BY)];
    let mut sb = ge_scalarmult(&b_point, &s_bytes);
    let ka = ge_scalarmult(&neg_a, &k);
    ge_add(&mut sb, &ka);
    let r_calc = ge_pack(&sb);

    let mut diff = 0u8;
    for i in 0..32 {
        diff |= r_calc[i] ^ r_bytes[i];
    }
    diff == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Isolated test-only Ed25519 public key (private key discarded at test generation time;
    /// NEVER matches `PUBLIC_KEY_K1_2026`).
    const TEST_PUBKEY: [u8; 32] = [
        0xe5, 0x7f, 0x7e, 0x1f, 0x74, 0x5f, 0xb3, 0x50, 0x08, 0xba, 0x40, 0x46, 0x5a, 0x61, 0xfa, 0xbe,
        0x3c, 0x19, 0x70, 0x17, 0xdd, 0x6c, 0xf0, 0x2c, 0xc3, 0xed, 0x5b, 0x43, 0x26, 0xcd, 0x4f, 0xb1,
    ];

    #[test]
    fn test_rfc8032_vector1_and_tamper() {
        let pub_key: [u8; 32] = [
            0xd7, 0x5a, 0x98, 0x01, 0x82, 0xb1, 0x0a, 0xb7, 0xd5, 0x4b, 0xfe, 0xd3, 0xc9, 0x64,
            0x07, 0x3a, 0x0e, 0xe1, 0x72, 0xf3, 0xda, 0xa6, 0x23, 0x25, 0xaf, 0x02, 0x1a, 0x68,
            0xf7, 0x07, 0x51, 0x1a,
        ];
        let sig: [u8; 64] = [
            0xe5, 0x56, 0x43, 0x00, 0xc3, 0x60, 0xac, 0x72, 0x90, 0x86, 0xe2, 0xcc, 0x80, 0x6e,
            0x82, 0x8a, 0x84, 0x87, 0x7f, 0x1e, 0xb8, 0xe5, 0xd9, 0x74, 0xd8, 0x73, 0xe0, 0x65,
            0x22, 0x49, 0x01, 0x55, 0x5f, 0xb8, 0x82, 0x15, 0x90, 0xa3, 0x3b, 0xac, 0xc6, 0x1e,
            0x39, 0x70, 0x1c, 0xf9, 0xb4, 0x6b, 0xd2, 0x5b, 0xf5, 0xf0, 0x59, 0x5b, 0xbe, 0x24,
            0x65, 0x51, 0x41, 0x43, 0x8e, 0x7a, 0x10, 0x0b,
        ];
        assert!(ed25519_verify(&pub_key, b"", &sig));
        assert!(!ed25519_verify(&pub_key, b"x", &sig));
    }

    #[test]
    fn test_isolated_fixture_token_normalization_and_expiry() {
        let valid_test_token = "LXEPH-1.eyJ2IjoxLCJraWQiOiJrMS0yMDI2IiwibGlkIjoiTElDLVRFU1QtMDAxIiwiaG9sZGVyIjoiVW5pdCBUZXN0IEZpeHR1cmUiLCJ0aWVyIjoiY29tbWVyY2lhbCIsImlzc3VlZF9hdCI6IjIwMjYtMTAtMjciLCJ1cGRhdGVzX3VudGlsIjpudWxsfQ.8mcKwP75ero9NKpO9YtvGR-d5aKRtIVe8zwJ8a69gAbgc00oQy5y1fIlHoEjDTN4QU-iGDlyj4cIp3NBWde4Cg";

        // 1. Test fixture token MUST NOT validate against production PUBLIC_KEY_K1_2026
        assert_eq!(
            verify_license_key(valid_test_token).status,
            LicenseStatus::InvalidSignature
        );

        // 2. Validates against isolated TEST_PUBKEY
        let info = verify_license_key_with_pubkey(valid_test_token, &TEST_PUBKEY);
        assert_eq!(info.status, LicenseStatus::Valid);
        assert!(info.is_commercial_valid);
        assert_eq!(info.license_id.as_deref(), Some("LIC-TEST-001"));
        assert_eq!(info.holder.as_deref(), Some("Unit Test Fixture"));

        // 3. Surrounding quotes and embedded line-wraps are normalized automatically
        let quoted_wrapped = format!(
            "  \"{}\n{}\"  \n",
            &valid_test_token[..60],
            &valid_test_token[60..]
        );
        let norm_info = verify_license_key_with_pubkey(&quoted_wrapped, &TEST_PUBKEY);
        assert_eq!(norm_info.status, LicenseStatus::Valid);

        // 4. Expired token (updates_until = 2025-10-27 < ENGINE_RELEASE_DATE)
        let expired_test_token = "LXEPH-1.eyJ2IjoxLCJraWQiOiJrMS0yMDI2IiwibGlkIjoiTElDLVRFU1QtRVhQIiwiaG9sZGVyIjoiRXhwaXJlZCBUZXN0IEZpeHR1cmUiLCJ0aWVyIjoiY29tbWVyY2lhbCIsImlzc3VlZF9hdCI6IjIwMjQtMTAtMjciLCJ1cGRhdGVzX3VudGlsIjoiMjAyNS0xMC0yNyJ9.HeIsLgwsZWFkI2XpU9Lsq-quHQOD_3p43kDyjdKRvimBx9CDwLRxGcQpTfH4OtiZmFqvJRUj0gm2PfVS40jmAA";
        let exp_info = verify_license_key_with_pubkey(expired_test_token, &TEST_PUBKEY);
        assert_eq!(exp_info.status, LicenseStatus::Expired);
        assert!(!exp_info.is_commercial_valid);

        // 5. Strict calendar date check
        assert!(is_valid_iso_date("2026-10-27"));
        assert!(is_valid_iso_date("2028-02-29"));
        assert!(!is_valid_iso_date("2027-02-29"));
        assert!(!is_valid_iso_date("2027-99-99"));
        assert!(!is_valid_iso_date("2027-13-01"));
    }
}
