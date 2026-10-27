//! `Engine`: immutable, thread-safe (`Send + Sync`) ephemeris engine with offline kernel loading
//! (`from_data_dir`, `from_bsp_paths`, `from_bytes`), operational `safe_interval` enforcement,
//! capabilities/coverage introspection, and fork-safe lazy `rayon::ThreadPool` batch evaluation.

use crate::almanac::{
    calc_lots, calc_panchanga, calc_shodashavarga, calc_varga, calc_vimshottari_dasha, chart_sect,
    map_zodiac_placement, planetary_hours, ChartSectInfo, HermeticLot, PanchangaReport,
    PlanetaryHoursReport, VargaPlacement, VimshottariDashaReport, ZodiacPlacement,
};
use crate::aspects::{
    calc_chart_aspects, calc_declination_aspects, calc_synastry_aspects, mutual_aspects_in_range,
    natal_transits_in_range, search_exact_aspect, void_of_course_moon_in_range, AspectHit,
    AspectPreset, DeclinationAspectHit, ExactAspectEvent, VoidOfCourseInterval,
};
use crate::astrometry::{
    apply_aberration_and_deflection_icrf, azalt, azalt_from_ecliptic, azalt_rev,
    calc_body_in_moment, eval_spk_target_vector_au_seeded, resolve_spk_target, rise_trans,
    time_equation_days, AtmosphereParams, CalcResult, HorizontalCoords, RiseTransEvent,
    RiseTransOptions, C_KM_S,
};
use crate::ayanamsa::{eval_ayanamsa, AyanamsaValue};
use crate::eclipses::{
    lunar_eclipse_at, lunar_eclipse_how_at_ut1, lunar_eclipses_in_range, next_local_solar_eclipse,
    next_lunar_eclipse, next_solar_eclipse, solar_eclipse_at, solar_eclipse_where_at_ut1,
    solar_eclipses_in_range, EclipseOptions, LunarEclipseInfo, LunarEclipseKind, SolarEclipseInfo,
    SolarEclipseKind,
};
use crate::error::EphemerisError;
use crate::events::{
    body_longitude_crossings_in_range, find_parans, human_design_epoch, lunar_phase_summary,
    lunar_phases_in_range, planetary_return, planetary_returns_in_range, search_body_longitude,
    search_declination_extrema_in_range, search_lunar_phase, search_stations_in_range,
    seasons_in_year, solar_terms_in_range, DeclinationExtremumEvent, LongitudeCrossingEvent,
    LunarPhaseEvent, LunarPhaseSummary, ParanEvent, ParanOptions, PlanetaryReturnEvent,
    SeasonsOfYear, SolarTermEvent, StationEvent,
};
use crate::frame::{
    cart_state_to_spherical, cart_to_spherical, mat_transpose, mat_vec, spherical_state_to_cart,
    vec_add, vec_dot, vec_norm, vec_scale, vec_sub, EclNutInfo, FrameTransform, SphericalState,
    Vec3,
};
use crate::heliacal::{heliacal_phenomena, HeliacalOptions, HeliacalResult};
use crate::houses::{
    calc_houses, calc_houses_in_moment, gauquelin_sectors, house_pos, houses_armc, GauquelinResult,
    HouseOptions, HousesResult,
};
use crate::lunar::{JD_INTERP_LILITH_SAFE_MAX, JD_INTERP_LILITH_SAFE_MIN};
use crate::manifest::{sha256_hex, DataManifest};
use crate::moment::Moment;
use crate::moshier::{MoshierAnalyticBackend, JD_MOSHIER_MAX, JD_MOSHIER_MIN};
use crate::nodes_apsides::{calc_nodes_apsides, NodesApsidesResult};
use crate::occultations::{
    local_lunar_occultations_in_range, lunar_occultations_in_range, next_local_lunar_occultation,
    next_lunar_occultation, LunarOccultationInfo,
};
use crate::phenomena::{calc_phenomena_in_moment, compute_apparent_magnitude_v, PhenomenaResult};
use crate::spk::{AniseSpkBackend, SpkSource};
use crate::stars::{
    calc_star_entry_in_moment, find_star_by_hip, find_star_by_name_or_bayer, list_stars,
    load_custom_stars_file, star_catalog_52, star_info, StarEntry, StarMetadata,
};
use crate::time::{DeltaTInfo, SiderealTimeInfo, TimeInput, TimeTable};
use crate::types::{
    AbMode, AyanamsaModel, BodyId, CalcOptions, HouseSystem, LunarPointId, NutationMode, Observer,
    Origin, Plane, PolarPolicy, Zodiac, AU_KM, JD_1850_01_01, JD_2100_01_01, MINOR_BODIES_68,
    SECONDS_PER_DAY,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// Configuration options for initializing an [`Engine`].
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct EngineConfig {
    /// Maximum number of worker threads in the engine's private, lazily initialized Rayon pool
    /// (`None` defaults to `std::thread::available_parallelism().min(8)`).
    pub max_threads: Option<usize>,
    /// Optional pinned time table version (`"embedded"` or a specific synced IERS release date)
    /// for 100% deterministic test reproducibility.
    pub pin_time_table_version: Option<String>,
    /// Whether to verify SHA-256 hashes of `.bsp` files against `manifest.json` on load.
    pub verify_manifest_sha256: bool,
}

/// Operational and raw coverage interval for a specific [`BodyId`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BodyCoverage {
    pub body: BodyId,
    pub target_id: i32,
    pub raw_start_jd_tdb: f64,
    pub raw_end_jd_tdb: f64,
    pub safe_start_jd_tt: f64,
    pub safe_end_jd_tt: f64,
}

/// Summary of loaded capabilities and data versions on an [`Engine`] (`R10.4`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EngineCapabilities {
    pub ephemeris_source: String,
    pub ephemeris_version: String,
    pub loaded_segment_count: usize,
    pub minor_bodies_loaded_count: usize,
    pub time_table_version: String,
    pub house_systems: Vec<char>,
    pub ayanamsa_models: Vec<String>,
}

/// Enriched body placement inside a [`ChartSnapshot`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChartBodySnapshot {
    pub body: BodyId,
    pub body_name: String,
    pub calc: CalcResult,
    pub placement: ZodiacPlacement,
    pub house_number: u8,
    pub house_pos_3d: f64,
    pub is_retrograde: bool,
}

/// Complete, one-call astrological and astronomical chart snapshot (`Wave H`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChartSnapshot {
    pub jd_ut1: f64,
    pub jd_tt: f64,
    pub observer: Observer,
    pub bodies: Vec<ChartBodySnapshot>,
    pub houses: HousesResult,
    pub aspects: Vec<AspectHit>,
    pub lunar_phase: LunarPhaseSummary,
    pub sect: ChartSectInfo,
    pub lots: Vec<HermeticLot>,
}

/// Immutable, thread-safe (`Send + Sync`) astronomical and astrological ephemeris engine.
pub struct Engine {
    spk: Option<AniseSpkBackend>,
    moshier: MoshierAnalyticBackend,
    custom_stars: Vec<StarEntry>,
    time_table: TimeTable,
    config: EngineConfig,
    manifest: Option<DataManifest>,
    loaded_kernel_names: Vec<String>,
    thread_pool: OnceLock<Result<rayon::ThreadPool, String>>,
}

static_assertions::assert_impl_all!(Engine: Send, Sync);

impl Default for Engine {
    /// Creates a zero-friction [`Engine`]: if the `LXEPHEMERIS_DATA_DIR` environment variable is set
    /// and contains valid `.bsp` kernels, loads SPK mode; otherwise initializes the fileless
    /// analytical Moshier engine (`Engine::new_moshier`).
    fn default() -> Self {
        Self::from_env_or_moshier(EngineConfig::default())
    }
}

impl Engine {
    /// Creates an [`Engine`] by checking `LXEPHEMERIS_DATA_DIR` first and falling back
    /// to the fileless analytical Moshier engine (`Engine::new_moshier`) only when `.bsp`
    /// files are absent (`MissingData`). Returns `Err` if `LXEPHEMERIS_DATA_DIR` points to
    /// a corrupt manifest or kernel.
    pub fn try_from_env_or_moshier(config: EngineConfig) -> Result<Self, EphemerisError> {
        if let Ok(env_dir) = std::env::var("LXEPHEMERIS_DATA_DIR") {
            let trimmed = env_dir.trim();
            if !trimmed.is_empty() {
                return Self::from_data_dir_or_moshier(trimmed, config);
            }
        }
        Ok(Self::new_moshier(config))
    }

    /// Creates an [`Engine`] by checking `LXEPHEMERIS_DATA_DIR` first and seamlessly
    /// falling back to the fileless analytical Moshier engine (`Engine::new_moshier`) if no `.bsp`
    /// files are present. Emits a diagnostic warning to `stderr` if `LXEPHEMERIS_DATA_DIR` contains
    /// a corrupt manifest or kernel.
    #[must_use]
    pub fn from_env_or_moshier(config: EngineConfig) -> Self {
        match Self::try_from_env_or_moshier(config.clone()) {
            Ok(eng) => eng,
            Err(err) => {
                use std::io::Write as _;
                let _ = writeln!(
                    std::io::stderr(),
                    "[lxephemeris] Warning: Failed to load SPK bundle from LXEPHEMERIS_DATA_DIR ({err}); falling back to analytical Moshier mode."
                );
                Self::new_moshier(config)
            }
        }
    }

    /// Creates a fileless analytical ("Moshier" / `SEFLG_MOSEPH`) [`Engine`] covering
    /// `[-3000 CE, +3000 CE]` (`JD [625307.5, 2816787.5]`) with **zero external `.bsp` files required**.
    #[must_use]
    pub fn new_moshier(config: EngineConfig) -> Self {
        crate::license::maybe_emit_production_license_notice(None);
        Self {
            spk: None,
            moshier: MoshierAnalyticBackend::new(),
            custom_stars: Vec::new(),
            time_table: TimeTable::embedded(),
            config,
            manifest: None,
            loaded_kernel_names: Vec::new(),
            thread_pool: OnceLock::new(),
        }
    }

    /// Loads an [`Engine`] from `data_dir` if `.bsp` kernels exist, or seamlessly falls back
    /// to the fileless analytical Moshier engine (`Engine::new_moshier`) when `.bsp` files are absent.
    pub fn from_data_dir_or_moshier(
        data_dir: impl AsRef<Path>,
        config: EngineConfig,
    ) -> Result<Self, EphemerisError> {
        let dir = data_dir.as_ref();
        match Self::from_data_dir(dir, config.clone()) {
            Ok(eng) => Ok(eng),
            Err(EphemerisError::MissingData { .. }) => {
                let mut eng = Self::new_moshier(config);
                if dir.exists() {
                    let sef_path = dir.join("sefstars.txt");
                    if sef_path.exists() {
                        let _ = eng.load_custom_stars(&sef_path);
                    }
                }
                Ok(eng)
            }
            Err(other) => Err(other),
        }
    }

    /// Loads an [`Engine`] from a directory containing `manifest.json` (or `.bsp` files) and optional
    /// synced IERS time tables. Never accesses the network.
    pub fn from_data_dir(
        data_dir: impl AsRef<Path>,
        config: EngineConfig,
    ) -> Result<Self, EphemerisError> {
        let dir = data_dir.as_ref();
        if !dir.exists() {
            return Err(EphemerisError::MissingData {
                detail: format!("Data directory does not exist: {}", dir.display()),
            });
        }

        let manifest_path = dir.join("manifest.json");
        let (bsp_paths, manifest_opt) = if manifest_path.exists() {
            let text = fs::read_to_string(&manifest_path).map_err(|e| EphemerisError::MissingData {
                detail: format!("Failed to read {}: {e}", manifest_path.display()),
            })?;
            let manifest: DataManifest =
                serde_json::from_str(&text).map_err(|e| EphemerisError::CorruptKernel {
                    path: manifest_path.display().to_string(),
                    reason: format!("Invalid manifest.json: {e}"),
                })?;
            manifest.validate_with_dir(dir)?;

            let p_path = dir.join(&manifest.planetary_kernel.filename);
            let m_path = dir.join(&manifest.minor_bodies_kernel.filename);

            if config.verify_manifest_sha256 {
                let p_bytes = fs::read(&p_path).map_err(|e| EphemerisError::MissingData {
                    detail: format!("Failed to read {}: {e}", p_path.display()),
                })?;
                if sha256_hex(&p_bytes) != manifest.planetary_kernel.sha256 {
                    return Err(EphemerisError::CorruptKernel {
                        path: p_path.display().to_string(),
                        reason: "SHA-256 mismatch against manifest.json".to_owned(),
                    });
                }
                let m_bytes = fs::read(&m_path).map_err(|e| EphemerisError::MissingData {
                    detail: format!("Failed to read {}: {e}", m_path.display()),
                })?;
                if sha256_hex(&m_bytes) != manifest.minor_bodies_kernel.sha256 {
                    return Err(EphemerisError::CorruptKernel {
                        path: m_path.display().to_string(),
                        reason: "SHA-256 mismatch against manifest.json".to_owned(),
                    });
                }
            }

            (vec![p_path, m_path], Some(manifest))
        } else {
            let mut paths: Vec<PathBuf> = Vec::new();
            for candidate in [
                "inpop21a.bsp",
                "inpop21a_TDB_m100_p100_spice.bsp",
                "de441.bsp",
                "de441_part-1.bsp",
                "de441_part-2.bsp",
                "de440.bsp",
                "de440s.bsp",
                "asteroids_68_1850_2100.bsp",
            ] {
                let p = dir.join(candidate);
                if p.exists() {
                    paths.push(p);
                }
            }
            if paths.is_empty() {
                return Err(EphemerisError::MissingData {
                    detail: format!("No .bsp kernels found in {}", dir.display()),
                });
            }
            (paths, None)
        };

        let loaded_kernel_names = bsp_paths
            .iter()
            .map(|p| {
                p.file_name()
                    .map_or_else(|| p.display().to_string(), |n| n.to_string_lossy().to_string())
            })
            .collect();

        let spk = AniseSpkBackend::from_bsp_paths(&bsp_paths)?;
        let time_table =
            TimeTable::load_or_embedded(dir, config.pin_time_table_version.as_deref())?;
        let custom_stars = if dir.join("sefstars.txt").exists() {
            load_custom_stars_file(dir.join("sefstars.txt")).unwrap_or_default()
        } else {
            Vec::new()
        };

        crate::license::maybe_emit_production_license_notice(None);
        Ok(Self {
            spk: Some(spk),
            moshier: MoshierAnalyticBackend::new(),
            custom_stars,
            time_table,
            config,
            manifest: manifest_opt,
            loaded_kernel_names,
            thread_pool: OnceLock::new(),
        })
    }

    /// Loads an [`Engine`] from explicit `.bsp` file paths using the embedded time table
    /// (or a pinned version).
    pub fn from_bsp_paths<P: AsRef<Path>>(
        paths: &[P],
        config: EngineConfig,
    ) -> Result<Self, EphemerisError> {
        let loaded_kernel_names = paths
            .iter()
            .map(|p| {
                p.as_ref().file_name().map_or_else(
                    || p.as_ref().display().to_string(),
                    |n| n.to_string_lossy().to_string(),
                )
            })
            .collect();
        let spk = AniseSpkBackend::from_bsp_paths(paths)?;
        let time_table = if let Some(first) = paths.first() {
            if let Some(parent) = first.as_ref().parent() {
                TimeTable::load_or_embedded(parent, config.pin_time_table_version.as_deref())?
            } else {
                TimeTable::embedded()
            }
        } else {
            TimeTable::embedded()
        };
        crate::license::maybe_emit_production_license_notice(None);
        Ok(Self {
            spk: Some(spk),
            moshier: MoshierAnalyticBackend::new(),
            custom_stars: Vec::new(),
            time_table,
            config,
            manifest: None,
            loaded_kernel_names,
            thread_pool: OnceLock::new(),
        })
    }

    /// Loads an [`Engine`] directly from in-memory `(name, bytes)` `.bsp` buffers.
    pub fn from_bytes(
        kernels: &[(String, Vec<u8>)],
        config: EngineConfig,
    ) -> Result<Self, EphemerisError> {
        let loaded_kernel_names = kernels.iter().map(|(n, _)| n.clone()).collect();
        let spk = AniseSpkBackend::from_named_bytes(kernels)?;
        let time_table = TimeTable::embedded();
        crate::license::maybe_emit_production_license_notice(None);
        Ok(Self {
            spk: Some(spk),
            moshier: MoshierAnalyticBackend::new(),
            custom_stars: Vec::new(),
            time_table,
            config,
            manifest: None,
            loaded_kernel_names,
            thread_pool: OnceLock::new(),
        })
    }

    /// Returns the active [`crate::license::LicenseInfo`] for this process (`PolyForm-Noncommercial-1.0.0`
    /// or verified `Ed25519` commercial license).
    #[must_use]
    pub fn license_info(&self) -> crate::license::LicenseInfo {
        crate::license::resolve_license_info(None)
    }

    /// Loads an [`Engine`] from `LXEPHEMERIS_DATA_DIR` (if set) or `./data`.
    pub fn from_env(config: EngineConfig) -> Result<Self, EphemerisError> {
        let dir = std::env::var("LXEPHEMERIS_DATA_DIR")
            .map_or_else(|_| PathBuf::from("data"), PathBuf::from);
        Self::from_data_dir(dir, config)
    }

    /// Loads additional custom stars from a `sefstars.txt` or JSON file into this [`Engine`].
    pub fn load_custom_stars(&mut self, path: impl AsRef<Path>) -> Result<usize, EphemerisError> {
        let loaded = load_custom_stars_file(path)?;
        let count = loaded.len();
        self.custom_stars.extend(loaded);
        Ok(count)
    }

    /// Appends parsed [`StarEntry`] records into this [`Engine`]'s custom star catalog.
    pub fn add_custom_stars(&mut self, stars: impl IntoIterator<Item = StarEntry>) -> usize {
        let before = self.custom_stars.len();
        self.custom_stars.extend(stars);
        self.custom_stars.len() - before
    }

    /// Builder helper to attach a custom star catalog (`sefstars.txt` or JSON) to this [`Engine`].
    pub fn with_custom_stars_file(
        mut self,
        path: impl AsRef<Path>,
    ) -> Result<Self, EphemerisError> {
        let _ = self.load_custom_stars(path)?;
        Ok(self)
    }

    /// Returns `true` if this engine is operating in pure fileless analytical ("Moshier") mode
    /// without external `.bsp` kernels.
    #[must_use]
    pub const fn is_moshier_mode(&self) -> bool {
        self.spk.is_none()
    }

    /// Returns the loaded [`DataManifest`] if `from_data_dir` found `manifest.json`.
    #[must_use]
    pub const fn manifest(&self) -> Option<&DataManifest> {
        self.manifest.as_ref()
    }

    /// Returns the active [`TimeTable`].
    #[must_use]
    pub const fn time_table(&self) -> &TimeTable {
        &self.time_table
    }

    /// Returns the active [`SpkSource`] (`AniseSpkBackend` if loaded, or `MoshierAnalyticBackend` in fileless mode).
    #[must_use]
    pub fn spk(&self) -> &dyn SpkSource {
        self.spk
            .as_ref()
            .map_or(&self.moshier as &dyn SpkSource, |s| s as &dyn SpkSource)
    }

    /// Returns the fileless analytical [`MoshierAnalyticBackend`] (`SEFLG_MOSEPH`).
    #[must_use]
    pub const fn moshier_spk(&self) -> &MoshierAnalyticBackend {
        &self.moshier
    }

    /// Returns a summary of loaded capabilities, body counts, and active time table version (`R10.4`).
    #[must_use]
    pub fn capabilities(&self) -> EngineCapabilities {
        let active_spk = self.spk();
        let minor_bodies_loaded_count = MINOR_BODIES_68
            .iter()
            .filter(|spec| active_spk.target_coverage(spec.mpc.standard_naif_spk_id()).is_ok())
            .count();
        let ephemeris_source = if self.spk.is_none() {
            "Moshier-DE440-Recurrent".to_owned()
        } else {
            let names_lower: Vec<String> = if let Some(m) = &self.manifest {
                vec![
                    m.planetary_kernel.filename.to_ascii_lowercase(),
                    m.minor_bodies_kernel.filename.to_ascii_lowercase(),
                ]
            } else {
                self.loaded_kernel_names
                    .iter()
                    .map(|n| n.to_ascii_lowercase())
                    .collect()
            };
            if names_lower.iter().any(|n| n.contains("inpop21")) {
                "INPOP21a".to_owned()
            } else if names_lower.iter().any(|n| n.contains("de441")) {
                "DE441".to_owned()
            } else if names_lower.iter().any(|n| n.contains("de440s")) {
                "DE440s".to_owned()
            } else if names_lower.iter().any(|n| n.contains("de440")) {
                "DE440".to_owned()
            } else {
                self.manifest.as_ref().map_or_else(
                    || "DE440s".to_owned(),
                    |m| {
                        format!(
                            "{} + {}",
                            m.planetary_kernel.filename, m.minor_bodies_kernel.filename
                        )
                    },
                )
            }
        };
        let ephemeris_version = self.manifest.as_ref().map_or_else(
            || format!("lxephemeris-{}", env!("CARGO_PKG_VERSION")),
            |m| {
                format!(
                    "lxephemeris-{}/schema-{}",
                    env!("CARGO_PKG_VERSION"),
                    m.schema_version
                )
            },
        );
        EngineCapabilities {
            ephemeris_source,
            ephemeris_version,
            loaded_segment_count: active_spk.segments().len(),
            minor_bodies_loaded_count,
            time_table_version: self.time_table.version().to_owned(),
            house_systems: HouseSystem::ALL.iter().map(|h| h.code()).collect(),
            ayanamsa_models: AyanamsaModel::ALL
                .iter()
                .map(|a| a.code().to_owned())
                .collect(),
        }
    }

    /// Returns the raw and operational `safe_interval` coverage for `body`.
    pub fn coverage(&self, body: BodyId) -> Result<BodyCoverage, EphemerisError> {
        let active_spk = self.spk();
        match body {
            BodyId::Planet(_) => {
                let (target_id, _, h_days) = resolve_spk_target(body)?;
                let cov = active_spk.target_coverage(target_id)?;
                let margin_days = (2.0 * h_days + 0.5).max(1.0);
                let safe_start = cov.start_jd_tdb + margin_days;
                let safe_end = cov.end_jd_tdb - margin_days;
                Ok(BodyCoverage {
                    body,
                    target_id,
                    raw_start_jd_tdb: cov.start_jd_tdb,
                    raw_end_jd_tdb: cov.end_jd_tdb,
                    safe_start_jd_tt: safe_start,
                    safe_end_jd_tt: safe_end,
                })
            }
            BodyId::MinorPlanet(_) => {
                let (target_id, _, _) = resolve_spk_target(body)?;
                let cov = active_spk.target_coverage(target_id)?;
                let (safe_start, safe_end) = if self.spk.is_none() {
                    (cov.start_jd_tdb + 4.0, cov.end_jd_tdb - 4.0)
                } else {
                    // Stitched minor bodies bundle has a 4-day pad around [1850-01-01, 2100-01-01)
                    (
                        (cov.start_jd_tdb + 4.0).max(JD_1850_01_01),
                        (cov.end_jd_tdb - 4.0).min(JD_2100_01_01),
                    )
                };
                Ok(BodyCoverage {
                    body,
                    target_id,
                    raw_start_jd_tdb: cov.start_jd_tdb,
                    raw_end_jd_tdb: cov.end_jd_tdb,
                    safe_start_jd_tt: safe_start,
                    safe_end_jd_tt: safe_end,
                })
            }
            BodyId::LunarPoint(lp) => {
                let cov = active_spk.target_coverage(301)?;
                let (safe_start, safe_end) = if matches!(
                    lp,
                    LunarPointId::InterpolatedLilith | LunarPointId::InterpolatedPerigee
                ) {
                    if self.spk.is_none() {
                        (cov.start_jd_tdb + 60.0, cov.end_jd_tdb - 60.0)
                    } else {
                        (
                            (cov.start_jd_tdb + 59.0).max(JD_INTERP_LILITH_SAFE_MIN),
                            (cov.end_jd_tdb - 60.0).min(JD_INTERP_LILITH_SAFE_MAX),
                        )
                    }
                } else {
                    (cov.start_jd_tdb + 1.0, cov.end_jd_tdb - 1.0)
                };
                Ok(BodyCoverage {
                    body,
                    target_id: 301,
                    raw_start_jd_tdb: cov.start_jd_tdb,
                    raw_end_jd_tdb: cov.end_jd_tdb,
                    safe_start_jd_tt: safe_start,
                    safe_end_jd_tt: safe_end,
                })
            }
            BodyId::Hypothetical(id) => {
                let cov = active_spk.target_coverage(399)?;
                Ok(BodyCoverage {
                    body,
                    target_id: -(id as i32 + 1),
                    raw_start_jd_tdb: cov.start_jd_tdb,
                    raw_end_jd_tdb: cov.end_jd_tdb,
                    safe_start_jd_tt: cov.start_jd_tdb,
                    safe_end_jd_tt: cov.end_jd_tdb - 1.0,
                })
            }
            BodyId::FixedStar(hip) => {
                if find_star_by_hip(hip).is_none()
                    && !self.custom_stars.iter().any(|s| s.hip == hip.0)
                {
                    return Err(EphemerisError::UnsupportedBody {
                        body: format!("FixedStar(HIP {})", hip.0),
                    });
                }
                let cov = active_spk.target_coverage(399)?;
                Ok(BodyCoverage {
                    body,
                    target_id: hip.0 as i32,
                    raw_start_jd_tdb: cov.start_jd_tdb,
                    raw_end_jd_tdb: cov.end_jd_tdb,
                    safe_start_jd_tt: cov.start_jd_tdb,
                    safe_end_jd_tt: cov.end_jd_tdb - 1.0,
                })
            }
            BodyId::EclNut => {
                let cov = active_spk.target_coverage(399)?;
                Ok(BodyCoverage {
                    body,
                    target_id: -1,
                    raw_start_jd_tdb: cov.start_jd_tdb,
                    raw_end_jd_tdb: cov.end_jd_tdb,
                    safe_start_jd_tt: cov.start_jd_tdb,
                    safe_end_jd_tt: cov.end_jd_tdb,
                })
            }
        }
    }

    /// Evaluates $\Delta T = \text{TT} - \text{UT1}$ and $\text{TT} - \text{TDB}$ at `epoch` (`R1.2`).
    #[must_use]
    pub fn delta_t(&self, epoch: TimeInput) -> DeltaTInfo {
        let (_, _, _, dt_info) = self.time_table.resolve_scales(epoch, None);
        dt_info
    }

    /// Evaluates Greenwich and Local Sidereal Time (`GMST`, `GAST`, `LAST`, `ARMC`) at `epoch` (`R1.4`).
    #[must_use]
    pub fn sidtime(
        &self,
        epoch: TimeInput,
        nutation: NutationMode,
        longitude_deg: f64,
    ) -> SiderealTimeInfo {
        self.time_table.sidtime(epoch, nutation, longitude_deg)
    }

    /// Evaluates ecliptic obliquity ($\varepsilon, \varepsilon_0$) and IAU 2000A nutation angles ($\Delta\psi, \Delta\varepsilon$) at `epoch` (`R2.5`).
    pub fn ecl_nut(&self, epoch: TimeInput) -> Result<EclNutInfo, EphemerisError> {
        let (_, tt, _, _) = self.time_table.resolve_scales(epoch, None);
        let ft = FrameTransform::at_tt(tt);
        Ok(ft.ecl_nut_info())
    }

    /// Returns [`StarMetadata`] for a fixed star looked up by canonical name, Bayer code, or `"hip:<id>"` (`R6.1`).
    pub fn star_info(&self, name_or_bayer: &str) -> Result<StarMetadata, EphemerisError> {
        if let Ok(meta) = star_info(name_or_bayer) {
            return Ok(meta);
        }
        let trimmed = name_or_bayer.trim().trim_start_matches(',');
        let prefix = trimmed.split(',').next().unwrap_or(trimmed).trim();
        if let Some(custom) = self.custom_stars.iter().find(|s| {
            s.name.eq_ignore_ascii_case(prefix)
                || s.bayer.eq_ignore_ascii_case(prefix)
                || s.name.eq_ignore_ascii_case(trimmed)
                || s.bayer.eq_ignore_ascii_case(trimmed)
        }) {
            return Ok(custom.to_metadata());
        }
        Err(EphemerisError::UnsupportedBody {
            body: format!("FixedStar({name_or_bayer})"),
        })
    }

    /// Returns [`StarMetadata`] for all stars in the embedded catalog (plus any loaded custom stars; `R6.1`).
    #[must_use]
    pub fn list_stars(&self) -> Vec<StarMetadata> {
        let mut out = list_stars();
        out.extend(self.custom_stars.iter().map(StarEntry::to_metadata));
        out
    }

    /// Resolves a [`StarEntry`] by case-insensitive canonical `name`, `bayer` code, or `"hip:<id>"`
    /// from the embedded catalog or loaded custom stars.
    pub fn resolve_star_entry(&self, name_or_bayer: &str) -> Result<StarEntry, EphemerisError> {
        if let Some(star) = find_star_by_name_or_bayer(name_or_bayer) {
            return Ok(star.clone());
        }
        let trimmed = name_or_bayer.trim().trim_start_matches(',');
        let prefix = trimmed.split(',').next().unwrap_or(trimmed).trim();
        if let Some(hip_s) = prefix.to_ascii_lowercase().strip_prefix("hip:") {
            if let Ok(hip_num) = hip_s.trim().parse::<u32>() {
                if let Some(custom) = self.custom_stars.iter().find(|s| s.hip == hip_num) {
                    return Ok(custom.clone());
                }
            }
        }
        if let Some(custom) = self.custom_stars.iter().find(|s| {
            s.name.eq_ignore_ascii_case(prefix)
                || s.bayer.eq_ignore_ascii_case(prefix)
                || s.name.eq_ignore_ascii_case(trimmed)
                || s.bayer.eq_ignore_ascii_case(trimmed)
        }) {
            return Ok(custom.clone());
        }
        Err(EphemerisError::UnsupportedBody {
            body: format!("FixedStar({name_or_bayer})"),
        })
    }

    /// Returns all [`StarEntry`] records in the embedded catalog plus any loaded custom stars.
    #[must_use]
    pub fn all_star_entries(&self) -> Vec<StarEntry> {
        let mut out = star_catalog_52().to_vec();
        out.extend(self.custom_stars.iter().cloned());
        out
    }

    /// Creates a per-epoch evaluation context [`Moment`] at `input` for an optional topocentric `observer`.
    pub fn at(
        &self,
        input: TimeInput,
        observer: Option<Observer>,
    ) -> Result<Moment<'_>, EphemerisError> {
        Moment::new(self.spk(), &self.time_table, input, observer)
    }

    /// Creates a per-epoch evaluation context [`Moment`] backed explicitly by the fileless
    /// analytical [`MoshierAnalyticBackend`] (`SEFLG_MOSEPH`).
    pub fn at_moshier(
        &self,
        input: TimeInput,
        observer: Option<Observer>,
    ) -> Result<Moment<'_>, EphemerisError> {
        Moment::new(&self.moshier, &self.time_table, input, observer)
    }

    /// Evaluates `body` inside an already constructed [`Moment`], enforcing `safe_interval`.
    pub fn calc_in_moment(
        &self,
        moment: &Moment<'_>,
        body: BodyId,
        options: CalcOptions,
    ) -> Result<CalcResult, EphemerisError> {
        let cov = self.coverage(body)?;
        let jd_tt = moment.center().tt.jd_f64();
        if jd_tt < cov.safe_start_jd_tt || jd_tt >= cov.safe_end_jd_tt {
            return Err(EphemerisError::OutOfCoverage {
                body: format!("{body:?}"),
                epoch_jd_tdb: jd_tt,
                start_jd_tdb: cov.safe_start_jd_tt,
                end_jd_tdb: cov.safe_end_jd_tt,
            });
        }
        if let BodyId::FixedStar(hip) = body {
            if find_star_by_hip(hip).is_none() {
                if let Some(custom) = self.custom_stars.iter().find(|s| s.hip == hip.0) {
                    let sid = moment.sidereal_offset_and_rate(&options)?;
                    return calc_star_entry_in_moment(moment, custom, &options, sid);
                }
            }
        }
        match body {
            BodyId::Planet(_)
                if options.include_speed && options.plane != Plane::Horizontal =>
            {
                calc_spk_body_analytical(moment, body, &options)
            }
            _ => calc_body_in_moment(moment, body, &options),
        }
    }

    /// Evaluates `body` inside `moment` using single-pass analytical Chebyshev velocity propagation
    /// (`calc_spk_body_fast`), bypassing the 4-point Cartesian stencil while maintaining `< 1e-6 deg/day` accuracy.
    pub fn calc_analytical_velocity_in_moment(
        &self,
        moment: &Moment<'_>,
        body: BodyId,
        options: &CalcOptions,
    ) -> Result<CalcResult, EphemerisError> {
        let cov = self.coverage(body)?;
        let jd_tt = moment.center().tt.jd_f64();
        if jd_tt < cov.safe_start_jd_tt || jd_tt >= cov.safe_end_jd_tt {
            return Err(EphemerisError::OutOfCoverage {
                body: format!("{body:?}"),
                epoch_jd_tdb: jd_tt,
                start_jd_tdb: cov.safe_start_jd_tt,
                end_jd_tdb: cov.safe_end_jd_tt,
            });
        }
        match body {
            BodyId::Planet(_) | BodyId::MinorPlanet(_) => {
                if options.plane == Plane::Horizontal {
                    calc_body_in_moment(moment, body, options)
                } else {
                    calc_spk_body_analytical(moment, body, options)
                }
            }
            _ => calc_body_in_moment(moment, body, options),
        }
    }

    /// Evaluates `body` at `input` with `options` using single-pass analytical Chebyshev velocity propagation.
    pub fn calc_analytical_velocity(
        &self,
        body: BodyId,
        input: TimeInput,
        options: CalcOptions,
    ) -> Result<CalcResult, EphemerisError> {
        let obs = match options.origin {
            Origin::Topocentric(o) => Some(o),
            _ => None,
        };
        let cov = self.coverage(body)?;
        let (_, tt, _, _) = self.time_table.resolve_scales(input, obs);
        let jd_tt = tt.jd_f64();
        if jd_tt < cov.safe_start_jd_tt || jd_tt >= cov.safe_end_jd_tt {
            return Err(EphemerisError::OutOfCoverage {
                body: format!("{body:?}"),
                epoch_jd_tdb: jd_tt,
                start_jd_tdb: cov.safe_start_jd_tt,
                end_jd_tdb: cov.safe_end_jd_tt,
            });
        }
        let moment = self.at(input, obs)?;
        self.calc_analytical_velocity_in_moment(&moment, body, &options)
    }

    /// Evaluates a single `body` at `input` with `options`.
    pub fn calc(
        &self,
        body: BodyId,
        input: TimeInput,
        options: CalcOptions,
    ) -> Result<CalcResult, EphemerisError> {
        let obs = match options.origin {
            Origin::Topocentric(o) => Some(o),
            _ => None,
        };
        let cov = self.coverage(body)?;
        let (_, tt, _, _) = self.time_table.resolve_scales(input, obs);
        let jd_tt = tt.jd_f64();
        if jd_tt < cov.safe_start_jd_tt || jd_tt >= cov.safe_end_jd_tt {
            return Err(EphemerisError::OutOfCoverage {
                body: format!("{body:?}"),
                epoch_jd_tdb: jd_tt,
                start_jd_tdb: cov.safe_start_jd_tt,
                end_jd_tdb: cov.safe_end_jd_tt,
            });
        }
        let moment = self.at(input, obs)?;
        self.calc_in_moment(&moment, body, options)
    }

    /// Evaluates `body` at `input` with `options` explicitly using the fileless analytical
    /// [`MoshierAnalyticBackend`] (`SEFLG_MOSEPH`, `[-3000 CE, +3000 CE]`), even if `.bsp` kernels are loaded.
    pub fn calc_moshier(
        &self,
        body: BodyId,
        input: TimeInput,
        options: CalcOptions,
    ) -> Result<CalcResult, EphemerisError> {
        let obs = match options.origin {
            Origin::Topocentric(o) => Some(o),
            _ => None,
        };
        let (_, tt, _, _) = self.time_table.resolve_scales(input, obs);
        let jd_tt = tt.jd_f64();
        let margin = if matches!(
            body,
            BodyId::LunarPoint(
                LunarPointId::InterpolatedLilith | LunarPointId::InterpolatedPerigee
            )
        ) {
            60.0
        } else {
            4.0
        };
        if jd_tt < JD_MOSHIER_MIN + margin || jd_tt >= JD_MOSHIER_MAX - margin {
            return Err(EphemerisError::OutOfCoverage {
                body: format!("{body:?}"),
                epoch_jd_tdb: jd_tt,
                start_jd_tdb: JD_MOSHIER_MIN + margin,
                end_jd_tdb: JD_MOSHIER_MAX - margin,
            });
        }
        let moment = self.at_moshier(input, obs)?;
        if let BodyId::FixedStar(hip) = body {
            if find_star_by_hip(hip).is_none() {
                if let Some(custom) = self.custom_stars.iter().find(|s| s.hip == hip.0) {
                    let sid = moment.sidereal_offset_and_rate(&options)?;
                    return calc_star_entry_in_moment(&moment, custom, &options, sid);
                }
            }
        }
        calc_body_in_moment(&moment, body, &options)
    }

    /// Evaluates a fixed star by case-insensitive canonical `name`, `bayer` designation, or `"hip:<id>"`.
    pub fn calc_star_by_name(
        &self,
        query: &str,
        input: TimeInput,
        options: CalcOptions,
    ) -> Result<CalcResult, EphemerisError> {
        if let Some(star) = find_star_by_name_or_bayer(query) {
            return self.calc(BodyId::FixedStar(star.hip_id()), input, options);
        }
        let trimmed = query.trim().trim_start_matches(',');
        let prefix = trimmed.split(',').next().unwrap_or(trimmed).trim();
        if let Some(custom) = self.custom_stars.iter().find(|s| {
            s.name.eq_ignore_ascii_case(prefix)
                || s.bayer.eq_ignore_ascii_case(prefix)
                || s.name.eq_ignore_ascii_case(trimmed)
                || s.bayer.eq_ignore_ascii_case(trimmed)
        }) {
            return self.calc(BodyId::FixedStar(custom.hip_id()), input, options);
        }
        Err(EphemerisError::UnsupportedBody {
            body: format!("FixedStar({query})"),
        })
    }

    /// Evaluates the requested [`AyanamsaModel`] at `input` with `nutation` mode.
    pub fn ayanamsa(
        &self,
        input: TimeInput,
        model: AyanamsaModel,
        nutation: NutationMode,
    ) -> Result<AyanamsaValue, EphemerisError> {
        let (_, tt, _, _) = self.time_table.resolve_scales(input, None);
        let jd_tt = tt.jd_f64();
        match eval_ayanamsa(self.spk(), model, jd_tt, nutation) {
            Ok(v) => Ok(v),
            Err(EphemerisError::OutOfCoverage { .. }) => {
                eval_ayanamsa(&self.moshier, model, jd_tt, nutation)
            }
            Err(e) => Err(e),
        }
    }

    ///Lazily initializes and returns the engine's private, bounded `rayon::ThreadPool`.
    fn get_thread_pool(&self) -> Result<&rayon::ThreadPool, EphemerisError> {
        let res = self.thread_pool.get_or_init(|| {
            let num_threads = self.config.max_threads.unwrap_or_else(|| {
                std::thread::available_parallelism()
                    .map(std::num::NonZeroUsize::get)
                    .unwrap_or(4)
                    .min(8)
            });
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads.max(1))
                .thread_name(|idx| format!("lxephemeris-worker-{idx}"))
                .build()
                .map_err(|e| format!("Failed to initialize Rayon thread pool: {e}"))
        });
        match res {
            Ok(pool) => Ok(pool),
            Err(msg) => Err(EphemerisError::InvalidInput {
                message: msg.clone(),
            }),
        }
    }

    /// Evaluates `bodies` across `epochs` in parallel using the engine's private Rayon pool
    /// (atomic mode: fails fast if any epoch/body fails).
    pub fn calc_batch(
        &self,
        bodies: &[BodyId],
        epochs: &[TimeInput],
        options: CalcOptions,
    ) -> Result<Vec<Vec<CalcResult>>, EphemerisError> {
        let obs = match options.origin {
            Origin::Topocentric(o) => Some(o),
            _ => None,
        };
        if epochs.len() <= 2 || cfg!(target_arch = "wasm32") {
            epochs
                .iter()
                .map(|&ep| {
                    let moment = self.at(ep, obs)?;
                    bodies
                        .iter()
                        .map(|&b| self.calc_in_moment(&moment, b, options))
                        .collect::<Result<Vec<_>, _>>()
                })
                .collect()
        } else if let Ok(pool) = self.get_thread_pool() {
            pool.install(|| {
                epochs
                    .par_iter()
                    .map(|&ep| {
                        let moment = self.at(ep, obs)?;
                        bodies
                            .iter()
                            .map(|&b| self.calc_in_moment(&moment, b, options))
                            .collect::<Result<Vec<_>, _>>()
                    })
                    .collect()
            })
        } else {
            epochs
                .iter()
                .map(|&ep| {
                    let moment = self.at(ep, obs)?;
                    bodies
                        .iter()
                        .map(|&b| self.calc_in_moment(&moment, b, options))
                        .collect::<Result<Vec<_>, _>>()
                })
                .collect()
        }
    }

    /// Evaluates `bodies` across `epochs` in parallel using the engine's private Rayon pool,
    /// returning per-item `Result<CalcResult, EphemerisError>` without aborting the batch on single-item errors.
    pub fn calc_batch_per_item(
        &self,
        bodies: &[BodyId],
        epochs: &[TimeInput],
        options: CalcOptions,
    ) -> Result<Vec<Vec<Result<CalcResult, EphemerisError>>>, EphemerisError> {
        let obs = match options.origin {
            Origin::Topocentric(o) => Some(o),
            _ => None,
        };
        let eval_epoch = |ep: TimeInput| -> Vec<Result<CalcResult, EphemerisError>> {
            match self.at(ep, obs) {
                Ok(moment) => bodies
                    .iter()
                    .map(|&b| self.calc_in_moment(&moment, b, options))
                    .collect(),
                Err(err) => bodies.iter().map(|_| Err(err.clone())).collect(),
            }
        };

        if epochs.len() <= 2 || cfg!(target_arch = "wasm32") {
            Ok(epochs.iter().map(|&ep| eval_epoch(ep)).collect())
        } else if let Ok(pool) = self.get_thread_pool() {
            Ok(pool.install(|| epochs.par_iter().map(|&ep| eval_epoch(ep)).collect()))
        } else {
            Ok(epochs.iter().map(|&ep| eval_epoch(ep)).collect())
        }
    }

    /// Computes horizontal coordinates (`Azimuth`, `Altitude`) from True Equatorial of Date `(ra_deg, dec_deg)` at `epoch` for `observer` (`R2.6`).
    #[must_use]
    pub fn azalt_equatorial(
        &self,
        observer: Observer,
        epoch: TimeInput,
        ra_deg: f64,
        dec_deg: f64,
        atmosphere: AtmosphereParams,
    ) -> HorizontalCoords {
        let (ut1, tt, _, _) = self.time_table.resolve_scales(epoch, Some(observer));
        let gast = crate::time::gast_rad(ut1, tt);
        azalt(observer, gast, ra_deg, dec_deg, atmosphere)
    }

    /// Computes horizontal coordinates (`Azimuth`, `Altitude`) from True Ecliptic of Date `(ecl_lon_deg, ecl_lat_deg)` at `epoch` for `observer` (`R2.6 & R7.3`).
    #[must_use]
    pub fn azalt_ecliptic(
        &self,
        observer: Observer,
        epoch: TimeInput,
        ecl_lon_deg: f64,
        ecl_lat_deg: f64,
        atmosphere: AtmosphereParams,
    ) -> HorizontalCoords {
        let (ut1, tt, _, _) = self.time_table.resolve_scales(epoch, Some(observer));
        let ft = FrameTransform::at_tt(tt);
        let gast = crate::time::gast_rad(ut1, tt);
        azalt_from_ecliptic(
            observer,
            gast,
            ft.true_eps_rad,
            ecl_lon_deg,
            ecl_lat_deg,
            atmosphere,
        )
    }

    /// Computes horizontal coordinates (`Azimuth`, `Altitude`) for `body` at `input` for `observer`.
    pub fn calc_azalt(
        &self,
        body: BodyId,
        input: TimeInput,
        observer: Observer,
        atmosphere: AtmosphereParams,
    ) -> Result<HorizontalCoords, EphemerisError> {
        let moment = self.at(input, Some(observer))?;
        let eq_opts = CalcOptions {
            origin: Origin::Topocentric(observer),
            plane: crate::types::Plane::Equatorial,
            equinox: crate::types::Equinox::OfDate,
            nutation: crate::types::NutationMode::True,
            aberration: crate::types::AbMode::Apparent,
            zodiac: crate::types::Zodiac::Tropical,
            include_speed: false,
        };
        let eq_res = self.calc_in_moment(&moment, body, eq_opts)?;
        Ok(azalt(
            observer,
            moment.center().gast_rad,
            eq_res.lon_deg,
            eq_res.lat_deg,
            atmosphere,
        ))
    }

    /// Finds the next rise, set, or meridian transit (`JD_UT1`) of `body` after `jd_ut1_start` at `observer`.
    pub fn rise_trans(
        &self,
        body: BodyId,
        observer: Observer,
        jd_ut1_start: f64,
        event: RiseTransEvent,
        options: RiseTransOptions,
    ) -> Result<Option<f64>, EphemerisError> {
        rise_trans(
            self.spk(),
            &self.time_table,
            body,
            observer,
            jd_ut1_start,
            event,
            options,
        )
    }

    /// Finds the next rise, set, or meridian transit (`JD_UT1`) of `body` after `jd_ut1_start` at `observer`
    /// explicitly using the fileless analytical [`MoshierAnalyticBackend`] (`SEFLG_MOSEPH`).
    pub fn rise_trans_moshier(
        &self,
        body: BodyId,
        observer: Observer,
        jd_ut1_start: f64,
        event: RiseTransEvent,
        options: RiseTransOptions,
    ) -> Result<Option<f64>, EphemerisError> {
        rise_trans(
            &self.moshier,
            &self.time_table,
            body,
            observer,
            jd_ut1_start,
            event,
            options,
        )
    }

    /// Computes annual heliacal rising and setting phenomena (`heliacal_rise_jd_ut1`, `morning_set_jd_ut1`,
    /// `evening_rise_jd_ut1`, `heliacal_set_jd_ut1`) for `body` (`FixedStar`, `Planet`, or `MinorPlanet`)
    /// using Bradley Schaefer's visibility model (`R6.3`).
    pub fn heliacal_phenomena(
        &self,
        body: BodyId,
        start_jd_ut1: f64,
        observer: Observer,
        options: HeliacalOptions,
    ) -> Result<HeliacalResult, EphemerisError> {
        heliacal_phenomena(
            self.spk(),
            &self.time_table,
            body,
            start_jd_ut1,
            observer,
            options,
        )
    }

    /// Computes annual heliacal rising and setting phenomena explicitly using the fileless analytical
    /// [`MoshierAnalyticBackend`] (`SEFLG_MOSEPH`).
    pub fn heliacal_phenomena_moshier(
        &self,
        body: BodyId,
        start_jd_ut1: f64,
        observer: Observer,
        options: HeliacalOptions,
    ) -> Result<HeliacalResult, EphemerisError> {
        heliacal_phenomena(
            &self.moshier,
            &self.time_table,
            body,
            start_jd_ut1,
            observer,
            options,
        )
    }

    /// Evaluates astrological house cusps and angles for `observer` at `input` with `options`.
    pub fn houses(
        &self,
        input: TimeInput,
        observer: Observer,
        options: HouseOptions,
    ) -> Result<HousesResult, EphemerisError> {
        match calc_houses(self.spk(), &self.time_table, input, observer, options) {
            Ok(res) => Ok(res),
            Err(EphemerisError::OutOfCoverage { .. }) => {
                calc_houses(&self.moshier, &self.time_table, input, observer, options)
            }
            Err(e) => Err(e),
        }
    }

    /// Evaluates astrological house cusps and angles inside an existing [`Moment`] for `observer`,
    /// reusing precomputed sidereal time (`GAST`) and obliquity at zero `nut00a` overhead.
    pub fn houses_in_moment(
        &self,
        moment: &Moment<'_>,
        observer: Observer,
        options: HouseOptions,
    ) -> Result<HousesResult, EphemerisError> {
        calc_houses_in_moment(moment, observer, options)
    }

    /// Evaluates astrological house cusps and angles across `epochs` in parallel using the engine's
    /// private Rayon pool (`R9.3`).
    pub fn houses_batch(
        &self,
        epochs: &[TimeInput],
        observer: Observer,
        options: HouseOptions,
    ) -> Result<Vec<HousesResult>, EphemerisError> {
        let spk = self.spk();
        if epochs.len() <= 2 || cfg!(target_arch = "wasm32") {
            epochs
                .iter()
                .map(|&ep| calc_houses(spk, &self.time_table, ep, observer, options))
                .collect()
        } else if let Ok(pool) = self.get_thread_pool() {
            pool.install(|| {
                epochs
                    .par_iter()
                    .map(|&ep| calc_houses(spk, &self.time_table, ep, observer, options))
                    .collect()
            })
        } else {
            epochs
                .iter()
                .map(|&ep| calc_houses(spk, &self.time_table, ep, observer, options))
                .collect()
        }
    }

    /// Evaluates astrological house cusps and angles directly from `armc_deg`, `latitude_deg`,
    /// `obliquity_deg`, and `sun_declination_deg`.
    pub fn houses_armc(
        &self,
        armc_deg: f64,
        latitude_deg: f64,
        obliquity_deg: f64,
        sun_declination_deg: f64,
        options: HouseOptions,
    ) -> Result<HousesResult, EphemerisError> {
        houses_armc(
            armc_deg,
            latitude_deg,
            obliquity_deg,
            sun_declination_deg,
            options,
        )
    }

    /// Evaluates the 36 Gauquelin sectors at `input` for `observer`.
    pub fn gauquelin_sectors(
        &self,
        input: TimeInput,
        observer: Observer,
        polar_policy: PolarPolicy,
    ) -> Result<GauquelinResult, EphemerisError> {
        let m = self.at(input, Some(observer))?;
        let armc_deg =
            (m.center().gast_rad.to_degrees() + observer.longitude_deg).rem_euclid(360.0);
        let obl_deg = m.center().frame.true_eps_rad.to_degrees();
        gauquelin_sectors(armc_deg, observer.latitude_deg, obl_deg, polar_policy)
    }

    /// Evaluates whether a Solar Eclipse occurs near `input` (within $\pm 0.6$ days of syzygy),
    /// returning global Besselian contacts and optional local topocentric circumstances for `observer`.
    pub fn solar_eclipse_at(
        &self,
        input: TimeInput,
        observer: Option<Observer>,
        options: EclipseOptions,
    ) -> Result<Option<SolarEclipseInfo>, EphemerisError> {
        let (_, tt, _, _) = self.time_table.resolve_scales(input, observer);
        solar_eclipse_at(self.spk(), &self.time_table, tt.jd_f64(), observer, options)
    }

    /// Searches forward (`backward == false`) or backward (`backward == true`) from `start_input`
    /// for the next Solar Eclipse.
    pub fn next_solar_eclipse(
        &self,
        start_input: TimeInput,
        backward: bool,
        observer: Option<Observer>,
        options: EclipseOptions,
    ) -> Result<SolarEclipseInfo, EphemerisError> {
        let (_, tt, _, _) = self.time_table.resolve_scales(start_input, observer);
        match next_solar_eclipse(
            self.spk(),
            &self.time_table,
            tt.jd_f64(),
            backward,
            observer,
            options,
        ) {
            Ok(info) => Ok(info),
            Err(EphemerisError::OutOfCoverage { .. }) => next_solar_eclipse(
                &self.moshier,
                &self.time_table,
                tt.jd_f64(),
                backward,
                observer,
                options,
            ),
            Err(e) => Err(e),
        }
    }

    /// Evaluates the instantaneous geographic position and attributes of a Solar Eclipse at `jd_ut1`.
    #[allow(clippy::type_complexity)]
    pub fn solar_eclipse_where_at_ut1(
        &self,
        jd_ut1: f64,
        options: EclipseOptions,
        use_moshier: bool,
    ) -> Result<Option<(SolarEclipseKind, bool, [f64; 10], [f64; 20])>, EphemerisError> {
        if use_moshier || self.is_moshier_mode() {
            return solar_eclipse_where_at_ut1(&self.moshier, &self.time_table, jd_ut1, options);
        }
        match solar_eclipse_where_at_ut1(self.spk(), &self.time_table, jd_ut1, options) {
            Ok(res) => Ok(res),
            Err(EphemerisError::OutOfCoverage { .. }) => {
                solar_eclipse_where_at_ut1(&self.moshier, &self.time_table, jd_ut1, options)
            }
            Err(e) => Err(e),
        }
    }

    /// Evaluates whether a Lunar Eclipse occurs near `input` (within $\pm 0.6$ days of syzygy).
    pub fn lunar_eclipse_at(
        &self,
        input: TimeInput,
        options: EclipseOptions,
    ) -> Result<Option<LunarEclipseInfo>, EphemerisError> {
        let (_, tt, _, _) = self.time_table.resolve_scales(input, None);
        lunar_eclipse_at(self.spk(), &self.time_table, tt.jd_f64(), options)
    }

    /// Searches forward (`backward == false`) or backward (`backward == true`) from `start_input`
    /// for the next Lunar Eclipse.
    pub fn next_lunar_eclipse(
        &self,
        start_input: TimeInput,
        backward: bool,
        options: EclipseOptions,
    ) -> Result<LunarEclipseInfo, EphemerisError> {
        let (_, tt, _, _) = self.time_table.resolve_scales(start_input, None);
        match next_lunar_eclipse(self.spk(), &self.time_table, tt.jd_f64(), backward, options) {
            Ok(info) => Ok(info),
            Err(EphemerisError::OutOfCoverage { .. }) => next_lunar_eclipse(
                &self.moshier,
                &self.time_table,
                tt.jd_f64(),
                backward,
                options,
            ),
            Err(e) => Err(e),
        }
    }

    /// Evaluates the instantaneous attributes of a Lunar Eclipse at `jd_ut1`.
    pub fn lunar_eclipse_how_at_ut1(
        &self,
        jd_ut1: f64,
        options: EclipseOptions,
        use_moshier: bool,
    ) -> Result<Option<(LunarEclipseKind, [f64; 20])>, EphemerisError> {
        if use_moshier || self.is_moshier_mode() {
            return lunar_eclipse_how_at_ut1(&self.moshier, &self.time_table, jd_ut1, options);
        }
        match lunar_eclipse_how_at_ut1(self.spk(), &self.time_table, jd_ut1, options) {
            Ok(res) => Ok(res),
            Err(EphemerisError::OutOfCoverage { .. }) => {
                lunar_eclipse_how_at_ut1(&self.moshier, &self.time_table, jd_ut1, options)
            }
            Err(e) => Err(e),
        }
    }

    /// Searches forward (`backward == false`) or backward (`backward == true`) from `start_input`
    /// for the next Solar Eclipse that is locally visible above the horizon at `observer` (`R8.2`).
    pub fn next_local_solar_eclipse(
        &self,
        start_input: TimeInput,
        observer: Observer,
        backward: bool,
        options: EclipseOptions,
    ) -> Result<SolarEclipseInfo, EphemerisError> {
        let (_, tt, _, _) = self.time_table.resolve_scales(start_input, Some(observer));
        match next_local_solar_eclipse(
            self.spk(),
            &self.time_table,
            tt.jd_f64(),
            observer,
            backward,
            options,
        ) {
            Ok(info) => Ok(info),
            Err(EphemerisError::OutOfCoverage { .. }) => next_local_solar_eclipse(
                &self.moshier,
                &self.time_table,
                tt.jd_f64(),
                observer,
                backward,
                options,
            ),
            Err(e) => Err(e),
        }
    }

    /// Finds all Solar Eclipses occurring within `[start_input, end_input]` (`R8.3`).
    pub fn solar_eclipses_in_range(
        &self,
        start_input: TimeInput,
        end_input: TimeInput,
        observer: Option<Observer>,
        options: EclipseOptions,
    ) -> Result<Vec<SolarEclipseInfo>, EphemerisError> {
        let (_, tt_start, _, _) = self.time_table.resolve_scales(start_input, observer);
        let (_, tt_end, _, _) = self.time_table.resolve_scales(end_input, observer);
        solar_eclipses_in_range(
            self.spk(),
            &self.time_table,
            tt_start.jd_f64(),
            tt_end.jd_f64(),
            observer,
            options,
        )
    }

    /// Finds all Lunar Eclipses occurring within `[start_input, end_input]` (`R8.3`).
    pub fn lunar_eclipses_in_range(
        &self,
        start_input: TimeInput,
        end_input: TimeInput,
        options: EclipseOptions,
    ) -> Result<Vec<LunarEclipseInfo>, EphemerisError> {
        let (_, tt_start, _, _) = self.time_table.resolve_scales(start_input, None);
        let (_, tt_end, _, _) = self.time_table.resolve_scales(end_input, None);
        lunar_eclipses_in_range(
            self.spk(),
            &self.time_table,
            tt_start.jd_f64(),
            tt_end.jd_f64(),
            options,
        )
    }

    /// Evaluates planetary, solar, lunar, or minor-body illumination and photometric phenomena
    /// (`phase_angle_deg`, `phase_illuminated`, `elongation_deg`, `apparent_diameter_deg`,
    /// `apparent_magnitude_v`, `horizontal_parallax_deg`), matching `swe_pheno` / `swe_pheno_ut`.
    pub fn phenomena(
        &self,
        body: BodyId,
        input: TimeInput,
        options: CalcOptions,
    ) -> Result<PhenomenaResult, EphemerisError> {
        let obs = match options.origin {
            Origin::Topocentric(o) => Some(o),
            _ => None,
        };
        let moment = self.at(input, obs)?;
        self.phenomena_in_moment(&moment, body, options)
    }

    /// Evaluates planetary, solar, lunar, or minor-body illumination and photometric phenomena
    /// explicitly using the fileless analytical [`MoshierAnalyticBackend`] (`SEFLG_MOSEPH`).
    pub fn phenomena_moshier(
        &self,
        body: BodyId,
        input: TimeInput,
        options: CalcOptions,
    ) -> Result<PhenomenaResult, EphemerisError> {
        let obs = match options.origin {
            Origin::Topocentric(o) => Some(o),
            _ => None,
        };
        let moment = self.at_moshier(input, obs)?;
        self.phenomena_in_moment(&moment, body, options)
    }

    /// Evaluates planetary, solar, lunar, or minor-body phenomena inside an already constructed [`Moment`].
    pub fn phenomena_in_moment(
        &self,
        moment: &Moment<'_>,
        body: BodyId,
        options: CalcOptions,
    ) -> Result<PhenomenaResult, EphemerisError> {
        if let BodyId::FixedStar(hip) = body {
            if let Some(star) = self
                .custom_stars
                .iter()
                .find(|s| s.hip == hip.0)
                .or_else(|| find_star_by_hip(hip))
            {
                return crate::phenomena::calc_star_entry_phenomena_in_moment(
                    moment, star, &options,
                );
            }
        }
        calc_phenomena_in_moment(moment, body, &options)
    }

    /// Evaluates planetary or lunar orbital nodes (ascending, descending) and apsides
    /// (perihelion/perigee, aphelion/apogee or second focal point), matching `swe_nod_aps` / `swe_nod_aps_ut`.
    pub fn nod_aps(
        &self,
        body: BodyId,
        input: TimeInput,
        options: CalcOptions,
        method_bits: u32,
    ) -> Result<NodesApsidesResult, EphemerisError> {
        calc_nodes_apsides(self.spk(), &self.time_table, body, input, options, method_bits)
    }

    /// Evaluates planetary or lunar orbital nodes and apsides explicitly using the fileless analytical
    /// [`MoshierAnalyticBackend`] (`SEFLG_MOSEPH`).
    pub fn nod_aps_moshier(
        &self,
        body: BodyId,
        input: TimeInput,
        options: CalcOptions,
        method_bits: u32,
    ) -> Result<NodesApsidesResult, EphemerisError> {
        calc_nodes_apsides(
            &self.moshier,
            &self.time_table,
            body,
            input,
            options,
            method_bits,
        )
    }

    /// Searches forward (`backward == false`) or backward (`backward == true`) from `start_input`
    /// for the next global Lunar Occultation of `target` (`FixedStar`, `Planet`, or `MinorPlanet`),
    /// matching `swe_lun_occult_when_glob`.
    pub fn next_lunar_occultation(
        &self,
        target: BodyId,
        start_input: TimeInput,
        backward: bool,
        max_search_days: f64,
    ) -> Result<LunarOccultationInfo, EphemerisError> {
        let (ut1, _, _, _) = self.time_table.resolve_scales(start_input, None);
        let start_jd = ut1.jd_f64();
        let search_on = |spk: &dyn SpkSource| -> Result<Option<LunarOccultationInfo>, EphemerisError> {
            if backward {
                let mut chunk_end = start_jd;
                let mut result = None;
                let chunk_size = 60.0_f64;
                while chunk_end > start_jd - max_search_days {
                    let chunk_start = (chunk_end - chunk_size).max(start_jd - max_search_days);
                    let events = lunar_occultations_in_range(
                        spk,
                        &self.time_table,
                        target,
                        chunk_start,
                        chunk_end,
                    )?;
                    if let Some(last) = events.into_iter().rfind(|e| e.jd_ut1_max < start_jd) {
                        result = Some(last);
                        break;
                    }
                    chunk_end = chunk_start;
                }
                Ok(result)
            } else {
                next_lunar_occultation(
                    spk,
                    &self.time_table,
                    target,
                    start_jd,
                    max_search_days,
                )
            }
        };
        let found = match search_on(self.spk()) {
            Ok(res) => res,
            Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                search_on(&self.moshier).map_err(|_| orig_err)?
            }
            Err(e) => return Err(e),
        };
        found.ok_or_else(|| EphemerisError::NonConvergence {
            solver: "next_lunar_occultation".to_owned(),
            iterations: (max_search_days * 4.0) as u32,
            residual: 0.0,
        })
    }

    /// Searches forward (`backward == false`) or backward (`backward == true`) from `start_input`
    /// for the next locally visible Lunar Occultation of `target` at `observer`,
    /// matching `swe_lun_occult_when_loc`.
    pub fn next_local_lunar_occultation(
        &self,
        target: BodyId,
        observer: Observer,
        start_input: TimeInput,
        backward: bool,
        max_search_days: f64,
    ) -> Result<LunarOccultationInfo, EphemerisError> {
        let (ut1, _, _, _) = self.time_table.resolve_scales(start_input, Some(observer));
        let start_jd = ut1.jd_f64();
        let search_on = |spk: &dyn SpkSource| -> Result<Option<LunarOccultationInfo>, EphemerisError> {
            if backward {
                let mut chunk_end = start_jd;
                let mut result = None;
                let chunk_size = 60.0_f64;
                while chunk_end > start_jd - max_search_days {
                    let chunk_start = (chunk_end - chunk_size).max(start_jd - max_search_days);
                    let events = local_lunar_occultations_in_range(
                        spk,
                        &self.time_table,
                        target,
                        chunk_start,
                        chunk_end,
                        observer,
                    )?;
                    if let Some(last) = events.into_iter().rfind(|e| e.jd_ut1_max < start_jd) {
                        result = Some(last);
                        break;
                    }
                    chunk_end = chunk_start;
                }
                Ok(result)
            } else {
                next_local_lunar_occultation(
                    spk,
                    &self.time_table,
                    target,
                    start_jd,
                    max_search_days,
                    observer,
                )
            }
        };
        let found = match search_on(self.spk()) {
            Ok(res) => res,
            Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                search_on(&self.moshier).map_err(|_| orig_err)?
            }
            Err(e) => return Err(e),
        };
        found.ok_or_else(|| EphemerisError::NonConvergence {
            solver: "next_local_lunar_occultation".to_owned(),
            iterations: (max_search_days * 4.0) as u32,
            residual: 0.0,
        })
    }

    /// Finds all global Lunar Occultations of `target` within `[start_input, end_input]`.
    pub fn lunar_occultations_in_range(
        &self,
        target: BodyId,
        start_input: TimeInput,
        end_input: TimeInput,
    ) -> Result<Vec<LunarOccultationInfo>, EphemerisError> {
        let (ut1_start, _, _, _) = self.time_table.resolve_scales(start_input, None);
        let (ut1_end, _, _, _) = self.time_table.resolve_scales(end_input, None);
        lunar_occultations_in_range(
            self.spk(),
            &self.time_table,
            target,
            ut1_start.jd_f64(),
            ut1_end.jd_f64(),
        )
    }

    /// Finds all locally visible Lunar Occultations of `target` at `observer` within `[start_input, end_input]`.
    pub fn local_lunar_occultations_in_range(
        &self,
        target: BodyId,
        observer: Observer,
        start_input: TimeInput,
        end_input: TimeInput,
    ) -> Result<Vec<LunarOccultationInfo>, EphemerisError> {
        let (ut1_start, _, _, _) = self.time_table.resolve_scales(start_input, Some(observer));
        let (ut1_end, _, _, _) = self.time_table.resolve_scales(end_input, Some(observer));
        local_lunar_occultations_in_range(
            self.spk(),
            &self.time_table,
            target,
            ut1_start.jd_f64(),
            ut1_end.jd_f64(),
            observer,
        )
    }

    /// Computes the 3D fractional mundane house position $H \in [1.0, 13.0)$ of a sky point
    /// `(ecl_lon_deg, ecl_lat_deg)` at `(armc_deg, latitude_deg, obliquity_deg)` for `system`,
    /// matching `swe_house_pos`.
    pub fn house_pos(
        &self,
        armc_deg: f64,
        latitude_deg: f64,
        obliquity_deg: f64,
        system: HouseSystem,
        ecl_lon_deg: f64,
        ecl_lat_deg: f64,
    ) -> Result<f64, EphemerisError> {
        house_pos(
            armc_deg,
            latitude_deg,
            obliquity_deg,
            system,
            ecl_lon_deg,
            ecl_lat_deg,
        )
    }

    /// Evaluates the Equation of Time ($E = \text{Apparent Solar Time} - \text{Mean Solar Time}$)
    /// in **days** at `input`, matching `swe_time_equ`.
    pub fn time_equation(&self, input: TimeInput) -> Result<f64, EphemerisError> {
        let moment = match self.at(input, None) {
            Ok(m) => m,
            Err(orig_err @ EphemerisError::OutOfCoverage { .. }) => {
                self.at_moshier(input, None).map_err(|_| orig_err)?
            }
            Err(e) => return Err(e),
        };
        time_equation_days(&moment)
    }

    /// Inverse horizontal coordinate transformation (`swe_azalt_rev`): converts geometric
    /// (`airless`) topocentric `(azimuth_north_east_deg, altitude_true_deg)` back into True Equatorial
    /// `(ra_deg, dec_deg)` or True Ecliptic `(lon_deg, lat_deg)` at `epoch` for `observer`.
    #[must_use]
    pub fn azalt_rev(
        &self,
        observer: Observer,
        epoch: TimeInput,
        azimuth_north_east_deg: f64,
        altitude_true_deg: f64,
        target_plane: Plane,
    ) -> [f64; 2] {
        let (ut1, tt, _, _) = self.time_table.resolve_scales(epoch, Some(observer));
        let ft = FrameTransform::at_tt(tt);
        let gast = crate::time::gast_rad(ut1, tt);
        azalt_rev(
            observer,
            gast,
            ft.true_eps_rad,
            azimuth_north_east_deg,
            altitude_true_deg,
            target_plane,
        )
    }

    /// Finds the first instant starting from `start_input` (moving forward if `forward == true`,
    /// or backward if `forward == false`, up to `max_search_days`) when `body`'s ecliptic longitude
    /// reaches `target_lon_deg.rem_euclid(360.0)` under `options` (`R9.1`, `R9.2`).
    pub fn search_body_longitude(
        &self,
        body: BodyId,
        start_input: TimeInput,
        target_lon_deg: f64,
        forward: bool,
        max_search_days: f64,
        options: &CalcOptions,
    ) -> Result<LongitudeCrossingEvent, EphemerisError> {
        search_body_longitude(
            self,
            body,
            start_input,
            target_lon_deg,
            forward,
            max_search_days,
            options,
        )
    }

    /// Finds **all** crossings of `target_lon_deg.rem_euclid(360.0)` by `body` in `[start_input, end_input]`
    /// in chronological order, including 3-pass retrograde triple returns (`R9.1`, `R9.2`).
    pub fn body_longitude_crossings_in_range(
        &self,
        body: BodyId,
        start_input: TimeInput,
        end_input: TimeInput,
        target_lon_deg: f64,
        options: &CalcOptions,
    ) -> Result<Vec<LongitudeCrossingEvent>, EphemerisError> {
        body_longitude_crossings_in_range(
            self,
            body,
            start_input,
            end_input,
            target_lon_deg,
            options,
        )
    }

    /// Finds all Chinese 24 Solar Terms (where **geocentric apparent tropical** solar longitude
    /// $\lambda_\odot = k \cdot 15^\circ$) in `[start_input, end_input]` in chronological order (`R12`).
    pub fn solar_terms_in_range(
        &self,
        start_input: TimeInput,
        end_input: TimeInput,
    ) -> Result<Vec<SolarTermEvent>, EphemerisError> {
        solar_terms_in_range(self, start_input, end_input)
    }

    /// Computes all 4 cardinal seasons (`0°`, `90°`, `180°`, `270°` geocentric apparent tropical solar longitude)
    /// for `year`.
    pub fn seasons_in_year(&self, year: i32) -> Result<SeasonsOfYear, EphemerisError> {
        seasons_in_year(self, year)
    }

    /// Computes the exact Human Design "Design" prenatal epoch (`R9.2`): evaluates the Sun's
    /// geocentric apparent tropical longitude $\lambda_{\text{birth}}$ at `birth_input` and solves
    /// backward (~88–89 days) for the exact instant when $\lambda_\odot(t_{\text{design}}) = (\lambda_{\text{birth}} - 88^\circ) \pmod{360^\circ}$.
    pub fn human_design_epoch(
        &self,
        birth_input: TimeInput,
    ) -> Result<LongitudeCrossingEvent, EphemerisError> {
        human_design_epoch(self, birth_input)
    }

    /// Finds the next (`forward == true`) or previous (`forward == false`) instant from `start_input`
    /// when geocentric apparent tropical lunar-solar elongation $(\lambda_{\text{Moon}} - \lambda_\odot) \pmod{360^\circ}$
    /// reaches `target_elongation_deg.rem_euclid(360.0)` (`R9.2`).
    pub fn search_lunar_phase(
        &self,
        start_input: TimeInput,
        target_elongation_deg: f64,
        forward: bool,
    ) -> Result<LunarPhaseEvent, EphemerisError> {
        search_lunar_phase(self, start_input, target_elongation_deg, forward)
    }

    /// Returns all 4 principal lunar phases (`NewMoon`, `FirstQuarter`, `FullMoon`, `LastQuarter`)
    /// occurring within `[start_input, end_input]` in chronological order (`R9.2`).
    pub fn lunar_phases_in_range(
        &self,
        start_input: TimeInput,
        end_input: TimeInput,
    ) -> Result<Vec<LunarPhaseEvent>, EphemerisError> {
        lunar_phases_in_range(self, start_input, end_input)
    }

    /// Computes the complete [`LunarPhaseSummary`] at `input` (`R0`, `R2.2`, `R9.2`).
    pub fn lunar_phase_summary(&self, input: TimeInput) -> Result<LunarPhaseSummary, EphemerisError> {
        lunar_phase_summary(self, input)
    }

    /// Finds all retrograde and direct stations ($\dot{\lambda}(t) = 0$) of `body` within `[start_input, end_input]`
    /// under `options` (`R9.2`).
    pub fn search_stations_in_range(
        &self,
        body: BodyId,
        start_input: TimeInput,
        end_input: TimeInput,
        options: &CalcOptions,
    ) -> Result<Vec<StationEvent>, EphemerisError> {
        search_stations_in_range(self, body, start_input, end_input, options)
    }

    /// Finds all northern and southern declination extrema ($\dot{\delta}(t) = 0$) of `body`
    /// within `[start_input, end_input]` (`R9.2`).
    pub fn search_declination_extrema_in_range(
        &self,
        body: BodyId,
        start_input: TimeInput,
        end_input: TimeInput,
        options: &CalcOptions,
    ) -> Result<Vec<DeclinationExtremumEvent>, EphemerisError> {
        search_declination_extrema_in_range(self, body, start_input, end_input, options)
    }

    /// Computes Bernadette Brady diurnal fixed-star parans (`R6.2`) for `observer` starting at `start_input`.
    pub fn find_parans(
        &self,
        start_input: TimeInput,
        observer: Observer,
        options: &ParanOptions,
    ) -> Result<Vec<ParanEvent>, EphemerisError> {
        find_parans(self, start_input, observer, options)
    }

    // =========================================================================
    // Wave H: Planetary Returns, Aspects, Transits, VoC Moon, Panchanga & Almanac
    // =========================================================================

    /// Finds the next exact planetary return (`Solar Return`, `Lunar Return`, or any planetary/minor-body return)
    /// after `search_after_input` for `body` relative to its natal position at `natal_input` under `options`
    /// (supports both Tropical and any of the 44 Sidereal / Precession-Corrected zodiac models).
    pub fn planetary_return(
        &self,
        body: BodyId,
        natal_input: TimeInput,
        search_after_input: TimeInput,
        options: CalcOptions,
    ) -> Result<PlanetaryReturnEvent, EphemerisError> {
        planetary_return(self, body, natal_input, search_after_input, options)
    }

    /// Finds all exact planetary returns of `body` within `[start_input, end_input]` relative to its
    /// natal position at `natal_input` under `options`, numbering retrograde triple-pass series (`pass_index: 1..=3`).
    pub fn planetary_returns_in_range(
        &self,
        body: BodyId,
        natal_input: TimeInput,
        start_input: TimeInput,
        end_input: TimeInput,
        options: CalcOptions,
    ) -> Result<Vec<PlanetaryReturnEvent>, EphemerisError> {
        planetary_returns_in_range(self, body, natal_input, start_input, end_input, options)
    }

    /// Evaluates static aspects among `bodies` at `input` using `preset` and `orb_scale`,
    /// including analytical `is_applying` (`true` when approaching exactness).
    pub fn calc_chart_aspects(
        &self,
        input: TimeInput,
        bodies: &[BodyId],
        preset: AspectPreset,
        orb_scale: f64,
        options: CalcOptions,
    ) -> Result<Vec<AspectHit>, EphemerisError> {
        calc_chart_aspects(self, input, bodies, preset, orb_scale, options)
    }

    /// Evaluates cross-chart (synastry or transit-to-natal) aspects between `bodies1` at `input1`
    /// and `bodies2` at `input2`.
    #[allow(clippy::too_many_arguments)]
    pub fn calc_synastry_aspects(
        &self,
        input1: TimeInput,
        bodies1: &[BodyId],
        input2: TimeInput,
        bodies2: &[BodyId],
        preset: AspectPreset,
        orb_scale: f64,
        options: CalcOptions,
    ) -> Result<Vec<AspectHit>, EphemerisError> {
        calc_synastry_aspects(
            self, input1, bodies1, input2, bodies2, preset, orb_scale, options,
        )
    }

    /// Evaluates equatorial declination aspects (`Parallel` and `Contraparallel`) among `bodies` at `input`.
    pub fn calc_declination_aspects(
        &self,
        input: TimeInput,
        bodies: &[BodyId],
        max_orb_deg: f64,
        options: CalcOptions,
    ) -> Result<Vec<DeclinationAspectHit>, EphemerisError> {
        calc_declination_aspects(self, input, bodies, max_orb_deg, options)
    }

    /// Finds the next (`forward == true`) or previous (`forward == false`) exact aspect of `target_angle_deg`
    /// between moving bodies `body1` and `body2` starting from `start_input`.
    pub fn search_exact_aspect(
        &self,
        body1: BodyId,
        body2: BodyId,
        target_angle_deg: f64,
        start_input: TimeInput,
        forward: bool,
        options: CalcOptions,
    ) -> Result<ExactAspectEvent, EphemerisError> {
        search_exact_aspect(
            self,
            body1,
            body2,
            target_angle_deg,
            start_input,
            forward,
            options,
        )
    }

    /// Finds all exact mutual aspects among `bodies` occurring within `[start_input, end_input]`,
    /// including triple passes across retrograde loops.
    pub fn mutual_aspects_in_range(
        &self,
        bodies: &[BodyId],
        preset: AspectPreset,
        start_input: TimeInput,
        end_input: TimeInput,
        options: CalcOptions,
    ) -> Result<Vec<ExactAspectEvent>, EphemerisError> {
        mutual_aspects_in_range(self, bodies, preset, start_input, end_input, options)
    }

    /// Finds all exact transiting aspects from `transiting_bodies` to `natal_bodies` (evaluated at `natal_input`)
    /// within `[start_input, end_input]`, including retrograde triple passes.
    #[allow(clippy::too_many_arguments)]
    pub fn natal_transits_in_range(
        &self,
        transiting_bodies: &[BodyId],
        natal_input: TimeInput,
        natal_bodies: &[BodyId],
        preset: AspectPreset,
        start_input: TimeInput,
        end_input: TimeInput,
        options: CalcOptions,
    ) -> Result<Vec<ExactAspectEvent>, EphemerisError> {
        natal_transits_in_range(
            self,
            transiting_bodies,
            natal_input,
            natal_bodies,
            preset,
            start_input,
            end_input,
            options,
        )
    }

    /// Computes all Void-of-Course Moon intervals within `[start_input, end_input]`.
    pub fn void_of_course_moon_in_range(
        &self,
        start_input: TimeInput,
        end_input: TimeInput,
        major_planets: Option<&[BodyId]>,
        options: CalcOptions,
    ) -> Result<Vec<VoidOfCourseInterval>, EphemerisError> {
        void_of_course_moon_in_range(self, start_input, end_input, major_planets, options)
    }

    /// Computes the complete Vedic Panchanga (`Vara`, `Tithi`, `Nakshatra + Pada`, `Yoga`, `Karana`,
    /// `Rahukalam`, `Yamagandam`, `Gulikakalam`, `Abhijit Muhurta`) at `input` for `observer` and `ayanamsa`.
    pub fn panchanga(
        &self,
        input: TimeInput,
        observer: Observer,
        ayanamsa: AyanamsaModel,
    ) -> Result<PanchangaReport, EphemerisError> {
        calc_panchanga(self, input, observer, ayanamsa)
    }

    /// Computes the 120-year Vedic Vimshottari Mahadasha and Antardasha schedule at `birth_input` for `ayanamsa`.
    pub fn vimshottari_dasha(
        &self,
        birth_input: TimeInput,
        ayanamsa: AyanamsaModel,
        year_length_days: Option<f64>,
    ) -> Result<VimshottariDashaReport, EphemerisError> {
        calc_vimshottari_dasha(self, birth_input, ayanamsa, year_length_days)
    }

    /// Evaluates a Parashara divisional chart (`Varga`, `D1..=D60`) placement for `sidereal_lon_deg`.
    pub fn calc_varga(
        &self,
        sidereal_lon_deg: f64,
        division: u16,
    ) -> Result<VargaPlacement, EphemerisError> {
        calc_varga(sidereal_lon_deg, division)
    }

    /// Evaluates all 16 classical Parashara Shodashavarga placements (`D1..D60`) for `sidereal_lon_deg`.
    pub fn calc_shodashavarga(
        &self,
        sidereal_lon_deg: f64,
    ) -> Result<Vec<VargaPlacement>, EphemerisError> {
        calc_shodashavarga(sidereal_lon_deg)
    }

    /// Evaluates traditional chart Sect (`Day` vs `Night` chart and `near_horizon` flag) at `input` for `observer`.
    pub fn chart_sect(
        &self,
        input: TimeInput,
        observer: Observer,
    ) -> Result<ChartSectInfo, EphemerisError> {
        chart_sect(self, input, observer)
    }

    /// Computes the 24 unequal Chaldean Planetary Hours anchored to local sunrise and sunset at `input` for `observer`.
    pub fn planetary_hours(
        &self,
        input: TimeInput,
        observer: Observer,
    ) -> Result<PlanetaryHoursReport, EphemerisError> {
        planetary_hours(self, input, observer)
    }

    /// Computes the 7 classical Hermetic / Arabic Lots (*Fortune*, *Spirit*, *Eros*, *Necessity*,
    /// *Courage*, *Victory*, *Nemesis*) at `input` for `observer`.
    pub fn calc_lots(
        &self,
        input: TimeInput,
        observer: Observer,
        house_options: HouseOptions,
    ) -> Result<Vec<HermeticLot>, EphemerisError> {
        calc_lots(self, input, observer, house_options)
    }

    /// Computes a complete one-call [`ChartSnapshot`] (all requested `bodies` with zodiac/nakshatra/HD gate
    /// placement, 2D and 3D house positions, 12 houses + 8 angles, static aspects, lunar phase summary,
    /// chart sect, and 7 Hermetic lots) at `input` for `observer`.
    #[allow(clippy::too_many_arguments)]
    pub fn chart_snapshot(
        &self,
        input: TimeInput,
        observer: Observer,
        bodies: &[BodyId],
        house_options: HouseOptions,
        calc_options: CalcOptions,
        aspect_preset: AspectPreset,
        orb_scale: f64,
    ) -> Result<ChartSnapshot, EphemerisError> {
        let mut opts = calc_options;
        opts.include_speed = true;
        let obs_opt = match opts.origin {
            Origin::Topocentric(o) => Some(o),
            _ => Some(observer),
        };
        let moment = self.at(input, obs_opt)?;
        let houses_res = self.houses_in_moment(&moment, observer, house_options)?;
        let armc_deg = houses_res.angles.armc_deg;
        let eps_deg = moment.center().frame.true_eps_rad.to_degrees();
        let calc_aya_deg = match opts.zodiac {
            Zodiac::Sidereal(m) => self.ayanamsa(input, m, opts.nutation)?.effective_deg,
            Zodiac::Tropical => 0.0,
        };
        let house_aya_deg = match house_options.zodiac {
            Zodiac::Sidereal(m) => {
                self.ayanamsa(input, m, house_options.nutation)?
                    .effective_deg
            }
            Zodiac::Tropical => 0.0,
        };
        let eff_sys = houses_res.effective_house_system;
        let is_3d_mundane = matches!(
            eff_sys,
            HouseSystem::Placidus
                | HouseSystem::PolichPage
                | HouseSystem::Regiomontanus
                | HouseSystem::Campanus
                | HouseSystem::Azimuthal
                | HouseSystem::Meridian
                | HouseSystem::Morinus
                | HouseSystem::Alcabitius
                | HouseSystem::Koch
        );

        let mut body_snapshots = Vec::with_capacity(bodies.len());
        let mut samples = Vec::with_capacity(bodies.len());
        for &body in bodies {
            let c = self.calc_in_moment(&moment, body, opts)?;
            let placement = map_zodiac_placement(c.lon_deg);
            let trop_lon_deg = (c.lon_deg + calc_aya_deg).rem_euclid(360.0);
            let house_lon_deg = (trop_lon_deg - house_aya_deg).rem_euclid(360.0);
            let h_num = crate::houses::house_of_lon(house_lon_deg, &houses_res.cusps_deg, 1.0e-9);
            let h_pos_3d = if is_3d_mundane {
                house_pos(
                    armc_deg,
                    observer.latitude_deg,
                    eps_deg,
                    eff_sys,
                    trop_lon_deg,
                    c.lat_deg,
                )
                .unwrap_or_else(|_| {
                    crate::houses::interp_cusps_lon(house_lon_deg, &houses_res.cusps_deg)
                })
            } else {
                crate::houses::interp_cusps_lon(house_lon_deg, &houses_res.cusps_deg)
            };
            samples.push(crate::aspects::BodyPositionSample {
                body,
                lon_deg: c.lon_deg,
                lat_deg: c.lat_deg,
                lon_speed_deg_day: c.lon_speed_deg_day,
                dec_deg: None,
                dec_speed_deg_day: None,
            });
            body_snapshots.push(ChartBodySnapshot {
                body,
                body_name: crate::events::body_display_name(body),
                is_retrograde: c.lon_speed_deg_day < 0.0,
                calc: c,
                placement,
                house_number: h_num,
                house_pos_3d: h_pos_3d,
            });
        }

        let defs = crate::aspects::aspect_definitions(aspect_preset);
        let aspects = crate::aspects::calc_aspects(&samples, &defs, orb_scale);
        let lunar_phase = self.lunar_phase_summary(input)?;
        let sect = self.chart_sect(input, observer)?;
        let lots = self.calc_lots(input, observer, house_options)?;

        Ok(ChartSnapshot {
            jd_ut1: moment.center().ut1.jd_f64(),
            jd_tt: moment.center().tt.jd_f64(),
            observer,
            bodies: body_snapshots,
            houses: houses_res,
            aspects,
            lunar_phase,
            sect,
            lots,
        })
    }
}

/// Resolves observer barycentric position (km), velocity (km/s), and center ID for `origin`.
fn resolve_origin_bary_state(
    ctx: &crate::moment::EpochContext,
    origin: Origin,
) -> (Vec3, Vec3, i32) {
    match origin {
        Origin::Geocentric => (
            ctx.earth_bary_km.pos_km,
            ctx.earth_bary_km.vel_km_s,
            399,
        ),
        Origin::Topocentric(obs) => {
            if ctx.observer == Some(obs) {
                (
                    ctx.observer_bary_km.pos_km,
                    ctx.observer_bary_km.vel_km_s,
                    399,
                )
            } else {
                let lon_rad = obs.longitude_deg.to_radians();
                let lat_rad = obs.latitude_deg.to_radians();
                let mut pv_m = [[0.0_f64; 3]; 2];
                sofars::astro::pvtob(
                    lon_rad,
                    lat_rad,
                    obs.altitude_m,
                    0.0,
                    0.0,
                    0.0,
                    ctx.gast_rad,
                    &mut pv_m,
                );
                let pos_true_eq_km = [
                    pv_m[0][0] / 1000.0,
                    pv_m[0][1] / 1000.0,
                    pv_m[0][2] / 1000.0,
                ];
                let vel_true_eq_km_s = [
                    pv_m[1][0] / 1000.0,
                    pv_m[1][1] / 1000.0,
                    pv_m[1][2] / 1000.0,
                ];
                let r_bpn_t = mat_transpose(&ctx.frame.r_bpn);
                let geo_pos_km = mat_vec(&r_bpn_t, &pos_true_eq_km);
                let geo_vel_km_s = mat_vec(&r_bpn_t, &vel_true_eq_km_s);
                (
                    vec_add(&ctx.earth_bary_km.pos_km, &geo_pos_km),
                    vec_add(&ctx.earth_bary_km.vel_km_s, &geo_vel_km_s),
                    399,
                )
            }
        }
        Origin::Heliocentric => (ctx.sun_bary_km.pos_km, ctx.sun_bary_km.vel_km_s, 10),
        Origin::Barycentric => ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0),
    }
}

/// Evaluates an SPK-backed body using single-pass analytical Chebyshev velocity propagation
/// (`calc_spk_body_fast`), avoiding multi-point stencil evaluations.
fn calc_spk_body_analytical(
    moment: &Moment<'_>,
    body: BodyId,
    options: &CalcOptions,
) -> Result<CalcResult, EphemerisError> {
    let (target_id, target_kind, _) = resolve_spk_target(body)?;
    let center = moment.center();
    let effective_origin = if options.plane == Plane::Horizontal
        && matches!(options.origin, Origin::Geocentric)
    {
        center.observer.map_or(Origin::Geocentric, Origin::Topocentric)
    } else {
        options.origin
    };
    let (p_obs_km, v_obs_kms, center_id) = resolve_origin_bary_state(center, effective_origin);

    // 1. Position and converged light-time via single-pass SPK lookup
    let (pos_0, tau_0) = eval_spk_target_vector_au_seeded(
        moment.spk(),
        center,
        target_id,
        options,
        None,
    )?;

    if !options.include_speed {
        let c = cart_to_spherical(&pos_0);
        let mut sph = SphericalState {
            lon_deg: c.lon_deg,
            lat_deg: c.lat_deg,
            dist_au: c.dist,
            lon_speed_deg_day: 0.0,
            lat_speed_deg_day: 0.0,
            dist_speed_au_day: 0.0,
        };
        let mut cart_pos = pos_0;
        let mut cart_vel = [0.0, 0.0, 0.0];
        if let Some((aya_deg, _)) = moment.sidereal_offset_and_rate(options)? {
            sph.lon_deg = (sph.lon_deg - aya_deg).rem_euclid(360.0);
            let (p_sid, v_sid) = spherical_state_to_cart(&sph, false);
            cart_pos = p_sid;
            cart_vel = v_sid;
        }
        return Ok(CalcResult {
            body,
            target_id,
            center_id,
            target_kind,
            jd_tt: center.tt.jd_f64(),
            jd_ut1: center.ut1.jd_f64(),
            jd_tdb: center.tdb.jd_f64(),
            delta_t_sec: center.delta_t_info.delta_t_sec,
            time_extrapolated: center.delta_t_info.time_extrapolated,
            time_table_version: center.delta_t_info.time_table_version.clone(),
            lon_deg: sph.lon_deg,
            lat_deg: sph.lat_deg,
            dist_au: sph.dist_au,
            dist_km: sph.dist_au * AU_KM,
            lon_speed_deg_day: 0.0,
            lat_speed_deg_day: 0.0,
            dist_speed_au_day: 0.0,
            dist_speed_km_day: 0.0,
            magnitude_v: None,
            cartesian_au: cart_pos,
            cartesian_vel_au_day: cart_vel,
            warnings: Vec::new(),
        });
    }

    // 2. Query target state at retarded epoch tdb - tau_0
    let st_targ = moment.spk().state_icrf_km_split(
        target_id,
        0,
        center.tdb.0.jd1,
        center.tdb.0.jd2 - tau_0,
    )?;

    let rho_vec_km = vec_sub(&st_targ.pos_km, &p_obs_km);
    let rho_km = vec_norm(&rho_vec_km);
    let e_hat = if rho_km > 1e-12 {
        vec_scale(&rho_vec_km, 1.0 / rho_km)
    } else {
        [0.0, 0.0, 0.0]
    };

    let v_rel_kms = vec_sub(&st_targ.vel_km_s, &v_obs_kms);
    let denom = C_KM_S + vec_dot(&e_hat, &st_targ.vel_km_s);
    let dtau_dt = if denom.abs() > 1e-12 {
        vec_dot(&e_hat, &v_rel_kms) / denom
    } else {
        0.0
    };

    // Astrometric velocity in km/s: drho/dt = v_targ * (1 - dtau/dt) - v_obs
    let v_astrom_kms = vec_sub(
        &vec_scale(&st_targ.vel_km_s, 1.0 - dtau_dt),
        &v_obs_kms,
    );

    // 3. Apparent position and velocity in ICRF (including deflection & stellar aberration rate)
    let km_s_to_au_day = SECONDS_PER_DAY / AU_KM;
    let r_astrom_au = vec_scale(&rho_vec_km, 1.0 / AU_KM);
    let v_astrom_au_day = vec_scale(&v_astrom_kms, km_s_to_au_day);

    let pos_icrf_corrected = apply_aberration_and_deflection_icrf(
        &r_astrom_au,
        &st_targ.pos_km,
        &p_obs_km,
        &v_obs_kms,
        center,
        target_id,
        effective_origin,
        options.aberration,
        tau_0,
    );

    let v_icrf_au_day = match options.aberration {
        AbMode::Geometric => {
            vec_scale(&v_rel_kms, km_s_to_au_day)
        }
        AbMode::Astrometric => {
            v_astrom_au_day
        }
        AbMode::Apparent | AbMode::NoAberration | AbMode::NoDeflection => {
            let is_earth_based = matches!(effective_origin, Origin::Geocentric | Origin::Topocentric(_));
            if is_earth_based && rho_km > 1e-12 {
                // Earth orbital acceleration: a_earth = -mu_sun / r_sun^3 * r_sun_vec
                let mu_sun = 1.327_124_400_419_393_8e11_f64; // km^3/s^2
                let r_se_km = vec_sub(&p_obs_km, &center.sun_bary_km.pos_km);
                let d_se_km = vec_norm(&r_se_km);
                let a_earth_kms2 = if d_se_km > 1e6 {
                    vec_scale(&r_se_km, -mu_sun / (d_se_km * d_se_km * d_se_km))
                } else {
                    [0.0, 0.0, 0.0]
                };

                // Topocentric centrifugal acceleration in ICRF: a_topo = -omega^2 * (r_topo - (r_topo . z_cip) * z_cip)
                let a_topo_kms2 = if let Origin::Topocentric(_) = effective_origin {
                    let r_topo_km = vec_sub(&p_obs_km, &center.earth_bary_km.pos_km);
                    let omega = 7.292_115_146_706_979e-5_f64; // rad/s
                    let omega2 = omega * omega;
                    let z_cip = center.frame.r_bpn[2];
                    let z_proj = vec_dot(&r_topo_km, &z_cip);
                    let r_perp_km = [
                        r_topo_km[0] - z_proj * z_cip[0],
                        r_topo_km[1] - z_proj * z_cip[1],
                        r_topo_km[2] - z_proj * z_cip[2],
                    ];
                    vec_scale(&r_perp_km, -omega2)
                } else {
                    [0.0, 0.0, 0.0]
                };

                let a_obs_kms2 = vec_add(&a_earth_kms2, &a_topo_kms2);

                // Differentiate the exact deflection + SOFA aberration pipeline across +/- 1e-4 days (8.64 s)
                let h_dt_days = 1.0e-4_f64;
                let h_dt_s = h_dt_days * SECONDS_PER_DAY;

                let r_astrom_p_au = vec_add(&r_astrom_au, &vec_scale(&v_astrom_au_day, h_dt_days));
                let r_astrom_m_au = vec_sub(&r_astrom_au, &vec_scale(&v_astrom_au_day, h_dt_days));

                let p_targ_p_km = vec_add(&st_targ.pos_km, &vec_scale(&st_targ.vel_km_s, h_dt_s * (1.0 - dtau_dt)));
                let p_targ_m_km = vec_sub(&st_targ.pos_km, &vec_scale(&st_targ.vel_km_s, h_dt_s * (1.0 - dtau_dt)));

                let p_obs_p_km = vec_add(&p_obs_km, &vec_scale(&v_obs_kms, h_dt_s));
                let p_obs_m_km = vec_sub(&p_obs_km, &vec_scale(&v_obs_kms, h_dt_s));

                let v_obs_p_kms = vec_add(&v_obs_kms, &vec_scale(&a_obs_kms2, h_dt_s));
                let v_obs_m_kms = vec_sub(&v_obs_kms, &vec_scale(&a_obs_kms2, h_dt_s));

                let pos_corr_p = apply_aberration_and_deflection_icrf(
                    &r_astrom_p_au,
                    &p_targ_p_km,
                    &p_obs_p_km,
                    &v_obs_p_kms,
                    center,
                    target_id,
                    effective_origin,
                    options.aberration,
                    tau_0,
                );
                let pos_corr_m = apply_aberration_and_deflection_icrf(
                    &r_astrom_m_au,
                    &p_targ_m_km,
                    &p_obs_m_km,
                    &v_obs_m_kms,
                    center,
                    target_id,
                    effective_origin,
                    options.aberration,
                    tau_0,
                );

                let inv_2h = 1.0 / (2.0 * h_dt_days);
                [
                    (pos_corr_p[0] - pos_corr_m[0]) * inv_2h,
                    (pos_corr_p[1] - pos_corr_m[1]) * inv_2h,
                    (pos_corr_p[2] - pos_corr_m[2]) * inv_2h,
                ]
            } else {
                v_astrom_au_day
            }
        }
    };

    // 4. Output frame rotation: v_frame = R * v_icrf + dR/dt * r_icrf_corrected
    let mut stencil_opts = *options;
    if options.plane == Plane::Ecliptic && matches!(options.zodiac, Zodiac::Sidereal(_)) {
        stencil_opts.nutation = NutationMode::Mean;
    }
    let r_mat = center.frame.icrf_to_frame_matrix(stencil_opts.plane, stencil_opts.equinox, stencil_opts.nutation);
    let dr_mat = center.frame.icrf_to_frame_matrix_rate(stencil_opts.plane, stencil_opts.equinox, stencil_opts.nutation);
    let v_rot_term = mat_vec(&dr_mat, &pos_icrf_corrected);
    let v_frame_term = mat_vec(&r_mat, &v_icrf_au_day);
    let vel_0 = vec_add(&v_frame_term, &v_rot_term);

    // 5. Convert to spherical coordinates
    let pos_for_rate = if stencil_opts.nutation != options.nutation {
        eval_spk_target_vector_au_seeded(moment.spk(), center, target_id, &stencil_opts, Some(tau_0))?.0
    } else {
        pos_0
    };
    let mut sph = cart_state_to_spherical(&pos_for_rate, &vel_0);
    if stencil_opts.nutation != options.nutation {
        let c0 = cart_to_spherical(&pos_0);
        sph.lon_deg = c0.lon_deg;
        sph.lat_deg = c0.lat_deg;
        sph.dist_au = c0.dist;
    }

    // 6. Sidereal offset & rate
    let mut cart_pos = pos_0;
    let mut cart_vel = vel_0;
    if let Some((aya_deg, aya_rate_deg_day)) = moment.sidereal_offset_and_rate(options)? {
        sph.lon_deg = (sph.lon_deg - aya_deg).rem_euclid(360.0);
        sph.lon_speed_deg_day -= aya_rate_deg_day;
        let (p_sid, v_sid) = spherical_state_to_cart(&sph, true);
        cart_pos = p_sid;
        cart_vel = v_sid;
    }

    // 7. Visual magnitude
    let magnitude_v = {
        let r_helio_au = vec_scale(
            &vec_sub(&st_targ.pos_km, &center.sun_bary_km.pos_km),
            1.0 / AU_KM,
        );
        let delta_obs_au = vec_scale(&rho_vec_km, 1.0 / AU_KM);
        let r_obs_sun_au = vec_scale(
            &vec_sub(&center.sun_bary_km.pos_km, &p_obs_km),
            1.0 / AU_KM,
        );
        compute_apparent_magnitude_v(
            body,
            &r_helio_au,
            &delta_obs_au,
            &r_obs_sun_au,
        )
    };

    Ok(CalcResult {
        body,
        target_id,
        center_id,
        target_kind,
        jd_tt: center.tt.jd_f64(),
        jd_ut1: center.ut1.jd_f64(),
        jd_tdb: center.tdb.jd_f64(),
        delta_t_sec: center.delta_t_info.delta_t_sec,
        time_extrapolated: center.delta_t_info.time_extrapolated,
        time_table_version: center.delta_t_info.time_table_version.clone(),
        lon_deg: sph.lon_deg,
        lat_deg: sph.lat_deg,
        dist_au: sph.dist_au,
        dist_km: sph.dist_au * AU_KM,
        lon_speed_deg_day: sph.lon_speed_deg_day,
        lat_speed_deg_day: sph.lat_speed_deg_day,
        dist_speed_au_day: sph.dist_speed_au_day,
        dist_speed_km_day: sph.dist_speed_au_day * AU_KM,
        magnitude_v,
        cartesian_au: cart_pos,
        cartesian_vel_au_day: cart_vel,
        warnings: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spk::AniseSpkBackend;
    use crate::time::{TimeInput, TimeTable};
    use crate::types::{
        AbMode, BodyId, CalcOptions, Equinox, NutationMode, Observer, Origin, Plane, PlanetId,
        Zodiac, JD_J2000,
    };
    use std::path::PathBuf;

    fn get_test_backend_and_time() -> Option<(AniseSpkBackend, TimeTable)> {
        let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(|p| p.parent())
            .map(|p| p.join("data/de440s.bsp"))?;
        if !p.exists() {
            return None;
        }
        let spk = AniseSpkBackend::from_bsp_paths(&[p]).ok()?;
        let time_table = TimeTable::embedded();
        Some((spk, time_table))
    }

    #[test]
    fn test_calc_spk_body_analytical_matches_stencil() {
        let Some((spk, time_table)) = get_test_backend_and_time() else {
            return;
        };
        let obs = Observer::new(8.54, 47.37, 400.0).expect("valid observer");

        for origin in [Origin::Geocentric, Origin::Topocentric(obs)] {
            for aberration in [AbMode::Apparent, AbMode::NoAberration, AbMode::NoDeflection] {
                let opts = CalcOptions {
                    origin,
                    plane: Plane::Ecliptic,
                    equinox: Equinox::OfDate,
                    nutation: NutationMode::True,
                    aberration,
                    zodiac: Zodiac::Tropical,
                    include_speed: true,
                };

                let moment = Moment::new(
                    &spk,
                    &time_table,
                    TimeInput::from_tt_jd(JD_J2000 + 100.0),
                    match origin {
                        Origin::Topocentric(o) => Some(o),
                        _ => None,
                    },
                )
                .expect("moment");

                // Test Mars and Jupiter
                for planet in [PlanetId::Mars, PlanetId::Jupiter] {
                    let body = BodyId::Planet(planet);
                    let res_ana =
                        calc_spk_body_analytical(&moment, body, &opts).expect("analytical");
                    let res_stencil = calc_body_in_moment(&moment, body, &opts).expect("stencil");

                    // Positions must match to < 1e-10 deg
                    assert!(
                        (res_ana.lon_deg - res_stencil.lon_deg).abs() < 1e-10,
                        "Position mismatch for {planet:?} under {aberration:?}: ana={}, stencil={}",
                        res_ana.lon_deg,
                        res_stencil.lon_deg
                    );

                    // Longitudinal velocity must match 4-point stencil to < 1e-6 deg/day (geocentric) or < 1e-5 deg/day (topocentric)
                    let max_d_speed = match origin {
                        Origin::Topocentric(_) => 1.0e-5,
                        _ => 1.0e-6,
                    };
                    let d_speed = (res_ana.lon_speed_deg_day - res_stencil.lon_speed_deg_day).abs();
                    assert!(
                        d_speed < max_d_speed,
                        "Velocity mismatch for {planet:?} under {aberration:?} ({origin:?}): ana={}, stencil={}, diff={d_speed:e} deg/day",
                        res_ana.lon_speed_deg_day,
                        res_stencil.lon_speed_deg_day
                    );

                    // 3D Cartesian velocity must match to < 5e-8 AU/day (geocentric) or < 1e-6 AU/day (topocentric)
                    let max_cart_diff = match origin {
                        Origin::Topocentric(_) => 1.0e-6,
                        _ => 5.0e-8,
                    };
                    let vx_diff =
                        (res_ana.cartesian_vel_au_day[0] - res_stencil.cartesian_vel_au_day[0]).abs();
                    let vy_diff =
                        (res_ana.cartesian_vel_au_day[1] - res_stencil.cartesian_vel_au_day[1]).abs();
                    let vz_diff =
                        (res_ana.cartesian_vel_au_day[2] - res_stencil.cartesian_vel_au_day[2]).abs();
                    assert!(
                        vx_diff < max_cart_diff && vy_diff < max_cart_diff && vz_diff < max_cart_diff,
                        "Cartesian velocity mismatch for {planet:?} under {aberration:?}: diff=[{vx_diff:e}, {vy_diff:e}, {vz_diff:e}] AU/day"
                    );
                }
            }
        }
    }

    #[test]
    fn test_topocentric_centrifugal_cip_acceleration() {
        let Some((spk, time_table)) = get_test_backend_and_time() else {
            return;
        };

        // Observer at equator: lat = 0 deg -> r_perp = R_earth ~ 6378 km
        let obs_eq = Observer::new(0.0, 0.0, 0.0).expect("equator observer");
        let m_eq = Moment::new(
            &spk,
            &time_table,
            TimeInput::from_tt_jd(JD_J2000),
            Some(obs_eq),
        )
        .expect("moment eq");
        let center_eq = m_eq.center();
        let r_topo_eq =
            vec_sub(&center_eq.observer_bary_km.pos_km, &center_eq.earth_bary_km.pos_km);
        let z_cip = center_eq.frame.r_bpn[2];
        let z_proj_eq = vec_dot(&r_topo_eq, &z_cip);
        let r_perp_eq = [
            r_topo_eq[0] - z_proj_eq * z_cip[0],
            r_topo_eq[1] - z_proj_eq * z_cip[1],
            r_topo_eq[2] - z_proj_eq * z_cip[2],
        ];
        let r_perp_norm = vec_norm(&r_perp_eq);
        assert!(
            (r_perp_norm - 6378.137).abs() < 1.0,
            "Equatorial r_perp must be ~6378 km, got {r_perp_norm}"
        );

        let omega = 7.292_115_146_706_979e-5_f64;
        let a_mag = omega * omega * r_perp_norm; // in km/s^2
        assert!(
            (a_mag - 3.39e-5).abs() < 1.0e-6,
            "Equatorial centrifugal acceleration ~3.39e-5 km/s^2, got {a_mag:e}"
        );

        // Observer at north pole: lat = 90 deg -> r_perp ~ 0
        let obs_pole = Observer::new(90.0, 0.0, 0.0).expect("pole observer");
        let m_pole = Moment::new(
            &spk,
            &time_table,
            TimeInput::from_tt_jd(JD_J2000),
            Some(obs_pole),
        )
        .expect("moment pole");
        let center_pole = m_pole.center();
        let r_topo_pole =
            vec_sub(&center_pole.observer_bary_km.pos_km, &center_pole.earth_bary_km.pos_km);
        let z_proj_pole = vec_dot(&r_topo_pole, &z_cip);
        let r_perp_pole = [
            r_topo_pole[0] - z_proj_pole * z_cip[0],
            r_topo_pole[1] - z_proj_pole * z_cip[1],
            r_topo_pole[2] - z_proj_pole * z_cip[2],
        ];
        let r_perp_pole_norm = vec_norm(&r_perp_pole);
        assert!(
            r_perp_pole_norm < 1.0,
            "Polar r_perp must be < 1 km, got {r_perp_pole_norm}"
        );
    }
}
