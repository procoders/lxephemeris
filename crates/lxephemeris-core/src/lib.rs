//! `lxephemeris-core`: Clean-room, thread-safe astronomical and astrological ephemeris engine.

#![forbid(unsafe_code)]

pub mod almanac;
pub mod aspects;
pub mod astrometry;
pub mod ayanamsa;
pub mod daf;
pub mod eclipses;
pub mod engine;
pub mod error;
pub mod events;
pub mod frame;
pub mod heliacal;
pub mod houses;
pub mod hypothetical;
pub mod license;
pub mod lunar;
pub mod manifest;
pub mod moment;
pub mod moshier;
pub mod nodes_apsides;
pub mod occultations;
pub mod phenomena;
pub mod spk;
pub mod stars;
#[cfg(feature = "swe-compat")]
pub mod swe_compat;
pub mod time;
pub mod types;

pub use almanac::{
    calc_lots, calc_lots_from_positions, calc_panchanga, calc_shodashavarga, calc_varga,
    calc_vimshottari_dasha, chart_sect, karana_name_for_index, map_zodiac_placement,
    normalize_lon_360, planetary_hours, resolve_vedic_day_sun_windows, varga_sanskrit_name,
    vimshottari_dasha, vimshottari_dasha_from_lon, ChartSectInfo, DashaPeriod, DashaSubPeriod,
    HermeticLot, KaranaInfo, NakshatraInfo, ObserverLocation, PanchangaReport, PlanetaryHourSlot,
    PlanetaryHoursReport, TimeWindow, TithiInfo, VaraInfo, VargaPlacement, VimshottariDashaReport,
    YogaInfo, ZodiacPlacement, CHALDEAN_ORDER_7, CHARA_KARANA_NAMES_7, DEFAULT_DASHAS_YEAR_DAYS,
    HD_MANDALA_GATES_64, HD_MANDALA_START_LON_DEG, NAKSHATRA_NAMES_27, TITHI_NAMES_30,
    VARA_LORDS_7, VARA_NAMES_7, VIMSHOTTARI_LORDS_9, VIMSHOTTARI_YEARS_9, YOGA_NAMES_27,
    ZODIAC_SIGN_LORDS, ZODIAC_SIGN_NAMES,
};
pub use aspects::{
    aspect_definitions, calc_aspects, calc_chart_aspects, calc_declination_aspects,
    calc_synastry_aspects, canonical_aspect_name, declination_aspects, mutual_aspects_in_range,
    natal_transits_in_range, search_exact_aspect, synastry_aspects, void_of_course_moon_in_range,
    AspectDefinition, AspectHit, AspectPreset, BodyPositionSample, DeclinationAspectHit,
    DeclinationAspectKind, ExactAspectEvent, VoidOfCourseInterval,
};
pub use astrometry::{
    azalt, azalt_from_ecliptic, azalt_rev, calc_body_in_moment, calc_spk_body_in_moment,
    equatorial_to_horizontal_vector_au, eval_spk_target_vector_au, resolve_spk_target, rise_trans,
    time_equation_days, AtmosphereParams, CalcResult, DiscMode, HorizontalCoords, RiseTransEvent,
    RiseTransOptions, BM_EARTH, BM_JUPITER, BM_SATURN, BM_SUN, C_AU_DAY, C_KM_S,
    LT_CONVERGENCE_DAYS, MAX_LT_ITER,
};
pub use ayanamsa::{
    epoch_model_anchor, eval_ayanamsa, eval_ayanamsa_in_moment, iae_general_precession_arcsec,
    star_model_anchor, AyanamsaValue,
};
pub use eclipses::{
    compute_lunar_saros, compute_solar_saros, lunar_eclipse_at, lunar_eclipses_in_range,
    next_local_solar_eclipse, next_lunar_eclipse, next_solar_eclipse, solar_eclipse_at,
    solar_eclipses_in_range, EclipseOptions, LocalSolarEclipse, LunarEclipseInfo, LunarEclipseKind,
    ShadowEnlargementRule, SolarEclipseInfo, SolarEclipseKind, K_MOON_MEAN, K_MOON_UMBRA, K_SUN,
    R_EARTH_EQ_KM, R_SUN_KM, WGS84_E2, WGS84_FLATTENING,
};
pub use engine::{
    BodyCoverage, ChartBodySnapshot, ChartSnapshot, Engine, EngineCapabilities, EngineConfig,
};
pub use error::{EphemerisError, ErrorCode};
pub use events::{
    body_display_name, body_longitude_crossings_in_range, characteristic_orbital_period_days,
    find_parans, human_design_epoch, lunar_phase_summary, lunar_phases_in_range, planetary_return,
    planetary_returns_in_range, search_body_longitude, search_declination_extrema_in_range,
    search_lunar_phase, search_stations_in_range, seasons_in_year, solar_terms_in_range,
    wrap_pm180, DeclinationExtremumEvent, DeclinationExtremumKind, LongitudeCrossingEvent,
    LunarPhaseEvent, LunarPhaseKind, LunarPhaseSummary, ParanAngle, ParanEvent, ParanOptions,
    PlanetaryReturnEvent, SeasonsOfYear, SolarTermEvent, StationEvent, StationKind,
    LUNAR_PHASE_8_NAMES, SOLAR_TERM_CHINESE_NAMES, SOLAR_TERM_ENGLISH_NAMES,
    SOLAR_TERM_PINYIN_NAMES,
};
pub use frame::{
    cart_state_to_spherical, cart_to_spherical, check_iau2006_validity, is_within_iau2006_domain,
    mat_identity, mat_mul, mat_transpose, mat_vec, rot_x, rot_y, rot_z, spherical_state_to_cart,
    stencil_4pt_derivative, EclNutInfo, FrameTransform, Mat3, SphericalCoords, SphericalState,
    Vec3, JD_IAU2006_MAX, JD_IAU2006_MIN,
};
pub use heliacal::{
    estimate_visual_magnitude, heliacal_phenomena, schaefer_critical_arcus_visionis_deg,
    HeliacalOptions, HeliacalResult,
};
pub use houses::{
    asc1, calc_houses, compute_angles, gauquelin_sectors, house_of_lon, house_pos, houses_armc,
    houses_armc_with_sidereal, is_cusp_ring_monotonic, porphyry_cusps, GauquelinResult,
    HouseAngles, HouseOptions, HousesResult,
};
pub use hypothetical::{
    calc_hypothetical_in_moment, eval_hypothetical_vector_au, hypothetical_elements,
    uranian_heliocentric_icrf_au, HypotheticalElements, GAUSS_K_DEG_DAY,
};
pub use license::{
    clear_commercial_license, ed25519_verify, is_valid_iso_date,
    maybe_emit_production_license_notice, register_commercial_license, resolve_license_info,
    verify_license_key, LicenseInfo, LicenseStatus, ENGINE_RELEASE_DATE, PUBLIC_KEY_K1_2026,
};
pub use lunar::{
    calc_lunar_point_in_moment, JD_INTERP_LILITH_SAFE_MAX, JD_INTERP_LILITH_SAFE_MIN,
    MEAN_MOON_DIST_AU, MU_EM_KM3_S2,
};
pub use manifest::{DataManifest, TimeSyncManifest};
pub use moment::{EpochContext, Moment, H_MOON_DAYS, H_PLANET_DAYS};
pub use moshier::{MoshierAnalyticBackend, JD_MOSHIER_MAX, JD_MOSHIER_MIN};
pub use nodes_apsides::{
    calc_nodes_apsides, NodesApsidesResult, SE_NODBIT_FOPOINT, SE_NODBIT_MEAN, SE_NODBIT_OSCU,
    SE_NODBIT_OSCU_BAR,
};
pub use occultations::{
    local_lunar_occultations_in_range, lunar_occultations_in_range, next_local_lunar_occultation,
    next_lunar_occultation, LunarOccultationInfo, OccultationKind,
};
pub use phenomena::{
    body_equatorial_radius_km, calc_phenomena, calc_phenomena_in_moment,
    compute_apparent_magnitude_v, PhenomenaResult, MINOR_BODIES_PHENO_68,
};
pub use spk::{AniseSpkBackend, SpkSource, SpkState, TargetCoverage};
pub use stars::{
    calc_star_entry_in_moment, calc_star_in_moment, eval_star_unit_vector, find_star_by_hip,
    find_star_by_name_or_bayer, list_stars, load_custom_stars_file, parse_sefstars_catalog,
    star_catalog_52, star_info, star_mean_ecliptic_lon_no_ab_at_ctx, StarEntry, StarMetadata,
    StarQualityWarning, DAS2R, JD_HIP_J1991_25, PARSEC_AU,
};
#[cfg(feature = "swe-compat")]
pub use swe_compat::{
    body_id_from_swe_ipl, calc_options_from_swe_flags, format_swe_state6, SweAdapter,
    SEFLG_BARYCTR, SEFLG_EQUATORIAL, SEFLG_HELCTR, SEFLG_ICRS, SEFLG_J2000, SEFLG_JPLEPH,
    SEFLG_MOSEPH, SEFLG_NOABERR, SEFLG_NOGDEFL, SEFLG_NONUT, SEFLG_RADIANS, SEFLG_SIDEREAL,
    SEFLG_SPEED, SEFLG_SPEED3, SEFLG_SWIEPH, SEFLG_TOPOCTR, SEFLG_TRUEPOS, SEFLG_XYZ, SE_ADMETOS,
    SE_APOLLON, SE_AST_OFFSET, SE_BIT_ASTRO_TWILIGHT, SE_BIT_CIVIL_TWILIGHT, SE_BIT_DISC_BOTTOM,
    SE_BIT_DISC_CENTER, SE_BIT_FIXED_DISC_SIZE, SE_BIT_GEOCTR_NO_ECL_LAT, SE_BIT_HINDU_RISING,
    SE_BIT_NAUTIC_TWILIGHT, SE_BIT_NO_REFRACTION, SE_CALC_ITRANSIT, SE_CALC_MTRANSIT, SE_CALC_RISE,
    SE_CALC_SET, SE_CERES, SE_CHIRON, SE_CUPIDO, SE_EARTH, SE_ECL2HOR, SE_ECL_1ST_VISIBLE,
    SE_ECL_2ND_VISIBLE, SE_ECL_3RD_VISIBLE, SE_ECL_4TH_VISIBLE, SE_ECL_ALLTYPES_LUNAR,
    SE_ECL_ALLTYPES_SOLAR, SE_ECL_ANNULAR, SE_ECL_ANNULAR_TOTAL, SE_ECL_CENTRAL, SE_ECL_HYBRID,
    SE_ECL_MAX_VISIBLE, SE_ECL_NONCENTRAL, SE_ECL_NUT, SE_ECL_OCC_BEG_DAYLIGHT,
    SE_ECL_OCC_END_DAYLIGHT, SE_ECL_ONE_TRY, SE_ECL_PARTIAL, SE_ECL_PENUMBRAL, SE_ECL_TOTAL,
    SE_ECL_VISIBLE, SE_EQU2HOR, SE_EVENING_FIRST, SE_FICT_OFFSET, SE_GREG_CAL, SE_HADES,
    SE_HELFLAG_HIGH_PRECISION, SE_HELIACAL_RISING, SE_HELIACAL_SETTING, SE_HOR2ECL, SE_HOR2EQU,
    SE_INTP_APOG, SE_INTP_PERG, SE_ISIS, SE_JUL_CAL, SE_JUNO, SE_JUPITER, SE_KRONOS, SE_MARS,
    SE_MEAN_APOG, SE_MEAN_NODE, SE_MERCURY, SE_MOON, SE_MORNING_LAST, SE_NEPTUNE, SE_OSCU_APOG,
    SE_PALLAS, SE_PHOLUS, SE_PLUTO, SE_POSEIDON, SE_SATURN, SE_SIDBIT_ECL_DATE, SE_SIDBIT_ECL_T0,
    SE_SIDBIT_NO_PREC_OFFSET, SE_SIDBIT_PREC_ORIG, SE_SIDBIT_SSY_PLANE, SE_SIDBIT_USER_UT,
    SE_SIDM_ALDEBARAN_15TAU, SE_SIDM_ARYABHATA, SE_SIDM_ARYABHATA_522, SE_SIDM_ARYABHATA_MSUN,
    SE_SIDM_B1950, SE_SIDM_BABYL_BRITTON, SE_SIDM_BABYL_ETPSC, SE_SIDM_BABYL_HUBER,
    SE_SIDM_BABYL_KUGLER1, SE_SIDM_BABYL_KUGLER2, SE_SIDM_BABYL_KUGLER3, SE_SIDM_DELUCE,
    SE_SIDM_DJWHAL_KHUL, SE_SIDM_FAGAN_BRADLEY, SE_SIDM_GALALIGN_MARDYKS, SE_SIDM_GALCENT_0SAG,
    SE_SIDM_GALCENT_COCHRANE, SE_SIDM_GALCENT_MULA_WILHELM, SE_SIDM_GALCENT_RGILBRAND,
    SE_SIDM_GALEQU_FIORENZA, SE_SIDM_GALEQ_IAU1958, SE_SIDM_GALEQ_MULA, SE_SIDM_GALEQ_TRUE,
    SE_SIDM_HIPPARCHOS, SE_SIDM_J1900, SE_SIDM_J2000, SE_SIDM_JN_BHASIN, SE_SIDM_KP_SENTHILATHIBAN,
    SE_SIDM_KRISHNAMURTI, SE_SIDM_KRISHNAMURTI_VP291, SE_SIDM_LAHIRI, SE_SIDM_LAHIRI_1940,
    SE_SIDM_LAHIRI_ICRC, SE_SIDM_LAHIRI_VP285, SE_SIDM_RAMAN, SE_SIDM_SASSANIAN, SE_SIDM_SS_CITRA,
    SE_SIDM_SS_REVATI, SE_SIDM_SURYASIDDHANTA, SE_SIDM_SURYASIDDHANTA_MSUN,
    SE_SIDM_TRUE_ALDEBARAN, SE_SIDM_TRUE_CITRA, SE_SIDM_TRUE_MULA, SE_SIDM_TRUE_PUSHYA,
    SE_SIDM_TRUE_REVATI, SE_SIDM_TRUE_SHERATAN, SE_SIDM_USER, SE_SIDM_USHASHASHI,
    SE_SIDM_VALENS_MOON, SE_SIDM_VEDIC_SHEORAN, SE_SIDM_YUKTESHWAR, SE_SUN, SE_TRUE_NODE,
    SE_URANUS, SE_VENUS, SE_VESTA, SE_VULKANUS, SE_ZEUS,
};
pub use time::{
    earth_rotation_angle_rad, gast_rad, gmst_rad, julday, julday_split, revjul, tdb_minus_tt_sec,
    tdb_to_tt, tt_to_tdb, CalendarDate, CalendarKind, DeltaTInfo, SiderealTimeInfo, SplitJd, Tdb,
    TimeInput, TimeTable, Tt, Ut1, Utc, EMBEDDED_TIME_TABLE_VERSION,
};
pub use types::{
    AbMode, AyanamsaModel, BodyId, CalcOptions, Equinox, HipId, HouseSystem, HouseWarning,
    HypotheticalId, LunarPointId, MinorBodyCategory, MinorBodySpec, MpcNumber, NutationMode,
    Observer, Origin, Plane, PlanetId, PolarCategory, PolarPolicy, TargetKind, Zodiac, AU_KM,
    JD_1850_01_01, JD_2100_01_01, JD_J2000, MINOR_BODIES_68, SECONDS_PER_DAY,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minor_bodies_registry_has_68_unique_sorted_entries() {
        assert_eq!(MINOR_BODIES_68.len(), 68);
        for pair in MINOR_BODIES_68.windows(2) {
            if let [a, b] = pair {
                assert!(
                    a.mpc < b.mpc,
                    "MPC numbers must be strictly sorted: {:?} vs {:?}",
                    a,
                    b
                );
            }
        }
    }

    #[test]
    fn test_house_systems_23_case_sensitive() {
        assert_eq!(HouseSystem::ALL.len(), 23);
        for hs in HouseSystem::ALL {
            let code = hs.code();
            let parsed = HouseSystem::from_code(code);
            assert_eq!(parsed, Ok(hs));
        }
        // Distinguish 'I' (Sunshine Treindl) and 'i' (Sunshine Makransky alt)
        assert_eq!(
            HouseSystem::from_code('I'),
            Ok(HouseSystem::SunshineTreindl)
        );
        assert_eq!(
            HouseSystem::from_code('i'),
            Ok(HouseSystem::SunshineMakransky)
        );
        // 'G' (Gauquelin) is not one of the 23 mandatory house system codes
        assert!(HouseSystem::from_code('G').is_err());
    }

    #[test]
    fn test_ayanamsa_models_12_strict_lookup() {
        assert_eq!(AyanamsaModel::ALL.len(), 12);
        for model in AyanamsaModel::ALL {
            assert_eq!(AyanamsaModel::from_code(model.code()), Ok(model));
        }
        // Unknown model must error, never fallback to Lahiri
        assert_eq!(
            AyanamsaModel::from_code("unknown_ayanamsa"),
            Err(EphemerisError::UnsupportedModel {
                model: "unknown_ayanamsa".to_owned()
            })
        );
    }

    #[test]
    fn test_planet_naif_mapping() {
        assert_eq!(PlanetId::Sun.naif_target_id(), 10);
        assert_eq!(PlanetId::Moon.naif_target_id(), 301);
        assert_eq!(PlanetId::Mercury.naif_target_id(), 199);
        assert_eq!(PlanetId::Venus.naif_target_id(), 299);
        assert_eq!(PlanetId::Earth.naif_target_id(), 399);
        assert_eq!(PlanetId::Mars.naif_target_id(), 4);
        assert_eq!(PlanetId::Jupiter.naif_target_id(), 5);
        assert_eq!(PlanetId::Saturn.naif_target_id(), 6);
        assert_eq!(PlanetId::Uranus.naif_target_id(), 7);
        assert_eq!(PlanetId::Neptune.naif_target_id(), 8);
        assert_eq!(PlanetId::Pluto.naif_target_id(), 9);
        assert_eq!(PlanetId::Mercury.target_kind(), TargetKind::BodyCenter);
        assert_eq!(PlanetId::Pluto.target_kind(), TargetKind::SystemBarycenter);
    }

    #[test]
    fn test_polar_policy_default_is_strict() {
        assert_eq!(PolarPolicy::default(), PolarPolicy::Strict);
    }

    #[test]
    fn test_corrupt_kernel_rejected_without_panic() {
        let bad_bytes = vec![0u8; 512];
        let res = AniseSpkBackend::from_named_bytes(&[("corrupt.bsp".to_owned(), bad_bytes)]);
        assert!(matches!(res, Err(EphemerisError::CorruptKernel { .. })));
    }

    #[test]
    fn test_de440s_and_asteroids_68_bundle_integration() {
        let workspace_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(|p| p.parent())
            .unwrap_or_else(|| std::path::Path::new("."));
        let de440s_path = workspace_dir.join("data/de440s.bsp");
        let ast_path = workspace_dir.join("data/asteroids_68_1850_2100.bsp");
        let manifest_path = workspace_dir.join("data/manifest.json");

        if !de440s_path.exists() || !ast_path.exists() || !manifest_path.exists() {
            return;
        }

        let manifest_str = std::fs::read_to_string(&manifest_path).expect("read manifest.json");
        let manifest: DataManifest =
            serde_json::from_str(&manifest_str).expect("parse manifest.json");
        assert_eq!(manifest.validate(), Ok(()));
        assert_eq!(manifest.minor_bodies_kernel.body_count, 68);

        let backend =
            AniseSpkBackend::from_bsp_paths(&[de440s_path, ast_path]).expect("load kernels");
        assert_eq!(backend.segments().len(), 14 + 68);

        // Verify all major planets at 1850, J2000, and 2100
        for p in [
            PlanetId::Sun,
            PlanetId::Moon,
            PlanetId::Earth,
            PlanetId::Mercury,
            PlanetId::Venus,
            PlanetId::Mars,
            PlanetId::Jupiter,
            PlanetId::Saturn,
            PlanetId::Uranus,
            PlanetId::Neptune,
            PlanetId::Pluto,
        ] {
            for jd in [JD_1850_01_01, JD_J2000, JD_2100_01_01] {
                let st = backend
                    .state_icrf_km(p.naif_target_id(), 0, jd)
                    .expect("planet state");
                assert!(st.pos_km[0].is_finite() && st.vel_km_s[0].is_finite());
            }
        }

        // Verify all 68 minor bodies at 1850, J2000, and 2100, and check OutOfCoverage outside bounds
        for entry in &manifest.minor_bodies_kernel.bodies {
            for jd in [JD_1850_01_01, JD_J2000, JD_2100_01_01] {
                let st = backend
                    .state_icrf_km(entry.naif_spk_id, 0, jd)
                    .expect("minor body state");
                assert!(st.pos_km[0].is_finite() && st.vel_km_s[0].is_finite());
            }
            let out_err = backend.state_icrf_km(entry.naif_spk_id, 0, 2_300_000.0);
            assert!(matches!(out_err, Err(EphemerisError::OutOfCoverage { .. })));
        }
    }
}
