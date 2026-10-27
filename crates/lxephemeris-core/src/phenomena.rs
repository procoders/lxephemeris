//! Planetary, lunar, and minor body phenomena (`swe_pheno` / `swe_pheno_ut` compatibility):
//! phase angle $\alpha$, illuminated fraction $k$, elongation $\psi$, apparent equatorial
//! disc diameter $d$, apparent visual magnitude $V$ (Mallama & Hilton 2018 / Hilton 2005 / Bowell 1989),
//! and equatorial horizontal parallax $\pi$.
//!
//! References:
//! - Mallama, A., & Hilton, J. L. (2018), "Computing apparent planetary magnitudes for The Astronomical Almanac",
//!   *Astronomy and Computing*, 25, 10–24.
//! - Hilton, J. L. (2005), "Improving the Photometric Magnitude Polynomials for the Major and Minor Planets",
//!   *Astronomical Journal*, 130, 2928–2935.
//! - Bowell, E., et al. (1989), "Application of Photometric Models to Asteroids", in *Asteroids II*, 524–556.
//! - Archinal, B. A., et al. (2018), "Report of the IAU Working Group on Cartographic Coordinates and
//!   Rotational Elements: 2015", *Celestial Mechanics and Dynamical Astronomy*, 130(22).

use crate::astrometry::calc_body_in_moment;
use crate::error::EphemerisError;
use crate::frame::{vec_dot, vec_norm};
use crate::moment::Moment;
use crate::spk::SpkSource;
use crate::stars::find_star_by_hip;
use crate::time::{TimeInput, TimeTable};
use crate::types::{
    AbMode, BodyId, CalcOptions, Equinox, MpcNumber, NutationMode, Origin, Plane, PlanetId,
    Zodiac, AU_KM,
};
use serde::{Deserialize, Serialize};

/// Earth equatorial radius in kilometers (IAU 2015 nominal / WGS84 $a$).
pub const R_EARTH_EQUATORIAL_KM: f64 = 6_378.136_6;

/// Visual and physical parameters for minor planets ($H$, $G$, nominal radius $R_{\text{km}}$).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MinorBodyPhenoData {
    /// Minor Planet Center (MPC) number.
    pub mpc: u32,
    /// Absolute visual magnitude $H$ in the IAU $(H, G)$ system.
    pub h_mag: f64,
    /// Slope parameter $G$ in the IAU $(H, G)$ system.
    pub g_slope: f64,
    /// Nominal physical equatorial/mean radius in kilometers.
    pub radius_km: f64,
}

/// Registry of visual and geometric properties for all 68 minor bodies in `MINOR_BODIES_68`.
///
/// Sorted by MPC number for $O(\log N)$ binary searching.
pub const MINOR_BODIES_PHENO_68: [MinorBodyPhenoData; 68] = [
    MinorBodyPhenoData { mpc: 1, h_mag: 3.34, g_slope: 0.12, radius_km: 469.73 },
    MinorBodyPhenoData { mpc: 2, h_mag: 4.13, g_slope: 0.11, radius_km: 256.00 },
    MinorBodyPhenoData { mpc: 3, h_mag: 5.33, g_slope: 0.32, radius_km: 127.00 },
    MinorBodyPhenoData { mpc: 4, h_mag: 3.20, g_slope: 0.32, radius_km: 262.70 },
    MinorBodyPhenoData { mpc: 5, h_mag: 6.95, g_slope: 0.15, radius_km: 53.35 },
    MinorBodyPhenoData { mpc: 6, h_mag: 5.62, g_slope: 0.24, radius_km: 92.59 },
    MinorBodyPhenoData { mpc: 7, h_mag: 5.70, g_slope: 0.15, radius_km: 99.92 },
    MinorBodyPhenoData { mpc: 8, h_mag: 6.62, g_slope: 0.28, radius_km: 73.75 },
    MinorBodyPhenoData { mpc: 9, h_mag: 6.18, g_slope: 0.17, radius_km: 95.00 },
    MinorBodyPhenoData { mpc: 10, h_mag: 5.43, g_slope: 0.15, radius_km: 217.00 },
    MinorBodyPhenoData { mpc: 11, h_mag: 6.73, g_slope: 0.15, radius_km: 71.44 },
    MinorBodyPhenoData { mpc: 12, h_mag: 7.30, g_slope: 0.22, radius_km: 57.54 },
    MinorBodyPhenoData { mpc: 13, h_mag: 6.92, g_slope: 0.15, radius_km: 101.32 },
    MinorBodyPhenoData { mpc: 14, h_mag: 6.55, g_slope: 0.15, radius_km: 76.00 },
    MinorBodyPhenoData { mpc: 15, h_mag: 5.42, g_slope: 0.23, radius_km: 115.84 },
    MinorBodyPhenoData { mpc: 16, h_mag: 5.90, g_slope: 0.20, radius_km: 113.00 },
    MinorBodyPhenoData { mpc: 18, h_mag: 6.35, g_slope: 0.25, radius_km: 69.80 },
    MinorBodyPhenoData { mpc: 19, h_mag: 7.49, g_slope: 0.10, radius_km: 100.00 },
    MinorBodyPhenoData { mpc: 26, h_mag: 7.48, g_slope: 0.15, radius_km: 47.40 },
    MinorBodyPhenoData { mpc: 30, h_mag: 7.57, g_slope: 0.15, radius_km: 46.39 },
    MinorBodyPhenoData { mpc: 37, h_mag: 7.29, g_slope: 0.24, radius_km: 54.17 },
    MinorBodyPhenoData { mpc: 40, h_mag: 6.57, g_slope: 0.15, radius_km: 55.63 },
    MinorBodyPhenoData { mpc: 42, h_mag: 7.69, g_slope: 0.15, radius_km: 55.50 },
    MinorBodyPhenoData { mpc: 55, h_mag: 7.82, g_slope: 0.15, radius_km: 42.40 },
    MinorBodyPhenoData { mpc: 78, h_mag: 8.39, g_slope: 0.08, radius_km: 60.30 },
    MinorBodyPhenoData { mpc: 80, h_mag: 8.06, g_slope: 0.15, radius_km: 34.28 },
    MinorBodyPhenoData { mpc: 100, h_mag: 7.73, g_slope: 0.15, radius_km: 42.87 },
    MinorBodyPhenoData { mpc: 105, h_mag: 8.75, g_slope: 0.10, radius_km: 47.43 },
    MinorBodyPhenoData { mpc: 114, h_mag: 8.45, g_slope: 0.15, radius_km: 47.09 },
    MinorBodyPhenoData { mpc: 128, h_mag: 7.68, g_slope: 0.15, radius_km: 81.26 },
    MinorBodyPhenoData { mpc: 227, h_mag: 9.13, g_slope: 0.15, radius_km: 61.95 },
    MinorBodyPhenoData { mpc: 251, h_mag: 9.95, g_slope: 0.15, radius_km: 13.75 },
    MinorBodyPhenoData { mpc: 259, h_mag: 7.80, g_slope: 0.15, radius_km: 87.16 },
    MinorBodyPhenoData { mpc: 344, h_mag: 8.28, g_slope: 0.15, radius_km: 62.09 },
    MinorBodyPhenoData { mpc: 375, h_mag: 7.48, g_slope: 0.27, radius_km: 108.00 },
    MinorBodyPhenoData { mpc: 399, h_mag: 8.91, g_slope: 0.15, radius_km: 19.78 },
    MinorBodyPhenoData { mpc: 407, h_mag: 9.16, g_slope: 0.15, radius_km: 47.53 },
    MinorBodyPhenoData { mpc: 433, h_mag: 10.40, g_slope: 0.46, radius_km: 8.42 },
    MinorBodyPhenoData { mpc: 447, h_mag: 9.05, g_slope: 0.15, radius_km: 42.61 },
    MinorBodyPhenoData { mpc: 471, h_mag: 6.34, g_slope: 0.37, radius_km: 74.06 },
    MinorBodyPhenoData { mpc: 490, h_mag: 8.66, g_slope: 0.15, radius_km: 59.40 },
    MinorBodyPhenoData { mpc: 504, h_mag: 10.12, g_slope: 0.15, radius_km: 15.22 },
    MinorBodyPhenoData { mpc: 763, h_mag: 12.67, g_slope: 0.15, radius_km: 3.50 },
    MinorBodyPhenoData { mpc: 944, h_mag: 10.55, g_slope: 0.15, radius_km: 19.00 },
    MinorBodyPhenoData { mpc: 1181, h_mag: 11.36, g_slope: 0.15, radius_km: 10.25 },
    MinorBodyPhenoData { mpc: 1221, h_mag: 17.37, g_slope: 0.15, radius_km: 0.50 },
    MinorBodyPhenoData { mpc: 1387, h_mag: 13.27, g_slope: 0.15, radius_km: 3.73 },
    MinorBodyPhenoData { mpc: 1388, h_mag: 11.14, g_slope: 0.15, radius_km: 10.68 },
    MinorBodyPhenoData { mpc: 1489, h_mag: 11.74, g_slope: 0.15, radius_km: 13.39 },
    MinorBodyPhenoData { mpc: 1566, h_mag: 16.53, g_slope: 0.15, radius_km: 0.50 },
    MinorBodyPhenoData { mpc: 1580, h_mag: 14.54, g_slope: 0.00, radius_km: 2.90 },
    MinorBodyPhenoData { mpc: 1585, h_mag: 10.30, g_slope: 0.15, radius_km: 27.64 },
    MinorBodyPhenoData { mpc: 1862, h_mag: 16.09, g_slope: 0.09, radius_km: 0.75 },
    MinorBodyPhenoData { mpc: 1930, h_mag: 11.11, g_slope: 0.15, radius_km: 17.22 },
    MinorBodyPhenoData { mpc: 2060, h_mag: 5.80, g_slope: 0.15, radius_km: 103.00 },
    MinorBodyPhenoData { mpc: 2101, h_mag: 18.69, g_slope: 0.15, radius_km: 0.30 },
    MinorBodyPhenoData { mpc: 3671, h_mag: 16.50, g_slope: 0.15, radius_km: 0.75 },
    MinorBodyPhenoData { mpc: 3782, h_mag: 12.80, g_slope: 0.15, radius_km: 2.96 },
    MinorBodyPhenoData { mpc: 3811, h_mag: 11.95, g_slope: 0.15, radius_km: 15.17 },
    MinorBodyPhenoData { mpc: 4386, h_mag: 12.93, g_slope: 0.15, radius_km: 7.83 },
    MinorBodyPhenoData { mpc: 5129, h_mag: 12.98, g_slope: 0.15, radius_km: 3.74 },
    MinorBodyPhenoData { mpc: 5145, h_mag: 7.00, g_slope: 0.15, radius_km: 95.00 },
    MinorBodyPhenoData { mpc: 6583, h_mag: 12.98, g_slope: 0.15, radius_km: 5.60 },
    MinorBodyPhenoData { mpc: 7066, h_mag: 9.57, g_slope: 0.15, radius_km: 25.62 },
    MinorBodyPhenoData { mpc: 90377, h_mag: 1.50, g_slope: 0.15, radius_km: 500.00 },
    MinorBodyPhenoData { mpc: 136108, h_mag: 0.20, g_slope: 0.15, radius_km: 816.00 },
    MinorBodyPhenoData { mpc: 136199, h_mag: -1.19, g_slope: 0.15, radius_km: 1163.00 },
    MinorBodyPhenoData { mpc: 136472, h_mag: -0.20, g_slope: 0.15, radius_km: 715.00 },
];

/// Looks up the `(H, G, radius_km)` physical parameters for a minor planet.
///
/// If `mpc` is in `MINOR_BODIES_PHENO_68`, the catalog parameters are returned.
/// For any other numbered minor body, standard default values are calculated:
/// $H = 6.50$, $G = 0.15$, and $R_{\text{km}} = \frac{664.5}{\sqrt{p_V}} 10^{-0.2 H}$ with $p_V = 0.10$.
#[must_use]
pub fn minor_body_pheno_lookup(mpc: MpcNumber) -> (f64, f64, f64) {
    if let Ok(idx) = MINOR_BODIES_PHENO_68.binary_search_by_key(&mpc.0, |d| d.mpc) {
        if let Some(entry) = MINOR_BODIES_PHENO_68.get(idx) {
            return (entry.h_mag, entry.g_slope, entry.radius_km);
        }
    }
    let default_h = 6.50;
    let default_g = 0.15;
    let default_pv: f64 = 0.10;
    let default_r_km = (664.5 / default_pv.sqrt()) * 10.0_f64.powf(-0.2 * default_h);
    (default_h, default_g, default_r_km)
}

/// Planetary and lunar phenomena result (`swe_pheno` / `swe_pheno_ut` 20-attribute compatible).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhenomenaResult {
    /// Phase angle $\alpha$ in degrees $\in [0^\circ, 180^\circ]$ (angle Sun–Target–Observer).
    pub phase_angle_deg: f64,
    /// Phase illuminated fraction $k = (1 + \cos\alpha)/2 \in [0, 1]$.
    pub phase_illuminated: f64,
    /// Elongation $\psi$ from the Sun in degrees $\in [0^\circ, 180^\circ]$ (angle Target–Observer–Sun).
    pub elongation_deg: f64,
    /// Apparent equatorial disc diameter $d = 2\arcsin(R_{\text{eq}} / \Delta_{\text{km}})$ in degrees.
    pub apparent_diameter_deg: f64,
    /// Apparent visual magnitude $V$.
    pub apparent_magnitude_v: f64,
    /// Equatorial horizontal parallax $\pi = \arcsin(R_\oplus / \Delta_{\text{km}})$ in degrees.
    pub horizontal_parallax_deg: f64,
}

impl PhenomenaResult {
    /// Formats the result as a 20-element `f64` array compatible with Swiss Ephemeris `swe_pheno`.
    ///
    /// The indices are mapped as:
    /// - `attr[0]`: `phase_angle_deg` ($\alpha^\circ$)
    /// - `attr[1]`: `phase_illuminated` ($k \in [0, 1]$)
    /// - `attr[2]`: `elongation_deg` ($\psi^\circ$)
    /// - `attr[3]`: `apparent_diameter_deg` ($d^\circ$)
    /// - `attr[4]`: `apparent_magnitude_v` ($V$)
    /// - `attr[5]`: `horizontal_parallax_deg` ($\pi^\circ$)
    /// - `attr[6..20]`: `0.0`
    #[must_use]
    pub fn to_swe_attr20(&self) -> [f64; 20] {
        let mut attr = [0.0; 20];
        attr[0] = self.phase_angle_deg;
        attr[1] = self.phase_illuminated;
        attr[2] = self.elongation_deg;
        attr[3] = self.apparent_diameter_deg;
        attr[4] = self.apparent_magnitude_v;
        attr[5] = self.horizontal_parallax_deg;
        attr
    }

    /// Apparent equatorial disc diameter in arcseconds (`apparent_diameter_deg * 3600.0`).
    #[must_use]
    pub fn apparent_diameter_arcsec(&self) -> f64 {
        self.apparent_diameter_deg * 3_600.0
    }

    /// Equatorial horizontal parallax in arcseconds (`horizontal_parallax_deg * 3600.0`).
    #[must_use]
    pub fn horizontal_parallax_arcsec(&self) -> f64 {
        self.horizontal_parallax_deg * 3_600.0
    }
}

/// Returns the IAU 2015 / WGCCRE (Archinal et al. 2018) nominal equatorial radius in kilometers.
///
/// Returns `None` for abstract mathematical points, reference barycenters, fixed stars,
/// or undefined bodies.
#[must_use]
pub fn body_equatorial_radius_km(body: BodyId) -> Option<f64> {
    match body {
        BodyId::Planet(p) => match p {
            PlanetId::Sun => Some(695_700.0),
            PlanetId::Moon => Some(1_737.40),
            PlanetId::Mercury => Some(2_439.70),
            PlanetId::Venus => Some(6_051.80),
            PlanetId::Earth => Some(R_EARTH_EQUATORIAL_KM),
            PlanetId::Mars => Some(3_396.19),
            PlanetId::Jupiter => Some(71_492.0),
            PlanetId::Saturn => Some(60_268.0),
            PlanetId::Uranus => Some(25_559.0),
            PlanetId::Neptune => Some(24_764.0),
            PlanetId::Pluto => Some(1_188.30),
            PlanetId::EarthMoonBarycenter | PlanetId::SolarSystemBarycenter => None,
        },
        BodyId::MinorPlanet(mpc) => {
            let (_, _, r_km) = minor_body_pheno_lookup(mpc);
            Some(r_km)
        }
        BodyId::FixedStar(_)
        | BodyId::LunarPoint(_)
        | BodyId::Hypothetical(_)
        | BodyId::EclNut => None,
    }
}

/// Evaluates the apparent visual magnitude $V$ using the IAU / Mallama & Hilton (2018)
/// and Hilton (2005) planetary phase curves, or the Bowell et al. (1989) IAU $(H, G)$ model.
///
/// # Arguments
/// - `body`: Target celestial object.
/// - `r_helio_au`: 3D Cartesian vector from the Sun to the target in AU.
/// - `delta_obs_au`: 3D Cartesian vector from the observer to the target in AU.
/// - `r_obs_sun_au`: 3D Cartesian vector from the observer to the Sun in AU.
#[must_use]
pub fn compute_apparent_magnitude_v(
    body: BodyId,
    r_helio_au: &[f64; 3],
    delta_obs_au: &[f64; 3],
    r_obs_sun_au: &[f64; 3],
) -> Option<f64> {
    let r_au = vec_norm(r_helio_au);
    let delta_au = vec_norm(delta_obs_au);
    let r_sun_au = vec_norm(r_obs_sun_au);

    // Phase angle alpha: angle Sun–Target–Observer
    // Vector Target -> Sun: -r_helio_au
    // Vector Target -> Observer: -delta_obs_au
    // Dot product: (-r_helio_au) . (-delta_obs_au) = r_helio_au . delta_obs_au
    let dot_phase = vec_dot(r_helio_au, delta_obs_au);
    let cos_alpha = if r_au > 1e-15 && delta_au > 1e-15 {
        (dot_phase / (r_au * delta_au)).clamp(-1.0, 1.0)
    } else {
        1.0
    };
    let alpha_deg = cos_alpha.acos().to_degrees();

    match body {
        BodyId::Planet(p) => match p {
            PlanetId::Sun => {
                let dist = if delta_au > 1e-12 {
                    delta_au
                } else if r_sun_au > 1e-12 {
                    r_sun_au
                } else {
                    1.0
                };
                Some(-26.74 + 5.0 * dist.max(1e-12).log10())
            }
            PlanetId::Moon => {
                let r_odot = if r_sun_au > 1e-12 { r_sun_au } else { 1.0 };
                let mean_dist_au = 384_400.0 / AU_KM;
                let log_term = 5.0 * (r_odot * (delta_au / mean_dist_au)).max(1e-12).log10();
                let mag = -12.74 + log_term + 0.026 * alpha_deg + 4.0e-9 * alpha_deg.powi(4);
                Some(mag)
            }
            PlanetId::Mercury => {
                let a = alpha_deg;
                let poly = 6.328e-2 * a
                    - 1.6336e-3 * a.powi(2)
                    + 3.3644e-5 * a.powi(3)
                    - 3.4265e-7 * a.powi(4)
                    + 1.6893e-9 * a.powi(5)
                    - 3.0334e-12 * a.powi(6);
                let mag = -0.613 + 5.0 * (r_au * delta_au).max(1e-12).log10() + poly;
                Some(mag)
            }
            PlanetId::Venus => {
                let a = alpha_deg;
                let delta_v = if a <= 163.7 {
                    -1.044e-3 * a + 3.687e-4 * a.powi(2) - 2.814e-6 * a.powi(3) + 8.938e-9 * a.powi(4)
                } else {
                    // Mallama & Hilton 2018 extreme-phase forward-scattering equation
                    // The constant 236.05828 replaces the baseline -4.384 term.
                    return Some(
                        5.0 * (r_au * delta_au).max(1e-12).log10()
                            + 236.05828
                            - 2.81914 * a
                            + 8.39034e-3 * a.powi(2),
                    );
                };
                let mag = -4.384 + 5.0 * (r_au * delta_au).max(1e-12).log10() + delta_v;
                Some(mag)
            }
            PlanetId::Earth => {
                let a = alpha_deg;
                let delta_v = -1.060e-3 * a + 2.054e-4 * a.powi(2);
                let mag = -3.99 + 5.0 * (r_au * delta_au).max(1e-12).log10() + delta_v;
                Some(mag)
            }
            PlanetId::Mars => {
                let a = alpha_deg;
                let poly = 0.02267 * a - 0.0001302 * a.powi(2);
                let mag = -1.601 + 5.0 * (r_au * delta_au).max(1e-12).log10() + poly;
                Some(mag)
            }
            PlanetId::Jupiter => {
                let a = alpha_deg;
                let poly = 0.03788 * a - 0.001325 * a.powi(2);
                let mag = -9.395 + 5.0 * (r_au * delta_au).max(1e-12).log10() + poly;
                Some(mag)
            }
            PlanetId::Saturn => {
                let a = alpha_deg;
                // Saturn north pole in ICRF (Archinal et al. 2018): alpha_0 = 40.589 deg, delta_0 = 83.537 deg
                let alpha_0 = 40.589_f64.to_radians();
                let delta_0 = 83.537_f64.to_radians();
                let p_pole = [
                    delta_0.cos() * alpha_0.cos(),
                    delta_0.cos() * alpha_0.sin(),
                    delta_0.sin(),
                ];
                let delta_hat = if delta_au > 1e-15 {
                    [
                        delta_obs_au[0] / delta_au,
                        delta_obs_au[1] / delta_au,
                        delta_obs_au[2] / delta_au,
                    ]
                } else {
                    [0.0, 0.0, 1.0]
                };
                // Ring tilt B satisfies sin(B) = -p_pole . delta_hat
                let sin_b = -(p_pole[0] * delta_hat[0] + p_pole[1] * delta_hat[1] + p_pole[2] * delta_hat[2]);
                let sin_b_abs = sin_b.abs();
                let mag = -8.914
                    + 5.0 * (r_au * delta_au).max(1e-12).log10()
                    - 2.60 * sin_b_abs
                    + 1.25 * sin_b.powi(2)
                    + 0.044 * a;
                Some(mag)
            }
            PlanetId::Uranus => {
                let a = alpha_deg;
                let mag = -7.110 + 5.0 * (r_au * delta_au).max(1e-12).log10() + 0.0028 * a;
                Some(mag)
            }
            PlanetId::Neptune => {
                let a = alpha_deg;
                let mag = -6.885 + 5.0 * (r_au * delta_au).max(1e-12).log10() + 0.0079 * a;
                Some(mag)
            }
            PlanetId::Pluto => {
                let a = alpha_deg;
                let mag = -1.00 + 5.0 * (r_au * delta_au).max(1e-12).log10() + 0.041 * a;
                Some(mag)
            }
            PlanetId::EarthMoonBarycenter | PlanetId::SolarSystemBarycenter => None,
        },
        BodyId::MinorPlanet(mpc) => {
            let (h_abs, g_slope, _) = minor_body_pheno_lookup(mpc);
            let alpha_rad = alpha_deg.to_radians();
            let tan_half = (0.5 * alpha_rad).tan().max(0.0);
            // Bowell et al. (1989) IAU (H, G) phase functions Phi_1 and Phi_2
            let phi1 = (-3.332 * tan_half.powf(0.631)).exp();
            let phi2 = (-1.862 * tan_half.powf(1.218)).exp();
            let phase_flux = ((1.0 - g_slope) * phi1 + g_slope * phi2).max(1e-6);
            let mag = h_abs + 5.0 * (r_au * delta_au).max(1e-12).log10() - 2.5 * phase_flux.log10();
            Some(mag)
        }
        BodyId::FixedStar(hip) => find_star_by_hip(hip).map(|s| s.vmag),
        BodyId::LunarPoint(_) | BodyId::Hypothetical(_) | BodyId::EclNut => None,
    }
}

/// Evaluates fixed-star phenomena for an explicit [`crate::stars::StarEntry`] (including custom stars)
/// inside an existing [`Moment`].
pub fn calc_star_entry_phenomena_in_moment(
    moment: &Moment<'_>,
    star: &crate::stars::StarEntry,
    options: &CalcOptions,
) -> Result<PhenomenaResult, EphemerisError> {
    let mut icrf_opts = *options;
    icrf_opts.plane = Plane::Equatorial;
    icrf_opts.equinox = Equinox::J2000;
    icrf_opts.nutation = NutationMode::Mean;
    icrf_opts.zodiac = Zodiac::Tropical;
    icrf_opts.include_speed = false;

    let sun_calc = calc_body_in_moment(moment, BodyId::Planet(PlanetId::Sun), &icrf_opts)?;
    let r_obs_sun_au = sun_calc.cartesian_au;
    let norm_sun = sun_calc.dist_au;

    let star_calc = crate::stars::calc_star_entry_in_moment(moment, star, &icrf_opts, None)?;
    let dot_elong = vec_dot(&star_calc.cartesian_au, &r_obs_sun_au);
    let norm_star = vec_norm(&star_calc.cartesian_au);
    let cos_elong = if norm_star > 0.0 && norm_sun > 0.0 {
        (dot_elong / (norm_star * norm_sun)).clamp(-1.0, 1.0)
    } else {
        1.0
    };
    let elongation_deg = cos_elong.acos().to_degrees();
    Ok(PhenomenaResult {
        phase_angle_deg: 0.0,
        phase_illuminated: 1.0,
        elongation_deg,
        apparent_diameter_deg: 0.0,
        apparent_magnitude_v: star.vmag,
        horizontal_parallax_deg: 0.0,
    })
}

/// Evaluates planetary, lunar, or minor body phenomena inside an existing [`Moment`].
///
/// Supports:
/// - `Planet` and `MinorPlanet`: phase angle, illuminated fraction, elongation, apparent diameter,
///   visual magnitude ($V$), and horizontal parallax.
/// - `FixedStar`: catalog visual magnitude, elongation from the Sun, $\alpha = 0^\circ$, $k = 1.0$.
/// - `LunarPoint` and `Hypothetical`: elongation from the Sun, zero diameter, and zero magnitude.
/// - `EclNut`: returns [`EphemerisError::UnsupportedBody`].
pub fn calc_phenomena_in_moment(
    moment: &Moment<'_>,
    body: BodyId,
    options: &CalcOptions,
) -> Result<PhenomenaResult, EphemerisError> {
    if body == BodyId::EclNut {
        return Err(EphemerisError::UnsupportedBody {
            body: "EclNut".to_owned(),
        });
    }

    // Set up ICRF options preserving observer origin and aberration mode for 3D vector evaluation
    let mut icrf_opts = *options;
    icrf_opts.plane = Plane::Equatorial;
    icrf_opts.equinox = Equinox::J2000;
    icrf_opts.nutation = NutationMode::Mean;
    icrf_opts.zodiac = Zodiac::Tropical;
    icrf_opts.include_speed = false;

    // Evaluate Sun position from observer in ICRF
    let sun_calc = calc_body_in_moment(moment, BodyId::Planet(PlanetId::Sun), &icrf_opts)?;
    let r_obs_sun_au = sun_calc.cartesian_au;
    let norm_sun = sun_calc.dist_au;

    // 1. Fixed stars
    if let BodyId::FixedStar(hip) = body {
        let star = find_star_by_hip(hip).ok_or_else(|| EphemerisError::UnsupportedBody {
            body: format!("FixedStar(HIP {})", hip.0),
        })?;
        return calc_star_entry_phenomena_in_moment(moment, star, options);
    }

    // 2. Synthetic lunar points
    if let BodyId::LunarPoint(_) = body {
        let lp_calc = calc_body_in_moment(moment, body, &icrf_opts)?;
        let dot_elong = vec_dot(&lp_calc.cartesian_au, &r_obs_sun_au);
        let norm_lp = vec_norm(&lp_calc.cartesian_au);
        let cos_elong = if norm_lp > 0.0 && norm_sun > 0.0 {
            (dot_elong / (norm_lp * norm_sun)).clamp(-1.0, 1.0)
        } else {
            1.0
        };
        let elongation_deg = cos_elong.acos().to_degrees();
        return Ok(PhenomenaResult {
            phase_angle_deg: 0.0,
            phase_illuminated: 1.0,
            elongation_deg,
            apparent_diameter_deg: 0.0,
            apparent_magnitude_v: 0.0,
            horizontal_parallax_deg: 0.0,
        });
    }

    // 3. Hypothetical bodies
    if let BodyId::Hypothetical(_) = body {
        let hyp_calc = calc_body_in_moment(moment, body, &icrf_opts)?;
        let dot_elong = vec_dot(&hyp_calc.cartesian_au, &r_obs_sun_au);
        let norm_hyp = vec_norm(&hyp_calc.cartesian_au);
        let cos_elong = if norm_hyp > 0.0 && norm_sun > 0.0 {
            (dot_elong / (norm_hyp * norm_sun)).clamp(-1.0, 1.0)
        } else {
            1.0
        };
        let elongation_deg = cos_elong.acos().to_degrees();
        return Ok(PhenomenaResult {
            phase_angle_deg: 0.0,
            phase_illuminated: 1.0,
            elongation_deg,
            apparent_diameter_deg: 0.0,
            apparent_magnitude_v: 0.0,
            horizontal_parallax_deg: 0.0,
        });
    }

    // 4. Physical bodies (Sun, Moon, Planets, Minor Planets)
    let targ_calc = calc_body_in_moment(moment, body, &icrf_opts)?;
    let delta_obs_au = targ_calc.cartesian_au;
    let delta_au = targ_calc.dist_au;
    let delta_km = targ_calc.dist_km;

    if body == BodyId::Planet(PlanetId::Sun) {
        let r_eq = body_equatorial_radius_km(body).unwrap_or(695_700.0);
        let apparent_diameter_deg = if delta_km > 0.0 {
            2.0 * (r_eq / delta_km).clamp(0.0, 1.0).asin().to_degrees()
        } else {
            0.0
        };
        let horizontal_parallax_deg = if delta_km > 0.0 {
            (R_EARTH_EQUATORIAL_KM / delta_km).clamp(0.0, 1.0).asin().to_degrees()
        } else {
            0.0
        };
        let apparent_magnitude_v = -26.74 + 5.0 * delta_au.max(1e-12).log10();
        return Ok(PhenomenaResult {
            phase_angle_deg: 0.0,
            phase_illuminated: 1.0,
            elongation_deg: 0.0,
            apparent_diameter_deg,
            apparent_magnitude_v,
            horizontal_parallax_deg,
        });
    }

    // Heliocentric vector of the target in consistent physical triangle (Astrometric / light-time retarded)
    let mut astrom_opts = icrf_opts;
    astrom_opts.aberration = AbMode::Astrometric;
    let targ_astrom = calc_body_in_moment(moment, body, &astrom_opts)?;
    let delta_phys_au = targ_astrom.cartesian_au;
    let delta_phys_dist_au = targ_astrom.dist_au;

    let mut helio_opts = astrom_opts;
    helio_opts.origin = Origin::Heliocentric;
    let helio_calc = calc_body_in_moment(moment, body, &helio_opts)?;
    let r_helio_au = helio_calc.cartesian_au;
    let r_helio_dist_au = helio_calc.dist_au;

    // Elongation from Sun (apparent angular separation between target and Sun as seen by observer)
    let dot_elong = vec_dot(&delta_obs_au, &r_obs_sun_au);
    let cos_elong = if delta_au > 0.0 && norm_sun > 0.0 {
        (dot_elong / (delta_au * norm_sun)).clamp(-1.0, 1.0)
    } else {
        1.0
    };
    let elongation_deg = cos_elong.acos().to_degrees();

    // Phase angle alpha: angle Sun–Target–Observer in consistent physical triangle
    let dot_phase = vec_dot(&r_helio_au, &delta_phys_au);
    let cos_phase = if r_helio_dist_au > 0.0 && delta_phys_dist_au > 0.0 {
        (dot_phase / (r_helio_dist_au * delta_phys_dist_au)).clamp(-1.0, 1.0)
    } else {
        1.0
    };
    let phase_angle_deg = cos_phase.acos().to_degrees();
    let phase_illuminated = ((1.0 + cos_phase) * 0.5).clamp(0.0, 1.0);

    // Apparent equatorial diameter
    let r_eq_km = body_equatorial_radius_km(body).unwrap_or(0.0);
    let apparent_diameter_deg = if r_eq_km > 0.0 && delta_km > 0.0 {
        2.0 * (r_eq_km / delta_km).clamp(0.0, 1.0).asin().to_degrees()
    } else {
        0.0
    };

    // Horizontal parallax
    let horizontal_parallax_deg = if delta_km > 0.0 {
        (R_EARTH_EQUATORIAL_KM / delta_km).clamp(0.0, 1.0).asin().to_degrees()
    } else {
        0.0
    };

    // Visual magnitude V from consistent physical triangle
    let apparent_magnitude_v = targ_calc
        .magnitude_v
        .or_else(|| compute_apparent_magnitude_v(body, &r_helio_au, &delta_phys_au, &r_obs_sun_au))
        .unwrap_or(0.0);

    Ok(PhenomenaResult {
        phase_angle_deg,
        phase_illuminated,
        elongation_deg,
        apparent_diameter_deg,
        apparent_magnitude_v,
        horizontal_parallax_deg,
    })
}

/// Evaluates planetary, lunar, or minor body phenomena at an arbitrary epoch and observer.
pub fn calc_phenomena(
    spk: &dyn SpkSource,
    time_table: &TimeTable,
    body: BodyId,
    input: TimeInput,
    options: CalcOptions,
) -> Result<PhenomenaResult, EphemerisError> {
    let observer = match options.origin {
        Origin::Topocentric(obs) => Some(obs),
        _ => None,
    };
    let moment = Moment::new(spk, time_table, input, observer)?;
    calc_phenomena_in_moment(&moment, body, &options)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_body_equatorial_radii() {
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Sun)), Some(695_700.0));
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Moon)), Some(1_737.40));
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Mercury)), Some(2_439.70));
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Venus)), Some(6_051.80));
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Earth)), Some(6_378.136_6));
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Mars)), Some(3_396.19));
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Jupiter)), Some(71_492.0));
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Saturn)), Some(60_268.0));
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Uranus)), Some(25_559.0));
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Neptune)), Some(24_764.0));
        assert_eq!(body_equatorial_radius_km(BodyId::Planet(PlanetId::Pluto)), Some(1_188.30));

        // Minor bodies
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(1))), Some(469.73));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(2))), Some(256.00));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(3))), Some(127.00));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(4))), Some(262.70));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(10))), Some(217.00));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(16))), Some(113.00));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(2060))), Some(103.00));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(5145))), Some(95.00));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(90377))), Some(500.00));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(136108))), Some(816.00));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(136199))), Some(1163.00));
        assert_eq!(body_equatorial_radius_km(BodyId::MinorPlanet(MpcNumber(136472))), Some(715.00));

        // Non-physical points
        assert_eq!(body_equatorial_radius_km(BodyId::EclNut), None);
        assert_eq!(
            body_equatorial_radius_km(BodyId::Planet(PlanetId::SolarSystemBarycenter)),
            None
        );
    }

    #[test]
    fn test_minor_body_lookup() {
        assert_eq!(MINOR_BODIES_PHENO_68.len(), 68);
        for pair in MINOR_BODIES_PHENO_68.windows(2) {
            if let [a, b] = pair {
                assert!(a.mpc < b.mpc, "Registry must be strictly sorted by MPC");
            }
        }

        let (h_ceres, g_ceres, r_ceres) = minor_body_pheno_lookup(MpcNumber(1));
        assert!((h_ceres - 3.34).abs() < 1e-4);
        assert!((g_ceres - 0.12).abs() < 1e-4);
        assert!((r_ceres - 469.73).abs() < 1e-4);

        // Unknown MPC number produces default H = 6.50, G = 0.15
        let (h_unk, g_unk, _) = minor_body_pheno_lookup(MpcNumber(999_999));
        assert!((h_unk - 6.50).abs() < 1e-4);
        assert!((g_unk - 0.15).abs() < 1e-4);
    }

    #[test]
    fn test_compute_apparent_magnitude_planets() {
        // Sun at 1 AU
        let sun_mag = compute_apparent_magnitude_v(
            BodyId::Planet(PlanetId::Sun),
            &[0.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
        );
        assert!(sun_mag.is_some());
        if let Some(m) = sun_mag {
            assert!((m - (-26.74)).abs() < 1e-3);
        }

        // Full Moon (alpha = 0) at mean distance 384,400 km
        let moon_delta_au = 384_400.0 / AU_KM;
        let moon_mag = compute_apparent_magnitude_v(
            BodyId::Planet(PlanetId::Moon),
            &[1.0, 0.0, 0.0],
            &[moon_delta_au, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
        );
        assert!(moon_mag.is_some());
        if let Some(m) = moon_mag {
            assert!((m - (-12.74)).abs() < 1e-2);
        }

        // Mercury at 1 AU heliocentric and 1 AU geocentric, alpha = 0
        let merc_mag = compute_apparent_magnitude_v(
            BodyId::Planet(PlanetId::Mercury),
            &[1.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
        );
        assert!(merc_mag.is_some());
        if let Some(m) = merc_mag {
            assert!((m - (-0.613)).abs() < 1e-3);
        }

        // Mars at 1 AU helio and 1 AU geo, alpha = 0
        let mars_mag = compute_apparent_magnitude_v(
            BodyId::Planet(PlanetId::Mars),
            &[1.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
        );
        assert!(mars_mag.is_some());
        if let Some(m) = mars_mag {
            assert!((m - (-1.601)).abs() < 1e-3);
        }

        // Jupiter at 1 AU helio and 1 AU geo, alpha = 0
        let jup_mag = compute_apparent_magnitude_v(
            BodyId::Planet(PlanetId::Jupiter),
            &[1.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
        );
        assert!(jup_mag.is_some());
        if let Some(m) = jup_mag {
            assert!((m - (-9.395)).abs() < 1e-3);
        }

        // Ceres at r = 1 AU, delta = 1 AU, alpha = 0
        let ceres_mag = compute_apparent_magnitude_v(
            BodyId::MinorPlanet(MpcNumber(1)),
            &[1.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0],
        );
        assert!(ceres_mag.is_some());
        if let Some(m) = ceres_mag {
            assert!((m - 3.34).abs() < 1e-2);
        }
    }

    #[test]
    fn test_phenomena_result_swe_attr20() {
        let res = PhenomenaResult {
            phase_angle_deg: 25.5,
            phase_illuminated: 0.951,
            elongation_deg: 45.2,
            apparent_diameter_deg: 0.533,
            apparent_magnitude_v: -26.74,
            horizontal_parallax_deg: 0.0024,
        };
        let attr = res.to_swe_attr20();
        assert_eq!(attr.len(), 20);
        assert_eq!(attr[0], 25.5);
        assert_eq!(attr[1], 0.951);
        assert_eq!(attr[2], 45.2);
        assert_eq!(attr[3], 0.533);
        assert_eq!(attr[4], -26.74);
        assert_eq!(attr[5], 0.0024);
        for &val in &attr[6..20] {
            assert_eq!(val, 0.0);
        }

        assert!((res.apparent_diameter_arcsec() - (0.533 * 3600.0)).abs() < 1e-9);
        assert!((res.horizontal_parallax_arcsec() - (0.0024 * 3600.0)).abs() < 1e-9);
    }

    #[test]
    fn test_phenomena_apparent_vs_astrometric_physical_triangle() {
        use crate::spk::AniseSpkBackend;
        use crate::types::JD_J2000;
        use std::path::PathBuf;

        let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("core parent")
            .parent()
            .expect("workspace root")
            .join("data/de440s.bsp");
        if !p.exists() {
            return;
        }
        let spk = AniseSpkBackend::from_bsp_paths(&[p]).expect("load de440s");
        let tt = TimeTable::embedded();

        let opts = CalcOptions {
            origin: Origin::Geocentric,
            plane: Plane::Ecliptic,
            equinox: Equinox::OfDate,
            nutation: NutationMode::True,
            aberration: AbMode::Apparent,
            zodiac: Zodiac::Tropical,
            include_speed: false,
        };

        let pheno = calc_phenomena(
            &spk,
            &tt,
            BodyId::Planet(PlanetId::Mars),
            TimeInput::from_tt_jd(JD_J2000),
            opts,
        )
        .expect("mars phenomena");

        assert!((0.0..=180.0).contains(&pheno.phase_angle_deg));
        assert!((0.0..=1.0).contains(&pheno.phase_illuminated));
        assert!((0.0..=180.0).contains(&pheno.elongation_deg));
        assert!(pheno.apparent_magnitude_v > -3.0 && pheno.apparent_magnitude_v < 3.0);
    }
}
