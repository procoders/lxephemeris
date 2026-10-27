//! Reference frames, IAU 2006/2000A Frame Bias + Precession + Nutation (`bp06`, `nut00a`, `num06a`),
//! long-term Vondrák 2011 precession fallback (`ltp`, `ltpecl`), obliquity (`obl06`),
//! 3D matrix/vector algebra, and 4th-order Cartesian stencil differentiation.

use crate::error::EphemerisError;
use crate::time::Tt;
use crate::types::{Equinox, NutationMode, Plane, JD_J2000};
use serde::{Deserialize, Serialize};
use sofars::coords::ecm06;
use sofars::pnp::{bp06, ltpb, ltpecl, numat, nut00a, nut00b, obl06};
use std::sync::OnceLock;

static R_ECL_J2000_CACHE: OnceLock<Mat3> = OnceLock::new();

#[must_use]
fn r_ecl_j2000_const() -> Mat3 {
    *R_ECL_J2000_CACHE.get_or_init(|| ecm06(JD_J2000, 0.0))
}

/// 3x3 rotation matrix (row-major).
pub type Mat3 = [[f64; 3]; 3];

/// 3D Cartesian vector.
pub type Vec3 = [f64; 3];

/// Lower bound of the IAU 2006 5th-order polynomial validity window (`1000-01-01 00:00 TT` Julian = `JD 2086302.5`).
pub const JD_IAU2006_MIN: f64 = 2_086_302.5;

/// Upper bound of the IAU 2006 5th-order polynomial validity window (`3000-01-01 00:00 TT` Gregorian).
pub const JD_IAU2006_MAX: f64 = 2_816_787.5;

/// Returns the 3x3 identity matrix.
#[must_use]
pub const fn mat_identity() -> Mat3 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

/// Returns the 3x3 zero matrix.
#[must_use]
pub const fn mat_zero() -> Mat3 {
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
}

/// Multiplies two 3x3 matrices `a * b`.
#[must_use]
pub const fn mat_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1],
            a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1],
            a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2],
        ],
        [
            a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0],
            a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1],
            a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2],
        ],
    ]
}

/// Multiplies a 3x3 matrix `m` by a 3D column vector `v`.
#[must_use]
pub fn mat_vec(m: &Mat3, v: &Vec3) -> Vec3 {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Transposes a 3x3 matrix.
#[must_use]
pub const fn mat_transpose(m: &Mat3) -> Mat3 {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

/// Rotation matrix $R_x(\theta)$ around the X-axis by `angle_rad` radians.
#[must_use]
pub fn rot_x(angle_rad: f64) -> Mat3 {
    let (s, c) = angle_rad.sin_cos();
    [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]]
}

/// Rotation matrix $R_y(\theta)$ around the Y-axis by `angle_rad` radians.
#[must_use]
pub fn rot_y(angle_rad: f64) -> Mat3 {
    let (s, c) = angle_rad.sin_cos();
    [[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]]
}

/// Rotation matrix $R_z(\theta)$ around the Z-axis by `angle_rad` radians.
#[must_use]
pub fn rot_z(angle_rad: f64) -> Mat3 {
    let (s, c) = angle_rad.sin_cos();
    [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]
}

/// Euclidean norm $\|\mathbf{v}\|_2$, scaled against intermediate overflow/underflow on extreme vectors.
#[must_use]
pub fn vec_norm(v: &Vec3) -> f64 {
    let (x, y, z) = (v[0], v[1], v[2]);
    if x.is_nan() || y.is_nan() || z.is_nan() {
        return f64::NAN;
    }
    let m = x.abs().max(y.abs()).max(z.abs());
    if m == 0.0 {
        0.0
    } else if m.is_infinite() {
        f64::INFINITY
    } else if !(1.0e-150..=1.0e150).contains(&m) {
        let (sx, sy, sz) = (x / m, y / m, z / m);
        m * (sx * sx + sy * sy + sz * sz).sqrt()
    } else {
        (x * x + y * y + z * z).sqrt()
    }
}

/// Dot product $\mathbf{a} \cdot \mathbf{b}$.
#[must_use]
pub fn vec_dot(a: &Vec3, b: &Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product $\mathbf{a} \times \mathbf{b}$.
#[must_use]
pub fn vec_cross(a: &Vec3, b: &Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Vector addition $\mathbf{a} + \mathbf{b}$.
#[must_use]
pub fn vec_add(a: &Vec3, b: &Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Vector subtraction $\mathbf{a} - \mathbf{b}$.
#[must_use]
pub fn vec_sub(a: &Vec3, b: &Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scalar multiplication $s \mathbf{v}$.
#[must_use]
pub fn vec_scale(v: &Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Unit vector $\mathbf{v} / \|\mathbf{v}\|_2$ (returns `[0.0; 3]` if norm is zero, `[NaN; 3]` if non-finite).
#[must_use]
pub fn vec_unit(v: &Vec3) -> Vec3 {
    let n = vec_norm(v);
    if !n.is_finite() {
        [f64::NAN, f64::NAN, f64::NAN]
    } else if n > 0.0 {
        [v[0] / n, v[1] / n, v[2] / n]
    } else {
        [0.0, 0.0, 0.0]
    }
}

/// Returns `true` if `jd_tt` lies within the `[1000 CE, 3000 CE]` domain of the IAU 2006
/// 5th-order precession polynomial (`sofars::pnp::bp06`).
#[must_use]
pub fn is_within_iau2006_domain(jd_tt: f64) -> bool {
    (JD_IAU2006_MIN..=JD_IAU2006_MAX).contains(&jd_tt)
}

/// Enforces that `jd_tt` is within the strict IAU 2006 polynomial validity window `[1000 CE, 3000 CE]`.
pub fn check_iau2006_validity(jd_tt: f64) -> Result<(), EphemerisError> {
    if is_within_iau2006_domain(jd_tt) {
        Ok(())
    } else {
        Err(EphemerisError::OutOfModelValidity {
            model: "IAU2006_Precession_bp06".to_owned(),
            epoch_jd_tt: jd_tt,
            min_jd_tt: JD_IAU2006_MIN,
            max_jd_tt: JD_IAU2006_MAX,
        })
    }
}

/// Ecliptic obliquity and IAU 2000A nutation angles (`SE_ECL_NUT` / `R2.5`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EclNutInfo {
    /// True obliquity of the ecliptic $\varepsilon = \varepsilon_0 + \Delta\varepsilon$ in degrees.
    pub true_obliquity_deg: f64,
    /// Mean obliquity of the ecliptic $\varepsilon_0$ in degrees.
    pub mean_obliquity_deg: f64,
    /// Nutation in longitude $\Delta\psi$ in degrees.
    pub nutation_lon_deg: f64,
    /// Nutation in obliquity $\Delta\varepsilon$ in degrees.
    pub nutation_obl_deg: f64,
}

/// Precomputed frame transformation matrices and nutation/obliquity angles at a single `TT` epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameTransform {
    /// Epoch in `TT`.
    pub tt: Tt,
    /// Mean obliquity of the ecliptic $\varepsilon_0$ in radians.
    pub eps0_rad: f64,
    /// Nutation in longitude $\Delta\psi$ in radians (IAU 2000A).
    pub dpsi_rad: f64,
    /// Nutation in obliquity $\Delta\varepsilon$ in radians (IAU 2000A).
    pub deps_rad: f64,
    /// True obliquity of the ecliptic $\varepsilon = \varepsilon_0 + \Delta\varepsilon$ in radians.
    pub true_eps_rad: f64,
    /// ICRF $\to$ Mean Equatorial of Date (`bp06` within `[1000, 3000]`, Vondrák 2011 `ltp` outside).
    pub r_bp: Mat3,
    /// ICRF $\to$ True Equatorial of Date (`num06a * r_bp`).
    pub r_bpn: Mat3,
    /// Time derivative $\dot{\mathbf{R}}_{\text{bpn}}(t)$ of ICRF $\to$ True Equatorial of Date (per day).
    pub r_bpn_dot: Mat3,
    /// ICRF $\to$ Mean Ecliptic of Date (`R_x(eps0) * r_bp`).
    pub r_ecl_mean: Mat3,
    /// ICRF $\to$ True Ecliptic of Date (`R_x(eps0 + deps) * r_bpn`).
    pub r_ecl_true: Mat3,
    /// ICRF $\to$ Mean Ecliptic J2000 (`ecm06(JD_J2000, 0.0)`).
    pub r_ecl_j2000: Mat3,
}

/// Evaluates the time derivative $(\dot{\Delta\psi}, \dot{\Delta\varepsilon})$ in radians/day
/// at `jd_tt` using centered finite difference on the 77-term IAU 2000B nutation series (`nut00b`).
///
/// Matches the full 1,365-term `nut00a` derivative to `< 1e-8 deg/day` at ~18x lower cost.
#[must_use]
pub fn nutation_rate_rad_day(jd_tt: f64) -> (f64, f64) {
    let h = 1.0e-3;
    let tt0 = Tt::from_jd(jd_tt).0;
    let tt_m = tt0.add_days(-h);
    let tt_p = tt0.add_days(h);
    let (psi_m, eps_m) = nut00b(tt_m.jd1, tt_m.jd2);
    let (psi_p, eps_p) = nut00b(tt_p.jd1, tt_p.jd2);
    ((psi_p - psi_m) / (2.0 * h), (eps_p - eps_m) / (2.0 * h))
}

/// Returns `a + s * b` for two 3x3 matrices.
#[must_use]
pub const fn mat_add_scaled(a: &Mat3, b: &Mat3, s: f64) -> Mat3 {
    [
        [
            a[0][0] + s * b[0][0],
            a[0][1] + s * b[0][1],
            a[0][2] + s * b[0][2],
        ],
        [
            a[1][0] + s * b[1][0],
            a[1][1] + s * b[1][1],
            a[1][2] + s * b[1][2],
        ],
        [
            a[2][0] + s * b[2][0],
            a[2][1] + s * b[2][1],
            a[2][2] + s * b[2][2],
        ],
    ]
}

/// Evaluates `(deps0_rad_day, dr_bp_day)` at `tt0` using a forward step of $10^{-3}\text{ d}$.
#[must_use]
pub fn precession_rate_per_day(tt0: Tt, eps0_0: f64, r_bp_0: &Mat3) -> (f64, Mat3) {
    let h = 1.0e-3;
    let tt_p = tt0.0.add_days(h);
    let jd_p = tt_p.to_f64();
    let (eps0_p, r_bp_p) = if is_within_iau2006_domain(jd_p) {
        let e0 = obl06(tt_p.jd1, tt_p.jd2);
        let (_, _, rbp) = bp06(tt_p.jd1, tt_p.jd2);
        (e0, rbp)
    } else {
        let epj = 2000.0 + (jd_p - JD_J2000) / 365.25;
        let rbp = ltpb(epj);
        let eq_pole = [rbp[2][0], rbp[2][1], rbp[2][2]];
        let ecl_pole = ltpecl(epj);
        let cos_eps = vec_dot(&eq_pole, &ecl_pole).clamp(-1.0, 1.0);
        (cos_eps.acos(), rbp)
    };
    let inv_h = 1.0 / h;
    let deps0_rate = (eps0_p - eps0_0) * inv_h;
    let dr_bp = [
        [
            (r_bp_p[0][0] - r_bp_0[0][0]) * inv_h,
            (r_bp_p[0][1] - r_bp_0[0][1]) * inv_h,
            (r_bp_p[0][2] - r_bp_0[0][2]) * inv_h,
        ],
        [
            (r_bp_p[1][0] - r_bp_0[1][0]) * inv_h,
            (r_bp_p[1][1] - r_bp_0[1][1]) * inv_h,
            (r_bp_p[1][2] - r_bp_0[1][2]) * inv_h,
        ],
        [
            (r_bp_p[2][0] - r_bp_0[2][0]) * inv_h,
            (r_bp_p[2][1] - r_bp_0[2][1]) * inv_h,
            (r_bp_p[2][2] - r_bp_0[2][2]) * inv_h,
        ],
    ];
    (deps0_rate, dr_bp)
}

impl FrameTransform {
    /// Computes all frame rotation matrices and nutation/obliquity angles at epoch `tt`.
    ///
    /// Within `[1000 CE, 3000 CE]`, uses IAU 2006/2000A (`bp06`, `obl06`, `nut00a`, `num06a`).
    /// Outside `[1000 CE, 3000 CE]`, automatically switches to Vondrák, Capitaine & Wallace (2011)
    /// long-term precession (`sofars::pnp::ltp` / `ltpecl`).
    #[must_use]
    pub fn at_tt(tt: Tt) -> Self {
        let (dpsi_rad, deps_rad) = nut00a(tt.0.jd1, tt.0.jd2);
        Self::at_tt_with_nutation(tt, dpsi_rad, deps_rad)
    }

    /// Computes all frame rotation matrices at epoch `tt` using pre-evaluated `(dpsi_rad, deps_rad)`.
    #[must_use]
    pub fn at_tt_with_nutation(tt: Tt, dpsi_rad: f64, deps_rad: f64) -> Self {
        let jd1 = tt.0.jd1;
        let jd2 = tt.0.jd2;
        let jd_tt = tt.jd_f64();

        let (eps0_rad, r_bp) = if is_within_iau2006_domain(jd_tt) {
            let eps0 = obl06(jd1, jd2);
            let (_, _, rbp) = bp06(jd1, jd2);
            (eps0, rbp)
        } else {
            // Julian epoch for Vondrák et al. (2011) long-term precession + frame bias (ltpb)
            let epj = 2000.0 + (jd_tt - JD_J2000) / 365.25;
            let rbp = ltpb(epj);
            let eq_pole = [rbp[2][0], rbp[2][1], rbp[2][2]];
            let ecl_pole = ltpecl(epj);
            let cos_eps = vec_dot(&eq_pole, &ecl_pole).clamp(-1.0, 1.0);
            (cos_eps.acos(), rbp)
        };

        Self::from_precomputed_bp(tt, eps0_rad, r_bp, dpsi_rad, deps_rad)
    }

    /// Computes all frame rotation matrices from pre-evaluated `(eps0_rad, r_bp, dpsi_rad, deps_rad)`.
    #[must_use]
    pub fn from_precomputed_bp(
        tt: Tt,
        eps0_rad: f64,
        r_bp: Mat3,
        dpsi_rad: f64,
        deps_rad: f64,
    ) -> Self {
        let true_eps_rad = eps0_rad + deps_rad;
        let r_nut = numat(eps0_rad, dpsi_rad, deps_rad);
        let r_bpn = mat_mul(&r_nut, &r_bp);
        let r_ecl_mean = mat_mul(&rot_x(eps0_rad), &r_bp);
        let r_ecl_true = mat_mul(&rot_x(true_eps_rad), &r_bpn);
        let r_ecl_j2000 = r_ecl_j2000_const();

        let h = 1.0e-4_f64;
        let jd_tt = tt.jd_f64();
        let (ddpsi, ddeps) = nutation_rate_rad_day(jd_tt);
        let (deps0, dr_bp) = precession_rate_per_day(tt, eps0_rad, &r_bp);

        let eps0_p = eps0_rad + h * deps0;
        let r_bp_p = mat_add_scaled(&r_bp, &dr_bp, h);
        let r_nut_p = numat(eps0_p, dpsi_rad + h * ddpsi, deps_rad + h * ddeps);
        let r_bpn_p = mat_mul(&r_nut_p, &r_bp_p);

        let eps0_m = eps0_rad - h * deps0;
        let r_bp_m = mat_add_scaled(&r_bp, &dr_bp, -h);
        let r_nut_m = numat(eps0_m, dpsi_rad - h * ddpsi, deps_rad - h * ddeps);
        let r_bpn_m = mat_mul(&r_nut_m, &r_bp_m);

        let inv_2h = 1.0 / (2.0 * h);
        let r_bpn_dot = [
            [
                (r_bpn_p[0][0] - r_bpn_m[0][0]) * inv_2h,
                (r_bpn_p[0][1] - r_bpn_m[0][1]) * inv_2h,
                (r_bpn_p[0][2] - r_bpn_m[0][2]) * inv_2h,
            ],
            [
                (r_bpn_p[1][0] - r_bpn_m[1][0]) * inv_2h,
                (r_bpn_p[1][1] - r_bpn_m[1][1]) * inv_2h,
                (r_bpn_p[1][2] - r_bpn_m[1][2]) * inv_2h,
            ],
            [
                (r_bpn_p[2][0] - r_bpn_m[2][0]) * inv_2h,
                (r_bpn_p[2][1] - r_bpn_m[2][1]) * inv_2h,
                (r_bpn_p[2][2] - r_bpn_m[2][2]) * inv_2h,
            ],
        ];

        Self {
            tt,
            eps0_rad,
            dpsi_rad,
            deps_rad,
            true_eps_rad,
            r_bp,
            r_bpn,
            r_bpn_dot,
            r_ecl_mean,
            r_ecl_true,
            r_ecl_j2000,
        }
    }

    /// Returns [`EclNutInfo`] (`true_obliquity_deg`, `mean_obliquity_deg`, `nutation_lon_deg`, `nutation_obl_deg`).
    #[must_use]
    pub fn ecl_nut_info(&self) -> EclNutInfo {
        EclNutInfo {
            true_obliquity_deg: self.true_eps_rad.to_degrees(),
            mean_obliquity_deg: self.eps0_rad.to_degrees(),
            nutation_lon_deg: self.dpsi_rad.to_degrees(),
            nutation_obl_deg: self.deps_rad.to_degrees(),
        }
    }

    /// Selects the 3x3 rotation matrix from ICRF to the requested `(Plane, Equinox, NutationMode)`.
    #[must_use]
    pub fn icrf_to_frame_matrix(
        &self,
        plane: Plane,
        equinox: Equinox,
        nutation: NutationMode,
    ) -> Mat3 {
        match (plane, equinox, nutation) {
            (Plane::Equatorial, Equinox::J2000, _) => mat_identity(),
            (Plane::Equatorial, Equinox::OfDate, NutationMode::Mean) => self.r_bp,
            (Plane::Equatorial, Equinox::OfDate, NutationMode::True) => self.r_bpn,
            (Plane::Ecliptic, Equinox::J2000, _) => self.r_ecl_j2000,
            (Plane::Ecliptic, Equinox::OfDate, NutationMode::Mean) => self.r_ecl_mean,
            (Plane::Ecliptic, Equinox::OfDate, NutationMode::True) => self.r_ecl_true,
            (Plane::Horizontal, _, _) => self.r_bpn,
        }
    }

    /// Evaluates the time derivative $\dot{\mathbf{R}}(t)$ (per day) of the 3x3 rotation matrix
    /// from ICRF to the requested `(Plane, Equinox, NutationMode)`.
    #[must_use]
    pub fn icrf_to_frame_matrix_rate(
        &self,
        plane: Plane,
        equinox: Equinox,
        nutation: NutationMode,
    ) -> Mat3 {
        if plane == Plane::Horizontal {
            return self.r_bpn_dot;
        }
        if equinox == Equinox::J2000 {
            return mat_zero();
        }
        if plane == Plane::Equatorial && nutation == NutationMode::True {
            return self.r_bpn_dot;
        }
        let h = 1.0e-4_f64;
        let jd_tt = self.tt.jd_f64();
        let (ddpsi, ddeps) = nutation_rate_rad_day(jd_tt);
        let (deps0, dr_bp) = precession_rate_per_day(self.tt, self.eps0_rad, &self.r_bp);

        let tt_p = Tt(self.tt.0.add_days(h));
        let tt_m = Tt(self.tt.0.add_days(-h));
        let ft_p = Self::from_precomputed_bp(
            tt_p,
            self.eps0_rad + h * deps0,
            mat_add_scaled(&self.r_bp, &dr_bp, h),
            self.dpsi_rad + h * ddpsi,
            self.deps_rad + h * ddeps,
        );
        let ft_m = Self::from_precomputed_bp(
            tt_m,
            self.eps0_rad - h * deps0,
            mat_add_scaled(&self.r_bp, &dr_bp, -h),
            self.dpsi_rad - h * ddpsi,
            self.deps_rad - h * ddeps,
        );
        let r_p = ft_p.icrf_to_frame_matrix(plane, equinox, nutation);
        let r_m = ft_m.icrf_to_frame_matrix(plane, equinox, nutation);
        let inv_2h = 1.0 / (2.0 * h);
        [
            [
                (r_p[0][0] - r_m[0][0]) * inv_2h,
                (r_p[0][1] - r_m[0][1]) * inv_2h,
                (r_p[0][2] - r_m[0][2]) * inv_2h,
            ],
            [
                (r_p[1][0] - r_m[1][0]) * inv_2h,
                (r_p[1][1] - r_m[1][1]) * inv_2h,
                (r_p[1][2] - r_m[1][2]) * inv_2h,
            ],
            [
                (r_p[2][0] - r_m[2][0]) * inv_2h,
                (r_p[2][1] - r_m[2][1]) * inv_2h,
                (r_p[2][2] - r_m[2][2]) * inv_2h,
            ],
        ]
    }
}

/// Spherical coordinates `(lon_deg, lat_deg, dist)` in degrees and distance units (AU or km).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SphericalCoords {
    /// Longitude or Right Ascension in degrees `[0, 360)`.
    pub lon_deg: f64,
    /// Latitude or Declination in degrees `[-90, +90]`.
    pub lat_deg: f64,
    /// Radial distance $r$.
    pub dist: f64,
}

/// Full spherical state `(lon_deg, lat_deg, dist, lon_speed_deg_day, lat_speed_deg_day, dist_speed_day)`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SphericalState {
    /// Longitude or Right Ascension in degrees `[0, 360)`.
    pub lon_deg: f64,
    /// Latitude or Declination in degrees `[-90, +90]`.
    pub lat_deg: f64,
    /// Radial distance $r$ in AU.
    pub dist_au: f64,
    /// Longitudinal angular velocity $\dot{\lambda}$ or $\dot{\alpha}$ in degrees/day.
    pub lon_speed_deg_day: f64,
    /// Latitudinal angular velocity $\dot{\beta}$ or $\dot{\delta}$ in degrees/day.
    pub lat_speed_deg_day: f64,
    /// Radial velocity $\dot{r}$ in AU/day.
    pub dist_speed_au_day: f64,
}

/// Converts a 3D Cartesian vector `u = [x, y, z]` into [`SphericalCoords`].
#[must_use]
pub fn cart_to_spherical(u: &Vec3) -> SphericalCoords {
    let (x, y, z) = (u[0], u[1], u[2]);
    let dist = vec_norm(u);
    if !dist.is_finite() {
        return SphericalCoords {
            lon_deg: f64::NAN,
            lat_deg: f64::NAN,
            dist: f64::NAN,
        };
    }
    let lon_deg = y.atan2(x).to_degrees().rem_euclid(360.0);
    let lat_deg = if dist > 0.0 {
        (z / dist).clamp(-1.0, 1.0).asin().to_degrees()
    } else {
        0.0
    };
    SphericalCoords {
        lon_deg,
        lat_deg,
        dist,
    }
}

/// Converts a 3D Cartesian position `u` (AU) and velocity `du` (AU/day) into [`SphericalState`],
/// handling the polar singularity $|\beta| \to 90^\circ$ gracefully without division by zero.
#[must_use]
pub fn cart_state_to_spherical(u: &Vec3, du: &Vec3) -> SphericalState {
    let (x, y, z) = (u[0], u[1], u[2]);
    let (dx, dy, dz) = (du[0], du[1], du[2]);
    let r2 = x * x + y * y + z * z;
    let r = r2.sqrt();
    let rho2 = x * x + y * y;
    let rho = rho2.sqrt();

    let lon_deg = y.atan2(x).to_degrees().rem_euclid(360.0);
    let lat_deg = if r > 0.0 {
        (z / r).clamp(-1.0, 1.0).asin().to_degrees()
    } else {
        0.0
    };

    let dist_speed_au_day = if r > 1e-15 {
        (x * dx + y * dy + z * dz) / r
    } else {
        0.0
    };

    // Guard against exact celestial/ecliptic poles (rho -> 0)
    let (lon_speed_deg_day, lat_speed_deg_day) = if rho > 1e-12 * r.max(1.0) && r2 > 1e-30 {
        let dlon_rad = (x * dy - y * dx) / rho2;
        let dlat_rad = (rho * dz - z * (x * dx + y * dy) / rho) / r2;
        (dlon_rad.to_degrees(), dlat_rad.to_degrees())
    } else if r > 1e-15 {
        let v_perp = (dx * dx + dy * dy).sqrt();
        let sign = if z >= 0.0 { -1.0 } else { 1.0 };
        (0.0, (sign * v_perp / r).to_degrees())
    } else {
        (0.0, 0.0)
    };

    SphericalState {
        lon_deg,
        lat_deg,
        dist_au: r,
        lon_speed_deg_day,
        lat_speed_deg_day,
        dist_speed_au_day,
    }
}

/// Evaluates the 4-point 4th-order symmetric finite-difference derivative:
/// $$\mathbf{v}(t_0) = \frac{\mathbf{r}(t_0 - 2h) - 8\mathbf{r}(t_0 - h) + 8\mathbf{r}(t_0 + h) - \mathbf{r}(t_0 + 2h)}{12 h}$$
#[must_use]
pub fn stencil_4pt_derivative(
    m2: &Vec3,
    m1: &Vec3,
    p1: &Vec3,
    p2: &Vec3,
    h_days: f64,
) -> Vec3 {
    let inv_12h = 1.0 / (12.0 * h_days);
    [
        (m2[0] - 8.0 * m1[0] + 8.0 * p1[0] - p2[0]) * inv_12h,
        (m2[1] - 8.0 * m1[1] + 8.0 * p1[1] - p2[1]) * inv_12h,
        (m2[2] - 8.0 * m1[2] + 8.0 * p1[2] - p2[2]) * inv_12h,
    ]
}

/// Converts a [`SphericalState`] into 3D Cartesian position (AU) and velocity (AU/day).
#[must_use]
pub fn spherical_state_to_cart(s: &SphericalState, include_speed: bool) -> (Vec3, Vec3) {
    let l_r = s.lon_deg.to_radians();
    let b_r = s.lat_deg.to_radians();
    let r = s.dist_au;
    let (sin_l, cos_l) = l_r.sin_cos();
    let (sin_b, cos_b) = b_r.sin_cos();
    let pos = [r * cos_b * cos_l, r * cos_b * sin_l, r * sin_b];
    let vel = if include_speed {
        let dl_r = s.lon_speed_deg_day.to_radians();
        let db_r = s.lat_speed_deg_day.to_radians();
        let dr = s.dist_speed_au_day;
        [
            dr * cos_b * cos_l - r * db_r * sin_b * cos_l - r * dl_r * cos_b * sin_l,
            dr * cos_b * sin_l - r * db_r * sin_b * sin_l + r * dl_r * cos_b * cos_l,
            dr * sin_b + r * db_r * cos_b,
        ]
    } else {
        [0.0, 0.0, 0.0]
    };
    (pos, vel)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_frame_transform_j2000_obliquity_and_orthogonality() {
        let ft = FrameTransform::at_tt(Tt::from_jd(JD_J2000));
        // Mean obliquity at J2000.0 is 84381.406" = 23.4392794444 deg
        let eps0_arcsec = ft.eps0_rad.to_degrees() * 3600.0;
        assert!((eps0_arcsec - 84_381.406).abs() < 1e-4);

        // Check orthogonality of R_ecl_true: R * R^T == I
        let prod = mat_mul(&ft.r_ecl_true, &mat_transpose(&ft.r_ecl_true));
        for (i, row) in prod.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((val - expected).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn test_iau2006_domain_guard_and_vondrak_fallback() {
        assert!(is_within_iau2006_domain(JD_J2000));
        assert!(check_iau2006_validity(JD_J2000).is_ok());
        // Year 500 CE is outside [1000 CE, 3000 CE]
        let jd_500 = 1_903_682.5;
        assert!(!is_within_iau2006_domain(jd_500));
        assert!(matches!(
            check_iau2006_validity(jd_500),
            Err(EphemerisError::OutOfModelValidity { .. })
        ));
        // FrameTransform::at_tt seamlessly uses Vondrák 2011 outside [1000, 3000]
        let ft_500 = FrameTransform::at_tt(Tt::from_jd(jd_500));
        assert!(ft_500.eps0_rad.to_degrees() > 23.0 && ft_500.eps0_rad.to_degrees() < 24.0);
    }

    #[test]
    fn test_4pt_stencil_exact_on_cubic_and_wrap_free() {
        // Test across 359.99 deg -> 0.01 deg wrap-around on unit circle:
        // r(t) = [cos(omega * t), sin(omega * t), 0], omega = 1.0 rad/day
        let omega = 1.0_f64;
        let h = 1e-3_f64;
        let pos_at = |t: f64| -> Vec3 { [(omega * t).cos(), (omega * t).sin(), 0.0] };
        let m2 = pos_at(-2.0 * h);
        let m1 = pos_at(-h);
        let r0 = pos_at(0.0);
        let p1 = pos_at(h);
        let p2 = pos_at(2.0 * h);
        let v0 = stencil_4pt_derivative(&m2, &m1, &p1, &p2, h);
        let sph = cart_state_to_spherical(&r0, &v0);
        assert!((sph.lon_speed_deg_day - omega.to_degrees()).abs() < 1e-11);
    }

    #[test]
    fn test_icrf_to_frame_matrix_rate() {
        let ft = FrameTransform::at_tt(Tt::from_jd(JD_J2000));
        // J2000 rate must be identically zero
        let r_rate_j2000 = ft.icrf_to_frame_matrix_rate(Plane::Ecliptic, Equinox::J2000, NutationMode::Mean);
        assert_eq!(r_rate_j2000, mat_zero());

        // OfDate rate must reflect precession rate (~50 arcsec/yr ~ 3.8e-7 rad/day)
        let r_rate_ecl = ft.icrf_to_frame_matrix_rate(Plane::Ecliptic, Equinox::OfDate, NutationMode::True);
        let mut max_rate = 0.0_f64;
        for row in &r_rate_ecl {
            for &val in row {
                assert!(val.is_finite());
                max_rate = max_rate.max(val.abs());
            }
        }
        assert!(max_rate > 1e-8 && max_rate < 1e-5, "Precession-nutation rate {max_rate} must be in [1e-8, 1e-5] rad/day");
    }

    #[test]
    fn test_icrf_to_frame_matrix_rate_horizontal_j2000() {
        let ft = FrameTransform::at_tt(Tt::from_jd(JD_J2000));
        // Plane::Horizontal with Equinox::J2000 must NOT be zero; it must equal r_bpn_dot
        let r_rate_hz = ft.icrf_to_frame_matrix_rate(Plane::Horizontal, Equinox::J2000, NutationMode::True);
        assert_eq!(r_rate_hz, ft.r_bpn_dot);
        assert_ne!(r_rate_hz, mat_zero());
    }
}
