//! Module to handle common constants, coordinate types, and utility functions shared across
//! GP parsing and SGP4 propagation.

// ------------------
// External Libraries
// ------------------
use std::f64::consts::PI;

// ------------------
// Internal Libraries
// ------------------

// -------
// Structs
// -------

/// World Geodetic System (WGS) parameters
///
/// This struct contains the important Earth parameters defined by different WGS
/// standards (e.g. WGS-72, WGS-84).
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
#[derive(Default, Clone, Copy)]
pub struct Wgs {
    /// Earth's standard gravitational parameter \[km^3/s^2\]
    pub mu: f64,

    /// Earth's equatorial radius \[km\]
    pub r_earth_eq: f64,

    /// Earth's J2 harmonic \[\]
    pub j2: f64,

    /// k2 constant \[Earth Radii^2\]
    pub k2: f64,

    /// Earth's J3 harmonic \[\]
    pub j3: f64,

    /// Earth's J4 harmonic \[\]
    pub j4: f64,

    /// k4 constant \[Earth Radii^4\]
    pub k4: f64,

    /// The square root of the standard gravitational parameter \[Earth radii^1.5 / min\]
    pub ke: f64,

    /// The inverse of ke \[min / Earth radii^1.5\]
    pub tumin: f64,
}

/// Satellite state vector
///
/// Position \[km\] and velocity \[km/s\] of a satellite in a specified coordinate frame.
#[derive(Default, Clone, Copy)]
pub struct StateVector {
    /// Position X component \[km\]
    pub r_x: f64,

    /// Position Y component \[km\]
    pub r_y: f64,

    /// Position Z component \[km\]
    pub r_z: f64,

    /// Velocity X component \[km/s\]
    pub v_x: f64,

    /// Velocity Y component \[km/s\]
    pub v_y: f64,

    /// Velocity Z component \[km/s\]
    pub v_z: f64,

    /// Coordinate frame of the state vector
    pub coordinate_frame: CoordinateFrame,
}

// -----
// Enums
// -----

/// Coordinate frames
///
/// Represents the coordinate frame used for a [`StateVector`].
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::CoordinateFrame;
///
/// let frame_teme = CoordinateFrame::TEME;
/// let frame_j2000 = CoordinateFrame::J2000;
/// ```
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateFrame {
    /// J2000, an Earth-centered inertial (ECI) coordinate frame
    #[default]
    J2000,

    /// True Equator Mean Equinox (TEME), an Earth-centered inertial (ECI) coordinate frame
    TEME,
}

// ---------
// Constants
// ---------

/// Fundamental and derived constants for WGS-72
///
/// Earth model parameters used as the default for TLE/GP processing with SGP4.
///
/// - `mu`: 398600.8 — standard gravitational parameter \[km^3 / s^2\]
/// - `r_earth_eq`: 6378.135 — Earth's equatorial radius \[km\]
/// - `j2`: 0.001082616 — second zonal harmonic (Earth's oblateness)
/// - `k2`: 0.000541308 — `0.5 * j2` \[Earth Radii^2\]
/// - `j3`: -0.00000253881 — third zonal harmonic (pear-shaped component)
/// - `j4`: -0.00000165597 — fourth zonal harmonic
/// - `k4`: 0.00000062098875 — `-3/8 * j4` \[Earth Radii^4\]
/// - `ke`: 0.07436691613317 — square root of `mu` \[Earth radii^1.5 / min\]
/// - `tumin`: 13.44683969695931 — inverse of `ke` \[min / Earth radii^1.5\]
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
pub const WGS72: Wgs = Wgs {
    mu: 398600.8,
    r_earth_eq: 6378.135,
    j2: 0.001082616,
    k2: 0.000541308,
    j3: -0.00000253881,
    j4: -0.00000165597,
    k4: 0.00000062098875,
    ke: 0.07436691613317,
    tumin: 13.44683969695931,
};

/// Fundamental and derived constants for WGS-84
///
/// Earth model parameters for the WGS-84 geodetic system.
///
/// - `mu`: 398600.5 — standard gravitational parameter \[km^3 / s^2\]
/// - `r_earth_eq`: 6378.137 — Earth's equatorial radius \[km\]
/// - `j2`: 0.00108262998905 — second zonal harmonic (Earth's oblateness)
/// - `k2`: 0.000541314994525 — `0.5 * j2` \[Earth Radii^2\]
/// - `j3`: -0.00000253215306 — third zonal harmonic (pear-shaped component)
/// - `j4`: -0.00000161098761 — fourth zonal harmonic
/// - `k4`: 0.0000006041203538 — `-3/8 * j4` \[Earth Radii^4\]
/// - `ke`: 0.07436685316871 — square root of `mu` \[Earth radii^1.5 / min\]
/// - `tumin`: 13.44685108204498 — inverse of `ke` \[min / Earth radii^1.5\]
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
pub const WGS84: Wgs = Wgs {
    mu: 398600.5,
    r_earth_eq: 6378.137,
    j2: 0.00108262998905,
    k2: 0.000541314994525,
    j3: -0.00000253215306,
    j4: -0.00000161098761,
    k4: 0.0000006041203538,
    ke: 0.07436685316871,
    tumin: 13.44685108204498,
};

// ---------
// Functions
// ---------

/// Convert an angle from degrees to radians
///
/// Multiplies the input angle by `pi / 180`.
///
/// # Arguments
/// * `theta` - The angle in degrees
///
/// # Returns
/// * `theta_rad` - The angle in radians
///
/// # Examples
/// ```rust
/// use std::f64::consts::PI;
/// use mako_sgp4::common::deg2rad;
///
/// // Define some angle in degrees
/// let theta = 90.0; // Degrees
///
/// // Convert the angle to radians
/// let theta_rad = deg2rad(theta);
///
/// // Assert the value is equal to the correct value
/// assert!((theta_rad - PI / 2.0).abs() < 1e-12);
/// ```
pub fn deg2rad(theta: f64) -> f64 {
    // Convert to radians
    let theta_rad = PI / 180. * theta;

    // Return theta in radians
    return theta_rad;
}

/// Calculate the orbital period from semi-major axis and gravitational parameter
///
/// Uses Kepler's third law: `T = 2 pi sqrt(a^3 / mu)`, returned in minutes.
///
/// # Arguments
/// * `a` - The semi-major axis \[km\]
/// * `mu` - The standard gravitational parameter \[km^3 / s^2\]
///
/// # Returns
/// * `period` - The period in minutes \[min\]
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::calc_period;
///
/// // Define the semi-major axis and the standard gravitational parameter
/// let a = 6378.137; // [km]
/// let mu = 398600.5; // [km^3 / s^2]
///
/// // Calculate the period
/// let period = calc_period(a, mu);
/// assert!(period > 0.0);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
pub fn calc_period(a: f64, mu: f64) -> f64 {
    // Calculate time period in minutes
    let period = 2. * PI * (a.powi(3) / mu).sqrt() / 60.;

    return period;
}

// ----------
// Unit Tests
// ----------
