// Module for propagating TLEs with SGP4

// ------------------
// External Libraries
// ------------------
use std::f64::consts::PI;
use std::fs;

// ------------------
// Internal Libraries
// ------------------
use crate::time::{dayofyr2utc, utc2jday, DateTime, Timezone};
use crate::common::{Wgs, WGS72, deg2rad, calc_period, StateVector, CoordinateFrame};

// -------
// Structs
// -------

/// Simplified General Perturbations 4 (SGP4) parameters
///
/// This struct contains the epoch parameters which are necessary to propagate the state vectors of a satellite with SGP4
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone)]
pub struct Sgp4 {
    /// WGS model
    pub wgs: Wgs,

    /// General Perturbation Element Set
    pub gp: GenPerturbElementSet,

    /// Julian date at epoch \[days\]
    pub jd0: f64,

    /// Fractional Julian date at epoch \[days\]
    pub jdfrac0: f64,

    /// Deep space satellite 
    pub deep_space: bool,

    /// Brouwer mean elements at epoch
    pub brouwer0: BrouwerMeanElements,

    /// Atmospheric Drag Parameters
    pub atm_params: AtmDragParams,

    /// Earth Zonal Harmonics Parameters
    pub zonal_params: EarthZonalParams,

    /// Solar 3rd Body Parameters
    pub solar_params: ThirdBodyParams,

    /// Lunar 3rd Body Parameters
    pub lunar_params: ThirdBodyParams,

    /// Account for whole day resonance effects of Earth's gravity
    pub whole_day_resonance: bool,

    /// Whole day resonance parameters of Earth's gravity
    pub whole_day_resonance_params: WholeDayResonanceParams,

    /// Account for half day resonance effects of Earth's gravity
    pub half_day_resonance: bool,

    /// Half day resonance parameters of Earth's gravity
    pub half_day_resonance_params: HalfDayResonanceParams,
}

/// General Perturbation Element Set for an Earth-orbiting satellite.
///
/// This struct represents the parsed contents of a standard General Perturbation Element Set.
/// The General Perturbation Element Set is the standard set of orbital elements used with the
/// SGP4 propagator. Elements are commonly distributed as Two-Line Element (TLE) text; other
/// formats (OMM, GP JSON/CSV) use the same fields but are not parsed by this crate yet.
///
/// References:
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
#[derive(Default, Clone)]
pub struct GenPerturbElementSet {
    /// Common name of the satellite (e.g., "ISS (ZARYA)")
    pub common_name: String,

    /// NORAD satellite catalog number
    pub satellite_catalog_number: i32,

    /// Classification (`U` = Unclassified, `C` = Classified, `S` = Secret)
    pub classification: char,

    /// International designator (launch year, launch number, piece)
    pub international_designator: String,

    /// Epoch UTC datetime
    pub epoch_datetime: DateTime,

    /// First time derivative of mean motion \[revs/day^2\]
    pub first_derivative_of_mean_motion: f64,

    /// Second time derivative of mean motion \[revs/day^3\]
    pub second_derivative_of_mean_motion: f64,

    /// B* drag term \[1/Earth radii\]
    pub bstar: f64,

    /// Ephemeris type (always zero)
    pub ephemeris_type: i32,

    /// Element set number
    pub element_set_number: i32,

    /// Orbital inclination \[degrees\]
    pub inclination: f64,

    /// Right ascension of the ascending node (RAAN) \[degrees\]
    pub right_ascension_of_ascending_node: f64,

    /// Orbital eccentricity \[\]
    pub eccentricity: f64,

    /// Argument of perigee \[degrees\]
    pub argument_of_perigee: f64,

    /// Mean anomaly \[degrees\]
    pub mean_anomaly: f64,

    /// Mean motion \[revs/day\]
    pub mean_motion: f64,

    /// Revolution number at epoch \[revs\]
    pub revolution_number_at_epoch: i64,
}

/// Brouwer Mean Orbital Elements
///
/// This struct contains the mean orbital elements of a TLE converted to Brouwer convention. TLEs report mean orbital elements
/// in Kozai convention.
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct BrouwerMeanElements {
    /// Orbital inclination \[rad\]
    pub i: f64,

    /// The cosine of the orbital inclination
    pub theta: f64,

    /// Right ascension of the ascending node (RAAN) \[rad\]
    pub raan: f64,

    /// Orbital eccentricity \[\]
    pub e: f64,

    /// The square root of 1 minus the orbital eccentricity squared \[\]
    pub beta: f64,

    /// Argument of perigee \[rad\]
    pub omega: f64,

    /// Mean anomaly \[rad\]
    pub m: f64,

    /// Mean motion \[rad/min\]
    pub n: f64,

    /// Semi-major axis \[Earth Radii\]
    pub a: f64,

    /// The orbital period \[mins\]
    pub period: f64,
}

/// Atmospheric Drag Effects
///
/// This struct contains the parameters necessary to account for the impacts of atmospheric drag on an orbit.
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct AtmDragParams {
    /// Perigee height \[km\]
    pub hp: f64,

    /// q0 parameter of power-law density function \[Earth Radii\]
    pub q0: f64,

    /// s parameter of power-law density function \[Earth Radii\]
    pub s: f64,

    /// Zeta constant \[1 / Earth Radii\]
    pub zeta: f64,

    /// Eta constant \[\]
    pub eta: f64,

    /// C1 constant \[\]
    pub c1: f64,

    /// C3 constant \[\]
    pub c3: f64,

    /// C4 constant \[\]
    pub c4: f64,

    /// C5 constant \[\]
    pub c5: f64,

    /// D2 constant \[\]
    pub d2: f64,

    /// D3 constant \[\]
    pub d3: f64,

    /// D4 constant \[\]
    pub d4: f64,
}

/// Earth Zonal Harmonics
///
/// This struct contains the parameters necessary to account for the impacts of Earth's zonal harmonics on an orbit.
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct EarthZonalParams {
    /// Rate of change of mean anomaly \[rad / min\]
    pub m_dot: f64,

    /// Rate of change of the argument of perigee \[rad / min\]
    pub omega_dot: f64,

    /// Rate of change of the right ascension of the ascending node \[rad / min\]
    pub raan_dot: f64,
}

/// Solar and Lunar 3rd Body Effects
///
/// This struct contains the parameters necessary to account for the impacts of the Sun and Moon on an orbit.
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct ThirdBodyParams {
    /// Third body orbital inclination cosine \[\]
    pub cos_i: f64,

    /// Third body orbital inclination sine \[\]
    pub sin_i: f64,

    /// Third body eccentricity \[\]
    pub e: f64,

    /// Third body mean motion \[rad/min\]
    pub n: f64,

    /// Third body argument of perigee cosine \[\]
    pub cos_omega: f64,

    /// Third body argument of perigee sine \[\]
    pub sin_omega: f64,

    /// Third body right ascension of the ascending node (RAAN) \[rad\]
    pub raan: f64,

    /// Third body mean anomaly \[rad\]
    pub m: f64,

    /// The square root of 1 minus the orbital eccentricity squared \[\]
    pub beta: f64,

    /// Third body perturbation coefficient \[rad/min\]
    pub c: f64,

    /// x1 constant
    pub x1: f64,

    /// x2 constant
    pub x2: f64,

    /// x3 constant
    pub x3: f64,

    /// x4 constant
    pub x4: f64,

    /// x5 constant
    pub x5: f64,

    /// x6 constant 
    pub x6: f64,

    /// x7 constant
    pub x7: f64,

    /// x8 constant
    pub x8: f64,

    /// z1 constant
    pub z1: f64,

    /// z2 constant
    pub z2: f64,

    /// z3 constant
    pub z3: f64,

    /// z11 constant
    pub z11: f64,

    /// z13 constant
    pub z13: f64,

    /// z21 constant
    pub z21: f64,

    /// z23 constant
    pub z23: f64,

    /// z22 constant
    pub z22: f64,

    /// z12 constant
    pub z12: f64,

    /// z31 constant
    pub z31: f64,

    /// z32 constant
    pub z32: f64,

    /// z33 constant
    pub z33: f64,

    /// Rate of change of the orbital eccentricity \[1 / min\]
    pub e_dot: f64,

    /// Rate of change of the orbital inclination \[rad / min\]
    pub i_dot: f64,

    /// Rate of change of the mean anomaly \[rad / min\]
    pub m_dot: f64,

    /// Rate of change of the argument of perigee \[rad / min\]
    pub omega_dot: f64,

    /// Rate of change of the right ascension of the ascending node \[rad / min\]
    pub raan_dot: f64,
}

/// Half day resonance effects of Earth's gravity
///
/// This struct contains the parameters necessary to account for the impacts of half day resonance effects on an orbit.
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct HalfDayResonanceParams {
    /// Greenwich sidereal time at epoch \[rad\]
    pub theta_g: f64,

    /// lam0 constant
    pub lam0: f64,
    
    /// lam0 rate of change
    pub lam0_dot: f64,

    /// d2201 constant
    pub d2201: f64,

    /// d2211 constant
    pub d2211: f64,

    /// d3210 constant
    pub d3210: f64,

    /// d3222 constant
    pub d3222: f64,

    /// d5220 constant
    pub d5220: f64,

    /// d5232 constant
    pub d5232: f64,

    /// d4422 constant
    pub d4422: f64,

    /// d5421 constant
    pub d5421: f64,

    /// d5433 constant
    pub d5433: f64,

    /// d4410 constant
    pub d4410: f64,
}

/// Whole day resonance effects of Earth's gravity
///
/// This struct contains the parameters necessary to account for the impacts of whole day resonance effects on an orbit.
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct WholeDayResonanceParams {
    /// Greenwich sidereal time at epoch \[rad\]
    pub theta_g: f64,

    /// lam0 constant
    pub lam0: f64,

    /// lam0 rate of change
    pub lam0_dot: f64,

    /// lam31 constant
    pub lam31: f64,

    /// lam22 constant
    pub lam22: f64,

    /// lam33 constant
    pub lam33: f64,

    /// delta1 constant
    pub delta1: f64,

    /// delta2 constant
    pub delta2: f64,

    /// delta3 constant
    pub delta3: f64,
}

// -----
// Enums
// -----

// ---------
// Constants
// ---------

/// A conversion from rev/day to rad/min
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
const XPDOTP: f64 = 1440.0 / (2.0 * PI);

/// The rotational velocity of the earth in rad/min
///
/// References:
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
const RPTIM: f64 =  4.37526908801129966e-3;

// ---------
// Functions
// ---------

/// Builds a [`Sgp4`] struct from the lines of a Two-Line Element (TLE) set.
///
/// Given the two required TLE lines (line 1 and line 2), and an optional
/// name line (line 0), this function parses the input into a [`Sgp4`] struct.
///
/// # Arguments
/// * `line1` - The first TLE data line (TLE line 1)
/// * `line2` - The second TLE data line (TLE line 2)
/// * `line0` - Optional name line (TLE line 0)
///
/// # Returns
/// * [`Sgp4`] - Struct containing the parsed SGP4 parameters.
///
/// # Examples
/// ```rust
/// // Define the TLE lines
/// let tle_line0 = "ISS (ZARYA)";
/// let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let tle_line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// 
/// // Parse the TLE lines into a SGP4 struct
/// let sgp4 = from_tle_lines(tle_line1, tle_line2, Some(tle_line0));
/// 
/// // Assert the TLE struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
/// ```
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
pub fn from_tle_lines(line1: &str, line2: &str, line0: Option<&str>) -> Sgp4 {
    // Create mutable General Perturbation Element Set struct
    let mut gp = GenPerturbElementSet::default();

    // Validate the TLE checksum
    if !tle_checksum(line1) || !tle_checksum(line2) {
        eprintln!("warning: TLE checksum failed; continuing with parse");
    }

    // Extract the common name of the satellite from line 0
    if let Some(name_line) = line0 {
        if name_line.len() < 1 || name_line.len() > 24 {
            panic!("TLE line 0 is invalid: name must be 1-24 characters, got {}", name_line.len());
        }
        gp.common_name = name_line.to_string();
    }
    
    // Parse through line 1 and populate TLE struct
    if line1.len() < 69 || line1.len() > 69 {
        panic!("TLE line 1 is invalid: must be 69 characters, got {}", line1.len());
    } else {
        // Line 1
        // Satellite catalog number
        gp.satellite_catalog_number = line1[2..7].trim().parse::<i32>().unwrap();

        // Classification
        gp.classification = line1[7..8].trim().parse::<char>().unwrap();

        // International designator
        if line1[9..17].trim().is_empty() {
            // Handle case where international designator is not present
            gp.international_designator = "".to_string();
        } else {
            gp.international_designator = line1[9..17].trim().to_string();
        }

        // Epoch year (last two numbers)
        let yr_two_digit = line1[18..20].trim().parse::<i32>().unwrap();
        let epoch_year: i32;
        if yr_two_digit < 57 {
            epoch_year = 2000 + yr_two_digit
        } else {
            epoch_year = 1900 + yr_two_digit
        }

        // Epoch day of year
        let epoch_day = line1[20..32].trim().parse::<f64>().unwrap();

        // Epoch UTC datetime
        let Some(epoch_datetime) = dayofyr2utc(epoch_year, epoch_day).ok() else {
            panic!("Error converting epoch day of year to UTC datetime: Epoch year: {}, Epoch day: {}", epoch_year, epoch_day);
        };
        gp.epoch_datetime = epoch_datetime;

        // 1st derivative of mean motion [revs/day^2]
        gp.first_derivative_of_mean_motion = line1[33..43].trim().parse::<f64>().unwrap() * 2.0;

        // 2nd derivative of mean motion [revs/days^3]
        // Account for - in 2nd derivative of mean motion
        if line1[44..45].parse::<char>().unwrap() == '-' {
            gp.second_derivative_of_mean_motion = format!("-0.{}", line1[45..50].trim()).parse::<f64>().unwrap() * 10.0_f64.powi(line1[50..52].parse::<i32>().unwrap()) * 6.0_f64;
        } else {
            gp.second_derivative_of_mean_motion = format!("0.{}", line1[45..50].trim()).parse::<f64>().unwrap() * 10.0_f64.powi(line1[50..52].parse::<i32>().unwrap()) * 6.0_f64;
        }

        // B* [1/Earth Radii]
        // Account for - in B* term
        if line1[53..54].parse::<char>().unwrap() == '-' {
            gp.bstar = format!("-0.{}", line1[54..59].trim()).parse::<f64>().unwrap() * 10.0_f64.powi(line1[59..61].parse::<i32>().unwrap());
        } else {
            gp.bstar = format!("0.{}", line1[54..59].trim()).parse::<f64>().unwrap() * 10.0_f64.powi(line1[59..61].parse::<i32>().unwrap());
        }

        // Ephemeris type
        if line1[62..63].trim().is_empty() {
            // Handle case where ephemeris type is not present
            gp.ephemeris_type = 0;
        } else {
            gp.ephemeris_type = line1[62..63].parse::<i32>().unwrap();
        }

        // Element set number
        gp.element_set_number = line1[64..68].trim().parse::<i32>().unwrap();
    }

    // Parse through line 2 and populate TLE struct
    if line2.len() < 69 || line2.len() > 69 {
        panic!("TLE line 2 is invalid: must be 69 characters, got {}", line2.len());
    } else {
        // Line 2
        // Inclination [degs]
        gp.inclination = line2[8..16].trim().parse::<f64>().unwrap();

        // Right ascension of ascending node [degs]
        gp.right_ascension_of_ascending_node = line2[17..25].trim().parse::<f64>().unwrap();

        // Eccentricity
        gp.eccentricity = format!("0.{}", line2[26..33].trim()).parse::<f64>().unwrap();

        // Argument of perigee [degs]
        gp.argument_of_perigee = line2[34..42].trim().parse::<f64>().unwrap();

        // Mean anomaly [degs]
        gp.mean_anomaly = line2[43..51].trim().parse::<f64>().unwrap();

        // Mean motion [revs/day]
        gp.mean_motion = line2[52..63].trim().parse::<f64>().unwrap();

        // Revolution number at epoch
        gp.revolution_number_at_epoch = line2[63..68].trim().parse::<i64>().unwrap();
    }

    // Initialize the SGP4 parameters
    let sgp4 = init_sgp4(&gp, None);

    return sgp4;
}

/// Builds a vector of [`Sgp4`] structs from a string containing Two-Line Element (TLE) sets.
///
/// This function parses a string containing one or more TLEs in either
/// 2-line or 3-line (name + 2 lines) format and returns all successfully
/// parsed entries.
///
/// # Arguments
/// * `tle_string` - A string containing one or more Two-Line Element (TLE) sets
///
/// # Returns
/// * `Vec<Sgp4>` - A vector containing all successfully parsed SGP4 parameters
///
/// # Examples
/// ```rust
/// // Define the TLE string
/// let tle_string = "ISS (ZARYA)\n1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921\n2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// 
/// // Parse the TLE string into a SGP4 struct
/// let sgp4s = from_tle_string(tle_string);
/// let sgp4 = &sgp4s[0];
/// 
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
/// ```
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
pub fn from_tle_string(tle_string: &str) -> Vec<Sgp4> {
    // Parse the string into lines, removing spaces
    let lines: Vec<&str> = tle_string
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .collect();

    // Create the tles vector
    let mut sgp4s = Vec::new();
    let mut i = 0;

    // Iterate through the lines building TLE structs when possible
    while i < lines.len() {
        // Find TLEs within the string, either 2 or 3 line entries
        if lines[i].starts_with('1') {
            // This is likely a 2 line entry, check length
            if i + 1 >= lines.len() {
                break;
            }
            // Check that next line starts with '2'
            if lines[i + 1].starts_with('2') {
                let sgp4 = from_tle_lines(lines[i], lines[i + 1], None);
                sgp4s.push(sgp4);
                i += 2;
            } else {
                i += 1;
            }
        } else {
            // This is likely a 3 line entry, check length
            if i + 2 >= lines.len() {
                break;
            }
            // Check that next line 2 lines starts with '1' and '2'
            if lines[i + 1].starts_with('1') && lines[i + 2].starts_with('2') {
                let sgp4 = from_tle_lines(lines[i + 1], lines[i + 2], Some(lines[i]));
                sgp4s.push(sgp4);
                i += 3;
            } else {
                i += 1;
            }
        }
    }
    // Return vector of SGP4 structs
    return sgp4s;
}

/// Builds a vector of [`Sgp4`] structs from a file containing Two-Line Element (TLE) sets.
///
/// This function parses a file containing one or more TLEs in either
/// 2-line or 3-line (name + 2 lines) format and returns all successfully
/// parsed entries.
///
/// # Arguments
/// * `file_path` - A path to a file containing one or more Two-Line Element (TLE) sets
///
/// # Returns
/// * `Vec<Sgp4>` - A vector containing all successfully parsed SGP4 structs.
///
/// # Examples
/// ```rust
/// // Define the TLE file path
/// let tle_file_path = "test/tle_file.txt";
/// 
/// // Parse the TLE file into a vector of SGP4 structs
/// let sgp4s = from_tle_file(tle_file_path);
/// let sgp4 = &sgp4s[12];
/// 
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
/// ```
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
pub fn from_tle_file(file_path: &str) -> Vec<Sgp4> {
    // Open the TLE file
    let tle_string = fs::read_to_string(file_path)
        .expect("Cannot read TLE file");
    
    // Parse tle string into a vector of SGP4 structs
    let sgp4s = from_tle_string(&tle_string);

    // Return the vector of SGP4 structs
    return sgp4s;
}

/// Calculate the checksum of the Two-Line Element (TLE) line.
///
/// Given a TLE line, calculate the checksum of that line. Follow the following rules: 
/// - Ignore alpha characters
/// - Sum digits 0-9 as integer values
/// - '-' is treated as 1
/// - Return checksum % 10
///
/// # Arguments
/// * `line` - The TLE line to calculate the checksum of
///
/// # Panics
/// * If the TLE line is invalid (must be 69 characters)
///
/// # Returns
/// * `checksum` - The checksum of the TLE line (integer 0-9)
///
/// # Examples
/// ```rust
/// // Define the TLE line
/// let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// 
/// // Calculate the checksum of the TLE line
/// let checksum = calc_checksum(tle_line1);
/// 
/// // Assert the checksum is correct
/// assert_eq!(checksum, 1);
/// ```
pub fn calc_checksum(line: &str) -> i32 {
    // Initialize checksum to 0
    let mut checksum = 0;

    // Loop through the line and calculate the checksum
    for c in line.chars().take(68) {
        match c {
            '0'..='9' => checksum += (c as u8 - b'0') as i32,
            '-' => checksum += 1,
            _ => {}
        }
    }

    // Calculate the checksum
    checksum = checksum % 10;

    // Return the checksum
    return checksum;
}

/// Check if the TLE line has been corrupted by running a checksum test.
///
/// Given a TLE line, check if the checksum of that line is valid.
///
/// # Arguments
/// * `line` - The TLE line to check the checksum of
///
/// # Returns
/// * `bool` - True if the checksum of the line is valid, false if otherwise
///
/// # Examples
/// ```rust
/// // Define the TLE line
/// let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// 
/// // Calculate the checksum of the TLE line
/// let checksum = tle_checksum(tle_line1);
/// 
/// // Assert the checksum is correct
/// assert_eq!(checksum, true);
/// ```
pub fn tle_checksum(line: &str) -> bool {
    // Calculate the checksum of the line
    let checksum = calc_checksum(line);

    // Compare the checksum to the last character of the line
    if checksum == line[68..69].parse::<i32>().unwrap() {
        return true;
    } else {
        return false;
    }
}

/// Build an [`Sgp4`] struct for state propagation from a [`GenPerturbElementSet`] struct
///
/// Given a [`GenPerturbElementSet`] struct, calculate the time-independent parameters necessary 
/// to propagate a satellite's states in time. These parameters include
/// - Brouwer mean orbital elements
/// - Atmospheric drag parameters
/// - Earth zonal harmonics parameters
/// - Solar and Lunar 3rd body effects
/// - Resonance effects of Earth's gravity
///
/// # Arguments
/// * `gp` - The General Perturbation Element Set parameters
/// * `wgs` - Optional, specify World Geodetic System (WGS) parameters (defaults to WGS-72, the standard for TLEs)
///
/// # Returns
/// * [`Sgp4`] - The time-independent parameters for the SGP4 propagator
///
/// # Examples
/// ```rust
/// // Define General Perturbation (GP) Element Set
/// let gp = GenPerturbElementSet::default();
/// 
/// // Define WGS model
/// let wgs = WGS72;
///
/// // Initialize the SGP4 propagator
/// let sgp4 = init_sgp4(&gp, Some(&wgs));
/// ```
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn init_sgp4(gp: &GenPerturbElementSet, wgs: Option<&Wgs>) -> Sgp4 {
    // Use WGS72 or custom WGS models if provided
    let wgs_sgp4 = if let Some(wgs_passed) = wgs { *wgs_passed } else { WGS72 };

    // Extract General Perturbation (GP) Element Set contents in proper units
    let i0 = deg2rad(gp.inclination); // [rad]
    let n0_kozai = gp.mean_motion / XPDOTP; // [rad/min]
    let e0 = gp.eccentricity; // []
    let omega0 = deg2rad(gp.argument_of_perigee); // [rad]
    let raan0 = deg2rad(gp.right_ascension_of_ascending_node); // [rad]
    let m0 = deg2rad(gp.mean_anomaly); // [rad]

    // Extract GP epoch in Julian day format
    let (jd0, jdfrac0) = utc2jday(&gp.epoch_datetime).unwrap();

    // Recover Brouwer mean motion from Kozai mean motion (mean motion in GP)
    let theta0 = i0.cos();
    let beta0 = (1. - e0.powi(2)).sqrt();
    let a1 = (wgs_sgp4.ke / n0_kozai).powf(2./3.);
    let delta1 = (3./2.) * (wgs_sgp4.k2 / a1.powf(2.)) * (3. * i0.cos().powf(2.) - 1.) / (1. - e0.powf(2.)).powf(3./2.);
    let a2 = a1 * (1. - (1./3.) * delta1 - delta1.powf(2.) - (134./81.) * delta1.powf(3.));
    let delta0 = (3./2.) * (wgs_sgp4.k2 / a2.powf(2.)) * (3. * i0.cos().powf(2.) - 1.) / (1. - e0.powf(2.)).powf(3./2.);
    let n0 = n0_kozai / (1. + delta0); // [rad/min]
    let a0 = (wgs_sgp4.ke / n0).powf(2./3.); // [Earth radii]
    let a0_km = a0 * wgs_sgp4.r_earth_eq; // [km]
    let period0 = calc_period(a0_km, wgs_sgp4.mu); // [min]

    // Store Brouwer mean elements
    let brouwer0 = BrouwerMeanElements {
        i: i0,
        theta: theta0,
        raan: raan0,
        e: e0,
        beta: beta0,
        omega: omega0,
        m: m0,
        n: n0,
        a: a0,
        period: period0,
    };

    // Initialize atmospheric drag parameters
    let atm_params = init_atm_effects(&wgs_sgp4, gp, &brouwer0);

    // Initialize Earth zonal harmonics parameters
    let zonal_params = init_zonal_effects(&wgs_sgp4, &brouwer0);

    // Check for deep space satellite
    let mut deep_space = false;
    if period0 >= 225. {
        deep_space = true;
    }

    // Lunar and solar gravity effects
    let (lunar_params, solar_params) = init_lunar_solar_effects(deep_space, jd0, jdfrac0, &brouwer0);

    // Earth gravity resonance effects (use Vallado criteria instead of Hoots)
    let mut whole_day_resonance = false;
    let mut half_day_resonance = false;
    let mut whole_day_resonance_params = WholeDayResonanceParams::default();
    let mut half_day_resonance_params = HalfDayResonanceParams::default();
    if (n0 > 0.0034906585) && (n0 < 0.0052359877) {
        whole_day_resonance = true;
        whole_day_resonance_params = init_earth_gravity_resonance_wholeday(jd0, jdfrac0, &brouwer0, &zonal_params, &lunar_params, &solar_params);
    }
    if (n0 >= 8.26e-3) && (n0 <= 9.24e-3) && (e0 >= 0.5) {
        half_day_resonance = true;
        half_day_resonance_params = init_earth_gravity_resonance_halfday(jd0, jdfrac0, &brouwer0, &zonal_params, &lunar_params, &solar_params);
    }

    // Construct SGP4 propagator
    let sgp4 = Sgp4 {
        wgs: wgs_sgp4,
        gp: gp.clone(),
        jd0: jd0,
        jdfrac0: jdfrac0,
        deep_space: deep_space,
        brouwer0: brouwer0,
        atm_params: atm_params,
        zonal_params: zonal_params,
        lunar_params: lunar_params,
        solar_params: solar_params,
        whole_day_resonance: whole_day_resonance,
        whole_day_resonance_params: whole_day_resonance_params,
        half_day_resonance: half_day_resonance,
        half_day_resonance_params: half_day_resonance_params,
    };

    // Propagate to epoch so initialization failures surface through the same checks as propagation
    let _ = sgp4_prop_delta(&sgp4, 0.0);

    return sgp4;
}

/// Initialize the atmospheric drag effects
///
/// # Arguments
/// * `wgs` - The WGS model
/// * `gp` - The General Perturbation (GP) Element Set
/// * `brouwer0` - The Brouwer mean elements at epoch
///
/// # Returns
/// * `AtmDragParams` - The atmospheric drag parameters
///
/// # Examples
/// ```rust
/// // Define WGS model
/// let wgs = WGS72;
///
/// // Initialize the atmospheric drag effects
/// let atm_params = init_atm_effects(&wgs, &gp, &brouwer0);
/// ```
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn init_atm_effects(wgs: &Wgs, gp: &GenPerturbElementSet, brouwer0: &BrouwerMeanElements) -> AtmDragParams {
    // Define initial constants
    let a30 = -wgs.j3; // [Earth Radii^3]
    let q0 = (120. + wgs.r_earth_eq) / wgs.r_earth_eq; // [Earth radii]

    // Determine parameter s based on perigee height
    let rp = brouwer0.a * (1. - brouwer0.e); // Radius of perigee [Earth Radii]
    let hp = (rp - 1.) * wgs.r_earth_eq; // Perigee height [km]
    
    let mut s: f64; // [Earth radii]
    if hp >= 156. {
        s = (78. + wgs.r_earth_eq) / wgs.r_earth_eq;
    } else if hp >= 98.{
        s = (hp - 78. + wgs.r_earth_eq) / wgs.r_earth_eq; // [Earth radii]
    } else {
        s = (20. + wgs.r_earth_eq) / wgs.r_earth_eq; // [Earth radii]
    }

    // Calculate atmospheric drag parameters
    let zeta = 1. / (brouwer0.a - s);
    let eta = brouwer0.a * brouwer0.e * zeta;
    let psisq = (1. - eta.powi(2)).abs(); // abs is used to handle the case when eta > 1 (sub-orbital / decayed orbits)
    
    let c2_1 = (q0 - s).powi(4) * zeta.powi(4) * brouwer0.n * psisq.powf(-3.5);
    let c2_2 = brouwer0.a * (1. + (3./2.) * eta.powi(2) + 4. * brouwer0.e * eta + brouwer0.e * eta.powi(3));
    let c2_3 = (3./2.) * (wgs.k2 * zeta / psisq) * (-(1./2.) + (3./2.) * brouwer0.theta.powi(2)) * (8. + 24. * eta.powi(2) + 3. * eta.powi(4));
    let c2 = c2_1 * (c2_2 + c2_3);
    
    let c1 = gp.bstar * c2;
    // Vallado drop C3 when eccentricity is too small (avoids /e blow-up)
    let c3 = if brouwer0.e > 1.0e-4 {
        ((q0 - s).powf(4.) * zeta.powf(5.) * a30 * brouwer0.n * brouwer0.i.sin()) / (wgs.k2 * brouwer0.e)
    } else {
        0.0
    };
    
    let c4_1 = 2. * brouwer0.n * (q0 - s).powi(4) * zeta.powi(4) * brouwer0.a * brouwer0.beta.powi(2) * psisq.powf(-3.5);
    let c4_2 = 2. * eta * (1. + brouwer0.e*eta) + 0.5 * brouwer0.e + 0.5 * eta.powi(3);
    let c4_3 = 2. * wgs.k2 * zeta / (brouwer0.a * psisq);
    let c4_4 = 3. * (1. - 3. * brouwer0.theta.powi(2)) * (1. + 3./2. * eta.powi(2) - 2. * brouwer0.e * eta - 0.5 * brouwer0.e * eta.powi(3));
    let c4_5 = 3./4. * (1. - brouwer0.theta.powi(2)) * (2. * eta.powi(2) - brouwer0.e * eta - brouwer0.e * eta.powi(3)) * (2. * brouwer0.omega).cos();
    let c4 = c4_1 * (c4_2 - c4_3 * (c4_4 + c4_5));
    
    let c5_1 = 2. * (q0 - s).powi(4) * zeta.powi(4) * brouwer0.a * brouwer0.beta.powi(2) * psisq.powf(-3.5);
    let c5_2 = 1. + 11./4. * eta * (eta + brouwer0.e) + brouwer0.e * eta.powi(3);
    let c5 = c5_1 * c5_2;
    
    let d2 = 4. * brouwer0.a * zeta * c1.powi(2);
    let d3 = 4./3. * brouwer0.a * zeta.powi(2) * (17. * brouwer0.a + s) * c1.powi(3);
    let d4 = 2./3. * brouwer0.a.powi(2) * zeta.powi(3) * (221. * brouwer0.a + 31. * s) * c1.powi(4);

    // Store atmospheric drag parameters
    let atm_params = AtmDragParams {
        hp: hp,
        q0: q0,
        s: s,
        zeta: zeta,
        eta: eta,
        c1: c1,
        c3: c3,
        c4: c4,
        c5: c5,
        d2: d2,
        d3: d3,
        d4: d4,
    };

    return atm_params;
}

/// Initialize the Earth zonal harmonics effects
///
/// # Arguments
/// * `wgs` - The WGS model
/// * `brouwer0` - The Brouwer mean elements at epoch
///
/// # Returns
/// * `EarthZonalParams` - The Earth zonal harmonics parameters
///
/// # Examples
/// ```rust
/// // Define Brouwer mean elements at epoch
/// let brouwer0 = BrouwerMeanElements::default();
///
/// // Initialize the Earth zonal harmonics effects
/// let zonal_params = init_zonal_effects(&WGS72, &brouwer0);
/// ```
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn init_zonal_effects(wgs: &Wgs, brouwer0: &BrouwerMeanElements) -> EarthZonalParams {
    // Calculate orbital element rates of change due to zonal harmonics
    let m_dot_1 = 3. * wgs.k2 * (-1. + 3. * brouwer0.theta.powi(2)) / (2. * brouwer0.a.powi(2) * brouwer0.beta.powi(3));
    let m_dot_2 = 3. * wgs.k2.powi(2) * (13. - 78. * brouwer0.theta.powi(2) + 137. * brouwer0.theta.powi(4)) / (16. * brouwer0.a.powi(4) * brouwer0.beta.powi(7));
    let m_dot = (m_dot_1 + m_dot_2) * brouwer0.n;

    let omega_dot_1 = -3. * wgs.k2 * (1. - 5. * brouwer0.theta.powi(2)) / (2. * brouwer0.a.powi(2) * brouwer0.beta.powi(4));
    let omega_dot_2 = 3. * wgs.k2.powi(2) * (7. - 114. * brouwer0.theta.powi(2) + 395. * brouwer0.theta.powi(4)) / (16. * brouwer0.a.powi(4) * brouwer0.beta.powi(8));
    let omega_dot_3 = 5. * wgs.k4 * (3. - 36. * brouwer0.theta.powi(2) + 49. * brouwer0.theta.powi(4)) / (4. * brouwer0.a.powi(4) * brouwer0.beta.powi(8));
    let omega_dot = (omega_dot_1 + omega_dot_2 + omega_dot_3) * brouwer0.n;

    let raan_dot_1 = -3. * wgs.k2 * brouwer0.theta / (brouwer0.a.powi(2) * brouwer0.beta.powi(4));
    let raan_dot_2 = 3. * wgs.k2.powi(2) * (4. * brouwer0.theta - 19. * brouwer0.theta.powi(3)) / (2. * brouwer0.a.powi(4) * brouwer0.beta.powi(8));
    let raan_dot_3 = 5. * wgs.k4 * brouwer0.theta * (3. - 7. * brouwer0.theta.powi(2)) / (2. * brouwer0.a.powi(4) * brouwer0.beta.powi(8));
    let raan_dot = (raan_dot_1 + raan_dot_2 + raan_dot_3) * brouwer0.n;

    // Store Earth zonal parameters
    let zonal_params = EarthZonalParams {
        m_dot: m_dot,
        omega_dot: omega_dot,
        raan_dot: raan_dot
    };

    return zonal_params;
}

/// Initialize the Lunar and Solar third body effects
///
/// # Arguments
/// * `deep_space` - Is this a deep space satellite
/// * `jd0` - The Julian date at epoch \[days\]
/// * `jdfrac0` - The fractional Julian date at epoch \[days\]
/// * `brouwer0` - The Brouwer mean elements at epoch
///
/// # Returns
/// * `(ThirdBodyParams, ThirdBodyParams)` - The Lunar and Solar third body parameters
///
/// # Examples
/// ```rust
/// // Define Brouwer mean elements at epoch
/// let brouwer0 = BrouwerMeanElements::default();
///
/// // Initialize the Lunar and Solar third body effects
/// let (lunar_params, solar_params) = init_lunar_solar_effects(true, jd0, jdfrac0, &brouwer0);
/// ```
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn init_lunar_solar_effects(deep_space: bool, jd0: f64, jdfrac0: f64, brouwer0: &BrouwerMeanElements) -> (ThirdBodyParams, ThirdBodyParams) {
    // Check if the satellite is not in deep space
    if !deep_space {
        return (ThirdBodyParams::default(), ThirdBodyParams::default());
    }

    // Lunar/Solar element epochs (12/31/1899 12:00:00 UTC) \[Julian date\]
    let epoch_sm = 2415020.0;
    
    // Lunar constants
    // Lunar eccentricity
    let e_m = 0.05490;

    // Lunar mean motion \[rad/min\] (Spacetrack/Vallado digits; extra precision drifts periodics)
    let n_m = 1.5835218e-4;

    // Lunar perturbation coefficient \[rad/min\]
    let c_m = 4.7968065e-7;

    // Lunar right ascension of the ascending node (RAAN) with respect to the ecliptic plane at epoch \[rad\]
    let raan_me0 = 4.5236020;

    // Lunar right ascension of the ascending node (RAAN) with respect to the ecliptic plane time rate of change at epoch \[rad/day\]
    let raan_me0_dot = -9.2422029e-4;

    // Lunar longitude of perigee with respect to the ecliptic plane at epoch \[rad\]
    let u_me0 = 5.8351514;

    // Lunar longitude of perigee with respect to the ecliptic plane time rate of change at epoch \[rad/day\]
    let u_me0_dot = 0.0019443680;

    // Lunar mean anomaly at epoch \[rad\]
    let m_m0 = 4.7199672;

    // Lunar mean anomaly time rate of change at epoch \[rad/day\]
    let m_m0_dot = 0.22997150;

    // Solar constants
    // Solar inclination sin and cos
    let sin_i_s = 0.39785416;
    let cos_i_s = 0.91744867;

    // Solar eccentricity
    let e_s = 0.01675;

    // Solar mean motion \[rad/min\]
    let n_s = 1.19459e-5;

    // Solar right ascension of the ascending node (RAAN) \[rad\]
    let raan_s = 0.0;

    // Solar argument of periapsis cos and sin
    let sin_omega_s = -0.98088458;
    let cos_omega_s = 0.1945905;

    // Solar perturbation coefficient \[rad/min\]
    let c_s = 2.9864797e-6;

    // Solar mean anomaly at epoch \[rad\]
    let m_s0 = 6.2565837;

    // Solar mean anomaly time rate of change at epoch \[rad/day\]
    let m_s0_dot = 0.017201977;

    // Find the difference in time between the Solar / Lunar epoch and the TLE epoch
    let delta_t = jd0 + jdfrac0 - epoch_sm;

    // Calculate the Lunar RAAN wrt to the ecliptic plane at TLE epoch
    let raan_me = (raan_me0 + raan_me0_dot * delta_t).rem_euclid(2.0 * PI);
    let sin_raan_me = raan_me.sin();
    let cos_raan_me = raan_me.cos();

    // Lunar equatorial inclination (Spacetrack/Hoots legacy linearization of the
    // ecliptic→equatorial transform — not the exact acos form)
    let cos_i_m = 0.91375164 - 0.03568096 * cos_raan_me;
    let sin_i_m = (1.0 - cos_i_m * cos_i_m).sqrt();

    // Calculate the Lunar longitude of perigee referred to the ecliptic
    let gamma_m = u_me0 + u_me0_dot * delta_t;

    // Lunar ascending-node trig in the equatorial frame (legacy SDP4)
    let sin_raan_m = 0.089683511 * sin_raan_me / sin_i_m;
    let cos_raan_m = (1.0 - sin_raan_m * sin_raan_m).sqrt();
    let raan_m = sin_raan_m.atan2(cos_raan_m);

    // Lunar argument of periapsis (legacy SDP4 zx formulation)
    let mut zx = 0.39785416 * sin_raan_me / sin_i_m;
    let zy = cos_raan_m * cos_raan_me + 0.91744867 * sin_raan_m * sin_raan_me;
    zx = zx.atan2(zy);
    let omega_m = gamma_m + zx - raan_me;
    let sin_omega_m = omega_m.sin();
    let cos_omega_m = omega_m.cos();

    // Calculate the Lunar mean anomaly \[rad\]
    let m_m = (m_m0 + m_m0_dot * delta_t - gamma_m).rem_euclid(2.0 * PI);

    // Calculate the Solar mean anomaly \[rad\]
    let m_s = (m_s0 + m_s0_dot * delta_t).rem_euclid(2.0 * PI);

    // Calculate the Lunar secular rates
    let lunar_params = calc_lunar_solar_secular_rates(cos_i_m, sin_i_m, e_m, n_m, cos_omega_m, sin_omega_m, raan_m, m_m, c_m, brouwer0);

    // Calculate the Solar secular rates
    let solar_params = calc_lunar_solar_secular_rates(cos_i_s, sin_i_s, e_s, n_s, cos_omega_s, sin_omega_s, raan_s, m_s, c_s, brouwer0);

    return (lunar_params, solar_params);
}

/// Calculate the secular rates of a third body's orbital elements
///
/// # Arguments
/// * `cos_i_x` - Cosine of the third-body inclination \[\]
/// * `sin_i_x` - Sine of the third-body inclination \[\]
/// * `e_x` - Third-body eccentricity \[\]
/// * `n_x` - Third-body mean motion \[rad/min\]
/// * `cos_omega_x` - Cosine of the third-body argument of perigee \[\]
/// * `sin_omega_x` - Sine of the third-body argument of perigee \[\]
/// * `raan_x` - Third-body right ascension of the ascending node (RAAN) \[rad\]
/// * `m_x` - Third-body mean anomaly at the satellite epoch \[rad\]
/// * `c_x` - Third-body perturbation coefficient \[rad/min\]
/// * `brouwer0` - Satellite Brouwer mean elements at epoch
///
/// # Returns
/// * `ThirdBodyParams` - Secular rates and frozen geometric coefficients (`x*`, `z*`)
///
/// # Examples
/// ```rust
/// let brouwer0 = BrouwerMeanElements::default();
/// let lunar_params = calc_lunar_solar_secular_rates(
///     cos_i_m, sin_i_m, e_m, n_m, cos_omega_m, sin_omega_m, raan_m, m_m, c_m, &brouwer0,
/// );
/// ```
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn calc_lunar_solar_secular_rates(cos_i_x: f64, sin_i_x: f64, e_x: f64, n_x: f64, cos_omega_x: f64, sin_omega_x: f64, raan_x: f64, m_x: f64, c_x: f64, brouwer0: &BrouwerMeanElements) -> ThirdBodyParams {
    // Precompute common quantities
    let cos_raan_diff = (brouwer0.raan - raan_x).cos();
    let sin_raan_diff = (brouwer0.raan - raan_x).sin();
    let cos_omega0 = brouwer0.omega.cos();
    let sin_omega0 = brouwer0.omega.sin();
    let cos_i0 = brouwer0.i.cos();
    let sin_i0 = brouwer0.i.sin();
    let beta_x = (1. - e_x.powi(2)).sqrt();

    // Calculate 3rd body constants
    let a1 = cos_omega_x * cos_raan_diff + sin_omega_x * cos_i_x * sin_raan_diff;
    let a3 = -sin_omega_x * cos_raan_diff + cos_omega_x * cos_i_x * sin_raan_diff;
    let a7 = -cos_omega_x * sin_raan_diff + sin_omega_x * cos_i_x * cos_raan_diff;
    let a8 = sin_omega_x * sin_i_x;
    let a9 = sin_omega_x * sin_raan_diff + cos_omega_x * cos_i_x * cos_raan_diff;
    let a10 = cos_omega_x * sin_i_x;
    let a2 = a7 * cos_i0 + a8 * sin_i0;
    let a4 = a9 * cos_i0 + a10 * sin_i0;
    let a5 = -a7 * sin_i0 + a8 * cos_i0;
    let a6 = -a9 * sin_i0 + a10 * cos_i0;

    let x1 = a1 * cos_omega0 + a2 * sin_omega0;
    let x2 = a3 * cos_omega0 + a4 * sin_omega0;
    let x3 = -a1 * sin_omega0 + a2 * cos_omega0;
    let x4 = -a3 * sin_omega0 + a4 * cos_omega0;
    let x5 = a5 * sin_omega0;
    let x6 = a6 * sin_omega0;
    let x7 = a5 * cos_omega0;
    let x8 = a6 * cos_omega0;
    
    let z31 = 12. * x1.powi(2) - 3. * x3.powi(2);
    let z32 = 24. * x1 * x2 - 6. * x3 * x4;
    let z33 = 12. * x2.powi(2) - 3. * x4.powi(2);
    let z1 = 6. * (a1.powi(2) + a2.powi(2)) + (1. + brouwer0.e.powi(2)) * z31;
    let z2 = 12. * (a1 * a3 + a2 * a4) + (1. + brouwer0.e.powi(2)) * z32;
    let z3 = 6. * (a3.powi(2) + a4.powi(2)) + (1. + brouwer0.e.powi(2)) * z33;
    let z11 = -6. * a1 * a5 + brouwer0.e.powi(2) * (-24. * x1 * x7 - 6. * x3 * x5);
    let z13 = -6. * a3 * a6 + brouwer0.e.powi(2) * (-24. * x2 * x8 - 6. * x4 * x6);
    let z21 = 6. * a2 * a5 + brouwer0.e.powi(2) * (24. * x1 * x5 - 6. * x3 * x7);
    let z23 = 6. * a4 * a6 + brouwer0.e.powi(2) * (24. * x2 * x6 - 6. * x4 * x8);
    let z22 = 6. * a4 * a5 + 6. * a2 * a6 + brouwer0.e.powi(2) * (24. * x2 * x5 + 24. * x1 * x6 - 6. * x4 * x7 - 6. * x3 * x8);
    let z12 = -6. * a1 * a6 - 6. * a3 * a5 - brouwer0.e.powi(2) * (24. * x2 * x7 + 24. * x1 * x8 + 6. * x3 * x6 + 6. * x4 * x5);

    // Calculate secular rates
    let e_x_dot = -15. * c_x * n_x * (brouwer0.e * brouwer0.beta / brouwer0.n) * (x1 * x3 + x2 * x4);
    
    let i_x_dot = (-c_x * n_x / (2. * brouwer0.n * brouwer0.beta)) * (z11 + z13);
    
    let m_x_dot = (-c_x * n_x / brouwer0.n) * (z1 + z3 - 14. - 6. * brouwer0.e.powi(2));
    
    let mut raan_x_dot = 0.;
    if brouwer0.i >= deg2rad(3.) {
        raan_x_dot = c_x * n_x / (2. * brouwer0.n * brouwer0.beta * sin_i0) * (z21 + z23);
    }

    let mut omega_x_dot = c_x * n_x * brouwer0.beta / brouwer0.n * (z31 + z33 - 6.);
    if brouwer0.i >= deg2rad(3.) {
        omega_x_dot = omega_x_dot - raan_x_dot * cos_i0;
    }

    // Store the 3rd body parameters
    let third_body_params = ThirdBodyParams {
        cos_i: cos_i_x,
        sin_i: sin_i_x,
        e: e_x,
        n: n_x,
        cos_omega: cos_omega_x,
        sin_omega: sin_omega_x,
        raan: raan_x,
        m: m_x,
        beta: beta_x,
        c: c_x,
        x1: x1,
        x2: x2,
        x3: x3,
        x4: x4,
        x5: x5,
        x6: x6,
        x7: x7,
        x8: x8,
        z1: z1,
        z2: z2,
        z3: z3,
        z11: z11,
        z13: z13,
        z21: z21,
        z23: z23,
        z22: z22,
        z12: z12,
        z31: z31,
        z32: z32,
        z33: z33,
        e_dot: e_x_dot,
        i_dot: i_x_dot,
        m_dot: m_x_dot,
        raan_dot: raan_x_dot,
        omega_dot: omega_x_dot,
    };

    // Return 3rd body parameters
    return third_body_params;
}

/// Initialize the half day resonance effects of Earth's gravity
///
/// # Arguments
/// * `jd0` - The Julian date at epoch \[days\]
/// * `jdfrac0` - The fractional Julian date at epoch \[days\]
/// * `brouwer0` - The Brouwer mean elements at epoch
/// * `zonal_params` - The Earth zonal harmonics parameters
/// * `lunar_params` - The Lunar third body parameters
/// * `solar_params` - The Solar third body parameters
///
/// # Returns
/// * `HalfDayResonanceParams` - The half day resonance parameters
///
/// # Examples
/// ```rust
/// // Define Brouwer mean elements at epoch
/// let brouwer0 = BrouwerMeanElements::default();
///
/// // Define Earth zonal harmonics parameters
/// let zonal_params = EarthZonalParams::default();
///
/// // Define Lunar third body parameters
/// let lunar_params = ThirdBodyParams::default();
///
/// // Define Solar third body parameters
/// let solar_params = ThirdBodyParams::default();
///
/// // Define Julian date at epoch
/// let jd0 = 2451545.0;
///
/// // Define fractional Julian date at epoch
/// let jdfrac0 = 0.0;
///
/// // Initialize half day resonance effects
/// let half_day_resonance_params = init_earth_gravity_resonance_halfday(jd0, jdfrac0, &brouwer0, &zonal_params, &lunar_params, &solar_params);
/// ```
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn init_earth_gravity_resonance_halfday(
    jd0: f64, 
    jdfrac0: f64, 
    brouwer0: &BrouwerMeanElements, 
    zonal_params: &EarthZonalParams, 
    lunar_params: &ThirdBodyParams, 
    solar_params: &ThirdBodyParams
) -> HalfDayResonanceParams {
    // Precompute common quantities
    let cos_i0 = brouwer0.i.cos();
    let sin_i0 = brouwer0.i.sin();

    // Define constants
    let c22s22 = 1.7891679e-6;
    let c32s32 = 3.7393792e-7;
    let c44s44 = 7.3636953e-9;
    let c52s52 = 1.1428639e-7;
    let c54s54 = 2.1765803e-9;

    // Calculate functions of inclination
    let f220 = (3./4.) * (1. + cos_i0).powi(2);
    let f221 = (3./2.) * sin_i0.powi(2);
    let f321 = (15./8.) * sin_i0 * (1. - 2. * cos_i0 - 3. * cos_i0.powi(2));
    let f322 = (-15./8.) * sin_i0 * (1. + 2. * cos_i0 - 3. * cos_i0.powi(2));
    let f441 = (105./4.) * sin_i0.powi(2) * (1. + cos_i0).powi(2);
    let f442 = (315./8.) * sin_i0.powi(4);
    let f522 = (315./32.) * (sin_i0.powi(3) - 2. * sin_i0.powi(3) * cos_i0 - 5. * sin_i0.powi(3) * cos_i0.powi(2) + sin_i0 * ((-2./3.) + (4./3.) * cos_i0 + 2. * cos_i0.powi(2)));
    let f523 = (105./16.) * sin_i0 * (1. + 2. * cos_i0 - 3. * cos_i0.powi(2) - (3./2.) * sin_i0.powi(2) * (1. + 2. * cos_i0 - 5. * cos_i0.powi(2)));
    let f542 = (945./32.) * sin_i0 * (2. - 8. * cos_i0 + cos_i0.powi(2) * (-12. + 8. * cos_i0 + 10. * cos_i0.powi(2)));
    let f543 = (945./32.) * sin_i0 * (cos_i0.powi(2) * (12. + 8. * cos_i0 - 10. * cos_i0.powi(2)) - 2. - 8. * cos_i0);

    // Calculate functions of eccentricity
    let g211: f64;
    let g201 = -0.306 - 0.44 * (brouwer0.e - 0.64);
    let g310: f64;
    let g322: f64;
    let g410: f64;
    let g422: f64;
    let g520: f64;
    let g521: f64;
    let g532: f64;
    let g533: f64;
    if brouwer0.e <= 0.65 {
        g211 = 3.616 - 13.247 * brouwer0.e + 16.29 * brouwer0.e.powi(2);
        g310 = -19.302 + 117.39 * brouwer0.e - 228.419 * brouwer0.e.powi(2) + 156.591 * brouwer0.e.powi(3);
        g322 = -18.9068 + 109.7927 * brouwer0.e - 214.6334 * brouwer0.e.powi(2) + 146.5816 * brouwer0.e.powi(3);
        g410 = -41.122 + 242.694 * brouwer0.e - 471.094 * brouwer0.e.powi(2) + 313.953 * brouwer0.e.powi(3);
        g422 = -146.407 + 841.88 * brouwer0.e - 1629.014 * brouwer0.e.powi(2) + 1083.435 * brouwer0.e.powi(3);
        g520 = -532.114 + 3017.977 * brouwer0.e - 5740.032 * brouwer0.e.powi(2) + 3708.276 * brouwer0.e.powi(3);
    } else {
        g211 = -72.099 + 331.819 * brouwer0.e - 508.738 * brouwer0.e.powi(2) + 266.724 * brouwer0.e.powi(3);
        g310 = -346.844 + 1582.851 * brouwer0.e - 2415.925 * brouwer0.e.powi(2) + 1246.113 * brouwer0.e.powi(3);
        g322 = -342.585 + 1554.908 * brouwer0.e - 2366.899 * brouwer0.e.powi(2) + 1215.972 * brouwer0.e.powi(3);
        g410 = -1052.797 + 4758.686 * brouwer0.e - 7193.992 * brouwer0.e.powi(2) + 3651.957 * brouwer0.e.powi(3);
        g422 = -3581.69 + 16178.11 * brouwer0.e - 24462.77 * brouwer0.e.powi(2) + 12422.52 * brouwer0.e.powi(3);
        if brouwer0.e < 0.715 {
            g520 = 1464.74 - 4664.75 * brouwer0.e + 3763.64 * brouwer0.e.powi(2);
        } else {
            g520 = -5149.66 + 29936.92 * brouwer0.e - 54087.36 * brouwer0.e.powi(2) + 31324.56 * brouwer0.e.powi(3);
        }
    }
    if brouwer0.e < 0.7 {
        g521 = -822.71072 + 4568.6173 * brouwer0.e - 8491.4146 * brouwer0.e.powi(2) + 5337.524 * brouwer0.e.powi(3);
        g532 = -853.666 + 4690.25 * brouwer0.e - 8624.77 * brouwer0.e.powi(2) + 5341.4 * brouwer0.e.powi(3);
        g533 = -919.2277 + 4988.61 * brouwer0.e - 9064.77 * brouwer0.e.powi(2) + 5542.21 * brouwer0.e.powi(3);
    } else {
        g521 = -51752.104 + 218913.95 * brouwer0.e - 309468.16 * brouwer0.e.powi(2) + 146349.42 * brouwer0.e.powi(3);
        g532 = -40023.88 + 170470.89 * brouwer0.e - 242699.48 * brouwer0.e.powi(2) + 115605.82 * brouwer0.e.powi(3);
        g533 = -37995.78 + 161616.52 * brouwer0.e - 229838.2 * brouwer0.e.powi(2) + 109377.94 * brouwer0.e.powi(3);
    }

    // Calculate the quadruples
    let d2201 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(2) * (c22s22 * f220 * g201);
    let d2211 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(2) * (c22s22 * f221 * g211);
    let d3210 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(3) * (c32s32 * f321 * g310);
    let d3222 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(3) * (c32s32 * f322 * g322);
    let d5220 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(5) * (c52s52 * f522 * g520);
    let d5232 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(5) * (c52s52 * f523 * g532);
    let d4422 = 6. * brouwer0.n.powi(2) / brouwer0.a.powi(4) * (c44s44 * f442 * g422); // 2x typo in Hoots et al 2004
    let d5421 = 6. * brouwer0.n.powi(2) / brouwer0.a.powi(5) * (c54s54 * f542 * g521); // 2x typo in Hoots et al 2004
    let d5433 = 6. * brouwer0.n.powi(2) / brouwer0.a.powi(5) * (c54s54 * f543 * g533); // Typo in Hoots et al 2004
    let d4410 = 6. * brouwer0.n.powi(2) / brouwer0.a.powi(4) * (c44s44 * f441 * g410); // 2x typo in Hoots et al 2004

    // Calculate the initial value for the auxilary variable lam0
    let theta_g = calc_theta_g(jd0, jdfrac0);
    let lam0 = (brouwer0.m + 2. * brouwer0.raan - 2. * theta_g).rem_euclid(2.0 * PI);
    let lam0_dot = zonal_params.m_dot + (lunar_params.m_dot + solar_params.m_dot) + 2. * zonal_params.raan_dot + 2. * (lunar_params.raan_dot + solar_params.raan_dot) - 2. * RPTIM;

    // Store resonance parameters
    let half_day_resonance_params = HalfDayResonanceParams {
        theta_g: theta_g,
        lam0: lam0,
        lam0_dot: lam0_dot,
        d2201: d2201,
        d2211: d2211,
        d3210: d3210,
        d3222: d3222,
        d5220: d5220,
        d5232: d5232,
        d4422: d4422,
        d5421: d5421,
        d5433: d5433,
        d4410: d4410,
    };

    return half_day_resonance_params;
}

/// Initialize the whole day resonance effects of Earth's gravity
///
/// # Arguments
/// * `jd0` - The Julian date at epoch \[days\]
/// * `jdfrac0` - The fractional Julian date at epoch \[days\]
/// * `brouwer0` - The Brouwer mean elements at epoch
/// * `zonal_params` - The Earth zonal harmonics parameters
/// * `lunar_params` - The Lunar third body parameters
/// * `solar_params` - The Solar third body parameters
///
/// # Returns
/// * `WholeDayResonanceParams` - The whole day resonance parameters
///
/// # Examples
/// ```rust
/// // Define Brouwer mean elements at epoch
/// let brouwer0 = BrouwerMeanElements::default();
///
/// // Define Earth zonal harmonics parameters
/// let zonal_params = EarthZonalParams::default();
///
/// // Define Lunar third body parameters
/// let lunar_params = ThirdBodyParams::default();
///
/// // Define Solar third body parameters
/// let solar_params = ThirdBodyParams::default();
///
/// // Define Julian date at epoch
/// let jd0 = 2451545.0;
///
/// // Define fractional Julian date at epoch
/// let jdfrac0 = 0.0;
///
/// // Initialize whole day resonance effects
/// let whole_day_resonance_params = init_earth_gravity_resonance_wholeday(jd0, jdfrac0, &brouwer0, &zonal_params, &lunar_params, &solar_params);
/// ```
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn init_earth_gravity_resonance_wholeday(
    jd0: f64, 
    jdfrac0: f64, 
    brouwer0: &BrouwerMeanElements, 
    zonal_params: &EarthZonalParams, 
    lunar_params: &ThirdBodyParams, 
    solar_params: &ThirdBodyParams
) -> WholeDayResonanceParams {
    // Precompute common quantities
    let cos_i0 = brouwer0.i.cos();
    let sin_i0 = brouwer0.i.sin();

    // Define constants
    let q31 = 2.1460748e-6;
    let q22 = 1.7891679e-6;
    let q33 = 2.2123015e-7;
    let lam31 = 0.13130908;
    let lam22 = 2.88431980;
    let lam33 = 0.37448087;

    // Calculate functions of inclination
    let f220 = (3./4.) * (1. + cos_i0).powi(2);
    let f311 = (15./16.) * sin_i0.powi(2) * (1. + 3. * cos_i0) - (3./4.) * (1. + cos_i0);
    let f330 = (15. / 8.) * (1. + cos_i0).powi(3);
    
    // Calculate functions of eccentricity
    let g200 = 1. - (5./2.) * brouwer0.e.powi(2) + (13. / 16.) * brouwer0.e.powi(4);
    let g310 = 1. + 2. * brouwer0.e.powi(2);
    let g300 = 1. - 6. * brouwer0.e.powi(2) + (423. / 64.) * brouwer0.e.powi(4);

    // Calculate coefficients of the resonance terms
    let delta1 = (3. * brouwer0.n.powi(2) / brouwer0.a.powi(3)) * f311 * g310 * q31;
    let delta2 = (6. * brouwer0.n.powi(2) / brouwer0.a.powi(2)) * f220 * g200 * q22;
    let delta3 = (9. * brouwer0.n.powi(2) / brouwer0.a.powi(3)) * f330 * g300 * q33;

    // Calculate the initial value for the auxilary variable lam0
    let theta_g = calc_theta_g(jd0, jdfrac0);
    let lam0 = brouwer0.m + brouwer0.raan + brouwer0.omega - theta_g;
    let lam0_dot_1 = zonal_params.m_dot + (lunar_params.m_dot + solar_params.m_dot) + zonal_params.raan_dot + (lunar_params.raan_dot + solar_params.raan_dot);
    let lam0_dot_2 = zonal_params.omega_dot + (lunar_params.omega_dot + solar_params.omega_dot) - RPTIM;
    let lam0_dot = lam0_dot_1 + lam0_dot_2;

    // Store resonance parameters
    let whole_day_resonance_params = WholeDayResonanceParams {
        theta_g: theta_g,
        lam0: lam0,
        lam0_dot:lam0_dot,
        lam31: lam31,
        lam22: lam22,
        lam33: lam33,
        delta1: delta1,
        delta2: delta2,
        delta3: delta3,
    };

    return whole_day_resonance_params;
}

/// Calculate Greenwich mean sidereal time (GMST) / longitude of Greenwich at a Julian date.
///
/// Used as θ_g when initializing 12 h / 24 h resonance terms.
///
/// # Arguments
/// * `jd0` - Julian day (integer part) \[days\]
/// * `jdfrac0` - Julian day fraction \[days\]
///
/// # Returns
/// * `theta_g` - GMST \[rad\], wrapped to \[0, 2π)
///
/// # Examples
/// ```rust
/// let theta_g = calc_theta_g(jd0, jdfrac0);
/// ```
///
/// References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
pub fn calc_theta_g(jd0: f64, jdfrac0: f64) -> f64 {
    // Calculate the Julian centuries since J2000.0
    let tut1 = (jd0 + jdfrac0 - 2451545.0) / 36525.0; // [centuries]

    // Calculate the Greenwich sidereal time in seconds
    let temp = -6.2e-6 * tut1.powi(3) + 0.093104 * tut1.powi(2) + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;  // [seconds]

    // Calculate the Greenwich sidereal time in radians
    let theta_g = (deg2rad(temp / 240.0)).rem_euclid(2.0 * PI); // [radians] 360/86400 = 1/240 degrees per second

    // Return the Greenwich sidereal time in radians
    return theta_g;
}

/// Propagate an initialized [`Sgp4`] model to a UTC [`DateTime`].
///
/// Converts `datetime` to Julian date, forms minutes since the TLE epoch, then
/// calls [`sgp4_prop_delta`]. Panics if `datetime` is not UTC or Julian conversion fails.
///
/// # Arguments
/// * `sgp4` - Initialized propagator
/// * `datetime` - Propagation epoch in UTC
///
/// # Returns
/// * `StateVector` - TEME position \[km\] and velocity \[km/s\]
///
/// # Examples
/// ```rust
/// let state_vector = sgp4_prop_datetime(&sgp4, &datetime);
/// ```
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn sgp4_prop_datetime(sgp4: &Sgp4, datetime: &DateTime) -> StateVector{
    // Convert datetime to Julian day format
    let (jd_prop, jdfrac_prop) = utc2jday(&datetime).unwrap();

    // Get minutes since epoch
    let delta_t = (jd_prop + jdfrac_prop - (sgp4.jd0 + sgp4.jdfrac0)) * 1440.;

    // Propagate the state vector
    return sgp4_prop_delta(sgp4, delta_t);
}

/// Propagate an initialized [`Sgp4`] model by `delta_t` minutes from the TLE epoch.
///
/// Applies secular drag/zonal updates (and deep-space lunar/solar + resonance when
/// applicable), long-period periodics, then the short-period Kepler / J₂ solution.
/// Panics on non-physical intermediate elements (e.g. mean motion ≤ 0, eccentricity
/// out of range, decayed radius), matching Vallado-style checks.
///
/// # Arguments
/// * `sgp4` - Initialized propagator
/// * `delta_t` - Minutes since epoch (may be negative)
///
/// # Returns
/// * `StateVector` - TEME position \[km\] and velocity \[km/s\]
///
/// # Examples
/// ```rust
/// let state_vector = sgp4_prop_delta(&sgp4, 360.0); // +6 hours
/// ```
///
/// References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn sgp4_prop_delta(sgp4: &Sgp4, delta_t: f64) -> StateVector{
    // Create mutable variables for the orbital elements
    let mut m: f64;
    let mut omega: f64;
    let mut raan: f64;
    let mut e = sgp4.brouwer0.e;
    let mut i = sgp4.brouwer0.i;
    let mut n = sgp4.brouwer0.n;

    // Account for Earth zonal gravity and partial atmospheric drag effects
    let m_df = sgp4.brouwer0.m + sgp4.brouwer0.n * delta_t + sgp4.zonal_params.m_dot * delta_t;
    let omega_df = sgp4.brouwer0.omega + sgp4.zonal_params.omega_dot * delta_t;
    let raan_df = sgp4.brouwer0.raan + sgp4.zonal_params.raan_dot * delta_t;

    // Neglect delta_omega and delta_m if deep space or perigee height is less than 220 km, or e ≤ 1e-4
    let delta_omega: f64;
    let delta_m: f64;
    if sgp4.deep_space || sgp4.atm_params.hp < 220. || sgp4.brouwer0.e <= 1.0e-4 {
        delta_omega = 0.;
        delta_m = 0.;
    } else {
        delta_omega = sgp4.gp.bstar * sgp4.atm_params.c3 * sgp4.brouwer0.omega.cos() * delta_t;
        delta_m = (-2./3.) * (sgp4.atm_params.q0 - sgp4.atm_params.s).powi(4) * sgp4.gp.bstar * sgp4.atm_params.zeta.powi(4) * (1. / (sgp4.brouwer0.e * sgp4.atm_params.eta)) * ((1. + sgp4.atm_params.eta*m_df.cos()).powi(3) - (1. + sgp4.atm_params.eta*sgp4.brouwer0.m.cos()).powi(3));
    }

    m = m_df + delta_omega + delta_m;
    omega = omega_df - delta_omega - delta_m;
    raan = raan_df - (21./2.) * (sgp4.brouwer0.n * sgp4.wgs.k2 * sgp4.brouwer0.theta / (sgp4.brouwer0.a.powi(2) * sgp4.brouwer0.beta.powi(2))) * sgp4.atm_params.c1 * delta_t.powi(2);

    // Account for Lunar and Solar third body effects
    if sgp4.deep_space {
        m += (sgp4.lunar_params.m_dot + sgp4.solar_params.m_dot) * delta_t;
        omega += (sgp4.lunar_params.omega_dot + sgp4.solar_params.omega_dot) * delta_t;
        raan += (sgp4.lunar_params.raan_dot + sgp4.solar_params.raan_dot) * delta_t;
        e += (sgp4.lunar_params.e_dot + sgp4.solar_params.e_dot) * delta_t;
        i += (sgp4.lunar_params.i_dot + sgp4.solar_params.i_dot) * delta_t;
    }

    // Account for the whole and half day resonance effects of Earth's gravity
    // Vallado dspace: ±720 min Euler-Maclaurin steps (works for negative tsince)
    if sgp4.half_day_resonance {
        let mut lami = sgp4.half_day_resonance_params.lam0;
        let mut ni = sgp4.brouwer0.n;
        let mut lami_dot: f64;
        let mut ni_dot: f64;
        let mut lami_ddot: f64;
        let mut ni_ddot: f64;
        let step = if delta_t >= 0.0 { 720.0 } else { -720.0 };
        let em_steps = (delta_t / step).floor() as i32;
        let t_em = delta_t - em_steps as f64 * step;

        (lami_dot, ni_dot, lami_ddot, ni_ddot) = half_day_euler_maclaurin_step(
            lami,
            ni,
            sgp4.brouwer0.omega,
            &sgp4.half_day_resonance_params,
        );

        for em_step in 0..em_steps {
            // h²/2 term uses |720|²; linear term uses signed step (Vallado step2 = 259200)
            lami += lami_dot * step + 0.5 * lami_ddot * 518400.;
            ni += ni_dot * step + 0.5 * ni_ddot * 518400.;

            let omegai = sgp4.brouwer0.omega
                + sgp4.zonal_params.omega_dot * (em_step + 1) as f64 * step;
            (lami_dot, ni_dot, lami_ddot, ni_ddot) = half_day_euler_maclaurin_step(
                lami,
                ni,
                omegai,
                &sgp4.half_day_resonance_params,
            );
        }

        lami = lami + (lami_dot * t_em) + (0.5 * lami_ddot * t_em.powi(2));
        ni = ni + (ni_dot * t_em) + (0.5 * ni_ddot * t_em.powi(2));

        let theta_t = (sgp4.half_day_resonance_params.theta_g + RPTIM * delta_t).rem_euclid(2.0 * PI);
        n = ni;
        m = lami - 2. * raan + 2. * theta_t;
    } else if sgp4.whole_day_resonance {
        let mut lami = sgp4.whole_day_resonance_params.lam0;
        let mut ni = sgp4.brouwer0.n;
        let mut lami_dot: f64;
        let mut ni_dot: f64;
        let mut lami_ddot: f64;
        let mut ni_ddot: f64;
        let step = if delta_t >= 0.0 { 720.0 } else { -720.0 };
        let em_steps = (delta_t / step).floor() as i32;
        let t_em = delta_t - em_steps as f64 * step;

        (lami_dot, ni_dot, lami_ddot, ni_ddot) =
            whole_day_euler_maclaurin_step(lami, ni, &sgp4.whole_day_resonance_params);

        for _ in 0..em_steps {
            lami += lami_dot * step + 0.5 * lami_ddot * 518400.;
            ni += ni_dot * step + 0.5 * ni_ddot * 518400.;

            (lami_dot, ni_dot, lami_ddot, ni_ddot) =
                whole_day_euler_maclaurin_step(lami, ni, &sgp4.whole_day_resonance_params);
        }

        lami = lami + (lami_dot * t_em) + (0.5 * lami_ddot * t_em.powi(2));
        ni = ni + (ni_dot * t_em) + (0.5 * ni_ddot * t_em.powi(2));

        let theta_t = (sgp4.whole_day_resonance_params.theta_g + RPTIM * delta_t).rem_euclid(2.0 * PI);
        n = ni;
        m = lami - raan - omega + theta_t;
    }

    // Mean motion must remain positive before the drag semi-major-axis update
    if n <= 0.0 {
        panic!("mean motion {} is less than or equal to zero", n);
    }

    // Account for remaining atmospheric drag effects
    let a: f64;
    let il_atm: f64;
    if sgp4.deep_space || sgp4.atm_params.hp < 220. {
        e += -sgp4.gp.bstar * (sgp4.atm_params.c4 * delta_t);
        let a_1 = 1. - sgp4.atm_params.c1 * delta_t; // Drop quadratic term, different from Hoots et al 2004
        a = (sgp4.wgs.ke / n).powf(2./3.) * a_1.powi(2);
        n = sgp4.wgs.ke / a.powf(1.5);
        let il_1 = 3./2. * sgp4.atm_params.c1 * delta_t.powi(2);
        il_atm = sgp4.brouwer0.n * (il_1);
    } else {
        e += -sgp4.gp.bstar * (sgp4.atm_params.c4 * delta_t + sgp4.atm_params.c5 * (m.sin() - sgp4.brouwer0.m.sin()));
        let a_1 = 1. - sgp4.atm_params.c1 * delta_t - sgp4.atm_params.d2 * delta_t.powi(2);
        let a_2 = - sgp4.atm_params.d3 * delta_t.powi(3) - sgp4.atm_params.d4 * delta_t.powi(4);
        a = (sgp4.wgs.ke / n).powf(2./3.) * (a_1 + a_2).powi(2);
        let il_1 = 3./2. * sgp4.atm_params.c1 * delta_t.powi(2);
        let il_2 = (sgp4.atm_params.d2 + 2. * sgp4.atm_params.c1.powi(2)) * delta_t.powi(3);
        let il_3 = 1./4. * (3. * sgp4.atm_params.d3 + 12. * sgp4.atm_params.c1 * sgp4.atm_params.d2 + 10. * sgp4.atm_params.c1.powi(3)) * delta_t.powi(4);
        let il_4 = 1./5. * (3. * sgp4.atm_params.d4 + 12. * sgp4.atm_params.c1 * sgp4.atm_params.d3 + 6. * sgp4.atm_params.d2.powi(2) + 30. * sgp4.atm_params.c1.powi(2) * sgp4.atm_params.d2 + 15. * sgp4.atm_params.c1.powi(4)) * delta_t.powi(5);
        il_atm = sgp4.brouwer0.n * (il_1 + il_2 + il_3 + il_4);
    }

    // Mean eccentricity check before the near-zero floor (Vallado)
    if e >= 1.0 || e < -0.001 {
        panic!("mean eccentricity {} is outside the range 0.0 to 1.0", e);
    }

    // Vallado eccentricity guard after atmospheric drag
    if e < 1.0e-6 {
        e = 1.0e-6;
    }

    // Account for long-period periodic effects of lunar and solar gravity
    if sgp4.deep_space {
        let m_m = sgp4.lunar_params.m + sgp4.lunar_params.n * delta_t;
        let m_s = sgp4.solar_params.m + sgp4.solar_params.n * delta_t;
        let f_m = m_m + 2. * sgp4.lunar_params.e * m_m.sin();
        let f_s = m_s + 2. * sgp4.solar_params.e * m_s.sin();
        let f2_m = 0.5 * f_m.sin().powi(2) - 0.25;
        let f2_s = 0.5 * f_s.sin().powi(2) - 0.25;
        let f3_m = -0.5 * f_m.sin() * f_m.cos();
        let f3_s = -0.5 * f_s.sin() * f_s.cos();
        let delta_e_m = -(30. * sgp4.brouwer0.beta * sgp4.lunar_params.c * sgp4.brouwer0.e / sgp4.brouwer0.n) * (f2_m * (sgp4.lunar_params.x2 * sgp4.lunar_params.x3 + sgp4.lunar_params.x1 * sgp4.lunar_params.x4) + f3_m * (sgp4.lunar_params.x2 * sgp4.lunar_params.x4 - sgp4.lunar_params.x1 * sgp4.lunar_params.x3));
        let delta_e_s = -(30. * sgp4.brouwer0.beta * sgp4.solar_params.c * sgp4.brouwer0.e / sgp4.brouwer0.n) * (f2_s * (sgp4.solar_params.x2 * sgp4.solar_params.x3 + sgp4.solar_params.x1 * sgp4.solar_params.x4) + f3_s * (sgp4.solar_params.x2 * sgp4.solar_params.x4 - sgp4.solar_params.x1 * sgp4.solar_params.x3));
        let delta_i_m = -(sgp4.lunar_params.c / sgp4.brouwer0.n / sgp4.brouwer0.beta) * (f2_m * sgp4.lunar_params.z12 + f3_m * (sgp4.lunar_params.z13 - sgp4.lunar_params.z11));
        let delta_i_s = -(sgp4.solar_params.c / sgp4.brouwer0.n / sgp4.brouwer0.beta) * (f2_s * sgp4.solar_params.z12 + f3_s * (sgp4.solar_params.z13 - sgp4.solar_params.z11));
        let delta_m_m = -(2. * sgp4.lunar_params.c / sgp4.brouwer0.n) * (f2_m * sgp4.lunar_params.z2 + f3_m * (sgp4.lunar_params.z3 - sgp4.lunar_params.z1) - 3. * sgp4.lunar_params.e * f_m.sin() * (7. + 3. * sgp4.brouwer0.e.powi(2)));
        let delta_m_s = -(2. * sgp4.solar_params.c / sgp4.brouwer0.n) * (f2_s * sgp4.solar_params.z2 + f3_s * (sgp4.solar_params.z3 - sgp4.solar_params.z1) - 3. * sgp4.solar_params.e * f_s.sin() * (7. + 3. * sgp4.brouwer0.e.powi(2)));
        let delta_raan_m = (sgp4.lunar_params.c / sgp4.brouwer0.n / sgp4.brouwer0.beta) * (f2_m * sgp4.lunar_params.z22 + f3_m * (sgp4.lunar_params.z23 - sgp4.lunar_params.z21));// / sgp4.lunar_params.i.sin();
        let delta_raan_s = (sgp4.solar_params.c / sgp4.brouwer0.n / sgp4.brouwer0.beta) * (f2_s * sgp4.solar_params.z22 + f3_s * (sgp4.solar_params.z23 - sgp4.solar_params.z21));// / sgp4.solar_params.i.sin();
        let delta_omega_m = (2. * sgp4.brouwer0.beta * sgp4.lunar_params.c / sgp4.brouwer0.n) * (f2_m * sgp4.lunar_params.z32 + f3_m * (sgp4.lunar_params.z33 - sgp4.lunar_params.z31) - 9. * sgp4.lunar_params.e * f_m.sin());// - delta_raan_m * sgp4.lunar_params.i.cos();
        let delta_omega_s = (2. * sgp4.brouwer0.beta * sgp4.solar_params.c / sgp4.brouwer0.n) * (f2_s * sgp4.solar_params.z32 + f3_s * (sgp4.solar_params.z33 - sgp4.solar_params.z31) - 9. * sgp4.solar_params.e * f_s.sin());// - delta_raan_s * sgp4.solar_params.i.cos();
        let delta_e_ls = delta_e_m + delta_e_s;
        let delta_i_ls = delta_i_m + delta_i_s;
        let delta_m_ls = delta_m_m + delta_m_s;
        let delta_raan_ls = delta_raan_m + delta_raan_s;
        let delta_omega_ls = delta_omega_m + delta_omega_s;

        e += delta_e_ls;
        if e < 0.0 || e > 1.0 {
            panic!("perturbed eccentricity {} is outside the range 0.0 to 1.0", e);
        }
        i += delta_i_ls;

        if i > 0.2 {
            // Notation is confusing in paper, delta_omega_ls stores more than just delta_omega
            // Same for delta_raan_ls, stores more than just delta_raan
            raan += delta_raan_ls / i.sin();
            omega += delta_omega_ls - delta_raan_ls * i.cos() / i.sin();
            m += delta_m_ls;
        } else {
            // Lyddane modification for inclinations below 0.2 rad (legacy)
            let alpha = i.sin() * raan.sin() + raan.cos() * delta_raan_ls + i.cos() * raan.sin() * delta_i_ls;
            let beta = i.sin() * raan.cos() - raan.sin() * delta_raan_ls + i.cos() * raan.cos() * delta_i_ls;
            let raan_old = raan;
            
            // Use Vallado's modification which is more numerically stable when i ~ 0
            let m_omega_raan = m + omega + delta_m_ls + delta_omega_ls + (i.cos() - delta_i_ls * i.sin()) * raan;
            
            // Calculate RAAN
            raan = alpha.atan2(beta);
            // Maintain RAAN continuity across atan2 branch cuts
            if (raan_old - raan).abs() > PI {
                if raan < raan_old {
                    raan += 2.0 * PI;
                } else {
                    raan -= 2.0 * PI;
                }
            }
            
            // Calculate omega
            m += delta_m_ls;
            omega = m_omega_raan - m - i.cos() * raan;
        }
    }

    // Vallado keep inclination in [0, π] after lunisolar periodics
    if sgp4.deep_space && i < 0.0 {
        i = -i;
        raan += PI;
        omega -= PI;
    }

    // Vallado mean-element recovery before long-period periodics
    m += il_atm;
    let mut lm = m + omega + raan;
    raan = raan.rem_euclid(2.0 * PI);
    omega = omega.rem_euclid(2.0 * PI);
    lm = lm.rem_euclid(2.0 * PI);
    m = (lm - omega - raan).rem_euclid(2.0 * PI);
    let il = m + omega + raan;

    // Account for long-period periodic effects of Earth's gravity
    let beta_update = (1. - e.powi(2)).sqrt();
    let a30 = -sgp4.wgs.j3; // [Earth Radii^3]
    let axn = e * omega.cos();
    let ill = a30 * i.sin() / (8. * sgp4.wgs.k2 * a * beta_update.powi(2)) * e * omega.cos() * (3. + 5. * i.cos()) / (1. + i.cos());
    let aynl = a30 * i.sin() / (4. * sgp4.wgs.k2 * a * beta_update.powi(2));
    let ilt = il + ill;
    let ayn = e * omega.sin() + aynl;

    // Account for short-period periodic effects of Earth's gravity (solve Kepler's equation)
    let u = (ilt - raan).rem_euclid(2.0 * PI);
    let mut e_omega = u.clone();
    let mut delta_e_omega: f64;

    // Newton-Raphson iteration to solve Kepler's equation (10 iterations max per Vallado)
    for _ in 0..10 {
        delta_e_omega = (u - ayn * e_omega.cos() + axn * e_omega.sin() - e_omega) / (1. - ayn * e_omega.sin() - axn * e_omega.cos());

        // Protect against oversized steps
        if delta_e_omega.abs() >= 0.95 {
            if delta_e_omega > 0.0 {
                delta_e_omega = 0.95;
            } else {
                delta_e_omega = -0.95;
            }
        }

        e_omega += delta_e_omega;
        
        // Verify convergence
        if delta_e_omega.abs() < 1e-12 {
            break;
        }
    }

    // Return position and velocity vectors in the TEME frame
    e = (axn.powi(2) + ayn.powi(2)).sqrt();
    let pl = a * (1. - e.powi(2));
    if pl < 0.0 {
        panic!("semilatus rectum {} is less than zero", pl);
    }
    let cos_ecc_anomaly = (axn * e_omega.cos() + ayn * e_omega.sin()) / e;
    let sin_ecc_anomaly = (axn * e_omega.sin() - ayn * e_omega.cos()) / e;
    let r = a * (1. - e * cos_ecc_anomaly);
    let r_dot = sgp4.wgs.ke * a.sqrt() * e * sin_ecc_anomaly / r;
    let r_f_dot = sgp4.wgs.ke * pl.sqrt() / r;
    let cos_u = a / r * (e_omega.cos() - axn + ayn * (e * sin_ecc_anomaly) / (1. + (1. - e.powi(2)).sqrt()));
    let sin_u = a / r * (e_omega.sin() - ayn - axn * (e * sin_ecc_anomaly) / (1. + (1. - e.powi(2)).sqrt()));
    let u = sin_u.atan2(cos_u);
    let delta_r = sgp4.wgs.k2 / (2. * pl) * (1. - i.cos().powi(2)) * (2. * u).cos();
    let delta_u = - sgp4.wgs.k2 / (4. * pl.powi(2)) * (7. * i.cos().powi(2) - 1.) * (2. * u).sin();
    let delta_raan = 3. * sgp4.wgs.k2 * i.cos() / (2. * pl.powi(2)) * (2. * u).sin();
    let delta_i = 3. * sgp4.wgs.k2 * i.cos() / (2. * pl.powi(2)) * i.sin() * (2. * u).cos();
    let delta_r_dot = - sgp4.wgs.k2 * n / pl * (1. - i.cos().powi(2)) * (2. * u).sin();
    let delta_r_f_dot = sgp4.wgs.k2 * n / pl * ((1. - i.cos().powi(2)) * (2. * u).cos() - 3./2. * (1. - 3. * i.cos().powi(2)));
    let rk = r * (1. - 3./2. * sgp4.wgs.k2 * (1. - e.powi(2)).sqrt() / pl.powi(2) * (3. * i.cos().powi(2) - 1.)) + delta_r;
    if rk < 1.0 {
        panic!("satellite has decayed: position radius {} Earth radii is less than 1.0", rk);
    }
    let uk = u + delta_u;
    let raan_k = raan + delta_raan;
    let i_k = i + delta_i;
    let r_dot_k = r_dot + delta_r_dot;
    let r_f_dot_k = r_f_dot + delta_r_f_dot;

    let mx = -raan_k.sin() * i_k.cos();
    let my = raan_k.cos() * i_k.cos();
    let mz = i_k.sin();

    let nx = raan_k.cos();
    let ny = raan_k.sin();
    let nz = 0.;

    let ux = mx * uk.sin() + nx * uk.cos();
    let uy = my * uk.sin() + ny * uk.cos();
    let uz = mz * uk.sin() + nz * uk.cos();

    let vx = mx * uk.cos() - nx * uk.sin();
    let vy = my * uk.cos() - ny * uk.sin();
    let vz = mz * uk.cos() - nz * uk.sin();

    let rx = rk * ux * sgp4.wgs.r_earth_eq;
    let ry = rk * uy * sgp4.wgs.r_earth_eq;
    let rz = rk * uz * sgp4.wgs.r_earth_eq;

    let r_dot_x = (r_dot_k * ux + r_f_dot_k * vx) * sgp4.wgs.r_earth_eq / 60.;
    let r_dot_y = (r_dot_k * uy + r_f_dot_k * vy) * sgp4.wgs.r_earth_eq / 60.;
    let r_dot_z = (r_dot_k * uz + r_f_dot_k * vz) * sgp4.wgs.r_earth_eq / 60.;

    let state_vector = StateVector {
        r_x: rx,
        r_y: ry,
        r_z: rz,
        v_x: r_dot_x,
        v_y: r_dot_y,
        v_z: r_dot_z,
        coordinate_frame: CoordinateFrame::TEME,
    };

    return state_vector;
}

/// Half day Euler-Maclaurin integration step
///
/// # Arguments
/// * `lami` - The auxilary variable at time i
/// * `ni` - The mean motion at time i
/// * `omegai` - The argument of perigee at time i
/// * `half_day_resonance_params` - The half day resonance parameters
///
/// # Returns
/// * `lami_dot` - The rate of change of the auxilary variable at time i+1
/// * `ni_dot` - The rate of change of the mean motion at time i+1
/// * `lami_ddot` - The 2nd derivative of the auxilary variable at time i+1
/// * `ni_ddot` - The 2nd derivative of the mean motion at time i+1
///
/// # Examples
/// ```rust
/// // Define auxilary variable at time i
/// let lami = 0.0;
///
/// // Define mean motion at time i
/// let ni = 0.0;
///
/// // Define half day resonance parameters
/// let half_day_resonance_params = HalfDayResonanceParams::default();
///
/// // Calculate the auxilary variable and mean motion at time i+1
/// let (lami_dot, ni_dot, lami_ddot, ni_ddot) = half_day_euler_maclaurin_step(lami, ni, omegai, &half_day_resonance_params);
/// ```
///
/// References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn half_day_euler_maclaurin_step(lami: f64, ni: f64, omegai: f64, half_day_resonance_params: &HalfDayResonanceParams) -> (f64, f64, f64, f64) {
    // Define constants
    let g22 = 5.7686396;
    let g32 = 0.95240898;
    let g44 = 1.8014998;
    let g52 = 1.0508330;
    let g54 = 4.4108898;

    // Calculate the rate of change of the auxilary variable
    let lami_dot = ni + half_day_resonance_params.lam0_dot;

    // Calculate the rate of change of the mean motion
    let ni_dot_2201 = half_day_resonance_params.d2201 * ((2. - 2. * 0.) * omegai + 2. / 2. * lami - g22).sin();
    let ni_dot_2211 = half_day_resonance_params.d2211 * ((2. - 2. * 1.) * omegai + 2. / 2. * lami - g22).sin();
    let ni_dot_3210 = half_day_resonance_params.d3210 * ((3. - 2. * 1.) * omegai + 2. / 2. * lami - g32).sin();
    let ni_dot_3222 = half_day_resonance_params.d3222 * ((3. - 2. * 2.) * omegai + 2. / 2. * lami - g32).sin();
    let ni_dot_5220 = half_day_resonance_params.d5220 * ((5. - 2. * 2.) * omegai + 2. / 2. * lami - g52).sin();
    let ni_dot_5232 = half_day_resonance_params.d5232 * ((5. - 2. * 3.) * omegai + 2. / 2. * lami - g52).sin();
    let ni_dot_4422 = half_day_resonance_params.d4422 * ((4. - 2. * 2.) * omegai + 4. / 2. * lami - g44).sin();
    let ni_dot_5421 = half_day_resonance_params.d5421 * ((5. - 2. * 2.) * omegai + 4. / 2. * lami - g54).sin();
    let ni_dot_5433 = half_day_resonance_params.d5433 * ((5. - 2. * 3.) * omegai + 4. / 2. * lami - g54).sin();
    let ni_dot_4410 = half_day_resonance_params.d4410 * ((4. - 2. * 1.) * omegai + 4. / 2. * lami - g44).sin();
    let ni_dot = ni_dot_2201 + ni_dot_2211 + ni_dot_3210 + ni_dot_3222 + ni_dot_5220 + ni_dot_5232 + ni_dot_4422 + ni_dot_5421 + ni_dot_5433 + ni_dot_4410;

    // Calculate the 2nd derivative of the auxilary variable
    let lami_ddot = ni_dot;

    // Calculate the 2nd derivative of the mean motion
    let ni_ddot_2201 = 2. / 2. * half_day_resonance_params.d2201 * ((2. - 2. * 0.) * omegai + 2. / 2. * lami - g22).cos();
    let ni_ddot_2211 = 2. / 2. * half_day_resonance_params.d2211 * ((2. - 2. * 1.) * omegai + 2. / 2. * lami - g22).cos();
    let ni_ddot_3210 = 2. / 2. * half_day_resonance_params.d3210 * ((3. - 2. * 1.) * omegai + 2. / 2. * lami - g32).cos();
    let ni_ddot_3222 = 2. / 2. * half_day_resonance_params.d3222 * ((3. - 2. * 2.) * omegai + 2. / 2. * lami - g32).cos();
    let ni_ddot_5220 = 2. / 2. * half_day_resonance_params.d5220 * ((5. - 2. * 2.) * omegai + 2. / 2. * lami - g52).cos();
    let ni_ddot_5232 = 2. / 2. * half_day_resonance_params.d5232 * ((5. - 2. * 3.) * omegai + 2. / 2. * lami - g52).cos();
    let ni_ddot_4422 = 4. / 2. * half_day_resonance_params.d4422 * ((4. - 2. * 2.) * omegai + 4. / 2. * lami - g44).cos();
    let ni_ddot_5421 = 4. / 2. * half_day_resonance_params.d5421 * ((5. - 2. * 2.) * omegai + 4. / 2. * lami - g54).cos();
    let ni_ddot_5433 = 4. / 2. * half_day_resonance_params.d5433 * ((5. - 2. * 3.) * omegai + 4. / 2. * lami - g54).cos();
    let ni_ddot_4410 = 4. / 2. * half_day_resonance_params.d4410 * ((4. - 2. * 1.) * omegai + 4. / 2. * lami - g44).cos();
    let ni_ddot = lami_dot * (ni_ddot_2201 + ni_ddot_2211 + ni_ddot_3210 + ni_ddot_3222 + ni_ddot_5220 + ni_ddot_5232 + ni_ddot_4422 + ni_ddot_5421 + ni_ddot_5433 + ni_ddot_4410);

    return (lami_dot, ni_dot, lami_ddot, ni_ddot);
}

/// Whole day Euler-Maclaurin integration step
///
/// # Arguments
/// * `lami` - The auxilary variable at time i
/// * `ni` - The mean motion at time i
/// * `whole_day_resonance_params` - The whole day resonance parameters
///
/// # Returns
/// * `lami_dot` - The rate of change of the auxilary variable at time i+1
/// * `ni_dot` - The rate of change of the mean motion at time i+1
/// * `lami_ddot` - The 2nd derivative of the auxilary variable at time i+1
/// * `ni_ddot` - The 2nd derivative of the mean motion at time i+1
///
/// # Examples
/// ```rust
/// // Define auxilary variable at time i
/// let lami = 0.0;
///
/// // Define mean motion at time i
/// let ni = 0.0;
///
/// // Define whole day resonance parameters
/// let whole_day_resonance_params = WholeDayResonanceParams::default();
///
/// // Calculate the auxilary variable and mean motion at time i+1
/// let (lami_dot, ni_dot, lami_ddot, ni_ddot) = whole_day_euler_maclaurin_step(lami, ni, &whole_day_resonance_params);
/// ```
///
/// References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn whole_day_euler_maclaurin_step(lami: f64, ni: f64, whole_day_resonance_params: &WholeDayResonanceParams) -> (f64, f64, f64, f64) {
    // Calculate the rate of change of the auxilary variable
    let lami_dot = ni + whole_day_resonance_params.lam0_dot;

    // Calculate the rate of change of the mean motion 
    let ni_dot_1 = whole_day_resonance_params.delta1 * (lami - whole_day_resonance_params.lam31).sin();
    let ni_dot_2 = whole_day_resonance_params.delta2 * (2. * (lami - whole_day_resonance_params.lam22)).sin();
    let ni_dot_3 = whole_day_resonance_params.delta3 * (3. * (lami - whole_day_resonance_params.lam33)).sin();
    let ni_dot = ni_dot_1 + ni_dot_2 + ni_dot_3;

    // Calculate the 2nd derivative of the auxilary variable
    let lami_ddot = ni_dot;

    // Calculate the 2nd derivative of the mean motion
    let ni_ddot_1 = whole_day_resonance_params.delta1 * (lami - whole_day_resonance_params.lam31).cos();
    let ni_ddot_2 = 2. * whole_day_resonance_params.delta2 * (2. * (lami - whole_day_resonance_params.lam22)).cos();
    let ni_ddot_3 = 3. * whole_day_resonance_params.delta3 * (3. * (lami - whole_day_resonance_params.lam33)).cos();
    let ni_ddot = lami_dot * (ni_ddot_1 + ni_ddot_2 + ni_ddot_3);

    return (lami_dot, ni_dot, lami_ddot, ni_ddot);
}

// ----------
// Unit Tests
// ----------
#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use serde::Deserialize;

    // -------------------------------------------------------
    // Structs for deserializing the Vallado test case TOML
    // -------------------------------------------------------

    #[derive(Deserialize)]
    struct ValladoCases {
        test: HashMap<String, ValladoCase>,
    }

    #[derive(Deserialize)]
    struct ValladoCase {
        name: String,
        tle: String,
        #[serde(default)]
        #[allow(dead_code)]
        start_mins_from_epoch: f64,
        #[serde(default)]
        #[allow(dead_code)]
        end_mins_from_epoch: f64,
        #[serde(default)]
        #[allow(dead_code)]
        delta_time_mins: f64,
        ephem: String,
        #[serde(default)]
        exception: bool,
    }

    /// One ephemeris row: time since epoch [min], TEME position [km], velocity [km/s].
    struct EphemRow {
        t_mins: f64,
        rx: f64,
        ry: f64,
        rz: f64,
        vx: f64,
        vy: f64,
        vz: f64,
    }

    /// Parse the whitespace-delimited Vallado ephem block into rows.
    /// Lines are `t_mins rx ry rz vx vy vz` with an optional trailing calendar stamp.
    fn parse_vallado_ephem(ephem: &str) -> Vec<EphemRow> {
        ephem
            .lines()
            .filter_map(|line| {
                let cols: Vec<&str> = line.split_whitespace().collect();
                if cols.len() < 7 {
                    return None;
                }
                Some(EphemRow {
                    t_mins: cols[0].parse().ok()?,
                    rx: cols[1].parse().ok()?,
                    ry: cols[2].parse().ok()?,
                    rz: cols[3].parse().ok()?,
                    vx: cols[4].parse().ok()?,
                    vy: cols[5].parse().ok()?,
                    vz: cols[6].parse().ok()?,
                })
            })
            .collect()
    }

    /// Sub-meter agreement with Vallado reference ephemerides [km] / [km/s].
    const VALLADO_STATE_TOL_KM: f64 = 1e-3;

    fn assert_state_near(key: &str, name: &str, t_mins: f64, got: &StateVector, row: &EphemRow) {
        let checks = [
            ("Position x", got.r_x, row.rx),
            ("Position y", got.r_y, row.ry),
            ("Position z", got.r_z, row.rz),
            ("Velocity x", got.v_x, row.vx),
            ("Velocity y", got.v_y, row.vy),
            ("Velocity z", got.v_z, row.vz),
        ];
        for (label, value, expected) in checks {
            assert!(
                (value - expected).abs() < VALLADO_STATE_TOL_KM,
                "{key}: {name}\n- {label} mismatch {value} vs {expected} (t_mins = {t_mins})"
            );
        }
    }

    #[test]
    fn test_checksum_calculation() {
        // Define the TLE line
        let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";

        // Calculate the checksum of the TLE line
        let checksum = calc_checksum(tle_line1);

        // Assert the checksum is correct
        assert_eq!(checksum, 1);
    }

    #[test]
    fn test_checksum_comparison() {
        // Define the TLE line
        let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
        let tle_line2 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2922";

        // Calculate the checksum of the TLE line
        let checksum = tle_checksum(tle_line1);
        let checksum2 = tle_checksum(tle_line2);

        // Assert the checksum is correct
        assert_eq!(checksum, true);
        assert_eq!(checksum2, false);
    }

    #[test]
    fn test_tle_parsing_from_lines() {
        // Define the TLE lines
        let tle_line0 = "ISS (ZARYA)";
        let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
        let tle_line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";

        // Parse the TLE lines into a TLE struct
        let sgp4 = from_tle_lines(tle_line1, tle_line2, Some(tle_line0));

        // Assert the TLE struct is correct
        assert_eq!(sgp4.gp.common_name, "ISS (ZARYA)");
        assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
        assert_eq!(sgp4.gp.classification, 'U');
        assert_eq!(sgp4.gp.international_designator, "98067A");
        assert_eq!(sgp4.gp.epoch_datetime.year, 2008);
        assert_eq!(sgp4.gp.epoch_datetime.month, 9);
        assert_eq!(sgp4.gp.epoch_datetime.day, 20);
        assert_eq!(sgp4.gp.epoch_datetime.hour, 12);
        assert_eq!(sgp4.gp.epoch_datetime.minute, 25);
        assert!((sgp4.gp.epoch_datetime.second - 40.104192).abs() < 1e-9);
        assert_eq!(sgp4.gp.first_derivative_of_mean_motion, -0.00004364);
        assert!((sgp4.gp.second_derivative_of_mean_motion + 6.0e-5).abs() < 1e-12);
        assert_eq!(sgp4.gp.bstar, -0.000011606);
        assert_eq!(sgp4.gp.ephemeris_type, 0);
        assert_eq!(sgp4.gp.element_set_number, 292);
        assert_eq!(sgp4.gp.inclination, 51.6416);
        assert_eq!(sgp4.gp.right_ascension_of_ascending_node, 247.4627);
        assert_eq!(sgp4.gp.eccentricity, 0.0006703);
        assert_eq!(sgp4.gp.argument_of_perigee, 130.536);
        assert_eq!(sgp4.gp.mean_anomaly, 325.0288);
        assert_eq!(sgp4.gp.mean_motion, 15.72125391);
        assert_eq!(sgp4.gp.revolution_number_at_epoch, 56353);
    }

    #[test]
    fn test_tle_parsing_from_string() {
        // Define the TLE string
        let tle_string = "ISS (ZARYA)\n1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921\n2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";

        // Parse the TLE string into a SGP4 struct
        let sgp4s = from_tle_string(tle_string);
        let sgp4 = &sgp4s[0];

        // Assert the TLE struct is correct
        assert_eq!(sgp4.gp.common_name, "ISS (ZARYA)");
        assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
        assert_eq!(sgp4.gp.classification, 'U');
        assert_eq!(sgp4.gp.international_designator, "98067A");
        assert_eq!(sgp4.gp.epoch_datetime.year, 2008);
        assert_eq!(sgp4.gp.epoch_datetime.month, 9);
        assert_eq!(sgp4.gp.epoch_datetime.day, 20);
        assert_eq!(sgp4.gp.epoch_datetime.hour, 12);
        assert_eq!(sgp4.gp.epoch_datetime.minute, 25);
        assert!((sgp4.gp.epoch_datetime.second - 40.104192).abs() < 1e-9);
        assert_eq!(sgp4.gp.first_derivative_of_mean_motion, -0.00004364);
        assert!((sgp4.gp.second_derivative_of_mean_motion + 6.0e-5).abs() < 1e-12);
        assert_eq!(sgp4.gp.bstar, -0.000011606);
        assert_eq!(sgp4.gp.ephemeris_type, 0);
        assert_eq!(sgp4.gp.element_set_number, 292);
        assert_eq!(sgp4.gp.inclination, 51.6416);
        assert_eq!(sgp4.gp.right_ascension_of_ascending_node, 247.4627);
        assert_eq!(sgp4.gp.eccentricity, 0.0006703);
        assert_eq!(sgp4.gp.argument_of_perigee, 130.536);
        assert_eq!(sgp4.gp.mean_anomaly, 325.0288);
        assert_eq!(sgp4.gp.mean_motion, 15.72125391);
        assert_eq!(sgp4.gp.revolution_number_at_epoch, 56353);
    }

    #[test]
    fn test_tle_parsing_from_file() {
        // Define the TLE file path
        let tle_file_path = "test/tle_file.txt";

        // Parse the TLE file into a SGP4 struct
        let sgp4s = from_tle_file(tle_file_path);
        let sgp4_count = sgp4s.len();
        let iss_sgp4 = &sgp4s[12];
        let hulianwang_sgp4 = &sgp4s[17];

        // Assert the TLE structs are correct
        assert_eq!(sgp4_count, 19);

        assert_eq!(iss_sgp4.gp.common_name, "ISS (ZARYA)");
        assert_eq!(iss_sgp4.gp.satellite_catalog_number, 25544);
        assert_eq!(iss_sgp4.gp.classification, 'U');
        assert_eq!(iss_sgp4.gp.international_designator, "98067A");
        assert_eq!(iss_sgp4.gp.epoch_datetime.year, 2008);
        assert_eq!(iss_sgp4.gp.epoch_datetime.month, 9);
        assert_eq!(iss_sgp4.gp.epoch_datetime.day, 20);
        assert_eq!(iss_sgp4.gp.epoch_datetime.hour, 12);
        assert_eq!(iss_sgp4.gp.epoch_datetime.minute, 25);
        assert!((iss_sgp4.gp.epoch_datetime.second - 40.104192).abs() < 1e-9);
        assert_eq!(iss_sgp4.gp.first_derivative_of_mean_motion, -0.00004364);
        assert!((iss_sgp4.gp.second_derivative_of_mean_motion + 6.0e-5).abs() < 1e-12);
        assert_eq!(iss_sgp4.gp.bstar, -0.000011606);
        assert_eq!(iss_sgp4.gp.ephemeris_type, 0);
        assert_eq!(iss_sgp4.gp.element_set_number, 292);
        assert_eq!(iss_sgp4.gp.inclination, 51.6416);
        assert_eq!(iss_sgp4.gp.right_ascension_of_ascending_node, 247.4627);
        assert_eq!(iss_sgp4.gp.eccentricity, 0.0006703);
        assert_eq!(iss_sgp4.gp.argument_of_perigee, 130.536);
        assert_eq!(iss_sgp4.gp.mean_anomaly, 325.0288);
        assert_eq!(iss_sgp4.gp.mean_motion, 15.72125391);
        assert_eq!(iss_sgp4.gp.revolution_number_at_epoch, 56353);

        assert_eq!(hulianwang_sgp4.gp.common_name, "HULIANWANG DIGUI-118");
        assert_eq!(hulianwang_sgp4.gp.satellite_catalog_number, 66957);
        assert_eq!(hulianwang_sgp4.gp.classification, 'U');
        assert_eq!(hulianwang_sgp4.gp.international_designator, "25287E");
        assert_eq!(hulianwang_sgp4.gp.epoch_datetime.year, 2025);
        assert_eq!(hulianwang_sgp4.gp.epoch_datetime.month, 12);
        assert_eq!(hulianwang_sgp4.gp.epoch_datetime.day, 12);
        assert_eq!(hulianwang_sgp4.gp.epoch_datetime.hour, 16);
        assert_eq!(hulianwang_sgp4.gp.epoch_datetime.minute, 47);
        assert!((hulianwang_sgp4.gp.epoch_datetime.second - 31.7748479989).abs() < 1e-9);
        assert_eq!(hulianwang_sgp4.gp.first_derivative_of_mean_motion, -0.00000302);
        assert_eq!(hulianwang_sgp4.gp.second_derivative_of_mean_motion, 0.0);
        assert!((hulianwang_sgp4.gp.bstar + 1.9373e-4).abs() < 1e-12);
        assert_eq!(hulianwang_sgp4.gp.ephemeris_type, 0);
        assert_eq!(hulianwang_sgp4.gp.element_set_number, 999);
        assert_eq!(hulianwang_sgp4.gp.inclination, 86.4945);
        assert_eq!(hulianwang_sgp4.gp.right_ascension_of_ascending_node, 346.1700);
        assert_eq!(hulianwang_sgp4.gp.eccentricity, 0.0007219);
        assert_eq!(hulianwang_sgp4.gp.argument_of_perigee, 190.5502);
        assert_eq!(hulianwang_sgp4.gp.mean_anomaly, 169.5507);
        assert_eq!(hulianwang_sgp4.gp.mean_motion, 13.69137019);
        assert_eq!(hulianwang_sgp4.gp.revolution_number_at_epoch, 52);
    }
    
    #[test]
    fn test_sgp4_vallado_cases() {
        let content = std::fs::read_to_string("test/vallado_cases.toml")
            .expect("could not read test/vallado_cases.toml");
        let cases: ValladoCases = toml::from_str(&content)
            .expect("could not parse test/vallado_cases.toml");

        let mut keys: Vec<&String> = cases.test.keys().collect();
        keys.sort();

        for key in keys {
            let case = &cases.test[key];

            if case.exception {
                let result = std::panic::catch_unwind(|| from_tle_string(&case.tle));
                assert!(
                    result.is_err(),
                    "case {key}: expected initialization to panic"
                );
                continue;
            }

            let sgp4s = from_tle_string(&case.tle);
            assert!(!sgp4s.is_empty(), "case {key}: TLE failed to parse");
            let sgp4 = &sgp4s[0];

            for row in parse_vallado_ephem(&case.ephem) {
                let state = sgp4_prop_delta(sgp4, row.t_mins);
                assert_state_near(key, &case.name, row.t_mins, &state, &row);
            }
        }
    }
}