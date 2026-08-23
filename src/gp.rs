//! Module to handle the input and processing of GP (General Perturbation)
//! elements. Element sets can be parsed from TLE (including Alpha-5 catalog
//! numbers) and OMM formats.

// ------------------
// External Libraries
// ------------------
use std::collections::HashMap;
use std::fs;

use csv::ReaderBuilder;
use roxmltree::{Document, Node};
use serde_json::Value;

// ------------------
// Internal Libraries
// ------------------
use crate::time::{dayofyr2utc, DateTime, Timezone};
use crate::sgp4::{init_sgp4, Sgp4};

// -------
// Structs
// -------

/// General Perturbation Element Set for an Earth-orbiting satellite.
///
/// This struct represents the parsed contents of a standard General Perturbation Element Set.
/// The General Perturbation Element Set is the standard set of orbital elements used with the
/// SGP4 propagator. Elements are commonly distributed as Two-Line Elements (TLEs) or Orbit
/// Mean-Elements Messages (OMMs).
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
#[derive(Default, Clone)]
pub struct GenPerturbElementSet {
    /// Common name of the satellite (e.g., "ISS (ZARYA)")
    pub common_name: String,

    /// NORAD satellite catalog number
    pub satellite_catalog_number: i32,

    /// Classification (`U` = Unclassified, `C` = Classified, `S` = Secret)
    pub classification: char,

    /// International designator in `YYYY-NNNP` form (launch year, launch number, piece)
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

// -----
// Enums
// -----

// ------
// Traits
// ------

/// Convert an OMM field value into a typed GP field
///
/// Used by KVN and XML OMM parsing to turn a present field string into String,
/// char, integer, floating-point, or DateTime values, and to supply defaults
/// when a field is missing or empty.
trait FromOmm: Sized {
    /// Default used when the OMM field is missing or empty
    fn omm_default() -> Self;

    /// Parse a present, non-empty OMM field value
    fn from_omm(field: &str, value: &str) -> Self;
}

impl FromOmm for String {
    fn omm_default() -> Self {
        String::new()
    }

    fn from_omm(_field: &str, value: &str) -> Self {
        value.to_string()
    }
}

impl FromOmm for char {
    fn omm_default() -> Self {
        'U'
    }

    fn from_omm(_field: &str, value: &str) -> Self {
        value.chars().next().unwrap_or('U')
    }
}

impl FromOmm for f64 {
    fn omm_default() -> Self {
        0.0
    }

    fn from_omm(field: &str, value: &str) -> Self {
        let normalized = value.replace('D', "E").replace('d', "E");
        normalized.parse::<f64>().unwrap_or_else(|_| {
            panic!("OMM field {} is not a valid number: {}", field, value);
        })
    }
}

impl FromOmm for i32 {
    fn omm_default() -> Self {
        0
    }

    fn from_omm(field: &str, value: &str) -> Self {
        value.parse::<i32>().unwrap_or_else(|_| {
            panic!("OMM field {} is not a valid integer: {}", field, value);
        })
    }
}

impl FromOmm for i64 {
    fn omm_default() -> Self {
        0
    }

    fn from_omm(field: &str, value: &str) -> Self {
        value.parse::<i64>().unwrap_or_else(|_| {
            panic!("OMM field {} is not a valid integer: {}", field, value);
        })
    }
}

impl FromOmm for DateTime {
    fn omm_default() -> Self {
        DateTime::default()
    }

    fn from_omm(_field: &str, value: &str) -> Self {
        parse_omm_epoch(value)
    }
}

// ---------
// Constants
// ---------

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
/// use mako_sgp4::gp::from_tle_lines;
///
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
            panic!(
                "TLE line 0 is invalid: name must be 1-24 characters, got {}",
                name_line.len()
            );
        }
        gp.common_name = name_line.to_string();
    }

    // Parse through line 1 and populate TLE struct
    if line1.len() < 69 || line1.len() > 69 {
        panic!(
            "TLE line 1 is invalid: must be 69 characters, got {}",
            line1.len()
        );
    } else {
        // Line 1
        // Satellite catalog number
        let catalog = line1[2..7].trim();
        let first_char = catalog.chars().next().unwrap();
        if first_char.is_ascii_digit() {
            // This is a classic numeric catalog number
            gp.satellite_catalog_number = catalog.parse::<i32>().unwrap();
        } else if first_char.is_ascii_alphabetic() {
            // This is an Alpha-5 catalog number
            gp.satellite_catalog_number =
                alpha5_digit(first_char).unwrap() * 10000 + catalog[1..].parse::<i32>().unwrap();
        } else {
            // This is an unexpected catalog number
            panic!(
                "TLE line 1 is invalid: catalog number is invalid: {}",
                catalog
            );
        }

        // Classification
        gp.classification = line1[7..8].trim().parse::<char>().unwrap();

        // International designator (COSPAR YYYY-NNNP)
        let intl_des = line1[9..17].trim();
        if intl_des.is_empty() {
            // Handle case where international designator is not present
            gp.international_designator = "".to_string();
        } else {
            // Expand two-digit launch year and insert the COSPAR dash
            let launch_year = tle_full_year(intl_des[0..2].parse::<i32>().unwrap());
            gp.international_designator = format!("{}-{}", launch_year, &intl_des[2..]);
        }

        // Epoch year (last two numbers)
        let yr_two_digit = line1[18..20].trim().parse::<i32>().unwrap();
        let epoch_year = tle_full_year(yr_two_digit);

        // Epoch day of year
        let epoch_day = line1[20..32].trim().parse::<f64>().unwrap();

        // Epoch UTC datetime
        let Some(epoch_datetime) = dayofyr2utc(epoch_year, epoch_day).ok() else {
            panic!(
                "Error converting epoch day of year to UTC datetime: Epoch year: {}, Epoch day: {}",
                epoch_year, epoch_day
            );
        };
        gp.epoch_datetime = epoch_datetime;

        // 1st derivative of mean motion [revs/day^2]
        gp.first_derivative_of_mean_motion = line1[33..43].trim().parse::<f64>().unwrap() * 2.0;

        // 2nd derivative of mean motion [revs/days^3]
        // Account for - in 2nd derivative of mean motion
        if line1[44..45].parse::<char>().unwrap() == '-' {
            gp.second_derivative_of_mean_motion = format!("-0.{}", line1[45..50].trim())
                .parse::<f64>()
                .unwrap()
                * 10.0_f64.powi(line1[50..52].parse::<i32>().unwrap())
                * 6.0_f64;
        } else {
            gp.second_derivative_of_mean_motion = format!("0.{}", line1[45..50].trim())
                .parse::<f64>()
                .unwrap()
                * 10.0_f64.powi(line1[50..52].parse::<i32>().unwrap())
                * 6.0_f64;
        }

        // B* [1/Earth Radii]
        // Account for - in B* term
        if line1[53..54].parse::<char>().unwrap() == '-' {
            gp.bstar = format!("-0.{}", line1[54..59].trim())
                .parse::<f64>()
                .unwrap()
                * 10.0_f64.powi(line1[59..61].parse::<i32>().unwrap());
        } else {
            gp.bstar = format!("0.{}", line1[54..59].trim())
                .parse::<f64>()
                .unwrap()
                * 10.0_f64.powi(line1[59..61].parse::<i32>().unwrap());
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
        panic!(
            "TLE line 2 is invalid: must be 69 characters, got {}",
            line2.len()
        );
    } else {
        // Line 2
        // Inclination [degs]
        gp.inclination = line2[8..16].trim().parse::<f64>().unwrap();

        // Right ascension of ascending node [degs]
        gp.right_ascension_of_ascending_node = line2[17..25].trim().parse::<f64>().unwrap();

        // Eccentricity
        gp.eccentricity = format!("0.{}", line2[26..33].trim())
            .parse::<f64>()
            .unwrap();

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
/// use mako_sgp4::gp::from_tle_string;
///
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
/// use mako_sgp4::gp::from_tle_file;
///
/// // Define the TLE file path
/// let tle_file_path = "test/tle_parsing_cases.txt";
///
/// // Parse the TLE file into a vector of SGP4 structs
/// let sgp4s = from_tle_file(tle_file_path);
/// let sgp4 = &sgp4s[0];
///
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
/// ```
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
pub fn from_tle_file(file_path: &str) -> Vec<Sgp4> {
    // Open the TLE file
    let tle_string = fs::read_to_string(file_path).expect("Cannot read TLE file");

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
/// use mako_sgp4::gp::calc_checksum;
///
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
/// use mako_sgp4::gp::tle_checksum;
///
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

/// Convert a two-digit TLE year to a four-digit Gregorian year.
///
/// TLE years follow the convention that 57-99 map to 1957-1999 and 00-56 map
/// to 2000-2056, matching the original SGP4 epoch encoding.
///
/// # Arguments
/// * `two_digit_year` - Year modulo 100 as stored in a TLE
///
/// # Returns
/// * Four-digit year in the range 1957-2056
fn tle_full_year(two_digit_year: i32) -> i32 {
    if two_digit_year < 57 {
        return 2000 + two_digit_year;
    } else {
        return 1900 + two_digit_year;
    }
}

/// Converts an Alpha-5 digit to an integer
///
/// Given a character, map it to the corresponding integer value as defined in the Alpha-5 standard.
///
/// # Arguments
/// * `c` - The character to convert
///
/// # Returns
/// * `Some(i)` - The integer value of the character
/// * `None` - If the character is not a valid Alpha-5 digit
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::alpha5_digit;
///
/// // Define the character
/// let c = 'A';
///
/// // Convert the character to an integer
/// let i = alpha5_digit(c);
///
/// // Assert the integer is correct
/// assert_eq!(i, Some(10));
/// ```
///
/// # References
/// - [Alpha-5 Standard](https://www.space-track.org/documentation#tle-alpha5)
pub fn alpha5_digit(c: char) -> Option<i32> {
    // Convert the character to uppercase to match the standard
    let c_upper = c.to_ascii_uppercase();

    // Match the character to the corresponding integer value
    match c_upper {
        'A'..='H' => Some((c_upper as i32) - ('A' as i32) + 10),
        'J'..='N' => Some((c_upper as i32) - ('J' as i32) + 18),
        'P'..='Z' => Some((c_upper as i32) - ('P' as i32) + 23),
        _ => None, // Standard skips 'I' and 'O'
    }
}

/// Builds a [`Sgp4`] struct from the lines of a single OMM KVN record.
///
/// Given the key-value lines of one Orbit Mean-Elements Message, this function
/// parses the input into a [`GenPerturbElementSet`] and initializes SGP4,
/// matching [`from_tle_lines`]. Missing numeric fields default to 0 and missing
/// string fields default to an empty string.
///
/// # Arguments
/// * `lines` - The KVN key-value lines of one OMM record
///
/// # Returns
/// * [`Sgp4`] - Struct containing the parsed SGP4 parameters
///
/// # Panics
/// * If a present field cannot be parsed as the requested type
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_kvn_lines;
///
/// // Define the OMM KVN lines
/// let omm = "\
/// OBJECT_NAME    = 2026-106A
/// OBJECT_ID      = 2026-106A
/// EPOCH          = 2026-06-14T15:07:48.259488
/// MEAN_MOTION    = 15.11169557
/// ECCENTRICITY   = .00147468
/// INCLINATION    = 97.5103
/// RA_OF_ASC_NODE = 247.7605
/// ARG_OF_PERICENTER = 169.6213
/// MEAN_ANOMALY   = 190.5325
/// NORAD_CAT_ID   = 69097
/// BSTAR          = .39221734E-3
/// MEAN_MOTION_DOT = .6535E-4
/// MEAN_MOTION_DDOT = 0
/// ";
/// let lines: Vec<&str> = omm.lines().map(str::trim).filter(|l| !l.is_empty()).collect();
///
/// // Parse the OMM lines into a SGP4 struct
/// let sgp4 = from_omm_kvn_lines(&lines);
///
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_kvn_lines(lines: &[&str]) -> Sgp4 {
    return sgp4_from_omm_lookup(|field| kvn_lookup(lines, field));
}

/// Builds a vector of [`Sgp4`] structs from a string containing Orbit Mean-Elements Message
/// (OMM) sets in CCSDS KVN (key-value notation) format.
///
/// This function parses a string containing one or more OMM records. A new record
/// starts at each CCSDS_OMM_VERS line. Each record is parsed by
/// [`from_omm_kvn_lines`].
///
/// # Arguments
/// * `omm_kvn_string` - A string containing one or more OMM records in KVN format
///
/// # Returns
/// * `Vec<Sgp4>` - A vector containing all successfully parsed SGP4 parameters
///
/// # Panics
/// * If a present OMM field cannot be parsed
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_kvn_string;
///
/// // Define the OMM KVN string
/// let omm = "\
/// CCSDS_OMM_VERS = 2.0
/// OBJECT_NAME    = 2026-106A
/// OBJECT_ID      = 2026-106A
/// EPOCH          = 2026-06-14T15:07:48.259488
/// MEAN_MOTION    = 15.11169557
/// ECCENTRICITY   = .00147468
/// INCLINATION    = 97.5103
/// RA_OF_ASC_NODE = 247.7605
/// ARG_OF_PERICENTER = 169.6213
/// MEAN_ANOMALY   = 190.5325
/// EPHEMERIS_TYPE = 0
/// CLASSIFICATION_TYPE = U
/// NORAD_CAT_ID   = 69097
/// ELEMENT_SET_NO = 999
/// REV_AT_EPOCH   = 459
/// BSTAR          = .39221734E-3
/// MEAN_MOTION_DOT = .6535E-4
/// MEAN_MOTION_DDOT = 0
/// ";
///
/// // Parse the OMM string into a SGP4 struct
/// let sgp4s = from_omm_kvn_string(omm);
/// let sgp4 = &sgp4s[0];
///
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_kvn_string(omm_kvn_string: &str) -> Vec<Sgp4> {
    // Split the string into individual OMM records
    let mut sgp4s = Vec::new();
    let mut current_lines: Vec<&str> = Vec::new();

    // Iterate through the lines, starting a new record at each CCSDS_OMM_VERS
    for line in omm_kvn_string.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        // A new CCSDS_OMM_VERS line starts the next concatenated record
        if trimmed.to_ascii_uppercase().starts_with("CCSDS_OMM_VERS") && !current_lines.is_empty() {
            let sgp4 = from_omm_kvn_lines(&current_lines);
            sgp4s.push(sgp4);
            current_lines.clear();
        }

        current_lines.push(trimmed);
    }

    // Parse the final record
    if !current_lines.is_empty() {
        let sgp4 = from_omm_kvn_lines(&current_lines);
        sgp4s.push(sgp4);
    }

    // Return vector of SGP4 structs
    return sgp4s;
}

/// Builds a vector of [`Sgp4`] structs from a file containing Orbit Mean-Elements Message
/// (OMM) sets in CCSDS KVN format.
///
/// This function parses a file containing one or more OMM records in KVN
/// format and returns all successfully parsed entries.
///
/// # Arguments
/// * `omm_kvn_file_path` - A path to a file containing one or more OMM records in KVN format
///
/// # Returns
/// * `Vec<Sgp4>` - A vector containing all successfully parsed SGP4 structs
///
/// # Panics
/// * If the file cannot be read
/// * If a present OMM field cannot be parsed
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_kvn_file;
///
/// // Define the OMM KVN file path
/// let omm_kvn_file_path = "test/omm_parsing_cases.txt";
///
/// // Parse the OMM file into a vector of SGP4 structs
/// let sgp4s = from_omm_kvn_file(omm_kvn_file_path);
/// let sgp4 = &sgp4s[0];
///
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_kvn_file(omm_kvn_file_path: &str) -> Vec<Sgp4> {
    // Open the OMM KVN file
    let omm_kvn_string = fs::read_to_string(omm_kvn_file_path).expect("Cannot read OMM KVN file");

    // Parse OMM string into a vector of SGP4 structs
    let sgp4s = from_omm_kvn_string(&omm_kvn_string);

    // Return the vector of SGP4 structs
    return sgp4s;
}

/// Parse a KVN field from OMM record lines
///
/// Search the provided lines for a key-value entry matching the given field
/// name. Missing or empty fields return a type-specific default: 0 for numbers
/// and an empty string for strings. Classification defaults to U.
///
/// # Arguments
/// * `lines` - The KVN lines of one OMM record
/// * `field` - The KVN keyword to look up
///
/// # Returns
/// * The parsed value as the requested type
///
/// # Panics
/// * If a present field cannot be parsed as the requested type
#[cfg(test)]
fn kvn_parse<T: FromOmm>(lines: &[&str], field: &str) -> T {
    return omm_typed_value(kvn_lookup(lines, field), field);
}

/// Look up a KVN keyword and return the cleaned value, if present
///
/// # Arguments
/// * `lines` - The KVN lines of one OMM record
/// * `field` - The KVN keyword to look up
///
/// # Returns
/// * `Some(value)` - The cleaned field value
/// * `None` - If the field is not present
fn kvn_lookup(lines: &[&str], field: &str) -> Option<String> {
    // Search each line for a matching KVN keyword
    for line in lines {
        let line = line.trim();

        // Skip empty lines and COMMENT lines
        if line.is_empty() || line.to_ascii_uppercase().starts_with("COMMENT") {
            continue;
        }

        // Split KEY = VALUE
        let Some((keyword, value)) = line.split_once('=') else {
            continue;
        };
        if !keyword.trim().eq_ignore_ascii_case(field) {
            continue;
        }

        return Some(clean_kvn_value(value));
    }

    return None;
}

/// Convert an optional OMM field string into type T
///
/// Missing or empty values use the type default.
///
/// # Arguments
/// * `value` - The raw field text, if the keyword was present
/// * `field` - The OMM keyword name, used in parse error messages
///
/// # Returns
/// * The parsed value as the requested type
///
/// # Panics
/// * If a present field cannot be parsed as the requested type
fn omm_typed_value<T: FromOmm>(value: Option<String>, field: &str) -> T {
    match value {
        Some(v) if !v.is_empty() => T::from_omm(field, &v),
        _ => T::omm_default(),
    }
}

/// Build an [`Sgp4`] struct from an OMM field lookup
///
/// Shared by KVN and XML parsers. Missing numeric fields default to 0 and
/// missing string fields default to an empty string. OMM MEAN_MOTION_DOT and
/// MEAN_MOTION_DDOT are the TLE-printed values (already divided by 2 and 6)
/// and are scaled so [`GenPerturbElementSet`] stores the true derivatives.
///
/// # Arguments
/// * `lookup` - Returns the text of one OMM keyword, if present
///
/// # Returns
/// * [`Sgp4`] - Struct containing the parsed SGP4 parameters
///
/// # Panics
/// * If a present field cannot be parsed as the requested type
fn sgp4_from_omm_lookup<F>(lookup: F) -> Sgp4
where
    F: Fn(&str) -> Option<String>,
{
    // Create mutable General Perturbation Element Set struct
    let mut gp = GenPerturbElementSet::default();

    // Common name of the satellite
    gp.common_name = omm_typed_value(lookup("OBJECT_NAME"), "OBJECT_NAME");

    // NORAD satellite catalog number
    gp.satellite_catalog_number = omm_typed_value(lookup("NORAD_CAT_ID"), "NORAD_CAT_ID");

    // Classification (default unclassified)
    gp.classification = omm_typed_value(lookup("CLASSIFICATION_TYPE"), "CLASSIFICATION_TYPE");

    // International designator
    gp.international_designator = omm_typed_value(lookup("OBJECT_ID"), "OBJECT_ID");

    // Epoch UTC datetime
    gp.epoch_datetime = omm_typed_value(lookup("EPOCH"), "EPOCH");

    // 1st derivative of mean motion [revs/day^2]
    // OMM stores the TLE-printed value, so multiply by 2
    gp.first_derivative_of_mean_motion =
        omm_typed_value::<f64>(lookup("MEAN_MOTION_DOT"), "MEAN_MOTION_DOT") * 2.0;

    // 2nd derivative of mean motion [revs/day^3]
    // OMM stores the TLE-printed value, so multiply by 6
    gp.second_derivative_of_mean_motion =
        omm_typed_value::<f64>(lookup("MEAN_MOTION_DDOT"), "MEAN_MOTION_DDOT") * 6.0;

    // BSTAR drag term [1/Earth radii]
    gp.bstar = omm_typed_value(lookup("BSTAR"), "BSTAR");

    // Ephemeris type
    gp.ephemeris_type = omm_typed_value(lookup("EPHEMERIS_TYPE"), "EPHEMERIS_TYPE");

    // Element set number
    gp.element_set_number = omm_typed_value(lookup("ELEMENT_SET_NO"), "ELEMENT_SET_NO");

    // Inclination [degs]
    gp.inclination = omm_typed_value(lookup("INCLINATION"), "INCLINATION");

    // Right ascension of ascending node [degs]
    gp.right_ascension_of_ascending_node =
        omm_typed_value(lookup("RA_OF_ASC_NODE"), "RA_OF_ASC_NODE");

    // Eccentricity
    gp.eccentricity = omm_typed_value(lookup("ECCENTRICITY"), "ECCENTRICITY");

    // Argument of perigee [degs]
    gp.argument_of_perigee = omm_typed_value(lookup("ARG_OF_PERICENTER"), "ARG_OF_PERICENTER");

    // Mean anomaly [degs]
    gp.mean_anomaly = omm_typed_value(lookup("MEAN_ANOMALY"), "MEAN_ANOMALY");

    // Mean motion [revs/day]
    gp.mean_motion = omm_typed_value(lookup("MEAN_MOTION"), "MEAN_MOTION");

    // Revolution number at epoch
    gp.revolution_number_at_epoch = omm_typed_value(lookup("REV_AT_EPOCH"), "REV_AT_EPOCH");

    // Initialize the SGP4 parameters
    let sgp4 = init_sgp4(&gp, None);

    return sgp4;
}

/// Clean a KVN value string
///
/// Trim whitespace, strip matching double quotes, and remove optional trailing
/// units written in square brackets.
///
/// # Arguments
/// * `value` - The raw value text to the right of the equals sign
///
/// # Returns
/// * The cleaned value string
fn clean_kvn_value(value: &str) -> String {
    let mut v = value.trim();

    // Strip matching double quotes
    if v.len() >= 2 && v.starts_with('"') && v.ends_with('"') {
        v = &v[1..v.len() - 1];
        v = v.trim();
    }

    // Strip optional units in square brackets
    if let Some(idx) = v.find('[') {
        v = v[..idx].trim();
    }

    return v.to_string();
}

/// Collect leaf OMM keyword values from one XML omm element
///
/// Wrapper tags (metadata, meanElements, tleParameters) and COMMENT elements
/// are skipped. Keys are stored in uppercase so lookup is case-insensitive.
///
/// # Arguments
/// * `omm` - The omm element to flatten
///
/// # Returns
/// * A map of OMM keyword to text value
fn xml_omm_fields(omm: Node) -> HashMap<String, String> {
    let mut fields = HashMap::new();

    for node in omm.descendants() {
        if !node.is_element() {
            continue;
        }

        let name = node.tag_name().name();

        // Skip the record root, comments, and wrapper elements with children
        if name.eq_ignore_ascii_case("omm") || name.eq_ignore_ascii_case("COMMENT") {
            continue;
        }
        if node.children().any(|child| child.is_element()) {
            continue;
        }

        let text = node.text().unwrap_or("").trim();
        fields.insert(name.to_ascii_uppercase(), text.to_string());
    }

    return fields;
}

/// Builds a vector of [`Sgp4`] structs from a string containing Orbit Mean-Elements Message
/// (OMM) sets in CCSDS XML format.
///
/// This function parses a string containing one or more OMM records. Each
/// omm element is flattened to keyword and value pairs by local tag name, then
/// initialized with the same GP mapping as [`from_omm_kvn_lines`]. COMMENT
/// elements and wrapper tags are ignored. Empty tags use the same defaults
/// as missing KVN fields.
///
/// # Arguments
/// * `omm_xml_string` - A string containing one or more OMM records in XML format
///
/// # Returns
/// * `Vec<Sgp4>` - A vector containing all successfully parsed SGP4 parameters
///
/// # Panics
/// * If the XML document cannot be parsed
/// * If a present OMM field cannot be parsed
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_xml_string;
///
/// // Define the OMM XML string
/// let omm = r#"<ndm>
/// <omm id="CCSDS_OMM_VERS" version="2.0">
/// <body><segment>
/// <metadata>
/// <OBJECT_NAME>2026-106A</OBJECT_NAME>
/// <OBJECT_ID>2026-106A</OBJECT_ID>
/// </metadata>
/// <data><meanElements>
/// <EPOCH>2026-06-14T15:07:48.259488</EPOCH>
/// <MEAN_MOTION>15.11169557</MEAN_MOTION>
/// <ECCENTRICITY>.00147468</ECCENTRICITY>
/// <INCLINATION>97.5103</INCLINATION>
/// <RA_OF_ASC_NODE>247.7605</RA_OF_ASC_NODE>
/// <ARG_OF_PERICENTER>169.6213</ARG_OF_PERICENTER>
/// <MEAN_ANOMALY>190.5325</MEAN_ANOMALY>
/// </meanElements><tleParameters>
/// <EPHEMERIS_TYPE>0</EPHEMERIS_TYPE>
/// <CLASSIFICATION_TYPE>U</CLASSIFICATION_TYPE>
/// <NORAD_CAT_ID>69097</NORAD_CAT_ID>
/// <ELEMENT_SET_NO>999</ELEMENT_SET_NO>
/// <REV_AT_EPOCH>459</REV_AT_EPOCH>
/// <BSTAR>.39221734E-3</BSTAR>
/// <MEAN_MOTION_DOT>.6535E-4</MEAN_MOTION_DOT>
/// <MEAN_MOTION_DDOT>0</MEAN_MOTION_DDOT>
/// </tleParameters></data>
/// </segment></body></omm></ndm>"#;
///
/// // Parse the OMM XML into a SGP4 struct
/// let sgp4s = from_omm_xml_string(omm);
/// let sgp4 = &sgp4s[0];
///
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_xml_string(omm_xml_string: &str) -> Vec<Sgp4> {
    // Parse the XML document
    let doc = Document::parse(omm_xml_string).unwrap_or_else(|err| {
        panic!("Cannot parse OMM XML: {}", err);
    });

    // Each omm element is one GP record
    let mut sgp4s = Vec::new();
    for node in doc.descendants() {
        if !node.is_element() {
            continue;
        }
        if !node.tag_name().name().eq_ignore_ascii_case("omm") {
            continue;
        }

        let fields = xml_omm_fields(node);
        let sgp4 = sgp4_from_omm_lookup(|field| fields.get(&field.to_ascii_uppercase()).cloned());
        sgp4s.push(sgp4);
    }

    // Return vector of SGP4 structs
    return sgp4s;
}

/// Builds a vector of [`Sgp4`] structs from a file containing Orbit Mean-Elements Message
/// (OMM) sets in CCSDS XML format.
///
/// This function parses a file containing one or more OMM records in XML
/// format and returns all successfully parsed entries.
///
/// # Arguments
/// * `omm_xml_file_path` - A path to a file containing one or more OMM records in XML format
///
/// # Returns
/// * `Vec<Sgp4>` - A vector containing all successfully parsed SGP4 structs
///
/// # Panics
/// * If the file cannot be read
/// * If the XML document cannot be parsed
/// * If a present OMM field cannot be parsed
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_xml_file;
///
/// // Define the OMM XML file path
/// let omm_xml_file_path = "test/omm_parsing_cases.xml";
///
/// // Parse the OMM file into a vector of SGP4 structs
/// let sgp4s = from_omm_xml_file(omm_xml_file_path);
/// let sgp4 = &sgp4s[0];
///
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_xml_file(omm_xml_file_path: &str) -> Vec<Sgp4> {
    // Open the OMM XML file
    let omm_xml_string = fs::read_to_string(omm_xml_file_path).expect("Cannot read OMM XML file");

    // Parse OMM string into a vector of SGP4 structs
    let sgp4s = from_omm_xml_string(&omm_xml_string);

    // Return the vector of SGP4 structs
    return sgp4s;
}

/// Build an [`Sgp4`] struct from one JSON OMM object
///
/// # Arguments
/// * `record` - One OMM record as a JSON object
///
/// # Returns
/// * [`Sgp4`] - Struct containing the parsed SGP4 parameters
///
/// # Panics
/// * If the JSON value is not an object
/// * If a present OMM field cannot be parsed
fn sgp4_from_json_record(record: &Value) -> Sgp4 {
    let fields = json_omm_fields(record);
    return sgp4_from_omm_lookup(|field| fields.get(&field.to_ascii_uppercase()).cloned());
}

/// Builds a vector of [`Sgp4`] structs from a string containing Orbit Mean-Elements Message
/// (OMM) sets in JSON format.
///
/// This function parses a string containing one OMM object or an array of OMM
/// objects. Each object is flattened to keyword and value pairs, then
/// initialized with the same GP mapping as [`from_omm_kvn_lines`]. JSON numbers
/// are converted to text. Null, missing, and empty string fields use the same
/// defaults as missing KVN fields.
///
/// # Arguments
/// * `omm_json_string` - A string containing one or more OMM records in JSON format
///
/// # Returns
/// * `Vec<Sgp4>` - A vector containing all successfully parsed SGP4 parameters
///
/// # Panics
/// * If the JSON document cannot be parsed
/// * If the top-level value is not an object or an array of objects
/// * If a present OMM field cannot be parsed
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_json_string;
///
/// // Define the OMM JSON string
/// let omm = r#"[
/// {
/// "OBJECT_NAME": "2026-106A",
/// "OBJECT_ID": "2026-106A",
/// "EPOCH": "2026-06-14T15:07:48.259488",
/// "MEAN_MOTION": 15.11169557,
/// "ECCENTRICITY": 0.00147468,
/// "INCLINATION": 97.5103,
/// "RA_OF_ASC_NODE": 247.7605,
/// "ARG_OF_PERICENTER": 169.6213,
/// "MEAN_ANOMALY": 190.5325,
/// "EPHEMERIS_TYPE": 0,
/// "CLASSIFICATION_TYPE": "U",
/// "NORAD_CAT_ID": 69097,
/// "ELEMENT_SET_NO": 999,
/// "REV_AT_EPOCH": 459,
/// "BSTAR": 0.00039221734,
/// "MEAN_MOTION_DOT": 6.535e-05,
/// "MEAN_MOTION_DDOT": 0
/// }
/// ]"#;
///
/// // Parse the OMM JSON into a SGP4 struct
/// let sgp4s = from_omm_json_string(omm);
/// let sgp4 = &sgp4s[0];
///
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_json_string(omm_json_string: &str) -> Vec<Sgp4> {
    // Parse the JSON document
    let value: Value = serde_json::from_str(omm_json_string).unwrap_or_else(|err| {
        panic!("Cannot parse OMM JSON: {}", err);
    });

    // A file may be one object or an array of objects
    let mut sgp4s = Vec::new();
    match &value {
        Value::Array(records) => {
            for record in records {
                sgp4s.push(sgp4_from_json_record(record));
            }
        }
        Value::Object(_) => {
            sgp4s.push(sgp4_from_json_record(&value));
        }
        _ => panic!("OMM JSON must be an object or an array of objects"),
    }

    // Return vector of SGP4 structs
    return sgp4s;
}

/// Collect OMM keyword values from one JSON object
///
/// Keys are stored in uppercase so lookup is case-insensitive. Strings are
/// used as-is, numbers are converted to text, and null becomes an empty
/// string. Nested arrays and objects are ignored.
///
/// # Arguments
/// * `record` - The JSON object to flatten
///
/// # Returns
/// * A map of OMM keyword to text value
///
/// # Panics
/// * If the JSON value is not an object
fn json_omm_fields(record: &Value) -> HashMap<String, String> {
    let Some(object) = record.as_object() else {
        panic!("OMM JSON record must be an object");
    };

    let mut fields = HashMap::new();
    for (key, value) in object {
        if key.eq_ignore_ascii_case("COMMENT") {
            continue;
        }

        // Convert JSON scalars to the same text form used by KVN and XML
        let text = match value {
            Value::Null => String::new(),
            Value::String(s) => s.trim().to_string(),
            Value::Number(n) => n.to_string(),
            Value::Bool(b) => b.to_string(),
            Value::Array(_) | Value::Object(_) => continue,
        };
        fields.insert(key.to_ascii_uppercase(), text);
    }

    return fields;
}

/// Builds a vector of [`Sgp4`] structs from a file containing Orbit Mean-Elements Message
/// (OMM) sets in JSON format.
///
/// This function parses a file containing one OMM object or an array of OMM
/// objects in JSON format and returns all successfully parsed entries.
///
/// # Arguments
/// * `omm_json_file_path` - A path to a file containing one or more OMM records in JSON format
///
/// # Returns
/// * `Vec<Sgp4>` - A vector containing all successfully parsed SGP4 structs
///
/// # Panics
/// * If the file cannot be read
/// * If the JSON document cannot be parsed
/// * If the top-level value is not an object or an array of objects
/// * If a present OMM field cannot be parsed
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_json_file;
///
/// // Define the OMM JSON file path
/// let omm_json_file_path = "test/omm_parsing_cases.json";
///
/// // Parse the OMM file into a vector of SGP4 structs
/// let sgp4s = from_omm_json_file(omm_json_file_path);
/// let sgp4 = &sgp4s[0];
///
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_json_file(omm_json_file_path: &str) -> Vec<Sgp4> {
    // Open the OMM JSON file
    let omm_json_string = fs::read_to_string(omm_json_file_path)
        .expect("Cannot read OMM JSON file");

    // Parse OMM string into a vector of SGP4 structs
    let sgp4s = from_omm_json_string(&omm_json_string);

    // Return the vector of SGP4 structs
    return sgp4s;
}

/// Build an [`Sgp4`] struct from one CSV OMM row
///
/// # Arguments
/// * `headers` - Uppercase OMM keywords from the CSV header row
/// * `record` - One OMM data row
///
/// # Returns
/// * [`Sgp4`] - Struct containing the parsed SGP4 parameters
///
/// # Panics
/// * If a present OMM field cannot be parsed
fn sgp4_from_csv_record(headers: &[String], record: &csv::StringRecord) -> Sgp4 {
    let fields = csv_omm_fields(headers, record);
    return sgp4_from_omm_lookup(|field| fields.get(&field.to_ascii_uppercase()).cloned());
}

/// Builds a vector of [`Sgp4`] structs from a string containing Orbit Mean-Elements Message
/// (OMM) sets in CSV format.
///
/// This function parses a string containing a header row of OMM keywords and
/// one data row per record. Each row is flattened to keyword and value pairs,
/// then initialized with the same GP mapping as [`from_omm_kvn_lines`]. Empty
/// cells use the same defaults as missing KVN fields.
///
/// # Arguments
/// * `omm_csv_string` - A string containing one or more OMM records in CSV format
///
/// # Returns
/// * `Vec<Sgp4>` - A vector containing all successfully parsed SGP4 parameters
///
/// # Panics
/// * If the CSV document cannot be parsed
/// * If a present OMM field cannot be parsed
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_csv_string;
///
/// // Define the OMM CSV string
/// let omm = "\
/// OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,EPHEMERIS_TYPE,CLASSIFICATION_TYPE,NORAD_CAT_ID,ELEMENT_SET_NO,REV_AT_EPOCH,BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT
/// 2026-106A,2026-106A,2026-06-14T15:07:48.259488,15.11169557,.00147468,97.5103,247.7605,169.6213,190.5325,0,U,69097,999,459,.39221734E-3,.6535E-4,0
/// ";
///
/// // Parse the OMM CSV into a SGP4 struct
/// let sgp4s = from_omm_csv_string(omm);
/// let sgp4 = &sgp4s[0];
///
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_csv_string(omm_csv_string: &str) -> Vec<Sgp4> {
    // Parse the CSV document, using the first row as OMM keywords
    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(omm_csv_string.as_bytes());

    let headers: Vec<String> = reader
        .headers()
        .unwrap_or_else(|err| {
            panic!("Cannot parse OMM CSV headers: {}", err);
        })
        .iter()
        .map(|header| header.trim().to_ascii_uppercase())
        .collect();

    // Each data row is one GP record
    let mut sgp4s = Vec::new();
    for result in reader.records() {
        let record = result.unwrap_or_else(|err| {
            panic!("Cannot parse OMM CSV: {}", err);
        });
        sgp4s.push(sgp4_from_csv_record(&headers, &record));
    }

    // Return vector of SGP4 structs
    return sgp4s;
}

/// Collect OMM keyword values from one CSV row
///
/// Keys are stored in uppercase so lookup is case-insensitive. COMMENT columns
/// are skipped. Empty cells become empty strings so missing-field defaults apply.
///
/// # Arguments
/// * `headers` - Uppercase OMM keywords from the CSV header row
/// * `record` - One OMM data row
///
/// # Returns
/// * A map of OMM keyword to text value
fn csv_omm_fields(headers: &[String], record: &csv::StringRecord) -> HashMap<String, String> {
    let mut fields = HashMap::new();

    for (index, header) in headers.iter().enumerate() {
        if header.is_empty() || header.eq_ignore_ascii_case("COMMENT") {
            continue;
        }

        let text = record.get(index).unwrap_or("").trim().to_string();
        fields.insert(header.clone(), text);
    }

    return fields;
}

/// Builds a vector of [`Sgp4`] structs from a file containing Orbit Mean-Elements Message
/// (OMM) sets in CSV format.
///
/// This function parses a file containing a header row of OMM keywords and
/// one data row per record and returns all successfully parsed entries.
///
/// # Arguments
/// * `omm_csv_file_path` - A path to a file containing one or more OMM records in CSV format
///
/// # Returns
/// * `Vec<Sgp4>` - A vector containing all successfully parsed SGP4 structs
///
/// # Panics
/// * If the file cannot be read
/// * If the CSV document cannot be parsed
/// * If a present OMM field cannot be parsed
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_csv_file;
///
/// // Define the OMM CSV file path
/// let omm_csv_file_path = "test/omm_parsing_cases.csv";
///
/// // Parse the OMM file into a vector of SGP4 structs
/// let sgp4s = from_omm_csv_file(omm_csv_file_path);
/// let sgp4 = &sgp4s[0];
///
/// // Assert the SGP4 struct is correct
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_csv_file(omm_csv_file_path: &str) -> Vec<Sgp4> {
    // Open the OMM CSV file
    let omm_csv_string = fs::read_to_string(omm_csv_file_path)
        .expect("Cannot read OMM CSV file");

    // Parse OMM string into a vector of SGP4 structs
    let sgp4s = from_omm_csv_string(&omm_csv_string);

    // Return the vector of SGP4 structs
    return sgp4s;
}

/// Parse an OMM EPOCH string into a UTC DateTime
///
/// Accepts an ISO-8601 date and time with optional fractional seconds and an
/// optional trailing Z. The date and time may be separated by T or a space.
///
/// # Arguments
/// * `epoch` - The EPOCH string from the OMM record
///
/// # Returns
/// * [`DateTime`] - The epoch as a UTC datetime
///
/// # Panics
/// * If the EPOCH string cannot be parsed
fn parse_omm_epoch(epoch: &str) -> DateTime {
    // Strip a trailing Z timezone marker
    let mut s = epoch.trim();
    if s.ends_with('Z') || s.ends_with('z') {
        s = &s[..s.len() - 1];
    }

    // Split date and time
    let (date, time) = if let Some((date, time)) = s.split_once('T') {
        (date, time)
    } else if let Some((date, time)) = s.split_once(' ') {
        (date, time)
    } else {
        panic!("OMM EPOCH is invalid: {}", epoch);
    };

    let date_parts: Vec<&str> = date.split('-').collect();
    let time_parts: Vec<&str> = time.split(':').collect();
    if date_parts.len() != 3 || time_parts.len() != 3 {
        panic!("OMM EPOCH is invalid: {}", epoch);
    }

    // Parse calendar date
    let year = date_parts[0].parse::<i32>().unwrap_or_else(|_| {
        panic!("OMM EPOCH year is invalid: {}", epoch);
    });
    let month = date_parts[1].parse::<i32>().unwrap_or_else(|_| {
        panic!("OMM EPOCH month is invalid: {}", epoch);
    });
    let day = date_parts[2].parse::<i32>().unwrap_or_else(|_| {
        panic!("OMM EPOCH day is invalid: {}", epoch);
    });

    // Parse clock time
    let hour = time_parts[0].parse::<i32>().unwrap_or_else(|_| {
        panic!("OMM EPOCH hour is invalid: {}", epoch);
    });
    let minute = time_parts[1].parse::<i32>().unwrap_or_else(|_| {
        panic!("OMM EPOCH minute is invalid: {}", epoch);
    });
    let second = f64::from_omm("EPOCH second", time_parts[2]);

    return DateTime {
        year,
        month,
        day,
        hour,
        minute,
        second,
        timezone: Timezone::UTC,
    };
}

// ----------
// Unit Tests
// ----------
#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;
    use std::collections::HashMap;
    use toml::from_str;

    #[test]
    fn test_tle_full_year() {
        assert_eq!(tle_full_year(0), 2000);
        assert_eq!(tle_full_year(25), 2025);
        assert_eq!(tle_full_year(56), 2056);
        assert_eq!(tle_full_year(57), 1957);
        assert_eq!(tle_full_year(98), 1998);
        assert_eq!(tle_full_year(99), 1999);
    }

    #[test]
    fn test_alpha5_digit() {
        assert_eq!(alpha5_digit('A'), Some(10));
        assert_eq!(alpha5_digit('H'), Some(17));
        assert_eq!(alpha5_digit('J'), Some(18));
        assert_eq!(alpha5_digit('N'), Some(22));
        assert_eq!(alpha5_digit('P'), Some(23));
        assert_eq!(alpha5_digit('Z'), Some(33));
        assert_eq!(alpha5_digit('a'), Some(10));
        assert_eq!(alpha5_digit('I'), None);
        assert_eq!(alpha5_digit('O'), None);
        assert_eq!(alpha5_digit('0'), None);
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

    // -------------------------------------------------------
    // Structs for deserializing the TLE / OMM parsing test TOML
    // -------------------------------------------------------

    #[derive(Deserialize)]
    struct ParsingCases {
        test: HashMap<String, ParsingCase>,
    }

    #[derive(Deserialize)]
    struct ParsingCase {
        name: String,
        #[serde(default)]
        tle: String,
        #[serde(default)]
        omm_kvn: String,
        #[serde(default)]
        exception: bool,
        #[serde(default)]
        common_name: String,
        #[serde(default)]
        satellite_catalog_number: i32,
        #[serde(default)]
        classification: String,
        #[serde(default)]
        international_designator: String,
        #[serde(default)]
        epoch_datetime_year: i32,
        #[serde(default)]
        epoch_datetime_month: i32,
        #[serde(default)]
        epoch_datetime_day: i32,
        #[serde(default)]
        epoch_datetime_hour: i32,
        #[serde(default)]
        epoch_datetime_minute: i32,
        #[serde(default)]
        epoch_datetime_second: f64,
        #[serde(default)]
        first_derivative_of_mean_motion: f64,
        #[serde(default)]
        second_derivative_of_mean_motion: f64,
        #[serde(default)]
        bstar: f64,
        #[serde(default)]
        ephemeris_type: i32,
        #[serde(default)]
        element_set_number: i32,
        #[serde(default)]
        inclination: f64,
        #[serde(default)]
        right_ascension_of_ascending_node: f64,
        #[serde(default)]
        eccentricity: f64,
        #[serde(default)]
        argument_of_perigee: f64,
        #[serde(default)]
        mean_anomaly: f64,
        #[serde(default)]
        mean_motion: f64,
        #[serde(default)]
        revolution_number_at_epoch: i64,
    }

    const TLE_PARSE_TOL: f64 = 1e-9;

    fn assert_near(key: &str, name: &str, source: &str, label: &str, value: f64, expected: f64) {
        assert!(
            (value - expected).abs() < TLE_PARSE_TOL,
            "{key} ({name}) via {source}: {label} mismatch {value} vs {expected}"
        );
    }

    fn assert_gp_matches(
        key: &str,
        name: &str,
        source: &str,
        gp: &GenPerturbElementSet,
        case: &ParsingCase,
    ) {
        assert_eq!(
            gp.common_name, case.common_name,
            "{key} ({name}) via {source}: common_name"
        );
        assert_eq!(
            gp.satellite_catalog_number, case.satellite_catalog_number,
            "{key} ({name}) via {source}: satellite_catalog_number"
        );
        assert_eq!(
            gp.classification.to_string(),
            case.classification,
            "{key} ({name}) via {source}: classification"
        );
        assert_eq!(
            gp.international_designator, case.international_designator,
            "{key} ({name}) via {source}: international_designator"
        );
        assert_eq!(
            gp.epoch_datetime.year, case.epoch_datetime_year,
            "{key} ({name}) via {source}: epoch year"
        );
        assert_eq!(
            gp.epoch_datetime.month, case.epoch_datetime_month,
            "{key} ({name}) via {source}: epoch month"
        );
        assert_eq!(
            gp.epoch_datetime.day, case.epoch_datetime_day,
            "{key} ({name}) via {source}: epoch day"
        );
        assert_eq!(
            gp.epoch_datetime.hour, case.epoch_datetime_hour,
            "{key} ({name}) via {source}: epoch hour"
        );
        assert_eq!(
            gp.epoch_datetime.minute, case.epoch_datetime_minute,
            "{key} ({name}) via {source}: epoch minute"
        );
        assert_near(
            key,
            name,
            source,
            "epoch second",
            gp.epoch_datetime.second,
            case.epoch_datetime_second,
        );
        assert_near(
            key,
            name,
            source,
            "n-dot",
            gp.first_derivative_of_mean_motion,
            case.first_derivative_of_mean_motion,
        );
        assert_near(
            key,
            name,
            source,
            "n-ddot",
            gp.second_derivative_of_mean_motion,
            case.second_derivative_of_mean_motion,
        );
        assert_near(key, name, source, "bstar", gp.bstar, case.bstar);
        assert_eq!(
            gp.ephemeris_type, case.ephemeris_type,
            "{key} ({name}) via {source}: ephemeris_type"
        );
        assert_eq!(
            gp.element_set_number, case.element_set_number,
            "{key} ({name}) via {source}: element_set_number"
        );
        assert_near(
            key,
            name,
            source,
            "inclination",
            gp.inclination,
            case.inclination,
        );
        assert_near(
            key,
            name,
            source,
            "raan",
            gp.right_ascension_of_ascending_node,
            case.right_ascension_of_ascending_node,
        );
        assert_near(
            key,
            name,
            source,
            "eccentricity",
            gp.eccentricity,
            case.eccentricity,
        );
        assert_near(
            key,
            name,
            source,
            "argp",
            gp.argument_of_perigee,
            case.argument_of_perigee,
        );
        assert_near(
            key,
            name,
            source,
            "mean anomaly",
            gp.mean_anomaly,
            case.mean_anomaly,
        );
        assert_near(
            key,
            name,
            source,
            "mean motion",
            gp.mean_motion,
            case.mean_motion,
        );
        assert_eq!(
            gp.revolution_number_at_epoch, case.revolution_number_at_epoch,
            "{key} ({name}) via {source}: revolution_number_at_epoch"
        );
    }

    fn sgp4_from_case_tle(tle: &str) -> Sgp4 {
        let lines: Vec<&str> = tle
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .collect();
        if lines[0].starts_with('1') {
            from_tle_lines(lines[0], lines[1], None)
        } else {
            from_tle_lines(lines[1], lines[2], Some(lines[0]))
        }
    }

    #[test]
    fn test_tle_parsing_cases() {
        let content = std::fs::read_to_string("test/tle_parsing_cases.toml")
            .expect("could not read test/tle_parsing_cases.toml");
        let cases: ParsingCases =
            from_str(&content).expect("could not parse test/tle_parsing_cases.toml");

        let from_file = from_tle_file("test/tle_parsing_cases.txt");
        let expected_file_count = cases.test.values().filter(|c| !c.exception).count();
        assert_eq!(
            from_file.len(),
            expected_file_count,
            "tle_parsing_cases.txt should contain one entry per non-error TOML case"
        );

        let mut keys: Vec<&String> = cases.test.keys().collect();
        keys.sort();

        for key in keys {
            let case = &cases.test[key];

            if case.exception {
                let lines_result = std::panic::catch_unwind(|| sgp4_from_case_tle(&case.tle));
                assert!(
                    lines_result.is_err(),
                    "case {key} ({}): expected from_tle_lines to panic",
                    case.name
                );
                let string_result = std::panic::catch_unwind(|| from_tle_string(&case.tle));
                assert!(
                    string_result.is_err(),
                    "case {key} ({}): expected from_tle_string to panic",
                    case.name
                );
                continue;
            }

            let from_lines = sgp4_from_case_tle(&case.tle);
            assert_gp_matches(key, &case.name, "from_tle_lines", &from_lines.gp, case);

            let from_string = from_tle_string(&case.tle);
            assert_eq!(
                from_string.len(),
                1,
                "case {key} ({}): expected one TLE from string",
                case.name
            );
            assert_gp_matches(key, &case.name, "from_tle_string", &from_string[0].gp, case);

            let file_match = from_file
                .iter()
                .find(|s| s.gp.satellite_catalog_number == case.satellite_catalog_number)
                .unwrap_or_else(|| {
                    panic!(
                        "case {key} ({}): catalog {} missing from tle_parsing_cases.txt parse",
                        case.name, case.satellite_catalog_number
                    )
                });
            assert_gp_matches(key, &case.name, "from_tle_file", &file_match.gp, case);
        }
    }

    fn sgp4_from_case_omm(omm_kvn: &str) -> Sgp4 {
        let lines: Vec<&str> = omm_kvn
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .collect();
        from_omm_kvn_lines(&lines)
    }

    #[test]
    fn test_omm_parsing_cases() {
        let content = std::fs::read_to_string("test/omm_parsing_cases.toml")
            .expect("could not read test/omm_parsing_cases.toml");
        let cases: ParsingCases =
            from_str(&content).expect("could not parse test/omm_parsing_cases.toml");

        let from_file = from_omm_kvn_file("test/omm_parsing_cases.txt");
        let from_xml_file = from_omm_xml_file("test/omm_parsing_cases.xml");
        let from_json_file = from_omm_json_file("test/omm_parsing_cases.json");
        let from_csv_file = from_omm_csv_file("test/omm_parsing_cases.csv");
        let expected_file_count = cases.test.values().filter(|c| !c.exception).count();
        assert_eq!(
            from_file.len(),
            expected_file_count,
            "omm_parsing_cases.txt should contain one entry per non-error TOML case"
        );
        assert_eq!(
            from_xml_file.len(),
            expected_file_count,
            "omm_parsing_cases.xml should contain one entry per non-error TOML case"
        );
        assert_eq!(
            from_json_file.len(),
            expected_file_count,
            "omm_parsing_cases.json should contain one entry per non-error TOML case"
        );
        assert_eq!(
            from_csv_file.len(),
            expected_file_count,
            "omm_parsing_cases.csv should contain one entry per non-error TOML case"
        );

        let mut keys: Vec<&String> = cases.test.keys().collect();
        keys.sort();

        for key in keys {
            let case = &cases.test[key];

            if case.exception {
                let lines_result = std::panic::catch_unwind(|| sgp4_from_case_omm(&case.omm_kvn));
                assert!(
                    lines_result.is_err(),
                    "case {key} ({}): expected from_omm_kvn_lines to panic",
                    case.name
                );
                let string_result = std::panic::catch_unwind(|| from_omm_kvn_string(&case.omm_kvn));
                assert!(
                    string_result.is_err(),
                    "case {key} ({}): expected from_omm_kvn_string to panic",
                    case.name
                );
                continue;
            }

            let from_lines = sgp4_from_case_omm(&case.omm_kvn);
            assert_gp_matches(key, &case.name, "from_omm_kvn_lines", &from_lines.gp, case);

            let from_string = from_omm_kvn_string(&case.omm_kvn);
            assert_eq!(
                from_string.len(),
                1,
                "case {key} ({}): expected one OMM from string",
                case.name
            );
            assert_gp_matches(
                key,
                &case.name,
                "from_omm_kvn_string",
                &from_string[0].gp,
                case,
            );

            let file_match = from_file
                .iter()
                .find(|s| s.gp.satellite_catalog_number == case.satellite_catalog_number)
                .unwrap_or_else(|| {
                    panic!(
                        "case {key} ({}): catalog {} missing from omm_parsing_cases.txt parse",
                        case.name, case.satellite_catalog_number
                    )
                });
            assert_gp_matches(key, &case.name, "from_omm_kvn_file", &file_match.gp, case);

            let xml_match = from_xml_file
                .iter()
                .find(|s| s.gp.satellite_catalog_number == case.satellite_catalog_number)
                .unwrap_or_else(|| {
                    panic!(
                        "case {key} ({}): catalog {} missing from omm_parsing_cases.xml parse",
                        case.name, case.satellite_catalog_number
                    )
                });
            assert_gp_matches(key, &case.name, "from_omm_xml_file", &xml_match.gp, case);

            let json_match = from_json_file
                .iter()
                .find(|s| s.gp.satellite_catalog_number == case.satellite_catalog_number)
                .unwrap_or_else(|| {
                    panic!(
                        "case {key} ({}): catalog {} missing from omm_parsing_cases.json parse",
                        case.name, case.satellite_catalog_number
                    )
                });
            assert_gp_matches(key, &case.name, "from_omm_json_file", &json_match.gp, case);

            let csv_match = from_csv_file
                .iter()
                .find(|s| s.gp.satellite_catalog_number == case.satellite_catalog_number)
                .unwrap_or_else(|| {
                    panic!(
                        "case {key} ({}): catalog {} missing from omm_parsing_cases.csv parse",
                        case.name, case.satellite_catalog_number
                    )
                });
            assert_gp_matches(key, &case.name, "from_omm_csv_file", &csv_match.gp, case);
        }
    }

    #[test]
    fn test_kvn_parse_defaults() {
        // Define an empty OMM record
        let lines: [&str; 0] = [];

        // Parse missing fields
        let name: String = kvn_parse(&lines, "OBJECT_NAME");
        let catalog: i32 = kvn_parse(&lines, "NORAD_CAT_ID");
        let rev: i64 = kvn_parse(&lines, "REV_AT_EPOCH");
        let motion: f64 = kvn_parse(&lines, "MEAN_MOTION");
        let class: char = kvn_parse(&lines, "CLASSIFICATION_TYPE");

        // Assert numeric fields default to 0, strings to empty, classification to U
        assert_eq!(name, "");
        assert_eq!(catalog, 0);
        assert_eq!(rev, 0);
        assert_eq!(motion, 0.0);
        assert_eq!(class, 'U');
    }
}
