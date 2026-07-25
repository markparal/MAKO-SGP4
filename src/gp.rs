//! Module to handle the input and processing of GP (General Perturbation)
//! elements. Element sets can be parsed from TLE (including Alpha-5 catalog
//! numbers) and OMM formats.

// ------------------
// External Libraries
// ------------------
use std::fs;

// ------------------
// Internal Libraries
// ------------------
use crate::time::{dayofyr2utc, DateTime};
use crate::sgp4::{init_sgp4, Sgp4};

// -------
// Structs
// -------

/// General Perturbation Element Set for an Earth-orbiting satellite.
///
/// This struct represents the parsed contents of a standard General Perturbation Element Set.
/// The General Perturbation Element Set is the standard set of orbital elements used with the
/// SGP4 propagator. Elements are commonly distributed as Two-Line Element (TLE) text; other
/// formats (OMM, GP JSON/CSV) use the same fields but are not parsed by this crate yet.
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

// -----
// Enums
// -----

// ---------
// Constants
// ---------

// ---------
// Functions
// ---------

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
        let catalog = line1[2..7].trim();
        let first_char = catalog.chars().next().unwrap();
        if first_char.is_ascii_digit() {
            // This is a classic numeric catalog number
            gp.satellite_catalog_number = catalog.parse::<i32>().unwrap();
        } else if first_char.is_ascii_alphabetic() {
            // This is an Alpha-5 catalog number
            gp.satellite_catalog_number = alpha5_digit(first_char).unwrap() * 10000 + catalog[1..].parse::<i32>().unwrap();
        } else {
            // This is an unexpected catalog number
            panic!("TLE line 1 is invalid: catalog number is invalid: {}", catalog);
        }

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

// ----------
// Unit Tests
// ----------
#[cfg(test)]
mod tests {
    use super::*;
    use toml::from_str;
    use std::collections::HashMap;
    use serde::Deserialize;

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
    // Structs for deserializing the TLE parsing test TOML
    // -------------------------------------------------------

    #[derive(Deserialize)]
    struct TleParsingCases {
        test: HashMap<String, TleParsingCase>,
    }

    #[derive(Deserialize)]
    struct TleParsingCase {
        name: String,
        tle: String,
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

    fn assert_gp_matches(key: &str, name: &str, source: &str, gp: &GenPerturbElementSet, case: &TleParsingCase) {
        assert_eq!(gp.common_name, case.common_name, "{key} ({name}) via {source}: common_name");
        assert_eq!(
            gp.satellite_catalog_number, case.satellite_catalog_number,
            "{key} ({name}) via {source}: satellite_catalog_number"
        );
        assert_eq!(
            gp.classification.to_string(), case.classification,
            "{key} ({name}) via {source}: classification"
        );
        assert_eq!(
            gp.international_designator, case.international_designator,
            "{key} ({name}) via {source}: international_designator"
        );
        assert_eq!(gp.epoch_datetime.year, case.epoch_datetime_year, "{key} ({name}) via {source}: epoch year");
        assert_eq!(gp.epoch_datetime.month, case.epoch_datetime_month, "{key} ({name}) via {source}: epoch month");
        assert_eq!(gp.epoch_datetime.day, case.epoch_datetime_day, "{key} ({name}) via {source}: epoch day");
        assert_eq!(gp.epoch_datetime.hour, case.epoch_datetime_hour, "{key} ({name}) via {source}: epoch hour");
        assert_eq!(gp.epoch_datetime.minute, case.epoch_datetime_minute, "{key} ({name}) via {source}: epoch minute");
        assert_near(key, name, source, "epoch second", gp.epoch_datetime.second, case.epoch_datetime_second);
        assert_near(
            key, name, source, "n-dot",
            gp.first_derivative_of_mean_motion, case.first_derivative_of_mean_motion,
        );
        assert_near(
            key, name, source, "n-ddot",
            gp.second_derivative_of_mean_motion, case.second_derivative_of_mean_motion,
        );
        assert_near(key, name, source, "bstar", gp.bstar, case.bstar);
        assert_eq!(gp.ephemeris_type, case.ephemeris_type, "{key} ({name}) via {source}: ephemeris_type");
        assert_eq!(gp.element_set_number, case.element_set_number, "{key} ({name}) via {source}: element_set_number");
        assert_near(key, name, source, "inclination", gp.inclination, case.inclination);
        assert_near(
            key, name, source, "raan",
            gp.right_ascension_of_ascending_node, case.right_ascension_of_ascending_node,
        );
        assert_near(key, name, source, "eccentricity", gp.eccentricity, case.eccentricity);
        assert_near(key, name, source, "argp", gp.argument_of_perigee, case.argument_of_perigee);
        assert_near(key, name, source, "mean anomaly", gp.mean_anomaly, case.mean_anomaly);
        assert_near(key, name, source, "mean motion", gp.mean_motion, case.mean_motion);
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
        let cases: TleParsingCases = from_str(&content)
            .expect("could not parse test/tle_parsing_cases.toml");

        let from_file = from_tle_file("test/tle_parsing_cases.txt");
        let expected_file_count = cases
            .test
            .values()
            .filter(|c| !c.exception)
            .count();
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
}