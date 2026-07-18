# MAKO-SGP4
**UNDER CONSTRUCTION**

A Rust crate to parse and propagate Two-Line Element (TLE) sets using Simplified Perturbations Models (SGP4 / SDP4).

## Motivations
I am pursuing this project for two reasons
1. To learn to write Rust code
2. To dive into the theory surrounding SGP4

## Status
- TLE parsing (optional name line + lines 1 and 2)
- Near-Earth SGP4 and deep-space SDP4 (Hoots-style formulation, with Vallado Rev 3 corrections where needed)
- Verification against Vallado Rev 3 cases in `test/vallado_cases.toml`

Known open item: case 21 (`e < 1e-4` drag-term drop) is not implemented yet.

## Plan / TODOs
- OMM / GP dataset handling
- Write a math spec
- Fit data to a TLE
- Conjunction screening

## Testing and Documentation
```bash
# To run the unit tests
cargo test

# To build the Rust Docs
cargo doc
```

The integration tests in `test/vallado_cases.toml` use reference TLEs and ephemerides reproduced from Vallado et al., "Revisiting Spacetrack Report #3: Rev 3," AIAA 2006-6753, 2006. Reproduced for algorithm verification purposes.

## Resources
- [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
- [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
- [Fundamentals of Astrodynamics Github Repository by Vallado et al](https://github.com/CelesTrak/fundamentals-of-astrodynamics?tab=readme-ov-file)
- [History of Analytical Orbit Modeling in the U. S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
