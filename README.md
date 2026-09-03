# MAKO-SGP4
**UNDER CONSTRUCTION**

A Rust crate to parse and propagate General Perturbation Element Sets (GPs) using Simplified Perturbations Models (SGP4 / SDP4). This code was implemented using the theory and equations found in ***History of Analytical Orbit Modeling in the U. S. Space Surveillance System*** by Hoots et al. Practical implementation adjustments were made based on ***Revisiting Spacetrack Report #3: Rev 3*** by Vallado et al.

The name MAKO-SGP4 pays tribute to the Shortfin Mako Shark, the fastest shark species. The speed and efficiency of the SGP4 propagator make this an apt name.

# Usage
In progress

## Accuracy
In progress

## Mathematical Specification
In progress

## Plan / TODOs
- Package as a crate
    - Return error instead of panic
    - Straighten out README
    - Add examples to docs
    - Write math spec
- Fit state data to GP
- Python wrapper

## Testing and Documentation
```bash
# To run the unit tests
cargo test

# To build the Rust Docs
cargo doc
```

## Resources
- [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
- [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
- [Fundamentals of Astrodynamics Github Repository by Vallado et al](https://github.com/CelesTrak/fundamentals-of-astrodynamics?tab=readme-ov-file)
- [History of Analytical Orbit Modeling in the U. S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
- [Space-Track](https://www.space-track.org/auth/login)
- [Celestrak](https://celestrak.org/)
