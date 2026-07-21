# ascent 0.1.1

## CRAN Resubmission Updates

* Expanded ASC-CFD acronym in `DESCRIPTION` to 'ASC-CFD' (Assemblage Shift Characterization - Community Functional Dynamics).
* Added missing `\value` tag in `plot.ascfcd_entities.Rd`.
* Refactored random number generation (RNG) in `asc_null()` to strictly comply with CRAN policies regarding `.GlobalEnv` and `.Random.seed` manipulation.

# ascent 0.1.0
## Initial CRAN release

### Core Framework
* `asc_paired()`: Temporal/paired functional restructuring (Layer 1-3).
* `asc_pairwise()`: Bidirectional spatial functional divergence across all community pairs.
* `asc_baseline()`: Absolute functional topology (CWM, FDis, FRic) for isolated communities.
* `asc_null()`: Hierarchical triad of null models (Structural/Curveball, Quantitative/SAD, Identity/Trait Shuffle).
* `asc_transitions()`: Species-level Functional Leverage decomposition.
* `asc_entities()`: Species clustering into discrete functional entities.
* `assess_functional_space()`: Diagnostic tool for evaluating PCoA quality before analysis.

### Architecture
* Centralized internal engines (`.build_functional_space()`, `.calc_core_topology()`) shared across all modules.
* S3 class system with `print()`, `summary()`, and `plot()` methods for all output types.
* PCoA quality metric ($Q = \sum\lambda^+ / \sum|\lambda|$) with automatic warning when $Q < 0.80$.

### Null Models
* Model A (Structural): Curveball algorithm preserving richness and prevalence.
* Model B (Quantitative): SAD reshuffle with fixed incidence. FRic reported as NA (hull invariant).
* Model C (Identity): Trait shuffle breaking species-trait association.

### Validation
* 28 `testthat` unit tests covering class structure, edge cases, orthogonal scenarios, null model isolation, and leverage additivity.
* 61 functional tests (standalone script) including cross-validation against manual computation.
