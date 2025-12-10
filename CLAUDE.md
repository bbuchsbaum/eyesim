# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

eyesim is an R package for analyzing eye-movement data, focusing on
computing fixation pattern similarity, density map generation, scanpath
analysis, and recognition memory studies. The package implements various
eye-tracking analysis methods including MultiMatch algorithm, CRQA, and
EMD-based similarity measures.

## Common Development Commands

### Building and Testing

``` r
# Install dependencies
devtools::install_deps()

# Run tests
devtools::test()
# Run a single test file
testthat::test_file("tests/testthat/test_similarity_emd.R")

# Check package (includes running tests)
devtools::check()

# Generate documentation from roxygen2 comments
devtools::document()

# Build package
devtools::build()

# Install locally for testing
devtools::install()
```

### Documentation

``` r
# Build pkgdown documentation site
pkgdown::build_site()

# Check test coverage
covr::package_coverage()
```

## Architecture

### Core Data Structures

- **eye_table**: Main S3 class for eye-tracking data (defined in
  R/eye_frame.R)
  - Contains fixations with x, y coordinates, onset, duration, subject
    info
  - Supports grouping by trial, subject, condition
- **eye_density**: S3 class for fixation density maps (R/density.R)
  - Grid-based representation of fixation patterns
  - Supports various density estimation methods
- **fixation_group**: Groups of fixations with metadata (R/fixations.R)

### Key Modules

- **Similarity Analysis** (R/similarity.R): Core similarity metrics
  including:
  - Template-based similarity
  - EMD (Earth Mover’s Distance) based similarity
  - Spatial and temporal similarity measures
- **MultiMatch** (R/multimatch.R): Implementation of MultiMatch scanpath
  comparison algorithm
  - Computes 5 similarity dimensions: vector, direction, length,
    position, duration
- **Density Estimation** (R/density.R): Various density map generation
  methods
  - density_by() for grouped density estimation
  - multiscale_density() for multi-resolution analysis
- **CRQA** (R/crqa.R): Cross-recurrence quantification analysis for
  scanpath dynamics

### Testing Approach

Tests use testthat framework with fixtures in tests/testthat/. Each
major function has corresponding test file. Tests focus on: - Correct
computation of similarity metrics - Proper handling of edge cases (empty
data, single fixations) - Consistency of density estimation methods -
Proper S3 method dispatch

### Package Dependencies

Key dependencies include: - **proxy**: Distance matrix computation -
**imager**: Image processing for density maps  
- **ggplot2**: Visualization - **dplyr/tidyr**: Data manipulation -
**Matrix**: Sparse matrix operations for efficiency
