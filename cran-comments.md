# Submission preparation

Preparation for a first submission of eyesim 0.3.0. No submission has been
performed. The current source includes the revisitation changes in
`f1cde6b1f14944e09333482de9d87e9c5a489657` and the CRAN-preparation changes.

# Current macOS validation (2026-10-29)

`devtools::check()` built the source archive and ran
`R CMD check --no-manual --as-cran`: **0 errors, 0 warnings, 0 notes**.
Installation, compiled-code checks, examples including `--run-donttest`,
tests, and vignette building and rebuilding passed. The environment was
macOS Sonoma 14.3 arm64, R 4.5.1, with Apple clang 15 and Homebrew Fortran
linkage. PDF manual generation was not rerun in this check.

Packaged tests passed 3,643 assertions with 39 skips and nine test warnings.
Four warnings concerned R dependencies built under R 4.5.2. Five concerned
optional Python discovery: sandbox access to the uv cache failed, so the
Python MultiMatch comparisons were skipped. These test warnings are separate
from the zero R CMD check warnings. Source-only research fixtures were
excluded from the archive.

* Source archive size: 1,449,829 bytes.
* SHA-256: `feb535914188693e05d44bd992107fa899ae900e9af75f25d3b0a5565ff115e0`.
* All changed package source and test files match the checked archive.
  This submission-preparation file is excluded from the package build.

The Windows and PDF-manual results below concern the earlier October 25
archive, before the revisitation changes. They do not validate the current
archive on Windows or establish a current PDF-manual result.

# Historical test environment (2026-10-25)

* Win-builder: Windows Server 2022 x64, R-devel 2026-10-24 r90762 ucrt,
  GCC and GNU Fortran 14.3.0; check dated 2026-10-25. Results supplied by
  the maintainer.
* macOS Sonoma 14.3, arm64, R 4.5.1.
* Apple clang 15, the default C++17 flags, and LAPACK/BLAS linkage.
* A temporary user Makevars selects the installed Homebrew Fortran runtime:
  this R installation's default `/opt/gfortran` runtime is absent. The
  package's Makevars is unchanged.
* All declared optional R dependencies were installed. CRAN-mode tests ran
  with `NOT_CRAN=false`; the slow and source-only research checks remain
  optional.

# Historical R CMD check results (2026-10-25)

## Windows R-devel (Win-builder)

**0 errors, 0 warnings, 1 note**. Installation, compiled-code checks,
examples, tests, rebuilt vignettes, and PDF and HTML manuals all passed.
Tests took 234 seconds; vignette rebuilding took 76 seconds.

The single incoming-feasibility note reports a new submission and possible
misspellings in DESCRIPTION. Dewhurst, Feng, Guibas, Rubner, Saenko and
Tomasi are cited authors' surnames. "et al." is the citation abbreviation,
and "scanpath" is an established eye-movement term. These are intentional.

## Local macOS

`R CMD build` followed by `R CMD check --as-cran` on the resulting tarball,
including PDF manual generation: **0 errors, 0 warnings, 2 notes**.

1. CRAN incoming feasibility reports **New submission**.
2. HTML manual validation was skipped because the local HTML Tidy executable
   is too old. PDF manual generation and all other documentation checks passed.

The CRAN-mode tests passed 2,337 assertions with 59 skips. Four test warnings
were traced to future, purrr, dplyr and tibble binaries built under R 4.5.2
being used by local R 4.5.1. These are separate from the zero R CMD check
warnings. The same upstream tests had no warnings in current-R CI.

All five vignettes built and rebuilt successfully. Normal examples and
`--run-donttest` examples passed. All six URLs found by urlchecker passed.

# Release changes

* Participant-level Wynn datasets and their help pages are no longer shipped.
  Package builds also exclude private copies and research input directories;
  the tutorials use synthetic data.
* The declared minimum R version now matches the ggplot2 dependency: R 4.1.0.
* A long revisitation example was wrapped to fit the PDF manual.
* Empty fixation groups now work when duration is inferred. Regression tests
  cover this behavior and the absence of the removed datasets.

# Historical artifact (2026-10-25)

* Source archive size: 1,443,893 bytes.
* SHA-256: `cbf65e3ea0fae2aa1babe345a5d7cf42ee22635ae0c7c7b3fd4b69f7fb22183f`.
* Based on `b6ccc1bc9967d32ab548db136a27d10dbaa0fd48` plus local changes.

# Before upload

Use an archive matching the intended release source. The current macOS
archive is identified in the October 29 validation section; the October 25
archive omits the subsequent revisitation changes. Refresh Windows R-devel
and PDF-manual validation for the intended release artifact before relying
on the historical results for submission.

The successful upstream Linux release/R-devel/oldrel, macOS release and
Windows release matrix below predates these local changes. It is supporting
evidence for the base revision rather than a check of this artifact:

<https://github.com/bbuchsbaum/eyesim/actions/runs/40187662008>

The upstream Linux extended suite passed 2,797 assertions with 38 skips;
Python parity, opt-in slow checks, and unavailable local research inputs
were among the skips. Its green status does not mean every research check ran.

CRAN submission guidance:
<https://cran.r-project.org/web/packages/submission_checklist.html>.
