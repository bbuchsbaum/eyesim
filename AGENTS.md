# Repository Guidelines

## Project Structure & Module Organization

- `R/`: core package functions grouped by topic (`similarity.R`,
  `multimatch.R`, `density.R`, `eye_frame.R`, etc.); `man/` is generated
  from roxygen, so never edit it directly.
- Tests live in `tests/testthat/` with fixtures in `test_data/`; mirror
  function names in `test_*.R` files.
- Data assets: packaged objects in `data/`, creation scripts in
  `data-raw/`; vignettes and styling in `vignettes/` (albersdown theme),
  site config in `_pkgdown.yml`, rendered site in `docs/`.
- Edit prose in `README.Rmd` (not `README.md`); pkgdown extras and
  assets live in `pkgdown/`.

## Build, Test, and Development Commands

- Install deps: `R -q -e "devtools::install_deps()"`.
- Run suite: `R -q -e "devtools::test()"`; target a file with
  `R -q -e "testthat::test_file('tests/testthat/test_similarity_emd.R')"`.
- Full check (tests, lint-style checks): `R -q -e "devtools::check()"`.
- Refresh docs/NAMESPACE after code changes:
  `R -q -e "devtools::document()"`.
- Build docs site locally: `R -q -e "pkgdown::build_site()"`.
- Optional: enable MultiMatch parity tests by installing the Python
  module once via
  `R -q -e "reticulate::py_install('multimatch_gaze', pip=TRUE)"`.

## Coding Style & Naming Conventions

- Use `<-`, 2-space indents, and `snake_case` for functions/objects; S3
  methods follow `method.class`.
- Favor tidyverse verbs and pipelines (`dplyr`, `purrr`, `%>%`) as in
  `eye_frame.R`; keep argument order consistent across functions (data,
  grouping, options).
- Document exports with roxygen2 blocks; annotate internal helpers with
  `@noRd`. Update examples to run quickly and deterministically.
- Keep file and object names descriptive but concise (e.g.,
  `sample_density_time`, `multi_match`); avoid long inline
  functions—extract helpers in the same module.

## Testing Guidelines

- Place new coverage in `tests/testthat/test_<feature>.R`; group
  expectations by scenario and seed randomness with
  [`set.seed()`](https://rdrr.io/r/base/Random.html).
- Prefer small, explicit fixtures (tibbles or `test_data/` snippets)
  over serialized binaries; assert both happy paths and edge cases
  (empty inputs, size mismatches).
- For Python-backed comparisons, follow the existing `skip_on_cran()`
  and
  [`reticulate::py_module_available()`](https://rstudio.github.io/reticulate/reference/py_module_available.html)
  guards so CI remains stable when dependencies are absent.
- Check coverage locally with `R -q -e "covr::package_coverage()"` when
  touching core similarity/density paths.

## Commit & Pull Request Guidelines

- Commits should be short, imperative, and scope-focused (e.g., “Tighten
  template_similarity bounds”, “Add density parity tests”); batch
  unrelated changes separately.
- PRs need a brief summary, command output for `devtools::test()` or
  `devtools::check()`, and linked issues when relevant. Note any new
  data assets or API changes.
- Include screenshots or pkgdown previews when plot outputs or vignettes
  change. Keep branches rebased to minimize conflicts.
