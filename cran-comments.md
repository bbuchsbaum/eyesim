## Submission

This is the first CRAN submission of eyesim (version 0.3.0).

## R CMD check results

0 errors | 0 warnings | 2 notes

* New submission.

* Installed package size is about 11 Mb. Of this, `libs` (about 4.7 Mb) is
  compiled RcppArmadillo code, and `data` (3.8 Mb; 3.5 Mb on disk, already
  xz-compressed) holds `wynn_study` and `wynn_test`. These are fixations from
  a published recognition-memory experiment, and the examples and vignettes
  depend on them. The source tarball is under 5 Mb.

## Test environments

* Local: Ubuntu 24.04, R 4.3.3, `R CMD check --as-cran`.

## Notes for the reviewer

* `install_multimatch()` installs a Python package through 'reticulate', so
  its example is wrapped in `\dontrun{}`. The function runs only when the
  user calls it.
* Functions that take a `seed` restore the user's random-number state on exit.
