# Reference scoring contracts (#53)

The subsystem separates map evaluation, reference construction, random
assignments, deterministic scoring, and fitted scaling/summaries. Existing APIs
retain their defaults and output shapes while using the same numerical helpers.

## Representations

* Density maps retain their existing x/y/z representation. Sampling preserves
  nearest-grid lookup, half-to-even ties, and edge clamping. Map normalization
  precedes an optional log(value + floor); floor is an explicit positive offset.
* A reference set stores unique crossed-factor keys, transformed source maps,
  grouped sums and counts. It does not duplicate four full maps per input row.
  A reference map accessor materializes one map for inspection. Log reference
  maps are score maps, not probability densities, and are never logged twice.
* Controls are a long table with case_id, control_id (donor input row, or NA for
  generated points), draw, stratum, is_case and x/y. Other columns belong to the
  case. These explicit assignments can be saved and rescored without RNG use.
* Mismatch assignments contain one replacement reference row per subject/trial
  and null draw. The same replacement applies to all fixations and control
  draws in that trial. Assignment construction and evaluation are separate APIs.
* Scores retain the strata columns and one numeric column per reference.
  Scaling is a separate fitted object with group keys, reference names, means,
  SDs and counts. Applying it never estimates new parameters.

## Scientific choices

Reference defaults normalize each input map to unit sum, then log with an
explicit offset, then aggregate. Density-space aggregation is also supported.
An explicit transform="log" after density-space aggregation supports log of
the arithmetic mean. Log-space aggregation rejects another log transform.
For each target subject/item: self_item uses its map; self_other_items averages
that subject's other valid maps; others_item averages other subjects' maps of
the target item; others_other_items first averages other items within each other
subject, then averages those subject means equally. NULL maps are missing, not
zeros. min_n counts independent averaging units (maps for self references,
subjects for other-subject references). Insufficient references produce NA
scores and a construction warning, preserving target rows. Duplicate keys and
incompatible grids are errors. Contributor counts are retained.

Observed controls are drawn within subject and optional strata, excluding the
case's entire trial and, when item metadata is available, all trials of the
case's item. The item argument defaults to the item column when present;
custom column names are explicit, and item=NULL preserves trial-only sampling.
Item-aware exclusion avoids attenuating matched scores in repeated-item designs.
Replacement is automatic only when the pool is too small;
empty pools are errors. Uniform controls use explicit bounds; map controls use
nonnegative grid masses with replacement and return grid coordinates. Sampling
without replacement distorts the intended map distribution and is rejected.
Atomic donor metadata is retained once as control_sources for replay and case
subsetting; assignment rows themselves carry case metadata. No control-only jitter.
New seeded sampling APIs restore the caller's RNG state. Legacy wrappers keep
their existing session-RNG behavior and candidate order.

Within-subject scaling is per reference across all case/control rows, separately
for each control draw. Thus adding control draws cannot change existing draws'
scaling. Missing scores remain missing; constant finite columns map to zero.
Contrasts require one case and at least one control per stratum; missing scores
propagate. Trial contrasts equally weight fixation contrasts. Draws remain
separate unless explicitly averaged by the caller.

Mismatch self_item replacement changes only that score; all replacement changes
the entire reference tuple. By default, controls from every trial of either
the case's item or the replacement item are excluded without redrawing,
including when replaying older trial-only assignments. The paired_observed and
contrasts tables use identical retained controls and the same scaling, fitted
on retained mismatched scores per draw (or the original observed fit when
restandardize=FALSE). Empty retained pools are errors; pools containing only
the two compared items cannot support the paired contrast, even with repeated
trials or more sampled controls. The fixed policy retains
all assignments for the original item-discrimination baseline; its null can be
negative under a genuine effect. The original observed contrasts remain
available for context and may be attenuated if supplied assignments include
own-item donors. Generated controls have no donor item. Neither policy alone
establishes exchangeability or inferential calibration. This API returns scores
and assignments, not automatic significance claims.

Overlap windows select positive interval overlap and default to clipped-duration
weighting. Explicit weights represent full-fixation mass and are multiplied by
overlap/duration. Explicit duration_weighted=FALSE retains unit weights with a
warning. Onset filtering and its unweighted eye_density default remain unchanged.
Times and durations must use
the same units; window intervals are half-open.

## Delivery and evidence

Extract shared grid, map preparation, matching and candidate sampling first;
migrate existing callers. Then add references and assignment/scoring APIs.
Validate against independent slow references on balanced and incomplete designs,
legacy behavior and RNG snapshots, deterministic replay, known-answer synthetic
scenarios, full package tests/checks and a runnable synthetic vignette. Bandwidth
selection, outcome models and publication are outside this change.

## Measured scale after review corrections

The reproducible benchmark in inst/benchmarks/reference-scoring.R now defaults
to 103 subjects x 120 items x 8 test fixations: 98,880 cases, 20 controls per
case, 2,076,480 assignment rows, 20x20 maps and 100 mismatch draws summarized by
trial. On local macOS arm64 / R 4.5.1 (2026-11-12), references took 3.7 s,
controls 23.8 s, scoring 5.5 s, and corrected paired mismatch draws 222.9 s
(255.9 s total; mismatch user CPU 148.2 s and system CPU 24.7 s). Package checks
ran concurrently; this is a local measurement, not an isolated hardware
benchmark. References occupied
41.2 MB and assignments 282.5 MB (object sizes, not peak process memory).
The earlier 12,360-case measurement used one fixation per trial and did not
establish realistic-scale performance.

A later revision (2026-11-13) reduces each null draw instead of rescoring it
from scratch. Matched case-minus-control differences are reduced once; each
draw removes only its excluded donors and rescores only the replaced
reference. Restandardization uses group moments shifted by the observed
centres, from which unchanged columns subtract the dropped donors. Outputs
matched the previous implementation to 1e-12 on 32 configurations and on the
full benchmark. With the same benchmark, again on a loaded machine (load
average about 20), mismatch draws took 64.5 s (52.4 s user CPU), controls
11.7 s, scoring 3.3 s and references 1.0 s. The CPU time meets the one-minute
target; wall time on an idle machine was not measured.
