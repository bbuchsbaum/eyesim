# Density Δ court: own-encoding gain in held-out retrieval gaze

Protocol `density-delta/1.0.0`. Written 2026-09-24, before any real-data
outcome was computed. The section above the `FROZEN-PROTOCOL-END` marker, the
script `gaze-weave-density-delta.R`, and the configuration returned by
`density_delta_config()` are hashed (SHA-256) by `density_delta_freeze()` into
`gaze-weave-density-delta-freeze/freeze-record.csv`. The real-data runners stop
unless the hashes still match. Results are appended below the marker.

## 1. Question and claim

The court asks whether, during the recognition retrieval period, a
participant's gaze on an image is better predicted when the participant's
**own** encoding gaze on that image is added to a model that already knows
(a) where **other participants** looked when encoding the same image and
(b) where **this participant** tends to look during retrieval in general.

A positive, attributable result supports participant-specific
encoding–retrieval correspondence across the retrieval period. It does not
establish memory-driven reinstatement (see section 9).

Precedent: Wynn, Ryan and Buchsbaum (2020, PNAS; PMC7084073) compared
retrieval gaze with the participant's own encoding gaze and with other
participants' encoding gaze of the same image, treating the other-participant
comparison as the control for image-driven similarity. This court applies the
same logic as a held-out predictive density comparison with fitted,
population-level mixture weights.

## 2. Data, cohort and folds (reused, not re-implemented)

- Inputs: `test_data/wynn_probe_delay/` (`study_fix_input_new.csv`,
  `testdelay_fix_input_matched.csv`, verified against `manifest.csv` SHA-256).
  The path may be supplied through `EYESIM_WYNN_DATA_DIR`. Inputs are read-only
  and Git-ignored; see that directory's README for the handling contract.
- Importer: `probe_delay_*` helpers (image rectangle x = 112–912,
  y = 84–684 translated to an 800 × 600 analysis frame; fixations off the image
  rectangle are dropped; fixations are clipped to the analysis window and kept
  only if the clipped duration exceeds 80 ms). Study window 0–2500 ms.
- Cohort: `full_recognition_select_cohort()` and
  `full_recognition_validate_design()` — 36 participants, 1295 old/lure
  participant–item pairs, 120 items, each pair with four complete study
  presentations (≥ 3 fixations each) and ≥ 3 retrieval fixations in 0–3000 ms.
- Folds: `full_recognition_fold_plan()` (seed 20260825): two participant folds
  × two item folds = four outer folds. For outer fold (P_f, I_f), evaluation
  trials are P_f × I_f and training trials are P_¬f × I_¬f, so no training
  trial shares a participant or an item with an evaluated trial. Every cohort
  trial is evaluated exactly once.
- Donor pool: `own_group_complete_study()` (all 46 participants' complete
  four-presentation study records). Donor matching, rotation and the
  wrong-item candidate plan come from `own_group_select_donors()` and
  `full_recognition_candidate_plan()`.

## 3. Estimand

For retrieval trial (participant s, item i) with fixations j = 1…n at screen
positions x_j with durations d_j (clipped to the window), total dwell
D = Σ_j d_j, define the duration-weighted mean log density

  ℓ_M(Y_si) = Σ_j (d_j / D) · log p_M(x_j),

that is, the time-average of log p_M over the dwelled retrieval time, in nats
per unit dwell (densities are per px² of the 800 × 600 screen). The primary
estimand is

  **Δ_si = ℓ_{group+bg+own}(Y_si) − ℓ_{group+bg}(Y_si)**, window 0–3000 ms,

and the primary parameter is the population mean of Δ_si over old and lure
trials of the cohort (each trial weighted equally).

Why per unit duration rather than a sum over fixations: (i) every trial then
contributes on the same scale regardless of how many fixations it contains,
so the crossed mean is not dominated by high-fixation-count trials or
participants; (ii) duration weighting matches the duration-weighted templates,
so the model is a density of where gaze *dwells*; (iii) the value is
invariant to how a long dwell is segmented into fixations. Δ is in nats and
equals the mean per-dwell log likelihood ratio; multiplying by the mean
fixation count gives an approximate per-trial log Bayes factor.

## 4. Components and templates

All kernels are isotropic Gaussians truncated to the screen rectangle and
divided by their on-screen mass
m(μ, h) = [Φ((800 − μx)/h) − Φ(−μx/h)]·[Φ((600 − μy)/h) − Φ(−μy/h)],
so every kernel, and hence every template density, integrates to one over the
screen. Uniform density U = 1/(800·600).

- **Own template O_si**: the participant's four study presentations of item i
  (their study image version). Each presentation's fixation weights are
  d/Σd (unit mass per presentation), and the four are averaged with weight
  1/4 each.
- **Group template G_si**: every other participant with four complete study
  presentations of the same image version, excluding the retrieval
  participant. Each donor's template is built as for O; donors are averaged
  with equal weight. For lure trials the version is the participant's study
  version (as in the own-group court). For newtest trials it is the version
  shown at retrieval.
- **Background template B_si**: the participant's own retrieval fixations,
  same window, from all their retrieval trials (old, lure, newtest) except
  item i. Each trial is normalised to unit mass and trials are averaged.
  B contains no fixation from the target trial.

## 5. Models and fitting (training folds only)

  p_M(x) = ε·U + (1 − ε)·Σ_k w_k f_k(x),   ε = 0.01 fixed,

with free weights w on the simplex. The base model M0 has components
{G(h_g), B(h_b), U}; the full model M1 adds O(h_o). The fixed ε·U floor keeps
every log density finite and is identical in both models.

For each outer fold, using only the training trials P_¬f × I_¬f:

1. For every (h_g, h_b) in the grid 20·√2^k px, k = 0…6 (20–160 px), fit
   the M0 weights by EM to maximise Σ_trials ℓ_M0 (relative tolerance 1e-10,
   at most 5000 iterations). Keep the maximising pair (the first in grid order
   on ties).
2. With (h_g, h_b) fixed, for every h_o in the same seven values, fit the M1
   weights; keep the maximising h_o.
3. Score every evaluation trial with both fitted models; Δ_si = ℓ_M1 − ℓ_M0.

Weights and bandwidths are population-level (one set per fold). There are no
per-trial or per-participant weights. Because M0 is nested in M1 only on the
training data, held-out Δ can be negative.

## 6. Inference

Crossed participant × item (pigeonhole) bootstrap,
`transport_exploratory_bootstrap_plan()`: 2000 draws, seed 20260924. Each
draw resamples participants and items independently, and each trial receives
weight = participant multiplicity × item multiplicity. Intervals are 95%
percentile intervals (type 8).

**Primary test**: H0: mean Δ ≤ 0 versus H1: mean Δ > 0 over the 1295 old/lure
trials, one-sided α = 0.025. Reject H0 if the 95% interval's lower bound is
> 0. Old and lure are reported separately as descriptive secondaries.

## 7. Nulls, attribution gate and tolerances

Each null is judged against a declared tolerance, not an expectation of exactly
zero. δ_tol = 0.01 nats.

1. **Null (1): simulated group + background gaze, no own contribution.**
   For each fold, fit M0 on the real training trials. For every cohort trial,
   keep the real fixation count and durations and redraw positions from that
   trial's fitted M0 mixture (ε·U, G_si at h_g, B_si at h_b, U; the fitted
   weights of the fold in which the trial is evaluated). Templates are
   unchanged. Rerun the full court (fold fitting, scoring, crossed test with
   499 bootstrap draws). 100 replicates.
   *Pass*: rejection rate at the primary test ≤ 0.05 and |mean Δ| ≤ 0.005.
   This null checks fold-fitting optimism, leakage and bootstrap calibration.
   It is lenient about one mechanism: the generator equals the fitted M0, so an
   own template cannot add image information that the group template lacks.
   The attribution gate and null (3) address that mechanism.
2. **Null (2): matched wrong-item own templates.** Replace O_si by the
   participant's own four-presentation template for a different item: the
   first non-target candidate in the frozen five-candidate plan
   (`full_recognition_candidate_plan()`, same item fold). G and B are
   unchanged. Refit and score on the same old/lure trials.
   *Pass*: 95% upper bound of mean Δ_wrong ≤ δ_tol **and** the paired contrast
   Δ_own − Δ_wrong has lower bound > 0.
3. **Null (3): matched pseudo-own templates on newtest trials.** Newtest
   trials of cohort participants on cohort items with ≥ 3 retrieval fixations.
   The participant never studied the image. The pseudo-own template is the
   four study presentations of one other participant: the first donor chosen
   by `own_group_select_donors()` for that participant, item and version. That
   donor is removed from G. Folds are the same participant × item folds.
   *Pass*: 95% upper bound of mean Δ_pseudo ≤ δ_tol **and** the unpaired
   crossed contrast Δ_own(old/lure) − Δ_pseudo(newtest) has lower bound > 0.
4. **Attribution gate (pseudo-own on the same old/lure trials).** Same
   construction as null (3), but on the 1295 old/lure trials: O is replaced by
   one other participant's four presentations of the participant's study
   version, and that donor is removed from G. *Pass*: paired
   Δ_own − Δ_pseudo lower bound > 0.

   This gate was added before the freeze after a synthetic smoke run. The
   generator had no own contribution at retrieval, but own study gaze was
   half image-driven and the group template was built from 16 donors
   (36 × 36 design, one replicate). There,
   raw Δ was positive (0.005 nats, rejecting): an own template adds a further
   sample of the image-driven density to a finite group template. A
   pseudo-own template from another participant carries the same extra image
   information but no participant specificity, so the paired difference
   isolates the participant-specific part.

### Verdict rule

- `supports_participant_specific_correspondence`: the primary test rejects and
  nulls (1)–(3) and the attribution gate all pass.
- `primary_positive_but_not_attributable`: the primary test rejects but at
  least one null or the attribution gate fails.
- `not_detected_at_protocol_sensitivity`: the primary test does not reject.

## 8. Pre-specified sensitivity (reported, not gates)

- **Delay-only window, 500–3000 ms**: Y and B rebuilt from delay-window
  fixations. Same cohort trials, less those with no delay fixation.
- **Donor support**: G replaced by the four-episode matched-donor template of
  the own-group court (`own_group_select_donors()`: presentations 1–4 from up
  to four distinct donors), equal in episode count to O.
- **Own bandwidth**: h_o fixed at 0.5× and 2× the training-selected value (two
  steps of √2 on the extended grid 10–320 px), weights refitted.

Consistency is judged by the sign and overlap of intervals, not by a second
significance test.

## 9. Simulation sweep (synthetic data only)

Generator (`density_delta_simulate()`): an 800 × 600 screen. Each item has
four Gaussian hotspots (sd 35 px). Each participant has a background centre
bias (sd 130 px). Each participant × item has an idiosyncratic pattern:
*low overlap* means two new spots, and *high overlap* means a
Dirichlet(0.3) reweighting of the item's own hotspots. Each study presentation
draws about 5 fixations from 45% image, 10% background and 50% idiosyncratic
gaze. Retrieval draws 3 + Poisson(2.5) fixations with own share
s_si = min(strength · m_s, 0.95), where m_s = exp(τz − τ²/2) is
participant heterogeneity. The remainder is 60% image and 40% background.
Group templates use 9 donors, the median real-data group support
(minimum 2; from a pre-freeze design audit that built templates only and
computed no outcome). Folds are 2 × 2 crossed.
Each replicate runs the full court with both own and pseudo-own templates.

Cells: strength ∈ {0, 0.05, 0.1, 0.2, 0.3} × five variants: base
(36 participants × 36 items, low overlap, τ = 0), high overlap,
heterogeneous (τ = 1), 24 × 36, and 12 × 24. There are 100 replicates per
cell with 499 bootstrap draws.

Three tests are reported per cell, each at the primary one-sided α:

- raw Δ;
- pseudo-own Δ;
- the paired own − pseudo contrast, which is the attribution test.

The rejection rate at strength 0 is the false-positive rate, and above 0 it is
power. These are **detection thresholds** under this generator, the smallest
own shares the design detects reliably. They are not ceilings on the true
effect, and they do not transfer to generators with other spatial structure.

## 10. Interpretation

- The claim licensed by a `supports_…` verdict is **participant-specific
  encoding–retrieval spatial correspondence** over the retrieval period. The
  participant's own encoding density predicts their retrieval gaze beyond
  other participants' encoding and beyond their own item-general retrieval
  habits.
- **Named alternative: participant × image preference.** A participant may
  look at the same idiosyncratically preferred regions of an image whenever
  they see it, with no memory retrieval involved. Retrieval here shows the
  studied or a similar image, so this court cannot separate reinstatement
  from a stable participant × image viewing preference. Newtest trials cannot
  resolve it either, because the participant never encoded the image.
- **A null is not absence of reactivation.** Failing to detect Δ > 0 means
  only that any own-specific component lies below the detection threshold of
  this design and model (section 9). It does not mean retrieval gaze lacks
  reinstatement. Reinstatement might be temporally local, sequence-based
  rather than density-based, or too small for 1295 trials.
- Old versus lure differences and behavioural (oldness) links are outside this
  protocol.

## 11. Running

```r
devtools::load_all()
source("inst/validation/gaze-weave-density-delta.R")
density_delta_freeze()           # writes the committed freeze record
run_density_delta_simulation()   # synthetic; writes gaze-weave-density-delta-simulation/
run_density_delta_null1()        # participant-linked; results dir (ignored)
run_density_delta_court()        # participant-linked; results dir (ignored)
```

Participant-linked outputs, including per-trial scores, fold fits and caches,
stay in `inst/validation/gaze-weave-density-delta-results/`, which is
Git-ignored and excluded from builds. Only aggregate statistics are copied
below. Any change to the frozen section, script or configuration after the
freeze is a deviation. It must be listed under Results with its reason, and
the run must be re-frozen under a new protocol version.

<!-- FROZEN-PROTOCOL-END -->

## Results

Pending.
