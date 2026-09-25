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

## Amendment 1 (protocol `density-delta/1.1.0`, 2026-09-24)

**Timing and reason, stated plainly.** An external plan review requested
this amendment. It arrived **after** the 1.0.0 court had been frozen
(commit `83442b0`), and after its null (1) and real-data analyses had been
run and their outcomes viewed. The 1.0.0 results are reported unchanged under
Results. The changes below respond to the review's points, not to the 1.0.0
outcome, but 1.1.0 is not blind to that outcome. Treat 1.1.0 as a
reviewer-mandated re-specification, not as an independent preregistration.
The 1.0.0 synthetic sweep was stopped before completion and superseded; only
the 1.1.0 sweep is reported. Everything above the `FROZEN-PROTOCOL-END`
marker is byte-identical to 1.0.0. Where the amendment conflicts with that
text, the amendment governs 1.1.0.

### A1. Observation model: the duration-weighted spatial score

- **Unit.** One observation is one retrieval fixation's screen position x_j.
  Fixation count n and total dwell D are conditioned on, not modelled.
  Durations enter only as within-trial weights. Fixation order and saccades
  are ignored.
- **Component choice.** The mixture component is latent **per fixation**:
  each fixation independently comes from ε·U or from one of the free
  components with probability w_k. It is not chosen per trial.
- **Dependence.** Fixations are treated as conditionally independent given
  the templates and weights. This is a working-independence device for
  defining the score, not a belief. Dependence within trials, participants
  and items is carried by the crossed participant × item bootstrap, which
  resamples whole participants and items. The score is therefore a
  **duration-weighted spatial score** (a weighted composite log density). It
  is **not** a joint scanpath likelihood and is not described as one.
- **Trial aggregation (primary).** ℓ_M(Y_si) = Σ_j (d_j/D)·log p_M(x_j), the
  duration-weighted mean. Justification as in section 3: equal trial weight
  in the crossed mean and invariance to how a dwell is segmented.
- **Trial aggregation (sensitivity, a different estimand).**
  ℓ^sum_M(Y_si) = Σ_j log p_M(x_j), the unweighted summed fixation log
  density. This weights trials by fixation count and ignores duration. It is
  run as `sensitivity_fixation_sum`, and fitting also uses the summed
  objective, so fitting and scoring always share one aggregation.
- **Nesting and equal nuisance.** M0 is M1 with w_own = 0. Both models use
  the same background template and the same group template. They also share
  the same (h_g, h_b), selected once by the M0 search over the same 49-pair
  grid on the same training rows, and the same EM (ε = 0.01, tolerance 1e-10,
  at most 5000 iterations). M1's only additional fitting is the own weight and
  the own bandwidth, a 7-value search. Held-out scoring penalises any
  overfitting that this adds.

### A2. Background under participant holdout: declared support split

Evaluation participants have no training rows, so B cannot come from the
training set. 1.1.0 adopts option **(b), a declared support/evaluation split
within each participant**. B_si is built from participant s's retrieval
trials (all probe types, same window) on items in the **other item fold** from
item i. Each trial is unit-normalised and trials are averaged. The 1.0.0
rule, "all other items", is withdrawn.

Consequences:

- No support trial of B_si is evaluated in the outer fold where (s, i) is
  evaluated.
- No support trial involves the target item or any candidate item. Candidate
  sets, and hence wrong-item templates, are drawn from the target's own item
  fold.
- The same rule applies to every row:
  - training rows use their other-fold trials, which are never evaluation
    rows of that fold;
  - evaluation rows;
  - nulls (2) and (3) and the attribution analysis;
  - the null (1) generator and model;
  - every synthetic replicate and test.
- Perturbing any evaluation row of a fold changes neither that fold's fit nor
  any other trial's background in that fold (unit-tested).

### A3. Null decision rules

The **decision rule** under test is the primary rule itself. The statistic is
the mean trial Δ. Its 95% crossed-bootstrap percentile interval (2000 draws
for real data, 499 in simulation) is computed, and H0 is rejected if the lower
bound is > 0. This is one-sided at α = 0.025.

- **Null (1).** R = **200** replicates. Each replicate reruns the whole
  procedure on freshly simulated retrieval gaze: the full (h_g, h_b) and h_o
  searches, EM weight fitting for every fold, scoring and the bootstrap
  decision. The false-positive rate (FPR) is the proportion of replicates in
  which the rule rejects.
  - *Tolerance*: FPR ≤ α + 2·MCSE(α) = 0.025 + 2·√(0.025·0.975/200) ≈ 0.047.
    Monte-Carlo uncertainty is reported as a Clopper–Pearson 95% interval.
    The null fails if the observed FPR exceeds the tolerance.
  - The 1.0.0 criterion |mean Δ| ≤ 0.005 is **withdrawn as a gate**. Under a
    null generated from the group + background model, E[Δ] equals minus a KL
    divergence (adding an unneeded component can only lose held-out
    likelihood in expectation), not zero. Mean Δ is reported descriptively.
- **Simulation sweep.** At strength 0 the FPR of each rule (raw Δ, pseudo Δ,
  own − pseudo) comes from R = 100 full-procedure replicates per cell, with
  Clopper–Pearson intervals (MCSE at α = 0.025 is about 0.016). Power cells
  use the same replicates and rule.
- **Nulls (2) and (3)** are single real datasets, so an FPR cannot be
  estimated for them. They are judged on the decision rule. A null *fails* if:
  - the primary rule rejects on the null Δ, or
  - its 95% upper bound exceeds δ_tol = 0.01 nats, or
  - the corresponding primary − null contrast (paired for null 2, unpaired
    for null 3) fails to reject.

  A null Δ below zero is expected (see the KL remark above) and is not a
  failure.
- **Attribution gate**: unchanged from 1.0.0 (paired own − pseudo, rule as
  above).

### A4. Verdict rule

The verdict rule is unchanged from 1.0.0, with the gates as redefined in A3.

<!-- AMENDMENT-1-END -->

## Results

All real-data statistics below are aggregates. Per-trial scores, fold fits
and logs are in the Git-ignored results directory: the top level for 1.0.0
and `v1.1.0/` for 1.1.0. Δ is in nats per unit of dwell. The interval is the
95% crossed participant × item bootstrap interval.

### Protocol 1.1.0 (governing)

**Null (1)**, with 200 full-procedure replicates:

- FPR = 0/200 (Clopper–Pearson 95% interval 0–0.018), against a tolerance of
  0.047. **Pass.**
- Mean Δ = −0.00028 (SD 0.00045). This is negative, as the KL argument
  predicts.
- Mean fitted own weight = 0.008.

**Real court.** All four folds converged. Fitted M0 weights are
background 0.88–0.95, group 0.05–0.12 and uniform ≈ 0. The M1 own weight is
0.04–0.07.

| Analysis | Trials | Δ | 95% interval | Rejects |
|---|---|---|---|---|
| **Primary, old + lure, 0–3000 ms** | 1295 | **0.0026** | **−0.0029, 0.0080** | **no** |
| Old only | 640 | 0.0054 | −0.0012, 0.0132 | no |
| Lure only | 655 | −0.0001 | −0.0074, 0.0076 | no |
| Null 2, wrong-item own | 1295 | −0.0021 | −0.0036, −0.0008 | no |
| Primary − null 2 (paired) | 1295 | 0.0047 | −0.0005, 0.0101 | no |
| Null 3, pseudo-own on newtest | 959 | −0.0002 | −0.0017, 0.0012 | no |
| Primary − null 3 (unpaired) | 2254 | 0.0028 | −0.0028, 0.0087 | no |
| Attribution: pseudo-own on old + lure | 1295 | 0.0005 | −0.0013, 0.0025 | no |
| Attribution: own − pseudo (paired) | 1295 | 0.0021 | −0.0035, 0.0080 | no |
| Sensitivity: delay only, 500–3000 ms | 1295 | 0.0023 | −0.0040, 0.0090 | no |
| Sensitivity: 4-episode donor support | 1295 | 0.0039 | −0.0027, 0.0107 | no |
| Sensitivity: own bandwidth × 0.5 | 1295 | 0.0011 | −0.0021, 0.0047 | no |
| Sensitivity: own bandwidth × 2 | 1295 | 0.0026 | −0.0018, 0.0083 | no |
| Sensitivity: fixation-sum aggregation | 1295 | 0.0256 | −0.0070, 0.0631 | no |

Gates:

- Null (1): pass.
- Nulls (2) and (3): within tolerance (neither rejects, and both upper bounds
  are ≤ 0.01).
- Primary rejects: fail.
- Both primary − null contrasts: fail.
- Attribution gate: fail.

**Verdict: `not_detected_at_protocol_sensitivity`.** Every estimate is
small and positive, and every interval includes zero. The sensitivity
analyses agree in sign and overlap. The fixation-sum score is on a per-trial
rather than a per-dwell scale, so it is roughly 5–10 times larger by
construction.

### Synthetic sweep (1.1.0 generator and rules)

The sweep used 100 replicates per cell, 499 bootstrap draws and 56 min of
wall time. It reports the rejection rate of the primary one-sided rule. FPR
is at strength 0, and power is above that. Full tables are in
`gaze-weave-density-delta-simulation/`.

| Variant | Raw Δ FPR | own − pseudo FPR | own − pseudo power, s = 0.05 / 0.10 / 0.20 |
|---|---|---|---|
| Base (36 × 36, low overlap) | 0.90 | 0.00 | 1.00 / 1.00 / 1.00 |
| High overlap | 0.65 | 0.00 | 0.11 / 0.63 / 1.00 |
| Heterogeneous (τ = 1) | 0.97 | 0.00 | 0.98 / 1.00 / 1.00 |
| 24 × 36 | 0.47 | 0.00 | 0.98 / 1.00 / 1.00 |
| 12 × 24 | 0.05 | 0.01 | 0.50 / 0.97 / 1.00 |

Raw Δ is **not** a valid test of participant-specific correspondence under
this generator. With no own contribution at retrieval, raw Δ is about
+0.008 nats and rejects in up to 97% of replicates, because an own study
template adds image information that a 9-donor group template lacks. The
pseudo-own attribution contrast holds its FPR at 0–0.01.

Mean own − pseudo Δ scales with strength:

- low overlap: about 0.03 nats at s = 0.05 and 0.07 at s = 0.10;
- high overlap: about 0.005 at s = 0.05 and 0.011 at s = 0.10.

These are detection thresholds under this generator. With low overlap the
design detects an own share of 5% of retrieval fixations. With high overlap
(the own pattern is a reweighting of the image's own hotspots) it needs a
share of about 10–20%. They are not ceilings on the real effect.

**Reading the real result against the sweep.** The observed paired
own − pseudo Δ is 0.0021 (upper bound 0.0080). That is below the mean
contrast at s = 0.05 under low overlap, and comparable to high-overlap
shares of about 0.05–0.10. So a distinct, low-overlap own pattern carrying
≥ 5% of retrieval dwell is unlikely under this generator. An own pattern that
reweights the same salient regions others fixate, at a share of about 10%,
is not excluded. None of this is evidence of absent reactivation (section 10).
Participant × image preference remains a named alternative for any future
positive result.

### Protocol 1.0.0 (superseded; reported for completeness)

- 1.0.0 was frozen at `83442b0` and run in full.
- Null (1), 100 replicates: FPR = 0 and mean Δ = −0.0003. It passed its
  1.0.0 criteria.
- Primary Δ = 0.0014 (interval −0.0036 to 0.0063), not rejected.
- Null 2: Δ = −0.0012 (−0.0022, −0.0005).
- Null 3: Δ = −0.0005 (−0.0019, 0.0007).
- Attribution own − pseudo: 0.0014 (−0.0038, 0.0066).
- All sensitivities included zero.
- Verdict `not_detected_at_protocol_sensitivity`.
- 1.0.0 differed from 1.1.0 mainly in its background (all other items rather
  than other-fold items). The 1.0.0 synthetic sweep was stopped unfinished
  when Amendment 1 arrived.
