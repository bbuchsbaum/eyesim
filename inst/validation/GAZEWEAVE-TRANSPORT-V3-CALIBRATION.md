# GazeWeave Transport v3 nested-calibration audit

Protocol version: 3.0.2
Gate: GW-18.7
Status: implemented and accepted on the response-blind fixture court

## Estimand and nesting

`gaze_transport_v3_cv()` holds out participant-item cells in outer folds. Each
outer-training set is split again by item for nuisance fitting and calibration.
The fold receipts record original row ids, encoded participant-item keys, the
warp receipt, and zero train/evaluation overlap for outer and inner folds. A
held-out source path is never used to fit its warp, temperature, or reliability
parameter.

Coverage remains the prospectively frozen `Beta(2, 2)` integral; it is not
retuned after outcomes. Warp hyperparameters remain fixed in the specification,
while warp parameters are fitted only from the applicable training cells.
Every candidate in a row receives that same fitted warp and solver policy.

## Candidates, episodes, and priors

The default implementation scores every permitted reference candidate in the
row's declared contrast stratum. `priorvar` supplies the actual positive design
weights; omitting it explicitly declares a uniform design. Candidate
probabilities are normalized over that exhaustive pool, and the endpoint is
computed directly as `log2(p_true / prior_true)`.

Separate presentation likelihoods receive equal prior weight after applying
one common calibrated energy scale. The candidate score at temperature `tau`
is `logmeanexp(episode_score / tau)`, so calibration occurs at the frozen
episode-likelihood level rather than as an approximate post-mixture rescaling.
All candidates use the intersection of valid presentation ids.

## Reliability

The only reliability covariate is response-blind duration-effective fixation
count. Reliability is `n_eff / (n_eff + kappa)` and shrinks the calibrated
posterior toward the declared candidate prior. Joint temperature/reliability
optimization contains `kappa = 0` as the exact temperature-only boundary. The
joint result is rejected whenever its unpenalized inner out-of-fold log loss is
worse than that boundary by more than `1e-10`.

Raw and coalesced fixation count, effective count, total duration, duration
concentration, duration entropy, normalized entropy, spatial dispersion, and
generic gaze quality are exported as diagnostics. Generic quality is explicitly
labelled a separate response-blind channel and is not item-specific evidence.

## Frozen synthetic court

The deterministic court has two participants, four items per participant,
four separate study presentations per item, a nonuniform design prior
proportional to item index, and two outer plus two inner folds.

| Check | Result |
|---|---:|
| Outer train/evaluation participant-item overlap | 0 in every fold |
| Inner train/evaluation participant-item overlap | 0 in every fold |
| Candidate count | 4 exhaustive candidates per row |
| Common presentation count | 4 per row |
| Maximum reliability log-loss worsening | -3.59e-6 |
| Held-out mean log loss | 0.1341 |
| Held-out expected calibration error | 0.0623 |
| Prior recovery under null scores | exact to 1e-12 |
| Null candidate-pool sizes checked | 2, 3, and 5 |

Candidate order, design priors, equal episode weights, warp parameters, solver
backend, and convergence are inspected in the executable fixture tests. Pool
size is exported with held-out mean information and log loss so the full court
can assess candidate-pool sensitivity without relabelling it as reinstatement.

No retrieval response, correctness, confidence, saliency, old/lure label, or
private Wynn outcome is read by this gate.
