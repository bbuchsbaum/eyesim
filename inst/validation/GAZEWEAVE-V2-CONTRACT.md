# GazeWeave v2 scientific and evaluation contract

Status: accepted for implementation
Date: 2026-08-15
Mote gate: `bd-01M03RQPJ6MACB9XJXDR588H41`

## Decision

GazeWeave will have two engines that answer different questions through one
candidate-evidence interface.

1. **Transport v2** is a symmetric, coverage-conditioned comparison of two
   duration-weighted gaze paths. It estimates how much mass can be matched and
   the conditional spatial and local-order fidelity of that match.
2. **Replay** is directional. It models a recall trajectory as a mixture of
   locally ordered replay of an encoding path and candidate-independent
   background gaze.

Neither engine is the default merely because it has been implemented. Replay
must earn default status for encoding-to-recall studies in the synthetic and
real-data courts. Transport remains available for comparisons without a
privileged direction.

The sole default inferential endpoint is held-out item information gain,

\[
G=\log_2\{p(k^*\mid Y,\mathcal C)/\pi_{k^*}\}.
\]

Coverage, selection, correspondence, spatial fidelity, chronology, warp,
candidate rank, calibration, and stability are diagnostics.

## Shared scientific objects

### Inputs

A fixation path contains finite `x`, `y`, `onset`, and non-negative `duration`.
Zero-duration events are discarded. Positive durations are normalized within
trial whenever they define mass. Coordinate units must be declared and agree
with the spatial bandwidth and screen specification.

Registration is learned without the scored participant/item unit and is
applied identically to every candidate. The default warp is isotropic
participant/session contraction about the known screen centre plus a
participant/session translation. Pair-specific translation and affine warps
are not defaults.

### Four distinct questions

Every result must keep these objects separate:

- **coverage**: how much recall gaze participates in the explanation;
- **selection**: which encoding and recall mass participates;
- **correspondence**: how selected mass or time points align;
- **evidence**: how much the candidate improves prediction of item identity.

An optimizer coupling is an alignment. It is not a posterior distribution.
Only quantities produced by a normalized probabilistic Replay model may be
called posterior probabilities, conditional on the fitted model.

### Shared engine result

A candidate engine returns a `gaze_engine_result` with:

- `engine`: stable engine identifier;
- `candidate_key`: candidate identity;
- `log_score`: a natural-log score for which larger is better;
- `diagnostics`: named, engine-specific audit quantities;
- `alignment`: an engine-specific alignment object;
- `convergence`: solver or forward-algorithm status;
- `provenance`: engine version, directionality, score semantics, duration
  semantics, and confirmation that candidate-invariant preprocessing was used.

Candidate scoring must never fit a warp, background model, temperature, or
other nuisance parameter from the pair being scored.

## Transport v2

Let `a` and `b` be unit duration masses and let

\[
\Pi=M\Gamma,\qquad |\Gamma|=1,
\]

where `M` is replay coverage. For each fixed coverage value, the selected
marginals must be dominated by the original measures:

\[
M\Gamma\mathbf 1\le a,\qquad
M\Gamma^\top\mathbf 1\le b.
\]

This prevents a fixation from contributing more mass than it contains. The
fixed-mass profile is evaluated over a preregistered grid that includes full
coverage. Coverage is selected or marginalized using a policy trained strictly
inside the outer fold. The implementation must not divide an unconstrained
objective by a fitted `M` near zero.

For normalized correspondence `Gamma`, define

\[
r_S=\sum_{ij}\Gamma_{ij}c(x_i,T(y_j))
\]

and

\[
r_T=\sum_{ii'jj'}\Gamma_{ij}\Gamma_{i'j'}
       (R^X_{ii'}-R^Y_{jj'})^2.
\]

The scientific fit at fixed coverage is

\[
Q_0(M,\Gamma)=M\{r_S+\lambda_T r_T\}.
\]

Correspondence regularization may be used to find a stable solution:

\[
Q_\epsilon=Q_0+\epsilon M
\mathrm{KL}(\Gamma\|\alpha\beta^\top),
\]

where `alpha` and `beta` are the selected unit marginals. This mutual
information term is computational smoothing: it penalizes concentrated
dependence and is not psychological uncertainty. Candidate evidence is based
on `Q_0` evaluated at the regularized solution, plus the cross-fitted coverage
policy; the regularizer is not presented as evidence. Stability as `epsilon`
decreases is required.

The default chronology is directed ordinal local order with a fixed
neighbourhood of the next two events. Physical-time chronology is an explicit
alternative, not silently combined with order chronology. Exact adjacent
same-location events are coalesced before ordinal relations are built; a
non-zero coalescing tolerance is a declared sensitivity parameter in visual
degrees. Returns after an intervening fixation are never coalesced.

Transport diagnostics are:

- `replay_coverage = M`;
- selected reference and source marginals;
- `spatial_residual` in native cost units and a distance summary in screen
  units;
- `local_order_error = r_T`, bounded to its declared relation scale;
- contraction and translation;
- objective and alignment variation across starts and continuation levels;
- a stationarity or conditional-gradient diagnostic when available.

## Directional Replay

Replay conditions on an encoding template `X` and models registered recall
gaze `Y`. Duration is represented on a fixed normalized gaze-time grid:

1. discard zero-duration events;
2. concatenate positive fixation intervals in duration mass time;
3. sample the piecewise-constant spatial trajectory at a fixed number `L` of
   bin midpoints.

Consequently, a long fixation occupies more observations, uniform temporal
dilation is irrelevant, and splitting or merging adjacent identical fixation
intervals leaves the sampled trajectory unchanged. `L` is fixed before model
training; convergence across a small declared resolution set is a required
sensitivity check.

For each grid point `t`, the state is

\[
z_t\in\{0,1,\ldots,m\},
\]

where `0` is background gaze and state `i` is encoding fixation `i`.

The transition family supports:

- staying at the current encoding state;
- moving forward by at most `k` encoding fixations;
- restarting at any encoding state under the normalized encoding-duration
  prior;
- entering, remaining in, and leaving the background state.

The default `k` is two. Transition probabilities are shared or hierarchically
pooled and estimated only from training trials. A restart permits global chunk
reordering but pays an explicit probability cost.

Replay spatial emissions are normalized bivariate Student densities with fixed
degrees of freedom and a cross-fitted scale. The background emission is a
low-capacity participant/group recall distribution learned only from training
data and shared by all candidates in an evaluation row. The same held-out warp
policy and warp distribution are used for every candidate.

Forward-backward inference is performed in log space. It returns a normalized
candidate log likelihood, state posterior conditional on the fitted model,
background/replay occupancy, restart expectations, encoding-state visitation,
and barycentric correspondence. Template length, spatial scale versus
background rate, restart versus local mismatch, and warp versus emission noise
are mandatory identifiability tests.

## Candidate probabilities and information

For candidate natural-log scores `s_k`, declared priors `pi_k`, and calibration
temperature `tau`,

\[
p(k\mid Y,\mathcal C)=
\frac{\pi_k\exp(s_k/\tau)}
     {\sum_l\pi_l\exp(s_l/\tau)}.
\]

Transport requires a single temperature fitted by inner-fold log loss. Replay
also fits one temperature from inner-held-out candidate scores within each
outer-training set. Its fixed weak Gaussian prior on log temperature is centred
at `tau = 1`, preventing separable small training sets from forcing a boundary
solution. If a contrast stratum cannot supply two candidate-bearing inner
folds, Replay records an explicit `tau = 1` fallback. Temperature, background,
warp, and engine parameters are frozen before outer-fold evaluation.

The primary trial score is

\[
\texttt{gaze_info_bits}=\log_2(p_{k^*}/\pi_{k^*}).
\]

The secondary odds score is

\[
\texttt{odds_bits}=\log_2
\left[
\frac{p_{k^*}/(1-p_{k^*})}
     {\pi_{k^*}/(1-\pi_{k^*})}
\right].
\]

`gaze_info_bits` is bounded above by `log2(1 / prior_true)` and unbounded below.
Its held-out mean is improvement in logarithmic score relative to the declared
prior. It is called mutual information only when the candidate prior and
calibrated posterior represent the target data-generating distribution.

Candidate tables retain total candidate count, prior, posterior, calibrated and
raw score, true-candidate indicator, rank with numerical tie tolerance, and
candidate-pool provenance. Comparisons across studies require the same target
candidate-sampling policy or an explicit adjustment.

## Cross-fitting and fair comparison

The outer held-out unit is the union of the declared participant and item unit.
No outer evaluation observation may contribute to:

- warp or background estimation;
- chronology, coverage, transition, spatial-noise, or temperature selection;
- learned MultiMatch or density-composite coefficients;
- candidate filtering based on outcomes.

The frozen-metric court compares untrained scores. The supervised court gives
Transport calibration, Replay calibration when justified, a nested ridge or
ranking MultiMatch composite, and a comparably supervised multiscale-density
model identical inner-fold access. Elastic looking-at-nothing matching is a
contraction-aware, order-free baseline where available.

## Invariance and failure contracts

Fast tests must cover:

- stable log-sum-exp candidate probabilities and declared non-uniform priors;
- correct information and odds identities;
- row-order and candidate-order invariance;
- no overlap of outer train and evaluation participant/item keys;
- uniform temporal dilation;
- splitting and merging adjacent identical fixation intervals;
- pixel/degree agreement after correct conversion;
- finite results for minimal paths and explicit errors for invalid paths.

Reference tests must compare optimized relation costs with direct four-index
oracles and forward likelihoods with exhaustive state enumeration on tiny
models. Stress tests cover fixation count, duration concentration, candidate
difficulty, near-zero/full coverage, solver starts, template length, and time
grid resolution.

## Frozen v1 baseline

The current files under `R/gaze_weave_*.R`, their tests, and
`inst/validation/results` define the v1 baseline. The existing report is
synthetic evidence only. It does not establish real-data validity, calibrated
probabilities, human interpretability, or universal superiority. V2 work must
retain a reproducible v1 path until an explicit compatibility decision is made
at GW-8.

## Non-goals

- No free pair-specific affine rescue.
- No claim that one number implies a parameter-free model.
- No interpretation of objective ledger components as independent evidence
  bits.
- No claim that a diffuse optimized transport coupling is a posterior.
- No default-engine decision based only on a generator tailored to that engine.
- No superiority claim based on comparing a trained GazeWeave model with
  untrained or multiplicity-penalized baselines.
