# Calibration controls for GazeWeave revision 2026.10

Controls how Transport and Replay turn held-out candidate scores into
candidate probabilities under revision \`"2026.10"\`. Pass the result as
\`calibration_control\` to \[gaze_transport_spec()\] or
\[gaze_replay_spec()\].

## Usage

``` r
gaze_calibration_control(
  method = c("evidence_scaled", "global"),
  gamma_bounds = c(0, 1),
  inverse_temperature_prior_sd = 3,
  stein_shrinkage = TRUE,
  typicality = NULL,
  typicality_sources = 24L,
  typicality_item_on = NULL
)
```

## Arguments

- method:

  \`"evidence_scaled"\` (default) fits \\T\\ and \\\gamma\\ as
  described. \`"global"\` reproduces the pre-A3 revision \`"2026.10"\`
  calibration: one temperature with a log-normal prior centred at one
  and, when the engine's \`reliability\` is \`"effective_fixations"\`,
  the \\\kappa\\ shrink. The pre-change Transport default used that
  shrink, so reproducing pre-A3 \`"2026.10"\` requires \`reliability =
  "effective_fixations"\` passed explicitly; with the (new) default
  \`reliability = "none"\` the global path fits a different model.
  Probabilities and bits then match exactly; the specification
  additionally carries this control. Rank and top-1 use the reference
  ranking described under Ranking, so they differ from the pre-change
  values only in rows where the fitted temperature or the \\\kappa\\
  mixture had reordered the candidates (multi-episode Transport, or a
  non-uniform prior).

- gamma_bounds:

  Bounds for \\\gamma\\; equal values fix it (for example \`c(0, 0)\`
  for one global temperature with the new prior).

- inverse_temperature_prior_sd:

  Standard deviation of the zero-centred Gaussian prior on the
  standardised inverse temperature \\\tilde\beta\\. \`Inf\` removes the
  penalty.

- stein_shrinkage:

  Logical. When \`TRUE\` (default), the fitted inverse temperature is
  multiplied by \\\max(0, 1 - p / LR)\\, where \\LR\\ is the inner
  likelihood-ratio statistic of the fit against the declared prior and
  \\p\\ its number of free parameters (see Details).

- typicality:

  \`"none"\`, \`"mean"\` (subtract the offset) or \`"standardized"\`
  (subtract the offset from the row-centred score and divide by the
  candidate's shrunk score standard deviation over the sources).
  \`NULL\` (the default) uses the engine default: \`"standardized"\` for
  Transport and \`"none"\` for Replay. On a centre-biased simulated null
  with six candidates per pool, \`"standardized"\` brought Transport's
  central-candidate argmax share per candidate from 0.40 to 0.16 (chance
  0.167, MC error 0.020) without lowering top-1 or AUC on signal data.
  For Replay it lowered the share only from 0.30 to about 0.25-0.27 (MC
  error 0.015), and it is therefore opt-in for Replay. The residual is
  not a heavy-tail effect. Under these nulls the HMM background state
  absorbs most recall fixations, so the candidates' likelihoods nearly
  coincide: in 57% of held-out rows every candidate lay within 0.01 nats
  of the others. The very large per-candidate SD ratio across rows comes
  from these near-degenerate pools, and the central bias from the argmax
  among near-ties: central candidates won 0.54 of the tie rows against a
  0.37 share of the candidates (few rows, so imprecise). An offset
  cannot remove it; a pre-declared tie tolerance for top-1 and AUC (see
  Ranking) is the appropriate treatment.

- typicality_sources:

  Maximum number of training recalls used per fold to estimate every
  candidate's offset.

- typicality_item_on:

  Columns defining an "item" for the other-item rule. \`NULL\` uses
  \`setdiff(match_on, contrast_on)\`, or \`match_on\` when that is
  empty.

## Value

A \`gaze_calibration_control\`.

## Details

\*\*Evidence-scaled temperature.\*\* Row \\i\\ is scored at temperature
\\T_i = T (\bar n / n_i)^\gamma\\, where \\n_i\\ is the row's evidence
count (Transport: duration-effective fixations of the recall; Replay:
the number of recall fixations the HMM observes) and \\\bar n\\ is the
geometric mean over the calibration rows. \\T\\ and \\\gamma\\ are
fitted jointly by minimising the candidate log loss of inner out-of-fold
rows; no held-out row enters the fit. \\\gamma = 0\\ is one global
temperature.

\*\*Ranking.\*\* Calibration changes probabilities and bits only. Under
revision \`"2026.10"\`, \`template_rank\` and \`top1_credit\` (and any
AUC computed from ranks) always use one fixed reference ranking: the log
posterior at inverse temperature one, \\\log \pi_k + r_k\\, where
\\r_k\\ is the engine's native candidate score after any typicality
offset (Transport: the log mean of the episode scores; Replay: the trial
log likelihood) and \\\pi_k\\ the declared prior (uniform unless
\`priorvar\` is given). This holds for both methods and every fitted
\\T\\, \\\gamma\\ or Stein factor, including a calibration that returns
the prior. The ranking cannot be read off the calibrated probabilities:
for multi-episode Transport the inverse temperature acts inside the log
mean over episodes, which ranks by the arithmetic episode mean as \\1/T
\to 0\\ and by the log mean at \\T = 1\\. Transport per-episode ranks
use the same rule. The reference score uses no label, so relabelling a
row permutes it but never changes it. Ties use a relative tolerance of
\\\sqrt{\epsilon}\\ on this log scale. Replay's trial log likelihoods
are often near-ties (see \`typicality\`), so a Replay top-1 or AUC
should be computed from the candidate tables' \`ranking_score\` with a
tie tolerance declared before the analysis (for example 0.01 nats),
splitting credit among tied candidates.

The fit penalises the standardised inverse temperature \\\tilde\beta =
\sigma / T\\ with a Gaussian prior centred at zero, where \\\sigma\\ is
the median within-row standard deviation of the calibration rows'
candidate scores (label-free). The prior therefore only ever pulls
toward the declared candidate prior (less confidence), never toward a
fixed temperature. The earlier prior on \\\log T\\, centred at \\T =
1\\, pulled toward overconfidence whenever the scores' natural scale
exceeded one nat (Replay's total trial log likelihood). The temperature
has no upper bound (\\T = \infty\\ returns the prior); the engine's
lower temperature bound still caps every row's confidence.

A calibration fitted on a few dozen inner rows is noisy, and under a
null every positive inverse temperature is overconfidence on held-out
rows. The fitted \\\tilde\beta\\ is therefore shrunk by the
positive-part Stein factor \\\max(0, 1 - p / LR)\\, where \\LR =
2N(L_0 - L)\\ is the inner likelihood-ratio statistic against the
declared prior (\\L_0\\: its log loss) and \\p\\ the number of free
parameters (two when \\\gamma\\ is fitted). Under a null \\LR\\ is about
\\\chi^2_p\\ and the factor is small or zero; with real evidence it is
close to one. A hard selection test (keep the fit only if it beats the
prior by AIC) was rejected: conditioning on passing inflates the
selected fit, which made weak signals overconfident.

The effective-fixation reliability shrink (\\\kappa\\) is not used by
\`"evidence_scaled"\`: it could only make sparse rows less confident,
never rich rows more confident, and with \\\gamma\\ it is a second,
weakly identified parameter on the same evidence axis.

\*\*Typicality offset.\*\* With \`typicality = TRUE\`, each candidate's
score has its typicality \\b_k\\ subtracted before the softmax. \\b_k\\
is the mean score of candidate \\k\\ against a fixed, seeded subsample
of at most \`typicality_sources\` recall paths from the fold's training
rows whose item (the \`typicality_item_on\` key) differs from candidate
\\k\\'s item. It is computed from the candidate's encodings and training
recalls only, so it never depends on the held-out recall or on which
candidate is labelled true, and relabelling a held-out row leaves every
candidate score unchanged. The subsample is drawn round-robin over
items. Scores are centred within each source before averaging; the
removed per-source constant is common to all candidates and cancels in
the softmax, but for Replay (total log likelihood) it would otherwise
dominate the sampling noise, because it grows with the recall's fixation
count. Each mean is shrunk toward the common mean by the empirical-Bayes
factor \\\tau^2 / (\tau^2 + s_k^2)\\, where \\s_k^2\\ is the sampling
variance of candidate \\k\\'s mean and \\\tau^2\\ the between-candidate
variance beyond sampling noise, so that offsets made mostly of noise do
not perturb the ranking. When a candidate has fewer than two other-item
sources, no offset is applied in that fold and the status is recorded.
Inner calibration rows use offsets estimated from their own inner
training rows. Offset-adjusted scores are ranking scores, not normalised
likelihoods. Two approximations: Transport offsets and scales are
estimated over each candidate's valid study episodes, while a scored row
uses the episodes valid for its whole pool; and Replay offsets always
use the training background protocol, also under \`background_support =
"held_out"\`. Subtracting a mean equalises the candidates' average
scores under generic gaze but not their score variances; the fit records
the per-candidate standard deviation over the typicality sources.
