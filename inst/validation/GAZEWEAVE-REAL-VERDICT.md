# GazeWeave real-control verdict

Final seed: `20260817`
Protocol version: 1 with Replay calibration implementation addendum B
Decision: six of seven implementation gates passed; the repeated-viewing gate
failed; there is **no real-data engine default** and superiority is not
supported.

## Court and gate result

The final court evaluated 8 participants (4 young, 4 older) and 8 images in
four Cartesian participant-by-item folds. Every model was trained after
excluding both the evaluated participants and evaluated items. Each held-out
path was identified among four within-participant encoding candidates. All
methods used the same folds, candidates, signal trials, and negative controls.

Six gates passed: blank-screen signal, local-order ablation, empirical null
error, null-safe registration, engine convergence, and two-way cross-fitting.
The repeated-viewing gate required both GazeWeave engines to have positive
mean item information and top-one credit above 0.25. Transport passed those
criteria; Replay had negative mean information and therefore failed the gate.

The machine decision is:

```text
advance = FALSE
provisional_default = no_real_data_default
strongest_baseline = density_ridge_registered
superiority_supported = FALSE
```

## Held-out signal

Chance top-one identification was 0.25 and chance log loss was
`log(4) = 1.386`.

| Task | Method | Information bits | Log loss | Top-one credit |
|---|---|---:|---:|---:|
| Repeated viewing | Registered density ridge | 0.467 | 1.062 | 0.547 |
| Repeated viewing | Transport v2 | 0.271 | 1.199 | 0.531 |
| Repeated viewing | MultiMatch ridge | 0.100 | 1.317 | 0.469 |
| Repeated viewing | Elastic ridge | -0.024 | 1.403 | 0.391 |
| Repeated viewing | Replay | -0.326 | 1.612 | 0.609 |
| Blank-screen imagery | Transport v2 | 0.038 | 1.360 | 0.391 |
| Blank-screen imagery | Registered density ridge | 0.032 | 1.364 | 0.313 |
| Blank-screen imagery | Elastic ridge | 0.002 | 1.385 | 0.297 |
| Blank-screen imagery | MultiMatch ridge | -0.005 | 1.390 | 0.297 |
| Blank-screen imagery | Replay | -0.197 | 1.523 | 0.359 |

Replay's rank behaviour and probability quality diverged. It had the highest
repeated-viewing top-one credit, but several confidently wrong predictions
dominated the proper log score. Inner-held-out calibration temperatures ranged
from 5.99 to 32.03 across the eight task-fold fits and were not boundary
solutions. The failure is therefore not numerical underflow; it is evidence
that calibration learned from the small outer-training domains did not
generalize reliably to held-out participants and items.

The crossed-bootstrap information interval included zero for every
blank-screen imagery method. The Replay interval was -0.859 to 0.392 bits;
Transport's was -0.105 to 0.173; registered density's was -0.040 to 0.110.
Using every deposited post-test fixation rather than the prespecified 3000-ms
window changed Replay to -0.153 bits, log loss 1.492, and top-one credit 0.359.

## Chronology and null behaviour

Registered density and elastic matching were invariant to order within
numerical precision. On repeated viewing, signal-minus-shuffled information was
+0.021 bits for Replay, -0.003 for Transport, and +0.238 for the learned
MultiMatch composite. Thus at least one GazeWeave engine responded in the
declared direction, but MultiMatch remained the clearest order-sensitive
comparator.

The empirical false-positive rate was 0.0508 for every calibrated method under
the frozen wrong-item and generic-gaze null. Registration did not rescue the
null: density top-one credit decreased from 0.302 raw to 0.294 registered. All
engine candidate fits converged, and every fold audit had zero participant and
item overlap.

## Alignment diagnostics and registration

Replay attributed 82.9% of repeated-viewing mass and 83.6% of blank-screen
imagery mass to replay, with the remainder assigned to background. Transport
coverage was 67.1% and 68.9%. Replay spatial residuals were 71.9 px for repeated
viewing and 122.8 px for imagery; Transport residuals were 73.1 and 105.3 px.

The learned blank-screen source-to-encoding scale was 1.46-1.83 across folds,
consistent with recall gaze occupying roughly 55-68% of the encoding
dispersion. Repeated-viewing scales were 0.87-1.10. These are held-out nuisance
parameters and useful diagnostics, but short-trial content and calibration
offset remain possible contributors; they are not separate inferential wins.

## Fair-baseline comparison and computation

Registered density was the strongest baseline. Pooling the two signal tasks,
its log loss was 1.213 versus 1.279 for Transport and 1.568 for Replay. Replay
minus registered density was -0.511 information bits (95% crossed-bootstrap
interval -1.478 to 0.123), -0.354 in log-loss advantage (-0.981 to 0.102), and
+0.055 in top-one credit (-0.094 to 0.211). No superiority rule passed.

Across all eight task-fold fits, Replay took 8.3 seconds, Transport v2 1970.1
seconds, and all baselines 111.9 seconds. Replay remains computationally
parsimonious; Transport remains expensive in the current R implementation.

## Bounded conclusion

The court supports three narrower statements:

1. Replay is a fast, directional, probabilistically explicit model with useful
   posterior alignments and strong rank performance in repeated viewing.
2. Transport v2 is the better-supported GazeWeave engine in this real court,
   but its symmetry and cost make it a deliberate model choice rather than a
   universal default.
3. Registered density remains a powerful and parsimonious item-identification
   baseline; MultiMatch remains the clearest diagnostic for order destruction.

It does **not** support saying that GazeWeave has better power than MultiMatch
or density, that Replay probabilities are reliably calibrated across held-out
participants and items, or that one engine should be silently selected for all
scientific questions.

The public API should therefore require an explicit invariance contract:

- choose Replay when the directional local-replay model is the hypothesis;
- choose Transport v2 for a symmetric partial-correspondence analysis;
- report registered density and a learned MultiMatch composite as fair
  comparators;
- treat `gaze_info_bits` as the primary proper-score endpoint and use rank,
  coverage, chronology, warp, and braid plots as diagnostics.

The next gate is independent external replication or substantially larger
training domains for calibration. Retuning the temperature against this frozen
court would invalidate its role as held-out evidence.
