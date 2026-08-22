# GazeWeave resolution-invariant replication verdict

Final seed: `20260818`
Protocol version: 2
Decision: six of seven gates passed; the blank-screen imagery gate failed;
there is **no real-data engine default** and superiority is not supported.

## What changed

Replay now scores candidates with mean log predictive density per normalized
duration bin. Its raw total HMM log likelihood remains an alignment diagnostic,
and inner-fold calibration estimates only a residual temperature on the mean
score. This makes the unit-centred temperature prior invariant to the chosen
duration-grid resolution.

Nothing else in the court changed. The models, thresholds, folds, candidate
pools, controls, and fair baselines were frozen before outcome scoring. The
replication used eight participants and eight items absent from the original
court, though both courts draw from the same deposited Wynn dataset.

## Gate result

Six gates passed: repeated viewing, local-order ablation, empirical null error,
null-safe registration, engine convergence, and two-way cross-fitting. The
blank-screen gate required at least one GazeWeave engine to have positive mean
information and top-one credit above 0.25. Top-one credit passed, but mean
information was -0.072 bits for Replay and -0.002 bits for Transport v2.

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
| Repeated viewing | Registered density ridge | 0.462 | 1.066 | 0.516 |
| Repeated viewing | Transport v2 | 0.384 | 1.120 | 0.531 |
| Repeated viewing | Replay | 0.229 | 1.227 | 0.453 |
| Repeated viewing | Elastic ridge | 0.015 | 1.376 | 0.203 |
| Repeated viewing | MultiMatch ridge | -0.672 | 1.852 | 0.422 |
| Blank-screen imagery | Elastic ridge | 0.022 | 1.371 | 0.422 |
| Blank-screen imagery | Transport v2 | -0.002 | 1.387 | 0.344 |
| Blank-screen imagery | Registered density ridge | -0.011 | 1.394 | 0.266 |
| Blank-screen imagery | MultiMatch ridge | -0.045 | 1.417 | 0.156 |
| Blank-screen imagery | Replay | -0.072 | 1.436 | 0.328 |

The correction repaired the specific repeated-viewing calibration pathology in
the original court: Replay changed from -0.326 to +0.229 mean information bits
on a disjoint participant-item sample, while retaining useful rank performance.
Its residual temperatures were on the intended per-bin scale rather than the
raw 48-bin total scale.

The harder imagery result did not become positive. Crossed-bootstrap 95%
information intervals were -0.271 to 0.079 bits for Replay, -0.177 to 0.135 for
Transport, and -0.090 to 0.099 for registered density. Using every deposited
post-test fixation yielded -0.082 Replay bits and did not alter the conclusion.

## Chronology, nulls, and diagnostics

Repeated-viewing signal exceeded shuffled-order information by 0.006 bits for
Replay and 0.014 for Transport. Registered density was invariant within
numerical precision. The empirical false-positive rate was 0.0508 for every
calibrated method, all engine candidate fits converged, and every participant
and item overlap receipt was zero.

Replay attributed 86.8% of repeated-viewing mass and 76.8% of imagery mass to
replay. Transport coverage was 68.7% and 69.2%. Replay spatial residuals were
69.9 px for repeated viewing and 111.2 px for imagery. These quantities explain
the fitted alignments; they are not additional inferential endpoints.

## Fair baselines and computation

Registered density remained strongest by pooled log loss: 1.230 versus 1.254
for Transport and 1.332 for Replay. Replay minus registered density was -0.147
information bits (95% crossed-bootstrap interval -0.340 to 0.012) and -0.102
in log-loss advantage (-0.241 to 0.010). No superiority rule passed.

Across all eight task-fold fits, Replay took 12.0 seconds, Transport v2 2516.6
seconds, and all baselines 141.0 seconds. Replay is computationally parsimonious
and its posterior alignment is easier to interpret, but runtime and mechanistic
clarity do not substitute for held-out item evidence.

## Bounded conclusion

The resolution-invariant correction is retained because it repairs an
algorithmic scale defect and transfers to a disjoint repeated-viewing control.
The package still requires explicit engine choice:

- choose Replay for a directional local-replay hypothesis and posterior braid;
- choose Transport v2 for symmetric partial correspondence;
- report registered density and a learned MultiMatch composite as fair
  predictive comparators; and
- use MultiMatch dimensions when the scientific question concerns a specific
  component such as direction or duration.

The evidence does not support a universal default, a GazeWeave power advantage,
or positive item-information reinstatement during blank-screen imagery in this
dataset. The next defensible gate is an external dataset or a prospectively
larger task-matched imagery study, not another post-result adjustment to this
court.
