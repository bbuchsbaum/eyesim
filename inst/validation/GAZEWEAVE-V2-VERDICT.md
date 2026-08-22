# GazeWeave v2 synthetic-court verdict

Final seed: `20260816`
Protocol version: 1 plus debug addenda A-C
Candidate pool: six items
Evaluation rows: 36 matched signal and 36 independent central-background null
Verdict: **advance to real-data controls; provisional directional default is Replay**

## Gate result

All seven frozen advancement gates passed: representation and stability, null
error, predictive signal, the fixed-density order control, contraction
recovery, null-registration safety, and engine convergence. Every engine score
was finite and declared converged. The fitted contraction scales were `0.7714`
and `0.7907` around the true value `0.78`.

All exact invariances had zero observed difference: fixation splitting and
merging, uniform temporal dilation, the declared detector coalescing variant,
and equivalent degree/pixel coordinates. Transport's two-start scientific
spread was `0.000128` against a `0.02` tolerance. Replay's 32/48/64-bin
contrastive likelihood spread was `0.0070` nats per bin against `0.15`.

## Predictive result

On the replay- and transport-favouring families, both GazeWeave engines were
near the six-candidate ceiling of `log2(6) = 2.585` bits. On the
fixed-density/order family:

| Method | Mean information | Log loss | Top-one credit |
|---|---:|---:|---:|
| Replay | 2.574 bits | 0.0078 | 1.00 |
| Transport v2 | 1.759 bits | 0.5728 | 1.00 |
| Nested MultiMatch | 2.574 bits | 0.0079 | 1.00 |
| Nested density | 0.000 bits | 1.7918 | 0.167 |

This establishes the intended order-versus-density distinction, but it does
**not** establish GazeWeave superiority over a fair learned MultiMatch
composite. Replay and MultiMatch were essentially tied in this particular
synthetic order control. At the empirical 4.17% null error, power was 1.00 for
MultiMatch, 0.722 for Transport, and 0.667 for Replay and density. The easy
spatial families saturated several methods, so these power values should not
be generalized.

Across the two directional families used by the default rule, pooled log loss
was `0.0039` for Replay and `0.2876` for Transport. Replay therefore satisfies
the prespecified noninferiority rule and wins the local-replay diagnostic.

## Parsimony and computation

Replay fitted and scored the full court in `1.85` seconds. Transport v2 took
`391.9` seconds using four candidate workers. The complete baseline court took
`18.8` seconds. Serialized fit sizes were 3.95 MB for Replay, 5.72 MB for
Transport, and 66.33 MB for the long-format baseline audit.

Replay calibration used two inner item folds inside each outer-training set and
a fixed weak prior centred on temperature one. This repaired a debug run in
which separable training candidates drove the unconstrained temperature to its
lower bound and produced anti-conservative null ties. The final log-domain
score implementation had no infinite information or log-loss values.

Replay also has the more economical scientific explanation for
encoding-to-recall: background occupancy, local continuation, restart, and
spatial emission noise are probabilistic parameters, and its braid is a state
posterior. Transport remains the appropriate symmetric comparison model and a
valuable sensitivity analysis, but its coupling diffuseness remains an
optimizer property and its fixed-mass continuation is much more expensive.

## Bounded conclusion

The synthetic evidence supports advancing Replay and Transport v2 to real
positive and negative controls. It supports Replay as the provisional default
for a directional encoding-to-recall API because Replay was stable, much
faster, probabilistically interpretable, and better than Transport in the
prespecified directional pool.

It does not show that GazeWeave has greater statistical power than MultiMatch,
that the candidate probabilities are calibrated in real data, or that the
synthetic effect sizes resemble looking-at-nothing experiments. Those are
GW-7 questions. Transport should not be removed: symmetry is scientifically
appropriate when neither path is privileged, and the final choice remains an
explicit invariance contract rather than a universal ranking.
