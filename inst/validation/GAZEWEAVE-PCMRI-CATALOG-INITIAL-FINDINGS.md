# pcmri catalog audit and initial reconstruction

Date: 2026-08-20
Protocol: `GAZEWEAVE-PCMRI-CATALOG-REPLICATION.md`
Status: density reconstruction, common-support comparison, and clustered
behavioral follow-up complete; all-four-presentation MultiMatch reconstruction
pending

## What the deployed report actually tests

The report labelled PROBE-WINDOW uses the full test-picture, mask, and delay
interval of about three seconds. This corresponds to the existing GazeWeave
`combined` interval, `[0, 3000)` ms, not the isolated 500-ms probe.

Its primary density quantity is matched Pearson-like map similarity minus
similarity to other studied images from the same participant. This is an
above-chance compatibility contrast. It is not the same estimand as held-out
five-candidate item information.

The report's metric-generation CSVs are not deployed and its embedded source
does not state the density bandwidth. Repository-local same-paradigm scripts
use duration-weighted density, `sigma = 80` px, an 80 by 60 grid, and exhaustive
within-participant other-image references. The reconstruction uses that
documented configuration.

## Audit of the rendered findings

The numerical tables support:

- positive own-template density reinstatement, about +0.028 matched minus
  nonmatch;
- lower density reinstatement for lures than old probes, about -0.013;
- a weak positive density-saliency slope;
- a positive old-trial response association for density difference in the
  catalog-compatible model, coefficient about +0.222 (`p = .034`);
- nominally positive chance differences for all five MultiMatch dimensions;
- no old-trial response association for any separate MultiMatch dimension;
- no response association for the report's combined MultiMatch shape score
  (`p = .537`);
- strong study-repeat density stability, about +0.114 above its baseline.

Several prose summaries contradict those current tables, apparently retaining
results from an earlier render. In particular, the page's concluding claims
that only position and duration reinstate, that vector and direction predict
recognition, and that shape predicts recognition do not follow from the
displayed probe-window tables.

## Sigma-80 reconstruction from the raw exports

The reconstruction used all 46 participants, all four study presentations in
each density template, and the combined retrieval interval. A participant-item
was retained only when all four study presentations existed. Of 2,056 mapped
old/lure trials, 2,014 had a finite density map after the existing fixation
filter; 970 finite old trials also had a nonzero response.

| Question | Public report | Reconstruction |
|---|---:|---:|
| Mean density difference | +0.028, `p < .001` | +0.0382, `p < .001` |
| Lure minus old | -0.013, `p = .035` | -0.0184, `p = .040` |
| Old saliency slope per 10% | +0.003, `p = .058` | +0.00684, `p = .003` |
| Lure by saliency | -0.002, `p = .364` | -0.00584, `p = .065` |
| Old response association | +0.222, `p = .034` | +0.1765, `p = .128` |

The spatial reinstatement result and old-lure ordering reproduce in sign and
rough scale. The old-trial behavioral association does not reproduce at the
conventional threshold under this reconstruction. This discrepancy could
reflect undeployed preprocessing/permutation details, the unspecified report
bandwidth, or ordinary sampling/model sensitivity; it must not be described as
an exact replication failure until the original metric-generation files are
available.

## Common-support method landscape

The completed GW-13 court provides a fair P4-to-combined-window comparison on
1,295 trials from 36 participants with the same five candidate identities for
every method and participant/item cross-fitting.

| Method | Mean information bits | Log loss | Mean rank | Top-one |
|---|---:|---:|---:|---:|
| Transport v2 | +0.00728 | 1.60439 | 2.851 | 0.241 |
| MultiMatch position | +0.00604 | 1.60525 | 2.842 | 0.242 |
| Density sigma 80 | +0.00594 | 1.60532 | 2.844 | 0.235 |
| Replay | +0.00316 | 1.60725 | 2.836 | 0.239 |
| Registered density ridge | +0.00242 | 1.60776 | 2.869 | 0.234 |
| MultiMatch ridge | -0.00083 | 1.61002 | 2.889 | 0.224 |

All values are close to zero and all methods are close to one another. Chance
top-one credit is 0.20. These descriptive values do not establish a power
advantage for GazeWeave, density, or MultiMatch.

In catalog-compatible old-trial response models using each method's held-out
information score, Replay was positive (+0.355 standardized log odds,
`p = .023`), while Transport (+0.094, `p = .544`), density sigma 80 (+0.131,
`p = .400`), registered density (+0.135, `p = .390`), MultiMatch position
(-0.120, `p = .479`), and the MultiMatch ridge composite (-0.218, `p = .193`)
were not. The models were singular and the Replay comparison was not a frozen
primary test.

## Replay behavioral follow-up

The follow-up trained a response model on participants and items disjoint from
each evaluation fold. Every method was compared with the same saliency-only
model. Improvement in held-out response log loss was converted to bits per
trial, with paired participant-by-item bootstrap intervals over 2,000 draws.

| Method | Held-out behavioral bits | Crossed 95% interval |
|---|---:|---:|
| Density sigma 80 | -0.00168 | [-0.00728, +0.00246] |
| Replay | -0.00388 | [-0.03912, +0.03161] |
| Registered density | -0.00773 | [-0.02262, +0.00324] |
| MultiMatch position | -0.00791 | [-0.02955, +0.00976] |
| Transport v2 | -0.00914 | [-0.02668, +0.00333] |
| MultiMatch composite | -0.01318 | [-0.04386, +0.01168] |

No method improved held-out prediction over saliency alone. Replay was not
better than density sigma 80: Replay minus density was -0.00220 bits, ordinary
95% interval [-0.03741, +0.03650], family interval
[-0.05035, +0.05671]. Every other Replay-versus-method interval also included
zero.

A separate crossed bootstrap of the all-trial standardized association reduced
Replay's coefficient to +0.304 with interval [-0.130, +0.778] and bootstrap
`p = .166`. Thus the earlier `p = .023` mixed-model result does not survive
cluster-aware uncertainty or out-of-fold behavioral prediction. The most
likely explanation is an unstable association in a highly imbalanced outcome:
the four evaluation folds contain 88%--92% old responses, and strict crossed
training leaves only 169--182 old trials per fold.

## Current conclusion

The public analysis strengthens the evidence that this dataset contains a
small but reliable spatial matched-minus-nonmatch signal when four study
presentations are pooled. It does not overturn the earlier finding that
item-identifying evidence during the combined recognition interval is weak.

GazeWeave is not presently more sensitive than density or MultiMatch at item
identification on this retrieval dataset. Replay's preliminary old-response
association did not survive its predictive and crossed-bootstrap follow-up and
therefore is not evidence of superior behavioral sensitivity.
