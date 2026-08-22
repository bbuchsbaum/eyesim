# Transport versus density: saliency and graded-response diagnostic

This post-court analysis asks whether canonical Transport shows the expected
increase in study-to-retrieval item information as stimulus saliency rises,
whether that increase exceeds simple density baselines on the same trials, and
whether either representation predicts graded recognition responses.

The frozen GW-18.10 result is unchanged. Transport and each density method use
the same 1,295 out-of-fold trials, candidates, participant folds, and item
folds. Trial-linked results remain local and Git-ignored.

## Estimands

The saliency diagnostic is the fitted change in `gaze_info_bits` from saliency
20 to 100. It is reported overall and separately for old and lure probes for:

- canonical Transport;
- raw density with sigma 80 and response-blind calibration;
- raw density with sigma 160 and response-blind calibration;
- the nested registered-density baseline.

Transport-minus-density differences are paired within trial. Uncertainty uses
the crossed participant/item bootstrap. A quality-adjusted slope controlling
effective fixation count and retained duration is reported as a sensitivity,
not substituted for the bare diagnostic.

The study records a 1--4 newness rating: 1 is least new and 4 is most new;
response 0 is missing. The model uses the equivalent `oldness = 5 - newness`
coding, so a positive gaze coefficient means greater oldness-confidence. All
four response categories are preserved.

For old and lure probes separately, proportional-odds models are trained on
rows sharing neither participants nor items with the evaluation fold. The base
contains saliency, log effective fixations, and retained duration. Held-out
ordinal log-loss gains compare Transport, each density score, their direct
predictions, and Transport-plus-density nested models. Positive bits mean the
second model predicts the four-category response better.

Run after loading the package:

```r
source("inst/validation/gaze-weave-transport-density-diagnostic.R")
run_gaze_weave_transport_density_diagnostic()
```

## Result

On old probes, Transport has the expected positive saliency diagnostic:
information rises by `0.03779` bits from saliency 20 to 100. Density sigma 80,
density sigma 160, and registered density instead have slopes `-0.00769`,
`-0.00909`, and `-0.00741` bits. The paired Transport-minus-density advantages
are about `0.045`--`0.047` bits, but their crossed 95% intervals include zero;
for sigma 80 the difference is `0.04548` (`-0.04990`, `0.13422`). The point
estimate favors Transport as the more diagnostically appropriate method, but
does not establish a reliable advantage. The five Transport cell means are
not strictly monotone (`0.0059`, `0.0282`, `0.0464`, `0.0158`, and `0.0603`
bits); saliency 80 contains a visible dip, so the positive linear trend should
not be described as a dose-response law.

On lure probes, Transport rises by `0.01534` bits and the density variants by
about `0.023`--`0.025` bits. None of those slopes or paired differences is
reliable. Across old and lure probes together, Transport's slope is `0.02630`
bits (`-0.05338`, `0.09634`) versus `0.00849` bits for density sigma 80
(`-0.04188`, `0.05637`). Quality adjustment barely changes the point estimates.

Transport does not reliably predict the 1--4 newness rating out of fold. Its
gain over the saliency-and-gaze-quality base is `0.00059` bits for old probes
(`-0.01163`, `0.01338`) and `0.00312` bits for lures (`-0.00892`, `0.01484`).
It does not reliably beat density sigma 80: the direct gains are `0.00193` bits
for old probes (`-0.01038`, `0.01538`) and `0.03141` bits for lures
(`-0.00998`, `0.08234`). The larger lure difference reflects density failure
in one outer fold, not a stable positive Transport effect. Adding Transport
beyond density sigma 80 also fails to improve held-out confidence prediction.

Thus Transport has a scientifically encouraging *point-estimate* saliency
signature for old probes, whereas the simple density methods do not. The
available participant/item support is insufficient to show that its slope is
reliably positive or reliably greater than density, and neither representation
predicts graded newness robustly.
