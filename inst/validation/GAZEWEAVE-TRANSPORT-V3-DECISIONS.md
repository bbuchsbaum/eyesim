# Transport v3 prospective decision table

Protocol: `gazeweave-transport-v3/3.0.2`
Frozen: 2026-08-20
Retrieval outcomes: sealed

| Gate | Evidence examined | Pass rule | Failure action |
|---|---|---|---|
| Estimator lock | Public v2 golden fixture | Every v2 summary/profile value agrees within `1e-8`; fixture inputs and expected files match their MD5 manifest | Stop; distinguish source drift from intended v3 change before implementation |
| Chronology | Analytic and synthetic paths of length 1-64 | Identical `<= 1e-10`, reversal `>= 0.95`, length span `<= 0.05`, dilation `<= 1e-10`, finite declared short-path behavior | Revise edge-normalized residual; do not inspect study or retrieval scores |
| Representation | Split/merge, duplicate locations, units, candidate order, batch size | Energy/coverage/selection within `1e-8`, unit conversion within `1e-6`, candidate ranks identical | Revise representation or backend; do not relax thresholds |
| Coverage | 12- versus 24-node quadrature and coherent partial replay | Energy change `<= 1e-3`, information change `<= 0.01` bits; low-coverage replay is intermediate; convenient-pair selection pays non-zero selection cost | Revise integration or selection model before calibration |
| Backend parity | Random, sparse, extreme, short, and public real-path fixtures | Energy `<= 1e-8`, mass `<= 1e-10`, barycentric RMSE `<= 1e-6` degrees, warp `<= 1e-10`, rank exact | Keep R reference authoritative and repair/fallback the optimized backend |
| Optimization | Entropic initialization and optional scientific polish | Polished scientific energy never worsens beyond `1e-10`; feasibility `<= 1e-8`; gap reported; structural-start ranks stable | Retain entropic backend and record polish as a no-go audit option |
| Calibration/null | Inner out-of-fold predictions and generated null rows | Reliability log loss no worse than temperature-only by `1e-10`; null information within `0.02` bits, top-one within `0.02`, ECE `<= 0.10` | Revise nested calibration without retrieval outcomes |
| Synthetic science | Frozen perturbation factorial | Fixed-density order destruction reduces v3 evidence; partial replay is intermediate; false-positive rate `<= 0.075`; finite/converged `>= 0.99` | Return to responsible algorithm gate with thresholds unchanged |
| Pair performance | Frozen 15-by-15 fixture | Median `<= 15` ms and speedup `>= 4x` versus the approximately 60 ms GW-14 reference | Do not promote as the exhaustive default; optimize without estimator drift |
| Fold performance | Frozen representative full fold | Speedup `>= 2x` beyond GW-14, bounded peak memory, all mandatory alignments converge | Do not open the study court |
| Study positive control | Response-blind P1-P2/P1-P3/P1-P4/adjacent comparisons | Intact mean information `> 0`, intact-control contrast `> 0`, and at least one crossed 95% interval lower bound `> 0` | Retain the negative result, keep retrieval sealed, publish bounded no-go/experimental role |
| Retrieval | Exact frozen 0-3000 ms score bundle, then outcomes | No post-open tuning; compare out-of-fold log loss on common support with crossed uncertainty; GLMMs diagnostic | Report the frozen result without rescue analyses |
| Product role | All prior ledgers | Promote only after invariance, distinctiveness, calibration, positive control, and runtime pass | Stable tie: explanatory symmetric engine; instability/runtime failure: experimental reference only |

All methods use identical outer folds, candidate pools, priors, registration
opportunities, common-support rows, and uncertainty resamples. `gaze_info_bits`
is the sole default inferential endpoint. Coverage and alignment fields remain
diagnostics.
