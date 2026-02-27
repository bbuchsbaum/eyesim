# Latent Transforms for Template Similarity

``` r
library(eyesim)
library(dplyr)
library(ggplot2)
```

## Why latent transforms?

[`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md)
compares density maps in raw screen space. But when encoding and
retrieval differ by a linear transform — different screen sizes,
calibration drift between sessions, or systematic participant
differences — raw correlation can underestimate true reinstatement.

Latent transforms address this by mapping both sides into a shared space
before computing similarity. You add them to
[`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md)
with two arguments; the core API stays unchanged.

## A concrete example

Suppose two devices record eye-movements on the same stimuli, but one
device scales coordinates differently. Let’s simulate this:

``` r
set.seed(1)

make_density <- function(vec) {
  fg <- fixation_group(
    x = vec[1:2] * 50 + 50,
    y = vec[3:4] * 50 + 50,
    onset = c(0, 100),
    duration = c(100, 100)
  )
  eye_density(fg, sigma = 30, xbounds = c(0, 100), ybounds = c(0, 100))
}

base_vecs <- replicate(8, rnorm(4), simplify = FALSE)
```

Now create a “source” set that is a scaled version of the reference — as
if recorded on a differently-calibrated device:

``` r
scale_factors <- c(2, 0.6, 1.7, 0.8)
source_vecs <- lapply(base_vecs, function(v) v * scale_factors)

ref_tab <- tibble(id = seq_along(base_vecs),
                  density = lapply(base_vecs, make_density))
source_tab <- tibble(id = seq_along(source_vecs),
                     density = lapply(source_vecs, make_density))
```

## Comparing raw vs. CORAL-transformed similarity

Without any transform, the scaling distortion reduces similarity:

``` r
raw <- template_similarity(ref_tab, source_tab,
                           match_on = "id",
                           permutations = 0,
                           method = "cosine")
cat("Raw mean similarity:", round(mean(raw$eye_sim), 3), "\n")
#> Raw mean similarity: 0.86
```

CORAL (CORrelation ALignment) whitens the source domain and re-colors it
with the reference covariance, correcting for the scaling difference:

``` r
coral_res <- template_similarity(
  ref_tab, source_tab,
  match_on = "id",
  permutations = 0,
  method = "cosine",
  similarity_transform = coral_transform,
  similarity_transform_args = list(comps = 4, shrink = 1e-6)
)
cat("CORAL mean similarity:", round(mean(coral_res$eye_sim), 3), "\n")
#> CORAL mean similarity: 0.642
```

![Paired comparison of raw vs. CORAL-transformed similarity for each
item. CORAL recovers higher similarity by correcting the scaling
distortion.](latent-transforms_files/figure-html/plot-comparison-1.png)

Paired comparison of raw vs. CORAL-transformed similarity for each item.
CORAL recovers higher similarity by correcting the scaling distortion.

CORAL recovers substantially higher similarity by correcting the
covariance mismatch between the two “devices.”

## Available transforms

| Transform                                                                                         | Supervised?               | Best for                                                 |
|:--------------------------------------------------------------------------------------------------|:--------------------------|:---------------------------------------------------------|
| [`latent_pca_transform()`](https://bbuchsbaum.github.io/eyesim/reference/latent_pca_transform.md) | No                        | Dimensionality reduction, mild noise smoothing           |
| [`coral_transform()`](https://bbuchsbaum.github.io/eyesim/reference/latent_pca_transform.md)      | No                        | Device/participant shifts (covariance-level differences) |
| [`cca_transform()`](https://bbuchsbaum.github.io/eyesim/reference/latent_pca_transform.md)        | Yes (needs matched pairs) | Item-level alignment when pairings are reliable          |

All are passed to
[`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md)
via the `similarity_transform` argument:

``` r
template_similarity(
  ref_tab, source_tab,
  match_on = "id",
  similarity_transform = cca_transform,
  similarity_transform_args = list(comps = 10, shrink = 0.01)
)
```

## Choosing a transform

- **PCA** is the safest default — dimensionality reduction with no
  assumptions about domain differences. Use when you want noise
  smoothing without domain adaptation.
- **CORAL** is unsupervised and assumes linear, covariance-level
  differences between domains. Good for calibration drift or different
  screen sizes.
- **CCA** is supervised and leverages matched pairs to find shared
  latent axes. Set `comps` modestly (5–15) and add `shrink` to stabilize
  small-N fits.

After the transform,
[`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md)
still uses your chosen `method` (cosine, pearson, fisherz, etc.) on the
transformed vectors.

## Notes and limitations

- Multiscale densities are supported only when all scales share the same
  grid size; otherwise latent transforms will error.
- Regularization (`shrink`) and component count (`comps`) affect
  stability on small-N data — tune as needed.
- CORAL and CCA are linear; they will not correct nonlinear spatial
  warps. Keep `comps` low to reduce overfitting, especially with few
  pairs.
- See
  [`vignette("eyesim")`](https://bbuchsbaum.github.io/eyesim/articles/eyesim.md)
  for the core
  [`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md)
  workflow without transforms.
