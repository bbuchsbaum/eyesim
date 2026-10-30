# Comparing Eye-Movement Patterns

``` r

library(eyesim)
library(dplyr)
library(ggplot2)
library(patchwork)
```

When you study an image, you move your eyes to different locations in a
particular sequence. If you later see that same image during a memory
test, do you look at the same places? This kind of “eye-movement
reinstatement” is a powerful marker of memory, and `eyesim` gives you
the tools to measure it.

This vignette walks you through the core workflow: representing
fixations, computing density maps, and measuring similarity between
fixation patterns across experimental conditions. It ends by scoring
revisits, meaning returns to earlier locations within a single viewing.

## How do you represent fixations?

A `fixation_group` holds a set of eye fixations — each defined by an x/y
screen position, an onset time (when the fixation began), and a duration
(how long the eye stayed there). You can create one directly:

``` r

fg <- fixation_group(
  x = c(-100, 0, 100),
  y = c(0, 100, 0),
  onset = c(0, 10, 60),
  duration = c(10, 50, 100)
)
```

![Three fixations. Point size reflects duration; color reflects onset
time (dark purple = early, green =
late).](eyesim_files/figure-html/plot-fixation-group-1.png)

Three fixations. Point size reflects duration; color reflects onset time
(dark purple = early, green = late).

Point area is proportional to how long each fixation lasted. Color
indicates when it occurred, from dark purple (early) to green (late).

Here is a more realistic group with 25 randomly placed fixations:

``` r

set.seed(42)
fg <- fixation_group(
  x = runif(25, 0, 100),
  y = runif(25, 0, 100),
  onset = cumsum(runif(25, 0, 100)),
  duration = runif(25, 50, 300)
)
```

![25 randomly placed
fixations.](eyesim_files/figure-html/plot-random-fixations-1.png)

25 randomly placed fixations.

## How do you visualize fixation density?

Individual points can be hard to interpret. Density maps show you where
fixations cluster. The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method supports
several display styles:

``` r

p1 <- plot(fg, typ = "contour", xlim = c(-10, 110), ylim = c(-10, 110), bandwidth = 35)
p2 <- plot(fg, typ = "raster", xlim = c(-10, 110), ylim = c(-10, 110), bandwidth = 35)
p3 <- plot(fg, typ = "filled_contour", xlim = c(-10, 110), ylim = c(-10, 110), bandwidth = 35)
p1 + p2 + p3
```

![Three density visualizations: contour, raster, and filled
contour.](eyesim_files/figure-html/density-types-1.png)

Three density visualizations: contour, raster, and filled contour.

The `bandwidth` parameter controls the smoothing level. Higher values
blur out fine detail and emphasize broad patterns:

``` r

p1 <- plot(fg, typ = "filled_contour", xlim = c(-10, 110), ylim = c(-10, 110), bandwidth = 20)
p2 <- plot(fg, typ = "filled_contour", xlim = c(-10, 110), ylim = c(-10, 110), bandwidth = 60)
p3 <- plot(fg, typ = "filled_contour", xlim = c(-10, 110), ylim = c(-10, 110), bandwidth = 100)
p1 + p2 + p3
```

![Effect of bandwidth: narrow (20), medium (60), and wide
(100).](eyesim_files/figure-html/density-bandwidth-1.png)

Effect of bandwidth: narrow (20), medium (60), and wide (100).

## How do you compare two fixation patterns?

To quantify how similar two fixation patterns are, convert them to
`eye_density` maps and call
[`similarity()`](https://bbuchsbaum.github.io/eyesim/reference/similarity.md).
Here we create two patterns that share roughly half their fixation
locations:

``` r

set.seed(123)
x_shared <- runif(12, 0, 100)
y_shared <- runif(12, 0, 100)

fg1 <- fixation_group(
  x = c(x_shared, runif(13, 0, 100)),
  y = c(y_shared, runif(13, 0, 100)),
  onset = cumsum(runif(25, 0, 100)),
  duration = runif(25, 50, 300)
)

fg2 <- fixation_group(
  x = c(x_shared, runif(13, 0, 100)),
  y = c(y_shared, runif(13, 0, 100)),
  onset = cumsum(runif(25, 0, 100)),
  duration = runif(25, 50, 300)
)
```

![Two fixation patterns sharing roughly half their
locations.](eyesim_files/figure-html/plot-pair-1.png)

Two fixation patterns sharing roughly half their locations.

Now convert to density maps and compute their similarity:

``` r

ed1 <- eye_density(fg1, sigma = 50, xbounds = c(0, 100), ybounds = c(0, 100))
ed2 <- eye_density(fg2, sigma = 50, xbounds = c(0, 100), ybounds = c(0, 100))
similarity(ed1, ed2)
#> [1] 0.7760392
```

The default metric is the Pearson correlation. Several alternatives are
available:

``` r

methods <- c("pearson", "spearman", "fisherz", "cosine", "l1", "jaccard")
if (requireNamespace("energy", quietly = TRUE)) {
  methods <- c(methods, "dcov")
}
results <- sapply(methods, function(m) similarity(ed1, ed2, method = m))
data.frame(method = methods, similarity = round(unlist(results), 4))
#>            method similarity
#> pearson   pearson     0.7760
#> spearman spearman     0.7714
#> fisherz   fisherz     1.0353
#> cosine     cosine     0.9934
#> l1             l1     0.9484
#> jaccard   jaccard     0.9869
#> dcov         dcov     0.7293
```

## How do you analyze a full experiment?

In a typical memory study, participants view images during encoding and
again during retrieval. You want to compare fixation patterns between
these phases for the same image and control for non-specific similarity.

Let’s simulate a small experiment: 3 participants, 20 images, encoding +
retrieval.

``` r

head(df)
#>          x        y     onset duration image    phase participant
#> 1 11.37817 50.54517  62.07862 197.0510  img1 encoding          s1
#> 2 68.42647 19.38365 235.32139 204.0351  img1 encoding          s1
#> 3 99.25088 63.69041 373.98848 263.6462  img1 encoding          s1
#> 4 53.49936 68.78001 539.99684 112.1470  img1 encoding          s1
#> 5 96.66141 64.01908 642.50973 233.1542  img1 encoding          s1
#> 6 67.14276 35.78854 693.41893 307.2941  img1 encoding          s1
cat("Rows:", nrow(df),
    " | Participants:", length(unique(df$participant)),
    " | Images:", length(unique(df$image)))
#> Rows: 783  | Participants: 3  | Images: 20
```

Wrap the raw data in an `eye_table`, which groups fixations by the
variables that define your experimental design:

``` r

eyetab <- eye_table("x", "y", "duration", "onset",
                    groupvar = c("participant", "phase", "image"),
                    data = df)
eyetab
#> Eye table: 120 groups, 783 total fixations
#>   origin: (640, 640)
#> # A tibble: 120 × 4
#>    participant phase    image fixgroup           
#>    <chr>       <chr>    <chr> <list>             
#>  1 s1          encoding img1  <fxtn_grp [10 × 6]>
#>  2 s1          encoding img10 <fxtn_grp [4 × 6]> 
#>  3 s1          encoding img11 <fxtn_grp [3 × 6]> 
#>  4 s1          encoding img12 <fxtn_grp [10 × 6]>
#>  5 s1          encoding img13 <fxtn_grp [10 × 6]>
#>  6 s1          encoding img14 <fxtn_grp [4 × 6]> 
#>  7 s1          encoding img15 <fxtn_grp [10 × 6]>
#>  8 s1          encoding img16 <fxtn_grp [7 × 6]> 
#>  9 s1          encoding img17 <fxtn_grp [3 × 6]> 
#> 10 s1          encoding img18 <fxtn_grp [4 × 6]> 
#> # ℹ 110 more rows
```

### Computing density maps by group

[`density_by()`](https://bbuchsbaum.github.io/eyesim/reference/density_by.md)
computes a density map for every combination of your grouping variables:

``` r

eyedens <- density_by(eyetab,
                      groups = c("phase", "image", "participant"),
                      sigma = 100,
                      xbounds = c(0, 100), ybounds = c(0, 100))
```

![Four density maps from the simulated
experiment.](eyesim_files/figure-html/plot-densities-1.png)

Four density maps from the simulated experiment.

``` r

eyedens
#> # A tibble: 120 × 5
#>    phase    image participant fixgroup            density       
#>    <chr>    <chr> <chr>       <list>              <list>        
#>  1 encoding img1  s1          <fxtn_grp [10 × 6]> <ey_dnsty [5]>
#>  2 encoding img1  s2          <fxtn_grp [8 × 6]>  <ey_dnsty [5]>
#>  3 encoding img1  s3          <fxtn_grp [8 × 6]>  <ey_dnsty [5]>
#>  4 encoding img10 s1          <fxtn_grp [4 × 6]>  <ey_dnsty [5]>
#>  5 encoding img10 s2          <fxtn_grp [4 × 6]>  <ey_dnsty [5]>
#>  6 encoding img10 s3          <fxtn_grp [8 × 6]>  <ey_dnsty [5]>
#>  7 encoding img11 s1          <fxtn_grp [3 × 6]>  <ey_dnsty [5]>
#>  8 encoding img11 s2          <fxtn_grp [3 × 6]>  <ey_dnsty [5]>
#>  9 encoding img11 s3          <fxtn_grp [9 × 6]>  <ey_dnsty [5]>
#> 10 encoding img12 s1          <fxtn_grp [10 × 6]> <ey_dnsty [5]>
#> # ℹ 110 more rows
```

### Template similarity analysis

[`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md)
compares each retrieval density map to its matched encoding density map.
It estimates a baseline by permuting image labels and reports the
corrected difference:

``` r

set.seed(1234)
enc_dens <- eyedens %>% filter(phase == "encoding")
ret_dens <- eyedens %>% filter(phase == "retrieval")

simres <- template_similarity(enc_dens, ret_dens,
                              match_on = "image",
                              method = "fisherz",
                              permutations = 50)
```

The result includes three key columns:

- `eye_sim` — raw similarity between matched pairs
- `perm_sim` — average similarity across non-matching pairs (baseline)
- `eye_sim_diff` — the corrected score (raw minus baseline)

Since our data is random, there should be no true reinstatement:

``` r

t.test(simres$eye_sim_diff)
#> 
#>  One Sample t-test
#> 
#> data:  simres$eye_sim_diff
#> t = 0.71493, df = 59, p-value = 0.4775
#> alternative hypothesis: true mean is not equal to 0
#> 95 percent confidence interval:
#>  -0.1111827  0.2347954
#> sample estimates:
#>  mean of x 
#> 0.06180634
```

![Distribution of raw, permuted, and corrected similarity
scores.](eyesim_files/figure-html/plot-similarity-1.png)

Distribution of raw, permuted, and corrected similarity scores.

As expected, the corrected similarity is centered near zero.

**A note on units.** With `method = "fisherz"`, all similarity values
are Fisher z scores (atanh of Pearson r). Convert back to correlations
with [`tanh()`](https://rdrr.io/r/base/Hyperbolic.html) if desired. With
`method = "pearson"` or `"spearman"`, values are correlations directly.

**A note on permutations.** The `permutations` argument is an upper
bound. If fewer non-matching candidates are available, `eyesim` uses all
of them exhaustively.

## What about multiscale analysis?

Choosing a single bandwidth is somewhat arbitrary. You can compute
density at multiple scales by passing a vector of sigma values:

``` r

eyedens_multi <- density_by(eyetab,
                            groups = c("phase", "image", "participant"),
                            sigma = c(25, 50, 100),
                            xbounds = c(0, 100), ybounds = c(0, 100))

enc_multi <- eyedens_multi %>% filter(phase == "encoding")
ret_multi <- eyedens_multi %>% filter(phase == "retrieval")

simres_multi <- template_similarity(enc_multi, ret_multi,
                                    match_on = "image",
                                    method = "fisherz",
                                    permutations = 50)

cat("Single-scale mean:", round(mean(simres$eye_sim_diff), 4), "\n")
#> Single-scale mean: 0.0618
cat("Multiscale mean:  ", round(mean(simres_multi$eye_sim_diff), 4), "\n")
#> Multiscale mean:   0.0351
```

Multiscale analysis provides a more robust similarity estimate by
averaging across spatial resolutions.

## How does similarity evolve over time?

Sometimes you want to know *when* during retrieval participants look at
previously-fixated locations.
[`sample_density_time()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density_time.md)
samples the encoding density at retrieval fixation locations at each
time point:

``` r

ret_eyetab <- eyetab %>% filter(phase == "retrieval")

enc_dens <- enc_dens %>%
  mutate(match_key = paste(participant, image, sep = "_"))
ret_matched <- ret_eyetab %>%
  mutate(match_key = paste(participant, image, sep = "_"))

temporal <- sample_density_time(
  template_tab = enc_dens,
  source_tab = ret_matched,
  match_on = "match_key",
  times = seq(0, 500, by = 50),
  time_bins = c(0, 250, 500),
  permutations = 10,
  permute_on = "participant"
)
```

``` r

temporal %>%
  select(participant, image, bin_1, bin_2,
         perm_bin_1, perm_bin_2, diff_bin_1, diff_bin_2) %>%
  head()
#> # A tibble: 6 × 8
#>   participant image    bin_1   bin_2 perm_bin_1 perm_bin_2 diff_bin_1 diff_bin_2
#>   <chr>       <chr>    <dbl>   <dbl>      <dbl>      <dbl>      <dbl>      <dbl>
#> 1 s1          img1  1.04 e-4 1.02e-4  0.0000966  0.0000992    7.85e-6    2.62e-6
#> 2 s1          img10 1.05 e-4 9.02e-5  0.000102   0.000101     2.83e-6   -1.06e-5
#> 3 s1          img11 9.87 e-5 1.02e-4  0.000101   0.000104    -2.43e-6   -2.07e-6
#> 4 s1          img12 1.000e-4 1.07e-4  0.0000935  0.000106     6.45e-6    9.48e-7
#> 5 s1          img13 9.71 e-5 1.00e-4  0.0000994  0.000100    -2.28e-6   -2.28e-7
#> 6 s1          img14 8.88 e-5 9.42e-5  0.0000979  0.0000977   -9.11e-6   -3.46e-6
```

![Encoding density sampled at retrieval fixation locations over time,
for six example trials.](eyesim_files/figure-html/plot-temporal-1.png)

Encoding density sampled at retrieval fixation locations over time, for
six example trials.

This reveals when during retrieval participants fixate locations that
mattered during encoding — a window into the time course of
memory-guided attention.

## How often does gaze return to earlier locations?

Similarity compares one viewing with another. Revisitation looks inside
a single viewing: how often does gaze return to a location it has
already fixated? Revisits during study have been linked to the
hippocampus and to later memory (Voss et al., 2011; Kragel et al.,
2021).
[`fixation_revisits()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_revisits.md)
scores spatial returns, a measure related to Kragel et al. (2021). It
does not implement Voss et al.’s backtracking-sequence classification:
that study counted sequences such as A-B-C-B-A but classified A-B-C-A as
“other.”

A fixation is a revisit when it lands within `radius` of an earlier
fixation at least `min_lag` fixations back. The default `min_lag = 2`
means a fixation never revisits the one just before it. By default, the
incoming step must be at least `radius`: the previous fixation must be
that far from the current one. This does not establish that gaze left a
region around the original location. With radius 10, positions 0, -9, 9
on one axis count as a return although both later positions remain
within 10 of the original position. Here gaze visits A, B and C, then
returns to A:

``` r

fg_abca <- fixation_group(x = c(10, 40, 40, 11), y = c(10, 10, 40, 9),
                          onset = c(0, 250, 500, 750), duration = rep(200, 4))
rv <- fixation_revisits(fg_abca, radius = 5)
rv
#> Revisits: 1 of 2 opportunities (4 fixations; radius 5, lag 2-Inf, departure required)
#> # A tibble: 4 × 10
#>   index onset opportunity revisit   lag match_index distance nearest_distance
#>   <int> <dbl> <lgl>       <lgl>   <int>       <int>    <dbl>            <dbl>
#> 1     1     0 FALSE       FALSE      NA          NA    NA               NA   
#> 2     2   250 FALSE       FALSE      NA          NA    NA               NA   
#> 3     3   500 TRUE        FALSE      NA          NA    42.4             42.4 
#> 4     4   750 TRUE        TRUE        3           1     1.41             1.41
#> # ℹ 2 more variables: match_distance <dbl>, expected <dbl>
```

Fixations 3 and 4 are *opportunities*: each has an earlier fixation far
enough back, and each arrives from at least `radius` away. Only fixation
4 is a revisit. It returns to fixation 1 (`match_index`), three
fixations back (`lag`). `match_distance` measures distance to this most
recent match; `nearest_distance` measures distance to the nearest
eligible fixation, which can be a different fixation. The original
`distance` column is retained as an alias for `nearest_distance`.
[`revisit_summary()`](https://bbuchsbaum.github.io/eyesim/reference/revisit_summary.md)
reduces the table to one row per viewing:

``` r

revisit_summary(rv)
#> # A tibble: 1 × 8
#>   n_fix opportunities revisits revisit_rate expected_revisits expected_rate
#>   <int>         <int>    <int>        <dbl>             <dbl>         <dbl>
#> 1     4             2        1          0.5                NA            NA
#> # ℹ 2 more variables: excess_revisits <dbl>, excess_rate <dbl>
```

Coordinates and onsets must be finite numeric values. Scoring rejects
invalid values instead of dropping fixations, because dropping a row
changes adjacency, fixation lags and denominators. Fixations are ordered
by onset, preserving input order among ties. Reference coordinates must
also be finite.

`radius` is in the units of your coordinates. To set it in degrees of
visual angle, convert with
[`visual_angle_to_px()`](https://bbuchsbaum.github.io/eyesim/reference/visual_angle_to_px.md)
for your display:

``` r

# 2 degrees on a 52 cm wide, 1920 px display viewed from 60 cm
visual_angle_to_px(2, distance = 60, screen_size = 52, screen_res = 1920)
#> [1] 77.33936
```

### Comparing revisits with a spatial null

Raw revisit counts partly restate how many fixations a viewer made and
how tightly they clustered. More fixations in a smaller area produce
more returns by chance. A spatial null estimates how many revisits to
expect from where people look. For each opportunity, it asks what
fraction of a set of reference locations would also have counted as a
revisit. With `require_departure = TRUE`, both the numerator and
denominator use only reference locations at least `radius` from the
preceding fixation.

The `eye_table` method scores every trial and builds the reference for
you. There are two choices:

- `null = "item"` uses fixations from other trials showing the same
  image (`item_var`). With `subject_var` given, all of the viewer’s own
  trials are left out. This keeps where the image draws gaze, but not
  the viewer’s own bias.
- `null = "subject"` uses the same participant’s other trials. This
  keeps the viewer’s spatial bias, but not the image content.

Score the encoding trials of the simulated experiment against the item
null:

``` r

enc_eyetab <- eyetab %>% filter(phase == "encoding")

revisits <- fixation_revisits(enc_eyetab, radius = 15,
                              null = "item", item_var = "image",
                              subject_var = "participant")

tibble::as_tibble(revisits) %>%
  select(participant, image, n_fix, opportunities, revisits,
         expected_revisits, excess_revisits) %>%
  head()
#> # A tibble: 6 × 7
#>   participant image n_fix opportunities revisits expected_revisits
#>   <chr>       <chr> <int>         <int>    <int>             <dbl>
#> 1 s1          img1     10             7        3             0.638
#> 2 s1          img10     4             2        0             0    
#> 3 s1          img11     3             0        0            NA    
#> 4 s1          img12    10             7        1             0.554
#> 5 s1          img13    10             8        0             4.30 
#> 6 s1          img14     4             2        0             0.160
#> # ℹ 1 more variable: excess_revisits <dbl>
```

Filter to a single phase before you score. Otherwise, a subject null
would include the same image from the other phase in the reference.

`excess_revisits` is the observed count minus the expected count, and
`excess_rate` is the same difference per opportunity. Trials with no
opportunities get `NA`, as do trials whose image no other participant
saw. An episode’s expectation also remains `NA` if even one opportunity
has no reference location satisfying the incoming-step condition.

These differences are descriptive corrections under the selected
reference distribution. They are not significance tests or a complete
generative model of scanpaths. Inference would need to account for
repeated observations and dependence introduced by shared reference
pools; a one-sample test treating all trial corrections as independent
does not do that.

![Observed against expected revisits for each encoding trial. Points
scatter around the identity line because the simulated fixations are
random.](eyesim_files/figure-html/plot-revisits-1.png)

Observed against expected revisits for each encoding trial. Points
scatter around the identity line because the simulated fixations are
random.

Set `keep_fixations = TRUE` to keep each trial’s fixation-level table in
a `revisit_fixations` column. Use `max_lag` to limit how far back a
return can reach. The `revisit_spec` attribute retains the radius, lag
bounds, distance metric, incoming-step rule, null choice and column
names even when only summaries are kept. Reference coverage is updated
as fixations enter and leave the lag window, without storing full
pairwise distance matrices. See
[`?fixation_revisits`](https://bbuchsbaum.github.io/eyesim/reference/fixation_revisits.md)
for the full definition and references.

## What’s next?

- [`vignette("Multimatch")`](https://bbuchsbaum.github.io/eyesim/articles/Multimatch.md)
  — compare scanpath structure across five dimensions
- [`vignette("RepetitiveSimilarity")`](https://bbuchsbaum.github.io/eyesim/articles/RepetitiveSimilarity.md)
  — analyze repeated viewing of the same stimulus
- [`vignette("latent-transforms")`](https://bbuchsbaum.github.io/eyesim/articles/latent-transforms.md)
  — domain adaptation for cross-device or cross-participant comparisons
- The package ships real data: `data(wynn_study)` and `data(wynn_test)`
  from a recognition memory experiment
