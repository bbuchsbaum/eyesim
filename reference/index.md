# Package index

## Density and sampling

Functions for computing density maps, sampling, and density matrix.

- [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  : eye_density
- [`get_density()`](https://bbuchsbaum.github.io/eyesim/reference/get_density.md)
  : get_density
- [`sample_density()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density.md)
  : sample_density
- [`density_matrix()`](https://bbuchsbaum.github.io/eyesim/reference/density_matrix.md)
  : Compute Density Matrix for a Given Object
- [`eye_density(`*`<fixation_group>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.fixation_group.md)
  : Compute a density map for a fixation group.
- [`gen_density()`](https://bbuchsbaum.github.io/eyesim/reference/gen_density.md)
  : This function creates a density object from the provided x, y, and z
  matrices. The density object is a list containing the x, y, and z
  values with a class attribute set to "density" and "list".
- [`as.data.frame(`*`<eye_density>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/as.data.frame.eye_density.md)
  : Convert an eye_density object to a data.frame.

## Spatial operations

Functions for manipulating spatial coordinates and transformations.

- [`coords()`](https://bbuchsbaum.github.io/eyesim/reference/coords.md)
  : extract coordinates
- [`rescale()`](https://bbuchsbaum.github.io/eyesim/reference/rescale.md)
  : rescale
- [`center()`](https://bbuchsbaum.github.io/eyesim/reference/center.md)
  : Center Eye-Movements in a New Coordinate System
- [`normalize()`](https://bbuchsbaum.github.io/eyesim/reference/normalize.md)
  : Normalize Eye-Movements to Unit Range
- [`match_scale()`](https://bbuchsbaum.github.io/eyesim/reference/match_scale.md)
  : Match Scaling Parameters for Fixation Data

## Fixation operations

Functions for working with fixation sequences.

- [`fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_group.md)
  : Create a Fixation Group Object
- [`fixation_entropy()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_entropy.md)
  : Compute Fixation Entropy
- [`fixation_entropy(`*`<default>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/fixation_entropy.default.md)
  [`fixation_entropy(`*`<eye_density>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/fixation_entropy.default.md)
  [`fixation_entropy(`*`<density>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/fixation_entropy.default.md)
  [`fixation_entropy(`*`<eye_density_multiscale>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/fixation_entropy.default.md)
  [`fixation_entropy(`*`<fixation_group>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/fixation_entropy.default.md)
  : Entropy of fixation patterns
- [`c(`*`<fixation_group>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/c.fixation_group.md)
  : Concatenate Fixation Groups
- [`rep_fixations()`](https://bbuchsbaum.github.io/eyesim/reference/rep_fixations.md)
  : rep_fixations
- [`sample_fixations()`](https://bbuchsbaum.github.io/eyesim/reference/sample_fixations.md)
  : sample_fixations
- [`fixation_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_similarity.md)
  : Fixation Similarity
- [`template_sample()`](https://bbuchsbaum.github.io/eyesim/reference/template_sample.md)
  : Sample density maps with coordinates derived from fixation groups.

## Similarity

Functions for computing similarity between objects.

- [`similarity()`](https://bbuchsbaum.github.io/eyesim/reference/similarity.md)
  : Compute Similarity Between Two Objects
- [`multi_match()`](https://bbuchsbaum.github.io/eyesim/reference/multi_match.md)
  : Compute MultiMatch Metrics for Scanpath Similarity
- [`scanpath_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/scanpath_similarity.md)
  : Scanpath Similarity
- [`similarity(`*`<scanpath>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/similarity.scanpath.md)
  : Compute Similarity Between Scanpaths
- [`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md)
  : template_similarity
- [`template_similarity_cv()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity_cv.md)
  : Cross-Fitted Template Similarity
- [`sample_density_time()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density_time.md)
  : Sample density maps at fixation locations over time
- [`repetitive_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/repetitive_similarity.md)
  : Repetitive Similarity Analysis for Density Maps
- [`install_multimatch()`](https://bbuchsbaum.github.io/eyesim/reference/install_multimatch.md)
  : Install Python multimatch_gaze Package

## GazeWeave

Cross-fitted probabilistic replay and registered temporal transport.

- [`gaze_weave_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_weave_cv.md)
  : Cross-fitted GazeWeave analysis
- [`gaze_replay_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_replay_cv.md)
  : Cross-fitted directional GazeWeave Replay analysis
- [`gaze_replay_spec()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_replay_spec.md)
  : Specify the directional GazeWeave replay model
- [`fit_gaze_replay_model()`](https://bbuchsbaum.github.io/eyesim/reference/fit_gaze_replay_model.md)
  : Fit a directional GazeWeave replay model on training trials
- [`gaze_replay_align()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_replay_align.md)
  : Align recall gaze to an encoding path with a fitted Replay model
- [`gaze_replay_align_episode()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_replay_align_episode.md)
  : Align recall gaze to a mixture of encoding presentations
- [`gaze_calibration_control()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_calibration_control.md)
  : Calibration controls for GazeWeave revision 2026.10
- [`gaze_transport_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_cv.md)
  : Nested-calibrated exhaustive-candidate Transport
- [`gaze_transport_spec()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_spec.md)
  : Specify edge-normalized GazeWeave Transport
- [`gaze_transport_prepare()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_prepare.md)
  : Prepare one immutable gaze path for repeated Transport scoring
- [`gaze_transport_align()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_align.md)
  : Align gaze paths with edge-normalized Transport
- [`gaze_transport_align_batch()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_align_batch.md)
  : Batch candidate alignments with Transport
- [`gaze_transport_result()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_result.md)
  : Extract one auditable Transport result
- [`gaze_gaussian_mixture()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_gaussian_mixture.md)
  : Declare a multiscale Gaussian spatial model for GazeWeave
- [`gaze_local_order()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_local_order.md)
  : Declare the local chronology model for GazeWeave
- [`gaze_order_neighbours()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_order_neighbours.md)
  : Declare fixed-neighbour ordinal chronology for Transport
- [`gaze_warp_none()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_warp_none.md)
  : Declare no geometric registration for GazeWeave
- [`gaze_warp_contraction()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_warp_contraction.md)
  : Declare cross-fitted contraction registration for GazeWeave
- [`gaze_screen()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_screen.md)
  : Declare screen geometry for GazeWeave
- [`tidy(`*`<gaze_replay_fit>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/tidy.gaze_replay_fit.md)
  : Tidy GazeWeave Replay results
- [`tidy(`*`<gaze_transport_fit>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/tidy.gaze_transport_fit.md)
  : Tidy a fitted Transport analysis
- [`autoplot(`*`<gaze_replay_alignment>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/autoplot.gaze_replay_alignment.md)
  : Plot a GazeWeave Replay alignment
- [`autoplot(`*`<gaze_replay_fit>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/autoplot.gaze_replay_fit.md)
  : Plot a fitted GazeWeave Replay analysis
- [`autoplot(`*`<gaze_transport_alignment>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/autoplot.gaze_transport_alignment.md)
  : Plot a Transport pair alignment
- [`autoplot(`*`<gaze_transport_fit>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/autoplot.gaze_transport_fit.md)
  : Plot one fitted Transport result

## GazeWeave comparators

Fair comparator court for GazeWeave (nested cross-validated baselines).

- [`gaze_baseline_spec()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_baseline_spec.md)
  : Specify the fair GazeWeave comparator court
- [`gaze_baseline_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_baseline_cv.md)
  : Nested cross-validated comparator court for GazeWeave
- [`gaze_baseline_availability()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_baseline_availability.md)
  : Report comparator availability
- [`tidy(`*`<gaze_baseline_fit>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/tidy.gaze_baseline_fit.md)
  : Tidy GazeWeave comparator results
- [`elastic_consensus_align()`](https://bbuchsbaum.github.io/eyesim/reference/elastic_consensus_align.md)
  : Elastic consensus matching for encoding and recall fixations

## Latent transforms

Domain-adaptation transforms for use with template_similarity().

- [`latent_pca_transform()`](https://bbuchsbaum.github.io/eyesim/reference/latent_pca_transform.md)
  [`contract_transform()`](https://bbuchsbaum.github.io/eyesim/reference/latent_pca_transform.md)
  [`affine_transform()`](https://bbuchsbaum.github.io/eyesim/reference/latent_pca_transform.md)
  [`coral_transform()`](https://bbuchsbaum.github.io/eyesim/reference/latent_pca_transform.md)
  [`cca_transform()`](https://bbuchsbaum.github.io/eyesim/reference/latent_pca_transform.md)
  : Latent-space transforms for template-based similarity

## Scanpath

Functions for creating and adding scanpaths.

- [`scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/scanpath.md)
  : Construct a Scanpath of a Fixation Group of Related Objects
- [`add_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/add_scanpath.md)
  : Add Scanpath to Dataset
- [`add_scanpath(`*`<data.frame>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/add_scanpath.data.frame.md)
  : Add Scanpath to a Data Frame
- [`add_scanpath(`*`<eye_table>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/add_scanpath.eye_table.md)
  : Add Scanpath to an Eye Table
- [`scanpath(`*`<fixation_group>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/scanpath.fixation_group.md)
  : Create a Scanpath for a Fixation Group
- [`calcangle()`](https://bbuchsbaum.github.io/eyesim/reference/calcangle.md)
  : Calculate the Angle Between Two Vectors
- [`cart2pol()`](https://bbuchsbaum.github.io/eyesim/reference/cart2pol.md)
  : Convert Cartesian Coordinates to Polar Coordinates
- [`template_multireg()`](https://bbuchsbaum.github.io/eyesim/reference/template_multireg.md)
  : Template Multiple Regression
- [`template_regression()`](https://bbuchsbaum.github.io/eyesim/reference/template_regression.md)
  : Template Regression

## Eye Table

Functions for working with eye tables.

- [`as_eye_table()`](https://bbuchsbaum.github.io/eyesim/reference/as_eye_table.md)
  : Reapply the 'eye_table' Class to an Object
- [`eye_table()`](https://bbuchsbaum.github.io/eyesim/reference/eye_table.md)
  : Construct an Eye-Movement Data Frame
- [`simulate_eye_table()`](https://bbuchsbaum.github.io/eyesim/reference/simulate_eye_table.md)
  : Generate a Simulated Eye-Movement Data Frame
- [`` `[`( ``*`<eye_table>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/sub-.eye_table.md)
  : Subset an 'eye_table' Object

## Plotting

Scanpath and density plots, and the shared eyesim plot theme.

- [`plot(`*`<fixation_group>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md)
  : Plot a fixation_group object
- [`plot(`*`<eye_density>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md)
  : Plot Eye Density
- [`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md)
  : Animate a Fixation Scanpath with gganimate
- [`theme_eyesim()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md)
  [`theme_eyesim_spatial()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md)
  : eyesim ggplot2 theme
- [`eyesim_colours()`](https://bbuchsbaum.github.io/eyesim/reference/eyesim_colours.md)
  : eyesim plot colours
- [`scale_fill_eyesim_density()`](https://bbuchsbaum.github.io/eyesim/reference/scale_fill_eyesim_density.md)
  [`scale_colour_eyesim_time()`](https://bbuchsbaum.github.io/eyesim/reference/scale_fill_eyesim_density.md)
  : eyesim colour scales
- [`element_text_wrap()`](https://bbuchsbaum.github.io/eyesim/reference/element_text_wrap.md)
  : Wrapping text element

## Density by Groups

Functions for density operations by groups.

- [`density_by()`](https://bbuchsbaum.github.io/eyesim/reference/density_by.md)
  : Calculate Eye Density by Groups
- [`suggest_sigma()`](https://bbuchsbaum.github.io/eyesim/reference/suggest_sigma.md)
  : Suggest Kernel Bandwidth for Density Estimation
- [`sample_density()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density.md)
  : sample_density

## Data

Example datasets included with the package.

- [`wynn_study`](https://bbuchsbaum.github.io/eyesim/reference/wynn_study.md)
  : Eye-tracking study data from Wynn et al.
- [`wynn_study_image`](https://bbuchsbaum.github.io/eyesim/reference/wynn_study_image.md)
  : Study image from Wynn et al.
- [`wynn_test`](https://bbuchsbaum.github.io/eyesim/reference/wynn_test.md)
  : Eye-tracking test data from Wynn et al.
