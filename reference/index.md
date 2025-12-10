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
- [`sample_density_time()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density_time.md)
  : Sample density maps at fixation locations over time
- [`repetitive_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/repetitive_similarity.md)
  : Repetitive Similarity Analysis for Density Maps
- [`install_multimatch()`](https://bbuchsbaum.github.io/eyesim/reference/install_multimatch.md)
  : Install Python multimatch_gaze Package

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
- [`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md)
  : Animate a Fixation Scanpath with gganimate
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
- [`plot(`*`<eye_density>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md)
  : Plot Eye Density
- [`plot(`*`<fixation_group>`*`)`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md)
  : Plot a fixation_group object

## Density by Groups

Functions for density operations by groups.

- [`density_by()`](https://bbuchsbaum.github.io/eyesim/reference/density_by.md)
  : Calculate Eye Density by Groups
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
