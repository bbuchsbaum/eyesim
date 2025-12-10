# Sample density maps with coordinates derived from fixation groups.

this function extracts the density for any arbitrary time point in a
trial. It simply extracts the value of the density map for the fixation
at time t. The fixations are taken from the \`fixgroup\` variable and
may be associated, for example, with an independent set of trials.

## Usage

``` r
template_sample(
  source_tab,
  template,
  fixgroup = "fixgroup",
  time = NULL,
  outcol = "sample_out"
)
```

## Arguments

- source_tab:

  the name of the table containing the density map and fixations

- template:

  the name of the template density variable

- fixgroup:

  the name of the fixation group supplying the spatiotemporal
  coordinates used to sample the template

- time:

  the time points used to extract coordinates from the
  \`fixation_group\`

- outcol:

  the name of the output variable
