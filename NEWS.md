# coiaf (development version)

## New vignettes

- New `vignette("example_real_data")` illustrates how to use the package on your
  own real data sets.
- New `vignette("sensitivity_analysis_discrete")` describes the sensitivity
  analysis for the discrete methods.
- New `vignette("sensitivity_analysis_continuous")` describes the sensitivity
  analysis for the continuous methods.
- New `vignette("example_coi_prediction")` explains the internal prediction
  methods.
- New `vignette("example_simulation")` explains how data is simulated.

## New features

- An `autoplot` and `plot` method have been written for simulated data (#22).
- Simulated data have been restructured and are now assigned a class `sim`
  (#22).
- The default sequence error threshold has been set to 1%.
- Our methods now can estimate the COI without grouping data points into
  buckets (#17).
- Within sample minor allele frequencies are now weighed by read coverage. If no
  coverage is supplied, the coverage is assumed to be uniform across all loci
  (#16.)
- New `theme_coiaf()` creates a custom theme for this package.

## Bug Fixes

- Midpoint calculation for buckets has been improved. Previously, if buckets
  consisted of a single value, we would combine the bucket with another bucket.
  We no longer do so and let buckets consist of only a single value (#19).
- In our estimation methods, each bucket is now weighed by the number of points
  in each bucket (#15).
- The Frequency Method is undefined for a COI less than two. When we suspect the
  COI should be less than two, we return `NaN`. We additionally return a note
  indicating that the COI may be one. If we are able to compute the COI, we also
  return the COI estimated by our methods (#14, #21, #23).

## Maintenance

- PLAF has been renamed to PLMAF to indicate the focus is the minor allele.
- WSAF has been renamed to WSMAF to indicate the focus is the minor allele.
- Fix partial argument match warnings.
- Update license year.
- The number of dependencies has been reduced.
- `{patchwork}` is used instead of `{ggpubr}` for combining plots.
- `{rlang}` is now used for messages, warnings, and errors.
- Superseded functions have been replaced with their alternatives (#6).
- Functions now exit implicitly and visibly (#8).
- Internal documentation has been improved.
- COI methods have been renamed for improved clarity. `"Method 1"` has been
  renamed to `"Variant Method"` and `"Method 2"` has been renamed to
  `"Frequency Method"` (#9).

# coiaf 0.1.0

- Initial public release.
