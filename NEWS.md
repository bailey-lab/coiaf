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

- New `theme_coiaf()` creates a custom theme for this package.

## Maintenance

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
