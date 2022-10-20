# compute_coi() args deprecated

    Code
      compute_coi(data, "sim", comparison = "end", use_bins = TRUE)
    Condition
      Warning:
      The `comparison` argument of `compute_coi()` is deprecated as of coiaf 0.2.0.
      i The comparison method will be fixed to "overall" in the next release.
      Warning:
      The `use_bins` argument of `compute_coi()` is deprecated as of coiaf 0.2.0.
      i The ability to use bins to estimate the COI will be dropped in the next release.
    Output
      $coi
      [1] 6
      
      $probability
       [1] 0.00579567 0.01084526 0.01921675 0.03726912 0.10980705 0.25458635
       [7] 0.07322100 0.04771043 0.03780684 0.03270041 0.02968371 0.02775688
      [13] 0.02646422 0.02556819 0.02493295 0.02447536 0.02414194 0.02389698
      [19] 0.02371590 0.02358143 0.02348125 0.02340642 0.02335043 0.02330847
      [25] 0.02327699
      

---

    Code
      compute_coi(data, "sim", distance = "abs_sum", use_bins = TRUE)
    Condition
      Warning:
      The `distance` argument of `compute_coi()` is deprecated as of coiaf 0.2.0.
      i The distance method will be fixed to "squared" in the next release.
      Warning:
      The `use_bins` argument of `compute_coi()` is deprecated as of coiaf 0.2.0.
      i The ability to use bins to estimate the COI will be dropped in the next release.
    Output
      $coi
      [1] 6
      
      $probability
       [1] 0.005794631 0.010843425 0.019213808 0.037264714 0.109809533 0.254663688
       [7] 0.073217458 0.047705762 0.037802413 0.032696253 0.029679763 0.027753088
      [13] 0.026460533 0.025564588 0.024929402 0.024471855 0.024138475 0.023893535
      [19] 0.023712471 0.023578019 0.023477845 0.023403024 0.023347034 0.023305078
      [25] 0.023273606
      

---

    The `use_bins` argument of `compute_coi()` is deprecated as of coiaf 0.2.0.
    i The ability to use bins to estimate the COI will be dropped in the next release.

---

    Code
      compute_coi(data, "sim", bin_size = 100, use_bins = TRUE)
    Condition
      Warning:
      The `use_bins` argument of `compute_coi()` is deprecated as of coiaf 0.2.0.
      i The ability to use bins to estimate the COI will be dropped in the next release.
      Warning:
      The `bin_size` argument of `compute_coi()` is deprecated as of coiaf 0.2.0.
      i The ability to use bins to estimate the COI will be dropped in the next release.
    Output
      $coi
      [1] 6
      
      $probability
       [1] 0.0003420419 0.0011977369 0.0037605936 0.0141457645 0.1228295740
       [6] 0.6605069038 0.0546083217 0.0231831616 0.0145569335 0.0108899679
      [11] 0.0089732816 0.0078460832 0.0071322633 0.0066574468 0.0063307298
      [16] 0.0061004773 0.0059353954 0.0058155500 0.0057277436 0.0056629740
      [21] 0.0056149565 0.0055792249 0.0055525613 0.0055326228 0.0055176898
      

# optimize_coi() args deprecated

    Code
      optimize_coi(data, "sim", distance = "abs_sum", use_bins = TRUE)
    Condition
      Warning:
      The `distance` argument of `optimize_coi()` is deprecated as of coiaf 0.2.0.
      i The distance method will be fixed to "squared" in the next release.
      Warning:
      The `use_bins` argument of `optimize_coi()` is deprecated as of coiaf 0.2.0.
      i The ability to use bins to estimate the COI will be dropped in the next release.
    Output
      [1] 5.6664

---

    The `use_bins` argument of `optimize_coi()` is deprecated as of coiaf 0.2.0.
    i The ability to use bins to estimate the COI will be dropped in the next release.

---

    Code
      optimize_coi(data, "sim", bin_size = 100, use_bins = TRUE)
    Condition
      Warning:
      The `use_bins` argument of `optimize_coi()` is deprecated as of coiaf 0.2.0.
      i The ability to use bins to estimate the COI will be dropped in the next release.
      Warning:
      The `bin_size` argument of `optimize_coi()` is deprecated as of coiaf 0.2.0.
      i The ability to use bins to estimate the COI will be dropped in the next release.
    Output
      [1] 5.6664

# compute_coi_regression() args deprecated

    The `distance` argument of `compute_coi_regression()` is deprecated as of coiaf 0.2.0.
    i The distance method will be fixed to "squared" in the next release.

# optimize_coi_regression() args deprecated

    The `distance` argument of `optimize_coi_regression()` is deprecated as of coiaf 0.2.0.
    i The distance method will be fixed to "squared" in the next release.

# bootstrap_ci() args deprecated

    The `use_bins` argument of `bootstrap_ci()` is deprecated as of coiaf 0.2.0.
    i The ability to use bins to estimate the COI will be dropped in the next release.

---

    Code
      bootstrap_ci(data, "sim", bin_size = 100, use_bins = TRUE, replicates = 10)
    Condition
      Warning:
      The `use_bins` argument of `bootstrap_ci()` is deprecated as of coiaf 0.2.0.
      i The ability to use bins to estimate the COI will be dropped in the next release.
      Warning:
      The `bin_size` argument of `bootstrap_ci()` is deprecated as of coiaf 0.2.0.
      i The ability to use bins to estimate the COI will be dropped in the next release.
      Warning in `norm.inter()`:
      extreme order statistics used as endpoints
    Output
      # A tibble: 1 x 6
          coi estimates       bias std.error conf.low conf.high
        <dbl> <list>         <dbl>     <dbl>    <dbl>     <dbl>
      1    25 <dbl [10 x 1]> -14.9      8.24        4        25
