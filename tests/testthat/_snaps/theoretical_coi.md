# coi_method argument is restricted

    Code
      theoretical_coi(1:5, seq(0, 0.5, 0.1))
    Output
      # A tibble: 6 x 6
        coi_1 coi_2 coi_3 coi_4 coi_5 plmaf
        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
      1     0  0     0    0     0       0  
      2     0  0.18  0.27 0.344 0.410   0.1
      3     0  0.32  0.48 0.589 0.672   0.2
      4     0  0.42  0.63 0.752 0.830   0.3
      5     0  0.48  0.72 0.845 0.912   0.4
      6     0  0.5   0.75 0.875 0.938   0.5

---

    Code
      theoretical_coi(1:5, seq(0, 0.5, 0.1), coi_method = "frequency")
    Output
      # A tibble: 6 x 6
        coi_1 coi_2   coi_3   coi_4   coi_5 plmaf
        <dbl> <dbl>   <dbl>   <dbl>   <dbl> <dbl>
      1   NaN NaN   NaN     NaN     NaN       0  
      2   NaN   0.5   0.367   0.291   0.244   0.1
      3   NaN   0.5   0.4     0.337   0.297   0.2
      4   NaN   0.5   0.433   0.388   0.359   0.3
      5   NaN   0.5   0.467   0.443   0.427   0.4
      6   NaN   0.5   0.5     0.5     0.5     0.5

---

    `coi_method` must be one of "variant" or "frequency", not "wrong method".

