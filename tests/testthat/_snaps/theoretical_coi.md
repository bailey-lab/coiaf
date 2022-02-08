# coi_method argument is restricted

    Code
      theoretical_coi(1:5, seq(0, 0.5, 0.1))
    Output
        coi_1 coi_2 coi_3  coi_4  coi_5 plaf
      1     0  0.00  0.00 0.0000 0.0000  0.0
      2     0  0.18  0.27 0.3438 0.4095  0.1
      3     0  0.32  0.48 0.5888 0.6720  0.2
      4     0  0.42  0.63 0.7518 0.8295  0.3
      5     0  0.48  0.72 0.8448 0.9120  0.4
      6     0  0.50  0.75 0.8750 0.9375  0.5

---

    Code
      theoretical_coi(1:5, seq(0, 0.5, 0.1), coi_method = "frequency")
    Output
        coi_1 coi_2     coi_3     coi_4     coi_5 plaf
      1   NaN   NaN       NaN       NaN       NaN  0.0
      2   NaN   0.5 0.3666667 0.2905759 0.2441758  0.1
      3   NaN   0.5 0.4000000 0.3369565 0.2971429  0.2
      4   NaN   0.5 0.4333333 0.3882682 0.3587342  0.3
      5   NaN   0.5 0.4666667 0.4431818 0.4273684  0.4
      6   NaN   0.5 0.5000000 0.5000000 0.5000000  0.5

---

    `coi_method` must be one of "variant" or "frequency", not "wrong method".

