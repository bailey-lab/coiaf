# cut midpoints work for a large data set (#18)

    Code
      compute_coi(tibble::as_tibble(data), "real", seq_error = data$seq_error,
      bin_size = data$bin_size, coi_method = data$coi_method, use_bins = TRUE)
    Output
      $coi
      [1] 1
      
      $probability
       [1] 0.7589392288 0.0970897028 0.0434247877 0.0246155376 0.0159032233
       [6] 0.0111576164 0.0082839487 0.0064095253 0.0051176318 0.0041885658
      [11] 0.0034974474 0.0029689706 0.0025555018 0.0022256993 0.0019582463
      [16] 0.0017382291 0.0015549565 0.0014005993 0.0012693162 0.0011566764
      [21] 0.0010592684 0.0009744307 0.0009000617 0.0008344838 0.0007763443
      
