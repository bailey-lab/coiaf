param_grid <- expand.grid(COI=2:6, epsilon = c(0.0,0.1))

repetitons <- 10

out <- lapply(seq_len(nrow(param_grid)), function(x) {

  out <- vapply(seq_len(repetitons), function(y) {

  sim1 <- sim_biallelic(COI = param_grid$COI[x], epsilon = param_grid$epsilon[x], PLAF = p)
  p_vec <- seq(0, 0.5, l=101)
  theory_cois <- theoretical_coi(seq(1, 7), p_vec)
  cut <- seq(0,0.5,0.01)
  df_grouped <- simulated_coi(sim1, seq_error = 0.01, cut, theory_cois)

  calc_coi <- compute_coi(seq(1, 7), df_grouped, cut, "overall", "squared")

  return(calc_coi)
  }, FUN.VALUE = character(1))

  return(out)

})

