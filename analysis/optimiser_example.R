
# remake this without an assert
single_theoretical_coi <- function(coi, plaf = seq(0, 0.5, l = 101), method = "1") {
  # Check inputs
  assert_vector(plaf)
  assert_bounded(plaf, left = 0, right = 0.5)
  assert_increasing(plaf)
  assert_single_string(method)
  assert_in(method, c("1", "2"))

  # Determine the curve
  if (method == "1"){
    curve <- 1 - plaf^coi - (1 - plaf)^coi
  } else if (method == "2"){
    curve <- (plaf - plaf^coi)/(1 - plaf^coi - (1 - plaf)^coi)
  }
  return(curve)
}

# Make some data
L <- 1e4
p <- stats::rbeta(L, 1, 5)
p[p > 0.5] <- 1 - p[p > 0.5]
sim <- sim_biallelic(COI = 3, PLAF = p, epsilon = 0)
processed_data <- coiaf:::process_simulated_coi(sim, seq_error = 0)

# Function to calculate our likelihood (i.e. the distance)
ll_function <- function(coi, processed_data) {

theory_coi <- single_theoretical_coi(coi, processed_data$midpoints, method = "1")
gap <- theory_coi - processed_data$m_variant
gap <- (gap * processed_data$bucket_size) / sum(processed_data$bucket_size)

# pick a distnace metric - am just doing sum of absolute error
gap <- sum(abs(gap))
return(gap)
}

# this is our optimiser. We pass in our parameters starting value (coi = 2 here)
# then our function that returns our value to be optimised
# then any other arguuments needed by ll_function (here the processed data)
# then the method - this is just so that our COI is bounded (might not be needed but good to be sure)
# then the lower and upper bounds
# then we have a control list for the optimiser
fit <- optim(par = 2,
             fn = ll_function,
             processed_data = processed_data,
             method = "L-BFGS-B",
             lower = 1, upper = 25,
             control = list(
               fnscale = 1, # this says this is a minimisation problem (i.e. reducing the value returned by our ll_function)
               ndeps = 1e-5 # this is something that controls out step sizes being taken in the optimiser.
             ))

# we get an object that returns us our estiamte of the parameter (coi)
coi <- fit$par
coi

# also other optimiser checks - convergence should be 0 if it has converged happily
fit$convergence

# for example if we remove our step size in ndeps (default is 1e-3 fyi)
fit <- optim(par = 2,
             fn = ll_function,
             processed_data = processed_data,
             method = "L-BFGS-B",
             lower = 1, upper = 25,
             control = list(
               fnscale = 1
             ))

# unlikely to converge
fit$convergence

# check out ?optim for more of this documentation
