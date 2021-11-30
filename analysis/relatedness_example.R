# relatedness attempt:
set.seed(1)

# Define number of loci, and distribution of minor allele frequencies
L <- 1e4
p <- rbeta(L, 1, 5)
p[p > 0.5] <- 1 - p[p > 0.5]
coverage <- 100
reps <- 20

cois <- 2:10
rels <- seq(0, 0.9, 0.1)
pars <- expand.grid(k = cois, r = rels)
results <- vector("list", nrow(pars))

for (i in seq_len(nrow(pars))) {
  message(i)
  k <- pars$k[i]
  r <- pars$r[i]

  dat <- lapply(seq_len(reps), sim_biallelic, coi = k, plaf = p, coverage = coverage, relatedness = r)

  disc_k <- lapply(dat, compute_coi, data_type = "sim", coi_method = "variant")
  disc_k <- unlist(lapply(disc_k, "[[", "coi"))

  disc_k_2 <- lapply(dat, compute_coi, data_type = "sim", coi_method = "frequency")
  disc_k_2 <- unlist(lapply(disc_k_2, "[[", "coi"))

  cont_k <- lapply(dat, optimize_coi, data_type = "sim", coi_method = "variant")
  cont_k_2 <- lapply(dat, optimize_coi, data_type = "sim", coi_method = "frequency")

  r_est <- unlist(cont_k_2) - unlist(cont_k)

  df <- data.frame(
    "coi" = k, "r" = r,
    "rep" = seq_len(reps), "r_est" = r_est,
    "coi_d" = unlist(disc_k), "coi_d_2" = unlist(disc_k_2),
    "coi_c" = unlist(cont_k), "coi_c_2" = unlist(cont_k_2)
  )

  results[[i]] <- df
}

# plot
plot <- do.call(rbind, results) %>%
  dplyr::mutate(coi_lab = paste0("COI = ", coi))

plot$coi_lab <- factor(plot$coi_lab,
  levels = c(
    "COI = 2", "COI = 3", "COI = 4", "COI = 5",
    "COI = 6", "COI = 7", "COI = 8", "COI = 9",
    "COI = 10"
  )
)

ggplot(aes(x = r, y = r_est), data = plot) +
  geom_point(color = "blue", alpha = 0.7, show.legend = FALSE) +
  geom_smooth(method = "lm", color = "blue") +
  geom_abline(color = "red", lty = 2, size = 1) +
  theme_classic() +
  theme(axis.line = element_line()) +
  facet_wrap(~coi_lab) +
  ylab("Method 2 Continous - Method 1 Continuous") +
  xlab("Relatedness")
