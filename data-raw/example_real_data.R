# Set the path and then access the data
path = "~/Desktop/Malaria/COI data/"
rmcl_wsafs   <- readRDS(paste0(path, "RMCL_wsafs_unique.rds"))
example_real_data <- rmcl_wsafs[1:2]

usethis::use_data(example_real_data, overwrite = TRUE)
