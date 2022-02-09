# Set the path and then access the data
path <- "~/Desktop/Malaria/COI data/"
wsmafs <- readRDS(paste0(path, "RMCL_wsafs_unique.rds"))
example_real_data <- wsmafs[[1]]
example_real_data <- example_real_data[11:20, 1:1000]

usethis::use_data(example_real_data, overwrite = TRUE)
