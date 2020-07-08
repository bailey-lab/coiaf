# Set the path and then access the data
path = "~/Desktop/Malaria/COI data/"
rmcl_wsafs   <- readRDS(paste0(path, "RMCL_wsafs_unique.rds"))
example_real_data <- rmcl_wsafs[1:2]
example_real_data$cat_region_1_vcf_0 <- example_real_data$cat_region_1_vcf_0[11:20, 1:1000]
example_real_data$cat_region_1_vcf_1 <- example_real_data$cat_region_1_vcf_1[11:20, 1:1000]

usethis::use_data(example_real_data, overwrite = TRUE)
