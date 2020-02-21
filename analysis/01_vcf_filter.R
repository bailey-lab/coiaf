
bcf_filter <- function(x){

  if (x < 10) {
    x <- paste0("0", x)
  }

  s <- paste0("/users/owatson1/data/owatson1/vcfs/pf3k_release5/SNP_INDEL_Pf3D7_",x)
  end <- "_v3.combined.filtered.vcf.gz"
  s <- paste0(s,end)

  t <- paste0("/users/owatson1/data/owatson1/vcfs/pf3k_release5_filtered/SNP_INDEL_Pf3D7_",x)
  t <- paste0(t,"_filtered", end);

  cmd <- paste0("module load bcftools/1.9; bcftools view -i 'FILTER=\"PASS\"' -m 2 -M 2 -v snps -O z \"",s, "\" > \"", t, "\"")

  system(command = cmd)

}

pars <- data.frame("x" = 1:14)
sopt <- list(time = '1:00:00', share = TRUE)

sjob <- rslurm::slurm_apply(bcf_filter, pars, jobname = 'test_filter',
                            nodes = 14, cpus_per_node = 1, slurm_options = sopt,
                            submit = TRUE)
res <- rslurm::get_slurm_out(sjob)
rslurm::print_job_status(sjob)
