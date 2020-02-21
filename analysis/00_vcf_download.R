ftp_address="ftp://ngs.sanger.ac.uk/production/pf3k/release_5/5.1/"

## Grab the filenames from the ftp address
filenames <- RCurl::getURL(ftp_address, ftp.use.epsv = FALSE, dirlistonly = TRUE) %>%
  strsplit("\n") %>%
  unlist()

## chromosome
chromosome <- 1

## look up the correct file for the chromosome specified
if(chromosome < 10){
  pattern <- paste0("_0",chromosome,".*gz$")
  chrom <- paste0("Pf3D7_0",chromosome,"_v3")
} else {
  pattern <- paste0("_",chromosome,".*gz$")
  chrom <- paste0("Pf3D7_",chromosome,"_v3")
}

## create full file path
dir.create("data-raw")
vcffile <- paste0(ftp_address,filenames[grep(pattern,filenames)])
tbi_vcffile <- paste0(vcffile, ".tbi")
download.file(vcffile, file.path("data-raw", basename(vcffile)))
download.file(tbi_vcffile, file.path("data-raw", basename(tbi_vcffile)))

