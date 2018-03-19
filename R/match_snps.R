## match the SNPs from the two major experiments


library(here)
source(here("R","prepare.R"))

if (!dir.exists(here("output"))) {
  dir.create(here("output"))
}

library(data.table) 

GWAS_height_summary <- data.table(read.table(
  dropbox("GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"), 
  header=TRUE, stringsAsFactors = FALSE))
head(GWAS_height_summary)

## This Height dats is  
## Wood AR, Esko T, Yang J, Vedantam S, Pers TH, Gustafsson S et al. (2014). 
## Defining the role of common variation in the genomic and biological architecture 
## of adult human height. Nature Genetics 11:1173-86. 

## Earlier Height
#
# wget http://portals.broadinstitute.org/collaboration/giant/images/b/b0/GIANT_Yang2012Nature_publicrelease_HapMapCeuFreq_Height.txt.gz
#  wget http://portals.broadinstitute.org/collaboration/giant/images/5/5e/GIANT_Yang2012Nature_publicrelease_HapMapCeuFreq_BMI.txt.gz
#
#



if (FALSE) {
  GWAS_height_summaryB <- data.table(read.table(
    dropbox("GIANT_Yang2012Nature_publicrelease_HapMapCeuFreq_Height.txt.gz" ), 
    header=TRUE, stringsAsFactors = FALSE))
  head(GWAS_height_summaryB)

  GWAS_BMI_summary <- data.table(read.table(
    dropbox("GIANT_Yang2012Nature_publicrelease_HapMapCeuFreq_BMI.txt.gz" ), 
    header=TRUE, stringsAsFactors = FALSE))
  head(GWAS_BMI_summary)

}

GWAS_scz_summary <- data.table(
  read.table(dropbox("scz2/ckqny.scz2snpres.gz"), header=TRUE, stringsAsFactors = FALSE)
  )

head(GWAS_scz_summary)

colnames(GWAS_height_summary)[1] <- "snpid"
setkey(GWAS_height_summary, "snpid")
setkey(GWAS_scz_summary, "snpid")
snp_both <- merge(GWAS_scz_summary, GWAS_height_summary)

snp_both <- snp_both[snp_both$a1==snp_both$Allele1]
snp_both <- snp_both[snp_both$a2==snp_both$Allele2]
snp_both[,Allele1 := NULL]
snp_both[,Allele2 := NULL]
snp_both[,N := NULL]
snp_both[,ngt := NULL]
colnames(snp_both) <- c("snpid", "chrom", "a1", "a2", "pos", "info", "or", "or.se", "p.scz", "freq", "b", "b.se", "p.height")

save(snp_both, file=here("output", "snp_match.rda"))

