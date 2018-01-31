## preparation
library(here)
source(here("R","prepare.R"))
## Load file to give locations for rs numbers
## https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP144.GRCh37.html
## load scz GWAS summary statistics
# 
# scz2.snp.results.txt.gz (header row and 9,444,230 SNPs)
# hg19chrc hg19 chromosome as character string (chr1-chr22, chrX)
# snpid rs ID of SNP
# a1 reference allele for OR (may not be minor allele)
# a2 alternate allele
# bp hg19 base pair position of SNP
# info imputation quality score
# or odds ratio in PGC GWAS data
# se standard error of ln(OR) in PGC GWAS data
# p p-value in PGC GWAS data
# ngt number of samples in which SNP directly genotyped



if (!dir.exists(here("output"))) {
  dir.create(here("output"))
}

if (file.exists(here("output", "GWAS_scz_summary.rda"))) {
  load(here("output", "GWAS_scz_summary.rda"))
} else {

  GWAS_scz_summary <- read.table(
    dropbox("scz2/ckqny.scz2snpres.gz"), 
    header=TRUE, stringsAsFactors = FALSE)
  head(GWAS_scz_summary)

  cat("There are", nrow(GWAS_scz_summary), "rows in the summary statistics.\n")

  GWAS_scz_summary$a1 <- NULL
  GWAS_scz_summary$a2 <- NULL
  GWAS_scz_summary$chr <- substring(GWAS_scz_summary$hg19chrc,4) 
  GWAS_scz_summary$chr[GWAS_scz_summary$chr=="X"] <- 22
  GWAS_scz_summary$chr <- as.numeric(GWAS_scz_summary$chr)
  GWAS_scz_summary$hg19chrc <- NULL
  save(GWAS_scz_summary, file=here("output", "GWAS_scz_summary.rda"))
}

if (FALSE) {   ## takes too long
  plot(density(GWAS_scz_summary$or))
  
  install.load("qqman")
  
  png(height=500, width=1200, file="scz.png")
  manhattan(GWAS_scz_summary, chr="chr", bp="bp", snp="snpid", p="p" )
  dev.off()
  
  qq(GWAS_scz_summary)
  
}



