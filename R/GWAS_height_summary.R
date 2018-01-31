## preparation
library(here)
source(here("R","prepare.R"))
## Load file to give locations for rs numbers
## https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP144.GRCh37.html
## load height GWAS summary statistics

# 2012-2015 Data File Description:
#   Each file consists of the following information for each SNP and its association to the specified trait based on meta-analysis in the respective publication. Significant digits for the p values, betas and standard errors are limited to two digits to further limit the possibility of identifiability.
# 
# MarkerName: The dbSNP name of the genetic marker
# Allele1: The first allele (hg19 + strand). Where the regression coefficients (betas) are provided, the first allele is the effect allele. Where betas are not provided (typically the 2010 data), the first allele is the trait-increasing allele.
# Allele2: The second allele (hg19 + strand)
# Freq.Allele1.HapMapCEU: The allele frequency of Allele1 in the HapMap CEU population
# b: beta
# SE: standard error
# p: p-value after meta-analysis using regression coefficients (beta and standard error), and after correction for inflation of test statistics using genomic control both at the individual study level and again after meta-analysis
# N: Number of observations


if (!dir.exists(here("output"))) {
  dir.create(here("output"))
}

if (file.exists(here("output", "GWAS_height_summary.rda"))) {
  load(here("output", "GWAS_height_summary.rda"))
} else {

  GWAS_height_summary <- read.table(
    dropbox("GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"), 
    header=TRUE, stringsAsFactors = FALSE)
  head(GWAS_height_summary)

  cat("There are", nrow(GWAS_height_summary), "rows in the summary statistics.\n")

  GWAS_height_summary <- GWAS_height_summary[substring(GWAS_height_summary$MarkerName, 1, 2)=="rs", ]

  install.load.bioc("SNPlocs.Hsapiens.dbSNP144.GRCh37")

  snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
  gr <- snpsById(snps, GWAS_height_summary$MarkerName, ifnotfound="warning")

  GWAS_height_summary <- GWAS_height_summary[GWAS_height_summary$MarkerName %in% gr$RefSNP_id,]

  cat("after leaving only those with positions, there are", 
     nrow(GWAS_height_summary), "rows in the summary statistics.\n")
  GWAS_height_summary <- cbind(GWAS_height_summary, gr)


  GWAS_height_summary$strand <- NULL
  GWAS_height_summary$RefSNP_id <- NULL
  GWAS_height_summary$alleles_as_ambig <- NULL
  GWAS_height_summary <- GWAS_height_summary[!(GWAS_height_summary$seqnames %in% c("MT","X","Y")),]
  GWAS_height_summary$seqnames <- as.numeric(GWAS_height_summary$seqnames)

  save(GWAS_height_summary, file=here("output", "GWAS_height_summary.rda"))
}

#install.load("ggplot2")
#ggplot(GWAS_height_summary, aes(x=b)) + geom_density() + xlim(-0.05, 0.05)
plot(density(GWAS_height_summary$b))

install.load("qqman")
manhattan(GWAS_height_summary, chr="seqnames", bp="pos", p="p", snp="MarkerName" )

qq(GWAS_height_summary)


