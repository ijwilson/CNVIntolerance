## preparation
library(here)
source(here("R","prepare.R"))
## Load file to give locations for rs numbers
## https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP144.GRCh37.html
install.load.bioc("SNPlocs.Hsapiens.dbSNP144.GRCh37")
## load height GWAS summary statistics
GWAS_height_summary <- read.table(
  dropbox("GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"), 
  header=TRUE, stringsAsFactors = FALSE)
head(GWAS_height_summary)

snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
#gr <- rsidsToGRanges(snps, GWAS_height_summary$MarkerName[1:20])

t100    <- system.time(gr <- snpsById(snps, GWAS_height_summary$MarkerName[1:100]))
t1000   <- system.time(gr <- snpsById(snps, GWAS_height_summary$MarkerName[1:1000]))
gr <- snpsById(snps, GWAS_height_summary$MarkerName, ifnotfound="warning")

GWAS_height_summary <- GWAS_height_summary[GWAS_height_summary$MarkerName %in% gr$RefSNP_id,]
table(GWAS_height_summary$MarkerName== gr$RefSNP_id)
GWAS_height_summary <- cbind(GWAS_height_summary, gr)
