## preparation
library(here)
source(here("R","prepare.R"))
## Load file to give locations for rs numbers
## https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP144.GRCh37.html
## load height GWAS summary statistics

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

lines(debninstall.load("qqman")
manhattan(GWAS_height_summary, chr="seqnames", bp="pos", p="p", snp="MarkerName" )

qq(GWAS_height_summary)


