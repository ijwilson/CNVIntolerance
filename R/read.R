## Read the GWAS

library(data.table)
## A list of ll associations
a <- fread("C:/Users/nijw/Dropbox/Intolerance/gwas_catalog_v1.0.1-associations_e90_r2017-12-04.tsv")

colnames(a)
table(a$CHR_ID)
hist(a$CHR_POS)
a$pos <- as.integer(a$CHR_POS)
