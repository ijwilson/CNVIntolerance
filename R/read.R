## Read the GWAS
library(here)
source(here("R","prepare.R"))
install.load("data.table")
## A list of ll associations
a <- fread(dropbox("gwas_catalog_v1.0.1-associations_e90_r2017-12-04.tsv"))
## there are no cnv studies
a$CNV <- NULL
a$pos <- as.integer(a$CHR_POS)
a$chr <- a$CHR_ID
a$chr[a$CHR_ID=="X"] <- 23
a$chr <- as.integer(a$chr)
table(a$chr)
colnames(a)
table(a$CHR_ID)
hist(a$CHR_POS)

cat("We have ", length(table(a$`DISEASE/TRAIT`)), "Different Diseases/Traits\n")
cat("From ",length(table(a$PUBMEDID)), "unique pubmed IDs\n")
#######
tb <- table(a$SNPS)
cat(length(table(a$SNPS)), "diffent significant snps with", sum(tb>1), "snps seen in more than one study and one snp seen in ",max(table(a$SNPS)),"studies\n")

##############################  Simple SNPs #######################
simple <- a$CHR_ID %in% c(1:23,"X")
simple2 <- !is.na(as.integer(a$CHR_POS))


a[!simple]

table(simple, simple2)

## OK lots of the SNPs reported are reported as multiple SNPs.
## Read the paper/information to see what to do about this

