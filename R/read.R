## Read in some data 
## GWAS
## https://www.ebi.ac.uk/gwas/
library(here)
source(here("R","prepare.R"))
install.load("data.table")
## A list of ll associations
a <- fread(dropbox("gwas_catalog_v1.0.1-associations_e90_r2017-12-04.tsv"))
a2 <-  fread(dropbox("gwas-catalog-associations_ontology-annotated.tsv"))
#a <- a2
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
##
## I believe that if multiple SNPs are reported this is because 
## th esuthors are reporting a haplotype.

## Also if the studies do not seem to have been completely entered in the database then 
# bits may be missing.  For example a recent paper on Breast Carcinoma just has rs numbers, but no chromosome
## ID or position.  About 1813 of the bad SNPs are like this.


## Ancestry
ancestry <- fread(dropbox("gwas_catalog-ancestry_r2017-12-11.tsv"))
ancestry$V13 <- NULL
colnames(ancestry) <- c( "STUDY_ACCCESSION",
                         "PUBMEDID",
                         "FIRST_AUTHOR",
                         "DATE",
                         "INITIAL_SAMPLE_DESCRIPTION",
                         "REPLICATION_SAMPLE_DESCRIPTION",
                         "STAGE",
                          "NUMBER_OF_INDIVDUALS", 
                         "BROAD_ANCESTRAL_CATEGORY",
                         "COUNTRY_OF_ORIGIN", 
                         "COUNTRY_OF_RECRUITMENT",
                         "ADDITONAL ANCESTRY DESCRIPTION")
hist(log10(ancestry$`NUMBER OF INDIVDUALS`))
sort(table(ancestry$BROAD_ANCESTRAL_CATEGORY), decreasing = TRUE)[1:20]
