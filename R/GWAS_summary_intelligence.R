## Genome-wide association meta-analysis of 78,308 individuals identifies new loci and genes influencing
## human intelligence
## Sniekers et al Nature Genetics volume 49, pages 1107–1112 (2017) doi:10.1038/ng.3869
#' https://www.nature.com/articles/ng.3869
#' 

library(here)
source(here("R","prepare.R"))

install.load("data.table", "qqman")

## neuroticism summary results from my own computer
a <- fread("D:/data/GWAS_Summary/Intelligence/Sniekers.txt")

head(a)

png(height=500, width=1200, file=here("figs","intelligence.png"))

manhattan(a, bp="position", chr="Chromosome", p="p_value", snp="rsid")

dev.off()

# can i get a merge?


