library(here)
source(here("R","prepare.R"))

install.load("data.table", "qqman")

## neuroticism summary results from my own computer
a <- fread("file:///D:/data/GWAS_Summary/Luciano_2017/SummaryStats.txt")

head(a)

png(height=500, width=1200, file="neuroticism.png")

manhattan(a, bp="BP", chr="CHR", p="p_value", snp="rsid")

dev.off()


## Lots with genome wide significance

