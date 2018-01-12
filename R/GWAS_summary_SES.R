library(here)
source(here("R","prepare.R"))

install.load("data.table", "qqman")

## income summary results from my own 
a <- fread("file:///D:/data/GWAS_Summary/Hill_CB_2016/Hill2016_UKB_Income_summary_results_21112016.txt")

head(a)

png(height=500, width=1200, file="SES.png")

manhattan(a, bp="Position", chr="Chromosome", p="P-value", snp="Markername")

dev.off()


#Nothing really that produces genome wide sifnificance.