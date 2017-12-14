## Height
library(data.table)
height <- fread("~/Dropbox/Intolerance/height_CNV_association_41467_2017_556_MOESM2_ESM.csv")

## Simple of Manhattan plot of BMI
plot(-log10(height$`Pvalue BMI`), col=height$CHR%%2+1)
## Simple plot of Height
plot(-log10(height$`Pvalue Height`), col=height$CHR%%2+1, main="Height")
